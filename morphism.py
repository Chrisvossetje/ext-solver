from dataclasses import dataclass
from typing import List, Self

import galois
from basis import Basis, BasisElement, BasisIndex, Grading, add_grade, comp_grade, sort_grades
from coalgebra import CoAlgebra, generate_tensored_moduled
import globals
from matrix import GradedMap, matrix_identity, one, reduce_to_pivots, zero
from comodule import CoModule
import numpy as np
from tensored import TensorIndex

@dataclass
class Morphism:
    domain: CoModule
    codomain: CoModule
    matrix: GradedMap
        
    def __matmul__(self, other: Self):
        assert self.domain == other.codomain, "Cannot compose morphisms with incompatible domains and codomains"
        
        for grade in other.matrix:
            if grade in self.matrix:
                other.matrix[grade] = self.matrix[grade] @ other.matrix[grade]
            else:
                if grade in self.codomain.basis:
                    other.matrix[grade] = zero(len(self.codomain.basis[grade]),len(other.domain.basis[grade]))
                else:
                    other.matrix[grade] = zero(0,len(other.domain.basis[grade]))
        return Morphism(other.domain, self.codomain, other.matrix)
    
    def __repr__(self):
        return ""
    
    def zero(domain: CoModule, codomain: CoModule):
        morph = {}
        for grade in domain.basis:
            target_length = 0
            if grade in codomain.basis:
                target_length += len(codomain.basis[grade])
            matrix = zero(target_length, len(domain.basis[grade]))
            morph[grade] = matrix
        return Morphism(domain, codomain, morph)
    
    def empty(domain: CoModule) -> Self:
        empty = domain.zero_module()
        morph = {}
        for grade in domain.basis:
            morph[grade] = zero(0, len(domain.basis[grade]))
        return Morphism(domain, empty, morph)

    def combine(f: Self, g: Self) -> tuple[Self, TensorIndex]:
        if globals.TEST:
            assert f.domain == g.domain
        
        # Matrix
        matrix: GradedMap = {}
        for grade in f.matrix:
            if grade in g.matrix:
                matrix[grade] = np.vstack((f.matrix[grade], g.matrix[grade]))
            else:
                matrix[grade] = f.matrix[grade]

        for grade in g.matrix:
            if grade not in matrix:
                matrix[grade] = g.matrix[grade]


        # Codomain Basis
        codomain = {}
        for grade in f.codomain.basis:
            if grade in g.codomain.basis:
                codomain[grade] = f.codomain.basis[grade] + g.codomain.basis[grade]
            else:
                codomain[grade] = f.codomain.basis[grade]

        for grade in g.codomain.basis:
            if grade not in codomain :
                codomain[grade] = g.codomain.basis[grade]
        
        # Codomain Coact 
        coaction: GradedMap = {}
        for grade in f.codomain.coaction:
            if grade in g.codomain.coaction:
                a = f.codomain.coaction[grade]
                b = g.codomain.coaction[grade]
                coaction[grade] = zero(a.shape[0] + b.shape[0], a.shape[1] + b.shape[1])
                coaction[grade][:a.shape[0], :a.shape[1]] = a
                coaction[grade][a.shape[0]:, a.shape[1]:] = b
            else:
                coaction[grade] = f.codomain.coaction[grade]

        for grade in g.codomain.coaction:
            if grade not in coaction:
                coaction[grade] = g.codomain.coaction[grade]
        

        # Tensored  
        tensored: TensorIndex = {}        
        for grade in f.codomain.tensored:
            if grade in g.codomain.tensored:
                temp = list(map(lambda x: (x[0], (x[1][0], len(f.codomain.basis[x[1][0]]) + x[1][1])), g.codomain.tensored[grade]))
                tensored[grade] = f.codomain.tensored[grade] + temp
            else:
                tensored[grade] = f.codomain.tensored[grade]
        for grade in g.codomain.tensored:
            if grade not in tensored:
                # TODO: Consider what happens when x[1][0] is not in f.codomain.basis
                def add_module_indices(x):
                    if x[1][0] in f.codomain.basis:
                        return (x[0], (x[1][0], len(f.codomain.basis[x[1][0]]) + x[1][1]))
                    else:
                        return x
                temp = list(map(add_module_indices, g.codomain.tensored[grade]))
                tensored[grade] = temp

        module = CoModule(f.codomain.coalgebra, codomain, coaction, tensored, None)
        return Morphism(f.domain, module, matrix)


    def verify(self):
        for grade in self.domain.basis:            
            if grade not in self.codomain.basis:
                assert self.matrix[grade].shape == (0, len(self.domain.basis[grade]))
            else:
                assert self.matrix[grade].shape == (len(self.codomain.basis[grade]), len(self.domain.basis[grade]))
        assert self.domain.coalgebra == self.codomain.coalgebra


    def __mat_mul__(self, other: Self):
        if globals.TEST:
            assert self.domain == other.codomain,\
                "Cannot compose morphisms with incompatible domains and codomains"
        for grade in other.matrix:
            other.matrix[grade] = self.matrix[grade] @ other.matrix[grade]
        return Morphism(other.domain, self.codomain, other.matrix)
    

    def structure_lines(self) -> List[tuple[BasisIndex,BasisIndex,int]]:
        # [(Domain F2 Generator BasisIndex, Codomain F2 Generator BasisIndex, h_i)]
        prims = []
        for prim_gr, prim_id in self.domain.primitive_indices():
            v = self.matrix[prim_gr][:, prim_id]
            for el_id in range(v.shape[0]):
                if v[el_id]:
                    gen = self.codomain.basis[prim_gr][el_id] 
                    if gen.generator:
                        primitive = self.domain.basis[prim_gr][prim_id].primitive
                        prim_generated_index = self.domain.basis[prim_gr][prim_id].generated_index
                        prim_gen = self.domain.find_generator(prim_generated_index)
                        prims.append((prim_gen, (prim_gr, el_id), primitive))
        return prims
    




    

def cokernel(F: Morphism) -> Morphism:
    M: GradedMap = {}
    for grade in F.codomain.basis:
        if grade in F.matrix:
            M[grade] = F.matrix[grade].left_null_space().row_reduce()
        else:
            M[grade] = matrix_identity(len(F.codomain.basis[grade]))
    Q_basis: Basis = {}
    pivots: dict[Grading, list[int]] = {}
    count = 0
    for grade in M:
        pivot = reduce_to_pivots(M[grade])
        if len(pivot) == 0:
            continue
        pivots[grade] = pivot
        Q_basis[grade] = [BasisElement(grade, str(count+n), False, None, -1) for n in pivots[grade]]
        count += len(pivot)

    alg = F.codomain.coalgebra
    Q_tensored, Q_moduled = generate_tensored_moduled(alg.basis, Q_basis)


    # Slightly slower for large coalgebras
    # def coact() -> GradedMap:
    #     coaction: GradedMap = {}

    #     for Q_grade in Q_basis:
    #         coaction[Q_grade] = zero(len(Q_tensored[Q_grade]) ,len(Q_basis[Q_grade]))
            
    #         for Q_index, F_index in enumerate(pivots[Q_grade]):
    #             F_coact_vec = F.codomain.coaction[Q_grade][:,F_index]
    #             for coact_id in range(F_coact_vec.shape[0]):
    #                 if F_coact_vec[coact_id]:
    #                     (alg_gr, alg_id), (mod_gr, mod_id) = F.codomain.tensored[Q_grade][coact_id]
    #                     target = M[mod_gr][:,mod_id]

    #                     for target_id in range(target.shape[0]):
    #                         if target[target_id]:
    #                             tensor_gr, tensor_id = Q_moduled[mod_gr][target_id][alg_gr][alg_id]
    #                             if globals.TEST:
    #                                 assert tensor_gr == Q_grade, "Resulting tensor grade not equal to original grade"
    #                             coaction[tensor_gr][tensor_id,Q_index] += target[target_id] * F_coact_vec[coact_id]
    #     return coaction
    

    def F_in_Q_map():
        F_in_Q_map = {}
        for grade in Q_basis:
            matrix = zero(len(F.codomain.basis[grade]), len(Q_basis[grade]))
            for Q_index, F_index in enumerate(pivots[grade]):
                matrix[F_index,Q_index] = one()
            F_in_Q_map[grade ]= matrix
        return F_in_Q_map

    def smart_coact() -> GradedMap:
        coaction: GradedMap = {}

        F_in_Q = F_in_Q_map()

        for Q_grade in Q_basis:
            coaction[Q_grade] = zero(len(Q_tensored[Q_grade]) ,len(Q_basis[Q_grade]))

            F_coact_vec = F.codomain.coaction[Q_grade] @ F_in_Q[Q_grade]
            for coact_id in range(F_coact_vec.shape[0]):
                if F_coact_vec[coact_id].any():
                    (alg_gr, alg_id), (mod_gr, mod_id) = F.codomain.tensored[Q_grade][coact_id]
                    target = M[mod_gr][:,mod_id]

                    for target_id in range(M[mod_gr].shape[0]):
                        if target[target_id]:
                            tensor_gr, tensor_id = Q_moduled[mod_gr][target_id][alg_gr][alg_id]
                            if globals.TEST:
                                assert tensor_gr == Q_grade, "Resulting tensor grade not equal to original grade"
                            coaction[tensor_gr][tensor_id] += target[target_id] * F_coact_vec[coact_id]
        return coaction


    Q = CoModule(F.codomain.coalgebra, Q_basis, smart_coact(), Q_tensored, Q_moduled)
    morph = Morphism(F.codomain, Q, M)
    morph.verify()
    return morph


def resolve(Q: CoModule, grade_limit: Grading) -> Morphism:
    growing_morphism: Morphism = Morphism.empty(Q)

    iteration = 0 # ascii for lowercase a
        
    grades = sort_grades(growing_morphism.matrix.keys(), grade_limit)
    prev_grade = 0

    while True:
        (Q_el, _, Q_index) = None, None, None
        for grade_id in range(prev_grade, len(grades)):
            grade = grades[grade_id]
            kernel = growing_morphism.matrix[grade].null_space()

            if kernel.shape[0] != 0:
                (Q_el, _, Q_index) = Q.lowest_graded_index_from_matrix(grade, kernel)
                prev_grade = grade_id
                break


        # No more elements in the kernel left
        if Q_el == None:
            break

        mapping_to_F: GradedMap = {}
        algebra_to_tensor = Q.moduled[Q_el.grading][Q_index]
        for a_gr in algebra_to_tensor:
            target_grade = add_grade(Q_el.grading, a_gr)
            if globals.TEST:
                assert len(algebra_to_tensor[a_gr]) == len(Q.coalgebra.basis[a_gr]),\
                        "Element is not represented by enough thingies"

            zero_matrix = zero(len(algebra_to_tensor[a_gr]) ,len(Q.basis[target_grade]))
            for a_id, (t_gr, t_id) in enumerate(algebra_to_tensor[a_gr]):
                if globals.TEST:
                    assert t_gr == target_grade, "Hmmm, grade should be equal"
                zero_matrix[a_id] = Q.coaction[t_gr][t_id]

            mapping_to_F[target_grade] = zero_matrix




        # F = CoModule.free_module(Q.coalgebra, Q_el.grading, iteration, chr(iteration + 97))
        F = CoModule.free_module_limit(Q.coalgebra, Q_el.grading, iteration, chr(iteration + 97), globals.ELEMENT_LIMIT)
        iteration += 1

        F_morphism = Morphism(Q, F, mapping_to_F)

        growing_morphism = Morphism.combine(growing_morphism, F_morphism)
    
    if globals.TEST:
        growing_morphism.verify()

    return growing_morphism

    