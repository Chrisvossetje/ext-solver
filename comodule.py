import copy
from dataclasses import dataclass
from typing import List, Self

import galois
from basis import Basis, BasisElement, BasisIndex, GradeZero, Grading, add_grade, comp_grade
from coalgebra import CoAlgebra 
import globals
from matrix import GradedMap, calculate_zeros_in_graded_map, zero
from tensored import ModuleIndex, TensorIndex, generate_tensored_moduled, verify_moduled_tensored, verify_tensored

@dataclass
class CoModule:
    coalgebra: CoAlgebra
    basis: Basis
    coaction: GradedMap
    
    tensored: TensorIndex
    moduled: ModuleIndex

    def __post_init__(self):
        if self.tensored == None:
            self.tensored, self.moduled = generate_tensored_moduled(self.coalgebra.basis, self.basis)
        if globals.TEST:
            if self.moduled != None:
                verify_moduled_tensored(self.tensored, self.moduled)
        
        if globals.TEST:
            self.test_coaction()
            verify_tensored(self.coalgebra.basis, self.basis, self.tensored)




    def __repr__(self) -> str:
        return str(self.dim()) + " Elements | " +  str(len(self.generators())) + " Generators"

    # Vector space dimension
    def dim(self) -> int:
        count = 0
        for grade in self.basis:
            count += len(self.basis[grade])
        return count

    def generate_luts(self):        
        self.tensored, self.moduled = generate_tensored_moduled(self.coalgebra.basis, self.basis)


    # TODO: Implement this for Cokernel module,
    # This should ALSO update the moduled thingy
    def reduce(self):        
        for grade in self.basis:
            m = self.coaction[grade]
            indices = []
            for i in range(m.shape[0]):
                if m[i].any():
                    indices.append(i)
            
            new_t = list(map(lambda x: self.tensored[grade][x], indices))
            self.tensored[grade] = new_t
            coact = zero(len(indices), len(self.basis[grade]))
            for new,old in enumerate(indices):
                coact[new] = m[old]

            self.coaction[grade] = coact

    def test_coaction(self):
        for grade in self.coaction:
            a_size = len(self.basis[grade])
            t_size = len(self.tensored[grade])
            assert self.coaction[grade].shape == (t_size, a_size)

        for grade in self.basis:
            rs, _ = self.coaction[grade].shape
            for r in range(rs):
                a, m = self.tensored[grade][r]
                if a[0] == GradeZero:
                    assert self.coaction[grade][r,m[1]] == 1,\
                        "coaction(x) does not contain 1⊕x"

    def tensor_str(self, tens_index):
        if globals.TEST:
            assert tens_index >= 0 and tens_index < self.dim() * self.coalgebra.dim(),\
                    "Tensor index is outside the range of this algebra"
        first = '{:25s} '.format(self.coalgebra.basis[tens_index % self.coalgebra.dim()].name)
        return first + " ⊗ " + self.basis[tens_index // self.coalgebra.dim()].name

    def zero_module(self):
        return CoModule(self.coalgebra, {}, {}, {}, {})
    

    def fp_module(coalgebra: CoAlgebra, grade: Grading = GradeZero):
        # Assumes there is 1 element in CoAlgebra AND it is in grade (0,0)
        basis = {(0,0): [BasisElement(grade, "F_p", False, None, 0)]}
        coaction = {
            (0,0): galois.GF(globals.FIELD)([[1]]), 
        }

        return CoModule(coalgebra, basis, coaction, None, None)
    
    def tensor_basis(self):
        tensor = {}
        for grade in self.basis:
            temp = []
            for a,m in self.tensored[grade]:
                a_el = self.coalgebra.basis[a[0]][a[1]]
                m_el = self.basis[m[0]][m[1]]
                temp.append(BasisElement.tensor(a_el, m_el))
            tensor[grade] = temp
        return tensor

    
    def free_module(coalgebra: CoAlgebra, grade: Grading, index: int, name: str) -> Self:
        tensored = {add_grade(tens_grade, grade):[(el[0], (add_grade(el[1][0],grade),el[1][1]))\
                            for el in coalgebra.tensored[tens_grade]]\
                                for tens_grade in coalgebra.tensored }

        basis = {add_grade(b_grade, grade):[BasisElement(add_grade(el.grading, grade),\
                                        el.name + " | "  + name,\
                                        el.generator, el.primitive, index) 
                            for el in coalgebra.basis[b_grade]]\
                                for b_grade in coalgebra.basis }
        
        coact = {add_grade(c_grade,grade):coalgebra.coaction[c_grade] 
                            for c_grade in coalgebra.coaction }

        module = CoModule(coalgebra, basis, coact, tensored, None)
        return module
    
    def free_module_limit(coalgebra: CoAlgebra, grade: Grading, index: int, name: str, limit: Grading) -> Self:
        tensored = {add_grade(tens_grade, grade):[(el[0], (add_grade(el[1][0],grade),el[1][1]))\
                            for el in coalgebra.tensored[tens_grade]]\
                                for tens_grade in coalgebra.tensored if comp_grade(tens_grade, limit)}

        basis = {add_grade(b_grade, grade):[BasisElement(add_grade(el.grading, grade),\
                                        el.name + " | "  + name,\
                                        el.generator, el.primitive, index) 
                            for el in coalgebra.basis[b_grade]]\
                                for b_grade in coalgebra.basis if comp_grade(b_grade, limit)}
        
        coact = {add_grade(c_grade,grade):coalgebra.coaction[c_grade] 
                            for c_grade in coalgebra.coaction if comp_grade(c_grade, limit) }

        module = CoModule(coalgebra, basis, coact, tensored, None)
        return module

    def symbol(self) -> str:
        if len(self.basis) == 0:
            return "0"
        As = ["A"]*(len(self.generators()))
        return "⊕".join(As)


    def find_generator(self, generated_index: int):
        for grade in self.basis:
            for id, el in enumerate(self.basis[grade]):
                if el.generator and el.generated_index == generated_index:
                    return (grade,id)
        assert False, "This grade + id does not exist :("


    def generators(self) -> List[BasisElement]:
        gens = []
        for grade in self.basis:
            gens += list(filter(lambda x: x.generator, self.basis[grade]))
        return gens
    
    
    def primitives(self) -> List[BasisElement]:
        prims = []
        for grade in self.basis:
            prims += list(filter(lambda x: x.primitive != None, self.basis[grade]))
        return prims
    
    
    def primitive_indices(self) -> List[BasisIndex]:
        prims = []
        for grade in self.basis:
            for p_id, p_el in enumerate(self.basis[grade]):
                if p_el.primitive != None:
                    prims.append((grade,p_id))
        return prims

    def lowest_graded_index_from_matrix(self, grade: Grading, matrix: galois.GF2) -> tuple[BasisElement, int, int]:
        rows, cols = matrix.shape
        
        row, col = (-1,-1)
        lowest_element = None

        for r in range(rows):
            for c in range(cols):
                if matrix[r,c]:
                    el = self.basis[grade][c]                    
                    # TODO: 
                    # BECAUSE IT IS ALL GRADED THIS SHOULD FREELY WORK
                    if lowest_element == None or el.grading[0] < el.grading[0]:
                        lowest_element = el
                        (row,col) = (r,c)
                        return (lowest_element, row, col)   
            
        
        if globals.TEST:
            assert lowest_element != None, "No vector with a grade has been found ?"
        
        return (lowest_element, row, col)    
    
    def calculate_zeros_in_coaction(self) -> tuple[str, bool]:
        return calculate_zeros_in_graded_map(self.coaction)

    def parse(filename: str, coalg: CoAlgebra, coalg_translate: dict[str, BasisIndex]) -> Self:
        with open(filename) as f:
            input = f.readlines()
        state = -1
        field = None
        basis = []
        coaction_lut = []
        field_constructor = None

        for i in input:
            i = i.strip()
            if len(i) == 0:
                continue
            if i[0] == "#":
                continue

            if i[0:7] == "- FIELD"[0:7]:
                if state != -1:
                    print("Field should be first to parse")
                    exit(1)
                state = 0

            elif i[0:7] == "- BASIS"[0:7]:
                if state != 0:
                    print("Previous state was not correct, expected field to be parsed first")
                    exit(1)
                state = 1
            
            elif i[0:7] == "- COACTION"[0:7]:
                if state != 1:
                    print("Previous state was not correct, expected basis to be parsed first")
                    exit(1)
                state = 2
            
            else:
                if state == 0:
                    if field != None:
                        print("Previous state was not correct, expected generators to be parsed first")
                        exit(1)
                    field = int(i)
                    field_constructor = galois.GF(field)

                elif state == 1:
                    a,b = i.split(":")
                    name = a.strip()
                    grade = b.strip()
                    t,v = grade.lstrip("(").rstrip(")").split(",")
                    basis.append((name,(int(t),int(v))))

                elif state == 2:
                    a,b = i.split(":")
                    name = a.strip()
                    tensors = b.split("+")
                    ts = []
                    for t in tensors:
                        t = t.strip()
                        if field == 2:
                            l,r = t.split("|")
                            ts.append((1, l.strip(),r.strip()))
                        else:
                            v,t = t.split("*")
                            l,r = t.split("|")
                            ts.append((v.strip(), l.strip(),r.strip()))
                    coaction_lut.append((name, ts))

        basis_dict: dict[str, BasisElement] = {}
        for b, gr in basis:
            if b in basis_dict:
                print("Name in basis appears twice")
                exit(1)
            basis_dict[b] = BasisElement(gr, b, False, None, 0)

        transformed = {}
        basis_translate = {}

        for b in basis_dict:
            el = basis_dict[b]
            gr = el.grading
            if gr not in transformed:
                transformed[el.grading] = []
            basis_translate[b] = (el.grading,len(transformed[el.grading]))
            transformed[el.grading].append(el)

        tensored, moduled = generate_tensored_moduled(transformed, transformed)
        coaction = {}
        for gr in transformed:
            coaction[gr] = zero(len(tensored[gr]), len(transformed[gr]))


        for b, ls in coaction_lut:
            gr, id = basis_translate[b]
            for v in ls:
                scalar, l, r = v
                l_gr, l_id = coalg_translate[l]
                r_gr, r_id = basis_translate[r]
                assert add_grade(l_gr,r_gr) == gr, "Grades are not homogenous"
                t_gr, t_id = moduled[r_gr][r_id][l_gr][l_id]
                assert t_gr == gr, "Tensored/Moduled is weird"
                coaction[gr][t_id,id] = field_constructor(scalar)

        return CoModule(coalg, transformed, coaction, tensored, moduled)