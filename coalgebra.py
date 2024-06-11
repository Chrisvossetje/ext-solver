from dataclasses import dataclass
from typing import List
import galois

import globals
from basis import Basis, BasisElement, BasisIndex, GradeZero, add_grade
from matrix import GradedMap, calculate_zeros_in_graded_map, zero
from tensored import TensorIndex, generate_tensored_moduled





@dataclass
class CoAlgebra:
    basis: Basis
    coaction: GradedMap
    tensored: TensorIndex
    field: int

    def __post_init__(self):
        self.set_primitives()
        if self.tensored == None:
            self.tensored, _ = generate_tensored_moduled(self, self.basis)
        self.reduce()

        if globals.TEST:
            assert self.calculate_zeros_in_coaction()[1], "Coaction has zero elements"
            self.test()
        assert globals.FIELD == self.field, "Field of imported coalgebra is NOT the same as the expected field in globals.py"

    def dim(self):
        count = 0
        for grade in self.basis:
            count += len(self.basis[grade])
        return count
    
    def test(self):
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



    def set_primitives(self):
        primitive_index = 0
        for grade in self.basis:
            for (index, el) in enumerate(self.basis[grade]):
                v = self.coaction[grade][:,index]
                count = 0
                for a in range(v.shape[0]):
                    if v[a]:
                        count += 1
                        if count > 2:
                            break
                if count == 2:
                    self.basis[grade][index].primitive = primitive_index
                    primitive_index += 1
    
    def calculate_zeros_in_coaction(self) ->  tuple[str, bool]:
        return calculate_zeros_in_graded_map(self.coaction)

    def parse(filename: str):
        with open(filename) as f:
            input = f.readlines()
        state = -1
        field = None
        basis = []
        generator = []
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
            
            elif i[0:7] == "- GENERATOR"[0:7]:
                if state != 1:
                    print("Previous state was not correct, expected basis to be parsed first")
                    exit(1)
                state = 2
            
            elif i[0:7] == "- COACTION"[0:7]:
                if state != 2:
                    print("Previous state was not correct, expected generators to be parsed first")
                    exit(1)
                state = 3
            
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
                    generator.append(i)

                elif state == 3:
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

        for g in generator:
            if g not in basis_dict:
                print("Generator name is not in basis")
                exit(1)
            basis_dict[g].generator = True

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
                l_gr, l_id = basis_translate[l]
                r_gr, r_id = basis_translate[r]
                assert add_grade(l_gr,r_gr) == gr, "Grades are not homogenous"
                t_gr, t_id = moduled[r_gr][r_id][l_gr][l_id]
                assert t_gr == gr, "Tensored/Moduled is weird"
                coaction[gr][t_id,id] = field_constructor(scalar)

        if state != 3:
            print("Coalgebra defintion is not complete")
            exit(1)

        return CoAlgebra(transformed, coaction, tensored, field), basis_translate
    
    def serialize(self, filename: str):
        with open(filename,'w') as f:
            f.write("- FIELD\n")
            f.write(str(self.field) + "\n")
            f.write("\n")
            f.write("- BASIS\n")
            for grade in self.basis:
                for el in self.basis[grade]:
                    f.write(basis_element_name(el) + " : " + str(el.grading) + "\n")
            f.write("\n")
            f.write("- GENERATOR\n")
            for grade in self.basis:
                for el in self.basis[grade]:
                    if el.generator:
                        f.write(basis_element_name(el) + "\n")
            f.write("\n")
            f.write("- COACTION\n")
            f.write(self.stringify_coaction())
            f.write("\n")

    def stringify_coaction(self) -> str:
        out = ""
        for grade in self.coaction:
            for i in range(self.coaction[grade].shape[1]):
                out += basis_element_name(self.basis[grade][i]) + " : "
                for j in range(self.coaction[grade].shape[0]):
                    if self.coaction[grade][j,i] != 0:
                        (algebra_grading, algebra_index), (module_grading,module_index) = self.tensored[grade][j]
                        if globals.FIELD != 2:
                            out += str(self.coaction[grade][j,i]) + "*" \
                                + basis_element_name(self.basis[algebra_grading][algebra_index]) + "|" \
                                + basis_element_name(self.basis[module_grading][module_index]) + " + "
                        else:
                            out += basis_element_name(self.basis[algebra_grading][algebra_index]) + "|" \
                                + basis_element_name(self.basis[module_grading][module_index]) + " + "
                out = out[:-3] + "\n"
        return out

    def __repr__(self):
        out = ""
        for list in self.basis.values():
            for b in list:
                out += str(b.name) + "\n"
        return out
    
def basis_element_name(b: BasisElement) -> str:
    return str(b.name).replace(" ","")
