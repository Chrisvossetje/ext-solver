import copy
from typing import List

import galois

from basis import Basis, BasisElement, Grading, add_grade
from coalgebra import CoAlgebra, generate_tensored_moduled
from comodule import CoModule
import globals
from matrix import zero


def a_0_dual_temp():
    basis = {
        (0,0): [BasisElement((0,0), "1", True)],
        (1,0): [BasisElement((1,0), "chi", False)],
    }
    coaction = {
        (0,0): galois.GF2([[1]]),
        (1,0): galois.GF2([[1],[1]])
    }
    tensored = {
        (0,0): [(((0,0),0), ((0,0),0))],
        (1,0): [(((1,0),0), ((0,0),0)),
                (((0,0),0), ((1,0),0))],
    }
    return CoAlgebra(basis, coaction, tensored)


def el_to_grade(ls):
    grade = 0
    for (index, power) in enumerate(ls):
        grade += power * (2**index - 1)
    return grade

def element(gen, power , length):
    el = [0]*length
    el[gen] = power
    return el


def mult_els(a,b):
    c = [0]*len(a)
    for i in range(len(a)):
        c[i] = a[i] + b[i]
    return c

def simple(chi_i, length):
    tensors = []
    for i in range(0, chi_i+1):
        two_to_the_i = 2**i
        tensors += [(element(chi_i-i, two_to_the_i, length), element(i, 1, length))]
    return tensors

def mult_tens(a,b):
    (a1,a2) = a
    (b1,b2) = b
    return (mult_els(a1,b1), mult_els(a2,b2))

def el_to_index(el, divisors):
    lut = [0,1] + divisors
    location = 0
    mult = 1
    for i in range(1, len(divisors) + 1):
        mult *= lut[i]
        location += el[i] * mult 
    return location

def mon_tensor_to_index(tens, divisors, dim):
    a,b = tens
    return el_to_index(a, divisors) + el_to_index(b, divisors)*dim





def prrr(xy):
    x = xy[0]
    y = xy[1]
    if x == 0:
        return "1"
    else:
        return "ξ" + to_sub_script(x) + to_super_script(y)

def els_to_string(els):
    out = ""
    for i,el in enumerate(els[1:]):
        if el != 0:
            out += prrr((i+1,el))
    if out == "":
        return "1"
    return out


def to_super_script(x):
    subs = ["⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"]
    out = ""
    while x > 0:
        out += subs[x % 10]
        x //= 10
    return out

def to_sub_script(x):
    subs = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]
    out = ""
    while x > 0:
        out += subs[x % 10]
        x //= 10
    return out


def basis_to_tensor_basis():
    pass



def create_a_n_dual(n: int, max_limit = globals.GRADE_LIMIT[0]):
    divisors = []
    for i in range(n+1,0, -1):
        divisors += [2**i]

    one = [0]*(n+2)
    one[0] = 1

    def recurse(list, index, amount):
        new_list = []
        for a in range(amount):
            one = [0]*(n+2)
            one[index] = a
            for ls in list:
                newls = copy.copy(ls)
                newls[index] = a
                if el_to_grade(newls) <= max_limit:
                    new_list.append(newls)
        return new_list

    els = [one]
    for i in range(0, n+1):
        els = recurse(els, i+1, divisors[i])

    # basis: Basis = {}
    # for el in els:
    #     gr = (el_to_grade(el), 0)
    #     gen = False
    #     if gr == (0,0):
    #         gen = True
        
        
    basis = [BasisElement((el_to_grade(el),0), els_to_string(el), False, None, 0) for el in els]
    dim = len(basis)
    coaction = zero(dim*dim,dim)
    basis[0].generator = True

    for el in els:
        prod = []
        for enum, power in enumerate(el):
            for a in range(power):
                if prod == []:
                    prod = simple(enum, n+2)
                else:
                    new_prod = []
                    simp = simple(enum, n+2)
                    for a in prod:
                        for b in simp:                            
                            new_prod += [mult_tens(a,b)]
                    prod = new_prod
        
        # filter
        def non_zero_el(p):
            a,b = p
            for i in range(1, n+1):
                if a[i] >= divisors[i-1]:
                    return False
                if b[i] >= divisors[i-1]:
                    return False
            return True
        prod = filter(non_zero_el, prod)
        
        # def filter_too_high_grade(p): 
        #     if el_to_grade(p) > max_limit:
        #         return False
        #     return True
        # prod = filter(filter_too_high_grade, prod)

        for a in prod:
            t_id = mon_tensor_to_index(a, divisors, dim)
            el_id = el_to_index(el, divisors)
            coaction[t_id, el_id] += galois.GF2(1)

    return basis, coaction


def init_from_basis(basis: List[BasisElement]) -> tuple[Basis, List[tuple[Grading, int]]]:
    graded_basis: dict[Grading, List[BasisElement]] = {}
    other_thing: List[tuple[Grading, int]] = []
    for e in basis:
        if not e.grading in graded_basis:
            graded_basis[e.grading] = []
        other_thing.append((e.grading ,len(graded_basis[e.grading])))
        graded_basis[e.grading].append(e)

    return graded_basis, other_thing


def init_from_hopfalgebra(basis: List[BasisElement], coaction: galois.GF2) -> CoAlgebra:
        B, basis_to_graded_basis = init_from_basis(basis) 
        tensored, moduled = generate_tensored_moduled(B, B)
        
        graded_coaction = {}
        for grade in B:
            graded_coaction[grade] = zero(len(tensored[grade]) , len(B[grade]))


        for col in range(coaction.shape[1]):
            (el_grade, el_id) = basis_to_graded_basis[col]            
            for row in range(coaction.shape[0]):
                if coaction[row,col]:
                    first_id = row // coaction.shape[1]
                    second_id = row % coaction.shape[1]
                    gr1, id1 = basis_to_graded_basis[first_id]
                    gr2, id2 = basis_to_graded_basis[second_id]
                    grade = add_grade(gr1,gr2)

                    finalgr, finalid = moduled[gr1][id1][gr2][id2]

                    graded_coaction[finalgr][finalid , el_id] = coaction[row,col]


        return CoAlgebra(B, graded_coaction, tensored, 2)

def create_real_a_n_dual(n: int):
    basis, coact = create_a_n_dual(n)
    return init_from_hopfalgebra(basis, coact)
