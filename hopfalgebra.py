from basis import Basis, BasisElement, BasisIndex, Grading, add_grade, comp_grade, mult_grade
from coalgebra import CoAlgebra
from matrix import GradedMap, one, to_galois_coeff, zero
from tensored import TensorIndex, generate_tensored_moduled


Monomial = tuple[int]
CoactionElement = list[tuple[int, Monomial, Monomial]]
Generator = tuple[str, Grading]


def createPolynomialHopfAlgebra(field: int, generators: list[Generator], coactions: list[CoactionElement], relations: list[Monomial], max_grading: Grading) -> CoAlgebra:
    # let n be the number of generators
    n = len(generators)

    # first we create the basis of the hopf algebra as an F_2 vectorspace

    queue: list[tuple[Monomial, int]] = []
    basis_information: dict[Monomial, CoactionElement] = {None: None}  # we add None to the dict since if a function returns none we can simply skip

    one_monomial = tuple([0] * n) # all exponents zero => 1
    basis_information[one_monomial] = [(1,one_monomial,one_monomial)]
    
    for i in range(n):
        queue.append((one_monomial,i))

    while len(queue) > 0:
        v, i = queue.pop(0)
        next = multiply_monomial_by_generator(v,i,relations)
        if next not in basis_information:
            if comp_grade(monomial_to_grade(next, generators), max_grading):
                basis_information[next] = multiply_coaction_elements(basis_information[v], coactions[i][1], relations, field)
                for i in range(n):
                    queue.append((next,i))
    basis_information.pop(None)

    # we construct the basis from the basis_information
    # we also construct a map from each monomial to its index in the basis
    basis: Basis = {}
    basis_index: dict[Monomial, BasisIndex] = {}
    for monomial in basis_information.keys():
        # check if basis_index has grading as key
        grading = monomial_to_grade(monomial, generators)
        basis_element = BasisElement(grading, monomial_to_string(monomial, generators), monomial == one_monomial, None, 0)
        if grading not in basis:
            basis[grading] = []
        basis_index[monomial] = (grading, len(basis[grading]))
        basis[grading].append(basis_element)



    
    ### old code
    # tensored: TensorIndex = {}
    # coaction: GradedMap = {}
    # for monomial, coaction_element in basis_information.items():
    #     grade, index = basis_index[monomial]
    #     if grade not in coaction:
    #         coaction[grade] = zero(0, len(basis[grade]))
    #     if grade not in tensored:
    #         tensored[grade] = []
    #     for a, b in enumerate(coaction_element): # we get an elemenet a|b
    #         a_grade, a_index = basis_index[a]
    #         b_grade, b_index = basis_index[b]
    #         if 
    #             row = zero(1, len(basis[grade]))
    #             row[0,i] = one()
    #             coaction[grade] = coaction[grade].vstack(row)

    # Module Grade + Index -> Algebra Grading + index -> Tensor Grading + index
    ModuleIndex = dict[Grading, list[dict[Grading, list[BasisIndex]]]]

    # Tensor Grade + Index -> (Algebra Grading + index, Module Grading + index)
    TensorIndex = dict[Grading, list[tuple[BasisIndex, BasisIndex]]]

    tensored, moduled = generate_tensored_moduled(basis, basis)

    # now we create the coaction map
    coaction: GradedMap = {}
    for grade in tensored:
        coaction[grade] = zero(len(tensored[grade]), len(basis[grade]))
    
    for monomial, coaction_element in basis_information.items():
        grade, index = basis_index[monomial]
        for c, a, b in coaction_element: # we get an elemenet c.a|b
            a_grade, a_index = basis_index[a]
            b_grade, b_index = basis_index[b]
            _, tensor_index = moduled[b_grade][b_index][a_grade][a_index]
            coaction[grade][tensor_index, index] = to_galois_coeff(c)


    return CoAlgebra(basis, coaction, tensored, field)


    


def monomial_to_grade(m: Monomial, generators: list[Generator]) -> Grading:
    total_grade = (0,0)
    for exp, (_, grade) in zip(m, generators):
        grade = mult_grade(grade, exp)
        total_grade = add_grade(total_grade, grade)
    return total_grade

def monomial_to_string(m: Monomial, generators: list[Generator]) -> str:
    name = ""
    for exp, (generator_name, _) in zip(m, generators):
        if exp != 0:
            name += generator_name + "^" + str(exp) + " "
    if name == "":
        return "1"
    return name

def multiply_tensor_terms(a: tuple[int, Monomial, Monomial], b: tuple[int, Monomial, Monomial], relations: list[Monomial], field: int) -> tuple[int, Monomial, Monomial]:
    a0, a1, a2 = a
    b0, b1, b2 = b
    c1 = multiply_monomials(a1,b1,relations)
    c2 = multiply_monomials(a2,b2,relations)
    if c1 == None or c2 == None:
        return None
    else:
        return ((a0*b0) % field,c1,c2)
    
def multiply_coaction_elements(a: CoactionElement, b: CoactionElement, relations: list[Monomial], field: int) -> CoactionElement:
    tensor_terms: dict[tuple[Monomial,Monomial],int] = {}
    for x in a:
        for y in b:
            c = multiply_tensor_terms(x,y,relations, field)
            if c != None:
                c0,c1,c2 = c
                if (c1,c2) in tensor_terms:
                    tensor_terms[(c1,c2)] += c0
                    tensor_terms[(c1,c2)] %= field
                else:
                    tensor_terms[(c1,c2)] = c0
    new = []
    for (c1,c2), c0 in tensor_terms.items():
        if c0 != 0:
            new.append((c0,c1,c2))
    return new



def increment_tuple(t: Monomial, index: int) -> Monomial:
    new = list(t)
    new[index] += 1
    return tuple(new)

def multiply_monomial_by_generator(m: Monomial, index: int, relations: list[Monomial]) -> Monomial:
    m = increment_tuple(m, index)
    m = mod_relations(m, relations)
    return m
    


def add_tuples(a: Monomial, b: Monomial) -> Monomial:
    new = []
    for x, y in zip(a,b):
        new.append(x + y)
    return tuple(new)

def multiply_monomials(a: Monomial, b: Monomial, relations: list[Monomial]) -> Monomial:
    c = add_tuples(a,b)
    c = mod_relations(c, relations)
    return c



def monomial_division_check(a: Monomial, b: Monomial) -> bool:
    for x, y in zip(a,b):
        if x < y:
            return False
    return True

def mod_relations(m: Monomial, relations: list[Monomial]) -> Monomial:
    for r in relations:
        if monomial_division_check(m,r):
            return None
    return m





def parse_monomial(name, generator_translate, size):
    els = name.split("*")
    mon = [0]*size
    for el in els:
        ls = el.split("^")
        
        if len(ls) == 1:
            expo = 1
        elif len(ls) == 2:
            expo = int(ls[1])
        else:
            print("Too many elements after splitting ^")
            exit(1)

        name = ls[0].strip()
        if name == '1':
            continue
        index = generator_translate[name]
        assert name in generator_translate, name + ": does not appear in the generators"
        mon[index] = expo

    return mon

def HopfAlgebraParse(filename: str):
    with open(filename) as f:
            input = f.readlines()
    state = -1
    field = None
    generators = []
    generator_translate = {}
    relations = []
    coaction = []

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

        elif i[0:7] == "- GENERATOR"[0:7]:
            if state != 0:
                print("Previous state was not correct, expected field to be parsed first")
                exit(1)
            state = 1
        
        elif i[0:7] == "- RELATIONS"[0:7]:
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
                # assert field == 2, "p != 2 field is not yet supported for processing hopfalgebra"

                
            elif state == 1:
                a,b = i.split(":")
                name = a.strip()
                grade = b.strip()
                t,v = grade.lstrip("(").rstrip(")").split(",")
                generator_translate[name] = len(generators)
                generators.append((name,(int(t),int(v))))

            elif state == 2:
                relations.append(parse_monomial(i, generator_translate, len(generators)))
                
            elif state == 3:
                a,b = i.split(":")
                name = a.strip()
                tensors = b.split("+")
                assert name == generators[len(coaction)][0],\
                    "Coaction should be given for the generators in the same order"
                ts = []
                for t in tensors:
                    t = t.strip()
                    if field == 2:
                        l,r = t.split("|")
                        lp = tuple(parse_monomial(l.strip(), generator_translate, len(generators)))
                        rp = tuple(parse_monomial(r.strip(), generator_translate, len(generators)))
                        ts.append((1, lp, rp))
                    else:
                        v,t = t.split("*")
                        l,r = t.split("|")
                        lp = tuple(parse_monomial(l.strip(), generator_translate, len(generators)))
                        rp = tuple(parse_monomial(r.strip(), generator_translate, len(generators)))
                        ts.append((int(v.strip()), lp, rp))
                coaction.append((len(coaction), ts))

    if state != 3:
        print("Coalgebra defintion is not complete")
        exit(1)

    return field, generators, relations, coaction

