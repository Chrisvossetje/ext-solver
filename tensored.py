
# Module Grading + index -> Algebra Grade + index -> Tensor Grade + index
from typing import List
from basis import Basis, BasisIndex, Grading, add_grade
import globals


# Module Grade + Index -> Algebra Grading + index -> Tensor Grading + index
ModuleIndex = dict[Grading, List[dict[Grading, List[BasisIndex]]]]



# Tensor Grade + Index -> (Algebra Grading + index, Module Grading + index)
TensorIndex = dict[Grading, List[tuple[BasisIndex, BasisIndex]]]



def generate_tensored_moduled(algebra: Basis, module: Basis) -> tuple[TensorIndex, ModuleIndex]: 
    tensored: TensorIndex = {}
    moduled: ModuleIndex = {}
    for grade_m in module:
        if grade_m not in moduled:
            moduled[grade_m] = []    

        for m_id, m_el in enumerate(module[grade_m]):
            alg_temp_to_total_index = {}

            for grade_a in algebra:
                grade = add_grade(grade_m, grade_a)
                if grade not in module:
                    continue
                if grade not in tensored:
                    tensored[grade] = []

                alg_temp_to_total_index[grade_a] = []
                for a_id, a_el in enumerate(algebra[grade_a]):
                    alg_temp_to_total_index[grade_a].append((grade,len(tensored[grade])))
                    tensored[grade].append(((a_el.grading, a_id), (m_el.grading, m_id)))
            moduled[grade_m].append(alg_temp_to_total_index)

    if globals.TEST:
        verify_moduled_tensored(tensored, moduled)
        verify_tensored(algebra, module, tensored)
    return (tensored, moduled)



# Verify bijective property of the two indexes
def verify_moduled_tensored(tensored: TensorIndex, moduled: ModuleIndex):
    # Tensored -> Moduled -> Tensored
    for t_grade in tensored:
        for t_id, t_el in enumerate(tensored[t_grade]):
            a,m = t_el
            m_gr, m_el = m
            a_gr, a_el = a
            compare_gr, compare_id = moduled[m_gr][m_el][a_gr][a_el]
            assert t_grade == compare_gr and add_grade(m_gr, a_gr) == t_grade,\
                "Tensored and module grading is not compatible"
            assert t_id == compare_id,\
                "Tensored and module index is not compatible"
    
    # Moduled -> Tensored -> Moduled
    for m_gr in moduled:
        for m_id, m_el in enumerate(moduled[m_gr]):
            alg = moduled[m_gr][m_id]
            for a_gr in alg:
                for a_id, a_el in enumerate(alg[a_gr]):
                    t_gr, t_id = alg[a_gr][a_id]
                    a,m = tensored[t_gr][t_id]
                    assert a == (a_gr,a_id) and add_grade(m_gr, a_gr) == t_gr,\
                        "Tensored and module grading is not compatible"
                    assert m == (m_gr, m_id),\
                        "Tensored and module index is not compatible"

def verify_tensored(algebra: Basis, module: Basis, tensored: TensorIndex):
    for t_grade in tensored:
        for _, t_el in enumerate(tensored[t_grade]):
            a,m = t_el
            m_gr, m_el = m
            a_gr, a_el = a
            assert a_gr in algebra, "algebra does not have this grading"
            assert a_el < len(algebra[a_gr]) and 0 <= a_el, "algebra does not have this grading"
            assert m_gr in module, "module does not have this grading"
            assert m_el < len(module[m_gr]) and 0 <= m_el, "algebra does not have this grading"


                
