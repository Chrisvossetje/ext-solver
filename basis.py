from dataclasses import dataclass
from typing import List, Self

Grading = tuple[int,int]
GradeZero = (0,0)


def add_grade(lhs: Grading, rhs: Grading) -> Grading:
    return (lhs[0] + rhs[0], lhs[1] + rhs[1])

def mult_grade(grading: Grading, n: int) -> Grading:
    return (grading[0] * n, grading[1] * n)

def comp_grade(lhs: Grading, rhs: Grading) -> bool:
    return lhs[0] <= rhs[0]

def sort_grades(grades: List[Grading], grade_limit: Grading) -> List[Grading]:
    
    return list(filter(lambda gr: comp_grade(gr, grade_limit), sorted(grades)))

@dataclass
class BasisElement:
    grading: Grading
    name: str
    generator: bool
    primitive: None | int
    generated_index: int

    def tensor(a: Self, b: Self) -> Self:
        grade = add_grade(a.grading, b.grading)
        name = a.name + " | " + b.name
        return BasisElement(grade,name, False, None, -1)
    
    def __repr__(self):
        return str(self.name)
    

Basis = dict[Grading, List[BasisElement]]
BasisIndex = tuple[Grading, int]


