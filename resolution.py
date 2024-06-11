

from dataclasses import dataclass
from typing import List
from basis import Grading, add_grade
from coalgebra import generate_tensored_moduled
from comodule import CoModule
import globals
from morphism import Morphism, cokernel, resolve

@dataclass
class Resolution:
    comodule: CoModule
    morphisms: List[Morphism]

    def __post_init__(self):
        if globals.TEST:
            self.verify()

    def verify(self):
        for i in range(1, len(self.morphisms)):
            assert self.morphisms[i].domain == self.morphisms[i-1].codomain,\
                "Domain and Codomains do not match up"
    

    def __str__(self):
        result = "Resolution" + ":\n"
        for i, morphism in enumerate(self.morphisms):
            result += "d_" + str(i) + " : " + str(morphism) + " | "
            for g in morphism.codomain.generators():
                result += str(g.grading[0] - i + 1) + ", "
            result += "\n"

        result += "\n"
        result += "0 --> F_2"
        for i in range(1, len(self.morphisms)):
            comodule_symbol = self.morphisms[i].codomain.symbol()
            result += " --> " + comodule_symbol
        return result

    def grading(self):
        gradings = []
        for i, morphism in enumerate(self.morphisms[1:]):
            temp = []
            for g in morphism.codomain.generators():
                temp.append(add_grade(g.grading, (-i,0)))
            gradings.append(temp)
        return gradings



def resolution(M: CoModule) -> Resolution:
    zero: Morphism = Morphism.zero(M.zero_module(), M)
    morphisms: List[Morphism] = [zero]

    for n in range(globals.FILTRATION_MAX):
        morph = morphisms[-1]
        print("Calculating cokernel  ", len(morphisms))
        coker = cokernel(morph)


        print("Calculating injection ", len(morphisms))
        injection_to_cofree = resolve(coker.codomain, add_grade(globals.GRADE_LIMIT, (n,0)))

        # "Forget" about the cokernel
        final = injection_to_cofree @ coker

        
        morphisms.append(final)

    return Resolution(M, morphisms)


