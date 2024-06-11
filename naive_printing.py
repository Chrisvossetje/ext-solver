from comodule import CoModule
from matrix import GradedMap, zero
from morphism import Morphism


# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS
# THIS FILE IS FOR BACKUP REASONSSS

def simple_coact(self: CoModule, filename: str):
        # els = []
        module_trans = {}
        count = 0
        for grade in self.basis:
            module_trans[grade] = []
            for _ in self.basis[grade]:
                module_trans[grade].append(count)
                count += 1

        alg_trans = {}
        count = 0
        for grade in self.coalgebra.basis:
            alg_trans[grade] = []
            for _ in self.coalgebra.basis[grade]:
                alg_trans[grade].append(count)
                count += 1
        
        mod_dim = self.dim()
        alg_dim = self.coalgebra.dim()
        
        coact = zero(mod_dim*alg_dim, mod_dim)
        for gr in self.coaction:
            rs,cs = self.coaction[gr].shape
            for r in range(rs):
                for c in range(cs):
                    if self.coaction[gr][r,c]:
                        a,m = self.tensored[gr][r]
                        a_gr, a_id = a
                        m_gr, m_id = m
                        m_count = module_trans[m_gr][m_id]
                        a_count = alg_trans[a_gr][a_id]
                        t_count = m_count*alg_dim + a_count
                        real_m = module_trans[gr][c]
                        coact[t_count,real_m] = self.coaction[gr][r,c]

        counts = []

        rs,cs = coact.shape

        for c in range(cs):
            count = 0
            for r in range(rs):
                if coact[r,c]:
                    count += 1
            counts.append(count) 
        
        counts_t = []

        for r in range(rs):
            count = 0
            for c in range(cs):
                if coact[r,c]:
                    count += 1
            counts_t.append(count) 

        with open(filename + ".txt", "w") as f:
            f.write(str(coact))
            
            f.write("\n")
            f.write("\n")
            f.write("Counts:\n")
            f.write(str(counts))
            
            f.write("\n")
            f.write("\n")
            f.write("Sorted:\n")
            counts.sort()
            f.write(str(counts))
            
            f.write("\n")
            f.write("\n")
            f.write("Counts:\n")
            f.write(str(counts_t))

            f.write("\n")
            f.write("\n")
            f.write("Sorted:\n")
            counts.sort()
            f.write(str(counts_t))

            f.write("\n")
            f.write("\n")
            f.write("Summed:\n")
            summed = sum(counts_t)
            f.write(str(summed))

        return coact
    
def simple_matrix(self: Morphism, filename: str):
    module_els = []
    for grade in self.domain.basis:
        for id, el in enumerate(self.domain.basis[grade]):
            module_els.append((grade,id))

    alg_els = []
    for grade in self.codomain.basis:
        for id, el in enumerate(self.codomain.basis[grade]):
            alg_els.append((grade,id))
    
    dom_dim = self.domain.dim()
    codom_dim = self.codomain.dim()
    
    matrix = zero(codom_dim, dom_dim)
    rs,cs = matrix.shape
    for r in range(rs):
        for c in range(cs):
            dom_gr, dom_id = module_els[c]
            codom_gr, codom_id = alg_els[r]
            
            if dom_gr == codom_gr:
                matrix[r,c] = self.matrix[dom_gr][codom_id,dom_id]


    counts = []

    for c in range(cs):
        count = 0
        for r in range(rs):
            if matrix[r,c]:
                count += 1
        counts.append(count) 
    
    counts_t = []

    for r in range(rs):
        count = 0
        for c in range(cs):
            if matrix[r,c]:
                count += 1
        counts_t.append(count) 

    with open(filename + ".txt", "w") as f:
        f.write(str(matrix))
        
        f.write("\n")
        f.write("\n")
        f.write("Counts:\n")
        f.write(str(counts))
        
        f.write("\n")
        f.write("\n")
        f.write("Sorted:\n")
        counts.sort()
        f.write(str(counts))

        f.write("\n")
        f.write("\n")
        f.write("Counts:\n")
        f.write(str(counts_t))

        f.write("\n")
        f.write("\n")
        f.write("Sorted:\n")
        counts.sort()
        f.write(str(counts_t))


        f.write("\n")
        f.write("\n")
        f.write("Summed:\n")
        summed = sum(counts_t)
        f.write(str(summed))
        
        
        f.write("\n")
        f.write("\n")
        f.write("Domain:\n")
        f.write(str(list(map(lambda x: self.domain.basis[x[0]][x[1]], module_els))))
        f.write("\n")
        f.write("\n")
        f.write("Codomain:\n")
        f.write(str(alg_els))
    
    return matrix


def print_stuff(domain: Basis, codomain: Basis, matrix: GradedMap, filename: str):
    if globals.PRINT == False:
        return
    
    with open("export/" + filename + ".txt", "w") as f:

        for grade in domain:
            f.write("\n\n\n\n------------------------------------------\n\n")
            f.write("Grade:\n" + str(grade))
            f.write("\n")
            f.write("\n")
            f.write("\n")
            f.write("Matrix:\n")
            f.write(str(matrix[grade]))
            f.write("\n")
            f.write("\n")
            f.write("Domain:\n")
            f.write(str(domain[grade]))
            f.write("\n")
            f.write("\n")

            if grade in codomain:
                f.write("Codomain:\n")
                f.write(str(codomain[grade]))
                f.write("\n")
                f.write("\n")