from basis import Grading
import globals
import galois

GradedMap = dict[Grading, galois.GF2]


def one():
    return galois.GF(globals.FIELD)(1)

def zero(rows: int, cols: int) -> galois.GF2:
    return galois.GF(globals.FIELD).Zeros((rows, cols))

def to_galois_coeff(coeff: int) -> galois.GF2:
    return galois.GF(globals.FIELD)(coeff)

def matrix_identity(dim):
    return galois.GF(globals.FIELD).Identity(dim)

def reduce_to_pivots(T: galois.GF2) -> list[int]:
    if T.shape[0] == 0:
        return []
    rank = 0
    pivot_indices = []
    for i in range(T.shape[1]):
        if T[rank,i] == galois.GF(globals.FIELD)(1):
            pivot_indices.append(i)
            rank += 1
            if rank == T.shape[0]:
                break
    return pivot_indices

def calculate_zeros_in_graded_map(graded_map) -> tuple[str, bool]:
    count = 0
    total = 0
    for grade in graded_map:
        m = graded_map[grade]
        for i in range(m.shape[0]):
            total += 1
            if not m[i].any():
                count += 1
    if total == 0:
        return "", True
    return str(total) + " | " + str(count) + " | " + str(count/total), count == 0