import cProfile
import pstats
from coalgebra import CoAlgebra
from comodule import CoModule
import globals
from hopfalgebra import HopfAlgebraParse, createPolynomialHopfAlgebra
from resolution import resolution, Resolution
from steenrod import a_0_dual_temp, create_a_n_dual, create_real_a_n_dual
from visualize import show, visualize_dots, visualize_structure_lines


def display_and_print_resolution(res: Resolution):
    print(res)
    visualize_dots(res)
    visualize_structure_lines(res)
    show()


def coalgebra_resolution(filename: str):
    print("Parsing Coalgebra file")
    print()
    A, _ = CoAlgebra.parse(filename)

    M = CoModule.fp_module(A)
    
    print("Creating a resolution")
    print()
    res = resolution(M)
    print()
    return res


def generated_poly_coalg_resolution(filename: str):
    print("Parsing hopfalgebra file")
    print()
    f,g,r,c = HopfAlgebraParse(filename)
    globals.FIELD = f

    print("Generating polynomial algebra")
    print()
    A = createPolynomialHopfAlgebra(f, g, c, r, globals.ELEMENT_LIMIT)

    M = CoModule.fp_module(A)
    print("Creating a resolution")
    print()
    res = resolution(M)
    return res

def specific_comodule_resolution(coalgebra_file: str, comodule_file: str):
    print("Parsing Coalgebra file")
    print()
    A, names_to_basis = CoAlgebra.parse(coalgebra_file)

    print("Parsing Comodule file")
    print()
    M = CoModule.parse(comodule_file, A, names_to_basis)
    
    print("Creating a resolution")
    print()
    res = resolution(M)
    print()
    return res

def export_hopfalgbera(import_filename: str, export_filename: str):
    print("Parsing hopfalgebra file")
    print()
    f,g,r,c = HopfAlgebraParse(import_filename)
    globals.FIELD = f

    print("Generating polynomial algebra")
    print()
    A = createPolynomialHopfAlgebra(f, g, c, r, globals.ELEMENT_LIMIT)

    A.serialize(export_filename)


if __name__ == "__main__":    
    # Here you should give filenames (with correct folder location) of you coalgebra / comodule definitions

    # # Do a resolution over F_p with a generated polynomial HopfAlgebra,
    # # No explicit comodule needs to be given
    # res = generated_poly_coalg_resolution("./examples/generating/gen_A(1).txt")

    # Do a resolution over F_p with a premade Coalgebra,
    # No explicit comodule needs to be given
    res = coalgebra_resolution("./examples/coalgebra/A(2).txt")

    # Do a resolution over a module M with a premade Coalgebra,
    # No explicit comodule needs to be given
    # res = specific_comodule_resolution("./examples/coalgebra/A(1).txt", "./examples/comodule/f2_mod.txt")


    # Display and print the actual resolution
    display_and_print_resolution(res)


    # # Export a hopfalgebra generator to a coalgebra file
    # export_hopfalgbera("./examples/generating/gen_A(1).txt","./examples/coalgebra/A(1).txt")
