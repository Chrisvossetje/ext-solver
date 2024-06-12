from resolution import Resolution
import matplotlib.pyplot as plt
import globals

def visualize_dots(resolution: Resolution):
    xs = []
    ys = []
    for y, arr in enumerate(resolution.grading()):
        for x in arr:
            xs.append(x[0])
            ys.append(y)
    
    plt.scatter(xs,ys)

def visualize_structure_lines(resolution: Resolution):
    COLORS = ['r-', 'b-', 'g:']
    for n in range(2, len(resolution.morphisms)):
        s = n-2
        M = resolution.morphisms[n]

        lines = M.structure_lines()

        for l in lines:
            dom = l[0]
            codom = l[1]
            prim = l[2]
            x_source = M.domain.basis[dom[0]][dom[1]].grading[0] - s
            x_target = M.codomain.basis[codom[0]][codom[1]].grading[0] - s - 1
            plt.plot((x_source,x_target), (s,s+1), COLORS[prim % len(COLORS)])


def show():
    plt.xlabel("t-s")
    plt.ylabel("s")
    plt.xticks(range(0, globals.GRADE_LIMIT[0]+1, 5))
    plt.yticks(range(0, globals.FILTRATION_MAX+1, 2))
    plt.show()