"""Compare the MATLAB/Sage o-graph layout against the Julia port.

Run under the Sage venv interpreter, not the system python:

    /Applications/SageMath-10-6.app/Contents/Frameworks/Sage.framework/\
Versions/Current/venv/bin/python3 julia/scripts/sage/sage_layout_parity.py

It executes the function definitions of the canonical MATLAB template
`Manifold/VirtualLink/allocate_pos.py` directly, so the comparison uses the
same code MATLAB ships to Sage rather than a reimplementation.

Part A isolates the minimum-bend MILP (`allocate_positions` first
`MixedIntegerLinearProgram`) across MIP backends.  That model is degenerate:
several distinct optimal bending vectors share the same objective value, so the
backend decides which one is drawn.

Part B pins the bending numbers, which removes Part A's choice, and prints the
resulting polylines in a normalised form (consecutive duplicates and collinear
interior points removed) so they can be diffed against the output of
`julia/scripts/sage/julia_layout_parity.jl`.

See `julia/PLOT_PIPELINE.md` for the recorded results and their interpretation.
"""

import os
import sys

import sage.all
from sage.all import flatten
from sage.knots.link import Link
from sage.numerical.mip import MixedIntegerLinearProgram

TEMPLATE = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "Manifold", "VirtualLink", "allocate_pos.py"))

# `allocate_pos.py` is a MATLAB string template: the driver below the function
# definitions references placeholders that MATLAB substitutes.  Only the
# definitions are needed here, so the driver is cut at its first statement.
DRIVER_MARKER = "comp = LINK._isolated_components()"

# (name, pd code, bending numbers pinned for part B)
CASES = [
    ("hopf", [[3, 2, 4, 1], [2, 3, 1, 4]], [2, 2, 4, 0]),
    ("trefoil", [[1, 2, 4, 3], [3, 4, 6, 5], [5, 6, 2, 1]], [1, 3, 0, 2, 0, 2]),
    ("borromean", [[1, 5, 2, 8], [5, 9, 6, 12], [9, 1, 10, 4],
                   [3, 6, 4, 7], [7, 10, 8, 11], [11, 2, 12, 3]],
     [1, 1, 1, 1, 2, 2, 0, 0, 4, 0, 0, 0]),
    ("whitehead", [[6, 1, 7, 2], [10, 7, 5, 8], [4, 5, 1, 6],
                   [2, 10, 3, 9], [8, 4, 9, 3]],
     [-1, -1, -1, -1, 0, -3, 0, 0, 3, 0]),
]

SOLVERS = ["Coin", "GLPK", "PPL"]


def load_template_functions():
    """Execute the definitions of the canonical MATLAB template.

    `SageWrapper.exec` runs the template inside `sage.all.__dict__`, so the
    template body resolves bare Sage names such as `sign` and `flatten` from
    that namespace.  Seeding the execution namespace the same way keeps this
    harness faithful to the MATLAB call path.
    """
    with open(TEMPLATE) as handle:
        source = handle.read()
    cut = source.index(DRIVER_MARKER)
    # `PRINT_BENDING_NUMBERS` is a MATLAB-substituted placeholder.
    definitions = source[:cut].replace("PRINT_BENDING_NUMBERS", "False")
    namespace = dict(vars(sage.all))
    exec(compile(definitions, TEMPLATE, "exec"), namespace)
    return namespace["allocate_positions"]


def bending_numbers(pd, solver):
    """Re-solve the first MILP of `allocate_positions` with a chosen backend.

    The model is transcribed from the template rather than called through it
    because `allocate_positions` does not expose a backend argument.
    """
    link = Link(pd)
    regions = sorted(link.regions(), key=len)
    edges = sorted(set(flatten(link.pd_code())))
    program = MixedIntegerLinearProgram(maximization=False, solver=solver)
    variables = program.new_variable(nonnegative=True, integer=True)

    def source(edge):
        if edge > 0:
            return variables[2 * edges.index(edge)]
        return variables[2 * edges.index(-edge) + 1]

    count = len(regions)
    for index, region in enumerate(regions):
        capacity = len(region) - 4 if index < count - 1 else len(region) + 4
        program.add_constraint(
            sum(source(-edge) - source(edge) for edge in region) == capacity)
    program.set_objective(program.sum(variables.values()))
    objective = program.solve()
    values = program.get_values(variables)
    bending = [int(round(float(values[2 * i] - values[2 * i + 1])))
               for i in range(len(edges))]
    return bending, float(objective)


def normalise(points):
    """Drop repeated points and collinear interior points.

    Sage emits one segment per bend and repeats the shared endpoint; the Julia
    router emits corner points only.  Normalising both sides makes the two
    representations directly comparable.
    """
    cleaned = []
    for point in points:
        if not cleaned or cleaned[-1] != point:
            cleaned.append(point)
    if len(cleaned) < 3:
        return cleaned
    reduced = [cleaned[0]]
    for previous, corner, following in zip(cleaned, cleaned[1:], cleaned[2:]):
        cross = ((corner[0] - previous[0]) * (following[1] - corner[1])
                 - (corner[1] - previous[1]) * (following[0] - corner[0]))
        if cross != 0:
            reduced.append(corner)
    reduced.append(cleaned[-1])
    return reduced


def polylines(coords, edge_ids):
    """Concatenate the per-bend segments of each edge in emission order."""
    order = []
    joined = {}
    for edge_id, segment in zip(edge_ids, coords):
        if edge_id not in joined:
            order.append(edge_id)
            joined[edge_id] = []
        joined[edge_id].extend([float(p[0]), float(p[1])] for p in segment)
    return [(edge_id, normalise([tuple(p) for p in joined[edge_id]]))
            for edge_id in order]


def main():
    allocate_positions = load_template_functions()

    print("== Part A: minimum-bend MILP across MIP backends ==")
    print("Equal objectives with different vectors mean the model is degenerate")
    print("and the backend, not the algorithm, selects the drawn embedding.\n")
    for name, pd, _ in CASES:
        for solver in SOLVERS:
            try:
                bending, objective = bending_numbers(pd, solver)
                print(f"{name:10s} {solver:5s} bending={bending} objective={objective}")
            except Exception as error:  # a backend may be absent in this Sage
                print(f"{name:10s} {solver:5s} ERROR {type(error).__name__}: {error}")
        print()

    print("== Part B: layout with the bending numbers pinned ==")
    print("Compare against julia/scripts/sage/julia_layout_parity.jl.\n")
    for name, pd, pinned in CASES:
        try:
            coords, edge_ids, crossings = allocate_positions(Link(pd), pinned)
        except Exception as error:
            print(f"{name:10s} ERROR {type(error).__name__}: {error}\n")
            continue
        print(f"-- {name} (bending pinned to {pinned})")
        positions = {key: tuple(float(x) for x in value)
                     for key, value in crossings.items()}
        print(f"   crossings: {positions}")
        for edge_id, points in polylines(coords, edge_ids):
            print(f"   edge {edge_id}: {points}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
