#!/usr/bin/env julia
#
# Julia counterpart of `sage_layout_parity.py`.
#
#     julia --project=julia julia/scripts/sage/julia_layout_parity.jl
#
# Part A prints this port's minimum-bend solution so it can be placed beside
# the Sage backends. Part B pins the bending numbers and prints the resulting
# polylines in the same normalised form the Python script emits, so the two
# outputs can be diffed directly.
#
# Recorded results and their interpretation: `julia/PLOT_PIPELINE.md`.

using Topology

"""Drop repeated points and collinear interior points.

Sage emits one segment per bend and repeats the shared endpoint, while
`_route_edges` emits corners only. Normalising both sides makes the two
representations comparable.
"""
function normalise(points)
    cleaned = ComplexF64[]
    for point in points
        (isempty(cleaned) || last(cleaned) != point) && push!(cleaned, point)
    end
    length(cleaned) < 3 && return cleaned
    reduced = ComplexF64[first(cleaned)]
    for index in 2:length(cleaned)-1
        previous, corner, following = cleaned[index-1], cleaned[index], cleaned[index+1]
        cross = real(corner-previous)*imag(following-corner) -
                imag(corner-previous)*real(following-corner)
        cross != 0 && push!(reduced, corner)
    end
    push!(reduced, last(cleaned))
    reduced
end

show_points(points) = "[" * join(("($(real(p)), $(imag(p)))" for p in points), ", ") * "]"

"""The part-A/B cases, built the way the Live Scripts build them.

The trefoil case of `sage_layout_parity.py` is absent on purpose: its PD matrix
does not survive a round trip through `set_data!`, so a bending vector indexed
to it cannot be reused here. See `PLOT_PIPELINE.md` IS4.
"""
function cases()
    hopf = VirtualLink()
    set_data!(hopf; gauss=[[1,-2],[-1,2]], orientation=[1,1])
    borromean = VirtualLink()
    set_data!(borromean; gauss=[[-1,6,-4,3],[-2,4,-5,1],[-3,5,-6,2]],
              orientation=[1,1,1,-1,-1,-1])
    whitehead = VirtualLink()
    set_data!(whitehead; gauss=[[1,-4,5,-3],[3,-1,2,-5,4,-2]],
              orientation=[-1,-1,-1,1,1])
    [("hopf", hopf, [2,2,4,0]),
     ("borromean", borromean, [1,1,1,1,2,2,0,0,4,0,0,0]),
     ("whitehead", whitehead, [-1,-1,-1,-1,0,-3,0,0,3,0])]
end

function main()
    println("== Part A: this port's minimum-bend solution (JuMP + GLPK) ==\n")
    for (name, link, _) in cases()
        pd = pd_code(link)
        println("$(rpad(name,10)) pd=$(pd)")
        println("$(rpad("",10)) bending=$(minimal_bending_numbers(pd))")
    end

    println("\n== Part B: layout with the bending numbers pinned ==")
    println("Compare against part B of sage_layout_parity.py.\n")
    for (name, link, pinned) in cases()
        layout = orthogonal_layout(link; bending_numbers=pinned)
        println("-- $name (bending pinned to $pinned)")
        println("   crossings: $([(real(p), imag(p)) for p in layout.crossings])")
        for edge in layout.edges
            println("   edge $(edge.id): $(show_points(normalise(edge.points)))")
        end
        println()
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
