using Topology
function run_C240623(; output_dir=nothing)
    samples=[
        ([-1,5,-2,1,-3,2,-4,3,-5,4],[1,1,1,1,1]),
        ([1,-2,3,-1,2,-3],[1,1,1]),
        ([1,-2,-3,-4,4,3,2,-1],[1,-1,1,1]),
        ([6,-6,2,3,-3,-4,-5,5,4,1,-1,-2],[1,1,-1,1,1,-1]),
    ]
    diagrams=URDiagram[]
    for (V,W) in samples
        d=URDiagram(); set_data_cvw!(d,1,V,W); validate_data(d); push!(diagrams,d)
    end
    # This assignment is part of the fourth Live Script cell; it must happen
    # before simplification (MATLAB indexing is one-based, as is Julia's).
    diagrams[end].W[6]=sym(:a); diagrams[end].W[5]=sym(:b)
    simplify!(diagrams[end])
    output_dir === nothing || (mkpath(output_dir);
        plot_svg(diagrams[end],joinpath(output_dir,"C240623URD.svg"));
        plot_urd_structure_svg(diagrams[end],joinpath(output_dir,"C240623URD_structure.svg")))
    diagrams
end

"""Port the three `D.trace` cells as exact shared-epsilon limit checks.

`trace` invokes the rank-one determinant formula (`trace3` here).  See
`run_C240610_trace` for why one Laurent perturbation is used pending a general
multivariate symbolic-limit implementation.
"""
function run_C240623_traces()
    epsilon=laurent_variable(:epsilon)
    samples=[
        ([-1,5,-2,1,-3,2,-4,3,-5,4],[1,1,1,1,1]),
        ([1,-2,3,-1,2,-3],[1,1,1]),
        ([1,-2,-3,-4,4,3,2,-1],[1,-1,1,1]),
    ]
    results=NamedTuple[]
    for (V,W) in samples
        d=set_data_cvw!(URDiagram(),1,V,[weight+epsilon for weight in W])
        value=trace3(d)
        push!(results,(rational=value,limit=limit_zero(value)))
    end
    results
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_C240623(output_dir=joinpath(@__DIR__,"..","artifacts","plots"))
end
