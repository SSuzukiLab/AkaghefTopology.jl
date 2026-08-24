using Topology
function run_URDforKnotCmpl(; output_dir=nothing)
    knots=Dict(
        "unknot"=>[1,-1],
        "trefoil"=>[-1,2,-3,1,-2,3],
        "figure-eight"=>[-1,2,-4,1,-3,4,-2,3],
    )
    diagrams=Dict{String,URDiagram}()
    for (name,k) in knots
        V=trans_knot_complement(k); W=repeat([1,-1,-1,1],outer=div(length(k),2))
        d=URDiagram(); set_data_cvw!(d,1,V,W); validate_data(d); diagrams[name]=d
    end
    output_dir === nothing || (mkpath(output_dir); plot_svg(diagrams["figure-eight"],joinpath(output_dir,"URDforKnotCmpl.svg")))
    diagrams
end

"""First manual D2 reduction sequence from `URDforKnotCmpl.mlx`.

This stops at the only self-contained manual state verified against the source.
The following Live Script cells are exploratory: `lim0` returns a value without
mutating `D2`, and later moves require the exact preceding cell state. The
intermediate vertex word is stored in the source Live Script output and is an
exact structural oracle.
"""
function run_URDforKnotCmpl_manual_prefix()
    d=URDiagram()
    set_data_cvw!(d,1,[-5,1,-4,4,3,2,-3,5,-2,-1],[-1,-1,-1,1,1])
    swap!(d,[-2,5]); reduct3!(d); swap!(d,[2,-2]); dilation!(d,2,6)
    swap!(d,[-2,4]); reduct3!(d)
    d
end

"""Replay the complete final exploratory D2 cell in `URDforKnotCmpl.mlx`.

The source cell is one ordered mutation sequence.  Its `lim0` call is omitted
from the mutations deliberately: MATLAB returns the limit but does not replace
`D2`, so the following swaps act on the same diagram.  Named snapshots make
that state boundary inspectable without treating stale embedded Live Script
output as an oracle.
"""
function run_URDforKnotCmpl_d2_exploratory_cell()
    d=run_URDforKnotCmpl_manual_prefix()
    prefix=(V=copy(d.V),C=d.C,W=copy(d.W))

    # `D2.lim0` is non-mutating in the MATLAB implementation.
    swap!(d,[2,-2]); dilation!(d,2,5); swap!(d,[1,5]); simplify2!(d)
    after_simplify=(V=copy(d.V),C=d.C,W=copy(d.W))

    swap!(d,[2,-3]); swap!(d,[2,5]); dilation!(d,3,5); swap!(d,[5,6])
    before_last_dilation=(V=copy(d.V),C=d.C,W=copy(d.W))
    dilation!(d,4,5); swap!(d,[1,5])
    (prefix=prefix,after_simplify=after_simplify,
     before_last_dilation=before_last_dilation,final=(V=copy(d.V),C=d.C,W=copy(d.W)))
end

"""Current-source numeric `simplify2; trace2` cell for the figure-eight case."""
function run_URDforKnotCmpl_figure_eight_trace()
    knot=[-1,2,-4,1,-3,4,-2,3]
    vertices=trans_knot_complement(knot)
    signs=repeat([1,-1,-1,1],outer=div(length(knot),2))
    diagram=URDiagram()
    # Match `setDataFromCVS(..., 1e-5)`: the current MATLAB implementation
    # runs this cell numerically, not with independent symbolic epsilons.
    set_data_cvw!(diagram,1,vertices,Float64.(signs).+1e-5)
    simplify2!(diagram)
    (diagram=diagram,trace=trace2(diagram))
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_URDforKnotCmpl(output_dir=joinpath(@__DIR__,"..","artifacts","plots"))
end
