using Topology
function run_C240610(; output_dir=nothing)
    d=URDiagram(); set_data_cvs!(d,1,[-5,1,-4,5,-3,4,-2,3,-1,2],ones(5))
    validate_data(d); simplify!(d)
    output_dir === nothing || (mkpath(output_dir);
        plot_svg(d,joinpath(output_dir,"C240610URD.svg"));
        plot_urd_structure_svg(d,joinpath(output_dir,"C240610URD_structure.svg")))
    d
end

"""Port the final `simplify2; trace2` cell with a common exact perturbation.

MATLAB's `setDataFromCVS` assigns a distinct `eps*` symbol to every edge and
`trace2` subsequently sends each one to zero.  The Julia symbolic surface does
not yet implement multivariate limits, so this executable oracle uses one
Laurent variable and takes its exact zero limit.  It is valid for this cell
only after the shared-perturbation result has been checked against the source
calculation.
"""
function run_C240610_trace()
    epsilon=laurent_variable(:epsilon)
    d=URDiagram()
    set_data_cvw!(d,1,[-5,1,-4,5,-3,4,-2,3,-1,2],fill(one(epsilon)+epsilon,5))
    validate_data(d); simplify2!(d)
    value=trace2(d)
    (diagram=d, rational=value, limit=limit_zero(value))
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_C240610(output_dir=joinpath(@__DIR__,"..","artifacts","plots"))
end
