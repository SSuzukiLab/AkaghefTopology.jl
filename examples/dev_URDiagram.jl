using Topology
function run_dev_urdiagram(; output_dir=nothing)
    invalid=URDiagram(); set_data_cvw!(invalid,1,[1,-2,-3,3,1,2],[1,-1,1])
    invalid_error = try validate_data(invalid); nothing catch e; e end
    d=URDiagram(); set_data_cvw!(d,1,[1,-2,-3,3,-1,2],[sym(:x),sym(:y),sym(:z)]); validate_data(d)
    swap!(d,[3,-3]); check_data(d,1/(1-sym(:z)),[1,-2,3,-3,-1,2],[sym(:x),sym(:y),sym(:z)/(1-sym(:z))])
    set_data_cvw!(d,1,[1,-2,-3,3,-1,2],[sym(:x),sym(:y),sym(:z)]); swap!(d,[-2,-3]); swap!(d,[-2,-3])
    output_dir === nothing || (mkpath(output_dir); plot_svg(d,joinpath(output_dir,"dev_URDiagram.svg")))
    d, invalid_error
end

"""Independent structural reduction cells from `dev_URDiagram.mlx`.

The Live Script displays these exploratory diagrams rather than asserting a
symbolic closed form.  This port preserves the source operation boundaries and
returns the vertex words, which remain stable before symbolic normalization.
"""
function run_dev_reduction_cells()
    reduct3_diagram=URDiagram()
    set_data_v!(reduct3_diagram,[1,2,-3,3,4,-1,-2,0,-4])
    reduct3!(reduct3_diagram)

    t=laurent_variable(:t)
    reduct5_diagram=URDiagram()
    set_data_cvw!(reduct5_diagram,1,[1,2,-3,3,-1,-2],[t,t^2,t^3])
    reduct5!(reduct5_diagram)
    after_first_reduct5=copy(reduct5_diagram.V)
    reduct5!(reduct5_diagram)

    simplify1_diagram=URDiagram()
    set_data_cvw!(simplify1_diagram,1,[1,-2,3,-1,2,-3],[-1,-1,-1])
    simplify1!(simplify1_diagram)
    (reduct3=copy(reduct3_diagram.V),reduct5_first=after_first_reduct5,
     reduct5_second=copy(reduct5_diagram.V),simplify1=copy(simplify1_diagram.V))
end

"""Port the `denom` loop using exact decimal-epsilon arithmetic.

The source supplies `1e-7` as a double.  Repeated UR reductions contain large
near-cancellations, so Float64 changes the N=4/5 rational reconstruction.
Representing that terminating decimal as a BigInt rational preserves the
source loop's intended algebra, then applies its `rat(...,1e-4)` acceptance
criterion only at the output boundary.
"""
function run_dev_denominator_loop(;maximum_order=15)
    maximum_order>=1 || throw(ArgumentError("maximum_order must be positive"))
    epsilon=big(1)//big(10)^7
    denominators=Float64[]
    traces=Any[]
    for order in 1:maximum_order
        vertices=fill(0,2order)
        for index in 1:order-1
            vertices[2index-1]=-index; vertices[2index+2]=index
        end
        vertices[[2,2order-1]]=[order,-order]
        diagram=URDiagram()
        set_data_cvw!(diagram,1,vertices,fill(big(1)+epsilon,order))
        simplify2!(diagram)
        trace=trace2(diagram); push!(traces,trace)
        approximation=rationalize(Int,Float64(trace),1e-4)
        # MATLAB records a non-unit large result as Inf.  Its implementation
        # is sign-sensitive; use magnitude so exact arithmetic cannot change
        # classification merely through a cancellation sign convention.
        push!(denominators,abs(numerator(approximation))!=1 && abs(Float64(trace))>100 ?
             Inf : Float64(denominator(approximation)))
    end
    (denominators=denominators,traces=traces)
end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_dev_urdiagram(output_dir=joinpath(@__DIR__,"..","artifacts","plots"))
end
