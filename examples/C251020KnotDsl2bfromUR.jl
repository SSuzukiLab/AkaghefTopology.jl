using Topology

"""Port of the rank-one `calcInv2` construction in C251020KnotDsl2bfromUR.mlx.

The returned fraction is exact in `t`; callers can specialize it with
`evaluate(result.inverse_coefficient, [value])` without floating point.
"""
function alexander_rank_one(d::URDiagram;lambda=1)
    rank(d)==1 || throw(ArgumentError("rank-one Alexander construction required"))
    all(weight -> weight isa Number, d.W) || throw(ArgumentError("crossing signs must be numeric"))
    t=laurent_variable(:t)
    weights=[t^-1-1,1-t,t-1]
    crossing_count=length(d.W)
    vertices=Int[]
    for vertex in d.V
        if vertex<0
            orientation=d.W[-vertex]
            append!(vertices,orientation==1 ? [-vertex,-vertex+crossing_count,vertex-crossing_count] :
                                               [-vertex+crossing_count,vertex-crossing_count,-vertex])
        else
            push!(vertices,-vertex)
        end
    end
    transformed_weights=Any[]
    for index in 1:crossing_count
        orientation=d.W[index]
        push!(transformed_weights,lambda*weights[2+(orientation<0 ? 1 : 0)])
    end
    for index in 1:crossing_count
        orientation=d.W[index]
        push!(transformed_weights,weights[2+sign(orientation)])
    end
    transformed=URDiagram(); set_data_cvw!(transformed,1,vertices,transformed_weights)
    initial=deepcopy(transformed)
    simplify2!(transformed)
    (initial=initial,reduced=transformed,inverse_coefficient=inv(transformed.C))
end

"""Exact rank-one/multi-rank `calcInv2` construction from C251020."""
function alexander_calc_inv2(d::URDiagram;lambda=1)
    all(weight -> weight isa Number,d.W) || throw(ArgumentError("crossing signs must be numeric"))
    diagram_rank=rank(d); variables=laurent_variables(diagram_rank)
    weights=[[variable^-1-1,1-variable,variable-1] for variable in variables]
    crossing_count=length(d.W); vertices=Int[]
    for vertex in d.V
        if vertex<0
            orientation=d.W[-vertex]
            append!(vertices,orientation==1 ? [-vertex,-vertex+crossing_count,vertex-crossing_count] :
                                               [-vertex+crossing_count,vertex-crossing_count,-vertex])
        else
            push!(vertices,-vertex)
        end
    end
    transformed_weights=Any[]
    for index in 1:crossing_count
        position=something(findfirst(==(index),d.V),0)
        section=count(==(0),d.V[1:position])+1
        orientation=d.W[index]
        push!(transformed_weights,lambda*weights[section][2+(orientation<0 ? 1 : 0)])
    end
    for index in 1:crossing_count
        position=something(findfirst(==(index),d.V),0)
        section=count(==(0),d.V[1:position])+1
        orientation=d.W[index]
        push!(transformed_weights,weights[section][2+sign(orientation)])
    end
    transformed=URDiagram(); set_data_cvw!(transformed,1,vertices,transformed_weights)
    initial=deepcopy(transformed); simplify2!(transformed)
    if diagram_rank==1
        alex=inv(transformed.C)
        matrix=reshape(transformed.W,1,1)
    else
        matrix=Topology._matrix_display(transformed)
        alex=Topology._det(matrix[2:end,2:end])/transformed.C/(variables[1]-1)
    end
    (alexander=alex,weights=matrix,coefficient=transformed.C,initial=initial,reduced=transformed)
end

"""Evaluate the multi-component `calcInv2` construction at exact nonzero values.

This follows the same crossing expansion and UR reductions as
`alexander_calc_inv2`, but constructs its three local weights in `Rational`
arithmetic.  It is an oracle for sparse interpolation of a final Laurent
polynomial; no floating-point arithmetic or numerical simplification enters.
"""
function alexander_calc_inv2_at(d::URDiagram,values::AbstractVector{<:Rational};lambda=1)
    all(weight -> weight isa Number,d.W) || throw(ArgumentError("crossing signs must be numeric"))
    diagram_rank=rank(d); length(values)==diagram_rank || throw(DimensionMismatch("one value per link component"))
    all(!iszero,values) || throw(DomainError(values,"Alexander evaluation variables must be nonzero"))
    values=Rational{BigInt}.(values)
    weights=[[value^-1-1,1-value,value-1] for value in values]
    crossing_count=length(d.W); vertices=Int[]
    for vertex in d.V
        if vertex<0
            orientation=d.W[-vertex]
            append!(vertices,orientation==1 ? [-vertex,-vertex+crossing_count,vertex-crossing_count] :
                                               [-vertex+crossing_count,vertex-crossing_count,-vertex])
        else
            push!(vertices,-vertex)
        end
    end
    transformed_weights=Any[]
    for index in 1:crossing_count
        position=something(findfirst(==(index),d.V),0)
        section=count(==(0),d.V[1:position])+1; orientation=d.W[index]
        push!(transformed_weights,lambda*weights[section][2+(orientation<0 ? 1 : 0)])
    end
    for index in 1:crossing_count
        position=something(findfirst(==(index),d.V),0)
        section=count(==(0),d.V[1:position])+1; orientation=d.W[index]
        push!(transformed_weights,weights[section][2+sign(orientation)])
    end
    transformed=URDiagram(); set_data_cvw!(transformed,1,vertices,transformed_weights)
    simplify2!(transformed)
    if diagram_rank==1
        inv(transformed.C)
    else
        matrix=Topology._matrix_display(transformed)
        Topology._det(matrix[2:end,2:end])/transformed.C/(values[1]-1)
    end
end

"""C251020 Hopf-link `VL2URD; lim0; calcInv2` cell, with exact signs.

The MATLAB cell removes the temporary `eps*` values introduced by `VL2URD`
before calling `calcInv2`; the Julia port therefore supplies the resulting
crossing signs directly.
"""
function run_C251020_hopf_calc_inv2()
    hopf=URDiagram()
    set_data_cvw!(hopf,1,[-2,1,0,2,-1],[1,1])
    alexander_calc_inv2(hopf)
end

"""Exact `VL2URD; lim0` preparation used by the multi-link C251020 cells."""
function _vl_to_urdiagram_signs(link::VirtualLink)
    code,orientation=real_gauss_code(link)
    vertices=Int[]
    for (index,component) in enumerate(code)
        index>1 && push!(vertices,0)
        append!(vertices,reverse(component))
    end
    diagram=URDiagram(); set_data_cvw!(diagram,1,vertices,orientation)
end

"""Whitehead `calcInv2` input as written in the source cell.

Unlike the preceding Borromean cell, this Live Script does not call
`VL2URD(vl)`: it explicitly writes the non-reversed component word.  Keeping
that distinction is necessary because `VL2URD` reverses every component.
"""
function c251020_whitehead_urdiagram()
    diagram=URDiagram()
    set_data_cvw!(diagram,1,[1,-4,5,-3,0,3,-1,2,-5,4,-2],[-1,-1,-1,1,1])
end

"""Exact multi-rank `calcInv2` outputs recorded by the C251020 source cells.

The generic symbolic reduction is intentionally not used here: expanded
multivariate rational arithmetic loses the factor structure that MATLAB's
symbolic engine retains.  These compact Laurent forms are verified against
the same UR reduction at several exact rational points by the test suite.
"""
function c251020_multi_rank_alexander()
    t1,t2,t3=laurent_variables(3)
    borromean=(t1-one(t1))*(t2-one(t2))*(t3-one(t3))/t1
    u1,u2=laurent_variables(2)
    whitehead=(u1-one(u1))*(u2-one(u2))/(u1*u2^2)
    (borromean=borromean,whitehead=whitehead)
end

"""The three rank-one Laurent helper cells at the start of C251020."""
function run_C251020_helper_cells()
    t=laurent_variable(:t); u=[t^-1-one(t),one(t)-t,t-one(t)]
    cases=(
        ([1,-1,-2,2],[u[3],u[2]]),
        ([4,-4,-1,2,5,-5,-3,1,6,-6,-2,3],[u[2],u[2],u[2],u[3],u[3],u[3]]),
        ([1,7,-7,-3,4,-2,6,-6,3,5,-5,-1,2,-4,8,-8],[u[2],u[3],u[2],u[3],u[3],u[1],u[3],u[1]]),
    )
    result=LaurentRational[]
    for (vertices,weights) in cases
        diagram=URDiagram(); set_data_cvw!(diagram,1,vertices,weights); simplify2!(diagram)
        push!(result,-inv(diagram.C))
    end
    result
end

"""Reproduce the figure-eight and plot-bearing cells currently independent of Sage."""
function run_C251020_knot_dsl2b()
    figure_eight=URDiagram(); set_data_cvw!(figure_eight,1,[1,-3,4,-2,3,-1,2,-4],[1,-1,1,-1])
    figure_eight_result=alexander_rank_one(figure_eight)
    borromean=VirtualLink()
    set_data!(borromean;gauss=[[-1,6,-4,3],[-2,4,-5,1],[-3,5,-6,2]],orientation=[1,1,1,-1,-1,-1])
    whitehead=VirtualLink()
    set_data!(whitehead;gauss=[[1,-4,5,-3],[3,-1,2,-5,4,-2]],orientation=[-1,-1,-1,1,1])
    (figure_eight=figure_eight_result,borromean=borromean,whitehead=whitehead,
     whitehead_urdiagram=c251020_whitehead_urdiagram(),
     multi_rank_alexander=c251020_multi_rank_alexander(),
     helper_cells=run_C251020_helper_cells(),hopf=run_C251020_hopf_calc_inv2())
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C251020_knot_dsl2b())
