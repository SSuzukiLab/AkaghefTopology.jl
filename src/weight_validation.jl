"""Validation and equation witnesses for integral normal o-graph weights.

The MATLAB source stores weights on `REdgeTable` rows, i.e. on real arcs
between real crossings.  The virtual Gauss presentation may contain virtual
crossing visits in between, so all public checks below use the stable
`real_edge_table` traversal order rather than raw Gauss positions.
"""

function _integer_edge_weights(v::VirtualLink)
    v.is_weighted || throw(ArgumentError("VirtualLink has no assigned integral weights"))
    edges=real_edge_table(v)
    values=Float64[edge.weight for edge in edges]
    any(value -> !isfinite(value),values) &&
        throw(ArgumentError("every real arc must have a finite weight"))
    all(isinteger,values) ||
        throw(ArgumentError("integral normal o-graph weights must be integer-valued"))
    BigInt[BigInt(value) for value in values]
end

"""Return the integer weight equation and its current residual.

The convention follows MST2 and MATLAB `calcWeight`: if `delta` is the
integer disk-boundary matrix and `cp2` is the doubled Euler cochain, a framed
integral normal o-graph satisfies

```text
2 * delta * edge_weights + cp2 == 0
```

`edge_weights` is ordered exactly as `real_edge_table(v)` and as the columns of
`disk_table(v).delta`.
"""
function weight_equation(v::VirtualLink)
    edge_weights=_integer_edge_weights(v)
    disks=disk_table(v)
    length(edge_weights)==size(disks.delta,2) ||
        throw(DimensionMismatch("weight count does not match real-edge count"))
    residual=2 .* BigInt.(disks.delta)*edge_weights+BigInt.(disks.cp2)
    (
        edge_weights=edge_weights,
        delta=disks.delta,
        cp2=disks.cp2,
        residual=residual,
    )
end

"""Return the integer residual of the framed weight equation."""
weight_residual(v::VirtualLink)=weight_equation(v).residual

"""Validate an integral normal o-graph weight assignment.

With `framed=false` this checks the integral edge-weight data contract.  With
`framed=true` it additionally requires the MST2 framing equation
`2*delta*w+cp2 == 0`.
"""
function validate_weight(v::VirtualLink; framed=false)
    equation=weight_equation(v)
    framed && any(!=(0),equation.residual) &&
        throw(ArgumentError("weights do not satisfy the framed integral equation; residual=$(equation.residual)"))
    true
end

"""Return whether the current assigned weights define a framed integral graph."""
function is_framed_weight(v::VirtualLink)
    try
        validate_weight(v;framed=true)
        true
    catch error
        error isa InterruptException && rethrow()
        false
    end
end

"""Return exact `dim_Q ker(delta) - rank_Q(H)` for a closed normal o-graph.

Checks `delta*H=0` and uses exact integer Smith ranks. This Betti number is
independent of Euler-cochain feasibility and is not a count of framings.
"""
function first_betti_number(v::VirtualLink)
    is_closed(v)[1] || throw(ArgumentError("first Betti number requires a closed normal o-graph"))
    delta=BigInt.(disk_table(v).delta)
    H=BigInt.(hmove_matrix(v))
    all(iszero,delta*H) || throw(ArgumentError("disk and H incidence do not form a cochain complex"))
    exact_rank(A)=begin
        S=AbstractAlgebra.snf(AbstractAlgebra.matrix(AbstractAlgebra.ZZ,A))
        count(i -> !iszero(S[i,i]),1:min(size(S)...))
    end
    size(delta,2)-exact_rank(delta)-exact_rank(H)
end
