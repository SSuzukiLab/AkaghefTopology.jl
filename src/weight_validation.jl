"""Validation and equation witnesses for integral normal o-graph weights.

The MATLAB source stores weights on `REdgeTable` rows, i.e. on real arcs
between real crossings.  The virtual Gauss presentation may contain virtual
crossing visits in between, so all public checks below use the stable
`real_edge_table` traversal order rather than raw Gauss positions.
"""

function _integer_edge_weights(v::VirtualLink)
    v.weight_mode==:integer || throw(ArgumentError("integral weights require integer mode"))
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
function _boundary_endpoint_identity(endpoint)
    endpoint===nothing && return nothing
    endpoint isa NamedTuple && hasproperty(endpoint,:strand) &&
        hasproperty(endpoint,:side) && hasproperty(endpoint,:lane) ||
        throw(ArgumentError("boundary endpoint metadata must expose strand, side, and lane"))
    (strand=Int(endpoint.strand),side=Symbol(endpoint.side),lane=Int(endpoint.lane))
end

function _boundary_identity(disks,row)
    (_boundary_endpoint_identity(disks.boundary_start[row]),
     _boundary_endpoint_identity(disks.boundary_end[row]))
end

function _boundary_open_rows(disks)
    findall(value -> !value,Bool.(disks.isCycle))
end

function _boundary_target_value(value,row)
    value isa Bool && throw(ArgumentError("boundary target for row $row must be an exact integer"))
    value isa Integer || throw(ArgumentError("boundary target for row $row must be an exact integer"))
    BigInt(value)
end

"""Resolve explicit targets for the open disk paths of `disks`.

The preferred mapping is a dictionary keyed by `(boundary_start, boundary_end)`;
the two endpoint values are the stable `(strand, side, lane)` identities exposed
by `disk_table` (extra endpoint metadata fields are ignored).  A row-aligned vector is accepted as a convenience for a
single exact `disk_table` result: closed rows must contain `nothing` or zero,
and every open row must contain an integer target.
"""
function _boundary_targets(disks,boundary)
    rows=size(disks.delta,1)
    length(disks.isCycle)==rows || throw(DimensionMismatch("disk cycle metadata is not row-aligned"))
    length(disks.boundary_start)==rows && length(disks.boundary_end)==rows ||
        throw(DimensionMismatch("disk boundary metadata is not row-aligned"))
    open_rows=_boundary_open_rows(disks)
    isempty(open_rows) && begin
        if boundary isa AbstractVector
            length(boundary)==rows || throw(DimensionMismatch("row-aligned boundary targets are required"))
            for row in 1:rows
                boundary[row]===nothing || _boundary_target_value(boundary[row],row)==0 ||
                    throw(ArgumentError("closed disk row $row cannot have a nonzero boundary target"))
            end
        elseif boundary isa AbstractDict
            isempty(boundary) || throw(ArgumentError("boundary mapping contains no open disk identity"))
        elseif boundary!==nothing
            throw(ArgumentError("boundary targets must be a row-aligned vector or dictionary"))
        end
        return zeros(BigInt,rows)
    end
    boundary===nothing && throw(ArgumentError("explicit boundary targets are required for open disk paths"))
    targets=zeros(BigInt,rows)
    seen=falses(rows)
    if boundary isa AbstractVector
        length(boundary)==rows || throw(DimensionMismatch("row-aligned boundary targets are required"))
        for row in 1:rows
            if row in open_rows
                boundary[row]===nothing && throw(ArgumentError("missing boundary target for open disk row $row"))
                targets[row]=_boundary_target_value(boundary[row],row); seen[row]=true
            else
                boundary[row]===nothing || _boundary_target_value(boundary[row],row)==0 ||
                    throw(ArgumentError("closed disk row $row cannot have a nonzero boundary target"))
            end
        end
    elseif boundary isa AbstractDict
        for (key,value) in pairs(boundary)
            row = if key isa Integer && !(key isa Bool)
                1 <= key <= rows || throw(BoundsError(disks.delta,(key,)))
                Int(key)
            elseif key isa Tuple && length(key)==2
                key_identity=(_boundary_endpoint_identity(key[1]),
                              _boundary_endpoint_identity(key[2]))
                matches=[row for row in open_rows if _boundary_identity(disks,row)==key_identity]
                length(matches)==1 || throw(ArgumentError("boundary identity does not identify exactly one open disk path"))
                only(matches)
            else
                throw(ArgumentError("boundary mapping keys must be disk rows or (boundary_start,boundary_end) identities"))
            end
            row in open_rows || throw(ArgumentError("boundary mapping refers to a closed disk row"))
            seen[row] && throw(ArgumentError("duplicate boundary target for open disk row $row"))
            targets[row]=_boundary_target_value(value,row); seen[row]=true
        end
    else
        throw(ArgumentError("boundary targets must be a row-aligned vector or dictionary"))
    end
    all(seen[open_rows]) || throw(ArgumentError("boundary targets must cover every open disk path"))
    BigInt[targets[row] for row in 1:rows]
end

"""Return exact passage contributions, without imposing boundary targets.

`measured[row]` is `2*delta[row,:]*w - Ndots[row]`.  For a closed disk row,
the framed residual is `measured[row] + 2`; for an open row the caller must
provide the relative boundary target against which `measured[row]` is compared.
"""
function boundary_weight_data(v::VirtualLink)
    edge_weights=_integer_edge_weights(v)
    disks=disk_table(v)
    length(edge_weights)==size(disks.delta,2) ||
        throw(DimensionMismatch("weight count does not match real-edge count"))
    delta=BigInt.(disks.delta)
    ndots=BigInt.(disks.Ndots)
    measured=2 .* delta*edge_weights-ndots
    closed_residual=measured+2 .* BigInt.(disks.isCycle)
    (
        edge_weights=edge_weights,
        delta=delta,
        Ndots=ndots,
        isCycle=copy(disks.isCycle),
        cp2=BigInt.(disks.cp2),
        measured=measured,
        closed_residual=closed_residual,
        boundary_start=copy(disks.boundary_start),
        boundary_end=copy(disks.boundary_end),
        boundary_identities=[_boundary_identity(disks,row) for row in axes(disks.delta,1)],
    )
end

function weight_equation(v::VirtualLink;boundary=nothing)
    data=boundary_weight_data(v)
    disks=disk_table(v)
    targets=_boundary_targets(disks,boundary)
    # Sign convention: cp2[i] = 2 - Ndots[i] = 2*c_P(MST2), so delta*w = -c_P(MST2).
    # MST3 §3.4 defines c_P with the opposite sign: c_P(MST3) = -c_P(MST2).
    # Thus its delta*w = c_P(MST3) is the SAME equation; do not negate w or delta.
    # Sources: MST2 §4.3 / Def.4.4; MST3 §§3.1,3.4. calculate_weight uses RHS -cp2.
    residual=data.measured+2 .* BigInt.(disks.isCycle)-targets
    (
        edge_weights=data.edge_weights,
        delta=disks.delta,
        cp2=disks.cp2,
        boundary_targets=targets,
        boundary_start=data.boundary_start,
        boundary_end=data.boundary_end,
        residual=residual,
    )
end

"""Return the integer residual of the framed weight equation."""
weight_residual(v::VirtualLink;boundary=nothing)=weight_equation(v;boundary=boundary).residual

"""Validate an integral normal o-graph weight assignment.

With `framed=false` this checks the integral edge-weight data contract.  With
`framed=true` it additionally requires the MST2 framing equation
`2*delta*w+cp2 == 0`.
"""
function validate_weight(v::VirtualLink; framed=false,boundary=nothing)
    if v.weight_mode==:z2
        framed && throw(ArgumentError("Z2 weights cannot be validated as integral framed weights; use validate_spin_weight"))
        _spin_edge_weights(v); return true
    end
    framed || (_integer_edge_weights(v); return true)
    equation=weight_equation(v;boundary=boundary)
    any(!=(0),equation.residual) &&
        throw(ArgumentError("weights do not satisfy the framed integral equation; residual=$(equation.residual)"))
    true
end

"""Return whether the current assigned weights define a closed framed integral graph.

Relative boundary equations are checked by `validate_weight`; a cut
presentation is never classified as a closed framed graph by this predicate.
"""
function is_framed_weight(v::VirtualLink;boundary=nothing)
    any(_is_cut_component, v.gauss) && return false
    try
        validate_weight(v;framed=true,boundary=boundary)
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
