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

"""
    boundary_weight_data(v::VirtualLink) -> NamedTuple

Return the measured integral weight data on the disk circuits and open paths.

# Arguments
Checked: W1 of [weight transport](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility);
`ArgumentError` identifies absent weights, wrong mode or invalid entries.
An edge/disk column-count mismatch raises `DimensionMismatch`.
The caller supplies coherent diagram data as required by [`disk_table`](@ref);
its errors propagate. No boundary targets are required for measurement.

# Returns
A `NamedTuple` with exact `edge_weights`, `delta`, `Ndots`, `cp2`, `measured`
and `closed_residual`; `isCycle`, `boundary_start`, `boundary_end` and
`boundary_identities` locate the rows in the relative-target model of §12.1.

# Guarantees
1. The returned measurements are the exact integral passage quantities of
   the [relative weight model](../specs/virtual-link-moves.md#121-supported-weighted-moves),
   using the disk table's edge and row coordinates. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Measured values do not choose independent boundary conditions. The returned
closed residual is not a framing certificate for open rows or a closedness check.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.1 — relative measurements and directed boundary targets.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); calculate_weight(v; assign=true);

julia> all(iszero, boundary_weight_data(v).closed_residual)
true
```

Role: invariant
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

"""
    weight_equation(v::VirtualLink; boundary=nothing) -> NamedTuple

Return the exact residual data for the integral o-graph weight equation.

# Arguments
`v` must have assigned finite integer real-arc weights in integer mode;
violations raise `ArgumentError` naming the mode or weight defect. Disk and
edge counts must agree or raise `DimensionMismatch`. Open disk paths require
explicit integer boundary targets accepted by the relative-boundary API;
missing, duplicate or invalid targets raise `ArgumentError`, a row-count
mismatch raises `DimensionMismatch`, and an out-of-range row raises
`BoundsError`. Diagram and disk-construction errors propagate. Geometric
framing interpretation additionally requires the hypotheses of
[weight transport §12](../specs/virtual-link-moves.md#12-weight-transport)
and the [closed-normal conditions](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3),
which this calculation does not certify.

# Returns
A `NamedTuple` containing `edge_weights`, disk-boundary matrix `delta`, doubled
Euler data `cp2`, explicit/default `boundary_targets`, directed endpoint
metadata `boundary_start` and `boundary_end`, and exact `BigInt` `residual`.

# Guarantees
1. The residual measures the assigned weights against the integral disk
   equation and supplied relative targets. It is computed with unbounded
   integer arithmetic from the returned discrete data.
   [proved/agent: by construction]
2. Closed disk rows use their fixed closed targets; open targets are not
   inferred from the measured residual. [proved/agent: by construction]
3. The diagram and supplied boundary data are not modified.
   [proved/agent: by construction]

# Not guaranteed
Nonzero residuals are returned rather than rejected. Zero residual alone
does not prove closedness, choose or compare framing classes, or establish
homeomorphism of represented manifolds. This is not the Z2 spin equation.
The geometric interpretation of arbitrary relative targets is not certified.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.0 — exact integral disk equation and the limit of a successful equation check on nonclosed presentations.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); calculate_weight(v; assign=true); equation = weight_equation(v);

julia> all(iszero, equation.residual) && equation.boundary_targets == zeros(BigInt, length(equation.residual))
true
```

See also [`validate_weight`](@ref), [`boundary_weight_data`](@ref), [`solve_integer_weights`](@ref).

Role: invariant
"""
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

"""
    weight_residual(v::VirtualLink; boundary=nothing) -> Vector{BigInt}

Return the exact residual of the assigned integral disk equation.

# Arguments
Checked: W1–W2 and the exceptions of `weight_equation`, as specified in
[weight transport](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
Invalid weights or targets raise `ArgumentError`, length failures
`DimensionMismatch`, and invalid target row keys `BoundsError`; diagram errors
propagate. Caller: supply the geometric hypotheses for a framing interpretation.

# Returns
The `residual` vector of `weight_equation`, in disk-table row order.

# Guarantees
1. Entries are the exact closed or relative integer-equation residuals for
   the assigned weights and prescribed targets. [proved/agent: by construction]
2. The diagram and targets are not modified. [proved/agent: by construction]

# Not guaranteed
A returned vector need not vanish. Zero residual alone does not establish
closed-normality, circuit identification, or equality of framing classes.

# Sources
cf: [virtual-link moves §12](../specs/virtual-link-moves.md#12-weight-transport) — residual and relative target semantics.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); all(iszero, weight_residual(v))
true
```

Role: invariant
"""
weight_residual(v::VirtualLink;boundary=nothing)=weight_equation(v;boundary=boundary).residual

"""
    validate_weight(v::VirtualLink; framed=false, boundary=nothing) -> Bool

Validate assigned integer or Z2 weights, optionally requiring the integral framing equation.

# Arguments
Checked: W1 and, with `framed=true`, W2 and zero residual of
[weight validation](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
Z2 mode instead checks assigned canonical entries and rejects `framed=true`.
Invalid mode, weights or residual raise `ArgumentError`; framing target
lengths and row keys can raise `DimensionMismatch` or `BoundsError`.
Caller: `framed` is Boolean; coherent diagram data and the geometric
closed-normal/circuit hypotheses are required for interpreting a framing.
Without a framing check, `boundary` is ignored.

# Returns
`true` on successful validation; failure raises rather than returning false.

# Guarantees
1. Basic validation certifies the applicable stored-weight domain; integral
   framed validation additionally certifies its disk-equation residual is zero.
   [proved/agent: by construction]
2. The diagram and boundary data are not modified. [proved/agent: by construction]

# Not guaranteed
Basic Z2 validation does not check the spin equation; use
`validate_spin_weight`. Integral validation does not certify every
closed-normal hypothesis or identify a geometric framing class.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — validation and interpretation boundaries.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]);

julia> validate_weight(v; framed=true)
true
```

Role: condition
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

"""
    is_framed_weight(v::VirtualLink; boundary=nothing) -> Bool

Test whether an uncut presentation passes the integral framed-weight equation.

# Arguments
Use W1–W2 of [weight transport](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
Checked: cut components give `false`; failures from integral validation give
`false`, except `InterruptException`, which propagates. The caller supplies
closed-normality and the source circuit correspondence for a geometric framing
interpretation; this predicate does not certify them.

# Returns
`true` exactly when no cut component is found and `validate_weight` succeeds
with `framed=true`; otherwise `false` for caught validation failures.

# Guarantees
1. A true result certifies the implemented integral equation on the supplied
   presentation, with the meaning and limits of §12. [proved/agent: by construction]
2. The diagram and boundary arguments are not modified. [proved/agent: by construction]

# Not guaranteed
This does not certify closed-normality or identify a framing class. Failure
does not distinguish a bad weight assignment from malformed diagram data.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12 — integral validation and geometric interpretation boundary.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); calculate_weight(v; assign=true);

julia> is_framed_weight(v)
true
```

Role: condition
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

"""
    first_betti_number(v::VirtualLink) -> Int

Return the rational cohomology dimension of the closed o-graph's disk/H complex.

# Arguments
Checked: W3 of [weight transport](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
`ArgumentError` reports a nonclosed input or failure of the cochain-complex
identity. Disk, H-incidence and exact-matrix errors propagate.
The caller supplies coherent diagram data and the geometric identification of these coordinates when
interpreting the answer as the first Betti number of a space.

# Returns
An `Int` giving the rational dimension of the kernel of the disk map modulo
the image of the H map.

# Guarantees
1. Under the checked complex identity, the dimension is exact and independent
   of solvability of the Euler-cochain equation. [proved/agent: by construction]
2. No weights or diagram data are changed. [proved/agent: by construction]

# Not guaranteed
This is not a count of framings. The rank computation alone does not identify
the complex with the cellular cohomology of a specified manifold.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.0 — exact disk/H ranks and their interpretation boundary.

# Examples
```jldoctest
julia> first_betti_number(lens_space_link(3,1))
0
```

Role: invariant
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
