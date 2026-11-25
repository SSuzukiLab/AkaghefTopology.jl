"""Canonical Z2 edge weights, stored as exactly representable 0/1 values."""
function _spin_edge_weights(v::VirtualLink)
    v.weight_mode==:z2 || throw(ArgumentError("spin weight APIs require :z2 mode"))
    v.is_weighted || throw(ArgumentError("no assigned Z2 weights"))
    values=[e.weight for e in real_edge_table(v)]
    all(x -> x in (0,1),values) || throw(ArgumentError("stored Z2 weights must be canonical 0/1"))
    Int.(values)
end

function _spin_target_value(value,row)
    value isa Bool && throw(ArgumentError("spin boundary target for row $row must be an exact integer in Z4"))
    value isa Integer || throw(ArgumentError("spin boundary target for row $row must be an exact integer in Z4"))
    mod(BigInt(value),4)
end

function _spin_boundary_targets(disks,boundary)
    rows=size(disks.delta,1)
    length(disks.isCycle)==rows || throw(DimensionMismatch("disk cycle metadata is not row-aligned"))
    length(disks.boundary_start)==rows && length(disks.boundary_end)==rows ||
        throw(DimensionMismatch("disk boundary metadata is not row-aligned"))
    open_rows=findall(!,Bool.(disks.isCycle))
    targets=fill(BigInt(2),rows)
    if isempty(open_rows)
        boundary===nothing && return targets
        if boundary isa AbstractVector
            length(boundary)==rows || throw(DimensionMismatch("row-aligned spin boundary targets are required"))
            for row in 1:rows
                boundary[row]===nothing && continue
                _spin_target_value(boundary[row],row)==2 ||
                    throw(ArgumentError("closed disk row $row has fixed spin target 2 in Z4"))
            end
            return targets
        elseif boundary isa AbstractDict
            isempty(boundary) || throw(ArgumentError("spin boundary mapping refers to a closed disk row"))
            return targets
        end
        throw(ArgumentError("spin boundary targets must be a row-aligned vector or dictionary"))
    end
    boundary===nothing && throw(ArgumentError("explicit Z4 spin boundary targets are required for open disk paths"))
    seen=falses(rows)
    if boundary isa AbstractVector
        length(boundary)==rows || throw(DimensionMismatch("row-aligned spin boundary targets are required"))
        for row in 1:rows
            if row in open_rows
                boundary[row]===nothing && throw(ArgumentError("missing spin boundary target for open disk row $row"))
                targets[row]=_spin_target_value(boundary[row],row); seen[row]=true
            elseif boundary[row]!==nothing
                _spin_target_value(boundary[row],row)==2 ||
                    throw(ArgumentError("closed disk row $row has fixed spin target 2 in Z4"))
            end
        end
    elseif boundary isa AbstractDict
        for (key,value) in pairs(boundary)
            row=if key isa Integer && !(key isa Bool)
                1<=key<=rows || throw(BoundsError(disks.delta,(key,)))
                Int(key)
            elseif key isa Tuple && length(key)==2
                identity=(_boundary_endpoint_identity(key[1]),_boundary_endpoint_identity(key[2]))
                matches=[candidate for candidate in open_rows
                         if _boundary_identity(disks,candidate)==identity]
                length(matches)==1 || throw(ArgumentError("spin boundary identity does not identify exactly one open disk path"))
                only(matches)
            else
                throw(ArgumentError("spin boundary mapping keys must be disk rows or (boundary_start,boundary_end) identities"))
            end
            row in open_rows || throw(ArgumentError("spin boundary mapping refers to a closed disk row"))
            seen[row] && throw(ArgumentError("duplicate spin boundary target for open disk row $row"))
            targets[row]=_spin_target_value(value,row); seen[row]=true
        end
    else
        throw(ArgumentError("spin boundary targets must be a row-aligned vector or dictionary"))
    end
    all(seen[open_rows]) || throw(ArgumentError("spin boundary targets must cover every open disk path"))
    targets
end

function _spin_weight_system(v::VirtualLink;boundary=nothing)
    v.weight_mode==:z2 || throw(ArgumentError("spin weight APIs require :z2 mode"))
    _validate_gauss(v.gauss,v.orientation)
    disks=disk_table(v)
    targets=_spin_boundary_targets(disks,boundary)
    ndots=BigInt.(disks.Ndots)
    parity=targets+ndots
    all(iseven,parity) || throw(ArgumentError("spin boundary target plus Ndots must be even before division by 2; parity=$(Int.(mod.(parity,2)))"))
    A=Int.(mod.(disks.delta,2))
    b=Int[mod(div(parity[row],2),2) for row in axes(disks.delta,1)]
    A,b,targets
end

"""
    spin_boundary_data(v::VirtualLink) -> NamedTuple

Measure the doubled spin boundary cochain of an assigned Z2-weighted presentation.

# Arguments
Checked: Z2 mode, assigned canonical weights, and alignment with disk columns;
violations raise `ArgumentError` or `DimensionMismatch`. Disk construction
errors propagate.
Caller: provide coherent Gauss data and interpret the
[relative boundary coordinates](../specs/spin-weights.md#relative-boundary-data)
using the intended directed circuit structure.

# Returns
A record with mode, modulus, edge weights, disk incidence, dot counts, cycle
flags, Euler data, measured residues, endpoints and directed boundary identities.

# Guarantees
1. The measurements are the doubled spin boundary cochain in the
   [repository coordinates](../specs/spin-weights.md#repository-coordinates),
   with values modulo four. [proved/agent: by construction]
2. The function chooses no prescribed boundary target and does not modify
   the input diagram. [proved/agent: by construction]

# Not guaranteed
Measurement alone does not establish the spin equation, closed-normality or
a geometric boundary trivialization. It is not a choice of boundary condition.

# Sources
cf: [spin specification](../specs/spin-weights.md#relative-boundary-data) — measured cochain and directed endpoints.

# Examples
```jldoctest
julia> v = VirtualLink(; weight_mode=:z2); set_data!(v; rgauss=[-2,-1,1,2,0], orientation=[1,1]); set_weight!(v, zeros(5));
 
julia> Tuple(findall(!, spin_boundary_data(v).isCycle))
(1, 3, 4)
```

Role: invariant
"""
function spin_boundary_data(v::VirtualLink)
    edge_weights=_spin_edge_weights(v)
    disks=disk_table(v)
    length(edge_weights)==size(disks.delta,2) ||
        throw(DimensionMismatch("weight count does not match real-edge count"))
    delta=BigInt.(disks.delta)
    ndots=BigInt.(disks.Ndots)
    measured4=mod.(2 .* delta*BigInt.(edge_weights)-ndots,4)
    (
        mode=:z2,
        modulus=4,
        edge_weights=edge_weights,
        delta=delta,
        Ndots=ndots,
        isCycle=copy(disks.isCycle),
        cp2=BigInt.(disks.cp2),
        measured4=measured4,
        boundary_start=copy(disks.boundary_start),
        boundary_end=copy(disks.boundary_end),
        boundary_identities=[_boundary_identity(disks,row) for row in axes(disks.delta,1)],
    )
end

"""
    spin_weight_equation(v::VirtualLink; boundary=nothing) -> NamedTuple

Return the assigned closed or relative spin equation and its exact residual.

# Arguments
Checked: the data and parity conditions S1, S3 and S4 of the
[spin specification](../specs/spin-weights.md#numbered-conditions), including
Gauss validity and weight alignment. Invalid values raise `ArgumentError`,
metadata or target lengths `DimensionMismatch`, and invalid row keys
`BoundsError`; other diagram errors propagate.
Caller: the closed geometric interpretation additionally requires the
closed-normal hypotheses and identification with the source circuits.

# Returns
For closed disk rows, a modulus-two record with `delta`, `rhs`,
`edge_weights` and `residual`. If any row is open, a modulus-four record
also includes `measured4`, `boundary_targets` and directed endpoints, as
specified under [equation records](../specs/spin-weights.md#equation-records).

# Guarantees
1. The residual is exactly that of the applicable
   [coordinate equation](../specs/spin-weights.md#definitions); zero is the
   algebraic spin condition on this presentation. [proved/agent: by construction]
2. Neither the diagram nor the supplied targets are modified.
   [proved/agent: by construction]

# Not guaranteed
Returning a record does not assert that its residual vanishes. The equation
does not certify closed-normality, identify equivalent spin structures, or
establish a geometric relative spin-extension theorem.

# Sources
cf:[[Mih23#^m23-blk-0020]] closed condition S; [spin specification](../specs/spin-weights.md#relative-boundary-data) — relative algebraic extension.

# Examples
```jldoctest
julia> v = VirtualLink(; weight_mode=:z2); set_data!(v; rgauss=[-2,-1,1,2,0], orientation=[1,1]); set_weight!(v, zeros(5));

julia> e = spin_weight_equation(v; boundary=[0,2,2,0]); (e.modulus, all(iszero, e.residual))
(4, true)
```

Role: invariant
"""
function spin_weight_equation(v::VirtualLink;boundary=nothing)
    A,b,targets=_spin_weight_system(v;boundary=boundary)
    w=_spin_edge_weights(v)
    length(w)==size(A,2) || throw(DimensionMismatch("one Z2 weight per real edge required"))
    disks=disk_table(v)
    if any(!,disks.isCycle)
        data=spin_boundary_data(v)
        return (mode=:z2,modulus=4,delta=data.delta,rhs=targets,edge_weights=w,
                measured4=data.measured4,boundary_targets=targets,
                boundary_start=data.boundary_start,boundary_end=data.boundary_end,
                residual=mod.(data.measured4-targets,4))
    end
    (mode=:z2,modulus=2,delta=A,rhs=b,edge_weights=w,residual=mod.(A*w-b,2))
end
"""
    validate_spin_weight(v::VirtualLink; boundary=nothing) -> Bool

Check the assigned Z2 weights against the closed or relative spin equation.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/spin-weights.md).

`v` and `boundary` must satisfy S1–S4 of the
[spin specification](../specs/spin-weights.md#numbered-conditions).
These conditions are checked by the supported equation interface.
`ArgumentError` identifies an
invalid mode, missing/noncanonical weights, missing or invalid boundary
targets, parity obstruction, or nonzero spin residual. Metadata/weight-count
mismatches raise `DimensionMismatch`; an out-of-range boundary row raises
`BoundsError`. Other diagram errors propagate. A closed-normal geometric
interpretation also requires N1, N2 and C1–C3 of
[the terminology specification §11](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3),
which this function does not fully certify.

# Returns
`true` when the supported spin equation has zero residual; failure raises
an exception instead of returning `false`.

# Guarantees
1. Success means every residual reported by `spin_weight_equation` vanishes.
   On closed normal data with the source's circuit interpretation, this checks
   [condition S](../specs/spin-weights.md#closed-spin-normal-o-graph)
   for the assigned weights. [proved/agent: by construction]
2. The diagram and supplied boundary targets are not modified.
   [proved/agent: by construction]

# Not guaranteed
This does not validate an integral framing, identify a spin structure up to
equivalence, or certify all closed-normal hypotheses. The open-path check
is the repository's relative boundary equation; the cited closed condition S
does not by itself prove a relative spin-extension theorem.

# Sources
cf:[[Mih23#^m23-blk-0020]] definition of spin normal o-graph, condition S — closed circuit condition; open relative boundary equations require additional interpretation.

# Examples
```jldoctest
julia> v = VirtualLink(; weight_mode=:z2); set_data!(v; rgauss=[1,2,-5,3,-2,-3,-4,4,-1,5,6,-7,7,-6], orientation=[1,1,1,-1,-1,-1,1]); solve_spin_weights(v; assign=true);

julia> validate_spin_weight(v)
true
```

Role: condition
"""
function validate_spin_weight(v::VirtualLink;boundary=nothing)
    all(iszero,spin_weight_equation(v;boundary=boundary).residual) || throw(ArgumentError("weights do not satisfy the spin equation"))
    true
end
"""
    is_spin_weight(v::VirtualLink; boundary=nothing) -> Bool

Test whether the assigned weights satisfy the supported closed or relative spin equation.

# Arguments
Checked: the conditions checked by `validate_spin_weight`; a cut component
without explicit `boundary` returns `false`. Validation failures return
`false`, while interruption propagates.
Caller: the geometric interpretation requires the closed-normal hypotheses
and circuit interpretation of the [spin specification](../specs/spin-weights.md#closed-spin-normal-o-graph).

# Returns
Whether validation succeeds, with the cut-input restriction above.

# Guarantees
1. A true result certifies the assigned equation of
   [S1–S5](../specs/spin-weights.md#numbered-conditions), within the applicable
   closed or relative case. [proved/agent: by construction]
2. The diagram and boundary data are not modified. [proved/agent: by construction]

# Not guaranteed
A false result does not distinguish invalid data from an equation obstruction.
The test does not identify a spin structure up to equivalence or establish a
relative spin-extension theorem for an open diagram.

# Sources
cf:[[Mih23#^m23-blk-0020]] condition S; see the closed and relative interpretation in the spin specification.

# Examples
```jldoctest
julia> v = VirtualLink(; weight_mode=:z2); set_data!(v; rgauss=[1,-1,0], orientation=[1]); set_weight!(v, [0,0,0]);

julia> is_spin_weight(v)
false
```

Role: condition
"""
function is_spin_weight(v::VirtualLink;boundary=nothing)
    any(_is_cut_component,v.gauss) && boundary===nothing && return false
    try
        validate_spin_weight(v;boundary=boundary)
    catch error
        error isa InterruptException && rethrow()
        false
    end
end

"""
    solve_spin_weight_system(A::AbstractMatrix, b::AbstractVector) -> NamedTuple

Solve an exact binary linear system and return its full affine solution space or obstruction.

# Arguments
Checked: one-based indexing, one right-hand side per row (`DimensionMismatch`),
and exact integer entries (`ArgumentError`). Entries, including booleans, are
interpreted modulo two as specified under
[binary solution records](../specs/spin-weights.md#binary-solution-records).
Caller: choose coordinates identifying this matrix with the intended spin equation.

# Returns
A versioned solution record with reduced `A`, `b`, `status`, optional
`particular`, column `kernel`, row certificate and optional `obstruction`.
`assignment_status` is `:not_requested`.

# Guarantees
1. `:spin_solution` provides a particular solution and a full binary kernel
   basis; `:inconsistent` provides a row-combination obstruction. The kernel
   is supplied in either case. [proved/agent: by construction]
2. Input arrays are not modified. [proved/agent: by construction]

# Not guaranteed
The matrix calculation alone does not identify a spin structure on a space,
select an integer lift, or quotient solutions by geometric equivalence.

# Sources
cf: [spin specification](../specs/spin-weights.md#binary-solution-records) — exact solution and witness semantics.

# Examples
```jldoctest
julia> r = solve_spin_weight_system([1 1 0; 0 1 1], [1,0]); (r.status, size(r.kernel), mod.(r.A*r.particular, 2) == r.b)
(:spin_solution, (3, 1), true)
```

Role: construct
"""
function solve_spin_weight_system(A::AbstractMatrix,b::AbstractVector)
    Base.require_one_based_indexing(A,b)
    size(A,1)==length(b) || throw(DimensionMismatch("one RHS per row required"))
    all(x -> x isa Integer,A) && all(x -> x isa Integer,b) || throw(ArgumentError("exact integer coefficients required"))
    A=Int.(mod.(A,2)); b=Int.(mod.(b,2)); m,n=size(A)
    R=hcat(A,b); U=Matrix{Int}(I,m,m); pivots=Int[]; r=0
    for col in 1:n
        row=findfirst(i -> R[i,col]==1,r+1:m)
        row===nothing && continue
        row+=r; r+=1
        R[[r,row],:]=R[[row,r],:]; U[[r,row],:]=U[[row,r],:]
        for i in 1:m
            i==r && continue
            if R[i,col]==1
                R[i,:]=mod.(R[i,:]+R[r,:],2); U[i,:]=mod.(U[i,:]+U[r,:],2)
            end
        end
        push!(pivots,col)
    end
    bad=findfirst(i -> all(iszero,R[i,1:n]) && R[i,n+1]==1,1:m)
    free=setdiff(collect(1:n),pivots); kernel=zeros(Int,n,length(free))
    for (j,col) in enumerate(free)
        kernel[col,j]=1
        for (i,pivot) in enumerate(pivots)
            kernel[pivot,j]=R[i,col]
        end
    end
    particular=bad===nothing ? zeros(Int,n) : nothing
    if particular!==nothing
        for (i,pivot) in enumerate(pivots)
            particular[pivot]=R[i,n+1]
        end
    end
    (schema="topology-spin-weight-solution-v1",mode=:z2,modulus=2,
     status=bad===nothing ? :spin_solution : :inconsistent,A=A,b=b,
     particular=particular,kernel=kernel,
     certificate=(rref=R,row_transform=U,pivots=pivots),
     obstruction=bad===nothing ? nothing : (row_combination=copy(U[bad,:]),rhs=1),
     assignment_status=:not_requested)
end

"""
    solve_spin_weights(v::VirtualLink; assign=false, boundary=nothing) -> NamedTuple

Solve the presentation's closed or relative binary spin equations.

# Arguments
Checked: `assign` must be Boolean; Z2 mode, Gauss validity, boundary target
coverage and parity satisfy the applicable
[S3–S4 conditions](../specs/spin-weights.md#numbered-conditions).
Violations raise `ArgumentError`, metadata length failures `DimensionMismatch`,
and invalid row keys `BoundsError`; other diagram errors propagate.
Caller: closed-normality and source-circuit correspondence are required for
the geometric spin interpretation. Prior weight assignment is unnecessary.

# Returns
The binary solution record with `boundary_targets` and `assignment_status`,
which is `:not_requested`, `:not_spin_feasible`, or `:assigned`.

# Guarantees
1. The solution space is exactly that of the
   [closed or relative coordinate equation](../specs/spin-weights.md#definitions).
   [proved/agent: by construction and exact binary solving]
2. With `assign=false` the diagram is unchanged. Successful requested
   assignment installs validated canonical Z2 weights and sets `is_weighted`;
   infeasibility leaves it unchanged. [proved/agent: by construction]

# Not guaranteed
The function does not choose an integer lift or certify all closed-normal
hypotheses. Relative equation solvability does not yet establish a geometric
boundary spin extension.

# Sources
cf:[[Mih23#^m23-blk-0020]] closed condition S; [relative model](../specs/spin-weights.md#relative-boundary-data).

# Examples
```jldoctest
julia> v = VirtualLink(; weight_mode=:z2); set_data!(v; rgauss=[-2,-1,1,2,0], orientation=[1,1]);

julia> r = solve_spin_weights(v; boundary=[0,2,2,0], assign=true); (r.status, r.assignment_status)
(:spin_solution, :assigned)
```

Role: convert
"""
function solve_spin_weights(v::VirtualLink;assign=false,boundary=nothing)
    assign isa Bool || throw(ArgumentError("assign must be Boolean"))
    A,b,targets=_spin_weight_system(v;boundary=boundary)
    result=solve_spin_weight_system(A,b)
    assignment=:not_requested
    if assign
        if result.particular===nothing
            assignment=:not_spin_feasible
        else
            target=deepcopy(v)
            set_weight!(target,result.particular)
            validate_spin_weight(target;boundary=boundary)
            v.weights=target.weights; v.is_weighted=true
            assignment=:assigned
        end
    end
    merge(result,(assignment_status=assignment,boundary_targets=targets))
end
