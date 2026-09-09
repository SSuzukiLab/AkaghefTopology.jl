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

"""Return Z2 edge weights together with the measured doubled spin boundary data.

For every disk row, `measured4` is the exact residue
`2*delta*z - Ndots (mod 4)`.  Open rows have endpoint identities in the same
directed `(boundary_start, boundary_end)` convention used by integer relative
weight data.  No boundary target is inferred from this measurement.
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

"""[Mih23, def:spin normal o-graph]: delta*w = 1 + Ndots/2 (mod 2).
Canonical paper-db Work: references/items/article/w-f2194518-7624-4b07-9eb5-3e04e0a7d036,
locator `m23-blk-0020` (source condition S).
Divide the even doubled Euler cochain BEFORE reduction mod 2.
This is a spin equation, not an integral framing test.
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
function validate_spin_weight(v::VirtualLink;boundary=nothing)
    all(iszero,spin_weight_equation(v;boundary=boundary).residual) || throw(ArgumentError("weights do not satisfy the spin equation"))
    true
end
"""Check spin equations; cut inputs require the explicit relative target map."""
function is_spin_weight(v::VirtualLink;boundary=nothing)
    any(_is_cut_component,v.gauss) && boundary===nothing && return false
    try
        validate_spin_weight(v;boundary=boundary)
    catch error
        error isa InterruptException && rethrow()
        false
    end
end

"""Exact elimination over F2, with full affine kernel and a row witness on failure."""
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

"""Solve spin disk equations, including explicit relative targets on cut paths.

Assignment preserves Z2 mode; never selects an integer lift.  For an open disk
row, `boundary` must provide a target in Z4 keyed by its directed endpoint
identity (or a row-aligned vector).
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
