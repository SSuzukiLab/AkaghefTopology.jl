"""Canonical Z2 edge weights, stored as exactly representable 0/1 values."""
function _spin_edge_weights(v::VirtualLink)
    v.weight_mode==:z2 || throw(ArgumentError("spin weight APIs require :z2 mode"))
    v.is_weighted || throw(ArgumentError("no assigned Z2 weights"))
    values=[e.weight for e in real_edge_table(v)]
    all(x -> x in (0,1),values) || throw(ArgumentError("stored Z2 weights must be canonical 0/1"))
    Int.(values)
end

function _spin_weight_system(v::VirtualLink;boundary=nothing)
    v.weight_mode==:z2 || throw(ArgumentError("spin weight APIs require :z2 mode"))
    _validate_gauss(v.gauss,v.orientation)
    disks=disk_table(v)
    any(!,disks.isCycle) && throw(ArgumentError("spin equations currently require closed disk circuits; relative spin boundary conventions are not defined"))
    boundary===nothing || throw(ArgumentError("spin boundary targets are not supported"))
    cp2=BigInt.(disks.cp2)
    all(iseven,cp2) || throw(ArgumentError("spin obstruction requires even doubled Euler cochain on every disk"))
    A=Int.(mod.(disks.delta,2)); b=Int.(mod.(-div.(cp2,2),2))
    A,b
end

"""[Mih23, def:spin normal o-graph]: delta*w = 1 + Ndots/2 (mod 2).
Canonical paper-db Work: references/items/article/w-f2194518-7624-4b07-9eb5-3e04e0a7d036,
locator `m23-blk-0020` (source condition S).
Divide the even doubled Euler cochain BEFORE reduction mod 2.
This is a spin equation, not an integral framing test.
"""
function spin_weight_equation(v::VirtualLink;boundary=nothing)
    A,b=_spin_weight_system(v;boundary=boundary)
    w=_spin_edge_weights(v)
    length(w)==size(A,2) || throw(DimensionMismatch("one Z2 weight per real edge required"))
    (mode=:z2,modulus=2,delta=A,rhs=b,edge_weights=w,residual=mod.(A*w-b,2))
end
function validate_spin_weight(v::VirtualLink;boundary=nothing)
    all(iszero,spin_weight_equation(v;boundary=boundary).residual) || throw(ArgumentError("weights do not satisfy the spin equation"))
    true
end
"""Check closed spin weight equations only; this does not certify C1-C3 normality."""
function is_spin_weight(v::VirtualLink)
    any(_is_cut_component,v.gauss) && return false
    try
        validate_spin_weight(v)
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

"""Solve closed disk spin equations. Assignment preserves Z2 mode; never selects an integer lift."""
function solve_spin_weights(v::VirtualLink;assign=false,boundary=nothing)
    assign isa Bool || throw(ArgumentError("assign must be Boolean"))
    A,b=_spin_weight_system(v;boundary=boundary)
    result=solve_spin_weight_system(A,b)
    assignment=:not_requested
    if assign
        if result.particular===nothing
            assignment=:not_spin_feasible
        else
            target=deepcopy(v)
            set_weight!(target,result.particular)
            validate_spin_weight(target)
            v.weights=target.weights; v.is_weighted=true
            assignment=:assigned
        end
    end
    merge(result,(assignment_status=assignment,))
end
