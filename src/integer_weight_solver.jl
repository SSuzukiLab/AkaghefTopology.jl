import AbstractAlgebra

const INTEGER_WEIGHT_SOLUTION_SCHEMA="topology-integer-weight-solution-v1"
const INTEGER_WEIGHT_LIFT_SCHEMA="topology-integer-weight-lift-v1"
const _INTEGER_WEIGHT_SAFE_MAX=big(2)^53-1

_weight_bigint_matrix(M)=[BigInt(M[i,j]) for i in 1:size(M,1), j in 1:size(M,2)]
_weight_exact_matrix(M)=AbstractAlgebra.matrix(AbstractAlgebra.ZZ,M)
_weight_certificate_require(condition,message)=condition ||
    throw(ArgumentError("invalid integer weight certificate: $message"))

function _weight_integer_inputs(A::AbstractMatrix,b::AbstractVector)
    Base.require_one_based_indexing(A,b)
    size(A,1)==length(b) || throw(DimensionMismatch("one right-hand side per equation is required"))
    all(x -> x isa Integer,A) && all(x -> x isa Integer,b) ||
        throw(ArgumentError("exact integer coefficients are required; floating inputs are not rounded"))
    BigInt.(A),BigInt.(b)
end

"""
    solve_integer_weight_system(A, b)

Solve `A*x=b` over the integers with an exact Smith certificate. `status` is
`:integer_solution`, `:rational_only`, or `:inconsistent`. `particular` is a
`Vector{BigInt}` only in the first case. `kernel` has one full integer-lattice
basis vector per COLUMN, so every integer solution is `particular + kernel*z`
for a unique integer vector `z`. The kernel is also returned for infeasible
inhomogeneous systems. `rational_particular` distinguishes a rational-only
solution from rational inconsistency; it is never substituted for an integer
solution.

The certificate stores `U*A*V=S`, with unimodular `U,V`, and `U*b`. A nonzero
entry of `U*b` in a zero Smith row proves rational inconsistency. Otherwise a
failed diagonal divisibility proves integer infeasibility. These statements
concern this matrix, not every presentation of an underlying manifold.

Uses the installed AbstractAlgebra `snf_with_transform`; no floating optimizer,
kernel denominator clearing, or numerical rank is involved.
"""
function solve_integer_weight_system(A::AbstractMatrix,b::AbstractVector)
    A,b=_weight_integer_inputs(A,b)
    smith,left,right=AbstractAlgebra.snf_with_transform(_weight_exact_matrix(A))
    S=_weight_bigint_matrix(smith); U=_weight_bigint_matrix(left); V=_weight_bigint_matrix(right)
    m,n=size(A)
    # Fix row units, retaining the full exact transform, for positive divisors.
    for i in 1:min(m,n)
        if S[i,i]<0
            S[i,:].*=-1; U[i,:].*=-1
        end
    end
    r=count(i -> S[i,i]!=0,1:min(m,n))
    c=U*b
    kernel=copy(V[:,r+1:n])
    inconsistent=findfirst(i -> c[i]!=0,r+1:m)
    indivisible=findfirst(i -> mod(c[i],S[i,i])!=0,1:r)
    particular=nothing; rational_particular=nothing; obstruction=nothing
    if inconsistent!==nothing
        row=r+inconsistent
        status=:inconsistent
        obstruction=(kind=:rational_inconsistency,row=row,
            row_combination=copy(U[row,:]),lhs_coefficients=vec(transpose(U[row,:])*A),
            divisor=big(0),rhs_value=c[row],remainder=c[row])
    else
        y=zeros(Rational{BigInt},n)
        for i in 1:r
            y[i]=c[i]//S[i,i]
        end
        rational_particular=V*y
        if indivisible!==nothing
            row=indivisible
            status=:rational_only
            obstruction=(kind=:divisibility,row=row,
                row_combination=copy(U[row,:]),lhs_coefficients=vec(transpose(U[row,:])*A),
                divisor=S[row,row],rhs_value=c[row],remainder=mod(c[row],S[row,row]))
        else
            status=:integer_solution
            particular=BigInt[numerator(x) for x in rational_particular]
        end
    end
    result=(schema=INTEGER_WEIGHT_SOLUTION_SCHEMA,status=status,A=A,b=b,
        particular=particular,kernel=kernel,rational_particular=rational_particular,
        certificate=(S=S,U=U,V=V,transformed_rhs=c,rank=r),obstruction=obstruction,
        assignment_status=:not_requested)
    verify_integer_weight_solution(result)
    result
end

"""
    verify_integer_weight_solution(result) -> true

Recheck a stored result using exact multiplication, determinants, divisibility,
and equations. This does NOT recompute a Smith form. The unimodularity of `V`
and the exact equality `kernel == V[:,rank+1:end]` certify the whole integer
kernel lattice, beyond rational independence or `A*kernel==0` alone. Throws
`ArgumentError` for an invalid certificate.
"""
function verify_integer_weight_solution(result)
    result.schema==INTEGER_WEIGHT_LIFT_SCHEMA && return _verify_integer_weight_lift(result)
    req=_weight_certificate_require
    req(result.schema==INTEGER_WEIGHT_SOLUTION_SCHEMA,"unknown schema")
    A,b=_weight_integer_inputs(result.A,result.b)
    m,n=size(A); cert=result.certificate
    S,U,V=cert.S,cert.U,cert.V
    req(size(S)==(m,n) && size(U)==(m,m) && size(V)==(n,n),"transform dimensions")
    req(all(x -> x isa Integer,Iterators.flatten((S,U,V))),"noninteger transform")
    req(U*A*V==S,"U*A*V is not S")
    req(abs(BigInt(AbstractAlgebra.det(_weight_exact_matrix(U))))==1,"U is not unimodular")
    req(abs(BigInt(AbstractAlgebra.det(_weight_exact_matrix(V))))==1,"V is not unimodular")
    req(all(i==j || S[i,j]==0 for i in 1:m,j in 1:n),"S is not diagonal")
    r=cert.rank
    req(r isa Integer && 0<=r<=min(m,n),"rank range")
    req(all(S[i,i]>0 for i in 1:r),"nonpositive pivot")
    req(all(S[i,i]==0 for i in r+1:min(m,n)),"nonzero trailing pivot")
    req(all(mod(S[i+1,i+1],S[i,i])==0 for i in 1:r-1),"Smith divisibility chain")
    c=U*b
    req(all(x -> x isa Integer,cert.transformed_rhs) && cert.transformed_rhs==c,"transformed right-hand side")
    req(all(x -> x isa Integer,result.kernel) && size(result.kernel)==(n,n-r) &&
        result.kernel==V[:,r+1:n],"full lattice basis")
    req(A*result.kernel==zeros(BigInt,m,n-r),"kernel equation")
    rational_feasible=all(c[i]==0 for i in r+1:m)
    integer_feasible=rational_feasible && all(mod(c[i],S[i,i])==0 for i in 1:r)
    expected=integer_feasible ? :integer_solution : rational_feasible ? :rational_only : :inconsistent
    req(result.status==expected,"feasibility status")
    if rational_feasible
        q=result.rational_particular
        req(q isa AbstractVector && length(q)==n,"rational solution dimensions")
        req(all(x -> x isa Union{Integer,Rational},q),"rational solution must be exact")
        req(A*q==b,"rational solution equation")
    else
        req(result.rational_particular===nothing,"inconsistent system has a rational solution")
    end
    if integer_feasible
        x=result.particular
        req(x isa AbstractVector && length(x)==n && all(y -> y isa Integer,x),"integer solution type")
        req(A*x==b,"integer solution equation")
        req(result.obstruction===nothing,"feasible system has an obstruction")
    else
        req(result.particular===nothing,"infeasible system has an integer solution")
        o=result.obstruction
        req(o!==nothing,"missing obstruction")
        row=o.row
        req(row isa Integer && 1<=row<=m,"obstruction row")
        req(all(x -> x isa Integer,o.row_combination) && all(x -> x isa Integer,o.lhs_coefficients) &&
            all(x -> x isa Integer,(o.divisor,o.rhs_value,o.remainder)),"noninteger obstruction")
        req(o.row_combination==U[row,:],"obstruction row combination")
        lhs=vec(transpose(U[row,:])*A)
        req(o.lhs_coefficients==lhs && o.rhs_value==c[row],"obstruction equation")
        if rational_feasible
            req(o.kind==:divisibility && 1<=row<=r,"divisibility obstruction kind")
            d=S[row,row]
            req(o.divisor==d && d>0,"obstruction divisor")
            req(all(mod(x,d)==0 for x in lhs),"divisor does not divide the row coefficients")
            req(o.remainder==mod(c[row],d) && o.remainder!=0,"divisibility obstruction remainder")
        else
            req(o.kind==:rational_inconsistency && row>r,"rational obstruction kind")
            req(all(iszero,lhs) && o.divisor==0 && c[row]!=0 && o.remainder==c[row],
                "rational inconsistency witness")
        end
    end
    true
end

function _current_integer_weight_system(v::VirtualLink;boundary=nothing)
    v.weight_mode==:integer || throw(ArgumentError("integer weight solving requires integer mode; choose an explicit integer-mode copy"))
    all(x -> x in (-1,0,1),v.orientation) || throw(ArgumentError("invalid crossing orientation"))
    _validate_gauss(v.gauss,v.orientation)
    code,orientation=real_gauss_code(v)
    visits=reduce(vcat,code;init=Int[])
    all(i -> count(==(i),visits)==1 && count(==(-i),visits)==1,eachindex(orientation)) ||
        throw(ArgumentError("one visit of each sign is required at each real crossing"))
    disks=disk_table(v)
    # Convert BEFORE multiplying, so machine-integer overflow is never part of
    # the solver's construction of the actual doubled disk equation.
    targets=_boundary_targets(disks,boundary)
    2 .* BigInt.(disks.delta),targets-BigInt.(disks.cp2)
end

function _assign_integer_weight_particular!(v::VirtualLink,result;boundary=nothing)
    verify_integer_weight_solution(result)
    A,b=_current_integer_weight_system(v;boundary=boundary)
    A==result.A && b==result.b || throw(ArgumentError("integer solution belongs to another disk equation"))
    result.status==:integer_solution || return :not_integer_feasible
    any(x -> abs(x)>_INTEGER_WEIGHT_SAFE_MAX,result.particular) && return :outside_safe_integer_range
    values=Float64[]
    for exact in result.particular
        value=try
            Float64(exact)
        catch error
            error isa Union{InexactError,OverflowError} || rethrow()
            return :not_float64_representable
        end
        isfinite(value) && BigInt(value)==exact || return :not_float64_representable
        push!(values,value)
    end
    # Stage the complete assignment; preserve all raw/virtual presentation data
    # and any prior weights if validation or storage conversion fails.
    target=deepcopy(v)
    set_weight!(target,values)
    validate_weight(target;framed=true,boundary=boundary)
    v.weights=target.weights
    v.is_weighted=target.is_weighted
    :assigned
end

"""
    solve_integer_weights(v::VirtualLink; assign=false, boundary=nothing)

Solve the current presentation's exact disk equations over the integers.
Closed rows use `2*delta*w+cp2=0`; open rows use the measured
`2*delta*w-Ndots` and require `boundary`, an explicit target mapping keyed by
stable endpoint identities (or a row-aligned target vector). Edge coordinates
are the `real_edge_table` traversal order, including the real arcs that pass
through virtual crossings. A solution of the disk equation on a nonclosed
presentation is not a closed framing claim.

No mutation occurs by default. With explicit `assign=true`, the result reports
`:assigned`, `:not_integer_feasible`, `:outside_safe_integer_range`, or
`:not_float64_representable` in `assignment_status`. In the
last two cases the input remains unchanged and the exact result is still
returned. The storage failure refers to this particular solution; it does not
assert that its affine lattice contains no other representable solution.
Assignment uses the same contiguous integer range `abs(w)<=2^53-1` as moves.
Assignment is initial weight selection, never move transport. This function
does not replace or call the legacy `calculate_weight` API.
"""
function solve_integer_weights(v::VirtualLink;assign=false,parity=nothing,boundary=nothing)
    assign isa Bool || throw(ArgumentError("assign must be an explicit Boolean"))
    parity===nothing || return solve_integer_weight_lift(v,parity;assign=assign,boundary=boundary)
    A,b=_current_integer_weight_system(v;boundary=boundary)
    result=solve_integer_weight_system(A,b)
    closed,conditions=is_closed(v)
    code,orientation=real_gauss_code(v)
    presentation=(gauss=deepcopy(v.gauss),orientation=copy(v.orientation),real_only=v.real_only,
        strand_types=copy(v.strand_types),real_gauss=code,real_orientation=orientation,
        closed_normal=closed,closed_conditions=conditions)
    assignment=assign ? _assign_integer_weight_particular!(v,result;boundary=boundary) : :not_requested
    merge(result,(presentation=presentation,boundary_targets=_boundary_targets(disk_table(v),boundary),
                  assignment_status=assignment))
end

"""
    solve_integer_weight_lift(A, b, parity)
    solve_integer_weight_lift(v::VirtualLink, parity; assign=false, boundary=nothing)

Solve the additional constraint `w mod 2 == parity`, with each specified
coordinate exactly 0 or 1. The transformed integer system is
`(2A)*t=b-A*parity`, `w=parity+2t`. Its complete solver result is stored as
`t_solution`; that Smith certificate concerns `t`, not the original `w`.
The returned `kernel=2*t_solution.kernel` parametrizes all solutions with THIS
parity and is generally not the full unrestricted integer kernel of `A`.

`:rational_only` means this parity has no integral lift even though the
underlying rational equation is consistent. It does not imply that every
parity is obstructed or that any manifold-level structure is absent. This is
an initial lifting calculation, not a CP transport rule.
"""
function solve_integer_weight_lift(A::AbstractMatrix,b::AbstractVector,parity::AbstractVector)
    A,b=_weight_integer_inputs(A,b)
    Base.require_one_based_indexing(parity)
    length(parity)==size(A,2) || throw(DimensionMismatch("one parity value per edge is required"))
    all(x -> x isa Integer && x in (0,1),parity) ||
        throw(ArgumentError("parity coordinates must be integers 0 or 1"))
    z=BigInt.(parity)
    t_solution=solve_integer_weight_system(2 .* A,b-A*z)
    particular=t_solution.particular===nothing ? nothing : z+2 .* t_solution.particular
    rational=t_solution.rational_particular===nothing ? nothing : z+2 .* t_solution.rational_particular
    result=(schema=INTEGER_WEIGHT_LIFT_SCHEMA,status=t_solution.status,A=A,b=b,parity=z,
        particular=particular,kernel=2 .* t_solution.kernel,rational_particular=rational,
        t_solution=t_solution,obstruction=t_solution.obstruction,assignment_status=:not_requested)
    verify_integer_weight_solution(result)
    result
end

function _verify_integer_weight_lift(result)
    req=_weight_certificate_require
    A,b=_weight_integer_inputs(result.A,result.b)
    z=result.parity
    req(length(z)==size(A,2) && all(x -> x isa Integer && x in (0,1),z),"lift parity")
    t=result.t_solution
    req(t.schema==INTEGER_WEIGHT_SOLUTION_SCHEMA,"lift must contain an ordinary integer certificate")
    verify_integer_weight_solution(t)
    req(t.A==2 .* A && t.b==b-A*z,"lift transformed equation")
    req(result.status==t.status && result.obstruction==t.obstruction,"lift status or obstruction")
    req(all(x -> x isa Integer,result.kernel) && result.kernel==2 .* t.kernel,"full fixed-parity lattice")
    expected=t.particular===nothing ? nothing : z+2 .* t.particular
    req(result.particular==expected,"lift particular reconstruction")
    rational=t.rational_particular===nothing ? nothing : z+2 .* t.rational_particular
    req(result.rational_particular==rational,"lift rational reconstruction")
    rational===nothing || req(all(x -> x isa Union{Integer,Rational},result.rational_particular),
                              "lift rational solution must be exact")
    if expected!==nothing
        req(all(x -> x isa Integer,result.particular),"noninteger lift particular")
        req(A*expected==b && mod.(expected,2)==z,"lift equation or parity")
    end
    true
end

function solve_integer_weight_lift(v::VirtualLink,parity::AbstractVector;assign=false,boundary=nothing)
    assign isa Bool || throw(ArgumentError("assign must be an explicit Boolean"))
    A,b=_current_integer_weight_system(v;boundary=boundary)
    result=solve_integer_weight_lift(A,b,parity)
    closed,conditions=is_closed(v)
    code,orientation=real_gauss_code(v)
    presentation=(gauss=deepcopy(v.gauss),orientation=copy(v.orientation),real_only=v.real_only,
        strand_types=copy(v.strand_types),real_gauss=code,real_orientation=orientation,
        closed_normal=closed,closed_conditions=conditions)
    assignment=assign ? _assign_integer_weight_particular!(v,result;boundary=boundary) : :not_requested
    merge(result,(presentation=presentation,boundary_targets=_boundary_targets(disk_table(v),boundary),
                  assignment_status=assignment))
end
