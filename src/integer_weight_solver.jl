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
    solve_integer_weight_system(A::AbstractMatrix, b::AbstractVector) -> NamedTuple

Solve an integer linear system with a certificate of its complete solution lattice or infeasibility.

# Arguments
Checked: W4 of [exact integer solving](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility):
one-based indexing, row/right-hand-side alignment and integer coefficients.
Misalignment raises `DimensionMismatch`, noninteger entries `ArgumentError`;
exact-algebra and certificate-check failures propagate.
Caller: identify the supplied coordinates with the intended geometric problem.

# Returns
A versioned record with `A`, `b`, `status`, `particular`, `kernel`,
`rational_particular`, `certificate`, `obstruction` and
`assignment_status=:not_requested`, with the meanings in §12.0.

# Guarantees
1. The status distinguishes integral solvability, rational-only solvability,
   and rational inconsistency. When integral, the particular solution and
   kernel columns parameterize every integer solution uniquely.
   [proved/agent: by construction]
2. The full integer kernel is returned even when the inhomogeneous equation
   is infeasible; an obstruction witnesses the reported failure.
   [proved/agent: by construction]
3. Input arrays are not modified. [proved/agent: by construction]

# Not guaranteed
This is a statement about the supplied matrix, not every presentation of a
manifold; the returned lattice is not a quotient by geometric equivalence.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — Smith solution records and their algebraic scope.

# Examples
```jldoctest
julia> r = solve_integer_weight_system(reshape([2], 1, 1), [1]); (r.status, r.particular === nothing)
(:rational_only, true)
```

Role: construct
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
    verify_integer_weight_solution(result) -> Bool

Verify the exact certificate and solution-space claims of an integer weight result.

# Arguments
Checked: the ordinary or fixed-parity schema, exact matrix data, transforms,
lattice basis, feasibility status, solution and obstruction of
[§12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
Failed certificate assertions raise `ArgumentError`; absent fields, invalid
shapes, indexing or exact-algebra failures may propagate other exceptions.
Caller: provide a structurally readable solver record; geometric coordinate
identification remains the caller's responsibility.

# Returns
`true` if all certificate assertions pass; failure throws instead of returning false.

# Guarantees
1. Success certifies the claimed complete integer solution lattice, or
   fixed-parity affine lattice for a lift record, and the reported feasibility
   obstruction. [proved/agent: by construction]
2. The result record is not modified. [proved/agent: by construction]

# Not guaranteed
The check does not authenticate a record's origin, its optional presentation
metadata, or any manifold-level interpretation. It is not a check of
assignment history.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — exact certificates and fixed-parity records.

# Examples
```jldoctest
julia> verify_integer_weight_solution(solve_integer_weight_system(reshape([2], 1, 1), [1]))
true
```

Role: condition
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
    solve_integer_weights(v::VirtualLink; assign=false, parity=nothing, boundary=nothing) -> NamedTuple

Solve the current presentation's exact integer disk equations, optionally with prescribed parity.

# Arguments
Checked: Boolean `assign`, integer mode, valid Gauss/orientation data, one visit
of each sign per real crossing, and applicable W2 boundary targets of
[§12.0–12.1](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
Invalid values raise `ArgumentError`; target lengths and row keys may raise
`DimensionMismatch` or `BoundsError`. A nonnothing `parity` uses
`solve_integer_weight_lift` and its checks.
Caller: supply the geometric hypotheses and circuit identification for a
framing interpretation. Existing assigned weights are not required.

# Returns
An ordinary or fixed-parity solution record, with presentation metadata,
normalized `boundary_targets`, and the assignment outcome of §12.0.

# Guarantees
1. The exact solution space concerns this presentation's closed or relative
   disk equation in real-edge coordinates. [proved/agent: by construction]
2. Default operation leaves `v` unchanged. Requested assignment installs the
   validated particular solution only when safe storage succeeds; infeasible
   or unrepresentable outcomes leave the diagram unchanged.
   [proved/agent: by construction]

# Not guaranteed
A nonclosed solution is not a closed framing claim. Storage failure for the
chosen particular solution does not exclude another representable point in
the affine lattice. Initial selection is not move transport.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — integer selection, storage and geometric scope.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); r = solve_integer_weights(v; assign=true); (r.status, r.assignment_status, is_framed_weight(v))
(:integer_solution, :assigned, true)
```

Role: convert
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
    solve_integer_weight_lift(A::AbstractMatrix, b::AbstractVector, parity::AbstractVector) -> NamedTuple

Solve an integer linear system with a fixed parity in every coordinate.

# Arguments
Checked: W4 of [§12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility),
one-based parity indexing and one integer binary parity entry per column.
Length failures raise `DimensionMismatch`; invalid coefficients or parity
raise `ArgumentError`; exact solver errors propagate.
Caller: interpret these matrix coordinates and the specified parity geometrically.

# Returns
A fixed-parity lift record with `A`, `b`, `parity`, feasibility `status`,
optional particular solutions, column `kernel`, the certified `t_solution`,
an optional obstruction and `assignment_status=:not_requested`.

# Guarantees
1. An integral result parameterizes exactly all solutions with the specified
   parity; the returned columns form a basis for its translation lattice.
   [proved/agent: by construction]
2. The obstruction and status concern this fixed parity, with the
   [lift semantics](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility).
   [proved/agent: by construction]
3. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
The kernel is generally not the unrestricted integer kernel. Failure does not
rule out another parity or a geometric structure in another presentation;
the rational particular need not have a meaningful integral parity.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — fixed-parity substitution and certificate meaning.

# Examples
```jldoctest
julia> r = solve_integer_weight_lift(reshape([1], 1, 1), [1], [0]); r.status
:rational_only
```

Role: construct
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

"""
    solve_integer_weight_lift(v::VirtualLink, parity::AbstractVector; assign=false, boundary=nothing) -> NamedTuple

Solve the presentation's integer disk equation with a fixed parity on its real edges.

# Arguments
Checked: the integer-mode diagram and boundary conditions of
`solve_integer_weights`, Boolean `assign`, and the exact binary parity vector
conditions of the matrix overload. Invalid data raise `ArgumentError`, length
failures `DimensionMismatch`, and invalid row keys `BoundsError`; other errors propagate.
Caller: use [real-edge coordinates](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility)
and supply the geometric hypotheses for a framing interpretation.

# Returns
A fixed-parity lift record with presentation metadata, normalized boundary
targets, and an assignment outcome.

# Guarantees
1. The record certifies exactly the integer disk-equation solutions having
   this parity. [proved/agent: by construction]
2. Default operation does not mutate the diagram. Requested assignment
   installs the validated particular solution only after safe storage succeeds;
   infeasibility or storage failure leaves it unchanged.
   [proved/agent: by construction]

# Not guaranteed
This is initial lifting, not CP transport. Failure for one parity is not a
manifold-level obstruction; unsafe storage of one particular solution does
not exclude another representable solution.

# Sources
cf: [virtual-link moves §12.0](../specs/virtual-link-moves.md#120-initial-integer-weights-and-diagnostic-compatibility) — fixed parity and assignment semantics.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); p = Int.(mod.([e.weight for e in real_edge_table(v)], 2));

julia> r = solve_integer_weight_lift(v, p); (r.status, mod.(r.particular, 2) == p)
(:integer_solution, true)
```

Role: convert
"""
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
