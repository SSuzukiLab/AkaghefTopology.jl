# cf: .research-scratch/PJ260922/astra-lanes/theory/REPORT-2.md §数学仕様と出典 — BJP Theorem 5.5 and Deloup (16); primary Work registration pending, unverified
import AbstractAlgebra
import AkaghefAlgebra

const _HCLinkRat = Rational{BigInt}
const _HCLinkTE = AkaghefAlgebra.TaftExact

"""
    HCLinkingInputError(reason, detail) -> HCLinkingInputError

Store the reason and explanation for a rejected linking-form computation.

# Arguments
`reason::Symbol` and `detail::String`, or values convertible to these types.
Conversion errors propagate. Construction accepts any reason tag and does not
throw the resulting exception.

# Returns
`HCLinkingInputError` with the supplied reason tag and explanatory text.
`showerror` displays both fields.

# Guarantees
1. Fields equal the supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. When thrown by a linking-form operation, the value reports a failed input
   or consistency check; it is diagnostic data, not a linking pairing on a
   3-manifold. [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
The constructor does not check that the stated failure occurred. An error
record supplies no certificate of a manifold's torsion linking form.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = HCLinkingInputError(:isotropy, "isotropy required"); (e.reason, e.detail)
(:isotropy, "isotropy required")
```

See also [`hc_linking_form`](@ref).

Role: construct
"""
struct HCLinkingInputError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::HCLinkingInputError) = print(io, "HCLinkingInputError(", e.reason, "): ", e.detail)
_hcl_require(ok, reason, detail) = ok || throw(HCLinkingInputError(reason, detail))
_hcl_mod1(q) = mod(q, one(q))
_hcl_bigmatrix(M) = [BigInt(M[i,j]) for i in axes(M,1), j in axes(M,2)]

function _hcl_smith(M::Matrix{BigInt})
    m,n = size(M)
    if m == 0 || n == 0
        return copy(M), Matrix{BigInt}(I,m,m), Matrix{BigInt}(I,n,n)
    end
    s,u,v = AbstractAlgebra.snf_with_transform(AbstractAlgebra.matrix(AbstractAlgebra.ZZ,M))
    S,U,V = _hcl_bigmatrix(s),_hcl_bigmatrix(u),_hcl_bigmatrix(v)
    for i in 1:min(m,n)
        if S[i,i] < 0
            S[i,:] .*= -1; U[i,:] .*= -1
        end
    end
    U*M*V == S || error("invalid Smith transform certificate")
    S,U,V
end

function _hcl_inv_integer(M::Matrix{BigInt})
    isempty(M) && return copy(M)
    R = inv(_HCLinkRat.(M))
    _hcl_require(all(x -> denominator(x)==1,R), :nonunimodular,
        "the reference pairing or Smith transform is not unimodular")
    numerator.(R)
end

"""
    HCLinkingForm(free_rank, torsion, matrix, generators, certificate) -> HCLinkingForm

Store an exact torsion-linking pairing record and its homology-coordinate data.

# Arguments
`free_rank::Int`, `torsion::Vector{BigInt}`,
`matrix::Matrix{Rational{BigInt}}`, `generators::Matrix{BigInt}` and
`certificate::NamedTuple`, or convertible values. Conversion errors propagate.
The constructor does not check shapes, symmetry, nondegeneracy or compatibility
with a surface presentation; these are the caller's responsibility.

# Returns
An `HCLinkingForm` retaining the supplied free rank, torsion orders, exact
pairing matrix, generator lifts and certificate fields.

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers are
   retained without mutation. [proved/agent: by construction]
2. The record can carry the pairing computed by `hc_linking_form`; constructing
   the record itself performs no linking-form calculation.
   [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
Arbitrary matrices need not represent the linking form of a closed oriented
3-manifold. The [linking-coordinate specification](../specs/heegaard-code.md#linking-generator-conventions)
states the Smith-generator and lift convention. Presentation-dependent arrays alone do
not certify isometry of linked groups.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> f = HCLinkingForm(0, BigInt[2], reshape(Rational{BigInt}[1//2],1,1), reshape(BigInt[1,0],2,1), NamedTuple()); (f.free_rank, Tuple(f.torsion), f.matrix[1,1])
(0, (2,), 1//2)
```

See also [`hc_linking_form`](@ref), [`linking_pairing`](@ref).

Role: construct
"""
struct HCLinkingForm
    free_rank::Int
    torsion::Vector{BigInt}
    matrix::Matrix{_HCLinkRat}
    generators::Matrix{BigInt}
    certificate::NamedTuple
end

function _hcl_lattice_form(A::Matrix{BigInt}, B::Matrix{BigInt}, omega::Matrix{BigInt};
                           boundary_side::Symbol=:alpha, reference_names=String[])
    boundary_side in (:alpha,:beta) || throw(ArgumentError("boundary_side must be :alpha or :beta"))
    n,g = size(A)
    size(B)==(n,g) && n==2g && size(omega)==(n,n) || throw(DimensionMismatch("two genus-g meridian systems required"))
    _hcl_require(omega == -transpose(omega), :reference_pairing, "surface form must be alternating")
    _hcl_inv_integer(omega)
    for (family,F) in ((:alpha,A),(:beta,B))
        S,_,_ = _hcl_smith(F)
        _hcl_require(all(i -> S[i,i]==1,1:g), :nonprimitive_meridians,
            "$family must span a primitive rank-g lattice")
        _hcl_require(iszero(transpose(F)*omega*F), :nonisotropic_meridians,
            "$family must be isotropic")
    end
    C = hcat(A,B)
    S,U,V = _hcl_smith(C)
    r = count(i -> S[i,i]!=0,1:n)
    indices = [i for i in 1:r if S[i,i]>1]
    d = BigInt[S[i,i] for i in indices]
    Uinv = _hcl_inv_integer(U)
    X = Uinv[:,indices]
    P,R = V[1:g,indices], V[g+1:n,indices]
    E,F = A*P, B*R
    E+F == X*Diagonal(d) || error("invalid meridian decomposition certificate")
    numer = transpose(E)*omega*X
    Q = _HCLinkRat[numer[i,j]//d[i] for i in eachindex(d), j in eachindex(d)]
    Q-transpose(Q) == transpose(X)*omega*X || error("invalid linking symmetry certificate")
    sign = boundary_side == :alpha ? 1 : -1
    L = _hcl_mod1.(sign .* Q)
    L == transpose(L) || error("linking form is not symmetric modulo integers")
    all(i -> all(j -> denominator(d[i]*L[i,j])==1 && denominator(d[j]*L[i,j])==1,
        eachindex(d)),eachindex(d)) || error("pairing does not descend to the torsion group")
    if !isempty(d)
        K = numerator.(Diagonal(d)*transpose(L))
        SK,_,_ = _hcl_smith(hcat(Matrix(Diagonal(d)),K))
        all(i -> SK[i,i]==1,eachindex(d)) || error("linking form is singular")
    end
    cert = (presentation=C,smith=S,left=U,right=V,alpha=A,beta=B,
        alpha_coefficients=P,beta_coefficients=R,alpha_parts=E,beta_parts=F,
        free_lifts=Uinv[:,r+1:n],surface_form=omega,reference_names=copy(reference_names),
        boundary_side=boundary_side)
    HCLinkingForm(n-r,d,L,X,cert)
end

"""
    hc_linking_form(hc::HeegaardCode; boundary_side::Symbol=:alpha) -> HCLinkingForm

Return the torsion linking-form record of the HC's meridian lattices.

# Arguments
`hc` satisfies the [Heegaard domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
The nonempty input is validated; meridian counts, self/family intersections,
connected active ribbon, reference-coordinate consistency and primitive
isotropic meridian ranks are checked with `HCLinkingInputError` naming the
condition. `boundary_side` is `:alpha` or `:beta`, otherwise `ArgumentError`.
Validation and exact-linear-algebra errors propagate. The caller supplies
the geometric handlebody interpretation of the meridian systems.

# Returns
`HCLinkingForm` containing free rank, torsion orders, a rational linking
matrix, integral generator lifts and a certificate of presentation, Smith,
meridian and reference-surface coordinates.

# Guarantees
1. The returned torsion pairing passes symmetry, denominator and
   nonsingularity checks on the reported finite group.
   [proved/agent: by construction]
2. Changing `boundary_side` negates the reported form modulo integers.
   [proved/agent: by construction]
3. The HC is not modified. [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
The [linking specification](../specs/heegaard-code.md#linking-generator-conventions)
states the coordinates and normalization. This record alone is not an independent proof identifying
the algebraic pairing with the linking form of an arbitrary claimed manifold.
Joining disks, combings and spin structures are not classified here.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#linking-generator-conventions — linking lifts, boundary orientation and identification limits.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> f = H.hc_linking_form(hc); (f.free_rank, Tuple(f.torsion), size(f.matrix))
(0, (), (0, 0))
```

Role: invariant
"""
function hc_linking_form(hc::HeegaardCode; boundary_side::Symbol=:alpha)
    g = hc.lc.genus
    if g==0 && isempty(hc.lc.words) && isempty(hc.roles)
        return _hcl_lattice_form(zeros(BigInt,0,0),zeros(BigInt,0,0),zeros(BigInt,0,0);
            boundary_side=boundary_side)
    end
    validate_diagram(hc)
    alphas,betas = hc_curve_names(hc,:alpha),hc_curve_names(hc,:beta)
    _hcl_require(length(alphas)==g && length(betas)==g, :meridian_count,
        "both families must have genus-many curves")
    owners = _hc_owners(hc.lc.words)
    adjacency = Dict(c=>Set{String}() for c in keys(hc.lc.words))
    for (plus,minus) in values(owners)
        _hcl_require(plus!=minus,:self_intersection,"curve $plus crosses itself")
        rp,rm = get(hc.roles,plus,:aux),get(hc.roles,minus,:aux)
        _hcl_require(!(rp==rm && rp in (:alpha,:beta)),:family_intersection,
            "$plus and $minus intersect within one meridian family")
        push!(adjacency[plus],minus); push!(adjacency[minus],plus)
    end
    active = Set(c for (c,w) in hc.lc.words if !isempty(w))
    reached = Set{String}(); pending = isempty(active) ? String[] : [first(active)]
    while !isempty(pending)
        c = pop!(pending)
        c in reached && continue
        push!(reached,c); append!(pending,setdiff(adjacency[c],reached))
    end
    _hcl_require(reached==active,:disconnected_ribbon,
        "a connected nonempty reference ribbon is required")
    names = vcat(["l$i" for i in 1:g],["m$i" for i in 1:g])
    function intersection(c,d)
        c==d && return big(0)
        sum((p==c && m==d ? big(-1) : m==c && p==d ? big(1) : big(0)
            for (p,m) in values(owners));init=big(0))
    end
    omega = BigInt[intersection(c,d) for c in names,d in names]
    invomega = _hcl_inv_integer(omega)
    coord(c) = invomega*BigInt[intersection(k,c) for k in names]
    A = hcat((coord(c) for c in alphas)...)
    B = hcat((coord(c) for c in betas)...)
    allnames = sort!(collect(keys(hc.lc.words)))
    coords = Dict(c=>coord(c) for c in allnames)
    for c in allnames,d in allnames
        _hcl_require(intersection(c,d)==dot(coords[c],omega*coords[d]),
            :coordinate_inconsistency,"reference coordinates disagree for $c and $d")
    end
    _hcl_lattice_form(A,B,omega;boundary_side=boundary_side,reference_names=names)
end

"""
    linking_pairing(form::HCLinkingForm, x::AbstractVector{<:Integer}, y::AbstractVector{<:Integer}) -> Rational{BigInt}

Evaluate the recorded torsion pairing on two integer coordinate vectors.

# Arguments
`x` and `y` each have one entry per element of `form.torsion`; otherwise `DimensionMismatch` is thrown. The caller supplies a consistent `HCLinkingForm`, normally from `hc_linking_form`. Matrix-size and conversion errors propagate; this accessor does not validate a hand-built form.

# Returns
The exact rational representative in `[0,1)` of the matrix pairing modulo the integers.

# Guarantees
1. The result evaluates the bilinear pairing recorded in `form.matrix` on the supplied generator coordinates. Its geometric meaning, when established, is the torsion linking pairing of the oriented manifold represented by the producer's HC. [proved/agent: matrix evaluation; geometric identification conditional on the producer]
2. Inputs are not modified. [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
This does not independently certify symmetry, nonsingularity, the generator orders or a manifold interpretation. It does not translate coordinates between different `HCLinkingForm` records.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md — reference-marked HC data and invariant interpretation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> f = H.HCLinkingForm(0, BigInt[3], reshape([big(2)//3],1,1), reshape(BigInt[1],1,1), (;)); H.linking_pairing(f, [3], [1])
0//1
```

Role: invariant
"""
function linking_pairing(form::HCLinkingForm,x::AbstractVector{<:Integer},y::AbstractVector{<:Integer})
    length(x)==length(y)==length(form.torsion) || throw(DimensionMismatch("one coordinate per torsion generator required"))
    _hcl_mod1(dot(BigInt.(x),form.matrix*BigInt.(y)))
end

function _hcl_square_parts(n::BigInt; max_squarefree::BigInt)
    n>0 || throw(ArgumentError("positive radicand required"))
    a,d = big(1),big(1)
    p = big(2)
    while p*p<=n
        s = isqrt(n)
        s*s==n && return a*s,d
        # If all primes <= the allowed squarefree part have been removed,
        # a nonsquare remainder necessarily exceeds the conductor budget.
        p<=max_squarefree || throw(ArgumentError("output radical exceeds max_conductor"))
        e = 0
        while mod(n,p)==0
            n=div(n,p); e+=1
        end
        a *= p^div(e,2)
        isodd(e) && (d*=p)
        d<=max_squarefree || throw(ArgumentError("output radical exceeds max_conductor"))
        p = p==2 ? big(3) : p+2
    end
    d*n<=max_squarefree || throw(ArgumentError("output radical exceeds max_conductor"))
    a,d*n
end

function _hcl_positive_sqrt(d::BigInt,root,c::Int)
    d==1 && return one(root)
    period = Int(4d)
    gauss4 = sum((root^Int(mod(BigInt(k)^2*(c÷period),c))
        for k in 0:period-1);init=zero(root))
    radical = gauss4*(1-root^(c÷4))*(big(1)//4)
    radical*radical==d || error("invalid positive radical certificate")
    radical
end

function _hcl_embed_cyclotomic(x::_HCLinkTE.Cyclotomic{M},root::_HCLinkTE.Cyclotomic{C}) where {M,C}
    C%M==0 || throw(ArgumentError("incompatible cyclotomic embedding"))
    C==M && return x
    sum((a*root^((j-1)*(C÷M)) for (j,a) in enumerate(_HCLinkTE.coefficients(x)));
        init=zero(root))
end

# Deloup Theorems 1–2: nonzero values are sqrt(|H¹(M;Z/N)|) times an eighth
# root of unity. This is exact field conversion, preserving the computed sum.
function _hcl_canonical_moo(value::_HCLinkTE.Cyclotomic{C},form::HCLinkingForm,N::BigInt;
                             max_conductor::Integer) where {C}
    coefficients = _HCLinkTE.coefficients(value)
    if all(iszero,coefficients[2:end])
        return AkaghefAlgebra.Cyclotomic(AkaghefAlgebra.RootSpec(2,1),coefficients[1])
    end
    kernel_order = prod((gcd(N,d) for d in form.torsion);init=big(1))
    radicand = (isodd(form.free_rank) ? N : big(1))*kernel_order
    square,d = _hcl_square_parts(radicand;
        max_squarefree=max(big(1),div(BigInt(max_conductor),4)))
    square *= N^div(form.free_rank,2)
    target = lcm(big(8),4d)
    common = lcm(BigInt(C),target)
    common<=max_conductor && common<=typemax(Int) ||
        throw(ArgumentError("canonical comparison field exceeds max_conductor"))
    t,c = Int(target),Int(common)
    root = _HCLinkTE.zeta(t); commonroot = _HCLinkTE.zeta(c)
    lifted = _hcl_embed_cyclotomic(value,commonroot)
    candidate = square*_hcl_positive_sqrt(d,root,t)
    phase = root^(t÷8)
    for _ in 0:7
        if _hcl_embed_cyclotomic(candidate,commonroot)==lifted
            return AkaghefAlgebra.Cyclotomic(candidate,AkaghefAlgebra.RootSpec(t,1))
        end
        candidate *= phase
    end
    error("MOO value has no exact eighth-root phase at its predicted magnitude")
end

"""
    Z_MOO(form::HCLinkingForm, N::Integer=2; max_terms=1_000_000, max_conductor=100_000) -> AkaghefAlgebra.Cyclotomic

Return the exact even-level MOO value computed from a torsion-linking-form record.

# Arguments
`form` must be a valid linking-form result with consistent torsion coordinates;
this method does not revalidate that geometric input. `N` must be a positive
even integer other than `Bool`. `max_terms > 0` and `max_conductor >= 2` are
required. These conditions and resource bounds are checked with `ArgumentError`
identifying the level, enumeration size or field bound. Pairing and exact
arithmetic errors propagate.

# Returns
An exact `AkaghefAlgebra.Cyclotomic` value. Its root metadata describes the
output field, rather than the requested MOO level. Zero is a valid result.

# Guarantees
1. The returned value uses exact arithmetic for the supplied form and level.
   [proved/agent: by construction]
2. Its intended interpretation is the MOO invariant of the oriented closed
   3-manifold whose linking form and first Betti number were supplied.
   [unverified: general source alignment under docs/specs/heegaard-code.md#moo-tv-normalization]
3. The input form is not modified. [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
The method does not certify the input as a manifold's linking form, choose a
spin refinement or support odd levels. The general identification in guarantee
2 remains limited by the unresolved source registration at file scope.
Resource limits can reject a value that would simplify into a smaller field.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> f = HCLinkingForm(0, BigInt[2], reshape(Rational{BigInt}[1//2],1,1), reshape(BigInt[1,0],2,1), NamedTuple()); Z_MOO(f, 2) == 0
true
```

See also [`hc_linking_form`](@ref).

Role: invariant
"""
function Z_MOO(form::HCLinkingForm,N::Integer=2;max_terms::Integer=1_000_000,
               max_conductor::Integer=100_000)
    !(N isa Bool) && N>0 && iseven(N) || throw(ArgumentError("N must be a positive even integer"))
    max_terms>0 && max_conductor>=2 || throw(ArgumentError("positive resource limits required"))
    order = prod(form.torsion;init=big(1))
    order<=max_terms || throw(ArgumentError("torsion enumeration exceeds max_terms"))
    Nbig = BigInt(N)
    # Remove the known square N^(2*floor(b1/2)) before trial division.
    radicand = (isodd(form.free_rank) ? Nbig : big(1))*order
    square,d = _hcl_square_parts(radicand;
        max_squarefree=max(big(1),div(BigInt(max_conductor),4)))
    square *= Nbig^div(form.free_rank,2)
    exponent = foldl(lcm,form.torsion;init=big(1))
    conductor = lcm(exponent,d==1 ? big(2) : 4d)
    conductor<=max_conductor && conductor<=typemax(Int) || throw(ArgumentError("output field exceeds max_conductor"))
    c = Int(conductor)
    root = _HCLinkTE.zeta(c)
    ranges = (0:Int(di)-1 for di in form.torsion)
    counts = Dict{Int,BigInt}()
    for x in Iterators.product(ranges...)
        phase = _hcl_mod1(div(Nbig,2)*linking_pairing(form,BigInt[x...],BigInt[x...]))
        power = conductor*phase
        denominator(power)==1 || error("field does not contain linking phase")
        k = Int(numerator(power))
        counts[k] = get(counts,k,big(0))+1
    end
    gauss = sum((root^k*v for (k,v) in counts);init=zero(root))
    radical = _hcl_positive_sqrt(d,root,c)
    value = gauss*radical*(square//order)
    _hcl_canonical_moo(value,form,Nbig;max_conductor=max_conductor)
end

"""
    Z_MOO(hc::HeegaardCode, N::Integer=2; boundary_side=:alpha, kwargs...) -> AkaghefAlgebra.Cyclotomic

Return the even-level MOO value of the linking form computed from an HC.

# Arguments
`hc` and `boundary_side` must satisfy `hc_linking_form`; its input and
geometry errors propagate. `N` and resource keywords follow the linking-form
overload of `Z_MOO`, including rejection of nonpositive or odd levels with
`ArgumentError`. Unsupported keywords raise `MethodError`.

# Returns
The exact cyclotomic value returned by `Z_MOO` for the computed linking form.

# Guarantees
1. The result equals evaluating the form overload on
   `hc_linking_form(hc; boundary_side)` at the same level and limits.
   [proved/agent: by construction]
2. The intended invariant belongs to the oriented closed 3-manifold of the
   valid HC and chosen boundary-side convention.
   [unverified: general linking-form identification and MOO normalization source alignment]
3. The input HC is not modified. [proved/agent: by construction]

# Conventions
Linking coordinates and state-sum normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#linking-generator-conventions).

# Not guaranteed
This wrapper does not supply an independent manifold or normalization
certificate; the linking-form overload's source and coordinate limits apply.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc, Dict("m1"=>:alpha, "l1"=>:beta), Dict(), Dict(), Dict(), nothing); Z_MOO(h, 2) == 1
true
```

See also [`hc_linking_form`](@ref).

Role: invariant
"""
Z_MOO(hc::HeegaardCode,N::Integer=2;boundary_side::Symbol=:alpha,kwargs...) =
    Z_MOO(hc_linking_form(hc;boundary_side=boundary_side),N;kwargs...)
