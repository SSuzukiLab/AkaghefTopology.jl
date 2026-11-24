# cf: .research-scratch/PJ260922/astra-lanes/theory/REPORT-2.md §数学仕様と出典 — BJP Theorem 5.5 and Deloup (16); primary Work registration pending, unverified
import AbstractAlgebra
import AkaghefAlgebra

const _HCLinkRat = Rational{BigInt}
const _HCLinkTE = AkaghefAlgebra.TaftExact

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
An exact torsion linking form in Smith generators. `torsion[i]` is the order
of generator i; `matrix[i,j]` is its pairing in Q/Z, represented in [0,1).
The columns of `generators` are lifts in the ordered l₁,…,l_g,m₁,…,m_g basis.
`free_rank` is b₁. `certificate` retains the integral presentation, both Smith
transforms, the meridian decompositions, and the actual surface form.
Generators and their matrix depend on the presentation; the linked group is
the invariant. Stored arrays follow the subsystem's immutable-by-convention API.
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
    hc_linking_form(hc::HeegaardCode; boundary_side=:alpha)

Compute Tor H₁ and its symmetric Q/Z linking form using exact integer Smith
transforms. The LC passage rule is used literally: + on c and − on d gives
c⋅d = −1. Coordinates come from intersections with the fixed l/m curves.

By default the encoded oriented surface is ∂Hα (outward normal first), so
`m*u=a+b` gives λ(u,v)=a⋅v/m. `boundary_side=:beta` selects the opposite
orientation. Reversing every passage sign also reverses the form. Individual
curve orientations and Smith generators do not affect the isometry class.

Requires genus-many embedded disjoint meridians in each family, with primitive
rank-g homology lattices. Checks passages, the connected capped ribbon surface,
reference basis, geometric disjointness, and the complete cut-system ranks.
Joining disks, weights, combing, framing, and spin structures are unused.
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

"""Evaluate a linking form on integer vectors in its reported torsion generators."""
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
    Z_MOO(hc_or_form, N::Integer=2; max_terms=1_000_000, max_conductor=100_000)

The MOO invariant at the positive primitive 2N-th root exp(πi/N), for positive
even N. Returns the existing `AkaghefAlgebra.Cyclotomic` exact value. Its root
metadata describes the *output field*, not the MOO level; there is no floating
approximation. For T=Tor H₁, the formula is
`N^(b₁/2)/sqrt(|T|) * sum(exp(+πi*N*λ(x,x)) for x in T)`.

Rational values use Q(ζ₂). Other values use Q(ζ_lcm(8,4d)), where d is the
squarefree part of the squared absolute value. Equal invariant values therefore
have equal field parameters and compare with ordinary `==`, even for different
torsion groups. Exact comparison with the original sum certifies the conversion.

The positive square root is embedded by G(4d)=2(1+i)√d after removing square
factors. Cost includes |T| summands, radical Gauss sums, and exact embeddings.
`max_conductor` bounds both construction and canonical comparison fields;
it may reject an answer that would simplify to a smaller field. Limits may be
raised explicitly. The zero Gauss sum is a valid result. Odd levels
require a separate definition and are rejected. No spin refinement is chosen.
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

Z_MOO(hc::HeegaardCode,N::Integer=2;boundary_side::Symbol=:alpha,kwargs...) =
    Z_MOO(hc_linking_form(hc;boundary_side=boundary_side),N;kwargs...)
