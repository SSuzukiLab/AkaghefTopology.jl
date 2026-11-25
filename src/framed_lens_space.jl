"""
    framed_lens_space_input(p::Integer, q::Integer) -> NamedTuple

Return the deterministic integer-weighted normal o-graph presentation of the framed lens-space family.

# Arguments
Checked: F1 of [framed-family conditions](../specs/virtual-link-examples.md#framed-family-conditions-and-stored-data).
Caller: F2; outside that range the construction guarantees do not apply.
Boolean parameters are rejected with `ArgumentError`; conversion to `Int`
must succeed. The ordering, positivity and coprimality check raises
`ArgumentError` naming the required coprime ordered domain. Overflow is not
generally checked. Internal consistency failures raise `ErrorException`.

# Returns
A named tuple with `schema`, input `p,q`, `rgauss`, `orientation`, `weights`,
`orientation_convention`, `edge_weight_convention`, `construction`,
branch-specific `parameters`, `search_used` and `linear_solver_used`.
The word/sign/weight fields are inputs for the
[integer-framed example route](../specs/virtual-link-examples.md#one-integer-framed-lens-space).
The metadata identifies the construction; parameter fields differ by branch.

# Guarantees
1. The record gives the deterministic presentation and integer weights of the linked construction, with both search/solver flags false. Inputs are not modified and each call allocates fresh word/sign/weight containers. [proved/agent: by construction]
2. On the representable domain the construction represents the integer-framed lens-space family in the repository's naming convention. [proved/agent: threads/PJ260829/docs/notes/2026-09-27-034949-deterministic-framed-lens-intake.md, inherited general construction proof]

# Conventions
[Lens-space names and orientation](../conventions/LensSpaceNames.md) and
[real-edge weight coordinates](../specs/virtual-link-data-model.md).

# Not guaranteed
The general construction proof is inherited from the linked record, not
owner-checked. The function performs no manifold recognition or framing
verification on its output. Minimal vertex count, a canonical framing,
Hopf degree and a translation to every external lens convention are not
established. Fixed-width overflow voids the construction guarantee.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-examples.md Evidence boundary — general agent-proof status and finite framed-family evidence.
cf:threads/PJ260829/docs/notes/2026-09-27-034949-deterministic-framed-lens-intake.md — deterministic family, proof provenance and finite checks.

# Examples
```jldoctest
julia> x = framed_lens_space_input(5,2);

julia> (x.parameters.r, x.parameters.A, length(x.orientation))
(4, 4, 20)
```

Role: construct
"""
function framed_lens_space_input(p::Integer,q::Integer)
    p isa Bool && throw(ArgumentError("p must be an integer, not Bool"))
    q isa Bool && throw(ArgumentError("q must be an integer, not Bool"))
    p=Int(p); q=Int(q)
    p>q>0 && gcd(p,q)==1 ||
        throw(ArgumentError("require coprime integers p>q>0"))

    if q==p-1
        word=vcat([-p],collect(-1:-1:-(p-1)),collect(1:p))
        weights=zeros(Int,2p)
        weights[1]=1
        weights[p+1]=-1
        weights[end]=-1
        return (
            schema="framed-lens-o-graph-v2",
            p=p,
            q=q,
            orientation_convention="local_all_positive_lens_space_link",
            rgauss=[word],
            orientation=ones(Int,p),
            weights=[weights],
            edge_weight_convention="edge entering visit; word read cyclically",
            construction="previous_q_minus_one_base",
            parameters=(vertices=p,),
            search_used=false,
            linear_solver_used=false,
        )
    end

    m=p+q
    d=gcd(2,p)
    modulus=p÷d
    r=mod(((q-1)÷d)*invmod((2q)÷d,modulus),modulus)
    numerator=2r*(p-q)-(p-q+1)
    mod(numerator,p)==0 || error("internal modular-division inconsistency")
    A=numerator÷p

    word,signs,augmented_weights=_framed_lens_patch(r,A,0)
    patch_vertices=length(signs)
    residues=[_framed_lens_residue(k*p,m) for k in m:-1:1]
    exterior=vcat(collect(p:m),-residues,collect(1:p-1))
    exterior_word=[sign(v)*(abs(v)+patch_vertices) for v in exterior]
    exterior_weights=[_framed_lens_exterior_weight(v,p,q,m,r,A) for v in exterior]
    augmented_weights[end]==0 && exterior_weights[1]==0 ||
        error("patch/cap edge convention mismatch")

    word,signs,weights=_framed_lens_mirror(
        vcat(word,exterior_word),
        vcat(signs,ones(Int,m)),
        vcat(augmented_weights[begin:end-1],exterior_weights),
    )
    (
        schema="framed-lens-o-graph-v2",
        p=p,
        q=q,
        orientation_convention="local_all_positive_lens_space_link",
        rgauss=[word],
        orientation=signs,
        weights=[weights],
        edge_weight_convention="edge entering visit; word read cyclically",
        construction="fixed_insertion_in_L_p_plus_q_p_then_P_r_and_mirror",
        parameters=(
            m=m,
            d=d,
            r_modulus=modulus,
            r=r,
            A=A,
            B=0,
            patch_vertices=patch_vertices,
            vertices=length(signs),
            congruence="2*q*r == q-1 (mod p)",
            raw_seed_orientation="opposite; mirror already applied",
        ),
        search_used=false,
        linear_solver_used=false,
    )
end

"""
    framed_lens_space_link(p::Integer, q::Integer) -> VirtualLink

Construct the weighted VirtualLink of the deterministic framed lens-space family.

# Arguments
Checked: F1 of [framed-family conditions](../specs/virtual-link-examples.md#framed-family-conditions-and-stored-data).
Caller: F2–F3; construction and exact stored-weight guarantees require these conditions.
The scalar input and internal consistency errors are those of `framed_lens_space_input(p,q)`.
VirtualLink data/weight assignment errors also propagate; construction does
not run an independent closed-normal or framing check.

# Returns
A fresh `VirtualLink` with the generated real-Gauss word, crossing signs and
assigned integer-mode weights, in the linked family's coordinates.

# Guarantees
1. Real-Gauss data and assigned weights equal the fields returned by `framed_lens_space_input(p,q)`, subject to the storage conversions. [proved/agent: by construction]
2. The represented framed lens-space interpretation has that construction's inherited general proof. [proved/agent: threads/PJ260829/docs/notes/2026-09-27-034949-deterministic-framed-lens-intake.md]
3. Inputs are unchanged; each call returns independent mutable link containers. [proved/agent: by construction]

# Conventions
[Lens-space names and orientation](../conventions/LensSpaceNames.md).

# Not guaranteed
No canonical framing or minimal presentation is supplied. Integer generation and floating weight storage must both be
representable; out-of-range data have no exact-weight guarantee.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-examples.md One integer-framed lens space and Evidence boundary.

# Examples
```jldoctest
julia> is_closed(framed_lens_space_link(5,2))
(true, (true, true, true))
```

Role: construct
"""
function framed_lens_space_link(p::Integer,q::Integer)
    input=framed_lens_space_input(p,q)
    link=VirtualLink()
    set_data!(link;rgauss=input.rgauss,orientation=input.orientation)
    set_weight!(link,input.weights)
end

"""
    framed_lens_space_series(maximum_p::Integer; minimum_p::Integer=2) -> Vector

Construct a finite parameter sweep of the deterministic framed lens-space family.

# Arguments
Checked: F4 of [framed-family conditions](../specs/virtual-link-examples.md#framed-family-conditions-and-stored-data).
Every generated pair inherits F1; the caller ensures F2–F3 for every row.
Boolean bounds raise `ArgumentError`; bounds must convert to `Int`.
The F4 interval check raises `ArgumentError` naming the required bounds
if violated. Each generated pair inherits `framed_lens_space_link`'s
arithmetic, storage and allocation obligations and errors.

# Returns
A vector of `(p,q,link)` records for coprime pairs in the requested interval,
ordered by increasing `p`, then `q`. Each link is a fresh framed-family value.

# Guarantees
1. Every admissible pair in the interval occurs once, with its `framed_lens_space_link(p,q)` presentation. [proved/agent: by construction]
2. Mutable link containers are independent between rows; inputs are unchanged. [proved/agent: by construction]

# Not guaranteed
The family is not deduplicated by lens-space homeomorphism or framing class.
It retains the constructor's inherited proof and representability limits;
the whole vector is materialized, without a time or memory budget.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-examples.md Finite parameter sweep — framed-family experiment inputs.

# Examples
```jldoctest
julia> Tuple((r.p,r.q) for r in framed_lens_space_series(3))
((2, 1), (3, 1), (3, 2))
```

Role: construct
"""
function framed_lens_space_series(maximum_p::Integer;minimum_p::Integer=2)
    maximum_p isa Bool && throw(ArgumentError("maximum_p must be an integer, not Bool"))
    minimum_p isa Bool && throw(ArgumentError("minimum_p must be an integer, not Bool"))
    maximum_p=Int(maximum_p); minimum_p=Int(minimum_p)
    2<=minimum_p<=maximum_p ||
        throw(ArgumentError("require 2 <= minimum_p <= maximum_p"))
    [(p=p,q=q,link=framed_lens_space_link(p,q))
     for p in minimum_p:maximum_p for q in 1:p-1 if gcd(p,q)==1]
end

_framed_lens_residue(a::Integer,modulus::Integer)=mod(a-1,modulus)+1

function _framed_lens_exterior_weight(v::Integer,p::Integer,q::Integer,
                                      m::Integer,r::Integer,A::Integer)
    v>0 && return v==1 ? -1 : 0
    t=-v
    t<=q && return -1
    t<p && return A+1-2r
    t<m && return A+1
    1
end

"""Return `(word, signs, augmented_weights)` for the fixed `P_r(A,B)` patch."""
function _framed_lens_patch(r::Integer,A::Integer,B::Integer)
    if r==0
        return [-1,2,-2,1],[-1,1],[B,0,1,A-B,0]
    elseif r==1
        return (
            [-1,-2,3,-4,2,4,5,-5,1,-3],
            [-1,-1,1,-1,1],
            [B,-1,A-B+1,-1,-1,0,0,1,0,0,0],
        )
    elseif iseven(r)
        h=r÷2
        word,signs,weights=_framed_lens_patch(h,2A-2B-4h+1,0)
        N=length(word); n=length(signs)
        a,b,c,d=n+1:n+4
        output_weights=zeros(Int,N+9)
        output_weights[4:N+4]=weights
        output_weights[1]=A
        output_weights[2]=1-2h
        output_weights[3]=-A+B+2h-1
        output_weights[N+5]=-A+B+2h-1
        return (
            vcat([a,b,-c],word,[d,-b,c,-d,-a]),
            vcat(signs,[1,-1,-1,-1]),
            output_weights,
        )
    end

    h=1-r
    word,signs,weights=_framed_lens_patch(h,-1,A-B+4h-3)
    N=length(word); n=length(signs)
    a,b,c,d=n+1:n+4
    output_weights=zeros(Int,N+9)
    output_weights[5:N+5]=weights
    output_weights[1:4]=[B-2h+2,2h-2,A-B+4h-2,1]
    (
        vcat([a,-b,-c,d],word,[-d,c,b,-a]),
        vcat(signs,fill(-1,4)),
        output_weights,
    )
end

function _framed_lens_mirror(word::Vector{Int},signs::Vector{Int},weights::Vector{Int})
    reverse(word),-signs,vcat(weights[1],reverse(weights[2:end]))
end
