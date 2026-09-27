"""Return the deterministic integer-framed normal o-graph input for `L(p,q)`.

The result follows the local all-positive `lens_space_link(p,q)` orientation
convention and stores one weight per Gauss visit (the weight of the entering
edge).  Generation uses exact integer arithmetic and a fixed recursive patch;
it performs no move search, manifold recognition, or weight solve.

The general construction is currently `proved/agent`, not owner-checked.  Its
finite verification record belongs to PJ260829.
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

"""Construct a fresh weighted `VirtualLink` for the deterministic framed `L(p,q)` example."""
function framed_lens_space_link(p::Integer,q::Integer)
    input=framed_lens_space_input(p,q)
    link=VirtualLink()
    set_data!(link;rgauss=input.rgauss,orientation=input.orientation)
    set_weight!(link,input.weights)
end

"""Materialize all coprime deterministic framed lens-space examples in a finite range.

Each row is `(p, q, link)` and owns a fresh `VirtualLink`, so experiments may
mutate one row without changing another.
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
