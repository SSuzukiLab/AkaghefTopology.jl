# Combing plus an optional second field, for b1 = 0.
# p1 and theta stay uncomputed here. Nothing is not the number zero.
# The Lescop evaluator is a consumer of this type and is not called from here.

"""
    CombingFraming(b1, presentation, combing, euler_representative, framing_status, second_field, mod2_status, spin_c, p1, theta, certificate) -> CombingFraming

Store a combing classification with optional framing and invariant data.

# Arguments
The fields accept `Int`, `Symbol`, `CMR3State`, optional integer Euler vector,
`Symbol`, optional `Vector{BigInt}` second field, `Symbol`, `nothing`,
optional `Rational{BigInt}` p1, optional four-tuple of `Rational{BigInt}` theta,
and `Dict{String,Any}` certificate, respectively. Convertible values are
accepted and conversion errors propagate. The raw constructor does not
check Betti number, state consistency or any equation.

# Returns
A `CombingFraming` retaining the supplied classification, optional fields
and certificate annotations.

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers and
   the state are retained without modification. [proved/agent: by construction]
2. The record can carry the combing and additional-structure data associated
   with a closed oriented 3-manifold; it does not construct vector fields.
   [proved/agent: by construction]

# Conventions
CMR3 framing records: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Arbitrary fields do not certify a framing, Spin-c structure, p1 or theta.
The [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the status and second-field convention. Use `combing_framing_from_cmr3` for its checked domain.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §10.3 — distinction between manifold and additional-structure claims.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> c = CombingFraming(0, :cmr3, CMR3State(:combed,nothing,nothing), nothing, :unjudged, nothing, :unjudged, nothing, nothing, nothing, Dict{String,Any}()); (c.b1, c.second_field, c.p1)
(0, nothing, nothing)
```

See also [`combing_framing_from_cmr3`](@ref).

Role: convert
"""
struct CombingFraming
    b1::Int
    presentation::Symbol
    combing::CMR3State
    euler_representative::Union{Nothing,Vector{Int}}
    framing_status::Symbol
    second_field::Union{Nothing,Vector{BigInt}}
    mod2_status::Symbol
    spin_c::Nothing
    p1::Union{Nothing,Rational{BigInt}}
    theta::Union{Nothing,NTuple{4,Rational{BigInt}}}
    certificate::Dict{String,Any}
end

function _combing_mod2_status(mod2::Union{Nothing,CMR3State})
    mod2 === nothing && return :unjudged
    mod2.kind == :mod2_solution && return :mod2_solution
    mod2.kind == :euler_obstructed && mod2.diagnostic === :mod2_unsolvable && return :mod2_unsolvable
    throw(ArgumentError("mod2 must be a mod-2 classification"))
end

function _combing_framing_fields(state::CMR3State, solution::Union{Nothing,AbstractVector{<:Integer}})
    state.kind in (:combed, :framing_extendable, :euler_obstructed) ||
        throw(ArgumentError("combing state must be combed, framing_extendable, or euler_obstructed"))
    if state.kind == :combed
        return :unjudged, nothing
    elseif state.kind == :euler_obstructed
        return :obstructed, nothing
    elseif solution === nothing
        return :extendable, nothing
    else
        return :framed, BigInt.(solution)
    end
end

"""
    combing_framing_from_cmr3(state::CMR3State; b1, c_p=nothing, solution=nothing, mod2=nothing) -> CombingFraming

Build a zero-Betti combing record from a classification state and optional second-field data.

# Arguments
`b1` is an integer required to equal zero. `state` must be a combing or integer
obstruction classification; `mod2`, when supplied, must be a modular
classification. Checked failures raise `ArgumentError` naming the unsupported
Betti number or state. A framed result requires a supplied Euler vector and
nonempty solution; failures raise `ArgumentError`. Vector conversions can fail.
The caller guarantees that a supplied solution solves the intended equation;
this operation receives no matrix and cannot check it.

# Returns
A `CombingFraming` with the state, optional copied Euler/second-field vectors,
separate mod-two status and certificate annotations. Spin-c, p1 and theta
are absent.

# Guarantees
1. The record reflects the supplied integer classification and keeps
   mod-two evidence separate from the second field.
   [proved/agent: by construction]
2. Inputs are not modified; the optional vectors are converted into new
   arrays. [proved/agent: by construction]

# Conventions
CMR3 framing records: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
This does not independently certify a vector field, solve the integer
equation, identify a Spin-c structure or evaluate p1/theta. The [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the geometric compatibility obligations.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §10.3 — additional structures require compatible geometric correspondence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> c = combing_framing_from_cmr3(cmr3_combed(); b1=0); (c.framing_status, c.second_field, c.spin_c)
(:unjudged, nothing, nothing)
```

Role: convert
"""
function combing_framing_from_cmr3(state::CMR3State;
        b1::Integer,
        c_p::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        solution::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        mod2::Union{Nothing,CMR3State}=nothing)
    b1 == 0 || throw(ArgumentError("CombingFraming is implemented for b1 = 0 only"))
    status, field = _combing_framing_fields(state, solution)
    representative = c_p === nothing ? nothing : Int.(c_p)
    if status == :framed
        representative === nothing && throw(ArgumentError("a framed combing needs the Euler representative"))
        isempty(field) && throw(ArgumentError("a framed combing needs the second field"))
    end
    CombingFraming(
        0, :cmr3, state, representative, status, field,
        _combing_mod2_status(mod2),
        nothing, nothing, nothing,
        Dict{String,Any}(
            "scope" => "b1=0",
            "presentation" => "cmr3",
            "spin_c" => "not_identified",
            "p1" => "uncomputed",
            "theta" => "uncomputed",
            "theta_order" => "(ell_2, s_ell, e, Theta)",
            "second_field_equation" => "delta * x = -c_P",
            "mod2_is_not_framing" => true,
        ),
    )
end

"""
    combing_framing_from_cmr3(reading::CMR3Reading, c_p::AbstractVector{<:Integer}; b1, mod2=nothing) -> CombingFraming

Build a zero-Betti combing record from an Euler-class reading and representative.

# Arguments
`reading` is a compatible result of `classify_euler_class`, and `c_p` is the
same Euler representative. Compatibility is the caller's responsibility.
`b1` and optional `mod2` follow the state overload; its `ArgumentError` and
conversion failures propagate.

# Returns
The `CombingFraming` obtained from the reading's state and integer solution,
with the supplied Euler vector and separate modular status.

# Guarantees
1. The supplied reading's state and solution determine the resulting record
   under the state overload's checks. [proved/agent: by construction]
2. The reading and representative are not modified.
   [proved/agent: by construction]

# Conventions
CMR3 framing records: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
The reading's modular-witness residue fields are not retained as a geometric
certificate. This wrapper does not recheck the integer equation or certify
the manifold/field correspondence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.5, 10.3 — cellular coordinates and additional-structure scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> r = classify_euler_class(reshape([1],1,1), [0]); c = combing_framing_from_cmr3(r, [0]; b1=0); (c.framing_status, Tuple(c.second_field))
(:framed, (0,))
```

Role: convert
"""
function combing_framing_from_cmr3(reading::CMR3Reading, c_p::AbstractVector{<:Integer};
        b1::Integer, mod2::Union{Nothing,CMR3State}=nothing)
    combing_framing_from_cmr3(reading.state; b1=b1, c_p=c_p,
        solution=reading.integer_solution, mod2=mod2)
end
