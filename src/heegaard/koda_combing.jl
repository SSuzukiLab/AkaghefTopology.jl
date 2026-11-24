# Combing plus an optional second field, for b1 = 0.
# p1 and theta stay uncomputed here. Nothing is not the number zero.
# The Lescop evaluator is a consumer of this type and is not called from here.

"""A combing on a closed oriented 3-manifold with b1 = 0, plus an optional second field.

`presentation` is `:cmr3`. The combing state is one of `:combed`,
`:framing_extendable`, and `:euler_obstructed`. A mod-2 solution is
`mod2_status`, not a framing. `framing_status` is `:unjudged`, `:extendable`,
`:framed`, or `:obstructed`. `:framed` means an integer second field for
`δx = -c_P` is stored. `spin_c` is `nothing` because the Euler class does not
identify a Spin^c structure. `p1` and `theta` are `nothing` while uncomputed.
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

"""Build the b1 = 0 combing record from one CMR3 state.

Pass `b1` explicitly. This constructor rejects every other first Betti number.
`solution` is the integer vector of `δx = -c_P` when the class vanishes.
`mod2` is stored separately and never becomes the second field.
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

function combing_framing_from_cmr3(reading::CMR3Reading, c_p::AbstractVector{<:Integer};
        b1::Integer, mod2::Union{Nothing,CMR3State}=nothing)
    combing_framing_from_cmr3(reading.state; b1=b1, c_p=c_p,
        solution=reading.integer_solution, mod2=mod2)
end
