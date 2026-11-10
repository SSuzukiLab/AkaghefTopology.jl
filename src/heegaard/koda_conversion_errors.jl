# Source: lib/algebra.jl/docs/260909_interface_design.md (値がないとき);
# docs/specs/heegaard-koda-conversion.md §2.2–2.3 (failure kinds).
"""Typed failures for the HC ↔ o-data conversion pipeline."""

abstract type KodaConversionFailure <: Exception end

struct KodaInputError <: KodaConversionFailure
    stage::Symbol
    message::String
    detail::Dict{String, Any}
end

struct KodaSpecificationPending <: KodaConversionFailure
    stage::Symbol
    topic::Symbol
    message::String
    detail::Dict{String, Any}
end

struct KodaImplementationError <: KodaConversionFailure
    stage::Symbol
    message::String
    detail::Dict{String, Any}
end

"""未計算: the definition does not exclude the input, but this encoding does not
produce a value. Not 未定義. See lib/algebra.jl/docs/260909_interface_design.md."""
struct KodaNotComputed <: KodaConversionFailure
    stage::Symbol
    reason::Symbol
    detail::Dict{String, Any}
end

"""Reserved for an exclusion written in the definition. Vertex-free singular
circles are not this case: see KodaNotComputed."""
struct KodaDomainExcluded <: KodaConversionFailure
    stage::Symbol
    condition::Symbol
    witness::Dict{String, Any}
end

function Base.showerror(io::IO, err::KodaDomainExcluded)
    print(io, "Koda domain excluded at ", err.stage, " (", err.condition, ")")
end

function Base.showerror(io::IO, err::KodaInputError)
    print(io, "Koda input error at ", err.stage, ": ", err.message)
end
function Base.showerror(io::IO, err::KodaSpecificationPending)
    print(io, "Koda specification pending at ", err.stage, " (", err.topic, "): ", err.message)
end
function Base.showerror(io::IO, err::KodaNotComputed)
    print(io, "Koda not computed at ", err.stage, " (", err.reason, ")")
end
function Base.showerror(io::IO, err::KodaImplementationError)
    print(io, "Koda implementation error at ", err.stage, ": ", err.message)
end

struct KodaConversionResult{T}
    value::T
    choices::Dict{String, Any}
    correspondence::Dict{String, Any}
    source_status::Symbol
end

function koda_not_computed(stage::Symbol, reason::Symbol, detail::Dict{String, Any}=Dict{String, Any}())
    throw(KodaNotComputed(stage, reason, detail))
end

function koda_partial(stage::Symbol, message::String; detail=Dict{String, Any}(), partial=nothing)
    throw(KodaImplementationError(stage, message, merge(Dict("partial" => partial), detail)))
end
