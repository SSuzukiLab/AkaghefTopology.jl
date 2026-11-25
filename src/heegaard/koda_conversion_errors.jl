# Source: lib/algebra.jl/docs/260909_interface_design.md (値がないとき);
# docs/specs/heegaard-koda-conversion.md §2.2–2.3 (failure kinds).
"""Typed failures for the HC ↔ o-data conversion pipeline."""

abstract type KodaConversionFailure <: Exception end

"""
    KodaInputError(stage, message, detail) -> KodaInputError

Store an input-failure diagnostic for a stage of the Koda conversion.

# Arguments
`stage::Symbol`, `message::String`, `detail::Dict{String,Any}`, or convertible
values. Conversion errors propagate. No stage name or diagnostic content is
validated, and constructing the exception does not throw it.

# Returns
`KodaInputError` with the stage, human-readable message and diagnostic data.
`showerror` displays the stage and message.

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed detail
   dictionary is retained and not modified. [proved/agent: by construction]
2. The record carries the input-error category of
   [the conversion result specification](../specs/heegaard-koda-conversion.md#23-返却する結果),
   rather than a converted mathematical object. [proved/agent: by construction]

# Not guaranteed
An arbitrary constructed record does not establish a violated mathematical
condition. Use the emitting operation's contract to interpret the failure.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.3 — conversion failure categories and diagnostic scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = KodaInputError(:surface_identification, "orientation required", Dict{String,Any}()); (e.stage, e.message)
(:surface_identification, "orientation required")
```

See also [`KodaSpecificationPending`](@ref), [`KodaImplementationError`](@ref).

Role: construct
"""
struct KodaInputError <: KodaConversionFailure
    stage::Symbol
    message::String
    detail::Dict{String, Any}
end

"""
    KodaSpecificationPending(stage, topic, message, detail) -> KodaSpecificationPending

Store an unresolved specification choice encountered by a Koda conversion.

# Arguments
`stage::Symbol`, `topic::Symbol`, `message::String`, and
`detail::Dict{String,Any}`, or convertible values. Conversion errors
propagate. No stage, topic or explanatory content is validated.

# Returns
`KodaSpecificationPending` with the stage, unresolved topic, message and
details. `showerror` displays the stage, topic and message. Construction
does not throw the exception.

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed detail
   dictionary is retained and not modified. [proved/agent: by construction]
2. This is the specification-pending result category of
   [the Koda conversion boundary](../specs/heegaard-koda-conversion.md#23-返却する結果).
   It supplies no converted surface or spine. [proved/agent: by construction]

# Not guaranteed
A pending choice does not imply that the mathematical construction is
undefined. Constructing this record does not check whether the named choice
is still unresolved or whether supplying it suffices to complete conversion.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.3 — unresolved specification choices as a distinct failure category.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> KodaSpecificationPending(:odata_to_hc, :surface_identification, "choose an identification", Dict{String,Any}()).topic
:surface_identification
```

See also [`odata_to_hc`](@ref), [`KodaInputError`](@ref).

Role: construct
"""
struct KodaSpecificationPending <: KodaConversionFailure
    stage::Symbol
    topic::Symbol
    message::String
    detail::Dict{String, Any}
end

"""
    KodaImplementationError(stage, message, detail) -> KodaImplementationError

Store a diagnostic for a Koda conversion stage that could not complete its construction.

# Arguments
`stage::Symbol`, `message::String`, and `detail::Dict{String,Any}`, or
convertible values. Conversion errors propagate. No diagnostic is checked,
and constructing the value does not throw it.

# Returns
`KodaImplementationError` with stage, explanatory message and detail data,
which may contain a partial result. `showerror` displays stage and message.

# Guarantees
1. Fields preserve supplied values after conversion; an already-typed detail
   dictionary is retained and not modified. [proved/agent: by construction]
2. This is a failure record in
   [the conversion result boundary](../specs/heegaard-koda-conversion.md#23-返却する結果),
   not a completed conversion or a mathematical obstruction theorem.
   [proved/agent: by construction]

# Not guaranteed
Arbitrary details, including a `partial` value, are not checked for validity
or completeness. The exception category does not prove that the input is
valid or that the desired mathematical object fails to exist.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.3 — incomplete constructions and optional partial results.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = KodaImplementationError(:hc_to_odata, "construction incomplete", Dict{String,Any}()); (e.stage, e.message)
(:hc_to_odata, "construction incomplete")
```

See also [`hc_to_odata`](@ref), [`KodaInputError`](@ref).

Role: construct
"""
struct KodaImplementationError <: KodaConversionFailure
    stage::Symbol
    message::String
    detail::Dict{String, Any}
end

"""
    KodaNotComputed(stage, reason, detail) -> KodaNotComputed

Store a Koda conversion failure for a value the current encoding has not computed.

# Arguments
`stage::Symbol`, `reason::Symbol` and `detail::Dict{String,Any}`, or
convertible values. Conversion errors propagate. Any stage/reason tags are
accepted; construction does not throw the resulting exception.

# Returns
A `KodaNotComputed` with the supplied stage, reason and diagnostic details.

# Guarantees
1. Fields preserve supplied values after conversion; a typed detail
   dictionary is retained without mutation. [proved/agent: by construction]
2. The record distinguishes unavailable computation from a definition's
   exclusion in [Koda status terminology](../specs/heegaard-koda-conversion.md#12-状態語).
   [proved/agent: by construction]

# Not guaranteed
An arbitrary reason tag does not prove that an input belongs to the
mathematical domain or that a proposed implementation would succeed.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = KodaNotComputed(:spine_to_odata, :vertex_free_odata_not_encoded, Dict{String,Any}()); (e.stage, e.reason)
(:spine_to_odata, :vertex_free_odata_not_encoded)
```

See also [`KodaDomainExcluded`](@ref).

Role: construct
"""
struct KodaNotComputed <: KodaConversionFailure
    stage::Symbol
    reason::Symbol
    detail::Dict{String, Any}
end

"""
    KodaDomainExcluded(stage, condition, witness) -> KodaDomainExcluded

Store a claimed exclusion by a mathematical condition of the Koda conversion domain.

# Arguments
`stage::Symbol`, `condition::Symbol` and `witness::Dict{String,Any}`, or
convertible values. Conversion errors propagate. The constructor does not
validate the claimed condition or throw the exception value.

# Returns
A `KodaDomainExcluded` with the supplied stage, condition and witness.

# Guarantees
1. Fields preserve supplied values after conversion; a typed witness
   dictionary is retained without mutation. [proved/agent: by construction]
2. The exception record can identify a domain exclusion under
   [Koda status terminology](../specs/heegaard-koda-conversion.md#12-状態語).
   [proved/agent: by construction]

# Not guaranteed
The record itself does not prove exclusion. Unsupported vertex-free encoding
is reported with `KodaNotComputed`, rather than establishing a domain exclusion.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = KodaDomainExcluded(:example, :condition, Dict{String,Any}()); (e.stage, e.condition)
(:example, :condition)
```

See also [`KodaNotComputed`](@ref).

Role: construct
"""
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

"""
    KodaConversionResult(value, choices::Dict{String,Any}, correspondence::Dict{String,Any}, source_status::Symbol) -> KodaConversionResult

Package a Koda conversion value with its choices, correspondence and evidence tag.

# Arguments
`value` may have any type `T`; the two dictionaries and the status symbol
must have the displayed types for this inferred constructor. A mismatched
signature raises `MethodError`. No mathematical validity or status is checked.

# Returns
`KodaConversionResult{T}` holding `value`, recorded construction `choices`,
source-to-result `correspondence`, and `source_status`. Their interpretation
is owned by [the result specification](../specs/heegaard-koda-conversion.md#23-返却する結果)
and the operation producing the value.

# Guarantees
1. All arguments are retained without copying or modification.
   [proved/agent: by construction]
2. The wrapper carries the supplied object and correspondence together; it
   does not change the mathematical object stored as `value`.
   [proved/agent: by construction]

# Not guaranteed
No correspondence, source-status claim or preservation of manifold, branching,
combing or framing is established by this constructor. A status symbol is
supplied metadata, not a validation result. Shared containers remain mutable.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§2.3, 6.4 — result packaging and the source-to-spine correspondence.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "HeegaardAbalone.jl"));

julia> e = HeegaardAbaloneExample.fixture(); r = H.hc_to_odata(e.diagram; joining_disk=e.disk); s = H.KodaConversionResult(r.value, r.choices, r.correspondence, r.source_status);

julia> (s.value isa Topology.VirtualLink, s.correspondence["source_object"], s.correspondence["image"])
(true, "HeegaardDiagram", "VirtualLink")
```

See also [`hc_to_odata`](@ref), [`heegaard_to_spine`](@ref).

Role: convert
"""
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
