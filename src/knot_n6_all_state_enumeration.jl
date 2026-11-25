"""Streaming enumeration of canonical signed one-component Gauss states.

This file deliberately does not import a database or a classical-realizability
oracle.  `enumerate_knot_gauss` is the repository's authoritative rooted
iterator; this lane visits it in its native deterministic order and emits a
state exactly when `canonicalize(diagram) == diagram`.  Consequently each
cyclic-root orbit has its repository representative emitted once, including
orbits with nontrivial rotational symmetry, without retaining either the raw
inputs or an accumulated canonical `Set`.

The callback receives a named tuple
`(diagram=..., serialization=..., planar=...)`.  `planar` is `missing` unless
the caller supplies `classical_canonical_set` (a set of repository canonical
serialization strings) or `planar_membership` (a predicate on such a string).
The supplied set/predicate is only queried; this file never writes a DB or
copies it.  Thus a direct N=6 classical census can be injected later without
changing the all-state enumeration route.

Completeness is a finite representation statement: every rooted signed
diagram is visited by `enumerate_knot_gauss`; every repository canonical state
has at least one rooted representative, and the minimum representative is
exactly the one satisfying the identity test.  The identity test is also an
exact periodic-symmetry deduplication rule, so no hash/set of all states is
needed.  This says nothing about classical realizability; `planar=true` is
only the result of the supplied external membership relation.
"""

using SHA

const KNOT_N6_ALL_STATE_ENUMERATION_VERSION = 1
const KNOT_GAUSS_EXPECTED_CANONICAL_COUNTS = (1, 2, 14, 164, 3388, 96776)
const _KNOT_N6_ALL_STATE_NEWLINE = UInt8[0x0a]

"""    KnotGaussCanonicalStreamStats(n, raw_total, raw_seen, canonical_emitted, planar_emitted, nonplanar_emitted, planar_unknown, elapsed_seconds, stopped, stop_reason, emitted_serialization_bytes, max_emitted_serialization_bytes, emitted_fingerprint) -> KnotGaussCanonicalStreamStats

Store the counts and audit data of a canonical signed-Gauss stream.

# Arguments
The fields and their meaning are defined in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).
K7 ([conditions](../specs/knot-gauss.md#numbered-conditions)) is the caller's
responsibility; the field constructor does not check counts, stopping state
or fingerprint. Julia conversions to the declared field types may fail.

# Returns
A `KnotGaussCanonicalStreamStats` with the crossing count `n`, rooted counts
`raw_total`, `raw_seen`, emitted counts `canonical_emitted`, `planar_emitted`,
`nonplanar_emitted`, `planar_unknown`, timing `elapsed_seconds`, stopping data
`stopped`, `stop_reason`, and byte/fingerprint fields
`emitted_serialization_bytes`, `max_emitted_serialization_bytes`,
`emitted_fingerprint`, as defined in the linked entry.

# Guarantees
1. Fields equal the supplied values after conversion; for K7 they summarize a finite stream of the repository's canonical representatives. [proved/agent: by construction]
2. Construction does not change any argument. [proved/agent: by construction]

# Conventions
[Stream fingerprint and membership](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No diagrams or population certificate are retained. Membership counts do not
establish classical realizability independently of the supplied membership
relation. Construction does not verify a digest or run an enumeration.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration records — canonical stream and audit definitions.

# Examples
```jldoctest
julia> s = stream_canonical_knot_gauss(0);

julia> (s.canonical_emitted, s.planar_unknown, s.stopped, s.stop_reason)
(1, 1, false, :complete)
```

See also [`stream_canonical_knot_gauss`](@ref).

Role: construct
"""
struct KnotGaussCanonicalStreamStats
    n::Int
    raw_total::Int
    raw_seen::Int
    canonical_emitted::Int
    planar_emitted::Int
    nonplanar_emitted::Int
    planar_unknown::Int
    elapsed_seconds::Float64
    stopped::Bool
    stop_reason::Symbol
    emitted_serialization_bytes::Int
    max_emitted_serialization_bytes::Int
    emitted_fingerprint::String
end

const KnotGaussCanonicalEnumerationStats = KnotGaussCanonicalStreamStats

function _knot_n6_stream_planar_lookup(serialization::String,
                                       diagram::KnotGaussDiagram,
                                       lookup)
    lookup === nothing && return missing
    lookup isa AbstractSet && return serialization in lookup
    lookup isa Function && return Bool(lookup(serialization))
    throw(ArgumentError(
        "planar_membership must be an AbstractSet of canonical strings or a predicate"))
end

function _knot_n6_stream_callback(callback, record)
    callback === nothing && return nothing
    # One-argument callbacks get the auditable record.  A callback explicitly
    # typed for `KnotGaussDiagram` is also accepted for diagram-only consumers;
    # the two-argument form is provided for `(diagram, planar)`, and the
    # three-argument form exposes all three record fields positionally.
    if applicable(callback, record)
        callback(record)
    elseif applicable(callback, record.diagram)
        callback(record.diagram)
    elseif applicable(callback, record.diagram, record.planar)
        callback(record.diagram, record.planar)
    elseif applicable(callback, record.diagram, record.serialization, record.planar)
        callback(record.diagram, record.serialization, record.planar)
    else
        throw(ArgumentError(
            "callback must accept record, diagram, (diagram, planar), or " *
            "(diagram, serialization, planar)"))
    end
    nothing
end

"""Serialize a diagram already known to be repository-canonical.

This is byte-for-byte the existing `canonical_serialization` format, but it
avoids canonicalizing the same emitted diagram a second time.  The caller is
the stream's identity check, so using it on a noncanonical diagram is not
part of the API contract.
"""
function _knot_n6_stream_serialization(d::KnotGaussDiagram)
    "KnotGaussDiagram;canon=v$(KNOT_GAUSS_CANON_VERSION);n=$(crossing_count(d));" *
        "gauss=$(join(d.gauss, ','));orientation=$(join(d.orientation, ','))"
end

"""
    stream_canonical_knot_gauss(n::Integer; callback=nothing, classical_canonical_set=nothing, planar_membership=nothing, max_raw=nothing, max_seconds=nothing) -> KnotGaussCanonicalStreamStats

Stream canonical signed Gauss representatives and return their counts and fingerprint.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
Range, callback type, membership exclusivity/type and positive limits are
checked with `ArgumentError` naming the failure; numeric conversion errors
propagate. Callback and membership exceptions propagate. The callback can
also be a second positional argument or first argument in do-block form.
The one-record callback receives `(diagram, serialization, planar)` as named
fields; `planar` is `missing` without a membership source and a Boolean
otherwise. Other callback shapes are specified by K10.

# Returns
A `KnotGaussCanonicalStreamStats` with all fields defined in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).
No diagram collection is retained in the result.

# Guarantees
1. Without early termination, exactly one representative of every canonical coordinate class at `n` crossings is emitted; each emitted diagram equals its canonicalization. [proved/agent: by construction]
2. The returned counts and fingerprint describe the emitted prefix. Membership counts equal the supplied lookup's answers, or unknown counts when none is supplied. [proved/agent: by construction]
3. A limit stops before another raw input is processed and records its stopping reason; completing the population records `:complete`. [proved/agent: by construction]
4. Input sets are only queried. The routine does not mutate an input container; callback/predicate side effects remain the caller's responsibility. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream order fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No classicality or knot-type decision is performed internally. The two
membership keyword names do not certify their supplied relation. A time limit
is not an interrupt deadline inside a raw item or callback. Early termination
does not give a complete census; callback results are ignored and cannot
request termination by their return value.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite canonical population and stream evidence.

# Examples
```jldoctest
julia> s = stream_canonical_knot_gauss(0);

julia> (s.canonical_emitted, s.planar_unknown, s.stopped, s.stop_reason)
(1, 1, false, :complete)
```

Role: construct
"""
function stream_canonical_knot_gauss(n::Integer;
                                    callback=nothing,
                                    classical_canonical_set=nothing,
                                    planar_membership=nothing,
                                    max_raw=nothing,
                                    max_seconds=nothing)
    n = Int(n)
    0 <= n <= 6 || throw(ArgumentError("n must lie in 0:6"))
    callback === nothing || callback isa Function ||
        throw(ArgumentError("callback must be a function or nothing"))
    classical_canonical_set === nothing || planar_membership === nothing ||
        throw(ArgumentError(
            "supply at most one of classical_canonical_set and planar_membership"))
    lookup = classical_canonical_set === nothing ? planar_membership :
        classical_canonical_set
    lookup === nothing || lookup isa AbstractSet || lookup isa Function ||
        throw(ArgumentError(
            "planar_membership must be an AbstractSet of canonical strings or a predicate"))

    raw_limit = if max_raw === nothing
        nothing
    else
        value = Int(max_raw)
        value > 0 || throw(ArgumentError("max_raw must be positive"))
        value
    end
    time_limit = if max_seconds === nothing
        nothing
    else
        value = Float64(max_seconds)
        isfinite(value) && value > 0 ||
            throw(ArgumentError("max_seconds must be finite and positive"))
        value
    end

    iterator = enumerate_knot_gauss(n; orientations=true)
    raw_total = length(iterator)
    raw_seen = 0
    canonical_emitted = 0
    planar_emitted = 0
    nonplanar_emitted = 0
    planar_unknown = 0
    emitted_bytes = 0
    max_emitted_bytes = 0
    started_ns = time_ns()
    stopped = false
    stop_reason = :complete
    fingerprint_context = SHA.SHA256_CTX()

    for diagram in iterator
        if raw_limit !== nothing && raw_seen >= raw_limit
            stopped = true
            stop_reason = :max_raw
            break
        end
        if time_limit !== nothing &&
           (time_ns() - started_ns) / 1e9 >= time_limit
            stopped = true
            stop_reason = :max_seconds
            break
        end

        raw_seen += 1
        canonicalize(diagram) == diagram || continue

        serialization = _knot_n6_stream_serialization(diagram)
        planar = _knot_n6_stream_planar_lookup(serialization, diagram, lookup)
        canonical_emitted += 1
        if ismissing(planar)
            planar_unknown += 1
        elseif planar
            planar_emitted += 1
        else
            nonplanar_emitted += 1
        end
        bytes = ncodeunits(serialization)
        emitted_bytes += bytes
        max_emitted_bytes = max(max_emitted_bytes, bytes)
        SHA.update!(fingerprint_context, codeunits(serialization))
        SHA.update!(fingerprint_context, _KNOT_N6_ALL_STATE_NEWLINE)
        _knot_n6_stream_callback(callback,
                                  (diagram=diagram,
                                   serialization=serialization,
                                   planar=planar))
    end

    elapsed_seconds = (time_ns() - started_ns) / 1e9
    fingerprint = bytes2hex(SHA.digest!(fingerprint_context))
    KnotGaussCanonicalStreamStats(n, raw_total, raw_seen, canonical_emitted,
        planar_emitted, nonplanar_emitted, planar_unknown, elapsed_seconds,
        stopped, stop_reason, emitted_bytes, max_emitted_bytes, fingerprint)
end

# Both forms support normal calls and Julia's `f(n) do record ... end` syntax.
"""
    stream_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) -> KnotGaussCanonicalStreamStats

Stream canonical signed Gauss representatives and return their counts and fingerprint.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
Range, callback type, membership exclusivity/type and positive limits are
checked with `ArgumentError` naming the failure; numeric conversion errors
propagate. Callback and membership exceptions propagate. The callback can
also be a second positional argument or first argument in do-block form.
The one-record callback receives `(diagram, serialization, planar)` as named
fields; `planar` is `missing` without a membership source and a Boolean
otherwise. Other callback shapes are specified by K10.

# Returns
A `KnotGaussCanonicalStreamStats` with all fields defined in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).
No diagram collection is retained in the result.

# Guarantees
1. Without early termination, exactly one representative of every canonical coordinate class at `n` crossings is emitted; each emitted diagram equals its canonicalization. [proved/agent: by construction]
2. The returned counts and fingerprint describe the emitted prefix. Membership counts equal the supplied lookup's answers, or unknown counts when none is supplied. [proved/agent: by construction]
3. A limit stops before another raw input is processed and records its stopping reason; completing the population records `:complete`. [proved/agent: by construction]
4. Input sets are only queried. The routine does not mutate an input container; callback/predicate side effects remain the caller's responsibility. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream order fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No classicality or knot-type decision is performed internally. The two
membership keyword names do not certify their supplied relation. A time limit
is not an interrupt deadline inside a raw item or callback. Early termination
does not give a complete census; callback results are ignored and cannot
request termination by their return value.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite canonical population and stream evidence.

# Examples
```jldoctest
julia> s = stream_canonical_knot_gauss(0, x -> nothing);

julia> (s.canonical_emitted, s.planar_unknown, s.stopped, s.stop_reason)
(1, 1, false, :complete)
```

Role: construct
"""
stream_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) =
    stream_canonical_knot_gauss(n; callback=callback, kwargs...)
"""
    stream_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) -> KnotGaussCanonicalStreamStats

Stream canonical signed Gauss representatives and return their counts and fingerprint.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
Range, callback type, membership exclusivity/type and positive limits are
checked with `ArgumentError` naming the failure; numeric conversion errors
propagate. Callback and membership exceptions propagate. The callback can
also be a second positional argument or first argument in do-block form.
The one-record callback receives `(diagram, serialization, planar)` as named
fields; `planar` is `missing` without a membership source and a Boolean
otherwise. Other callback shapes are specified by K10.

# Returns
A `KnotGaussCanonicalStreamStats` with all fields defined in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).
No diagram collection is retained in the result.

# Guarantees
1. Without early termination, exactly one representative of every canonical coordinate class at `n` crossings is emitted; each emitted diagram equals its canonicalization. [proved/agent: by construction]
2. The returned counts and fingerprint describe the emitted prefix. Membership counts equal the supplied lookup's answers, or unknown counts when none is supplied. [proved/agent: by construction]
3. A limit stops before another raw input is processed and records its stopping reason; completing the population records `:complete`. [proved/agent: by construction]
4. Input sets are only queried. The routine does not mutate an input container; callback/predicate side effects remain the caller's responsibility. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream order fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No classicality or knot-type decision is performed internally. The two
membership keyword names do not certify their supplied relation. A time limit
is not an interrupt deadline inside a raw item or callback. Early termination
does not give a complete census; callback results are ignored and cannot
request termination by their return value.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite canonical population and stream evidence.

# Examples
```jldoctest
julia> s = stream_canonical_knot_gauss(x -> nothing, 0);

julia> (s.canonical_emitted, s.planar_unknown, s.stopped, s.stop_reason)
(1, 1, false, :complete)
```

Role: construct
"""
stream_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) =
    stream_canonical_knot_gauss(n; callback=callback, kwargs...)

"""    enumerate_canonical_knot_gauss(n::Integer; callback=nothing, classical_canonical_set=nothing, planar_membership=nothing, max_raw=nothing, max_seconds=nothing) -> KnotGaussCanonicalStreamStats

Enumerate canonical signed Gauss states through the canonical-stream interface.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
The same arguments and checks as `stream_canonical_knot_gauss` are used.
The callback may also be supplied as a second positional argument or first
argument in Julia's do-block form. Invalid range, callback type, membership
combination or limits raise `ArgumentError` describing the failure; conversion,
lookup and callback exceptions propagate.

# Returns
`KnotGaussCanonicalStreamStats` with the field meanings in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).

# Guarantees
1. This entry point has the same emitted records, stopping behavior and result semantics as `stream_canonical_knot_gauss` for the same arguments. [proved/agent: by construction]
2. Without a limit or exception, every canonical signed state for `n` is emitted exactly once. With a limit, counts describe the completed prefix. [proved/agent: by construction]
3. Input sets are only queried; user callbacks and predicates may perform their own side effects. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
The result is not a materialized collection of diagrams. Membership labels
are supplied externally and do not prove classicality. A time limit is not
an interruption deadline inside one item or callback. No knot-type quotient
or census is implied.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite representation population and stream semantics.

# Examples
```jldoctest
julia> s = enumerate_canonical_knot_gauss(0);

julia> (s.canonical_emitted, s.planar_unknown, s.stop_reason)
(1, 1, :complete)
```

See also [`stream_canonical_knot_gauss`](@ref).

Role: construct
"""
enumerate_canonical_knot_gauss(n::Integer; kwargs...) =
    stream_canonical_knot_gauss(n; kwargs...)
"""
    enumerate_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) -> KnotGaussCanonicalStreamStats

Enumerate canonical signed Gauss states through the canonical-stream interface.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
The same arguments and checks as `stream_canonical_knot_gauss` are used.
The callback may also be supplied as a second positional argument or first
argument in Julia's do-block form. Invalid range, callback type, membership
combination or limits raise `ArgumentError` describing the failure; conversion,
lookup and callback exceptions propagate.

# Returns
`KnotGaussCanonicalStreamStats` with the field meanings in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).

# Guarantees
1. This entry point has the same emitted records, stopping behavior and result semantics as `stream_canonical_knot_gauss` for the same arguments. [proved/agent: by construction]
2. Without a limit or exception, every canonical signed state for `n` is emitted exactly once. With a limit, counts describe the completed prefix. [proved/agent: by construction]
3. Input sets are only queried; user callbacks and predicates may perform their own side effects. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
The result is not a materialized collection of diagrams. Membership labels
are supplied externally and do not prove classicality. A time limit is not
an interruption deadline inside one item or callback. No knot-type quotient
or census is implied.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite representation population and stream semantics.

# Examples
```jldoctest
julia> s = enumerate_canonical_knot_gauss(0, x -> nothing);

julia> (s.canonical_emitted, s.planar_unknown, s.stop_reason)
(1, 1, :complete)
```

See also [`stream_canonical_knot_gauss`](@ref).

Role: construct
"""
enumerate_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) =
    stream_canonical_knot_gauss(n, callback; kwargs...)
"""
    enumerate_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) -> KnotGaussCanonicalStreamStats

Enumerate canonical signed Gauss states through the canonical-stream interface.

# Arguments
K4 and K10–K12 ([conditions](../specs/knot-gauss.md#numbered-conditions)) apply.
The same arguments and checks as `stream_canonical_knot_gauss` are used.
The callback may also be supplied as a second positional argument or first
argument in Julia's do-block form. Invalid range, callback type, membership
combination or limits raise `ArgumentError` describing the failure; conversion,
lookup and callback exceptions propagate.

# Returns
`KnotGaussCanonicalStreamStats` with the field meanings in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).

# Guarantees
1. This entry point has the same emitted records, stopping behavior and result semantics as `stream_canonical_knot_gauss` for the same arguments. [proved/agent: by construction]
2. Without a limit or exception, every canonical signed state for `n` is emitted exactly once. With a limit, counts describe the completed prefix. [proved/agent: by construction]
3. Input sets are only queried; user callbacks and predicates may perform their own side effects. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[stream fingerprint](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
The result is not a materialized collection of diagrams. Membership labels
are supplied externally and do not prove classicality. A time limit is not
an interruption deadline inside one item or callback. No knot-type quotient
or census is implied.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Canonical coordinates and Enumeration records — finite representation population and stream semantics.

# Examples
```jldoctest
julia> s = enumerate_canonical_knot_gauss(x -> nothing, 0);

julia> (s.canonical_emitted, s.planar_unknown, s.stop_reason)
(1, 1, :complete)
```

See also [`stream_canonical_knot_gauss`](@ref).

Role: construct
"""
enumerate_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) =
    stream_canonical_knot_gauss(callback, n; kwargs...)
