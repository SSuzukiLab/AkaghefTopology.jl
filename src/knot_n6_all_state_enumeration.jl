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

"""Compact result returned by `stream_canonical_knot_gauss`.

No diagram or canonical collection is retained in this object.  The
`emitted_fingerprint` is SHA-256 of emitted canonical serializations in
emission order, each followed by one newline; it is a stream audit value, not
the per-state `canonical_hash` value.
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

"""Emit repository-canonical signed Gauss states without retaining a census.

Arguments and keywords:

* `n` is bounded by the repository iterator to `0:6`.
* `callback` may receive the complete one-argument record, a typed
  diagram-only argument, `(diagram, planar)`, or
  `(diagram, serialization, planar)`.  It is optional when only
  counts/fingerprints are wanted.
* `classical_canonical_set` is an externally-owned set of canonical
  serialization strings.  `planar_membership` is an alternative predicate on
  those strings; supplying both is an error.
* `max_raw` and `max_seconds` stop the stream before the next raw item.  They
  are intended for bounded N=6 benchmarks and are not needed for an exact
  run through N=5.

The repository identity convention is cyclic root shift followed by
first-occurrence relabeling; mirror and traversal reversal remain distinct.
The source iterator itself is never collected, and the function retains only
scalar counters plus the SHA-256 context.
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
stream_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) =
    stream_canonical_knot_gauss(n; callback=callback, kwargs...)
stream_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) =
    stream_canonical_knot_gauss(n; callback=callback, kwargs...)

"""Compatibility spelling for callers treating this as an enumerator."""
enumerate_canonical_knot_gauss(n::Integer; kwargs...) =
    stream_canonical_knot_gauss(n; kwargs...)
enumerate_canonical_knot_gauss(n::Integer, callback::Function; kwargs...) =
    stream_canonical_knot_gauss(n, callback; kwargs...)
enumerate_canonical_knot_gauss(callback::Function, n::Integer; kwargs...) =
    stream_canonical_knot_gauss(callback, n; kwargs...)
