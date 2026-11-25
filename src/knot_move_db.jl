"""Area-2 adapter between signed Gauss diagrams and the verified move store."""

# The store is a child of the caller when a script includes it first, and a
# sibling once this file is included from KnotGauss.
if isdefined(@__MODULE__, :MoveGraphStore) &&
        getfield(@__MODULE__, :MoveGraphStore) isa Module &&
        parentmodule(getfield(@__MODULE__, :MoveGraphStore)) === @__MODULE__
    @eval using .MoveGraphStore
else
    @eval using ..MoveGraphStore
end

const KNOT_MOVE_ADAPTER_VERSION = 1

function _knot_payload(d::KnotGaussDiagram)
    Dict("gauss" => collect(d.gauss), "orientation" => collect(d.orientation))
end

"""
    knot_canonical_state(d::KnotGaussDiagram) -> CanonicalState

Convert a signed Gauss diagram to its canonical move-store state.

# Arguments
`d` satisfies K1 ([Gauss conditions](../specs/knot-gauss.md#numbered-conditions)).
The diagram constructor establishes that representation condition; no
classicality check is performed here. Store-state construction checks its
supported version and metadata under M1
([store conditions](../specs/move-graph-store.md#numbered-conditions)).

# Returns
A `CanonicalState` with canonical text, hash and versioned state ID, area,
representation, complexity and payload defined by the
[Knot adapter](../specs/move-graph-store.md#knot-adapter).

# Guarantees
1. The state describes the canonical coordinate class of `d` and its payload agrees with the stored canonical serialization. [proved/agent: by construction]
2. Cyclic root shifts and consistent crossing relabellings give the same state identity. [proved/agent: by construction]
3. The input is unchanged and no database is written. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[store identity](../specs/move-graph-store.md#textual-state-identity).

# Not guaranteed
The state ID is not a complete knot-type invariant or a collision-free
mathematical identifier. Classical realizability and move equivalence are
separate from canonical coordinate identity.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Knot adapter — canonical diagram-to-storage mapping.

# Examples
```jldoctest
julia> d = KnotGaussDiagram(Int[], Int[]); s = knot_canonical_state(d);

julia> s.serialization == Topology.KnotGauss.canonical_serialization(d)
true
```

See also [`persist_reidemeister_transition!`](@ref).

Role: convert
""" 
function knot_canonical_state(d::KnotGaussDiagram)
    c = canonicalize(d)
    CanonicalState(canonical_serialization(c); area=AREA_2, representation="KnotGauss",
                   payload=_knot_payload(c), complexity=Float64(crossing_count(c)))
end

function _candidate_metadata(before, after, candidate)
    Dict("adapter_version" => KNOT_MOVE_ADAPTER_VERSION,
         "theorem" => "Reidemeister local witness with independent source-bound verifier",
         "move_type" => candidate.type,
         "direction" => String(candidate.direction),
         "variant" => String(candidate.variant),
         "vertices" => collect(candidate.vertices),
         "edge_positions" => [collect(edge) for edge in candidate.edge_positions],
         "local_before" => collect(candidate.before),
         "local_after" => collect(candidate.after),
         "source_serialization" => canonical_serialization(before),
         "target_serialization" => canonical_serialization(after))
end

"""
    persist_reidemeister_transition!(store::MoveGraphStoreDB, before::KnotGaussDiagram, candidate::ReidemeisterCandidate) -> Int

Persist a verified signed-Gauss Reidemeister transition and its local witness and return the stored transition identifier.

# Arguments
`before` and `candidate` satisfy K1 and K3
([Gauss conditions](../specs/knot-gauss.md#numbered-conditions)); stale source
or failed independent verification raises `ArgumentError`. Local application
checks and errors are those of `apply_reidemeister`. `store` is open and
writable with the supported schema and satisfies M3
([store conditions](../specs/move-graph-store.md#numbered-conditions)).
Store validation and SQLite failures propagate. A manually malformed witness
may raise indexing errors before a validation message.

# Returns
The integer transition identifier returned by the move store. Source and
target states and evidence use the [Knot adapter](../specs/move-graph-store.md#knot-adapter).

# Guarantees
1. Before persistence, the candidate is applied to its source and the resulting transition passes `verify_reidemeister_transition`. [proved/agent: by construction]
2. Successful return refers to a stored transition between the canonical source and target, carrying the local candidate, both serializations and adapter/verifier versions. [proved/agent: by construction]
3. The store is modified; `before` and `candidate` are unchanged. [proved/agent: by construction]

# Conventions
[Knot adapter direction labels](../specs/move-graph-store.md#knot-adapter).

# Not guaranteed
No global knot-type or classicality decision is performed. The source
identification of the local move rules retains the Gauss specification's
evidence boundary. No additional all-or-nothing guarantee is made beyond the
underlying store operation, and repeated calls need not create a new state.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Knot adapter — persisted local-move evidence.

# Examples
```jldoctest
julia> store = open_store(":memory:"); d = KnotGaussDiagram(Int[], Int[]);

julia> id = try
           persist_reidemeister_transition!(store, d, first(movable(d, "R1"; param=:insert)))
       finally
           close_store!(store)
       end;

julia> id isa Integer
true
```

See also [`expand_reidemeister_neighbors!`](@ref).

Role: move
""" 
function persist_reidemeister_transition!(store::MoveGraphStore.MoveGraphStoreDB,
                                          before::KnotGaussDiagram,
                                          candidate::ReidemeisterCandidate)
    candidate.source == before || throw(ArgumentError("candidate is stale for before diagram"))
    after = apply_reidemeister(before, candidate)
    verify_reidemeister_transition(before, after, candidate) ||
        throw(ArgumentError("independent Reidemeister verification failed"))
    source = knot_canonical_state(before)
    target = knot_canonical_state(after)
    metadata = _candidate_metadata(before, after, candidate)
    verifier = (s, t, m) -> verify_reidemeister_transition(before, after, candidate)
    direction = candidate.direction == :insert ? :forward :
                candidate.direction == :delete ? :inverse : :forward
    MoveGraphStore.insert_verified_transition!(store, source, target, candidate.type;
        direction=direction,
        location=metadata["edge_positions"], parameters=metadata,
        engine_version="Topology/$(KNOT_MOVE_ADAPTER_VERSION)",
        verifier_version="Reidemeister/$(KNOT_MOVE_ADAPTER_VERSION)",
        verifier=verifier, verified=true,
        certificate_paths=[metadata])
end

"""
    expand_reidemeister_neighbors!(store::MoveGraphStoreDB, before::KnotGaussDiagram; store_cap::Integer, max_crossings::Integer=typemax(Int)) -> NamedTuple

Persist a bounded selection of verified one-step signed-Gauss neighbors in the move store and return attempted, persisted and rejected counts with the Boolean capped flag.

# Arguments
`store`, `store_cap` and `max_crossings` satisfy M6
([store conditions](../specs/move-graph-store.md#numbered-conditions)). Negative
bounds raise `ArgumentError` naming the bound. `before` satisfies K1
([Gauss conditions](../specs/knot-gauss.md#numbered-conditions)). The local
witness checks are those of `apply_reidemeister` and
`persist_reidemeister_transition!`. Application and persistence exceptions
for individual candidates are counted as rejected attempts; failures outside
those guarded operations can still propagate.

# Returns
A named tuple `(attempted, persisted, rejected, capped)` with integer counts
of candidates visited, successful persistence calls and caught failures, and
a Boolean saying the cap stopped further eligible work. `persisted` counts
calls in this invocation, not the number of newly inserted database rows.

# Guarantees
1. `persisted <= store_cap`; every successful persistence call uses a target with crossing count at most `max_crossings`, and no target canonical state ID is attempted for persistence twice in this invocation. [proved/agent: by construction]
2. Successful writes carry the adapter's source-bound local witness and verification metadata. Their mathematical scope is the finite local transition described by the adapter specification. [proved/agent: by construction]
3. The store is modified by successful writes; `before` is unchanged. Earlier writes remain if a later attempt fails. [proved/agent: by construction]

# Not guaranteed
`capped=false` does not mean that no candidate was rejected. Counts need not
sum to `attempted`, since filtered or duplicate targets are not failures.
The result does not certify all knot-equivalent states or a complete move
neighborhood after a cap or failure. Database exceptions can be summarized
as rejections rather than propagated individually.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Knot adapter — stored finite transition evidence and its source boundaries.

# Examples
```jldoctest
julia> store = open_store(":memory:");

julia> bounded = try
           expand_reidemeister_neighbors!(store, KnotGaussDiagram(Int[], Int[]); store_cap=1, max_crossings=1)
       finally
           close_store!(store)
       end;

julia> (bounded.persisted <= 1, bounded.capped || bounded.persisted < 1)
(true, true)
```

See also [`persist_reidemeister_transition!`](@ref).

Role: move
""" 
function expand_reidemeister_neighbors!(store::MoveGraphStore.MoveGraphStoreDB,
                                        before::KnotGaussDiagram; store_cap::Integer,
                                        max_crossings::Integer=typemax(Int))
    store_cap >= 0 || throw(ArgumentError("store_cap must be nonnegative"))
    max_crossings >= 0 || throw(ArgumentError("max_crossings must be nonnegative"))
    attempted = persisted = rejected = 0
    seen = Set{String}()
    for typ in ("R1", "R2", "R3")
        for candidate in movable(before, typ)
            attempted += 1
            after = try
                apply_reidemeister(before, candidate)
            catch
                rejected += 1
                continue
            end
            crossing_count(after) <= max_crossings || continue
            sid = knot_canonical_state(after).state_id
            sid in seen && continue
            push!(seen, sid)
            persisted >= store_cap && return (attempted=attempted, persisted=persisted,
                                              rejected=rejected, capped=true)
            try
                persist_reidemeister_transition!(store, before, candidate)
                persisted += 1
            catch
                rejected += 1
            end
        end
    end
    (attempted=attempted, persisted=persisted, rejected=rejected, capped=false)
end
