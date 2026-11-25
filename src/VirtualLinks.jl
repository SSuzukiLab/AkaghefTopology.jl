"""O-graphs and links: the `VirtualLink` presentation, its moves, weights, and search database."""
module VirtualLinks

using LinearAlgebra
using JuMP
using GLPK
using Graphs
using ..MoveGraphStore
import ..Algebra: SparseFactor, FiniteHopfAlgebra, TaftKind, _contract, _factor_with_variables
import ..URDiagrams: URDiagram, URAlgebra, URMSTExpression, set_data_cvs!, trace_value
import ..MoveAPI: movable, move, move!,
    search_state, put_search_state!, load_search_state,
    begin_search_run!, finish_search_run!, expand_search_state!,
    search_neighbors, query_search_paths, verify_search_store
import AkaghefAlgebra
import AkaghefAlgebra: simplify_invariant, evaluate_invariant

include("virtual_link.jl")
include("framed_lens_space.jl")
include("circuit_b.jl")
include("integer_weight_solver.jl")
include("ograph_state.jl")
include("ograph_colours.jl")
include("ograph_colour_normal.jl")
include("ograph_face_cycles.jl")
include("weight_validation.jl")
include("spin_weights.jl")
include("layer8_planarization.jl")
include("virtual_planarization.jl")
include("orthogonal_layout.jl")
include("layout_invariants.jl")
include("helpers.jl")
include("moves.jl")
include("ograph_move_verifier.jl")
include("ograph_tetrahedra.jl")
include("weighted_ps.jl")
include("weighted_move_verifier.jl")
include("spin_cp.jl")
include("hopf_link.jl")
include("taft_invariant.jl")
include("ur_mst_link.jl")
include("ograph_gauss_diagram.jl")
include("ograph_move_db.jl")
include("ograph_search_db.jl")
include("ograph_path_constraints.jl")
include("ograph_path_lifting.jl")
include("ograph_invariants.jl")

# VirtualLink methods of the MoveAPI search names. The ograph_search_* functions stay.
search_state(link::VirtualLink) = ograph_search_state(link)
"""
    put_search_state!(store::MoveGraphStoreDB, link::VirtualLink) -> CanonicalState

Store an unweighted virtual-link search vertex through the common move API.

# Arguments
`store` and `link` meet the checked preconditions of
[`put_ograph_search_state!`](@ref). Schema, representation and persistence
errors propagate unchanged; weighted inputs are projected to unweighted data.

# Returns
The `CanonicalState` of the stored presentation, as in the o-graph-specific
entry point, including its identity and serialization.

# Guarantees
1. This method has the result and database effects of
   `put_ograph_search_state!(store, link)`; the link is unchanged.
   [proved/agent: by construction]
2. The identity is that of an [unweighted search presentation](../specs/ograph-search-db.md#contract),
   not a manifold equivalence class. [proved/agent: by construction]

# Not guaranteed
Weights are not saved, existing metadata are not replayed, and the two storage
steps are not one transaction. This method does not certify closed-normality.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — unweighted exact presentation identity.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); s = Topology.MoveAPI.put_search_state!(store, vl_example("s3_abalone"));

julia> load_ograph_search_state(store, s.state_id).is_weighted
false

julia> close_store!(store);
```

See also [`put_ograph_search_state!`](@ref), [`load_ograph_search_state`](@ref).

Role: other
"""
put_search_state!(store::MoveGraphStoreDB, link::VirtualLink) =
    put_ograph_search_state!(store, link)
"""
    load_search_state(store::MoveGraphStoreDB, state_id::AbstractString) -> Union{Nothing,VirtualLink}

Load an unweighted virtual-link search presentation through the common move API.

# Arguments
`store` must have readable o-graph search tables. An absent joined state
returns `nothing`; schema, JSON and reconstruction errors of
[`load_ograph_search_state`](@ref) propagate unchanged.

# Returns
A fresh unweighted integer-mode `VirtualLink` reconstructed from the stored
presentation, or `nothing` when absent.

# Guarantees
1. This method returns `load_ograph_search_state(store, state_id)` and leaves
   database rows unchanged. Mutating the result does not update those rows.
   [proved/agent: by construction]
2. It recovers the [search presentation](../specs/ograph-search-db.md#contract),
   including its raw diagram fields, with weights removed.
   [proved/agent: by construction]

# Not guaranteed
This two-argument method selects the o-graph storage format. It does not
replay fingerprints, validate closed-normality or recover a previous framing.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — exact unweighted presentation storage.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); s = put_ograph_search_state!(store, vl_example("s3_abalone"));

julia> Topology.MoveAPI.load_search_state(store, s.state_id).is_weighted
false

julia> close_store!(store);
```

See also [`load_ograph_search_state`](@ref), [`put_ograph_search_state!`](@ref).

Role: access
"""
load_search_state(store::MoveGraphStoreDB, state_id::AbstractString) =
    load_ograph_search_state(store, state_id)
"""
    begin_search_run!(store::MoveGraphStoreDB, script_path::AbstractString; label::AbstractString, scope=Dict{String,Any}()) -> String

Register an o-graph search run through the common move API and return its identifier.

# Arguments
Use an initialized search `store`, readable existing `script_path`, text
`label` and JSON-compatible `scope`, as required by
[`begin_ograph_search_run!`](@ref). Its checked schema, absent-file,
serialization and database errors propagate unchanged.

# Returns
A `String` identifying the new running search record.

# Guarantees
1. This method has exactly the result and writes of `begin_ograph_search_run!`
   with these arguments, including the script path and content digest snapshot.
   [proved/agent: by construction]
2. It records computational provenance for the
   [unweighted search workflow](../specs/ograph-search-db.md#workflow).
   It does not execute a search. [proved/agent: by construction]

# Not guaranteed
The scope is metadata, not an enforced bound or a mathematical certificate.
The script need not actually be running; timestamp-based identity need not
distinguish identical registrations at the same timestamp.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — search-run provenance.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl");

julia> Topology.MoveAPI.begin_search_run!(store, path; label="literal-02-smoke"); Topology.MoveAPI.begin_search_run!(store, path; label="literal-02-reuse");

julia> verify_ograph_search_store(store).runs
2

julia> close_store!(store);
```

See also [`begin_ograph_search_run!`](@ref), [`ograph_search_run`](@ref).

Role: other
"""
function begin_search_run!(store::MoveGraphStoreDB, script_path::AbstractString;
                           label::AbstractString, scope=Dict{String,Any}())
    begin_ograph_search_run!(store, script_path; label, scope)
end
"""
    finish_search_run!(store::MoveGraphStoreDB, run_id::AbstractString; status::AbstractString, counts=Dict{String,Any}()) -> Union{Nothing,NamedTuple}

Record a terminal o-graph search outcome and return its run summary, or nothing for an absent run.

# Arguments
Use the search store, terminal status and JSON-compatible counts accepted by
[`finish_ograph_search_run!`](@ref). Invalid status throws `ArgumentError`;
serialization and database errors propagate. The run may be absent or finished.

# Returns
The `ograph_search_run` summary, or `nothing` for an absent run. It contains
provenance, status, reported counts and linked-expansion totals.

# Guarantees
1. The result and database effects equal those of `finish_ograph_search_run!`
   on these arguments. Only a running named run is updated.
   [proved/agent: by construction]
2. Terminal status records a caller-reported outcome within the
   [search workflow](../specs/ograph-search-db.md#workflow), not a theorem
   about all move paths. [proved/agent: by construction]

# Not guaranteed
It neither cancels executing jobs nor validates claimed counts or completeness.
Already terminal runs retain their stored outcome.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Exact H-insertion feasibility — finite-scope terminal run statuses.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke");

julia> result = Topology.MoveAPI.finish_search_run!(store, r; status="completed-none"); result.candidates == result.successes
true

julia> close_store!(store);
```

See also [`finish_ograph_search_run!`](@ref), [`ograph_search_run`](@ref).

Role: other
"""
function finish_search_run!(store::MoveGraphStoreDB, run_id::AbstractString;
                            status::AbstractString, counts=Dict{String,Any}())
    finish_ograph_search_run!(store, run_id; status, counts)
end
"""
    expand_search_state!(store::MoveGraphStoreDB, run_id::AbstractString, source::VirtualLink, move_kind::AbstractString; selector=Dict{String,Any}()) -> Vector{NamedTuple}

Store one selected unweighted move neighborhood and return its operation outcome records through the common move API.

# Arguments
`store` must have the search schema, `run_id` must exist, and `source` must
be coherent. `move_kind` and `selector` follow [`ograph_search_candidates`](@ref).
Unsupported kinds, H, a previously failed expansion or another running owner's
unfinished expansion throw `ArgumentError`. Enumeration errors propagate after
recording failure; candidate application errors are recorded in returned outcomes.

# Returns
Ordered operation rows for the expansion, including operation and endpoint
IDs, move kind/direction, engine version, candidate JSON, outcome and error JSON.
Rejected or errored operations have no target state.

# Guarantees
1. This method has the same result and writes as `expand_ograph_search_state!`.
   Completed matching expansions are reused and linked to this run. The source
   link is unchanged. [proved/agent: by construction]
2. Operations concern the [unweighted search graph](../specs/ograph-search-db.md#contract);
   outcomes distinguish successful applications, rejected candidates and errors.
   [proved/agent: by construction]

# Not guaranteed
An expansion can contain failures; completion does not mean every candidate
succeeded. No weighted admissibility or closed-normal certificate is implied.
Search and provenance rows can remain after an interrupted or failed expansion.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract and §Workflow — reusable selected neighborhoods and recorded outcomes.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke");

julia> edges = Topology.MoveAPI.expand_search_state!(store, r, vl_example("s2xs1"), "02"; selector=Dict("direction" => "forward"));

julia> (!isempty(edges), all(String(e.outcome) == "success" for e in edges))
(true, true)

julia> close_store!(store);
```

See also [`expand_ograph_search_state!`](@ref), [`ograph_search_neighbors`](@ref).

Role: other
"""
function expand_search_state!(store::MoveGraphStoreDB, run_id::AbstractString,
                              source::VirtualLink, move_kind::AbstractString;
                              selector=Dict{String,Any}())
    expand_ograph_search_state!(store, run_id, source, move_kind; selector)
end
"""
    search_neighbors(store::MoveGraphStoreDB, state_id::AbstractString; kinds=("MP", "02")) -> Vector{NamedTuple}

Read stored outgoing o-graph moves through the common move API.

# Arguments
Use a readable search store and a state ID. `kinds` must contain only MP,
literal 02 or PS; unsupported kinds throw `ArgumentError`. Database errors
propagate. Empty kinds and absent sources return empty vectors.

# Returns
The rows of [`ograph_search_neighbors`](@ref), carrying operation ID,
endpoint IDs, move kind, direction, engine version and candidate JSON.

# Guarantees
1. The result equals `ograph_search_neighbors(store, state_id; kinds)`:
   only successful stored operations, sorted by move kind and operation ID.
   Database rows are unchanged. [proved/agent: by construction]
2. These are directed edges of the [unweighted presentation graph](../specs/ograph-search-db.md#contract),
   with no new mathematical certification. [proved/agent: by construction]

# Not guaranteed
The call does not generate neighbors, replay edges or certify manifold
equivalence. Empty output does not prove there are no possible moves.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — the persisted directed geometric graph.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> all(String(e.move_kind) == "02" for e in Topology.MoveAPI.search_neighbors(store, s.state_id; kinds=("02",)))
true

julia> close_store!(store);
```

See also [`ograph_search_neighbors`](@ref), [`query_ograph_search_paths`](@ref).

Role: access
"""
function search_neighbors(store::MoveGraphStoreDB, state_id::AbstractString; kinds=("MP", "02"))
    ograph_search_neighbors(store, state_id; kinds)
end
"""
    query_search_paths(store::MoveGraphStoreDB, source_state_id::AbstractString; kwargs...) -> Union{Nothing,Vector}

Read bounded directed o-graph search paths through the common move API.

# Arguments
`store` must have the search schema. Supply exactly one of `target_state_id`
and `target_geometry_hash`. Keyword defaults are `kinds=("MP","02")`,
`kind_pattern=nothing`, `min_depth=1`, `max_depth=12`, `limit=1000`.
Depth bounds must be ordered and nonnegative, and the limit positive;
violations throw `ArgumentError`. A supplied move-kind pattern sets both
depths to its length. Otherwise `kinds` selects MP, literal 02 or PS; unsupported
selected kinds throw `ArgumentError`. SQL and stored-data errors propagate.

# Returns
A vector of paths, each a vector of step records with IDs, move kind,
direction, engine version and parsed candidate data. No matches gives an empty
vector. Empty `kinds` without a pattern instead returns `nothing`.

# Guarantees
1. The result equals `query_ograph_search_paths(store, source_state_id; kwargs...)`:
   paths use successful stored edges and are ordered by depth then operation-ID
   sequence, up to the limit. Database rows are unchanged.
   [proved/agent: by construction]
2. The [search-path restriction](../specs/ograph-search-db.md#exact-h-insertion-feasibility)
   applies: states do not recur except that the source may be the final endpoint.
   These are unweighted presentation paths. [proved/agent: by construction]

# Not guaranteed
This does not expand the graph, certify weighted admissibility, or identify
manifolds by a geometry hash. No bounded stored path is not a nonexistence proof.
An empty pattern is not promised to form a valid SQL query in this implementation.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow and §Exact H-insertion feasibility — path targets, geometric selection and final-source loops.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); edges = expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> target = String(first(edges).target_state_id); path = query_ograph_search_path(store, s.state_id; target_state_id=target, kinds=("02",), max_depth=1);

julia> Topology.MoveAPI.query_search_paths(store, s.state_id; target_state_id=target, kinds=("02",), max_depth=1, limit=10) == [path]
true

julia> close_store!(store);
```

See also [`query_ograph_search_paths`](@ref), [`replay_weighted_ograph_path`](@ref).

Role: access
"""
function query_search_paths(store::MoveGraphStoreDB, source_state_id::AbstractString; kwargs...)
    query_ograph_search_paths(store, source_state_id; kwargs...)
end
"""
    verify_search_store(store::MoveGraphStoreDB) -> NamedTuple

Return o-graph search-schema and record counts after checking storage integrity through the common move API.

# Arguments
`store` must be open with readable search tables. As in
[`verify_ograph_search_store`](@ref), missing or unsupported schema metadata
and foreign-key failures throw `ArgumentError`; database and parsing errors
can occur first.

# Returns
The search-store summary: integer `schema_version`, `states`, `runs`,
`expansions`, `operations`, `successful_edges` and `weighted_trials` counts.

# Guarantees
1. The result equals `verify_ograph_search_store(store)` and the database is
   unchanged. Successful return confirms the current version and absence of
   reported foreign-key violations. [proved/agent: by construction]
2. Counts refer to [search records](../specs/ograph-search-db.md#workflow),
   not distinct spaces or mathematically verified equivalences.
   [proved/agent: by construction]

# Not guaranteed
For this method, verification does not replay transitions, recompute hashes,
check weighted trials or establish neighborhood completeness. The generic
API's broad description does not supply those checks here.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — record families of the o-graph search store.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl");

julia> begin_ograph_search_run!(store, path; label="literal-02-smoke"); begin_ograph_search_run!(store, path; label="literal-02-reuse");

julia> Topology.MoveAPI.verify_search_store(store).runs
2

julia> close_store!(store);
```

See also [`verify_ograph_search_store`](@ref), [`replay_weighted_ograph_path`](@ref).

Role: condition
"""
verify_search_store(store::MoveGraphStoreDB) = verify_ograph_search_store(store)

end
