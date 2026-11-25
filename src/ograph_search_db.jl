using Dates
using JSON
using SHA

const OGRAPH_SEARCH_SCHEMA_VERSION = 1
const OGRAPH_SEARCH_STATE_VERSION = "ograph-search-state/1"
const OGRAPH_SEARCH_ENUMERATOR_VERSION = "ograph-search-enumerator/1"
const OGRAPH_SEARCH_ENGINE_VERSION = "ograph-search-engine/1"

function _ograph_search_schema_path()
    normpath(joinpath(@__DIR__, "..", "db", "ograph_search_schema.sql"))
end

"""
    migrate_ograph_search_store!(store::MoveGraphStoreDB) -> MoveGraphStoreDB

Install the current unweighted o-graph search schema in an open move-graph store.

# Arguments
`store` must be an open writable store. SQL, schema-file and parsing errors
propagate. A stored search-schema version other than the supported version
throws `ArgumentError` identifying the unsupported schema.

# Returns
The same `MoveGraphStoreDB` object, ready for the
[search workflow](../specs/ograph-search-db.md#workflow).

# Guarantees
1. On successful return the generic move-graph migration has completed and
   the search-schema version is the supported current version. [proved/agent: by construction]
2. Repeating this call on the current schema retains existing search records.
   [proved/agent: by construction]
3. Migration creates storage for unweighted presentation searches; it
   performs no move search or mathematical verification. [proved/agent: by construction]

# Not guaranteed
This is not an automatic conversion of arbitrary future or incompatible
schemas. An error after an earlier migration step need not roll back every
preceding step of the call. Existing records are not replayed or repaired.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — search-store setup and its presentation-level scope.

# Examples
```jldoctest
julia> store = open_store();

julia> migrate_ograph_search_store!(store) === migrate_ograph_search_store!(store)
true

julia> close_store!(store);
```

See also [`begin_ograph_search_run!`](@ref), [`verify_ograph_search_store`](@ref).

Role: other
"""
function migrate_ograph_search_store!(store::MoveGraphStoreDB)
    MoveGraphStore.migrate!(store)
    MoveGraphStore._transaction(store.db) do
        MoveGraphStore._execute_schema!(store.db, read(_ograph_search_schema_path(), String))
    end
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph search schema version"))
    store
end

function _require_ograph_search_store(store::MoveGraphStoreDB)
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("o-graph search schema is absent; call migrate_ograph_search_store!"))
    store
end

function _ograph_search_unweighted(link::VirtualLink)
    projected = VirtualLink(weight_mode=:integer)
    kwargs = (orientation=copy(link.orientation), strand_types=copy(link.strand_types))
    if link.real_only
        set_data!(projected; rgauss=deepcopy(link.gauss), kwargs...)
    else
        set_data!(projected; gauss=deepcopy(link.gauss), kwargs...)
    end
    projected
end

function _ograph_search_state_data(link::VirtualLink)
    projected = _ograph_search_unweighted(link)
    real_gauss, real_orientation = real_gauss_code(projected)
    raw = Dict(
        "schema" => OGRAPH_SEARCH_STATE_VERSION,
        "gauss" => projected.gauss,
        "orientation" => projected.orientation,
        "real_only" => projected.real_only,
        "strand_types" => string.(projected.strand_types),
    )
    geometry = Dict(
        "schema" => "ograph-search-geometry/1",
        "real_gauss" => real_gauss,
        "real_orientation" => real_orientation,
    )
    serialization = _ograph_canonical_json(raw)
    geometry_serialization = _ograph_canonical_json(geometry)
    (projected=projected, raw=raw, serialization=serialization,
     raw_hash=_ograph_sha256(serialization),
     geometry_hash=_ograph_sha256(geometry_serialization),
     real_gauss=real_gauss, real_orientation=real_orientation)
end

"""
    ograph_search_state(link::VirtualLink) -> CanonicalState

Return the unweighted exact-presentation descriptor used as a geometric search vertex.

# Arguments
`link` must be coherent under the
[VirtualLink data model](../specs/virtual-link-data-model.md). Reconstruction
and canonical JSON errors propagate. Weighted inputs and either coefficient
mode are accepted; this function does not require closed-normality.

# Returns
A `CanonicalState` carrying identity hashes, schema version, raw serialization,
raw payload, Area 3, representation `"o-graph-search"`, and real-vertex
complexity. The payload retains raw Gauss data, orientation, real-only status
and strand types, with weights omitted.

# Guarantees
1. The result uses the exact unweighted presentation identity defined by the search specification. [proved/agent: by construction]
2. Changing assigned weights alone leaves that identity unchanged. [proved/agent: by construction]
3. Neither the link nor a database is modified. [proved/agent: by construction]

# Not guaranteed
The descriptor is not canonical under re-rooting, relabelling or moves and
does not identify a manifold. Weight validity and transport require separate
replay; this function does not persist the descriptor.

# Sources
cf:docs/specs/ograph-search-db.md#contract — unweighted exact presentations and the separate geometry lookup.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); a = ograph_search_state(v);

julia> set_weight!(v, [1,-1,-1,0,1,0,-1,0]);

julia> ograph_search_state(v).state_id == a.state_id
true
```

Role: convert
"""
function ograph_search_state(link::VirtualLink)
    data = _ograph_search_state_data(link)
    CanonicalState(data.serialization; area=AREA_3,
        representation="o-graph-search", payload=data.raw,
        canonical_version=OGRAPH_SEARCH_SCHEMA_VERSION,
        complexity=length(data.real_orientation))
end

"""
    put_ograph_search_state!(store::MoveGraphStoreDB, link::VirtualLink) -> CanonicalState

Persist the unweighted exact presentation used as an o-graph search vertex.

# Arguments
`store` must be open and initialized by [`migrate_ograph_search_store!`](@ref).
An absent or unsupported schema is rejected; database errors can precede the
schema diagnostic. `link` must satisfy the
[VirtualLink data model](../specs/virtual-link-data-model.md); reconstruction
and serialization errors propagate. Weighted inputs are accepted and projected.

# Returns
The `CanonicalState` returned by `ograph_search_state(link)`, with identity,
serialization, representation metadata and real-vertex complexity for the
[search presentation](../specs/ograph-search-db.md#contract).

# Guarantees
1. Successful return registers the generic state and its search lookup row.
   Repeating the insertion reuses their identities. The link is unchanged.
   [proved/agent: by construction]
2. Input weights do not affect the search identity. Raw diagram data remain
   distinct from the real-Gauss geometry lookup, neither of which identifies
   a manifold. [proved/agent: by construction]

# Not guaranteed
This does not certify closed-normality, insert edges or preserve weights for
later loading. It does not verify or repair preexisting search metadata, and
the generic and search-specific writes are not one atomic transaction.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — exact unweighted state identity and secondary geometry lookup.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s3_abalone");

julia> s = put_ograph_search_state!(store, v);

julia> load_ograph_search_state(store, s.state_id).is_weighted
false

julia> close_store!(store);
```

See also [`load_ograph_search_state`](@ref), [`expand_ograph_search_state!`](@ref).

Role: other
"""
function put_ograph_search_state!(store::MoveGraphStoreDB, link::VirtualLink)
    _require_ograph_search_store(store)
    data = _ograph_search_state_data(link)
    state = put_state!(store, ograph_search_state(link))
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_states(
          state_id,state_version,raw_hash,geometry_hash,raw_gauss_json,
          orientation_json,real_gauss_json,real_orientation_json,real_vertex_count)
        VALUES (?,?,?,?,?,?,?,?,?)
    """, (state.state_id, OGRAPH_SEARCH_STATE_VERSION, data.raw_hash,
          data.geometry_hash, _ograph_canonical_json(data.projected.gauss),
          _ograph_canonical_json(data.projected.orientation),
          _ograph_canonical_json(data.real_gauss),
          _ograph_canonical_json(data.real_orientation),
          length(data.real_orientation)))
    state
end

"""
    load_ograph_search_state(store::MoveGraphStoreDB, state_id::AbstractString) -> Union{Nothing,VirtualLink}

Reconstruct the unweighted presentation stored under an o-graph search-state identity.

# Arguments
`store` must have the search schema and readable rows. `state_id` is an exact
stored identity. An absent joined search-state row returns `nothing`;
missing tables, malformed JSON or invalid reconstruction data propagate
database, parsing, conversion or `set_data!` errors.

# Returns
A fresh unweighted `VirtualLink` in integer weight mode, or `nothing`.
The returned Gauss data, vertex signs, strand types and `real_only` view
come from the stored search serialization.

# Guarantees
1. Rows written by [`put_ograph_search_state!`](@ref) reconstruct their
   unweighted presentation, including raw virtual visits and cut markers.
   [proved/agent: by construction]
2. Modifying the returned link does not update the database. Loading does
   not modify stored rows. [proved/agent: by construction]
3. This recovers the [search presentation](../specs/ograph-search-db.md#contract),
   not a manifold identity or the weights of an earlier weighted input.
   [proved/agent: by construction]

# Not guaranteed
The call does not replay the stored hash, verify the entire database or
certify closed-normality. Manually corrupted rows may fail or reconstruct
data that disagree with their claimed identity.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — exact raw presentations with weights removed.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store);

julia> s = put_ograph_search_state!(store, vl_example("s3_abalone"));

julia> load_ograph_search_state(store, s.state_id).is_weighted
false

julia> close_store!(store);
```

See also [`put_ograph_search_state!`](@ref), [`verify_ograph_search_store`](@ref).

Role: access
"""
function load_ograph_search_state(store::MoveGraphStoreDB, state_id::AbstractString)
    row = MoveGraphStore._one(store.db, """
        SELECT state.serialization FROM states AS state
        JOIN ograph_search_states AS search_state USING(state_id)
        WHERE state_id=?
    """, (String(state_id),))
    row === nothing && return nothing
    data = JSON.parse(String(row.serialization))
    link = VirtualLink()
    gauss = [Int.(component) for component in data["gauss"]]
    orientation = Int.(data["orientation"])
    strand_types = only.(String.(data["strand_types"]))
    if Bool(data["real_only"])
        set_data!(link; rgauss=gauss, orientation=orientation,
                  strand_types=strand_types)
    else
        set_data!(link; gauss=gauss, orientation=orientation,
                  strand_types=strand_types)
    end
    link
end

function _ograph_search_script_snapshot(script_path::AbstractString)
    path = abspath(String(script_path))
    isfile(path) || throw(ArgumentError("provenance script does not exist: $path"))
    path, bytes2hex(SHA.sha256(read(path)))
end

"""
    begin_ograph_search_run!(store::MoveGraphStoreDB, script_path::AbstractString; label::AbstractString, scope=Dict{String,Any}()) -> String

Register a running search with the executing script's provenance snapshot and return the new run identifier.

# Arguments
`store` must be open and migrated with [`migrate_ograph_search_store!`](@ref).
An absent or unsupported schema is rejected; an absent metadata table can
raise a database error before the schema diagnostic. `script_path` must name
an existing readable file, otherwise `ArgumentError` identifies the absent
path or the file-read error propagates. `label` is stored as text and `scope`
must be accepted by [`ograph_canonical_json`](@ref). Database and serialization
errors propagate.

# Returns
A `String` run identifier for the newly inserted search-run record.

# Guarantees
1. The stored run has status `"running"`, the supplied label and serialized
   scope, empty reported counts, and the absolute script path and digest
   of the bytes read by this call. [proved/agent: by construction]
2. The run records computational provenance for the
   [unweighted search graph](../specs/ograph-search-db.md#contract).
   Registration neither starts a search nor establishes a relation between
   represented spaces. [proved/agent: by construction]

# Not guaranteed
The script need not actually be executing and its contents are not executed
by this call. Run identities depend on the current time; uniqueness for
identical calls within the same timestamp is not promised. Scope is recorded
metadata, not an enforced budget or a mathematical source citation.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — script provenance and the unweighted search scope.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl");

julia> begin_ograph_search_run!(store, path; label="literal-02-smoke"); begin_ograph_search_run!(store, path; label="literal-02-reuse");

julia> verify_ograph_search_store(store).runs
2

julia> close_store!(store);
```

See also [`ograph_search_run`](@ref), [`finish_ograph_search_run!`](@ref).

Role: other
"""
function begin_ograph_search_run!(store::MoveGraphStoreDB,
                                  script_path::AbstractString;
                                  label::AbstractString,
                                  scope=Dict{String,Any}())
    _require_ograph_search_store(store)
    path, digest = _ograph_search_script_snapshot(script_path)
    started = string(Dates.now(Dates.UTC))
    run_id = "ograph-search-run-" * _ograph_sha256(
        join((replace(path, '\\' => '/'), digest, String(label), started), '\0'))
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO ograph_search_runs(
          run_id,script_path,script_sha256,label,scope_json,status,counts_json,
          started_at_utc)
        VALUES (?,?,?,?,?,'running','{}',?)
    """, (run_id, replace(path, '\\' => '/'), digest, String(label),
          _ograph_canonical_json(scope), started))
    run_id
end


"""
    finish_ograph_search_run!(store::MoveGraphStoreDB, run_id::AbstractString; status::AbstractString, counts=Dict{String,Any}()) -> Union{Nothing,NamedTuple}

Record a terminal outcome for a running o-graph search and return its summary.

# Arguments
`store` must have the search schema. `status` must be `"completed-found"`,
`"completed-none"`, `"budget-exhausted"` or `"failed"`; another value throws
`ArgumentError` identifying an invalid terminal status. `counts` must be
accepted by [`ograph_canonical_json`](@ref). Database and serialization errors
propagate. `run_id` may name an absent or already finished run.

# Returns
The [`ograph_search_run`](@ref) summary for `run_id`, or `nothing` when absent.
Its `counts` field is reported metadata, separate from the expansion totals.

# Guarantees
1. If the named run is running, its status, reported counts and finish time
   are updated. Other runs are unchanged. [proved/agent: by construction]
2. An already terminal run is not overwritten; an absent run is not created.
   [proved/agent: by construction]
3. The status records the caller's outcome for the
   [unweighted search workflow](../specs/ograph-search-db.md#workflow).
   It does not independently prove existence or nonexistence of a move path.
   [proved/agent: by construction]

# Not guaranteed
Counts and completion claims are not validated against a search budget.
This call does not cancel jobs, expand states or certify geometric equivalence.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — search-run records and independently stored expansions.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke");

julia> summary = finish_ograph_search_run!(store, r; status="completed-none");

julia> summary.candidates == summary.successes
true

julia> close_store!(store);
```

See also [`begin_ograph_search_run!`](@ref), [`ograph_search_run`](@ref).

Role: other
"""
function finish_ograph_search_run!(store::MoveGraphStoreDB,
                                   run_id::AbstractString;
                                   status::AbstractString,
                                   counts=Dict{String,Any}())
    status in ("completed-found", "completed-none", "budget-exhausted", "failed") ||
        throw(ArgumentError("invalid terminal o-graph search status"))
    MoveGraphStore._dbexec(store.db, """
        UPDATE ograph_search_runs
        SET status=?,counts_json=?,finished_at_utc=?
        WHERE run_id=? AND status='running'
    """, (String(status), _ograph_canonical_json(counts),
          string(Dates.now(Dates.UTC)), String(run_id)))
    ograph_search_run(store, run_id)
end

"""
    ograph_search_run(store::MoveGraphStoreDB, run_id::AbstractString) -> Union{Nothing,NamedTuple}

Read provenance, reported outcome and expansion totals for an o-graph search run.

# Arguments
`store` must contain readable search tables. An absent `run_id` returns
`nothing`. Missing tables, malformed stored JSON and database errors propagate.

# Returns
`nothing`, or a `NamedTuple` with `run_id`, script path and digest, `label`,
parsed `scope`, `status`, parsed reported `counts`, and integer totals
`expansions`, `candidates`, `successes`, `rejects`, `errors` for linked expansions.

# Guarantees
1. Totals include reused expansions linked to this run. Reported `counts`
   remain separate from those totals. Stored rows are unchanged.
   [proved/agent: by construction]
2. The result describes a finite execution record in the
   [search workflow](../specs/ograph-search-db.md#workflow), not a proof that
   all presentations of a manifold have been searched.
   [proved/agent: by construction]

# Not guaranteed
The script is not reread, stored outcomes are not replayed, and reported
scope or counts are not independently validated. Separate operations can
reach the same target; `successes` is not a distinct-target count.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — run and expansion records.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke");

julia> summary = ograph_search_run(store, r); summary.candidates == summary.successes
true

julia> close_store!(store);
```

See also [`finish_ograph_search_run!`](@ref), [`verify_ograph_search_store`](@ref).

Role: access
"""
function ograph_search_run(store::MoveGraphStoreDB, run_id::AbstractString)
    row = MoveGraphStore._one(store.db, """
        SELECT run.*,
               COUNT(DISTINCT link.expansion_id) AS expansions,
               COALESCE(SUM(expansion.candidate_count),0) AS candidates,
               COALESCE(SUM(expansion.success_count),0) AS successes,
               COALESCE(SUM(expansion.reject_count),0) AS rejects,
               COALESCE(SUM(expansion.error_count),0) AS errors
        FROM ograph_search_runs AS run
        LEFT JOIN ograph_search_run_expansions AS link ON link.run_id=run.run_id
        LEFT JOIN ograph_search_expansions AS expansion
          ON expansion.expansion_id=link.expansion_id
        WHERE run.run_id=? GROUP BY run.run_id
    """, (String(run_id),))
    row === nothing && return nothing
    (run_id=String(row.run_id), script_path=String(row.script_path),
     script_sha256=String(row.script_sha256), label=String(row.label),
     scope=JSON.parse(String(row.scope_json)), status=String(row.status),
     counts=JSON.parse(String(row.counts_json)), expansions=Int(row.expansions),
     candidates=Int(row.candidates), successes=Int(row.successes),
     rejects=Int(row.rejects), errors=Int(row.errors))
end

function _ograph_search_candidate_data(candidate::MoveCandidate)
    Dict(
        "type" => candidate.type,
        "parameter" => candidate.parameter,
        "vertices" => candidate.vertices,
        "edges" => [collect(edge) for edge in candidate.edges],
        "direction" => String(candidate.direction),
        "locations" => [collect(location) for location in candidate.locations],
        "ps_split_k" => candidate.ps_split_k === nothing ? nothing :
                        string(candidate.ps_split_k),
    )
end

_ograph_search_candidate_key(candidate::MoveCandidate) =
    _ograph_canonical_json(_ograph_search_candidate_data(candidate))

function _ograph_search_unique_candidates(candidates)
    seen = Set{String}()
    result = MoveCandidate[]
    for candidate in candidates
        key = _ograph_search_candidate_key(candidate)
        key in seen && continue
        push!(seen, key)
        push!(result, candidate)
    end
    result
end

function _ograph_search_mp_candidates(link::VirtualLink, selector)
    arities = haskey(selector, "vertex_arities") ?
        sort!(unique(Int.(selector["vertex_arities"]))) : [2, 3]
    all(arity -> arity in (2, 3), arities) ||
        throw(ArgumentError("MP vertex_arities must be drawn from [2,3]"))
    result = MoveCandidate[]
    n = length(link.orientation)
    for arity in arities,
        vertices in Iterators.product(ntuple(_ -> 1:n, arity)...)
        length(unique(vertices)) == arity || continue
        append!(result, movable(link, "MP"; vertices=collect(vertices)))
    end
    _ograph_search_unique_candidates(result)
end

function _ograph_search_02_candidates(link::VirtualLink, selector)
    direction = get(selector, "direction", "both")
    direction in ("forward", "inverse", "both") ||
        throw(ArgumentError("02 direction must be forward, inverse, or both"))
    result = MoveCandidate[]
    if direction in ("forward", "both")
        code, _ = real_gauss_code(link)
        for word in code, index in eachindex(word)
            isempty(word) && continue
            vertices = [word[index], word[mod1(index + 1, length(word))]]
            append!(result, movable(link, "02"; param="dir", vertices=vertices))
        end
    end
    if direction in ("inverse", "both")
        n = length(link.orientation)
        for first_vertex in 1:n, second_vertex in 1:n
            first_vertex == second_vertex && continue
            append!(result, movable(link, "02"; param="inv",
                                    vertices=[first_vertex, second_vertex]))
        end
    end
    _ograph_search_unique_candidates(result)
end

function _ograph_search_ps_candidates(link::VirtualLink, selector)
    parameters = haskey(selector, "parameters") ?
        String.(selector["parameters"]) : ["1", "2", "3", "4"]
    all(parameter -> parameter in ("1", "2", "3", "4"), parameters) ||
        throw(ArgumentError("PS parameters must be drawn from 1,2,3,4"))
    direction = get(selector, "direction", "both")
    direction in ("forward", "inverse", "both") ||
        throw(ArgumentError("PS direction must be forward, inverse, or both"))
    result = MoveCandidate[]
    if direction in ("inverse", "both")
        n = length(link.orientation)
        for parameter in parameters, first_vertex in 1:n, second_vertex in 1:n
            first_vertex == second_vertex && continue
            append!(result, movable(link, "PS"; param=parameter,
                                    vertices=[first_vertex, second_vertex]))
        end
    end
    if direction in ("forward", "both")
        code, _ = real_gauss_code(link)
        locations = [(component, position) for (component, word) in enumerate(code)
                     for position in eachindex(word)]
        for parameter in parameters, first_location in locations,
            second_location in locations
            first_location == second_location && continue
            append!(result, movable(link, "PS"; param=parameter,
                                    edge_occurrences=[first_location, second_location]))
        end
    end
    _ograph_search_unique_candidates(result)
end

"""
    ograph_search_candidates(link::VirtualLink, move_kind::AbstractString; selector=Dict{String,Any}()) -> Vector{MoveCandidate}

Enumerate selected unweighted local move candidates on a virtual-link presentation.

# Arguments
`link` must satisfy the [VirtualLink data model](../specs/virtual-link-data-model.md).
Reconstruction errors propagate. `move_kind` is `"MP"`, `"02"` or `"PS"`;
unsupported kinds throw `ArgumentError`, with a replay-only diagnostic for H.
`selector` is a key-value collection with string-convertible keys. MP accepts
`vertex_arities`; 02 accepts `direction`; PS accepts `parameters` and
`direction`. Invalid recognized values throw `ArgumentError` naming the
selector. Unrecognized keys are ignored. Closed-normality is not required or
certified. Defaults and selector meanings are fixed in the
[selector specification](../specs/ograph-search-db.md#selectors-and-query-boundaries).

# Returns
A `Vector{MoveCandidate}` for the selected move family. Candidates have the
meaning of the [move API](../specs/virtual-link-moves.md) on the
[unweighted search projection](../specs/ograph-search-db.md#contract).

# Guarantees
1. Candidates are independent of input weights and have the requested move
   kind. Duplicate candidate records are removed. The link is unchanged.
   [proved/agent: by construction]
2. Every candidate is a locally matched unweighted operation; applying it can
   still fail. This call neither executes moves nor inserts graph edges.
   [proved/agent: by construction]

# Not guaranteed
This is not a closed-normal move certificate, a weighted admissibility test,
or a complete search of manifold equivalences. PS and literal 02 remain
distinct. A selector can deliberately restrict the neighborhood.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Contract — MP, literal 02 and PS candidates on unweighted states.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> candidates = ograph_search_candidates(vl_example("s2xs1"), "02"; selector=Dict("direction" => "forward"));

julia> all(c.type == "02" for c in candidates)
true
```

See also [`expand_ograph_search_state!`](@ref), [`verify_ograph_move`](@ref).

Role: access
"""
function ograph_search_candidates(link::VirtualLink, move_kind::AbstractString;
                                  selector=Dict{String,Any}())
    projected = _ograph_search_unweighted(link)
    kind = String(move_kind)
    normalized = Dict(String(key) => value for (key, value) in pairs(selector))
    kind == "MP" && return _ograph_search_mp_candidates(projected, normalized)
    kind == "02" && return _ograph_search_02_candidates(projected, normalized)
    kind == "PS" && return _ograph_search_ps_candidates(projected, normalized)
    kind == "H" && throw(ArgumentError(
        "H is not a geometric search edge; insert it only during weighted replay"))
    throw(ArgumentError("o-graph search supports MP, literal 02, and PS"))
end

function _ograph_search_link_expansion!(store, run_id, expansion_id, reused)
    ordinal = Int(MoveGraphStore._one(store.db, """
        SELECT COALESCE(MAX(ordinal),0)+1 AS ordinal
        FROM ograph_search_run_expansions WHERE run_id=?
    """, (String(run_id),)).ordinal)
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_run_expansions(
          run_id,expansion_id,ordinal,reused) VALUES (?,?,?,?)
    """, (String(run_id), String(expansion_id), ordinal, reused ? 1 : 0))
end

function _ograph_search_operations(store, expansion_id)
    MoveGraphStore._rows(store.db, """
        SELECT operation.* FROM ograph_search_expansion_operations AS link
        JOIN ograph_search_operations AS operation USING(operation_id)
        WHERE link.expansion_id=? ORDER BY link.ordinal
    """, (String(expansion_id),))
end

"""
    expand_ograph_search_state!(store::MoveGraphStoreDB, run_id::AbstractString, source::VirtualLink, move_kind::AbstractString; selector=Dict{String,Any}()) -> Vector{NamedTuple}

Store a finite unweighted o-graph move neighborhood and return its operation outcome records.

# Arguments
- `store`: an open migrated search database; database errors propagate.
- `run_id`: an existing run in that database, required by the foreign key.
- `source`: a coherent presentation accepted by [`ograph_search_state`](@ref).
  Its assigned weights are discarded in the stored search copy.
- `move_kind`: `"MP"`, `"02"` or `"PS"`; `ArgumentError` rejects other kinds
  and identifies `"H"` as replay-only.
- `selector`: the finite scope passed to [`ograph_search_candidates`](@ref).
  Invalid selectors propagate that function's errors. Reusing a failed
  expansion or resuming one owned by another running run throws `ArgumentError`.

# Returns
The stored operation rows in candidate order, with `operation_id`,
`source_state_id`, optional `target_state_id`, `move_kind`, `direction`,
`engine_version`, `candidate_json`, `outcome` and optional `error_json`.
Outcomes distinguish successful, rejected and unexpectedly failed applications.

# Guarantees
1. Successful return records every candidate of the selected finite scope, or reuses its already completed expansion. [proved/agent: by construction]
2. Successful rows are directed edges between the exact unweighted presentations specified by the search database contract. [proved/agent: by construction]
3. The database and run-expansion association are updated; `source` is not modified. [proved/agent: by construction]

# Not guaranteed
Completion of this selected neighborhood is not completeness of a move class
or a manifold-equivalence proof. Failed candidate applications are returned as
rows, not necessarily thrown. Enumeration failure can leave a stored failed
expansion, and source/run records can survive a later failure.

# Sources
cf:docs/specs/ograph-search-db.md#contract — presentation-level state identity, finite expansions and reuse.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); run = begin_ograph_search_run!(store, path; label="literal-02-smoke");

julia> v = VirtualLink(); set_data!(v; gauss=[[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]], orientation=[1,1,1,1,0,0]);

julia> edges = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> (!isempty(edges), all(String(e.outcome) == "success" for e in edges))
(true, true)

julia> close_store!(store);
```

Role: other
"""
function expand_ograph_search_state!(store::MoveGraphStoreDB,
                                     run_id::AbstractString,
                                     source::VirtualLink,
                                     move_kind::AbstractString;
                                     selector=Dict{String,Any}())
    kind = String(move_kind)
    kind == "H" && throw(ArgumentError(
        "H is replay-only and cannot enter the unweighted graph"))
    kind in ("MP", "02", "PS") || throw(ArgumentError(
        "o-graph search supports MP, literal 02, and PS"))
    source_link = _ograph_search_unweighted(source)
    source_state = put_ograph_search_state!(store, source_link)
    selector_json = _ograph_canonical_json(
        Dict(String(key) => value for (key, value) in pairs(selector)))
    expansion_id = "ograph-search-expansion-" * _ograph_sha256(join((
        source_state.state_id, kind, OGRAPH_SEARCH_ENUMERATOR_VERSION,
        selector_json), '\0'))
    existing = MoveGraphStore._one(store.db,
        "SELECT * FROM ograph_search_expansions WHERE expansion_id=?",
        (expansion_id,))
    if existing !== nothing && String(existing.status) == "completed"
        _ograph_search_link_expansion!(store, run_id, expansion_id, true)
        return _ograph_search_operations(store, expansion_id)
    elseif existing !== nothing
        String(existing.status)=="building" || throw(ArgumentError(
            "existing search expansion failed: $expansion_id"))
        active=MoveGraphStore._one(store.db,"""
            SELECT COUNT(*) AS n FROM ograph_search_run_expansions AS link
            JOIN ograph_search_runs AS run USING(run_id)
            WHERE link.expansion_id=? AND run.status='running' AND run.run_id<>?
        """,(expansion_id,String(run_id)))
        Int(active.n)==0 || throw(ArgumentError(
            "search expansion is owned by another running run: $expansion_id"))
        # A terminated run can leave the row created before enumeration. The
        # actual edge batch is transactional; resume without discarding rows.
    end

    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_expansions(
          expansion_id,source_state_id,move_kind,enumerator_version,selector_json,
          status,candidate_count,success_count,reject_count,error_count)
        VALUES (?,?,?,?,?,'building',0,0,0,0)
    """, (expansion_id, source_state.state_id, kind,
          OGRAPH_SEARCH_ENUMERATOR_VERSION, selector_json))
    _ograph_search_link_expansion!(store, run_id, expansion_id, false)

    candidates = try
        ograph_search_candidates(source_link, kind; selector=selector)
    catch error_value
        error_json = _ograph_canonical_json(Dict(
            "type" => string(typeof(error_value)),
            "message" => sprint(showerror, error_value)))
        MoveGraphStore._dbexec(store.db, """
            UPDATE ograph_search_expansions
            SET status='failed',error_count=1,error_json=? WHERE expansion_id=?
        """, (error_json, expansion_id))
        rethrow()
    end

    success_count = 0
    reject_count = 0
    error_count = 0
    MoveGraphStore._transaction(store.db) do
        for (ordinal, candidate) in enumerate(candidates)
            candidate_data = _ograph_search_candidate_data(candidate)
            candidate_json = _ograph_canonical_json(candidate_data)
            operation_id = "ograph-search-operation-" * _ograph_sha256(join((
                source_state.state_id, OGRAPH_SEARCH_ENGINE_VERSION,
                candidate_json), '\0'))
            target_state_id = nothing
            outcome = "success"
            error_json = nothing
            try
                target = move(source_link, candidate)
                target_state_id = put_ograph_search_state!(store, target).state_id
                success_count += 1
            catch error_value
                expected = error_value isa Union{ArgumentError,DimensionMismatch,BoundsError}
                outcome = expected ? "rejected" : "error"
                error_json = _ograph_canonical_json(Dict(
                    "type" => string(typeof(error_value)),
                    "message" => sprint(showerror, error_value)))
                expected ? (reject_count += 1) : (error_count += 1)
            end
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_search_operations(
                  operation_id,source_state_id,target_state_id,move_kind,direction,
                  engine_version,candidate_json,outcome,error_json)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, (operation_id, source_state.state_id,
                  target_state_id === nothing ? missing : target_state_id,
                  kind, String(candidate.direction), OGRAPH_SEARCH_ENGINE_VERSION,
                  candidate_json, outcome,
                  error_json === nothing ? missing : error_json))
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_search_expansion_operations(
                  expansion_id,ordinal,operation_id) VALUES (?,?,?)
            """, (expansion_id, ordinal, operation_id))
        end
        MoveGraphStore._dbexec(store.db, """
            UPDATE ograph_search_expansions
            SET status='completed',candidate_count=?,success_count=?,reject_count=?,
                error_count=?,error_json=NULL WHERE expansion_id=?
        """, (length(candidates), success_count, reject_count, error_count,
              expansion_id))
    end
    _ograph_search_operations(store, expansion_id)
end

"""
    ograph_search_neighbors(store::MoveGraphStoreDB, state_id::AbstractString; kinds=("MP", "02")) -> Vector{NamedTuple}

Read the stored successful outgoing moves from one unweighted search state.

# Arguments
`store` must have readable search tables. `kinds` contains `"MP"`, `"02"`
or `"PS"`; other values throw `ArgumentError` naming the supported kinds.
An empty collection returns an empty vector without reading the database.
An absent source identity has no neighbors. Database errors propagate.

# Returns
A vector of rows with `operation_id`, source and target state IDs,
`move_kind`, `direction`, `engine_version` and `candidate_json`. The latter
is serialized candidate data, not a `MoveCandidate` instance.

# Guarantees
1. Only stored successful operations of allowed kinds with the exact source
   ID are returned, sorted by move kind and operation ID. Duplicate kinds in
   the query have no effect. Stored rows are unchanged.
   [proved/agent: by construction]
2. These are directed edges of the
   [unweighted search graph](../specs/ograph-search-db.md#contract), not
   weighted moves or independently certified manifold equivalences.
   [proved/agent: by construction]

# Not guaranteed
No neighborhood is generated and no edge is replayed. An empty result does
not establish that no move exists. Multiple operations may share a target.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — persisted successful operations.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke");

julia> expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> all(String(e.move_kind) == "02" for e in ograph_search_neighbors(store, s.state_id; kinds=("02",)))
true

julia> close_store!(store);
```

See also [`expand_ograph_search_state!`](@ref), [`query_ograph_search_path`](@ref).

Role: access
"""
function ograph_search_neighbors(store::MoveGraphStoreDB,
                                 state_id::AbstractString;
                                 kinds=("MP", "02"))
    allowed = unique(String.(collect(kinds)))
    isempty(allowed) && return NamedTuple[]
    all(kind -> kind in ("MP", "02", "PS"), allowed) ||
        throw(ArgumentError("search-neighbor kinds must be MP, 02, or PS"))
    placeholders = join(fill("?", length(allowed)), ",")
    MoveGraphStore._rows(store.db, """
        SELECT * FROM successful_ograph_search_edges
        WHERE source_state_id=? AND move_kind IN ($placeholders)
        ORDER BY move_kind,operation_id
    """, Tuple(Any[String(state_id); allowed]))
end

function _ograph_search_path_rows(store, operation_ids)
    result = NamedTuple[]
    for operation_id in operation_ids
        row = MoveGraphStore._one(store.db, """
            SELECT * FROM successful_ograph_search_edges WHERE operation_id=?
        """, (String(operation_id),))
        row === nothing && throw(ArgumentError("search path operation is absent"))
        push!(result, (operation_id=String(row.operation_id),
                       source_state_id=String(row.source_state_id),
                       target_state_id=String(row.target_state_id),
                       move_kind=String(row.move_kind),
                       direction=String(row.direction),
                       engine_version=String(row.engine_version),
                       candidate=JSON.parse(String(row.candidate_json))))
    end
    result
end

"""
    query_ograph_search_paths(store::MoveGraphStoreDB, source_state_id::AbstractString; target_state_id=nothing, target_geometry_hash=nothing, kinds=("MP","02"), kind_pattern=nothing, min_depth::Integer=1, max_depth::Integer=12, limit::Integer=1000) -> Union{Nothing,Vector{Vector{NamedTuple}}}

Read bounded directed paths in the stored unweighted presentation graph.

# Arguments
- `store`: an open migrated search store; database/decoding errors propagate.
- `source_state_id`: exact stored source identity.
- Exactly one of `target_state_id` and `target_geometry_hash` must be supplied;
  otherwise `ArgumentError` identifies the target requirement.
- `kinds`: admitted MP, literal 02 or PS move names. Without a pattern,
  unsupported kinds throw `ArgumentError`, and an empty collection returns `nothing`.
- `kind_pattern`: optional sequence of those move names; it overrides the
  depth bounds after their initial validation and fixes the path length and
  move at each step. Unsupported entries throw `ArgumentError`.
- `min_depth`, `max_depth`, `limit`: initially require nonnegative minimum,
  maximum at least minimum, and positive limit. Violations throw `ArgumentError`;
  the limit must convert to `Int`.

# Returns
A vector of paths, ordered by length then the stored operation-ID sequence.
Each path is a vector of rows containing `operation_id`, `source_state_id`,
`target_state_id`, `move_kind`, `direction`, `engine_version` and decoded
`candidate`. No matches give an empty vector; empty kinds have the separate
`nothing` return noted above.

# Guarantees
1. Every returned step is a stored successful directed unweighted operation and its endpoints compose with adjacent steps. [proved/agent: by construction]
2. Paths obey the target and selected finite scope, have no repeated state except a possible terminal return to the source, and number at most `limit`. [proved/agent: by construction]
3. The database is not modified. [proved/agent: by construction]

# Not guaranteed
These are paths in the stored graph, not an exhaustive search of all moves or
all manifold equivalences. Stored operations are not replayed here; weights
and framings require separate replay. Reaching a geometry hash does not
certify an induced homology marking. An empty pattern is not a supported
zero-step query and can produce a database syntax error.

# Sources
cf:docs/specs/ograph-search-db.md — geometric path selection, source-return loops and separate weighted replay.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); sid = put_ograph_search_state!(store, v).state_id;

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward")); target = String(first(rows).target_state_id);

julia> path = query_ograph_search_path(store, sid; target_state_id=target, kinds=("02",), max_depth=1);

julia> query_ograph_search_paths(store, sid; target_state_id=target, kinds=("02",), max_depth=1, limit=10) == [path]
true

julia> close_store!(store);
```

Role: access
"""
function query_ograph_search_paths(store::MoveGraphStoreDB,
                                   source_state_id::AbstractString;
                                   target_state_id=nothing,
                                   target_geometry_hash=nothing,
                                   kinds=("MP", "02"),
                                   kind_pattern=nothing,
                                   min_depth::Integer=1,
                                   max_depth::Integer=12,
                                   limit::Integer=1000)
    (target_state_id === nothing) == (target_geometry_hash === nothing) &&
        throw(ArgumentError("specify exactly one path target"))
    min_depth >= 0 && max_depth >= min_depth ||
        throw(ArgumentError("invalid path depth bounds"))
    limit > 0 || throw(ArgumentError("path query limit must be positive"))
    allowed = unique(String.(collect(kinds)))
    pattern = kind_pattern === nothing ? nothing : String.(collect(kind_pattern))
    if pattern !== nothing
        all(kind -> kind in ("MP", "02", "PS"), pattern) ||
            throw(ArgumentError("path pattern contains an unsupported move"))
        max_depth = length(pattern)
        min_depth = length(pattern)
    else
        all(kind -> kind in ("MP", "02", "PS"), allowed) ||
            throw(ArgumentError("path kinds must be MP, 02, or PS"))
        isempty(allowed) && return nothing
    end
    edge_filter = if pattern === nothing
        "edge.move_kind IN (" * join(fill("?", length(allowed)), ",") * ")"
    else
        clauses = ["WHEN $(index - 1) THEN ?" for index in eachindex(pattern)]
        "edge.move_kind = CASE walk.depth " * join(clauses, " ") * " END"
    end
    goal_filter = target_state_id === nothing ?
        "target.geometry_hash=?" : "walk.state_id=?"
    sql = """
        WITH RECURSIVE walk(state_id,depth,operation_ids,visited) AS (
          SELECT ?,0,'','|' || ? || '|'
          UNION ALL
          SELECT edge.target_state_id,walk.depth+1,
                 CASE WHEN walk.operation_ids=''
                      THEN edge.operation_id
                      ELSE walk.operation_ids || ',' || edge.operation_id END,
                 walk.visited || edge.target_state_id || '|'
          FROM walk
          JOIN successful_ograph_search_edges AS edge
            ON edge.source_state_id=walk.state_id
          WHERE walk.depth < ? AND $edge_filter
            AND (walk.depth=0 OR walk.state_id<>?)
            AND (instr(walk.visited,'|' || edge.target_state_id || '|')=0
                 OR (edge.target_state_id=? AND walk.depth+1>=?))
        )
        SELECT walk.operation_ids,walk.depth FROM walk
        JOIN ograph_search_states AS target ON target.state_id=walk.state_id
        WHERE walk.depth>=? AND $goal_filter
        ORDER BY walk.depth,walk.operation_ids LIMIT ?
    """
    params = Any[String(source_state_id), String(source_state_id), max_depth]
    append!(params, pattern === nothing ? allowed : pattern)
    append!(params, [String(source_state_id),String(source_state_id),min_depth])
    push!(params, min_depth)
    push!(params, target_state_id === nothing ? String(target_geometry_hash) :
          String(target_state_id))
    push!(params, Int(limit))
    rows = MoveGraphStore._rows(store.db, sql, Tuple(params))
    [_ograph_search_path_rows(store,
        isempty(String(row.operation_ids)) ? String[] :
        split(String(row.operation_ids), ',')) for row in rows]
end

"""
    query_ograph_search_path(store::MoveGraphStoreDB, source_state_id::AbstractString; kwargs...) -> Union{Nothing,Vector{NamedTuple}}

Select one shortest stored directed move path to a specified search target.

# Arguments
`store` must contain the search schema. Supply exactly one of `target_state_id`
and `target_geometry_hash`, and the keyword options of
[`query_ograph_search_paths`](@ref): `kinds`, `kind_pattern`, `min_depth`,
`max_depth`, and optionally `limit`. Target and bound violations throw
`ArgumentError` with the corresponding diagnostic; SQL and stored-data errors
propagate. Defaults are depths 1 through 12 and kinds MP and literal 02.
An empty `kinds` collection without a pattern is unsupported here: the plural
query returns `nothing`, and this wrapper raises `MethodError` on it.

# Returns
`nothing` when no matching path is found, or a vector of step records with
operation ID, endpoint IDs, move kind, direction, engine version and parsed
candidate data. A permitted zero-length path is an empty vector, not `nothing`.

# Guarantees
1. The selected path uses stored successful edges and is first in the plural
   query's depth/operation-ID order. A pattern, when supplied, fixes its length
   and move kinds. The database is unchanged. [proved/agent: by construction]
2. The path is a walk in the [unweighted search graph](../specs/ograph-search-db.md#workflow),
   subject to that query's simple-path restriction and final-source return.
   It is not yet a weighted move word. [proved/agent: by construction]

# Not guaranteed
No new states are explored, no edges are replayed and no manifold equivalence
is certified. Failure to find a stored bounded path is not a nonexistence proof.
A geometry-hash target is not a manifold-identity target.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow and §Exact H-insertion feasibility — geometric path selection before weighted replay.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); edges = expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, s.state_id; target_state_id=String(first(edges).target_state_id), kinds=("02",), max_depth=1);

julia> (length(path), only(path).move_kind)
(1, "02")

julia> close_store!(store);
```

See also [`replay_weighted_ograph_path`](@ref), [`solve_ograph_path_weights`](@ref).

Role: access
"""
function query_ograph_search_path(store::MoveGraphStoreDB,
                                  source_state_id::AbstractString; kwargs...)
    paths = query_ograph_search_paths(store, source_state_id; limit=1, kwargs...)
    isempty(paths) ? nothing : first(paths)
end

function _ograph_search_match_candidate(candidate::MoveCandidate, data)
    candidate.type == String(data["type"]) || return false
    candidate.parameter == String(data["parameter"]) || return false
    candidate.vertices == Int.(data["vertices"]) || return false
    [collect(edge) for edge in candidate.edges] ==
        [Int.(edge) for edge in data["edges"]] || return false
    String(candidate.direction) == String(data["direction"]) || return false
    [collect(location) for location in candidate.locations] ==
        [Int.(location) for location in data["locations"]] || return false
    true
end

function _ograph_search_weighted_candidate(link, data)
    kind = String(data["type"])
    parameter = String(data["parameter"])
    vertices = Int.(data["vertices"])
    candidates = if kind == "PS" && String(data["direction"]) == "forward"
        locations = [Tuple(Int.(location)) for location in data["locations"]]
        movable(link, kind; param=parameter, edge_occurrences=locations, strict=true)
    else
        movable(link, kind; param=parameter, vertices=vertices, strict=true)
    end
    matching = filter(candidate -> _ograph_search_match_candidate(candidate, data),
                      candidates)
    length(matching) == 1 ? only(matching) : nothing
end

"""
    replay_weighted_ograph_path(store::MoveGraphStoreDB, weighted_source::VirtualLink, path; h_moves=NamedTuple[]) -> NamedTuple

Replay a selected geometric move word on assigned weights with optional H insertions.

# Arguments
- `weighted_source`: must have assigned weights, otherwise `ArgumentError`.
  Its search identity must match the first path source when nonempty;
  mismatch throws `ArgumentError`. Valid integer framing equations are a
  caller precondition for a framed interpretation, not a global check here.
- `path`: coherent geometric path rows. The caller is responsible for their
  provenance and endpoint correspondence; each supplied candidate must be
  available on the current weighted presentation.
- `h_moves`: records with `after_step`, `vertex`, `amount`; the step must be
  between zero and the path length, otherwise `ArgumentError`. Integer
  conversions and candidate-enumeration errors can propagate.
- `store`: retained for the public search API; this method does not read or
  modify its database.

# Returns
`(valid, failure_kind, failure_step, message, target, states)`. Success has
`valid=true`, a final `VirtualLink` target and `nothing` failure fields.
A reported failure has no target and retains snapshots up to that failure.
`states` contains the initial state and a snapshot after each successful move.

# Guarantees
1. A successful result is obtained by applying the specified weighted candidates and H insertions in order to a copy of the source. [proved/agent: by construction]
2. The result has the weighted path interpretation of the linked specification, subject to the caller's framing and endpoint preconditions. [proved/agent: by construction]
3. The source and path are not modified; snapshots have independent mutable data. [proved/agent: by construction]

# Not guaranteed
Success does not independently verify disk equations at every stage, compare
each result with the row's stored target ID, authenticate operation IDs, or
certify a prescribed final weight vector. Candidate enumeration can throw
instead of producing a failure record; numeric weight storage limitations apply.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — weighted replay after geometric selection and required endpoint/framed-residual checks.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); set_weight!(v, [1,-1,-1,0,1,0,-1,0]); sid = put_ograph_search_state!(store, v).state_id;

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, sid; target_state_id=String(first(rows).target_state_id), kinds=("02",), max_depth=1);

julia> replay_weighted_ograph_path(store, v, path).valid
true

julia> close_store!(store);
```

Role: move
"""
function replay_weighted_ograph_path(store::MoveGraphStoreDB,
                                     weighted_source::VirtualLink,
                                     path;
                                     h_moves=NamedTuple[])
    weighted_source.is_weighted || throw(ArgumentError(
        "weighted replay requires a weighted source"))
    steps = collect(path)
    source_state = ograph_search_state(weighted_source)
    if !isempty(steps) && String(first(steps).source_state_id) != source_state.state_id
        throw(ArgumentError("weighted source geometry does not match path source"))
    end
    grouped = Dict{Int,Vector{Any}}()
    for hmove in h_moves
        step = Int(hmove.after_step)
        0 <= step <= length(steps) || throw(ArgumentError("H after_step is outside path"))
        push!(get!(grouped, step, Any[]), hmove)
    end
    current = deepcopy(weighted_source)
    states = VirtualLink[deepcopy(current)]
    function apply_h!(step)
        for hmove in get(grouped, step, Any[])
            try
                move!(current, "H"; param=Int(hmove.amount),
                      vertices=[Int(hmove.vertex)], strict=true)
            catch error_value
                return (valid=false, failure_kind="H", failure_step=step,
                        message=sprint(showerror, error_value), target=nothing,
                        states=states)
            end
            push!(states, deepcopy(current))
        end
        nothing
    end
    failure = apply_h!(0)
    failure === nothing || return failure
    for (index, step) in enumerate(steps)
        candidate = _ograph_search_weighted_candidate(current, step.candidate)
        candidate === nothing && return (
            valid=false, failure_kind=String(step.move_kind), failure_step=index,
            message="unweighted path edge is not weight-admissible",
            target=nothing, states=states)
        try
            move!(current, candidate)
        catch error_value
            return (valid=false, failure_kind=String(step.move_kind),
                    failure_step=index, message=sprint(showerror, error_value),
                    target=nothing, states=states)
        end
        push!(states, deepcopy(current))
        failure = apply_h!(index)
        failure === nothing || return failure
    end
    (valid=true, failure_kind=nothing, failure_step=nothing, message=nothing,
     target=current, states=states)
end

function _ograph_search_weight_payload(link::VirtualLink)
    [[isinteger(value) ? Int(round(value)) : value for value in component]
     for component in link.weights]
end

"""
    record_weighted_ograph_path_trial!(store::MoveGraphStoreDB, run_id::AbstractString, weighted_source::VirtualLink, path; h_moves=NamedTuple[]) -> NamedTuple

Store a geometric move word's weighted replay outcome and return its trial identifier, path fingerprint and replay result.

# Arguments
- `store`, `run_id`: an open migrated search store and an existing run.
  Database errors propagate.
- `weighted_source`: assigned weights are required, otherwise `ArgumentError`.
- `path`, `h_moves`: satisfy [`replay_weighted_ograph_path`](@ref); its domain
  errors and representation/conversion failures propagate. The caller supplies
  authentic path records and integer-convertible H data.

# Returns
`(trial_id, path_fingerprint, replay)`, with string identifiers and the current
replay result, including any reported failure.

# Guarantees
1. Successful return stores an admissible or rejected weighted trial separately from the unweighted search graph. [proved/agent: by construction]
2. Repeating the same run/path/source-weight/H identity increments its observation count. [proved/agent: by construction]
3. The database is modified; source, path and H records are not modified. [proved/agent: by construction]

# Not guaranteed
On identity conflict the original stored outcome is retained, even though
the returned replay is freshly computed. Source insertion and trial recording
are not one transaction. Replay has the framing and endpoint limits of the
linked operation; persistence adds no mathematical certificate.

# Sources
cf:docs/specs/ograph-search-db.md — separate weighted trials and unweighted graph identity.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store); v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); set_weight!(v, [1,-1,-1,0,1,0,-1,0]); sid = put_ograph_search_state!(store, v).state_id;

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke"); rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, sid; target_state_id=String(first(rows).target_state_id), kinds=("02",), max_depth=1); trial = record_weighted_ograph_path_trial!(store, run, v, path);

julia> (trial.replay.valid, verify_ograph_search_store(store).weighted_trials)
(true, 1)

julia> close_store!(store);
```

Role: other
"""
function record_weighted_ograph_path_trial!(store::MoveGraphStoreDB,
                                            run_id::AbstractString,
                                            weighted_source::VirtualLink,
                                            path;
                                            h_moves=NamedTuple[])
    weighted_source.is_weighted || throw(ArgumentError(
        "weighted path trial requires a weighted source"))
    steps = collect(path)
    source_state = put_ograph_search_state!(store, weighted_source)
    operation_ids = String[String(step.operation_id) for step in steps]
    operations_json = _ograph_canonical_json(operation_ids)
    source_weights_json = _ograph_canonical_json(
        _ograph_search_weight_payload(weighted_source))
    h_data = [Dict("after_step" => Int(move.after_step),
                   "vertex" => Int(move.vertex),
                   "amount" => Int(move.amount)) for move in h_moves]
    h_json = _ograph_canonical_json(h_data)
    path_fingerprint = _ograph_sha256(operations_json)
    trial_id = "ograph-search-weighted-trial-" * _ograph_sha256(join((
        String(run_id), path_fingerprint, source_weights_json, h_json), '\0'))
    replay = replay_weighted_ograph_path(store, weighted_source, steps;
                                         h_moves=h_moves)
    outcome = replay.valid ? "admissible" : "rejected"
    target_weights_json = replay.valid ? _ograph_canonical_json(
        _ograph_search_weight_payload(replay.target)) : nothing
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO ograph_search_weighted_trials(
          trial_id,run_id,source_state_id,path_fingerprint,operation_ids_json,
          source_weights_json,h_moves_json,outcome,target_weights_json,
          failure_kind,failure_step,message,observed_count)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,1)
        ON CONFLICT(run_id,path_fingerprint,source_weights_json,h_moves_json)
        DO UPDATE SET observed_count=observed_count+1
    """, (trial_id, String(run_id), source_state.state_id, path_fingerprint,
          operations_json, source_weights_json, h_json, outcome,
          target_weights_json === nothing ? missing : target_weights_json,
          replay.failure_kind === nothing ? missing : String(replay.failure_kind),
          replay.failure_step === nothing ? missing : Int(replay.failure_step),
          replay.message === nothing ? missing : String(replay.message)))
    (trial_id=trial_id, path_fingerprint=path_fingerprint, replay=replay)
end

"""
    verify_ograph_search_store(store::MoveGraphStoreDB) -> NamedTuple

Return the schema version and search-record counts after checking the version and foreign-key integrity.

# Arguments
`store` must be open with readable search tables. Missing version metadata,
an unsupported version or a foreign-key violation throws `ArgumentError`
identifying that failure. Missing tables or malformed version text can raise
database or parsing errors first.

# Returns
A `NamedTuple` of integers: `schema_version`, `states`, `runs`, `expansions`,
`operations`, `successful_edges` and `weighted_trials`. Each count concerns
the corresponding [search record family](../specs/ograph-search-db.md#workflow).

# Guarantees
1. Successful return means the stored version matches this implementation
   and the database's foreign-key check reported no violations. The database
   is unchanged. [proved/agent: by construction]
2. Counts describe stored rows and successful directed edges, not distinct
   manifolds, complete move classes or independently verified calculations.
   [proved/agent: by construction]

# Not guaranteed
This function does not recompute identities, replay moves or weighted trials,
validate expansion completeness, or count constraint trials. The name
`verify` does not imply those stronger checks.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Workflow — separately stored runs, expansions, operations and weighted trials.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); begin_ograph_search_run!(store, path; label="literal-02-smoke"); begin_ograph_search_run!(store, path; label="literal-02-reuse");

julia> verify_ograph_search_store(store).runs
2

julia> close_store!(store);
```

See also [`ograph_search_run`](@ref), [`replay_weighted_ograph_path`](@ref).

Role: condition
"""
function verify_ograph_search_store(store::MoveGraphStoreDB)
    version = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    version !== nothing || throw(ArgumentError("o-graph search schema is absent"))
    parse(Int, String(version.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph search schema version"))
    isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check")) ||
        throw(ArgumentError("o-graph search foreign-key check failed"))
    count(table) = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM $table").n)
    (schema_version=OGRAPH_SEARCH_SCHEMA_VERSION,
     states=count("ograph_search_states"),
     runs=count("ograph_search_runs"),
     expansions=count("ograph_search_expansions"),
     operations=count("ograph_search_operations"),
     successful_edges=count("successful_ograph_search_edges"),
     weighted_trials=count("ograph_search_weighted_trials"))
end
