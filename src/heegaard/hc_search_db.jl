# cf: lib/topology.jl/docs/specs/module-structure.md §3.1 — hc_search_* mirrors ograph_search_* on MoveGraphStoreDB
# cf: lib/topology.jl/src/ograph_search_db.jl — run, state, and transition tables over the same store
# cf: lib/topology.jl/src/heegaard/hc_compare_canonical.jl — canonical key, or the literal writing when that key is not used

import JSON
import Dates
import ..MoveGraphStore

const HC_SEARCH_SCHEMA_VERSION = 1
const HC_SEARCH_STATE_VERSION = "hc-search-state/1"
const HC_SEARCH_ENGINE_VERSION = "hc-search-engine/1"
const _HC_SEARCH_IDENTITIES = ("canonical", "literal")
const _HC_SEARCH_TERMINAL = ("completed", "budget-exhausted", "failed")

function _hc_search_schema_path()
    normpath(joinpath(@__DIR__, "..", "..", "db", "hc_search_schema.sql"))
end

"""
    migrate_hc_search_store!(store::MoveGraphStoreDB) -> MoveGraphStoreDB

Install the HC-search tables and views in a move-graph store.

# Arguments
`store` must be an open writable store compatible with the base
`MoveGraphStore` schema. The installed HC-search schema file must be readable.
Database, file and schema errors propagate; an absent or unsupported
HC-search schema version raises `ArgumentError` after the migration attempt.

# Returns
The same `store`, ready for the HC-search record operations.

# Guarantees
1. Successful return means the HC-search schema-version check passes.
   [proved/agent: by construction]
2. The store is modified to install the base and HC-search schema objects;
   this operation constructs database structure, not a mathematical move
   witness. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
This does not validate stored diagrams or prove the mathematical claims of
search records. Base migration and HC migration are separate operations;
failure is not a promise that the entire store remains unchanged.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:");

julia> migrate_hc_search_store!(store) === store
true

julia> Topology.close_store!(store);
```

See also [`begin_hc_search_run!`](@ref), [`put_hc_search_state!`](@ref).

Role: other
"""
function migrate_hc_search_store!(store::MoveGraphStore.MoveGraphStoreDB)
    MoveGraphStore.migrate!(store)
    MoveGraphStore._transaction(store.db) do
        MoveGraphStore._execute_schema!(store.db, read(_hc_search_schema_path(), String))
    end
    _require_hc_search_store(store)
end

function _require_hc_search_store(store::MoveGraphStore.MoveGraphStoreDB)
    present = MoveGraphStore._one(store.db, """
        SELECT 1 AS present FROM sqlite_master
        WHERE type='table' AND name='hc_search_schema_meta'
    """)
    present === nothing && throw(ArgumentError(
        "HC search schema is absent; call migrate_hc_search_store!"))
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM hc_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == HC_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported HC search schema version"))
    store
end

_hc_identity(identity::AbstractString) = identity in _HC_SEARCH_IDENTITIES ?
    String(identity) : throw(ArgumentError("identity must be canonical or literal"))
_hc_identity(identity::Symbol) = _hc_identity(String(identity))

function _hc_stable_json(value)
    value === nothing && return "null"
    value isa Bool && return value ? "true" : "false"
    if value isa Integer
        integer = BigInt(value)
        typemin(Int64) <= integer <= typemax(Int64) ||
            throw(ArgumentError("search JSON integer does not fit Int64"))
        return string(Int64(integer))
    end
    value isa AbstractString && return JSON.json(String(value))
    value isa Symbol && return JSON.json(String(value))
    value isa Pair && return _hc_stable_json([first(value), last(value)])
    if value isa NamedTuple
        return _hc_stable_json(Dict{String,Any}(String(key) => getfield(value, key) for key in keys(value)))
    end
    if value isa AbstractDict
        parts = String[]
        for key in sort!(collect(keys(value)); by=string)
            push!(parts, JSON.json(string(key)) * ":" * _hc_stable_json(value[key]))
        end
        return "{" * join(parts, ",") * "}"
    end
    if value isa Union{AbstractVector,Tuple}
        return "[" * join((_hc_stable_json(item) for item in value), ",") * "]"
    end
    throw(ArgumentError("not a search JSON value: $(typeof(value))"))
end

function _hc_pair_list(dict)
    [[string(key), dict[key]] for key in sort!(collect(keys(dict)); by=string)]
end

function _hc_writing(hc::HeegaardCode)
    lc = hc.lc
    order = hc_component_names(lc)
    Dict{String,Any}(
        "schema" => "hc-search-writing/1",
        "genus" => lc.genus,
        "words" => [Dict{String,Any}("name" => name, "word" => copy(lc.words[name]))
                    for name in order if haskey(lc.words, name)],
        "places" => _hc_pair_list(lc.places),
        "lc_disks" => _hc_pair_list(lc.disks),
        "lc_points" => _hc_pair_list(lc.points),
        "lc_weights" => _hc_pair_list(lc.weights),
        "roles" => [Dict{String,Any}("name" => name, "role" => String(hc.roles[name]))
                    for name in order if haskey(hc.roles, name)],
        "disks" => _hc_pair_list(hc.disks),
        "points" => _hc_pair_list(hc.points),
        "weights" => _hc_pair_list(hc.weights),
    )
end

function _hc_disk_dict(items)
    Dict{String,Vector{String}}(String(item[1]) => String.(item[2]) for item in items)
end

function _hc_string_dict(items)
    Dict{String,String}(String(item[1]) => String(item[2]) for item in items)
end

function _hc_int_dict(items)
    Dict{String,Int}(String(item[1]) => Int(item[2]) for item in items)
end

function _hc_from_writing(data)
    String(data["schema"]) == "hc-search-writing/1" ||
        throw(ArgumentError("unknown HC search writing"))
    words = Dict{String,Vector{Int}}(
        String(item["name"]) => Int.(item["word"]) for item in data["words"])
    lc = LickorishCode(
        Int(data["genus"]), words,
        _hc_string_dict(data["places"]),
        _hc_disk_dict(data["lc_disks"]),
        _hc_string_dict(data["lc_points"]),
        _hc_int_dict(data["lc_weights"]),
        nothing)
    roles = Dict{String,Symbol}(
        String(item["name"]) => Symbol(String(item["role"])) for item in data["roles"])
    hc = HeegaardCode(lc, roles,
        _hc_disk_dict(data["disks"]),
        _hc_string_dict(data["points"]),
        _hc_int_dict(data["weights"]),
        nothing)
    validate_hc_input(hc)
    hc
end

function _hc_form_data(form::HCCanonical)
    Dict{String,Any}(
        "genus" => form.genus,
        "k_words" => form.k_words,
        "alpha_words" => form.alpha_words,
        "beta_words" => form.beta_words,
        "other_words" => form.other_words,
        "places" => form.places,
        "lc_disks" => form.lc_disks,
        "hc_disks" => form.hc_disks,
        "lc_points" => form.lc_points,
        "hc_points" => form.hc_points,
        "lc_weights" => form.lc_weights,
        "hc_weights" => form.hc_weights,
    )
end

"""
    hc_search_state(hc::HeegaardCode; identity=:canonical) -> CanonicalState

Return the move-graph state record for an HC under a selected identity mode.

# Arguments
`identity` is `:canonical`, `:literal` or the corresponding string; other
values raise `ArgumentError`. Canonical mode requires the
[canonical comparison domain](../specs/heegaard-code.md#comparison-records)
and can raise `HCNotComputed`; writing and serialization errors propagate.

# Returns
`CanonicalState` with versioned identity data, the serialized HC writing
payload, three-dimensional area, `hc-search` representation and genus complexity.

# Guarantees
1. Canonical mode keys the canonical comparison record; literal mode keys
   the stored writing. [proved/agent: by construction]
2. The HC is not modified and no database write occurs.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
The key is not a manifold homeomorphism invariant. Canonical identity has
the limits of `hc_canonical`; serialization versions are part of identity.
The [search-record specification](../specs/heegaard-code.md#hc-search-records) states the payload and identity schema.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — scope of canonical identity.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> H.hc_search_state(hc).state_id == H.hc_search_state(hc).state_id
true
```

Role: convert
"""
function hc_search_state(hc::HeegaardCode; identity=:canonical)
    mode = _hc_identity(identity)
    writing = _hc_writing(hc)
    body = if mode == "literal"
        Dict{String,Any}("mode" => "literal", "writing" => writing)
    else
        Dict{String,Any}("mode" => "canonical", "form" => _hc_form_data(hc_canonical(hc)))
    end
    MoveGraphStore.CanonicalState(_hc_stable_json(body); area=MoveGraphStore.AREA_3,
        representation="hc-search", payload=_hc_stable_json(writing),
        canonical_version=HC_SEARCH_SCHEMA_VERSION,
        complexity=Float64(hc.lc.genus))
end

function _store_hc_state!(store, state, hc, mode)
    MoveGraphStore.put_state!(store, state)
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO hc_search_states(
          state_id, state_version, identity_mode, raw_json)
        VALUES (?,?,?,?)
    """, (state.state_id, HC_SEARCH_STATE_VERSION, mode, _hc_stable_json(_hc_writing(hc))))
    state
end

"""
    put_hc_search_state!(store::MoveGraphStoreDB, hc::HeegaardCode; identity=:canonical) -> CanonicalState

Store an HC writing under its selected search identity.

# Arguments
`store` must be open, writable and migrated by `migrate_hc_search_store!`;
an absent or unsupported HC-search schema raises `ArgumentError`.
`identity` is `:canonical`, `:literal`, or the corresponding string; other
symbol or string values raise `ArgumentError` naming the accepted modes.
`hc` must satisfy the input obligations of `hc_search_state`; canonicalization,
serialization and database errors propagate.

# Returns
The `CanonicalState` computed for this input, including its state ID,
serialized identity, writing payload, schema version and genus complexity.

# Guarantees
1. Successful return stores the state and an HC-writing record under its ID.
   Reusing an existing ID retains its previously stored HC writing.
   [proved/agent: by construction]
2. Canonical identity uses [HC comparison records](../specs/heegaard-code.md#comparison-records);
   literal identity uses the stored writing representation. The identity is
   a key for the chosen HC comparison, not a classification of 3-manifolds.
   [proved/agent: by construction]
3. The database is modified; `hc` is not modified.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
The returned writing need not be the representative retained for an already
present ID. This call does not validate a move or certify manifold equivalence,
and a failed database operation need not roll back every earlier write.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §比較用記録 — HC identities and correspondence scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> h = HeegaardCode(lc, Dict("l1" => :beta, "m1" => :alpha), Dict(), Dict(), Dict(), nothing);

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store); s = put_hc_search_state!(store, h);

julia> hc_search_state(load_hc_search_state(store, s.state_id)).state_id == s.state_id
true

julia> Topology.close_store!(store);
```

See also [`hc_search_state`](@ref), [`load_hc_search_state`](@ref).

Role: other
"""
function put_hc_search_state!(store::MoveGraphStore.MoveGraphStoreDB, hc::HeegaardCode; identity=:canonical)
    _require_hc_search_store(store)
    mode = _hc_identity(identity)
    _store_hc_state!(store, hc_search_state(hc; identity=mode), hc, mode)
end

"""
    load_hc_search_state(store::MoveGraphStoreDB, state_id::AbstractString) -> Union{Nothing,HeegaardCode}

Reconstruct an HC from the stored writing of a search state.

# Arguments
`store` must be open and migrated; absent or unsupported HC-search schemas
raise `ArgumentError`. An unknown `state_id` is accepted. Stored JSON must
have the supported writing schema and valid HC input: schema errors raise
`ArgumentError`, HC validation errors propagate, and malformed fields can
raise JSON, lookup or conversion errors.

# Returns
`nothing` when absent; otherwise a newly reconstructed `HeegaardCode` with
the saved LC words, roles and annotations, and no source file paths.

# Guarantees
1. A present state is reconstructed from its recorded writing and passes
   `validate_hc_input` before return. [proved/agent: by construction]
2. The store is not modified. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
This does not reauthenticate a state ID, recompute its search identity or
certify a manifold. If several writings share an identity, the stored writing
is returned, not a requested representative or a canonical drawing.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> h = HeegaardCode(lc, Dict("l1" => :beta, "m1" => :alpha), Dict(), Dict(), Dict(), nothing);

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store); s = put_hc_search_state!(store, h; identity=:canonical);

julia> hc_search_state(load_hc_search_state(store, s.state_id); identity=:canonical).state_id == s.state_id
true

julia> Topology.close_store!(store);
```

See also [`put_hc_search_state!`](@ref), [`hc_search_state`](@ref).

Role: access
"""
function load_hc_search_state(store::MoveGraphStore.MoveGraphStoreDB, state_id::AbstractString)
    _require_hc_search_store(store)
    row = MoveGraphStore._one(store.db, """
        SELECT raw_json FROM hc_search_states WHERE state_id=?
    """, (String(state_id),))
    row === nothing && return nothing
    _hc_from_writing(JSON.parse(String(row.raw_json)))
end

function _hc_script_snapshot(script_path::AbstractString)
    path = abspath(String(script_path))
    isfile(path) || throw(ArgumentError("provenance script does not exist: $path"))
    path, bytes2hex(sha256(read(path)))
end

"""
    begin_hc_search_run!(store::MoveGraphStoreDB, script_path::AbstractString; label::AbstractString, scope=Dict{String,Any}()) -> String

Create a running HC-search record with a snapshot of its provenance script's identity.

# Arguments
`store` must be open and migrated by `migrate_hc_search_store!`; an absent or
unsupported HC schema throws `ArgumentError` explaining the schema problem.
`script_path` must name an existing readable file; absence throws
`ArgumentError` naming the path, and read errors propagate. `label` is a
required string. `scope` must be representable by the search JSON format;
unsupported values or integers outside `Int64` throw `ArgumentError`.
Database errors propagate.

# Returns
The new run's identifier as a `String`, for `hc_search_run`, expansion and
`finish_hc_search_run!`.

# Guarantees
1. On success, the database contains a run with status `"running"`, the
   supplied label and scope, and the script's absolute path and content
   digest at this call. [proved/agent: by construction]
2. The store is modified; the script and scope are not modified and the
   script is not executed. [proved/agent: by construction]
3. The record identifies an experiment in the
   [HC search interface](../specs/module-structure.md#31-moveapi);
   it certifies no move, equivalence or manifold invariant.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
No search is started, no script bytes are archived, and no run is finished
automatically. The identifier is not a reproducible mathematical state key.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store);

julia> script = joinpath(pkgdir(Topology), "test", "heegaard", "hc_search.jl");

julia> id = begin_hc_search_run!(store, script; label="round-trip", scope=Dict{String,Any}("identity" => "canonical"));

julia> r = hc_search_run(store, id); (r.label, r.scope["identity"])
("round-trip", "canonical")

julia> Topology.close_store!(store);
```

See also [`hc_search_run`](@ref), [`finish_hc_search_run!`](@ref).

Role: other
"""
function begin_hc_search_run!(store::MoveGraphStore.MoveGraphStoreDB, script_path::AbstractString;
                              label::AbstractString, scope=Dict{String,Any}())
    _require_hc_search_store(store)
    path, digest = _hc_script_snapshot(script_path)
    started = string(Dates.now(Dates.UTC))
    run_id = "hc-search-run-" * bytes2hex(sha256(join(
        (replace(path, '\\' => '/'), digest, String(label), started), '\0')))
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO hc_search_runs(
          run_id, script_path, script_sha256, label, scope_json, status, counts_json,
          started_at_utc)
        VALUES (?,?,?,?,?,'running','{}',?)
    """, (run_id, replace(path, '\\' => '/'), digest, String(label),
          _hc_stable_json(scope), started))
    run_id
end

"""
    hc_search_run(store::MoveGraphStoreDB, run_id::AbstractString) -> Union{Nothing,NamedTuple}

Return the recorded status and operation counts of an HC search run.

# Arguments
`store` must be open and migrated by `migrate_hc_search_store!`.
Schema checks raise `ArgumentError` for an absent or unsupported HC-search
schema; database errors propagate. `run_id` is an opaque stored identifier;
an unknown identifier is accepted.

# Returns
`nothing` if absent; otherwise a named tuple containing `run_id`,
`script_path`, `script_sha256`, `label`, decoded `scope` and `counts`,
`status`, and the integer counts `operations`, `successes`, `rejects`, `errors`.

# Guarantees
1. The returned status and supplied counts describe the stored run; operation
   counts describe its currently stored operation rows. [proved/agent: by construction]
2. Neither the store nor its recorded HC states are modified.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
A recorded success is not a mathematical proof that its move preserves the
represented manifold. The distinction is owned by
[move witnesses](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store);

julia> script = joinpath(dirname(pathof(Topology)), "..", "test", "heegaard", "hc_search.jl");

julia> id = begin_hc_search_run!(store, script; label="round-trip");

julia> hc_search_run(store, id).label
"round-trip"

julia> Topology.close_store!(store);
```

See also [`begin_hc_search_run!`](@ref), [`finish_hc_search_run!`](@ref).

Role: access
"""
function hc_search_run(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString)
    _require_hc_search_store(store)
    row = MoveGraphStore._one(store.db, """
        SELECT run.*,
          (SELECT COUNT(*) FROM hc_search_operations WHERE run_id=run.run_id) AS operations,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='success') AS successes,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='rejected') AS rejects,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='error') AS errors
        FROM hc_search_runs AS run WHERE run.run_id=?
    """, (String(run_id),))
    row === nothing && return nothing
    (run_id=String(row.run_id), script_path=String(row.script_path),
     script_sha256=String(row.script_sha256), label=String(row.label),
     scope=JSON.parse(String(row.scope_json)), status=String(row.status),
     counts=JSON.parse(String(row.counts_json)),
     operations=Int(row.operations), successes=Int(row.successes),
     rejects=Int(row.rejects), errors=Int(row.errors))
end

"""
    finish_hc_search_run!(store::MoveGraphStoreDB, run_id::AbstractString; status::AbstractString, counts=Dict{String,Any}()) -> Union{Nothing,NamedTuple}

Finish a running HC-search record and return its stored summary.

# Arguments
`store` must be an open migrated HC-search store. `status` is `"completed"`,
`"budget-exhausted"` or `"failed"`; other values throw `ArgumentError`
identifying an invalid terminal status. `counts` must fit the search JSON
format, with unsupported values raising `ArgumentError`. Database and schema
errors propagate. An absent `run_id` is accepted and yields `nothing`.

# Returns
`nothing` for an absent run, otherwise the named-tuple summary from
`hc_search_run`: run ID, script path and digest, label, scope, status, counts,
and recorded operation/success/rejection/error counts.

# Guarantees
1. A currently running record is updated with the requested terminal status,
   counts and finish time. An already terminal record is left unchanged.
   [proved/agent: by construction]
2. The store can be modified; the supplied counts are not modified.
   [proved/agent: by construction]
3. This records the caller's experiment status in the
   [search interface](../specs/module-structure.md#31-moveapi), without
   validating a mathematical search result. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
The supplied counts are not checked against operations. Marking a run
completed does not prove exhaustive search or any manifold equivalence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store);

julia> id = begin_hc_search_run!(store, joinpath(pkgdir(Topology), "test", "heegaard", "hc_search.jl"); label="round-trip");

julia> finish_hc_search_run!(store, id; status="completed").status
"completed"

julia> Topology.close_store!(store);
```

See also [`begin_hc_search_run!`](@ref), [`hc_search_run`](@ref).

Role: other
"""
function finish_hc_search_run!(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString;
                               status::AbstractString, counts=Dict{String,Any}())
    status in _HC_SEARCH_TERMINAL || throw(ArgumentError("invalid terminal HC search status"))
    MoveGraphStore._dbexec(store.db, """
        UPDATE hc_search_runs
        SET status=?, counts_json=?, finished_at_utc=?
        WHERE run_id=? AND status='running'
    """, (String(status), _hc_stable_json(counts), string(Dates.now(Dates.UTC)), String(run_id)))
    hc_search_run(store, run_id)
end

function _hc_jsonable(value)
    value === nothing && return nothing
    value isa Bool && return value
    if value isa Integer
        return typemin(Int) <= value <= typemax(Int) ? Int(value) : string(value)
    end
    value isa AbstractString && return String(value)
    value isa Symbol && return String(value)
    value isa Rational && return string(value)
    value isa Pair && return [_hc_jsonable(first(value)), _hc_jsonable(last(value))]
    value isa NamedTuple && return Dict{String,Any}(String(key) => _hc_jsonable(getfield(value, key)) for key in keys(value))
    value isa AbstractDict && return Dict{String,Any}(string(key) => _hc_jsonable(item) for (key, item) in value)
    value isa Union{AbstractVector,Tuple} && return [_hc_jsonable(item) for item in value]
    string(value)
end

function _hc_support_brief(value)
    value isa AbstractDict || return _hc_jsonable(value)
    brief = Dict{String,Any}()
    for key in ("kind", "core", "curve", "side")
        haskey(value, key) || continue
        item = value[key]
        brief[key] = item isa AbstractDict ? _hc_support_brief(item) : _hc_jsonable(item)
    end
    brief
end

function _hc_step_brief(step::LCMoveWitness)
    data = Dict{String,Any}(
        "row_id" => step.row_id,
        "locations" => [collect(Int, location) for location in step.locations],
        "before" => step.before_frags,
        "after" => step.after_frags,
    )
    for key in ("q", "curve", "side", "acted", "acted_components", "disk")
        haskey(step.extra, key) || continue
        data[key] = _hc_jsonable(step.extra[key])
    end
    haskey(step.extra, "support") && (data["support"] = _hc_support_brief(step.extra["support"]))
    data
end

function _hc_candidate_data(witness::LCMoveWitness)
    data = Dict{String,Any}(
        "row_id" => witness.row_id,
        "names" => witness.names,
        "locations" => [collect(Int, location) for location in witness.locations],
        "vertices" => witness.vertices,
        "before" => witness.before_frags,
        "after" => witness.after_frags,
        "genus_delta" => witness.genus_delta,
    )
    for key in ("q", "support", "acted", "acted_components", "disk", "curve", "side",
                "k_identification", "name", "target", "curves", "band")
        haskey(witness.extra, key) || continue
        value = witness.extra[key]
        data[key] = key == "support" ? _hc_support_brief(value) : _hc_jsonable(value)
    end
    # Prepared twists and stabilizations keep the choice inside the sequence, not on the outer witness.
    if haskey(witness.extra, "preparation")
        data["preparation"] = _hc_jsonable(witness.extra["preparation"])
    end
    if haskey(witness.extra, "sequence")
        steps = witness.extra["sequence"]
        data["sequence_length"] = length(steps)
        data["sequence"] = [_hc_step_brief(step) for step in steps if step isa LCMoveWitness]
    end
    if haskey(witness.extra, "created_components")
        data["created_components"] = _hc_jsonable(witness.extra["created_components"])
    end
    attachment = get(witness.extra, "attachment", nothing)
    if attachment isa AbstractDict
        brief = Dict{String,Any}()
        for key in ("requested_disk", "actual_face", "reassigned")
            haskey(attachment, key) && (brief[key] = _hc_jsonable(attachment[key]))
        end
        data["attachment"] = brief
    end
    data
end

_hc_candidate_json(witness::LCMoveWitness) = _hc_stable_json(_hc_candidate_data(witness))

"""
    hc_search_candidates(hc::HeegaardCode) -> Vector{LCMoveWitness}

Return distinct recorded move candidates for an HC search state.

# Arguments
`hc` must be accepted by the registered move enumeration in
[the HC move domain](../specs/heegaard-moves.md).
Validation, candidate-enumeration and serialization errors propagate.

# Returns
A vector of `LCMoveWitness` candidates, unique under their stable serialized
selection data and retaining first-occurrence order.

# Guarantees
1. Each returned witness comes from the input's registered move enumeration;
   candidates with equal serialized selections appear only once.
   [proved/agent: by construction]
2. The input is not modified and no candidate is applied.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
Candidate enumeration is not an exhaustive list of geometric Heegaard moves.
A candidate need not apply successfully or yield a new state or a distinct
manifold; execution and inverse checks are separate.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §6 — candidate enumeration and move witnesses.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> length(H.hc_search_candidates(hc))
9
```

Role: access
"""
function hc_search_candidates(hc::HeegaardCode)
    seen = Set{String}()
    result = LCMoveWitness[]
    for witness in movable(hc)
        key = _hc_candidate_json(witness)
        key in seen && continue
        push!(seen, key)
        push!(result, witness)
    end
    result
end

function _hc_err_json(err)
    message = sprint(showerror, err)
    length(message) > 800 && (message = first(message, 800))
    _hc_stable_json(Dict{String,Any}("type" => string(typeof(err)), "message" => message))
end

function _hc_known_state(store, state_id)
    MoveGraphStore._one(store.db,
        "SELECT 1 AS present FROM hc_search_states WHERE state_id=?", (state_id,)) !== nothing
end

function _hc_bundle_text(hc::HeegaardCode)
    try
        lc_text(hc.lc) * "\n--- hc ---\n" * hc_text(hc; lc_filename="diagram.lc", encoding=:explicit)
    catch err
        "text-failed: " * sprint(showerror, err)
    end
end

function _hc_snapshot(hc::HeegaardCode)
    homology = try
        value = hc_first_homology(hc)
        Dict{String,Any}(
            "ok" => true,
            "free_rank" => value.free_rank,
            "torsion" => copy(value.torsion),
            "invariant_factors" => copy(value.invariant_factors),
        )
    catch err
        Dict{String,Any}("ok" => false, "error" => sprint(showerror, err))
    end
    linking = try
        form = hc_linking_form(hc)
        rows, cols = size(form.matrix)
        entries = String[]
        for col in 1:cols, row in 1:rows
            push!(entries, string(form.matrix[row, col]))
        end
        Dict{String,Any}(
            "ok" => true,
            "free_rank" => form.free_rank,
            "torsion" => [Int(order) for order in form.torsion],
            "matrix" => entries,
        )
    catch err
        Dict{String,Any}("ok" => false, "error" => sprint(showerror, err))
    end
    (homology=homology, linking=linking)
end

"""
    hc_search_same_homology(left::AbstractDict, right::AbstractDict) -> Bool

Compare two successful HC-search homology summaries by their abelian-group data.

# Arguments
`left` and `right` are string-keyed homology dictionaries from
the HC search snapshot. Missing or false `"ok"` returns `false`.
When both are successful, the caller supplies `"free_rank"` and `"torsion"`
in the normalized form of `hc_first_homology`; absent keys raise `KeyError`,
and non-Boolean `"ok"` values can raise a Boolean-context `TypeError`.
Hand-written, unnormalized torsion lists have no group-isomorphism guarantee.

# Returns
`Bool`, requiring both successful flags and equal free rank and torsion data.

# Guarantees
1. On normalized successful summaries, `true` means the recorded finitely
   generated abelian groups are isomorphic; Smith factors equal to one do
   not affect the result. [proved/agent: by construction]
2. Neither dictionary is modified. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
Failure to compute either summary returns `false` without establishing that
the groups differ. Equality does not identify the represented 3-manifolds.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — first homology as the abelianization of the presentation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> h = HeegaardCode(lc, Dict("l1" => :beta, "m1" => :alpha), Dict(), Dict(), Dict(), nothing);

julia> hom = hc_first_homology(h); d = Dict("ok" => true, "free_rank" => hom.free_rank, "torsion" => hom.torsion);

julia> hc_search_same_homology(d, d)
true
```

See also [`hc_first_homology`](@ref), [`hc_search_explore`](@ref).

Role: condition
"""
function hc_search_same_homology(left::AbstractDict, right::AbstractDict)
    # `invariant_factors` in this code keeps the Smith 1s, so it changes with the presentation.
    # The group is the free rank together with the torsion coefficients above 1.
    get(left, "ok", false) && get(right, "ok", false) || return false
    left["free_rank"] == right["free_rank"] && left["torsion"] == right["torsion"]
end

function _hc_rational_matrix(entries)
    n = isqrt(length(entries))
    n * n == length(entries) || throw(ArgumentError("linking matrix is not square"))
    values = [parse(Rational{BigInt}, String(entry)) for entry in entries]
    reshape(values, n, n)
end

"""
    hc_search_same_linking(left::AbstractDict, right::AbstractDict) -> Bool

Compare two stored linking-form snapshots under the search's supported coordinate tests.

# Arguments
Each successful snapshot supplies `ok`, `free_rank`, `torsion` and serialized
rational `matrix` entries from the HC-search invariant snapshot format.
Missing successful-record fields, invalid rational strings and nonsquare
matrix arrays raise lookup, parse or `ArgumentError` exceptions.
The caller supplies valid positive torsion orders and consistent dimensions.
An absent or false `ok` returns false.

# Returns
A Boolean indicating that the available snapshots pass the supported
linking-form comparison.

# Guarantees
1. A true result has matching free rank and torsion orders and passes the
   implemented matrix/isometry comparison. [proved/agent: by construction]
2. Neither dictionary is modified. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
This is not a general finite-linking-form isometry algorithm: multiple
factors and cyclic orders above 10000 use literal matrix equality. False can
mean unavailable data or unsupported coordinate equivalence, not different
manifolds. The snapshot schema is specified in [HC search records](../specs/heegaard-code.md#hc-search-records).

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md — linking-form coordinate and search interpretation boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> s = Dict("ok"=>true,"free_rank"=>0,"torsion"=>Int[]); H.hc_search_same_linking(s,s)
true
```

Role: condition
"""
function hc_search_same_linking(left::AbstractDict, right::AbstractDict)
    get(left, "ok", false) && get(right, "ok", false) || return false
    left["free_rank"] == right["free_rank"] || return false
    torsion = Int.(left["torsion"])
    torsion == Int.(right["torsion"]) || return false
    isempty(torsion) && return true
    left_matrix = _hc_rational_matrix(left["matrix"])
    right_matrix = _hc_rational_matrix(right["matrix"])
    if length(torsion) == 1
        order = torsion[1]
        order > 10_000 && return left_matrix == right_matrix
        return any(unit -> gcd(unit, order) == 1 &&
            isinteger(left_matrix[1, 1] - right_matrix[1, 1] * unit * unit), 1:(order - 1))
    end
    left_matrix == right_matrix
end

function _hc_note!(problems, kind, row_id, depth, detail, on_problem)
    key = (kind, row_id)
    if !haskey(problems, key) || depth < problems[key]["depth"]
        record = Dict{String,Any}(
            "kind" => kind,
            "row_id" => row_id,
            "depth" => depth,
        )
        for (name, value) in detail
            record[String(name)] = value
        end
        problems[key] = record
        on_problem === nothing || on_problem(record)
    end
end

function _hc_invariant_notes!(problems, start_snap, snap, row_id, depth, detail, on_problem)
    if !snap.homology["ok"]
        _hc_note!(problems, "homology_error", row_id, depth,
            merge(detail, Dict("value" => snap.homology)), on_problem)
    elseif start_snap.homology["ok"] && !hc_search_same_homology(start_snap.homology, snap.homology)
        _hc_note!(problems, "homology", row_id, depth,
            merge(detail, Dict("value" => snap.homology, "start" => start_snap.homology)), on_problem)
    end
    if !snap.linking["ok"]
        _hc_note!(problems, "linking_error", row_id, depth,
            merge(detail, Dict("value" => snap.linking)), on_problem)
    elseif start_snap.linking["ok"] && !hc_search_same_linking(start_snap.linking, snap.linking)
        _hc_note!(problems, "linking", row_id, depth,
            merge(detail, Dict("value" => snap.linking, "start" => start_snap.linking)), on_problem)
    end
end

"""
    expand_hc_search_state!(store::MoveGraphStoreDB, run_id::AbstractString, source::HeegaardCode; depth=0, identity=:canonical, room=typemax(Int), check_inverse=true, check_invariants=false, start_snap=nothing, problems=nothing, on_problem=nothing) -> NamedTuple

Store attempted outgoing HC-move transitions from one search state.

# Arguments
`store` must be open, writable and migrated for HC search; the schema is checked.
`run_id` must be usable by the store's run records. `identity` is `:canonical`
or `:literal`, checked with `ArgumentError`. `depth` and `room` convert to `Int`;
callers supply meaningful nonnegative search limits. `source` must be accepted
by the [HC move domain](../specs/heegaard-moves.md).
Optional invariant checking needs `start_snap` and a mutable `problems` collection;
`on_problem`, when supplied, must accept the recorded problem.

# Returns
A named tuple with `records` (operation, move-row, outcome, inverse status,
target ID, target diagram and target depth), `created` (new states), `stopped`
(capacity stopped enumeration), `source_state_id` and `candidate_count`.

# Guarantees
1. Returned records describe attempted transitions from the stored source;
   an `ok` inverse status means the chosen state identity returned to the
   source identity. [proved/agent: by construction]
2. At most `room` new target states are admitted when `room` is nonnegative;
   existing target states do not consume capacity. [proved/agent: by construction]
3. The store and supplied problem collection can be modified. No mutation
   of `source` is requested by this entry. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
State-key agreement is not a proof of manifold equivalence beyond the move's
own witness. Failed operations may already have written states. Callbacks,
SQL writes and candidate construction can throw; there is no all-or-nothing
transaction. At the capacity boundary a candidate can have been computed
without its operation being recorded.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §6 — move candidates, witnesses and partial domains.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> store = Topology.open_store(":memory:"); H.migrate_hc_search_store!(store); hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> try H.expand_hc_search_state!(store, "unused", hc; identity=:unknown) catch e; e isa ArgumentError end
true

julia> Topology.close_store!(store);
```

Role: other
"""
function expand_hc_search_state!(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString,
                                 source::HeegaardCode;
                                 depth::Integer=0, identity=:canonical, room::Integer=typemax(Int),
                                 check_inverse::Bool=true, check_invariants::Bool=false,
                                 start_snap=nothing, problems=nothing, on_problem=nothing)
    _require_hc_search_store(store)
    mode = _hc_identity(identity)
    source_depth = Int(depth)
    source_state = put_hc_search_state!(store, source; identity=mode)
    records = NamedTuple[]
    created = 0
    stopped = false
    candidates = hc_search_candidates(source)
    for witness in candidates
        candidate = _hc_candidate_data(witness)
        candidate_json = _hc_stable_json(candidate)
        outcome = "error"
        inverse_status = "unchecked"
        target = nothing
        target_state_id = nothing
        error_json = nothing
        local result
        try
            result = move(source, witness)
            target = result.output
            target isa HeegaardCode || throw(LCMoveError(:expansion, "move did not return an HC",
                Dict("row" => witness.row_id, "type" => string(typeof(target)))))
            target_state = hc_search_state(target; identity=mode)
            known = _hc_known_state(store, target_state.state_id)
            if !known && room - created <= 0
                stopped = true
                target = nothing
                break
            end
            _store_hc_state!(store, target_state, target, mode)
            known || (created += 1)
            target_state_id = target_state.state_id
            outcome = "success"
            if check_inverse
                try
                    back = move(target, result.inverse)
                    back.output isa HeegaardCode || error("inverse did not return an HC")
                    back_state = hc_search_state(back.output; identity=mode)
                    inverse_status = back_state.state_id == source_state.state_id ? "ok" : "mismatch"
                    if inverse_status == "mismatch" && problems !== nothing
                        _hc_note!(problems, "inverse_mismatch", witness.row_id, source_depth + 1,
                            Dict("source_text" => _hc_bundle_text(source),
                                 "target_text" => _hc_bundle_text(target),
                                 "back_text" => _hc_bundle_text(back.output),
                                 "witness" => candidate,
                                 "source_state_id" => source_state.state_id,
                                 "back_state_id" => back_state.state_id), on_problem)
                    end
                catch err
                    inverse_status = "error"
                    error_json = _hc_err_json(err)
                    problems === nothing || _hc_note!(problems,
                        err isa HCNotComputed ? "key" : "inverse_error",
                        witness.row_id, source_depth + 1,
                        Dict("source_text" => _hc_bundle_text(source),
                             "target_text" => _hc_bundle_text(target),
                             "witness" => candidate,
                             "value" => JSON.parse(error_json)), on_problem)
                end
            end
            if check_invariants && start_snap !== nothing && problems !== nothing && !known
                snap = _hc_snapshot(target)
                _hc_invariant_notes!(problems, start_snap, snap, witness.row_id, source_depth + 1,
                    Dict("source_text" => _hc_bundle_text(source),
                         "target_text" => _hc_bundle_text(target),
                         "witness" => candidate), on_problem)
            end
        catch err
            err isa InterruptException && rethrow()
            target = nothing
            target_state_id = nothing
            outcome = err isa LCMoveError ? "rejected" : "error"
            inverse_status = "unchecked"
            error_json = _hc_err_json(err)
            kind = if err isa HCNotComputed
                "key"
            elseif err isa LCMoveError
                "move_rejected"
            else
                "move_error"
            end
            problems === nothing || _hc_note!(problems, kind, witness.row_id, source_depth + 1,
                Dict("source_text" => _hc_bundle_text(source),
                     "witness" => candidate,
                     "value" => JSON.parse(error_json)), on_problem)
        end
        operation_id = "hc-op-" * bytes2hex(sha256(join(
            (String(run_id), source_state.state_id, HC_SEARCH_ENGINE_VERSION, candidate_json), '\0')))
        MoveGraphStore._dbexec(store.db, """
            INSERT OR IGNORE INTO hc_search_operations(
              operation_id, run_id, source_state_id, target_state_id, depth, row_id,
              engine_version, candidate_json, outcome, inverse_status, error_json)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, (operation_id, String(run_id), source_state.state_id,
              target_state_id === nothing ? missing : target_state_id,
              source_depth + 1, witness.row_id, HC_SEARCH_ENGINE_VERSION, candidate_json,
              outcome, inverse_status, error_json === nothing ? missing : error_json))
        push!(records, (
            operation_id=operation_id, row_id=witness.row_id, outcome=outcome,
            inverse_status=inverse_status, target_state_id=target_state_id,
            target=target, depth=source_depth + 1,
        ))
        stopped && break
    end
    (records=records, created=created, stopped=stopped, source_state_id=source_state.state_id,
     candidate_count=length(candidates))
end

"""
    hc_search_neighbors(store::MoveGraphStoreDB, state_id::AbstractString; row_ids=nothing) -> Vector

Return stored successful outgoing edges of an HC-search state.

# Arguments
`store` must be open and migrated by `migrate_hc_search_store!`; schema
checks raise `ArgumentError` and database errors propagate. `state_id` is
the stored source identifier. Optional `row_ids` is an iterable of move-row
names convertible to strings; conversion errors propagate. An empty filter
is accepted and selects no edges.

# Returns
A vector of named tuples with `operation_id`, `run_id`, `source_state_id`,
`target_state_id`, `depth`, `row_id`, `inverse_status`, and decoded `candidate`.

# Guarantees
1. The result consists of the stored successful edges for that source and
   filter, ordered by depth, row ID and operation ID.
   [proved/agent: by construction]
2. The query does not create neighbors or modify the store.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
The result is not the set of all mathematically possible moves, and successful
rows do not independently certify manifold equivalence. See
[move witnesses](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store);

julia> isempty(hc_search_neighbors(store, "absent"; row_ids=String[]))
true

julia> Topology.close_store!(store);
```

See also [`expand_hc_search_state!`](@ref), [`query_hc_search_paths`](@ref).

Role: access
"""
function hc_search_neighbors(store::MoveGraphStore.MoveGraphStoreDB, state_id::AbstractString; row_ids=nothing)
    _require_hc_search_store(store)
    if row_ids === nothing
        rows = MoveGraphStore._rows(store.db, """
            SELECT * FROM successful_hc_search_edges
            WHERE source_state_id=? ORDER BY depth, row_id, operation_id
        """, (String(state_id),))
    else
        allowed = unique(String.(collect(row_ids)))
        isempty(allowed) && return NamedTuple[]
        marks = join(fill("?", length(allowed)), ",")
        rows = MoveGraphStore._rows(store.db, """
            SELECT * FROM successful_hc_search_edges
            WHERE source_state_id=? AND row_id IN ($marks)
            ORDER BY depth, row_id, operation_id
        """, Tuple(Any[String(state_id); allowed]))
    end
    [(operation_id=String(row.operation_id), run_id=String(row.run_id),
      source_state_id=String(row.source_state_id), target_state_id=String(row.target_state_id),
      depth=Int(row.depth), row_id=String(row.row_id),
      inverse_status=String(row.inverse_status),
      candidate=JSON.parse(String(row.candidate_json))) for row in rows]
end

"""
    query_hc_search_paths(store::MoveGraphStoreDB, source_state_id::AbstractString; target_state_id, min_depth=1, max_depth=12, limit=100, row_ids=nothing) -> Vector{NamedTuple}

Return bounded simple paths through recorded successful HC-search transitions.

# Arguments
`store` must be open and migrated by `migrate_hc_search_store!`; schema checks
raise `ArgumentError`, and database errors propagate. Source and target are
stored state identifiers. Checked: `0 <= min_depth <= max_depth` and
`limit > 0`; `ArgumentError` identifies invalid depth bounds or a nonpositive
limit. `row_ids`, when supplied, is an iterable of move-row names convertible
to strings. The caller uses identifiers produced by the store.

# Returns
A vector of `(depth, steps)` named tuples. Each step has `operation_id`,
`source_state_id`, `target_state_id`, stored `depth`, `row_id`,
`inverse_status`, and decoded `candidate`. Path depth counts its edges;
a step's stored depth belongs to the recorded exploration.

# Guarantees
1. Returned paths connect the requested identifiers using successful edges
   allowed by the filter, with no repeated state, within the depth bounds.
   Results are ordered by path depth and their serialized operation-ID
   sequence, truncated to `limit`. [proved/agent: by construction]
2. An empty row filter returns no paths. The query does not modify the store.
   [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
A missing path is not evidence of inequivalent manifolds. A recorded path
inherits the [move witness limitations](../specs/heegaard-moves.md#comparison-and-move-witnesses).
This query does not generate moves or replay their witnesses.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> store = Topology.open_store(":memory:"); migrate_hc_search_store!(store);

julia> paths = query_hc_search_paths(store, "absent"; target_state_id="absent", min_depth=0); (length(paths), paths[1].depth, isempty(paths[1].steps))
(1, 0, true)

julia> Topology.close_store!(store);
```

See also [`hc_search_neighbors`](@ref).

Role: access
"""
function query_hc_search_paths(store::MoveGraphStore.MoveGraphStoreDB, source_state_id::AbstractString;
                               target_state_id::AbstractString, min_depth::Integer=1,
                               max_depth::Integer=12, limit::Integer=100, row_ids=nothing)
    _require_hc_search_store(store)
    min_depth >= 0 && max_depth >= min_depth || throw(ArgumentError("invalid path depth bounds"))
    limit > 0 || throw(ArgumentError("path query limit must be positive"))
    allowed = row_ids === nothing ? nothing : unique(String.(collect(row_ids)))
    allowed !== nothing && isempty(allowed) && return NamedTuple[]
    edge_filter = allowed === nothing ? "1=1" :
        "edge.row_id IN (" * join(fill("?", length(allowed)), ",") * ")"
    sql = """
        WITH RECURSIVE walk(state_id, depth, operation_ids, visited) AS (
          SELECT ?, 0, '', '|' || ? || '|'
          UNION ALL
          SELECT edge.target_state_id, walk.depth + 1,
                 CASE WHEN walk.operation_ids='' THEN edge.operation_id
                      ELSE walk.operation_ids || ',' || edge.operation_id END,
                 walk.visited || edge.target_state_id || '|'
          FROM walk
          JOIN successful_hc_search_edges AS edge ON edge.source_state_id=walk.state_id
          WHERE walk.depth < ? AND $edge_filter
            AND instr(walk.visited, '|' || edge.target_state_id || '|')=0
        )
        SELECT operation_ids, depth FROM walk
        WHERE state_id=? AND depth>=? AND depth<=?
        ORDER BY depth, operation_ids LIMIT ?
    """
    params = Any[String(source_state_id), String(source_state_id), Int(max_depth)]
    allowed === nothing || append!(params, allowed)
    append!(params, [String(target_state_id), Int(min_depth), Int(max_depth), Int(limit)])
    rows = MoveGraphStore._rows(store.db, sql, Tuple(params))
    paths = NamedTuple[]
    for row in rows
        ids = isempty(String(row.operation_ids)) ? String[] : split(String(row.operation_ids), ',')
        steps = NamedTuple[]
        for operation_id in ids
            step = MoveGraphStore._one(store.db, """
                SELECT * FROM successful_hc_search_edges WHERE operation_id=?
            """, (String(operation_id),))
            step === nothing && throw(ArgumentError("search path operation is absent"))
            push!(steps, (
                operation_id=String(step.operation_id),
                source_state_id=String(step.source_state_id),
                target_state_id=String(step.target_state_id),
                depth=Int(step.depth), row_id=String(step.row_id),
                inverse_status=String(step.inverse_status),
                candidate=JSON.parse(String(step.candidate_json)),
            ))
        end
        push!(paths, (depth=Int(row.depth), steps=steps))
    end
    paths
end

"""
    hc_search_explore(root::HeegaardCode, store::MoveGraphStoreDB, script_path::AbstractString; label::AbstractString, identity=:canonical, max_depth=3, max_states=2000, check_invariants=true, progress=false, on_problem=nothing) -> NamedTuple

Explore and record a bounded HC move graph from a root diagram.

# Arguments
`store` must be open and writable; the HC-search schema is migrated.
`script_path` must satisfy `begin_hc_search_run!` and `label` names the run.
`identity` is a supported canonical/literal mode, and nonnegative depth and
positive state bounds are checked with `ArgumentError`. The root must satisfy
[the HC move domain](../specs/heegaard-moves.md). Optional `on_problem` accepts
a recorded problem; database, script, callback and root-state errors propagate.

# Returns
A named tuple with run identity and record, `status`, counts, reached
`states`, `transitions`, `depth`, initial invariant snapshot and `problems`.

# Guarantees
1. The run records attempted moves at source depths below `max_depth` and
   uses one identity mode throughout. [proved/agent: by construction]
2. Invariant discrepancies are reported as problems, not silently converted
   into successful equivalence proofs. [proved/agent: by construction]
3. The store is modified, optional progress is printed and callbacks may run;
   the root is not modified by this entry. [proved/agent: by construction]

# Conventions
HC search-record format: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#hc-search-records).

# Not guaranteed
This is not a complete manifold-equivalence search. Completion means the
implemented bounded exploration completed, not that every geometric move was
tried. A prepopulated store and shared state identities affect capacity
accounting; there is no transaction rolling back an entire failed run.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md#comparison-and-move-witnesses — search records and equivalence evidence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> store = Topology.open_store(":memory:");

julia> try H.hc_search_explore(hc,store,"unused"; label="example",max_depth=-1) catch e; e isa ArgumentError end
true

julia> Topology.close_store!(store);
```

Role: other
"""
function hc_search_explore(root::HeegaardCode, store::MoveGraphStore.MoveGraphStoreDB,
                           script_path::AbstractString;
                           label::AbstractString, identity=:canonical,
                           max_depth::Integer=3, max_states::Integer=2000,
                           check_invariants::Bool=true, progress::Bool=false,
                           on_problem=nothing)
    mode = _hc_identity(identity)
    max_depth >= 0 && max_states >= 1 || throw(ArgumentError("invalid search bounds"))
    migrate_hc_search_store!(store)
    scope = Dict{String,Any}(
        "identity" => mode, "max_depth" => Int(max_depth), "max_states" => Int(max_states),
        "label" => String(label),
    )
    run_id = begin_hc_search_run!(store, script_path; label=String(label), scope=scope)
    problems = Dict{Tuple{String,String},Dict{String,Any}}()
    status = "completed"
    counts = Dict{String,Any}()
    try
        root_state = put_hc_search_state!(store, root; identity=mode)
        start_snap = check_invariants ? _hc_snapshot(root) : nothing
        if start_snap !== nothing
            _hc_invariant_notes!(problems, start_snap, start_snap, "(start)", 0,
                Dict("source_text" => _hc_bundle_text(root), "witness" => nothing), on_problem)
        end
        queue = [(root, 0)]
        depths = Dict(root_state.state_id => 0)
        created_total = 1
        transitions = 0
        while !isempty(queue) && status == "completed"
            hc, depth = popfirst!(queue)
            depth >= max_depth && continue
            if progress
                println("hc-search label=$label depth=$depth states=$created_total queue=$(length(queue)) transitions=$transitions")
                flush(stdout)
            end
            room = max_states - created_total
            expansion = expand_hc_search_state!(store, run_id, hc;
                depth=depth, identity=mode, room=room, check_inverse=true,
                check_invariants=check_invariants, start_snap=start_snap,
                problems=problems, on_problem=on_problem)
            transitions += length(expansion.records)
            for record in expansion.records
                record.outcome == "success" || continue
                record.target_state_id in keys(depths) && continue
                depths[record.target_state_id] = record.depth
                created_total += 1
                record.depth < max_depth && push!(queue, (record.target, record.depth))
            end
            expansion.stopped && (status = "budget-exhausted")
        end
        counts = Dict{String,Any}(
            "states" => created_total,
            "transitions" => transitions,
            "max_depth_reached" => isempty(depths) ? 0 : maximum(values(depths)),
            "identity" => mode,
            "problems" => length(problems),
        )
        finished = finish_hc_search_run!(store, run_id; status=status, counts=counts)
        return (
            run_id=run_id, status=status, identity=mode, counts=counts, run=finished,
            states=created_total, transitions=transitions,
            depth=counts["max_depth_reached"],
            start=start_snap, problems=collect(values(problems)),
        )
    catch err
        err isa InterruptException && rethrow()
        counts["error"] = sprint(showerror, err)
        try
            finish_hc_search_run!(store, run_id; status="failed", counts=counts)
        catch
        end
        rethrow()
    end
end
