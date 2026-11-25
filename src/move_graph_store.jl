module MoveGraphStore

using SQLite
using DBInterface
using JSON
using SHA
using Dates

include("domain_taxonomy.jl")
using .DomainTaxonomy

export CanonicalState, MoveType, MoveGraphStoreDB, Store,
       AREA_2, AREA_2_5, AREA_3,
       open_store, close_store!, migrate!, state_id, state_hash,
       put_state!, get_state, register_move_type!, get_move_type,
       insert_verified_transition!, record_transition!, add_certificate!,
       add_invariant!, neighbors, bfs_certificate_path, certificate_path,
       import_knotinfo!, knotinfo_records, state_from_serialization,
       ClassificationRun,
       register_classification_run!, put_state_classification!,
       record_state_classification!, get_classification_run,
       classification_run, get_state_classification, state_classification,
       classification_rows, classification_counts

"""
    CanonicalState(serialization::AbstractString; area, representation, payload=nothing, canonical_version=1, complexity=nothing, state_id=nothing, state_hash=nothing) -> CanonicalState

Package a caller-supplied canonical serialization as a state of a move graph.

# Arguments
- `serialization` is caller-canonicalized text for the chosen representation;
  no parser or mathematical canonicalizer is called.
- The keyword constructor checks positive integer `canonical_version`, an
  area accepted by [`normalize_area`](@ref), and any supplied identifier or
  hash against the supplied text. `ArgumentError` names a nonpositive version,
  invalid area, or mismatched identifier/hash. Field conversion errors propagate.
- `representation` is normalized with [`normalize_representation`](@ref).
  `complexity` is absent or convertible to `Float64`; `payload` is arbitrary.
- The eight-field constructor checks only identifier/hash consistency.
  Its caller supplies valid domain, version and representation metadata.

# Returns
`CanonicalState` with the serialization, its content hash, version-sensitive
identifier, canonical version, payload, area, representation and optional
complexity. The storage API belongs to
[MoveGraphStore](../specs/move-graph-store.md#definitions).

# Guarantees
1. Both constructors reject inconsistent identifiers and hashes; neither
   changes the serialization's text. [proved/agent: by construction]
2. The payload is retained without copying or modification. Other inputs
   are not modified. [proved/agent: by construction]
3. The result identifies the supplied textual state; it identifies a
   mathematical equivalence class only to the extent that the caller's
   canonicalization does so. [proved/agent: by construction]

# Not guaranteed
No diagram validity, canonicality, collision-free hashing, or equivalence of
represented knots, manifolds or other spaces is certified. A payload need
not agree with the serialization. The eight-field constructor does not
provide the keyword constructor's metadata validation.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, CanonicalState — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> s = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss", payload=[1, 2]);

julia> store = open_store(); saved = put_state!(store, s);

julia> saved.state_id == s.state_id
true

julia> close_store!(store);
```

See also [`state_from_serialization`](@ref), [`put_state!`](@ref).

Role: construct
"""
struct CanonicalState
    state_id::String
    state_hash::String
    canonical_version::Int
    serialization::String
    payload
    area::Float64
    representation::String
    complexity::Union{Nothing,Float64}

    function CanonicalState(state_id::String, state_hash::String,
                            canonical_version::Int, serialization::String,
                            payload, area::Float64, representation::String,
                            complexity::Union{Nothing,Float64})
        expected_hash = bytes2hex(SHA.sha256(codeunits(serialization)))
        state_hash == expected_hash || throw(ArgumentError(
            "state_hash does not match serialization"))
        expected_id = bytes2hex(SHA.sha256(
            codeunits("v$(canonical_version)\0$(serialization)")))
        state_id == expected_id || throw(ArgumentError(
            "state_id does not match canonical serialization"))
        new(state_id, state_hash, canonical_version, serialization, payload,
            area, representation, complexity)
    end
end

function CanonicalState(serialization::AbstractString; area, representation,
                        payload=nothing, canonical_version::Integer=1,
                        complexity=nothing, state_id=nothing, state_hash=nothing)
    v = Int(canonical_version)
    v > 0 || throw(ArgumentError("canonical_version must be positive"))
    text = String(serialization)
    expected_hash = bytes2hex(SHA.sha256(codeunits(text)))
    h = state_hash === nothing ? expected_hash : String(state_hash)
    h == expected_hash || throw(ArgumentError("state_hash does not match serialization"))
    expected_id = bytes2hex(SHA.sha256(codeunits("v$(v)\0$(text)")))
    sid = state_id === nothing ? expected_id : String(state_id)
    sid == expected_id || throw(ArgumentError("state_id does not match canonical serialization"))
    c = complexity === nothing ? nothing : Float64(complexity)
    CanonicalState(sid, h, v, text, payload, normalize_area(area), normalize_representation(representation), c)
end

"""
    state_from_serialization(serialization::AbstractString; area, representation, payload=nothing, canonical_version=1, complexity=nothing, state_id=nothing, state_hash=nothing) -> CanonicalState

Package a canonical serialization as a move graph state with checked identifiers.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

Arguments and checks are those of the keyword constructor of
[`CanonicalState`](@ref): the caller supplies canonical text and coherent
metadata; positive version, accepted domain and any supplied identifier/hash
are checked. `ArgumentError` identifies a version, domain or digest mismatch;
field conversion errors propagate.

# Returns
The corresponding `CanonicalState`, with normalized representation label,
content hash, version-sensitive identifier and supplied optional payload
and complexity.

# Guarantees
1. The result has the postconditions of `CanonicalState(serialization; kwargs...)`.
   [proved/agent: by construction]
2. Inputs are not modified; the payload is shared without copying.
   [proved/agent: by construction]

# Not guaranteed
Serialization is not parsed or canonicalized. No validity or equivalence of
represented mathematical objects is established by computing the identifiers.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, state_from_serialization — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> b = Topology.MoveGraphStore.state_from_serialization("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> a.state_id == b.state_id
true
```

Role: construct
"""
state_from_serialization(serialization; kwargs...) = CanonicalState(serialization; kwargs...)

"""
    state_id(s::CanonicalState) -> String

Return the version-sensitive identifier of a stored textual state.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is a `CanonicalState`; no additional checks are performed. Its
serialization and identifiers obey the constructor's checks, while the
caller supplies their mathematical interpretation.

# Returns
The `String` in `s.state_id`.

# Guarantees
1. The result is the state's recorded identifier, including its canonical
   version in the identifier scheme. [proved/agent: by construction]
2. The state and its payload are not modified. [proved/agent: by construction]

# Not guaranteed
An identifier is not a proof of mathematical equivalence, and hashing is
not guaranteed collision-free. Domain and payload are not part of the
identifier's consistency check.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, state_id — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); s = open_store();

julia> state_id(put_state!(s, a)) == state_id(a)
true

julia> close_store!(s);
```

See also [`state_hash`](@ref).

Role: access
"""
state_id(s::CanonicalState) = s.state_id

"""
    state_hash(s::CanonicalState) -> String

Return the content hash of a state's serialization.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is a `CanonicalState`; this accessor performs no new validation.
The constructor checks consistency of this hash with the serialization.

# Returns
The `String` in `s.state_hash`.

# Guarantees
1. The result is the recorded hash of the serialization text; it does not
   depend on the canonical version or payload. [proved/agent: by construction]
2. The state and payload are not modified. [proved/agent: by construction]

# Not guaranteed
Hash equality does not certify mathematical equivalence or even collision-free
textual identity. Use [`state_id`](@ref) when the canonical version matters.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, state_hash — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); s = open_store(); put_state!(s, a);

julia> state_hash(get_state(s, a.state_id)) == state_hash(a)
true

julia> close_store!(s);
```

Role: access
"""
state_hash(s::CanonicalState) = s.state_hash

"""
    MoveType(move_type, area, inverse_type, inverse_direction) -> MoveType

Store the declared domain and inverse metadata of a move name.

# Arguments
Arguments must have or convert to the field types `String`, `Float64`,
`String`, and `Symbol`, respectively. Conversion errors propagate.
This constructor performs no taxonomy or inverse-move validation; use
[`register_move_type!`](@ref) for checked registration in a store.

# Returns
`MoveType` with the supplied name, domain, inverse name and inverse direction,
as metadata for [MoveGraphStore](../specs/move-graph-store.md#definitions).

# Guarantees
1. Fields preserve the supplied values after conversion.
   [proved/agent: by construction]
2. Arguments are not modified and no database is changed.
   [proved/agent: by construction]

# Not guaranteed
The record does not establish that either name denotes a valid move or
that the two moves are inverses on any mathematical structure.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, MoveType — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> m = MoveType("MP", Topology.AREA_3, "MP", :inverse);

julia> store = open_store(); registered = register_move_type!(store, "MP"; area=Topology.AREA_3);

julia> m.area == registered.area == Topology.AREA_3
true

julia> close_store!(store);
```

Role: construct
"""
struct MoveType
    move_type::String
    area::Float64
    inverse_type::String
    inverse_direction::Symbol
end

"""
    MoveGraphStoreDB(database_or_path) -> MoveGraphStoreDB

Return a handle to a persistent or in-memory move graph store.

# Arguments
A `SQLite.DB` is wrapped as supplied, without opening, migrating or checking
it; its caller is responsible for a usable connection and schema.
An `AbstractString` path uses [`open_store`](@ref), including schema migration
and its errors. The path `":memory:"` selects a new in-memory database.

# Returns
`MoveGraphStoreDB` whose `db::SQLite.DB` field is the database connection.
`Store` is an alias of this type. See the
[MoveGraphStore module](../specs/move-graph-store.md#definitions).

# Guarantees
1. Wrapping a connection retains the same connection without modifying it.
   [proved/agent: by construction]
2. Successful path construction has the postconditions of `open_store(path)`.
   [proved/agent: by construction]

# Not guaranteed
The handle does not certify the mathematical correctness of stored states,
invariants, or supplied verification evidence. Wrapping a connection does
not initialize a schema.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, MoveGraphStoreDB — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> store = MoveGraphStoreDB(":memory:");

julia> s = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> put_state!(store, s).state_id == s.state_id
true

julia> close_store!(store);
```

See also [`close_store!`](@ref).

Role: construct
"""
struct MoveGraphStoreDB
    db::SQLite.DB
end
const Store = MoveGraphStoreDB

"""
    ClassificationRun(method_version::String, canonical_version::Int, max_n::Int, barrier::Int, census_fingerprint::String, classification_fingerprint::String, counts, conventions, created_at_utc::String) -> ClassificationRun

Package the provenance metadata of a classification dataset.

# Arguments
The fields have the meaning of M4 and M13 in the
[store specification](../specs/move-graph-store.md#numbered-conditions).
The raw constructor converts declared field types but does not validate those
conditions; conversion failures propagate. Use `register_classification_run!`
to validate and persist a run. `counts` and `conventions` are unrestricted values.

# Returns
A `ClassificationRun` with the named fields in
[Schema and classifications](../specs/move-graph-store.md#schema-and-classifications).

# Guarantees
1. Fields contain the supplied values after conversion to their declared types. [proved/agent: by construction]
2. Construction does not modify arguments or a database. Mutable field values are shared with their arguments. [proved/agent: by construction]

# Not guaranteed
A raw record does not certify valid bounds, a reproducible census, correct
classification or agreement of its counts and fingerprints with any dataset.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md Schema and classifications — the metadata record.

# Examples
```jldoctest
julia> s = open_store(); r = register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> r isa ClassificationRun
true

julia> close_store!(s);
```

Role: construct
"""
struct ClassificationRun
    method_version::String
    canonical_version::Int
    max_n::Int
    barrier::Int
    census_fingerprint::String
    classification_fingerprint::String
    counts
    conventions
    created_at_utc::String
end

# Builders that need a reproducible artifact pass an explicit release timestamp.
_classification_now_utc() = string(Dates.now(Dates.UTC))

function _dbexec(db, sql, params=())
    DBInterface.execute(db, sql, params)
end
function _rows(db, sql, params=())
    # SQLite.Row is a forward-only view into the live statement.  Collecting
    # Row objects first leaves them pointing at the reset statement (and thus
    # reads as `missing`).  Materialize each row before advancing the cursor.
    DBInterface.execute(db, sql, params) do query
        names = Tuple(query.names)
        [NamedTuple{names}(Tuple(row[i] for i in 1:length(names))) for row in query]
    end
end
_one(db, sql, params=()) = (r = _rows(db, sql, params); isempty(r) ? nothing : first(r))
_val(row, key) = row === nothing ? nothing : row[key]

function _sql_path()
    normpath(joinpath(@__DIR__, "..", "db", "knot_move_schema.sql"))
end

function _execute_schema!(db, sql::AbstractString)
    # Split ordinary statements while keeping the semicolon inside CREATE
    # TRIGGER ... BEGIN ... END; intact.
    statement = IOBuffer()
    in_trigger = false
    for line in split(String(sql), '\n')
        stripped = strip(line)
        startswith(uppercase(stripped), "CREATE TRIGGER") && (in_trigger = true)
        print(statement, line, '\n')
        done = in_trigger ? endswith(stripped, "END;") : endswith(stripped, ";")
        if done
            text = strip(String(take!(statement)))
            endswith(text, ";") && (text = text[1:prevind(text, lastindex(text))])
            isempty(text) || _dbexec(db, text)
            in_trigger = false
        end
    end
    isempty(strip(String(take!(statement)))) ||
        throw(ArgumentError("schema contains an unterminated SQL statement"))
    db
end

function _existing_schema_version(db)
    has_meta = _one(db, """
        SELECT 1 AS present FROM sqlite_master
        WHERE type='table' AND name='schema_meta'
    """)
    if has_meta === nothing
        user_tables = _rows(db, """
            SELECT name FROM sqlite_master
            WHERE type='table' AND name NOT LIKE 'sqlite_%'
        """)
        isempty(user_tables) || throw(ArgumentError(
            "database has tables but no schema_meta; refusing to guess its schema"))
        return nothing
    end
    row = _one(db, "SELECT value FROM schema_meta WHERE key='schema_version'")
    row === nothing && throw(ArgumentError("schema_meta has no schema_version"))
    try
        parse(Int, String(row.value))
    catch
        throw(ArgumentError("schema_meta.schema_version is not an integer"))
    end
end

function _validate_legacy_evidence(row, raw::AbstractString)
    evidence = try
        JSON.parse(String(raw))
    catch
        throw(ArgumentError(
            "cannot migrate v1: transition $(row.id) has malformed certificate JSON"))
    end
    evidence isa AbstractDict || throw(ArgumentError(
        "cannot migrate v1: transition $(row.id) certificate must be a JSON object"))
    required = ("adapter_version", "move_type", "source_serialization",
                "target_serialization", "theorem")
    all(haskey(evidence, key) for key in required) || throw(ArgumentError(
        "cannot migrate v1: transition $(row.id) certificate is missing required fields"))
    evidence["adapter_version"] == 1 || throw(ArgumentError(
        "cannot migrate v1: transition $(row.id) has an unsupported adapter_version"))
    String(evidence["move_type"]) == String(row.move_type) || throw(ArgumentError(
        "cannot migrate v1: transition $(row.id) move_type does not match its certificate"))
    String(evidence["source_serialization"]) == String(row.source_serialization) ||
        throw(ArgumentError(
            "cannot migrate v1: transition $(row.id) source does not match its certificate"))
    String(evidence["target_serialization"]) == String(row.target_serialization) ||
        throw(ArgumentError(
            "cannot migrate v1: transition $(row.id) target does not match its certificate"))
    !isempty(strip(String(evidence["theorem"]))) || throw(ArgumentError(
        "cannot migrate v1: transition $(row.id) certificate theorem is empty"))
    true
end

function _migrate_v1_to_v2!(db)
    columns = Set(String(row.name) for row in _rows(db, "PRAGMA table_info(transitions)"))
    "verification_evidence" in columns && throw(ArgumentError(
        "cannot migrate v1: transitions already has verification_evidence"))

    transition_count = Int(_one(db, "SELECT COUNT(*) AS n FROM transitions").n)
    certificate_count = Int(_one(db, "SELECT COUNT(*) AS n FROM certificate_paths").n)
    distinct_certified = Int(_one(db,
        "SELECT COUNT(DISTINCT transition_id) AS n FROM certificate_paths").n)
    orphan_count = Int(_one(db, """
        SELECT COUNT(*) AS n FROM certificate_paths c
        WHERE NOT EXISTS (SELECT 1 FROM transitions t WHERE t.id=c.transition_id)
    """).n)
    transition_count == certificate_count == distinct_certified ||
        throw(ArgumentError(
            "cannot migrate v1: every transition requires exactly one certificate"))
    orphan_count == 0 || throw(ArgumentError(
        "cannot migrate v1: orphan certificate paths are present"))

    rows = _rows(db, """
        SELECT t.id, t.move_type, t.parameters,
               source.serialization AS source_serialization,
               target.serialization AS target_serialization,
               certificate.path
        FROM transitions t
        JOIN states source ON source.state_id=t.source_state_id
        JOIN states target ON target.state_id=t.target_state_id
        JOIN certificate_paths certificate ON certificate.transition_id=t.id
        ORDER BY t.id
    """)
    length(rows) == transition_count || throw(ArgumentError(
        "cannot migrate v1: transition/certificate join is incomplete"))
    for row in rows
        path = String(row.path)
        isempty(strip(path)) && throw(ArgumentError(
            "cannot migrate v1: transition $(row.id) has an empty certificate"))
        _validate_legacy_evidence(row, path)
    end

    _dbexec(db, "ALTER TABLE transitions ADD COLUMN verification_evidence TEXT")
    _dbexec(db, """
        UPDATE transitions
        SET verification_evidence = (
            SELECT certificate.path FROM certificate_paths certificate
            WHERE certificate.transition_id=transitions.id
        )
    """)
    nothing
end

function _assert_v2_postconditions(db)
    version = _existing_schema_version(db)
    version == 2 || throw(ArgumentError("schema migration did not reach version 2"))
    columns = Set(String(row.name) for row in _rows(db, "PRAGMA table_info(transitions)"))
    "verification_evidence" in columns || throw(ArgumentError(
        "schema v2 is missing transitions.verification_evidence"))
    missing = Int(_one(db, """
        SELECT COUNT(*) AS n FROM transitions
        WHERE verification_evidence IS NULL
           OR length(trim(verification_evidence))=0
    """).n)
    missing == 0 || throw(ArgumentError(
        "schema v2 contains verified transitions without evidence"))
    isempty(_rows(db, "PRAGMA foreign_key_check")) || throw(ArgumentError(
        "schema v2 foreign-key check failed"))
    true
end

"""
    migrate!(store::MoveGraphStoreDB) -> MoveGraphStoreDB

Initialize or migrate the common move graph store to schema version two.

# Arguments
`store` must hold an open, writable SQLite connection outside a transaction.
Supported inputs are an empty database and schema versions one and two.
Legacy transitions must have unambiguous certificate data accepted by the
legacy evidence checks. `ArgumentError` identifies an unsupported or malformed
schema version, unversioned populated database, invalid legacy certificate,
missing evidence or failed foreign-key check. Database/I/O errors propagate.
See [MoveGraphStore](../specs/move-graph-store.md#definitions) for the module boundary.

# Returns
The same `MoveGraphStoreDB` handle.

# Guarantees
1. Success records schema version two, with a transition evidence column,
   nonempty evidence on existing transitions and no foreign-key violations.
   [proved/agent: by construction]
2. Compatible version-one certificates are retained and copied to transition
   evidence; existing move parameters are preserved. [proved/agent: by construction]
3. Schema/data changes in the migration transaction roll back on failure;
   the connection's foreign-key setting is enabled outside that transaction.
   [proved/agent: by construction]

# Not guaranteed
Migration does not replay mathematical certificates or repair unsupported
schemas. It is not a general schema-diff validator for arbitrary modified
databases. Existing caller transactions are not supported.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql schema version two — persistent common-store relations and constraints.

# Examples
```jldoctest
julia> s = open_store(); migrate!(s); migrate!(s);

julia> Topology.MoveGraphStore._one(s.db, "SELECT value FROM schema_meta WHERE key='schema_version'").value == "2"
true

julia> close_store!(s);
```

Role: other
"""
function migrate!(store::MoveGraphStoreDB)
    sql = read(_sql_path(), String)
    _dbexec(store.db, "PRAGMA foreign_keys = ON")
    _transaction(store.db) do
        version = _existing_schema_version(store.db)
        version === nothing || version <= 2 ||
            throw(ArgumentError("unsupported schema version $version"))
        if version === nothing
            _execute_schema!(store.db, sql)
        elseif version == 1
            _migrate_v1_to_v2!(store.db)
            _execute_schema!(store.db, sql)
            _dbexec(store.db,
                    "UPDATE schema_meta SET value='2' WHERE key='schema_version'")
        else
            _execute_schema!(store.db, sql)
        end
        _assert_v2_postconditions(store.db)
    end
    store
end

"""
    open_store(path::AbstractString=":memory:"; migrate=true, busy_timeout_ms::Integer=10000) -> MoveGraphStoreDB

Open an in-memory or on-disk move graph store.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`path` is a SQLite database path; the default selects a fresh in-memory
database. `busy_timeout_ms` is a nonnegative integer, checked with
`ArgumentError` describing a negative timeout. `migrate` is a Boolean:
when true, [`migrate!`](@ref) initializes/checks the supported schema.
File, SQLite and migration errors propagate. A writable path and compatible
database are the caller's responsibility.

# Returns
A `MoveGraphStoreDB` holding the open connection. The caller closes it
with [`close_store!`](@ref).

# Guarantees
1. The successful default call returns a store ready for the common-store
   API, with the postconditions of `migrate!`. [proved/agent: by construction]
2. Opening an on-disk path can create or migrate that database; opening
   `":memory:"` creates a separate in-memory store. [proved/agent: by construction]

# Not guaranteed
With `migrate=false`, schema and foreign-key enforcement are not initialized
by this function. Opening a store does not certify its mathematical contents.
A migration failure does not return a usable handle to the caller.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql — the common-store schema initialized by default.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> put_state!(s, a).state_id == a.state_id
true

julia> close_store!(s);
```

Role: construct
"""
function open_store(path::AbstractString=":memory:"; migrate=true,busy_timeout_ms::Integer=10_000)
    busy_timeout_ms>=0 || throw(ArgumentError("busy timeout must be nonnegative"))
    s = MoveGraphStoreDB(SQLite.DB(String(path)))
    SQLite.busy_timeout(s.db,busy_timeout_ms)
    migrate && migrate!(s)
    s
end
MoveGraphStoreDB(path::AbstractString) = open_store(path)
"""
    close_store!(s::MoveGraphStoreDB) -> Nothing

Close the connection of a move graph store.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is a handle returned by [`open_store`](@ref) or wrapping a SQLite
connection. SQLite connection/finalization errors propagate.

# Returns
`nothing`.

# Guarantees
1. On success the connection is closed; further operations requiring it
   cannot use it as an open store. [proved/agent: by construction]
2. Any other handles sharing `s.db` share the closed connection.
   [proved/agent: by construction]

# Not guaranteed
This does not delete an on-disk database, validate its mathematical contents,
or preserve an in-memory database for later reopening.

# Sources
cf:lib/topology.jl/src/move_graph_store.jl, close_store! — persistent record identity and connection semantics; no mathematical evidence is verified by construction.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> put_state!(s, a).state_id == a.state_id
true

julia> close_store!(s);
```

Role: other
"""
close_store!(s::MoveGraphStoreDB) = SQLite.close(s.db)

function _payload(x)
    x === nothing && return nothing
    x isa AbstractVector{UInt8} && return x
    x isa AbstractString && return String(x)
    JSON.json(x)
end

"""
    put_state!(s::MoveGraphStoreDB, state::CanonicalState) -> CanonicalState

Insert or update the payload of a textual state in the common move graph store.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` must be open and writable with the common-store schema. `state` supplies
serialization, identifiers and metadata. Its payload must be absent, text,
bytes, or JSON-serializable. Conversion, encoding and database errors propagate;
SQL constraints reject invalid version/domain metadata. No mathematical
canonicalization or validity check is performed.

# Returns
The stored `CanonicalState` read back with [`get_state`](@ref). Structured
payloads are returned as encoded JSON text, not as their original Julia type.

# Guarantees
1. A new state is stored with its supplied metadata. Re-inserting the same
   hash/version/serialization updates payload and complexity, retaining the
   existing area's and representation's values. [proved/agent: by construction]
2. A successful result has consistent identifier/hash data as checked by
   `get_state`. [proved/agent: by construction]
3. The database is modified; `state` and its payload are not modified.
   [proved/agent: by construction]

# Not guaranteed
An identical serialization is not checked to denote the same mathematical
object across domains or representations. The returned payload need not
have the input's Julia type. A failure while reading back is not a promise
to roll back an insertion performed outside a caller transaction.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql states — textual identities and persistence constraints.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss", payload=[1, 2]);

julia> put_state!(s, a).state_id == a.state_id
true

julia> close_store!(s);
```

Role: other
"""
function put_state!(s::MoveGraphStoreDB, state::CanonicalState)
    _dbexec(s.db, """
        INSERT INTO states(state_id,state_hash,canonical_version,serialization,payload,area,representation,complexity)
        VALUES (?,?,?,?,?,?,?,?)
        ON CONFLICT(state_hash,canonical_version,serialization) DO UPDATE SET
          payload=excluded.payload, complexity=excluded.complexity
    """, (state.state_id, state.state_hash, state.canonical_version, state.serialization,
           _payload(state.payload), state.area, state.representation, state.complexity))
    row = _one(s.db, "SELECT state_id FROM states WHERE state_hash=? AND canonical_version=? AND serialization=?",
               (state.state_hash, state.canonical_version, state.serialization))
    get_state(s, String(row.state_id))
end

"""
    get_state(s::MoveGraphStoreDB, sid::AbstractString) -> Union{Nothing,CanonicalState}

Read a stored textual state and check its identifier and content hash.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema; `sid` is the exact state identifier.
Missing identifiers are permitted. Invalid stored identifier/hash data raise
`ArgumentError` describing the mismatch. Database and field conversion
errors propagate.

# Returns
A `CanonicalState` containing the stored serialization, identifiers, version,
payload, area, representation and optional complexity, or `nothing` if absent.
Payload JSON text is returned as stored, without decoding it into the
original Julia object.

# Guarantees
1. A returned state's identifier and hash agree with its serialization and
   canonical version. [proved/agent: by construction]
2. The database is not modified. [proved/agent: by construction]

# Not guaranteed
No mathematical object is parsed or canonicalized, and no diagram validity
or agreement between payload and serialization is checked. Reconstruction
does not apply the keyword constructor's metadata normalization.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql states — persisted state representation and metadata.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a);

julia> get_state(s, a.state_id).state_id == a.state_id
true

julia> close_store!(s);
```

Role: access
"""
function get_state(s::MoveGraphStoreDB, sid::AbstractString)
    r = _one(s.db, "SELECT * FROM states WHERE state_id=?", (String(sid),))
    r === nothing && return nothing
    payload = r.payload === missing ? nothing : r.payload
    complexity = (r.complexity === nothing || r.complexity === missing) ? nothing : Float64(r.complexity)
    CanonicalState(String(r.state_id), String(r.state_hash), Int(r.canonical_version),
                   String(r.serialization), payload, Float64(r.area), String(r.representation), complexity)
end

function _classification_text(value, field::AbstractString)
    text = strip(String(value))
    isempty(text) && throw(ArgumentError("$field must not be empty"))
    text
end

function _classification_planar(value)
    value isa Bool && return value
    value isa Integer && value in (0, 1) && return Bool(value)
    throw(ArgumentError("planar must be Bool or SQLite 0/1"))
end

function _classification_status(value)
    text = lowercase(strip(String(value)))
    text in ("classified", "unresolved", "nonclassical") ||
        throw(ArgumentError("classification_status must be classified, unresolved, or nonclassical"))
    text
end

function _classification_knot_type(value)
    if value === nothing || value === missing
        return nothing
    end
    text = String(value)
    isempty(strip(text)) && throw(ArgumentError("knot_type must not be empty"))
    text
end

function _validate_classification_fields(planar, knot_type, classification_status)
    p = _classification_planar(planar)
    status = _classification_status(classification_status)
    type = _classification_knot_type(knot_type)
    if !p
        type === nothing && status == "nonclassical" ||
            throw(ArgumentError("planar=false requires knot_type=NULL and status=nonclassical"))
    elseif status == "classified"
        type === nothing && throw(ArgumentError("classified planar state requires knot_type"))
    elseif status == "unresolved"
        type === nothing || throw(ArgumentError("unresolved state requires knot_type=NULL"))
    else
        throw(ArgumentError("planar=true cannot have status=nonclassical"))
    end
    p, type, status
end

function _classification_run_from_row(row)
    row === nothing && return nothing
    ClassificationRun(
        String(row.method_version), Int(row.canonical_version), Int(row.max_n),
        Int(row.barrier), String(row.census_fingerprint),
        String(row.classification_fingerprint), JSON.parse(String(row.counts_json)),
        JSON.parse(String(row.conventions_json)), String(row.created_at_utc))
end

"""
    get_classification_run(s::MoveGraphStoreDB, method_version::AbstractString) -> Union{Nothing,ClassificationRun}

Read the stored metadata of one classification run.

# Arguments
Requires M4 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

`s` is an open store with the classification schema. `method_version`
is the exact stored method key; no whitespace normalization is performed.
Database errors, field conversion errors and invalid stored JSON propagate.


# Returns
A `ClassificationRun` with method and canonical versions, census bound,
barrier, census and classification fingerprints, decoded counts and
conventions, and creation time; `nothing` if the key is absent.

# Guarantees
1. The metadata are those persisted for the requested method, with counts
   and conventions decoded from stored JSON. [proved/agent: by construction]
2. The store is not modified. [proved/agent: by construction]

# Not guaranteed
Fingerprints, counts, conventions and timestamp are not revalidated.
Stored metadata do not certify a complete census or correct classification.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql classification_runs — dataset-level classification metadata.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> r = get_classification_run(s, "fixture-v2"); (r.max_n, r.barrier)
(7, 9)

julia> close_store!(s);
```

Role: access
"""
function get_classification_run(s::MoveGraphStoreDB, method_version::AbstractString)
    row = _one(s.db, "SELECT * FROM classification_runs WHERE method_version=?",
               (String(method_version),))
    _classification_run_from_row(row)
end

"""
    classification_run(s::MoveGraphStoreDB, method_version::AbstractString) -> Union{Nothing,ClassificationRun}

Read the stored metadata of one classification run.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is an open store with the classification schema. `method_version`
is the exact stored method key; no whitespace normalization is performed.
Database errors, field conversion errors and invalid stored JSON propagate.
This accepts the calls of [`get_classification_run`](@ref).

# Returns
A `ClassificationRun` with method and canonical versions, census bound,
barrier, census and classification fingerprints, decoded counts and
conventions, and creation time; `nothing` if the key is absent.

# Guarantees
1. The metadata are those persisted for the requested method, with counts
   and conventions decoded from stored JSON. [proved/agent: by construction]
2. The store is not modified. [proved/agent: by construction]

# Not guaranteed
Fingerprints, counts, conventions and timestamp are not revalidated.
Stored metadata do not certify a complete census or correct classification.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql classification_runs — dataset-level classification metadata.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> r = classification_run(s, "fixture-v2"); (r.max_n, r.barrier)
(7, 9)

julia> close_store!(s);
```

Role: access
"""
classification_run(args...; kwargs...) = get_classification_run(args...; kwargs...)

"""
    register_classification_run!(s::MoveGraphStoreDB, method_version::AbstractString; canonical_version::Integer, max_n::Integer, barrier::Integer, census_fingerprint::AbstractString, classification_fingerprint::AbstractString, counts=Dict{String,Any}(), conventions=Dict{String,Any}(), created_at_utc=string(Dates.now(Dates.UTC))) -> ClassificationRun

Register immutable provenance for a finite classification dataset.

# Arguments
Requires M4 and M13 ([conditions](../specs/move-graph-store.md#numbered-conditions)).
`s` is open, writable, migrated and outside a transaction. Text and bound
checks raise `ArgumentError` naming the invalid field. A conflicting existing
method raises `ArgumentError` about different metadata. Integer conversion,
JSON and database failures propagate. Truth of metadata is caller-supplied.

# Returns
The persisted `ClassificationRun`, including decoded `counts` and
`conventions`, as defined in [Schema and classifications](../specs/move-graph-store.md#schema-and-classifications).

# Guarantees
1. A new method key stores the supplied metadata. An identical repeated registration returns the existing row and retains its original timestamp. [proved/agent: by construction]
2. A conflicting row is not overwritten. Database writes commit together or roll back on failure; input containers are not modified. [proved/agent: by construction]

# Not guaranteed
Counts, fingerprints and timestamp are not recomputed or authenticated.
This does not run a census, decide a knot type or prove completeness.
JSON text equality can distinguish differently encoded equivalent metadata.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md M4 and M13 — classification provenance and immutable keys.

# Examples
```jldoctest
julia> s = open_store();

julia> r = register_classification_run!(s, "fixture-v2"; canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> (r.max_n, r.barrier)
(7, 9)

julia> close_store!(s);
```

Role: other
"""
function register_classification_run!(s::MoveGraphStoreDB,
                                      method_version::AbstractString;
                                      canonical_version::Integer,
                                      max_n::Integer,
                                      barrier::Integer,
                                      census_fingerprint::AbstractString,
                                      classification_fingerprint::AbstractString,
                                      counts=Dict{String,Any}(),
                                      conventions=Dict{String,Any}(),
                                      created_at_utc::AbstractString=_classification_now_utc())
    method = _classification_text(method_version, "method_version")
    canon = Int(canonical_version)
    maximum_n = Int(max_n)
    barrier_n = Int(barrier)
    canon > 0 || throw(ArgumentError("canonical_version must be positive"))
    maximum_n >= 0 || throw(ArgumentError("max_n must be nonnegative"))
    barrier_n >= maximum_n || throw(ArgumentError("require barrier >= max_n"))
    census = _classification_text(census_fingerprint, "census_fingerprint")
    classification = _classification_text(classification_fingerprint,
                                           "classification_fingerprint")
    counts_json = JSON.json(counts)
    conventions_json = JSON.json(conventions)
    created = _classification_text(created_at_utc, "created_at_utc")
    _transaction(s.db) do
        existing = _one(s.db,
                        "SELECT * FROM classification_runs WHERE method_version=?",
                        (method,))
        if existing === nothing
            _dbexec(s.db, """
                INSERT INTO classification_runs(
                  method_version,canonical_version,max_n,barrier,census_fingerprint,
                  classification_fingerprint,counts_json,conventions_json,created_at_utc)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, (method, canon, maximum_n, barrier_n, census, classification,
                  counts_json, conventions_json, created))
        else
            expected = (canon, maximum_n, barrier_n, census, classification,
                        counts_json, conventions_json)
            actual = (Int(existing.canonical_version), Int(existing.max_n),
                      Int(existing.barrier), String(existing.census_fingerprint),
                      String(existing.classification_fingerprint),
                      String(existing.counts_json), String(existing.conventions_json))
            actual == expected || throw(ArgumentError(
                "classification run $method already exists with different metadata"))
        end
    end
    get_classification_run(s, method)
end

"""
    register_classification_run!(s::MoveGraphStoreDB; method_version, kwargs...) -> ClassificationRun

Register immutable provenance for a finite classification dataset.

# Arguments
Requires M4 and M13 ([conditions](../specs/move-graph-store.md#numbered-conditions)).
`s` is open, writable, migrated and outside a transaction. Text and bound
checks raise `ArgumentError` naming the invalid field. A conflicting existing
method raises `ArgumentError` about different metadata. Integer conversion,
JSON and database failures propagate. Truth of metadata is caller-supplied.

# Returns
The persisted `ClassificationRun`, including decoded `counts` and
`conventions`, as defined in [Schema and classifications](../specs/move-graph-store.md#schema-and-classifications).

# Guarantees
1. A new method key stores the supplied metadata. An identical repeated registration returns the existing row and retains its original timestamp. [proved/agent: by construction]
2. A conflicting row is not overwritten. Database writes commit together or roll back on failure; input containers are not modified. [proved/agent: by construction]

# Not guaranteed
Counts, fingerprints and timestamp are not recomputed or authenticated.
This does not run a census, decide a knot type or prove completeness.
JSON text equality can distinguish differently encoded equivalent metadata.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md M4 and M13 — classification provenance and immutable keys.

# Examples
```jldoctest
julia> s = open_store();

julia> r = register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> (r.max_n, r.barrier)
(7, 9)

julia> close_store!(s);
```

Role: other
"""
register_classification_run!(s::MoveGraphStoreDB; method_version, kwargs...) =
    register_classification_run!(s, method_version; kwargs...)

function _classification_state_id(state)
    state isa CanonicalState && return state.state_id
    String(state)
end

function _classification_row(row)
    row === nothing && return nothing
    (state_id=String(row.state_id), planar=Int(row.planar) == 1,
     knot_type=(row.knot_type === nothing || row.knot_type === missing) ? nothing :
               String(row.knot_type),
     classification_status=String(row.classification_status),
     method_version=String(row.method_version))
end

"""
    put_state_classification!(s::MoveGraphStoreDB, state; planar, knot_type=nothing, classification_status, method_version::AbstractString) -> NamedTuple

Record an immutable classification assertion for a state under a registered method and return the stored classification record.

# Arguments
Requires M4 and M14 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

The method version may also be supplied as the third positional argument. `s` is open and writable outside a transaction; `state`
is a stored `CanonicalState` or identifier. The method must already be
registered and the state must belong to `AREA_2`. The planarity, knot-type
and status fields must be mutually consistent. These checks and conflicting
updates raise `ArgumentError` identifying the missing record, domain mismatch,
invalid field or conflicting values. Database/conversion errors propagate.
The caller supplies the mathematical truth of the assertion.

# Returns
A `NamedTuple` with `state_id::String`, `planar::Bool`, optional
`knot_type::String`, `classification_status::String`, and `method_version::String`.

# Guarantees
1. Success stores the requested assertion under the state and method key;
   an identical repeat preserves the existing row. [proved/agent: by construction]
2. Conflicting values cannot replace that assertion through this API.
   [proved/agent: by construction]
3. The store is modified; the supplied state is not modified. Database
   changes of a failed transaction roll back. [proved/agent: by construction]

# Not guaranteed
This does not decide planarity, identify a knot, or verify a classification
certificate. Agreement of field values is a storage consistency check.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql state_classifications — assertions indexed by state and classification method.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a);

julia> put_state_classification!(s, a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2").knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: other
"""
function put_state_classification!(s::MoveGraphStoreDB, state;
                                   planar,
                                   knot_type=nothing,
                                   classification_status,
                                   method_version::AbstractString)
    p, type, status = _validate_classification_fields(planar, knot_type,
                                                      classification_status)
    method = _classification_text(method_version, "method_version")
    sid = _classification_state_id(state)
    _transaction(s.db) do
        get_classification_run(s, method) === nothing && throw(ArgumentError(
            "classification run $method does not exist"))
        state_record = get_state(s, sid)
        state_record === nothing && throw(ArgumentError(
            "state $sid does not exist"))
        state_record.area == AREA_2 || throw(ArgumentError(
            "state classifications require Area 2 states"))
        existing = _one(s.db, """
            SELECT * FROM state_classifications
            WHERE method_version=? AND state_id=?
        """, (method, sid))
        if existing === nothing
            _dbexec(s.db, """
                INSERT INTO state_classifications(
                  state_id,planar,knot_type,classification_status,method_version)
                VALUES (?,?,?,?,?)
            """, (sid, p ? 1 : 0, type, status, method))
        else
            _classification_row(existing) ==
                (state_id=sid, planar=p, knot_type=type,
                 classification_status=status, method_version=method) ||
                throw(ArgumentError("state classification already exists with different values"))
        end
    end
    get_state_classification(s, sid; method_version=method)
end

"""
    put_state_classification!(s::MoveGraphStoreDB, state, method_version::AbstractString; planar, knot_type=nothing, classification_status) -> NamedTuple

Record an immutable classification assertion for a state under a registered method and return the stored classification record.

# Arguments
Requires M4 and M14 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

The method version may also be supplied as the third positional argument. `s` is open and writable outside a transaction; `state`
is a stored `CanonicalState` or identifier. The method must already be
registered and the state must belong to `AREA_2`. The planarity, knot-type
and status fields must be mutually consistent. These checks and conflicting
updates raise `ArgumentError` identifying the missing record, domain mismatch,
invalid field or conflicting values. Database/conversion errors propagate.
The caller supplies the mathematical truth of the assertion.

# Returns
A `NamedTuple` with `state_id::String`, `planar::Bool`, optional
`knot_type::String`, `classification_status::String`, and `method_version::String`.

# Guarantees
1. Success stores the requested assertion under the state and method key;
   an identical repeat preserves the existing row. [proved/agent: by construction]
2. Conflicting values cannot replace that assertion through this API.
   [proved/agent: by construction]
3. The store is modified; the supplied state is not modified. Database
   changes of a failed transaction roll back. [proved/agent: by construction]

# Not guaranteed
This does not decide planarity, identify a knot, or verify a classification
certificate. Agreement of field values is a storage consistency check.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql state_classifications — assertions indexed by state and classification method.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a);

julia> put_state_classification!(s, a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2").knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: other
"""
put_state_classification!(s::MoveGraphStoreDB, state, method_version::AbstractString;
                          kwargs...) =
    put_state_classification!(s, state; method_version=method_version, kwargs...)
"""
    record_state_classification!(s::MoveGraphStoreDB, state; planar, knot_type=nothing, classification_status, method_version::AbstractString) -> NamedTuple

Record an immutable classification assertion for a state under a registered method and return the stored classification record.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

Calls are those of [`put_state_classification!`](@ref), including the positional
method-version form. `s` is open and writable outside a transaction; `state`
is a stored `CanonicalState` or identifier. The method must already be
registered and the state must belong to `AREA_2`. The planarity, knot-type
and status fields must be mutually consistent. These checks and conflicting
updates raise `ArgumentError` identifying the missing record, domain mismatch,
invalid field or conflicting values. Database/conversion errors propagate.
The caller supplies the mathematical truth of the assertion.

# Returns
A `NamedTuple` with `state_id::String`, `planar::Bool`, optional
`knot_type::String`, `classification_status::String`, and `method_version::String`.

# Guarantees
1. Success stores the requested assertion under the state and method key;
   an identical repeat preserves the existing row. [proved/agent: by construction]
2. Conflicting values cannot replace that assertion through this API.
   [proved/agent: by construction]
3. The store is modified; the supplied state is not modified. Database
   changes of a failed transaction roll back. [proved/agent: by construction]

# Not guaranteed
This does not decide planarity, identify a knot, or verify a classification
certificate. Agreement of field values is a storage consistency check.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql state_classifications — assertions indexed by state and classification method.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a);

julia> record_state_classification!(s, a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2").knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: other
"""
record_state_classification!(args...; kwargs...) =
    put_state_classification!(args...; kwargs...)

"""
    get_state_classification(s::MoveGraphStoreDB, state; method_version=nothing) -> Union{Nothing,NamedTuple}

Read a persisted classification assertion for a state.

# Arguments
Requires M4 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

`s` is open with the common-store schema; `state` is a `CanonicalState` or
identifier convertible to `String`. Supply `method_version` to select an
exact stored method key. Omitting it is accepted only when at most one
assertion exists for that state; otherwise `ArgumentError` requests the
method version. Database/conversion errors propagate. 

# Returns
`nothing` if absent, otherwise a `NamedTuple` with state identifier,
planarity flag, optional knot-type name, status and method version.

# Guarantees
1. The result contains the stored fields for the selected assertion.
   [proved/agent: by construction]
2. The store and supplied state are not modified. [proved/agent: by construction]

# Not guaranteed
Stored classification or planarity claims are not recomputed or verified.
Missing data do not establish that the state is unclassifiable.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql state_classifications — stored method-specific assertions.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a); record_state_classification!(s, a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2");

julia> get_state_classification(s, a).knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: access
"""
function get_state_classification(s::MoveGraphStoreDB, state;
                                  method_version=nothing)
    sid = _classification_state_id(state)
    row = if method_version === nothing
        rows = _rows(s.db, """
            SELECT * FROM state_classifications
            WHERE state_id=? ORDER BY method_version
        """, (sid,))
        isempty(rows) && return nothing
        length(rows) == 1 || throw(ArgumentError(
            "method_version is required when a state has multiple classifications"))
        first(rows)
    else
        _one(s.db, """
            SELECT * FROM state_classifications
            WHERE method_version=? AND state_id=?
        """, (String(method_version), sid))
    end
    _classification_row(row)
end

"""
    state_classification(s::MoveGraphStoreDB, state; method_version=nothing) -> Union{Nothing,NamedTuple}

Read a persisted classification assertion for a state.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema; `state` is a `CanonicalState` or
identifier convertible to `String`. Supply `method_version` to select an
exact stored method key. Omitting it is accepted only when at most one
assertion exists for that state; otherwise `ArgumentError` requests the
method version. Database/conversion errors propagate. Calls are those of
[`get_state_classification`](@ref).

# Returns
`nothing` if absent, otherwise a `NamedTuple` with state identifier,
planarity flag, optional knot-type name, status and method version.

# Guarantees
1. The result contains the stored fields for the selected assertion.
   [proved/agent: by construction]
2. The store and supplied state are not modified. [proved/agent: by construction]

# Not guaranteed
Stored classification or planarity claims are not recomputed or verified.
Missing data do not establish that the state is unclassifiable.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql state_classifications — stored method-specific assertions.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a); record_state_classification!(s, a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2");

julia> state_classification(s, a).knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: access
"""
state_classification(args...; kwargs...) =
    get_state_classification(args...; kwargs...)

"""
    classification_rows(s::MoveGraphStoreDB, method_version::AbstractString) -> Vector

Read the stored classification assertions for one method in state-ID order.

# Arguments
`s` is open with the common classification schema. The exact method string
selects M4 assertions ([conditions](../specs/move-graph-store.md#numbered-conditions)).
Method existence is not checked. Database and field conversion errors
propagate. Valid stored fields are a caller obligation for meaningful results.

# Returns
A vector of records with `state_id`, `planar`, optional `knot_type`,
`classification_status` and `method_version`, as defined in
[Schema and classifications](../specs/move-graph-store.md#schema-and-classifications).
An absent method or a method with no assertions gives an empty vector.

# Guarantees
1. Results contain exactly the stored assertions with the selected method key, sorted by state ID. [proved/agent: by construction]
2. The store and arguments are not modified. [proved/agent: by construction]

# Not guaranteed
Stored assertions are not independently verified knot classifications.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md Schema and classifications — method-specific assertions.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s,a); put_state_classification!(s,a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2");

julia> only(classification_rows(s, "fixture-v2")).knot_type == "future-6_1"
true

julia> close_store!(s);
```

Role: access
"""
function classification_rows(s::MoveGraphStoreDB, method_version::AbstractString)
    rows = _rows(s.db, """
        SELECT * FROM state_classifications
        WHERE method_version=? ORDER BY state_id
    """, (String(method_version),))
    _classification_row.(rows)
end

"""
    classification_counts(s::MoveGraphStoreDB, method_version::AbstractString) -> Dict{String,Any}

Count the stored classification assertions for one method.

# Arguments
`s` is open with the common classification schema. The exact method string
selects M4 assertions ([conditions](../specs/move-graph-store.md#numbered-conditions)).
Method existence is not checked. Database and field conversion errors
propagate. Valid stored fields are a caller obligation for meaningful results.

# Returns
A dictionary reporting `total`, `planar_true`, `planar_false`, `classified`,
`unresolved`, `nonclassical` and `knot_type_counts`; the last is a
string-to-integer dictionary for non-null knot names. These count the
[classification assertions](../specs/move-graph-store.md#schema-and-classifications).

# Guarantees
1. Counts describe the stored assertions with the selected method key, independently of the counts declared in the run's metadata. [proved/agent: by construction]
2. The store and arguments are not modified. [proved/agent: by construction]

# Not guaranteed
Counts do not independently certify planarity, knot type, or a complete census.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md Schema and classifications — method-specific assertions.

# Examples
```jldoctest
julia> s = open_store(); register_classification_run!(s; method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9, census_fingerprint="census-fixture", classification_fingerprint="class-fixture");

julia> a = CanonicalState("fixture-area2"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s,a); put_state_classification!(s,a; planar=true, knot_type="future-6_1", classification_status="classified", method_version="fixture-v2");

julia> classification_counts(s, "fixture-v2")["total"]
1

julia> close_store!(s);
```

Role: access
"""
function classification_counts(s::MoveGraphStoreDB, method_version::AbstractString)
    method = String(method_version)
    total_row = _one(s.db, """
        SELECT COUNT(*) AS total,
               SUM(planar = 1) AS planar_true,
               SUM(planar = 0) AS planar_false
        FROM state_classifications WHERE method_version=?
    """, (method,))
    total = Int(total_row.total)
    planar_true = total_row.planar_true === nothing ? 0 : Int(total_row.planar_true)
    planar_false = total_row.planar_false === nothing ? 0 : Int(total_row.planar_false)
    status = Dict{String,Int}(String(row.classification_status) => Int(row.count)
                              for row in _rows(s.db, """
        SELECT classification_status, COUNT(*) AS count
        FROM state_classifications WHERE method_version=?
        GROUP BY classification_status
    """, (method,)))
    for key in ("classified", "unresolved", "nonclassical")
        get!(status, key, 0)
    end
    types = Dict{String,Int}(String(row.knot_type) => Int(row.count)
                             for row in _rows(s.db, """
        SELECT knot_type, COUNT(*) AS count
        FROM state_classifications
        WHERE method_version=? AND knot_type IS NOT NULL
        GROUP BY knot_type ORDER BY knot_type
    """, (method,)))
    Dict{String,Any}(
        "total" => total, "planar_true" => planar_true,
        "planar_false" => planar_false,
        "classified" => status["classified"],
        "unresolved" => status["unresolved"],
        "nonclassical" => status["nonclassical"],
        "knot_type_counts" => types)
end

"""
    register_move_type!(s::MoveGraphStoreDB, move_type; area=nothing, inverse_type=move_type, inverse_direction=:inverse) -> MoveType

Register a recognized move name and its declared inverse metadata.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open and writable with the common-store schema. `move_type` and
`inverse_type` convert to `String`. Omitted `area` is obtained from
[`move_area`](@ref); a supplied area must agree. The direction must convert
to `:forward` or `:inverse` ignoring case. `ArgumentError` identifies an
unknown name, mismatched area or invalid direction. Conversion and database
errors propagate. The inverse name itself is not validated.

# Returns
The stored `MoveType` with name, domain and declared inverse metadata.

# Guarantees
1. The exact move-name text is registered; a repeat updates its metadata.
   [proved/agent: by construction]
2. The accepted name belongs to the supplied or inferred domain according
   to the taxonomy. [proved/agent: by construction]
3. The store is modified; arguments are not modified.
   [proved/agent: by construction]

# Not guaranteed
No move is applied or checked on diagrams, and the declared inverse is not
proved to be an inverse. Lookup keys retain their supplied case/whitespace
even though name validation ignores those differences.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql move_types — registration of named moves.

# Examples
```jldoctest
julia> s = open_store();

julia> register_move_type!(s, "MP"; area=Topology.AREA_3).area == Topology.AREA_3
true

julia> close_store!(s);
```

Role: other
"""
function register_move_type!(s::MoveGraphStoreDB, move_type; area=nothing,
                             inverse_type=move_type, inverse_direction=:inverse)
    a = area === nothing ? move_area(move_type) : normalize_area(area)
    validate_topological_move(move_type, a)
    d = Symbol(lowercase(String(inverse_direction)))
    d in (:forward, :inverse) || throw(ArgumentError("inverse_direction must be :forward or :inverse"))
    _dbexec(s.db, """
      INSERT INTO move_types(move_type,area,inverse_type,inverse_direction) VALUES (?,?,?,?)
      ON CONFLICT(move_type) DO UPDATE SET area=excluded.area,
        inverse_type=excluded.inverse_type,inverse_direction=excluded.inverse_direction
    """, (String(move_type), a, String(inverse_type), String(d)))
    get_move_type(s, move_type)
end

"""
    get_move_type(s::MoveGraphStoreDB, move_type) -> Union{Nothing,MoveType}

Read the registered metadata of a move name.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` must be open with the common-store schema; `move_type` converts to
`String`. Lookup uses that exact text, without the taxonomy's case or
whitespace normalization. Connection and conversion errors propagate.

# Returns
The stored `MoveType` containing its name, area, inverse name and inverse
direction, or `nothing` if the exact name is absent.

# Guarantees
1. Fields reproduce the stored metadata after conversion to the `MoveType`
   field types. [proved/agent: by construction]
2. The database is not modified. [proved/agent: by construction]

# Not guaranteed
Lookup does not recheck taxonomy membership or verify an inverse move on
diagrams or spaces. An unregistered alias is not resolved automatically.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql move_types — move registration metadata.

# Examples
```jldoctest
julia> s = open_store(); register_move_type!(s, "MP"; area=Topology.AREA_3);

julia> get_move_type(s, "MP").area == Topology.AREA_3
true

julia> close_store!(s);
```

Role: access
"""
function get_move_type(s::MoveGraphStoreDB, move_type)
    r = _one(s.db, "SELECT * FROM move_types WHERE move_type=?", (String(move_type),))
    r === nothing && return nothing
    MoveType(String(r.move_type), Float64(r.area), String(r.inverse_type), Symbol(r.inverse_direction))
end

function _ensure_move!(s, move_type, area)
    m = get_move_type(s, move_type)
    m === nothing ? register_move_type!(s, move_type; area=area) :
        (validate_topological_move(m.move_type, area); m.area == normalize_area(area) ||
            throw(ArgumentError("registered move type belongs to a different area")); m)
end

function _verification_result(x)
    x === true && return true
    x === false && return false
    x isa NamedTuple && haskey(x, :verified) && return Bool(x.verified)
    x isa AbstractDict && (haskey(x, :verified) || haskey(x, "verified")) &&
        return Bool(get(x, haskey(x, :verified) ? :verified : "verified", false))
    hasproperty(x, :verified) && return Bool(getproperty(x, :verified))
    false
end

function _run_verifier(verifier, source, target, move)
    verifier === nothing && return false
    # `applicable` distinguishes an unsupported arity from a MethodError
    # raised by the verifier's own body.  Never catch the latter: a verifier
    # failure must not be retried at a weaker arity and accidentally certify.
    result = if applicable(verifier, source, target, move)
        verifier(source, target, move)
    elseif applicable(verifier, source, target)
        verifier(source, target)
    elseif applicable(verifier)
        verifier()
    else
        return false
    end
    _verification_result(result)
end

function _json_or_nothing(x)
    x === nothing && return nothing
    x isa AbstractString ? String(x) : JSON.json(x)
end

function _transaction(action::Function, db)
    _dbexec(db, "BEGIN IMMEDIATE")
    try
        result = action()
        _dbexec(db, "COMMIT")
        result
    catch err
        _dbexec(db, "ROLLBACK")
        Base.throw(err)
    end
end

"""
    insert_verified_transition!(s::MoveGraphStoreDB, source::CanonicalState, target::CanonicalState, move_type; kwargs...) -> Int

Persist a directed move-graph edge authorized by caller-supplied verification and return the stored transition identifier.

# Arguments
Requires M3 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

Keywords are `direction=:forward`, `location=nothing`, `parameters=nothing`,
`engine_version="unknown"`, `verifier_version="unknown"`, `verifier=nothing`,
`verification_result=nothing`, `verified=false`, `certificate_paths=String[]`,
and `evidence=nothing`.
`s` must be open and writable outside a transaction; states must have the
same area and the move name must belong to that area. `direction` must be
`:forward` or `:inverse`. Supply `verified=true` and either a verifier callback
or an accepted verification result; a flag alone is rejected. These checks
raise `ArgumentError` describing the domain, direction or verification failure.
A callback may accept `(source, target, move)`, `(source, target)`, or no
arguments, in that preference order; callback errors propagate without retry.
Optional location, parameters, evidence and certificate data must be encodable
as text/JSON. Empty persisted evidence and other constraint violations raise
`SQLite.SQLiteException`. Correct mathematical verification is the caller's
responsibility.

# Returns
The stored transition identifier as `Int`, including for a duplicate edge.

# Guarantees
1. Success stores the endpoints and an edge accepted by the supplied
   verification, with nonempty stored evidence. [proved/agent: by construction]
2. A duplicate with the same endpoints, move, direction and encoded parameters
   retains existing edge metadata; endpoint payloads/complexities can update
   and certificate associations can be added. [proved/agent: by construction]
3. Database changes belonging to this operation commit together or roll back
   together on failure. The API does not modify input state objects; effects
   of a caller's callback are outside that promise. [proved/agent: by construction]

# Not guaranteed
The store trusts the supplied verification and does not independently prove
equivalence of represented spaces. A callback returning true can admit an
invalid mathematical move. An inverse direction label does not add a reverse
edge. Callback side effects are not rolled back.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql transitions and evidence triggers — stored, verifier-authorized directed edges.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> insert_verified_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true) isa Integer
true

julia> close_store!(s);
```

Role: other
""" 
function insert_verified_transition!(s::MoveGraphStoreDB, source::CanonicalState,
                                     target::CanonicalState, move_type;
                                     direction=:forward, location=nothing, parameters=nothing,
                                     engine_version="unknown", verifier_version="unknown",
                                     verifier=nothing, verification_result=nothing,
                                     verified=false, certificate_paths=String[], evidence=nothing)
    source.area == target.area || throw(ArgumentError("source and target areas differ"))
    move_name = String(move_type)
    stored_move = get_move_type(s, move_name)
    m = if stored_move === nothing
        area = normalize_area(source.area)
        validate_topological_move(move_name, area)
        MoveType(move_name, area, move_name, :inverse)
    else
        validate_topological_move(stored_move.move_type, source.area)
        stored_move.area == normalize_area(source.area) || throw(ArgumentError(
            "registered move type belongs to a different area"))
        stored_move
    end
    d = Symbol(lowercase(String(direction)))
    d in (:forward, :inverse) || throw(ArgumentError("direction must be :forward or :inverse"))
    allowed = verifier === nothing ? _verification_result(verification_result) :
              _run_verifier(verifier, source, target, m)
    allowed || throw(ArgumentError("transition verifier did not certify this edge"))
    verified === true || throw(ArgumentError("verified=true is not evidence; pass a verifier callback/result"))
    stored_parameters = parameters === nothing ? "" : _json_or_nothing(parameters)
    stored_evidence = _json_or_nothing(evidence === nothing ? verification_result : evidence)
    # A callback is itself the independently executed verifier witness when
    # the caller did not provide a richer structured result.  Keep a compact
    # nonempty marker so the SQL evidence trigger cannot be bypassed by the
    # verified flag alone.
    stored_evidence === nothing && verifier !== nothing &&
        (stored_evidence = JSON.json(Dict("verifier_callback" => true)))
    tid_ref = Ref{Any}(nothing)
    _transaction(s.db) do
        # Registration is part of the same transaction as the edge.  A failed
        # verifier above therefore cannot leave a new move type behind.
        _ensure_move!(s, move_name, source.area)
        put_state!(s, source); put_state!(s, target)
        _dbexec(s.db, """
          INSERT INTO transitions(source_state_id,target_state_id,move_type,direction,location,parameters,engine_version,verifier_version,verification_evidence,verified)
          VALUES (?,?,?,?,?,?,?,?,?,1)
          ON CONFLICT(source_state_id,target_state_id,move_type,direction,parameters) DO NOTHING
        """, (source.state_id,target.state_id,move_name,String(d),_json_or_nothing(location),
               stored_parameters,String(engine_version),String(verifier_version),stored_evidence))
        row = _one(s.db, "SELECT * FROM transitions WHERE source_state_id=? AND target_state_id=? AND move_type=? AND direction=? AND parameters IS ?",
                    (source.state_id,target.state_id,move_name,String(d),stored_parameters))
        tid = Int(row.id)
        tid_ref[] = tid
        for path in certificate_paths
            add_certificate!(s, tid, path)
        tid
    end
    tid_ref[]
end
end
"""
    record_transition!(s::MoveGraphStoreDB, source::CanonicalState, target::CanonicalState, move_type; kwargs...) -> Int

Persist a directed move-graph edge authorized by caller-supplied verification and return the stored transition identifier.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

Calls and keyword defaults are those of [`insert_verified_transition!`](@ref).
`s` must be open and writable outside a transaction; states must have the
same area and the move name must belong to that area. `direction` must be
`:forward` or `:inverse`. Supply `verified=true` and either a verifier callback
or an accepted verification result; a flag alone is rejected. These checks
raise `ArgumentError` describing the domain, direction or verification failure.
A callback may accept `(source, target, move)`, `(source, target)`, or no
arguments, in that preference order; callback errors propagate without retry.
Optional location, parameters, evidence and certificate data must be encodable
as text/JSON. Empty persisted evidence and other constraint violations raise
`SQLite.SQLiteException`. Correct mathematical verification is the caller's
responsibility.

# Returns
The stored transition identifier as `Int`, including for a duplicate edge.

# Guarantees
1. Success stores the endpoints and an edge accepted by the supplied
   verification, with nonempty stored evidence. [proved/agent: by construction]
2. A duplicate with the same endpoints, move, direction and encoded parameters
   retains existing edge metadata; endpoint payloads/complexities can update
   and certificate associations can be added. [proved/agent: by construction]
3. Database changes belonging to this operation commit together or roll back
   together on failure. The API does not modify input state objects; effects
   of a caller's callback are outside that promise. [proved/agent: by construction]

# Not guaranteed
The store trusts the supplied verification and does not independently prove
equivalence of represented spaces. A callback returning true can admit an
invalid mathematical move. An inverse direction label does not add a reverse
edge. Callback side effects are not rolled back.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql transitions and evidence triggers — stored, verifier-authorized directed edges.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> record_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true) isa Integer
true

julia> close_store!(s);
```

Role: other
"""
record_transition!(args...; kwargs...) = insert_verified_transition!(args...; kwargs...)

"""
    add_certificate!(s::MoveGraphStoreDB, transition_id::Integer, path) -> Int

Attach a certificate reference or certificate data to a stored transition and return its transition identifier.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` must be open with the common-store schema; `transition_id` must identify
an existing edge. The foreign-key check raises `SQLite.SQLiteException`
for an absent edge when foreign keys are enabled, as by [`open_store`](@ref).
`path` is text or JSON-serializable data. Conversion/serialization failures
propagate; a null path is ignored by the insertion's conflict policy.

# Returns
The supplied identifier converted to `Int`, including for an existing
association.

# Guarantees
1. A non-null, successfully encoded association is present after success;
   a duplicate adds no duplicate row. [proved/agent: by construction]
2. The store is modified; supplied certificate data are not modified.
   [proved/agent: by construction]

# Not guaranteed
No path is opened and no certificate is checked. An association alone does
not prove that its transition preserves any mathematical object.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql certificate_paths — associations with stored transitions, not a proof checker.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> id = record_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true, certificate_paths=["cert/a-b"]);

julia> add_certificate!(s, id, "cert/a-b") == id && length(Topology.MoveGraphStore.neighbors(s, a.state_id)) == 1
true

julia> close_store!(s);
```

Role: other
"""
function add_certificate!(s::MoveGraphStoreDB, transition_id::Integer, path)
    _dbexec(s.db, "INSERT OR IGNORE INTO certificate_paths(transition_id,path) VALUES (?,?)",
            (Int(transition_id), _json_or_nothing(path)))
    Int(transition_id)
end

"""
    add_invariant!(s::MoveGraphStoreDB, state, name, value; source=nothing) -> SQLite.Query

Attach a caller-supplied invariant value to a stored state and return the insertion's `SQLite.Query`.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema. `state` is a `CanonicalState` or
an identifier convertible to `String` and must already be stored; a missing
state raises `SQLite.SQLiteException` when foreign keys are enabled.
`name` and optional `source` convert to `String`; `value` is text or
JSON-serializable data. Conversion/encoding errors propagate. The caller
supplies the invariant's meaning and correctness.

# Returns
The `SQLite.Query` for the insertion, not the stored invariant value or a
mathematical verification result.

# Guarantees
1. A non-null encoded value is associated with the named state and invariant;
   an existing identical association retains its existing source.
   [proved/agent: by construction]
2. The database is modified; the supplied state and value are not modified.
   [proved/agent: by construction]

# Not guaranteed
The invariant is not computed or checked. Conflicting values of the same
named invariant may coexist. A null value is ignored by the insertion's
conflict policy; the returned query does not certify that a row was added.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql invariants — stored assertions associated with states.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); put_state!(s, a);

julia> add_invariant!(s, a, "fixture", 1);

julia> get_state(s, a.state_id).state_id == a.state_id
true

julia> close_store!(s);
```

Role: other
"""
function add_invariant!(s::MoveGraphStoreDB, state, name, value; source=nothing)
    sid = state isa CanonicalState ? state.state_id : String(state)
    _dbexec(s.db, "INSERT OR IGNORE INTO invariants(state_id,name,value,source) VALUES (?,?,?,?)",
            (sid,String(name),_json_or_nothing(value),source === nothing ? nothing : String(source)))
end

"""
    neighbors(s::MoveGraphStoreDB, sid; direction=:outgoing) -> Vector

Return the stored directed transitions incident to a state.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema; `sid` converts to a state identifier
string. `direction=:incoming` selects edges targeting it, `:both` selects
either endpoint, and every other value selects outgoing edges. The selector
is not validated. No state-membership check is performed. Database and
conversion errors propagate.

# Returns
A vector of materialized row `NamedTuple`s, in transition-identifier order.
Each row contains the edge identifier, endpoint identifiers, move, direction,
location, parameters, engine/verifier versions, evidence and verified flag.
An isolated or absent identifier gives an empty vector.

# Guarantees
1. Exactly the stored rows matching the selected incidence are returned;
   `:both` returns a self-loop only once. [proved/agent: by construction]
2. The store is not modified. [proved/agent: by construction]

# Not guaranteed
Rows are transitions, not distinct neighboring states. Parallel moves can
reach the same target. Evidence is not replayed and mathematical validity
is not checked by this query.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql transitions — directed edges and their metadata.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> record_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true);

julia> length(Topology.MoveGraphStore.neighbors(s, a.state_id))
1

julia> close_store!(s);
```

Role: access
"""
function neighbors(s::MoveGraphStoreDB, sid; direction=:outgoing)
    id = String(sid)
    q = direction == :incoming ? "SELECT t.* FROM transitions t WHERE target_state_id=? ORDER BY t.id" :
        direction == :both ? "SELECT t.* FROM transitions t WHERE source_state_id=? OR target_state_id=? ORDER BY t.id" :
        "SELECT t.* FROM transitions t WHERE source_state_id=? ORDER BY t.id"
    direction == :both ? _rows(s.db,q,(id,id)) : _rows(s.db,q,(id,))
end

"""
    bfs_certificate_path(s::MoveGraphStoreDB, source, target; max_depth=typemax(Int)) -> Union{Nothing,Vector{String}}

Return the state identifiers along a shortest directed path in the stored move graph.

# Arguments
Requires M15 ([conditions](../specs/move-graph-store.md#numbered-conditions)).

`s` is open with the common-store schema. `source` and `target` convert to
identifier strings. Neither identifier is checked for membership in `states`.
Use a nonnegative integer `max_depth` to bound the number of edges; this
precondition is the caller's responsibility. Connection and conversion
errors propagate. 

# Returns
A `Vector{String}` containing both endpoints, or `nothing` when no path
is found within the bound. Equal identifiers return a one-element vector
even if the state is absent from the database.

# Guarantees
1. For distinct endpoints and a nonnegative integer bound, the returned
   sequence follows stored directed edges and has minimum edge count among
   such paths within the bound. [proved/agent: by construction]
2. The store and arguments are not modified. [proved/agent: by construction]

# Not guaranteed
The result does not replay or verify edges. A mathematical equivalence
follows only if stored transitions are independently valid. `nothing`
does not establish inequivalence. Reverse traversal requires stored reverse
edges; certificate contents and edge identifiers are not returned.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql transitions — the directed graph whose stored edges are searched.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> record_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true);

julia> bfs_certificate_path(s, a.state_id, b.state_id) == [a.state_id, b.state_id]
true

julia> close_store!(s);
```

Role: access
""" 
function bfs_certificate_path(s::MoveGraphStoreDB, source, target; max_depth=typemax(Int))
    src, dst = String(source), String(target)
    src == dst && return [src]
    queue = [src]; parent = Dict{String,Union{Nothing,String}}(src => nothing); depth = Dict(src => 0)
    head = 1
    while head <= length(queue)
        current = queue[head]; head += 1
        depth[current] >= max_depth && continue
        for row in neighbors(s, current)
            nxt = String(row.target_state_id)
            haskey(parent,nxt) && continue
            parent[nxt] = current; depth[nxt] = depth[current] + 1; push!(queue,nxt)
            if nxt == dst
                head = length(queue) + 1
                break
            end
        end
    end
    haskey(parent,dst) || return nothing
    path = String[]; current = dst
    while current !== nothing
        push!(path,current); current = parent[current]
    end
    reverse!(path)
end
"""
    certificate_path(s::MoveGraphStoreDB, source, target; max_depth=typemax(Int)) -> Union{Nothing,Vector{String}}

Return the state identifiers along a shortest directed path in the stored move graph.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema. `source` and `target` convert to
identifier strings. Neither identifier is checked for membership in `states`.
Use a nonnegative integer `max_depth` to bound the number of edges; this
precondition is the caller's responsibility. Connection and conversion
errors propagate. Calls are those of [`bfs_certificate_path`](@ref).

# Returns
A `Vector{String}` containing both endpoints, or `nothing` when no path
is found within the bound. Equal identifiers return a one-element vector
even if the state is absent from the database.

# Guarantees
1. For distinct endpoints and a nonnegative integer bound, the returned
   sequence follows stored directed edges and has minimum edge count among
   such paths within the bound. [proved/agent: by construction]
2. The store and arguments are not modified. [proved/agent: by construction]

# Not guaranteed
The result does not replay or verify edges. A mathematical equivalence
follows only if stored transitions are independently valid. `nothing`
does not establish inequivalence. Reverse traversal requires stored reverse
edges; certificate contents and edge identifiers are not returned.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql transitions — the directed graph whose stored edges are searched.

# Examples
```jldoctest
julia> s = open_store(); a = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss"); b = CanonicalState("state-b"; area=Topology.AREA_2, representation="KnotGauss");

julia> record_transition!(s, a, b, "R3"; verifier=(args...)->true, verified=true);

julia> certificate_path(s, a.state_id, b.state_id) == [a.state_id, b.state_id]
true

julia> close_store!(s);
```

Role: access
"""
certificate_path(args...; kwargs...) = bfs_certificate_path(args...; kwargs...)

"""
    import_knotinfo!(s::MoveGraphStoreDB, json_path::AbstractString) -> Int

Import the first sheet of a KnotInfo snapshot as reference records and named knot metadata.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is an open, writable store. `json_path` must name a readable snapshot
in the repository's KnotInfo JSON format, with source metadata, a first
sheet, headers and indexed rows. The caller supplies that format; missing
keys, malformed JSON, field conversion, I/O and database errors propagate.
The supplied snapshot digest is trusted rather than recomputed.

# Returns
The number of rows in the snapshot's first sheet as `Int`, including rows
that replace previously imported records.

# Guarantees
1. On success, first-sheet records are stored under the supplied snapshot
   digest and row indices; named knot metadata are updated from those rows.
   [proved/agent: by construction]
2. No diagram state is identified from a knot name or created for those
   reference records. [proved/agent: by construction]
3. The database is modified; the input file is only read.
   [proved/agent: by construction]

# Not guaranteed
No snapshot authenticity, digest correctness, knot invariant or mathematical
classification is verified. Additional sheets are ignored. The import is not
atomic; an error can leave a partial import. Reimporting an existing digest
does not refresh the original stored raw snapshot metadata.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql knotinfo_snapshots, knotinfo_raw_records and knot_types — reference data separate from diagram states.

# Examples
```jldoctest
julia> path = normpath(joinpath(dirname(pathof(Topology)), "..", "..", "..", "threads", "PJ260830", "data", "knotinfo-n6.json")); s = open_store();

julia> Topology.MoveGraphStore.import_knotinfo!(s, path)
8

julia> close_store!(s);
```

Role: other
"""
function import_knotinfo!(s::MoveGraphStoreDB, json_path::AbstractString)
    raw = read(json_path, String); data = JSON.parse(raw)
    sha = String(data["source"]["snapshot_sha256"])
    sheet = first(data["sheets"]); headers = sheet["headers"]
    _dbexec(s.db, "INSERT OR IGNORE INTO knotinfo_snapshots(snapshot_sha256,source_url,schema_version,raw_json) VALUES (?,?,?,?)",
            (sha,get(data["source"],"url",nothing),get(data,"schema_version",nothing),raw))
    snapshot = _one(s.db,"SELECT id FROM knotinfo_snapshots WHERE snapshot_sha256=?",(sha,)); snapid=Int(snapshot.id)
    for rec in sheet["rows"]
        vals=rec["values"]; name=String(vals[1]); rowidx=Int(rec["source_row_index"])
        _dbexec(s.db,"INSERT OR REPLACE INTO knotinfo_raw_records(snapshot_id,source_row_index,crossing_number,name,headers_json,values_json) VALUES (?,?,?,?,?,?)",
                (snapid,rowidx,Int(rec["crossing_number"]),name,JSON.json(headers),JSON.json(vals)))
        invariants=Dict{String,Any}()
        for (h,v) in zip(headers,vals)
            invariants[String(h)] = v
        end
        _dbexec(s.db,"INSERT OR REPLACE INTO knot_types(name,area,knotinfo_snapshot_id,knotinfo_row_index,invariants_json) VALUES (?,?,?,?,?)",
                (name,AREA_2,snapid,rowidx,JSON.json(invariants)))
    end
    length(sheet["rows"])
end

"""
    knotinfo_records(s::MoveGraphStoreDB; snapshot_sha256=nothing) -> Vector

Read imported KnotInfo reference rows from the store.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/move-graph-store.md).

`s` is open with the common-store schema. Omitted `snapshot_sha256` selects
all snapshots; otherwise it converts to `String` and selects that exact
stored digest. Database and conversion errors propagate. No digest or
record-content validation is performed.

# Returns
A vector of row `NamedTuple`s ordered by source row index, with record and
snapshot identifiers, source row index, crossing number, knot name, and
JSON text for headers and values. An absent snapshot gives an empty vector.

# Guarantees
1. The result contains the imported rows matching the requested snapshot
   filter. [proved/agent: by construction]
2. The store is not modified; JSON columns remain encoded text.
   [proved/agent: by construction]

# Not guaranteed
Rows are reference metadata, not verified diagram identifications. Ordering
of equal row indices across different snapshots is unspecified.

# Sources
cf:lib/topology.jl/db/knot_move_schema.sql knotinfo_raw_records — imported snapshot rows.

# Examples
```jldoctest
julia> path = normpath(joinpath(dirname(pathof(Topology)), "..", "..", "..", "threads", "PJ260830", "data", "knotinfo-n6.json")); s = open_store(); Topology.MoveGraphStore.import_knotinfo!(s, path);

julia> length(Topology.MoveGraphStore.knotinfo_records(s))
8

julia> close_store!(s);
```

Role: access
"""
function knotinfo_records(s::MoveGraphStoreDB; snapshot_sha256=nothing)
    q = snapshot_sha256 === nothing ? "SELECT * FROM knotinfo_raw_records ORDER BY source_row_index" :
        "SELECT r.* FROM knotinfo_raw_records r JOIN knotinfo_snapshots s ON s.id=r.snapshot_id WHERE s.snapshot_sha256=? ORDER BY r.source_row_index"
    snapshot_sha256 === nothing ? _rows(s.db,q) : _rows(s.db,q,(String(snapshot_sha256),))
end

end
