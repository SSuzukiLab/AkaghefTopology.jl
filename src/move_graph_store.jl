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

state_from_serialization(serialization; kwargs...) = CanonicalState(serialization; kwargs...)
state_id(s::CanonicalState) = s.state_id
state_hash(s::CanonicalState) = s.state_hash

struct MoveType
    move_type::String
    area::Float64
    inverse_type::String
    inverse_direction::Symbol
end

struct MoveGraphStoreDB
    db::SQLite.DB
end
const Store = MoveGraphStoreDB

"""Dataset-level metadata for one reproducible classification method."""
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

function open_store(path::AbstractString=":memory:"; migrate=true)
    s = MoveGraphStoreDB(SQLite.DB(String(path)))
    migrate && migrate!(s)
    s
end
MoveGraphStoreDB(path::AbstractString) = open_store(path)
close_store!(s::MoveGraphStoreDB) = SQLite.close(s.db)

function _payload(x)
    x === nothing && return nothing
    x isa AbstractVector{UInt8} && return x
    x isa AbstractString && return String(x)
    JSON.json(x)
end

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

"""Read one dataset-level classification run by its method version."""
function get_classification_run(s::MoveGraphStoreDB, method_version::AbstractString)
    row = _one(s.db, "SELECT * FROM classification_runs WHERE method_version=?",
               (String(method_version),))
    _classification_run_from_row(row)
end

classification_run(args...; kwargs...) = get_classification_run(args...; kwargs...)

"""Insert one immutable metadata row for a reproducible classification run."""
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

"""Persist one state classification with Julia-side and SQL-side checks."""
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

put_state_classification!(s::MoveGraphStoreDB, state, method_version::AbstractString;
                          kwargs...) =
    put_state_classification!(s, state; method_version=method_version, kwargs...)
record_state_classification!(args...; kwargs...) =
    put_state_classification!(args...; kwargs...)

"""Read one persisted state classification."""
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

state_classification(args...; kwargs...) =
    get_state_classification(args...; kwargs...)

"""Read all classifications for one method version in state-id order."""
function classification_rows(s::MoveGraphStoreDB, method_version::AbstractString)
    rows = _rows(s.db, """
        SELECT * FROM state_classifications
        WHERE method_version=? ORDER BY state_id
    """, (String(method_version),))
    _classification_row.(rows)
end

"""Return SQL-aggregated status and knot-type totals for one run."""
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

"Atomically persist states and a verifier-authorized, verified edge." 
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
record_transition!(args...; kwargs...) = insert_verified_transition!(args...; kwargs...)

function add_certificate!(s::MoveGraphStoreDB, transition_id::Integer, path)
    _dbexec(s.db, "INSERT OR IGNORE INTO certificate_paths(transition_id,path) VALUES (?,?)",
            (Int(transition_id), _json_or_nothing(path)))
    Int(transition_id)
end

function add_invariant!(s::MoveGraphStoreDB, state, name, value; source=nothing)
    sid = state isa CanonicalState ? state.state_id : String(state)
    _dbexec(s.db, "INSERT OR IGNORE INTO invariants(state_id,name,value,source) VALUES (?,?,?,?)",
            (sid,String(name),_json_or_nothing(value),source === nothing ? nothing : String(source)))
end

function neighbors(s::MoveGraphStoreDB, sid; direction=:outgoing)
    id = String(sid)
    q = direction == :incoming ? "SELECT t.* FROM transitions t WHERE target_state_id=? ORDER BY t.id" :
        direction == :both ? "SELECT t.* FROM transitions t WHERE source_state_id=? OR target_state_id=? ORDER BY t.id" :
        "SELECT t.* FROM transitions t WHERE source_state_id=? ORDER BY t.id"
    direction == :both ? _rows(s.db,q,(id,id)) : _rows(s.db,q,(id,))
end

"Return the state IDs on a shortest verified path, including both endpoints." 
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
certificate_path(args...; kwargs...) = bfs_certificate_path(args...; kwargs...)

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

function knotinfo_records(s::MoveGraphStoreDB; snapshot_sha256=nothing)
    q = snapshot_sha256 === nothing ? "SELECT * FROM knotinfo_raw_records ORDER BY source_row_index" :
        "SELECT r.* FROM knotinfo_raw_records r JOIN knotinfo_snapshots s ON s.id=r.snapshot_id WHERE s.snapshot_sha256=? ORDER BY r.source_row_index"
    snapshot_sha256 === nothing ? _rows(s.db,q) : _rows(s.db,q,(String(snapshot_sha256),))
end

end
