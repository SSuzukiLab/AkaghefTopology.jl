"""Area-3 persistence adapter for closed normal o-graph presentations.

The generic move graph remains schema v2.  This adapter adds a separately
versioned extension for closedness evidence, source provenance, relation
profiles, finite campaign bounds, and ordered relation certificates.
"""

using JSON
using SHA
using .MoveGraphStore

const OGRAPH_STORE_SCHEMA_VERSION = 2
const OGRAPH_IDENTITY_DOMAIN = "area3/closed-normal-ograph/presentation-v1"

function _ograph_schema_path()
    normpath(joinpath(@__DIR__, "..", "db", "ograph_relation_schema.sql"))
end

function _ograph_canonical_json(value)
    value === nothing && return "null"
    value === missing && throw(ArgumentError("missing is not canonical JSON"))
    value isa Bool && return value ? "true" : "false"
    if value isa Integer
        typemin(Int64) <= big(value) <= typemax(Int64) || throw(ArgumentError(
            "canonical JSON integers must fit signed 64-bit replay"))
        return string(Int64(value))
    end
    if value isa AbstractFloat
        isfinite(value) || throw(ArgumentError("non-finite numbers are not canonical JSON"))
        return JSON.json(value)
    end
    value isa Number && throw(ArgumentError(
        "unsupported canonical JSON number $(typeof(value))"))
    value isa AbstractString && return JSON.json(String(value))
    if value isa NamedTuple
        return _ograph_canonical_json(Dict(String(key) => getfield(value, key)
                                           for key in keys(value)))
    end
    if value isa AbstractDict
        entries = [(key isa AbstractString || key isa Symbol) ?
                   (String(key), item) : throw(ArgumentError(
                       "canonical JSON object keys must be strings or symbols"))
                   for (key, item) in value]
        length(unique(first.(entries))) == length(entries) || throw(ArgumentError(
            "canonical JSON object keys collide after string conversion"))
        sort!(entries; by=first)
        return "{" * join((JSON.json(key) * ":" * _ograph_canonical_json(item)
                            for (key, item) in entries), ",") * "}"
    end
    if value isa Tuple || value isa AbstractVector
        return "[" * join((_ograph_canonical_json(item) for item in value), ",") * "]"
    end
    throw(ArgumentError("unsupported canonical JSON value $(typeof(value))"))
end

ograph_canonical_json(value) = _ograph_canonical_json(value)

_ograph_sha256(text::AbstractString) =
    bytes2hex(SHA.sha256(codeunits(String(text))))

function _ograph_existing_schema_version(db)
    meta = MoveGraphStore._one(db, """
        SELECT 1 AS present FROM sqlite_master
        WHERE type='table' AND name='ograph_schema_meta'
    """)
    if meta === nothing
        partial = MoveGraphStore._rows(db, """
            SELECT name FROM sqlite_master
            WHERE type='table' AND name LIKE 'ograph_%'
        """)
        isempty(partial) || throw(ArgumentError(
            "o-graph extension tables exist without ograph_schema_meta"))
        return nothing
    end
    row = MoveGraphStore._one(db,
        "SELECT value FROM ograph_schema_meta WHERE key='schema_version'")
    row === nothing && throw(ArgumentError(
        "ograph_schema_meta has no schema_version"))
    try
        parse(Int, String(row.value))
    catch
        throw(ArgumentError("ograph_schema_meta.schema_version is not an integer"))
    end
end

function _assert_ograph_schema(store::MoveGraphStoreDB)
    _ograph_existing_schema_version(store.db) == OGRAPH_STORE_SCHEMA_VERSION ||
        throw(ArgumentError(
            "o-graph schema migration did not reach version $(OGRAPH_STORE_SCHEMA_VERSION)"))
    columns = Set(String(row.name) for row in
                  MoveGraphStore._rows(store.db, "PRAGMA table_info(transitions)"))
    "relation_profile_id" in columns || throw(ArgumentError(
        "o-graph store is missing transitions.relation_profile_id"))
    isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check")) ||
        throw(ArgumentError("o-graph store foreign-key check failed"))
    true
end

"""Install the o-graph extension after the generic v2 store is valid."""
function migrate_ograph_store!(store::MoveGraphStoreDB)
    MoveGraphStore.migrate!(store)
    version = _ograph_existing_schema_version(store.db)
    version === nothing || version <= OGRAPH_STORE_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph schema version $version"))
    sql = read(_ograph_schema_path(), String)
    MoveGraphStore._transaction(store.db) do
        columns = Set(String(row.name) for row in
                      MoveGraphStore._rows(store.db, "PRAGMA table_info(transitions)"))
        "relation_profile_id" in columns || MoveGraphStore._dbexec(store.db,
            "ALTER TABLE transitions ADD COLUMN relation_profile_id TEXT")
        if version == 1
            # Append-only guards are a contract change.  Replay the complete
            # legacy extension before freezing its rows, and roll the version
            # update back if any evidence cannot be reconstructed exactly.
            MoveGraphStore._dbexec(store.db, """
                UPDATE ograph_schema_meta SET value='2' WHERE key='schema_version'
            """)
            verify_ograph_store(store)
        end
        MoveGraphStore._execute_schema!(store.db, sql)
        _assert_ograph_schema(store)
    end
    store
end

struct OGraphProvenance
    source_id::String
    source_kind::String
    locator::String
    snapshot_sha256::Union{Nothing,String}
    metadata_json::String
end

function OGraphProvenance(source_kind::AbstractString, locator::AbstractString;
                          snapshot_sha256=nothing, metadata=Dict{String,Any}())
    kind = strip(String(source_kind))
    location = strip(String(locator))
    isempty(kind) && throw(ArgumentError("source_kind must not be empty"))
    isempty(location) && throw(ArgumentError("source locator must not be empty"))
    snapshot = snapshot_sha256 === nothing ? nothing : lowercase(strip(String(snapshot_sha256)))
    snapshot === nothing || occursin(r"^[0-9a-f]{64}$", snapshot) ||
        throw(ArgumentError("snapshot_sha256 must be a 64-character hexadecimal digest"))
    metadata_json = _ograph_canonical_json(metadata)
    identity = _ograph_canonical_json(Dict(
        "source_kind" => kind,
        "locator" => location,
        "snapshot_sha256" => snapshot,
        "metadata" => JSON.parse(metadata_json),
    ))
    OGraphProvenance("ograph-source-" * _ograph_sha256(identity), kind,
                     location, snapshot, metadata_json)
end

function _closedness_dict(witness::ClosednessWitness)
    Dict(
        "canonical_hash" => witness.canonical_hash,
        "verifier_version" => witness.verifier_version,
        "criteria" => collect(witness.criteria),
        "c1_oriented_circuits" => witness.c1_oriented_circuits,
        "c2_trivalent_components" => witness.c2_trivalent_components,
        "c3_method_circuits" => witness.c3_method_circuits,
        "true_vertices" => witness.true_vertices,
        "method_circuit_lengths" => collect(witness.method_circuit_lengths),
    )
end

function _ograph_state_payload(graph::ClosedNormalOGraph)
    _ograph_canonical_json(Dict(
        "identity_domain" => OGRAPH_IDENTITY_DOMAIN,
        "schema" => OGRAPH_PRESENTATION_SCHEMA,
        "canonicalization" => OGRAPH_PRESENTATION_CANONICALIZATION,
        "gauss" => collect(graph.state.gauss),
        "vertex_signs" => collect(graph.state.vertex_signs),
        "convention_id" => graph.state.convention_id,
        "closedness" => _closedness_dict(graph.witness),
    ))
end

"""Wrap one certified presentation in the generic domain-qualified state."""
function ograph_canonical_state(graph::ClosedNormalOGraph)
    CanonicalState(ograph_canonical_serialization(graph.state);
        area=AREA_3, representation="o-graph",
        payload=_ograph_state_payload(graph),
        canonical_version=OGRAPH_PRESENTATION_CANON_VERSION,
        complexity=Float64(true_vertex_count(graph.state)))
end

_ograph_nullable_string(value) =
    value === nothing || value === missing ? nothing : String(value)

function _ograph_payload_text(value)
    value === nothing && return nothing
    value === missing && return nothing
    value isa AbstractVector{UInt8} && return String(value)
    String(value)
end

function _same_canonical_state(actual::CanonicalState, expected::CanonicalState)
    actual.state_id == expected.state_id &&
    actual.state_hash == expected.state_hash &&
    actual.canonical_version == expected.canonical_version &&
    actual.serialization == expected.serialization &&
    _ograph_payload_text(actual.payload) == _ograph_payload_text(expected.payload) &&
    actual.area == expected.area &&
    actual.representation == expected.representation &&
    actual.complexity == expected.complexity
end

function _assert_source_row(row, source::OGraphProvenance)
    metadata = try
        JSON.parse(source.metadata_json)
    catch error_value
        throw(ArgumentError("o-graph source metadata is not valid JSON: $(error_value)"))
    end
    canonical = OGraphProvenance(source.source_kind, source.locator;
        snapshot_sha256=source.snapshot_sha256, metadata=metadata)
    provided_identity = (source.source_id, source.source_kind, source.locator,
                         source.snapshot_sha256, source.metadata_json)
    derived_identity = (canonical.source_id, canonical.source_kind,
                        canonical.locator, canonical.snapshot_sha256,
                        canonical.metadata_json)
    provided_identity == derived_identity || throw(ArgumentError(
        "o-graph source identity is not derived from its canonical evidence"))
    row === nothing && throw(ArgumentError("failed to persist o-graph source"))
    actual = (String(row.source_id), String(row.source_kind), String(row.locator),
              _ograph_nullable_string(row.snapshot_sha256),
              String(row.metadata_json))
    expected = (source.source_id, source.source_kind, source.locator,
                source.snapshot_sha256, source.metadata_json)
    actual == expected || throw(ArgumentError(
        "o-graph source id already exists with different evidence"))
    true
end

"""Persist a certified state and any number of provenance aliases atomically."""
function put_closed_ograph!(store::MoveGraphStoreDB, graph::ClosedNormalOGraph;
                            provenance::AbstractVector{<:OGraphProvenance}=OGraphProvenance[])
    _assert_ograph_schema(store)
    state = ograph_canonical_state(graph)
    evidence = _ograph_canonical_json(_closedness_dict(graph.witness))
    MoveGraphStore._transaction(store.db) do
        persisted = get_state(store, state.state_id)
        if persisted === nothing
            persisted = put_state!(store, state)
        else
            _same_canonical_state(persisted, state) || throw(ArgumentError(
                "canonical o-graph state already exists with different base data"))
        end
        MoveGraphStore._dbexec(store.db, """
            INSERT OR IGNORE INTO ograph_states(
              state_id,identity_domain,presentation_hash,convention_id,
              closedness_verifier_version,closedness_evidence)
            VALUES (?,?,?,?,?,?)
        """, (state.state_id, OGRAPH_IDENTITY_DOMAIN,
                ograph_canonical_hash(graph.state), graph.state.convention_id,
                graph.witness.verifier_version, evidence))
        row = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_states WHERE state_id=?", (state.state_id,))
        expected = (OGRAPH_IDENTITY_DOMAIN, ograph_canonical_hash(graph.state),
                    graph.state.convention_id, graph.witness.verifier_version, evidence)
        actual = (String(row.identity_domain), String(row.presentation_hash),
                  String(row.convention_id), String(row.closedness_verifier_version),
                  String(row.closedness_evidence))
        actual == expected || throw(ArgumentError(
            "o-graph state already exists with different evidence"))
        for source in provenance
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_sources(
                  source_id,source_kind,locator,snapshot_sha256,metadata_json)
                VALUES (?,?,?,?,?)
            """, (source.source_id, source.source_kind, source.locator,
                    source.snapshot_sha256, source.metadata_json))
            stored = MoveGraphStore._one(store.db,
                "SELECT * FROM ograph_sources WHERE source_id=?", (source.source_id,))
            _assert_source_row(stored, source)
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_state_sources(state_id,source_id)
                VALUES (?,?)
            """, (state.state_id, source.source_id))
        end
        state
    end
end

struct OGraphRelationProfile
    profile_id::String
    profile_name::String
    profile_version::Int
    parent_profile_id::Union{Nothing,String}
    moves::Tuple{Vararg{String}}
    moves_complete::Bool
    evidence_status::String
    conventions_json::String
    profile_fingerprint::String
    created_at_utc::String
end

function _expected_profile_moves(name::String)
    name == "presentation" && return ()
    name == "combed" && return ("MP", "PS")
    name == "manifold" && return ("CP", "MP", "PS")
    throw(ArgumentError(
        "weighted framed/spin profiles require a separate typed state adapter"))
end

function OGraphRelationProfile(profile_name::AbstractString;
        profile_version::Integer=1, parent_profile_id=nothing,
        moves_complete::Bool=false, evidence_status::AbstractString="unimplemented",
        conventions=Dict{String,Any}(), created_at_utc::AbstractString)
    name = lowercase(strip(String(profile_name)))
    version = Int(profile_version)
    version > 0 || throw(ArgumentError("profile_version must be positive"))
    moves = _expected_profile_moves(name)
    parent = parent_profile_id === nothing ? nothing : String(parent_profile_id)
    status = String(evidence_status)
    if name == "presentation"
        moves_complete || throw(ArgumentError("presentation profile is identity-only and complete"))
        status == "identity-only" || throw(ArgumentError(
            "presentation profile evidence_status must be identity-only"))
        parent === nothing || throw(ArgumentError("presentation profile cannot have a parent"))
    else
        parent === nothing && throw(ArgumentError("non-presentation profile requires a parent"))
        status in ("unimplemented", "verified") || throw(ArgumentError(
            "relation profile evidence_status must be unimplemented or verified"))
        moves_complete == (status == "verified") || throw(ArgumentError(
            "verified relation profiles must be complete and incomplete profiles unimplemented"))
    end
    conventions_json = _ograph_canonical_json(conventions)
    fingerprint_input = _ograph_canonical_json(Dict(
        "profile_name" => name,
        "profile_version" => version,
        "parent_profile_id" => parent,
        "moves" => collect(moves),
        "moves_complete" => moves_complete,
        "evidence_status" => status,
        "conventions" => JSON.parse(conventions_json),
    ))
    fingerprint = _ograph_sha256(fingerprint_input)
    OGraphRelationProfile("ograph-profile-" * fingerprint, name, version, parent,
        moves, moves_complete, status, conventions_json, fingerprint,
        String(created_at_utc))
end

function _ograph_profile_moves(store::MoveGraphStoreDB, profile_id::AbstractString)
    Tuple(String(row.move_type) for row in MoveGraphStore._rows(store.db, """
        SELECT move_type FROM ograph_relation_profile_moves
        WHERE profile_id=? ORDER BY move_type
    """, (String(profile_id),)))
end


function _assert_profile_row(store::MoveGraphStoreDB, row,
                             profile::OGraphRelationProfile)
    conventions = try
        JSON.parse(profile.conventions_json)
    catch error_value
        throw(ArgumentError(
            "relation profile conventions are not valid JSON: $(error_value)"))
    end
    canonical = OGraphRelationProfile(profile.profile_name;
        profile_version=profile.profile_version,
        parent_profile_id=profile.parent_profile_id,
        moves_complete=profile.moves_complete,
        evidence_status=profile.evidence_status,
        conventions=conventions,
        created_at_utc=profile.created_at_utc)
    provided_identity = (profile.profile_id, profile.profile_name,
        profile.profile_version, profile.parent_profile_id, profile.moves,
        profile.moves_complete, profile.evidence_status, profile.conventions_json,
        profile.profile_fingerprint, profile.created_at_utc)
    derived_identity = (canonical.profile_id, canonical.profile_name,
        canonical.profile_version, canonical.parent_profile_id, canonical.moves,
        canonical.moves_complete, canonical.evidence_status,
        canonical.conventions_json, canonical.profile_fingerprint,
        canonical.created_at_utc)
    provided_identity == derived_identity || throw(ArgumentError(
        "relation profile identity is not derived from its canonical semantics"))
    row === nothing && throw(ArgumentError("relation profile is not registered"))
    parent = _ograph_nullable_string(row.parent_profile_id)
    actual = (String(row.profile_id), String(row.profile_name),
              Int(row.profile_version), parent,
              String(row.profile_fingerprint), Int(row.moves_complete),
              String(row.evidence_status), String(row.conventions_json),
              String(row.created_at_utc))
    expected = (profile.profile_id, profile.profile_name, profile.profile_version,
                profile.parent_profile_id, profile.profile_fingerprint,
                profile.moves_complete ? 1 : 0, profile.evidence_status,
                profile.conventions_json, profile.created_at_utc)
    actual == expected || throw(ArgumentError(
        "relation profile id already exists with different metadata"))
    _ograph_profile_moves(store, profile.profile_id) == profile.moves ||
        throw(ArgumentError("relation profile move set does not match its fingerprint"))
    if profile.profile_name != "presentation"
        parent = MoveGraphStore._one(store.db, """
            SELECT profile_name FROM ograph_relation_profiles WHERE profile_id=?
        """, (profile.parent_profile_id,))
        expected_parent = profile.profile_name == "combed" ? "presentation" : "combed"
        parent !== nothing && String(parent.profile_name) == expected_parent ||
            throw(ArgumentError(
                "$(profile.profile_name) profile requires a $expected_parent parent"))
    end
    for move_type in profile.moves
        move = get_move_type(store, move_type)
        move !== nothing && move.area == AREA_3 &&
            move.inverse_type == move_type && move.inverse_direction == :inverse ||
            throw(ArgumentError("move type $move_type has incompatible stored semantics"))
    end
    true
end

function register_ograph_profile!(store::MoveGraphStoreDB,
                                  profile::OGraphRelationProfile)
    _assert_ograph_schema(store)
    MoveGraphStore._transaction(store.db) do
        row = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_relation_profiles WHERE profile_id=?",
            (profile.profile_id,))
        if row !== nothing
            _assert_profile_row(store, row, profile)
            return profile
        end
        for move_type in profile.moves
            existing_move = get_move_type(store, move_type)
            if existing_move === nothing
                register_move_type!(store, move_type; area=AREA_3,
                                    inverse_type=move_type,
                                    inverse_direction=:inverse)
            else
                existing_move.area == AREA_3 &&
                    existing_move.inverse_type == move_type &&
                    existing_move.inverse_direction == :inverse ||
                    throw(ArgumentError("move type $move_type has incompatible stored semantics"))
            end
        end
        MoveGraphStore._dbexec(store.db, """
            INSERT INTO ograph_relation_profiles(
              profile_id,profile_name,profile_version,parent_profile_id,
              profile_fingerprint,moves_complete,evidence_status,
              conventions_json,created_at_utc)
            VALUES (?,?,?,?,?,?,?,?,?)
        """, (profile.profile_id, profile.profile_name, profile.profile_version,
                profile.parent_profile_id, profile.profile_fingerprint,
                profile.moves_complete ? 1 : 0, profile.evidence_status,
                profile.conventions_json, profile.created_at_utc))
        for move_type in profile.moves
            MoveGraphStore._dbexec(store.db, """
                INSERT INTO ograph_relation_profile_moves(profile_id,move_type)
                VALUES (?,?)
            """, (profile.profile_id, move_type))
        end
        row = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_relation_profiles WHERE profile_id=?",
            (profile.profile_id,))
        _assert_profile_row(store, row, profile)
        profile
    end
end

function _assert_campaign_row(store::MoveGraphStoreDB, row, expected, ids)
    row === nothing && throw(ArgumentError("campaign is not registered"))
    actual = (campaign_id=String(row.campaign_id),
              profile_id=String(row.profile_id),
              n_input=Int(row.n_input), b_path=Int(row.b_path),
              input_count=Int(row.input_count), caps_json=String(row.caps_json),
              input_fingerprint=String(row.input_fingerprint),
              campaign_fingerprint=String(row.campaign_fingerprint),
              engine_version=String(row.engine_version),
              verifier_version=String(row.verifier_version),
              status=String(row.status), stop_reason=String(row.stop_reason),
              counts_json=String(row.counts_json),
              created_at_utc=String(row.created_at_utc))
    actual == expected || throw(ArgumentError(
        "campaign id already exists with different metadata"))
    inputs = MoveGraphStore._rows(store.db, """
        SELECT ordinal,state_id FROM ograph_relation_campaign_inputs
        WHERE campaign_id=? ORDER BY ordinal
    """, (expected.campaign_id,))
    actual_inputs = [(Int(row.ordinal), String(row.state_id)) for row in inputs]
    expected_inputs = [(ordinal - 1, state_id) for (ordinal, state_id) in enumerate(ids)]
    actual_inputs == expected_inputs || throw(ArgumentError(
        "campaign input rows do not exactly match their fingerprint"))
    true
end

"""Register a bounded campaign.  N_input and B_path are true-vertex bounds."""
function register_ograph_campaign!(store::MoveGraphStoreDB, campaign_id::AbstractString,
        profile::OGraphRelationProfile, input_states::AbstractVector{<:CanonicalState};
        n_input::Integer, b_path::Integer, caps,
        engine_version::AbstractString, verifier_version::AbstractString,
        status::AbstractString, stop_reason::AbstractString, counts,
        created_at_utc::AbstractString)
    _assert_ograph_schema(store)
    n_bound, path_bound = Int(n_input), Int(b_path)
    n_bound >= 0 || throw(ArgumentError("N_input must be nonnegative"))
    path_bound >= n_bound || throw(ArgumentError("B_path must be at least N_input"))
    ids = sort(unique(String(state.state_id) for state in input_states))
    length(ids) == length(input_states) || throw(ArgumentError(
        "campaign inputs must have distinct canonical state IDs"))
    for state in input_states
        state.area == AREA_3 || throw(ArgumentError("campaign inputs must be Area 3 states"))
        state.complexity === nothing && throw(ArgumentError(
            "campaign input is missing true-vertex complexity"))
        isfinite(state.complexity) && isinteger(state.complexity) || throw(ArgumentError(
            "campaign input true-vertex complexity must be a finite integer"))
        state.complexity <= n_bound || throw(ArgumentError(
            "campaign input exceeds N_input"))
    end
    caps_json = _ograph_canonical_json(caps)
    counts_json = _ograph_canonical_json(counts)
    input_fingerprint = _ograph_sha256(_ograph_canonical_json(ids))
    identifier = strip(String(campaign_id))
    isempty(identifier) && throw(ArgumentError("campaign_id must not be empty"))
    campaign_data = Dict(
        "campaign_id" => identifier,
        "profile_id" => profile.profile_id,
        "N_input" => n_bound,
        "B_path" => path_bound,
        "input_state_ids" => ids,
        "caps" => JSON.parse(caps_json),
        "input_fingerprint" => input_fingerprint,
        "engine_version" => String(engine_version),
        "verifier_version" => String(verifier_version),
        "status" => String(status),
        "stop_reason" => String(stop_reason),
        "counts" => JSON.parse(counts_json),
    )
    campaign_fingerprint = _ograph_sha256(_ograph_canonical_json(campaign_data))
    expected = (campaign_id=identifier, profile_id=profile.profile_id,
                n_input=n_bound, b_path=path_bound, input_count=length(ids),
                caps_json=caps_json, input_fingerprint=input_fingerprint,
                campaign_fingerprint=campaign_fingerprint,
                engine_version=String(engine_version),
                verifier_version=String(verifier_version), status=String(status),
                stop_reason=String(stop_reason), counts_json=counts_json,
                created_at_utc=String(created_at_utc))
    MoveGraphStore._transaction(store.db) do
        profile_row = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_relation_profiles WHERE profile_id=?",
            (profile.profile_id,))
        _assert_profile_row(store, profile_row, profile)
        existing = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_relation_campaigns WHERE campaign_id=?", (identifier,))
        if existing !== nothing
            _assert_campaign_row(store, existing, expected, ids)
            return (campaign_id=identifier, campaign_fingerprint=campaign_fingerprint,
                    input_fingerprint=input_fingerprint)
        end
        for state in input_states
            stored_state = get_state(store, state.state_id)
            stored_state !== nothing && _same_canonical_state(stored_state, state) ||
                throw(ArgumentError("campaign input is not the exact persisted base state"))
            MoveGraphStore._one(store.db,
                "SELECT state_id FROM ograph_states WHERE state_id=?", (state.state_id,)) !== nothing ||
                throw(ArgumentError("campaign input is not a registered o-graph state"))
        end
        MoveGraphStore._dbexec(store.db, """
            INSERT INTO ograph_relation_campaigns(
              campaign_id,profile_id,n_input,b_path,input_count,caps_json,
              input_fingerprint,campaign_fingerprint,engine_version,
              verifier_version,status,stop_reason,counts_json,created_at_utc)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (identifier, profile.profile_id, n_bound, path_bound, length(ids),
                caps_json, input_fingerprint, campaign_fingerprint,
                String(engine_version), String(verifier_version), String(status),
                String(stop_reason), counts_json, String(created_at_utc)))
        stored = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_relation_campaigns WHERE campaign_id=?", (identifier,))
        for (ordinal, state_id_value) in enumerate(ids)
            MoveGraphStore._dbexec(store.db, """
                INSERT INTO ograph_relation_campaign_inputs(
                  campaign_id,ordinal,state_id) VALUES (?,?,?)
            """, (identifier, ordinal - 1, state_id_value))
        end
        _assert_campaign_row(store, stored, expected, ids)
        (campaign_id=identifier, campaign_fingerprint=campaign_fingerprint,
         input_fingerprint=input_fingerprint)
    end
end

function _decode_ograph_state_payload(raw)
    text = raw isa AbstractVector{UInt8} ? String(raw) : String(raw)
    payload = JSON.parse(text)
    payload["identity_domain"] == OGRAPH_IDENTITY_DOMAIN || throw(ArgumentError(
        "stored o-graph payload has the wrong identity domain"))
    OGraphState(Int.(payload["gauss"]), Int.(payload["vertex_signs"]);
                convention_id=String(payload["convention_id"]))
end

"""Recompute every v1 o-graph pilot identity and finite evidence record.

Schema v1 deliberately has no public Area-3 transition writer.  Until a typed,
source-bound MP/PS/CP verifier can replay each step (including Method-B evidence
for nonlocal PS), this verifier rejects every stored o-graph transition or
relation certificate instead of treating a generic callback as mathematical
evidence.
"""
function verify_ograph_store(store::MoveGraphStoreDB)
    _assert_ograph_schema(store)
    base_version = MoveGraphStore._existing_schema_version(store.db)
    base_version == 2 || throw(ArgumentError("o-graph extension requires base schema v2"))
    integrity = String(MoveGraphStore._one(store.db, "PRAGMA integrity_check").integrity_check)
    integrity == "ok" || throw(ArgumentError("SQLite integrity_check failed: $integrity"))
    isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check")) ||
        throw(ArgumentError("SQLite foreign_key_check failed"))

    state_rows = MoveGraphStore._rows(store.db, """
        SELECT base.*, ext.identity_domain, ext.presentation_hash,
               ext.convention_id, ext.closedness_verifier_version,
               ext.closedness_evidence
        FROM states base JOIN ograph_states ext ON ext.state_id=base.state_id
        ORDER BY base.state_id
    """)
    base_ograph_count = Int(MoveGraphStore._one(store.db, """
        SELECT COUNT(*) AS n FROM states
        WHERE area=3.0 AND representation='o-graph'
    """).n)
    base_ograph_count == length(state_rows) || throw(ArgumentError(
        "Area-3 o-graph base states and extension rows are not one-to-one"))
    for row in state_rows
        payload_text = _ograph_payload_text(row.payload)
        payload_text === nothing && throw(ArgumentError(
            "stored o-graph state is missing its replay payload"))
        graph = certify_closed_ograph(_decode_ograph_state_payload(payload_text))
        expected = ograph_canonical_state(graph)
        actual = CanonicalState(String(row.state_id), String(row.state_hash),
            Int(row.canonical_version), String(row.serialization), payload_text,
            Float64(row.area), String(row.representation),
            row.complexity === missing || row.complexity === nothing ? nothing :
                Float64(row.complexity))
        _same_canonical_state(actual, expected) || throw(ArgumentError(
            "stored o-graph base state failed exact replay"))
        extension_actual = (String(row.identity_domain), String(row.presentation_hash),
            String(row.convention_id), String(row.closedness_verifier_version),
            String(row.closedness_evidence))
        extension_expected = (OGRAPH_IDENTITY_DOMAIN,
            ograph_canonical_hash(graph.state), graph.state.convention_id,
            graph.witness.verifier_version,
            _ograph_canonical_json(_closedness_dict(graph.witness)))
        extension_actual == extension_expected || throw(ArgumentError(
            "stored o-graph extension evidence failed exact replay"))
    end

    source_rows = MoveGraphStore._rows(store.db,
        "SELECT * FROM ograph_sources ORDER BY source_id")
    for row in source_rows
        metadata_json = String(row.metadata_json)
        metadata = JSON.parse(metadata_json)
        _ograph_canonical_json(metadata) == metadata_json || throw(ArgumentError(
            "stored o-graph source metadata is not canonical JSON"))
        source = OGraphProvenance(String(row.source_kind), String(row.locator);
            snapshot_sha256=_ograph_nullable_string(row.snapshot_sha256),
            metadata=metadata)
        _assert_source_row(row, source)
    end
    orphan_sources = Int(MoveGraphStore._one(store.db, """
        SELECT COUNT(*) AS n FROM ograph_sources source
        WHERE NOT EXISTS (
          SELECT 1 FROM ograph_state_sources link WHERE link.source_id=source.source_id
        )
    """).n)
    orphan_sources == 0 || throw(ArgumentError(
        "o-graph source provenance contains unassociated rows"))

    profile_rows = MoveGraphStore._rows(store.db,
        "SELECT * FROM ograph_relation_profiles ORDER BY profile_id")
    profiles = Dict{String,OGraphRelationProfile}()
    for row in profile_rows
        conventions_json = String(row.conventions_json)
        conventions = JSON.parse(conventions_json)
        _ograph_canonical_json(conventions) == conventions_json || throw(ArgumentError(
            "relation profile conventions are not canonical JSON"))
        profile = OGraphRelationProfile(String(row.profile_name);
            profile_version=Int(row.profile_version),
            parent_profile_id=_ograph_nullable_string(row.parent_profile_id),
            moves_complete=Int(row.moves_complete) == 1,
            evidence_status=String(row.evidence_status),
            conventions=conventions, created_at_utc=String(row.created_at_utc))
        _assert_profile_row(store, row, profile)
        profiles[profile.profile_id] = profile
    end

    campaigns = MoveGraphStore._rows(store.db,
        "SELECT * FROM ograph_relation_campaigns ORDER BY campaign_id")
    for campaign in campaigns
        profile_id = String(campaign.profile_id)
        haskey(profiles, profile_id) || throw(ArgumentError(
            "campaign references an unreplayable relation profile"))
        inputs = MoveGraphStore._rows(store.db, """
            SELECT input.ordinal,input.state_id,base.complexity
            FROM ograph_relation_campaign_inputs input
            JOIN ograph_states ext ON ext.state_id=input.state_id
            JOIN states base ON base.state_id=input.state_id
            WHERE input.campaign_id=? ORDER BY input.ordinal
        """, (String(campaign.campaign_id),))
        ids = String[input.state_id for input in inputs]
        expected_ordinals = collect(0:length(ids)-1)
        Int[input.ordinal for input in inputs] == expected_ordinals ||
            throw(ArgumentError("campaign input ordinals are not contiguous from zero"))
        ids == sort(unique(ids)) || throw(ArgumentError(
            "campaign inputs are not in deterministic state-id order"))
        length(ids) == Int(campaign.input_count) || throw(ArgumentError(
            "campaign input_count does not match its input rows"))
        for input in inputs
            complexity = Float64(input.complexity)
            isfinite(complexity) && isinteger(complexity) || throw(ArgumentError(
                "campaign input complexity is not a finite integer"))
            complexity <= Int(campaign.n_input) || throw(ArgumentError(
                "campaign input exceeds N_input"))
        end
        Int(campaign.b_path) >= Int(campaign.n_input) || throw(ArgumentError(
            "campaign B_path is smaller than N_input"))
        caps_json, counts_json = String(campaign.caps_json), String(campaign.counts_json)
        caps, counts = JSON.parse(caps_json), JSON.parse(counts_json)
        _ograph_canonical_json(caps) == caps_json || throw(ArgumentError(
            "campaign caps are not canonical JSON"))
        _ograph_canonical_json(counts) == counts_json || throw(ArgumentError(
            "campaign counts are not canonical JSON"))
        input_fingerprint = _ograph_sha256(_ograph_canonical_json(ids))
        campaign_data = Dict(
            "campaign_id" => String(campaign.campaign_id),
            "profile_id" => profile_id,
            "N_input" => Int(campaign.n_input),
            "B_path" => Int(campaign.b_path),
            "input_state_ids" => ids,
            "caps" => caps,
            "input_fingerprint" => input_fingerprint,
            "engine_version" => String(campaign.engine_version),
            "verifier_version" => String(campaign.verifier_version),
            "status" => String(campaign.status),
            "stop_reason" => String(campaign.stop_reason),
            "counts" => counts,
        )
        expected = (campaign_id=String(campaign.campaign_id), profile_id=profile_id,
            n_input=Int(campaign.n_input), b_path=Int(campaign.b_path),
            input_count=length(ids), caps_json=caps_json,
            input_fingerprint=input_fingerprint,
            campaign_fingerprint=_ograph_sha256(_ograph_canonical_json(campaign_data)),
            engine_version=String(campaign.engine_version),
            verifier_version=String(campaign.verifier_version),
            status=String(campaign.status), stop_reason=String(campaign.stop_reason),
            counts_json=counts_json, created_at_utc=String(campaign.created_at_utc))
        _assert_campaign_row(store, campaign, expected, ids)
        haskey(counts, "input_states") && Int(counts["input_states"]) != length(ids) &&
            throw(ArgumentError("campaign declared input-state count is false"))
    end

    ograph_transitions = Int(MoveGraphStore._one(store.db, """
        SELECT COUNT(*) AS n FROM transitions transition
        WHERE transition.relation_profile_id IS NOT NULL
           OR EXISTS (SELECT 1 FROM ograph_states
                      WHERE state_id=transition.source_state_id)
           OR EXISTS (SELECT 1 FROM ograph_states
                      WHERE state_id=transition.target_state_id)
    """).n)
    certificates = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM ograph_relation_certificates").n)
    certificate_steps = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM ograph_relation_certificate_steps").n)
    ograph_transitions == 0 || throw(ArgumentError(
        "Area-3 transition replay is not implemented in o-graph schema v1"))
    certificates == 0 && certificate_steps == 0 || throw(ArgumentError(
        "Area-3 certificate replay is not implemented in o-graph schema v1"))
    for campaign in campaigns
        counts = JSON.parse(String(campaign.counts_json))
        haskey(counts, "transitions") && Int(counts["transitions"]) != 0 &&
            throw(ArgumentError("campaign declares unsupported o-graph transitions"))
        haskey(counts, "certificates") && Int(counts["certificates"]) != 0 &&
            throw(ArgumentError("campaign declares unsupported o-graph certificates"))
    end

    (
        base_schema_version=base_version,
        ograph_schema_version=OGRAPH_STORE_SCHEMA_VERSION,
        integrity=integrity,
        states=length(state_rows),
        sources=length(source_rows),
        profiles=length(profile_rows),
        campaigns=length(campaigns),
        transitions=ograph_transitions,
        certificates=certificates,
    )
end
