using Test
using Topology
using SQLite
using DBInterface

function _ograph_db_link(code, signs)
    link = VirtualLink()
    set_data!(link; rgauss=[Int.(code)], orientation=Int.(signs))
end

const OGRAPH_TEST_RELEASE = "2026-09-04T00:00:00Z"

@testset "o-graph store schema and finite pilot" begin
    store = open_store()
    migrate_ograph_store!(store)
    migrate_ograph_store!(store)
    @test Topology._ograph_existing_schema_version(store.db) == 2
    transition_columns = Set(String(row.name) for row in
        MoveGraphStore._rows(store.db, "PRAGMA table_info(transitions)"))
    @test "relation_profile_id" in transition_columns
    @test ograph_canonical_json(Dict("b" => 2, "a" => 1)) ==
          ograph_canonical_json(Dict("a" => 1, "b" => 2))
    @test_throws ArgumentError ograph_canonical_json(Dict{Any,Int}(:a => 1, "a" => 2))
    @test_throws ArgumentError ograph_canonical_json(NaN)
    @test_throws ArgumentError ograph_canonical_json(Inf)
    @test ograph_canonical_json(typemax(Int64)) == string(typemax(Int64))
    @test_throws ArgumentError ograph_canonical_json(big(typemax(Int64)) + 1)

    source_file = OGraphProvenance(
        "local-source",
        "threads/PJ260825/sources/MST1/heisenberg-double-04-07.tex:318-336";
        snapshot_sha256=repeat("a", 64), metadata=Dict("claim" => "closed normal o-graph"))
    source_constructor = OGraphProvenance(
        "implementation",
        "lib/topology.jl/src/virtual_link.jl:35-53";
        metadata=Dict("constructor" => "lens_space_link(2,1)"))
    @test source_file.source_id == OGraphProvenance(
        "local-source",
        "threads/PJ260825/sources/MST1/heisenberg-double-04-07.tex:318-336";
        snapshot_sha256=repeat("a", 64), metadata=Dict("claim" => "closed normal o-graph")).source_id
    @test_throws ArgumentError OGraphProvenance("source", "path";
                                               snapshot_sha256="not-a-sha")

    s3 = certify_closed_ograph(_ograph_db_link([1, -1], [1]))
    l21 = certify_closed_ograph(lens_space_link(2, 1))
    s2xs1 = certify_closed_ograph(_ograph_db_link(
        [-1, -2, 3, 4, -4, -3, 2, 1], [1, 1, 1, 1]))
    fresh_forged_source = OGraphProvenance(
        "not-derived", "local-source", "fresh-forgery:1", nothing, "{}")
    @test_throws ArgumentError put_closed_ograph!(store, s3;
        provenance=[fresh_forged_source])
    @test isempty(MoveGraphStore._rows(store.db, "SELECT * FROM states"))
    s3_state = put_closed_ograph!(store, s3; provenance=[source_file])
    l21_state = put_closed_ograph!(store, l21;
                                  provenance=[source_file, source_constructor])
    s2_state = put_closed_ograph!(store, s2xs1; provenance=[source_constructor])
    @test put_closed_ograph!(store, s3; provenance=[source_file]).state_id == s3_state.state_id
    conflicting_source = OGraphProvenance(
        source_file.source_id, "changed-kind", source_file.locator,
        source_file.snapshot_sha256, source_file.metadata_json)
    @test_throws ArgumentError put_closed_ograph!(store, s3;
                                                   provenance=[conflicting_source])
    @test length(MoveGraphStore._rows(store.db, "SELECT * FROM ograph_states")) == 3
    @test length(MoveGraphStore._rows(store.db, "SELECT * FROM ograph_sources")) == 2
    @test all(state.area == AREA_3 for state in (s3_state, l21_state, s2_state))
    @test all(state.representation == "o-graph" for state in (s3_state, l21_state, s2_state))
    @test startswith(s3_state.serialization, "OGraphState;schema=akaghef.ograph/1;")

    presentation = OGraphRelationProfile("presentation";
        moves_complete=true, evidence_status="identity-only",
        conventions=Dict("canonicalization" => OGRAPH_PRESENTATION_CANONICALIZATION),
        created_at_utc=OGRAPH_TEST_RELEASE)
    combed = OGraphRelationProfile("combed";
        parent_profile_id=presentation.profile_id,
        moves_complete=false, evidence_status="unimplemented",
        conventions=Dict("generators" => ["primary MP", "PS with Method-B witness"]),
        created_at_utc=OGRAPH_TEST_RELEASE)
    manifold = OGraphRelationProfile("manifold";
        parent_profile_id=combed.profile_id,
        moves_complete=false, evidence_status="unimplemented",
        conventions=Dict("generators" => ["primary MP", "PS", "CP"]),
        created_at_utc=OGRAPH_TEST_RELEASE)
    @test length(unique([presentation.profile_id, combed.profile_id,
                         manifold.profile_id])) == 3
    @test_throws ArgumentError OGraphRelationProfile("presentation";
        moves_complete=false, evidence_status="unimplemented",
        created_at_utc=OGRAPH_TEST_RELEASE)
    @test_throws ArgumentError OGraphRelationProfile("spin";
        parent_profile_id=manifold.profile_id, moves_complete=false,
        evidence_status="unimplemented", created_at_utc=OGRAPH_TEST_RELEASE)
    fresh_forged_profile = OGraphRelationProfile(
        "not-derived", presentation.profile_name, presentation.profile_version,
        presentation.parent_profile_id, presentation.moves,
        presentation.moves_complete, presentation.evidence_status,
        presentation.conventions_json, presentation.profile_fingerprint,
        presentation.created_at_utc)
    @test_throws ArgumentError register_ograph_profile!(store,
        fresh_forged_profile)
    for profile in (presentation, combed, manifold)
        register_ograph_profile!(store, profile)
        register_ograph_profile!(store, profile)
    end
    @test length(MoveGraphStore._rows(store.db,
        "SELECT * FROM ograph_relation_profiles")) == 3
    @test Set(String(row.move_type) for row in MoveGraphStore._rows(store.db,
        "SELECT move_type FROM ograph_relation_profile_moves")) == Set(["MP", "PS", "CP"])

    campaign = register_ograph_campaign!(store,
        "ograph-presentation-pilot-v1", presentation,
        [s3_state, l21_state, s2_state];
        n_input=4, b_path=4,
        caps=Dict("max_depth" => 0, "transition_budget" => 0,
                  "meaning" => "identity-and-closedness pilot; no move search"),
        engine_version="Topology/ograph-state-1",
        verifier_version=OGRAPH_CLOSEDNESS_VERIFIER_VERSION,
        status="completed", stop_reason="identity-only profile completed",
        counts=Dict("input_states" => 3, "transitions" => 0,
                    "certificates" => 0),
        created_at_utc=OGRAPH_TEST_RELEASE)
    @test campaign == register_ograph_campaign!(store,
        "ograph-presentation-pilot-v1", presentation,
        [s3_state, l21_state, s2_state];
        n_input=4, b_path=4,
        caps=Dict("max_depth" => 0, "transition_budget" => 0,
                  "meaning" => "identity-and-closedness pilot; no move search"),
        engine_version="Topology/ograph-state-1",
        verifier_version=OGRAPH_CLOSEDNESS_VERIFIER_VERSION,
        status="completed", stop_reason="identity-only profile completed",
        counts=Dict("input_states" => 3, "transitions" => 0,
                    "certificates" => 0),
        created_at_utc=OGRAPH_TEST_RELEASE)
    @test_throws ArgumentError register_ograph_campaign!(store,
        "bad-bounds", presentation, [s3_state]; n_input=4, b_path=3,
        caps=Dict(), engine_version="test", verifier_version="test",
        status="failed", stop_reason="test", counts=Dict(),
        created_at_utc=OGRAPH_TEST_RELEASE)

    # The extension independently blocks an Area-3 edge without a relation
    # profile; a generic `(args...)->true` callback cannot bypass this SQL gate.
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
        INSERT INTO transitions(
          source_state_id,target_state_id,move_type,direction,parameters,
          engine_version,verifier_version,verification_evidence,verified)
        VALUES (?,?,?,?,?,?,?,?,1)
    """, (s3_state.state_id, l21_state.state_id, "MP", "forward", "{}",
            "test", "test", "{}"))

    # Area 3 is broader than the o-graph representation.  The extension does
    # not seize unrelated Area-3 transitions merely because their area is 3.
    other_a = CanonicalState("other-area3-a"; area=AREA_3, representation="Heegaard")
    other_b = CanonicalState("other-area3-b"; area=AREA_3, representation="Heegaard")
    @test insert_verified_transition!(store, other_a, other_b, "MP";
        verifier=(args...)->true, verified=true) isa Integer

    summary = verify_ograph_store(store)
    @test summary == (
        base_schema_version=2,
        ograph_schema_version=2,
        integrity="ok",
        states=3,
        sources=2,
        profiles=3,
        campaigns=1,
        transitions=0,
        certificates=0,
    )

    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "UPDATE ograph_states SET closedness_evidence='{}' WHERE state_id=?",
        (s3_state.state_id,))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "UPDATE states SET complexity=99 WHERE state_id=?", (s3_state.state_id,))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "DELETE FROM states WHERE state_id=?", (s3_state.state_id,))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "UPDATE ograph_sources SET locator='changed' WHERE source_id=?",
        (source_file.source_id,))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "DELETE FROM ograph_state_sources WHERE state_id=? AND source_id=?",
        (s3_state.state_id, source_file.source_id))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "DELETE FROM ograph_relation_campaign_inputs WHERE campaign_id=?",
        (campaign.campaign_id,))
    @test verify_ograph_store(store).states == 3

    # Schema v2 refuses the row at its SQL boundary; a nominal `verified=1`
    # cannot enter before a typed source-bound move verifier exists.
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
        INSERT INTO transitions(
          source_state_id,target_state_id,move_type,direction,parameters,
          engine_version,verifier_version,verification_evidence,verified,
          relation_profile_id)
        VALUES (?,?,?,?,?,?,?,?,1,?)
    """, (s3_state.state_id, l21_state.state_id, "MP", "forward", "{}",
            "test", "unimplemented-test", "{\"synthetic\":true}",
            combed.profile_id))
    close_store!(store)
end

@testset "o-graph extension v1 to v2 replay gate" begin
    valid = open_store()
    migrate_ograph_store!(valid)
    DBInterface.execute(valid.db,
        "UPDATE ograph_schema_meta SET value='1' WHERE key='schema_version'")
    migrate_ograph_store!(valid)
    @test Topology._ograph_existing_schema_version(valid.db) == 2
    @test verify_ograph_store(valid).states == 0
    close_store!(valid)

    invalid = open_store()
    migrate_ograph_store!(invalid)
    DBInterface.execute(invalid.db,
        "INSERT INTO ograph_sources VALUES (?,?,?,?,?)",
        ("orphan", "test", "test:1", missing, "{}"))
    DBInterface.execute(invalid.db,
        "UPDATE ograph_schema_meta SET value='1' WHERE key='schema_version'")
    @test_throws ArgumentError migrate_ograph_store!(invalid)
    @test Topology._ograph_existing_schema_version(invalid.db) == 1
    close_store!(invalid)
end

@testset "o-graph extension future-version rejection" begin
    store = open_store()
    DBInterface.execute(store.db,
        "CREATE TABLE ograph_schema_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    DBInterface.execute(store.db,
        "INSERT INTO ograph_schema_meta(key,value) VALUES ('schema_version','3')")
    before = Set(String(row.name) for row in
        MoveGraphStore._rows(store.db, "PRAGMA table_info(transitions)"))
    @test !("relation_profile_id" in before)
    @test_throws ArgumentError migrate_ograph_store!(store)
    after = Set(String(row.name) for row in
        MoveGraphStore._rows(store.db, "PRAGMA table_info(transitions)"))
    @test after == before
    close_store!(store)
end
