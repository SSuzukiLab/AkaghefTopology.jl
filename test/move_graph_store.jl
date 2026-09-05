using Test
using SQLite
using DBInterface
using JSON

include(joinpath(@__DIR__, "..", "src", "move_graph_store.jl"))
using .MoveGraphStore

const FIXTURE = normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260830",
                                  "data", "knotinfo-n6.json"))

struct DualVerifier
    calls::Base.RefValue{Int}
end
function (v::DualVerifier)(source, target, move)
    v.calls[] += 1
    1 + "internal MethodError"
end
(::DualVerifier)() = true

function legacy_v1_store(; fault=:none, parameters=:certificate)
    db = SQLite.DB()
    store = MoveGraphStoreDB(db)
    for sql in (
        "CREATE TABLE schema_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)",
        "CREATE TABLE states (id INTEGER PRIMARY KEY, state_id TEXT NOT NULL UNIQUE, state_hash TEXT NOT NULL, canonical_version INTEGER NOT NULL CHECK (canonical_version > 0), serialization TEXT NOT NULL, payload BLOB, area REAL NOT NULL CHECK (area IN (2.0, 2.5, 3.0)), representation TEXT NOT NULL, complexity REAL, UNIQUE(state_hash, canonical_version, serialization))",
        "CREATE TABLE move_types (move_type TEXT PRIMARY KEY, area REAL NOT NULL CHECK (area IN (2.0, 2.5, 3.0)), inverse_type TEXT NOT NULL, inverse_direction TEXT NOT NULL CHECK (inverse_direction IN ('forward', 'inverse')), UNIQUE(move_type, area))",
        "CREATE TABLE transitions (id INTEGER PRIMARY KEY, source_state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE, target_state_id TEXT NOT NULL REFERENCES states(state_id) ON DELETE CASCADE, move_type TEXT NOT NULL REFERENCES move_types(move_type), direction TEXT NOT NULL CHECK (direction IN ('forward', 'inverse')), location TEXT, parameters TEXT, engine_version TEXT NOT NULL, verifier_version TEXT NOT NULL, verified INTEGER NOT NULL CHECK (verified = 1), UNIQUE(source_state_id, target_state_id, move_type, direction, parameters))",
        "CREATE TABLE certificate_paths (id INTEGER PRIMARY KEY, transition_id INTEGER NOT NULL REFERENCES transitions(id) ON DELETE CASCADE, path TEXT NOT NULL, UNIQUE(transition_id, path))",
    )
        DBInterface.execute(db, sql)
    end
    DBInterface.execute(db,
        "INSERT INTO schema_meta(key,value) VALUES ('schema_version','1')")
    source = CanonicalState("legacy-source"; area=AREA_2, representation="KnotGauss")
    target = CanonicalState("legacy-target"; area=AREA_2, representation="KnotGauss")
    for state in (source, target)
        DBInterface.execute(db, """
            INSERT INTO states(state_id,state_hash,canonical_version,serialization,payload,area,representation,complexity)
            VALUES (?,?,?,?,?,?,?,?)
        """, (state.state_id, state.state_hash, state.canonical_version,
                state.serialization, nothing, state.area, state.representation, nothing))
    end
    DBInterface.execute(db,
        "INSERT INTO move_types(move_type,area,inverse_type,inverse_direction) VALUES ('R1',2.0,'R1','inverse')")
    evidence = Dict(
        "adapter_version" => 1,
        "move_type" => fault == :move_mismatch ? "R2" : "R1",
        "source_serialization" => fault == :source_mismatch ? "wrong-source" : source.serialization,
        "target_serialization" => target.serialization,
        "theorem" => "legacy independently verified Reidemeister witness",
    )
    encoded = fault == :malformed ? "not-json" : JSON.json(evidence)
    stored_parameters = parameters === :certificate ? encoded :
                        parameters === :null ? nothing :
                        parameters === :distinct ? "{\"move_parameter\":1}" :
                        throw(ArgumentError("unsupported fixture parameters mode"))
    DBInterface.execute(db, """
        INSERT INTO transitions(source_state_id,target_state_id,move_type,direction,parameters,engine_version,verifier_version,verified)
        VALUES (?,?,?,?,?,?,?,1)
    """, (source.state_id, target.state_id, "R1", "forward", stored_parameters,
            "Topology/1", "Reidemeister/1"))
    transition_id = Int(first(DBInterface.execute(db,
        "SELECT id FROM transitions")).id)
    if fault != :missing
        DBInterface.execute(db,
            "INSERT INTO certificate_paths(transition_id,path) VALUES (?,?)",
            (transition_id, encoded))
    end
    if fault == :multiple
        second = JSON.json(merge(evidence, Dict("theorem" => "second witness")))
        DBInterface.execute(db,
            "INSERT INTO certificate_paths(transition_id,path) VALUES (?,?)",
            (transition_id, second))
    end
    store, encoded
end

@testset "move graph v1 to v2 migration" begin
    store, evidence = legacy_v1_store()
    migrate!(store)
    @test MoveGraphStore._one(store.db,
        "SELECT value FROM schema_meta WHERE key='schema_version'").value == "2"
    @test MoveGraphStore._one(store.db,
        "SELECT verification_evidence FROM transitions").verification_evidence == evidence
    @test length(MoveGraphStore._rows(store.db, "SELECT * FROM transitions")) == 1
    @test length(MoveGraphStore._rows(store.db, "SELECT * FROM certificate_paths")) == 1
    @test isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check"))
    migrate!(store)
    @test MoveGraphStore._one(store.db,
        "SELECT verification_evidence FROM transitions").verification_evidence == evidence
    close_store!(store)

    for parameters in (:null, :distinct)
        compatible, compatible_evidence = legacy_v1_store(; parameters=parameters)
        migrate!(compatible)
        row = MoveGraphStore._one(compatible.db,
            "SELECT parameters,verification_evidence FROM transitions")
        @test row.verification_evidence == compatible_evidence
        @test parameters === :null ? row.parameters === missing :
              String(row.parameters) == "{\"move_parameter\":1}"
        close_store!(compatible)
    end

    for fault in (:missing, :multiple, :malformed, :move_mismatch, :source_mismatch)
        bad, _ = legacy_v1_store(; fault=fault)
        @test_throws ArgumentError migrate!(bad)
        @test MoveGraphStore._one(bad.db,
            "SELECT value FROM schema_meta WHERE key='schema_version'").value == "1"
        columns = Set(String(row.name) for row in
                      MoveGraphStore._rows(bad.db, "PRAGMA table_info(transitions)"))
        @test !("verification_evidence" in columns)
        close_store!(bad)
    end

    future = MoveGraphStoreDB(SQLite.DB())
    DBInterface.execute(future.db,
        "CREATE TABLE schema_meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    DBInterface.execute(future.db,
        "INSERT INTO schema_meta(key,value) VALUES ('schema_version','3')")
    DBInterface.execute(future.db, "CREATE TABLE marker (value TEXT)")
    schema_rows() = [(name=String(row.name),
                      sql=row.sql === missing ? nothing : String(row.sql))
                     for row in MoveGraphStore._rows(future.db,
                         "SELECT name,sql FROM sqlite_master ORDER BY name")]
    schema_before = schema_rows()
    @test_throws ArgumentError migrate!(future)
    @test schema_rows() == schema_before
    close_store!(future)

    unknown = MoveGraphStoreDB(SQLite.DB())
    DBInterface.execute(unknown.db, "CREATE TABLE marker (value TEXT)")
    @test_throws ArgumentError migrate!(unknown)
    @test length(MoveGraphStore._rows(unknown.db,
        "SELECT name FROM sqlite_master WHERE type='table' AND name='marker'")) == 1
    close_store!(unknown)
end

@testset "move graph store" begin
    store = open_store()
    migrate!(store)
    @test first(MoveGraphStore._rows(store.db, "SELECT value FROM schema_meta WHERE key='schema_version'")).value == "2"

    @test AREA_2 == 2.0
    @test AREA_2_5 == 2.5
    @test AREA_3 == 3.0
    @test length(unique([AREA_2, AREA_2_5, AREA_3])) == 3
    @test_throws ArgumentError register_move_type!(store, "representation_transform"; area=AREA_2)
    # Move families belong to Area 3, while representation/object names never
    # become graph edges.
    for move_type in ("MP", "PS", "CP")
        @test register_move_type!(store, move_type; area=AREA_3).area == AREA_3
        @test_throws ArgumentError register_move_type!(store, move_type; area=AREA_2)
        @test_throws ArgumentError register_move_type!(store, move_type; area=AREA_2_5)
    end
    for non_move in ("o-data", "odata", "3-manifold")
        @test_throws ArgumentError register_move_type!(store, non_move; area=AREA_3)
    end

    a = CanonicalState("state-a"; area=AREA_2, representation="KnotGauss", payload=[1, 2])
    @test_throws ArgumentError CanonicalState("state-a"; area=AREA_2, representation="KnotGauss", state_id="different-id")
    @test_throws ArgumentError CanonicalState("state-a"; area=AREA_2, representation="KnotGauss", state_hash="different-hash")
    @test_throws ArgumentError CanonicalState("bad-id", "bad-hash", 1, "state-a", nothing, AREA_2, "knotgauss", nothing)
    b = CanonicalState("state-b"; area=AREA_2, representation="KnotGauss")
    @test put_state!(store, a).state_id == a.state_id
    @test length(MoveGraphStore._rows(store.db, "SELECT id FROM states")) == 1
    DBInterface.execute(store.db, "UPDATE states SET state_hash='corrupt' WHERE state_id=?", (a.state_id,))
    @test_throws ArgumentError get_state(store, a.state_id)
    DBInterface.execute(store.db, "UPDATE states SET state_hash=? WHERE state_id=?", (state_hash(a), a.state_id))

    @test_throws ArgumentError insert_verified_transition!(store, a, b, "R3"; verified=true)
    @test_throws ArgumentError insert_verified_transition!(store, a, CanonicalState("c"; area=AREA_3, representation="o-data"), "R3";
                                                           verifier=(args...)->true, verified=true)
    unregistered_target = CanonicalState("unregistered-target";
        area=AREA_2, representation="KnotGauss")
    @test_throws ArgumentError insert_verified_transition!(store, a,
        unregistered_target, "R2"; verifier=(args...)->false, verified=true)
    @test get_move_type(store, "R2") === nothing
    @test get_state(store, unregistered_target.state_id) === nothing

    edge = insert_verified_transition!(store, a, b, "R3"; direction=:forward,
                                       location="source:1",
                                       engine_version="engine/1",
                                       verifier_version="verifier/1",
                                       evidence=Dict("certificate" => "strong"),
                                       verifier=(args...)->true, verified=true,
                                       certificate_paths=["cert/a-b"])
    @test edge == insert_verified_transition!(store, a, b, "R3"; direction=:forward,
                                               verifier=(args...)->true, verified=true)
    retained = MoveGraphStore._one(store.db,
        "SELECT location,engine_version,verifier_version,verification_evidence FROM transitions WHERE id=?",
        (edge,))
    @test String(retained.location) == "source:1"
    @test String(retained.engine_version) == "engine/1"
    @test String(retained.verifier_version) == "verifier/1"
    @test JSON.parse(String(retained.verification_evidence)) ==
          Dict("certificate" => "strong")
    @test length(neighbors(store, a.state_id)) == 1
    @test bfs_certificate_path(store, a.state_id, b.state_id) == [a.state_id, b.state_id]
    inverse = insert_verified_transition!(store, b, a, "R3"; direction=:inverse,
                                          verification_result=(verified=true,), verified=true)
    @test MoveGraphStore._one(store.db, "SELECT verification_evidence FROM transitions WHERE id=?", (edge,)).verification_evidence !== nothing
    @test inverse != edge
    @test length(neighbors(store, b.state_id)) == 1

    # An internal MethodError from an applicable 3-argument verifier must not
    # be mistaken for an unsupported arity and retried through the permissive
    # zero-argument fallback.
    bad_target = CanonicalState("bad-target"; area=AREA_2, representation="KnotGauss")
    calls = Ref(0)
    @test_throws MethodError insert_verified_transition!(store, a, bad_target, "R3";
                                                         verifier=DualVerifier(calls), verified=true)
    @test calls[] == 1
    @test get_state(store, bad_target.state_id) === nothing
    @test isempty(MoveGraphStore._rows(store.db, "SELECT * FROM transitions WHERE target_state_id=?",
                                       (bad_target.state_id,)))

    # The SQL layer independently rejects an unverified edge, even when all
    # foreign keys and other edge fields are otherwise valid.
    sql_target = CanonicalState("sql-target"; area=AREA_2, representation="KnotGauss")
    put_state!(store, sql_target)
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "INSERT INTO transitions(source_state_id,target_state_id,move_type,direction,parameters,engine_version,verifier_version,verified) VALUES (?,?,?,?,?,?,?,1)",
        (a.state_id, sql_target.state_id, "R3", "forward", "", "test", "test"))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "INSERT INTO transitions(source_state_id,target_state_id,move_type,direction,parameters,engine_version,verifier_version,verification_evidence,verified) VALUES (?,?,?,?,?,?,?,?,1)",
        (a.state_id, sql_target.state_id, "R3", "forward", "x", "test", "test", ""))

    # Failed work inside the transaction must not leave either endpoint/edge.
    x = CanonicalState("rollback-x"; area=AREA_2, representation="KnotGauss")
    y = CanonicalState("rollback-y"; area=AREA_2, representation="KnotGauss")
    @test_throws ErrorException MoveGraphStore._transaction(() -> begin
        put_state!(store, x)
        put_state!(store, y)
        error("intentional rollback")
    end, store.db)
    @test get_state(store, x.state_id) === nothing
    @test get_state(store, y.state_id) === nothing

    # Area 2.5 and Area 3 remain separate from Area 2.
    m25 = CanonicalState("mixed"; area=AREA_2_5, representation="mixed/lo-graph")
    m25b = CanonicalState("mixed-b"; area=AREA_2_5, representation="mixed/lo-graph")
    o3 = CanonicalState("odata"; area=AREA_3, representation="o-data")
    o3b = CanonicalState("odata-b"; area=AREA_3, representation="Heegaard")
    @test register_move_type!(store, "mixed"; area=AREA_2_5).area == AREA_2_5
    @test register_move_type!(store, "MP"; area=AREA_3).area == AREA_3
    @test_throws ArgumentError insert_verified_transition!(store, m25, m25b, "R3";
                                                           verifier=(args...)->true, verified=true)
    @test insert_verified_transition!(store, m25, m25b, "mixed";
                                      verifier=(args...)->true, verified=true) isa Integer
    @test insert_verified_transition!(store, o3, o3b, "MP";
                                      verifier=(args...)->true, verified=true) isa Integer

    imported = import_knotinfo!(store, FIXTURE)
    @test imported == 8
    @test length(knotinfo_records(store)) == 8
    names = [String(r.name) for r in knotinfo_records(store)]
    @test sort(names) == sort(["0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"])
    @test isempty(MoveGraphStore._rows(store.db, "SELECT name FROM knot_types WHERE state_id IS NOT NULL"))
    @test first(MoveGraphStore._rows(store.db, "SELECT snapshot_sha256 FROM knotinfo_snapshots")).snapshot_sha256 ==
          "578be92fc1d37b1351c8273c3beab00654913ea009ca751a71e99bac2209751f"

    # The database must independently enforce the common-store domain boundary
    # and retain verifier evidence, not merely the verified flag.
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "INSERT INTO transitions(source_state_id,target_state_id,move_type,direction,parameters,engine_version,verifier_version,verified) VALUES (?,?,?,?,?,?,?,1)",
        (a.state_id, m25b.state_id, "R3", "forward", "", "test", "test"))
    evidenced = insert_verified_transition!(store, a, CanonicalState("evidence-target"; area=AREA_2, representation="KnotGauss"), "R3";
        verifier=(args...)->true, verified=true, verification_result=(verified=true, certificate="proof-1"), evidence=(certificate="proof-1", score=1))
    @test MoveGraphStore._one(store.db, "SELECT verification_evidence FROM transitions WHERE id=?", (evidenced,)).verification_evidence !== nothing
    close_store!(store)
end
