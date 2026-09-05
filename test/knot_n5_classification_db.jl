using Test
using SQLite
using DBInterface
using JSON
using SHA
using Topology

const N5_DB_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const N5_DB_JSON = joinpath(N5_DB_REPO_ROOT, "threads", "PJ260830", "data",
                            "knot-n5-classification.json")
const N5_DB_PATH = joinpath(N5_DB_REPO_ROOT, "references", "items", "database",
                            "knotinfo", "knot-n5-classification.sqlite3")

@testset "classification store v2 API" begin
    store = Topology.open_store()
    try
        schema_version = first(Topology.MoveGraphStore._rows(
            store.db, "SELECT value FROM schema_meta WHERE key='schema_version'"))
        @test schema_version.value == "2"

        run = Topology.register_classification_run!(store;
            method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9,
            census_fingerprint="census-fixture", classification_fingerprint="class-fixture",
            counts=Dict("total" => 2), conventions=Dict("knot_type" => "arbitrary"),
            created_at_utc="2026-08-31T00:00:00Z")
        @test run isa Topology.ClassificationRun
        @test run.max_n == 7 && run.barrier == 9
        @test Topology.get_classification_run(store, "fixture-v2").method_version ==
              "fixture-v2"
        @test Topology.register_classification_run!(store;
            method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9,
            census_fingerprint="census-fixture", classification_fingerprint="class-fixture",
            counts=Dict("total" => 2), conventions=Dict("knot_type" => "arbitrary"),
            created_at_utc="different-created-at").method_version == "fixture-v2"
        @test_throws ArgumentError Topology.register_classification_run!(store;
            method_version="fixture-v2", canonical_version=1, max_n=7, barrier=9,
            census_fingerprint="changed", classification_fingerprint="class-fixture",
            counts=Dict("total" => 2), conventions=Dict("knot_type" => "arbitrary"))
        @test_throws ArgumentError Topology.register_classification_run!(store;
            method_version="bad-bound", canonical_version=1, max_n=8, barrier=7,
            census_fingerprint="x", classification_fingerprint="y")

        area2 = Topology.CanonicalState("fixture-area2";
            area=Topology.AREA_2, representation="KnotGauss")
        area25 = Topology.CanonicalState("fixture-area25";
            area=Topology.AREA_2_5, representation="mixed/lo-graph")
        area3 = Topology.CanonicalState("fixture-area3";
            area=Topology.AREA_3, representation="o-data")
        Topology.put_state!(store, area2)
        Topology.put_state!(store, area25)
        Topology.put_state!(store, area3)

        @test Topology.put_state_classification!(store, area2;
            planar=true, knot_type="future-6_1", classification_status="classified",
            method_version="fixture-v2").knot_type == "future-6_1"
        # The Julia boundary rejects Area 2.5/3 before SQL can persist a row.
        @test_throws ArgumentError Topology.put_state_classification!(store, area25;
            planar=false, classification_status="nonclassical", method_version="fixture-v2")
        @test_throws ArgumentError Topology.put_state_classification!(store, area3;
            planar=false, classification_status="nonclassical", method_version="fixture-v2")

        # SQL independently enforces the Area-2 boundary and field consistency.
        @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
            INSERT INTO state_classifications(state_id,planar,knot_type,classification_status,method_version)
            VALUES (?,?,?,?,?)
        """, (area25.state_id, 0, nothing, "nonclassical", "fixture-v2"))
        @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
            INSERT INTO state_classifications(state_id,planar,knot_type,classification_status,method_version)
            VALUES (?,?,?,?,?)
        """, (area2.state_id, 0, "not-null-for-nonclassical", "nonclassical", "fixture-v2"))
        @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
            INSERT INTO state_classifications(state_id,planar,knot_type,classification_status,method_version)
            VALUES (?,?,?,?,?)
        """, (area2.state_id, 1, "", "classified", "fixture-v2"))

        @test Topology.get_state_classification(store, area2.state_id).knot_type ==
              "future-6_1"
        @test Topology.classification_counts(store, "fixture-v2") == Dict{String,Any}(
            "total" => 1, "planar_true" => 1, "planar_false" => 0,
            "classified" => 1, "unresolved" => 0, "nonclassical" => 0,
            "knot_type_counts" => Dict("future-6_1" => 1))

        @test Topology.put_state_classification!(store, area2;
            planar=true, knot_type="future-6_1", classification_status="classified",
            method_version="fixture-v2").state_id == area2.state_id
        @test_throws ArgumentError Topology.put_state_classification!(store, area2;
            planar=true, knot_type="future-6_1", classification_status="unresolved",
            method_version="fixture-v2")

        # A second method makes an implicit lookup ambiguous rather than
        # silently selecting an arbitrary run.
        Topology.register_classification_run!(store;
            method_version="fixture-v2-b", canonical_version=1, max_n=0, barrier=0,
            census_fingerprint="census-fixture-b", classification_fingerprint="class-fixture-b")
        Topology.put_state_classification!(store, area2;
            planar=true, knot_type="other", classification_status="classified",
            method_version="fixture-v2-b")
        @test_throws ArgumentError Topology.get_state_classification(store, area2.state_id)

        columns = String.(getproperty.(Topology.MoveGraphStore._rows(
            store.db, "PRAGMA table_info(state_classifications)"), :name))
        @test columns == ["state_id", "planar", "knot_type", "classification_status", "method_version"]
        forbidden = Set(["genus", "face", "euler", "rotation", "face_cycles",
                         "certificate", "path", "parent", "alexander", "jones"])
        @test isempty(intersect(Set(lowercase.(columns)), forbidden))
        @test first(Topology.MoveGraphStore._rows(
            store.db, "SELECT COUNT(*) AS n FROM transitions")).n == 0
    finally
        Topology.close_store!(store)
    end
end

@testset "generated N5 classification DB summary" begin
    if isfile(N5_DB_JSON) && isfile(N5_DB_PATH)
        payload = JSON.parsefile(N5_DB_JSON)
        counts = payload["result"]["counts"]
        @test payload["schema_version"] == 2
        @test counts["canonical_state_total"] == 100345
        @test counts["planar_true"] == 6361
        @test counts["planar_false"] == 93984
        @test counts["classified"] == 6361
        @test counts["nonclassical"] == 93984
        @test counts["unresolved"] == 0
        @test counts["transitions"] == 0
        expected_types = Dict(
            "0_1" => 5973, "3_1^+" => 174, "3_1^-" => 174, "4_1" => 36,
            "5_1^+" => 1, "5_1^-" => 1, "5_2^+" => 1, "5_2^-" => 1)
        expected_raw_types = Dict(
            "0_1" => 58213, "3_1^+" => 1710, "3_1^-" => 1710, "4_1" => 348,
            "5_1^+" => 2, "5_1^-" => 2, "5_2^+" => 10, "5_2^-" => 10)
        @test counts["canonical_knot_type_counts"] == expected_types
        @test counts["raw_rooted_classical_counts"] == expected_raw_types
        db = SQLite.DB(N5_DB_PATH)
        try
            q(sql) = first(Topology.MoveGraphStore._rows(db, sql)).n
            @test q("SELECT COUNT(*) AS n FROM states") == 100345
            @test q("SELECT COUNT(*) AS n FROM state_classifications") == 100345
            @test q("SELECT SUM(planar=1) AS n FROM state_classifications") == 6361
            @test q("SELECT SUM(planar=0) AS n FROM state_classifications") == 93984
            @test q("SELECT COUNT(*) AS n FROM transitions") == 0
            @test q("SELECT COUNT(*) AS n FROM classification_runs") == 1
            @test q("SELECT COUNT(DISTINCT knot_type) AS n FROM state_classifications WHERE knot_type IS NOT NULL") == 8
        finally
            SQLite.close(db)
        end
        @test bytes2hex(SHA.sha256(read(N5_DB_PATH))) == payload["database"]["sha256"]
        @test filesize(N5_DB_PATH) == payload["database"]["bytes"]
    else
        @test true  # Full generated artifacts are checked when the builder runs.
    end
end
