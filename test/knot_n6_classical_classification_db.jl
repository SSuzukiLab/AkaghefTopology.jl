using Test
using JSON
using SQLite
using DBInterface

include(joinpath(@__DIR__, "..", "scripts",
                 "build_knot_n6_classical_classification_db.jl"))

function _n6cc_fixture_count(db, sql)
    rows = MoveGraphStore._rows(db, sql)
    Int(first(rows).n)
end

@testset "N6 classical-subset DB fixture" begin
    prepared = _n6cc_run_in_memory(3; progress=false, census=false)
    classification = prepared.classification
    bridge = prepared.bridge
    counts = prepared.counts

    @test length(classification.states) == 73
    @test classification.state_counts == [1, 2, 10, 60]
    @test classification.raw_rooted_classical_counts_by_n == [1, 4, 32, 336]
    @test isempty(classification.unresolved_components)
    @test isempty(bridge.seed_components)
    @test isempty(bridge.resolved_components)
    @test isempty(bridge.unresolved_components)
    @test counts["canonical_state_total"] == 73
    @test counts["planar_true"] == 73
    @test counts["planar_false"] == 0
    @test counts["classified"] == 73
    @test counts["unresolved"] == 0
    @test counts["transitions"] == 0
    @test counts["omitted_nonclassical_canonical_total"] == 108

    fixture_dir = mktempdir()
    fixture_db = joinpath(fixture_dir, "fixture.sqlite3")
    try
        _n6cc_persist_db(fixture_db, classification, bridge, prepared.labels,
                         counts, prepared.census, "fixture-snapshot-sha256")
        store = open_store(fixture_db)
        try
            db = store.db
            @test _n6cc_fixture_count(db, "SELECT COUNT(*) AS n FROM states") == 73
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM state_classifications") == 73
            @test _n6cc_fixture_count(db,
                "SELECT SUM(planar=1) AS n FROM state_classifications") == 73
            @test _n6cc_fixture_count(db,
                "SELECT SUM(planar=0) AS n FROM state_classifications") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM transitions") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM certificate_paths") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM invariants") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM states WHERE payload IS NOT NULL") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM classification_runs") == 1

            columns = String.(getproperty.(MoveGraphStore._rows(
                db, "PRAGMA table_info(state_classifications)"), :name))
            @test columns == ["state_id", "planar", "knot_type",
                              "classification_status", "method_version"]
            @test MoveGraphStore.classification_counts(
                store, N6CC_METHOD_VERSION) == Dict{String,Any}(
                    "total" => 73, "planar_true" => 73, "planar_false" => 0,
                    "classified" => 73, "unresolved" => 0,
                    "nonclassical" => 0,
                    "knot_type_counts" => Dict("0_1" => 71,
                                                  "3_1^+" => 1,
                                                  "3_1^-" => 1))
            run = MoveGraphStore.get_classification_run(
                store, N6CC_METHOD_VERSION)
            @test run.max_n == 3
            @test run.barrier == 7
            @test run.conventions["scope"] ==
                  "classical canonical subset; one-component diagrams; N<=6"
        finally
            close_store!(store)
        end
    finally
        rm(fixture_dir; recursive=true, force=true)
    end
end

@testset "production constants and omitted population" begin
    @test N6CC_METHOD_VERSION == "knot-n6-classical-classification-v1"
    @test N6CC_METHOD_RELEASED_AT_UTC == "2026-08-31T16:16:30Z"
    @test N6CC_ALL_CANONICAL_STATE_COUNTS ==
          [1, 2, 14, 164, 3388, 96776, 3548876]
    @test sum(N6CC_ALL_CANONICAL_STATE_COUNTS) == 3649221
end

@testset "generated N6 classical-subset artifacts" begin
    if isfile(N6CC_JSON_OUTPUT) && isfile(N6CC_DB_OUTPUT)
        payload = JSON.parsefile(N6CC_JSON_OUTPUT)
        counts = payload["result"]["counts"]
        @test payload["schema_version"] == 2
        @test payload["scope"] == "classical canonical subset"
        @test payload["method_version"] == N6CC_METHOD_VERSION
        @test counts["canonical_state_total"] == 77989
        @test counts["planar_true"] == 77989
        @test counts["planar_false"] == 0
        @test counts["classified"] + counts["unresolved"] == 77989
        @test counts["transitions"] == 0
        @test counts["omitted_nonclassical_canonical_total"] == 3571232
        @test length(payload["input"]["snapshot_sha256"]) == 64
        @test length(payload["input"]["classical_census_raw_fingerprint_sha256"]) == 64
        @test length(payload["input"]["classical_census_canonical_fingerprint_sha256"]) == 64
        @test length(payload["input"]["barrier6_graph_fingerprint_sha256"]) == 64
        @test length(payload["input"]["barrier7_bridge_fingerprint_sha256"]) == 64

        db = SQLite.DB(N6CC_DB_OUTPUT)
        try
            @test _n6cc_fixture_count(db, "SELECT COUNT(*) AS n FROM states") == 77989
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM state_classifications") == 77989
            @test _n6cc_fixture_count(db,
                "SELECT SUM(planar=1) AS n FROM state_classifications") == 77989
            @test _n6cc_fixture_count(db,
                "SELECT SUM(planar=0) AS n FROM state_classifications") == 0
            @test _n6cc_fixture_count(db,
                "SELECT SUM(classification_status='classified') AS n FROM state_classifications") == counts["classified"]
            @test _n6cc_fixture_count(db,
                "SELECT SUM(classification_status='unresolved') AS n FROM state_classifications") == counts["unresolved"]
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM transitions") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM certificate_paths") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM invariants") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM states WHERE payload IS NOT NULL") == 0
            @test _n6cc_fixture_count(db,
                "SELECT COUNT(*) AS n FROM classification_runs") == 1
            row = first(MoveGraphStore._rows(db,
                "SELECT max_n,barrier FROM classification_runs"))
            @test (Int(row.max_n), Int(row.barrier)) == (6, 7)
            columns = String.(getproperty.(MoveGraphStore._rows(
                db, "PRAGMA table_info(state_classifications)"), :name))
            @test columns == ["state_id", "planar", "knot_type",
                              "classification_status", "method_version"]
        finally
            SQLite.close(db)
        end
        @test bytes2hex(SHA.sha256(read(N6CC_DB_OUTPUT))) ==
              payload["database"]["sha256"]
        @test filesize(N6CC_DB_OUTPUT) == payload["database"]["bytes"]
    else
        @test true  # Full artifacts are checked after the production build.
    end
end
