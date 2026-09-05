using Test
using JSON
using SHA
using SQLite
using DBInterface

include(joinpath(@__DIR__, "..", "scripts",
                 "build_knot_n6_full_classification_db.jl"))

function _n6full_test_count(db, sql)
    rows = MoveGraphStore._rows(db, sql)
    isempty(rows) && throw(ErrorException("count query returned no row"))
    value = first(rows).n
    value === nothing ? 0 : Int(value)
end

function _n6full_test_state_rows(path)
    store = MoveGraphStore.open_store(path)
    try
        rows = MoveGraphStore._rows(store.db, """
            SELECT state_id,state_hash,canonical_version,serialization,payload,
                   area,representation,complexity
            FROM states ORDER BY id
        """)
        return [(String(row.state_id), String(row.state_hash),
                 Int(row.canonical_version), String(row.serialization),
                 row.payload === missing ? nothing : row.payload,
                 Float64(row.area), String(row.representation),
                 Float64(row.complexity)) for row in rows]
    finally
        MoveGraphStore.close_store!(store)
    end
end

function _n6full_test_classification_rows(path, method)
    store = MoveGraphStore.open_store(path)
    try
        rows = MoveGraphStore._rows(store.db, """
            SELECT state_id,planar,knot_type,classification_status,method_version
            FROM state_classifications WHERE method_version=? ORDER BY state_id
        """, (method,))
        return [(String(row.state_id), Int(row.planar),
                 row.knot_type === missing ? nothing : String(row.knot_type),
                 String(row.classification_status), String(row.method_version))
                for row in rows]
    finally
        MoveGraphStore.close_store!(store)
    end
end

function _n6full_test_production_gate_inputs(; unresolved=copy(
        N6FULL_BARRIER8_UNRESOLVED_COMPONENTS),
        barrier8_status=N6FULL_BARRIER8_STATUS,
        barrier8_unresolved=copy(unresolved),
        barrier8_fingerprint=N6FULL_BARRIER8_FINGERPRINT)
    classification = KnotN6Classification(
        6, 6, Int[], Int[], 0, 0, KnotGaussDiagram[],
        Dict{KnotGaussDiagram,Int}(), KnotN6Edge[], Vector{Vector{Int}}(),
        Int[], Int[], Dict{Int,String}(), Dict{String,Int}(), Dict{String,Int}(),
        0, KnotN6ComponentInvariant[], KnotN6KnownRepresentative[], Int[],
        0, 0, 0, 0, 0, 0, String[], N6FULL_BARRIER6_FINGERPRINT, 0.0)
    bridge7 = KnotN6TargetedBridgeResult(
        7, classification.fingerprint, Int[], 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, Dict{Int,Vector{String}}(),
        Dict{Int,Union{Nothing,Int}}(), Dict{Int,String}(), Dict{String,Int}(),
        Dict{String,Int}(), Vector{Vector{Int}}(), copy(unresolved),
        N6FULL_BARRIER7_FINGERPRINT, 0.0)
    barrier8 = KnotN6TargetedBarrier8Result(
        8, classification.fingerprint, bridge7.fingerprint, copy(unresolved),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        Dict{Int,Vector{String}}(), Dict{Int,Union{Nothing,Int}}(),
        Dict{Int,String}(), Dict{Int,String}(), Dict{String,Int}(),
        Dict{String,Int}(), Vector{Vector{Int}}(), copy(barrier8_unresolved),
        String(barrier8_status), "", String(barrier8_fingerprint), 0.0)
    classification, bridge7, barrier8
end

@testset "N6 full canonical DB fixture" begin
    classification = build_knot_n6_classification(3; barrier=3, progress=false)
    workspace_ref = Ref{Any}(nothing)
    bridge7 = run_knot_n6_targeted_bridge(
        classification; progress=false,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    barrier8 = run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; progress=false)
    prepared = _n6full_prepared(classification, bridge7; barrier8=barrier8)
    @test classification.state_counts == [1, 2, 10, 60]
    @test classification.raw_rooted_classical_counts_by_n == [1, 4, 32, 336]
    @test length(prepared.classical_map) == 73
    @test length(prepared.labels) == 3
    @test isempty(classification.unresolved_components)
    @test isempty(bridge7.seed_components)
    @test isempty(bridge7.resolved_components)
    @test barrier8.barrier == 8
    @test isempty(barrier8.resolved_components)
    @test isempty(barrier8.unresolved_components)

    fixture_dir = mktempdir()
    fixture_db = joinpath(fixture_dir, "fixture.sqlite3")
    method = "knot-n6-full-classification-fixture-v1"
    try
        result = _n6full_persist_stream(fixture_db, classification, bridge7;
            barrier8=prepared.barrier8, classical_map=prepared.classical_map,
            labels=prepared.labels, max_n=3, progress=false,
            method_version=method, snapshot_sha256="fixture-snapshot-sha256")
        counts = result.counts
        @test counts["canonical_state_counts"] == [1, 2, 14, 164]
        @test counts["canonical_state_total"] == 181
        @test counts["planar_true"] == 73
        @test counts["planar_false"] == 108
        @test counts["classified"] == 73
        @test counts["unresolved"] == 0
        @test counts["nonclassical"] == 108
        @test counts["transitions"] == 0
        @test length(result.census_fingerprint) == 64
        @test length(result.classification_fingerprint) == 64
        @test _n6full_verify_db(fixture_db, counts; method_version=method)

        # A fault after n=1 has committed must leave n=0,1 available for a
        # later invocation.  Resume must skip that prefix and reproduce the
        # clean-build rows and fingerprints exactly.
        resume_db = joinpath(fixture_dir, "resume.sqlite3")
        @test_throws ErrorException _n6full_persist_stream(resume_db,
            classification, bridge7; barrier8=prepared.barrier8,
            classical_map=prepared.classical_map, labels=prepared.labels,
            max_n=3, progress=false, method_version=method,
            snapshot_sha256="fixture-snapshot-sha256", stop_after_n=1)
        checkpoint_store = MoveGraphStore.open_store(resume_db)
        try
            @test _n6full_test_count(checkpoint_store.db,
                "SELECT COUNT(*) AS n FROM states WHERE complexity <= 1") == 3
            @test _n6full_test_count(checkpoint_store.db,
                "SELECT COUNT(*) AS n FROM states WHERE complexity >= 2") == 0
            checkpoint_run = MoveGraphStore.get_classification_run(
                checkpoint_store, method)
            @test checkpoint_run.counts["completed_n"] == 1
            @test length(checkpoint_run.counts["stream_n_stats"]) == 2
        finally
            MoveGraphStore.close_store!(checkpoint_store)
        end
        resumed = _n6full_persist_stream(resume_db, classification, bridge7;
            barrier8=prepared.barrier8, classical_map=prepared.classical_map,
            labels=prepared.labels, max_n=3, progress=false,
            method_version=method, snapshot_sha256="fixture-snapshot-sha256")
        @test resumed.census_fingerprint == result.census_fingerprint
        @test resumed.classification_fingerprint == result.classification_fingerprint
        for key in ("canonical_state_counts", "canonical_state_total",
                    "planar_true", "planar_false", "classified", "unresolved",
                    "nonclassical", "transitions", "stream_canonical_fingerprint_sha256",
                    "stored_state_total", "stored_classification_total")
            @test resumed.counts[key] == counts[key]
        end
        @test _n6full_test_state_rows(resume_db) ==
              _n6full_test_state_rows(fixture_db)
        @test _n6full_test_classification_rows(resume_db, method) ==
              _n6full_test_classification_rows(fixture_db, method)
        for n in 0:3
            resume_fingerprint_store = MoveGraphStore.open_store(resume_db)
            clean_fingerprint_store = MoveGraphStore.open_store(fixture_db)
            try
                resume_fingerprint = _n6full_db_serialization_fingerprint(
                    resume_fingerprint_store.db; n=n)
                clean_fingerprint = _n6full_db_serialization_fingerprint(
                    clean_fingerprint_store.db; n=n)
                @test resume_fingerprint == clean_fingerprint
            finally
                MoveGraphStore.close_store!(resume_fingerprint_store)
                MoveGraphStore.close_store!(clean_fingerprint_store)
            end
        end
        @test _n6full_verify_db(resume_db, resumed.counts; method_version=method)

        # Corrupting a committed per-n fingerprint is a malformed checkpoint,
        # not a signal to silently rebuild it.
        malformed_db = joinpath(fixture_dir, "malformed.sqlite3")
        @test_throws ErrorException _n6full_persist_stream(malformed_db,
            classification, bridge7; barrier8=prepared.barrier8,
            classical_map=prepared.classical_map, labels=prepared.labels,
            max_n=3, progress=false, method_version=method,
            snapshot_sha256="fixture-snapshot-sha256", stop_after_n=1)
        malformed_store = MoveGraphStore.open_store(malformed_db)
        try
            malformed_run = MoveGraphStore.get_classification_run(
                malformed_store, method)
            malformed_counts = deepcopy(malformed_run.counts)
            malformed_counts["stream_n_stats"][1]["emitted_fingerprint"] =
                repeat("0", 64)
            MoveGraphStore._dbexec(malformed_store.db, """
                UPDATE classification_runs SET counts_json=?
                WHERE method_version=?
            """, (JSON.json(malformed_counts), method))
        finally
            MoveGraphStore.close_store!(malformed_store)
        end
        @test_throws ErrorException _n6full_persist_stream(malformed_db,
            classification, bridge7; barrier8=prepared.barrier8,
            classical_map=prepared.classical_map, labels=prepared.labels,
            max_n=3, progress=false, method_version=method,
            snapshot_sha256="fixture-snapshot-sha256")

        db_sha, db_bytes = _n6full_sha256_file(fixture_db)
        fixture_bytes = read(fixture_db)
        @test (db_sha, db_bytes) ==
              (bytes2hex(SHA.sha256(fixture_bytes)), length(fixture_bytes))

        empty_path = joinpath(fixture_dir, "empty.bin")
        write(empty_path, UInt8[])
        @test _n6full_sha256_file(empty_path) ==
              (bytes2hex(SHA.sha256(UInt8[])), 0)
        boundary_bytes = UInt8[mod(i, 256) for i in
                               0:(N6FULL_HASH_CHUNK_BYTES - 1)]
        boundary_path = joinpath(fixture_dir, "boundary.bin")
        write(boundary_path, boundary_bytes)
        @test _n6full_sha256_file(boundary_path) ==
              (bytes2hex(SHA.sha256(boundary_bytes)),
               N6FULL_HASH_CHUNK_BYTES)
        boundary_plus_bytes = vcat(boundary_bytes, UInt8[0xa5])
        boundary_plus_path = joinpath(fixture_dir, "boundary-plus-one.bin")
        write(boundary_plus_path, boundary_plus_bytes)
        @test _n6full_sha256_file(boundary_plus_path) ==
              (bytes2hex(SHA.sha256(boundary_plus_bytes)),
               N6FULL_HASH_CHUNK_BYTES + 1)
        payload = _n6full_payload(prepared, result, "fixture-snapshot-sha256",
                                   1, db_sha, db_bytes, result.elapsed_seconds)
        payload_path = joinpath(fixture_dir, "fixture.json")
        _n6full_write_json(payload_path, payload)
        parsed_payload = JSON.parsefile(payload_path)
        @test parsed_payload["scope"] == "full canonical census"
        @test parsed_payload["result"]["counts"]["canonical_state_total"] == 181
        @test parsed_payload["input"]["barrier8_final_fingerprint_sha256"] ==
              barrier8.fingerprint
        @test parsed_payload["result"]["barrier8_summary"]["barrier"] == 8

        store = MoveGraphStore.open_store(fixture_db)
        try
            db = store.db
            @test _n6full_test_count(db, "SELECT COUNT(*) AS n FROM states") == 181
            @test _n6full_test_count(db,
                "SELECT COUNT(*) AS n FROM state_classifications") == 181
            @test _n6full_test_count(db,
                "SELECT SUM(planar=1) AS n FROM state_classifications") == 73
            @test _n6full_test_count(db,
                "SELECT SUM(planar=0) AS n FROM state_classifications") == 108
            @test _n6full_test_count(db,
                "SELECT SUM(classification_status='classified') AS n FROM state_classifications") == 73
            @test _n6full_test_count(db,
                "SELECT SUM(classification_status='nonclassical') AS n FROM state_classifications") == 108
            for table in ("transitions", "certificate_paths", "invariants",
                          "knotinfo_snapshots", "knotinfo_raw_records", "knot_types",
                          "move_types")
                @test _n6full_test_count(db, "SELECT COUNT(*) AS n FROM $table") == 0
            end
            @test _n6full_test_count(db,
                "SELECT COUNT(*) AS n FROM states WHERE payload IS NOT NULL") == 0
            columns = String.(getproperty.(MoveGraphStore._rows(
                db, "PRAGMA table_info(state_classifications)"), :name))
            @test columns == ["state_id", "planar", "knot_type",
                              "classification_status", "method_version"]
            run = MoveGraphStore.get_classification_run(store, method)
            @test run.max_n == 3
            @test run.barrier == 8
            @test run.census_fingerprint == result.census_fingerprint
            @test run.classification_fingerprint == result.classification_fingerprint
            @test run.counts["knotinfo_snapshot_sha256"] ==
                  "fixture-snapshot-sha256"
            @test run.conventions["scope"] ==
                  "full canonical signed one-component Gauss census; N<=6"

            first_row = first(MoveGraphStore._rows(db,
                "SELECT state_id, state_hash, serialization FROM states ORDER BY id"))
            state = MoveGraphStore.get_state(store, String(first_row.state_id))
            @test state !== nothing
            @test state.state_hash == String(first_row.state_hash)
            @test state.state_hash ==
                  bytes2hex(SHA.sha256(codeunits(state.serialization)))
            @test state.state_id ==
                  bytes2hex(SHA.sha256(codeunits(
                      "v$(state.canonical_version) $(state.serialization)")))
            row_classification = MoveGraphStore.get_state_classification(
                store, state; method_version=method)
            @test row_classification !== nothing
            @test row_classification.method_version == method

        finally
            MoveGraphStore.close_store!(store)
        end

        # The schema trigger is part of the Area-2 boundary: an Area-3 state
        # may exist, but it cannot receive this classification.
        trigger_store = MoveGraphStore.open_store(":memory:")
        try
            MoveGraphStore.register_classification_run!(trigger_store;
                method_version="trigger-fixture", canonical_version=1,
                max_n=0, barrier=8, census_fingerprint="trigger-census",
                classification_fingerprint="trigger-classification")
            area3 = CanonicalState("fixture-area3";
                area=AREA_3, representation="KnotGauss", complexity=0.0)
            MoveGraphStore.put_state!(trigger_store, area3)
            @test_throws Exception MoveGraphStore._dbexec(trigger_store.db, """
                INSERT INTO state_classifications(
                  state_id,planar,knot_type,classification_status,method_version)
                VALUES (?,?,?,?,?)
            """, (area3.state_id, 0, nothing, "nonclassical",
                  "trigger-fixture"))
        finally
            MoveGraphStore.close_store!(trigger_store)
        end
    finally
        GC.gc()
        rm(fixture_dir; recursive=true, force=true)
    end
end

@testset "N6 full production gate and exact constants" begin
    @test N6FULL_CANONICAL_STATE_COUNTS ==
          [1, 2, 14, 164, 3388, 96776, 3548876]
    @test sum(N6FULL_CANONICAL_STATE_COUNTS) == N6FULL_FULL_CANONICAL_TOTAL
    @test N6FULL_NONCLASSICAL_CANONICAL_TOTAL ==
          N6FULL_FULL_CANONICAL_TOTAL - N6FULL_CLASSICAL_CANONICAL_TOTAL
    @test_throws ArgumentError build_knot_n6_full_classification_db(progress=false)

    # This exercises the real --resume argument/production route without
    # starting the multi-million-state barrier chain.
    dispatch = _n6full_cli_dispatch([
        "--quiet", "--resume", "references/items/database/work.sqlite3"])
    @test dispatch.mode == :production_full
    @test dispatch.progress == false
    @test dispatch.resume_db_path ==
          "references/items/database/work.sqlite3"
    dispatch_equals = _n6full_cli_dispatch([
        "--resume=references/items/database/work.sqlite3"])
    @test dispatch_equals.mode == :production_full
    @test dispatch_equals.progress == true
    @test dispatch_equals.resume_db_path ==
          "references/items/database/work.sqlite3"

    gate_classification, gate_bridge7, gate_barrier8 =
        _n6full_test_production_gate_inputs()
    @test _n6full_validate_production_barrier8!(
        gate_classification, gate_bridge7, gate_barrier8) === nothing
    _, resolved_bridge7, resolved_barrier8 =
        _n6full_test_production_gate_inputs(
            barrier8_status="barrier8_resolved_all_seed_components",
            barrier8_unresolved=Int[])
    @test_throws ArgumentError _n6full_validate_production_barrier8!(
        gate_classification, resolved_bridge7, resolved_barrier8)
    _, fingerprint_bridge7, bad_fingerprint_barrier8 =
        _n6full_test_production_gate_inputs(
            barrier8_fingerprint=repeat("0", 64))
    @test_throws ArgumentError _n6full_validate_production_barrier8!(
        gate_classification, fingerprint_bridge7, bad_fingerprint_barrier8)
end
