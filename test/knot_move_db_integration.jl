using Test
using JSON
using SQLite
using Topology

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const KNOTINFO = joinpath(REPO_ROOT, "threads", "PJ260830", "data", "knotinfo-n6.json")
const GENERATED = joinpath(REPO_ROOT, "threads", "PJ260830", "data", "knot-n6-move-db.json")
const GENERATED_DB = joinpath(REPO_ROOT, "references", "items", "database", "knotinfo", "knot-n6-moves.sqlite3")

@testset "Knot move DB package integration" begin
    @test KnotGaussDiagram isa DataType
    @test ReidemeisterCandidate isa DataType
    @test any(m -> occursin("VirtualLink", string(m.sig)), collect(methods(movable)))
    @test any(m -> occursin("KnotGaussDiagram", string(m.sig)), collect(methods(movable)))

    diag_dir = mktempdir()
    store = MoveGraphStore.open_store(joinpath(diag_dir, "diagnostic.sqlite3"))
    empty = KnotGaussDiagram(Int[], Int[])
    insert_candidate = first(movable(empty, "R1"; param=:insert))
    one = apply_reidemeister(empty, insert_candidate)
    e0 = persist_reidemeister_transition!(store, empty, insert_candidate)
    delete_candidate = first(movable(one, "R1"; param=:delete))
    e1 = persist_reidemeister_transition!(store, one, delete_candidate)
    @test e0 isa Integer && e1 isa Integer
    s0, s1 = knot_canonical_state(empty), knot_canonical_state(one)
    @test MoveGraphStore.bfs_certificate_path(store, s0.state_id, s1.state_id) == [s0.state_id, s1.state_id]
    @test length(MoveGraphStore.neighbors(store, s1.state_id)) == 1
    @test_throws ArgumentError MoveGraphStore.insert_verified_transition!(store, s0, s1, "R1"; verified=true)

    bounded = expand_reidemeister_neighbors!(store, empty; store_cap=1, max_crossings=1)
    @test bounded.persisted <= 1
    @test bounded.capped || bounded.persisted < 1

    @test MoveGraphStore.import_knotinfo!(store, KNOTINFO) == 8
    @test length(MoveGraphStore.knotinfo_records(store)) == 8
    @test isempty(MoveGraphStore._rows(store.db, "SELECT name FROM knot_types WHERE state_id IS NOT NULL"))

    MoveGraphStore.close_store!(store)

    if isfile(GENERATED)
        payload = JSON.parsefile(GENERATED)
        @test payload["schema_version"] == 1
        @test payload["counts"]["seed_representations"] == 281
        @test payload["counts"]["unique_seed_states"] == 215
        @test payload["counts"]["seed_dedupes"] == 66
        @test payload["counts"]["candidate_count"] > 12
        @test payload["counts"]["persisted_edge_count"] > 0
        @test payload["counts"]["rejection_count"] == 0
        @test payload["counts"]["error_counts"] == Dict{String,Any}()
        @test haskey(payload["database"], "sha256")
        @test haskey(payload["runtime"], "command")
        @test isfile(GENERATED_DB)
        db = SQLite.DB(GENERATED_DB)
        try
            state_count = Int(first(DBInterface.execute(db, "SELECT count(*) FROM states"))[1])
            edge_count = Int(first(DBInterface.execute(db, "SELECT count(*) FROM transitions"))[1])
            unlinked_knotinfo = Int(first(DBInterface.execute(db, "SELECT count(*) FROM knot_types WHERE state_id IS NULL"))[1])
            @test state_count == payload["counts"]["state_count"]
            @test edge_count == payload["counts"]["persisted_edge_count"]
            @test unlinked_knotinfo == 8
        finally
            SQLite.close(db)
        end
    end
end
