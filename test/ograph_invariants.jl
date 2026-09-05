using Test
using Topology
using SQLite
using DBInterface

function _invariant_graph(code, signs)
    link = VirtualLink()
    set_data!(link; rgauss=[Int.(code)], orientation=Int.(signs))
    certify_closed_ograph(link)
end

const OGRAPH_INVARIANT_TEST_RELEASE = "2026-09-04T00:00:00Z"

@testset "deterministic o-graph pi1 presentation" begin
    s3 = _invariant_graph([1, -1], [1])
    presentation = ograph_pi1_presentation(s3)
    @test presentation.generator_count == 2
    @test presentation.relators == ((1, -2), (2, 1, -1))
    @test presentation == ograph_pi1_presentation(
        _invariant_graph([-1, 1], [1]))
    @test occursin("Koda/MST", presentation.evidence_json)

    l21 = ograph_pi1_presentation(certify_closed_ograph(lens_space_link(2, 1)))
    @test l21.generator_count == 4
    @test l21.relators == ((1, -2), (2, -3),
                           (4, 1, -1), (3, 2, -4))
    @test_throws ArgumentError OGraphPi1Presentation("synthetic", 1, [[0]])
    @test_throws ArgumentError OGraphPi1Presentation("synthetic", 1, [[2]])

    crossing_free = certify_closed_ograph(OGraphState(Int[], Int[]))
    @test crossing_free.witness.criteria == (true, true, true)
    @test crossing_free.witness.method_circuit_lengths == (0,)
    crossing_free_pi1 = ograph_pi1_presentation(crossing_free)
    @test crossing_free_pi1.generator_count == 1
    @test isempty(crossing_free_pi1.relators)
    crossing_free_h1 = ograph_h1(crossing_free_pi1)
    @test (crossing_free_h1.free_rank,
           crossing_free_h1.torsion_invariant_factors,
           crossing_free_h1.h1_finite_order) == (1, (), nothing)
end

@testset "exact H1 from signed-word presentations" begin
    trivial = ograph_h1(OGraphPi1Presentation("empty", 0, ()))
    @test trivial.free_rank == 0
    @test isempty(trivial.torsion_invariant_factors)
    @test trivial.torsion_subgroup_order == 1
    @test trivial.h1_finite_order == 1

    free = ograph_h1(OGraphPi1Presentation("free", 3, ()))
    @test free.free_rank == 3
    @test isempty(free.torsion_invariant_factors)
    @test free.torsion_subgroup_order == 1
    @test free.h1_finite_order === nothing

    cyclic = ograph_h1(OGraphPi1Presentation("cyclic", 1, [fill(1, 12)]))
    @test cyclic.free_rank == 0
    @test cyclic.torsion_invariant_factors == (BigInt(12),)
    @test cyclic.h1_finite_order == 12

    diagonal = ograph_h1(OGraphPi1Presentation(
        "diagonal", 2, [fill(1, 2), fill(2, 12)]))
    @test diagonal.torsion_invariant_factors == (BigInt(2), BigInt(12))
    @test diagonal.h1_finite_order == 24

    recombined = ograph_h1(OGraphPi1Presentation(
        "recombined", 2, [fill(1, 6), fill(2, 10)]))
    @test recombined.torsion_invariant_factors == (BigInt(2), BigInt(30))
    @test recombined.h1_finite_order == 60

    identified_generators = ograph_h1(OGraphPi1Presentation(
        "identified", 2, [[1, -2]]))
    @test identified_generators.free_rank == 1
    @test isempty(identified_generators.torsion_invariant_factors)

    unit = ograph_h1(OGraphPi1Presentation("unit", 1, [[1]]))
    @test unit.free_rank == 0
    @test isempty(unit.torsion_invariant_factors)
    @test unit.h1_finite_order == 1

    huge = big(typemax(Int64)) + 42
    huge_result = Topology._ograph_h1_from_matrix("huge", reshape([huge], 1, 1))
    @test huge_result.torsion_invariant_factors == (huge,)
    @test occursin(string(huge), huge_result.signature)
end

@testset "known closed-o-graph H1 regression" begin
    s3 = ograph_h1(ograph_pi1_presentation(
        _invariant_graph([1, -1], [1])))
    @test (s3.free_rank, s3.torsion_invariant_factors, s3.h1_finite_order) ==
          (0, (), BigInt(1))

    l21 = ograph_h1(ograph_pi1_presentation(
        certify_closed_ograph(lens_space_link(2, 1))))
    @test (l21.free_rank, l21.torsion_invariant_factors, l21.h1_finite_order) ==
          (0, (BigInt(2),), BigInt(2))

    s2xs1 = ograph_h1(ograph_pi1_presentation(_invariant_graph(
        [-1, -2, 3, 4, -4, -3, 2, 1], [1, 1, 1, 1])))
    @test (s2xs1.free_rank, s2xs1.torsion_invariant_factors,
           s2xs1.h1_finite_order) == (1, (), nothing)

    # This locally labelled Poincare-sphere regression has H1=0 just like S3.
    # No S3 identification is inferred; group triviality remains a separate,
    # proof-bearing enrichment status.
    poincare = ograph_h1(ograph_pi1_presentation(_invariant_graph(
        [4, -5, 3, -4, 2, -3, 1, -2, 5, -1], ones(Int, 5))))
    @test (poincare.free_rank, poincare.torsion_invariant_factors,
           poincare.h1_finite_order) == (0, (), BigInt(1))
end

@testset "append-only o-graph invariant store" begin
    store = open_store()
    migrate_ograph_invariant_store!(store)
    migrate_ograph_invariant_store!(store)
    @test Topology._ograph_invariant_existing_schema_version(store.db) == 1

    graphs = [
        _invariant_graph([1, -1], [1]),
        certify_closed_ograph(lens_space_link(2, 1)),
        _invariant_graph([-1, -2, 3, 4, -4, -3, 2, 1], [1, 1, 1, 1]),
    ]
    for graph in graphs
        put_closed_ograph!(store, graph)
        pair = put_ograph_invariants!(store, graph;
            created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE)
        @test put_ograph_pi1_presentation!(store, graph;
            created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE) == pair.presentation
        @test put_ograph_h1!(store, pair.presentation;
            created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE) == pair.h1
    end

    original = ograph_pi1_presentation(first(graphs))
    forged = OGraphPi1Presentation(
        original.presentation_id, original.state_id, original.presentation_hash,
        1, (), original.algorithm_id, original.algorithm_version,
        original.convention_id, original.serialization, original.evidence_json)
    @test_throws ArgumentError put_ograph_h1!(store, forged;
        created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE)
    summary = verify_ograph_invariant_store(store)
    @test summary == (
        ograph_schema_version=2,
        invariant_schema_version=1,
        states=3,
        presentations=3,
        h1_results=3,
        backend_runs=0,
        identifications=0,
    )
    @test_throws ArgumentError put_ograph_pi1_presentation!(store, first(graphs);
        created_at_utc="different-release")
    first_result = String(first(MoveGraphStore._rows(store.db,
        "SELECT h1_result_id FROM ograph_h1_results")).h1_result_id)
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db,
        "UPDATE ograph_h1_results SET free_rank=99 WHERE h1_result_id=?",
        (first_result,))
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
        INSERT INTO ograph_h1_results(
          h1_result_id,presentation_id,algorithm_id,algorithm_version,
          backend,backend_version,free_rank,
          torsion_invariant_factors_json,torsion_subgroup_order_decimal,
          h1_finite_order_decimal,h1_signature,result_fingerprint,
          evidence_json,created_at_utc)
        SELECT 'bad-decimal',presentation_id,algorithm_id,algorithm_version,
          backend,backend_version || '-bad-decimal',free_rank,
          torsion_invariant_factors_json,'1junk',h1_finite_order_decimal,
          h1_signature,
          'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb',
          evidence_json,created_at_utc
        FROM ograph_h1_results LIMIT 1
    """)
    @test_throws SQLite.SQLiteException DBInterface.execute(store.db, """
        INSERT INTO ograph_group_backend_runs(
          backend_run_id,presentation_id,backend,backend_version,operation,
          options_json,triviality_status,finiteness_status,pi1_order_decimal,
          external_group_id,result_json,evidence_json,transcript_sha256,created_at_utc)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, ("fake-run", first(ograph_pi1_presentation.(graphs)).presentation_id,
            "fake", "1", "triviality", "{}", "PROVEN_TRIVIAL",
            "PROVEN_FINITE", "1", missing, "{}", "{}", missing,
            OGRAPH_INVARIANT_TEST_RELEASE))
    @test verify_ograph_invariant_store(store).h1_results == 3
    close_store!(store)
end


@testset "o-graph H1 replay preserves stored backend provenance" begin
    store = open_store()
    migrate_ograph_invariant_store!(store)
    graph = _invariant_graph([1, -1], [1])
    put_closed_ograph!(store, graph)
    presentation = put_ograph_pi1_presentation!(store, graph;
        created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE)
    put_ograph_h1!(store, presentation;
        created_at_utc=OGRAPH_INVARIANT_TEST_RELEASE)

    old_result = Topology._ograph_h1_from_matrix(
        presentation.presentation_id,
        Topology._ograph_exponent_sum_matrix(presentation);
        backend_version="0.44.0")
    DBInterface.execute(store.db, """
        INSERT INTO ograph_h1_results(
          h1_result_id,presentation_id,algorithm_id,algorithm_version,
          backend,backend_version,free_rank,
          torsion_invariant_factors_json,torsion_subgroup_order_decimal,
          h1_finite_order_decimal,h1_signature,result_fingerprint,
          evidence_json,created_at_utc)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """, (old_result.result_id, old_result.presentation_id,
            old_result.algorithm_id, old_result.algorithm_version,
            old_result.backend, old_result.backend_version, old_result.free_rank,
            ograph_canonical_json(string.(old_result.torsion_invariant_factors)),
            string(old_result.torsion_subgroup_order),
            old_result.h1_finite_order === nothing ? nothing :
                string(old_result.h1_finite_order),
            old_result.signature, old_result.result_fingerprint,
            old_result.evidence_json, OGRAPH_INVARIANT_TEST_RELEASE))
    @test verify_ograph_invariant_store(store).h1_results == 2
    close_store!(store)
end

@testset "o-graph invariant future-version rejection" begin
    store = open_store()
    migrate_ograph_store!(store)
    DBInterface.execute(store.db, """
        CREATE TABLE ograph_invariant_schema_meta (
          key TEXT PRIMARY KEY, value TEXT NOT NULL)
    """)
    DBInterface.execute(store.db, """
        INSERT INTO ograph_invariant_schema_meta VALUES ('schema_version','2')
    """)
    @test_throws ArgumentError migrate_ograph_invariant_store!(store)
    @test isempty(MoveGraphStore._rows(store.db, """
        SELECT name FROM sqlite_master
        WHERE type='table' AND name='ograph_pi1_presentations'
    """))
    close_store!(store)
end
