using Test
using Topology
using Graphs

function _test_ograph_link(code, signs)
    link = VirtualLink()
    set_data!(link; rgauss=[Int.(code)], orientation=Int.(signs))
end

@testset "immutable o-graph presentation state" begin
    code = [-2, -1, 1, 2]
    signs = [1, -1]
    state = OGraphState(code, signs)
    code[1] = 99
    signs[1] = -1
    @test state.gauss == (-2, -1, 1, 2)
    @test state.vertex_signs == (1, -1)
    @test !Base.ismutabletype(typeof(state))

    rotated = OGraphState([1, 2, -2, -1], [1, -1])
    relabelled = OGraphState([-1, -2, 2, 1], [-1, 1])
    @test state != rotated
    @test canonicalize_ograph(state) == canonicalize_ograph(rotated)
    @test canonicalize_ograph(state) == canonicalize_ograph(relabelled)
    @test canonicalize_ograph(canonicalize_ograph(state)) == canonicalize_ograph(state)
    @test ograph_canonical_hash(state) == ograph_canonical_hash(rotated) ==
          ograph_canonical_hash(relabelled)
    @test startswith(ograph_canonical_serialization(state),
        "OGraphState;schema=akaghef.ograph/1;canonicalization=ograph-presentation/1;")

    changed_sign = OGraphState([-2, -1, 1, 2], [-1, -1])
    @test ograph_canonical_hash(state) != ograph_canonical_hash(changed_sign)
    @test ograph_canonical_hash(state) != ograph_canonical_hash(
        OGraphState(state.gauss, state.vertex_signs; convention_id="other-convention/1"))

    @test_throws ArgumentError OGraphState([1, 1], [1])
    @test_throws ArgumentError OGraphState([1, -1], [0])
    @test_throws ArgumentError OGraphState([1, -1, 2], [1])
    @test_throws ArgumentError OGraphState([1, -1], [1]; convention_id="")

    weighted = _test_ograph_link([1, -1], [1])
    set_weight!(weighted, [[-2, 1]])
    @test_throws ArgumentError ograph_state(weighted)
    multi = VirtualLink()
    set_data!(multi; rgauss=[[1, -1], [2, -2]], orientation=[1, 1])
    @test_throws ArgumentError ograph_state(multi)
end

@testset "closed normal o-graph certification" begin
    graph_a = ograph_graph_a(_test_ograph_link([1, -1], [1]))
    @test graph_a isa Graphs.SimpleGraph
    @test nv(graph_a) == 4
    @test ne(graph_a) == 4
    @test length(connected_components(graph_a)) == 1

    s3 = certify_closed_ograph(_test_ograph_link([1, -1], [1]))
    @test s3.witness.criteria == (true, true, true)
    @test s3.witness.true_vertices == 1
    @test s3.witness.c3_method_circuits == 2
    @test s3.witness.method_circuit_lengths == (5, 1)

    lens = certify_closed_ograph(lens_space_link(2, 1))
    @test lens.witness.criteria == (true, true, true)
    @test lens.witness.true_vertices == 2
    @test lens.witness.c3_method_circuits == 3
    @test lens.witness.method_circuit_lengths == (8, 3, 1)

    s2xs1 = certify_closed_ograph(_test_ograph_link(
        [-1, -2, 3, 4, -4, -3, 2, 1], [1, 1, 1, 1]))
    @test s2xs1.witness.criteria == (true, true, true)
    @test s2xs1.witness.true_vertices == 4
    @test s2xs1.witness.c3_method_circuits == 5
    @test s2xs1.witness.method_circuit_lengths == (8, 6, 1, 8, 1)

    c3_failure = OGraphState([1, 2, -1, -2], [1, 1])
    report = closedness_witness(c3_failure)
    @test report.criteria == (true, true, false)
    @test report.c3_method_circuits == 1
    @test_throws ArgumentError certify_closed_ograph(c3_failure)

    stale = ClosednessWitness(
        "wrong-hash", s3.witness.verifier_version, s3.witness.criteria,
        s3.witness.c1_oriented_circuits, s3.witness.c2_trivalent_components,
        s3.witness.c3_method_circuits, s3.witness.true_vertices,
        s3.witness.method_circuit_lengths)
    @test_throws ArgumentError ClosedNormalOGraph(s3.state, stale)
end
