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
    # The compatibility projection drops two loops and merges the duplicated
    # 3--4 edge retained by the authoritative C2 multigraph.
    @test ne(graph_a) == 3
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

@testset "BP Figure 1.2 C2 multigraph" begin
    one = c2_graph([[1, -1]], [1])
    @test one.node_count == 4
    @test one.edges == [(1, 1), (2, 2), (3, 4),
                        (2, 4), (3, 1), (4, 3)]
    @test one.origins == [(1, 1, 1), (1, 1, 2), (1, 1, 3),
                          (1, 2, 1), (1, 2, 2), (1, 2, 3)]

    # One cyclic word exercising +/+ (same), +/- (mixed), -/- (same),
    # and -/+ (mixed), including the wrap-around pair.
    four_pairs = c2_graph([[1, 2, -1, -2]], [1, 1])
    @test four_pairs.node_count == 8
    @test four_pairs.edges == [
        (1, 5), (2, 8), (3, 7),
        (5, 1), (6, 2), (7, 4),
        (2, 6), (3, 5), (4, 8),
        (6, 4), (7, 1), (8, 3),
    ]
    @test four_pairs.origins == [
        (1, occurrence, lane)
        for occurrence in 1:4 for lane in 1:3
    ]

    degree = zeros(Int, four_pairs.node_count)
    for (left, right) in four_pairs.edges
        degree[left] += 1
        degree[right] += 1
    end
    @test all(==(3), degree)
    @test length(four_pairs.edges) == 6length([1, 1])

    multi = c2_graph([[1, -1], [2, -2]], [1, -1])
    @test multi.node_count == 8
    @test length(multi.edges) == 12
    @test last(multi.origins) == (2, 2, 3)

    state = OGraphState([1, 2, -1, -2], [1, 1])
    from_state = c2_graph(state)
    @test from_state.node_count == four_pairs.node_count
    @test from_state.edges == four_pairs.edges
    @test from_state.origins == four_pairs.origins

    weighted = _test_ograph_link([1, -1], [1])
    set_weight!(weighted, [[-2, 1]])
    saved_weights = deepcopy(weighted.weights)
    @test c2_graph(weighted).edges == one.edges
    @test weighted.is_weighted
    @test weighted.weights == saved_weights

    crossing_free = VirtualLink()
    set_data!(crossing_free; rgauss=[Int[]], orientation=Int[])
    @test nv(ograph_graph_a(crossing_free)) == 0
    @test_throws ArgumentError c2_graph(crossing_free)

    virtual = VirtualLink()
    set_data!(virtual; gauss=[[1, -1]], orientation=[0])
    @test_throws ArgumentError c2_graph(virtual)
    @test is_closed(virtual) == (false, (true, false, false))
    cut = VirtualLink()
    set_data!(cut; rgauss=[[1, -1]], orientation=[1], cut=[true])
    @test_throws ArgumentError c2_graph(cut)
    @test is_closed(cut) == (false, (true, false, false))

    @test_throws ArgumentError c2_graph(Vector{Vector{Int}}(), Int[])
    @test_throws ArgumentError c2_graph([Int[]], [1])
    @test_throws ArgumentError c2_graph([[1, -1]], [0])
    @test_throws ArgumentError c2_graph([[1, 1]], [1])
    @test_throws ArgumentError c2_graph([[1, -1, 2, -2]], [1])
    @test_throws DimensionMismatch C2Graph(4, [(1, 2)], NTuple{3,Int}[])
    @test_throws ArgumentError C2Graph(4, [(1, 5)], [(1, 1, 1)])
end
