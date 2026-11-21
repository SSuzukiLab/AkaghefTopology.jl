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

# N1 [BP97, §1.4], N2 [BP97, §1.4], C1 [BP97, §1.4], C2 [BP97, §1.4], C3 [BP97, §1.4]
@testset "N1 N2 C1 C2 C3 of a closed normal o-graph" begin
    function _argument_error(body)
        try
            body()
            nothing
        catch error
            error
        end
    end

    empty_state = OGraphState(Int[], Int[])
    empty_witness = closedness_witness(empty_state)
    @test empty_witness.n1 == false
    @test empty_witness.n2_enforcement == "OGraphState"
    @test empty_witness.c1_enforcement == "ograph_state"
    # MST numeric contract read by the move verifier.  Not an N1 certificate.
    @test empty_witness.criteria == (true, true, true)
    @test empty_witness.c1_oriented_circuits == 1
    @test empty_witness.method_circuit_lengths == (0,)
    n1_error = _argument_error(() -> certify_closed_ograph(empty_state))
    @test n1_error isa ArgumentError
    @test occursin("N1", n1_error.msg)
    @test occursin("[BP97, §1.4]", n1_error.msg)
    direct = _argument_error(() -> ClosedNormalOGraph(empty_state, empty_witness))
    @test direct isa ArgumentError
    @test occursin("N1", direct.msg)
    @test occursin("[BP97, §1.4]", direct.msg)

    over_pair = _argument_error(() -> OGraphState([1, 1], [1]))
    @test over_pair isa ArgumentError
    @test occursin("N1", over_pair.msg)
    @test occursin("[BP97, §1.4]", over_pair.msg)
    @test !occursin("N2", over_pair.msg)

    odd_length = _argument_error(() -> OGraphState([1, -1, 2], [1]))
    @test odd_length isa ArgumentError
    @test occursin("N1", odd_length.msg)
    @test occursin("[BP97, §1.4]", odd_length.msg)

    missing_vertex = _argument_error(() -> OGraphState([1, 2, -2, 2], [1, 1]))
    @test missing_vertex isa ArgumentError
    @test occursin("N1", missing_vertex.msg)
    @test occursin("[BP97, §1.4]", missing_vertex.msg)

    # N2 has no representable failing word.  The site string is the check.
    hopf = closedness_witness(OGraphState([1, -1], [1]))
    @test hopf.n1
    @test hopf.n2_enforcement == Topology.N2_ENFORCEMENT_SITE
    @test hopf.c1_enforcement == Topology.C1_ENFORCEMENT_SITE

    multi = VirtualLink()
    set_data!(multi; rgauss=[[1, -1], [2, -2]], orientation=[1, 1])
    c1_error = _argument_error(() -> ograph_state(multi))
    @test c1_error isa ArgumentError
    @test occursin("C1", c1_error.msg)
    @test occursin("[BP97, §1.4]", c1_error.msg)

    crossing_free_link = VirtualLink()
    set_data!(crossing_free_link; rgauss=[Int[]], orientation=Int[])
    @test ograph_state(crossing_free_link) == empty_state
    @test is_closed(crossing_free_link) == (false, (true, true, false))

    c2_failure = OGraphState([1, 2, -1, -2], [1, -1])
    c2_report = closedness_witness(c2_failure)
    @test c2_report.n1
    @test c2_report.criteria == (true, false, true)
    @test c2_report.c2_trivalent_components == 2
    @test c2_report.c3_method_circuits == 3
    @test c2_report.true_vertices == 2
    c2_error = _argument_error(() -> certify_closed_ograph(c2_failure))
    @test c2_error isa ArgumentError
    @test occursin("fails C2 [BP97, §1.4]", c2_error.msg)
    @test !occursin("C3", c2_error.msg)

    c3_failure = OGraphState([1, 2, -1, -2], [1, 1])
    c3_report = closedness_witness(c3_failure)
    @test c3_report.criteria == (true, true, false)
    @test c3_report.c3_method_circuits == 1
    c3_error = _argument_error(() -> certify_closed_ograph(c3_failure))
    @test c3_error isa ArgumentError
    @test occursin("fails C3 [BP97, §1.4]", c3_error.msg)
    @test !occursin("C2", c3_error.msg)

    both = OGraphState([1, -1, -2, 2], [1, -1])
    both_report = closedness_witness(both)
    @test both_report.criteria == (true, false, false)
    both_error = _argument_error(() -> certify_closed_ograph(both))
    @test both_error isa ArgumentError
    @test occursin("fails C2, C3 [BP97, §1.4]", both_error.msg)

    empty_c2 = _argument_error(() -> c2_graph(Vector{Vector{Int}}(), Int[]))
    @test empty_c2 isa ArgumentError
    @test occursin("N1", empty_c2.msg)
    @test occursin("[BP97, §1.4]", empty_c2.msg)
    same_sign_c2 = _argument_error(() -> c2_graph([[1, 1]], [1]))
    @test same_sign_c2 isa ArgumentError
    @test occursin("N1", same_sign_c2.msg)
    @test occursin("[BP97, §1.4]", same_sign_c2.msg)

    isdefined(@__MODULE__, :vl_example) ||
        include(joinpath(@__DIR__, "../examples/VLExample.jl"))
    for name in ("s3_abalone", "lens_2_1", "lens_3_1_mst2_example23",
                 "poincare", "quaternion", "s2xs1", "s2xs1_2vertex")
        certified = certify_closed_ograph(vl_example(name))
        @test certified isa ClosedNormalOGraph
        @test certified.witness.n1
        @test certified.witness.criteria == (true, true, true)
        @test certified.witness.n2_enforcement == "OGraphState"
        @test certified.witness.c1_enforcement == "ograph_state"
    end

    circuit = Topology.mst_crossing_free_circuit()
    @test circuit isa Topology.MSTCrossingFreeCircuit
    @test circuit.state == empty_state
    @test circuit.witness == empty_witness
    nonempty = OGraphState([1, -1], [1])
    refused = _argument_error(() -> Topology.MSTCrossingFreeCircuit(
        nonempty, closedness_witness(nonempty)))
    @test refused isa ArgumentError
    @test occursin("no vertices", refused.msg)
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
