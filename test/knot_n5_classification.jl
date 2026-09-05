using Test

include(joinpath(@__DIR__, "..", "src", "knot_n5_classification.jl"))

@testset "bounded N<=5 classical Gauss classification" begin
    result = build_knot_n5_classification()

    # These values are recomputed by the run from the rooted iterator and
    # Carter predicate; they are assertions of the finite contract, not input
    # used to skip or manufacture any states.
    @test result.state_counts == [1, 2, 10, 60, 536, 5752]
    @test length(result.states) == 6361
    @test Set(keys(result.canonical_knot_type_counts)) ==
          Set(KNOT_N5_STRICT_KNOT_TYPES)
    @test sum(values(result.canonical_knot_type_counts)) == 6361
    @test result.canonical_knot_type_counts == Dict(
        "0_1" => 5973, "3_1^+" => 174, "3_1^-" => 174, "4_1" => 36,
        "5_1^+" => 1, "5_1^-" => 1, "5_2^+" => 1, "5_2^-" => 1)
    @test Set(keys(result.raw_rooted_classical_counts)) ==
          Set(KNOT_N5_STRICT_KNOT_TYPES)
    @test result.raw_rooted_total == 995573
    @test result.raw_rooted_classical_total == 62005
    @test sum(values(result.raw_rooted_classical_counts)) == 62005
    @test result.raw_rooted_classical_counts == Dict(
        "0_1" => 58213, "3_1^+" => 1710, "3_1^-" => 1710, "4_1" => 348,
        "5_1^+" => 2, "5_1^-" => 2, "5_2^+" => 10, "5_2^-" => 10)
    @test result.raw_rooted_classical_unresolved_count == 0
    @test length(result.component_roots) == 8
    @test length(result.component_invariants) == 8
    @test isempty(result.unresolved_components)
    @test result.path_replay_failures == 0
    @test result.path_replay_count > 0
    @test isempty(result.edge_invariant_failures)
    @test result.edge_invariant_checks == length(result.edges)

    # All retained R1/R2 edges point from the higher-crossing state to the
    # lower-crossing state.  Their undirected registration therefore contains
    # the inverse insertion edges too, including insertions in both directions
    # whenever both endpoints lie in the barrier.
    @test all(edge -> edge.move_type == :R3 ||
                    crossing_count(result.states[edge.source]) >
                    crossing_count(result.states[edge.target]), result.edges)
    @test result.direct_transition_count == 23830
    @test result.duplicate_edge_transition_count == 3176
    @test length(result.edges) == 20654
    @test result.fingerprint ==
          "35c4baa4d1800266b5334e64fc85476a49f9b600347e551176ff6563ee861047"

    known = result.known_representatives
    @test length(known) == 10
    @test [representative.label for representative in known] == [
        "0_1", "0_1_mirror", "3_1", "3_1_mirror", "4_1", "4_1_mirror",
        "5_1", "5_1_mirror", "5_2", "5_2_mirror"]
    @test all(representative -> representative.state_index >= 1 &&
                              representative.component >= 1,
              known)
    @test known[findfirst(x -> x.label == "3_1", known)].component !=
          known[findfirst(x -> x.label == "3_1_mirror", known)].component
    @test known[findfirst(x -> x.label == "4_1", known)].component ==
          known[findfirst(x -> x.label == "4_1_mirror", known)].component

    # The strict quotient keeps traversal reversal/mirror out of the
    # canonicalizer.  (Some individual knots, such as 4_1, can nevertheless
    # be connected to their mirror by the move graph.)
    trefoil = known[findfirst(x -> x.label == "3_1", known)].diagram
    trefoil_mirror = known[findfirst(x -> x.label == "3_1_mirror", known)].diagram
    @test canonicalize(trefoil) != canonicalize(trefoil_mirror)
    @test canonicalize(trefoil) == result.states[known[findfirst(
        x -> x.label == "3_1", known)].state_index]

    first_edge = first(result.edges)
    raw_target = apply_reidemeister(result.states[first_edge.source],
                                    first_edge.candidate)
    @test verify_reidemeister_transition(result.states[first_edge.source],
                                         raw_target, first_edge.candidate)
    @test canonicalize(raw_target) == result.states[first_edge.target]
end
