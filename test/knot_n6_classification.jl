using Test

include(joinpath(@__DIR__, "..", "src", "knot_n6_classification.jl"))

@testset "direct classical generator contract" begin
    expected_raw = [1, 4, 32, 336]
    expected_canonical = [1, 2, 10, 60]
    for n in 0:3
        census = direct_classical_gauss_enumeration(n)
        @test length(census.raw_keys) == expected_raw[n + 1]
        canonical = Set(canonical_serialization(
            classical_gauss_key_to_diagram(key, n)) for key in census.raw_keys)
        @test length(canonical) == expected_canonical[n + 1]
    end
end

@testset "barrier-3 N6 classifier fixture" begin
    result = build_knot_n6_classification(3; barrier=3, progress=false)
    @test result.state_counts == [1, 2, 10, 60]
    @test result.raw_rooted_classical_counts_by_n == [1, 4, 32, 336]
    @test result.raw_rooted_total == 1013
    @test result.raw_rooted_classical_total == 373
    @test length(result.states) == 73
    @test length(result.component_roots) == 3
    @test isempty(result.unresolved_components)
    @test knot_n6_summary(result)["unresolved_state_count"] == 0
    @test result.raw_rooted_classical_unresolved_count == 0
    @test result.spanning_tree_replay_failures == 0
    @test result.spanning_tree_edge_replay_count == 70
    @test result.edge_invariant_sample_count == 76
    @test isempty(result.edge_invariant_failures)
    @test result.nonclassical_target_count == 0
    @test result.component_type_labels == Dict(1 => "0_1", 2 => "3_1^+", 3 => "3_1^-")
    @test result.canonical_knot_type_counts ==
          Dict("0_1" => 71, "3_1^+" => 1, "3_1^-" => 1)
    @test result.raw_rooted_classical_type_counts ==
          Dict("0_1" => 369, "3_1^+" => 2, "3_1^-" => 2)
    @test sum(values(result.canonical_knot_type_counts)) == 73
    @test sum(values(result.raw_rooted_classical_type_counts)) == 373
    @test result.fingerprint ==
          "65d97988cb17fce54308a77eb61829d0d49a017f28ac181a81816b6215408be0"
    @test all(edge -> begin
        raw = apply_reidemeister(result.states[edge.source], edge.candidate)
        verify_reidemeister_transition(result.states[edge.source], raw, edge.candidate) &&
            canonicalize(raw) == result.states[edge.target]
    end, result.edges)
end
