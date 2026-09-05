using Test
using Topology

@testset "N6 package API integration" begin
    direct = direct_classical_gauss_enumeration(2)
    @test direct isa ClassicalGaussEnumerationResult
    @test direct.generated_decorations == 32
    @test length(direct.raw_keys) == 32
    key = first(direct.raw_keys)
    code, orientation = classical_gauss_key_to_code(key, 2)
    @test classical_gauss_key_to_diagram(key, 2).gauss == Tuple(code)
    @test classical_gauss_key_to_diagram(key, 2).orientation == Tuple(orientation)

    streamed = stream_canonical_knot_gauss(2)
    @test streamed isa KnotGaussCanonicalStreamStats
    @test streamed.raw_total == 48
    @test streamed.raw_seen == 48
    @test streamed.canonical_emitted == 14
    @test streamed.stop_reason === :complete
    @test enumerate_canonical_knot_gauss(1).canonical_emitted == 2

    classification = build_knot_n6_classification(
        3; barrier=3, progress=false)
    @test classification isa KnotN6Classification
    @test classification.state_counts == [1, 2, 10, 60]
    @test length(classification.states) == 73
    @test isempty(classification.unresolved_components)
    @test knot_n6_summary(classification)["total_states"] == 73

    bridge7 = run_knot_n6_targeted_bridge(
        classification; progress=false)
    @test bridge7 isa KnotN6TargetedBridgeResult
    @test isempty(bridge7.seed_components)
    @test bridge7.n7_state_count == 0
    @test knot_n6_targeted_bridge_summary(bridge7)["n7_state_count"] == 0

    barrier8 = build_knot_n6_targeted_barrier8(
        classification; progress=false)
    @test barrier8 isa KnotN6TargetedBarrier8Result
    @test isempty(barrier8.seed_components)
    @test barrier8.n7_state_count == 0
    @test barrier8.n8_state_count == 0
    @test barrier8.status == "barrier8_resolved_all_seed_components"
    @test knot_n6_targeted_barrier8_summary(
        barrier8)["n8_state_count"] == 0
end
