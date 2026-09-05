using Test

include(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_bridge.jl"))

@testset "targeted bridge input and empty fixture" begin
    classification = build_knot_n6_classification(3; barrier=3, progress=false)
    result = run_knot_n6_targeted_bridge(classification; progress=false)
    summary = knot_n6_targeted_bridge_summary(result)
    @test result.barrier == 7
    @test isempty(result.seed_components)
    @test result.seed_state_count == 0
    @test result.n7_state_count == 0
    @test result.candidate_count == 0
    @test result.adopted_transition_count == 0
    @test result.duplicate_transition_count == 0
    @test result.reused_n6_candidate_count == 0
    @test result.excluded_over_barrier_count == 0
    @test result.n7_processed_count == 0
    @test result.max_frontier == 0
    @test result.transition_verification_failures == 0
    @test result.n6_membership_failures == 0
    @test result.classicality_checks == 0
    @test result.classicality_failures == 0
    @test isempty(result.resolved_components)
    @test isempty(result.unresolved_components)
    @test summary["barrier"] == 7
    @test summary["n7_state_count"] == 0
    @test !haskey(summary, "parent")
    @test !haskey(summary, "path")
    @test length(result.fingerprint) == 64
    @test_throws ArgumentError run_knot_n6_targeted_bridge(
        classification; barrier=8, progress=false)
end

@testset "canonical N7 identity is deduplicated" begin
    uf = _N6TBUnionFind()
    adjacency = Vector{Vector{Int}}()
    index = Dict{KnotGaussDiagram,Int}()
    states = KnotGaussDiagram[]
    nodes = Int[]
    local_by_node = Dict{Int,Int}()
    parents = Int[]
    distances = Int[]
    diagram = KnotGaussDiagram([1, -1], [1])
    first_node, first_new = _n6tb_ensure_n7!(
        diagram, index, states, nodes, local_by_node, parents, distances,
        uf, adjacency)
    second_node, second_new = _n6tb_ensure_n7!(
        canonicalize(diagram), index, states, nodes, local_by_node, parents,
        distances, uf, adjacency)
    @test first_new
    @test !second_new
    @test first_node == second_node
    @test length(states) == 1
    @test local_by_node[first_node] == 1
end

