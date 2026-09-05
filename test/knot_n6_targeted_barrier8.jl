using Test

include(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_barrier8.jl"))

@testset "targeted barrier-8 empty fixture" begin
    classification = build_knot_n6_classification(
        3; barrier=3, progress=false)
    workspace_ref = Ref{Any}(nothing)
    bridge7 = run_knot_n6_targeted_bridge(
        classification; progress=false,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    result = run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; progress=false)
    summary = knot_n6_targeted_barrier8_summary(result)

    @test result.barrier == 8
    @test isempty(result.seed_components)
    @test result.seed_state_count == 0
    @test result.imported_n7_state_count == 0
    @test result.n7_state_count == 0
    @test result.n8_state_count == 0
    @test result.candidate_count == 0
    @test result.adopted_transition_count == 0
    @test result.duplicate_transition_count == 0
    @test result.imported_graph_edge_count == 0
    @test result.reused_n6_candidate_count == 0
    @test result.excluded_over_barrier_count == 0
    @test result.n7_processed_count == 0
    @test result.n8_processed_count == 0
    @test result.max_frontier == 0
    @test result.transition_verification_failures == 0
    @test result.n6_membership_failures == 0
    @test result.classicality_checks == 0
    @test result.classicality_failures == 0
    @test isempty(result.connected_known_labels)
    @test isempty(result.resolved_components)
    @test isempty(result.unresolved_components)
    @test result.status == "barrier8_resolved_all_seed_components"
    @test length(result.fingerprint) == 64
    @test summary["barrier"] == 8
    @test summary["n7_state_count"] == 0
    @test summary["n8_state_count"] == 0
    @test !haskey(summary, "path")
    @test !haskey(summary, "parent")
    @test_throws ArgumentError run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; barrier=7, progress=false)
end
@testset "targeted barrier-8 canonical identities" begin
    uf = _N6B8UnionFind()
    adjacency = Vector{Vector{Int}}()
    n7_index = Dict{KnotGaussDiagram,Int}()
    n7_states = KnotGaussDiagram[]
    n7_nodes = Int[]
    n7_local_by_node = Dict{Int,Int}()
    diagram = KnotGaussDiagram([1, -1], [1])
    first_node, first_new = _n6b8_ensure_n7!(
        diagram, n7_index, n7_states, n7_nodes, n7_local_by_node,
        uf, adjacency)
    second_node, second_new = _n6b8_ensure_n7!(
        canonicalize(diagram), n7_index, n7_states, n7_nodes, n7_local_by_node,
        uf, adjacency)
    @test first_new
    @test !second_new
    @test first_node == second_node
    @test length(n7_states) == 1

    n8_index = Dict{KnotGaussDiagram,Int}()
    n8_states = KnotGaussDiagram[]
    n8_nodes = Int[]
    n8_local_by_node = Dict{Int,Int}()
    n8_first, n8_is_new = _n6b8_ensure_n8!(
        diagram, n8_index, n8_states, n8_nodes, n8_local_by_node,
        uf, adjacency)
    n8_second, n8_second_new = _n6b8_ensure_n8!(
        canonicalize(diagram), n8_index, n8_states, n8_nodes,
        n8_local_by_node, uf, adjacency)
    @test n8_is_new
    @test !n8_second_new
    @test n8_first == n8_second
    @test length(n8_states) == 1
end
