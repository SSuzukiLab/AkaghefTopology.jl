using Test
using Topology

function _cut_method_b_fixture()
    source = VirtualLink()
    set_data!(source;
        rgauss=[[-1, 2, -2, 3, -3, 1, 0]],
        orientation=[1, 1, 1])
    source
end

@testset "Method-B uses open cut paths" begin
    source = _cut_method_b_fixture()
    candidate = only(movable(source, "PS";
        param="4", edge_occurrences=[(1, 1), (1, 1)]))
    condition = ps_global_condition(source, candidate)
    @test condition.valid
    @test condition.reason == :same_method_b_circuit

    witness = something(condition.witness)
    @test witness.paper_type == "PS-IV"
    @test witness.evaluated_state_side == :source_raw
    @test [gap.target_edge_id for gap in witness.gap_boundaries] == [1, 1]

    steps = Topology._method_b_steps(source)
    @test length(steps) == 3length(real_edge_table(source))
    boundary_limit = 6length(source.orientation)
    @test any(step -> step.from_node > boundary_limit || step.to_node > boundary_limit,
              steps)

    trace = witness.circuit_trace
    @test !isempty(trace)
    @test first(trace).from_node > boundary_limit
    @test last(trace).to_node > boundary_limit
    @test first(trace).from_node != last(trace).to_node
    @test all(trace[index].to_node == trace[index + 1].from_node
              for index in 1:length(trace)-1)
    @test any(step -> step.edge_id == witness.right_edge_id &&
                      step.lane == witness.right_lane, trace)
end
