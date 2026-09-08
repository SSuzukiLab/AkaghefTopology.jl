using Test
using Topology

@testset "PS preserves empty interval identity through full deletion" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[0],[0]],orientation=Int[])
    for parameter in 1:4, locations in ([(1,1),(1,1)],[(1,1),(2,1)])
        candidate=only(movable(source,"PS";param=string(parameter),edge_occurrences=locations))
        forward=move_boundary_transport(source,candidate)
        inverse=only(movable(forward.target,"PS";param=string(parameter),vertices=[1,2]))
        backward=move_boundary_transport(forward.target,inverse)
        @test backward.target.gauss==source.gauss
        @test isempty(backward.target.orientation)
        @test length(real_edge_table(backward.target))==2
        # Same lane on the same interval is connected; lanes from disjoint
        # intervals are not. This checks the actual global query separately
        # from locally applicable fragment substitution.
        expected=parameter==4 && locations[1]==locations[2]
        @test ps_global_condition(source,candidate).valid==expected
        @test ps_global_condition(forward.target,inverse).valid==expected
        back=Dict(row.source=>row.target for row in backward.endpoints)
        @test all(back[row.target]==row.source for row in forward.endpoints)
    end
end

function _multi_cut_method_b_fixture()
    source = VirtualLink()
    set_data!(source;
        rgauss=[[-1, 2, -2, 3, -3, 1, 0],
                [-4, 5, -5, 6, -6, 4, 0]],
        orientation=ones(Int, 6))
    source
end

@testset "Method-B authenticates a PS patch across open components" begin
    source = _multi_cut_method_b_fixture()
    candidate = only(movable(source, "PS";
        param="4", edge_occurrences=[(1, 1), (1, 1)]))

    reason, frozen, code, orientation =
        Topology._check_open_ps_candidate_source(source, candidate)
    @test reason == :ok
    @test length(code) == 2
    @test all(last(word) == 0 for word in code)
    @test orientation == ones(Int, 6)
    @test frozen.locations == [(1, 1), (1, 1)]

    condition = ps_global_condition(source, candidate)
    @test condition.valid
    @test condition.reason == :same_method_b_circuit
    witness = something(condition.witness)
    @test [gap.source_component for gap in witness.gap_boundaries] == [1, 1]
    @test [gap.target_edge_id for gap in witness.gap_boundaries] == [1, 1]

    # The second open strand remains in the lane graph, but cannot be used to
    # authenticate a selected subpath on the first strand.
    steps = Topology._method_b_steps(source)
    @test any(step -> step.component == 2, steps)
    @test all(any(step -> step == selected, steps)
              for selected in witness.circuit_trace)
end

@testset "Method-B preserves an additional marker-only cut interval" begin
    source = VirtualLink()
    set_data!(source;
        rgauss=[[-1, 2, -2, 3, -3, 1, 0], [0]],
        orientation=ones(Int, 3))
    candidate = only(movable(source, "PS";
        param="4", edge_occurrences=[(1, 1), (1, 1)]))
    reason, _, code, _ =
        Topology._check_open_ps_candidate_source(source, candidate)
    @test reason == :ok
    @test code[2] == [0]
    @test ps_global_condition(source, candidate).valid
end

@testset "Open PS source authentication rejects marker wrapping" begin
    source = _multi_cut_method_b_fixture()
    data = make_move_data()
    candidate = MoveCandidate("PS", "4", [1, 2], NTuple{2,Int}[], :inverse,
        data["PS-R"][4], only(data["PS-L"]), [(1, 5), (1, 5)])
    reason, _, _, _ = Topology._check_open_ps_candidate_source(source, candidate)
    @test reason == :source_mismatch
    @test !ps_global_condition(source, candidate).valid
end
