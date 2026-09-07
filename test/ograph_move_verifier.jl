using Test
using Topology

function d2_inverse_fixture()
    source = VirtualLink()
    set_data!(source;
        rgauss=[[-3, 1, 3, 2, -2, -1, 4, 5, 6, -6, -5, -4]],
        orientation=ones(Int, 6))
    candidate = only(movable(source, "MP";
        param="D2", vertices=[1, 2, 3], strict=true))
    expected = OGraphState([2, -1, 1, -2, 3, 4, 5, -5, -4, -3],
                           ones(Int, 5))
    source, candidate, expected
end

@testset "source-grounded o-graph move verifier" begin
    source, candidate, expected = d2_inverse_fixture()
    @test candidate.locations == [(1, 1), (1, 3), (1, 5)]

    report = verify_ograph_move(source, candidate;
        expected_raw_target=expected,
        expected_canonical_target=canonicalize_ograph(expected))
    @test report.valid
    @test report.reason == :verified
    witness = something(report.witness)
    @test witness.move_type == "MP"
    @test witness.parameter == "D2"
    @test witness.direction == :inverse
    @test witness.target_raw == expected
    @test witness.source_to_target_vertices == [1, 2, 0, 3, 4, 5]
    @test [(item.source_position, item.target_position,
            item.source_entry, item.target_entry)
           for item in witness.outside_occurrences] ==
          [(7, 5, 4, 3), (8, 6, 5, 4), (9, 7, 6, 5),
           (10, 8, -6, -5), (11, 9, -5, -4), (12, 10, -4, -3)]
    @test witness.outside_adjacencies ==
          [(7, 8, 5, 6), (8, 9, 6, 7), (9, 10, 7, 8),
           (10, 11, 8, 9), (11, 12, 9, 10)]
    @test witness.method_b === nothing
    @test witness.source_raw_sha256 == ograph_raw_hash(witness.source_raw)
    @test witness.target_raw_sha256 == ograph_raw_hash(witness.target_raw)
    @test replay_ograph_move(witness).valid
    encoded = ograph_move_witness_data(witness)
    @test ograph_move_witness_from_data(encoded) isa OGraphMoveWitness
    @test replay_ograph_move(encoded).valid

    unknown_field = (; encoded..., unsupported=true)
    @test replay_ograph_move(unknown_field).reason == :stored_witness_invalid
    wrong_hash = merge(encoded, (source_raw_sha256=repeat("0", 64),))
    @test replay_ograph_move(wrong_hash).reason == :witness_mismatch
    wrong_type = merge(encoded, (vertices=Any[true, 2, 3],))
    @test replay_ograph_move(wrong_type).reason == :stored_witness_invalid

    wrong_target = OGraphState([1, -1], [1])
    @test verify_ograph_move(source, candidate;
        expected_raw_target=wrong_target).reason == :raw_target_mismatch

    forged = MoveCandidate("MP", candidate.parameter, candidate.vertices,
        candidate.edges, candidate.direction, candidate.before, candidate.after,
        [(1, 2), candidate.locations[2], candidate.locations[3]])
    @test verify_ograph_move(source, forged).reason == :source_mismatch

    tampered = deepcopy(witness)
    tampered.source_to_target_vertices[4] = 99
    tampered_report = replay_ograph_move(tampered)
    @test !tampered_report.valid
    @test tampered_report.reason == :witness_mismatch

    unaudited = VirtualLink()
    set_data!(unaudited;
        rgauss=[[-1, -2, 3, 4, -4, -3, 2, 1]],
        orientation=[-1, 1, -1, 1])
    a2 = only(movable(unaudited, "MP";
        param="A2", vertices=[2, 1], strict=true))
    @test verify_ograph_move(unaudited, a2).reason ==
          :mp_parameter_not_source_audited
end

@testset "Method-B PS-III occurrence/lane witness" begin
    source = VirtualLink()
    set_data!(source;
        rgauss=[[-4, -3, -2, -1, 1, 2, 3, 4]],
        orientation=[1, -1, -1, -1])
    candidate = only(movable(source, "PS"; vertices=[1, 2], strict=true))
    @test candidate.parameter == "3"
    @test candidate.locations == [(1, 5), (1, 3)]

    report = verify_ograph_move(source, candidate)
    @test report.valid
    witness = something(report.witness)
    method = something(witness.method_b)
    @test method.paper_type == "PS-III"
    @test method.evaluated_state_side == :target_raw
    @test method.evaluated_state_sha256 == witness.target_raw_sha256
    @test (method.left_edge_id, method.left_lane) == (3, 1)
    @test (method.right_edge_id, method.right_lane) == (3, 2)
    @test [gap.fragment_index for gap in method.gap_boundaries] == [1, 2]
    @test [gap.source_fragment_positions for gap in method.gap_boundaries] ==
          [[5, 6], [3, 4]]
    @test [(gap.source_predecessor_position, gap.source_successor_position,
            gap.source_predecessor_entry, gap.source_successor_entry)
           for gap in method.gap_boundaries] == [(2, 7, -3, 3), (2, 7, -3, 3)]
    @test [(gap.target_component, gap.target_tail_occurrence,
            gap.target_head_occurrence, gap.target_tail_entry,
            gap.target_head_entry, gap.target_edge_id)
           for gap in method.gap_boundaries] ==
          [(1, 2, 3, -1, 1, 3), (1, 2, 3, -1, 1, 3)]
    @test any(step -> step.edge_id == method.right_edge_id &&
                      step.lane == method.right_lane,
              method.circuit_trace)
    @test method.circuit_id == minimum(
        min(step.from_node, step.to_node) for step in method.circuit_trace)
    @test replay_ograph_move(witness).valid
    @test replay_ograph_move(ograph_move_witness_data(witness)).valid

    encoded_ps = ograph_move_witness_data(witness)
    wrong_method = merge(encoded_ps.method_b,
                         (evaluated_state_sha256=repeat("0", 64),))
    wrong_side = merge(encoded_ps, (method_b=wrong_method,))
    @test replay_ograph_move(wrong_side).reason == :witness_mismatch

    tampered = deepcopy(witness)
    something(tampered.method_b).circuit_trace[1] =
        OGraphMethodBLaneStep(99, 1, 1, 2, 1, 1, 2)
    tampered_report = replay_ograph_move(tampered)
    @test !tampered_report.valid
    @test tampered_report.reason == :witness_mismatch

    tiny = VirtualLink()
    set_data!(tiny; rgauss=[[1, 2, -1, -2]], orientation=[1, -1])
    data = make_move_data()
    row1 = MoveCandidate("PS", "1", [1, 2], NTuple{2,Int}[], :inverse,
                         data["PS-R"][1], only(data["PS-L"]), [(1, 1), (1, 3)])
    @test verify_ograph_move(tiny, row1).reason ==
          :ps_parameter_not_source_audited
end
