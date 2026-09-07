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

    a2_source = VirtualLink()
    set_data!(a2_source;
        rgauss=[[-1, -2, 3, 4, -4, -3, 2, 1]],
        orientation=[-1, 1, -1, 1])
    a2 = only(movable(a2_source, "MP";
        param="A2", vertices=[2, 1], strict=true))
    @test verify_ograph_move(a2_source, a2).valid
end

@testset "GKS26 CP source gate and boundary-rooted certificate" begin
    gate = Topology.ograph_cp_catalogue_source_gate()
    @test gate.valid
    @test gate.source_rule_version == Topology.OGRAPH_CP_SOURCE_RULE_VERSION
    @test !gate.geometric_authority_verified
    @test !gate.certifiable
    tampered_catalogue = make_move_data()
    tampered_catalogue["CP-R"] = [MovePattern("inv", [[1, 2, -5, 3, -2,
        -3, -4, 4, -1, 5]], [1, 1, 1, -1, 1])]
    @test !Topology._cp_catalogue_source_gate(tampered_catalogue).valid

    recipe_audit = Topology.ograph_cp_negation_recipe_audit()
    @test recipe_audit.valid
    @test recipe_audit.reason == :conditional_recipe_conjugacy_verified
    @test recipe_audit.local_frame_cases == 8
    @test recipe_audit.face_pairing_cases == 16
    @test recipe_audit.gks_ordered_face_correspondence_cases == 16
    @test recipe_audit.tetrahedron_vertex_permutation == [1, 0, 3, 2]
    @test recipe_audit.germ_permutation == [4, 3, 2, 1, 6, 5]
    @test recipe_audit.lane_permutation == [2, 1, 3]
    @test recipe_audit.tetrahedron_permutation_even
    @test recipe_audit.local_orientation_sign_reversed
    @test recipe_audit.conditional_orientation_reversal
    @test recipe_audit.gks_figure_2_2_mapping_transcribed
    @test recipe_audit.gks_ordered_face_correspondence_verified
    @test recipe_audit.gks_figure_2_2_index_overlay_reviewed
    @test recipe_audit.recipe_source_authority_verified
    @test recipe_audit.cp_endpoint_preservation_proved
    @test recipe_audit.certifiable
    semantic_gate = Topology._cp_semantic_authority_gate()
    @test semantic_gate.valid
    @test semantic_gate.catalogue_correspondence_verified
    @test semantic_gate.recipe_naturality_verified
    @test semantic_gate.source_theorem_applied
    @test semantic_gate.certifiable

    certificate = Topology.ograph_cp_boundary_certificate()
    audit = Topology.audit_ograph_cp_boundary_certificate(certificate)
    @test audit.valid
    @test audit.left_move_count == 14
    @test audit.right_move_count == 6
    @test audit.max_vertices == 6
    @test audit.multi_circuit_left == collect(7:13)
    @test audit.multi_circuit_right == collect(1:5)
    @test audit.boundary_roots_preserved
    @test !audit.geometric_authority_verified
    @test !audit.certifiable
    @test !audit.semantic_moves_verified

    cyclically_rerooted = deepcopy(certificate)
    open_word = cyclically_rerooted.left_steps[10].target.open_word
    open_word[:] = circshift(open_word, 1)
    bad = Topology.audit_ograph_cp_boundary_certificate(cyclically_rerooted)
    @test !bad.valid
    @test bad.reason == :source_transcription_mismatch
end

@testset "source-grounded CP candidate and inverse replay" begin
    source = VirtualLink()
    set_data!(source; rgauss=[[-1, 1, -2, 2, -3, 3]],
              orientation=[-1, -1, -1])
    candidate = only(movable(source, "CP"; vertices=[1, 2, 3], strict=true))
    @test candidate.parameter == "dir"
    @test candidate.direction == :forward
    @test candidate.locations == [(1, 2)]
    expected_canonical = OGraphState(
        [-1, -2, -3, 3, -4, 5, 4, 1, -5, 2], [1, 1, -1, 1, -1])
    audit = Topology.audit_ograph_cp_candidate(source, candidate)
    @test audit.structural_checks_pass
    @test audit.geometric_authority_verified
    @test audit.certifiable
    @test audit.reason == :cp_source_geometry_verified
    @test canonicalize_ograph(something(audit.source_rule_target)) == expected_canonical
    @test audit.production_target == audit.source_rule_target
    @test audit.source_to_target_vertices == [1, 2, 3]
    report = verify_ograph_move(source, candidate;
        expected_canonical_target=expected_canonical)
    @test report.valid
    @test report.reason == :verified
    @test report.witness !== nothing
    @test replay_ograph_move(something(report.witness)).valid

    inverse_source = VirtualLink()
    raw_target = something(audit.source_rule_target)
    set_data!(inverse_source; rgauss=[collect(raw_target.gauss)],
              orientation=collect(raw_target.vertex_signs))
    inverse = only(movable(inverse_source, "CP";
        vertices=[1, 2, 3, 4, 5], strict=true))
    inverse_audit = Topology.audit_ograph_cp_candidate(inverse_source, inverse)
    @test inverse_audit.structural_checks_pass
    @test canonicalize_ograph(something(inverse_audit.source_rule_target)) ==
          canonicalize_ograph(ograph_state(source))
    @test inverse_audit.source_to_target_vertices == [1, 2, 3, 0, 0]
    @test inverse_audit.certifiable
    inverse_report = verify_ograph_move(inverse_source, inverse;
        expected_canonical_target=canonicalize_ograph(ograph_state(source)))
    @test inverse_report.valid
    @test replay_ograph_move(something(inverse_report.witness)).valid

    # Relabel so the source-defined removed roles 4--5 are concrete labels 1--2.
    # Production must retain roles 1--3 rather than sorting global labels.
    labels = [4, 5, 3, 1, 2]
    relabel_word(word) = [sign(entry) * labels[abs(entry)] for entry in word]
    relabel_signs = zeros(Int, 5)
    for old in 1:5
        relabel_signs[labels[old]] = raw_target.vertex_signs[old]
    end
    relabelled = VirtualLink()
    set_data!(relabelled; rgauss=[relabel_word(raw_target.gauss)],
              orientation=relabel_signs)
    inverse_roles = only(movable(relabelled, "CP"; vertices=labels, strict=true))
    relabelled_audit = Topology.audit_ograph_cp_candidate(relabelled, inverse_roles)
    @test relabelled_audit.structural_checks_pass
    @test relabelled_audit.certifiable
    @test relabelled_audit.source_to_target_vertices == [0, 0, 1, 2, 3]
    relabelled_report = verify_ograph_move(relabelled, inverse_roles)
    @test relabelled_report.valid
    @test replay_ograph_move(something(relabelled_report.witness)).valid

    # Swapping source-defined roles 4 and 5 while reusing the frozen location is
    # not a permissible global-label permutation and must fail at source match.
    wrong_roles = MoveCandidate(
        inverse_roles.type, inverse_roles.parameter,
        vcat(inverse_roles.vertices[1:3], reverse(inverse_roles.vertices[4:5])),
        copy(inverse_roles.edges), inverse_roles.direction,
        inverse_roles.before, inverse_roles.after, copy(inverse_roles.locations))
    negative = Topology.audit_ograph_cp_candidate(relabelled, wrong_roles)
    @test !negative.structural_checks_pass
    @test !negative.certifiable
    @test negative.reason != :cp_source_geometry_verified
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
          :source_not_closed
end
