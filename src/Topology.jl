module Topology

using LinearAlgebra
using JuMP
using GLPK
using Graphs
import AkaghefAlgebra

export URDiagram, VirtualLink, lens_space_rgauss, lens_space_link, SymExpr, sym,
       LaurentPolynomial, LaurentRational, laurent_variable, laurent_variables, evaluate, limit_zero,
       set_data_cvw!, set_data_cvs!, set_data_v!, set_data_cm!, set_data!, set_weight!,
       set_weight_by_crossing!,
       validate_data, check_data, swap!, delete_edge!, add_edges!, dilation!,
       compose_dilation!, put_zero!, reduct3!, reduct4!, reduct5!, reduct6!, reduct7!,
       reduct8!, simplify1!, simplify2!, simplify2_history!, simplify!,
       trace_value, trace1, trace2, trace3, standard_matrix_adjustment, standard_matrix_trace,
       rank, odata, weight_expression, unicode_art, trans_knot_complement,
       gauss_code, real_gauss_code, real_edge_table, RealEdge, pd_code, edge_directions, regions,
       head_map, is_closed, strand_cuts, set_cuts!,
       hmove_matrix, disk_table, disk_tetrahedra, calculate_weight, weight_equation, weight_residual,
       validate_weight, is_framed_weight, first_betti_number, boundary_weight_data,
       circuit_b_boundary, circuit_b_at, plot_circuit_b_svg,
       mirror_manifold!, cyclic_shift!,
       MoveCandidate, movable, move, move!, mp_move!, move_boundary_transport,
       connected_sum!, disjoint_sum!, knot_complement!, vl_to_urdiagram,
       ur_invariant, ur_prime_invariant, EdgeLayout, LinkLayout,
       minimal_bending_numbers, orthogonal_layout, plot_svg, plot_urd_structure_svg, bernstein, bezier,
       LayoutDiagnostic, check_layout, is_valid_layout, failing_rules, LAYOUT_RULES,
       planarize_virtual_gauss,
       MovePattern, make_move_data, SmallObject
export OGraphState, ClosednessWitness, ClosedNormalOGraph,
       OGRAPH_PRESENTATION_SCHEMA, OGRAPH_PRESENTATION_CANON_VERSION,
       OGRAPH_PRESENTATION_CANONICALIZATION, OGRAPH_DEFAULT_CONVENTION,
       OGRAPH_CLOSEDNESS_VERIFIER_VERSION, C2_GRAPH_SPEC_VERSION,
       C2Graph, c2_graph, true_vertex_count,
       canonicalize_ograph, ograph_canonical_serialization,
       ograph_canonical_hash, ograph_state, closedness_witness,
       certify_closed_ograph, ograph_graph_a,
       OGRAPH_MOVE_WITNESS_SCHEMA, OGRAPH_MOVE_VERIFIER_VERSION,
       OGraphOutsideOccurrence, OGraphCanonicalizationWitness,
       OGraphMethodBLaneStep, OGraphPSGapBoundary, OGraphPSMethodBWitness,
       OGraphMoveWitness, OGraphMoveVerification,
       ograph_raw_hash, verify_ograph_move, ps_global_condition,
       ograph_move_witness_data, ograph_move_witness_from_data,
       replay_ograph_move
export WeightedMoveEvidence, WeightedMoveVerification,
       WEIGHTED_MOVE_EVIDENCE_SCHEMA, WEIGHTED_MOVE_VERIFIER_VERSION,
       WEIGHTED_MOVE_SOURCE_RULE_VERSION,
       verify_weighted_move, weighted_move_evidence_data,
       weighted_move_evidence_from_data, replay_weighted_move, weighted_raw_hash,
       verify_weighted_relative_result, replay_weighted_relative_result
export weighted_ps_candidate
export spin_cp_move, spin_cp_transport
export weight_mode, set_weight_mode!, spin_weight_equation, validate_spin_weight, is_spin_weight,
       spin_boundary_data, solve_spin_weight_system, solve_spin_weights
export solve_integer_weight_system, solve_integer_weights,
       solve_integer_weight_lift, verify_integer_weight_solution
export planarize_virtual_gauss_layer8
export SparseFactor, FiniteHopfAlgebra, uqsl2_borel_small, sweedler_algebra,
       exterior_algebra, hopf_invariant, label_tensor, contract_network,
       dense_tensor, hopf_tensors, validate_tensor_expression,
       evaluate_tensor_expression, verify_hopf,
       is_involutory, is_unimodular
export legacy_hopf_tensors
export FiniteAlgebra, AlgebraElement, basis_element, group_algebra, cyclic_group_algebra,
       kac_paljutkin_algebra, dual_algebra, hopf_pairing, heisenberg_double,
       drinfeld_double

export KnotGaussDiagram, KnotGaussEnumeration, ReidemeisterCandidate,
       canonicalize, canonical_serialization, canonical_hash, crossing_count,
       apply_reidemeister, verify_reidemeister_transition,
       ClassicalGaussCertificate, ClassicalGaussImmersion,
       classical_gauss_certificate, classical_gauss_realizable,
       classical_gauss_realizable_fast,
       ExactKnotLaurentPolynomial, KnotPD, pd_matrix, pd_from_gauss, pd_to_gauss,
       pd_from_matrix, verify_pd, polynomial_terms, alexander_polynomial,
       jones_polynomial, invert_knot_variable, knot_polynomial_string,
       MoveGraphStore, CanonicalState, MoveType, MoveGraphStoreDB, Store,
       AREA_2, AREA_2_5, AREA_3, open_store, close_store!, migrate!,
       state_id, state_hash, put_state!, get_state, register_move_type!,
       get_move_type, insert_verified_transition!, record_transition!,
       add_certificate!, add_invariant!, neighbors, bfs_certificate_path,
       certificate_path, ClassificationRun,
       register_classification_run!, put_state_classification!,
       record_state_classification!, get_classification_run,
       classification_run, get_state_classification, state_classification,
       classification_rows, classification_counts,
       OGraphProvenance, OGraphRelationProfile,
       OGRAPH_STORE_SCHEMA_VERSION, OGRAPH_IDENTITY_DOMAIN,
       migrate_ograph_store!, ograph_canonical_state, put_closed_ograph!,
       OGRAPH_TRANSITION_SCHEMA_VERSION, OGRAPH_TRANSITION_EVIDENCE_VERSION,
       migrate_ograph_transition_store!, insert_ograph_verified_transition!,
       register_ograph_profile!, register_ograph_campaign!,
       verify_ograph_store, ograph_canonical_json,
       OGraphPi1Presentation, OGraphH1Result,
       OGRAPH_PI1_PRESENTATION_SCHEMA, OGRAPH_PI1_ALGORITHM,
       OGRAPH_PI1_ALGORITHM_VERSION, OGRAPH_PI1_CONVENTION,
       OGRAPH_H1_SCHEMA, OGRAPH_H1_ALGORITHM, OGRAPH_H1_ALGORITHM_VERSION,
       OGRAPH_INVARIANT_STORE_SCHEMA_VERSION,
       ograph_pi1_presentation, ograph_h1,
       migrate_ograph_invariant_store!, put_ograph_pi1_presentation!,
       put_ograph_h1!, put_ograph_invariants!,
       verify_ograph_invariant_store,
       knot_canonical_state, persist_reidemeister_transition!,
       expand_reidemeister_neighbors!,
       KnotN5Edge, KnotN5KnownRepresentative, KnotN5ComponentInvariant,
       KnotN5Classification, KNOT_N5_CLASSIFICATION_VERSION,
       KNOT_N5_EXPECTED_STATE_COUNTS, KNOT_N5_KNOWN_NAMES,
       KNOT_N5_STRICT_KNOT_TYPES, build_knot_n5_classification,
       classify_knot_n5, knot_n5_summary
export ClassicalGaussEnumerationResult, CLASSICAL_GAUSS_ENUMERATION_VERSION,
       rooted_unsigned_gauss_count, foreach_rooted_unsigned_gauss_word,
       enumerate_rooted_unsigned_gauss_words,
       classical_gauss_key_to_code, classical_gauss_key_to_diagram,
       direct_classical_gauss_enumeration,
       KnotN6Classification, KNOT_N6_CLASSIFICATION_VERSION, KNOT_N6_MAX_N,
       build_knot_n6_classification, classify_knot_n6, knot_n6_summary,
       KnotN6TargetedBridgeResult, KNOT_N6_TARGETED_BRIDGE_VERSION,
       KNOT_N6_TARGETED_BRIDGE_BARRIER, run_knot_n6_targeted_bridge,
       build_knot_n6_targeted_bridge, targeted_knot_n6_bridge,
       knot_n6_targeted_bridge_summary,
       KnotN6TargetedBarrier8Result, KNOT_N6_TARGETED_BARRIER8_VERSION,
       KNOT_N6_TARGETED_BARRIER8, run_knot_n6_targeted_barrier8,
       build_knot_n6_targeted_barrier8, targeted_knot_n6_barrier8,
       knot_n6_targeted_barrier8_summary,
       KnotGaussCanonicalStreamStats, KNOT_N6_ALL_STATE_ENUMERATION_VERSION,
       stream_canonical_knot_gauss, enumerate_canonical_knot_gauss

include("symbolic.jl")
include("laurent.jl")
include("ur_diagram.jl")
include("virtual_link.jl")
include("circuit_b.jl")
include("integer_weight_solver.jl")
include("ograph_state.jl")
include("weight_validation.jl")
include("spin_weights.jl")
include("layer8_planarization.jl")
include("virtual_planarization.jl")
include("orthogonal_layout.jl")
include("layout_invariants.jl")
include("helpers.jl")
include("moves.jl")
include("ograph_move_verifier.jl")
include("ograph_tetrahedra.jl")
include("weighted_ps.jl")
include("weighted_move_verifier.jl")
include("spin_cp.jl")
include("hopf_invariant.jl")
include("taft_invariant.jl")
include("ur_mst.jl")
include("tensor_network.jl")
include("hopf_properties.jl")
include("finite_algebra.jl")
include("svg_plot.jl")
include("knot_gauss.jl")
include("reidemeister.jl")
include("classical_gauss_realizability.jl")
include("knot_classical_invariants.jl")
include("classical_gauss_enumeration.jl")
include("knot_n6_classification.jl")
include("knot_n6_targeted_bridge.jl")
include("knot_n6_targeted_barrier8.jl")
include("knot_n6_all_state_enumeration.jl")
include("move_graph_store.jl")
using .MoveGraphStore
include("ograph_move_db.jl")
include("ograph_invariants.jl")
include("knot_move_db.jl")
include("knot_n5_classification.jl")

# Fixture-friendly package overload.  The source-level production builder is
# zero-argument and reconstructs max_n=6; this overload keeps the transient
# bridge workspace private while allowing a small package smoke test.
function build_knot_n6_targeted_barrier8(
        classification::KnotN6Classification; progress::Bool=true)
    workspace_ref = Ref{Any}(nothing)
    bridge7 = run_knot_n6_targeted_bridge(
        classification; barrier=7, progress=progress,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; barrier=8,
        progress=progress)
end

end
