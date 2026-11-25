module Topology

using LinearAlgebra
using JuMP
using GLPK
using Graphs
import AkaghefAlgebra
import AbstractAlgebra

export URDiagram, VirtualLink, HeegaardDiagram, lens_space_rgauss, lens_space_link,
       framed_lens_space_input, framed_lens_space_link, framed_lens_space_series,
       SymExpr, sym,
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
       switch_crossings, switch_crossings!, reflect_diagram, reflect_diagram!,
       reverse_strands, reverse_strands!, mst_sign_change, mst_sign_change!,
       mst_reflect_reverse, mst_reflect_reverse!,
       MoveCandidate, movable, move, move!, mp_move!, move_boundary_transport,
       connected_sum!, disjoint_sum!, knot_complement!, vl_to_urdiagram,
       ur_invariant, ur_prime_invariant, EdgeLayout, LinkLayout,
       minimal_bending_numbers, orthogonal_layout, plot_svg, plot_urd_structure_svg, bernstein, bezier,
       LayoutDiagnostic, check_layout, is_valid_layout, failing_rules, LAYOUT_RULES,
       gauss_diagram_svg, plot_gauss_diagram_svg,
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
       OGRAPH_COLOUR_CONVENTION_ID, OGRAPH_FACE_CYCLE_STATUS,
       ColouredOGraph, DiscGerm, EdgeColour, OGraphColourCheck,
       NormalEdgeOrientation, ClosedNormalColouring,
       fig25_tau, fig25_leg_direction, fig25_leg_pairing,
       fig25_gauged_colour, fig25_frame_reversal,
       edge_gluing_from_colour, colour_from_gluing, edge_disc_germs,
       validate_coloured_ograph, canonical_coloured_ograph,
       coloured_ograph_serialization, coloured_ograph_hash,
       prop321_colour, normal_coloured_ograph, validate_normal_colouring,
       normal_orientations, prop321_orientations, coloured_gc_display,
       FaceCycle, ColouredPolyhedron, ColourHomology, DualFacePairing,
       SpineBoundaryComponent, SpineBoundary, DiscPartitionComparison,
       face_cycles, cellular_homology, dual_face_pairings, dual_edge_orbits,
       sorted_germ_partition, disk_table_corners, compare_disc_partition,
       thickening_boundary,
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
       OGRAPH_SEARCH_SCHEMA_VERSION, OGRAPH_SEARCH_STATE_VERSION,
       OGRAPH_SEARCH_ENUMERATOR_VERSION, OGRAPH_SEARCH_ENGINE_VERSION,
       migrate_ograph_search_store!, ograph_search_state,
       put_ograph_search_state!, load_ograph_search_state,
       begin_ograph_search_run!, finish_ograph_search_run!, ograph_search_run,
       ograph_search_candidates, expand_ograph_search_state!,
       ograph_search_neighbors, query_ograph_search_path,
       query_ograph_search_paths,
       replay_weighted_ograph_path, record_weighted_ograph_path_trial!,
       verify_ograph_search_store,
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

include("Algebra.jl")
include("MoveAPI.jl")
include("Invariants.jl")
include("move_graph_store.jl")
using .MoveGraphStore
include("URDiagrams.jl")
include("KnotGauss.jl")
include("VirtualLinks.jl")

# Bind each presentation name on Topology, then re-export what that module
# exports. Heegaard, included below, reads VirtualLink and the move helpers
# through parentmodule, so this runs first. The same object is not imported twice.
let modules = (:Algebra, :MoveAPI, :Invariants, :URDiagrams, :KnotGauss, :VirtualLinks)
    skip(name::Symbol) = begin
        name === :eval && return true
        name === :include && return true
        Base.isidentifier(name) || return true
        text = String(name)
        isempty(text) || text[1] == '#'
    end
    for modname in modules
        child = getfield(@__MODULE__, modname)
        seen = Set{Symbol}()
        for name in names(child; all=true)
            push!(seen, name)
        end
        for name in names(child)
            push!(seen, name)
        end
        for name in seen
            skip(name) && continue
            isdefined(child, name) || continue
            if isdefined(@__MODULE__, name)
                getfield(@__MODULE__, name) === getfield(child, name) ||
                    error("Topology name clash: $(modname).$(name)")
            else
                @eval import .$modname: $name
            end
        end
        for name in names(child)
            skip(name) && continue
            isdefined(@__MODULE__, name) || continue
            @eval export $name
        end
    end
end

include("svg_plot.jl")

export ograph_move_weight_relation, solve_ograph_path_weights,
       record_ograph_path_constraint_trial!, query_ograph_path_constraint_trials,
       ograph_path_endpoint_permutations,
       lift_ograph_search_path!

# Fixture-friendly package overload.  The source-level production builder is
# zero-argument and reconstructs max_n=6; this overload keeps the transient
# bridge workspace private while allowing a small package smoke test.
"""
    build_knot_n6_targeted_barrier8(classification::KnotN6Classification; progress::Bool=true) -> KnotN6TargetedBarrier8Result

Return a targeted barrier-eight extension of an existing bounded knot classification.

# Arguments
`classification` must be a coherent result of
[`build_knot_n6_classification`](@ref) with `max_n` in `0:6`.
`progress` controls progress messages. An out-of-range `max_n` throws
`ArgumentError`; inconsistent indices, known labels, transitions or closure
data propagate the targeted-search errors. No independent validation of all
manually assembled classification fields is performed.

# Returns
`KnotN6TargetedBarrier8Result` with barrier, input fingerprints, seed and
resolved/unresolved component data, known-label connections, search and
verification counts, elapsed time, and `status`/`resource_reason` fields.
Use [`knot_n6_targeted_barrier8_summary`](@ref) for a compact view.

# Guarantees
1. The result describes finite verified-move connections extending the input's
   unresolved components through the targeted barrier-seven and barrier-eight
   searches. Its mathematical scope is that finite diagram graph. [proved/agent: by construction]
2. `barrier == 8`. Reported resolved labels come from connections to known
   input components; unresolved components are preserved as unresolved. [proved/agent: by construction]
3. The supplied classification is not modified by this operation. [proved/agent: by construction]

# Not guaranteed
This is not an eight-crossing census or an unbounded knot-equivalence
decision. Resource termination can produce a partial result; inspect
`status` and `resource_reason`. No full paths or parent arrays are returned.
The [module scope](../specs/module-structure.md#1-構成) does not itself
certify the local-move or classical-realizability conventions.

# Sources
cf: lib/topology.jl/src/reidemeister.jl verify_reidemeister_transition — local transition evidence underlying the finite graph.

# Examples
```jldoctest
julia> c = build_knot_n6_classification(3; barrier=3, progress=false);

julia> r = build_knot_n6_targeted_barrier8(c; progress=false);

julia> (isempty(r.seed_components), r.n7_state_count, r.n8_state_count, r.status)
(true, 0, 0, "barrier8_resolved_all_seed_components")
```

See also [`knot_n6_targeted_barrier8_summary`](@ref).

Role: convert
"""
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

export HeegaardCodes, Heegaard
include("heegaard/Heegaard.jl")
# Old module name, kept as an alias (docs/specs/module-structure.md §2).
const Heegaard = HeegaardCodes
# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/docs/DESIGN.md §判断 — 公開 HeegaardDiagram は解決済み HC（struct HeegaardCode）と同じ型。standard は α の支持を既存の m に束縛する入力略記であり、別型ではない
# cf: lib/topology.jl/src/heegaard/hc_compare.jl::compare_hc — same word, same homology, and move-equivalence
# cf: lib/topology.jl/docs/specs/heegaard-diagram.md — one resolved diagram; these names are the existing operations
using .HeegaardCodes: HeegaardCode, draw_heegaard, compare_hc, validate_diagram,
    hc_intersection, hc_first_homology, hc_to_odata, is_standard
const HeegaardDiagram = HeegaardCode
move(diagram::HeegaardDiagram, witness::Heegaard.LCMoveWitness) =
    Heegaard.move(diagram, witness)
movable(diagram::HeegaardDiagram, row::AbstractDict; kwargs...) =
    Heegaard.movable(diagram, row; kwargs...)
movable(diagram::HeegaardDiagram, row_id::AbstractString; kwargs...) =
    Heegaard.movable(diagram, row_id; kwargs...)
export draw_heegaard, compare_hc, validate_diagram, is_standard,
    hc_intersection, hc_first_homology, hc_to_odata
include("kuperberg_invariant.jl")
export Z_Kuperberg, KuperbergTensorExpression

# Manifold operations. Generics: Invariants.mirror, Invariants.connected_sum.
# MST25 Proposition A.1 sends Γ× to (M^op, -v^op). mst_sign_change is that
# checked operation. mirror_manifold! does the same array action and does not
# check that the input is a closed normal o-graph.
mirror(x::HeegaardCodes.HeegaardCode) = HeegaardCodes.hc_mirror(x)
mirror(x::VirtualLink) = mst_sign_change(x)
connected_sum(x::HeegaardCodes.HeegaardCode, y::HeegaardCodes.HeegaardCode; kwargs...) =
    HeegaardCodes.hc_connected_sum(x, y; kwargs...)
connected_sum(x::VirtualLink, y=(1, 2); transport_weights::Bool=true) =
    connected_sum!(deepcopy(x), y; transport_weights)

# Search database. The presentation type is the argument after the store
# (docs/specs/module-structure.md §3.1). The o-graph calls without that type
# stay in VirtualLinks. ograph_h1 takes an OGraphPi1Presentation, and no
# VirtualLink linking form exists, so those two have no method here.
load_search_state(store::MoveGraphStoreDB, ::Type{HeegaardCode}, state_id::AbstractString) =
    HeegaardCodes.load_hc_search_state(store, state_id)
load_search_state(store::MoveGraphStoreDB, ::Type{VirtualLink}, state_id::AbstractString) =
    load_ograph_search_state(store, state_id)
begin_search_run!(store::MoveGraphStoreDB, ::Type{HeegaardCode}, script_path::AbstractString; kwargs...) =
    HeegaardCodes.begin_hc_search_run!(store, script_path; kwargs...)
begin_search_run!(store::MoveGraphStoreDB, ::Type{VirtualLink}, script_path::AbstractString; kwargs...) =
    begin_ograph_search_run!(store, script_path; kwargs...)
finish_search_run!(store::MoveGraphStoreDB, ::Type{HeegaardCode}, run_id::AbstractString; kwargs...) =
    HeegaardCodes.finish_hc_search_run!(store, run_id; kwargs...)
finish_search_run!(store::MoveGraphStoreDB, ::Type{VirtualLink}, run_id::AbstractString; kwargs...) =
    finish_ograph_search_run!(store, run_id; kwargs...)
"""
    search_neighbors(store::MoveGraphStoreDB, ::Type{HeegaardCode}, state_id::AbstractString; row_ids=nothing) -> Vector

Read the stored outgoing edges of an HC presentation in the search graph.

# Arguments
`store` satisfies M7 of the [search-store conditions](../specs/move-graph-store.md#numbered-conditions);
`state_id` is the source presentation's stored identifier.
Checked: the HC schema must be installed and supported; `ArgumentError` names
the absent schema or unsupported version. Database and string-conversion errors
propagate. `row_ids` is absent or an iterable of move-row names.
Caller's responsibility: trust the stored certificates only to the extent
they have been independently verified.

# Returns
Records containing `operation_id`, `run_id`, source and target state IDs,
`depth`, `row_id`, `inverse_status`, and decoded `candidate` data.

# Guarantees
1. Only stored successful outgoing HC edges matching the requested row names
   are returned, ordered by depth, row name and operation ID. [proved/agent: by construction]
2. Database records are not modified. [proved/agent: by construction]

# Not guaranteed
This does not enumerate unstored moves or recheck equivalence of manifolds.
An empty result is compatible with undiscovered neighbors.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — representation-selected search interface.

# Examples
```jldoctest
julia> s = open_store(); Topology.HeegaardCodes.migrate_hc_search_store!(s);

julia> isempty(search_neighbors(s, Topology.HeegaardCodes.HeegaardCode, "absent"; row_ids=String[]))
true

julia> close_store!(s);
```

Role: access
"""
search_neighbors(store::MoveGraphStoreDB, ::Type{HeegaardCode}, state_id::AbstractString; kwargs...) =
    HeegaardCodes.hc_search_neighbors(store, state_id; kwargs...)
"""
    search_neighbors(store::MoveGraphStoreDB, ::Type{VirtualLink}, state_id::AbstractString; kinds=("MP", "02")) -> Vector

Read the stored outgoing edges of an o-graph presentation in the search graph.

# Arguments
Use the store, identifier and kind conditions of
[`ograph_search_neighbors`](@ref) and the [search specification](../specs/ograph-search-db.md#contract).
The store requires M7 of the [store conditions](../specs/move-graph-store.md#numbered-conditions).
Checked: unsupported move kinds raise `ArgumentError`; schema and database
errors propagate. The caller supplies a readable store and interprets its
stored evidence independently of this query.

# Returns
The neighbor records of `ograph_search_neighbors`, including operation and
endpoint IDs, kind, direction, engine version and candidate JSON.

# Guarantees
1. The returned rows equal `ograph_search_neighbors(store,state_id; kinds)`.
   These are stored successful directed operations, not newly generated moves. [proved/agent: by construction]
2. The store is not modified. [proved/agent: by construction]

# Not guaranteed
The call does not replay transitions or establish manifold equivalence.

# Sources
cf:lib/topology.jl/docs/specs/ograph-search-db.md §Contract — the persisted presentation graph.

# Examples
```jldoctest
julia> s = open_store(); migrate_ograph_search_store!(s);

julia> isempty(search_neighbors(s, VirtualLink, "absent"; kinds=String[]))
true

julia> close_store!(s);
```

Role: access
"""
search_neighbors(store::MoveGraphStoreDB, ::Type{VirtualLink}, state_id::AbstractString; kwargs...) =
    ograph_search_neighbors(store, state_id; kwargs...)
query_search_paths(store::MoveGraphStoreDB, ::Type{HeegaardCode}, source_state_id::AbstractString; kwargs...) =
    HeegaardCodes.query_hc_search_paths(store, source_state_id; kwargs...)
query_search_paths(store::MoveGraphStoreDB, ::Type{VirtualLink}, source_state_id::AbstractString; kwargs...) =
    query_ograph_search_paths(store, source_state_id; kwargs...)
put_search_state!(store::MoveGraphStoreDB, hc::HeegaardCode; kwargs...) =
    HeegaardCodes.put_hc_search_state!(store, hc; kwargs...)
"""
    search_state(hc::HeegaardCode; identity=:canonical) -> CanonicalState

Construct the stored presentation identity of a Heegaard code.

# Arguments
`hc` satisfies the [HC data specification](../specs/heegaard-code.md).
Checked: `identity` selects `:canonical` or `:literal`; an unsupported name
raises `ArgumentError`. The caller supplies valid HC data; writing or
canonicalization failures propagate, and invalid data has no certified
presentation-identity interpretation.

# Returns
A `CanonicalState` with the serialized identity, its hash and versioned ID,
HC writing as payload, genus as complexity, and representation `hc-search`.

# Guarantees
1. The result equals `Topology.HeegaardCodes.hc_search_state(hc; identity)`
   under the [common search interface](../specs/module-structure.md#31-moveapi). [proved/agent: by construction]
2. The input is not modified and no database row is inserted. [proved/agent: by construction]

# Not guaranteed
Identity of stored presentations is not a classification of represented
manifolds. Canonicalization has the scope of `hc_canonical`.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — HC search-state adapter.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H._lc_parse_text("genus: 1\\nl1 : -1\\nm1 : +1\\n"); h = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> search_state(h; identity=:literal).representation
"hc-search"
```

Role: convert
"""
search_state(hc::HeegaardCode; kwargs...) =
    HeegaardCodes.hc_search_state(hc; kwargs...)
expand_search_state!(store::MoveGraphStoreDB, run_id::AbstractString, hc::HeegaardCode; kwargs...) =
    HeegaardCodes.expand_hc_search_state!(store, run_id, hc; kwargs...)
first_homology(x::HeegaardCode) = HeegaardCodes.hc_first_homology(x)
linking_form(x::HeegaardCode; kwargs...) = HeegaardCodes.hc_linking_form(x; kwargs...)

end
