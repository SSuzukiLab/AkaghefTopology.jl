# Source scope, index gates, and replay contract:
# docs/specs/ograph-move-verifier.md
const OGRAPH_MOVE_WITNESS_SCHEMA = "akaghef.ograph-move-witness/1"
const OGRAPH_MOVE_VERIFIER_VERSION = "ograph-move-source-replay/3"
const OGRAPH_CP_SOURCE_RULE_VERSION = "gks26-cp-edata-N/2"
const OGRAPH_CP_BOUNDARY_CERTIFICATE_SCHEMA = "akaghef.ograph-cp-boundary-factorization/1"

"One E-datum with the boundary-rooted component stored linearly from a to b."
struct OGraphBoundaryEData
    open_word::Vector{Int}
    closed_components::Vector{Vector{Int}}
    vertex_signs::Vector{Int}
end

"One source-transcribed arrow in the GKS26 CP reduction."
struct OGraphCPBoundaryStep
    move_name::String
    direction::Symbol
    vertices::Vector{Int}
    target::OGraphBoundaryEData
end

"The two boundary-rooted reductions to the common GKS26 terminal E-datum."
struct OGraphCPBoundaryCertificate
    schema::String
    source_rule_version::String
    paper_left::OGraphBoundaryEData
    paper_right::OGraphBoundaryEData
    right_relabel::Vector{Int}
    right_initial::OGraphBoundaryEData
    left_steps::Vector{OGraphCPBoundaryStep}
    right_steps::Vector{OGraphCPBoundaryStep}
end

"One unchanged Gauss occurrence tracked across the raw local rewrite."
struct OGraphOutsideOccurrence
    source_component::Int
    source_position::Int
    target_component::Int
    target_position::Int
    source_entry::Int
    target_entry::Int
end

"The exact root and label map used by presentation canonicalization."
struct OGraphCanonicalizationWitness
    root_occurrence::Int
    raw_to_canonical_labels::Vector{Int}
end

"One oriented Method-B lane segment, retaining its GC edge occurrence."
struct OGraphMethodBLaneStep
    edge_id::Int
    component::Int
    tail_occurrence::Int
    head_occurrence::Int
    lane::Int
    from_node::Int
    to_node::Int
end

"One inverse-PS source patch gap transported to an edge occurrence of the raw target."
struct OGraphPSGapBoundary
    fragment_index::Int
    source_component::Int
    source_fragment_positions::Vector{Int}
    source_predecessor_position::Int
    source_successor_position::Int
    source_predecessor_entry::Int
    source_successor_entry::Int
    target_component::Int
    target_tail_occurrence::Int
    target_head_occurrence::Int
    target_tail_entry::Int
    target_head_entry::Int
    target_edge_id::Int
end

"Source-grounded global circuit evidence for a pure sliding candidate."
struct OGraphPSMethodBWitness
    paper_type::String
    parameter::String
    evaluated_state_side::Symbol
    evaluated_state_sha256::String
    gap_boundaries::Vector{OGraphPSGapBoundary}
    left_edge_id::Int
    right_edge_id::Int
    left_lane::Int
    right_lane::Int
    circuit_id::Int
    circuit_trace::Vector{OGraphMethodBLaneStep}
end

"Replayable evidence for one concrete raw move followed by canonicalization."
struct OGraphMoveWitness
    schema::String
    verifier_version::String
    move_type::String
    parameter::String
    direction::Symbol
    vertices::Vector{Int}
    source_locations::Vector{Tuple{Int,Int}}
    target_locations::Vector{Tuple{Int,Int}}
    source_raw::OGraphState
    target_raw::OGraphState
    source_canonical::OGraphState
    target_canonical::OGraphState
    source_raw_sha256::String
    target_raw_sha256::String
    source_closedness::ClosednessWitness
    target_closedness::ClosednessWitness
    source_to_target_vertices::Vector{Int}
    outside_occurrences::Vector{OGraphOutsideOccurrence}
    outside_adjacencies::Vector{NTuple{4,Int}}
    source_canonicalization::OGraphCanonicalizationWitness
    target_canonicalization::OGraphCanonicalizationWitness
    method_b::Union{Nothing,OGraphPSMethodBWitness}
end

struct OGraphMoveVerification
    valid::Bool
    reason::Symbol
    witness::Union{Nothing,OGraphMoveWitness}
    error::Union{Nothing,String}
end

_ograph_move_invalid(reason::Symbol; error=nothing) =
    OGraphMoveVerification(false, reason, nothing,
                           error === nothing ? nothing : sprint(showerror, error))

function _ograph_raw_serialization(state::OGraphState)
    "OGraphRaw;schema=$(OGRAPH_PRESENTATION_SCHEMA);" *
        "convention-bytes=$(ncodeunits(state.convention_id)):$(state.convention_id);" *
        "n=$(true_vertex_count(state));gauss=$(join(state.gauss, ','));" *
        "vertex_signs=$(join(state.vertex_signs, ','))"
end

ograph_raw_hash(state::OGraphState) =
    bytes2hex(SHA.sha256(codeunits(_ograph_raw_serialization(state))))

_cp_edata(open_word, signs, closed=Vector{Vector{Int}}()) =
    OGraphBoundaryEData(Int.(open_word), [Int.(word) for word in closed], Int.(signs))

function _cp_negate(edata::OGraphBoundaryEData)
    OGraphBoundaryEData(-edata.open_word, [-word for word in edata.closed_components],
                        -edata.vertex_signs)
end

function _cp_relabel(edata::OGraphBoundaryEData, labels::AbstractVector{<:Integer})
    relabel(word) = [sign(entry) * labels[abs(entry)] for entry in word]
    OGraphBoundaryEData(relabel(edata.open_word), relabel.(edata.closed_components),
                        edata.vertex_signs[sortperm(Int.(labels))])
end

function _same_boundary_edata(left::OGraphBoundaryEData,
                              right::OGraphBoundaryEData)
    left.open_word == right.open_word &&
        left.closed_components == right.closed_components &&
        left.vertex_signs == right.vertex_signs
end

function _valid_boundary_edata(edata::OGraphBoundaryEData)
    n = length(edata.vertex_signs)
    n > 0 || return false
    all(sign -> sign in (-1, 1), edata.vertex_signs) || return false
    words = vcat([edata.open_word], edata.closed_components)
    all(word -> !isempty(word), words) || return false
    entries = reduce(vcat, words; init=Int[])
    length(entries) == 2n || return false
    sort(entries) == sort(vcat(collect(1:n), -collect(1:n)))
end

function _gks26_cp_paper_endpoints()
    left = _cp_edata([-1, 2, -2, 3, -3, 1], [1, 1, 1])
    right = _cp_edata([-1, -2, 5, -3, 2, 3, 4, -4, 1, -5],
                      [-1, -1, -1, 1, 1])
    left, right
end

function _cp_expected_catalogue_patterns()
    paper_left, paper_right = _gks26_cp_paper_endpoints()
    runtime_left = _cp_negate(paper_left)
    runtime_right = _cp_negate(paper_right)
    MovePattern("dir", [runtime_left.open_word], runtime_left.vertex_signs),
    MovePattern("inv", [runtime_right.open_word], runtime_right.vertex_signs)
end

function _cp_catalogue_source_gate(data=make_move_data())
    expected_left, expected_right = _cp_expected_catalogue_patterns()
    left = get(data, "CP-L", MovePattern[])
    right = get(data, "CP-R", MovePattern[])
    left_match = length(left) == 1 && _same_move_pattern(only(left), expected_left)
    right_match = length(right) == 1 && _same_move_pattern(only(right), expected_right)
    (valid=left_match && right_match,
     source_rule_version=OGRAPH_CP_SOURCE_RULE_VERSION,
     left_match=left_match, right_match=right_match,
     geometric_authority_verified=false,
     certifiable=false)
end

"Check the finite CP catalogue correspondence to GKS26 E-data; this is not geometric authority."
ograph_cp_catalogue_source_gate() = _cp_catalogue_source_gate()

# The four ordered face triples are derived from the exact production SOT used
# by `disk_table`, not copied into the verifier as a second magic table.
function _cp_recipe_germs(sign::Int)
    sign in (-1, 1) || throw(ArgumentError("crossing sign must be +/-1"))
    slots = _disk_node_dictionary()[sign]
    ntuple(face -> Tuple(slots[(3face - 2):(3face)]), 4)
end

# GKS26 Figure 2.2: repository face index -> omitted ordered-tetrahedron
# vertex.  The two rows are crossing signs, not tetrahedron types (the source
# assigns type - to a positive crossing and type + to a negative crossing).
const _CP_GKS_FACE_OMITTED_VERTEX = Dict(
    1 => Int[2, 0, 3, 1],
   -1 => Int[3, 1, 2, 0],
)

function _cp_recipe_edge(table, germ::Int)
    Tuple(face for face in 0:3 if !(germ in table[face + 1]))
end

function _cp_recipe_permutations4()
    [Int[a, b, c, d] for a in 0:3 for b in 0:3 for c in 0:3 for d in 0:3
     if length(Set((a, b, c, d))) == 4]
end

function _cp_recipe_face_map(left_sign::Int, left_face::Int,
                             right_sign::Int, right_face::Int)
    left = _cp_recipe_germs(left_sign)
    right = _cp_recipe_germs(right_sign)
    candidates = Vector{Vector{Int}}()
    for permutation in _cp_recipe_permutations4()
        permutation[left_face + 1] == right_face || continue
        all(Set(permutation[vertex + 1]
                for vertex in _cp_recipe_edge(left, left_germ)) ==
            Set(_cp_recipe_edge(right, right_germ))
            for (left_germ, right_germ) in
                zip(left[left_face + 1], right[right_face + 1])) || continue
        push!(candidates, permutation)
    end
    length(candidates) == 1 || return nothing
    only(candidates)
end

function _cp_permutation_even(permutation::AbstractVector{<:Integer})
    inversions = sum(permutation[left] > permutation[right]
                     for left in eachindex(permutation)
                     for right in (left + 1):lastindex(permutation))
    iseven(inversions)
end

_cp_tetrahedron_type(crossing_sign::Int) = -crossing_sign

function _cp_order_preserving_face_map(left_omitted::Int, right_omitted::Int)
    source_vertices = [vertex for vertex in 0:3 if vertex != left_omitted]
    target_vertices = [vertex for vertex in 0:3 if vertex != right_omitted]
    permutation = fill(-1, 4)
    permutation[left_omitted + 1] = right_omitted
    for (source, target) in zip(source_vertices, target_vertices)
        permutation[source + 1] = target
    end
    permutation
end

function _cp_gks_face_map(left_sign::Int, left_face::Int,
                          right_sign::Int, right_face::Int)
    left_labels = _CP_GKS_FACE_OMITTED_VERTEX[left_sign]
    right_labels = _CP_GKS_FACE_OMITTED_VERTEX[right_sign]
    right_inverse = [findfirst(==(vertex), right_labels) - 1 for vertex in 0:3]
    ordered = _cp_order_preserving_face_map(left_labels[left_face + 1],
                                            right_labels[right_face + 1])
    [right_inverse[ordered[left_labels[vertex + 1] + 1] + 1]
     for vertex in 0:3]
end

"""
    ograph_cp_negation_recipe_audit()

Check the finite conjugacy of the current disk-germ dual-triangulation recipe
under `N(G, eps)=(-G,-eps)`. The face triples come from the production
`disk_table` SOT, and every induced face map is compared with GKS26 Figure 2.2's
unique order-preserving gluing. Together with the reviewed source index overlay
and the source CP theorem, this is the semantic authority gate.
"""
function ograph_cp_negation_recipe_audit()
    tetrahedron_vertices = Int[1, 0, 3, 2] # p=(0 1)(2 3), indexed by 0:3
    germ_labels = Int[4, 3, 2, 1, 6, 5]    # mu=(1 4)(2 3)(5 6)
    lanes = Int[2, 1, 3]
    source_tables = Dict(sign => _cp_recipe_germs(sign) for sign in (-1, 1))
    all(sort(vcat((collect(face) for face in table)...)) ==
        sort(repeat(collect(1:6), 2)) for table in values(source_tables)) ||
        return (valid=false, reason=:production_germ_incidence_invalid,
                certifiable=false)
    all(length(unique(_cp_recipe_edge(table, germ) for germ in 1:6)) == 6
        for table in values(source_tables)) ||
        return (valid=false, reason=:production_germ_edge_gauge_invalid,
                certifiable=false)
    tetrahedron_permutation_even = _cp_permutation_even(tetrahedron_vertices)
    tetrahedron_permutation_even || return (
        valid=false, reason=:tetrahedron_permutation_not_even,
        certifiable=false)
    local_orientation_sign_reversed = all(
        _cp_tetrahedron_type(-sign) == -_cp_tetrahedron_type(sign)
        for sign in (-1, 1))
    local_orientation_sign_reversed || return (
        valid=false, reason=:tetrahedron_type_not_reversed,
        certifiable=false)
    local_frames = 0
    for vertex_sign in (-1, 1), face in 0:3
        source = source_tables[vertex_sign][face + 1]
        target = source_tables[-vertex_sign][tetrahedron_vertices[face + 1] + 1]
        [germ_labels[germ] for germ in source] == collect(target)[lanes] ||
            return (valid=false, reason=:local_frame_conjugacy_failed,
                    certifiable=false)
        local_frames += 1
    end

    face_pairing_cases = 0
    gks_face_correspondence_cases = 0
    for left_sign in (-1, 1), right_sign in (-1, 1),
        tail_positive in (false, true), head_positive in (false, true)
        left_face = 2 + Int(tail_positive)
        right_face = Int(head_positive)
        negated_left_face = 2 + Int(!tail_positive)
        negated_right_face = Int(!head_positive)
        negated_left_face == tetrahedron_vertices[left_face + 1] ||
            return (valid=false, reason=:tail_face_transport_failed,
                    certifiable=false)
        negated_right_face == tetrahedron_vertices[right_face + 1] ||
            return (valid=false, reason=:head_face_transport_failed,
                    certifiable=false)
        original = _cp_recipe_face_map(left_sign, left_face,
                                       right_sign, right_face)
        negated = _cp_recipe_face_map(-left_sign, negated_left_face,
                                      -right_sign, negated_right_face)
        original === nothing && return (
            valid=false, reason=:original_face_map_not_unique, certifiable=false)
        negated === nothing && return (
            valid=false, reason=:negated_face_map_not_unique, certifiable=false)
        all(negated[tetrahedron_vertices[vertex + 1] + 1] ==
            tetrahedron_vertices[original[vertex + 1] + 1]
            for vertex in 0:3) || return (
                valid=false, reason=:face_pairing_conjugacy_failed,
                certifiable=false)
        original == _cp_gks_face_map(left_sign, left_face,
                                     right_sign, right_face) || return (
            valid=false, reason=:gks_ordered_face_correspondence_failed,
            certifiable=false)
        gks_face_correspondence_cases += 1
        face_pairing_cases += 1
    end
    (
        valid=true,
        reason=:conditional_recipe_conjugacy_verified,
        local_frame_cases=local_frames,
        face_pairing_cases=face_pairing_cases,
        gks_ordered_face_correspondence_cases=gks_face_correspondence_cases,
        tetrahedron_vertex_permutation=tetrahedron_vertices,
        germ_permutation=germ_labels,
        lane_permutation=lanes,
        tetrahedron_permutation_even=tetrahedron_permutation_even,
        local_orientation_sign_reversed=local_orientation_sign_reversed,
        conditional_orientation_reversal=true,
        gks_figure_2_2_mapping_transcribed=true,
        gks_ordered_face_correspondence_verified=true,
        gks_figure_2_2_index_overlay_reviewed=true,
        recipe_source_authority_verified=true,
        cp_endpoint_preservation_proved=true,
        certifiable=true,
    )
end

function _cp_semantic_authority_gate(data=make_move_data())
    catalogue = _cp_catalogue_source_gate(data)
    recipe = ograph_cp_negation_recipe_audit()
    valid = catalogue.valid && recipe.valid &&
        recipe.gks_figure_2_2_index_overlay_reviewed &&
        recipe.recipe_source_authority_verified &&
        recipe.cp_endpoint_preservation_proved && recipe.certifiable
    (
        valid=valid,
        source_rule_version=OGRAPH_CP_SOURCE_RULE_VERSION,
        catalogue_correspondence_verified=catalogue.valid,
        recipe_naturality_verified=recipe.valid,
        source_theorem_applied=true,
        certifiable=valid,
    )
end

_cp_step(name, direction, vertices, target) =
    OGraphCPBoundaryStep(String(name), Symbol(direction), Int.(vertices), target)

function _gks26_cp_boundary_certificate()
    paper_left, paper_right = _gks26_cp_paper_endpoints()
    left_steps = OGraphCPBoundaryStep[
        _cp_step("C3", :forward, [2, 1],
            _cp_edata([2, -4, -1, -2, 3, -3, 4, 1], [-1, 1, 1, 1])),
        _cp_step("C3", :forward, [2, 3],
            _cp_edata([5, 2, -4, -1, 3, -5, -2, -3, 4, 1], [-1, -1, 1, 1, 1])),
        _cp_step("ps-I", :inverse, [5, 2],
            _cp_edata([-3, -1, 2, -2, 3, 1], [-1, 1, 1])),
        _cp_step("C1", :forward, [1, 2],
            _cp_edata([-3, 2, -4, -1, -2, 3, 1, 4], [1, 1, 1, -1])),
        _cp_step("ps-III", :inverse, [1, 4],
            _cp_edata([-2, 1, -1, 2], [1, 1])),
        _cp_step("C3", :forward, [1, 2],
            _cp_edata([-1, -2, 3, 1, 2, -3], [-1, 1, 1])),
        _cp_step("E2", :forward, [1, 2],
            _cp_edata([2, 3, 4, -3], [1, -1, 1, -1],
                      [[-1, -2, -4], [1]])),
        _cp_step("D1", :forward, [3, 4],
            _cp_edata([2, 5, 4, -5, 3], [1, -1, -1, -1, 1],
                      [[-1, -2, -3, -4], [1]])),
        _cp_step("B2", :forward, [4, 1],
            _cp_edata([2, 5, -4, 6, -1, -5, 3], [-1, -1, -1, 1, 1, -1],
                      [[-2, -3, -6], [1, 4]])),
        _cp_step("bps-2", :inverse, [1, 4],
            _cp_edata([1, 3, -3, 2], [-1, -1, 1, -1],
                      [[-1, -2, -4], [4]])),
        _cp_step("E1", :forward, [1, 3],
            _cp_edata([5, 1, -2, -4, 3, 2], [-1, -1, 1, -1, 1],
                      [[-3, -1, -5], [4]])),
        _cp_step("ps-IV", :inverse, [5, 1],
            _cp_edata([-1, -3, 2, 1], [-1, 1, -1], [[-2], [3]])),
        _cp_step("C1", :forward, [3, 2],
            _cp_edata([-1, 2, -4, 1], [-1, 1, 1, -1],
                      [[-3, -2], [3, 4]])),
        _cp_step("bps-2", :inverse, [3, 4],
            _cp_edata([-1, 2, -2, 1], [-1, 1])),
    ]
    right_initial = _cp_edata([-1, -2, 3, -4, 2, 4, 5, -5, 1, -3],
                              [-1, -1, 1, -1, 1])
    right_steps = OGraphCPBoundaryStep[
        _cp_step("F1", :forward, [4, 5],
            _cp_edata([-1, -2, 3, 5, 1, -3], [-1, -1, 1, -1, 1, 1],
                      [[-5, -4, -6], [4, 2, 6]])),
        _cp_step("ps-III", :inverse, [6, 4],
            _cp_edata([-1, -2, 3, 4, 1, -3], [-1, -1, 1, 1], [[-4], [2]])),
        _cp_step("A1", :forward, [3, 4],
            _cp_edata([-1, -2, 5, 1, 3, -5, 4], [-1, -1, 1, -1, 1],
                      [[-3, -4], [2]])),
        _cp_step("bps-1", :inverse, [4, 3],
            _cp_edata([-1, -2, 3, 1], [-1, -1, 1], [[-3], [2]])),
        _cp_step("C1", :forward, [2, 3],
            _cp_edata([-1, 3, -4, 1], [-1, 1, 1, -1],
                      [[-2, -3], [2, 4]])),
        _cp_step("bps-2", :inverse, [2, 4],
            _cp_edata([-1, 2, -2, 1], [-1, 1])),
    ]
    # GKS relabels only the RHS before its six-step reduction:
    # old paper labels 1,2,3,4,5 become 1,2,4,5,3.
    OGraphCPBoundaryCertificate(OGRAPH_CP_BOUNDARY_CERTIFICATE_SCHEMA,
        OGRAPH_CP_SOURCE_RULE_VERSION, paper_left, paper_right,
        [1, 2, 4, 5, 3], right_initial, left_steps, right_steps)
end

"Return the boundary-rooted GKS26 14+6 CP reduction certificate."
ograph_cp_boundary_certificate() = _gks26_cp_boundary_certificate()

function _same_cp_step(left::OGraphCPBoundaryStep, right::OGraphCPBoundaryStep)
    left.move_name == right.move_name && left.direction == right.direction &&
        left.vertices == right.vertices && _same_boundary_edata(left.target, right.target)
end

function _same_cp_certificate(left::OGraphCPBoundaryCertificate,
                              right::OGraphCPBoundaryCertificate)
    left.schema == right.schema && left.source_rule_version == right.source_rule_version &&
        _same_boundary_edata(left.paper_left, right.paper_left) &&
        _same_boundary_edata(left.paper_right, right.paper_right) &&
        left.right_relabel == right.right_relabel &&
        _same_boundary_edata(left.right_initial, right.right_initial) &&
        length(left.left_steps) == length(right.left_steps) &&
        length(left.right_steps) == length(right.right_steps) &&
        all(_same_cp_step.(left.left_steps, right.left_steps)) &&
        all(_same_cp_step.(left.right_steps, right.right_steps))
end

"Audit exact source transcription and boundary roots; move semantics remain separate."
function audit_ograph_cp_boundary_certificate(certificate::OGraphCPBoundaryCertificate)
    expected = _gks26_cp_boundary_certificate()
    certificate.schema == OGRAPH_CP_BOUNDARY_CERTIFICATE_SCHEMA ||
        return (valid=false, reason=:schema_mismatch, semantic_moves_verified=false)
    certificate.source_rule_version == OGRAPH_CP_SOURCE_RULE_VERSION ||
        return (valid=false, reason=:source_rule_mismatch, semantic_moves_verified=false)
    states = vcat([certificate.paper_left, certificate.paper_right,
                   certificate.right_initial],
                  [step.target for step in certificate.left_steps],
                  [step.target for step in certificate.right_steps])
    all(_valid_boundary_edata, states) ||
        return (valid=false, reason=:malformed_edata, semantic_moves_verified=false)
    _same_cp_certificate(certificate, expected) ||
        return (valid=false, reason=:source_transcription_mismatch,
                semantic_moves_verified=false)
    _same_boundary_edata(_cp_relabel(certificate.paper_right,
                                    certificate.right_relabel),
                         certificate.right_initial) ||
        return (valid=false, reason=:right_relabel_mismatch,
                semantic_moves_verified=false)
    left_terminal = last(certificate.left_steps).target
    right_terminal = last(certificate.right_steps).target
    _same_boundary_edata(left_terminal, right_terminal) ||
        return (valid=false, reason=:terminal_mismatch, semantic_moves_verified=false)
    max_vertices = maximum(length(state.vertex_signs) for state in states)
    multi_circuit_left = [index for (index, step) in enumerate(certificate.left_steps)
                          if 1 + length(step.target.closed_components) > 1]
    multi_circuit_right = [index for (index, step) in enumerate(certificate.right_steps)
                           if 1 + length(step.target.closed_components) > 1]
    (valid=true, reason=:source_transcription_verified,
     source_rule_version=certificate.source_rule_version,
     left_move_count=length(certificate.left_steps),
     right_move_count=length(certificate.right_steps),
     max_vertices=max_vertices,
     multi_circuit_left=multi_circuit_left,
     multi_circuit_right=multi_circuit_right,
     boundary_roots_preserved=true,
     geometric_authority_verified=false,
     certifiable=false,
     semantic_moves_verified=false)
end

function _canonicalization_witness(state::OGraphState)
    canonical = canonicalize_ograph(state)
    isempty(state.gauss) && return OGraphCanonicalizationWitness(1, Int[])
    for root in eachindex(state.gauss)
        _ograph_relabel_from_root(state, root) == canonical || continue
        labels = zeros(Int, true_vertex_count(state))
        next_label = 1
        for offset in 0:length(state.gauss)-1
            old = abs(state.gauss[mod1(root + offset, length(state.gauss))])
            if labels[old] == 0
                labels[old] = next_label
                next_label += 1
            end
        end
        return OGraphCanonicalizationWitness(root, labels)
    end
    error("canonicalization witness could not be reconstructed")
end

function _expected_pattern_candidate(candidate::MoveCandidate)
    move_type = _normalize_move_type(candidate.type)
    move_type in ("MP", "PS", "CP", "S") || return nothing
    candidate.before === nothing && return nothing
    candidate.after === nothing && return nothing
    data = make_move_data()
    move_type == "CP" && !_cp_catalogue_source_gate(data).valid && return nothing
    move_type=="S" && !(only(data["S-R"]).gauss==[[1,2,-3,-4,4,3,-2,-1]] &&
                         only(data["S-R"]).orientation==[-1,-1,-1,-1]) && return nothing
    before_table, after_table, direction = try
        _pattern_pairs(move_type, length(candidate.vertices), data)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    direction == candidate.direction || return nothing
    before = findfirst(pattern -> pattern.parameter == candidate.parameter ||
                       (move_type=="PS" && direction==:forward && isempty(pattern.parameter)), before_table)
    before === nothing && return nothing
    canonical_before = before_table[before]
    after = findfirst(pattern -> _same_move_pattern(pattern, candidate.after),
                      _after_patterns(canonical_before, after_table, candidate.parameter))
    after === nothing && return nothing
    _same_move_pattern(candidate.before, canonical_before) || return nothing
    canonical_after = _after_patterns(canonical_before, after_table,
                                      candidate.parameter)[after]
    _same_move_pattern(candidate.after, canonical_after) || return nothing
    MoveCandidate(move_type, candidate.parameter, copy(candidate.vertices),
                  NTuple{2,Int}[], direction, canonical_before, canonical_after,
                  copy(candidate.locations))
end

function _check_candidate_source(source::VirtualLink, candidate::MoveCandidate)
    source.is_weighted && return (:weighted_unsupported, nothing, nothing, nothing)
    isempty(candidate.edges) || return (:edge_selected_unsupported, nothing, nothing, nothing)
    expected = _expected_pattern_candidate(candidate)
    expected === nothing && return (:catalogue_mismatch, nothing, nothing, nothing)
    code, orientation = real_gauss_code(source)
    all(!=(0), source.orientation) && all(word -> all(!=(0), word), source.gauss) ||
        return (:non_ograph_source, nothing, nothing, nothing)
    length(code) == 1 || return (:source_not_one_circuit, nothing, nothing, nothing)
    vertices = expected.vertices
    length(vertices) == length(something(expected.before).orientation) ||
        return (:source_mismatch, nothing, nothing, nothing)
    length(unique(vertices)) == length(vertices) ||
        return (:source_mismatch, nothing, nothing, nothing)
    all(vertex -> 1 <= vertex <= length(orientation), vertices) ||
        return (:source_mismatch, nothing, nothing, nothing)
    orientation[vertices] == something(expected.before).orientation ||
        return (:source_mismatch, nothing, nothing, nothing)
    mapped_before = [_map_fragment(fragment, vertices)
                     for fragment in something(expected.before).gauss]
    length(expected.locations) == length(mapped_before) ||
        return (:source_mismatch, nothing, nothing, nothing)
    for (fragment, (component, start)) in zip(mapped_before, expected.locations)
        1 <= component <= length(code) || return (:source_mismatch, nothing, nothing, nothing)
        1 <= start <= length(code[component]) ||
            return (:source_mismatch, nothing, nothing, nothing)
        (isempty(fragment) || start in _cyclic_occurrences(code[component], fragment)) ||
            return (:source_mismatch, nothing, nothing, nothing)
    end
    _locations_compatible(expected.locations, mapped_before, code) ||
        return (:source_mismatch, nothing, nothing, nothing)
    (:ok, expected, code, orientation)
end

function _target_fragment_locations(target::VirtualLink, candidate::MoveCandidate,
                                    after_vertices::AbstractVector)
    code, _ = real_gauss_code(target)
    locations = Tuple{Int,Int}[]
    for fragment in something(candidate.after).gauss
        mapped = _map_fragment(fragment, after_vertices)
        isempty(mapped) && (push!(locations, (0, 0)); continue)
        matches = [(component, start) for (component, word) in enumerate(code)
                   for start in _cyclic_occurrences(word, mapped)]
        length(matches) == 1 || return nothing
        push!(locations, only(matches))
    end
    locations
end

function _outside_witness(source_code, target_code, candidate::MoveCandidate,
                          mapped_before, old_to_new)
    consumed = Set{Tuple{Int,Int}}()
    for (fragment, (component, start)) in zip(mapped_before, candidate.locations)
        for position in _interval_positions(start, length(fragment),
                                            length(source_code[component]))
            push!(consumed, (component, position))
        end
    end
    occurrences = OGraphOutsideOccurrence[]
    source_to_target = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
    for (component, word) in enumerate(source_code), position in eachindex(word)
        (component, position) in consumed && continue
        source_entry = word[position]
        haskey(old_to_new, abs(source_entry)) || return nothing
        target_entry = sign(source_entry) * old_to_new[abs(source_entry)]
        matches = [(target_component, target_position)
                   for (target_component, target_word) in enumerate(target_code)
                   for target_position in findall(==(target_entry), target_word)]
        length(matches) == 1 || return nothing
        target_component, target_position = only(matches)
        push!(occurrences, OGraphOutsideOccurrence(
            component, position, target_component, target_position,
            source_entry, target_entry))
        source_to_target[(component, position)] = (target_component, target_position)
    end
    adjacencies = NTuple{4,Int}[]
    for (component, word) in enumerate(source_code), position in eachindex(word)
        next_position = mod1(position + 1, length(word))
        candidate.type in ("PS","S") && candidate.direction==:forward &&
            (component,next_position) in candidate.locations && continue
        haskey(source_to_target, (component, position)) || continue
        haskey(source_to_target, (component, next_position)) || continue
        target_component, target_position = source_to_target[(component, position)]
        next_component, next_target_position = source_to_target[(component, next_position)]
        target_component == next_component || return nothing
        mod1(target_position + 1, length(target_code[target_component])) == next_target_position ||
            return nothing
        push!(adjacencies, (position, next_position, target_position, next_target_position))
    end
    occurrences, adjacencies
end

const _METHOD_B_DP = (1, 2, 6, 3, 4, 2, 1, 3, 5, 5, 4, 6)
const _METHOD_B_DM = (1, 2, 3, 3, 4, 5, 1, 6, 5, 2, 4, 6)

function _method_b_steps(source::VirtualLink)
    code, orientation = real_gauss_code(source)
    edges = real_edge_table(source)
    steps = OGraphMethodBLaneStep[]
    for edge in edges
        source_crossing, target_crossing = edge.crossing
        source_slots = (1:3) .+ 6 .+ 3 * (source_crossing > 0)
        target_slots = (1:3) .+ 3 * (target_crossing > 0)
        source_dictionary = orientation[abs(source_crossing)] > 0 ? _METHOD_B_DP : _METHOD_B_DM
        target_dictionary = orientation[abs(target_crossing)] > 0 ? _METHOD_B_DP : _METHOD_B_DM
        source_nodes = 6 * (abs(source_crossing) - 1) .+
            [source_dictionary[index] for index in source_slots]
        target_nodes = 6 * (abs(target_crossing) - 1) .+
            [target_dictionary[index] for index in target_slots]
        word = code[edge.strand]
        head = findfirst(==(target_crossing), word)
        head === nothing && error("real edge head occurrence was lost")
        tail = mod1(head - 1, length(word))
        word[tail] == source_crossing || error("real edge tail occurrence was lost")
        for lane in 1:3
            from_node, to_node = lane == 3 ?
                (target_nodes[lane], source_nodes[lane]) :
                (source_nodes[lane], target_nodes[lane])
            push!(steps, OGraphMethodBLaneStep(edge.id, edge.strand, tail, head,
                                               lane, from_node, to_node))
        end
    end
    steps
end

function _method_b_circuit_trace(steps, selected::OGraphMethodBLaneStep)
    outgoing = Dict{Int,OGraphMethodBLaneStep}()
    for step in steps
        haskey(outgoing, step.from_node) && return nothing
        outgoing[step.from_node] = step
    end
    trace = OGraphMethodBLaneStep[]
    node = selected.from_node
    while true
        haskey(outgoing, node) || return nothing
        step = outgoing[node]
        push!(trace, step)
        node = step.to_node
        node == selected.from_node && break
        length(trace) <= length(steps) || return nothing
    end
    trace
end

function _ps_target_gap_boundaries(source_code, target::VirtualLink, target_code,
                                   candidate::MoveCandidate, mapped_before, old_to_new)
    consumed = Set{Tuple{Int,Int}}()
    fragment_positions = Vector{Vector{Int}}()
    for (fragment, (component, start)) in zip(mapped_before, candidate.locations)
        positions = [mod1(start + offset, length(source_code[component]))
                     for offset in 0:length(fragment)-1]
        push!(fragment_positions, positions)
        foreach(position -> push!(consumed, (component, position)), positions)
    end
    steps = _method_b_steps(target)
    boundaries = OGraphPSGapBoundary[]
    for (index, ((component, start), positions)) in
            enumerate(zip(candidate.locations, fragment_positions))
        n = length(source_code[component])
        length(positions) < n || return nothing
        predecessor = mod1(start - 1, n)
        while (component, predecessor) in consumed
            predecessor = mod1(predecessor - 1, n)
        end
        successor = mod1(start + length(positions), n)
        while (component, successor) in consumed
            successor = mod1(successor + 1, n)
        end
        predecessor_entry = source_code[component][predecessor]
        successor_entry = source_code[component][successor]
        haskey(old_to_new, abs(predecessor_entry)) || return nothing
        haskey(old_to_new, abs(successor_entry)) || return nothing
        target_tail_entry = sign(predecessor_entry) * old_to_new[abs(predecessor_entry)]
        target_head_entry = sign(successor_entry) * old_to_new[abs(successor_entry)]
        target_matches = Tuple{Int,Int,Int}[]
        for (target_component, word) in enumerate(target_code), tail in eachindex(word)
            head = mod1(tail + 1, length(word))
            word[tail] == target_tail_entry && word[head] == target_head_entry || continue
            edge_ids = unique(step.edge_id for step in steps
                              if step.component == target_component &&
                                 step.tail_occurrence == tail &&
                                 step.head_occurrence == head)
            length(edge_ids) == 1 || return nothing
            push!(target_matches, (target_component, tail, only(edge_ids)))
        end
        length(target_matches) == 1 || return nothing
        target_component, tail, edge_id = only(target_matches)
        head = mod1(tail + 1, length(target_code[target_component]))
        push!(boundaries, OGraphPSGapBoundary(
            index, component, positions, predecessor, successor,
            predecessor_entry, successor_entry, target_component, tail, head,
            target_tail_entry, target_head_entry, edge_id))
    end
    boundaries
end

function _ps_method_b_witness(target::VirtualLink, target_state::OGraphState,
                              candidate::MoveCandidate,
                              boundaries::Vector{OGraphPSGapBoundary})
    # MST25 definition 0083 and figure 0084. Row II swaps catalogue fragments.
    rules = Dict("1" => ("PS-I",3,2), "2" => ("PS-II",3,1),
                 "3" => ("PS-III",1,2), "4" => ("PS-IV",3,3))
    haskey(rules, candidate.parameter) || return nothing
    length(boundaries) == 2 || return nothing
    steps = _method_b_steps(target)
    selected_edges = [boundary.target_edge_id for boundary in boundaries]
    candidate.parameter=="2" && reverse!(selected_edges)
    paper_type, left_lane, right_lane = rules[candidate.parameter]
    left_matches = [step for step in steps
                    if step.edge_id == selected_edges[1] && step.lane == left_lane]
    right_matches = [step for step in steps
                     if step.edge_id == selected_edges[2] && step.lane == right_lane]
    length(left_matches) == 1 && length(right_matches) == 1 || return nothing
    left, right = only(left_matches), only(right_matches)
    trace = _method_b_circuit_trace(steps, left)
    trace === nothing && return nothing
    any(step -> step.edge_id == right.edge_id && step.lane == right.lane, trace) || return nothing
    circuit_id = minimum(min(step.from_node, step.to_node) for step in trace)
    OGraphPSMethodBWitness(paper_type, candidate.parameter,
                           candidate.direction==:forward ? :source_raw : :target_raw,
                           ograph_raw_hash(target_state), copy(boundaries),
                           left.edge_id, right.edge_id, left_lane, right_lane,
                           circuit_id, trace)
end

function _ps_source_gap_boundaries(source::VirtualLink, candidate::MoveCandidate)
    code,_=real_gauss_code(source)
    steps=_method_b_steps(source)
    result=OGraphPSGapBoundary[]
    for (i,(component,head)) in enumerate(candidate.locations)
        word=code[component]; tail=mod1(head-1,length(word))
        ids=unique(s.edge_id for s in steps if s.component==component &&
                   s.tail_occurrence==tail && s.head_occurrence==head)
        length(ids)==1 || return nothing
        push!(result,OGraphPSGapBoundary(i,component,Int[],tail,head,
            word[tail],word[head],component,tail,head,word[tail],word[head],only(ids)))
    end
    result
end

"Inspect the actual PS lane circuits; unlike closed-state validation this also handles a nonclosed normal diagram."
function ps_global_condition(source::VirtualLink,candidate::MoveCandidate)
    reason,frozen,code,orientation=_check_candidate_source(source,candidate)
    reason==:ok && frozen.type=="PS" || return (valid=false,reason=reason,witness=nothing)
    lhs=source
    boundaries=if frozen.direction==:forward
        _ps_source_gap_boundaries(source,frozen)
    else
        lhs=deepcopy(source); _apply_pattern_move!(lhs,frozen)
        mapping,_,_=_vertex_maps(orientation,frozen.vertices,frozen.before,frozen.after)
        target_code,_=real_gauss_code(lhs)
        _ps_target_gap_boundaries(code,lhs,target_code,frozen,
            [_map_fragment(f,frozen.vertices) for f in frozen.before.gauss],mapping)
    end
    boundaries===nothing && return (valid=false,reason=:ps_gap_transport_failed,witness=nothing)
    witness=_ps_method_b_witness(lhs,ograph_state(lhs),frozen,boundaries)
    (valid=witness!==nothing,reason=witness===nothing ? :different_method_b_circuits : :same_method_b_circuit,
     witness=witness)
end

function _coerce_expected_state(value)
    value === nothing && return nothing
    value isa OGraphState && return value
    value isa VirtualLink && return ograph_state(value)
    throw(ArgumentError("expected endpoint must be an OGraphState or VirtualLink"))
end

function _cp_source_vertex_maps(orientation, candidate::MoveCandidate)
    before = something(candidate.before)
    after = something(candidate.after)
    selected = candidate.vertices
    old_count = length(orientation)
    if candidate.direction == :forward
        length(selected) == 3 && before.vertex_count == 3 && after.vertex_count == 5 ||
            throw(ArgumentError("source CP forward roles must be 3-to-5"))
        old_to_new = Dict(index => index for index in 1:old_count)
        after_vertices = vcat(selected, collect(old_count + 1:old_count + 2))
        new_orientation = vcat(copy(orientation), zeros(Int, 2))
    elseif candidate.direction == :inverse
        length(selected) == 5 && before.vertex_count == 5 && after.vertex_count == 3 ||
            throw(ArgumentError("source CP inverse roles must be 5-to-3"))
        # GKS26 labels the five local roles.  The inverse retains roles 1--3
        # and removes roles 4--5; numeric magnitude is not a geometric role.
        removed = Set(selected[4:5])
        remaining = [vertex for vertex in 1:old_count if !(vertex in removed)]
        old_to_new = Dict(vertex => index for (index, vertex) in enumerate(remaining))
        after_vertices = [old_to_new[vertex] for vertex in selected[1:3]]
        new_orientation = orientation[remaining]
    else
        throw(ArgumentError("unsupported CP direction $(candidate.direction)"))
    end
    for local_vertex in eachindex(after.orientation)
        new_orientation[after_vertices[local_vertex]] = after.orientation[local_vertex]
    end
    old_to_new, after_vertices, new_orientation
end

function _apply_cp_source_rule(source::VirtualLink, candidate::MoveCandidate)
    before = something(candidate.before)
    after = something(candidate.after)
    code, orientation = real_gauss_code(source)
    code = [_clean_component(component) for component in code]
    weights = [zeros(Float64, length(component)) for component in code]
    old_to_new, after_vertices, new_orientation =
        _cp_source_vertex_maps(orientation, candidate)
    mapped_before = [_map_fragment(fragment, candidate.vertices) for fragment in before.gauss]
    mapped_after = [_map_fragment(fragment, after_vertices) for fragment in after.gauss]
    converted = length(before.gauss) == length(after.gauss) &&
                before.capacity == after.capacity ?
        _simple_pattern_rewrite(code, weights, candidate, mapped_before, mapped_after,
                                old_to_new, new_orientation; weighted=false) : nothing
    converted === nothing &&
        (converted = _reconnect_pattern(code, weights, candidate, mapped_before, mapped_after,
                                        old_to_new, new_orientation; weighted=false))
    result, _, result_orientation = converted
    target = deepcopy(source)
    set_data!(target; rgauss=result, orientation=result_orientation)
    target, old_to_new, after_vertices
end

"""
    audit_ograph_cp_candidate(source, candidate)

Independently apply the source-transcribed CP role rule and compare it with the
production catalogue rewrite. A successful structural comparison is certifying
only when the separate `N(G, eps)=(-G,-eps)` recipe/naturality gate also passes.
"""
function _audit_ograph_cp_candidate(source::VirtualLink, candidate::MoveCandidate)
    reason, frozen, _, _ = _check_candidate_source(source, candidate)
    reason == :ok || return (
        structural_checks_pass=false,
        geometric_authority_verified=false,
        certifiable=false,
        reason=reason,
        source_rule_target=nothing,
        production_target=nothing,
        source_to_target_vertices=Int[],
        error=nothing,
    )
    frozen.type == "CP" || return (
        structural_checks_pass=false,
        geometric_authority_verified=false,
        certifiable=false,
        reason=:not_cp_candidate,
        source_rule_target=nothing,
        production_target=nothing,
        source_to_target_vertices=Int[],
        error=nothing,
    )
    cp_result = try
        _apply_cp_source_rule(source, frozen)
    catch err
        return (
            structural_checks_pass=false,
            geometric_authority_verified=false,
            certifiable=false,
            reason=:cp_source_rule_application_failed,
            source_rule_target=nothing,
            production_target=nothing,
            source_to_target_vertices=Int[],
            error=sprint(showerror, err),
        )
    end
    source_rule_link, old_to_new, _ = cp_result
    source_rule_target = try
        ograph_state(source_rule_link)
    catch err
        return (
            structural_checks_pass=false,
            geometric_authority_verified=false,
            certifiable=false,
            reason=:cp_source_rule_target_invalid,
            source_rule_target=nothing,
            production_target=nothing,
            source_to_target_vertices=Int[],
            error=sprint(showerror, err),
        )
    end
    production_link = deepcopy(source)
    try
        _apply_pattern_move!(production_link, frozen)
    catch err
        err isa InterruptException && rethrow()
        return (
            structural_checks_pass=false,
            geometric_authority_verified=false,
            certifiable=false,
            reason=:raw_application_failed,
            source_rule_target=source_rule_target,
            production_target=nothing,
            source_to_target_vertices=Int[],
            error=sprint(showerror, err),
        )
    end
    production_target = try
        ograph_state(production_link)
    catch err
        return (
            structural_checks_pass=false,
            geometric_authority_verified=false,
            certifiable=false,
            reason=:invalid_raw_target,
            source_rule_target=source_rule_target,
            production_target=nothing,
            source_to_target_vertices=Int[],
            error=sprint(showerror, err),
        )
    end
    vertex_map = [get(old_to_new, vertex, 0)
                  for vertex in 1:true_vertex_count(ograph_state(source))]
    production_target == source_rule_target || return (
        structural_checks_pass=false,
        geometric_authority_verified=false,
        certifiable=false,
        reason=:cp_production_semantics_mismatch,
        source_rule_target=source_rule_target,
        production_target=production_target,
        source_to_target_vertices=vertex_map,
        error=nothing,
    )
    authority = _cp_semantic_authority_gate()
    (
        structural_checks_pass=true,
        geometric_authority_verified=authority.valid,
        certifiable=authority.certifiable,
        reason=authority.valid ? :cp_source_geometry_verified :
                                 :cp_geometric_authority_unverified,
        source_rule_target=source_rule_target,
        production_target=production_target,
        source_to_target_vertices=vertex_map,
        error=nothing,
    )
end

function audit_ograph_cp_candidate(source::VirtualLink, candidate::MoveCandidate)
    try
        _audit_ograph_cp_candidate(source, candidate)
    catch err
        err isa InterruptException && rethrow()
        (
            structural_checks_pass=false,
            geometric_authority_verified=false,
            certifiable=false,
            reason=:malformed_audit_input,
            source_rule_target=nothing,
            production_target=nothing,
            source_to_target_vertices=Int[],
            error=sprint(showerror, err),
        )
    end
end

function _verify_ograph_move(source::VirtualLink, candidate::MoveCandidate;
                             expected_raw_target=nothing,
                             expected_canonical_target=nothing,legacy_mp_labels=false)
    reason, frozen, source_code, source_orientation =
        _check_candidate_source(source, candidate)
    reason == :ok || return _ograph_move_invalid(reason)
    source_state = try ograph_state(source) catch err
        return _ograph_move_invalid(:invalid_source_state; error=err)
    end
    source_closedness = try closedness_witness(source_state) catch err
        return _ograph_move_invalid(:source_closedness_failed; error=err)
    end
    all(source_closedness.criteria) || return _ograph_move_invalid(:source_not_closed)

    before = something(frozen.before)
    after = something(frozen.after)
    cp_source_target = nothing
    if frozen.type == "CP"
        cp_result = try _apply_cp_source_rule(source, frozen) catch err
            return _ograph_move_invalid(:cp_source_rule_application_failed; error=err)
        end
        cp_source_target, old_to_new, after_vertices = cp_result
    else
        old_to_new, after_vertices, _ =
            _vertex_maps(source_orientation, frozen.vertices, before, after;
                         role_ordered=frozen.type=="MP" && !legacy_mp_labels)
    end
    mapped_before = [_map_fragment(fragment, frozen.vertices) for fragment in before.gauss]
    target = deepcopy(source)
    try
        _apply_pattern_move!(target, frozen;legacy_mp_labels=legacy_mp_labels)
    catch err
        err isa InterruptException && rethrow()
        return _ograph_move_invalid(:raw_application_failed; error=err)
    end
    if frozen.type == "CP"
        production_state = try ograph_state(target) catch err
            return _ograph_move_invalid(:invalid_raw_target; error=err)
        end
        source_rule_state = try ograph_state(something(cp_source_target)) catch err
            return _ograph_move_invalid(:cp_source_rule_target_invalid; error=err)
        end
        production_state == source_rule_state ||
            return _ograph_move_invalid(:cp_production_semantics_mismatch)
        authority = try _cp_semantic_authority_gate() catch err
            return _ograph_move_invalid(:cp_geometric_authority_check_failed; error=err)
        end
        authority.valid ||
            return _ograph_move_invalid(:cp_geometric_authority_unverified)
    end
    target_state = try ograph_state(target) catch err
        return _ograph_move_invalid(:invalid_raw_target; error=err)
    end
    target_closedness = try closedness_witness(target_state) catch err
        return _ograph_move_invalid(:target_closedness_failed; error=err)
    end
    all(target_closedness.criteria) || return _ograph_move_invalid(:target_not_closed)

    raw_expected = try _coerce_expected_state(expected_raw_target) catch err
        return _ograph_move_invalid(:invalid_expected_raw_target; error=err)
    end
    raw_expected === nothing || raw_expected == target_state ||
        return _ograph_move_invalid(:raw_target_mismatch)
    canonical_target = canonicalize_ograph(target_state)
    canonical_expected = try _coerce_expected_state(expected_canonical_target) catch err
        return _ograph_move_invalid(:invalid_expected_canonical_target; error=err)
    end
    canonical_expected === nothing || canonicalize_ograph(canonical_expected) == canonical_target ||
        return _ograph_move_invalid(:canonical_target_mismatch)

    target_locations = _target_fragment_locations(target, frozen, after_vertices)
    target_locations === nothing && return _ograph_move_invalid(:target_patch_mismatch)
    target_code, _ = real_gauss_code(target)
    outside = _outside_witness(source_code, target_code, frozen,
                               mapped_before, old_to_new)
    outside === nothing && return _ograph_move_invalid(:outside_region_mismatch)
    outside_occurrences, outside_adjacencies = outside
    vertex_map = [get(old_to_new, vertex, 0) for vertex in eachindex(source_orientation)]

    gap_boundaries = frozen.type == "PS" ?
        (frozen.direction==:forward ? _ps_source_gap_boundaries(source,frozen) :
         _ps_target_gap_boundaries(source_code, target, target_code, frozen,
                                  mapped_before, old_to_new)) : nothing
    frozen.type == "PS" && gap_boundaries === nothing &&
        return _ograph_move_invalid(:ps_target_gap_transport_failed)
    method_b = frozen.type == "PS" ?
        (frozen.direction==:forward ? _ps_method_b_witness(source,source_state,frozen,gap_boundaries) :
         _ps_method_b_witness(target, target_state, frozen, gap_boundaries)) : nothing
    frozen.type == "PS" && method_b === nothing &&
        return _ograph_move_invalid(:ps_method_b_condition_failed)

    source_canonical = canonicalize_ograph(source_state)
    witness = OGraphMoveWitness(
        OGRAPH_MOVE_WITNESS_SCHEMA, OGRAPH_MOVE_VERIFIER_VERSION,
        frozen.type, frozen.parameter, frozen.direction, copy(frozen.vertices),
        copy(frozen.locations), target_locations, source_state, target_state,
        source_canonical, canonical_target, ograph_raw_hash(source_state),
        ograph_raw_hash(target_state), source_closedness, target_closedness,
        vertex_map, outside_occurrences, outside_adjacencies,
        _canonicalization_witness(source_state),
        _canonicalization_witness(target_state), method_b)
    OGraphMoveVerification(true, :verified, witness, nothing)
end

"""Verify a frozen MP/PS/CP candidate against its actual raw source and target.

The move generator is not consulted. Catalogue row, concrete source fragments,
locations, raw rewrite, outside occurrences, C1--C3 endpoints, and presentation
canonicalization are recomputed. Malformed evidence fails closed as a report;
only process interrupts are rethrown. PS-I--IV require the prescribed
Method-B lane connection on the pre-insertion state. CP additionally requires the GKS source
catalogue correspondence, the production-SOT recipe/naturality audit, and the
independent role-wise rewrite before it enters the same closedness, outside-
region, canonical-target, persistence, and replay checks as MP/PS.
"""
function verify_ograph_move(source::VirtualLink, candidate::MoveCandidate;
                            expected_raw_target=nothing,
                            expected_canonical_target=nothing)
    try
        _verify_ograph_move(source, candidate;
            expected_raw_target=expected_raw_target,
            expected_canonical_target=expected_canonical_target)
    catch err
        err isa InterruptException && rethrow()
        _ograph_move_invalid(:malformed_verification_input; error=err)
    end
end

function _link_from_ograph_state(state::OGraphState)
    link = VirtualLink()
    set_data!(link; rgauss=[collect(state.gauss)],
              orientation=collect(state.vertex_signs))
end

function _candidate_from_witness(witness::OGraphMoveWitness)
    data = make_move_data()
    before_table, after_table, direction =
        _pattern_pairs(witness.move_type, length(witness.vertices), data)
    before = only(pattern for pattern in before_table
                  if pattern.parameter == witness.parameter ||
                     (witness.move_type=="PS" && direction==:forward && isempty(pattern.parameter)))
    after = only(_after_patterns(before, after_table, witness.parameter))
    MoveCandidate(witness.move_type, witness.parameter, copy(witness.vertices),
                  NTuple{2,Int}[], direction, before, after,
                  copy(witness.source_locations))
end

function ograph_move_witness_data(witness::OGraphMoveWitness)
    state_data(state) = (gauss=collect(state.gauss),
                         vertex_signs=collect(state.vertex_signs),
                         convention_id=state.convention_id)
    closedness_data(value) = (canonical_hash=value.canonical_hash,
        verifier_version=value.verifier_version, criteria=collect(value.criteria),
        c1_oriented_circuits=value.c1_oriented_circuits,
        c2_trivalent_components=value.c2_trivalent_components,
        c3_method_circuits=value.c3_method_circuits,
        true_vertices=value.true_vertices,
        method_circuit_lengths=collect(value.method_circuit_lengths))
    outside_data(value) = (source_component=value.source_component,
        source_position=value.source_position, target_component=value.target_component,
        target_position=value.target_position, source_entry=value.source_entry,
        target_entry=value.target_entry)
    canonicalization_data(value) = (root_occurrence=value.root_occurrence,
        raw_to_canonical_labels=copy(value.raw_to_canonical_labels))
    step_data(value) = (edge_id=value.edge_id, component=value.component,
        tail_occurrence=value.tail_occurrence, head_occurrence=value.head_occurrence,
        lane=value.lane, from_node=value.from_node, to_node=value.to_node)
    gap_data(value) = (fragment_index=value.fragment_index,
        source_component=value.source_component,
        source_fragment_positions=copy(value.source_fragment_positions),
        source_predecessor_position=value.source_predecessor_position,
        source_successor_position=value.source_successor_position,
        source_predecessor_entry=value.source_predecessor_entry,
        source_successor_entry=value.source_successor_entry,
        target_component=value.target_component,
        target_tail_occurrence=value.target_tail_occurrence,
        target_head_occurrence=value.target_head_occurrence,
        target_tail_entry=value.target_tail_entry,
        target_head_entry=value.target_head_entry,
        target_edge_id=value.target_edge_id)
    method_data(value) = value === nothing ? nothing : (
        paper_type=value.paper_type, parameter=value.parameter,
        evaluated_state_side=String(value.evaluated_state_side),
        evaluated_state_sha256=value.evaluated_state_sha256,
        gap_boundaries=gap_data.(value.gap_boundaries),
        left_edge_id=value.left_edge_id, right_edge_id=value.right_edge_id,
        left_lane=value.left_lane, right_lane=value.right_lane,
        circuit_id=value.circuit_id, circuit_trace=step_data.(value.circuit_trace))
    (schema=witness.schema, verifier_version=witness.verifier_version,
     move_type=witness.move_type, parameter=witness.parameter,
     direction=String(witness.direction), vertices=copy(witness.vertices),
     source_locations=collect.(witness.source_locations),
     target_locations=collect.(witness.target_locations),
     source_raw=state_data(witness.source_raw), target_raw=state_data(witness.target_raw),
     source_canonical=state_data(witness.source_canonical),
     target_canonical=state_data(witness.target_canonical),
     source_raw_sha256=witness.source_raw_sha256,
     target_raw_sha256=witness.target_raw_sha256,
     source_closedness=closedness_data(witness.source_closedness),
     target_closedness=closedness_data(witness.target_closedness),
     source_to_target_vertices=copy(witness.source_to_target_vertices),
     outside_occurrences=outside_data.(witness.outside_occurrences),
     outside_adjacencies=collect.(witness.outside_adjacencies),
     source_canonicalization=canonicalization_data(witness.source_canonicalization),
     target_canonicalization=canonicalization_data(witness.target_canonicalization),
     method_b=method_data(witness.method_b))
end

_move_data_keys(value) = Set(String(key) for key in keys(value))
function _move_data_field(value, key::String)
    value isa NamedTuple && return getproperty(value, Symbol(key))
    value isa AbstractDict || throw(ArgumentError("move witness object expected"))
    haskey(value, key) && return value[key]
    haskey(value, Symbol(key)) && return value[Symbol(key)]
    throw(ArgumentError("move witness field $key is missing"))
end

function _move_data_object(value, fields)
    value isa Union{NamedTuple,AbstractDict} ||
        throw(ArgumentError("move witness object expected"))
    actual = _move_data_keys(value)
    expected = Set(String.(fields))
    actual == expected || throw(ArgumentError(
        "move witness fields differ: missing=$(sort!(collect(setdiff(expected, actual)))) " *
        "extra=$(sort!(collect(setdiff(actual, expected))))"))
    value
end

_move_data_string(value, label) = value isa AbstractString ? String(value) :
    throw(ArgumentError("$label must be a string"))
_move_data_int(value, label) = value isa Integer && !(value isa Bool) ? Int(value) :
    throw(ArgumentError("$label must be an integer"))
_move_data_bool(value, label) = value isa Bool ? value :
    throw(ArgumentError("$label must be Boolean"))
function _move_data_vector(value, length_, item, label)
    value isa AbstractVector || throw(ArgumentError("$label must be an array"))
    length_ === nothing || length(value) == length_ ||
        throw(ArgumentError("$label has the wrong length"))
    [item(entry, "$label[$index]") for (index, entry) in enumerate(value)]
end
_move_data_int_vector(value, label; length_=nothing) =
    _move_data_vector(value, length_, _move_data_int, label)

function _move_data_state(value, label)
    fields = ("gauss", "vertex_signs", "convention_id")
    _move_data_object(value, fields)
    OGraphState(
        _move_data_int_vector(_move_data_field(value, "gauss"), "$label.gauss"),
        _move_data_int_vector(_move_data_field(value, "vertex_signs"),
                              "$label.vertex_signs");
        convention_id=_move_data_string(_move_data_field(value, "convention_id"),
                                        "$label.convention_id"))
end

function _move_data_closedness(value, label)
    fields = ("canonical_hash", "verifier_version", "criteria",
              "c1_oriented_circuits", "c2_trivalent_components",
              "c3_method_circuits", "true_vertices", "method_circuit_lengths")
    _move_data_object(value, fields)
    criteria = _move_data_vector(_move_data_field(value, "criteria"), 3,
                                 _move_data_bool, "$label.criteria")
    lengths = _move_data_int_vector(_move_data_field(value, "method_circuit_lengths"),
                                    "$label.method_circuit_lengths")
    ClosednessWitness(
        _move_data_string(_move_data_field(value, "canonical_hash"),
                          "$label.canonical_hash"),
        _move_data_string(_move_data_field(value, "verifier_version"),
                          "$label.verifier_version"), Tuple(criteria),
        _move_data_int(_move_data_field(value, "c1_oriented_circuits"),
                       "$label.c1_oriented_circuits"),
        _move_data_int(_move_data_field(value, "c2_trivalent_components"),
                       "$label.c2_trivalent_components"),
        _move_data_int(_move_data_field(value, "c3_method_circuits"),
                       "$label.c3_method_circuits"),
        _move_data_int(_move_data_field(value, "true_vertices"),
                       "$label.true_vertices"), Tuple(lengths))
end

function _move_data_outside(value, label)
    fields = ("source_component", "source_position", "target_component",
              "target_position", "source_entry", "target_entry")
    _move_data_object(value, fields)
    values = [_move_data_int(_move_data_field(value, field), "$label.$field")
              for field in fields]
    OGraphOutsideOccurrence(values...)
end

function _move_data_canonicalization(value, label)
    _move_data_object(value, ("root_occurrence", "raw_to_canonical_labels"))
    OGraphCanonicalizationWitness(
        _move_data_int(_move_data_field(value, "root_occurrence"),
                       "$label.root_occurrence"),
        _move_data_int_vector(_move_data_field(value, "raw_to_canonical_labels"),
                              "$label.raw_to_canonical_labels"))
end

function _move_data_step(value, label)
    fields = ("edge_id", "component", "tail_occurrence", "head_occurrence",
              "lane", "from_node", "to_node")
    _move_data_object(value, fields)
    values = [_move_data_int(_move_data_field(value, field), "$label.$field")
              for field in fields]
    OGraphMethodBLaneStep(values...)
end

function _move_data_gap(value, label)
    fields = ("fragment_index", "source_component", "source_fragment_positions",
              "source_predecessor_position", "source_successor_position",
              "source_predecessor_entry", "source_successor_entry",
              "target_component", "target_tail_occurrence", "target_head_occurrence",
              "target_tail_entry", "target_head_entry", "target_edge_id")
    _move_data_object(value, fields)
    OGraphPSGapBoundary(
        _move_data_int(_move_data_field(value, "fragment_index"),
                       "$label.fragment_index"),
        _move_data_int(_move_data_field(value, "source_component"),
                       "$label.source_component"),
        _move_data_int_vector(_move_data_field(value, "source_fragment_positions"),
                              "$label.source_fragment_positions"),
        (_move_data_int(_move_data_field(value, field), "$label.$field")
         for field in fields[4:end])...)
end

function _move_data_method(value, label)
    value === nothing && return nothing
    fields = ("paper_type", "parameter", "evaluated_state_side",
              "evaluated_state_sha256", "gap_boundaries", "left_edge_id",
              "right_edge_id", "left_lane", "right_lane", "circuit_id",
              "circuit_trace")
    _move_data_object(value, fields)
    side_string = _move_data_string(_move_data_field(value, "evaluated_state_side"),
                                    "$label.evaluated_state_side")
    side_string in ("target_raw","source_raw") ||
        throw(ArgumentError("$label.evaluated_state_side is unsupported"))
    gaps = _move_data_vector(_move_data_field(value, "gap_boundaries"), 2,
        (entry, item_label) -> _move_data_gap(entry, item_label),
        "$label.gap_boundaries")
    trace = _move_data_vector(_move_data_field(value, "circuit_trace"), nothing,
        (entry, item_label) -> _move_data_step(entry, item_label),
        "$label.circuit_trace")
    OGraphPSMethodBWitness(
        _move_data_string(_move_data_field(value, "paper_type"), "$label.paper_type"),
        _move_data_string(_move_data_field(value, "parameter"), "$label.parameter"),
        Symbol(side_string),
        _move_data_string(_move_data_field(value, "evaluated_state_sha256"),
                          "$label.evaluated_state_sha256"),
        gaps,
        (_move_data_int(_move_data_field(value, field), "$label.$field")
         for field in ("left_edge_id", "right_edge_id", "left_lane", "right_lane",
                       "circuit_id"))..., trace)
end

"Strictly decode canonical JSON-like move evidence into a typed witness."
function ograph_move_witness_from_data(value)
    fields = ("schema", "verifier_version", "move_type", "parameter", "direction",
              "vertices", "source_locations", "target_locations", "source_raw",
              "target_raw", "source_canonical", "target_canonical",
              "source_raw_sha256", "target_raw_sha256", "source_closedness",
              "target_closedness", "source_to_target_vertices",
              "outside_occurrences", "outside_adjacencies",
              "source_canonicalization", "target_canonicalization", "method_b")
    _move_data_object(value, fields)
    schema = _move_data_string(_move_data_field(value, "schema"), "schema")
    schema == OGRAPH_MOVE_WITNESS_SCHEMA ||
        throw(ArgumentError("unsupported move witness schema $schema"))
    verifier = _move_data_string(_move_data_field(value, "verifier_version"),
                                 "verifier_version")
    verifier in (OGRAPH_MOVE_VERIFIER_VERSION,"ograph-move-source-replay/2") ||
        throw(ArgumentError("unsupported move verifier version $verifier"))
    direction_string = _move_data_string(_move_data_field(value, "direction"), "direction")
    direction_string in ("forward", "inverse") ||
        throw(ArgumentError("unsupported move direction $direction_string"))
    pair_vector(field, label) = _move_data_vector(_move_data_field(value, field), nothing,
        (entry, item_label) -> Tuple(_move_data_int_vector(entry, item_label; length_=2)), label)
    outside = _move_data_vector(_move_data_field(value, "outside_occurrences"), nothing,
        (entry, label) -> _move_data_outside(entry, label), "outside_occurrences")
    adjacencies = _move_data_vector(_move_data_field(value, "outside_adjacencies"), nothing,
        (entry, label) -> Tuple(_move_data_int_vector(entry, label; length_=4)),
        "outside_adjacencies")
    OGraphMoveWitness(
        schema, verifier,
        _move_data_string(_move_data_field(value, "move_type"), "move_type"),
        _move_data_string(_move_data_field(value, "parameter"), "parameter"),
        Symbol(direction_string),
        _move_data_int_vector(_move_data_field(value, "vertices"), "vertices"),
        pair_vector("source_locations", "source_locations"),
        pair_vector("target_locations", "target_locations"),
        _move_data_state(_move_data_field(value, "source_raw"), "source_raw"),
        _move_data_state(_move_data_field(value, "target_raw"), "target_raw"),
        _move_data_state(_move_data_field(value, "source_canonical"), "source_canonical"),
        _move_data_state(_move_data_field(value, "target_canonical"), "target_canonical"),
        _move_data_string(_move_data_field(value, "source_raw_sha256"), "source_raw_sha256"),
        _move_data_string(_move_data_field(value, "target_raw_sha256"), "target_raw_sha256"),
        _move_data_closedness(_move_data_field(value, "source_closedness"),
                              "source_closedness"),
        _move_data_closedness(_move_data_field(value, "target_closedness"),
                              "target_closedness"),
        _move_data_int_vector(_move_data_field(value, "source_to_target_vertices"),
                              "source_to_target_vertices"),
        outside, adjacencies,
        _move_data_canonicalization(_move_data_field(value, "source_canonicalization"),
                                    "source_canonicalization"),
        _move_data_canonicalization(_move_data_field(value, "target_canonicalization"),
                                    "target_canonicalization"),
        _move_data_method(_move_data_field(value, "method_b"), "method_b"))
end

"Recompute a stored witness without using `movable` or a stored success flag."
function replay_ograph_move(witness::OGraphMoveWitness)
    source = _link_from_ograph_state(witness.source_raw)
    candidate = try _candidate_from_witness(witness) catch err
        return _ograph_move_invalid(:stored_candidate_invalid; error=err)
    end
    replay = verify_ograph_move(source, candidate;
        expected_raw_target=witness.target_raw,
        expected_canonical_target=witness.target_canonical)
    if witness.verifier_version=="ograph-move-source-replay/2" &&
       witness.move_type=="MP" && witness.parameter=="D2"
        # v2 numbered retained patch roles by sorted source labels. First
        # verify the current source-role move gives the SAME canonical target;
        # only then reproduce the v2 raw numbering for exact historical replay.
        modern=verify_ograph_move(source,candidate;
            expected_canonical_target=witness.target_canonical)
        modern.valid || return modern
        replay=_verify_ograph_move(source,candidate;expected_raw_target=witness.target_raw,
            expected_canonical_target=witness.target_canonical,legacy_mp_labels=true)
    end
    replay.valid || return replay
    merge(ograph_move_witness_data(something(replay.witness)),
          (verifier_version=witness.verifier_version,)) == ograph_move_witness_data(witness) ||
        return _ograph_move_invalid(:witness_mismatch)
    replay
end


function replay_ograph_move(value::Union{NamedTuple,AbstractDict})
    witness = try ograph_move_witness_from_data(value) catch err
        return _ograph_move_invalid(:stored_witness_invalid; error=err)
    end
    replay_ograph_move(witness)
end
