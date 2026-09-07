# Source scope, index gates, and replay contract:
# docs/specs/ograph-move-verifier.md
const OGRAPH_MOVE_WITNESS_SCHEMA = "akaghef.ograph-move-witness/1"
const OGRAPH_MOVE_VERIFIER_VERSION = "ograph-move-source-replay/1"

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
    move_type in ("MP", "PS") || return nothing
    candidate.before === nothing && return nothing
    candidate.after === nothing && return nothing
    data = make_move_data()
    before_table, after_table, direction = try
        _pattern_pairs(move_type, length(candidate.vertices), data)
    catch err
        err isa ArgumentError || rethrow()
        return nothing
    end
    direction == candidate.direction || return nothing
    before = findfirst(pattern -> pattern.parameter == candidate.parameter, before_table)
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
        start in _cyclic_occurrences(code[component], fragment) ||
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
    # Parameter 3 is the only current template row independently overlaid on
    # MST25 v3 Fig. 2.15. Rows 1, 2, and 4 remain fail-closed below.
    rules = Dict("3" => ("PS-III", 1, 2))
    haskey(rules, candidate.parameter) || return nothing
    before = something(candidate.before)
    length(before.gauss) == 2 && all(fragment -> length(fragment) == 2, before.gauss) ||
        return nothing
    length(boundaries) == 2 || return nothing
    steps = _method_b_steps(target)
    selected_edges = [boundary.target_edge_id for boundary in boundaries]
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
    OGraphPSMethodBWitness(paper_type, candidate.parameter, :target_raw,
                           ograph_raw_hash(target_state), copy(boundaries),
                           left.edge_id, right.edge_id, left_lane, right_lane,
                           circuit_id, trace)
end

function _coerce_expected_state(value)
    value === nothing && return nothing
    value isa OGraphState && return value
    value isa VirtualLink && return ograph_state(value)
    throw(ArgumentError("expected endpoint must be an OGraphState or VirtualLink"))
end

function _verify_ograph_move(source::VirtualLink, candidate::MoveCandidate;
                             expected_raw_target=nothing,
                             expected_canonical_target=nothing)
    reason, frozen, source_code, source_orientation =
        _check_candidate_source(source, candidate)
    reason == :ok || return _ograph_move_invalid(reason)
    frozen.type == "MP" && frozen.parameter != "D2" &&
        return _ograph_move_invalid(:mp_parameter_not_source_audited)
    frozen.type == "PS" && frozen.parameter != "3" &&
        return _ograph_move_invalid(:ps_parameter_not_source_audited)
    source_state = try ograph_state(source) catch err
        return _ograph_move_invalid(:invalid_source_state; error=err)
    end
    source_closedness = try closedness_witness(source_state) catch err
        return _ograph_move_invalid(:source_closedness_failed; error=err)
    end
    all(source_closedness.criteria) || return _ograph_move_invalid(:source_not_closed)

    before = something(frozen.before)
    after = something(frozen.after)
    old_to_new, after_vertices, _ =
        _vertex_maps(source_orientation, frozen.vertices, before, after)
    mapped_before = [_map_fragment(fragment, frozen.vertices) for fragment in before.gauss]
    target = deepcopy(source)
    try
        _apply_pattern_move!(target, frozen)
    catch err
        err isa InterruptException && rethrow()
        return _ograph_move_invalid(:raw_application_failed; error=err)
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
        _ps_target_gap_boundaries(source_code, target, target_code, frozen,
                                  mapped_before, old_to_new) : nothing
    frozen.type == "PS" && gap_boundaries === nothing &&
        return _ograph_move_invalid(:ps_target_gap_transport_failed)
    method_b = frozen.type == "PS" ?
        _ps_method_b_witness(target, target_state, frozen, gap_boundaries) : nothing
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

"""Verify a frozen MP/PS candidate against its actual raw source and target.

The move generator is not consulted. Catalogue row, concrete source fragments,
locations, raw rewrite, outside occurrences, C1--C3 endpoints, and presentation
canonicalization are recomputed. Malformed evidence fails closed as a report;
only process interrupts are rethrown. PS is limited to the independently
overlaid parameter-3/PS-III row.
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
                  if pattern.parameter == witness.parameter)
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
    side_string == "target_raw" ||
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
    verifier == OGRAPH_MOVE_VERIFIER_VERSION ||
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
    replay.valid || return replay
    ograph_move_witness_data(something(replay.witness)) ==
        ograph_move_witness_data(witness) ||
        return _ograph_move_invalid(:witness_mismatch)
    replay
end


function replay_ograph_move(value::Union{NamedTuple,AbstractDict})
    witness = try ograph_move_witness_from_data(value) catch err
        return _ograph_move_invalid(:stored_witness_invalid; error=err)
    end
    replay_ograph_move(witness)
end
