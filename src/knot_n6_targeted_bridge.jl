"""Targeted barrier-7 bridge for the unresolved barrier-6 components.

This is intentionally not an N<=7 census.  The input is one in-memory
`KnotN6Classification` result.  Only states in its unresolved components are
seeded, and only insertion targets with at most seven crossings are added.
At seven crossings the frontier is expanded by R1/R2 deletion and R3.  A
target at most six crossings must already be in the input state index; a
seven-crossing target is checked by the Carter predicate before it is added.

The barrier-6 graph is reused as follows: every state in an unresolved
barrier-6 component is a seed, so no N<=6 edge needs to be rediscovered.
Insertion candidates whose target is still at most six crossings are checked
against the existing edge graph and counted as reused.  A seven-crossing
state is represented once by its canonical identity.  The transient closure
graph has unresolved-component nodes, canonical seven-crossing nodes, and
known-component boundary nodes.  Known nodes are never expanded.  Thus every
boundary candidate is inspected, while full paths and parent arrays are
discarded before the result is returned.
"""

if !isdefined(@__MODULE__, :KnotN6Classification)
    include(joinpath(@__DIR__, "knot_n6_classification.jl"))
end

using SHA

const KNOT_N6_TARGETED_BRIDGE_VERSION = 1
const KNOT_N6_TARGETED_BRIDGE_BARRIER = 7

"Compact transient result of one targeted barrier-7 closure." 
struct KnotN6TargetedBridgeResult
    barrier::Int
    input_fingerprint::String
    seed_components::Vector{Int}
    seed_state_count::Int
    n7_state_count::Int
    candidate_count::Int
    adopted_transition_count::Int
    duplicate_transition_count::Int
    reused_n6_candidate_count::Int
    excluded_over_barrier_count::Int
    n7_processed_count::Int
    max_frontier::Int
    transition_verification_failures::Int
    n6_membership_failures::Int
    classicality_checks::Int
    classicality_failures::Int
    connected_known_labels::Dict{Int,Vector{String}}
    known_reach_distance::Dict{Int,Union{Nothing,Int}}
    resolved_components::Dict{Int,String}
    resolved_component_type_counts::Dict{String,Int}
    resolved_state_type_counts::Dict{String,Int}
    closure_groups::Vector{Vector{Int}}
    unresolved_components::Vector{Int}
    fingerprint::String
    elapsed_seconds::Float64
end

Base.length(result::KnotN6TargetedBridgeResult) = result.n7_state_count

mutable struct _N6TBUnionFind
    parent::Vector{Int}
    size::Vector{Int}
end

_N6TBUnionFind() = _N6TBUnionFind(Int[], Int[])

function _n6tb_new_node!(uf::_N6TBUnionFind,
                         adjacency::Vector{Vector{Int}})
    node = length(uf.parent) + 1
    push!(uf.parent, node)
    push!(uf.size, 1)
    push!(adjacency, Int[])
    node
end

function _n6tb_find(uf::_N6TBUnionFind, node::Int)
    root = node
    while uf.parent[root] != root
        root = uf.parent[root]
    end
    while uf.parent[node] != node
        next = uf.parent[node]
        uf.parent[node] = root
        node = next
    end
    root
end

function _n6tb_union!(uf::_N6TBUnionFind, first::Int, second::Int)
    left = _n6tb_find(uf, first)
    right = _n6tb_find(uf, second)
    left == right && return left
    uf.size[left] < uf.size[right] && ((left, right) = (right, left))
    uf.parent[right] = left
    uf.size[left] += uf.size[right]
    left
end

mutable struct _N6TBStats
    candidate_count::Int
    adopted_transition_count::Int
    duplicate_transition_count::Int
    reused_n6_candidate_count::Int
    excluded_over_barrier_count::Int
    n7_processed_count::Int
    max_frontier::Int
    transition_verification_failures::Int
    n6_membership_failures::Int
    classicality_checks::Int
    classicality_failures::Int
    transition_keys::Set{Tuple{Symbol,Int,Int,Symbol}}
    graph_keys::Set{Tuple{Int,Int}}
end

function _N6TBStats()
    _N6TBStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               Set{Tuple{Symbol,Int,Int,Symbol}}(), Set{Tuple{Int,Int}}())
end

function _n6tb_known_label_name(label::AbstractString)
    first(split(String(label), '^'; limit=2))
end

"Reject both distinct knot names and an accidental strict mirror merge." 
function _n6tb_check_known_conflict(uf::_N6TBUnionFind,
                                    known_nodes::Dict{Int,Int},
                                    classification::KnotN6Classification)
    entries = sort!(collect(known_nodes); by=first)
    for left_index in eachindex(entries)
        left_component, left_node = entries[left_index]
        left_label = classification.component_type_labels[left_component]
        for right_index in left_index+1:length(entries)
            right_component, right_node = entries[right_index]
            _n6tb_find(uf, left_node) == _n6tb_find(uf, right_node) || continue
            right_label = classification.component_type_labels[right_component]
            left_label == right_label && continue
            left_name = _n6tb_known_label_name(left_label)
            right_name = _n6tb_known_label_name(right_label)
            if left_name != right_name
                throw(ErrorException(
                    "barrier-7 closure joins distinct known knot names " *
                    "$left_label and $right_label"))
            end
            throw(ErrorException(
                "barrier-7 closure joins strict mirror labels " *
                "$left_label and $right_label"))
        end
    end
    nothing
end

function _n6tb_add_transition!(stats::_N6TBStats,
                               uf::_N6TBUnionFind,
                               adjacency::Vector{Vector{Int}},
                               known_nodes::Dict{Int,Int},
                               classification::KnotN6Classification,
                               key::Tuple{Symbol,Int,Int,Symbol},
                               first_node::Int, second_node::Int)
    if key in stats.transition_keys
        stats.duplicate_transition_count += 1
        return false
    end
    push!(stats.transition_keys, key)
    stats.adopted_transition_count += 1
    pair = first_node <= second_node ? (first_node, second_node) :
        (second_node, first_node)
    if !(pair in stats.graph_keys)
        push!(stats.graph_keys, pair)
        push!(adjacency[first_node], second_node)
        first_node == second_node || push!(adjacency[second_node], first_node)
        _n6tb_union!(uf, first_node, second_node)
        _n6tb_check_known_conflict(uf, known_nodes, classification)
    end
    true
end

function _n6tb_known_node!(component::Int,
                           known_nodes::Dict{Int,Int},
                           uf::_N6TBUnionFind,
                           adjacency::Vector{Vector{Int}})
    get!(known_nodes, component) do
        _n6tb_new_node!(uf, adjacency)
    end
end

function _n6tb_existing_edge(classification::KnotN6Classification,
                             source::Int, target::Int, move_type::Symbol)
    any(classification.adjacency[source]) do edge_id
        edge = classification.edges[edge_id]
        edge.move_type == move_type &&
            ((edge.source == source && edge.target == target) ||
             (edge.source == target && edge.target == source))
    end
end

function _n6tb_checked_apply(source::KnotGaussDiagram,
                             candidate::ReidemeisterCandidate,
                             stats::_N6TBStats,
                             context::AbstractString)
    raw_target = try
        apply_reidemeister(source, candidate)
    catch error
        stats.transition_verification_failures += 1
        throw(ErrorException("barrier-7 application failed ($context): " *
                             sprint(showerror, error)))
    end
    if !verify_reidemeister_transition(source, raw_target, candidate)
        stats.transition_verification_failures += 1
        throw(ErrorException("barrier-7 verifier rejected transition ($context)"))
    end
    raw_target
end

function _n6tb_ensure_n7!(diagram::KnotGaussDiagram,
                          n7_index::Dict{KnotGaussDiagram,Int},
                          n7_states::Vector{KnotGaussDiagram},
                          n7_nodes::Vector{Int},
                          n7_local_by_node::Dict{Int,Int},
                          n7_parent::Vector{Int},
                          n7_distance::Vector{Int},
                          uf::_N6TBUnionFind,
                          adjacency::Vector{Vector{Int}})
    canonical = canonicalize(diagram)
    existing = get(n7_index, canonical, 0)
    existing != 0 && return n7_nodes[existing], false
    node = _n6tb_new_node!(uf, adjacency)
    push!(n7_states, canonical)
    n7_index[canonical] = length(n7_states)
    push!(n7_nodes, node)
    n7_local_by_node[node] = length(n7_states)
    push!(n7_parent, 0)
    push!(n7_distance, -1)
    node, true
end

function _n6tb_classical_n7!(diagram::KnotGaussDiagram, stats::_N6TBStats)
    stats.classicality_checks += 1
    if !classical_gauss_realizable_fast(diagram)
        stats.classicality_failures += 1
        throw(ErrorException("barrier-7 frontier generated a nonclassical state"))
    end
    nothing
end

function _n6tb_component_connections(seed_components::Vector{Int},
                                     seed_nodes::Dict{Int,Int},
                                     known_nodes::Dict{Int,Int},
                                     uf::_N6TBUnionFind,
                                     classification::KnotN6Classification)
    connections = Dict{Int,Vector{String}}()
    for component in seed_components
        root = _n6tb_find(uf, seed_nodes[component])
        labels = String[]
        for (known_component, node) in sort!(collect(known_nodes); by=first)
            _n6tb_find(uf, node) == root || continue
            push!(labels, classification.component_type_labels[known_component])
        end
        connections[component] = sort!(unique(labels))
    end
    connections
end

function _n6tb_shortest_known_distance(seed_node::Int,
                                       known_nodes::Dict{Int,Int},
                                       adjacency::Vector{Vector{Int}})
    isempty(known_nodes) && return nothing
    distances = fill(-1, length(adjacency))
    queue = Int[seed_node]
    distances[seed_node] = 0
    head = 1
    best = nothing
    while head <= length(queue)
        current = queue[head]
        head += 1
        best !== nothing && distances[current] >= best && continue
        for (component, node) in known_nodes
            node == current || continue
            best = best === nothing ? distances[current] :
                min(best, distances[current])
        end
        for neighbor in adjacency[current]
            distances[neighbor] == -1 || continue
            distances[neighbor] = distances[current] + 1
            push!(queue, neighbor)
        end
    end
    best
end

function _n6tb_closure_groups(seed_components::Vector{Int},
                              seed_nodes::Dict{Int,Int},
                              uf::_N6TBUnionFind)
    groups = Dict{Int,Vector{Int}}()
    for component in seed_components
        root = _n6tb_find(uf, seed_nodes[component])
        push!(get!(groups, root, Int[]), component)
    end
    sort!([sort!(members) for members in values(groups)]; by=first)
end

function _n6tb_fingerprint(classification::KnotN6Classification,
                           n7_states::Vector{KnotGaussDiagram},
                           graph_keys::Set{Tuple{Int,Int}},
                           connections::Dict{Int,Vector{String}})
    lines = String[
        "KnotN6TargetedBridge;version=$(KNOT_N6_TARGETED_BRIDGE_VERSION);" *
        "barrier=$(KNOT_N6_TARGETED_BRIDGE_BARRIER);input=$(classification.fingerprint)"]
    append!(lines, "state7|" * canonical_serialization(diagram)
                    for diagram in sort!(copy(n7_states), by=canonical_serialization))
    for pair in sort!(collect(graph_keys))
        push!(lines, "graph|$(pair[1])|$(pair[2])")
    end
    for component in sort!(collect(keys(connections)))
        push!(lines, "connection|$component|$(join(connections[component], ','))")
    end
    bytes2hex(SHA.sha256(codeunits(join(lines, '\n'))))
end

"""Run the unresolved-component-only barrier-7 closure.

`classification` must be produced by `build_knot_n6_classification`; a
smaller result is accepted for a fast fixture and naturally has no N7 seeds.
The production script requires `classification.max_n == 6`.
"""
function run_knot_n6_targeted_bridge(
        classification::KnotN6Classification;
        barrier::Integer=KNOT_N6_TARGETED_BRIDGE_BARRIER,
        progress::Bool=true,
        _workspace_callback=nothing)
    Int(barrier) == KNOT_N6_TARGETED_BRIDGE_BARRIER ||
        throw(ArgumentError("targeted bridge is fixed to barrier=7"))
    0 <= classification.max_n <= KNOT_N6_MAX_N ||
        throw(ArgumentError("classification max_n is outside the N6 lane"))
    started = time_ns()
    seed_components = sort!(copy(classification.unresolved_components))
    seed_component_set = Set(seed_components)
    seed_state_indices = Int[
        index for index in eachindex(classification.states)
        if classification.components[index] in seed_component_set]
    sort!(seed_state_indices)

    uf = _N6TBUnionFind()
    closure_adjacency = Vector{Vector{Int}}()
    seed_nodes = Dict{Int,Int}()
    for component in seed_components
        seed_nodes[component] = _n6tb_new_node!(uf, closure_adjacency)
    end
    known_nodes = Dict{Int,Int}()
    stats = _N6TBStats()
    n7_index = Dict{KnotGaussDiagram,Int}()
    n7_states = KnotGaussDiagram[]
    n7_nodes = Int[]
    n7_local_by_node = Dict{Int,Int}()
    n7_parent = Int[]
    n7_distance = Int[]
    queue = Int[]
    queue_head = 1

    # Seed only insertion targets that can lie at the requested barrier.  A
    # target at most six crossings is already an N6 edge; checking it here
    # makes the reuse contract explicit without adding a duplicate bridge edge.
    for source_index in seed_state_indices
        source = classification.states[source_index]
        source_component = classification.components[source_index]
        for move_name in ("R1", "R2")
            move_type = Symbol(move_name)
            candidates = movable(source, move_name; param=:insert)
            for candidate in candidates
                stats.candidate_count += 1
                delta = move_type == :R1 ? 1 : 2
                target_crossings = crossing_count(source) + delta
                if target_crossings > KNOT_N6_TARGETED_BRIDGE_BARRIER
                    stats.excluded_over_barrier_count += 1
                    continue
                end
                raw_target = _n6tb_checked_apply(
                    source, candidate, stats,
                    "seed state=$source_index move=$move_name")
                target_crossings == crossing_count(raw_target) || throw(ErrorException(
                    "barrier-7 insertion changed crossing count unexpectedly"))
                canonical_target = canonicalize(raw_target)
                if target_crossings <= classification.max_n
                    target_index = get(classification.state_index, canonical_target, 0)
                    if target_index == 0
                        stats.n6_membership_failures += 1
                        throw(ErrorException(
                            "reused N<=6 insertion target is absent from input census"))
                    end
                    classification.components[target_index] == source_component ||
                        throw(ErrorException(
                            "reused N<=6 insertion crosses barrier-6 components"))
                    _n6tb_existing_edge(classification, source_index, target_index,
                                         move_type) || throw(ErrorException(
                        "existing barrier-6 edge was not found for reusable insertion"))
                    stats.reused_n6_candidate_count += 1
                    continue
                end
                target_crossings == 7 || throw(ErrorException(
                    "targeted bridge retained a target above seven crossings"))
                node, is_new = _n6tb_ensure_n7!(
                    canonical_target, n7_index, n7_states, n7_nodes,
                    n7_local_by_node, n7_parent, n7_distance, uf,
                    closure_adjacency)
                is_new && _n6tb_classical_n7!(canonical_target, stats)
                key = (:seed_to_n7, source_index, node, move_type)
                _n6tb_add_transition!(stats, uf, closure_adjacency, known_nodes,
                                      classification, key, seed_nodes[source_component], node)
                if is_new
                    push!(queue, node)
                    local_index = n7_local_by_node[node]
                    n7_parent[local_index] = seed_nodes[source_component]
                    n7_distance[local_index] = 1
                    stats.max_frontier = max(stats.max_frontier,
                                             length(queue) - queue_head + 1)
                end
            end
        end
    end

    # Expand only seven-crossing frontier states.  At this boundary all
    # deletion/R3 candidates are processed, including candidates that already
    # touch a known component; known N6 vertices themselves are never expanded.
    while queue_head <= length(queue)
        current_node = queue[queue_head]
        queue_head += 1
        local_index = get(n7_local_by_node, current_node, 0)
        local_index == 0 && error("internal N7 node-index inconsistency")
        source = n7_states[local_index]
        stats.n7_processed_count += 1
        for (move_name, direction) in (("R1", :delete), ("R2", :delete),
                                       ("R3", nothing))
            candidates = direction === nothing ? movable(source, move_name) :
                movable(source, move_name; param=direction)
            for candidate in candidates
                stats.candidate_count += 1
                raw_target = _n6tb_checked_apply(
                    source, candidate, stats,
                    "N7 node=$current_node move=$move_name")
                target_crossings = crossing_count(raw_target)
                target_crossings <= 7 || throw(ErrorException(
                    "barrier-7 deletion/R3 target exceeds seven crossings"))
                canonical_target = canonicalize(raw_target)
                if target_crossings <= classification.max_n
                    target_index = get(classification.state_index, canonical_target, 0)
                    if target_index == 0
                        stats.n6_membership_failures += 1
                        throw(ErrorException(
                            "N<=6 target from N7 frontier is absent from input census"))
                    end
                    target_component = classification.components[target_index]
                    target_node = if haskey(seed_nodes, target_component)
                        seed_nodes[target_component]
                    elseif haskey(classification.component_type_labels, target_component)
                        _n6tb_known_node!(target_component, known_nodes, uf,
                                          closure_adjacency)
                    else
                        stats.n6_membership_failures += 1
                        throw(ErrorException(
                            "N<=6 target belongs to an unclassified component"))
                    end
                    key = (:n7_to_n6, current_node, target_index,
                           Symbol(move_name))
                    _n6tb_add_transition!(stats, uf, closure_adjacency, known_nodes,
                                          classification, key, current_node, target_node)
                    continue
                end
                target_crossings == 7 || throw(ErrorException(
                    "frontier retained a target above seven crossings"))
                target_node, is_new = _n6tb_ensure_n7!(
                    canonical_target, n7_index, n7_states, n7_nodes,
                    n7_local_by_node, n7_parent, n7_distance, uf,
                    closure_adjacency)
                is_new && _n6tb_classical_n7!(canonical_target, stats)
                low, high = min(current_node, target_node),
                            max(current_node, target_node)
                key = (:n7_to_n7, low, high, Symbol(move_name))
                _n6tb_add_transition!(stats, uf, closure_adjacency, known_nodes,
                                      classification, key, current_node, target_node)
                if is_new
                    push!(queue, target_node)
                    target_local_index = n7_local_by_node[target_node]
                    n7_parent[target_local_index] = current_node
                    n7_distance[target_local_index] =
                        n7_distance[local_index] < 0 ? 1 : n7_distance[local_index] + 1
                    stats.max_frontier = max(stats.max_frontier,
                                             length(queue) - queue_head + 1)
                end
            end
        end
        progress && stats.n7_processed_count % 100 == 0 && println(stderr,
            "[knot-n6-bridge] processed=$(stats.n7_processed_count) " *
            "n7=$(length(n7_states)) candidates=$(stats.candidate_count) " *
            "adopted=$(stats.adopted_transition_count) frontier=" *
            "$(length(queue) - queue_head + 1)")
    end

    connections = _n6tb_component_connections(
        seed_components, seed_nodes, known_nodes, uf, classification)
    resolved = Dict{Int,String}()
    unresolved = Int[]
    for component in seed_components
        labels = connections[component]
        if isempty(labels)
            push!(unresolved, component)
        else
            length(labels) == 1 || throw(ErrorException(
                "barrier-7 closure has multiple known labels for component $component"))
            resolved[component] = only(labels)
        end
    end
    resolved_component_type_counts = Dict{String,Int}()
    resolved_state_type_counts = Dict{String,Int}()
    for (component, label) in resolved
        resolved_component_type_counts[label] =
            get(resolved_component_type_counts, label, 0) + 1
        resolved_state_type_counts[label] =
            get(resolved_state_type_counts, label, 0) +
            count(==(component), classification.components)
    end
    distances = Dict{Int,Union{Nothing,Int}}(
        component => _n6tb_shortest_known_distance(
            seed_nodes[component], known_nodes, closure_adjacency)
        for component in seed_components)
    groups = _n6tb_closure_groups(seed_components, seed_nodes, uf)
    fingerprint = _n6tb_fingerprint(classification, n7_states,
                                     stats.graph_keys, connections)
    elapsed = (time_ns() - started) / 1e9
    if _workspace_callback !== nothing
        # Internal bridge8 consumers receive only the transient graph needed
        # to seed the next barrier.  Full candidates, paths, and parent arrays
        # remain private to this call and are not exposed by the public result.
        _workspace_callback((
            seed_components=copy(seed_components),
            unresolved_components=copy(unresolved),
            seed_nodes=copy(seed_nodes),
            known_nodes=copy(known_nodes),
            n7_states=copy(n7_states),
            n7_nodes=copy(n7_nodes),
            n7_local_by_node=copy(n7_local_by_node),
            closure_adjacency=[copy(neighbors) for neighbors in closure_adjacency],
            uf_parent=copy(uf.parent),
            uf_size=copy(uf.size),
            connections=Dict(component => copy(labels)
                             for (component, labels) in connections),
            component_type_labels=copy(classification.component_type_labels)))
    end
    # Parent arrays, queue, adjacency, state-index maps, and full candidate
    # witnesses are deliberately local temporaries and are not returned.
    KnotN6TargetedBridgeResult(
        KNOT_N6_TARGETED_BRIDGE_BARRIER, classification.fingerprint,
        seed_components, length(seed_state_indices), length(n7_states),
        stats.candidate_count, stats.adopted_transition_count,
        stats.duplicate_transition_count, stats.reused_n6_candidate_count,
        stats.excluded_over_barrier_count, stats.n7_processed_count,
        stats.max_frontier, stats.transition_verification_failures,
        stats.n6_membership_failures, stats.classicality_checks,
        stats.classicality_failures, connections, distances, resolved,
        resolved_component_type_counts, resolved_state_type_counts, groups,
        unresolved, fingerprint, elapsed)
end

build_knot_n6_targeted_bridge(args...; kwargs...) =
    run_knot_n6_targeted_bridge(args...; kwargs...)

targeted_knot_n6_bridge(args...; kwargs...) =
    run_knot_n6_targeted_bridge(args...; kwargs...)

"Return only compact closure statistics; no paths or frontier states." 
function knot_n6_targeted_bridge_summary(result::KnotN6TargetedBridgeResult)
    Dict{String,Any}(
        "version" => KNOT_N6_TARGETED_BRIDGE_VERSION,
        "barrier" => result.barrier,
        "input_fingerprint" => result.input_fingerprint,
        "seed_components" => copy(result.seed_components),
        "seed_state_count" => result.seed_state_count,
        "n7_state_count" => result.n7_state_count,
        "candidate_count" => result.candidate_count,
        "adopted_transition_count" => result.adopted_transition_count,
        "duplicate_transition_count" => result.duplicate_transition_count,
        "reused_n6_candidate_count" => result.reused_n6_candidate_count,
        "excluded_over_barrier_count" => result.excluded_over_barrier_count,
        "n7_processed_count" => result.n7_processed_count,
        "max_frontier" => result.max_frontier,
        "transition_verification_failures" => result.transition_verification_failures,
        "n6_membership_failures" => result.n6_membership_failures,
        "classicality_checks" => result.classicality_checks,
        "classicality_failures" => result.classicality_failures,
        "connected_known_labels" => Dict(string(component) => copy(labels)
                                          for (component, labels) in
                                          sort!(collect(result.connected_known_labels);
                                                by=first)),
        "known_reach_distance" => Dict(string(component) => distance
                                        for (component, distance) in
                                        sort!(collect(result.known_reach_distance);
                                              by=first)),
        "resolved_components" => Dict(string(component) => label
                                       for (component, label) in
                                       sort!(collect(result.resolved_components);
                                             by=first)),
        "resolved_component_type_counts" =>
            copy(result.resolved_component_type_counts),
        "resolved_state_type_counts" => copy(result.resolved_state_type_counts),
        "closure_groups" => [copy(group) for group in result.closure_groups],
        "unresolved_components" => copy(result.unresolved_components),
        "fingerprint_sha256" => result.fingerprint,
        "elapsed_seconds" => result.elapsed_seconds)
end
