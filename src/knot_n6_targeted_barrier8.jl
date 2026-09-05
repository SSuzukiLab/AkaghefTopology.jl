"""Targeted barrier-8 closure for the unresolved barrier-7 components.

The production entry point reconstructs the N<=6 classification once and the
barrier-7 bridge once, capturing only the bridge's transient canonical N7
states and closure graph.  It then seeds the ten components still unresolved
at barrier 7.  Existing N<=6 edges and the complete N7 bridge closure are
imported; only N8 states and transitions are discovered on demand.

For a seed at most six crossings, every R1/R2 insertion whose target has at
most eight crossings is considered.  At seven crossings, R1/R2 insertions to
N8 and all R1/R2 deletions/R3 moves are considered.  At eight crossings, only
R1/R2 deletions and R3 are considered.  Targets above eight are excluded.
Every adopted transition is applied and independently checked by
`verify_reidemeister_transition`; every new N7/N8 canonical state is checked
by the Carter realizability predicate.  N<=6 targets must be present in the
original 77989-state index.  Known components are boundary nodes only and
are never expanded.

No path, parent array, candidate witness, Carter certificate, or invariant is
returned or persisted.  A resource limit produces a compact
`resource_limited` result instead of silently declaring the remaining seeds
unresolved.
"""

if !isdefined(@__MODULE__, :KnotN6TargetedBridgeResult)
    include(joinpath(@__DIR__, "knot_n6_targeted_bridge.jl"))
end

using SHA

const KNOT_N6_TARGETED_BARRIER8_VERSION = 1
const KNOT_N6_TARGETED_BARRIER8 = 8
const KNOT_N6_TARGETED_BARRIER8_STATE_LIMIT = 500_000
const KNOT_N6_TARGETED_BARRIER8_TIME_LIMIT_NS = UInt64(1_800_000_000_000)
const KNOT_N6_TARGETED_BARRIER8_RSS_LIMIT = UInt64(12) * UInt64(1024)^3

"Compact result of one targeted barrier-8 closure."
struct KnotN6TargetedBarrier8Result
    barrier::Int
    n6_input_fingerprint::String
    barrier7_fingerprint::String
    seed_components::Vector{Int}
    seed_state_count::Int
    imported_n7_state_count::Int
    n7_state_count::Int
    n8_state_count::Int
    candidate_count::Int
    adopted_transition_count::Int
    duplicate_transition_count::Int
    imported_graph_edge_count::Int
    reused_n6_candidate_count::Int
    excluded_over_barrier_count::Int
    n7_processed_count::Int
    n8_processed_count::Int
    max_frontier::Int
    transition_verification_failures::Int
    n6_membership_failures::Int
    classicality_checks::Int
    classicality_failures::Int
    connected_known_labels::Dict{Int,Vector{String}}
    known_reach_distance::Dict{Int,Union{Nothing,Int}}
    resolved_components::Dict{Int,String}
    component_status::Dict{Int,String}
    resolved_component_type_counts::Dict{String,Int}
    resolved_state_type_counts::Dict{String,Int}
    closure_groups::Vector{Vector{Int}}
    unresolved_components::Vector{Int}
    status::String
    resource_reason::String
    fingerprint::String
    elapsed_seconds::Float64
end

Base.length(result::KnotN6TargetedBarrier8Result) =
    result.n7_state_count + result.n8_state_count

mutable struct _N6B8UnionFind
    parent::Vector{Int}
    size::Vector{Int}
end

_N6B8UnionFind() = _N6B8UnionFind(Int[], Int[])

function _n6b8_new_node!(uf::_N6B8UnionFind,
                         adjacency::Vector{Vector{Int}})
    node = length(uf.parent) + 1
    push!(uf.parent, node)
    push!(uf.size, 1)
    push!(adjacency, Int[])
    node
end

function _n6b8_find(uf::_N6B8UnionFind, node::Int)
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

function _n6b8_union!(uf::_N6B8UnionFind, first::Int, second::Int)
    left = _n6b8_find(uf, first)
    right = _n6b8_find(uf, second)
    left == right && return left
    uf.size[left] < uf.size[right] && ((left, right) = (right, left))
    uf.parent[right] = left
    uf.size[left] += uf.size[right]
    left
end

function _n6b8_old_find(parent::Vector{Int}, node::Int)
    root = node
    while parent[root] != root
        root = parent[root]
    end
    root
end

mutable struct _N6B8Stats
    candidate_count::Int
    adopted_transition_count::Int
    duplicate_transition_count::Int
    imported_graph_edge_count::Int
    reused_n6_candidate_count::Int
    excluded_over_barrier_count::Int
    n7_processed_count::Int
    n8_processed_count::Int
    max_frontier::Int
    transition_verification_failures::Int
    n6_membership_failures::Int
    classicality_checks::Int
    classicality_failures::Int
    transition_keys::Set{Tuple{Symbol,Int,Int,Symbol}}
    graph_keys::Set{Tuple{Int,Int}}
    resource_limited::Bool
    resource_reason::String
end

function _N6B8Stats()
    _N6B8Stats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               Set{Tuple{Symbol,Int,Int,Symbol}}(), Set{Tuple{Int,Int}}(),
               false, "")
end

struct _N6B8ResourceLimit <: Exception
    reason::String
end

Base.showerror(io::IO, error::_N6B8ResourceLimit) =
    print(io, "targeted barrier-8 resource limit: ", error.reason)

function _n6b8_check_limits!(stats::_N6B8Stats, started::UInt64,
                             n7_count::Int, n8_count::Int)
    total = n7_count + n8_count
    total > KNOT_N6_TARGETED_BARRIER8_STATE_LIMIT &&
        throw(_N6B8ResourceLimit("canonical N7+N8 state limit reached"))
    elapsed = UInt64(time_ns()) - started
    elapsed >= KNOT_N6_TARGETED_BARRIER8_TIME_LIMIT_NS &&
        throw(_N6B8ResourceLimit("30-minute barrier-8 time limit reached"))
    rss = try
        UInt64(Sys.maxrss())
    catch
        UInt64(0)
    end
    rss >= KNOT_N6_TARGETED_BARRIER8_RSS_LIMIT &&
        throw(_N6B8ResourceLimit("RSS limit reached"))
    nothing
end

function _n6b8_check_new_state_limit!(stats::_N6B8Stats, started::UInt64,
                                      n7_count::Int, n8_count::Int)
    n7_count + n8_count + 1 > KNOT_N6_TARGETED_BARRIER8_STATE_LIMIT &&
        throw(_N6B8ResourceLimit("canonical N7+N8 state limit would be exceeded"))
    _n6b8_check_limits!(stats, started, n7_count, n8_count)
end

function _n6b8_known_label_name(label::AbstractString)
    first(split(String(label), '^'; limit=2))
end

function _n6b8_check_known_conflict(uf::_N6B8UnionFind,
                                    known_nodes::Dict{Int,Int},
                                    classification::KnotN6Classification)
    entries = sort!(collect(known_nodes); by=first)
    for left_index in eachindex(entries)
        left_component, left_node = entries[left_index]
        left_label = classification.component_type_labels[left_component]
        for right_index in left_index+1:length(entries)
            right_component, right_node = entries[right_index]
            _n6b8_find(uf, left_node) == _n6b8_find(uf, right_node) || continue
            right_label = classification.component_type_labels[right_component]
            left_label == right_label && continue
            if _n6b8_known_label_name(left_label) !=
               _n6b8_known_label_name(right_label)
                throw(ErrorException(
                    "barrier-8 closure joins distinct known knot names " *
                    "$left_label and $right_label"))
            end
            throw(ErrorException(
                "barrier-8 closure joins strict mirror labels " *
                "$left_label and $right_label"))
        end
    end
    nothing
end

function _n6b8_add_graph_edge!(stats::_N6B8Stats,
                               uf::_N6B8UnionFind,
                               adjacency::Vector{Vector{Int}},
                               known_nodes::Dict{Int,Int},
                               classification::KnotN6Classification,
                               first_node::Int, second_node::Int;
                               imported::Bool=false)
    pair = first_node <= second_node ? (first_node, second_node) :
        (second_node, first_node)
    pair in stats.graph_keys && return false
    push!(stats.graph_keys, pair)
    push!(adjacency[first_node], second_node)
    first_node == second_node || push!(adjacency[second_node], first_node)
    _n6b8_union!(uf, first_node, second_node)
    imported && (stats.imported_graph_edge_count += 1)
    _n6b8_check_known_conflict(uf, known_nodes, classification)
    true
end

function _n6b8_add_transition!(stats::_N6B8Stats,
                               uf::_N6B8UnionFind,
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
    _n6b8_add_graph_edge!(stats, uf, adjacency, known_nodes, classification,
                          first_node, second_node)
    true
end

function _n6b8_checked_apply(source::KnotGaussDiagram,
                             candidate::ReidemeisterCandidate,
                             stats::_N6B8Stats,
                             context::AbstractString)
    raw_target = try
        apply_reidemeister(source, candidate)
    catch error
        stats.transition_verification_failures += 1
        throw(ErrorException("barrier-8 application failed ($context): " *
                             sprint(showerror, error)))
    end
    if !verify_reidemeister_transition(source, raw_target, candidate)
        stats.transition_verification_failures += 1
        throw(ErrorException("barrier-8 verifier rejected transition ($context)"))
    end
    raw_target
end

function _n6b8_new_known_node!(component::Int,
                               known_nodes::Dict{Int,Int},
                               uf::_N6B8UnionFind,
                               adjacency::Vector{Vector{Int}})
    get!(known_nodes, component) do
        _n6b8_new_node!(uf, adjacency)
    end
end

function _n6b8_ensure_n7!(diagram::KnotGaussDiagram,
                          n7_index::Dict{KnotGaussDiagram,Int},
                          n7_states::Vector{KnotGaussDiagram},
                          n7_nodes::Vector{Int},
                          n7_local_by_node::Dict{Int,Int},
                          uf::_N6B8UnionFind,
                          adjacency::Vector{Vector{Int}})
    canonical = canonicalize(diagram)
    existing = get(n7_index, canonical, 0)
    existing != 0 && return n7_nodes[existing], false
    node = _n6b8_new_node!(uf, adjacency)
    push!(n7_states, canonical)
    n7_index[canonical] = length(n7_states)
    push!(n7_nodes, node)
    n7_local_by_node[node] = length(n7_states)
    node, true
end

function _n6b8_ensure_n8!(diagram::KnotGaussDiagram,
                          n8_index::Dict{KnotGaussDiagram,Int},
                          n8_states::Vector{KnotGaussDiagram},
                          n8_nodes::Vector{Int},
                          n8_local_by_node::Dict{Int,Int},
                          uf::_N6B8UnionFind,
                          adjacency::Vector{Vector{Int}})
    canonical = canonicalize(diagram)
    existing = get(n8_index, canonical, 0)
    existing != 0 && return n8_nodes[existing], false
    node = _n6b8_new_node!(uf, adjacency)
    push!(n8_states, canonical)
    n8_index[canonical] = length(n8_states)
    push!(n8_nodes, node)
    n8_local_by_node[node] = length(n8_states)
    node, true
end

function _n6b8_classical_new!(diagram::KnotGaussDiagram, stats::_N6B8Stats)
    stats.classicality_checks += 1
    if !classical_gauss_realizable_fast(diagram)
        stats.classicality_failures += 1
        throw(ErrorException("barrier-8 frontier generated a nonclassical state"))
    end
    nothing
end

function _n6b8_connect_known_labels!(n7_node::Int,
                                     labels::Vector{String},
                                     known_component_by_label::Dict{String,Int},
                                     known_nodes::Dict{Int,Int},
                                     uf::_N6B8UnionFind,
                                     adjacency::Vector{Vector{Int}},
                                     stats::_N6B8Stats,
                                     classification::KnotN6Classification)
    for label in labels
        component = get(known_component_by_label, label, 0)
        component == 0 && throw(ErrorException(
            "barrier-8 boundary refers to an unknown label $label"))
        known_node = _n6b8_new_known_node!(component, known_nodes, uf, adjacency)
        _n6b8_add_graph_edge!(stats, uf, adjacency, known_nodes, classification,
                              n7_node, known_node)
    end
    nothing
end

function _n6b8_seed_state_indices(classification::KnotN6Classification,
                                  seed_components::Vector{Int})
    wanted = Set(seed_components)
    indices = Int[index for index in eachindex(classification.states)
                  if classification.components[index] in wanted]
    sort!(indices)
end

function _n6b8_existing_edge(classification::KnotN6Classification,
                             source::Int, target::Int, move_type::Symbol)
    any(classification.adjacency[source]) do edge_id
        edge = classification.edges[edge_id]
        edge.move_type == move_type &&
            ((edge.source == source && edge.target == target) ||
             (edge.source == target && edge.target == source))
    end
end

function _n6b8_component_connections(seed_components::Vector{Int},
                                     seed_nodes::Dict{Int,Int},
                                     known_nodes::Dict{Int,Int},
                                     uf::_N6B8UnionFind,
                                     classification::KnotN6Classification)
    result = Dict{Int,Vector{String}}()
    for component in seed_components
        root = _n6b8_find(uf, seed_nodes[component])
        labels = String[]
        for (known_component, node) in sort!(collect(known_nodes); by=first)
            _n6b8_find(uf, node) == root || continue
            push!(labels, classification.component_type_labels[known_component])
        end
        result[component] = sort!(unique(labels))
    end
    result
end

function _n6b8_shortest_known_distance(seed_node::Int,
                                       known_nodes::Dict{Int,Int},
                                       adjacency::Vector{Vector{Int}})
    isempty(known_nodes) && return nothing
    distance = fill(-1, length(adjacency))
    queue = Int[seed_node]
    distance[seed_node] = 0
    head = 1
    best = nothing
    while head <= length(queue)
        current = queue[head]
        head += 1
        best !== nothing && distance[current] >= best && continue
        any(node == current for node in values(known_nodes)) &&
            (best = best === nothing ? distance[current] :
             min(best, distance[current]))
        for neighbor in adjacency[current]
            distance[neighbor] == -1 || continue
            distance[neighbor] = distance[current] + 1
            push!(queue, neighbor)
        end
    end
    best
end

function _n6b8_closure_groups(seed_components::Vector{Int},
                              seed_nodes::Dict{Int,Int},
                              uf::_N6B8UnionFind)
    groups = Dict{Int,Vector{Int}}()
    for component in seed_components
        root = _n6b8_find(uf, seed_nodes[component])
        push!(get!(groups, root, Int[]), component)
    end
    sort!([sort!(members) for members in values(groups)]; by=first)
end

function _n6b8_fingerprint(classification::KnotN6Classification,
                           bridge7::KnotN6TargetedBridgeResult,
                           n7_states::Vector{KnotGaussDiagram},
                           n8_states::Vector{KnotGaussDiagram},
                           graph_keys::Set{Tuple{Int,Int}},
                           component_status::Dict{Int,String})
    lines = String[
        "KnotN6TargetedBarrier8;version=$(KNOT_N6_TARGETED_BARRIER8_VERSION);" *
        "barrier=8;n6=$(classification.fingerprint);barrier7=$(bridge7.fingerprint)"]
    append!(lines, "state7|" * canonical_serialization(diagram)
                    for diagram in sort!(copy(n7_states), by=canonical_serialization))
    append!(lines, "state8|" * canonical_serialization(diagram)
                    for diagram in sort!(copy(n8_states), by=canonical_serialization))
    for pair in sort!(collect(graph_keys))
        push!(lines, "graph|$(pair[1])|$(pair[2])")
    end
    for component in sort!(collect(keys(component_status)))
        push!(lines, "component|$component|$(component_status[component])")
    end
    bytes2hex(SHA.sha256(codeunits(join(lines, '\n'))))
end

function _n6b8_import_bridge7_workspace!(
        classification::KnotN6Classification,
        bridge7::KnotN6TargetedBridgeResult,
        workspace,
        seed_components::Vector{Int},
        uf::_N6B8UnionFind,
        adjacency::Vector{Vector{Int}},
        seed_nodes::Dict{Int,Int},
        n7_index::Dict{KnotGaussDiagram,Int},
        n7_states::Vector{KnotGaussDiagram},
        n7_nodes::Vector{Int},
        n7_local_by_node::Dict{Int,Int},
        n7_boundary_labels::Dict{KnotGaussDiagram,Vector{String}},
        stats::_N6B8Stats,
        queue_kind::Vector{UInt8},
        queue_nodes::Vector{Int})
    workspace === nothing && throw(ErrorException(
        "barrier-7 workspace callback did not return a workspace"))
    Set(seed_components) == Set(bridge7.unresolved_components) ||
        throw(ErrorException("barrier-8 seed set disagrees with barrier-7 result"))
    old_parent = workspace.uf_parent
    old_seed_nodes = workspace.seed_nodes
    old_known_nodes = workspace.known_nodes
    old_adjacency = workspace.closure_adjacency
    old_n7_states = workspace.n7_states
    old_n7_nodes = workspace.n7_nodes
    old_n7_local_by_node = workspace.n7_local_by_node
    roots = Set(_n6b8_old_find(old_parent, old_seed_nodes[component])
                for component in seed_components)

    known_old_nodes = Set(values(old_known_nodes))
    for component in seed_components
        root = _n6b8_old_find(old_parent, old_seed_nodes[component])
        any(_n6b8_old_find(old_parent, node) == root for node in known_old_nodes) &&
            throw(ErrorException(
                "barrier-7 unresolved seed is already connected to a known node"))
    end

    old_n7_diagram_by_node = Dict{Int,KnotGaussDiagram}(
        old_n7_nodes[index] => old_n7_states[index]
        for index in eachindex(old_n7_states))
    for index in eachindex(old_n7_states)
        diagram = canonicalize(old_n7_states[index])
        old_root = _n6b8_old_find(old_parent, old_n7_nodes[index])
        labels = sort!(unique(String[
            workspace.component_type_labels[component]
            for (component, node) in old_known_nodes
            if _n6b8_old_find(old_parent, node) == old_root]))
        n7_boundary_labels[diagram] = labels
        old_n7_diagram_by_node[old_n7_nodes[index]] = diagram
        old_root in roots && isempty(labels) || old_root in roots || continue
    end

    # Map old seed vertices and every N7 vertex in the selected unresolved
    # closure into the new transient graph.  The ordering is canonical to
    # keep the eventual fingerprint independent of dictionary iteration.
    old_node_to_local = Dict{Int,Int}(
        old_seed_nodes[component] => seed_nodes[component]
        for component in seed_components)
    selected_old_n7 = sort!([old_n7_nodes[index]
                             for index in eachindex(old_n7_states)
                             if _n6b8_old_find(old_parent, old_n7_nodes[index]) in roots];
                            by=node -> canonical_serialization(
                                old_n7_diagram_by_node[node]))
    for old_node in selected_old_n7
        diagram = old_n7_diagram_by_node[old_node]
        isempty(n7_boundary_labels[diagram]) || throw(ErrorException(
            "selected unresolved barrier-7 N7 state has a known boundary"))
        local_node, is_new = _n6b8_ensure_n7!(
            diagram, n7_index, n7_states, n7_nodes, n7_local_by_node,
            uf, adjacency)
        is_new || error("duplicate imported barrier-7 N7 identity")
        old_node_to_local[old_node] = local_node
    end

    for old_first in eachindex(old_adjacency)
        for old_second in old_adjacency[old_first]
            old_first < old_second || continue
            selected_first = haskey(old_node_to_local, old_first)
            selected_second = haskey(old_node_to_local, old_second)
            selected_first || selected_second || continue
            (old_first in known_old_nodes || old_second in known_old_nodes) &&
                throw(ErrorException(
                    "selected barrier-7 unresolved closure touches a known node"))
            selected_first && selected_second || throw(ErrorException(
                "barrier-7 selected closure has an unmapped boundary vertex"))
            _n6b8_add_graph_edge!(stats, uf, adjacency, Dict{Int,Int}(),
                                  classification, old_node_to_local[old_first],
                                  old_node_to_local[old_second]; imported=true)
        end
    end

    # The graph-edge import above must not run known-label conflict checks with
    # an empty map as a substitute for the real map.  No known node can lie in
    # an unresolved root, so this is safe and preserves the compact import.
    for local_node in n7_nodes
        push!(queue_kind, UInt8(1))
        push!(queue_nodes, local_node)
    end
    stats.max_frontier = max(stats.max_frontier, length(queue_nodes))
    length(selected_old_n7), old_node_to_local, old_n7_diagram_by_node
end

function _n6b8_component_type_counts(classification::KnotN6Classification,
                                     seed_components::Vector{Int},
                                     resolved::Dict{Int,String})
    component_counts = Dict{String,Int}()
    state_counts = Dict{String,Int}()
    for (component, label) in resolved
        component_counts[label] = get(component_counts, label, 0) + 1
        state_counts[label] = get(state_counts, label, 0) +
            count(==(component), classification.components)
    end
    component_counts, state_counts
end

"""Run the unresolved barrier-7 components through barrier 8.

The three-argument method is useful for focused fixtures.  The production
script obtains `classification`, `bridge7`, and `workspace` in one run and
then calls this method.  No barrier above eight is accepted.
"""
function run_knot_n6_targeted_barrier8(
        classification::KnotN6Classification,
        bridge7::KnotN6TargetedBridgeResult,
        workspace;
        barrier::Integer=KNOT_N6_TARGETED_BARRIER8,
        progress::Bool=true)
    Int(barrier) == KNOT_N6_TARGETED_BARRIER8 ||
        throw(ArgumentError("targeted barrier8 closure is fixed to barrier=8"))
    0 <= classification.max_n <= 6 || throw(ArgumentError(
        "classification max_n must lie in 0:6"))
    bridge7.barrier == 7 || throw(ArgumentError(
        "barrier8 input must be a barrier-7 bridge result"))
    started = UInt64(time_ns())
    seed_components = sort!(copy(bridge7.unresolved_components))
    seed_state_indices = _n6b8_seed_state_indices(classification, seed_components)
    uf = _N6B8UnionFind()
    adjacency = Vector{Vector{Int}}()
    seed_nodes = Dict{Int,Int}()
    for component in seed_components
        seed_nodes[component] = _n6b8_new_node!(uf, adjacency)
    end
    known_nodes = Dict{Int,Int}()
    known_component_by_label = Dict{String,Int}(
        label => component for (component, label) in
        classification.component_type_labels)
    length(known_component_by_label) ==
        length(classification.component_type_labels) || throw(ErrorException(
            "classification has duplicate known component labels"))
    stats = _N6B8Stats()
    n7_index = Dict{KnotGaussDiagram,Int}()
    n7_states = KnotGaussDiagram[]
    n7_nodes = Int[]
    n7_local_by_node = Dict{Int,Int}()
    n7_imported_local_nodes = Set{Int}()
    n7_boundary_labels = Dict{KnotGaussDiagram,Vector{String}}()
    n8_index = Dict{KnotGaussDiagram,Int}()
    n8_states = KnotGaussDiagram[]
    n8_nodes = Int[]
    n8_local_by_node = Dict{Int,Int}()
    queue_kind = UInt8[]
    queue_nodes = Int[]
    queue_head = 1

    imported_n7_count, old_node_to_local, old_n7_diagram_by_node =
        _n6b8_import_bridge7_workspace!(
            classification, bridge7, workspace, seed_components, uf, adjacency,
            seed_nodes, n7_index, n7_states, n7_nodes, n7_local_by_node,
            n7_boundary_labels, stats, queue_kind, queue_nodes)
    union!(n7_imported_local_nodes, n7_nodes)

    function enqueue!(kind::UInt8, node::Int)
        push!(queue_kind, kind)
        push!(queue_nodes, node)
        stats.max_frontier = max(stats.max_frontier,
                                 length(queue_nodes) - queue_head + 1)
    end

    # Queue the imported N7 states exactly once.  The import helper populated
    # the arrays and queue directly so that it can preserve canonical ordering.
    # New states discovered below use this local enqueue helper.

    function known_boundary_labels(diagram::KnotGaussDiagram)
        get(n7_boundary_labels, canonicalize(diagram), String[])
    end

    function add_n7_target!(source_node::Int, source_kind::Symbol,
                            source_index::Int, raw_target::KnotGaussDiagram,
                            move_type::Symbol, transition_key_kind::Symbol;
                            require_bridge7::Bool=false,
                            enqueue_new::Bool=true)
        canonical_target = canonicalize(raw_target)
        labels = known_boundary_labels(canonical_target)
        if require_bridge7 && !haskey(n7_boundary_labels, canonical_target)
            throw(ErrorException(
                "barrier-7 N7 target was not reproduced in the barrier-8 import"))
        end
        _n6b8_check_new_state_limit!(stats, started,
                                     length(n7_states), length(n8_states))
        target_node, is_new = _n6b8_ensure_n7!(
            canonical_target, n7_index, n7_states, n7_nodes,
            n7_local_by_node, uf, adjacency)
        key = (transition_key_kind, source_index, target_node, move_type)
        _n6b8_add_transition!(stats, uf, adjacency, known_nodes, classification,
                              key, source_node, target_node)
        isempty(labels) || _n6b8_connect_known_labels!(
            target_node, labels, known_component_by_label, known_nodes, uf,
            adjacency, stats, classification)
        if is_new
            _n6b8_classical_new!(canonical_target, stats)
            enqueue_new && isempty(labels) && enqueue!(UInt8(1), target_node)
        end
        target_node, is_new
    end

    function add_n8_target!(source_node::Int, source_index::Int,
                            raw_target::KnotGaussDiagram, move_type::Symbol,
                            transition_key_kind::Symbol)
        canonical_target = canonicalize(raw_target)
        _n6b8_check_new_state_limit!(stats, started,
                                     length(n7_states), length(n8_states))
        target_node, is_new = _n6b8_ensure_n8!(
            canonical_target, n8_index, n8_states, n8_nodes,
            n8_local_by_node, uf, adjacency)
        key = (transition_key_kind, source_index, target_node, move_type)
        _n6b8_add_transition!(stats, uf, adjacency, known_nodes, classification,
                              key, source_node, target_node)
        if is_new
            _n6b8_classical_new!(canonical_target, stats)
            enqueue!(UInt8(2), target_node)
        end
        target_node, is_new
    end

    function connect_n6_target!(source_index::Int, source::KnotGaussDiagram,
                               raw_target::KnotGaussDiagram, move_type::Symbol)
        canonical_target = canonicalize(raw_target)
        target_index = get(classification.state_index, canonical_target, 0)
        if target_index == 0
            stats.n6_membership_failures += 1
            throw(ErrorException(
                "N<=6 target from barrier-8 transition is absent from input census"))
        end
        target_component = classification.components[target_index]
        target_node = if haskey(seed_nodes, target_component)
            seed_nodes[target_component]
        elseif haskey(classification.component_type_labels, target_component)
            _n6b8_new_known_node!(target_component, known_nodes, uf, adjacency)
        else
            stats.n6_membership_failures += 1
            throw(ErrorException("N<=6 target belongs to an unclassified component"))
        end
        key = (:to_n6, source_index, target_index, move_type)
        _n6b8_add_transition!(stats, uf, adjacency, known_nodes, classification,
                              key, source_index == 0 ? source_index :
                              get(n7_local_by_node, source_index, source_index),
                              target_node)
        target_index, target_component
    end

    # Seed all N<=6 states.  The source node in the transient closure graph is
    # the unresolved component node, while transition keys retain the actual
    # canonical N6 state index for deduplication and accounting.
    try
        for source_index in seed_state_indices
            source = classification.states[source_index]
            source_component = classification.components[source_index]
            source_node = seed_nodes[source_component]
            for move_name in ("R1", "R2")
                move_type = Symbol(move_name)
                candidates = movable(source, move_name; param=:insert)
                for candidate in candidates
                    stats.candidate_count += 1
                    delta = move_type == :R1 ? 1 : 2
                    target_crossings = crossing_count(source) + delta
                    if target_crossings > KNOT_N6_TARGETED_BARRIER8
                        stats.excluded_over_barrier_count += 1
                        continue
                    raw_target = _n6b8_checked_apply(
                        source, candidate, stats,
                        "seed state=$source_index move=$move_name")
                    canonical_target = canonicalize(raw_target)
                    if target_crossings <= 6
                        target_index = get(classification.state_index,
                                            canonical_target, 0)
                        target_index == 0 && throw(ErrorException(
                            "reusable N<=6 insertion target is absent"))
                        classification.components[target_index] == source_component ||
                            throw(ErrorException(
                                "reusable N<=6 insertion crosses components"))
                        _n6b8_existing_edge(classification, source_index,
                                             target_index, move_type) ||
                            throw(ErrorException(
                                "existing N<=6 edge missing for reusable insertion"))
                        stats.reused_n6_candidate_count += 1
                    elseif target_crossings == 7
                        haskey(n7_boundary_labels, canonical_target) ||
                            throw(ErrorException(
                                "direct N6 insertion did not reproduce barrier-7 N7 state"))
                        add_n7_target!(source_node, :n6, source_index, raw_target,
                                       move_type, :n6_to_n7;
                                       require_bridge7=true, enqueue_new=false)
                    else
                        target_crossings == 8 || throw(ErrorException(
                            "N6 insertion target exceeds barrier8"))
                        add_n8_target!(source_node, source_index, raw_target,
                                       move_type, :n6_to_n8)
                    end
                end
            end
        end

        while queue_head <= length(queue_nodes)
            kind = queue_kind[queue_head]
            current_node = queue_nodes[queue_head]
            queue_head += 1
            if kind == UInt8(1)
                local_index = get(n7_local_by_node, current_node, 0)
                local_index == 0 && error("N7 queue node-index inconsistency")
                source = n7_states[local_index]
                stats.n7_processed_count += 1

                # N7 insertion to N8, then deletion/R3 at the N7 boundary.
                for candidate in movable(source, "R1"; param=:insert)
                    stats.candidate_count += 1
                    raw_target = _n6b8_checked_apply(
                        source, candidate, stats,
                        "N7 node=$current_node move=R1 insert")
                    add_n8_target!(current_node, local_index, raw_target, :R1,
                                   :n7_to_n8)
                end
                stats.candidate_count += length(movable(source, "R2";
                                                        param=:insert))
                stats.excluded_over_barrier_count += length(movable(
                    source, "R2"; param=:insert))
                for (move_name, direction) in (("R1", :delete),
                                               ("R2", :delete),
                                               ("R3", nothing))
                    candidates = direction === nothing ? movable(source, move_name) :
                        movable(source, move_name; param=direction)
                    for candidate in candidates
                        stats.candidate_count += 1
                        raw_target = _n6b8_checked_apply(
                            source, candidate, stats,
                            "N7 node=$current_node move=$move_name")
                        target_crossings = crossing_count(raw_target)
                        target_crossings <= 7 || throw(ErrorException(
                            "N7 deletion/R3 target exceeds seven crossings"))
                        canonical_target = canonicalize(raw_target)
                        if target_crossings <= 6
                            target_index = get(classification.state_index,
                                                canonical_target, 0)
                            target_index == 0 && throw(ErrorException(
                                "N<=6 target from N7 is absent from input census"))
                            target_component = classification.components[target_index]
                            target_node = if haskey(seed_nodes, target_component)
                                seed_nodes[target_component]
                            elseif haskey(classification.component_type_labels,
                                          target_component)
                                _n6b8_new_known_node!(target_component, known_nodes,
                                                      uf, adjacency)
                            else
                                throw(ErrorException(
                                    "N<=6 target belongs to an unknown component"))
                            end
                            _n6b8_add_transition!(
                                stats, uf, adjacency, known_nodes, classification,
                                (:n7_to_n6, current_node, target_index,
                                 Symbol(move_name)), current_node, target_node)
                        else
                            add_n7_target!(
                                current_node, :n7, current_node, raw_target,
                                Symbol(move_name), :n7_to_n7;
                                require_bridge7=current_node in n7_imported_local_nodes,
                                enqueue_new=current_node in n7_imported_local_nodes)
                        end
                    end
                end
            elseif kind == UInt8(2)
                local_index = get(n8_local_by_node, current_node, 0)
                local_index == 0 && error("N8 queue node-index inconsistency")
                source = n8_states[local_index]
                stats.n8_processed_count += 1
                for (move_name, direction) in (("R1", :delete),
                                               ("R2", :delete),
                                               ("R3", nothing))
                    candidates = direction === nothing ? movable(source, move_name) :
                        movable(source, move_name; param=direction)
                    for candidate in candidates
                        stats.candidate_count += 1
                        raw_target = _n6b8_checked_apply(
                            source, candidate, stats,
                            "N8 node=$current_node move=$move_name")
                        target_crossings = crossing_count(raw_target)
                        target_crossings <= 8 || throw(ErrorException(
                            "N8 deletion/R3 target exceeds eight crossings"))
                        canonical_target = canonicalize(raw_target)
                        if target_crossings <= 6
                            target_index = get(classification.state_index,
                                                canonical_target, 0)
                            target_index == 0 && throw(ErrorException(
                                "N<=6 target from N8 is absent from input census"))
                            target_component = classification.components[target_index]
                            target_node = if haskey(seed_nodes, target_component)
                                seed_nodes[target_component]
                            elseif haskey(classification.component_type_labels,
                                          target_component)
                                _n6b8_new_known_node!(target_component, known_nodes,
                                                      uf, adjacency)
                            else
                                throw(ErrorException(
                                    "N<=6 target belongs to an unknown component"))
                            end
                            _n6b8_add_transition!(
                                stats, uf, adjacency, known_nodes, classification,
                                (:n8_to_n6, current_node, target_index,
                                 Symbol(move_name)), current_node, target_node)
                        elseif target_crossings == 7
                            add_n7_target!(
                                current_node, :n8, current_node, raw_target,
                                Symbol(move_name), :n8_to_n7;
                                require_bridge7=false, enqueue_new=true)
                        else
                            target_crossings == 8 || throw(ErrorException(
                                "N8 target exceeds barrier8"))
                            add_n8_target!(
                                current_node, local_index, raw_target,
                                Symbol(move_name), :n8_to_n8)
                        end
                    end
                end
            else
                error("unknown barrier-8 queue kind $kind")
            end
            total_processed = stats.n7_processed_count + stats.n8_processed_count
            total_processed % 1000 == 0 &&
                _n6b8_check_limits!(stats, started, length(n7_states),
                                    length(n8_states))
            progress && total_processed % 10_000 == 0 && println(stderr,
                "[knot-n6-barrier8] processed=$(total_processed) " *
                "n7=$(length(n7_states)) n8=$(length(n8_states)) " *
                "candidates=$(stats.candidate_count) frontier=" *
                "$(length(queue_nodes) - queue_head + 1)")
        end
    end
    catch error
        if error isa _N6B8ResourceLimit
            stats.resource_limited = true
            stats.resource_reason = error.reason
        else
            rethrow()
        end
    end

    connections = _n6b8_component_connections(
        seed_components, seed_nodes, known_nodes, uf, classification)
    resolved = Dict{Int,String}()
    unresolved = Int[]
    component_status = Dict{Int,String}()
    for component in seed_components
        labels = connections[component]
        if isempty(labels)
            push!(unresolved, component)
            component_status[component] = stats.resource_limited ?
                "resource_limited" : "barrier8_unresolved"
        else
            length(labels) == 1 || throw(ErrorException(
                "barrier-8 closure has multiple known labels for component $component"))
            resolved[component] = only(labels)
            component_status[component] = only(labels)
        end
    end
    resolved_component_type_counts, resolved_state_type_counts =
        _n6b8_component_type_counts(classification, seed_components, resolved)
    distances = Dict{Int,Union{Nothing,Int}}(
        component => _n6b8_shortest_known_distance(
            seed_nodes[component], known_nodes, adjacency)
        for component in seed_components)
    groups = _n6b8_closure_groups(seed_components, seed_nodes, uf)
    status = stats.resource_limited ? "resource_limited" :
        isempty(unresolved) ? "barrier8_resolved_all_seed_components" :
        "barrier8_unresolved"
    fingerprint = _n6b8_fingerprint(
        classification, bridge7, n7_states, n8_states, stats.graph_keys,
        component_status)
    elapsed = (UInt64(time_ns()) - started) / 1e9
    KnotN6TargetedBarrier8Result(
        KNOT_N6_TARGETED_BARRIER8, classification.fingerprint,
        bridge7.fingerprint, seed_components, length(seed_state_indices),
        imported_n7_count, length(n7_states), length(n8_states),
        stats.candidate_count, stats.adopted_transition_count,
        stats.duplicate_transition_count, stats.imported_graph_edge_count,
        stats.reused_n6_candidate_count, stats.excluded_over_barrier_count,
        stats.n7_processed_count, stats.n8_processed_count, stats.max_frontier,
        stats.transition_verification_failures, stats.n6_membership_failures,
        stats.classicality_checks, stats.classicality_failures, connections,
        distances, resolved, component_status,
        resolved_component_type_counts, resolved_state_type_counts, groups,
        unresolved, status, stats.resource_reason, fingerprint, elapsed)
end

function build_knot_n6_targeted_barrier8(; progress::Bool=true)
    classification = build_knot_n6_classification(6; barrier=6,
                                                   progress=progress)
    workspace_ref = Ref{Any}(nothing)
    bridge7 = run_knot_n6_targeted_bridge(
        classification; barrier=7, progress=progress,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; barrier=8,
        progress=progress)
end

targeted_knot_n6_barrier8(args...; kwargs...) =
    run_knot_n6_targeted_barrier8(args...; kwargs...)

"Return compact barrier-8 statistics without states, paths, or witnesses."
function knot_n6_targeted_barrier8_summary(
        result::KnotN6TargetedBarrier8Result)
    Dict{String,Any}(
        "version" => KNOT_N6_TARGETED_BARRIER8_VERSION,
        "barrier" => result.barrier,
        "n6_input_fingerprint" => result.n6_input_fingerprint,
        "barrier7_fingerprint" => result.barrier7_fingerprint,
        "seed_components" => copy(result.seed_components),
        "seed_state_count" => result.seed_state_count,
        "imported_n7_state_count" => result.imported_n7_state_count,
        "n7_state_count" => result.n7_state_count,
        "n8_state_count" => result.n8_state_count,
        "candidate_count" => result.candidate_count,
        "adopted_transition_count" => result.adopted_transition_count,
        "duplicate_transition_count" => result.duplicate_transition_count,
        "imported_graph_edge_count" => result.imported_graph_edge_count,
        "reused_n6_candidate_count" => result.reused_n6_candidate_count,
        "excluded_over_barrier_count" => result.excluded_over_barrier_count,
        "n7_processed_count" => result.n7_processed_count,
        "n8_processed_count" => result.n8_processed_count,
        "max_frontier" => result.max_frontier,
        "transition_verification_failures" =>
            result.transition_verification_failures,
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
        "component_status" => Dict(string(component) => status
                                    for (component, status) in
                                    sort!(collect(result.component_status);
                                          by=first)),
        "resolved_component_type_counts" =>
            copy(result.resolved_component_type_counts),
        "resolved_state_type_counts" => copy(result.resolved_state_type_counts),
        "closure_groups" => [copy(group) for group in result.closure_groups],
        "unresolved_components" => copy(result.unresolved_components),
        "status" => result.status,
        "resource_reason" => result.resource_reason,
        "fingerprint_sha256" => result.fingerprint,
        "elapsed_seconds" => result.elapsed_seconds)
end
