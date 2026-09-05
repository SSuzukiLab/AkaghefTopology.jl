"""Finite barrier-6 classification of classical signed one-component Gauss data.

This lane is deliberately independent of the N5 classifier and of the
persistent move database.  The direct classical generator supplies every
rooted classical input through six crossings; this file canonicalizes those
inputs, enumerates only deletion/R3 witnesses, and registers their reverse
insertions through the resulting undirected graph.  An insertion with both
endpoints at most six crossings is the inverse of a deletion on its higher
endpoint, so no separate insertion sweep (and no N5-scale insertion audit) is
needed.

All candidates that become graph edges are applied and checked with
`verify_reidemeister_transition`.  Paths are not retained: for a component
with a known KnotInfo representative, one BFS spanning tree is built and each
tree edge is replayed exactly once.  Alexander/Jones values are calculated at
one representative per component and on a deterministic edge sample only.
Unrepresented components remain `unresolved`; they are not assigned a new
knot type or treated as evidence that a higher barrier is unnecessary.
"""

if !isdefined(@__MODULE__, :KnotGaussDiagram)
    include(joinpath(@__DIR__, "knot_gauss.jl"))
end
if !isdefined(@__MODULE__, :ReidemeisterCandidate)
    include(joinpath(@__DIR__, "reidemeister.jl"))
end
if !isdefined(@__MODULE__, :ClassicalGaussCertificate)
    include(joinpath(@__DIR__, "classical_gauss_realizability.jl"))
end
if !isdefined(@__MODULE__, :ExactKnotLaurentPolynomial)
    include(joinpath(@__DIR__, "knot_classical_invariants.jl"))
end
if !isdefined(@__MODULE__, :direct_classical_gauss_enumeration)
    include(joinpath(@__DIR__, "classical_gauss_enumeration.jl"))
end

using SHA
using JSON

const KNOT_N6_CLASSIFICATION_VERSION = 1
const KNOT_N6_KNOWN_NAMES = ("0_1", "3_1", "4_1", "5_1", "5_2",
                            "6_1", "6_2", "6_3")
const KNOT_N6_MAX_N = 6

"One retained, independently verified edge in the bounded graph." 
struct KnotN6Edge
    source::Int
    target::Int
    move_type::Symbol
    candidate::ReidemeisterCandidate
end

"One local KnotInfo PD representative and its traversal mirror." 
struct KnotN6KnownRepresentative
    name::String
    label::String
    mirror::Bool
    diagram::KnotGaussDiagram
    state_index::Int
    component::Int
    alexander::ExactKnotLaurentPolynomial
    jones::ExactKnotLaurentPolynomial
end

"Invariant values retained for one transient representative per component." 
struct KnotN6ComponentInvariant
    component::Int
    representative::Int
    known::Bool
    alexander::ExactKnotLaurentPolynomial
    jones::ExactKnotLaurentPolynomial
end

"Result of one bounded N<=6 classification run." 
struct KnotN6Classification
    max_n::Int
    barrier::Int
    state_counts::Vector{Int}
    raw_rooted_classical_counts_by_n::Vector{Int}
    raw_rooted_total::Int
    raw_rooted_classical_total::Int
    states::Vector{KnotGaussDiagram}
    state_index::Dict{KnotGaussDiagram,Int}
    edges::Vector{KnotN6Edge}
    adjacency::Vector{Vector{Int}}
    components::Vector{Int}
    component_roots::Vector{Int}
    component_type_labels::Dict{Int,String}
    canonical_knot_type_counts::Dict{String,Int}
    raw_rooted_classical_type_counts::Dict{String,Int}
    raw_rooted_classical_unresolved_count::Int
    component_invariants::Vector{KnotN6ComponentInvariant}
    known_representatives::Vector{KnotN6KnownRepresentative}
    unresolved_components::Vector{Int}
    direct_transition_count::Int
    duplicate_edge_transition_count::Int
    nonclassical_target_count::Int
    spanning_tree_edge_replay_count::Int
    spanning_tree_replay_failures::Int
    edge_invariant_sample_count::Int
    edge_invariant_failures::Vector{String}
    fingerprint::String
    elapsed_seconds::Float64
end

Base.length(result::KnotN6Classification) = length(result.states)

function _n6_check_bound(max_n::Integer, barrier::Integer)
    m = Int(max_n)
    b = Int(barrier)
    0 <= m <= KNOT_N6_MAX_N || throw(ArgumentError("max_n must lie in 0:6"))
    m <= b <= KNOT_N6_MAX_N ||
        throw(ArgumentError("require max_n <= barrier <= 6"))
    m, b
end

function _n6_state_key(d::KnotGaussDiagram)
    canonical_serialization(d)
end

function _n6_state_key_with_n(d::KnotGaussDiagram)
    "n=$(crossing_count(d))|$(_n6_state_key(d))"
end

"Enumerate direct raw classical keys once and fold them to canonical states." 
function _n6_collect_states(max_n::Int; progress=nothing)
    states = KnotGaussDiagram[]
    state_index = Dict{KnotGaussDiagram,Int}()
    state_counts = zeros(Int, max_n + 1)
    raw_classical_by_n = zeros(Int, max_n + 1)
    raw_classical_multiplicity = Dict{KnotGaussDiagram,Int}()
    raw_total = 0
    raw_classical_total = 0

    for n in 0:max_n
        result = direct_classical_gauss_enumeration(n;
            progress=(words, spheres, generated) -> begin
                progress === nothing || progress(n, words, spheres, generated)
                nothing
            end)
        raw_total += rooted_normalized_signed_gc_count(n) * 2^n
        keys = sort!(collect(result.raw_keys))
        raw_classical_by_n[n + 1] = length(keys)
        raw_classical_total += length(keys)
        for key in keys
            diagram = classical_gauss_key_to_diagram(key, n)
            canonical = canonicalize(diagram)
            raw_classical_multiplicity[canonical] =
                get(raw_classical_multiplicity, canonical, 0) + 1
            if !haskey(state_index, canonical)
                push!(states, canonical)
                state_index[canonical] = length(states)
                state_counts[n + 1] += 1
            end
        end
    end
    states, state_index, state_counts, raw_classical_by_n,
    raw_classical_multiplicity, raw_total, raw_classical_total
end

function _n6_canonical_target(target::KnotGaussDiagram,
                              state_index::Dict{KnotGaussDiagram,Int})
    canonical = canonicalize(target)
    get(state_index, canonical, 0), canonical
end

function _n6_register_edge!(edges::Vector{KnotN6Edge},
                            adjacency::Vector{Vector{Int}},
                            edge_index::Dict{Tuple{Int,Int},Int},
                            source::Int, target::Int, move_type::Symbol,
                            candidate::ReidemeisterCandidate)
    pair = source <= target ? (source, target) : (target, source)
    haskey(edge_index, pair) && return false
    push!(edges, KnotN6Edge(source, target, move_type, candidate))
    edge_id = length(edges)
    edge_index[pair] = edge_id
    push!(adjacency[source], edge_id)
    source == target || push!(adjacency[target], edge_id)
    true
end

"""Build all deletion/R3 edges and verify every direct transition.

R1/R2 deletion is enumerated only at the higher endpoint.  Registering that
edge as undirected includes its insertion inverse whenever the lower and
higher endpoints both lie within the barrier.  R3 is self-inverse.
"""
function _n6_build_graph(states::Vector{KnotGaussDiagram},
                          state_index::Dict{KnotGaussDiagram,Int};
                          progress=nothing)
    adjacency = [Int[] for _ in eachindex(states)]
    edges = KnotN6Edge[]
    edge_index = Dict{Tuple{Int,Int},Int}()
    direct_count = 0
    duplicate_count = 0
    nonclassical_count = 0
    for source_index in eachindex(states)
        source = states[source_index]
        for (move_type, direction) in (("R1", :delete), ("R2", :delete),
                                       ("R3", nothing))
            candidates = direction === nothing ? movable(source, move_type) :
                movable(source, move_type; param=direction)
            for candidate in candidates
                direct_count += 1
                raw_target = try
                    apply_reidemeister(source, candidate)
                catch error
                    throw(ErrorException("N6 $move_type application failed at state " *
                                         "$source_index: " * sprint(showerror, error)))
                end
                verify_reidemeister_transition(source, raw_target, candidate) ||
                    throw(ErrorException("N6 $move_type verifier rejected state " *
                                         "$source_index"))
                target_index, canonical_target = _n6_canonical_target(raw_target,
                                                                       state_index)
                if target_index == 0
                    # A Reidemeister move of a classical diagram must remain
                    # in the direct classical census.  A missing target is
                    # therefore a generator or move-consistency defect,
                    # irrespective of a separate Carter predicate; never
                    # silently discard it as a nonclassical target.
                    throw(ErrorException(
                        "N6 classical target is absent from direct census " *
                        "(source=$source_index, move=$move_type)"))
                end
                _n6_register_edge!(edges, adjacency, edge_index, source_index,
                                   target_index, Symbol(move_type), candidate) ||
                    (duplicate_count += 1)
                states[target_index] == canonical_target || throw(ErrorException(
                    "N6 state-index canonical target mismatch"))
            end
        end
        progress === nothing || progress(source_index, length(states), direct_count,
                                         length(edges))
    end
    (edges=edges, adjacency=adjacency, direct_count=direct_count,
     duplicate_count=duplicate_count, nonclassical_count=nonclassical_count)
end

function _n6_components(adjacency::Vector{Vector{Int}},
                        edges::Vector{KnotN6Edge})
    component = zeros(Int, length(adjacency))
    roots = Int[]
    queue = Int[]
    for start in eachindex(adjacency)
        component[start] == 0 || continue
        component_id = length(roots) + 1
        push!(roots, start)
        empty!(queue)
        push!(queue, start)
        component[start] = component_id
        head = 1
        while head <= length(queue)
            current = queue[head]
            head += 1
            for edge_id in adjacency[current]
                edge = edges[edge_id]
                neighbor = edge.source == current ? edge.target : edge.source
                component[neighbor] == 0 || continue
                component[neighbor] = component_id
                push!(queue, neighbor)
            end
        end
    end
    component, roots
end

function _n6_pd_rows(value)
    isempty(strip(String(value))) && return Any[]
    value isa AbstractString ? JSON.parse(String(value)) : value
end

function _n6_load_knotinfo_rows(path::AbstractString)
    payload = JSON.parsefile(path)
    sheets = get(payload, "sheets", Any[])
    isempty(sheets) && throw(ArgumentError("KnotInfo snapshot has no sheet"))
    sheet = first(sheets)
    headers = String.(sheet["headers"])
    indices = Dict(header => i for (i, header) in enumerate(headers))
    haskey(indices, "name") || throw(ArgumentError("KnotInfo rows have no name field"))
    haskey(indices, "pd_notation") ||
        throw(ArgumentError("KnotInfo rows have no pd_notation field"))
    rows = Dict{String,Any}()
    for row in sheet["rows"]
        values = row["values"]
        rows[String(values[indices["name"]])] = row
    end
    rows, indices
end

function _n6_default_snapshot()
    normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260830", "data",
                      "knotinfo-n6.json"))
end

function _n6_known_representatives(state_index::Dict{KnotGaussDiagram,Int},
                                   components::Vector{Int};
                                   max_n::Int,
                                   snapshot_path::AbstractString)
    rows, indices = _n6_load_knotinfo_rows(snapshot_path)
    result = KnotN6KnownRepresentative[]
    for name in KNOT_N6_KNOWN_NAMES
        haskey(rows, name) || throw(ArgumentError("KnotInfo snapshot lacks $name"))
        row = rows[name]
        diagram = pd_to_gauss(_n6_pd_rows(row["values"][indices["pd_notation"]]))
        crossing_count(diagram) <= max_n || continue
        for mirror_flag in (false, true)
            selected = mirror_flag ? mirror(diagram) : diagram
            index = get(state_index, canonicalize(selected), 0)
            index == 0 && throw(ErrorException(
                "known representative $name is absent from N<=$max_n classical states"))
            label = mirror_flag ? "$(name)^-" : "$(name)^+"
            push!(result, KnotN6KnownRepresentative(
                name, label, mirror_flag, selected, index, components[index],
                alexander_polynomial(selected), jones_polynomial(selected)))
        end
    end
    result
end

"Assign labels only after graph components have been computed." 
function _n6_component_type_labels(known::Vector{KnotN6KnownRepresentative})
    by_component = Dict{Int,Dict{String,Vector{KnotN6KnownRepresentative}}}()
    for representative in known
        by_name = get!(by_component, representative.component,
                       Dict{String,Vector{KnotN6KnownRepresentative}}())
        push!(get!(by_name, representative.name, KnotN6KnownRepresentative[]),
              representative)
    end
    labels = Dict{Int,String}()
    for component in sort!(collect(keys(by_component)))
        by_name = by_component[component]
        length(by_name) == 1 || throw(ErrorException(
            "barrier-6 graph merges distinct KnotInfo names in component $component"))
        name = first(keys(by_name))
        representatives = by_name[name]
        components_seen = Set(rep.component for rep in representatives)
        length(components_seen) == 1 || throw(ErrorException(
            "internal component label inconsistency for $name"))
        # A base label is used when the local representative and its mirror
        # are connected; otherwise the local-PD side and mirror remain strict
        # +/-.  This discovers amphicheiral behavior from connectivity instead
        # of assuming a knot table's chiral count.
        labels[component] = if length(representatives) == 2
            name
        else
            only(representatives).label
        end
    end
    labels
end

function _n6_fixed_counts(labels::Dict{Int,String})
    Dict{String,Int}(label => 0 for label in sort!(collect(values(labels))))
end

function _n6_type_counts(states::Vector{KnotGaussDiagram},
                         state_index::Dict{KnotGaussDiagram,Int},
                         components::Vector{Int},
                         component_labels::Dict{Int,String},
                         raw_classical_multiplicity::Dict{KnotGaussDiagram,Int})
    canonical = _n6_fixed_counts(component_labels)
    for component in components
        haskey(component_labels, component) || continue
        canonical[component_labels[component]] += 1
    end
    raw = _n6_fixed_counts(component_labels)
    unresolved = 0
    for (diagram, multiplicity) in raw_classical_multiplicity
        index = state_index[diagram]
        component = components[index]
        if haskey(component_labels, component)
            raw[component_labels[component]] += multiplicity
        else
            unresolved += multiplicity
        end
    end
    canonical, raw, unresolved
end

function _n6_component_invariants(states::Vector{KnotGaussDiagram},
                                  roots::Vector{Int},
                                  known_components::Set{Int})
    result = KnotN6ComponentInvariant[]
    for (component, representative) in enumerate(roots)
        diagram = states[representative]
        push!(result, KnotN6ComponentInvariant(
            component, representative, component in known_components,
            alexander_polynomial(diagram), jones_polynomial(diagram)))
    end
    result
end

function _n6_invariant_pair(diagram::KnotGaussDiagram,
                            cache::Dict{KnotGaussDiagram,
                                        Tuple{ExactKnotLaurentPolynomial,
                                              ExactKnotLaurentPolynomial}})
    get!(cache, diagram) do
        (alexander_polynomial(diagram), jones_polynomial(diagram))
    end
end

"Audit only a deterministic bounded edge sample, not all state invariants." 
function _n6_audit_edge_sample(states::Vector{KnotGaussDiagram},
                               edges::Vector{KnotN6Edge}; sample_limit::Int=128)
    isempty(edges) && return 0, String[]
    stride = max(1, cld(length(edges), sample_limit))
    sample_ids = collect(1:stride:length(edges))
    cache = Dict{KnotGaussDiagram,
                 Tuple{ExactKnotLaurentPolynomial,
                       ExactKnotLaurentPolynomial}}()
    failures = String[]
    for edge_id in sample_ids
        edge = edges[edge_id]
        before = _n6_invariant_pair(states[edge.source], cache)
        after = _n6_invariant_pair(states[edge.target], cache)
        before[1] == after[1] || push!(failures,
            "edge $edge_id $(edge.move_type): Alexander changed")
        before[2] == after[2] || push!(failures,
            "edge $edge_id $(edge.move_type): Jones changed")
    end
    length(sample_ids), failures
end

function _n6_neighbor(edge::KnotN6Edge, current::Int)
    edge.source == current ? edge.target : edge.source
end

"Find and replay one edge in a requested orientation." 
function _n6_replay_edge(result::KnotN6Classification, edge_id::Int,
                         from::Int, to::Int)
    edge = result.edges[edge_id]
    source = result.states[from]
    if edge.source == from && edge.target == to
        raw = apply_reidemeister(source, edge.candidate)
        verify_reidemeister_transition(source, raw, edge.candidate) ||
            throw(ErrorException("stored N6 edge failed replay verification"))
        canonicalize(raw) == result.states[to] ||
            throw(ErrorException("stored N6 edge replay reached wrong state"))
        return
    end

    move_name = String(edge.move_type)
    candidates = edge.move_type == :R3 ? movable(source, move_name) :
        movable(source, move_name; param=:insert)
    for candidate in candidates
        raw = try
            apply_reidemeister(source, candidate)
        catch
            continue
        end
        verify_reidemeister_transition(source, raw, candidate) || continue
        canonicalize(raw) == result.states[to] || continue
        return
    end
    throw(ErrorException("could not replay reverse N6 edge $edge_id"))
end

"Replay each edge of one known component's spanning tree exactly once." 
function _n6_replay_known_spanning_tree(result::KnotN6Classification,
                                        root::Int, component_id::Int)
    parent = zeros(Int, length(result.states))
    parent_edge = zeros(Int, length(result.states))
    parent[root] = root
    queue = Int[root]
    head = 1
    while head <= length(queue)
        current = queue[head]
        head += 1
        result.components[current] == component_id || continue
        for edge_id in result.adjacency[current]
            neighbor = _n6_neighbor(result.edges[edge_id], current)
            result.components[neighbor] == component_id || continue
            parent[neighbor] == 0 || continue
            parent[neighbor] = current
            parent_edge[neighbor] = edge_id
            push!(queue, neighbor)
        end
    end
    replays = 0
    failures = 0
    for state in eachindex(result.states)
        result.components[state] == component_id || continue
        state == root && continue
        if parent[state] == 0
            failures += 1
            continue
        end
        try
            _n6_replay_edge(result, parent_edge[state], state, parent[state])
        catch
            failures += 1
        end
        replays += 1
    end
    replays, failures
end

function _n6_replay_known_trees(result::KnotN6Classification,
                                roots::Vector{Int}, known::Vector{KnotN6KnownRepresentative})
    roots_by_component = Dict{Int,Int}()
    for representative in known
        get!(roots_by_component, representative.component,
             representative.state_index)
    end
    total_replays = 0
    total_failures = 0
    for component in sort!(collect(keys(roots_by_component)))
        replays, failures = _n6_replay_known_spanning_tree(
            result, roots_by_component[component], component)
        total_replays += replays
        total_failures += failures
    end
    total_replays, total_failures
end

function _n6_unresolved_components(components::Vector{Int},
                                   component_labels::Dict{Int,String})
    sort!(collect(setdiff(Set(components), Set(keys(component_labels)))))
end

function _n6_fingerprint(states::Vector{KnotGaussDiagram},
                         edges::Vector{KnotN6Edge})
    context = SHA.SHA256_CTX()
    for diagram in states
        SHA.update!(context, codeunits("state|$(_n6_state_key_with_n(diagram))\n"))
    end
    pairs = [(min(edge.source, edge.target), max(edge.source, edge.target),
              edge.move_type) for edge in edges]
    sort!(pairs)
    for (source, target, move_type) in pairs
        SHA.update!(context, codeunits("edge|$(source)|$(target)|$(move_type)\n"))
    end
    bytes2hex(SHA.digest!(context))
end

function _n6_validate_known_invariants(result::KnotN6Classification)
    by_component = Dict{Int,Vector{KnotN6KnownRepresentative}}()
    for representative in result.known_representatives
        push!(get!(by_component, representative.component,
                   KnotN6KnownRepresentative[]), representative)
    end
    for invariant in result.component_invariants
        haskey(by_component, invariant.component) || continue
        representatives = by_component[invariant.component]
        any(representative -> representative.alexander == invariant.alexander &&
                              representative.jones == invariant.jones,
            representatives) || throw(ErrorException(
                "component invariant disagrees with known representative"))
    end
    true
end

"""Build the complete direct-classical barrier-6 graph in memory."""
function build_knot_n6_classification(max_n::Integer=6;
                                      barrier::Integer=max_n,
                                      snapshot_path::AbstractString=_n6_default_snapshot(),
                                      progress=true)
    max_n, barrier = _n6_check_bound(max_n, barrier)
    started = time_ns()
    progress_callback = progress ?
        (n, words, spheres, generated) -> begin
            words % 1000 == 0 && println(stderr,
                "[knot-n6] census n=$(n) unsigned_words=$(words) " *
                "spherical_immersions=$(spheres) generated=$(generated)")
            nothing
        end : nothing
    states, state_index, state_counts, raw_by_n, raw_multiplicity,
    raw_total, raw_classical_total = _n6_collect_states(max_n;
                                                         progress=progress_callback)
    graph_progress = progress ?
        (index, total, direct_count, edge_count) -> begin
            index % 1000 == 0 && println(stderr,
                "[knot-n6] graph states=$(index)/$(total) direct=$(direct_count) " *
                "edges=$(edge_count)")
            nothing
        end : nothing
    graph = _n6_build_graph(states, state_index; progress=graph_progress)
    graph.nonclassical_count == 0 || throw(ErrorException(
        "N6 graph retained nonclassical targets: $(graph.nonclassical_count)"))
    components, roots = _n6_components(graph.adjacency, graph.edges)
    known = _n6_known_representatives(state_index, components;
                                      max_n=max_n, snapshot_path=snapshot_path)
    component_labels = _n6_component_type_labels(known)
    canonical_type_counts, raw_type_counts, raw_unresolved = _n6_type_counts(
        states, state_index, components, component_labels, raw_multiplicity)
    known_components = Set(rep.component for rep in known)
    component_invariants = _n6_component_invariants(states, roots, known_components)
    unresolved = _n6_unresolved_components(components, component_labels)

    # Construct a provisional result for the tree replay, then fill its
    # verification counters without retaining parent arrays or paths.
    provisional = KnotN6Classification(
        max_n, barrier, state_counts, raw_by_n, raw_total, raw_classical_total,
        states, state_index, graph.edges, graph.adjacency, components, roots,
        component_labels, canonical_type_counts, raw_type_counts, raw_unresolved,
        component_invariants, known, unresolved, graph.direct_count,
        graph.duplicate_count, graph.nonclassical_count, 0, 0, 0, String[], "", 0.0)
    tree_replays, tree_failures = _n6_replay_known_trees(provisional, roots, known)
    tree_failures == 0 || throw(ErrorException(
        "N6 spanning-tree replay failed for $tree_failures edges"))
    sample_count, invariant_failures = _n6_audit_edge_sample(states, graph.edges)
    isempty(invariant_failures) || throw(ErrorException(
        "N6 sampled Alexander/Jones audit failed: " * join(invariant_failures, "; ")))
    _n6_validate_known_invariants(provisional)
    fingerprint = _n6_fingerprint(states, graph.edges)
    elapsed = (time_ns() - started) / 1e9
    KnotN6Classification(
        max_n, barrier, state_counts, raw_by_n, raw_total, raw_classical_total,
        states, state_index, graph.edges, graph.adjacency, components, roots,
        component_labels, canonical_type_counts, raw_type_counts, raw_unresolved,
        component_invariants, known, unresolved, graph.direct_count,
        graph.duplicate_count, graph.nonclassical_count, tree_replays,
        tree_failures, sample_count, invariant_failures, fingerprint, elapsed)
end

classify_knot_n6(args...; kwargs...) =
    build_knot_n6_classification(args...; kwargs...)

"Return a compact summary; no Carter data, paths, or state invariants." 
function knot_n6_summary(result::KnotN6Classification)
    known_components = sort!(unique(rep.component for rep in result.known_representatives))
    unresolved_state_count = sum(count(==(component), result.components)
                                 for component in result.unresolved_components;
                                 init=0)
    Dict{String,Any}(
        "version" => KNOT_N6_CLASSIFICATION_VERSION,
        "max_n" => result.max_n,
        "barrier" => result.barrier,
        "state_counts" => copy(result.state_counts),
        "raw_rooted_classical_counts_by_n" =>
            copy(result.raw_rooted_classical_counts_by_n),
        "raw_rooted_total" => result.raw_rooted_total,
        "raw_rooted_classical_total" => result.raw_rooted_classical_total,
        "total_states" => length(result.states),
        "edge_count" => length(result.edges),
        "direct_transition_count" => result.direct_transition_count,
        "duplicate_edge_transition_count" =>
            result.duplicate_edge_transition_count,
        "nonclassical_target_count" => result.nonclassical_target_count,
        "component_count" => length(result.component_roots),
        "component_type_labels" => Dict(string(component) => label
                                         for (component, label) in
                                         sort!(collect(result.component_type_labels))),
        "canonical_knot_type_counts" => copy(result.canonical_knot_type_counts),
        "raw_rooted_classical_type_counts" =>
            copy(result.raw_rooted_classical_type_counts),
        "raw_rooted_classical_unresolved_count" =>
            result.raw_rooted_classical_unresolved_count,
        "known_representative_count" => length(result.known_representatives),
        "known_component_count" => length(known_components),
        "known_reachable_state_count" => length(result.states) - unresolved_state_count,
        "unresolved_state_count" => unresolved_state_count,
        "unresolved_component_count" => length(result.unresolved_components),
        "unresolved_components" => copy(result.unresolved_components),
        "spanning_tree_edge_replay_count" => result.spanning_tree_edge_replay_count,
        "spanning_tree_replay_failures" => result.spanning_tree_replay_failures,
        "edge_invariant_sample_count" => result.edge_invariant_sample_count,
        "edge_invariant_failures" => copy(result.edge_invariant_failures),
        "fingerprint_sha256" => result.fingerprint,
        "elapsed_seconds" => result.elapsed_seconds)
end
