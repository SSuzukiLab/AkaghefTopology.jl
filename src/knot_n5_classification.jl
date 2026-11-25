"""Finite barrier-`N` classical signed-Gauss move classification.

This file is deliberately a standalone research lane.  It includes the four
small KnotGauss dependencies when they are not already present in the caller,
but it does not modify `Topology.jl` or persist a database.  The default
`barrier=5` run enumerates every Carter-genus-zero canonical state through
five crossings and builds the induced R1/R2/R3 graph.

The graph is generated from deletions and R3 moves only.  Every R1/R2
insertion whose two endpoints have crossing number at most the barrier has a
deletion inverse on its higher-crossing endpoint.  Since every canonical
higher endpoint is visited and all local deletion witnesses are enumerated,
registering each deletion edge as undirected is exactly the complete
insertion closure within the barrier.  This is a completeness argument for
the finite graph, not an isotopy claim outside the enumerated barrier.
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

using SHA
using JSON

const KNOT_N5_CLASSIFICATION_VERSION = 1
const KNOT_N5_EXPECTED_STATE_COUNTS = (1, 2, 10, 60, 536, 5752)
const KNOT_N5_KNOWN_NAMES = ("0_1", "3_1", "4_1", "5_1", "5_2")
const KNOT_N5_STRICT_KNOT_TYPES =
    ("0_1", "3_1^+", "3_1^-", "4_1", "5_1^+", "5_1^-", "5_2^+", "5_2^-")

"""
    KnotN5Edge(source, target, move_type, candidate) -> KnotN5Edge

Store a directed Reidemeister witness between indexed states of a bounded graph.

# Arguments
The field meanings are owned by [Result field meanings](../specs/knot-classification.md#result-field-meanings).
Q6 ([conditions](../specs/knot-classification.md#numbered-conditions)) is the
caller's responsibility. The field constructor does not check run consistency,
witness validity, counts or labels. Julia converts arguments to their declared
field types; incompatible conversions can raise `MethodError`, `InexactError`
or other conversion exceptions.

# Returns
A `KnotN5Edge` with fields `source`, `target`, `move_type`, `candidate`.
Their individual meanings are given in the linked result-field entry.

# Guarantees
1. The returned fields equal the supplied values after conversion. Under Q6 they describe the finite graph or subordinate record specified by the linked definition, with exactly the enclosing run's population and bounds. [proved/agent: by construction]
2. Construction leaves the arguments unchanged; vectors, dictionaries and contained records already of the declared types remain shared with the result. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[classification labels](../specs/knot-classification.md#evidence-and-boundaries).

# Not guaranteed
The field constructor is not a validator or a classification computation.
A stored witness, resolved label or zero failure count is not independently
verified here. Distinct bounded components need not be distinct knot types.
The general classicality and knot-polynomial source identifications retain
the evidence boundaries stated in the specification. Mutable fields can
invalidate Q6 after construction.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Definitions and Q6 — finite graph and result semantics; no unrestricted knot classification claim.

# Examples
```jldoctest
julia> d = KnotGaussDiagram(Int[], Int[]);

julia> c = first(movable(d, "R1"; param=:insert));

julia> e = KnotN5Edge(1, 2, :R1, c);

julia> verify_reidemeister_transition(d, apply_reidemeister(d, e.candidate), e.candidate)
true

```

Role: convert
"""
struct KnotN5Edge
    source::Int
    target::Int
    move_type::Symbol
    candidate::ReidemeisterCandidate
end

"""
    KnotN5KnownRepresentative(name, label, mirror, diagram, state_index, component, alexander, jones) -> KnotN5KnownRepresentative

Store a named knot-diagram representative and its location in a finite move graph.

# Arguments
The field meanings are owned by [Result field meanings](../specs/knot-classification.md#result-field-meanings).
Q6 ([conditions](../specs/knot-classification.md#numbered-conditions)) is the
caller's responsibility. The field constructor does not check run consistency,
witness validity, counts or labels. Julia converts arguments to their declared
field types; incompatible conversions can raise `MethodError`, `InexactError`
or other conversion exceptions.

# Returns
A `KnotN5KnownRepresentative` with fields `name`, `label`, `mirror`, `diagram`, `state_index`, `component`, `alexander`, `jones`.
Their individual meanings are given in the linked result-field entry.

# Guarantees
1. The returned fields equal the supplied values after conversion. Under Q6 they describe the finite graph or subordinate record specified by the linked definition, with exactly the enclosing run's population and bounds. [proved/agent: by construction]
2. Construction leaves the arguments unchanged; vectors, dictionaries and contained records already of the declared types remain shared with the result. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[classification labels](../specs/knot-classification.md#evidence-and-boundaries).

# Not guaranteed
The field constructor is not a validator or a classification computation.
A stored witness, resolved label or zero failure count is not independently
verified here. Distinct bounded components need not be distinct knot types.
The general classicality and knot-polynomial source identifications retain
the evidence boundaries stated in the specification. Mutable fields can
invalidate Q6 after construction.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Definitions and Q6 — finite graph and result semantics; no unrestricted knot classification claim.

# Examples
```jldoctest
julia> r = build_knot_n5_classification(0);

julia> first(r.known_representatives).label
"0_1"

```

Role: convert
"""
struct KnotN5KnownRepresentative
    name::String
    label::String
    mirror::Bool
    diagram::KnotGaussDiagram
    state_index::Int
    component::Int
    alexander::ExactKnotLaurentPolynomial
    jones::ExactKnotLaurentPolynomial
end

"""
    KnotN5ComponentInvariant(component, representative, alexander, jones) -> KnotN5ComponentInvariant

Associate two polynomial records with a representative of a finite move component.

# Arguments
The field meanings are owned by [Result field meanings](../specs/knot-classification.md#result-field-meanings).
Q6 ([conditions](../specs/knot-classification.md#numbered-conditions)) is the
caller's responsibility. The field constructor does not check run consistency,
witness validity, counts or labels. Julia converts arguments to their declared
field types; incompatible conversions can raise `MethodError`, `InexactError`
or other conversion exceptions.

# Returns
A `KnotN5ComponentInvariant` with fields `component`, `representative`, `alexander`, `jones`.
Their individual meanings are given in the linked result-field entry.

# Guarantees
1. The returned fields equal the supplied values after conversion. Under Q6 they describe the finite graph or subordinate record specified by the linked definition, with exactly the enclosing run's population and bounds. [proved/agent: by construction]
2. Construction leaves the arguments unchanged; vectors, dictionaries and contained records already of the declared types remain shared with the result. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[classification labels](../specs/knot-classification.md#evidence-and-boundaries).

# Not guaranteed
The field constructor is not a validator or a classification computation.
A stored witness, resolved label or zero failure count is not independently
verified here. Distinct bounded components need not be distinct knot types.
The general classicality and knot-polynomial source identifications retain
the evidence boundaries stated in the specification. Mutable fields can
invalidate Q6 after construction.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Definitions and Q6 — finite graph and result semantics; no unrestricted knot classification claim.

# Examples
```jldoctest
julia> p = alexander_polynomial(pd_to_gauss(Int[]));

julia> r = KnotN5ComponentInvariant(1, 1, p, p);

julia> knot_polynomial_string(r.alexander)
"1"

```

Role: convert
"""
struct KnotN5ComponentInvariant
    component::Int
    representative::Int
    alexander::ExactKnotLaurentPolynomial
    jones::ExactKnotLaurentPolynomial
end

"""
    KnotN5Classification(max_n, barrier, state_counts, expected_state_counts, canonical_knot_type_counts, raw_rooted_classical_counts, raw_rooted_total, raw_rooted_classical_total, raw_rooted_classical_unresolved_count, states, state_index, edges, adjacency, components, component_roots, component_invariants, known_representatives, unresolved_components, direct_transition_count, duplicate_edge_transition_count, nonclassical_target_count, path_replay_count, path_replay_failures, edge_invariant_checks, edge_invariant_failures, fingerprint, elapsed_seconds) -> KnotN5Classification

Store a bounded signed-Gauss move graph and its classification diagnostics.

# Arguments
The field meanings are owned by [Result field meanings](../specs/knot-classification.md#result-field-meanings).
Q6 ([conditions](../specs/knot-classification.md#numbered-conditions)) is the
caller's responsibility. The field constructor does not check run consistency,
witness validity, counts or labels. Julia converts arguments to their declared
field types; incompatible conversions can raise `MethodError`, `InexactError`
or other conversion exceptions.

# Returns
A `KnotN5Classification` with fields `max_n`, `barrier`, `state_counts`, `expected_state_counts`, `canonical_knot_type_counts`, `raw_rooted_classical_counts`, `raw_rooted_total`, `raw_rooted_classical_total`, `raw_rooted_classical_unresolved_count`, `states`, `state_index`, `edges`, `adjacency`, `components`, `component_roots`, `component_invariants`, `known_representatives`, `unresolved_components`, `direct_transition_count`, `duplicate_edge_transition_count`, `nonclassical_target_count`, `path_replay_count`, `path_replay_failures`, `edge_invariant_checks`, `edge_invariant_failures`, `fingerprint`, `elapsed_seconds`.
Their individual meanings are given in the linked result-field entry.

# Guarantees
1. The returned fields equal the supplied values after conversion. Under Q6 they describe the finite graph or subordinate record specified by the linked definition, with exactly the enclosing run's population and bounds. [proved/agent: by construction]
2. Construction leaves the arguments unchanged; vectors, dictionaries and contained records already of the declared types remain shared with the result. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[classification labels](../specs/knot-classification.md#evidence-and-boundaries).

# Not guaranteed
The field constructor is not a validator or a classification computation.
A stored witness, resolved label or zero failure count is not independently
verified here. Distinct bounded components need not be distinct knot types.
The general classicality and knot-polynomial source identifications retain
the evidence boundaries stated in the specification. Mutable fields can
invalidate Q6 after construction.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Definitions and Q6 — finite graph and result semantics; no unrestricted knot classification claim.

# Examples
```jldoctest
julia> r = build_knot_n5_classification(0);

julia> Tuple(r.state_counts)
(1,)

```

Role: convert
"""
struct KnotN5Classification
    max_n::Int
    barrier::Int
    state_counts::Vector{Int}
    expected_state_counts::Vector{Int}
    canonical_knot_type_counts::Dict{String,Int}
    raw_rooted_classical_counts::Dict{String,Int}
    raw_rooted_total::Int
    raw_rooted_classical_total::Int
    raw_rooted_classical_unresolved_count::Int
    states::Vector{KnotGaussDiagram}
    state_index::Dict{KnotGaussDiagram,Int}
    edges::Vector{KnotN5Edge}
    adjacency::Vector{Vector{Int}}
    components::Vector{Int}
    component_roots::Vector{Int}
    component_invariants::Vector{KnotN5ComponentInvariant}
    known_representatives::Vector{KnotN5KnownRepresentative}
    unresolved_components::Vector{Int}
    direct_transition_count::Int
    duplicate_edge_transition_count::Int
    nonclassical_target_count::Int
    path_replay_count::Int
    path_replay_failures::Int
    edge_invariant_checks::Int
    edge_invariant_failures::Vector{String}
    fingerprint::String
    elapsed_seconds::Float64
end

Base.length(result::KnotN5Classification) = length(result.states)

"""Canonical key used in fingerprints and diagnostics."""
function _n5_state_key(d::KnotGaussDiagram)
    canonical_serialization(d)
end

function _n5_state_key_with_n(d::KnotGaussDiagram)
    "n=$(crossing_count(d))|$(_n5_state_key(d))"
end

function _n5_check_bound(max_n::Integer, barrier::Integer)
    m = Int(max_n)
    b = Int(barrier)
    0 <= m <= 5 || throw(ArgumentError("max_n must lie in 0:5"))
    m <= b <= 5 || throw(ArgumentError("require max_n <= barrier <= 5"))
    m, b
end

"""Enumerate and retain every Carter-genus-zero canonical state."""
function _n5_collect_states(max_n::Int)
    states = KnotGaussDiagram[]
    state_index = Dict{KnotGaussDiagram,Int}()
    counts = zeros(Int, max_n + 1)
    classical_cache = Dict{KnotGaussDiagram,Bool}()
    raw_classical_multiplicity = Dict{KnotGaussDiagram,Int}()
    raw_total = 0
    raw_classical_total = 0

    for n in 0:max_n
        for diagram in enumerate_knot_gauss(n; orientations=true)
            raw_total += 1
            canonical = canonicalize(diagram)
            # The source is classical by this branch.  Cache the canonical
            # decision so the same state is never re-certified repeatedly.
            is_classical = get!(classical_cache, canonical) do
                classical_gauss_realizable_fast(canonical)
            end
            is_classical || continue
            raw_classical_total += 1
            raw_classical_multiplicity[canonical] =
                get(raw_classical_multiplicity, canonical, 0) + 1
            if !haskey(state_index, canonical)
                push!(states, canonical)
                state_index[canonical] = length(states)
                counts[n + 1] += 1
            end
        end
    end
    states, state_index, counts, raw_classical_multiplicity, raw_total,
    raw_classical_total
end

function _n5_canonical_target(target::KnotGaussDiagram,
                              state_index::Dict{KnotGaussDiagram,Int},
                              classical_cache::Dict{KnotGaussDiagram,Bool})
    canonical = canonicalize(target)
    is_classical = get!(classical_cache, canonical) do
        classical_gauss_realizable_fast(canonical)
    end
    is_classical, canonical, get(state_index, canonical, 0)
end

"""Add an already verified transition to the simple undirected graph."""
function _n5_register_edge!(edges::Vector{KnotN5Edge},
                            adjacency::Vector{Vector{Int}},
                            edge_index::Dict{Tuple{Int,Int},Int},
                            source::Int, target::Int, move_type::Symbol,
                            candidate::ReidemeisterCandidate)
    pair = source <= target ? (source, target) : (target, source)
    haskey(edge_index, pair) && return false
    push!(edges, KnotN5Edge(source, target, move_type, candidate))
    edge_id = length(edges)
    edge_index[pair] = edge_id
    push!(adjacency[source], edge_id)
    source == target || push!(adjacency[target], edge_id)
    true
end

"""Build all deletion/R3 edges, with application and independent verification."""
function _n5_build_graph(states::Vector{KnotGaussDiagram},
                          state_index::Dict{KnotGaussDiagram,Int})
    adjacency = [Int[] for _ in eachindex(states)]
    edges = KnotN5Edge[]
    edge_index = Dict{Tuple{Int,Int},Int}()
    classical_cache = Dict{KnotGaussDiagram,Bool}(d => true for d in states)
    direct_count = 0
    duplicate_count = 0
    nonclassical_count = 0

    # Deletion edges are enough for R1/R2 insertions as explained in the file
    # docstring: every bounded insertion is the inverse of one enumerated
    # deletion on its canonical higher-crossing endpoint.
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
                    throw(ErrorException("N5 $move_type application failed at state " *
                                         "$source_index: " * sprint(showerror, error)))
                end
                verify_reidemeister_transition(source, raw_target, candidate) ||
                    throw(ErrorException("N5 $move_type transition verifier rejected " *
                                         "state $source_index"))
                classical, canonical_target, target_index =
                    _n5_canonical_target(raw_target, state_index, classical_cache)
                classical || begin
                    nonclassical_count += 1
                    continue
                end
                target_index == 0 && throw(ErrorException(
                    "verified classical target is outside the N<=barrier state set"))
                _n5_register_edge!(edges, adjacency, edge_index, source_index,
                                   target_index, Symbol(move_type), candidate) ||
                    (duplicate_count += 1)
                # Keep this assertion close to registration: canonical target
                # must really be the target state selected by the graph map.
                states[target_index] == canonical_target ||
                    throw(ErrorException("state-index canonical target mismatch"))
            end
        end
    end
    (edges=edges, adjacency=adjacency, direct_count=direct_count,
     duplicate_count=duplicate_count, nonclassical_count=nonclassical_count)
end

"""Connected components of the simple undirected retained graph."""
function _n5_components(adjacency::Vector{Vector{Int}}, edges::Vector{KnotN5Edge})
    count_states = length(adjacency)
    component = zeros(Int, count_states)
    roots = Int[]
    queue = Int[]
    for start in 1:count_states
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

function _n5_index_for_state(d::KnotGaussDiagram,
                             state_index::Dict{KnotGaussDiagram,Int})
    canonical = canonicalize(d)
    get(state_index, canonical, 0)
end

function _n5_pd_rows(value)
    isempty(strip(String(value))) && return Any[]
    value isa AbstractString && return JSON.parse(String(value))
    value
end

"""Read the five local KnotInfo PD rows without using their invariants."""
function _n5_load_knotinfo_rows(path::AbstractString)
    payload = JSON.parsefile(path)
    sheets = get(payload, "sheets", Any[])
    isempty(sheets) && throw(ArgumentError("KnotInfo snapshot has no sheet"))
    sheet = first(sheets)
    headers = String.(sheet["headers"])
    header_index = Dict(header => i for (i, header) in enumerate(headers))
    haskey(header_index, "name") || throw(ArgumentError("KnotInfo row has no name field"))
    rows = Dict{String,Any}()
    for row in sheet["rows"]
        values = row["values"]
        name = String(values[header_index["name"]])
        rows[name] = row
    end
    rows, header_index
end

function _n5_default_snapshot()
    normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260830", "data",
                      "knotinfo-n6.json"))
end

"""Construct base and mirrored KnotInfo representatives from PD notation."""
function _n5_known_representatives(state_index::Dict{KnotGaussDiagram,Int},
                                   component::Vector{Int};
                                   max_n::Int=5,
                                   snapshot_path::AbstractString=_n5_default_snapshot())
    rows, header_index = _n5_load_knotinfo_rows(snapshot_path)
    haskey(header_index, "pd_notation") ||
        throw(ArgumentError("KnotInfo snapshot has no pd_notation field"))
    representatives = KnotN5KnownRepresentative[]
    for name in KNOT_N5_KNOWN_NAMES
        haskey(rows, name) || throw(ArgumentError("KnotInfo snapshot lacks $name"))
        row = rows[name]
        diagram = pd_to_gauss(_n5_pd_rows(row["values"][header_index["pd_notation"]]))
        # Focused max_n<5 runs intentionally have only the seeds that fit in
        # their state set; the full barrier=5 run includes all ten labels.
        crossing_count(diagram) <= max_n || continue
        for mirror_flag in (false, true)
            selected = mirror_flag ? mirror(diagram) : diagram
            index = _n5_index_for_state(selected, state_index)
            index == 0 && throw(ErrorException(
                "known representative $name is absent from the classical state set"))
            label = mirror_flag ? "$(name)_mirror" : name
            push!(representatives,
                  KnotN5KnownRepresentative(
                      name, label, mirror_flag, selected, index, component[index],
                      alexander_polynomial(selected), jones_polynomial(selected)))
        end
    end
    representatives
end

function _n5_strict_knot_type(representative::KnotN5KnownRepresentative)
    representative.name in ("0_1", "4_1") && return representative.name
    "$(representative.name)^$(representative.mirror ? '-' : '+')"
end

"""Map move components to strict local-PD/mirror labels."""
function _n5_component_knot_types(known::Vector{KnotN5KnownRepresentative})
    result = Dict{Int,String}()
    for representative in known
        type = _n5_strict_knot_type(representative)
        previous = get(result, representative.component, nothing)
        previous === nothing || previous == type || throw(ErrorException(
            "known representatives disagree on component $(representative.component)"))
        result[representative.component] = type
    end
    result
end

function _n5_fixed_type_counts()
    Dict{String,Int}(type => 0 for type in KNOT_N5_STRICT_KNOT_TYPES)
end

"""Aggregate canonical and raw rooted counts after move components are known."""
function _n5_knot_type_counts(states::Vector{KnotGaussDiagram},
                              state_index::Dict{KnotGaussDiagram,Int},
                              components::Vector{Int},
                              known::Vector{KnotN5KnownRepresentative},
                              raw_classical_multiplicity::Dict{KnotGaussDiagram,Int})
    component_types = _n5_component_knot_types(known)
    canonical_counts = _n5_fixed_type_counts()
    for index in eachindex(states)
        type = get(component_types, components[index], nothing)
        type === nothing || (canonical_counts[type] += 1)
    end
    raw_counts = _n5_fixed_type_counts()
    unresolved = 0
    for (state, multiplicity) in raw_classical_multiplicity
        type = get(component_types, components[state_index[state]], nothing)
        if type === nothing
            unresolved += multiplicity
        else
            raw_counts[type] += multiplicity
        end
    end
    canonical_counts, raw_counts, unresolved
end

function _n5_component_invariants(states::Vector{KnotGaussDiagram},
                                  roots::Vector{Int})
    result = KnotN5ComponentInvariant[]
    for (component_id, state_index) in enumerate(roots)
        diagram = states[state_index]
        push!(result, KnotN5ComponentInvariant(component_id, state_index,
                                                alexander_polynomial(diagram),
                                                jones_polynomial(diagram)))
    end
    result
end

function _n5_invariant_pair(d::KnotGaussDiagram,
                            cache::Dict{KnotGaussDiagram,
                                        Tuple{ExactKnotLaurentPolynomial,
                                              ExactKnotLaurentPolynomial}})
    get!(cache, d) do
        (alexander_polynomial(d), jones_polynomial(d))
    end
end

"""Audit Alexander/Jones preservation on every retained edge.

The cache is local to this audit and is intentionally not a field of the
result.  Component representatives are retained; per-state values are not.
"""
function _n5_audit_edge_invariants(states::Vector{KnotGaussDiagram},
                                   edges::Vector{KnotN5Edge})
    cache = Dict{KnotGaussDiagram,
                 Tuple{ExactKnotLaurentPolynomial,
                       ExactKnotLaurentPolynomial}}()
    failures = String[]
    checks = 0
    for (edge_id, edge) in enumerate(edges)
        before = states[edge.source]
        after = states[edge.target]
        before_pair = _n5_invariant_pair(before, cache)
        after_pair = _n5_invariant_pair(after, cache)
        checks += 1
        before_pair[1] == after_pair[1] || push!(failures,
            "edge $edge_id $(edge.move_type): Alexander changed")
        before_pair[2] == after_pair[2] || push!(failures,
            "edge $edge_id $(edge.move_type): Jones changed")
    end
    checks, failures
end

"""Return the neighbor endpoint of an edge from one endpoint."""
function _n5_neighbor(edge::KnotN5Edge, current::Int)
    edge.source == current ? edge.target : edge.source
end

"""Find and replay the edge witness in a requested orientation."""
function _n5_replay_edge(result::KnotN5Classification, edge_id::Int,
                         from::Int, to::Int,
                         reverse_cache::Dict{Tuple{Int,Int},
                                             Tuple{ReidemeisterCandidate,
                                                   KnotGaussDiagram}})
    edge = result.edges[edge_id]
    source = result.states[from]
    if edge.source == from && edge.target == to
        candidate = edge.candidate
        raw = apply_reidemeister(source, candidate)
        verify_reidemeister_transition(source, raw, candidate) ||
            throw(ErrorException("stored direct edge failed replay verification"))
        canonicalize(raw) == result.states[to] ||
            throw(ErrorException("stored direct edge replay reached wrong state"))
        return raw
    end

    key = (from, to)
    if haskey(reverse_cache, key)
        candidate, expected_raw = reverse_cache[key]
        raw = apply_reidemeister(source, candidate)
        verify_reidemeister_transition(source, raw, candidate) ||
            throw(ErrorException("cached reverse edge failed replay verification"))
        canonicalize(raw) == result.states[to] ||
            throw(ErrorException("cached reverse edge reached wrong state"))
        return raw
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
        reverse_cache[key] = (candidate, raw)
        return raw
    end
    throw(ErrorException("could not find reverse witness for edge $edge_id " *
                         "from state $from to $to"))
end

"""Replay one parent-pointer tree for all states in one known component."""
function _n5_replay_known_component(result::KnotN5Classification,
                                    root::Int,
                                    component_id::Int)
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
            neighbor = _n5_neighbor(result.edges[edge_id], current)
            result.components[neighbor] == component_id || continue
            parent[neighbor] == 0 || continue
            parent[neighbor] = current
            parent_edge[neighbor] = edge_id
            push!(queue, neighbor)
        end
    end

    reverse_cache = Dict{Tuple{Int,Int},Tuple{ReidemeisterCandidate,
                                                KnotGaussDiagram}}()
    replay_count = 0
    failures = 0
    for state in eachindex(result.states)
        result.components[state] == component_id || continue
        parent[state] == 0 && (failures += 1; continue)
        current = state
        guard = 0
        while current != root
            guard += 1
            guard <= length(result.states) ||
                throw(ErrorException("parent-pointer cycle during N5 replay"))
            previous = parent[current]
            edge_id = parent_edge[current]
            try
                _n5_replay_edge(result, edge_id, current, previous, reverse_cache)
                replay_count += 1
            catch
                failures += 1
            end
            current = previous
        end
    end
    replay_count, failures
end

"""Replay every state in every component containing a known representative."""
function _n5_replay_known_paths(result::KnotN5Classification)
    roots_by_component = Dict{Int,Int}()
    for representative in result.known_representatives
        get!(roots_by_component, representative.component,
             representative.state_index)
    end
    total_replays = 0
    total_failures = 0
    for (component_id, root) in sort!(collect(roots_by_component))
        replays, failures = _n5_replay_known_component(result, root, component_id)
        total_replays += replays
        total_failures += failures
    end
    total_replays, total_failures
end

function _n5_unresolved_components(component::Vector{Int},
                                   known::Vector{KnotN5KnownRepresentative})
    known_components = Set(rep.component for rep in known)
    sort!(collect(setdiff(Set(component), known_components)))
end

"""Stable SHA-256 fingerprint of canonical states and retained edge pairs."""
function _n5_fingerprint(states::Vector{KnotGaussDiagram},
                         edges::Vector{KnotN5Edge})
    lines = String[]
    append!(lines, "state|$(_n5_state_key_with_n(d))\n" for d in states)
    edge_lines = String[]
    for edge in edges
        pair = edge.source <= edge.target ? (edge.source, edge.target) :
            (edge.target, edge.source)
        push!(edge_lines, "edge|$(pair[1])|$(pair[2])|$(edge.move_type)\n")
    end
    append!(lines, sort!(edge_lines))
    payload = join(lines)
    bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(payload))))
end

"""
    build_knot_n5_classification(max_n::Integer=5; barrier::Integer=max_n, snapshot_path=_n5_default_snapshot()) -> KnotN5Classification

Construct the finite move graph of signed Gauss diagrams accepted by the Carter genus-zero test.

# Arguments
`max_n` and `barrier` must satisfy Q7 ([conditions](../specs/knot-classification.md#numbered-conditions)); the range is checked and `ArgumentError` names the violated bound. Integer conversion errors can propagate.
`snapshot_path` must satisfy Q9 on that page. It defaults to the repository's
`threads/PJ260830/data/knotinfo-n6.json`. Missing files, invalid JSON and PD
conversion failures propagate; missing names or columns raise `ArgumentError`.
Consistency, witness, polynomial or replay failures raise `ErrorException`
with the failed audit identified. Snapshot name identification is the caller's
responsibility; without it a connected label is only a snapshot label.

# Returns
A `KnotN5Classification` with the fields described individually in
[Result field meanings](../specs/knot-classification.md#result-field-meanings).

# Guarantees
1. States form the canonical population through `max_n` retained by the Carter test; edges are replay-verified transitions between those states and components describe connectivity in that finite graph. [proved/agent: by construction]
2. Resolved labels are attached through known representatives in the same component. Components with no such representative remain unresolved. [proved/agent: by construction]
3. Successful return reports no failures in the path replay and edge-polynomial checks. [proved/agent: by construction]
4. Inputs and the snapshot file are not modified. The result is in memory; no database is written. [proved/agent: by construction]

# Conventions
[Canonical coordinates](../specs/knot-gauss.md#canonical-coordinates) and
[classification labels](../specs/knot-classification.md#evidence-and-boundaries).

# Not guaranteed
A larger `barrier` alone does not enumerate states above `max_n`. Distinct
bounded components need not be distinct knot types; no unrestricted isotopy
decision follows. Classicality, polynomial and seed-name source identifications
retain the specification's evidence boundaries. The finite graph does not
retain complete paths or certify the external snapshot's naming convention.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Bounded graph and Classification record — finite population and connectivity meanings.

# Examples
```jldoctest
julia> r = build_knot_n5_classification(0);

julia> Tuple(r.state_counts)
(1,)
```

Role: construct
"""
function build_knot_n5_classification(max_n::Integer=5;
                                      barrier::Integer=max_n,
                                      snapshot_path::AbstractString=_n5_default_snapshot())
    max_n, barrier = _n5_check_bound(max_n, barrier)
    started = time_ns()
    states, state_index, state_counts, raw_classical_multiplicity, raw_total,
    raw_classical_total = _n5_collect_states(max_n)
    expected_counts = collect(KNOT_N5_EXPECTED_STATE_COUNTS[1:max_n + 1])
    state_counts == expected_counts || throw(ErrorException(
        "classical canonical state counts mismatch: got $(state_counts), " *
        "expected $(expected_counts)"))

    graph = _n5_build_graph(states, state_index)
    components, roots = _n5_components(graph.adjacency, graph.edges)
    known = _n5_known_representatives(state_index, components;
                                      max_n=max_n,
                                      snapshot_path=snapshot_path)
    canonical_type_counts, raw_type_counts, raw_unresolved =
        _n5_knot_type_counts(states, state_index, components, known,
                             raw_classical_multiplicity)
    component_invariants = _n5_component_invariants(states, roots)
    edge_checks, invariant_failures =
        _n5_audit_edge_invariants(states, graph.edges)
    isempty(invariant_failures) || throw(ErrorException(
        "Alexander/Jones edge-preservation audit failed: " *
        join(invariant_failures, "; ")))

    unresolved = _n5_unresolved_components(components, known)
    # Paths are reconstructed through transient parent pointers only.  They
    # are not fields in the result and are never persisted.
    provisional = KnotN5Classification(
        max_n, barrier, state_counts, expected_counts, canonical_type_counts,
        raw_type_counts, raw_total, raw_classical_total, raw_unresolved,
        states, state_index,
        graph.edges, graph.adjacency, components, roots, component_invariants,
        known, unresolved, graph.direct_count, graph.duplicate_count,
        graph.nonclassical_count, 0, 0, edge_checks, invariant_failures, "",
        0.0)
    path_count, path_failures = _n5_replay_known_paths(provisional)
    path_failures == 0 || throw(ErrorException(
        "known-representative path replay failed for $path_failures transitions"))
    fingerprint = _n5_fingerprint(states, graph.edges)
    elapsed = (time_ns() - started) / 1e9
    KnotN5Classification(
        max_n, barrier, state_counts, expected_counts, canonical_type_counts,
        raw_type_counts, raw_total, raw_classical_total, raw_unresolved,
        states, state_index,
        graph.edges, graph.adjacency, components, roots, component_invariants,
        known, unresolved, graph.direct_count, graph.duplicate_count,
        graph.nonclassical_count, path_count, path_failures, edge_checks,
        invariant_failures, fingerprint, elapsed)
end

"""
    classify_knot_n5(max_n::Integer=5; barrier::Integer=max_n, snapshot_path=_n5_default_snapshot()) -> KnotN5Classification

Return a bounded classification of classical signed Gauss representations by verified diagram moves.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/knot-classification.md).

`max_n` and `barrier` must convert to `Int` and satisfy `0 <= max_n <= barrier <= 5`;
otherwise `ArgumentError` reports the violated range. `snapshot_path` is a
readable local KnotInfo JSON snapshot accepted by
[`build_knot_n5_classification`](@ref). Snapshot parsing, missing data and
internal consistency failures propagate. Unsupported argument forms fail
by dispatch in that entry point.

# Returns
`KnotN5Classification` with canonical `states`, their index, verified `edges`,
adjacency and component data; known representatives and component invariants;
raw and canonical counts; unresolved components; replay and invariant-audit
outcomes; a graph fingerprint and elapsed time.

# Guarantees
1. The returned graph concerns the classical signed representations with
   at most `max_n` crossings. Edges are retained local Reidemeister moves;
   components give finite move-connectivity evidence among those states. [proved/agent: by construction]
2. The result has the same semantics as `build_knot_n5_classification`.
   Known type labels come from connected snapshot representatives; other
   components remain unresolved. [proved/agent: by construction]
3. Successful return has no reported path-replay or edge-invariant failures. [proved/agent: by construction]
4. Arguments and the input snapshot are not modified. [proved/agent: by construction]

# Not guaranteed
Setting `barrier > max_n` does not add states above `max_n`. Distinct finite
components are not thereby proved to be distinct knot types. The result
does not retain full paths, certify the snapshot's naming conventions or
establish an unbounded equivalence decision.

# Sources
cf: lib/topology.jl/src/reidemeister.jl verify_reidemeister_transition — independent local checks on retained graph edges.
cf: lib/topology.jl/src/classical_gauss_realizability.jl classical_gauss_certificate — classical-realizability predicate used to restrict the diagram state set.

# Examples
```jldoctest
julia> r = classify_knot_n5(0);

julia> Tuple(r.state_counts)
(1,)
```

See also [`build_knot_n5_classification`](@ref), [`knot_n5_summary`](@ref).

Role: construct
"""
classify_knot_n5(args...; kwargs...) =
    build_knot_n5_classification(args...; kwargs...)

"""
    knot_n5_summary(result::KnotN5Classification) -> Dict{String,Any}

Return a compact dictionary summary of a finite knot-move result.

# Arguments
`result` satisfies Q6 ([conditions](../specs/knot-classification.md#numbered-conditions))
for the intended interpretation. This is the caller's responsibility; the
summary does not replay or validate a run. Inconsistent manually constructed
records can produce inconsistent counts rather than a validation error.

# Returns
A new `Dict{String,Any}` with keys
`version`, `max_n`, `barrier`, `state_counts`, `canonical_knot_type_counts`, `raw_rooted_total`, `raw_rooted_classical_total`, `raw_rooted_classical_counts`, `raw_rooted_classical_unresolved_count`, `total_states`, `edge_count`, `direct_transition_count`, `duplicate_edge_transition_count`, `component_count`, `known_representative_count`, `known_component_count`, `known_reachable_state_count`, `unresolved_component_count`, `unresolved_components`, `path_replay_count`, `path_replay_failures`, `edge_invariant_checks`, `edge_invariant_failures`, `fingerprint_sha256`, `elapsed_seconds`.
Same-named result fields have their meanings in
[Result field meanings](../specs/knot-classification.md#result-field-meanings).
`version` identifies the summary's implementation schema;
`fingerprint_sha256` is the stored `fingerprint`. Where present,
`total_states`, `edge_count`, `component_count`, `known_representative_count`,
`known_component_count` and `unresolved_component_count` count the corresponding
collections. Unresolved state count counts states in listed unresolved
components; known-reachable state count is its complement in the population.
Component-keyed output dictionaries use decimal string keys.

# Guarantees
1. The summary reports the supplied finite result's data and derived collection counts; it performs no new mathematical classification. [proved/agent: by construction]
2. The input is unchanged. Returned vectors and dictionaries are copied, including the nested label lists and closure groups when present. [proved/agent: by construction]

# Not guaranteed
The dictionary is not a self-contained replay certificate and contains no
full path or Carter surface data. It does not strengthen the source result's
classicality, knot-type or completeness evidence. No cross-version schema
compatibility is promised.

# Sources
cf:lib/topology.jl/docs/specs/knot-classification.md, Result field meanings — finite-run summary interpretation.

# Examples
```jldoctest
julia> r = build_knot_n5_classification(0);

julia> knot_n5_summary(r)["total_states"]
1

```

Role: convert
"""
function knot_n5_summary(result::KnotN5Classification)
    known_components = sort!(unique(rep.component for rep in result.known_representatives))
    unresolved_state_count = sum(count(==(component_id), result.components)
                                 for component_id in result.unresolved_components;
                                 init=0)
    Dict{String,Any}(
        "version" => KNOT_N5_CLASSIFICATION_VERSION,
        "max_n" => result.max_n,
        "barrier" => result.barrier,
        "state_counts" => copy(result.state_counts),
        "canonical_knot_type_counts" => copy(result.canonical_knot_type_counts),
        "raw_rooted_total" => result.raw_rooted_total,
        "raw_rooted_classical_total" => result.raw_rooted_classical_total,
        "raw_rooted_classical_counts" => copy(result.raw_rooted_classical_counts),
        "raw_rooted_classical_unresolved_count" =>
            result.raw_rooted_classical_unresolved_count,
        "total_states" => length(result.states),
        "edge_count" => length(result.edges),
        "direct_transition_count" => result.direct_transition_count,
        "duplicate_edge_transition_count" => result.duplicate_edge_transition_count,
        "component_count" => length(result.component_roots),
        "known_representative_count" => length(result.known_representatives),
        "known_component_count" => length(known_components),
        "known_reachable_state_count" => length(result.states) - unresolved_state_count,
        "unresolved_component_count" => length(result.unresolved_components),
        "unresolved_components" => copy(result.unresolved_components),
        "path_replay_count" => result.path_replay_count,
        "path_replay_failures" => result.path_replay_failures,
        "edge_invariant_checks" => result.edge_invariant_checks,
        "edge_invariant_failures" => copy(result.edge_invariant_failures),
        "fingerprint_sha256" => result.fingerprint,
        "elapsed_seconds" => result.elapsed_seconds,
    )
end
