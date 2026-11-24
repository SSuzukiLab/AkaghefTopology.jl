# cf:[[BP97]] §3.1 Cor. 3.1.7, §3.2 Prop. 3.2.1 Figs. 3.4–3.6

"""Coloured o-graph for the converse of Prop. 3.2.1.

Slots `1:4` run counter-clockwise. `edge_darts[e] = (tail, head)` names the
two ends and does not orient the colour. `colours[e]` is the one Z₃ element
on the undirected edge. P2 of a candidate arrow is compared with that element.
The opposite arrow needs the negated element, stored as its own tuple.
This is not the shared coloured type; `coloured_ograph` adapts that type.
"""
struct BranchColourGraph
    nvertices::Int
    over_slots::Vector{NTuple{2,Int}}
    edge_darts::Vector{NTuple{2,Int}}
    colours::Vector{Int}
end

"""O0–O3 for Def. 2.1.4. O2 is the combinatorial normal crossing only."""
struct OGraphConditionReport
    O0::Bool
    O1::Bool
    O2::Bool
    O3::Bool
    failures::Vector{String}
end

"""One straight-through orientation. Each circuit is a cycle of outgoing darts."""
struct P1Orientation
    circuits::Vector{Vector{Int}}
end

"""Oriented normal diagram, one word per straight-through circuit.

`vertex_indices` follows Figs. 3.4 and 3.5 for this orientation. Reversing
every circuit keeps the indices. Reversing one circuit of a crossing where
the two strands lie on different circuits exchanges `+1` and `−1`.
"""
struct NormalDiagram
    components::Vector{Vector{Int}}
    vertex_indices::Vector{Int}
end

"""Normal diagrams recovered from a coloured o-graph, plus the disc gap.

`states` keeps the diagrams with one circuit, because `OGraphState` is a
one-circuit presentation. `diagrams` keeps every compatible orientation,
including several circuits. `discs == :unread_fig25` means the four T's of
Fig. 2.5 (left) have not been matched to the four ports, so disc orientations
are not constructed. Edge directions and vertex indices do not use that matching.
"""
struct CompatibleBranching
    states::Vector{OGraphState}
    diagrams::Vector{NormalDiagram}
    discs::Symbol
    provisional::Bool
end

const FIG25_PORT_ASSIGNMENT = :unread
const PRO_REVIEW_MOMENT1 = "threads/PJ261120/sources/pro/branching-20261122-000145-answer.md"
const PRO_REVIEW_MOMENT1_STATUS = :read
const PRO_REVIEW_MOMENT3 = "threads/PJ261120/sources/pro/branching-20261122-013056.md"
const PRO_REVIEW_MOMENT3_STATUS = :read

function conditions_passed(report::OGraphConditionReport)
    report.O0 && report.O1 && report.O2 && report.O3
end

"""P2. Above is an over-passage and below is an under-passage.

P2a both ends the same → 0, P2b first below and second above → 1,
P2c first above and second below → −1. Swapping the ends swaps P2b and P2c.
"""
function p2_colour(first_above::Bool, second_above::Bool)::Int
    first_above == second_above && return 0
    return first_above ? -1 : 1
end

"""Representative in `{-1,0,1}` of an integer class in Z₃. `2` is `−1`."""
function colour_mod3(value::Integer)::Int
    reduced = mod(Int(value), 3)
    reduced == 2 ? -1 : reduced
end

dart_id(vertex::Int, slot::Int) = (vertex - 1) * 4 + slot
dart_vertex(dart::Int) = (dart - 1) ÷ 4 + 1
dart_slot(dart::Int) = mod1(dart, 4)

function opposite_slot(slot::Int)
    shifted = slot + 2
    return shifted > 4 ? slot - 2 : shifted
end

opposite_dart(dart::Int) = dart_id(dart_vertex(dart), opposite_slot(dart_slot(dart)))

function is_over(graph::BranchColourGraph, dart::Int)
    slot = dart_slot(dart)
    left, right = graph.over_slots[dart_vertex(dart)]
    slot == left || slot == right
end

"""Port slot in the counter-clockwise frame of Figs. 3.4 and 3.5.

Sign `+1` is Fig. 3.4: slots `(−i, +o, −o, +i)`.
Sign `−1` is Fig. 3.5: slots `(−i, +i, −o, +o)`.
"""
function port_slot(vertex_sign::Int, port::Symbol)
    vertex_sign in (1, -1) || throw(ArgumentError(
        "IDX: vertex index must be the Fig. 3.4 or Fig. 3.5 picture"))
    if vertex_sign == 1
        port === :under_in && return 1
        port === :over_out && return 2
        port === :under_out && return 3
        port === :over_in && return 4
    else
        port === :under_in && return 1
        port === :over_in && return 2
        port === :under_out && return 3
        port === :over_out && return 4
    end
    throw(ArgumentError("unknown port $port"))
end

"""Forward half of Prop. 3.2.1, from a one-circuit normal o-graph.

The reference direction of each edge is the step `g[k] → g[k+1]`.
"""
function coloured_ograph_from_normal(state::OGraphState)
    nvertices = length(state.vertex_signs)
    gauss = collect(state.gauss)
    nedges = length(gauss)
    edge_darts = Vector{NTuple{2,Int}}(undef, nedges)
    colours = Vector{Int}(undef, nedges)
    for index in 1:nedges
        departure = gauss[index]
        arrival = gauss[mod1(index + 1, nedges)]
        departure_vertex = abs(departure)
        arrival_vertex = abs(arrival)
        tail_port = departure > 0 ? :over_out : :under_out
        head_port = arrival > 0 ? :over_in : :under_in
        tail = dart_id(departure_vertex, port_slot(state.vertex_signs[departure_vertex], tail_port))
        head = dart_id(arrival_vertex, port_slot(state.vertex_signs[arrival_vertex], head_port))
        edge_darts[index] = (tail, head)
        colours[index] = p2_colour(departure > 0, arrival > 0)
    end
    over_slots = fill((2, 4), nvertices)
    BranchColourGraph(nvertices, over_slots, edge_darts, colours)
end

function _mate_map(graph::BranchColourGraph)
    mate = Dict{Int,Int}()
    for (tail, head) in graph.edge_darts
        mate[tail] = head
        mate[head] = tail
    end
    mate
end

function _opposite_pair(slots::NTuple{2,Int})
    left, right = slots
    pair = (min(left, right), max(left, right))
    pair == (1, 3) || pair == (2, 4)
end

"""Def. 2.1.4 as O0–O3. Failure strings name the condition id."""
function ograph_conditions(graph::BranchColourGraph)
    failures = String[]
    nvertices = graph.nvertices
    if nvertices < 1
        push!(failures, "O1: an o-graph has at least one vertex")
    end
    nedges = length(graph.edge_darts)
    if length(graph.colours) != nedges
        push!(failures, "O3: each edge needs one colour")
    end
    if length(graph.over_slots) != nvertices
        push!(failures, "O2: each vertex needs one over pair")
    end
    degree = zeros(Int, max(nvertices, 0))
    seen = Dict{Int,Int}()
    for (edge_index, (tail, head)) in enumerate(graph.edge_darts)
        for dart in (tail, head)
            if dart < 1 || (nvertices > 0 && dart > 4 * nvertices)
                push!(failures, "O0: edge $edge_index meets a dart outside the four ports")
                continue
            end
            seen[dart] = get(seen, dart, 0) + 1
            vertex = dart_vertex(dart)
            if 1 <= vertex <= nvertices
                degree[vertex] += 1
            end
        end
        if edge_index <= length(graph.colours) && !(graph.colours[edge_index] isa Integer)
            push!(failures, "O3: edge $edge_index has no Z3 colour")
        end
    end
    if nvertices >= 1
        for dart in 1:(4 * nvertices)
            count = get(seen, dart, 0)
            count == 1 ||
                push!(failures, "O0: port dart $dart occurs $count times")
        end
        for vertex in 1:nvertices
            degree[vertex] == 4 ||
                push!(failures, "O0: vertex $vertex is not quadrivalent")
        end
        if length(graph.over_slots) == nvertices
            for vertex in 1:nvertices
                _opposite_pair(graph.over_slots[vertex]) ||
                    push!(failures, "O2: vertex $vertex over pair is not an opposite pair")
            end
        end
        if nvertices > 0 && all(==(4), degree) &&
                !any(text -> startswith(text, "O0:"), failures)
            stack = [1]
            visited = falses(nvertices)
            nvertices >= 1 && (visited[1] = true)
            while !isempty(stack)
                vertex = pop!(stack)
                for (tail, head) in graph.edge_darts
                    for dart in (tail, head)
                        dart_vertex(dart) == vertex || continue
                        other = dart_vertex(dart == tail ? head : tail)
                        if 1 <= other <= nvertices && !visited[other]
                            visited[other] = true
                            push!(stack, other)
                        end
                    end
                end
            end
            all(visited) || push!(failures, "O0: the graph is not connected")
        end
    end
    named(prefix) = any(failure -> startswith(failure, prefix), failures)
    OGraphConditionReport(!named("O0:"), !named("O1:"), !named("O2:"), !named("O3:"),
                          failures)
end

function _outgoing_cycles(graph::BranchColourGraph)
    mate = _mate_map(graph)
    ndarts = 4 * graph.nvertices
    seen = falses(ndarts)
    cycles = Vector{Vector{Int}}()
    for start in 1:ndarts
        seen[start] && continue
        cycle = Int[]
        dart = start
        while !seen[dart]
            seen[dart] = true
            push!(cycle, dart)
            haskey(mate, dart) || throw(ArgumentError(
                "O0: dart $dart is not incident to an edge"))
            dart = opposite_dart(mate[dart])
        end
        push!(cycles, cycle)
    end
    cycles
end

function _direction_pairs(graph::BranchColourGraph)
    cycles = _outgoing_cycles(graph)
    by_set = Dict{Set{Int},Vector{Int}}(Set(cycle) => cycle for cycle in cycles)
    pairs = Vector{NTuple{2,Vector{Int}}}()
    used = Set{Set{Int}}()
    for cycle in cycles
        key = Set(cycle)
        key in used && continue
        reverse_key = Set(opposite_dart(dart) for dart in cycle)
        reverse_key == key && throw(ArgumentError(
            "O2: a straight-through walk meets both ends of one strand"))
        haskey(by_set, reverse_key) || throw(ArgumentError(
            "O2: the opposite direction of a walk is missing"))
        push!(pairs, (cycle, by_set[reverse_key]))
        push!(used, key)
        push!(used, reverse_key)
    end
    pairs
end

"""P1 orientations: both directions of every straight-through circuit."""
function p1_orientations(graph::BranchColourGraph)
    pairs = _direction_pairs(graph)
    orientations = P1Orientation[]
    ncircuits = length(pairs)
    for mask in 0:(2^ncircuits - 1)
        circuits = Vector{Int}[]
        for index in 1:ncircuits
            bit = (mask >> (index - 1)) & 1
            push!(circuits, bit == 0 ? pairs[index][1] : pairs[index][2])
        end
        push!(orientations, P1Orientation(circuits))
    end
    orientations
end

function _outgoing_set(orientation::P1Orientation)
    outgoing = Set{Int}()
    for circuit in orientation.circuits
        for dart in circuit
            push!(outgoing, dart)
        end
    end
    outgoing
end

"""Whether this P1 orientation induces the stored undirected colours.

Def. 2.1.4 puts one Z₃ element on the edge. P2 reads that element from
the orientation's first and second ends. The opposite orientation reads
the negated value, so it matches the negated colour tuple, not the same one.
`edge_darts` only names the two ends; it does not orient the colour.
"""
function satisfies_p2(graph::BranchColourGraph, orientation::P1Orientation)
    outgoing = _outgoing_set(orientation)
    for (edge_index, (tail, head)) in enumerate(graph.edge_darts)
        colour = colour_mod3(graph.colours[edge_index])
        tail_out = tail in outgoing
        head_out = head in outgoing
        if tail_out && !head_out
            p2_colour(is_over(graph, tail), is_over(graph, head)) == colour ||
                return false
        elseif head_out && !tail_out
            p2_colour(is_over(graph, head), is_over(graph, tail)) == colour ||
                return false
        else
            return false
        end
    end
    true
end

function compatible_orientations(graph::BranchColourGraph)
    filter(orientation -> satisfies_p2(graph, orientation), p1_orientations(graph))
end

admits_compatible_branching(graph::BranchColourGraph) =
    !isempty(compatible_orientations(graph))

"""Index +1 when the port counter-clockwise from the incoming under-passage
is the outgoing over-passage (Fig. 3.4). Index −1 when that port is the
incoming over-passage (Fig. 3.5).
"""
function vertex_index(graph::BranchColourGraph, vertex::Int, outgoing::Set{Int})
    over_left, over_right = graph.over_slots[vertex]
    over = (over_left, over_right)
    incoming_under = Int[]
    for slot in 1:4
        (slot == over_left || slot == over_right) && continue
        dart_id(vertex, slot) in outgoing && continue
        push!(incoming_under, slot)
    end
    length(incoming_under) == 1 || throw(ArgumentError(
        "P1: vertex $vertex has $(length(incoming_under)) incoming under-passages"))
    next_slot = mod1(only(incoming_under) + 1, 4)
    next_slot in over || throw(ArgumentError(
        "O2: vertex $vertex is not an alternating normal crossing"))
    return dart_id(vertex, next_slot) in outgoing ? 1 : -1
end

function recovered_normal_diagrams(graph::BranchColourGraph)
    diagrams = NormalDiagram[]
    for orientation in compatible_orientations(graph)
        outgoing = _outgoing_set(orientation)
        components = Vector{Vector{Int}}()
        for circuit in orientation.circuits
            push!(components, Int[is_over(graph, dart) ? dart_vertex(dart) : -dart_vertex(dart)
                                  for dart in circuit])
        end
        signs = Int[vertex_index(graph, vertex, outgoing) for vertex in 1:graph.nvertices]
        push!(diagrams, NormalDiagram(components, signs))
    end
    diagrams
end

function _state_from_diagram(diagram::NormalDiagram, convention_id::AbstractString)
    length(diagram.components) == 1 || return nothing
    OGraphState(only(diagram.components), diagram.vertex_indices; convention_id=convention_id)
end

"""One-circuit normal o-graphs whose P2 colours are the stored colours.

Both directions are returned when both fit in `OGraphState`. Canonicalization
is left to the caller: it does not identify a global reversal.
"""
function recovered_normal_states(graph::BranchColourGraph;
                                 convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION)
    states = OGraphState[]
    for diagram in recovered_normal_diagrams(graph)
        state = _state_from_diagram(diagram, convention_id)
        state === nothing && continue
        push!(states, state)
    end
    states
end

function recovered_canonical_states(graph::BranchColourGraph; kwargs...)
    unique(canonicalize_ograph, recovered_normal_states(graph; kwargs...))
end

"""Prop. 3.2.1 data that this lane constructs.

`states` is empty when every compatible orientation has two or more circuits.
Those orientations stay in `diagrams`. Disc orientations stay unread until the
Fig. 2.5 port chart is fixed. `provisional` is true for that reason even when
a normal o-graph is recovered.
"""
function compatible_branching(graph::BranchColourGraph;
                              convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION)
    diagrams = recovered_normal_diagrams(graph)
    states = OGraphState[]
    for diagram in diagrams
        state = _state_from_diagram(diagram, convention_id)
        state === nothing && continue
        push!(states, state)
    end
    CompatibleBranching(states, diagrams, :unread_fig25, true)
end

function shift_colour(graph::BranchColourGraph, edge::Int)
    colours = copy(graph.colours)
    reduced = colour_mod3(colours[edge])
    colours[edge] = reduced == 1 ? -1 : reduced + 1
    BranchColourGraph(graph.nvertices, copy(graph.over_slots), graph.edge_darts, colours)
end

function negate_colours(graph::BranchColourGraph)
    BranchColourGraph(graph.nvertices, copy(graph.over_slots), graph.edge_darts,
                      [colour_mod3(-colour) for colour in graph.colours])
end

"""Count reference-direction colour tuples that satisfy P2 for some P1 orientation.

The census is the set `{-1,0,1}` to the power of the edge count. It is
refused above 10 edges.
"""
function branching_colour_census(graph::BranchColourGraph)
    nedges = length(graph.edge_darts)
    nedges <= 10 || throw(ArgumentError(
        "colour census stops at 10 edges; this graph has $nedges"))
    orientations = p1_orientations(graph)
    values = (-1, 0, 1)
    branching = 0
    total = 0
    for index in CartesianIndices(ntuple(_ -> 3, nedges))
        colours = [values[index[axis]] for axis in 1:nedges]
        coloured = BranchColourGraph(graph.nvertices, graph.over_slots,
                                      graph.edge_darts, colours)
        any(orientation -> satisfies_p2(coloured, orientation), orientations) &&
            (branching += 1)
        total += 1
    end
    (total=total, branching=branching, nonbranching=total - branching)
end

"""Read a `ColouredOGraph` without using its normal-o-graph constructor.

The shared type's O0–O3 checks must already pass. Half-edges are renumbered
to the counter-clockwise slots of this file. The stored colour is the
undirected Z₃ element and is not given a direction.
"""
function coloured_ograph_from_core(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "coloured o-graph fails " * join((check.message for check in failed), "; ")))
    nvertices = length(graph.vertex_halfedges)
    halfedge_dart = Dict{Int,Int}()
    over_slots = Vector{NTuple{2,Int}}(undef, nvertices)
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for (slot, halfedge) in enumerate(cyclic)
            halfedge_dart[halfedge] = dart_id(vertex, slot)
        end
        positions = Tuple(Int(position) for position in graph.over_positions[vertex])
        length(positions) == 2 || throw(ArgumentError(
            "O2: vertex $vertex over pair does not have two positions"))
        over_slots[vertex] = (positions[1], positions[2])
    end
    edge_darts = NTuple{2,Int}[]
    colours = Int[]
    seen = Set{Int}()
    for record in graph.edge_colours
        record.halfedge in seen && continue
        other = graph.mate[record.halfedge]
        push!(seen, record.halfedge)
        push!(seen, other)
        push!(edge_darts, (halfedge_dart[record.halfedge], halfedge_dart[other]))
        push!(colours, colour_mod3(record.colour))
    end
    BranchColourGraph(nvertices, over_slots, edge_darts, colours)
end

function _adapt_core_coloured(value)
    if isdefined(Topology, :ColouredOGraph) && value isa ColouredOGraph
        return coloured_ograph_from_core(value)
    end
    throw(ArgumentError(
        "O0: pass an OGraphState, a BranchColourGraph, or a ColouredOGraph"))
end

"""Adapt a normal presentation, this lane's graph, or a later core coloured type."""
function coloured_ograph(value)
    value isa BranchColourGraph && return value
    value isa OGraphState && return coloured_ograph_from_normal(value)
    _adapt_core_coloured(value)
end

"""One-circuit presentation read in the opposite direction.

Both strands reverse, so the vertex indices stay. A diagram with two or
more circuits is not an `OGraphState`.
"""
function reversed_normal_presentation(state::OGraphState)
    OGraphState(reverse(collect(state.gauss)), collect(state.vertex_signs);
                convention_id=state.convention_id)
end
