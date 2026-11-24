# cf:[[BP97#^bp97-sec-0566]] Proposition 3.2.1 — colours and the orientation witness of a closed normal o-graph; adjustment: above = positive passage joins the repository Gauss code to the over branch; unverified

"""One directed edge in the Gauss-word orientation of a normal o-graph.

`tail` is the first end and `head` is the second end. The colour stays on the
undirected edge. This direction is not a field of `ColouredOGraph`.
"""
struct NormalEdgeOrientation
    tail::Int
    head::Int
end

"""A closed normal o-graph converted by Proposition 3.2.1.

`colours` follows `real_edge_table`: index `k` is the edge whose head lies at
word position `k`. `state` is the presentation that was converted. The stored
colour is the value Proposition 3.2.1 assigns to `edges`. Figure 2.5 does not
say that the second end is the right-hand T, so this conversion does not place
the head on the right and does not negate the colour.

`witness` is the existing C1–C3 replay. It is not an N1 or N2 proof.
Face cycles are `face_cycles`, under `OGRAPH_FACE_CYCLE_STATUS`.
"""
struct ClosedNormalColouring
    graph::ColouredOGraph
    edges::Tuple{Vararg{NormalEdgeOrientation}}
    colours::Tuple{Vararg{Int}}
    witness::ClosednessWitness
    state::OGraphState
end

"""Colour in `{-1, 0, 1}` from the above/below marks of the two ends.

Above is the over branch. In this adapter a positive Gauss passage is that
branch. Both marks equal gives 0. First below and second above gives 1.
First above and second below gives −1.
"""
function prop321_colour(first_above::Bool, second_above::Bool)
    first = first_above ? 0 : 1
    second = second_above ? 0 : 1
    _colour_mod3(first - second)
end

function _normal_port_slot(vertex_sign::Int, port::Symbol)
    vertex_sign == 1 && port === :mi && return 1
    vertex_sign == 1 && port === :po && return 2
    vertex_sign == 1 && port === :mo && return 3
    vertex_sign == 1 && port === :pi && return 4
    vertex_sign == -1 && port === :mi && return 1
    vertex_sign == -1 && port === :pi && return 2
    vertex_sign == -1 && port === :mo && return 3
    vertex_sign == -1 && port === :po && return 4
    throw(ArgumentError("vertex sign and port do not name a Fig. 3.4 or Fig. 3.5 slot"))
end

_normal_halfedge(vertex::Int, slot::Int) = (vertex - 1) * 4 + slot

function _normal_pieces(state::OGraphState)
    n = true_vertex_count(state)
    n >= 1 || throw(ArgumentError("N1 failed: an o-graph has at least one vertex"))
    link = _ograph_virtual_link(state)
    table = real_edge_table(link)
    length(table) == length(state.gauss) || throw(ArgumentError(
        "real edge table does not have one edge per Gauss entry"))
    vertices = Tuple(ntuple(slot -> _normal_halfedge(vertex, slot), 4) for vertex in 1:n)
    overs = Tuple((2, 4) for _ in 1:n)
    mate = zeros(Int, 4n)
    edges = NormalEdgeOrientation[]
    colours = Int[]
    for edge in table
        tail_label, head_label = edge.crossing
        (tail_label == 0 || head_label == 0) && throw(ArgumentError(
            "a closed normal colouring has no cut endpoint"))
        tail_vertex = abs(tail_label)
        head_vertex = abs(head_label)
        tail_port = tail_label > 0 ? :po : :mo
        head_port = head_label > 0 ? :pi : :mi
        tail = _normal_halfedge(tail_vertex, _normal_port_slot(state.vertex_signs[tail_vertex], tail_port))
        head = _normal_halfedge(head_vertex, _normal_port_slot(state.vertex_signs[head_vertex], head_port))
        (tail == head || mate[tail] != 0 || mate[head] != 0) && throw(ArgumentError(
            "Gauss passages do not give four distinct half-edges at each vertex"))
        mate[tail] = head
        mate[head] = tail
        push!(edges, NormalEdgeOrientation(tail, head))
        push!(colours, prop321_colour(tail_label > 0, head_label > 0))
    end
    all(!iszero, mate) || throw(ArgumentError("a half-edge was left unpaired"))
    records = Tuple(EdgeColour(edge.tail, colour) for (edge, colour) in zip(edges, colours))
    graph = ColouredOGraph(vertices, Tuple(mate), overs, records)
    graph, Tuple(edges), Tuple(colours)
end

"""Convert a closed normal presentation. The name `coloured_ograph` belongs to the branching adapter."""
function normal_coloured_ograph(state::OGraphState)
    true_vertex_count(state) >= 1 || throw(ArgumentError(
        "N1 failed: an o-graph has at least one vertex"))
    witness = closedness_witness(state)
    all(witness.criteria) || throw(ArgumentError(
        "C1-C3 failed under $(witness.verifier_version): $(witness.criteria)"))
    graph, edges, colours = _normal_pieces(state)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "converted o-graph failed Definition 2.1.4: " *
        join((check.message for check in failed), "; ")))
    colouring = ClosedNormalColouring(graph, edges, colours, witness, state)
    checks = validate_normal_colouring(colouring)
    failed_normal = filter(check -> !check.passed, checks)
    isempty(failed_normal) || throw(ArgumentError(
        "converted o-graph failed the normal colour checks: " *
        join((check.message for check in failed_normal), "; ")))
    colouring
end

normal_coloured_ograph(closed::ClosedNormalOGraph) = normal_coloured_ograph(closed.state)

function normal_coloured_ograph(link::VirtualLink)
    normal_coloured_ograph(ograph_state(link))
end

# cf:[[BP97#^bp97-sec-0536]] N1 and N2 — opposite-edge directions, checked apart from closedness_witness; unverified

function _check_n1(graph::ColouredOGraph)
    reasons = String[]
    if isempty(graph.vertex_halfedges)
        push!(reasons, "an o-graph has at least one vertex")
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        length(cyclic) == 4 || push!(reasons, "vertex $vertex is not a normal crossing")
        detail = _o2_detail(cyclic, graph.over_positions[vertex])
        isempty(detail) || push!(reasons, "vertex $vertex has $detail")
    end
    isempty(reasons) && return _colour_check("N1", true,
        "at least one vertex, each a marked normal crossing")
    _colour_check("N1", false, join(reasons, "; "))
end

function _above_halfedges(graph::ColouredOGraph)
    above = falses(length(graph.mate))
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        overs = Set(graph.over_positions[vertex])
        for (position, halfedge) in enumerate(cyclic)
            1 <= halfedge <= length(above) || continue
            above[halfedge] = position in overs
        end
    end
    above
end

function _orientation_is_pairing(graph::ColouredOGraph, edges)
    seen = Set{Tuple{Int,Int}}()
    for edge in edges
        tail, head = edge.tail, edge.head
        if tail < 1 || head < 1 || tail > length(graph.mate) || head > length(graph.mate)
            return false
        end
        if graph.mate[tail] != head || graph.mate[head] != tail || tail == head
            return false
        end
        push!(seen, minmax(tail, head))
    end
    length(seen) == length(edges) || return false
    for halfedge in 1:length(graph.mate)
        other = graph.mate[halfedge]
        halfedge < other || continue
        minmax(halfedge, other) in seen || return false
    end
    true
end

function _directions_match(graph::ColouredOGraph, edges)
    _orientation_is_pairing(graph, edges) || return false
    tail_at = falses(length(graph.mate))
    for edge in edges
        tail_at[edge.tail] = true
    end
    for cyclic in graph.vertex_halfedges
        length(cyclic) == 4 || return false
        for pair in ((1, 3), (2, 4))
            count(position -> tail_at[cyclic[position]], pair) == 1 || return false
        end
    end
    true
end

function _check_n2(graph::ColouredOGraph, edges)
    _directions_match(graph, edges) && return _colour_check("N2", true,
        "directions of opposite edges match through the vertices")
    _colour_check("N2", false, "directions of opposite edges do not match through the vertices")
end

function _check_p1(graph::ColouredOGraph, edges)
    _directions_match(graph, edges) && return _colour_check("P1", true,
        "edge orientations match across marked vertices")
    _colour_check("P1", false, "edge orientations do not match across marked vertices")
end

function _check_p2(graph::ColouredOGraph, edges)
    recorded, reasons = _recorded_colours(graph)
    isempty(reasons) || return _colour_check("P2", false, join(reasons, "; "))
    above = _above_halfedges(graph)
    mismatches = String[]
    for edge in edges
        key = minmax(edge.tail, edge.head)
        if !haskey(recorded, key)
            push!(mismatches, "edge $key has no colour")
            continue
        end
        expected = prop321_colour(above[edge.tail], above[edge.head])
        expected == recorded[key] || push!(mismatches, "edge $key colour does not match its ends")
    end
    isempty(mismatches) && return _colour_check("P2", true,
        "every edge colour agrees with its first and second ends")
    _colour_check("P2", false, join(mismatches, "; "))
end

function _check_closed(witness::ClosednessWitness, index::Int, id::String)
    passed = witness.criteria[index]
    passed && return _colour_check(id, true, "closedness witness criterion $(id)")
    _colour_check(id, false, "closedness witness rejected $(id)")
end

"""Return N1, N2, C1, C2, C3, P1, and P2 for one converted colouring.

N2 and P1 use the same opposite-edge matching. C1–C3 are copied from
`witness.criteria` and are not recomputed here.
"""
function validate_normal_colouring(colouring::ClosedNormalColouring)
    graph = colouring.graph
    edges = colouring.edges
    (
        _check_n1(graph),
        _check_n2(graph, edges),
        _check_closed(colouring.witness, 1, "C1"),
        _check_closed(colouring.witness, 2, "C2"),
        _check_closed(colouring.witness, 3, "C3"),
        _check_p1(graph, edges),
        _check_p2(graph, edges),
    )
end

const _NORMAL_ORIENTATION_MAX_EDGES = 12

function _paired_edges(graph::ColouredOGraph)
    edges = Tuple{Int,Int}[]
    for halfedge in 1:length(graph.mate)
        other = graph.mate[halfedge]
        halfedge < other && push!(edges, (halfedge, other))
    end
    edges
end

"""N2-compatible directions of an undirected coloured o-graph.

Each result orients every edge. An undirected graph may admit none, one, or
several. This list is not a count of oriented branchings of the discs.
"""
function normal_orientations(graph::ColouredOGraph)
    checks = validate_coloured_ograph(graph)
    passed = Dict(check.id => check.passed for check in checks)
    (passed["O0"] && passed["O1"] && passed["O2"]) || throw(ArgumentError(
        "normal orientations require O0, O1, and O2"))
    edges = _paired_edges(graph)
    m = length(edges)
    m <= _NORMAL_ORIENTATION_MAX_EDGES || throw(ArgumentError(
        "normal orientation search supports at most $(_NORMAL_ORIENTATION_MAX_EDGES) edges"))
    found = Vector{NTuple{m,NormalEdgeOrientation}}()
    for mask in 0:(1 << m) - 1
        oriented = Vector{NormalEdgeOrientation}(undef, m)
        for index in 1:m
            left, right = edges[index]
            if (mask >> (index - 1)) & 1 == 1
                oriented[index] = NormalEdgeOrientation(left, right)
            else
                oriented[index] = NormalEdgeOrientation(right, left)
            end
        end
        _directions_match(graph, oriented) || continue
        push!(found, Tuple(oriented))
    end
    found
end

"""Orientations from `normal_orientations` whose P2 value equals the stored colour."""
function prop321_orientations(graph::ColouredOGraph)
    filter(edges -> _check_p2(graph, edges).passed, normal_orientations(graph))
end

function coloured_gc_display(colouring::ClosedNormalColouring)
    state = colouring.state
    length(colouring.colours) == length(state.gauss) || throw(ArgumentError(
        "coloured GC display needs one colour per Gauss entry"))
    passages = join(map(label -> label > 0 ? "+$(label)" : string(label), state.gauss), ",")
    signs = join(map(sign -> sign > 0 ? "+" : "−", state.vertex_signs), ",")
    colours = join(map(string, colouring.colours), ",")
    "coloured GC [[[$(passages)]];[$(signs)];[$(colours)]]"
end
