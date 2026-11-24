# cf:[[BP97#^bp97-sec-0549]] Definition 2.1.4 and Fig. 2.5 — O0–O3 and the edge chart (1 2 3)^α; adjustment: upright legs (upper-right, up, lower-left), upside-down legs (lower-left, down, upper-right); the colour is the undirected Z_3 label and is not negated when the drawing is reversed; unverified

"""Coloured o-graphs in the sense of Definition 2.1.4.

The general object is a half-edge multigraph. Loops and parallel edges stay.
It has no vertex sign, no edge direction, and no Gauss word.

Each vertex stores its incident half-edges in cyclic order. On a quadrivalent
vertex the opposite positions are `(1, 3)` and `(2, 4)`. `over_positions[v]`
lists the local positions marked over. Exactly one of those two pairs must be
marked.

Each edge carries one colour, an element of Z_3, stored on either of its
half-edges. Integers are read modulo 3, with the representative in
`{-1, 0, 1}`: `2`, `-2`, and `4` are colours, not malformed data. The printed
Fig. 2.5 rule, once an end has been drawn on the right, glues right leg `i` to
left leg `τ^α(i)`. Reversing the drawing reverses the transverse frame by
`ρ = (1 3)` and keeps `α`. Freezing the leg numbers and only reversing the
bijection produces `-α`; that is a different calculation, and it is not how a
colour is stored. Status of this transport: draft.

`bp97-fig2.5-rl-tau/draft1` names the chart. Figure 2.5 (left) specifies the
local polyhedron. The over pair is the upright T and the under pair is the
upside-down T, and leg 2 meets the opposite half-edge. Figure 2.7 fixes the
local embedding once disc D has the usual orientation of the horizontal plane
and the screw is left-handed. Which side sheet of that embedding is leg 1 is
not registered onto the counter-clockwise slot. Candidate B is the calibration
derived from Proposition 3.2.1 together with the §2.1 reconstruction and the
§3.2 construction, under the draft dependencies on the convention page.
Candidate A on the closed presentation `s2xs1_2vertex` gives five discs and
Euler characteristic 3. `face_cycles` builds the discs. The default
calibration is `:b`. `OGRAPH_FACE_CYCLE_STATUS` records that derivation.
The status stays draft.

The reading check is `threads/PJ261120/sources/pro/core-20261121-235637.md`,
harvested at `threads/PJ261120/sources/pro/core-20261121-235637-answer.md`.
The colour sits on the undirected edge, as Definition 2.1.4 states. The frame
reversal `ρ = (1 3)` and the frozen-label inverse `-α` are inferences from
Figure 2.5. The twelve pairings are derived from the four panels of Figure 2.5
(right) and the gluing rule `τ^α`. The PDF does not print those twelve rows.
The condition-id check is `threads/PJ261120/sources/pro/core-20261122-012651.md`.
The slot-transcription check is `threads/PJ261120/sources/pro/core-20261122-033946.md`.
It leaves both charts of leg 1 open. The later calibration check is
`threads/PJ261120/sources/pro/core-20261122-045101.md`. Candidate B is the
calibration derived from Proposition 3.2.1, and the status stays draft.
"""

using SHA

const OGRAPH_COLOUR_CONVENTION_ID = "bp97-fig2.5-rl-tau/draft1"
const OGRAPH_FACE_CYCLE_STATUS = "prop321-calibration-b/draft"
const _COLOUR_CANON_MAX_VERTICES = 6

struct DiscGerm
    halfedge::Int
    leg::Int

    function DiscGerm(halfedge::Integer, leg::Integer)
        halfedge > 0 || throw(ArgumentError("disc-germ half-edge must be positive"))
        leg in 1:3 || throw(ArgumentError("disc-germ leg must be 1, 2, or 3"))
        new(Int(halfedge), Int(leg))
    end
end

struct EdgeColour
    halfedge::Int
    colour::Int
end

struct OGraphColourCheck
    id::String
    passed::Bool
    message::String
end

struct ColouredOGraph
    vertex_halfedges::Tuple{Vararg{Tuple{Vararg{Int}}}}
    mate::Tuple{Vararg{Int}}
    over_positions::Tuple{Vararg{Tuple{Vararg{Int}}}}
    edge_colours::Tuple{Vararg{EdgeColour}}
    convention_id::String

    function ColouredOGraph(vertex_halfedges, mate, over_positions, edge_colours;
                            convention_id::AbstractString=OGRAPH_COLOUR_CONVENTION_ID)
        convention = strip(String(convention_id))
        isempty(convention) && throw(ArgumentError("convention_id must not be empty"))
        vertices = Tuple(Tuple(Int.(collect(cyclic))) for cyclic in vertex_halfedges)
        partners = Tuple(Int.(collect(mate)))
        overs = Tuple(Tuple(Int.(collect(positions))) for positions in over_positions)
        stored = Tuple(_edge_colour(record) for record in edge_colours)
        length(overs) == length(vertices) || throw(ArgumentError(
            "over_positions must contain one entry per vertex"))
        new(vertices, partners, overs, stored, convention)
    end
end

_edge_colour(record::EdgeColour) = record

function _edge_colour(record::Tuple{Integer,Integer})
    EdgeColour(Int(record[1]), Int(record[2]))
end

function _edge_colour(record)
    throw(ArgumentError("an edge colour must be an EdgeColour or a (halfedge, colour) pair"))
end

function _colour_mod3(value::Integer)
    reduced = mod(Int(value), 3)
    reduced == 2 ? -1 : reduced
end

"""Return τ^α(leg), where τ = (1 2 3) and τ(1) = 2, τ(2) = 3, τ(3) = 1."""
function fig25_tau(alpha::Integer, leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    shift = mod(Int(alpha), 3)
    shift == 0 && return Int(leg)
    shift == 1 && return leg == 3 ? 1 : Int(leg) + 1
    return leg == 1 ? 3 : Int(leg) - 1
end

"""Direction of one numbered leg. `true` is right-way-up; `false` is upside-down.

The components are `(1, 1)` upper-right, `(0, 1)` up, `(-1, -1)` lower-left,
and `(0, -1)` down. This is the locked Fig. 2.5 (right) chart, still draft.
"""
function fig25_leg_direction(upright::Bool, leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    upright && return ((1, 1), (0, 1), (-1, -1))[leg]
    return ((-1, -1), (0, -1), (1, 1))[leg]
end

"""Three physical leg pairs for one drawing, ordered by the right leg 1, 2, 3.

Each pair is `(direction of the right leg, direction of the left leg it meets)`.
The left leg number is τ^α of the right leg number. `left_up` and `right_up`
say whether those T's are right-way-up in this drawing.
"""
function fig25_leg_pairing(left_up::Bool, right_up::Bool, alpha::Integer)
    colour = _colour_mod3(alpha)
    ntuple(leg -> (fig25_leg_direction(right_up, leg),
                   fig25_leg_direction(left_up, fig25_tau(colour, leg))), 3)
end

"""Left leg numbers glued to right legs 1, 2, 3 in one drawing of colour `alpha`."""
function edge_gluing_from_colour(alpha::Integer)
    colour = _colour_mod3(alpha)
    ntuple(leg -> fig25_tau(colour, leg), 3)
end

"""Recover the representative in `{-1, 0, 1}` from left legs glued to right legs 1, 2, 3."""
function colour_from_gluing(left_legs::NTuple{3,Integer})
    images = ntuple(leg -> Int(left_legs[leg]), 3)
    for alpha in (-1, 0, 1)
        edge_gluing_from_colour(alpha) == images && return alpha
    end
    throw(ArgumentError("gluing is not a power of (1 2 3) under the Fig. 2.5 chart"))
end

"""Colour read after numbering shifts τ^u on the right and τ^v on the left.

The new label of the physical leg whose old label was `i` is τ^u(i) on the right
and τ^v(i) on the left. The resulting colour is α + v − u, reduced into {-1, 0, 1}.
The stored colour of an o-graph uses the unshifted chart.
"""
function fig25_gauged_colour(alpha::Integer, right_shift::Integer, left_shift::Integer)
    _colour_mod3(Int(alpha) + Int(left_shift) - Int(right_shift))
end

"""Transverse frame reversal ρ = (1 3). Leg 2, the stem, stays put."""
function fig25_frame_reversal(leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    leg == 2 ? 2 : 4 - Int(leg)
end

"""Disc-germ pairs for one drawing that places `right_halfedge` on the right.

Legs use that drawing's printed numbers. The opposite drawing reverses both
frames by `fig25_frame_reversal` and keeps the same colour.
"""
function edge_disc_germs(colour::Integer, right_halfedge::Integer, mate::Tuple{Vararg{Int}})
    1 <= right_halfedge <= length(mate) || throw(ArgumentError(
        "drawing right end is outside the half-edge set"))
    left = mate[right_halfedge]
    alpha = _colour_mod3(colour)
    ntuple(leg -> (DiscGerm(right_halfedge, leg),
                   DiscGerm(left, fig25_tau(alpha, leg))), 3)
end

function _colour_check(id::AbstractString, passed::Bool, detail::AbstractString)
    message = passed ? "$id passed: $detail" : "$id failed: $detail"
    occursin(id, message) || error("colour check message lost its condition id")
    OGraphColourCheck(String(id), passed, message)
end

function _listed_halfedges(graph::ColouredOGraph)
    listed = Int[]
    for cyclic in graph.vertex_halfedges
        append!(listed, cyclic)
    end
    listed
end

function _colour_edges(graph::ColouredOGraph)
    edges = Tuple{Int,Int}[]
    reasons = String[]
    n = length(graph.mate)
    for (halfedge, other) in enumerate(graph.mate)
        if other < 1 || other > n || graph.mate[other] != halfedge || other == halfedge
            push!(reasons, "half-edge $halfedge is not paired")
            continue
        end
        halfedge < other && push!(edges, (halfedge, other))
    end
    edges, reasons
end

function _colour_components(graph::ColouredOGraph, edges::Vector{Tuple{Int,Int}})
    n = length(graph.vertex_halfedges)
    parent = collect(1:n)
    function find(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        vertex
    end
    owner = Dict{Int,Int}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for halfedge in cyclic
            owner[halfedge] = vertex
        end
    end
    for (left, right) in edges
        haskey(owner, left) || continue
        haskey(owner, right) || continue
        a = find(owner[left])
        b = find(owner[right])
        a == b || (parent[b] = a)
    end
    n == 0 ? 0 : length(unique(find(vertex) for vertex in 1:n))
end

function _check_o0(graph::ColouredOGraph)
    reasons = String[]
    listed = _listed_halfedges(graph)
    n = length(graph.mate)
    if sort(listed) != collect(1:n) || !allunique(listed)
        push!(reasons, "half-edges are not exactly the set 1:$n")
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        length(cyclic) == 4 || push!(reasons, "vertex $vertex is not quadrivalent")
    end
    edges, pair_reasons = _colour_edges(graph)
    append!(reasons, unique(pair_reasons))
    if isempty(reasons) && _colour_components(graph, edges) != 1
        push!(reasons, "graph is not connected")
    end
    isempty(reasons) && return _colour_check("O0", true, "connected finite quadrivalent graph")
    _colour_check("O0", false, join(reasons, "; "))
end

function _check_o1(graph::ColouredOGraph)
    if length(graph.vertex_halfedges) >= 1
        return _colour_check("O1", true, "at least one vertex")
    end
    _colour_check("O1", false, "an o-graph has at least one vertex")
end

function _o2_detail(cyclic::Tuple, positions::Tuple)
    length(cyclic) == 4 || return "missing an over pair"
    isempty(positions) && return "missing an over pair"
    any(position -> position < 1 || position > length(cyclic), positions) &&
        return "missing an over pair"
    marked = Set(positions)
    (marked == Set((1, 3)) || marked == Set((2, 4))) && return ""
    marked == Set((1, 2, 3, 4)) && return "more than one over pair"
    length(marked) == 2 && return "over branches are not opposite"
    return "missing an over pair"
end

function _check_o2(graph::ColouredOGraph)
    reasons = String[]
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        detail = _o2_detail(cyclic, graph.over_positions[vertex])
        isempty(detail) || push!(reasons, "vertex $vertex has $detail")
    end
    isempty(reasons) && return _colour_check("O2", true,
        "every vertex has one marked opposite over pair")
    _colour_check("O2", false, join(reasons, "; "))
end

function _recorded_colours(graph::ColouredOGraph)
    recorded = Dict{Tuple{Int,Int},Int}()
    reasons = String[]
    for record in graph.edge_colours
        right = record.halfedge
        if right < 1 || right > length(graph.mate)
            push!(reasons, "colour record does not select an edge")
            continue
        end
        left = graph.mate[right]
        if left < 1 || left > length(graph.mate) || graph.mate[left] != right || left == right
            push!(reasons, "colour record does not select an edge")
            continue
        end
        key = minmax(left, right)
        reduced = _colour_mod3(record.colour)
        if haskey(recorded, key) && recorded[key] != reduced
            push!(reasons, "edge $key has more than one colour")
        else
            recorded[key] = reduced
        end
    end
    recorded, reasons
end

function _check_o3(graph::ColouredOGraph)
    edges, _ = _colour_edges(graph)
    recorded, reasons = _recorded_colours(graph)
    for edge in edges
        haskey(recorded, edge) || push!(reasons, "edge $edge has no colour")
    end
    isempty(reasons) && return _colour_check("O3", true, "every edge has one colour in Z_3")
    _colour_check("O3", false, join(reasons, "; "))
end

"""Return the O0, O1, O2, and O3 checks, in that order.

O0 is the preamble of Definition 2.1.4: a connected finite quadrivalent graph.
O1, O2, and O3 are items 1, 2, and 3. A colour may be written as any integer;
it is read in Z_3. Two records on one edge that reduce to the same value are
one colour. Two records that reduce to different values fail O3. A failing
`message` contains its condition id.
"""
function validate_coloured_ograph(graph::ColouredOGraph)
    (_check_o0(graph), _check_o1(graph), _check_o2(graph), _check_o3(graph))
end

function _require_coloured_ograph(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "canonical form requires O0-O3: " * join((check.message for check in failed), "; ")))
    length(graph.vertex_halfedges) <= _COLOUR_CANON_MAX_VERTICES || throw(ArgumentError(
        "coloured o-graph canonical search supports at most $(_COLOUR_CANON_MAX_VERTICES) vertices"))
    graph
end

function _rotated_over_code(positions::Tuple, rotation::Integer)
    moved = Set(mod1(position - rotation, 4) for position in positions)
    moved == Set((1, 3)) && return 1
    moved == Set((2, 4)) && return 2
    throw(ArgumentError("over pair is not one opposite pair"))
end

function _permutations(n::Integer)
    results = Vector{Int}[]
    used = falses(n)
    current = Int[]
    function walk()
        if length(current) == n
            push!(results, copy(current))
            return
        end
        for vertex in 1:n
            used[vertex] && continue
            used[vertex] = true
            push!(current, vertex)
            walk()
            pop!(current)
            used[vertex] = false
        end
    end
    walk()
    results
end

function _rotations(n::Integer)
    results = Vector{Int}[]
    current = zeros(Int, n)
    function walk(index)
        if index > n
            push!(results, copy(current))
            return
        end
        for rotation in 0:3
            current[index] = rotation
            walk(index + 1)
        end
    end
    walk(1)
    results
end

function _canonical_candidate(graph::ColouredOGraph, permutation::Vector{Int},
                              rotations::Vector{Int}, colour_of_edge::Dict{Tuple{Int,Int},Int})
    n = length(permutation)
    new_of_old = Dict{Int,Int}()
    new_over = Int[]
    next_id = 1
    for new_vertex in 1:n
        old_vertex = permutation[new_vertex]
        cyclic = graph.vertex_halfedges[old_vertex]
        rotation = rotations[new_vertex]
        rotated = (cyclic[rotation + 1:end]..., cyclic[1:rotation]...)
        for halfedge in rotated
            new_of_old[halfedge] = next_id
            next_id += 1
        end
        push!(new_over, _rotated_over_code(graph.over_positions[old_vertex], rotation))
    end
    mate = zeros(Int, 4n)
    colour_by_end = Dict{Int,Int}()
    for (left, right) in keys(colour_of_edge)
        new_left = new_of_old[left]
        new_right = new_of_old[right]
        mate[new_left] = new_right
        mate[new_right] = new_left
        smaller = min(new_left, new_right)
        colour_by_end[smaller] = colour_of_edge[(left, right)]
    end
    ordered_colours = Int[]
    for halfedge in 1:4n
        halfedge < mate[halfedge] || continue
        push!(ordered_colours, colour_by_end[halfedge])
    end
    (graph.convention_id, Tuple(new_over), Tuple(mate), Tuple(ordered_colours))
end

"""Relabel vertices and cyclic origins. Do not quotient by C or MP.

The representative stores vertices `1:n` and half-edges `4(v - 1) + 1:4v`.
Each colour is the representative in `{-1, 0, 1}` and is written on the
smaller half-edge id. Exchanging the ends does not change the colour.
Graphs above the search bound are rejected rather than given a partial
representative. Cyclic order may be rotated. It is not reversed.
"""
function canonical_coloured_ograph(graph::ColouredOGraph)
    _require_coloured_ograph(graph)
    n = length(graph.vertex_halfedges)
    colour_of_edge, _ = _recorded_colours(graph)
    best = nothing
    for permutation in _permutations(n)
        for rotations in _rotations(n)
            candidate = _canonical_candidate(graph, permutation, rotations, colour_of_edge)
            if best === nothing || candidate < best
                best = candidate
            end
        end
    end
    convention, over_codes, mate, colours = best
    vertices = Tuple(Tuple((4 * (vertex - 1) + 1):(4 * vertex)) for vertex in 1:n)
    overs = Tuple(code == 1 ? (1, 3) : (2, 4) for code in over_codes)
    records = EdgeColour[]
    colour_index = 1
    for halfedge in 1:length(mate)
        halfedge < mate[halfedge] || continue
        push!(records, EdgeColour(halfedge, colours[colour_index]))
        colour_index += 1
    end
    ColouredOGraph(vertices, mate, overs, records; convention_id=convention)
end

function coloured_ograph_serialization(graph::ColouredOGraph)
    canonical = canonical_coloured_ograph(graph)
    colours = join((string(record.halfedge, ":", record.colour) for record in canonical.edge_colours), ",")
    mates = join(canonical.mate, ",")
    overs = join((join(positions, "+") for positions in canonical.over_positions), ",")
    "ColouredOGraph;convention=$(canonical.convention_id);n=$(length(canonical.vertex_halfedges));" *
        "over_positions=$overs;mate=$mates;colours=$colours"
end

function coloured_ograph_hash(graph::ColouredOGraph)
    bytes2hex(SHA.sha256(codeunits(coloured_ograph_serialization(graph))))
end
