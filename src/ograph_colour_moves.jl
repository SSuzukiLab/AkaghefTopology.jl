module OGraphColourMoves

# cf:[[BP97#^bp97-sec-0549]] Figure 2.4, Definition 2.1.4, Theorem 2.1.5 — C, MP, and the subsegment sum; unverified
# cf:[[BP97#^bp97-sec-0551]] Proposition 2.1.6 — C preserves the oriented spine and MP translates the oriented Matveev–Piergallini move; unverified
# cf:[[BP97#^bp97-sec-0574]] Example 3.4.1 and Figure 3.10 — one MP followed by three C moves; unverified

using Topology

export ColourMoveGraph, CMove, MPForward, MPInverse,
       MOVE_C, MOVE_MP, MOVE_SUB,
       C_ARM_DELTA, MP_DELETED_COLOUR, MP_COLOUR_TM, MP_COLOUR_MD, MP_COLOUR_TD,
       colour_mod3, subsegment_sum,
       from_coloured, to_coloured,
       applicable_c, applicable_mp, c_move, mp_move,
       move_sanity

"""Condition id of the Figure 2.4 C-move."""
const MOVE_C = "C"

"""Condition id of the Figure 2.4 MP-move."""
const MOVE_MP = "MP"

"""Condition id of the subsegment-sum sentence after Theorem 2.1.5."""
const MOVE_SUB = "SUB"

"""Printed arm contribution of one forward C. The inverse subtracts the same integer."""
const C_ARM_DELTA = 1

"""Printed colour of the edge deleted by the positive MP-move."""
const MP_DELETED_COLOUR = 0

"""Printed colour of the new edge from the top vertex to the middle vertex."""
const MP_COLOUR_TM = -1

"""Printed colour of the new edge from the middle vertex to the bottom vertex."""
const MP_COLOUR_MD = 1

"""Printed colour of the new edge from the top vertex to the bottom vertex."""
const MP_COLOUR_TD = 0

"""One Z₃ colour, stored as the same representative as `_colour_mod3`."""
function colour_mod3(value::Integer)
    reduced = mod(Int(value), 3)
    reduced == 2 ? -1 : reduced
end

"""Sum of localized subsegment colours, reduced into {-1, 0, 1}.

Theorem 2.1.5 is followed by the sentence that an edge may carry one colour,
or several localized colours whose sum is that colour. The stored datum is the sum.
An empty list sums to 0. A blank arm in a figure is that omission, not a stored list.
"""
function subsegment_sum(parts)
    parts isa AbstractVector || throw(ArgumentError("SUB: subsegment colours are a list of integers"))
    total = 0
    for part in parts
        part isa Integer || throw(ArgumentError("SUB: a subsegment colour must be an integer"))
        total += Int(part)
    end
    colour_mod3(total)
end

"""Coloured o-graph used by the moves. Half-edges are `1:length(mate)`.

`cyclic[v]` is the counter-clockwise order. Opposite positions are `(1, 3)` and
`(2, 4)`. `over[v]` is one of those pairs, sorted. `colour` keys are sorted
half-edge pairs. One colour is stored per edge.
"""
struct ColourMoveGraph
    cyclic::Vector{NTuple{4,Int}}
    mate::Vector{Int}
    over::Vector{NTuple{2,Int}}
    colour::Dict{Tuple{Int,Int},Int}
    convention_id::String
end

# The generic `==` is identity once a struct holds a vector or a dict.
function Base.:(==)(left::ColourMoveGraph, right::ColourMoveGraph)
    left.cyclic == right.cyclic &&
        left.mate == right.mate &&
        left.over == right.over &&
        left.colour == right.colour &&
        left.convention_id == right.convention_id
end

function Base.hash(graph::ColourMoveGraph, h::UInt)
    h = hash(graph.cyclic, h)
    h = hash(graph.mate, h)
    h = hash(graph.over, h)
    h = hash(graph.colour, h)
    hash(graph.convention_id, h)
end

"""One placement of the Figure 2.4 C-move.

`frame` is the west slot. The next slot is south, the opposite slot is east,
and the remaining slot is north. `delta` is `+1` forward and `-1` inverse.
`east` and `south` are the slots whose half-edges are exchanged.
"""
struct CMove
    id::String
    vertex::Int
    frame::Int
    direction::Symbol
    delta::Int
    east::Int
    south::Int
end

"""Positive Figure 2.4 MP-move: two vertices to three.

`vertex_a` carries the central edge in its east slot and `vertex_b` carries
the same edge in its west slot. Both ends are over. The deleted colour is 0.
"""
struct MPForward
    id::String
    direction::Symbol
    vertex_a::Int
    vertex_b::Int
    frame_a::Int
    frame_b::Int
    deleted_colour::Int
    colour_tm::Int
    colour_md::Int
    colour_td::Int
end

"""Inverse Figure 2.4 MP-move: three vertices to two.

The colour-0 edge runs from the east slot of `vertex_t` to the east slot of
`vertex_d`. Both of those ends are under. The middle vertex is `vertex_m`.
"""
struct MPInverse
    id::String
    direction::Symbol
    vertex_t::Int
    vertex_m::Int
    vertex_d::Int
    frame_t::Int
    frame_m::Int
    frame_d::Int
    deleted_colour::Int
    colour_tm::Int
    colour_md::Int
    colour_td::Int
end

function _pair(left::Int, right::Int)
    left < right ? (left, right) : (right, left)
end

"""Slot `delta` steps counter-clockwise from `frame`. `frame` itself is step 0."""
function _slot(frame::Int, delta::Int)
    mod(frame - 1 + delta, 4) + 1
end

"""Page compass in the stored counter-clockwise order.

The repository's next slot is counter-clockwise. On the page, the
counter-clockwise successor of the west arm is the south arm. Step 1 is south,
step 2 is east, and step 3 is north. The opposite assignment is the reflection.
"""
function _compass(frame::Int)
    (
        west = frame,
        south = _slot(frame, 1),
        east = _slot(frame, 2),
        north = _slot(frame, 3),
    )
end

function _mate_dict(mate::Vector{Int})
    Dict{Int,Int}(index => mate[index] for index in eachindex(mate))
end

function _bump!(colour::Dict{Tuple{Int,Int},Int}, mate::Dict{Int,Int},
                halfedge::Int, delta::Int, id::String)
    delta == 0 && return colour
    haskey(mate, halfedge) || throw(ArgumentError("$id: half-edge $halfedge has no mate"))
    other = mate[halfedge]
    key = minmax(halfedge, other)
    haskey(colour, key) || throw(ArgumentError("$id: edge $key has no colour"))
    colour[key] = colour_mod3(colour[key] + delta)
    colour
end

function _drop_edge!(mate::Dict{Int,Int}, colour::Dict{Tuple{Int,Int},Int},
                     halfedge::Int, id::String)
    haskey(mate, halfedge) || throw(ArgumentError("$id: half-edge $halfedge has no mate"))
    other = mate[halfedge]
    delete!(mate, halfedge)
    delete!(mate, other)
    delete!(colour, minmax(halfedge, other))
    mate
end

function _add_edge!(mate::Dict{Int,Int}, colour::Dict{Tuple{Int,Int},Int},
                    left::Int, right::Int, value::Int)
    mate[left] = right
    mate[right] = left
    colour[_pair(left, right)] = colour_mod3(value)
    mate
end

function _require_coloured(graph::ColourMoveGraph, id::String)
    coloured = to_coloured(graph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(coloured))
    isempty(failed) || throw(ArgumentError(
        "$id: " * join((check.message for check in failed), "; ")))
    graph
end

function _assemble(cyclic::Vector{NTuple{4,Int}}, mate::Dict{Int,Int},
                   over::Vector{NTuple{2,Int}}, colour::Dict{Tuple{Int,Int},Int},
                   convention_id::String, id::String)
    used = Int[]
    for tuple_cyclic in cyclic
        append!(used, tuple_cyclic)
    end
    sort!(used)
    length(unique(used)) == length(used) ||
        throw(ArgumentError("$id: a half-edge is repeated"))
    remap = Dict(old => new for (new, old) in enumerate(used))
    new_cyclic = NTuple{4,Int}[ntuple(index -> remap[tuple_cyclic[index]], 4) for tuple_cyclic in cyclic]
    new_mate = zeros(Int, length(used))
    new_colour = Dict{Tuple{Int,Int},Int}()
    for halfedge in used
        haskey(mate, halfedge) || throw(ArgumentError("$id: half-edge $halfedge has no mate"))
        other = mate[halfedge]
        haskey(remap, other) || throw(ArgumentError("$id: mate $other leaves the graph"))
        new_mate[remap[halfedge]] = remap[other]
        halfedge < other || continue
        key = minmax(halfedge, other)
        haskey(colour, key) || throw(ArgumentError("$id: edge $key has no colour"))
        new_colour[_pair(remap[halfedge], remap[other])] = colour_mod3(colour[key])
    end
    graph = ColourMoveGraph(new_cyclic, new_mate, over, new_colour, convention_id)
    _require_coloured(graph, id)
end

"""Read a validated coloured o-graph. The convention id is kept."""
function from_coloured(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(join((check.message for check in failed), "; ")))
    cyclic = NTuple{4,Int}[]
    over = NTuple{2,Int}[]
    for (vertex, tuple_cyclic) in enumerate(graph.vertex_halfedges)
        length(tuple_cyclic) == 4 ||
            throw(ArgumentError("O0: vertex $vertex is not quadrivalent"))
        push!(cyclic, ntuple(index -> Int(tuple_cyclic[index]), 4))
        positions = sort(collect(Int.(graph.over_positions[vertex])))
        length(positions) == 2 || throw(ArgumentError("O2: vertex $vertex has no over pair"))
        push!(over, (positions[1], positions[2]))
    end
    mate = Int[graph.mate...]
    colour = Dict{Tuple{Int,Int},Int}()
    for record in graph.edge_colours
        halfedge = Int(record.halfedge)
        other = mate[halfedge]
        colour[minmax(halfedge, other)] = colour_mod3(record.colour)
    end
    ColourMoveGraph(cyclic, mate, over, colour, graph.convention_id)
end

"""Write one colour on the smaller half-edge, in the stored convention."""
function to_coloured(graph::ColourMoveGraph)
    records = EdgeColour[]
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < other || continue
        key = minmax(halfedge, other)
        haskey(graph.colour, key) ||
            throw(ArgumentError("O3: edge $key has no colour"))
        push!(records, EdgeColour(halfedge, graph.colour[key]))
    end
    ColouredOGraph(graph.cyclic, graph.mate, graph.over, records;
                   convention_id=graph.convention_id)
end

"""Placements of C. Two forward frames and two inverse frames meet a fixed over pair.

Each frame is a rotation of the printed picture. No frame is a reflection:
every forward witness has delta `+1`.
"""
function applicable_c(graph::ColourMoveGraph)
    moves = CMove[]
    for vertex in eachindex(graph.cyclic)
        over = Set(graph.over[vertex])
        for frame in 1:4
            compass = _compass(frame)
            if over == Set((compass.west, compass.east))
                push!(moves, CMove(MOVE_C, vertex, frame, :forward, C_ARM_DELTA,
                                    compass.east, compass.south))
            elseif over == Set((compass.north, compass.south))
                push!(moves, CMove(MOVE_C, vertex, frame, :inverse, -C_ARM_DELTA,
                                    compass.east, compass.south))
            end
        end
    end
    moves
end

function c_move(graph::ColourMoveGraph, move::CMove)
    move.id == MOVE_C || throw(ArgumentError("C: witness id is $(move.id)"))
    move.direction in (:forward, :inverse) ||
        throw(ArgumentError("C: direction $(move.direction) is not forward or inverse"))
    1 <= move.vertex <= length(graph.cyclic) ||
        throw(ArgumentError("C: vertex $(move.vertex) is outside the graph"))
    1 <= move.frame <= 4 || throw(ArgumentError("C: frame $(move.frame) is outside 1:4"))
    compass = _compass(move.frame)
    (move.east, move.south) == (compass.east, compass.south) ||
        throw(ArgumentError("C: witness slots do not match frame $(move.frame)"))
    cyclic = Int[graph.cyclic[move.vertex]...]
    over = Set(graph.over[move.vertex])
    if move.direction == :forward
        over == Set((compass.west, compass.east)) || throw(ArgumentError(
            "C: forward frame $(move.frame) at vertex $(move.vertex) does not match the over pair"))
        delta = C_ARM_DELTA
        marked = _pair(compass.north, compass.south)
    else
        over == Set((compass.north, compass.south)) || throw(ArgumentError(
            "C: inverse frame $(move.frame) at vertex $(move.vertex) does not match the over pair"))
        delta = -C_ARM_DELTA
        marked = _pair(compass.west, compass.east)
    end
    move.delta == delta || throw(ArgumentError("C: witness delta $(move.delta) is not $delta"))
    colour = copy(graph.colour)
    mate = _mate_dict(graph.mate)
    for slot in 1:4
        _bump!(colour, mate, cyclic[slot], delta, MOVE_C)
    end
    cyclic[compass.east], cyclic[compass.south] = cyclic[compass.south], cyclic[compass.east]
    new_cyclic = copy(graph.cyclic)
    new_cyclic[move.vertex] = ntuple(index -> cyclic[index], 4)
    new_over = copy(graph.over)
    new_over[move.vertex] = marked
    _assemble(new_cyclic, mate, new_over, colour, graph.convention_id, MOVE_C)
end

function _mp_ports(graph::ColourMoveGraph, vertex::Int, frame::Int)
    1 <= vertex <= length(graph.cyclic) ||
        throw(ArgumentError("MP: vertex $vertex is outside the graph"))
    1 <= frame <= 4 || throw(ArgumentError("MP: frame $frame is outside 1:4"))
    compass = _compass(frame)
    cyclic = graph.cyclic[vertex]
    (
        west = cyclic[compass.west],
        south = cyclic[compass.south],
        east = cyclic[compass.east],
        north = cyclic[compass.north],
        compass = compass,
    )
end

function _mp_forward_matches(graph::ColourMoveGraph)
    moves = MPForward[]
    n = length(graph.cyclic)
    for vertex_a in 1:n, vertex_b in 1:n
        vertex_a == vertex_b && continue
        for frame_a in 1:4, frame_b in 1:4
            compass_a = _compass(frame_a)
            compass_b = _compass(frame_b)
            Set(graph.over[vertex_a]) == Set((compass_a.west, compass_a.east)) || continue
            Set(graph.over[vertex_b]) == Set((compass_b.west, compass_b.east)) || continue
            east = graph.cyclic[vertex_a][compass_a.east]
            west = graph.cyclic[vertex_b][compass_b.west]
            graph.mate[east] == west || continue
            colour_mod3(graph.colour[minmax(east, west)]) == MP_DELETED_COLOUR || continue
            push!(moves, MPForward(MOVE_MP, :forward, vertex_a, vertex_b, frame_a, frame_b,
                                    MP_DELETED_COLOUR, MP_COLOUR_TM, MP_COLOUR_MD, MP_COLOUR_TD))
        end
    end
    moves
end

function _frame_south_north(cyclic::NTuple{4,Int}, south_halfedge::Int, north_halfedge::Int)
    south_slot = findfirst(==(south_halfedge), cyclic)
    north_slot = findfirst(==(north_halfedge), cyclic)
    south_slot === nothing && return nothing
    north_slot === nothing && return nothing
    frame = _slot(south_slot, -1)
    _slot(frame, 1) == south_slot || return nothing
    _slot(frame, 3) == north_slot || return nothing
    frame
end

function _mp_inverse_matches(graph::ColourMoveGraph)
    moves = MPInverse[]
    n = length(graph.cyclic)
    for vertex_t in 1:n, vertex_d in 1:n
        vertex_t == vertex_d && continue
        for frame_t in 1:4, frame_d in 1:4
            compass_t = _compass(frame_t)
            compass_d = _compass(frame_d)
            Set(graph.over[vertex_t]) == Set((compass_t.north, compass_t.south)) || continue
            Set(graph.over[vertex_d]) == Set((compass_d.north, compass_d.south)) || continue
            top = graph.cyclic[vertex_t][compass_t.east]
            bottom = graph.cyclic[vertex_d][compass_d.east]
            graph.mate[top] == bottom || continue
            graph.colour[minmax(top, bottom)] == MP_COLOUR_TD || continue
            south_t = graph.cyclic[vertex_t][compass_t.south]
            north_d = graph.cyclic[vertex_d][compass_d.north]
            haskey(graph.colour, minmax(south_t, graph.mate[south_t])) || continue
            graph.colour[minmax(south_t, graph.mate[south_t])] == MP_COLOUR_TM || continue
            graph.colour[minmax(north_d, graph.mate[north_d])] == MP_COLOUR_MD || continue
            middle_north = graph.mate[south_t]
            middle_south = graph.mate[north_d]
            owner = Dict{Int,Int}()
            for (vertex, cyclic) in enumerate(graph.cyclic)
                for halfedge in cyclic
                    owner[halfedge] = vertex
                end
            end
            vertex_m = get(owner, middle_north, 0)
            vertex_m == 0 && continue
            (vertex_m == vertex_t || vertex_m == vertex_d) && continue
            owner[middle_south] == vertex_m || continue
            frame_m = _frame_south_north(graph.cyclic[vertex_m], middle_south, middle_north)
            frame_m === nothing && continue
            compass_m = _compass(frame_m)
            Set(graph.over[vertex_m]) == Set((compass_m.west, compass_m.east)) || continue
            push!(moves, MPInverse(MOVE_MP, :inverse, vertex_t, vertex_m, vertex_d,
                                    frame_t, frame_m, frame_d,
                                    MP_DELETED_COLOUR, MP_COLOUR_TM, MP_COLOUR_MD, MP_COLOUR_TD))
        end
    end
    moves
end

"""Forward and inverse placements of MP. Reflections are not generated."""
function applicable_mp(graph::ColourMoveGraph)
    moves = Union{MPForward,MPInverse}[]
    append!(moves, _mp_forward_matches(graph))
    append!(moves, _mp_inverse_matches(graph))
    moves
end

function _kept(graph::ColourMoveGraph, dropped::Set{Int})
    cyclic = NTuple{4,Int}[]
    over = NTuple{2,Int}[]
    for vertex in eachindex(graph.cyclic)
        vertex in dropped && continue
        push!(cyclic, graph.cyclic[vertex])
        push!(over, graph.over[vertex])
    end
    cyclic, over
end

function mp_move(graph::ColourMoveGraph, move::MPForward)
    move.id == MOVE_MP || throw(ArgumentError("MP: witness id is $(move.id)"))
    move.direction == :forward || throw(ArgumentError("MP: forward witness has the wrong direction"))
    (move.deleted_colour, move.colour_tm, move.colour_md, move.colour_td) ==
        (MP_DELETED_COLOUR, MP_COLOUR_TM, MP_COLOUR_MD, MP_COLOUR_TD) ||
        throw(ArgumentError("MP: witness colours are not the printed MP colours"))
    move.vertex_a == move.vertex_b &&
        throw(ArgumentError("MP: the central edge has two distinct ends"))
    ports_a = _mp_ports(graph, move.vertex_a, move.frame_a)
    ports_b = _mp_ports(graph, move.vertex_b, move.frame_b)
    Set(graph.over[move.vertex_a]) == Set((ports_a.compass.west, ports_a.compass.east)) ||
        throw(ArgumentError("MP: vertex $(move.vertex_a) is not horizontal over in frame $(move.frame_a)"))
    Set(graph.over[move.vertex_b]) == Set((ports_b.compass.west, ports_b.compass.east)) ||
        throw(ArgumentError("MP: vertex $(move.vertex_b) is not horizontal over in frame $(move.frame_b)"))
    graph.mate[ports_a.east] == ports_b.west ||
        throw(ArgumentError("MP: the east slot of A is not mated to the west slot of B"))
    graph.colour[minmax(ports_a.east, ports_b.west)] == MP_DELETED_COLOUR ||
        throw(ArgumentError("MP: the deleted edge is not colour 0"))
    mate = _mate_dict(graph.mate)
    colour = copy(graph.colour)
    _bump!(colour, mate, ports_b.north, 1, MOVE_MP)
    _bump!(colour, mate, ports_b.south, -1, MOVE_MP)
    _drop_edge!(mate, colour, ports_a.east, MOVE_MP)
    base = length(graph.mate)
    tm_t, tm_m = base + 1, base + 2
    md_m, md_d = base + 3, base + 4
    td_t, td_d = base + 5, base + 6
    _add_edge!(mate, colour, tm_t, tm_m, MP_COLOUR_TM)
    _add_edge!(mate, colour, md_m, md_d, MP_COLOUR_MD)
    _add_edge!(mate, colour, td_t, td_d, MP_COLOUR_TD)
    cyclic, over = _kept(graph, Set((move.vertex_a, move.vertex_b)))
    # Counter-clockwise from west: west, south, east, north.
    push!(cyclic, (ports_a.north, tm_t, td_t, ports_b.north))
    push!(over, (2, 4))
    push!(cyclic, (ports_a.west, md_m, ports_b.east, tm_m))
    push!(over, (1, 3))
    push!(cyclic, (ports_a.south, ports_b.south, td_d, md_d))
    push!(over, (2, 4))
    _assemble(cyclic, mate, over, colour, graph.convention_id, MOVE_MP)
end

function mp_move(graph::ColourMoveGraph, move::MPInverse)
    move.id == MOVE_MP || throw(ArgumentError("MP: witness id is $(move.id)"))
    move.direction == :inverse || throw(ArgumentError("MP: inverse witness has the wrong direction"))
    (move.deleted_colour, move.colour_tm, move.colour_md, move.colour_td) ==
        (MP_DELETED_COLOUR, MP_COLOUR_TM, MP_COLOUR_MD, MP_COLOUR_TD) ||
        throw(ArgumentError("MP: witness colours are not the printed MP colours"))
    length(Set((move.vertex_t, move.vertex_m, move.vertex_d))) == 3 ||
        throw(ArgumentError("MP: the inverse uses three distinct vertices"))
    ports_t = _mp_ports(graph, move.vertex_t, move.frame_t)
    ports_m = _mp_ports(graph, move.vertex_m, move.frame_m)
    ports_d = _mp_ports(graph, move.vertex_d, move.frame_d)
    Set(graph.over[move.vertex_t]) == Set((ports_t.compass.north, ports_t.compass.south)) ||
        throw(ArgumentError("MP: the top vertex is not vertical over"))
    Set(graph.over[move.vertex_d]) == Set((ports_d.compass.north, ports_d.compass.south)) ||
        throw(ArgumentError("MP: the bottom vertex is not vertical over"))
    Set(graph.over[move.vertex_m]) == Set((ports_m.compass.west, ports_m.compass.east)) ||
        throw(ArgumentError("MP: the middle vertex is not horizontal over"))
    graph.mate[ports_t.east] == ports_d.east ||
        throw(ArgumentError("MP: the top and bottom east slots are not one edge"))
    graph.colour[minmax(ports_t.east, ports_d.east)] == MP_COLOUR_TD ||
        throw(ArgumentError("MP: the top-bottom edge is not colour 0"))
    graph.mate[ports_t.south] == ports_m.north ||
        throw(ArgumentError("MP: the top south slot is not mated to the middle north slot"))
    graph.colour[minmax(ports_t.south, ports_m.north)] == MP_COLOUR_TM ||
        throw(ArgumentError("MP: the top-middle edge is not colour -1"))
    graph.mate[ports_d.north] == ports_m.south ||
        throw(ArgumentError("MP: the bottom north slot is not mated to the middle south slot"))
    graph.colour[minmax(ports_d.north, ports_m.south)] == MP_COLOUR_MD ||
        throw(ArgumentError("MP: the middle-bottom edge is not colour 1"))
    mate = _mate_dict(graph.mate)
    colour = copy(graph.colour)
    _bump!(colour, mate, ports_t.north, -1, MOVE_MP)
    _bump!(colour, mate, ports_d.south, 1, MOVE_MP)
    _drop_edge!(mate, colour, ports_t.south, MOVE_MP)
    _drop_edge!(mate, colour, ports_d.north, MOVE_MP)
    _drop_edge!(mate, colour, ports_t.east, MOVE_MP)
    base = length(graph.mate)
    central_a, central_b = base + 1, base + 2
    _add_edge!(mate, colour, central_a, central_b, MP_DELETED_COLOUR)
    cyclic, over = _kept(graph, Set((move.vertex_t, move.vertex_m, move.vertex_d)))
    push!(cyclic, (ports_m.west, ports_d.west, central_a, ports_t.west))
    push!(over, (1, 3))
    push!(cyclic, (central_b, ports_d.south, ports_m.east, ports_t.north))
    push!(over, (1, 3))
    _assemble(cyclic, mate, over, colour, graph.convention_id, MOVE_MP)
end

"""Counts named in the 2026-11-22 sanity rule, under calibration `:b` unless another is passed.

Face cycles, homology, and the boundary use the draft candidate B. A failing
count throws. The throw is a fault in the object being checked.
"""
function move_sanity(graph::ColouredOGraph; calibration::Symbol=:b)
    poly = face_cycles(graph; calibration)
    boundary = thickening_boundary(graph; calibration)
    homology = cellular_homology(graph; calibration)
    pairings = dual_face_pairings(graph; calibration)
    orbits = dual_edge_orbits(graph; calibration)
    n = poly.vertex_count
    problems = String[]
    poly.edge_count == 2 * n || push!(problems, "edges $(poly.edge_count)")
    poly.halfedge_count == 4 * n || push!(problems, "half-edges $(poly.halfedge_count)")
    poly.germ_count == 6 * n || push!(problems, "germs $(poly.germ_count)")
    poly.boundary_length == 6 * n || push!(problems, "boundary length $(poly.boundary_length)")
    poly.euler == poly.disc_count - n || push!(problems, "euler $(poly.euler)")
    boundary.euler_sum == 2 * poly.euler ||
        push!(problems, "boundary euler $(boundary.euler_sum) against $(poly.euler)")
    for (index, component) in enumerate(boundary.components)
        component.closed_surface || push!(problems, "component $index is not a closed surface")
        component.orientable || push!(problems, "component $index is not orientable")
        component.euler <= 2 || push!(problems, "component $index has euler $(component.euler)")
    end
    length(pairings) == 2 * n || push!(problems, "face pairs $(length(pairings))")
    all(pairing -> pairing.orientation_reversing, pairings) ||
        push!(problems, "a dual face pair does not reverse orientation")
    length(orbits) == poly.disc_count ||
        push!(problems, "dual edges $(length(orbits)) against discs $(poly.disc_count)")
    report = (
        n = n,
        edges = poly.edge_count,
        halfedges = poly.halfedge_count,
        germs = poly.germ_count,
        discs = poly.disc_count,
        euler = poly.euler,
        boundary_length = poly.boundary_length,
        component_euler = [component.euler for component in boundary.components],
        boundary_euler = boundary.euler_sum,
        h1_free = homology.free_rank,
        h1_torsion = homology.torsion_factors,
        dual_edges = length(orbits),
        dual_pairs = length(pairings),
        calibration = calibration,
        status = OGRAPH_FACE_CYCLE_STATUS,
    )
    isempty(problems) || throw(ArgumentError(
        "sanity failed: " * join(problems, "; ") * " :: " * string(report)))
    report
end

end
