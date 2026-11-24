# cf:[[BP97#^bp97-sec-0549]] Fig. 2.5 — discs pasted after the edge gluings; adjustment: leg 1 is candidate B, derived from Proposition 3.2.1, draft
# cf:[[BP97#^bp97-sec-0566]] Proposition 3.2.1 — the P2 colouring represents the same spine; candidate A is refuted by the closed-spine count

import AbstractAlgebra

"""One face cycle of the standard polyhedron.

`germs[k]` is the sheet end where step `k` leaves an edge. `corners[k]` is the
vertex disc germ `(vertex, slot_lo, slot_hi)` where that step arrives.
`edge_letters[k]` is the signed index of the singular edge. A positive letter
leaves the smaller half-edge id. That sign is a chain orientation, not a
branching.
"""
struct FaceCycle
    germs::Vector{DiscGerm}
    corners::Vector{NTuple{3,Int}}
    edge_letters::Vector{Int}
end

"""Cell counts of the standard polyhedron built from a coloured o-graph.

`calibration` is `:b`, the chart derived from Proposition 3.2.1, or `:a`, the
chart that proposition refutes. The default is `:b`.
"""
struct ColouredPolyhedron
    vertex_count::Int
    edge_count::Int
    disc_count::Int
    euler::Int
    halfedge_count::Int
    germ_count::Int
    boundary_length::Int
    cycles::Vector{FaceCycle}
    calibration::Symbol
    status::String
end

struct ColourHomology
    free_rank::Int
    torsion_factors::Vector{BigInt}
end

struct DualFacePairing
    source_vertex::Int
    source_slot::Int
    target_vertex::Int
    target_slot::Int
    orientation_reversing::Bool
end

"""One boundary component of a regular neighbourhood of the polyhedron.

The complex is one connected piece of the links of the tetrahedron corners.
Each corner contributes a triangle, and each dual face gluing pairs a side.
`euler` is `vertices - edges + faces` of that piece. `closed_surface` means
every side is paired. `orientable` means every face gluing in the piece
reverses orientation.
"""
struct SpineBoundaryComponent
    euler::Int
    face_count::Int
    edge_count::Int
    vertex_count::Int
    orientable::Bool
    closed_surface::Bool
end

struct SpineBoundary
    components::Vector{SpineBoundaryComponent}
    euler_sum::Int
end

"""Germ-by-germ comparison of one calibration with `disk_table`.

`ours` and `pipeline` are partitions: each disc is its corners
`(vertex, slot_lo, slot_hi)`, sorted, and the discs are sorted. Equality
includes repeated corners. The pipeline corners come from
`disk_table_corners`.
"""
struct DiscPartitionComparison
    equal::Bool
    ours::Vector{Vector{NTuple{3,Int}}}
    pipeline::Vector{Vector{NTuple{3,Int}}}
    germ_count::Int
    unique_germs::Int
end

function _face_owner(graph::ColouredOGraph)
    owner = Dict{Int,Tuple{Int,Int}}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for (slot, halfedge) in enumerate(cyclic)
            owner[Int(halfedge)] = (vertex, slot)
        end
    end
    owner
end

function _face_colours(graph::ColouredOGraph)
    colour = Dict{Tuple{Int,Int},Int}()
    for record in graph.edge_colours
        other = Int(graph.mate[record.halfedge])
        colour[minmax(Int(record.halfedge), other)] = _colour_mod3(record.colour)
    end
    colour
end

function _face_edges(graph::ColouredOGraph)
    ids = zeros(Int, length(graph.mate))
    edges = Tuple{Int,Int}[]
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        push!(edges, (halfedge, Int(other)))
        ids[halfedge] = length(edges)
        ids[Int(other)] = length(edges)
    end
    edges, ids
end

"""Leg chart. `:a` puts leg 1 on the counter-clockwise side of an over slot.
`:b` swaps legs 1 and 3. Leg 2 stays the opposite half-edge.
"""
function _leg_chart(graph::ColouredOGraph, calibration::Symbol)
    calibration in (:a, :b) || throw(ArgumentError(
        "leg calibration must be :a or :b"))
    charts = Vector{NTuple{3,Int}}(undef, length(graph.mate))
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        over = graph.over_positions[vertex]
        for (slot, halfedge) in enumerate(cyclic)
            opposite = Int(cyclic[mod1(slot + 2, 4)])
            ccw = Int(cyclic[mod1(slot + 1, 4)])
            cw = Int(cyclic[mod1(slot - 1, 4)])
            upright = slot in over
            adjacent = upright ? (ccw, cw) : (cw, ccw)
            calibration == :b && (adjacent = (adjacent[2], adjacent[1]))
            charts[Int(halfedge)] = (adjacent[1], opposite, adjacent[2])
        end
    end
    charts
end

"Transport of one sheet across an edge: ρ∘τ^α, with the stored colour."
function _sheet_transport(colour::Integer, leg::Integer)
    fig25_frame_reversal(fig25_tau(Int(colour), leg))
end

function _corner(owner, halfedge::Int, partner::Int)
    vertex, slot_a = owner[halfedge]
    other_vertex, slot_b = owner[partner]
    other_vertex == vertex || throw(ArgumentError(
        "a vertex corner joins half-edges at two vertices"))
    lo, hi = minmax(slot_a, slot_b)
    (vertex, lo, hi)
end

"""Partition the 6n sheet germs into face cycles.

The default calibration is B. Pass `calibration=:a` only to reproduce the
refuted chart.
"""
function face_cycles(graph::ColouredOGraph; calibration::Symbol=:b)
    _require_coloured_ograph_shape(graph)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    _, edge_id = _face_edges(graph)
    n = length(graph.mate)
    edge_count = count(i -> i < graph.mate[i], 1:n)
    used = falses(n, 3)
    cycles = FaceCycle[]
    for start_halfedge in 1:n, start_leg in 1:3
        used[start_halfedge, start_leg] && continue
        germs = DiscGerm[]
        corners = NTuple{3,Int}[]
        letters = Int[]
        halfedge, leg = start_halfedge, start_leg
        while true
            used[halfedge, leg] && throw(ArgumentError(
                "sheet end ($halfedge, $leg) is reused"))
            used[halfedge, leg] = true
            other = Int(graph.mate[halfedge])
            alpha = get(colours, minmax(halfedge, other), nothing)
            alpha === nothing && throw(ArgumentError(
                "edge ($halfedge, $other) has no colour"))
            other_leg = _sheet_transport(alpha, leg)
            used[other, other_leg] && throw(ArgumentError(
                "arrival sheet end ($other, $other_leg) is reused"))
            used[other, other_leg] = true
            partner = charts[other][other_leg]
            back = findfirst(==(other), charts[partner])
            back === nothing && throw(ArgumentError(
                "half-edge $partner has no leg returning to $other"))
            push!(germs, DiscGerm(halfedge, leg))
            push!(corners, _corner(owner, other, partner))
            letter = edge_id[halfedge]
            push!(letters, halfedge < other ? letter : -letter)
            if (partner, back) == (start_halfedge, start_leg)
                break
            end
            halfedge, leg = partner, back
        end
        push!(cycles, FaceCycle(germs, corners, letters))
    end
    all(used) || throw(ArgumentError("a sheet end was not used"))
    vertex_count = length(graph.vertex_halfedges)
    disc_count = length(cycles)
    germ_count = sum(length, (cycle.germs for cycle in cycles))
    germ_count == 6 * vertex_count || throw(ArgumentError(
        "face cycles cover $germ_count germs, not $(6 * vertex_count)"))
    ColouredPolyhedron(
        vertex_count, edge_count, disc_count,
        vertex_count - edge_count + disc_count,
        n, germ_count, germ_count, cycles, calibration, OGRAPH_FACE_CYCLE_STATUS)
end

function _require_coloured_ograph_shape(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "face cycles require O0-O3: " * join((check.message for check in failed), "; ")))
    graph
end

function _spanning_tree(edges::Vector{Tuple{Int,Int}}, owner, vertex_count::Int)
    parent = collect(1:vertex_count)
    function find(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        vertex
    end
    tree = Int[]
    for (index, (left, right)) in enumerate(edges)
        a = find(owner[left][1])
        b = find(owner[right][1])
        a == b && continue
        parent[b] = a
        push!(tree, index)
    end
    length(tree) == vertex_count - 1 || throw(ArgumentError(
        "the 1-skeleton is not connected"))
    tree
end

"""Cellular H₁ of the polyhedron.

A spanning tree of the singular graph is contracted. The chain orientation is
the sign stored on `FaceCycle.edge_letters`.
"""
function cellular_homology(graph::ColouredOGraph; calibration::Symbol=:b)
    poly = face_cycles(graph; calibration)
    owner = _face_owner(graph)
    edges, _ = _face_edges(graph)
    tree = _spanning_tree(edges, owner, poly.vertex_count)
    generator_of = zeros(Int, poly.edge_count)
    generator_count = 0
    tree_set = Set(tree)
    for edge in 1:poly.edge_count
        edge in tree_set && continue
        generator_count += 1
        generator_of[edge] = generator_count
    end
    rows = length(poly.cycles)
    matrix = zeros(BigInt, rows, generator_count)
    for (row, cycle) in enumerate(poly.cycles)
        for letter in cycle.edge_letters
            generator = generator_of[abs(letter)]
            generator == 0 && continue
            matrix[row, generator] += sign(letter)
        end
    end
    diagonal = BigInt[]
    if rows > 0 && generator_count > 0
        smith = AbstractAlgebra.snf(AbstractAlgebra.matrix(AbstractAlgebra.ZZ, matrix))
        diagonal = [abs(BigInt(smith[index, index])) for index in 1:min(size(smith)...)]
    end
    nonzero = filter(!iszero, diagonal)
    all(value > 0 for value in nonzero) || throw(ArgumentError(
        "Smith diagonal must be nonnegative"))
    all(nonzero[index + 1] % nonzero[index] == 0 for index in 1:max(length(nonzero) - 1, 0)) ||
        throw(ArgumentError("Smith invariant factors are not divisibility ordered"))
    free_rank = generator_count - length(nonzero)
    free_rank >= 0 || throw(ArgumentError("Smith rank exceeds the generator count"))
    factors = BigInt[value for value in nonzero if value > 1]
    ColourHomology(free_rank, factors)
end

function _perm_sign(values)
    sign = 1
    sequence = collect(Int.(values))
    for left in 1:length(sequence) - 1, right in left + 1:length(sequence)
        sequence[left] == sequence[right] && return 0
        sequence[left] > sequence[right] && (sign = -sign)
    end
    sign
end

"""Induced orientation of the face opposite slot `omitted`.

The positive orientation of the tetrahedron is the cyclic order of the four
slots, begun at either slot of the over pair. The two over slots are
opposite, so the two beginnings differ by a double transposition and are the
same orientation. In that order, the face opposite the i-th slot carries
(−1)^i times the other three slots. The result is three slot indices. It
does not read half-edge ids.
"""
function _positive_face_order(over, omitted::Int)
    slots = Tuple(sort(Int.(over)))
    length(slots) == 2 && slots[2] - slots[1] == 2 ||
        throw(ArgumentError("the over pair must be opposite slots"))
    positive = ntuple(step -> mod1(slots[1] + step - 1, 4), 4)
    index = findfirst(==(omitted), positive)
    index === nothing && throw(ArgumentError("slot $omitted is outside 1:4"))
    others = [positive[step] for step in 1:4 if step != index]
    isodd(index) ? (others[1], others[2], others[3]) : (others[1], others[3], others[2])
end

function _order_sign(reference::NTuple{3,Int}, image::NTuple{3,Int})
    position = Dict(reference[index] => index for index in 1:3)
    _perm_sign(ntuple(index -> position[image[index]], 3))
end

"""Face identifications of the dual ideal triangulation.

Each record sends the three slots of the face opposite `source_slot` to the
three slots of the face opposite `target_slot`. The reversal flag compares
the induced orders of those faces. The half-edge comparison only selects one
end of the undirected edge. `dual_face_pairings` lists the smaller
`(vertex, slot)` first, and exchanging the ends keeps the flag.
"""
function _slot_gluings(graph::ColouredOGraph, calibration::Symbol)
    _require_coloured_ograph_shape(graph)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    gluings = NamedTuple[]
    seen = Set{Tuple{Int,Int}}()
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        source_vertex, source_slot = owner[halfedge]
        target_vertex, target_slot = owner[Int(other)]
        alpha = colours[minmax(halfedge, Int(other))]
        images = Dict{Int,Int}()
        for leg in 1:3
            neighbour = charts[halfedge][leg]
            _, neighbour_slot = owner[neighbour]
            image_leg = _sheet_transport(alpha, leg)
            image_halfedge = charts[Int(other)][image_leg]
            _, image_slot = owner[image_halfedge]
            images[neighbour_slot] = image_slot
        end
        source_positive = _positive_face_order(graph.over_positions[source_vertex], source_slot)
        image_positive = ntuple(index -> images[source_positive[index]], 3)
        target_positive = _positive_face_order(graph.over_positions[target_vertex], target_slot)
        target_vertices = Set(slot for slot in 1:4 if slot != target_slot)
        Set(values(images)) == target_vertices && length(images) == 3 ||
            throw(ArgumentError(
                "face ($source_vertex, $source_slot) does not match three slots"))
        reversing = _order_sign(target_positive, image_positive) == -1
        push!(gluings, (; source_vertex, source_slot, target_vertex, target_slot,
            images, reversing))
        push!(seen, (source_vertex, source_slot))
        push!(seen, (target_vertex, target_slot))
    end
    vertex_count = length(graph.vertex_halfedges)
    expected = Set((vertex, slot) for vertex in 1:vertex_count for slot in 1:4)
    seen == expected || throw(ArgumentError("a tetrahedron face is unpaired"))
    length(gluings) == 2 * vertex_count || throw(ArgumentError(
        "a quadrivalent graph has 2n face pairs"))
    gluings
end

"""Face pairings of the dual ideal triangulation.

Local vertices of each tetrahedron are the four slots. There are `n`
tetrahedra and `2n` pairs.
"""
function dual_face_pairings(graph::ColouredOGraph; calibration::Symbol=:b)
    rows = DualFacePairing[]
    for gluing in _slot_gluings(graph, calibration)
        source = (gluing.source_vertex, gluing.source_slot)
        target = (gluing.target_vertex, gluing.target_slot)
        # The reversal flag is a property of the unordered pair. The record
        # lists the end with the smaller (vertex, slot) first.
        if target < source
            source, target = target, source
        end
        push!(rows, DualFacePairing(source[1], source[2], target[1], target[2],
            gluing.reversing))
    end
    rows
end

"""Orbits of the 6n vertex germs under the sheet identifications.

The count is the number of edges of the dual ideal triangulation.
"""
function dual_edge_orbits(graph::ColouredOGraph; calibration::Symbol=:b)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    parent = Dict{NTuple{3,Int},NTuple{3,Int}}()
    function find(germ)
        while parent[germ] != germ
            parent[germ] = parent[parent[germ]]
            germ = parent[germ]
        end
        germ
    end
    function unite(left, right)
        left_root = find(left)
        right_root = find(right)
        left_root == right_root || (parent[right_root] = left_root)
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for left in 1:3, right in left + 1:4
            germ = (vertex, left, right)
            parent[germ] = germ
        end
    end
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        alpha = colours[minmax(halfedge, Int(other))]
        for leg in 1:3
            tail_partner = charts[halfedge][leg]
            other_leg = _sheet_transport(alpha, leg)
            head_partner = charts[Int(other)][other_leg]
            unite(_corner(owner, halfedge, tail_partner),
                  _corner(owner, Int(other), head_partner))
        end
    end
    orbits = Dict{NTuple{3,Int},Vector{NTuple{3,Int}}}()
    for germ in keys(parent)
        push!(get!(orbits, find(germ), NTuple{3,Int}[]), germ)
    end
    collect(values(orbits))
end

"""Discs as a sorted partition of corners `(vertex, slot_lo, slot_hi)`.

Repeated corners stay repeated. Cyclic order inside a disc is forgotten.
"""
function sorted_germ_partition(discs)
    parts = Vector{NTuple{3,Int}}[sort!(collect(disc)) for disc in discs]
    sort!(parts)
    parts
end

function _port_order(sign::Integer)
    sign == 1 && return (:mi, :po, :mo, :pi)
    sign == -1 && return (:mi, :pi, :mo, :po)
    throw(ArgumentError("vertex sign must be +1 or -1"))
end

const _DISK_PORT_LANES = Dict(
    :mi => (1, 2, 3),
    :pi => (4, 5, 6),
    :mo => (7, 8, 9),
    :po => (10, 11, 12),
)

"""Vertex germs of each `disk_table` disc, as `(vertex, slot_lo, slot_hi)`.

The node-to-port pairs come from `_disk_node_dictionary`. The port-to-slot
pairs come from the Figure 3.4 and Figure 3.5 order at `state.vertex_signs`.
"""
function disk_table_corners(link::VirtualLink, state::OGraphState)
    table = disk_table(link)
    dictionary = _disk_node_dictionary()
    slot_of_port = Dict{Tuple{Int,Symbol},Int}()
    node_ports = Dict{Tuple{Int,Int},NTuple{2,Symbol}}()
    for (vertex, sign) in enumerate(state.vertex_signs)
        for (slot, port) in enumerate(_port_order(sign))
            slot_of_port[(vertex, port)] = slot
        end
        slots = dictionary[sign]
        owners = Dict{Int,Vector{Symbol}}()
        for (port, indices) in _DISK_PORT_LANES
            for index in indices
                push!(get!(owners, slots[index], Symbol[]), port)
            end
        end
        for (node, ports) in owners
            unique_ports = unique(ports)
            length(unique_ports) == 2 || throw(ArgumentError(
                "disk node $node at vertex $vertex meets ports $ports"))
            a, b = sort(collect(unique_ports))
            node_ports[(vertex, node)] = (a, b)
        end
    end
    discs = Vector{NTuple{3,Int}}[]
    for path in table.paths
        germs = NTuple{3,Int}[]
        for node in path
            vertex = (node - 1) ÷ 6 + 1
            local_node = mod1(node, 6)
            port_a, port_b = node_ports[(vertex, local_node)]
            slot_a = slot_of_port[(vertex, port_a)]
            slot_b = slot_of_port[(vertex, port_b)]
            push!(germs, (vertex, min(slot_a, slot_b), max(slot_a, slot_b)))
        end
        push!(discs, germs)
    end
    discs
end

"""Exact partition comparison with `disk_table` on `link`.

`link` is the registered presentation and `colouring.state` supplies the
vertex signs of the Figure 3.4 and Figure 3.5 slot order. The comparison
does not choose the calibration.
"""
function compare_disc_partition(link::VirtualLink, colouring::ClosedNormalColouring;
                                calibration::Symbol=:b)
    poly = face_cycles(colouring.graph; calibration)
    ours = [collect(cycle.corners) for cycle in poly.cycles]
    pipeline = disk_table_corners(link, colouring.state)
    left = sorted_germ_partition(ours)
    right = sorted_germ_partition(pipeline)
    flat = collect(Iterators.flatten(left))
    DiscPartitionComparison(left == right, left, right, length(flat), length(unique(flat)))
end

"""Boundary of a regular neighbourhood of the polyhedron.

One triangle per tetrahedron corner. Its three sides are the three incident
faces, in the induced orientation of the face opposite that corner. Each face
gluing pairs a side with the image side. The vertices of the surface are the
cycles of that rotation system. Every side is paired, so each component is a
closed surface.
"""
function thickening_boundary(graph::ColouredOGraph; calibration::Symbol=:b)
    gluings = _slot_gluings(graph, calibration)
    partner = Dict{NTuple{3,Int},NTuple{3,Int}}()
    reversing_of = Dict{NTuple{3,Int},Bool}()
    for gluing in gluings
        for (source_slot, target_slot) in gluing.images
            left = (gluing.source_vertex, source_slot, gluing.source_slot)
            right = (gluing.target_vertex, target_slot, gluing.target_slot)
            partner[left] = right
            partner[right] = left
            reversing_of[left] = gluing.reversing
            reversing_of[right] = gluing.reversing
        end
    end
    vertex_count = length(graph.vertex_halfedges)
    triangles = Vector{NTuple{3,Int}}[]
    for vertex in 1:vertex_count, slot in 1:4
        order = _positive_face_order(graph.over_positions[vertex], slot)
        push!(triangles, [(vertex, slot, order[index]) for index in 1:3])
    end
    next_of = Dict{NTuple{3,Int},NTuple{3,Int}}()
    owner = Dict{NTuple{3,Int},Int}()
    for (index, triangle) in enumerate(triangles)
        for position in 1:3
            side = triangle[position]
            haskey(partner, side) || throw(ArgumentError("link side $side is unpaired"))
            partner[partner[side]] == side ||
                throw(ArgumentError("link side pairing is not an involution"))
            next_of[side] = triangle[mod1(position + 1, 3)]
            owner[side] = index
        end
    end
    parent = collect(1:length(triangles))
    function find(index)
        while parent[index] != index
            parent[index] = parent[parent[index]]
            index = parent[index]
        end
        index
    end
    for side in keys(partner)
        haskey(owner, side) || throw(ArgumentError("link side $side is not on a corner"))
        left = find(owner[side])
        right = find(owner[partner[side]])
        left == right || (parent[right] = left)
    end
    seen = Set{NTuple{3,Int}}()
    cycle_of = Dict{NTuple{3,Int},Int}()
    cycle_count = 0
    for side in keys(next_of)
        side in seen && continue
        cycle_count += 1
        current = side
        while true
            current in seen && throw(ArgumentError("link vertex walk is not a cycle"))
            push!(seen, current)
            cycle_of[current] = cycle_count
            current = next_of[partner[current]]
            current == side && break
        end
    end
    groups = Dict{Int,Vector{Int}}()
    for index in 1:length(triangles)
        push!(get!(groups, find(index), Int[]), index)
    end
    components = SpineBoundaryComponent[]
    for indices in values(groups)
        sides = NTuple{3,Int}[]
        orientable = true
        for index in indices, side in triangles[index]
            push!(sides, side)
            reversing_of[side] || (orientable = false)
        end
        faces = length(indices)
        edge_count = length(sides) ÷ 2
        vertices = length(unique(cycle_of[side] for side in sides))
        push!(components, SpineBoundaryComponent(
            vertices - edge_count + faces, faces, edge_count, vertices, orientable, true))
    end
    sort!(components, by=component ->
        (component.euler, component.face_count, component.edge_count, component.vertex_count))
    SpineBoundary(components, sum(component.euler for component in components; init=0))
end

