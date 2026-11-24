# cf:[[BP97#^bp97-sec-0551]] §2.2 Figs. 2.8–2.13 Prop. 2.2.1 — Θ₀ from Fig. 2.12, its Z₂ strand colours, and the Fig. 2.8 ribbon; adjustment: PD slot 1 (−i) is the south arm of that cross; unverified

"""Boundary of a coloured o-graph by the Section 2.2 ribbon.

`reading_status` is `:provisional`. The strand permutations and the four
trivalent vertices are read from Figures 2.12 and 2.9. Which arm of that
cross is the incoming under-passage is not settled by those figures.
"""
module OGraphColourBoundary

export BoundarySurface, BoundaryComponent,
       boundary_from_gauss, boundary_from_coloured, prop321_colour

const READING_STATUS = :provisional

struct BoundaryComponent
    theta_vertices::Int
    theta_edges::Int
    ribbon_cycles::Int
    euler::Int
    # Lengths of one directed traversal of each unoriented boundary circle.
    orbit_degrees::Vector{Int}
end

struct BoundarySurface
    reading::Symbol
    components::Vector{BoundaryComponent}
end

# Page-left-to-right stub order when the arm points down the page.
# Stubs: N0 west-north, N1 mid-north, N2 east-north, and likewise S, W, E.
const DOWN_ORDER = (
    S = (:S0, :S1, :S2),
    N = (:N2, :N1, :N0),
    E = (:E2, :E1, :E0),
    W = (:W0, :W1, :W2),
)

const STUBS = (:N0, :N1, :N2, :S0, :S1, :S2, :W0, :W1, :W2, :E0, :E1, :E2)
const STUB_INDEX = Dict(stub => i for (i, stub) in enumerate(STUBS))

# Θ vertex that owns each stub, then the counterclockwise order at that vertex.
const OWNER = (
    N0 = :W, N1 = :N, N2 = :E,
    S0 = :W, S1 = :S, S2 = :E,
    W0 = :N, W1 = :W, W2 = :S,
    E0 = :N, E1 = :E, E2 = :S,
)
const ROTATION = (
    N = (:N1, :W0, :E0),
    S = (:E2, :W2, :S1),
    W = (:N0, :W1, :S0),
    E = (:N2, :S2, :E1),
)
const THETA_VERTICES = (:N, :E, :S, :W)

# perm[i] is the bottom position, 0-based, of the strand that starts at the
# top in position i-1. z2 is indexed by that bottom position.
# Same-type rows are the top row of Fig. 2.12, keyed by (both-above ? +1 : -1) * colour.
# Mixed rows are the bottom row, drawn with the above end at the top.
const SAME_ROW = Dict(
    0 => ((0, 1, 2), (0, 0, 0)),
    1 => ((2, 0, 1), (1, 1, 0)),
    -1 => ((1, 2, 0), (0, 1, 1)),
)
const MIXED_ROW = Dict(
    0 => ((2, 1, 0), (1, 1, 1)),
    1 => ((0, 2, 1), (1, 0, 0)),
    -1 => ((1, 0, 2), (0, 0, 1)),
)

function _mod3(value::Integer)
    reduced = mod(Int(value), 3)
    reduced == 2 ? -1 : reduced
end

"""Colour in `{-1,0,1}` from the two passage marks and the edge orientation.

`first_above` is the tail. Above means an over-passage.
"""
function prop321_colour(first_above::Bool, second_above::Bool)
    first_above == second_above && return 0
    first_above && return -1
    return 1
end

# cf:[[BP97#^bp97-sec-0553]] Fig. 2.12 bottom row — the undrawn mixed pattern, by a half-turn of the drawn row; the Z₂ mark stays on the strand; reconstructed
function _rotate180(perm::NTuple{3,Int}, z2::NTuple{3,Int})
    inverse = ntuple(position -> findfirst(==(position - 1), perm) - 1, 3)
    new_perm = ntuple(j -> 2 - inverse[3 - (j - 1)], 3)
    # Old top 2-(j-1) becomes new bottom j-1, and it still carries z2 at its old bottom.
    new_z2 = ntuple(j -> z2[perm[4 - j] + 1], 3)
    new_perm, new_z2
end

function _strand_table(first_above::Bool, second_above::Bool, colour::Int)
    colour = _mod3(colour)
    if first_above == second_above
        sign = first_above ? 1 : -1
        return SAME_ROW[sign * colour]
    end
    perm, z2 = MIXED_ROW[colour]
    if first_above
        return perm, z2
    end
    _rotate180(perm, z2)
end

function _ports(vertex_sign::Integer)
    vertex_sign == 1 && return (:mi, :po, :mo, :pi)
    vertex_sign == -1 && return (:mi, :pi, :mo, :po)
    throw(ArgumentError("vertex sign must be +1 or -1"))
end

const SLOT_DIRECTION = (:S, :E, :N, :W)

function _above(port::Symbol)
    port === :po || port === :pi
end

"""One oriented edge of a coloured o-graph, in the plain input of this lane."""
struct PlainEdge
    colour::Int
    tail_vertex::Int
    tail_port::Symbol
    head_vertex::Int
    head_port::Symbol
end

struct PlainGraph
    signs::Vector{Int}
    edges::Vector{PlainEdge}
end

function boundary_from_gauss(words, vertex_signs; alter_edge::Union{Nothing,Integer}=nothing)
    signs = Int[sign == 1 || sign == -1 ? Int(sign) :
                throw(ArgumentError("vertex sign must be +1 or -1"))
                for sign in vertex_signs]
    n = length(signs)
    edges = PlainEdge[]
    for word in words
        m = length(word)
        m == 0 && continue
        for index in 1:m
            tail = Int(word[index])
            head = Int(word[mod1(index + 1, m)])
            tail_vertex = abs(tail)
            head_vertex = abs(head)
            (1 <= tail_vertex <= n && 1 <= head_vertex <= n) ||
                throw(ArgumentError("Gauss label out of range"))
            tail_port = tail > 0 ? :po : :mo
            head_port = head > 0 ? :pi : :mi
            colour = prop321_colour(tail > 0, head > 0)
            push!(edges, PlainEdge(colour, tail_vertex, tail_port, head_vertex, head_port))
        end
    end
    if alter_edge !== nothing
        index = Int(alter_edge)
        1 <= index <= length(edges) || throw(ArgumentError("alter_edge is out of range"))
        edge = edges[index]
        edges[index] = PlainEdge(_mod3(edge.colour + 1), edge.tail_vertex, edge.tail_port,
                                  edge.head_vertex, edge.head_port)
    end
    _boundary(PlainGraph(signs, edges))
end

function boundary_from_coloured(graph; alter_edge::Union{Nothing,Integer}=nothing)
    mate = graph.mate
    owners = Dict{Int,Tuple{Int,Int}}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for (slot, halfedge) in enumerate(cyclic)
            owners[Int(halfedge)] = (vertex, slot)
        end
    end
    colour_at = Dict{Int,Int}()
    for record in graph.edge_colours
        colour_at[Int(record.halfedge)] = _mod3(record.colour)
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        length(cyclic) == 4 || throw(ArgumentError("vertex $vertex is not quadrivalent"))
        overs = Tuple(sort(Int.(graph.over_positions[vertex])))
        overs == (1, 3) || overs == (2, 4) ||
            throw(ArgumentError("over marks at vertex $vertex are not an opposite pair"))
    end
    # A coloured o-graph has no edge direction. The smaller half-edge is the tail.
    direct = Tuple{Int,Symbol,Bool,Int,Symbol,Bool,Int}[]
    for (halfedge, other) in enumerate(mate)
        halfedge < Int(other) || continue
        tv, tslot = owners[halfedge]
        hv, hslot = owners[Int(other)]
        tail_over = tslot in graph.over_positions[tv]
        head_over = hslot in graph.over_positions[hv]
        colour = get(colour_at, halfedge, get(colour_at, Int(other), nothing))
        colour === nothing && throw(ArgumentError("edge ($halfedge, $other) has no colour"))
        push!(direct, (tv, SLOT_DIRECTION[tslot], tail_over, hv, SLOT_DIRECTION[hslot], head_over, colour))
    end
    if alter_edge !== nothing
        index = Int(alter_edge)
        1 <= index <= length(direct) || throw(ArgumentError("alter_edge is out of range"))
        row = direct[index]
        direct[index] = (row[1], row[2], row[3], row[4], row[5], row[6], _mod3(row[7] + 1))
    end
    _boundary_directed(length(graph.vertex_halfedges), direct)
end

function _direction(sign::Integer, port::Symbol)
    slot = findfirst(==(port), _ports(sign))
    slot === nothing && throw(ArgumentError("port $port is not on this vertex"))
    SLOT_DIRECTION[slot]
end

function _boundary(graph::PlainGraph)
    rows = Tuple{Int,Symbol,Bool,Int,Symbol,Bool,Int}[]
    for edge in graph.edges
        push!(rows, (
            edge.tail_vertex, _direction(graph.signs[edge.tail_vertex], edge.tail_port),
            _above(edge.tail_port),
            edge.head_vertex, _direction(graph.signs[edge.head_vertex], edge.head_port),
            _above(edge.head_port),
            edge.colour,
        ))
    end
    _boundary_directed(length(graph.signs), rows)
end

function _boundary_directed(n_vertices::Integer, rows)
    n = Int(n_vertices)
    nstubs = 12 * n
    partner = zeros(Int, nstubs)
    twist = zeros(Int, nstubs)
    function stub_id(vertex, name)
        (vertex - 1) * 12 + STUB_INDEX[name]
    end
    for (tail_v, tail_dir, tail_over, head_v, head_dir, head_over, colour) in rows
        perm, z2 = _strand_table(tail_over, head_over, colour)
        tail_order = DOWN_ORDER[tail_dir]
        head_order = DOWN_ORDER[head_dir]
        for top in 0:2
            bottom = perm[top + 1]
            # The head sits at the bottom of the panel, so left and right swap.
            head_name = head_order[3 - bottom]
            tail_name = tail_order[top + 1]
            a = stub_id(tail_v, tail_name)
            b = stub_id(head_v, head_name)
            (partner[a] == 0 && partner[b] == 0) ||
                throw(ArgumentError("stub $tail_name or $head_name is used twice"))
            partner[a] = b
            partner[b] = a
            twist[a] = z2[bottom + 1]
            twist[b] = z2[bottom + 1]
        end
    end
    all(!=(0), partner) || throw(ArgumentError("a stub was not paired"))

    ccw = zeros(Int, nstubs)
    cw = zeros(Int, nstubs)
    theta_of_stub = zeros(Int, nstubs)
    for vertex in 1:n
        for (local_index, name) in enumerate(THETA_VERTICES)
            theta_id = (vertex - 1) * 4 + local_index
            order = ROTATION[name]
            ids = ntuple(i -> stub_id(vertex, order[i]), 3)
            for i in 1:3
                ccw[ids[i]] = ids[mod1(i + 1, 3)]
                cw[ids[i]] = ids[mod1(i - 1, 3)]
                theta_of_stub[ids[i]] = theta_id
            end
        end
    end

    image = zeros(Int, nstubs * 2)
    for stub in 1:nstubs, side in 0:1
        arrival = partner[stub]
        arrival_side = twist[stub] == 1 ? 1 - side : side
        nxt = arrival_side == 0 ? ccw[arrival] : cw[arrival]
        image[(stub - 1) * 2 + side + 1] = (nxt - 1) * 2 + arrival_side + 1
    end
    sort(image) == collect(1:nstubs * 2) ||
        throw(ArgumentError("the ribbon boundary is not a permutation of the sides"))

    # cf:[[BP97#^bp97-sec-0551]] paragraph before Fig. 2.8 — one disc on each
    # boundary circle of the ribbon; adjustment: R(h, s) = (partner(h), 1 ⊕ s ⊕ twist(h))
    # pairs the two directions, using that the two turns at a vertex are inverses
    # and that the Z₂ mark agrees on the two ends of a strand; reconstructed
    # R² = 1 and (R T)² = 1, so each geometric circle is one reversal orbit.
    for stub in 1:nstubs, side in 0:1
        forward = image[(stub - 1) * 2 + side + 1]
        a = (forward - 1) ÷ 2 + 1
        b = (forward - 1) % 2
        a, b = partner[a], xor(1, b, twist[a])
        stepped = image[(a - 1) * 2 + b + 1]
        c = (stepped - 1) ÷ 2 + 1
        d = (stepped - 1) % 2
        c, d = partner[c], xor(1, d, twist[c])
        (c == stub && d == side) ||
            throw(ArgumentError("reversal does not invert the side-walk"))
    end
    counted = falses(nstubs, 2)
    cycle_theta = Vector{Int}[]
    for start in 1:nstubs, side in 0:1
        counted[start, side + 1] && continue
        cursor = start
        cursor_side = side
        visited = Int[]
        states = Tuple{Int,Int}[]
        while true
            length(states) <= 2 * nstubs ||
                throw(ArgumentError("a side-walk did not close"))
            push!(states, (cursor, cursor_side))
            push!(visited, theta_of_stub[cursor])
            arrival = partner[cursor]
            arrival_side = twist[cursor] == 1 ? 1 - cursor_side : cursor_side
            cursor = arrival_side == 0 ? ccw[arrival] : cw[arrival]
            cursor_side = arrival_side
            cursor == start && cursor_side == side && break
        end
        for (stub, side_bit) in states
            counted[stub, side_bit + 1] = true
            reverse_stub = partner[stub]
            reverse_side = xor(1, side_bit, twist[stub])
            counted[reverse_stub, reverse_side + 1] = true
        end
        push!(cycle_theta, visited)
    end
    all(counted) || throw(ArgumentError("a side-state was left out of the reversal orbits"))

    parent = collect(1:4n)
    function find(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        vertex
    end
    for stub in 1:nstubs
        a = find(theta_of_stub[stub])
        b = find(theta_of_stub[partner[stub]])
        a == b || (parent[b] = a)
    end
    groups = Dict{Int,Vector{Int}}()
    for vertex in 1:4n
        push!(get!(groups, find(vertex), Int[]), vertex)
    end
    edge_count = Dict{Int,Int}()
    for stub in 1:nstubs
        stub < partner[stub] || continue
        root = find(theta_of_stub[stub])
        edge_count[root] = get(edge_count, root, 0) + 1
    end
    cycle_count = Dict{Int,Int}()
    degree_lists = Dict{Int,Vector{Int}}()
    for visited in cycle_theta
        root = find(visited[1])
        cycle_count[root] = get(cycle_count, root, 0) + 1
        push!(get!(degree_lists, root, Int[]), length(visited))
    end
    components = BoundaryComponent[]
    for (root, vertices) in sort!(collect(groups); by=first)
        v = length(vertices)
        e = get(edge_count, root, 0)
        b = get(cycle_count, root, 0)
        degs = sort!(get(degree_lists, root, Int[]))
        sum(degs) == 2 * e ||
            throw(ArgumentError("boundary degrees sum to $(sum(degs)) on a component with $e edges"))
        push!(components, BoundaryComponent(v, e, b, v - e + b, degs))
    end
    BoundarySurface(READING_STATUS, components)
end

end
