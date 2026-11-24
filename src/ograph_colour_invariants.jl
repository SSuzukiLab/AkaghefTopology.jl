# cf:[[BP97#^bp97-sec-0703]] §10.1 — cellular chains of the spine; adjustment: the discs are face_cycles under prop321-calibration-b/draft, and this presentation only drops a spanning tree from those letters; reconstructed

"""π₁ presentation of the polyhedron whose discs are `face_cycles`.

The leg chart and the disc partition belong to `face_cycles`. The default
calibration is B, status `prop321-calibration-b/draft`. This function does
not build a second chart. A positive letter is the sign `face_cycles` stored
on the disc. A spanning tree of the singular graph is deleted from the words.
"""
struct ColourPi1
    generator_count::Int
    relators::Vector{Vector{Int}}
    tree_edges::Vector{Int}
    edge_count::Int
end

function _presentation_edges(graph::ColouredOGraph)
    edges = Tuple{Int,Int}[]
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        push!(edges, (halfedge, Int(other)))
    end
    edges
end

function _presentation_owner(graph::ColouredOGraph)
    owner = Dict{Int,Int}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for halfedge in cyclic
            owner[Int(halfedge)] = vertex
        end
    end
    owner
end

function _presentation_tree(edges, owner, vertex_count::Int)
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
        a = find(owner[left])
        b = find(owner[right])
        a == b && continue
        parent[b] = a
        push!(tree, index)
    end
    length(tree) == vertex_count - 1 || throw(ArgumentError(
        "the 1-skeleton is not connected"))
    tree
end

"""Presentation of `face_cycles(graph)`.

The call uses the default calibration. The singular edges are ordered by the
smaller half-edge id. The contracted tree is Kruskal in that order, the same
order `cellular_homology` contracts. Abelianization of these relators is that
homology.
"""
function coloured_pi1(graph::ColouredOGraph)
    poly = face_cycles(graph)
    edges = _presentation_edges(graph)
    tree = _presentation_tree(edges, _presentation_owner(graph), poly.vertex_count)
    generator_of = zeros(Int, poly.edge_count)
    generator_count = 0
    tree_set = Set(tree)
    for edge in 1:poly.edge_count
        edge in tree_set && continue
        generator_count += 1
        generator_of[edge] = generator_count
    end
    relators = Vector{Int}[]
    for cycle in poly.cycles
        relator = Int[]
        for letter in cycle.edge_letters
            generator = generator_of[abs(letter)]
            generator == 0 && continue
            push!(relator, sign(letter) * generator)
        end
        push!(relators, relator)
    end
    ColourPi1(generator_count, relators, tree, poly.edge_count)
end
