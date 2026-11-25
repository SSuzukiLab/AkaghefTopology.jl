"""Oriented combinatorial ribbon graph and boundary cycles."""

_lexmin(values) = first(sort(collect(values); by=repr))

function _intersection_map(code::CombinatorialHeegaardCode)
    Dict(point.id => point for point in code.intersections)
end

function _curve_vector(code::CombinatorialHeegaardCode)
    vcat(code.alpha, code.beta)
end

function _port(family::Symbol, direction::Symbol)
    family == :alpha ? (direction == :out ? :alpha_out : :alpha_in) :
        (direction == :out ? :beta_out : :beta_in)
end

"""
    ribbon_graph(object) -> RibbonGraph

Construct the oriented ribbon graph of a legacy curve arrangement.

# Arguments
`object` is a legacy combinatorial code, decorated code or cellular diagram. The caller supplies consistent incidence in the [legacy representation](../specs/heegaard-terminology.md#legacy-cellular-records). Unknown intersections, duplicate ports and missing outgoing darts raise `ArgumentError`; unsupported objects also raise `ArgumentError`. This operation does not run every validation level.

# Returns
A `RibbonGraph` with vertices, directed curve edges, paired darts, local rotations and the boundary-successor permutation.

# Guarantees
1. On valid input, the ribbon graph represents the oriented regular neighborhood of the curve arrangement, in the [legacy port convention](../specs/heegaard-terminology.md#legacy-cut-indices). [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Empty curves have no graph edges. The ribbon graph alone does not supply a complement partition, complete cut systems or a 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cut-indices — ribbon ports and boundary numbering.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> g = H.ribbon_graph(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)])); (length(g.vertices), length(g.edges), length(g.darts))
(1, 2, 4)
```

Role: convert
"""
function ribbon_graph(object)
    code = _undecorated(object)
    point_map = _intersection_map(code)
    vertices = sort(collect(keys(point_map)))
    edges = RibbonEdge[]
    darts = RibbonDart[]
    next_edge = 1
    for curve in _curve_vector(code)
        m = length(curve.visits)
        m == 0 && continue
        for k in 1:m
            from = curve.visits[k]
            to = curve.visits[mod1(k + 1, m)]
            haskey(point_map, from) || throw(ArgumentError("curve $(curve.label) references unknown intersection $from"))
            haskey(point_map, to) || throw(ArgumentError("curve $(curve.label) references unknown intersection $to"))
            from_port = _port(curve.family, :out)
            to_port = _port(curve.family, :in)
            push!(edges, RibbonEdge(next_edge, curve.family, curve.id, from, to, from_port, to_port))
            push!(darts, RibbonDart(2next_edge - 1, next_edge, from, to, from_port, to_port))
            push!(darts, RibbonDart(2next_edge, next_edge, to, from, to_port, from_port))
            next_edge += 1
        end
    end

    rotation = Dict{Int,Vector{Symbol}}()
    for point in code.intersections
        rotation[point.id] = point.sign == 1 ?
            Symbol[:alpha_out, :beta_out, :alpha_in, :beta_in] :
            Symbol[:alpha_out, :beta_in, :alpha_in, :beta_out]
    end

    outgoing = Dict{Tuple{Int,Symbol},Int}()
    for dart in darts
        key = (dart.origin, dart.origin_port)
        haskey(outgoing, key) && throw(ArgumentError("duplicate ribbon port at $(key[1])/$(key[2])"))
        outgoing[key] = dart.id
    end
    next_dart = zeros(Int, length(darts))
    for dart in darts
        rot = get(rotation, dart.target, Symbol[])
        isempty(rot) && throw(ArgumentError("ribbon dart targets unknown vertex $(dart.target)"))
        idx = findfirst(==(dart.target_port), rot)
        idx === nothing && throw(ArgumentError("missing target port $(dart.target_port) at $(dart.target)"))
        port = rot[mod1(idx - 1, length(rot))]
        next_dart[dart.id] = get(outgoing, (dart.target, port), 0)
        next_dart[dart.id] == 0 && throw(ArgumentError("missing outgoing ribbon dart at $(dart.target)/$port"))
    end
    RibbonGraph(vertices, edges, darts, rotation, next_dart)
end

function boundary_cycles(graph::RibbonGraph)
    visited = falses(length(graph.darts))
    cycles = Vector{Vector{Int}}()
    for dart in graph.darts
        visited[dart.id] && continue
        cycle = Int[]
        current = dart.id
        while !visited[current]
            visited[current] = true
            push!(cycle, current)
            current = graph.next_dart[current]
        end
        current == first(cycle) || throw(ArgumentError("ribbon boundary permutation is not a disjoint cycle"))
        push!(cycles, cycle)
    end
    cycles
end

"""
    boundary_cycles(object) -> Vector{Vector{Int}}

Return the oriented boundary circuits of a legacy curve arrangement's ribbon neighbourhood.

# Arguments
`object` is a `CombinatorialHeegaardCode`, `DecoratedHeegaardCode` or
`CellularHeegaardDiagram`; alternatively supply a `RibbonGraph` directly.
The caller must supply consistent incidence and a valid boundary permutation
as described in [legacy cellular records](../specs/heegaard-terminology.md#legacy-cellular-records).
A traversal that joins an already visited circuit throws `ArgumentError`
explaining that the boundary permutation is not disjoint cycles. Invalid dart
indices can raise `BoundsError`; incidence errors from `ribbon_graph` propagate.

# Returns
Fresh vectors of dart IDs, in the
[legacy boundary order](../specs/heegaard-terminology.md#legacy-cellular-records).

# Guarantees
1. For a valid ribbon graph, every dart occurs once among the circuits, and
   each circuit follows the graph's boundary successor cyclically.
   [proved/agent: by construction]
2. The circuits bound the oriented ribbon neighbourhood; identifying which
   circuits bound one complementary region requires a separate decoration.
   [proved/agent: by construction]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Circuits alone do not determine a closed surface or Heegaard splitting.
This operation does not fully validate a hand-built graph or use LC face numbering.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy cellular 記録 — ribbon boundary circuits and numbering.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = ribbon_graph(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]));

julia> sort(vcat(boundary_cycles(g)...)) == collect(1:length(g.darts))
true
```

See also [`ribbon_graph`](@ref), [`RegionPartition`](@ref), [`complement_regions`](@ref).

Role: access
"""
boundary_cycles(object) = boundary_cycles(ribbon_graph(object))

function _dart_map(graph::RibbonGraph)
    Dict(dart.id => dart for dart in graph.darts)
end

function _cycle_signature(graph::RibbonGraph, cycle::AbstractVector{<:Integer};
                          intersection_map=identity, curve_map=identity)
    darts = _dart_map(graph)
    raw = [(intersection_map(darts[id].origin),
            darts[id].origin_port,
            curve_map(darts[id].edge),
            intersection_map(darts[id].target),
            darts[id].target_port) for id in cycle]
    isempty(raw) && return raw
    rotations = [vcat(raw[i:end], raw[1:i-1]) for i in eachindex(raw)]
    _lexmin(rotations)
end

"Equality of cycles up to cyclic rotation."
function _cyclic_equal(first_cycle, second_cycle)
    length(first_cycle) == length(second_cycle) || return false
    isempty(first_cycle) && return true
    any(all(first_cycle[index] == second_cycle[mod1(index + shift, length(second_cycle))]
            for index in eachindex(first_cycle))
        for shift in 0:length(second_cycle)-1)
end
