"""Oriented combinatorial ribbon graph and boundary cycles."""

_lexmin(values) = first(sort(collect(values); by=repr))

function _intersection_map(code::HeegaardCode)
    Dict(point.id => point for point in code.intersections)
end

function _curve_vector(code::HeegaardCode)
    vcat(code.alpha, code.beta)
end

function _port(family::Symbol, direction::Symbol)
    family == :alpha ? (direction == :out ? :alpha_out : :alpha_in) :
        (direction == :out ? :beta_out : :beta_in)
end

"""Build the oriented ribbon graph using the fixed sign convention.

The four local ports at a positive crossing are ordered
`alpha_out, beta_out, alpha_in, beta_in` counter-clockwise.  At a negative
crossing the beta ports are exchanged.  The boundary permutation uses the
clockwise (right-turn) predecessor in this order.
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
