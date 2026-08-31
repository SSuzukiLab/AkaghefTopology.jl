"""A bounded fixed-endpoint route search for small alpha-cut fixtures.

This is deliberately a prototype. It searches a finite visibility graph in
the already chosen alpha-cut chart; it does not construct a general planar
embedding and it does not move punctures or apply a radial cleanup pass.
"""

function _router_anchor(anchors, index)
    anchors === nothing && return nothing
    value = anchors isa AbstractDict ? get(anchors, index, nothing) :
        (index <= length(anchors) ? anchors[index] : nothing)
    value === nothing && return nothing
    coordinates = collect(value)
    length(coordinates) == 2 || throw(ArgumentError("router anchor $index must be 2D"))
    point = (Float64(coordinates[1]), Float64(coordinates[2]))
    all(isfinite, point) || throw(ArgumentError("router anchor $index is not finite"))
    point
end

function _router_segment_clear(a, b, punctures)
    _same_point(a, b) && return false
    all(_segment_distance_squared(puncture.center, a, b) >=
        (puncture.radius - 1e-10)^2 for puncture in punctures)
end

function _router_path_clear(points, punctures, routes)
    length(points) >= 2 || return false
    all(_router_segment_clear(points[index], points[index + 1], punctures)
        for index in 1:length(points)-1) || return false
    candidate = Route(:beta, 0, 0, nothing, nothing, points)
    _self_crosses(candidate) && return false
    for route in routes
        for first in _route_segments(Route[candidate]), second in _route_segments(Route[route])
            _segments_cross(first.a, first.b, second.a, second.b) && return false
        end
    end
    !any(_route_enters_puncture(candidate, puncture) for puncture in punctures)
end

function _router_bounds(punctures, ports; margin=1.5)
    points = vcat([puncture.center for puncture in punctures],
                  [port.point for port in ports])
    isempty(points) && throw(ArgumentError("fixed-endpoint router requires punctures or ports"))
    (minimum(point[1] for point in points) - margin,
     maximum(point[1] for point in points) + margin,
     minimum(point[2] for point in points) - margin,
     maximum(point[2] for point in points) + margin)
end

function _router_grid_points(punctures, ports; grid_step, margin)
    minimum_x, maximum_x, minimum_y, maximum_y = _router_bounds(punctures, ports; margin=margin)
    points = NTuple{2,Float64}[]
    for x in minimum_x:grid_step:maximum_x, y in minimum_y:grid_step:maximum_y
        point = (Float64(x), Float64(y))
        any(_segment_distance_squared(puncture.center, point, point) <
            (puncture.radius + 1e-8)^2 for puncture in punctures) || push!(points, point)
    end
    # A fixed endpoint on a puncture boundary needs a short, radial collar
    # before it can enter the visibility grid.  The collar is not a radial
    # postprocess: it is part of the finite candidate graph and is validated
    # like every other segment.
    collar = max(0.18, 1.5 * grid_step)
    puncture_map = Dict(puncture.id => puncture for puncture in punctures)
    for port in ports
        puncture = get(puncture_map, port.boundary_id, nothing)
        puncture === nothing && continue
        dx = port.point[1] - puncture.center[1]
        dy = port.point[2] - puncture.center[2]
        norm = hypot(dx, dy)
        norm > 0 || continue
        point = (puncture.center[1] + (puncture.radius + collar) * dx / norm,
                 puncture.center[2] + (puncture.radius + collar) * dy / norm)
        any(_segment_distance_squared(other.center, point, point) <
            (other.radius + 1e-8)^2 for other in punctures) || push!(points, point)
    end
    points
end

function _router_neighbors(points, punctures, routes; max_edge_length=Inf)
    neighbors = [Tuple{Int,Float64}[] for _ in points]
    for first in eachindex(points)
        for second in first+1:length(points)
            hypot(points[second][1] - points[first][1],
                  points[second][2] - points[first][2]) <= max_edge_length || continue
            _router_segment_clear(points[first], points[second], punctures) || continue
            candidate = Route(:beta, 0, 0, nothing, nothing,
                              [points[first], points[second]])
            _router_path_clear(candidate.points, punctures, routes) || continue
            distance = hypot(points[second][1] - points[first][1],
                             points[second][2] - points[first][2])
            push!(neighbors[first], (second, distance))
            push!(neighbors[second], (first, distance))
        end
    end
    neighbors
end

function _router_shortest_path(points, start_id, finish_id, punctures, routes;
                               anchor=nothing, bend_weight=0.0, grid_step=0.18,
                               max_expansions=25000)
    span = maximum(maximum(point[axis] for point in points) -
                   minimum(point[axis] for point in points) for axis in 1:2)
    neighbors = _router_neighbors(points, punctures, routes;
                                  max_edge_length=min(span, max(3.1 * grid_step, 1.25 * grid_step + 0.18)))
    distances = fill(Inf, length(points))
    bends = fill(typemax(Int), length(points))
    previous = fill(0, length(points))
    visited = falses(length(points))
    distances[start_id] = 0.0
    bends[start_id] = 0
    expansions = 0
    while true
        unvisited = findall(!, visited)
        isempty(unvisited) && break
        current = unvisited[argmin((distances[index], bends[index], index)
                                    for index in unvisited)]
        isfinite(distances[current]) || break
        current == finish_id && break
        visited[current] = true
        expansions += 1
        expansions <= max_expansions || return nothing
        for (next, edge_length) in neighbors[current]
            visited[next] && continue
            candidate_distance = distances[current] + edge_length
            anchor_penalty = anchor === nothing ? 0.0 :
                0.04 * hypot(points[next][1] - anchor[1], points[next][2] - anchor[2])
            candidate_bends = bends[current] + 1
            candidate_distance += anchor_penalty + bend_weight * candidate_bends
            if (candidate_distance, candidate_bends, next) <
               (distances[next], bends[next], next)
                distances[next] = candidate_distance
                bends[next] = candidate_bends
                previous[next] = current
            end
        end
    end
    isfinite(distances[finish_id]) || return nothing
    ids = Int[]
    current = finish_id
    while current != 0
        push!(ids, current)
        current == start_id && break
        current = previous[current]
    end
    last(ids) == start_id || return nothing
    reverse!([points[index] for index in ids])
end

function _router_candidate(points, wire, ports, punctures, routes; anchor=nothing,
                           bend_weight=0.0, grid_step=0.18, max_expansions=25000)
    port_map = Dict(port.id => port for port in ports)
    start = get(port_map, wire.from_port, nothing)
    finish = get(port_map, wire.to_port, nothing)
    (start === nothing || finish === nothing) &&
        throw(ArgumentError("wiring references an unknown endpoint"))
    augmented = NTuple{2,Float64}[start.point, finish.point]
    append!(augmented, points)
    anchor === nothing || push!(augmented, anchor)
    path = _router_shortest_path(augmented, 1, 2, punctures, routes;
                                 anchor=anchor, bend_weight=bend_weight,
                                 grid_step=grid_step,
                                 max_expansions=max_expansions)
    path === nothing && return nothing
    _router_path_clear(path, punctures, routes) ? path : nothing
end

function _router_wire_order(wiring::WiringProblem, ports)
    port_map = Dict(port.id => port for port in ports)
    sortperm(wiring.wires; by=wire -> begin
        start = port_map[wire.from_port].point
        finish = port_map[wire.to_port].point
        (-hypot(start[1] - finish[1], start[2] - finish[2]), wire.beta_id,
         wire.from_port, wire.to_port)
    end)
end

function _feasible_fixed_endpoint_routes(H::HeegaardDiagram,
                                         punctures::Vector{Puncture},
                                         ports::Vector{Port},
                                         wiring::WiringProblem;
                                         anchors=nothing, grid_step=nothing,
                                         margin=1.5, max_attempts=12,
                                         max_expansions=25000)
    wiring.diagram_hash == canonical_hash(H) ||
        throw(ArgumentError("wiring problem does not belong to the diagram"))
    isempty(wiring.wires) && return Route[]
    length(unique(port.id for port in ports)) == length(ports) ||
        throw(ArgumentError("router requires unique port identifiers"))
    step = grid_step === nothing ? min(0.18, max(minimum(p.radius for p in punctures) / 2, 0.06)) : Float64(grid_step)
    step > 0 || throw(ArgumentError("grid_step must be positive"))
    base_points = _router_grid_points(punctures, ports; grid_step=step, margin=margin)
    base_order = _router_wire_order(wiring, ports)
    segment_ids = Dict{Int,Int}()
    per_beta = Dict{Int,Int}()
    for (wire_index, wire) in enumerate(wiring.wires)
        segment_ids[wire_index] = get!(per_beta, wire.beta_id, 0) + 1
        per_beta[wire.beta_id] = segment_ids[wire_index]
    end
    routes = Route[]
    for attempt in 1:max_attempts
        empty!(routes)
        shift = mod(fld(attempt - 1, 2), length(base_order))
        order = circshift(base_order, shift)
        iseven(attempt) && reverse!(order)
        bend_weight = 0.01 * ((attempt - 1) % 4)
        success = true
        for wire_index in order
            wire = wiring.wires[wire_index]
            anchor = _router_anchor(anchors, wire_index)
            path = _router_candidate(base_points, wire, ports, punctures, routes;
                                     anchor=anchor, bend_weight=bend_weight,
                                     grid_step=step,
                                     max_expansions=max_expansions)
            if path === nothing
                success = false
                break
            end
            push!(routes, Route(:beta, wire.beta_id, segment_ids[wire_index],
                                wire.from_port, wire.to_port, path))
        end
        success || continue
        sort!(routes; by=route -> findfirst(wire_index ->
            wiring.wires[wire_index].from_port == route.start_port &&
            wiring.wires[wire_index].to_port == route.end_port, order))
        candidate = HeegaardLayout(canonical_hash(H), :alpha_cut, Dict{String,Any}(),
                                   punctures, ports, routes, _route_crossings(routes),
                                   Dict{String,Any}(),
                                   LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
        report = validate_layout(H, candidate; expected_hash=canonical_hash(H))
        report.valid && return routes
    end
    throw(ArgumentError("bounded fixed-endpoint router found no verified solution after $max_attempts attempts"))
end

feasible_fixed_endpoint_routes(H::HeegaardDiagram, punctures, ports, wiring; kwargs...) =
    _feasible_fixed_endpoint_routes(H, Vector{Puncture}(punctures), Vector{Port}(ports), wiring; kwargs...)
