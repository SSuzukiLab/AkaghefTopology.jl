"""Constraint-preserving mechanical relaxation of verified alpha-cut routes."""

_vadd(a, b) = (a[1] + b[1], a[2] + b[2])
_vsub(a, b) = (a[1] - b[1], a[2] - b[2])
_vscale(value, a) = (value * a[1], value * a[2])
_vnorm(a) = hypot(a[1], a[2])

function _unit(a)
    length = _vnorm(a)
    length <= eps(Float64) ? (0.0, 0.0) : _vscale(inv(length), a)
end

function _densify_route(points; max_spacing)
    length(points) <= 2 && return copy(points)
    dense = NTuple{2,Float64}[points[1]]
    for segment in 1:length(points)-1
        first, second = points[segment], points[segment + 1]
        # The endpoint guard segments are deliberately retained exactly.  They
        # leave/enter the puncture boundary along an already certified PL edge.
        pieces = segment in (1, length(points) - 1) ? 1 :
                 max(1, ceil(Int, _vnorm(_vsub(second, first)) / max_spacing))
        for part in 1:pieces
            parameter = part / pieces
            push!(dense, ((1 - parameter) * first[1] + parameter * second[1],
                          (1 - parameter) * first[2] + parameter * second[2]))
        end
    end
    dense
end

function _closest_on_segment(point, first, second)
    direction = _vsub(second, first)
    denominator = direction[1]^2 + direction[2]^2
    denominator <= eps(Float64) && return first
    parameter = clamp(((point[1] - first[1]) * direction[1] +
                       (point[2] - first[2]) * direction[2]) / denominator, 0.0, 1.0)
    _vadd(first, _vscale(parameter, direction))
end

function _segment_distance(first_a, first_b, second_a, second_b)
    _segments_cross(first_a, first_b, second_a, second_b) && return 0.0
    minimum((sqrt(_segment_distance_squared(first_a, second_a, second_b)),
             sqrt(_segment_distance_squared(first_b, second_a, second_b)),
             sqrt(_segment_distance_squared(second_a, first_a, first_b)),
             sqrt(_segment_distance_squared(second_b, first_a, first_b))))
end

function _closest_on_segment_with_parameter(point, first, second)
    direction = _vsub(second, first)
    denominator = direction[1]^2 + direction[2]^2
    denominator <= eps(Float64) && return first, 0.0
    parameter = clamp(((point[1] - first[1]) * direction[1] +
                       (point[2] - first[2]) * direction[2]) / denominator, 0.0, 1.0)
    _vadd(first, _vscale(parameter, direction)), parameter
end

function _closest_segment_pair(first_a, first_b, second_a, second_b)
    candidates = Tuple{Float64,NTuple{2,Float64},NTuple{2,Float64},Float64,Float64}[]
    projected, parameter = _closest_on_segment_with_parameter(first_a, second_a, second_b)
    push!(candidates, (_vnorm(_vsub(first_a, projected)), first_a, projected, 0.0, parameter))
    projected, parameter = _closest_on_segment_with_parameter(first_b, second_a, second_b)
    push!(candidates, (_vnorm(_vsub(first_b, projected)), first_b, projected, 1.0, parameter))
    projected, parameter = _closest_on_segment_with_parameter(second_a, first_a, first_b)
    push!(candidates, (_vnorm(_vsub(projected, second_a)), projected, second_a, parameter, 0.0))
    projected, parameter = _closest_on_segment_with_parameter(second_b, first_a, first_b)
    push!(candidates, (_vnorm(_vsub(projected, second_b)), projected, second_b, parameter, 1.0))
    candidates[argmin(first.(candidates))]
end

function _clearance_barrier(distance, clearance, strength)
    distance >= clearance && return 0.0
    0.5 * strength * (clearance - max(distance, 0.0))^2
end

function _clearance_force(distance, clearance, strength)
    distance >= clearance && return 0.0
    strength * (clearance - max(distance, 0.0))
end

function _multiscale_density(point, route_index, point_routes, spread_range)
    density = 0.0
    gradient = (0.0, 0.0)
    # The widest kernel reaches across a puncture corridor.  Without this term
    # the force only reacts after neighboring routes have already formed a
    # visually collapsed bundle.
    scales = (0.5 * spread_range, 1.5 * spread_range, 4.0 * spread_range)
    scale_weights = (0.15, 0.35, 0.50)
    for other_route in eachindex(point_routes)
        other_route == route_index && continue
        points = point_routes[other_route]
        for segment in 1:length(points)-1
            closest = _closest_on_segment(point, points[segment], points[segment + 1])
            difference = _vsub(point, closest)
            squared_distance = difference[1]^2 + difference[2]^2
            segment_length = _vnorm(_vsub(points[segment + 1], points[segment]))
            for (sigma, weight) in zip(scales, scale_weights)
                kernel = weight * segment_length / sigma *
                         exp(-squared_distance / (2sigma^2))
                density += kernel
                gradient = _vadd(gradient, _vscale(-kernel / sigma^2, difference))
            end
        end
    end
    density, gradient
end

function _polyline_clearance(point_routes, punctures)
    minimum_clearance = Inf
    for first_route in eachindex(point_routes)
        first_points = point_routes[first_route]
        for first_segment in 1:length(first_points)-1
            minimum_clearance = min(minimum_clearance,
                _vnorm(_vsub(first_points[first_segment + 1], first_points[first_segment])) / 2)
            for second_route in first_route:length(point_routes)
                second_points = point_routes[second_route]
                start_segment = second_route == first_route ? first_segment + 2 : 1
                for second_segment in start_segment:length(second_points)-1
                    minimum_clearance = min(minimum_clearance,
                        _segment_distance(first_points[first_segment], first_points[first_segment + 1],
                                          second_points[second_segment], second_points[second_segment + 1]))
                end
            end
            # Endpoint guard segments are fixed and already certified against
            # their incident puncture, so only movable interior segments enter
            # the global step bound.
            if first_segment != 1 && first_segment != length(first_points) - 1
                for puncture in punctures
                    margin = sqrt(_segment_distance_squared(
                        puncture.center, first_points[first_segment],
                        first_points[first_segment + 1])) - puncture.radius
                    minimum_clearance = min(minimum_clearance, max(margin, 0.0))
                end
            end
        end
    end
    isfinite(minimum_clearance) ? minimum_clearance : 1.0
end

function _route_pair_clearance(point_routes)
    minimum_clearance = Inf
    for first_route in eachindex(point_routes)
        first_points = point_routes[first_route]
        for second_route in first_route+1:length(point_routes)
            second_points = point_routes[second_route]
            for first_segment in 1:length(first_points)-1,
                second_segment in 1:length(second_points)-1
                minimum_clearance = min(minimum_clearance,
                    _segment_distance(first_points[first_segment],
                                      first_points[first_segment + 1],
                                      second_points[second_segment],
                                      second_points[second_segment + 1]))
            end
        end
    end
    isfinite(minimum_clearance) ? minimum_clearance : Inf
end

_polyline_length(points) = sum(_vnorm(_vsub(points[index + 1], points[index]))
                               for index in 1:length(points)-1; init=0.0)

function _point_in_polygon(point, polygon)
    inside = false
    previous = length(polygon)
    for index in eachindex(polygon)
        first, second = polygon[index], polygon[previous]
        if (first[2] > point[2]) != (second[2] > point[2]) &&
           point[1] < (second[1] - first[1]) * (point[2] - first[2]) /
                      (second[2] - first[2]) + first[1]
            inside = !inside
        end
        previous = index
    end
    inside
end

function _shortcut_region_clear(point_routes, route_index, first, last, punctures)
    polygon = point_routes[route_index][first:last]
    length(polygon) >= 3 || return true
    for puncture in punctures
        _point_in_polygon(puncture.center, polygon) && return false
        for index in eachindex(polygon)
            next = index == length(polygon) ? 1 : index + 1
            sqrt(_segment_distance_squared(puncture.center, polygon[index], polygon[next])) <=
                puncture.radius && return false
        end
    end
    for other_route in eachindex(point_routes)
        for (index, point) in enumerate(point_routes[other_route])
            other_route == route_index && first <= index <= last && continue
            _point_in_polygon(point, polygon) && return false
        end
    end
    true
end

function _global_shortcut_polylines(H, layout, point_routes; minimum_gain=1e-3)
    current = deepcopy(point_routes)
    before = sum(_polyline_length, current)
    accepted = 0
    while true
        changed = false
        for route_index in eachindex(current)
            points = current[route_index]
            maximum_span = length(points) - 3
            maximum_span < 3 && continue
            spans = unique(vcat(maximum_span, reverse([2^power for power in 2:10
                                                       if 2^power <= maximum_span])))
            for span in spans, first in 2:length(points)-span-1
                last = first + span
                old_length = _polyline_length(points[first:last])
                chord_length = _vnorm(_vsub(points[last], points[first]))
                old_length - chord_length > minimum_gain || continue
                _shortcut_region_clear(current, route_index, first, last,
                                       layout.punctures) || continue
                candidate_routes = deepcopy(current)
                candidate_routes[route_index] = vcat(points[1:first], points[last:end])
                _validated_polyline_candidate(H, layout, candidate_routes) === nothing && continue
                current = candidate_routes
                accepted += 1
                changed = true
                break
            end
            changed && break
        end
        changed || break
    end
    current, (accepted=accepted, before=before, after=sum(_polyline_length, current))
end

function _reduce_cut_word(tokens)
    reduced = Int[]
    for token in tokens
        !isempty(reduced) && last(reduced) == -token ? pop!(reduced) : push!(reduced, token)
    end
    reduced
end

function _path_cut_word(points, punctures, radii)
    center_y = sum(puncture.center[2] for puncture in punctures) / length(punctures)
    tokens = Int[]
    for segment in 1:length(points)-1
        first, second = points[segment], points[segment + 1]
        delta_x = second[1] - first[1]
        abs(delta_x) <= 1e-12 && continue
        events = Tuple{Float64,Int}[]
        for (index, puncture) in enumerate(punctures)
            parameter = (puncture.center[1] - first[1]) / delta_x
            1e-10 < parameter < 1 - 1e-10 || continue
            crossing_y = first[2] + parameter * (second[2] - first[2])
            outward = puncture.center[2] >= center_y ? 1.0 : -1.0
            cut_start = puncture.center[2] + outward * radii[index]
            outward * (crossing_y - cut_start) > 1e-10 || continue
            push!(events, (parameter, delta_x > 0 ? index : -index))
        end
        sort!(events; by=event -> event[1])
        append!(tokens, last.(events))
    end
    _reduce_cut_word(tokens)
end

function _visibility_radii(start_point, end_point, punctures, clearance)
    [begin
        endpoint_margin = min(_vnorm(_vsub(start_point, puncture.center)) - puncture.radius,
                              _vnorm(_vsub(end_point, puncture.center)) - puncture.radius)
        puncture.radius + min(clearance, 0.8 * max(endpoint_margin, 0.0))
    end for puncture in punctures]
end

function _visible_disk_segment(first, second, punctures, radii)
    for (puncture, radius) in zip(punctures, radii)
        sqrt(_segment_distance_squared(puncture.center, first, second)) >= radius - 1e-8 ||
            return false
    end
    true
end

function _visible_route_segment(first, second, barrier_routes, minimum_clearance)
    for points in barrier_routes, segment in 1:length(points)-1
        _segments_cross(first, second, points[segment], points[segment + 1]) && return false
        _segment_distance(first, second, points[segment], points[segment + 1]) >=
            minimum_clearance || return false
    end
    true
end

function _barrier_navigation_nodes(barrier_routes, offset)
    nodes = NTuple{2,Float64}[]
    for points in barrier_routes
        length(points) < 2 && continue
        for index in eachindex(points)
            previous = points[max(index - 1, 1)]
            following = points[min(index + 1, length(points))]
            tangent = _vsub(following, previous)
            norm = _vnorm(tangent)
            norm <= 1e-12 && continue
            normal = (-tangent[2] / norm, tangent[1] / norm)
            push!(nodes, _vadd(points[index], _vscale(offset, normal)))
            push!(nodes, _vsub(points[index], _vscale(offset, normal)))
        end
    end
    nodes
end

function _homotopic_visibility_path(points, punctures, barrier_routes;
                                    clearance, route_clearance, samples=32)
    length(points) <= 2 && return copy(points)
    start_point, end_point = first(points), last(points)
    radii = _visibility_radii(start_point, end_point, punctures, clearance)
    target_word = _path_cut_word(points, punctures, radii)
    nodes = NTuple{2,Float64}[start_point, end_point]
    for (puncture, radius) in zip(punctures, radii)
        vertex_radius = radius / cospi(1 / samples) + 1e-6
        for sample in 0:samples-1
            angle = 2π * sample / samples
            push!(nodes, (puncture.center[1] + vertex_radius * cos(angle),
                          puncture.center[2] + vertex_radius * sin(angle)))
        end
    end
    append!(nodes, _barrier_navigation_nodes(barrier_routes,
                                              max(0.6 * clearance, 3e-3)))
    adjacency = [Tuple{Int,Float64,Vector{Int}}[] for _ in nodes]
    for first_index in eachindex(nodes), second_index in first_index+1:length(nodes)
        _visible_disk_segment(nodes[first_index], nodes[second_index], punctures, radii) ||
            continue
        _visible_route_segment(nodes[first_index], nodes[second_index], barrier_routes,
                               route_clearance) ||
            continue
        weight = _vnorm(_vsub(nodes[second_index], nodes[first_index]))
        forward_word = _path_cut_word([nodes[first_index], nodes[second_index]], punctures, radii)
        push!(adjacency[first_index], (second_index, weight, forward_word))
        push!(adjacency[second_index], (first_index, weight, -reverse(forward_word)))
    end
    frontier = [(distance=0.0, node=1, word=Int[], path=Int[1])]
    best = Dict{Tuple{Int,String},Float64}((1, "") => 0.0)
    while !isempty(frontier)
        frontier_index = argmin(state.distance for state in frontier)
        state = popat!(frontier, frontier_index)
        state.distance > get(best, (state.node, join(state.word, ',')), Inf) + 1e-12 && continue
        state.node == 2 && state.word == target_word && return nodes[state.path]
        for (neighbor, edge_length, edge_word) in adjacency[state.node]
            word = _reduce_cut_word(vcat(state.word, edge_word))
            length(word) <= length(target_word) || continue
            word == target_word[1:length(word)] || continue
            distance = state.distance + edge_length
            key = (neighbor, join(word, ','))
            distance + 1e-12 < get(best, key, Inf) || continue
            best[key] = distance
            push!(frontier, (distance=distance, node=neighbor, word=word,
                             path=vcat(state.path, neighbor)))
        end
    end
    copy(points)
end

function _geodesically_tighten_polylines(H, layout, point_routes;
                                         clearance, route_clearance,
                                         minimum_relative_gain=0.0)
    current = deepcopy(point_routes)
    accepted = 0
    before = sum(_polyline_length, current)
    attempts = Dict{String,Any}[]
    for route_index in eachindex(current)
        points = current[route_index]
        length(points) <= 4 && continue
        interior = points[2:end-1]
        barriers = [current[index] for index in eachindex(current) if index != route_index]
        geodesic = _homotopic_visibility_path(
            interior, layout.punctures, barriers; clearance, route_clearance)
        old_length = _polyline_length(interior)
        new_length = _polyline_length(geodesic)
        route = layout.routes[route_index]
        detail = Dict{String,Any}(
            "family" => String(route.family),
            "curve_id" => route.curve_id,
            "segment_id" => route.segment_id,
            "old_length" => old_length,
            "candidate_length" => new_length,
            "status" => "not-shorter")
        relative_gain = (old_length - new_length) / max(old_length, eps())
        detail["relative_gain"] = relative_gain
        if !(new_length + 1e-3 < old_length &&
             relative_gain >= minimum_relative_gain)
            push!(attempts, detail)
            continue
        end
        proposal = deepcopy(current)
        proposal[route_index] = vcat(points[1:1], geodesic, points[end:end])
        if _validated_polyline_candidate(H, layout, proposal) === nothing
            detail["status"] = "validator-rejected"
            push!(attempts, detail)
            continue
        end
        current = proposal
        accepted += 1
        detail["status"] = "accepted"
        push!(attempts, detail)
    end
    current, (accepted=accepted, before=before, after=sum(_polyline_length, current),
              attempts=attempts)
end

function _mechanical_energy(point_routes, punctures;
                            tension, bending, route_clearance, route_repulsion,
                            puncture_clearance, puncture_repulsion,
                            density_weight, spread_range)
    energy = 0.0
    for points in point_routes
        for index in 1:length(points)-1
            delta = _vsub(points[index + 1], points[index])
            energy += tension * _vnorm(delta)
        end
        for index in 3:length(points)-2
            second = _vadd(_vsub(points[index - 1], _vscale(2.0, points[index])), points[index + 1])
            energy += bending * (second[1]^2 + second[2]^2)
        end
    end
    for first_route in eachindex(point_routes), first_segment in 1:length(point_routes[first_route])-1
        first_points = point_routes[first_route]
        for second_route in first_route:length(point_routes)
            second_points = point_routes[second_route]
            start_segment = second_route == first_route ? first_segment + 4 : 1
            for second_segment in start_segment:length(second_points)-1
                distance = _segment_distance(first_points[first_segment], first_points[first_segment + 1],
                                             second_points[second_segment], second_points[second_segment + 1])
                energy += _clearance_barrier(distance, route_clearance, route_repulsion)
            end
        end
    end
    for points in point_routes, segment in 2:length(points)-2, puncture in punctures
        margin = sqrt(_segment_distance_squared(puncture.center,
                      points[segment], points[segment + 1])) - puncture.radius
        energy += _clearance_barrier(margin, puncture_clearance, puncture_repulsion)
    end
    for route_index in eachindex(point_routes), index in 3:length(point_routes[route_index])-2
        density, _ = _multiscale_density(point_routes[route_index][index], route_index,
                                         point_routes, spread_range)
        energy += density_weight * density^3 / 3
    end
    energy
end

function _mechanical_forces(point_routes, punctures;
                            tension, bending, route_clearance, route_repulsion,
                            puncture_clearance, puncture_repulsion,
                            density_weight, spread_range)
    forces = [[(0.0, 0.0) for _ in points] for points in point_routes]
    for route_index in eachindex(point_routes)
        points = point_routes[route_index]
        for index in 2:length(points)-1
            incoming = _vsub(points[index], points[index - 1])
            outgoing = _vsub(points[index + 1], points[index])
            force = _vscale(tension, _vsub(_unit(outgoing), _unit(incoming)))
            if 3 <= index <= length(points)-2
                fourth = _vadd(
                    _vadd(points[index - 2], _vscale(-4.0, points[index - 1])),
                    _vadd(_vscale(6.0, points[index]),
                          _vadd(_vscale(-4.0, points[index + 1]), points[index + 2])))
                force = _vadd(force, _vscale(-bending, fourth))
            end
            density, density_gradient = _multiscale_density(
                points[index], route_index, point_routes, spread_range)
            force = _vadd(force,
                _vscale(-density_weight * density^2, density_gradient))
            forces[route_index][index] = force
        end
    end
    for first_route in eachindex(point_routes), first_segment in 1:length(point_routes[first_route])-1
        first_points = point_routes[first_route]
        for second_route in first_route:length(point_routes)
            second_points = point_routes[second_route]
            start_segment = second_route == first_route ? first_segment + 4 : 1
            for second_segment in start_segment:length(second_points)-1
                distance, first_point, second_point, first_parameter, second_parameter =
                    _closest_segment_pair(first_points[first_segment], first_points[first_segment + 1],
                                          second_points[second_segment], second_points[second_segment + 1])
                magnitude = _clearance_force(distance, route_clearance, route_repulsion)
                magnitude == 0.0 && continue
                difference = _vsub(first_point, second_point)
                direction = distance <= 1e-12 ?
                    _unit((-(second_points[second_segment + 1][2] - second_points[second_segment][2]),
                            second_points[second_segment + 1][1] - second_points[second_segment][1])) :
                    _vscale(inv(distance), difference)
                first_force = _vscale(magnitude, direction)
                second_force = _vscale(-magnitude, direction)
                forces[first_route][first_segment] = _vadd(
                    forces[first_route][first_segment], _vscale(1 - first_parameter, first_force))
                forces[first_route][first_segment + 1] = _vadd(
                    forces[first_route][first_segment + 1], _vscale(first_parameter, first_force))
                forces[second_route][second_segment] = _vadd(
                    forces[second_route][second_segment], _vscale(1 - second_parameter, second_force))
                forces[second_route][second_segment + 1] = _vadd(
                    forces[second_route][second_segment + 1], _vscale(second_parameter, second_force))
            end
        end
    end
    for route_index in eachindex(point_routes)
        points = point_routes[route_index]
        for segment in 2:length(points)-2, puncture in punctures
            closest, parameter = _closest_on_segment_with_parameter(
                puncture.center, points[segment], points[segment + 1])
            radial = _vsub(closest, puncture.center)
            margin = _vnorm(radial) - puncture.radius
            magnitude = _clearance_force(margin, puncture_clearance, puncture_repulsion)
            magnitude == 0.0 && continue
            force = _vscale(magnitude, _unit(radial))
            forces[route_index][segment] = _vadd(
                forces[route_index][segment], _vscale(1 - parameter, force))
            forces[route_index][segment + 1] = _vadd(
                forces[route_index][segment + 1], _vscale(parameter, force))
        end
    end
    forces
end

function _layout_from_polylines(layout, point_routes)
    routes = [Route(route.family, route.curve_id, route.segment_id,
                    route.start_port, route.end_port, point_routes[index])
              for (index, route) in enumerate(layout.routes)]
    HeegaardLayout(layout.diagram_hash, layout.view, deepcopy(layout.chart),
                   deepcopy(layout.punctures), deepcopy(layout.ports), routes,
                   Tuple{Int,Int}[], deepcopy(layout.backend_trace),
                   LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
end

function _validated_polyline_candidate(H, layout, point_routes)
    provisional = _layout_from_polylines(layout, point_routes)
    report = validate_layout(H, provisional; expected_hash=layout.diagram_hash)
    report.valid || return nothing
    HeegaardLayout(provisional.diagram_hash, provisional.view, provisional.chart,
                   provisional.punctures, provisional.ports, provisional.routes,
                   _route_crossings(provisional.routes), provisional.backend_trace, report)
end

function _resample_polyline(points, count)
    count <= 2 && return [first(points), last(points)]
    cumulative = zeros(Float64, length(points))
    for index in 2:length(points)
        cumulative[index] = cumulative[index - 1] +
                            _vnorm(_vsub(points[index], points[index - 1]))
    end
    cumulative[end] <= eps(Float64) && return fill(first(points), count)
    result = Tuple{Float64,Float64}[]
    segment = 1
    for target in range(0.0, cumulative[end]; length=count)
        while segment < length(points) - 1 && cumulative[segment + 1] < target
            segment += 1
        end
        denominator = cumulative[segment + 1] - cumulative[segment]
        fraction = denominator <= eps(Float64) ? 0.0 :
                   (target - cumulative[segment]) / denominator
        push!(result, _vadd(points[segment],
                            _vscale(fraction, _vsub(points[segment + 1], points[segment]))))
    end
    result
end

function _redistribute_route(points)
    length(points) <= 4 && return points
    middle = _resample_polyline(points[2:end-1], length(points) - 2)
    vcat(points[1:1], middle, points[end:end])
end

function _relax_polylines(H, layout, point_routes;
                          iterations, time_step, tension, bending,
                          route_clearance, route_repulsion,
                          puncture_clearance, puncture_repulsion,
                          density_weight, spread_range,
                          safety_fraction, convergence_tolerance)
    current = deepcopy(point_routes)
    energy = _mechanical_energy(current, layout.punctures;
        tension, bending, route_clearance, route_repulsion,
        puncture_clearance, puncture_repulsion, density_weight, spread_range)
    history = Float64[energy]
    accepted_steps = 0
    converged = false
    minimum_clearance = _polyline_clearance(current, layout.punctures)
    for _ in 1:iterations
        forces = _mechanical_forces(current, layout.punctures;
            tension, bending, route_clearance, route_repulsion,
            puncture_clearance, puncture_repulsion, density_weight, spread_range)
        force_norm = maximum((_vnorm(force) for route in forces for force in route); init=0.0)
        force_norm <= convergence_tolerance && (converged = true; break)
        clearance = _polyline_clearance(current, layout.punctures)
        minimum_clearance = min(minimum_clearance, clearance)
        safe_step = min(time_step, safety_fraction * clearance / max(force_norm, eps(Float64)))
        safe_step <= eps(Float64) && break
        accepted = false
        trial_step = safe_step
        for _ in 1:20
            proposal = deepcopy(current)
            for route_index in eachindex(proposal), index in 3:length(proposal[route_index])-2
                proposal[route_index][index] = _vadd(
                    current[route_index][index], _vscale(trial_step, forces[route_index][index]))
            end
            proposal = [_redistribute_route(points) for points in proposal]
            candidate = _validated_polyline_candidate(H, layout, proposal)
            proposal_energy = candidate === nothing ? Inf :
                _mechanical_energy(proposal, layout.punctures;
                    tension, bending, route_clearance, route_repulsion,
                    puncture_clearance, puncture_repulsion,
                    density_weight, spread_range)
            if candidate !== nothing && proposal_energy <= energy + 1e-12
                current, energy = proposal, proposal_energy
                push!(history, energy)
                accepted_steps += 1
                accepted = true
                maximum(trial_step * _vnorm(force) for route in forces for force in route) <=
                    convergence_tolerance && (converged = true)
                break
            end
            trial_step /= 2
        end
        (!accepted || converged) && break
    end
    current, (history=history, accepted_steps=accepted_steps,
              converged=converged, minimum_clearance=minimum_clearance)
end

function _point_line_distance(point, first, second)
    sqrt(_segment_distance_squared(point, first, second))
end

function _rdp_indices(points, tolerance)
    length(points) <= 2 && return collect(eachindex(points))
    distances = [_point_line_distance(points[index], first(points), last(points))
                 for index in 2:length(points)-1]
    isempty(distances) && return [1, length(points)]
    maximum_distance, offset = findmax(distances)
    maximum_distance <= tolerance && return [1, length(points)]
    split = offset + 1
    left = _rdp_indices(points[1:split], tolerance)
    right = _rdp_indices(points[split:end], tolerance)
    vcat(left[1:end-1], [split - 1 + index for index in right])
end

function _simplify_polylines(H, layout, point_routes; tolerance)
    original_count = sum(length, point_routes)
    trial_tolerance = tolerance
    for _ in 1:20
        proposal = [begin
            if length(points) <= 4
                copy(points)
            else
                interior = points[2:end-1]
                selected = _rdp_indices(interior, trial_tolerance)
                vcat(points[1:1], interior[selected], points[end:end])
            end
        end for points in point_routes]
        candidate = _validated_polyline_candidate(H, layout, proposal)
        candidate !== nothing &&
            return proposal, original_count - sum(length, proposal)
        trial_tolerance /= 2
    end
    deepcopy(point_routes), 0
end

function _vertex_tangent(points, index)
    index == 1 && return _unit(_vsub(points[2], points[1]))
    index == length(points) && return _unit(_vsub(points[end], points[end - 1]))
    incoming = _unit(_vsub(points[index], points[index - 1]))
    outgoing = _unit(_vsub(points[index + 1], points[index]))
    tangent = _unit(_vadd(incoming, outgoing))
    _vnorm(tangent) <= eps(Float64) ? outgoing : tangent
end

function _bezier_route(route, points; handle_scale, flatness, max_depth)
    tangents = [_vertex_tangent(points, index) for index in eachindex(points)]
    if length(points) >= 4
        tangents[2] = _unit(_vsub(points[2], points[1]))
        tangents[end - 1] = _unit(_vsub(points[end], points[end - 1]))
    end
    handles = zeros(Float64, length(points))
    for index in 2:length(points)-1
        handles[index] = handle_scale * min(_vnorm(_vsub(points[index], points[index - 1])),
                                             _vnorm(_vsub(points[index + 1], points[index])))
    end
    segments = CubicBezierSegment[]
    for index in 1:length(points)-1
        first, second = points[index], points[index + 1]
        if index == 1 || index == length(points)-1
            delta = _vsub(second, first)
            push!(segments, CubicBezierSegment((first, _vadd(first, _vscale(1 / 3, delta)),
                                                _vadd(first, _vscale(2 / 3, delta)), second)))
        else
            push!(segments, CubicBezierSegment((
                first, _vadd(first, _vscale(handles[index], tangents[index])),
                _vsub(second, _vscale(handles[index + 1], tangents[index + 1])), second)))
        end
    end
    curve = BezierChain(segments; flatness, max_depth)
    Route(route.family, route.curve_id, route.segment_id,
          route.start_port, route.end_port, _curve_points(curve); curve)
end

function _fit_bezier_layout(H, layout, point_routes;
                            handle_scale, flatness, max_depth,
                            minimum_route_separation)
    scale = handle_scale
    diagnostic = ""
    for _ in 1:20
        routes = [_bezier_route(route, point_routes[index];
                    handle_scale=scale, flatness, max_depth)
                  for (index, route) in enumerate(layout.routes)]
        provisional = HeegaardLayout(layout.diagram_hash, layout.view, deepcopy(layout.chart),
            deepcopy(layout.punctures), deepcopy(layout.ports), routes,
            Tuple{Int,Int}[], deepcopy(layout.backend_trace),
            LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
        report = validate_layout(H, provisional; expected_hash=layout.diagram_hash)
        endpoint_tangents = all(begin
            chain = route.curve
            start_port = only(port for port in layout.ports if port.id == route.start_port)
            end_port = only(port for port in layout.ports if port.id == route.end_port)
            start_puncture = only(puncture for puncture in layout.punctures
                                  if puncture.id == start_port.boundary_id)
            end_puncture = only(puncture for puncture in layout.punctures
                                if puncture.id == end_port.boundary_id)
            start_normal = _vsub(start_port.point, start_puncture.center)
            end_normal = _vsub(end_port.point, end_puncture.center)
            start_tangent = _vsub(first(chain.segments).p1, first(chain.segments).p0)
            reverse_end_tangent = _vsub(last(chain.segments).p2, last(chain.segments).p3)
            start_normal[1] * start_tangent[1] + start_normal[2] * start_tangent[2] >= -1e-12 &&
            end_normal[1] * reverse_end_tangent[1] + end_normal[2] * reverse_end_tangent[2] >= -1e-12
        end for route in routes)
        route_separation = _route_pair_clearance([route.points for route in routes])
        diagnostic = "issues=$(join(string.(getfield.(report.issues,:code)),',')); endpoint_tangents=$endpoint_tangents; separation=$route_separation; required=$minimum_route_separation"
        if report.valid && endpoint_tangents &&
           route_separation >= minimum_route_separation
            accepted = HeegaardLayout(provisional.diagram_hash, provisional.view, provisional.chart,
                provisional.punctures, provisional.ports, routes, _route_crossings(routes),
                provisional.backend_trace, report)
            minimum_puncture_margin = Inf
            for route in routes
                points = _route_points(route)
                for puncture in layout.punctures, index in 1:length(points)-1
                    minimum_puncture_margin = min(minimum_puncture_margin,
                        sqrt(_segment_distance_squared(puncture.center,
                            points[index], points[index + 1])) - puncture.radius)
                end
            end
            return accepted, scale, (
                minimum_separation=route_separation,
                minimum_puncture_margin=minimum_puncture_margin,
                tolerance=flatness,
                endpoint_tangents=true)
        end
        scale /= 2
    end
    throw(ArgumentError("mechanical relaxation could not fit a validated Bezier chain: " * diagnostic))
end

"""Relax a verified PL layout by constrained tension/bending/repulsion dynamics.

The exact endpoint guard edges remain fixed.  Every accepted time step is
bounded by the current route/puncture clearance and is revalidated.  The
relaxed dense polyline is simplified by a fail-closed shortcut proposal and
then fitted by a G1 cubic Bezier chain. The final adaptive flattening is passed
through the ordinary layout validator before the result is returned.
"""
function mechanically_relax_layout(H::HeegaardDiagram, layout::HeegaardLayout;
        max_spacing=0.10, iterations=120, time_step=0.04,
        tension=1.0, bending=0.05,
        route_clearance=0.22, route_repulsion=8.0,
        minimum_route_separation=0.004,
        puncture_clearance=0.18, puncture_repulsion=10.0,
        density_weight=0.08, spread_range=0.45,
        repulsion_distance=nothing, repulsion_strength=nothing,
        safety_fraction=0.12,
        convergence_tolerance=1e-5, simplify_tolerance=0.04,
        handle_scale=1 / 3, flatness=1e-4, max_depth=20)
    layout.view == :alpha_cut || throw(ArgumentError("mechanical relaxation requires alpha-cut view"))
    validate_layout(H, layout; expected_hash=layout.diagram_hash).valid ||
        throw(ArgumentError("mechanical relaxation requires a verified feasible layout"))
    all(route.curve isa PolylineCurve for route in layout.routes) ||
        throw(ArgumentError("mechanical relaxation requires PL input routes"))
    max_spacing > 0 || throw(ArgumentError("max_spacing must be positive"))
    iterations >= 1 || throw(ArgumentError("iterations must be positive"))
    0 < safety_fraction < 0.5 || throw(ArgumentError("safety_fraction must be in (0,1/2)"))
    flatness > 0 || throw(ArgumentError("flatness must be positive"))
    repulsion_distance === nothing ||
        (route_clearance = puncture_clearance = Float64(repulsion_distance))
    repulsion_strength === nothing ||
        (route_repulsion = puncture_repulsion = Float64(repulsion_strength))
    route_clearance > 0 || throw(ArgumentError("route_clearance must be positive"))
    minimum_route_separation > 0 ||
        throw(ArgumentError("minimum_route_separation must be positive"))
    puncture_clearance > 0 || throw(ArgumentError("puncture_clearance must be positive"))
    route_repulsion > 0 || throw(ArgumentError("route_repulsion must be positive"))
    puncture_repulsion > 0 || throw(ArgumentError("puncture_repulsion must be positive"))
    density_weight >= 0 || throw(ArgumentError("density_weight must be nonnegative"))
    spread_range > 0 || throw(ArgumentError("spread_range must be positive"))

    source_points = [_route_points(route) for route in layout.routes]
    shortened, shortcuts = _global_shortcut_polylines(H, layout, source_points)
    geodesic, geodesics = _geodesically_tighten_polylines(
        H, layout, shortened; clearance=puncture_clearance,
        route_clearance=minimum_route_separation)
    dense = [_densify_route(points; max_spacing) for points in geodesic]
    relaxed, relaxation = _relax_polylines(H, layout, dense;
        iterations, time_step, tension, bending, route_clearance,
        route_repulsion, puncture_clearance, puncture_repulsion,
        density_weight, spread_range,
        safety_fraction, convergence_tolerance)
    presimplified, first_removed = _simplify_polylines(
        H, layout, relaxed; tolerance=simplify_tolerance)
    tightened, final_geodesics = _geodesically_tighten_polylines(
        H, layout, presimplified; clearance=puncture_clearance,
        route_clearance=minimum_route_separation,
        minimum_relative_gain=0.05)
    simplified, second_removed = _simplify_polylines(
        H, layout, tightened; tolerance=simplify_tolerance)
    removed = first_removed + second_removed
    result, accepted_handle_scale, certificate = _fit_bezier_layout(H, layout, simplified;
        handle_scale, flatness, max_depth, minimum_route_separation)

    trace = deepcopy(result.backend_trace)
    trace["smoothing"] = Dict{String,Any}(
        "method" => "global-shortcut-mechanical-relaxation-v2",
        "shortcut_count" => shortcuts.accepted,
        "geodesic_route_count" => geodesics.accepted,
        "length_before_shortcut" => shortcuts.before,
        "length_after_shortcut" => shortcuts.after,
        "length_after_geodesic" => geodesics.after,
        "geodesic_attempts" => geodesics.attempts,
        "length_after_relaxation" => sum(_polyline_length, relaxed),
        "final_geodesic_route_count" => final_geodesics.accepted,
        "length_after_final_geodesic" => final_geodesics.after,
        "final_geodesic_attempts" => final_geodesics.attempts,
        "dense_points" => sum(length, dense),
        "relaxed_points" => sum(length, relaxed),
        "simplified_points" => sum(length, simplified),
        "removed_points" => removed,
        "bezier_segments" => sum(length(route.curve.segments) for route in result.routes),
        "control_points" => sum(3length(route.curve.segments) + 1 for route in result.routes),
        "accepted_steps" => relaxation.accepted_steps,
        "converged" => relaxation.converged,
        "energy_history" => relaxation.history,
        "initial_energy" => first(relaxation.history),
        "final_energy" => last(relaxation.history),
        "minimum_step_clearance" => relaxation.minimum_clearance,
        "minimum_flattened_separation" => certificate.minimum_separation,
        "minimum_flattened_puncture_margin" => certificate.minimum_puncture_margin,
        "flatness_bound" => certificate.tolerance,
        "endpoint_tangents_certified" => certificate.endpoint_tangents,
        "handle_scale" => accepted_handle_scale,
        "parameters" => Dict{String,Any}(
            "max_spacing" => max_spacing, "iterations" => iterations,
            "time_step" => time_step, "tension" => tension, "bending" => bending,
            "route_clearance" => route_clearance,
            "minimum_route_separation" => minimum_route_separation,
            "route_repulsion" => route_repulsion,
            "puncture_clearance" => puncture_clearance,
            "puncture_repulsion" => puncture_repulsion,
            "density_weight" => density_weight,
            "spread_range" => spread_range,
            "safety_fraction" => safety_fraction,
            "convergence_tolerance" => convergence_tolerance,
            "simplify_tolerance" => simplify_tolerance))
    haskey(trace, "pipeline") &&
        (trace["pipeline"]["smoothing"] = "global-shortcut-mechanical-relaxation-v2")
    HeegaardLayout(result.diagram_hash, result.view, result.chart,
        result.punctures, result.ports, result.routes, result.crossings,
        trace, result.report)
end

function smooth_layout(H::HeegaardDiagram, layout::HeegaardLayout;
                       step_fraction=nothing, clearance_fraction=nothing,
                       puncture_fraction=nothing, kwargs...)
    mechanically_relax_layout(H, layout; kwargs...)
end
