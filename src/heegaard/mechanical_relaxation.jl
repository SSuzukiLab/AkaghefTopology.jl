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
    for _ in 1:8
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

function _mechanical_energy(point_routes, punctures;
                            tension, bending, repulsion_distance, repulsion_strength)
    energy = 0.0
    for points in point_routes
        for index in 1:length(points)-1
            delta = _vsub(points[index + 1], points[index])
            energy += tension * _vnorm(delta)
        end
        for index in 2:length(points)-1
            second = _vadd(_vsub(points[index - 1], _vscale(2.0, points[index])), points[index + 1])
            energy += bending * (second[1]^2 + second[2]^2)
        end
    end
    for route_index in eachindex(point_routes), point_index in 3:length(point_routes[route_index])-2
        point = point_routes[route_index][point_index]
        for other_route in eachindex(point_routes)
            other = point_routes[other_route]
            for segment in 1:length(other)-1
                other_route == route_index && abs(segment - point_index) <= 2 && continue
                distance = sqrt(_segment_distance_squared(point, other[segment], other[segment + 1]))
                distance < repulsion_distance &&
                    (energy += repulsion_strength * (repulsion_distance - distance)^2)
            end
        end
        for puncture in punctures
            margin = _vnorm(_vsub(point, puncture.center)) - puncture.radius
            margin < repulsion_distance &&
                (energy += repulsion_strength * (repulsion_distance - margin)^2)
        end
    end
    energy
end

function _mechanical_forces(point_routes, punctures;
                            tension, bending, repulsion_distance, repulsion_strength)
    forces = [[(0.0, 0.0) for _ in points] for points in point_routes]
    for route_index in eachindex(point_routes)
        points = point_routes[route_index]
        for index in 3:length(points)-2
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
            point = points[index]
            for other_route in eachindex(point_routes)
                other = point_routes[other_route]
                for segment in 1:length(other)-1
                    other_route == route_index && abs(segment - index) <= 2 && continue
                    closest = _closest_on_segment(point, other[segment], other[segment + 1])
                    difference = _vsub(point, closest)
                    distance = _vnorm(difference)
                    if distance < repulsion_distance
                        direction = distance <= 1e-12 ?
                            _unit((- (other[segment + 1][2] - other[segment][2]),
                                    other[segment + 1][1] - other[segment][1])) :
                            _vscale(inv(distance), difference)
                        magnitude = repulsion_strength *
                                    (repulsion_distance - max(distance, 1e-12)) /
                                    repulsion_distance
                        force = _vadd(force, _vscale(magnitude, direction))
                    end
                end
            end
            for puncture in punctures
                radial = _vsub(point, puncture.center)
                distance = _vnorm(radial)
                margin = distance - puncture.radius
                if margin < repulsion_distance
                    magnitude = repulsion_strength *
                                (repulsion_distance - margin) / repulsion_distance
                    force = _vadd(force, _vscale(magnitude, _unit(radial)))
                end
            end
            forces[route_index][index] = force
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
                          repulsion_distance, repulsion_strength,
                          safety_fraction, convergence_tolerance)
    current = deepcopy(point_routes)
    energy = _mechanical_energy(current, layout.punctures;
        tension, bending, repulsion_distance, repulsion_strength)
    history = Float64[energy]
    accepted_steps = 0
    converged = false
    minimum_clearance = _polyline_clearance(current, layout.punctures)
    for _ in 1:iterations
        forces = _mechanical_forces(current, layout.punctures;
            tension, bending, repulsion_distance, repulsion_strength)
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
                    tension, bending, repulsion_distance, repulsion_strength)
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
                            handle_scale, flatness, max_depth)
    scale = handle_scale
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
        if report.valid && endpoint_tangents
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
                minimum_separation=_polyline_clearance(
                    [route.points for route in routes], layout.punctures),
                minimum_puncture_margin=minimum_puncture_margin,
                tolerance=flatness,
                endpoint_tangents=true)
        end
        scale /= 2
    end
    throw(ArgumentError("mechanical relaxation could not fit a validated Bezier chain"))
end

"""Relax a verified PL layout by constrained tension/bending/repulsion dynamics.

The exact endpoint guard edges remain fixed.  Every accepted time step is
bounded by the current route/puncture clearance and is revalidated.  The
relaxed dense polyline is simplified by a fail-closed shortcut proposal and
then fitted by a G1 cubic Bezier chain. The final adaptive flattening is passed
through the ordinary layout validator before the result is returned.
"""
function mechanically_relax_layout(H::HeegaardDiagram, layout::HeegaardLayout;
        max_spacing=0.12, iterations=60, time_step=0.08,
        tension=1.0, bending=0.05, repulsion_distance=0.16,
        repulsion_strength=0.08, safety_fraction=0.18,
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

    source_points = [_route_points(route) for route in layout.routes]
    shortened, shortcuts = _global_shortcut_polylines(H, layout, source_points)
    dense = [_densify_route(points; max_spacing) for points in shortened]
    relaxed, relaxation = _relax_polylines(H, layout, dense;
        iterations, time_step, tension, bending, repulsion_distance,
        repulsion_strength, safety_fraction, convergence_tolerance)
    simplified, removed = _simplify_polylines(H, layout, relaxed;
                                               tolerance=simplify_tolerance)
    result, accepted_handle_scale, certificate = _fit_bezier_layout(H, layout, simplified;
        handle_scale, flatness, max_depth)

    trace = deepcopy(result.backend_trace)
    trace["smoothing"] = Dict{String,Any}(
        "method" => "global-shortcut-mechanical-relaxation-v2",
        "shortcut_count" => shortcuts.accepted,
        "length_before_shortcut" => shortcuts.before,
        "length_after_shortcut" => shortcuts.after,
        "length_after_relaxation" => sum(_polyline_length, relaxed),
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
            "repulsion_distance" => repulsion_distance,
            "repulsion_strength" => repulsion_strength,
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
