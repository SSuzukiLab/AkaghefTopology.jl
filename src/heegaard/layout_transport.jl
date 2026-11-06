"""Ambient-isotopy transport and validation-preserving route smoothing.

The terminal ribbon map is first drawn without prescribed vertex coordinates.
Puncture disks are then carried to the requested chart positions by a finite
composition of compactly supported bi-Lipschitz maps.  A translation step has

    F(x) = x + v b(norm(x-c)),

where `b` is one on an inner disk and decreases linearly to zero on the support
annulus.  The step condition `norm(v) < support_radius - inner_radius` makes
`F_t = id + t(F-id)` bi-Lipschitz for every `0 <= t <= 1`.  Radial changes of a
puncture radius use a strictly increasing piecewise-linear radial coordinate.
Together these maps are an orientation-preserving ambient isotopy of the
one-point compactified alpha-cut chart.

Routes are split at the affine-piece boundaries of the square bump/radial
maps, so each segment has its exact piecewise-affine image. Every stage is
accepted only after the full layout validator succeeds. The trace records
every accepted isotopy stage and its Lipschitz/Jacobian bound.
"""

_point_distance(a, b) = hypot(a[1] - b[1], a[2] - b[2])

function _map_layout_geometry(punctures, ports, routes, point_map;
                              radius_updates=Dict{Int,Float64}(),
                              mapped_routes=nothing)
    mapped_punctures = [Puncture(
        puncture.id, puncture.pair, puncture.side, point_map(puncture.center),
        get(radius_updates, puncture.id, puncture.radius), puncture.orientation,
    ) for puncture in punctures]
    mapped_ports = [Port(
        port.id, port.intersection_id, port.boundary_id, port.beta_id,
        port.parameter, point_map(port.point),
    ) for port in ports]
    result_routes = mapped_routes === nothing ? [Route(
        route.family, route.curve_id, route.segment_id,
        route.start_port, route.end_port, point_map.(route.points),
    ) for route in routes] : mapped_routes
    mapped_punctures, mapped_ports, result_routes
end


function _compress_collinear_points(points)
    result = NTuple{2,Float64}[]
    for point in points
        (isempty(result) || !_same_point(last(result), point)) && push!(result, point)
        while length(result) >= 3
            first_point, middle, last_point = result[end - 2], result[end - 1], result[end]
            scale = max(_point_distance(first_point, middle) *
                        _point_distance(middle, last_point), 1.0)
            abs(_segment_orientation(first_point, middle, last_point)) <= 1e-12scale || break
            _between(first_point, last_point, middle) || break
            result[end - 1] = last_point
            pop!(result)
        end
    end
    result
end

function _map_routes_on_line_arrangement(routes, point_map, lines;
                                         split_predicate=point -> true)
    [begin
        points = NTuple{2,Float64}[]
        for index in 1:length(route.points)-1
            first_point, second_point = route.points[index], route.points[index + 1]
            parameters = Float64[0.0, 1.0]
            for (a, b, c) in lines
                first_value = a * first_point[1] + b * first_point[2] + c
                second_value = a * second_point[1] + b * second_point[2] + c
                denominator = second_value - first_value
                abs(denominator) <= 1e-14 && continue
                parameter = -first_value / denominator
                if 1e-12 < parameter < 1 - 1e-12
                    crossing = (
                        first_point[1] + parameter * (second_point[1] - first_point[1]),
                        first_point[2] + parameter * (second_point[2] - first_point[2]),
                    )
                    split_predicate(crossing) && push!(parameters, parameter)
                end
            end
            sort!(unique!(parameters))
            source_piece = [(
                first_point[1] + parameter * (second_point[1] - first_point[1]),
                first_point[2] + parameter * (second_point[2] - first_point[2]),
            ) for parameter in parameters]
            mapped_piece = point_map.(source_piece)
            isempty(points) ? append!(points, mapped_piece) : append!(points, mapped_piece[2:end])
        end
        filtered = _compress_collinear_points(points)
        Route(route.family, route.curve_id, route.segment_id,
              route.start_port, route.end_port, filtered)
    end for route in routes]
end

function _piecewise_affine_bump_routes(routes, point_map, center,
                                       inner_radius, support_radius)
    # The L-infinity bump changes affine formula only on these lines.  After
    # splitting at every crossing, endpoint mapping gives the exact image of
    # each route segment, rather than a chordal approximation.
    lines = NTuple{3,Float64}[
        (1.0, 0.0, -(center[1] - inner_radius)),
        (1.0, 0.0, -(center[1] + inner_radius)),
        (0.0, 1.0, -(center[2] - inner_radius)),
        (0.0, 1.0, -(center[2] + inner_radius)),
        (1.0, 0.0, -(center[1] - support_radius)),
        (1.0, 0.0, -(center[1] + support_radius)),
        (0.0, 1.0, -(center[2] - support_radius)),
        (0.0, 1.0, -(center[2] + support_radius)),
        (1.0, -1.0, -center[1] + center[2]),
        (1.0, 1.0, -center[1] - center[2]),
    ]
    inside_support = point -> max(abs(point[1] - center[1]),
                                  abs(point[2] - center[2])) <= support_radius + 1e-9
    _map_routes_on_line_arrangement(
        routes, point_map, lines; split_predicate=inside_support)
end

function _line_through(first_point, second_point)
    a = first_point[2] - second_point[2]
    b = second_point[1] - first_point[1]
    c = -(a * first_point[1] + b * first_point[2])
    scale = hypot(a, b)
    scale > 1e-14 || throw(ArgumentError("cannot form a line through coincident points"))
    (a / scale, b / scale, c / scale)
end

function _triangle_orientation(triangle)
    _segment_orientation(triangle[1], triangle[2], triangle[3])
end

function _pl_radial_square_map(center, old_radius, new_radius, support_radius)
    inner = [(-old_radius, old_radius), (old_radius, old_radius),
             (old_radius, -old_radius), (-old_radius, -old_radius)]
    outer = [(-support_radius, support_radius), (support_radius, support_radius),
             (support_radius, -support_radius), (-support_radius, -support_radius)]
    scale = new_radius / old_radius
    mapped_inner = [(scale * point[1], scale * point[2]) for point in inner]
    quadrilaterals = [
        (inner[1], inner[2], outer[2], outer[1]),
        (inner[2], inner[3], outer[3], outer[2]),
        (inner[3], inner[4], outer[4], outer[3]),
        (inner[4], inner[1], outer[1], outer[4]),
    ]
    mapped_quadrilaterals = [
        (mapped_inner[1], mapped_inner[2], outer[2], outer[1]),
        (mapped_inner[2], mapped_inner[3], outer[3], outer[2]),
        (mapped_inner[3], mapped_inner[4], outer[4], outer[3]),
        (mapped_inner[4], mapped_inner[1], outer[1], outer[4]),
    ]
    triangles = NTuple{3,NTuple{2,Float64}}[]
    mapped_triangles = NTuple{3,NTuple{2,Float64}}[]
    for (quad, mapped_quad) in zip(quadrilaterals, mapped_quadrilaterals)
        for indices in ((1, 2, 3), (1, 3, 4))
            triangle = (quad[indices[1]], quad[indices[2]], quad[indices[3]])
            mapped_triangle = (mapped_quad[indices[1]], mapped_quad[indices[2]],
                               mapped_quad[indices[3]])
            if _triangle_orientation(triangle) < 0
                triangle = (triangle[1], triangle[3], triangle[2])
                mapped_triangle = (mapped_triangle[1], mapped_triangle[3], mapped_triangle[2])
            end
            push!(triangles, triangle)
            push!(mapped_triangles, mapped_triangle)
        end
    end
    determinant_ratios = [_triangle_orientation(mapped_triangle) /
                          _triangle_orientation(triangle)
                          for (triangle, mapped_triangle) in zip(triangles, mapped_triangles)]
    minimum(determinant_ratios) > 0 || throw(ArgumentError(
        "radial square-annulus triangulation reverses orientation"))
    function local_map(local_point)
        max(abs(local_point[1]), abs(local_point[2])) <= old_radius + 1e-12 &&
            return (scale * local_point[1], scale * local_point[2])
        max(abs(local_point[1]), abs(local_point[2])) >= support_radius - 1e-12 &&
            return local_point
        for (triangle, mapped_triangle) in zip(triangles, mapped_triangles)
            denominator = _triangle_orientation(triangle)
            first_weight = _segment_orientation(triangle[2], triangle[3], local_point) / denominator
            second_weight = _segment_orientation(triangle[3], triangle[1], local_point) / denominator
            third_weight = 1 - first_weight - second_weight
            min(first_weight, second_weight, third_weight) >= -1e-9 || continue
            return (
                first_weight * mapped_triangle[1][1] +
                    second_weight * mapped_triangle[2][1] +
                    third_weight * mapped_triangle[3][1],
                first_weight * mapped_triangle[1][2] +
                    second_weight * mapped_triangle[2][2] +
                    third_weight * mapped_triangle[3][2],
            )
        end
        throw(ArgumentError("point in radial support was not covered by its annulus triangulation"))
    end
    point_map = point -> begin
        local_point = (point[1] - center[1], point[2] - center[2])
        mapped = local_map(local_point)
        (center[1] + mapped[1], center[2] + mapped[2])
    end
    lines = NTuple{3,Float64}[]
    for triangle in triangles
        for (first_point, second_point) in ((triangle[1], triangle[2]),
                                            (triangle[2], triangle[3]),
                                            (triangle[3], triangle[1]))
            global_first = (center[1] + first_point[1], center[2] + first_point[2])
            global_second = (center[1] + second_point[1], center[2] + second_point[2])
            push!(lines, _line_through(global_first, global_second))
        end
    end
    point_map, unique(lines), minimum(determinant_ratios)
end


function _require_valid_geometry(H, diagram_hash, punctures, ports, routes, stage)
    report = _layout_candidate_report(H, diagram_hash, punctures, ports, routes)
    report.valid || throw(ArgumentError(
        "$stage failed strict layout validation: $(_layout_issue_summary(report))"))
    report
end

function _radial_transport_step(H, diagram_hash, punctures, ports, routes,
                                puncture_id, new_radius, stage)
    index = findfirst(puncture -> puncture.id == puncture_id, punctures)
    index === nothing && throw(ArgumentError("unknown puncture $puncture_id"))
    selected = punctures[index]
    old_radius = selected.radius
    isapprox(old_radius, new_radius; atol=1e-12) &&
        return punctures, ports, routes, Dict{String,Any}(
            "stage" => stage, "puncture_id" => puncture_id,
            "old_radius" => old_radius, "new_radius" => new_radius,
            "status" => "identity",
        )
    new_radius > 0 || throw(ArgumentError("puncture radius must remain positive"))
    other_clearance = minimum((
        max(abs(selected.center[1] - other.center[1]),
            abs(selected.center[2] - other.center[2])) - other.radius
        for other in punctures if other.id != puncture_id
    ); init=Inf)
    if !isfinite(other_clearance)
        other_clearance = max(old_radius, new_radius) + 2.0
    end
    support_radius = 0.8other_clearance
    support_radius > max(old_radius, new_radius) + 1e-8 || throw(ArgumentError(
        "$stage cannot change puncture $puncture_id from radius $old_radius to $new_radius: " *
        "available radial support is $support_radius"))
    center = selected.center
    point_map, arrangement_lines, minimum_annulus_determinant =
        _pl_radial_square_map(center, old_radius, new_radius, support_radius)
    inside_support = point -> max(abs(point[1] - center[1]),
                                  abs(point[2] - center[2])) <= support_radius + 1e-9
    mapped_routes = _map_routes_on_line_arrangement(
        routes, point_map, arrangement_lines; split_predicate=inside_support)
    mapped = _map_layout_geometry(
        punctures, ports, routes, point_map;
        radius_updates=Dict(puncture_id => Float64(new_radius)),
        mapped_routes=mapped_routes,
    )
    report = _require_valid_geometry(H, diagram_hash, mapped..., stage)
    minimum_jacobian = min((new_radius / old_radius)^2,
                           minimum_annulus_determinant, 1.0)
    mapped..., Dict{String,Any}(
        "stage" => stage,
        "puncture_id" => puncture_id,
        "old_radius" => old_radius,
        "new_radius" => new_radius,
        "support_radius" => support_radius,
        "minimum_jacobian_determinant" => minimum_jacobian,
        "annulus_triangles" => 8,
        "route_image" => "exact-piecewise-affine",
        "validation_checks" => Dict(string(key) => value for (key, value) in report.checks),
        "status" => "accepted",
    )
end

function _segment_avoids_disks(first_point, second_point, obstacles)
    all(sqrt(_segment_distance_squared(obstacle.center, first_point, second_point)) >=
        obstacle.radius - 1e-9 for obstacle in obstacles)
end

function _visibility_path(start, goal, obstacles)
    _segment_avoids_disks(start, goal, obstacles) && return [start, goal]
    sides = 24
    nodes = NTuple{2,Float64}[start, goal]
    for obstacle in obstacles
        waypoint_radius = obstacle.radius / cos(π / sides) + 0.025
        for index in 0:sides-1
            angle = 2π * index / sides
            push!(nodes, (obstacle.center[1] + waypoint_radius * cos(angle),
                          obstacle.center[2] + waypoint_radius * sin(angle)))
        end
    end
    node_count = length(nodes)
    adjacency = [Tuple{Int,Float64}[] for _ in 1:node_count]
    for first_index in 1:node_count-1, second_index in first_index+1:node_count
        _segment_avoids_disks(nodes[first_index], nodes[second_index], obstacles) || continue
        distance = _point_distance(nodes[first_index], nodes[second_index])
        distance <= 1e-12 && continue
        push!(adjacency[first_index], (second_index, distance))
        push!(adjacency[second_index], (first_index, distance))
    end
    distances = fill(Inf, node_count)
    predecessors = fill(0, node_count)
    visited = falses(node_count)
    distances[1] = 0.0
    for _ in 1:node_count
        current = 0
        best = Inf
        for node in 1:node_count
            !visited[node] && distances[node] < best || continue
            current, best = node, distances[node]
        end
        current == 0 && break
        current == 2 && break
        visited[current] = true
        for (neighbor, weight) in adjacency[current]
            candidate = distances[current] + weight
            candidate + 1e-12 < distances[neighbor] || continue
            distances[neighbor] = candidate
            predecessors[neighbor] = current
        end
    end
    isfinite(distances[2]) || throw(ArgumentError(
        "no verified puncture-center path exists through the disk complement"))
    indices = Int[2]
    while first(indices) != 1
        predecessor = predecessors[first(indices)]
        predecessor != 0 || throw(ArgumentError("visibility path reconstruction failed"))
        pushfirst!(indices, predecessor)
    end
    nodes[indices]
end

function _translation_step_map(center, displacement, inner_radius, support_radius)
    gap = support_radius - inner_radius
    lipschitz = hypot(displacement...) / gap
    lipschitz < 1 || throw(ArgumentError(
        "translation bump is not injective: Lipschitz bound $lipschitz is not below one"))
    function point_map(point)
        radius = max(abs(point[1] - center[1]), abs(point[2] - center[2]))
        weight = radius <= inner_radius ? 1.0 :
                 radius < support_radius ? (support_radius - radius) / gap : 0.0
        (point[1] + weight * displacement[1],
         point[2] + weight * displacement[2])
    end
    point_map, lipschitz
end

function _move_puncture_along_path(H, diagram_hash, punctures, ports, routes,
                                   puncture_id, goal; phase)
    selected_index = findfirst(puncture -> puncture.id == puncture_id, punctures)
    selected_index === nothing && throw(ArgumentError("unknown puncture $puncture_id"))
    selected = punctures[selected_index]
    nearest = minimum((_point_distance(selected.center, other.center)-other.radius-selected.radius
                       for other in punctures if other.id!=puncture_id);init=1.0)
    path_margin = min(max(0.12, 2selected.radius), 0.25nearest)
    obstacles = [(center=other.center,
                  radius=other.radius + selected.radius + path_margin)
                 for other in punctures if other.id != puncture_id]
    path = _visibility_path(selected.center, goal, obstacles)
    stage_records = Dict{String,Any}[]
    for waypoint_index in 2:length(path)
        waypoint = path[waypoint_index]
        while _point_distance(punctures[selected_index].center, waypoint) > 1e-9
            selected = punctures[selected_index]
            center = selected.center
            other_clearance = minimum((
                max(abs(center[1] - other.center[1]),
                    abs(center[2] - other.center[2])) - other.radius
                for other in punctures if other.id != puncture_id
            ); init=Inf)
            isfinite(other_clearance) || (other_clearance = _point_distance(center, waypoint) + 2.0)
            free_gap = other_clearance - selected.radius
            free_gap > 1e-8 || throw(ArgumentError(
                "$phase puncture $puncture_id has no positive isotopy support clearance"))
            # Keep the map constant on a collar outside the puncture.  This
            # makes the PL chord from a port to the first affected route
            # vertex stay outside the moving disk, rather than merely being
            # an arbitrarily close approximation to the curved image.
            inner_radius = selected.radius + 0.35free_gap
            support_radius = selected.radius + 0.82free_gap
            available = support_radius - inner_radius
            remaining = _point_distance(center, waypoint)
            direction = ((waypoint[1] - center[1]) / remaining,
                         (waypoint[2] - center[2]) / remaining)
            step_length = min(remaining, 0.72available)
            mapped = nothing
            report = nothing
            last_report = nothing
            lipschitz = 0.0
            rejected_lengths = Float64[]
            while step_length >= min(remaining, available) * 1e-7
                displacement = (step_length * direction[1], step_length * direction[2])
                point_map, lipschitz = _translation_step_map(
                    center, displacement, inner_radius, support_radius)
                candidate_routes = _piecewise_affine_bump_routes(
                    routes, point_map, center, inner_radius, support_radius)
                candidate = _map_layout_geometry(
                    punctures, ports, routes, point_map;
                    mapped_routes=candidate_routes)
                candidate_report = _layout_candidate_report(
                    H, diagram_hash, candidate...)
                last_report = candidate_report
                if candidate_report.valid
                    mapped, report = candidate, candidate_report
                    break
                end
                push!(rejected_lengths, step_length)
                step_length /= 2
            end
            mapped === nothing && throw(ArgumentError(
                "$phase puncture $puncture_id waypoint $waypoint_index cannot be PL-straightened " *
                "after an injective bump step without violating strict layout validation: " *
                _layout_issue_summary(last_report)))
            punctures, ports, routes = mapped
            selected_index = findfirst(puncture -> puncture.id == puncture_id, punctures)
            push!(stage_records, Dict{String,Any}(
                "waypoint" => waypoint_index,
                "step_length" => step_length,
                "inner_radius" => inner_radius,
                "support_radius" => support_radius,
                "bump_lipschitz" => lipschitz,
                "rejected_step_lengths" => rejected_lengths,
                "validation_checks" => Dict(string(key) => value for (key, value) in report.checks),
            ))
        end
    end
    punctures, ports, routes, Dict{String,Any}(
        "phase" => phase,
        "puncture_id" => puncture_id,
        "path" => [collect(point) for point in path],
        "homeomorphism_steps" => stage_records,
        "step_count" => length(stage_records),
    )
end

function _transport_to_placement(H, diagram_hash, punctures, ports, routes, target_placement)
    targets = Dict(puncture.id => puncture for puncture in target_placement)
    Set(keys(targets)) == Set(puncture.id for puncture in punctures) ||
        throw(ArgumentError("target placement must contain every realized puncture exactly once"))
    initial_report = _require_valid_geometry(
        H, diagram_hash, punctures, ports, routes, "pre-transport realization")
    records = Dict{String,Any}[]

    pair_distances = Float64[]
    for first_index in 1:length(punctures)-1, second_index in first_index+1:length(punctures)
        push!(pair_distances, _point_distance(punctures[first_index].center,
                                             punctures[second_index].center))
    end
    for first_index in 1:length(target_placement)-1, second_index in first_index+1:length(target_placement)
        push!(pair_distances, _point_distance(target_placement[first_index].center,
                                             target_placement[second_index].center))
    end
    shrink_radius = min(0.04, minimum(pair_distances; init=1.0) / 12)
    for puncture_id in sort(collect(keys(targets)))
        index = findfirst(puncture -> puncture.id == puncture_id, punctures)
        desired = min(punctures[index].radius, shrink_radius)
        punctures, ports, routes, record = _radial_transport_step(
            H, diagram_hash, punctures, ports, routes, puncture_id, desired,
            "shrink-before-transport",
        )
        push!(records, record)
    end

    all_centers = [puncture.center for puncture in punctures]
    append!(all_centers, [puncture.center for puncture in target_placement])
    xmin, xmax = extrema(first.(all_centers))
    ymin, ymax = extrema(last.(all_centers))
    span = max(xmax - xmin, ymax - ymin, 1.0)
    staging_x = xmax + span + 3.0
    ordered_ids = sort(collect(keys(targets)))
    staging = Dict(puncture_id =>
        (staging_x + 0.8(index - 1), ymax + span + 2.0 + 0.9(index - 1))
        for (index, puncture_id) in enumerate(ordered_ids))

    remaining = Set(ordered_ids)
    staged = Set{Int}()
    while !isempty(remaining)
        movable = nothing
        for puncture_id in sort(collect(remaining))
            target = targets[puncture_id]
            destination_clear = all(begin
                other.id == puncture_id ||
                    _point_distance(target.center, other.center) >
                    target.radius + other.radius + 0.08
            end for other in punctures)
            destination_clear || continue
            movable = puncture_id
            break
        end
        if movable === nothing
            candidates = sort(collect(setdiff(remaining, staged)))
            isempty(candidates) && throw(ArgumentError(
                "target puncture placement has an unresolved occupancy cycle after staging"))
            puncture_id = first(candidates)
            punctures, ports, routes, record = _move_puncture_along_path(
                H, diagram_hash, punctures, ports, routes,
                puncture_id, staging[puncture_id]; phase="occupancy-cycle-to-staging",
            )
            push!(staged, puncture_id)
            push!(records, record)
            continue
        end
        puncture_id = movable
        phase = puncture_id in staged ? "staging-to-target" : "direct-to-target"
        punctures, ports, routes, record = _move_puncture_along_path(
            H, diagram_hash, punctures, ports, routes,
            puncture_id, targets[puncture_id].center; phase=phase,
        )
        delete!(remaining, puncture_id)
        push!(records, record)
    end
    for puncture_id in ordered_ids
        punctures, ports, routes, record = _radial_transport_step(
            H, diagram_hash, punctures, ports, routes,
            puncture_id, targets[puncture_id].radius, "expand-at-target",
        )
        push!(records, record)
    end

    final_puncture_map = Dict(puncture.id => puncture for puncture in punctures)
    for puncture_id in ordered_ids
        actual, target = final_puncture_map[puncture_id], targets[puncture_id]
        _point_distance(actual.center, target.center) <= 1e-8 || throw(ArgumentError(
            "puncture $puncture_id did not reach its prescribed target center"))
        isapprox(actual.radius, target.radius; atol=1e-10) || throw(ArgumentError(
            "puncture $puncture_id did not reach its prescribed target radius"))
    end
    final_report = _require_valid_geometry(
        H, diagram_hash, punctures, ports, routes, "post-transport grid realization")
    trace = Dict{String,Any}(
        "algorithm" => "compactly-supported-bilipschitz-ambient-isotopy",
        "chart" => "sphere-at-infinity",
        "orientation_preserving" => true,
        "puncture_count" => length(punctures),
        "shrink_radius" => shrink_radius,
        "staging_positions" => Dict(string(key) => collect(value) for (key, value) in staging),
        "stages" => records,
        "homeomorphism_step_count" => sum(get(record, "step_count", 0) for record in records),
        "initial_checks" => Dict(string(key) => value for (key, value) in initial_report.checks),
        "final_checks" => Dict(string(key) => value for (key, value) in final_report.checks),
    )
    punctures, ports, routes, trace
end

function _layout_candidate_report(H, hash, punctures, ports, routes)
    candidate=HeegaardLayout(hash,:alpha_cut,Dict{String,Any}(),punctures,ports,routes,
        _route_crossings(routes),Dict{String,Any}(),
        LayoutValidationReport(false,ValidationIssue[],Dict{Symbol,Bool}()))
    validate_layout(H,candidate;expected_hash=hash)
end
_layout_issue_summary(report)=join([string(i.code) for i in report.issues],", ")
