"""Semantic layouts.  Geometry is derived data and never enters the core hash."""

function _point_for_port(ports, id)
    port = findfirst(port -> port.id == id, ports)
    port === nothing && throw(ArgumentError("unknown port id $id"))
    ports[port].point
end

function _gauss_position(curve_index, curve_count, position, length_word, y)
    x = curve_count == 1 ? 0.0 : 3.0 * (curve_index - 1) / (curve_count - 1) - 1.5
    angle = 2π * (position - 1) / max(length_word, 1)
    (x + 0.8cos(angle), y + 0.8sin(angle))
end

function gauss_layout(H::HeegaardDiagram; kwargs...)
    code = H.code.code
    alpha_positions = Dict{Int,NTuple{2,Float64}}()
    beta_positions = Dict{Int,NTuple{2,Float64}}()
    for (index, curve) in enumerate(code.alpha), (position, point) in enumerate(curve.visits)
        alpha_positions[point] = _gauss_position(index, max(length(code.alpha), 1), position,
                                                 length(curve.visits), 1.5)
    end
    for (index, curve) in enumerate(code.beta), (position, point) in enumerate(curve.visits)
        beta_positions[point] = _gauss_position(index, max(length(code.beta), 1), position,
                                                length(curve.visits), -1.5)
    end
    routes = Route[]
    for (index, point) in enumerate(code.intersections)
        push!(routes, Route(:chord, point.id, index, nothing, nothing,
                            [alpha_positions[point.id], beta_positions[point.id]]))
    end
    report = LayoutValidationReport(true, ValidationIssue[], Dict(:semantic=>true, :geometric=>true))
    HeegaardLayout(canonical_hash(H), :gauss,
        Dict{String,Any}("chart"=>"two-row-circles", "width"=>4.0, "height"=>4.0),
        Puncture[], Port[], routes, Tuple{Int,Int}[],
        Dict{String,Any}("backend"=>"pure-julia", "seed"=>nothing), report)
end

function _segment_orientation(a, b, c)
    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
end

function _between(a, b, c)
    min(a[1], b[1]) - 1e-10 <= c[1] <= max(a[1], b[1]) + 1e-10 &&
    min(a[2], b[2]) - 1e-10 <= c[2] <= max(a[2], b[2]) + 1e-10
end

function _segments_cross(a, b, c, d)
    o1, o2 = _segment_orientation(a, b, c), _segment_orientation(a, b, d)
    o3, o4 = _segment_orientation(c, d, a), _segment_orientation(c, d, b)
    if o1 == 0 && _between(a, b, c) || o2 == 0 && _between(a, b, d) ||
       o3 == 0 && _between(c, d, a) || o4 == 0 && _between(c, d, b)
        return true
    end
    (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0)
end

function _same_point(a, b)
    isapprox(a[1], b[1]; atol=1e-10) && isapprox(a[2], b[2]; atol=1e-10)
end

function _route_crossings(routes::Vector{Route})
    result = Tuple{Int,Int}[]
    for i in 1:length(routes), j in i+1:length(routes)
        first_route, second_route = routes[i], routes[j]
        first_route.family == :alpha && continue
        second_route.family == :alpha && continue
        for first_segment in 1:length(first_route.points)-1,
            second_segment in 1:length(second_route.points)-1
            a, b = first_route.points[first_segment], first_route.points[first_segment+1]
            c, d = second_route.points[second_segment], second_route.points[second_segment+1]
            shared_endpoint = _same_point(a, c) || _same_point(a, d) ||
                              _same_point(b, c) || _same_point(b, d)
            shared_endpoint && continue
            _segments_cross(a, b, c, d) && (push!(result, (i, j)); break)
        end
    end
    unique(result)
end

function _alpha_cut_routes(H::HeegaardDiagram, cut, punctures, ports)
    point_to_port = Dict(port.intersection_id => port.id for port in ports)
    routes = Route[]
    for curve in H.code.code.beta
        visits = curve.visits
        isempty(visits) && continue
        for index in eachindex(visits)
            current = visits[index]
            next = visits[mod1(index + 1, length(visits))]
            start_port, end_port = point_to_port[current], point_to_port[next]
            start_point, end_point = _point_for_port(ports, start_port), _point_for_port(ports, end_port)
            if start_port == end_port
                control = (start_point[1] + 0.3, start_point[2] + 0.3)
                points = [start_point, control, start_point]
            else
                points = [start_point, end_point]
            end
            push!(routes, Route(:beta, curve.id, index, start_port, end_port, points))
        end
    end
    routes
end

function validate_layout(H::HeegaardDiagram, layout::HeegaardLayout)
    issues = ValidationIssue[]
    checks = Dict{Symbol,Bool}()
    layout.diagram_hash == canonical_hash(H) || _add_issue!(issues, :layout, :hash_mismatch,
        "layout does not belong to this core diagram", "diagram_hash")
    checks[:hash] = !any(issue -> issue.code == :hash_mismatch, issues)
    length(layout.punctures) == 2length(H.code.code.alpha) || _add_issue!(issues, :layout, :puncture_count,
        "alpha-cut layout must contain two punctures per alpha curve", "punctures")
    port_ids = [port.id for port in layout.ports]
    length(port_ids) == length(unique(port_ids)) || _add_issue!(issues, :layout, :duplicate_port,
        "port identifiers must be unique", "ports")
    point_ids = [point.id for point in H.code.code.intersections]
    sort([port.intersection_id for port in layout.ports]) == sort(point_ids) ||
        _add_issue!(issues, :layout, :port_incidence, "ports must cover intersections exactly once", "ports")
    endpoint_use = Dict{Int,Int}()
    for route in layout.routes
        route.family == :beta || continue
        route.start_port === nothing || (endpoint_use[route.start_port] = get(endpoint_use, route.start_port, 0) + 1)
        route.end_port === nothing || (endpoint_use[route.end_port] = get(endpoint_use, route.end_port, 0) + 1)
        all(isfinite(point[1]) && isfinite(point[2]) for point in route.points) ||
            _add_issue!(issues, :layout, :nonfinite_coordinate, "route contains a non-finite coordinate", "routes")
        length(route.points) >= 2 || _add_issue!(issues, :layout, :degenerate_route, "route has fewer than two points", "routes")
    end
    all(value == 2 for value in values(endpoint_use)) ||
        _add_issue!(issues, :layout, :endpoint_use, "each cut port must be used by two beta arc endpoints", "routes")
    isempty(layout.crossings) || _add_issue!(issues, :layout, :geometric_crossing,
        "accepted alpha-cut layouts cannot contain non-endpoint beta crossings", "crossings")
    checks[:semantic] = isempty(issues)
    checks[:geometric] = isempty(layout.crossings)
    LayoutValidationReport(isempty(issues), issues, checks)
end

function layout_heegaard(H::HeegaardDiagram; backend=:catlab, seed=nothing, placement=nothing, kwargs...)
    backend in (:catlab, :pure_julia, :fallback) || throw(ArgumentError("unknown layout backend $backend"))
    cut = cut_along_alpha(H)
    punctures = placement === nothing ? default_puncture_placement(cut) : placement
    ports = place_ports(cut, punctures)
    routes = _alpha_cut_routes(H, cut, punctures, ports)
    crossings = _route_crossings(routes)
    trace = Dict{String,Any}(
        "requested_backend" => String(backend),
        "backend" => backend == :catlab ? "pure-julia-fallback" : String(backend),
        "seed" => seed,
        "catlab_used" => false,
    )
    provisional = HeegaardLayout(canonical_hash(H), :alpha_cut,
        Dict{String,Any}("chart"=>"sphere-at-infinity", "outer_halo"=>true),
        punctures, ports, routes, crossings, trace,
        LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
    report = validate_layout(H, provisional)
    HeegaardLayout(provisional.diagram_hash, provisional.view, provisional.chart,
                   provisional.punctures, provisional.ports, provisional.routes,
                   provisional.crossings, provisional.backend_trace, report)
end
