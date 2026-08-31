"""Semantic layouts.  Geometry is derived data and never enters the core hash."""

function _curve_points(curve::PolylineCurve)
    copy(curve.points)
end

function _midpoint(first, second)
    ((first[1] + second[1]) / 2, (first[2] + second[2]) / 2)
end

function _cubic_flat_enough(segment::CubicBezierSegment)
    chord = (segment.p0, segment.p3)
    max(_segment_distance_squared(segment.p1, chord[1], chord[2]),
        _segment_distance_squared(segment.p2, chord[1], chord[2])) <= 0.0
end

function _flatten_cubic(segment::CubicBezierSegment, tolerance, max_depth, depth=0)
    chord_distance = max(sqrt(_segment_distance_squared(segment.p1, segment.p0, segment.p3)),
                         sqrt(_segment_distance_squared(segment.p2, segment.p0, segment.p3)))
    (depth >= max_depth || chord_distance <= tolerance) && return [segment.p0, segment.p3]
    q01 = _midpoint(segment.p0, segment.p1)
    q12 = _midpoint(segment.p1, segment.p2)
    q23 = _midpoint(segment.p2, segment.p3)
    r012 = _midpoint(q01, q12)
    r123 = _midpoint(q12, q23)
    center = _midpoint(r012, r123)
    left = CubicBezierSegment((segment.p0, q01, r012, center))
    right = CubicBezierSegment((center, r123, q23, segment.p3))
    vcat(_flatten_cubic(left, tolerance, max_depth, depth + 1)[1:end-1],
         _flatten_cubic(right, tolerance, max_depth, depth + 1))
end

function _curve_points(curve::BezierChain)
    result = NTuple{2,Float64}[]
    for (index, segment) in enumerate(curve.segments)
        flattened = _flatten_cubic(segment, curve.flatness, curve.max_depth)
        index == 1 ? append!(result, flattened) : append!(result, flattened[2:end])
    end
    result
end

_flatten_curve(curve::RouteCurve) = _curve_points(curve)

function _curve_points(curve::RouteCurve)
    throw(ArgumentError("unsupported route curve representation $(typeof(curve))"))
end

function _route_points(route::Route)
    _curve_points(route.curve)
end

function _line_curve(first, second)
    CubicBezierSegment((first, first, second, second))
end

function _route_corner_step(points, route_index, vertex_index, routes, punctures;
                            step_fraction, clearance_fraction, puncture_fraction)
    vertex = points[vertex_index]
    previous = points[vertex_index - 1]
    next = points[vertex_index + 1]
    previous_length = hypot(vertex[1] - previous[1], vertex[2] - previous[2])
    next_length = hypot(next[1] - vertex[1], next[2] - vertex[2])
    local_bound = step_fraction * min(previous_length, next_length)

    other_clearance = Inf
    for (other_index, other_route) in enumerate(routes)
        other_index == route_index && continue
        other_points = _route_points(other_route)
        for segment_index in 1:length(other_points)-1
            other_clearance = min(other_clearance,
                sqrt(_segment_distance_squared(vertex, other_points[segment_index],
                                               other_points[segment_index + 1])))
        end
    end

    puncture_margin = Inf
    for puncture in punctures
        distance = hypot(vertex[1] - puncture.center[1], vertex[2] - puncture.center[2])
        puncture_margin = min(puncture_margin, distance - puncture.radius)
    end
    min(local_bound, clearance_fraction * other_clearance,
        puncture_fraction * max(puncture_margin, 0.0))
end

function _mechanically_smooth_route(route::Route, route_index, routes, punctures;
                                    step_fraction, clearance_fraction,
                                    puncture_fraction, flatness, max_depth)
    route.curve isa PolylineCurve || throw(ArgumentError(
        "mechanical smoothing accepts only a feasible polyline route"))
    points = _route_points(route)
    segments = CubicBezierSegment[]
    current = first(points)
    steps = Float64[]
    smoothed_corners = 0
    for vertex_index in 2:length(points)-1
        vertex = points[vertex_index]
        step = _route_corner_step(points, route_index, vertex_index, routes, punctures;
                                   step_fraction, clearance_fraction, puncture_fraction)
        push!(steps, step)
        previous = points[vertex_index - 1]
        next = points[vertex_index + 1]
        turn = (vertex[1] - previous[1]) * (next[2] - vertex[2]) -
               (vertex[2] - previous[2]) * (next[1] - vertex[1])
        scale = max(hypot(vertex[1] - previous[1], vertex[2] - previous[2]) *
                    hypot(next[1] - vertex[1], next[2] - vertex[2]), 1.0)
        if step <= eps(Float64) || abs(turn) <= 1e-12 * scale
            push!(segments, _line_curve(current, vertex))
            current = vertex
            continue
        end
        incoming_length = hypot(vertex[1] - previous[1], vertex[2] - previous[2])
        outgoing_length = hypot(next[1] - vertex[1], next[2] - vertex[2])
        incoming = (vertex[1] + step * (previous[1] - vertex[1]) / incoming_length,
                    vertex[2] + step * (previous[2] - vertex[2]) / incoming_length)
        outgoing = (vertex[1] + step * (next[1] - vertex[1]) / outgoing_length,
                    vertex[2] + step * (next[2] - vertex[2]) / outgoing_length)
        push!(segments, _line_curve(current, incoming))
        push!(segments, CubicBezierSegment((incoming, vertex, vertex, outgoing)))
        current = outgoing
        smoothed_corners += 1
    end
    push!(segments, _line_curve(current, last(points)))
    curve = BezierChain(segments; flatness, max_depth)
    flattened = _curve_points(curve)
    smoothed = Route(route.family, route.curve_id, route.segment_id,
                     route.start_port, route.end_port, flattened; curve)
    smoothed, (steps=steps, smoothed_corners=smoothed_corners,
               max_step=isempty(steps) ? 0.0 : maximum(steps),
               min_positive_step=isempty(filter(step -> step > 0.0, steps)) ? 0.0 :
                   minimum(filter(step -> step > 0.0, steps)))
end

"""Mechanically round a valid feasible layout inside a clearance tube.

The input must already be a valid alpha-cut PL layout.  Each corner step is
bounded by incident-edge length, distance to other routes, and puncture
margin.  The returned layout is accepted only after flattening its actual
Bezier representation and re-running the full layout validator.
"""
function smooth_layout(H::HeegaardDiagram, layout::HeegaardLayout;
                       step_fraction=0.18, clearance_fraction=0.25,
                       puncture_fraction=0.5, flatness=1e-4, max_depth=20)
    layout.view == :alpha_cut || throw(ArgumentError("mechanical smoothing requires an alpha-cut layout"))
    0.0 < step_fraction <= 1 / 3 || throw(ArgumentError("step_fraction must be in (0, 1/3]"))
    0.0 < clearance_fraction <= 0.5 || throw(ArgumentError("clearance_fraction must be in (0, 1/2]"))
    0.0 < puncture_fraction <= 0.5 || throw(ArgumentError("puncture_fraction must be in (0, 1/2]"))
    input_report = validate_layout(H, layout)
    input_report.valid || throw(ArgumentError("mechanical smoothing requires a valid feasible layout"))
    all(route.curve isa PolylineCurve for route in layout.routes) ||
        throw(ArgumentError("mechanical smoothing requires PL input routes"))

    routes = Route[]
    route_metrics = Any[]
    for (route_index, route) in enumerate(layout.routes)
        smoothed, metrics = _mechanically_smooth_route(
            route, route_index, layout.routes, layout.punctures;
            step_fraction, clearance_fraction, puncture_fraction, flatness, max_depth)
        push!(routes, smoothed)
        push!(route_metrics, Dict{String,Any}(
            "curve_id" => route.curve_id, "segment_id" => route.segment_id,
            "smoothed_corners" => metrics.smoothed_corners,
            "max_step" => metrics.max_step,
            "min_positive_step" => metrics.min_positive_step,
            "flattened_points" => length(smoothed.points)))
    end
    trace = deepcopy(layout.backend_trace)
    trace["smoothing"] = Dict{String,Any}(
        "method" => "mechanical-clearance-tube",
        "step_fraction" => Float64(step_fraction),
        "clearance_fraction" => Float64(clearance_fraction),
        "puncture_fraction" => Float64(puncture_fraction),
        "flatness" => Float64(flatness), "max_depth" => Int(max_depth),
        "route_metrics" => route_metrics)
    haskey(trace, "pipeline") &&
        (trace["pipeline"]["smoothing"] = "mechanical-clearance-tube")
    candidate = HeegaardLayout(layout.diagram_hash, layout.view, deepcopy(layout.chart),
                               deepcopy(layout.punctures), deepcopy(layout.ports), routes,
                               Tuple{Int,Int}[], trace,
                               LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
    report = validate_layout(H, candidate; expected_hash=layout.diagram_hash)
    report.valid || throw(ArgumentError(
        "mechanical smoothing failed closed: " *
        join(["$(issue.code): $(issue.message)" for issue in report.issues], "; ")))
    HeegaardLayout(candidate.diagram_hash, candidate.view, candidate.chart,
                   candidate.punctures, candidate.ports, candidate.routes,
                   _route_crossings(candidate.routes), candidate.backend_trace, report)
end

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
    o1 * o2 < 0 && o3 * o4 < 0
end

function _same_point(a, b)
    isapprox(a[1], b[1]; atol=1e-10) && isapprox(a[2], b[2]; atol=1e-10)
end

function _route_crossings(routes::Vector{Route})
    segments = _route_segments(routes)
    route_pairs = Set{Tuple{Int,Int}}()
    for (first_id, second_id) in _spatial_segment_pairs(segments)
        first, second = segments[first_id], segments[second_id]
        first.route == second.route && continue
        routes[first.route].family == :alpha && continue
        routes[second.route].family == :alpha && continue
        pair = minmax(first.route, second.route)
        pair in route_pairs && continue
        shared_endpoint = _same_point(first.a, second.a) || _same_point(first.a, second.b) ||
                          _same_point(first.b, second.a) || _same_point(first.b, second.b)
        shared_endpoint && continue
        _segments_cross(first.a, first.b, second.a, second.b) && push!(route_pairs, pair)
    end
    sort!(collect(route_pairs))
end

function _self_crosses(route::Route)
    segments = _route_segments(Route[route])
    for (first_id, second_id) in _spatial_segment_pairs(segments)
        first, second = segments[first_id], segments[second_id]
        abs(first.segment - second.segment) <= 1 && continue
        _segments_cross(first.a, first.b, second.a, second.b) && return true
    end
    false
end

function _route_segments(routes)
    segments = NamedTuple[]
    for (route_id, route) in enumerate(routes)
        points = _route_points(route)
        for segment_id in 1:length(points)-1
            a, b = points[segment_id], points[segment_id + 1]
            push!(segments, (route=route_id, segment=segment_id, a=a, b=b,
                             xmin=min(a[1], b[1]), xmax=max(a[1], b[1]),
                             ymin=min(a[2], b[2]), ymax=max(a[2], b[2])))
        end
    end
    segments
end

"""Bounding-box broad phase; every returned pair still receives the exact segment test."""
function _spatial_segment_pairs(segments)
    isempty(segments) && return Tuple{Int,Int}[]
    xmin = minimum(segment.xmin for segment in segments)
    xmax = maximum(segment.xmax for segment in segments)
    ymin = minimum(segment.ymin for segment in segments)
    ymax = maximum(segment.ymax for segment in segments)
    span = max(xmax - xmin, ymax - ymin, 1.0)
    cell_size = max(span / max(ceil(Int, sqrt(length(segments))), 1), 1e-4)
    cells = Dict{Tuple{Int,Int},Vector{Int}}()
    for (segment_id, segment) in enumerate(segments)
        first_x = floor(Int, (segment.xmin - xmin) / cell_size)
        last_x = floor(Int, (segment.xmax - xmin) / cell_size)
        first_y = floor(Int, (segment.ymin - ymin) / cell_size)
        last_y = floor(Int, (segment.ymax - ymin) / cell_size)
        for x in first_x:last_x, y in first_y:last_y
            push!(get!(cells, (x, y), Int[]), segment_id)
        end
    end
    pairs = Set{Tuple{Int,Int}}()
    for ids in values(cells)
        length(ids) < 2 && continue
        for first_index in 1:length(ids)-1, second_index in first_index+1:length(ids)
            first_id, second_id = minmax(ids[first_index], ids[second_index])
            first, second = segments[first_id], segments[second_id]
            overlaps = !(first.xmax < second.xmin || second.xmax < first.xmin ||
                         first.ymax < second.ymin || second.ymax < first.ymin)
            overlaps && push!(pairs, (first_id, second_id))
        end
    end
    collect(pairs)
end

function _segment_distance_squared(point, first, second)
    dx, dy = second[1] - first[1], second[2] - first[2]
    denominator = dx * dx + dy * dy
    denominator == 0 && return (point[1] - first[1])^2 + (point[2] - first[2])^2
    parameter = clamp(((point[1] - first[1]) * dx + (point[2] - first[2]) * dy) /
                      denominator, 0.0, 1.0)
    closest = (first[1] + parameter * dx, first[2] + parameter * dy)
    (point[1] - closest[1])^2 + (point[2] - closest[2])^2
end

function _route_enters_puncture(route::Route, puncture::Puncture; tolerance=1e-9)
    threshold = max(puncture.radius - tolerance, 0.0)^2
    points = _route_points(route)
    any(_segment_distance_squared(puncture.center, points[index],
                                  points[index + 1]) < threshold
        for index in 1:length(points)-1)
end

function _alpha_cut_routes(H::HeegaardDiagram, cut, punctures, ports)
    point_to_port = _halfedge_port_map(cut, ports)
    routes = Route[]
    for curve in H.code.code.beta
        visits = curve.visits
        isempty(visits) && continue
        for index in eachindex(visits)
            current = visits[index]
            next = visits[mod1(index + 1, length(visits))]
            start_port = point_to_port[(current, :out)]
            end_port = point_to_port[(next, :in)]
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

function _cyclic_equal(first_cycle, second_cycle)
    length(first_cycle) == length(second_cycle) || return false
    isempty(first_cycle) && return true
    any(all(first_cycle[index] == second_cycle[mod1(index + shift, length(second_cycle))]
            for index in eachindex(first_cycle))
        for shift in 0:length(second_cycle)-1)
end

function _terminal_ribbon_report(H::HeegaardDiagram, cut::CutSurface,
                                 ports::Vector{Port}, wiring::WiringProblem)
    port_map = Dict(port.id => port for port in ports)
    rotations = Dict(boundary.id =>
        [port.id for port in sort([port for port in ports
                                   if port.boundary_id == boundary.id];
                                  by=port -> (-port.parameter, -port.id))]
        for boundary in cut.boundaries)
    origin_port = Dict{Int,Int}()
    origin_dart = Dict{Int,Int}()
    for (edge, wire) in enumerate(wiring.wires)
        forward, reverse = 2edge - 1, 2edge
        origin_port[forward], origin_port[reverse] = wire.from_port, wire.to_port
        origin_dart[wire.from_port], origin_dart[wire.to_port] = forward, reverse
    end
    next_dart = Dict{Int,Int}()
    for dart in sort(collect(keys(origin_port)))
        reverse = isodd(dart) ? dart + 1 : dart - 1
        target_port = origin_port[reverse]
        boundary = port_map[target_port].boundary_id
        rotation = rotations[boundary]
        target_index = findfirst(==(target_port), rotation)
        target_index === nothing && throw(ArgumentError("terminal rotation omits port $target_port"))
        predecessor = rotation[mod1(target_index - 1, length(rotation))]
        next_dart[dart] = origin_dart[predecessor]
    end
    seen = Set{Int}()
    faces = Vector{Vector{Int}}()
    for dart in sort(collect(keys(next_dart)))
        dart in seen && continue
        face = Int[]
        current = dart
        while !(current in seen)
            push!(seen, current)
            push!(face, current)
            current = next_dart[current]
        end
        current == first(face) || throw(ArgumentError("terminal ribbon permutation is not cyclic"))
        push!(faces, face)
    end
    vertex_count = length(rotations)
    edge_count = length(wiring.wires)
    face_count = length(faces)
    euler = vertex_count - edge_count + face_count
    genus = (2 - euler) ÷ 2
    graph = ribbon_graph(H)
    original_cycles = boundary_cycles(graph)
    alpha_edge_count = sum(length(curve.visits) for curve in H.code.code.alpha)
    original_beta_darts = Set(vcat(
        [2(alpha_edge_count + edge) - 1 for edge in 1:edge_count],
        [2(alpha_edge_count + edge) for edge in 1:edge_count],
    ))
    filtered_original = [[dart for dart in cycle if dart in original_beta_darts]
                         for cycle in original_cycles]
    any(isempty, filtered_original) &&
        throw(ArgumentError("a ribbon boundary cycle contains no beta dart"))
    face_to_boundary = Int[]
    for face in faces
        mapped = [isodd(dart) ? 2(alpha_edge_count + cld(dart, 2)) - 1 :
                                2(alpha_edge_count + cld(dart, 2))
                  for dart in face]
        match = findfirst(cycle -> _cyclic_equal(mapped, cycle), filtered_original)
        match === nothing &&
            throw(ArgumentError("terminal face does not match an original ribbon boundary cycle"))
        push!(face_to_boundary, match)
    end
    sort(face_to_boundary) == collect(eachindex(filtered_original)) ||
        throw(ArgumentError("terminal faces do not bijectively match ribbon boundary cycles"))
    (planar=euler == 2, vertices=vertex_count, edges=edge_count,
     faces=face_count, euler=euler, genus=genus,
     face_lengths=sort(length.(faces)), face_cycles=faces,
     face_to_boundary_cycle=face_to_boundary,
     boundary_cycle_bijection=true)
end

"""Run Catlab's wiring-diagram layout on an ordered, non-topological view.

The Catlab object deliberately contains one unary box per ordered beta arc and
one outer input/output port for that arc.  It is an adapter for initial route
anchors only: the returned positions are never used to choose punctures,
validate intersections, or compute a surface invariant.
"""
function _catlab_wiring_layout(problem::WiringProblem)
    catlab = Base.require(Main, :Catlab)
    theories = getproperty(catlab, :Theories)
    wiring = getproperty(catlab, :WiringDiagrams)
    graphics = getproperty(catlab, :Graphics)
    category = getproperty(theories, :FreeSymmetricMonoidalCategory)
    Ob = getproperty(theories, :Ob)
    Hom = getproperty(theories, :Hom)
    WiringDiagram = getproperty(wiring, :WiringDiagram)
    add_box! = getproperty(wiring, Symbol("add_box!"))
    add_wire! = getproperty(wiring, Symbol("add_wire!"))
    input_id = getproperty(wiring, :input_id)
    output_id = getproperty(wiring, :output_id)
    input_ports = getproperty(wiring, :input_ports)
    output_ports = getproperty(wiring, :output_ports)
    layout_diagram = getproperty(graphics, :layout_diagram)
    layouts = getproperty(graphics, :WiringDiagramLayouts)
    position = getproperty(layouts, :position)

    count = length(problem.wires)
    port_type = Base.invokelatest(Ob, category, :heegaard_port)
    diagram = Base.invokelatest(WiringDiagram, fill(port_type, count), fill(port_type, count))
    for (index, wire) in enumerate(problem.wires)
        box = Base.invokelatest(add_box!, diagram, Base.invokelatest(wiring.Box,
            Base.invokelatest(Hom, Symbol("beta$(wire.beta_id)_$(index)"), port_type, port_type),
            [port_type], [port_type]))
        Base.invokelatest(add_wire!, diagram,
                          Pair((Base.invokelatest(input_id, diagram), index), (box, 1)))
        Base.invokelatest(add_wire!, diagram,
                          Pair((box, 1), (Base.invokelatest(output_id, diagram), index)))
    end
    laid_out = Base.invokelatest(layout_diagram, category, diagram; outer_ports_layout=:fixed)
    in_layout = Base.invokelatest(input_ports, laid_out)
    out_layout = Base.invokelatest(output_ports, laid_out)
    # Do not broadcast over Catlab's SVector here.  On Julia 1.12 the
    # StaticArrays `size` method can be newer than the world in which this
    # adapter was compiled; indexed access through invokelatest keeps the
    # optional backend boundary world-age safe.
    point = p -> (Float64(Base.invokelatest(getindex, p, 1)),
                  Float64(Base.invokelatest(getindex, p, 2)))
    inputs = [point(Base.invokelatest(position, port)) for port in in_layout]
    outputs = [point(Base.invokelatest(position, port)) for port in out_layout]
    box_layout = Base.invokelatest(wiring.boxes, laid_out)
    boxes = [point(Base.invokelatest(position, box.value)) for box in box_layout]
    (inputs=inputs, outputs=outputs, boxes=boxes, box_count=count,
     wire_count=length(Base.invokelatest(wiring.wires, laid_out)))
end

function _catlab_route_points(route::Route, anchor)
    isempty(route.points) && return route.points
    start_point, end_point = first(route.points), last(route.points)
    cat_x, cat_y = anchor
    offset = 0.12 * tanh(cat_y)
    if _same_point(start_point, end_point)
        control = (start_point[1] + 0.3 + 0.04 * tanh(cat_x),
                   start_point[2] + 0.3 + offset)
        return [start_point, control, end_point]
    end
    midpoint = ((start_point[1] + end_point[1]) / 2,
                (start_point[2] + end_point[2]) / 2 + offset)
    [start_point, midpoint, end_point]
end

function _catlab_route_anchors(problem::WiringProblem)
    # Catlab is an optional late-loaded adapter.  The adapter function may have
    # been compiled before Catlab installs its generated `Ob` methods, so enter
    # the optional backend through the current world explicitly.
    result = Base.invokelatest(_catlab_wiring_layout, problem)
    anchors = [((result.inputs[index][1] + result.outputs[index][1]) / 2,
                (result.inputs[index][2] + result.outputs[index][2]) / 2)
               for index in eachindex(problem.wires)]
    result, anchors
end

function _route_certificate_points(path::AbstractString, diagram_hash::AbstractString)
    isfile(path) || throw(ArgumentError("route certificate does not exist: $path"))
    parsed = JSON3.read(read(path, String))
    encoded = if parsed isa AbstractVector
        parsed
    else
        schema = string(_spec_field(parsed, :schema, ""))
        schema == "akaghef.heegaard.layout/1" ||
            throw(ArgumentError("route certificate must be a route array or .hgd.layout.json"))
        stored_hash = string(_spec_field(parsed, :diagram_hash, ""))
        stored_hash == diagram_hash ||
            throw(ArgumentError("route certificate core hash does not match the diagram"))
        [_spec_field(route, :points, Any[])
         for route in _as_vector(_spec_field(parsed, :routes, Any[]))]
    end
    [[begin
        coordinates = collect(point)
        length(coordinates) == 2 || throw(ArgumentError("route certificate points must be 2D"))
        value = (Float64(coordinates[1]), Float64(coordinates[2]))
        all(isfinite, value) || throw(ArgumentError("route certificate contains a non-finite point"))
        value
    end for point in _as_vector(route)] for route in _as_vector(encoded)]
end

function _apply_route_certificate(routes, ports, path::AbstractString, diagram_hash::AbstractString)
    points = _route_certificate_points(path, diagram_hash)
    length(points) == length(routes) || throw(ArgumentError(
        "route certificate has $(length(points)) routes; expected $(length(routes))"))
    port_map = Dict(port.id => port for port in ports)
    [begin
        length(points[index]) >= 2 || throw(ArgumentError("route certificate route $index is degenerate"))
        route.start_port === nothing && throw(ArgumentError("certified beta route lacks a start port"))
        route.end_port === nothing && throw(ArgumentError("certified beta route lacks an end port"))
        certified = copy(points[index])
        certified[1] = port_map[route.start_port].point
        certified[end] = port_map[route.end_port].point
        Route(route.family, route.curve_id, route.segment_id,
              route.start_port, route.end_port, certified)
    end for (index, route) in enumerate(routes)]
end

function validate_layout(H::HeegaardDiagram, layout::HeegaardLayout; expected_hash=nothing)
    issues = ValidationIssue[]
    checks = Dict{Symbol,Bool}()
    core_hash = expected_hash === nothing ? canonical_hash(H) : string(expected_hash)
    layout.diagram_hash == core_hash || _add_issue!(issues, :layout, :hash_mismatch,
        "layout does not belong to this core diagram", "diagram_hash")
    checks[:hash] = !any(issue -> issue.code == :hash_mismatch, issues)
    length(layout.punctures) == 2length(H.code.code.alpha) || _add_issue!(issues, :layout, :puncture_count,
        "alpha-cut layout must contain two punctures per alpha curve", "punctures")
    expected_boundary_keys = Set((pair, side)
        for pair in 1:length(H.code.code.alpha), side in (:plus, :minus))
    actual_boundary_keys = Set((puncture.pair, puncture.side) for puncture in layout.punctures)
    actual_boundary_keys == expected_boundary_keys || _add_issue!(issues, :layout, :puncture_pairs,
        "layout must contain exactly one plus/minus puncture for every alpha curve", "punctures")
    length(Set(puncture.id for puncture in layout.punctures)) == length(layout.punctures) ||
        _add_issue!(issues, :layout, :duplicate_puncture,
                    "puncture identifiers must be unique", "punctures")
    port_ids = [port.id for port in layout.ports]
    length(port_ids) == length(unique(port_ids)) || _add_issue!(issues, :layout, :duplicate_port,
        "port identifiers must be unique", "ports")
    code = H.code.code
    points = Dict(point.id => point for point in code.intersections)
    point_ids = sort(collect(keys(points)))
    sort([port.intersection_id for port in layout.ports]) == sort(vcat(point_ids, point_ids)) ||
        _add_issue!(issues, :layout, :port_incidence, "ports must contain two cut half-edges per intersection", "ports")
    punctures = Dict(puncture.id => puncture for puncture in layout.punctures)
    alpha_pairs = Dict(curve.id => pair for (pair, curve) in enumerate(code.alpha))
    port_roles = Dict{Int,Symbol}()
    for port in layout.ports
        point = get(points, port.intersection_id, nothing)
        puncture = get(punctures, port.boundary_id, nothing)
        if point === nothing || puncture === nothing || get(alpha_pairs, point.alpha, 0) != puncture.pair ||
           port.beta_id != point.beta
            _add_issue!(issues, :layout, :port_boundary, "port boundary/beta incidence does not match its intersection", "ports")
            continue
        end
        port_roles[port.id] = (puncture.side == :plus) == (point.sign == 1) ? :out : :in
        abs(hypot(port.point[1] - puncture.center[1],
                  port.point[2] - puncture.center[2]) - puncture.radius) <= 1e-8 ||
            _add_issue!(issues, :layout, :port_geometry,
                        "port $(port.id) is not on its puncture boundary", "ports")
    end
    for point_id in point_ids
        ids = [port.id for port in layout.ports if port.intersection_id == point_id]
        roles = sort([get(port_roles, id, :invalid) for id in ids]; by=string)
        roles == [:in, :out] || _add_issue!(issues, :layout, :halfedge_roles,
            "intersection $point_id must have one beta_in and one beta_out cut port", "ports")
    end
    if length(punctures) == 2length(code.alpha)
        cut_for_ports = cut_along_alpha(H; diagram_hash=core_hash)
        try
            expected_ports = Dict(port.id => port
                for port in place_ports(cut_for_ports, layout.punctures))
            for port in layout.ports
                expected = get(expected_ports, port.id, nothing)
                expected !== nothing &&
                    (port.intersection_id, port.boundary_id, port.beta_id) ==
                    (expected.intersection_id, expected.boundary_id, expected.beta_id) &&
                    isapprox(port.parameter, expected.parameter; atol=1e-12) ||
                    _add_issue!(issues, :layout, :port_order,
                                "port $(port.id) does not preserve its boundary cyclic order", "ports")
            end
        catch
            _add_issue!(issues, :layout, :port_order,
                        "port cyclic order could not be reconstructed from puncture placement", "ports")
        end
    end
    endpoint_use = Dict{Int,Int}()
    actual_routes = Set{Tuple{Int,Int,Int,Int}}()
    port_map = Dict(port.id => port for port in layout.ports)
    curve_geometry_valid = true
    for route in layout.routes
        if route.family != :beta
            _add_issue!(issues, :layout, :route_family,
                        "alpha-cut route family must be beta", "routes")
            continue
        end
        route.start_port === nothing || get(port_roles, route.start_port, :invalid) == :out ||
            _add_issue!(issues, :layout, :route_start_role, "beta route must start at a beta_out port", "routes")
        route.end_port === nothing || get(port_roles, route.end_port, :invalid) == :in ||
            _add_issue!(issues, :layout, :route_end_role, "beta route must end at a beta_in port", "routes")
        route.start_port === nothing || (endpoint_use[route.start_port] = get(endpoint_use, route.start_port, 0) + 1)
        route.end_port === nothing || (endpoint_use[route.end_port] = get(endpoint_use, route.end_port, 0) + 1)
        if route.start_port !== nothing && route.end_port !== nothing
            push!(actual_routes, (route.curve_id, route.segment_id,
                                  route.start_port, route.end_port))
        end
        route_points = try
            _route_points(route)
        catch error
            curve_geometry_valid = false
            _add_issue!(issues, :layout, :unsupported_curve,
                        "route curve cannot be flattened: $(sprint(showerror, error))", "routes")
            nothing
        end
        route_points === nothing && continue
        route_points == route.points || begin
            curve_geometry_valid = false
            _add_issue!(issues, :layout, :curve_points_mismatch,
                        "route flattened cache does not match its curve representation", "routes")
        end
        if route.start_port !== nothing && route.end_port !== nothing
            start = get(port_map, route.start_port, nothing)
            finish = get(port_map, route.end_port, nothing)
            start === nothing || isempty(route_points) || first(route_points) == start.point ||
                _add_issue!(issues, :layout, :route_start_geometry,
                            "route curve does not start exactly at its declared port", "routes")
            finish === nothing || isempty(route_points) || last(route_points) == finish.point ||
                _add_issue!(issues, :layout, :route_end_geometry,
                            "route curve does not end exactly at its declared port", "routes")
        end
        all(isfinite(point[1]) && isfinite(point[2]) for point in route_points) ||
            _add_issue!(issues, :layout, :nonfinite_coordinate, "route contains a non-finite coordinate", "routes")
        length(route_points) >= 2 || _add_issue!(issues, :layout, :degenerate_route, "route has fewer than two points", "routes")
        any(_same_point(route_points[index], route_points[index + 1])
            for index in 1:length(route_points)-1) &&
            _add_issue!(issues, :layout, :degenerate_segment,
                        "route contains a zero-length segment", "routes")
        _self_crosses(route) && _add_issue!(issues, :layout, :self_crossing,
            "route has a non-adjacent self-crossing", "routes")
        any(_route_enters_puncture(route, puncture) for puncture in layout.punctures) &&
            _add_issue!(issues, :layout, :puncture_intrusion,
                        "route enters a puncture interior", "routes")
    end
    cut = cut_along_alpha(H; diagram_hash=core_hash)
    halfedges = _halfedge_port_map(cut, layout.ports)
    expected_routes = Set{Tuple{Int,Int,Int,Int}}()
    for curve in code.beta, segment_id in eachindex(curve.visits)
        isempty(curve.visits) && continue
        current = curve.visits[segment_id]
        next = curve.visits[mod1(segment_id + 1, length(curve.visits))]
        push!(expected_routes, (curve.id, segment_id,
                                halfedges[(current, :out)], halfedges[(next, :in)]))
    end
    actual_routes == expected_routes || _add_issue!(issues, :layout, :beta_successor,
        "routes do not realize the beta cyclic successor relation", "routes")
    try
        terminal_report = _terminal_ribbon_report(
            H, cut, layout.ports, build_wiring_problem(H, cut, layout.ports))
        terminal_report.planar || _add_issue!(issues, :layout, :terminal_nonplanar,
            "cut terminal rotation system has genus $(terminal_report.genus)", "ports")
    catch error
        _add_issue!(issues, :layout, :terminal_embedding,
                    "terminal ribbon embedding failed: $(sprint(showerror, error))", "ports")
    end
    sort(collect(keys(endpoint_use))) == sort(port_ids) && all(value == 1 for value in values(endpoint_use)) ||
        _add_issue!(issues, :layout, :endpoint_use, "each beta half-edge cut port must be used exactly once", "routes")
    computed_crossings = try
        _route_crossings(layout.routes)
    catch error
        curve_geometry_valid = false
        _add_issue!(issues, :layout, :unsupported_curve,
                    "route crossings cannot be computed: $(sprint(showerror, error))", "routes")
        Tuple{Int,Int}[]
    end
    layout.crossings == computed_crossings || _add_issue!(issues, :layout, :crossing_mismatch,
        "stored crossing list does not match route geometry", "crossings")
    isempty(computed_crossings) || _add_issue!(issues, :layout, :geometric_crossing,
        "accepted alpha-cut layouts cannot contain non-endpoint beta crossings", "crossings")
    checks[:semantic] = isempty(issues)
    checks[:curve_representation] = curve_geometry_valid
    checks[:geometric] = curve_geometry_valid && isempty(computed_crossings) &&
                         !any(issue -> issue.code in (:self_crossing, :puncture_intrusion,
                                                      :degenerate_segment), issues)
    LayoutValidationReport(isempty(issues), issues, checks)
end

function layout_heegaard(H::HeegaardDiagram; backend=:catlab, fallback=:pure_julia,
                         seed=nothing, placement=nothing, route_certificate=nothing,
                         router=:auto, anchors=nothing,
                         grid_step=nothing, margin=1.5, max_attempts=12,
                         max_expansions=25000,
                         kwargs...)
    backend in (:catlab, :pure_julia, :fallback) || throw(ArgumentError("unknown layout backend $backend"))
    fallback in (:pure_julia, :error, :none) || throw(ArgumentError("unknown Catlab fallback policy $fallback"))
    router in (:auto, :none, :feasible) || throw(ArgumentError("unknown route policy $router"))
    router == :feasible && route_certificate !== nothing &&
        throw(ArgumentError("router=:feasible cannot be combined with route_certificate"))
    route_policy = router == :auto ?
        (route_certificate === nothing ? :feasible : :certificate) : router
    diagram_hash = canonical_hash(H)
    cut = cut_along_alpha(H; diagram_hash=diagram_hash)
    punctures = placement === nothing ? default_puncture_placement(cut) : placement
    ports = place_ports(cut, punctures)
    wiring = build_wiring_problem(H, cut, ports)
    router_anchors = anchors
    terminal_report = _terminal_ribbon_report(H, cut, ports, wiring)
    terminal_report.planar || throw(ArgumentError(
        "the fixed-puncture terminal ribbon map is not planar; no feasible alpha-cut route exists"))
    routes = _alpha_cut_routes(H, cut, punctures, ports)
    catlab_error = nothing
    catlab_result = nothing
    if backend == :catlab
        try
            catlab_result, catlab_anchors = _catlab_route_anchors(wiring)
            router_anchors === nothing && (router_anchors = catlab_anchors)
            if route_policy == :none
                routes = [Route(route.family, route.curve_id, route.segment_id,
                                route.start_port, route.end_port,
                                _catlab_route_points(route, catlab_anchors[index]))
                          for (index, route) in enumerate(routes)]
            end
        catch error
            catlab_error = sprint(showerror, error)
            fallback in (:error, :none) && throw(ArgumentError(
                "Catlab layout backend failed; use backend=:pure_julia or fallback=:pure_julia. " *
                catlab_error))
        end
    end
    post_router = "none"
    certificate_trace = nothing
    if route_policy == :feasible
        routes = feasible_fixed_endpoint_routes(H, punctures, ports, wiring;
                                                anchors=router_anchors,
                                                grid_step=grid_step, margin=margin,
                                                max_attempts=max_attempts,
                                                max_expansions=max_expansions)
        post_router = "bounded-fixed-endpoint-feasible"
    elseif route_certificate !== nothing
        certificate_path = string(route_certificate)
        routes = _apply_route_certificate(routes, ports, certificate_path, diagram_hash)
        post_router = "certified-fixed-boundary-pl"
        certificate_trace = Dict{String,Any}(
            "file" => basename(certificate_path),
            "sha256" => bytes2hex(sha256(read(certificate_path))),
        )
    end
    crossings = _route_crossings(routes)
    trace = Dict{String,Any}(
        "requested_backend" => String(backend),
        "backend" => catlab_result === nothing ? (backend == :catlab ? "pure-julia-fallback" : String(backend)) : "catlab",
        "seed" => seed,
        "catlab_used" => catlab_result !== nothing,
        "fallback_policy" => String(fallback),
        "provenance" => source_manifest(H),
        "post_router" => post_router,
        "pipeline" => Dict{String,Any}(
            "puncture_placement" => "fixed-2-by-g",
            "ordered_ports" => true,
            "seed_layout" => catlab_result === nothing ? "pure-julia" : "catlab",
            "fixed_endpoint_fitting" => true,
            "feasible_router" => post_router,
            "smoothing" => "not-run",
        ),
        "feasible_router" => route_policy == :feasible,
        "router_anchor_count" => router_anchors === nothing ? 0 : length(router_anchors),
        "router_limits" => Dict{String,Any}(
            "grid_step" => grid_step,
            "margin" => margin,
            "max_attempts" => max_attempts,
            "max_expansions" => max_expansions,
        ),
        "terminal_embedding" => Dict{String,Any}(
            "planar" => terminal_report.planar,
            "vertices" => terminal_report.vertices,
            "edges" => terminal_report.edges,
            "faces" => terminal_report.faces,
            "euler" => terminal_report.euler,
            "genus" => terminal_report.genus,
            "face_lengths" => terminal_report.face_lengths,
            "boundary_cycle_bijection" => terminal_report.boundary_cycle_bijection,
        ),
    )
    catlab_error === nothing || (trace["catlab_error"] = catlab_error)
    certificate_trace === nothing || (trace["route_certificate"] = certificate_trace)
    catlab_result === nothing || (trace["catlab_layout"] = Dict{String,Any}(
        "box_count" => catlab_result.box_count,
        "wire_count" => catlab_result.wire_count,
        "ordered_port_count" => length(catlab_result.inputs) + length(catlab_result.outputs),
    ))
    provisional = HeegaardLayout(diagram_hash, :alpha_cut,
        Dict{String,Any}("chart"=>"sphere-at-infinity", "outer_halo"=>true),
        punctures, ports, routes, crossings, trace,
        LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
    report = validate_layout(H, provisional; expected_hash=diagram_hash)
    if !report.valid
        codes = join(sort(unique(string(issue.code) for issue in report.issues)), ", ")
        throw(ArgumentError(
            "fixed-endpoint seed is not a feasible Heegaard layout ($codes); " *
            "no verified topological router result is available"))
    end
    if route_policy == :none
        trace["post_router"] = "verified-seed-identity"
        trace["pipeline"]["feasible_router"] = "verified-seed-identity"
    end
    HeegaardLayout(provisional.diagram_hash, provisional.view, provisional.chart,
                   provisional.punctures, provisional.ports, provisional.routes,
                   provisional.crossings, provisional.backend_trace, report)
end
