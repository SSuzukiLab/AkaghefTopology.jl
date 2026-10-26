"""Checked geometric realization of the terminal planar ribbon map.

Tutte's theorem applies to simple 3-connected plane graphs, while terminal
Heegaard maps may have loops, parallel edges, bridges, and cut vertices.  For
each connected terminal component we therefore take a flag-level barycentric
refinement: every barycentric incidence edge is subdivided and every resulting
flag hexagon is coned.  The result is a simple triangulated sphere, hence a
3-connected planar graph.  A convex triangular outer face and uniform positive
barycentric weights then give a crossing-free straight-line drawing.

An original terminal edge is recovered as

    vertex -- vertex/edge incidence -- edge center -- incidence -- vertex.

Small disks truncate the original vertices.  Their intersection with the
incident first segments gives the puncture ports in the prescribed rotation.
"""

_reverse_terminal_dart(dart::Int) = isodd(dart) ? dart + 1 : dart - 1

function _terminal_component_refinement(terminal, component_vertices)
    vertex_set = Set(component_vertices)
    darts = sort([dart for dart in keys(terminal.origin_vertex)
                  if terminal.origin_vertex[dart] in vertex_set])
    isempty(darts) && return (
        trivial=true,
        terminal_vertices=sort(collect(vertex_set)),
        terminal_edges=Int[],
        terminal_faces=Int[],
    )

    dart_set = Set(darts)
    terminal_edges = sort(unique(cld(dart, 2) for dart in darts))
    face_of = Dict{Int,Int}()
    for (face_id, face) in enumerate(terminal.face_cycles), dart in face
        face_of[dart] = face_id
    end
    all(haskey(face_of, dart) for dart in darts) ||
        throw(ArgumentError("terminal face permutation omits a component dart"))
    terminal_faces = sort(unique(face_of[dart] for dart in darts))

    node_of = Dict{Any,Int}()
    labels = Any[]
    function node(key)
        get!(node_of, key) do
            push!(labels, key)
            length(labels)
        end
    end

    vertex_nodes = Dict(vertex => node((:vertex, vertex)) for vertex in sort(collect(vertex_set)))
    edge_nodes = Dict(edge => node((:edge, edge)) for edge in terminal_edges)
    face_nodes = Dict(face => node((:face, face)) for face in terminal_faces)
    ve_nodes = Dict(dart => node((:vertex_edge, dart)) for dart in darts)
    ef_nodes = Dict(dart => node((:edge_face, dart)) for dart in darts)
    vf_nodes = Dict(dart => node((:vertex_face, dart)) for dart in darts)

    triangles = NTuple{3,Int}[]
    for dart in darts
        reverse = _reverse_terminal_dart(dart)
        reverse in dart_set || throw(ArgumentError("terminal component is not closed under edge reversal"))
        edge = cld(dart, 2)
        vertex = terminal.origin_vertex[dart]
        right_corner = terminal.next_dart[reverse]
        right_corner in dart_set ||
            throw(ArgumentError("terminal component is not closed under its face permutation"))

        left_hexagon = (
            vertex_nodes[vertex],
            ve_nodes[dart],
            edge_nodes[edge],
            ef_nodes[dart],
            face_nodes[face_of[dart]],
            vf_nodes[dart],
        )
        right_hexagon = (
            vertex_nodes[vertex],
            vf_nodes[right_corner],
            face_nodes[face_of[reverse]],
            ef_nodes[reverse],
            edge_nodes[edge],
            ve_nodes[dart],
        )
        for (side, hexagon) in ((:left, left_hexagon), (:right, right_hexagon))
            center = node((:flag, side, dart))
            for index in eachindex(hexagon)
                push!(triangles, (center, hexagon[index], hexagon[mod1(index + 1, 6)]))
            end
        end
    end

    edge_multiplicity = Dict{Tuple{Int,Int},Int}()
    for triangle in triangles
        for (first, second) in ((triangle[1], triangle[2]),
                                (triangle[2], triangle[3]),
                                (triangle[3], triangle[1]))
            edge = minmax(first, second)
            first == second && throw(ArgumentError("flag refinement produced a loop"))
            edge_multiplicity[edge] = get(edge_multiplicity, edge, 0) + 1
        end
    end
    all(==(2), values(edge_multiplicity)) ||
        throw(ArgumentError("flag refinement is not a closed triangulated surface"))
    edges = sort(collect(keys(edge_multiplicity)))
    node_count = length(labels)
    length(edges) == 3node_count - 6 ||
        throw(ArgumentError("flag refinement does not satisfy E=3V-6"))
    length(triangles) == 2node_count - 4 ||
        throw(ArgumentError("flag refinement does not satisfy F=2V-4"))

    adjacency = [Int[] for _ in 1:node_count]
    for (first, second) in edges
        push!(adjacency[first], second)
        push!(adjacency[second], first)
    end
    foreach(sort!, adjacency)
    terminal_node_set = Set(values(vertex_nodes))
    outer_node = hasproperty(terminal,:outer_faces) ?
        face_nodes[terminal.outer_faces[minimum(component_vertices)]] : nothing
    outer_index = findfirst(triangle -> all(!(vertex in terminal_node_set) for vertex in triangle) &&
                            (outer_node===nothing || outer_node in triangle),
                            triangles)
    outer_index === nothing && throw(ArgumentError("flag refinement has no auxiliary outer face"))

    (
        trivial=false,
        terminal_vertices=sort(collect(vertex_set)),
        terminal_edges=terminal_edges,
        terminal_faces=terminal_faces,
        labels=labels,
        triangles=triangles,
        edges=edges,
        adjacency=adjacency,
        vertex_nodes=vertex_nodes,
        edge_nodes=edge_nodes,
        ve_nodes=ve_nodes,
        outer_face=triangles[outer_index],
    )
end

function _tutte_outer_triangle(reflected; terminal_count)
    # This is chart data, not a boundary condition for the punctures.  Its
    # size depends only on the terminal component being realized, so a later
    # fixed-puncture transport cannot feed back into the Tutte solve.
    center = (0.0, 0.0)
    scale = max(1.0, sqrt(Float64(terminal_count) / 2))
    triangle = [
        (center[1] - 2scale, center[2] - 1.5scale),
        (center[1] + 2scale, center[2] - 1.5scale),
        (center[1], center[2] + 2scale),
    ]
    reflected ? [triangle[1], triangle[3], triangle[2]] : triangle
end

function _tutte_embedding(refinement, fixed_terminal; reflected=false)
    refinement.trivial && throw(ArgumentError("an isolated terminal vertex has no Tutte refinement"))
    pins = Dict{Int,NTuple{2,Float64}}(
        refinement.vertex_nodes[vertex] => point for (vertex, point) in fixed_terminal)
    outer_points = _tutte_outer_triangle(reflected;
                                         terminal_count=length(refinement.terminal_vertices))
    for (node, point) in zip(refinement.outer_face, outer_points)
        haskey(pins, node) && throw(ArgumentError("terminal pin collides with the auxiliary outer face"))
        pins[node] = point
    end

    node_count = length(refinement.labels)
    free = [node for node in 1:node_count if !haskey(pins, node)]
    free_index = Dict(node => index for (index, node) in enumerate(free))
    laplacian = zeros(Float64, length(free), length(free))
    x_rhs = zeros(Float64, length(free))
    y_rhs = zeros(Float64, length(free))
    for node in free
        row = free_index[node]
        neighbors = refinement.adjacency[node]
        laplacian[row, row] = length(neighbors)
        for neighbor in neighbors
            if haskey(free_index, neighbor)
                laplacian[row, free_index[neighbor]] -= 1
            else
                point = pins[neighbor]
                x_rhs[row] += point[1]
                y_rhs[row] += point[2]
            end
        end
    end
    x_solution = laplacian \ x_rhs
    y_solution = laplacian \ y_rhs
    coordinates = Dict{Int,NTuple{2,Float64}}(pins)
    for node in free
        index = free_index[node]
        coordinates[node] = (x_solution[index], y_solution[index])
    end
    all(all(isfinite, point) for point in values(coordinates)) ||
        throw(ArgumentError("Tutte solve produced a non-finite coordinate"))
    coordinates
end

function _component_rotation_ok(terminal, refinement, coordinates)
    for vertex in refinement.terminal_vertices
        prescribed_clockwise = terminal.rotations[vertex]
        length(prescribed_clockwise) <= 2 && continue
        center = coordinates[refinement.vertex_nodes[vertex]]
        actual_ccw = sort(copy(prescribed_clockwise); by=port -> begin
            dart = terminal.origin_dart[port]
            point = coordinates[refinement.ve_nodes[dart]]
            atan(point[2] - center[2], point[1] - center[1])
        end)
        _cyclic_equal(actual_ccw, reverse(prescribed_clockwise)) || return false
    end
    true
end

function _refinement_crossings(refinement, coordinates)
    crossings = Tuple{Int,Int}[]
    for first_index in 1:length(refinement.edges)-1
        first = refinement.edges[first_index]
        for second_index in first_index+1:length(refinement.edges)
            second = refinement.edges[second_index]
            isempty(intersect(first, second)) || continue
            _segments_cross(coordinates[first[1]], coordinates[first[2]],
                            coordinates[second[1]], coordinates[second[2]]) &&
                push!(crossings, (first_index, second_index))
        end
    end
    crossings
end

function _solve_terminal_coordinates(terminal; fixed_placement=nothing)
    fixed = fixed_placement === nothing ? Dict{Int,NTuple{2,Float64}}() :
        Dict(puncture.id => puncture.center for puncture in fixed_placement)
    refinements = [_terminal_component_refinement(terminal, component)
                   for component in terminal.component_vertices]
    component_coordinates = Vector{Dict{Int,NTuple{2,Float64}}}(undef, length(refinements))

    for (index, refinement) in enumerate(refinements)
        if refinement.trivial
            point = fixed_placement === nothing ? (0.0, 0.0) : fixed[only(refinement.terminal_vertices)]
            component_coordinates[index] = Dict(1 => point)
            continue
        end
        terminal_pins = Dict(vertex => fixed[vertex]
                             for vertex in refinement.terminal_vertices if haskey(fixed, vertex))
        fixed_placement === nothing || length(terminal_pins) == length(refinement.terminal_vertices) ||
            throw(ArgumentError("fixed placement omits a terminal component vertex"))
        selected = nothing
        for reflected in (false, true)
            candidate = _tutte_embedding(refinement, terminal_pins; reflected=reflected)
            _component_rotation_ok(terminal, refinement, candidate) || continue
            selected = candidate
            break
        end
        selected === nothing &&
            throw(ArgumentError("Tutte embedding does not realize the prescribed terminal rotations"))
        component_coordinates[index] = selected
    end

    if fixed_placement === nothing
        nontrivial = [(index, refinement) for (index, refinement) in enumerate(refinements)
                      if !refinement.trivial]
        if !isempty(nontrivial)
            widths = Float64[]
            heights = Float64[]
            bounds = Dict{Int,NTuple{4,Float64}}()
            for (index, _) in nontrivial
                points = collect(values(component_coordinates[index]))
                xmin, xmax = extrema(first.(points))
                ymin, ymax = extrema(last.(points))
                bounds[index] = (xmin, xmax, ymin, ymax)
                push!(widths, xmax - xmin)
                push!(heights, ymax - ymin)
            end
            cell_width = maximum(widths) + 4.0
            cell_height = maximum(heights) + 4.0
            columns = ceil(Int, sqrt(length(refinements)))
            for (index, refinement) in enumerate(refinements)
                if refinement.trivial
                    column = mod(index - 1, columns)
                    row = fld(index - 1, columns)
                    component_coordinates[index][1] = (column * cell_width, -row * cell_height)
                    continue
                end
                xmin, _, ymin, _ = bounds[index]
                column = mod(index - 1, columns)
                row = fld(index - 1, columns)
                shift = (column * cell_width - xmin, -row * cell_height - ymin)
                for node in keys(component_coordinates[index])
                    point = component_coordinates[index][node]
                    component_coordinates[index][node] = (point[1] + shift[1], point[2] + shift[2])
                end
            end
        end
    end

    vertex_points = Dict{Int,NTuple{2,Float64}}()
    edge_points = Dict{Int,NTuple{2,Float64}}()
    incidence_points = Dict{Int,NTuple{2,Float64}}()
    refinement_crossings = 0
    total_nodes = 0
    total_edges = 0
    total_triangles = 0
    for (index, refinement) in enumerate(refinements)
        coordinates = component_coordinates[index]
        if refinement.trivial
            vertex_points[only(refinement.terminal_vertices)] = coordinates[1]
            continue
        end
        for (vertex, node) in refinement.vertex_nodes
            vertex_points[vertex] = coordinates[node]
        end
        for (edge, node) in refinement.edge_nodes
            edge_points[edge] = coordinates[node]
        end
        for (dart, node) in refinement.ve_nodes
            incidence_points[dart] = coordinates[node]
        end
        refinement_crossings += length(_refinement_crossings(refinement, coordinates))
        total_nodes += length(refinement.labels)
        total_edges += length(refinement.edges)
        total_triangles += length(refinement.triangles)
    end
    (
        vertex_points=vertex_points,
        edge_points=edge_points,
        incidence_points=incidence_points,
        refinement_crossings=refinement_crossings,
        refinement_nodes=total_nodes,
        refinement_edges=total_edges,
        refinement_triangles=total_triangles,
    )
end

function _computed_terminal_punctures(cut, terminal, geometry)
    skeleton = NamedTuple{(:first, :second, :incident),
                          Tuple{NTuple{2,Float64},NTuple{2,Float64},Set{Int}}}[]
    for edge in 1:terminal.edges
        forward, reverse = 2edge - 1, 2edge
        first_vertex = terminal.origin_vertex[forward]
        second_vertex = terminal.origin_vertex[reverse]
        first_center = geometry.vertex_points[first_vertex]
        second_center = geometry.vertex_points[second_vertex]
        first_incidence = geometry.incidence_points[forward]
        second_incidence = geometry.incidence_points[reverse]
        edge_center = geometry.edge_points[edge]
        push!(skeleton, (first=first_center, second=first_incidence,
                         incident=Set([first_vertex])))
        push!(skeleton, (first=first_incidence, second=edge_center, incident=Set{Int}()))
        push!(skeleton, (first=edge_center, second=second_incidence, incident=Set{Int}()))
        push!(skeleton, (first=second_incidence, second=second_center,
                         incident=Set([second_vertex])))
    end

    punctures = Puncture[]
    for boundary in cut.boundaries
        center = geometry.vertex_points[boundary.id]
        clearance = Inf
        for other in cut.boundaries
            other.id == boundary.id && continue
            point = geometry.vertex_points[other.id]
            clearance = min(clearance, hypot(point[1] - center[1], point[2] - center[2]) / 2)
        end
        for segment in skeleton
            boundary.id in segment.incident && continue
            clearance = min(clearance,
                sqrt(_segment_distance_squared(center, segment.first, segment.second)))
        end
        for port in terminal.rotations[boundary.id]
            dart = terminal.origin_dart[port]
            point = geometry.incidence_points[dart]
            clearance = min(clearance, hypot(point[1] - center[1], point[2] - center[2]))
        end
        isfinite(clearance) || (clearance = 1.0)
        clearance > 1e-8 ||
            throw(ArgumentError("Tutte drawing has no positive puncture clearance at boundary $(boundary.id)"))
        push!(punctures, Puncture(boundary.id, boundary.pair, boundary.side,
                                  center, 0.18clearance, boundary.orientation))
    end
    punctures
end

function _ports_on_terminal_rays(base_ports, terminal, geometry, punctures)
    puncture_map = Dict(puncture.id => puncture for puncture in punctures)
    [begin
        puncture = puncture_map[port.boundary_id]
        dart = terminal.origin_dart[port.id]
        incidence = geometry.incidence_points[dart]
        dx = incidence[1] - puncture.center[1]
        dy = incidence[2] - puncture.center[2]
        length = hypot(dx, dy)
        length > 1e-10 ||
            throw(ArgumentError("terminal incidence ray collapses at port $(port.id)"))
        point = (puncture.center[1] + puncture.radius * dx / length,
                 puncture.center[2] + puncture.radius * dy / length)
        Port(port.id, port.intersection_id, port.boundary_id, port.beta_id,
             port.parameter, point)
    end for port in base_ports]
end

function _routes_on_terminal_refinement(template_routes, terminal, geometry, ports)
    port_map = Dict(port.id => port for port in ports)
    length(template_routes) == terminal.edges ||
        throw(ArgumentError("terminal edge count does not match beta successor routes"))
    [begin
        forward, reverse = 2edge - 1, 2edge
        route.start_port == terminal.origin_port[forward] ||
            throw(ArgumentError("terminal forward dart does not match beta successor route $edge"))
        route.end_port == terminal.origin_port[reverse] ||
            throw(ArgumentError("terminal reverse dart does not match beta successor route $edge"))
        Route(route.family, route.curve_id, route.segment_id,
              route.start_port, route.end_port,
              [port_map[route.start_port].point,
               geometry.incidence_points[forward],
               geometry.edge_points[edge],
               geometry.incidence_points[reverse],
               port_map[route.end_port].point])
    end for (edge, route) in enumerate(template_routes)]
end

function _realize_terminal_map(H, cut, base_ports, terminal; fixed_placement=nothing)
    geometry = _solve_terminal_coordinates(terminal; fixed_placement=fixed_placement)
    geometry.refinement_crossings == 0 ||
        throw(ArgumentError("Tutte refinement has $(geometry.refinement_crossings) nonincident edge crossings"))
    punctures = fixed_placement === nothing ?
        _computed_terminal_punctures(cut, terminal, geometry) : collect(fixed_placement)
    ports = _ports_on_terminal_rays(base_ports, terminal, geometry, punctures)
    template_routes = _alpha_cut_routes(H, cut, punctures, ports)
    routes = _routes_on_terminal_refinement(template_routes, terminal, geometry, ports)
    trace = Dict{String,Any}(
        "algorithm" => "tutte-barycentric-flag-refinement",
        "components" => terminal.components,
        "refinement_nodes" => geometry.refinement_nodes,
        "refinement_edges" => geometry.refinement_edges,
        "refinement_triangles" => geometry.refinement_triangles,
        "refinement_crossings" => geometry.refinement_crossings,
    )
    punctures, ports, routes, trace
end
