"""The deterministic alpha-cut chart and ordered-port data."""

function cut_along_alpha(H::HeegaardDiagram; diagram_hash=canonical_hash(H))
    compression_boundary_report(H).alpha_complete || throw(ArgumentError(
        "alpha-cut coordinates require a complete alpha cut system; use the explicit auxiliary marked cuts otherwise"))
    code = H.code.code
    boundaries = CutBoundary[]
    orders = Dict{Int,Vector{Int}}()
    beta_orders = Dict{Int,Vector{Int}}(curve.id => copy(curve.visits) for curve in code.beta)
    signs = Dict(point.id => point.sign for point in code.intersections)
    betas = Dict(point.id => point.beta for point in code.intersections)
    for (pair, curve) in enumerate(code.alpha)
        push!(boundaries, CutBoundary(2pair - 1, pair, :plus, curve.id, 1))
        push!(boundaries, CutBoundary(2pair, pair, :minus, curve.id, -1))
        orders[curve.id] = copy(curve.visits)
    end
    CutSurface(string(diagram_hash), length(code.alpha), boundaries, orders, beta_orders, signs, betas)
end

function default_puncture_placement(cut::CutSurface)
    g = cut.genus
    [Puncture(boundary.id, boundary.pair, boundary.side,
              (Float64(boundary.pair), boundary.side == :plus ? 1.0 : -1.0),
              0.28, boundary.orientation) for boundary in cut.boundaries]
end

function _placement_map(placement)
    Dict(p.id => p for p in placement)
end

function place_ports(cut::CutSurface, placement=default_puncture_placement(cut))
    punctures = _placement_map(placement)
    length(punctures) == length(cut.boundaries) || throw(ArgumentError("placement must contain every cut boundary"))
    ports = Port[]
    next_port_id = 1
    for (pair, plus_boundary_id) in enumerate(1:2:length(cut.boundaries))
        alpha_curve_id = cut.boundaries[plus_boundary_id].curve_id
        visits = get(cut.alpha_orders, alpha_curve_id, Int[])
        # Cutting an alpha curve splits every alpha/beta intersection into two
        # boundary half-edges.  The minus boundary has the opposite induced
        # orientation, hence its reversed cyclic order.
        for (points, boundary_id) in ((visits, 2pair - 1), (reverse(visits), 2pair))
            for (rank, point) in enumerate(points)
                puncture = punctures[boundary_id]
                parameter = isempty(points) ? 0.5 : (rank - 0.5) / length(points)
                angle = 2π * parameter
                coordinate = (puncture.center[1] + puncture.radius * cos(angle),
                              puncture.center[2] + puncture.radius * sin(angle))
                push!(ports, Port(next_port_id, point, boundary_id, cut.intersection_betas[point],
                                  parameter, coordinate))
                next_port_id += 1
            end
        end
    end
    sort!(ports; by=port -> (port.boundary_id, port.parameter, port.intersection_id))
    ports
end

function _halfedge_port_map(cut::CutSurface, ports::Vector{Port})
    boundaries = Dict(boundary.id => boundary for boundary in cut.boundaries)
    result = Dict{Tuple{Int,Symbol},Int}()
    for port in ports
        boundary = get(boundaries, port.boundary_id, nothing)
        boundary === nothing && throw(ArgumentError("port $(port.id) references unknown cut boundary $(port.boundary_id)"))
        sign = get(cut.intersection_signs, port.intersection_id, 0)
        sign in (-1, 1) || throw(ArgumentError("port $(port.id) references an unknown intersection sign"))
        # For a positive intersection beta_out exits on alpha+ and beta_in on
        # alpha-.  A negative intersection exchanges the two sides.
        role = (boundary.side == :plus) == (sign == 1) ? :out : :in
        key = (port.intersection_id, role)
        haskey(result, key) && throw(ArgumentError("duplicate $role half-edge port for intersection $(port.intersection_id)"))
        result[key] = port.id
    end
    result
end

place_ports(H::HeegaardDiagram, cut::CutSurface,
            placement=default_puncture_placement(cut)) = place_ports(cut, placement)

function build_wiring_problem(cut::CutSurface, ports::Vector{Port})
    port_map = _halfedge_port_map(cut, ports)
    wires = NamedTuple{(:beta_id, :from_port, :to_port),Tuple{Int,Int,Int}}[]
    for (beta_id, order) in cut.beta_orders
        # A generic CutSurface retains beta cyclic words as ordered wires.
        isempty(order) && continue
        for index in eachindex(order)
            from = get(port_map, (order[index], :out), 0)
            to = get(port_map, (order[mod1(index + 1, length(order))], :in), 0)
            from == 0 || to == 0 || push!(wires, (beta_id=beta_id, from_port=from, to_port=to))
        end
    end
    WiringProblem(cut.diagram_hash, sort(collect(keys(cut.alpha_orders))), ports, wires)
end

function build_wiring_problem(H::HeegaardDiagram, cut::CutSurface, ports::Vector{Port})
    port_map = _halfedge_port_map(cut, ports)
    wires = NamedTuple{(:beta_id, :from_port, :to_port),Tuple{Int,Int,Int}}[]
    for curve in H.code.code.beta
        isempty(curve.visits) && continue
        for index in eachindex(curve.visits)
            from = port_map[(curve.visits[index], :out)]
            to = port_map[(curve.visits[mod1(index + 1, length(curve.visits))], :in)]
            push!(wires, (beta_id=curve.id, from_port=from, to_port=to))
        end
    end
    WiringProblem(cut.diagram_hash, [curve.id for curve in H.code.code.beta], ports, wires)
end
