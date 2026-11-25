"""The deterministic alpha-cut chart and ordered-port data."""

"""
    cut_along_alpha(H::CellularHeegaardDiagram; diagram_hash=canonical_hash(H)) -> CutSurface

Return the paired boundaries and beta incidences of a legacy Heegaard surface cut along alpha.

# Arguments
`H` is a certified cellular diagram whose alpha curves form the complete
cut system of [the Heegaard domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
The compression-boundary report checks completeness; failure throws
`ArgumentError` explaining that explicit auxiliary cuts are required.
`diagram_hash` is converted with `string`; by default it is computed by
`canonical_hash`, whose errors propagate. An override is trusted as an identity
label and is not checked against `H`.

# Returns
`CutSurface` with the identity label, genus, paired boundary records,
copied alpha and beta visit orders, crossing signs and crossing-to-beta map.
Their numbering is specified in [legacy cut conventions](../specs/heegaard-terminology.md#legacy-cut-indices).

# Guarantees
1. For a valid complete alpha system, the result describes its cut surface
   and the beta endpoints on the resulting paired boundaries.
   [proved/agent: by construction]
2. There are twice as many boundaries as alpha curves. The input is not
   modified and its visit vectors are copied. [proved/agent: by construction]

# Not guaranteed
This supplies no planar placement or embedded beta arcs. A supplied identity
label is not a certificate, and later edits of the input are not synchronized
with the returned cut.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 and lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — complete alpha-cut surface and its paired boundaries.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> c = cut_along_alpha(d); (length(default_puncture_placement(c)), length(place_ports(c)))
(2, 2)
```

See also [`default_puncture_placement`](@ref), [`place_ports`](@ref), [`build_wiring_problem`](@ref).

Role: convert
"""
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

"""
    default_puncture_placement(cut::CutSurface) -> Vector{Puncture}

Return the default circular display placements for a legacy alpha-cut surface's boundaries.

# Arguments
`cut` normally comes from `cut_along_alpha`. The caller must supply consistent
boundary records under [legacy cut conventions](../specs/heegaard-terminology.md#legacy-cut-indices);
no boundary, genus or pair-consistency check is performed here.

# Returns
A fresh `Vector{Puncture}` in stored boundary order, carrying each boundary's
ID, pair, side and orientation with its default center and radius.

# Guarantees
1. Each stored cut boundary has one placement in the
   [legacy working chart](../specs/heegaard-terminology.md#legacy-cut-indices).
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not place beta arcs, certify an embedding, or identify the legacy
chart with the current LC planar convention. Malformed boundary pairs can
produce coincident punctures without an error.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 — legacy puncture working chart.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> length(default_puncture_placement(cut_along_alpha(d)))
2
```

See also [`cut_along_alpha`](@ref), [`place_ports`](@ref).

Role: draw
"""
function default_puncture_placement(cut::CutSurface)
    g = cut.genus
    [Puncture(boundary.id, boundary.pair, boundary.side,
              (Float64(boundary.pair), boundary.side == :plus ? 1.0 : -1.0),
              0.28, boundary.orientation) for boundary in cut.boundaries]
end

function _placement_map(placement)
    Dict(p.id => p for p in placement)
end

"""
    place_ports(cut::CutSurface, placement=default_puncture_placement(cut)) -> Vector{Port}

Place the beta endpoints on the paired boundaries of a legacy alpha-cut surface.

# Arguments
`cut` must have the consistent paired boundaries and incidence dictionaries
returned by `cut_along_alpha`. `placement` must provide one `Puncture` for each
boundary ID. A different number of distinct placement IDs raises
`ArgumentError`; an absent required boundary or beta incidence raises `KeyError`.
Other consistency conditions are caller obligations.

# Returns
A fresh `Vector{Port}` with endpoint IDs, intersection and beta IDs, boundary
parameters and planar coordinates, ordered by boundary and parameter.

# Guarantees
1. For a consistent cut, each alpha visit supplies its two boundary endpoints
   under [the legacy cut convention](../specs/heegaard-terminology.md#legacy-cut-indices).
   These are endpoints of the cut beta curves, not beta arcs.
   [proved/agent: by construction]
2. Coordinates lie on the supplied puncture circles up to floating-point
   rounding. The inputs are not modified. [proved/agent: by construction]

# Not guaranteed
This does not validate a complete cut system, disjoint placement circles,
positive radii, finite coordinates, or a planar embedding of the beta arcs.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 — boundary endpoint interpretation and placement convention.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> length(place_ports(cut_along_alpha(d)))
2
```

See also [`cut_along_alpha`](@ref), [`default_puncture_placement`](@ref), [`build_wiring_problem`](@ref).

Role: draw
"""
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

"""
    place_ports(H::CellularHeegaardDiagram, cut::CutSurface, placement=default_puncture_placement(cut)) -> Vector{Port}

Place the beta endpoints of a supplied legacy cut in a diagram-compatible call form.

# Arguments
`H` is a cellular diagram. The obligations and errors for `cut` and `placement`
are those of `place_ports(cut, placement)`. Their correspondence with `H` is
not checked.

# Returns
The `Vector{Port}` returned by `place_ports(cut, placement)`.

# Guarantees
1. For a consistent cut, the result records its cut beta endpoints under
   [the legacy cut convention](../specs/heegaard-terminology.md#legacy-cut-indices).
   [proved/agent: by construction]
2. Changing `H` alone does not change the result. No input is modified.
   [proved/agent: by construction]

# Not guaranteed
The result is not certified to belong to `H`; the caller must supply its cut.
No embedding of the beta arcs is constructed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 — legacy cut endpoint interpretation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> length(place_ports(d, cut_along_alpha(d)))
2
```

See also [`cut_along_alpha`](@ref), [`build_wiring_problem`](@ref).

Role: draw
"""
place_ports(H::HeegaardDiagram, cut::CutSurface,
            placement=default_puncture_placement(cut)) = place_ports(cut, placement)

"""
    build_wiring_problem(cut::CutSurface, ports::Vector{Port}) -> WiringProblem

Return the beta-arc connections recorded by a legacy alpha-cut surface and its ports.

# Arguments
`cut` and `ports` use the [legacy cut conventions](../specs/heegaard-terminology.md#legacy-cut-indices).
Unknown boundaries, invalid intersection signs and duplicate half-edge
assignments throw `ArgumentError` naming the problem. The caller must supply
all required ports for a complete wiring: a missing endpoint causes that
wire to be omitted. No source identity or placement check is made.

# Returns
`WiringProblem` with the cut's `diagram_hash`, the sorted keys of
`cut.alpha_orders` in `boxes`, the supplied `ports` and `wires`
recording beta IDs with directed endpoint port IDs.

# Guarantees
1. Each wire connects the indicated consecutive beta visits under the
   [half-edge convention](../specs/heegaard-terminology.md#legacy-cut-indices).
   With complete consistent ports this is the beta wiring of the cut surface.
   [proved/agent: by construction]
2. Inputs are not modified; the result shares the port vector.
   [proved/agent: by construction]

# Not guaranteed
Incomplete ports can silently omit wires. Wire order is not canonical.
The result supplies no embedded paths or validation of a 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 — beta half-edge matching.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> c = cut_along_alpha(d); w = build_wiring_problem(c, place_ports(c));

julia> w.wires == [(beta_id=1, from_port=1, to_port=2)]
true
```

See also [`cut_along_alpha`](@ref), [`place_ports`](@ref).

Role: convert
"""
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

"""
    build_wiring_problem(H::CellularHeegaardDiagram, cut::CutSurface, ports::Vector{Port}) -> WiringProblem

Return beta-arc connections using the cellular diagram's visits and supplied alpha-cut ports.

# Arguments
`H`, `cut` and `ports` must describe the same diagram under the
[legacy cut conventions](../specs/heegaard-terminology.md#legacy-cut-indices);
this correspondence is the caller's responsibility. Unknown boundaries,
invalid signs or duplicate half-edge ports throw `ArgumentError`. A missing
endpoint for a diagram visit throws `KeyError`. Hashes and coordinates are
not compared.

# Returns
`WiringProblem` with the cut's `diagram_hash`, the supplied `ports` and
directed beta `wires`. In this overload, the field named `boxes`
contains the diagram's beta curve IDs in their stored order.

# Guarantees
1. For mutually consistent inputs, the wires represent consecutive visits
   along the diagram's beta curves on the alpha-cut surface.
   [proved/agent: by construction]
2. Inputs are not modified; the result shares the port vector.
   [proved/agent: by construction]

# Not guaranteed
The `boxes` field does not give alpha IDs in this overload. Consistency
of the diagram and cut is not checked, and no geometric wiring is constructed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §切断境界と port の番号 — beta half-edge matching.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> c = cut_along_alpha(d);

julia> build_wiring_problem(d, c, place_ports(c)).wires == [(beta_id=1, from_port=1, to_port=2)]
true
```

See also [`cut_along_alpha`](@ref), [`place_ports`](@ref).

Role: convert
"""
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
