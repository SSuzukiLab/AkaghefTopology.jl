"""Small semantic SVG/TikZ renderers with no topology dependency on CairoMakie."""

_svg_escape(value) = replace(string(value), "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\""=>"&quot;")
_fmt(value) = string(round(Float64(value); digits=6))

function _svg_point(point)
    "$( _fmt(point[1]) ),$( _fmt(-point[2]) )"
end

function _layout_json(layout::HeegaardLayout)
    Dict{String,Any}(
        "schema" => "akaghef.heegaard.layout/1",
        "diagram_hash" => layout.diagram_hash,
        "view" => String(layout.view),
        "chart" => deepcopy(layout.chart),
        "punctures" => [Dict("id"=>p.id, "pair"=>p.pair, "side"=>String(p.side),
                             "center"=>[p.center...], "radius"=>p.radius, "orientation"=>p.orientation)
                        for p in layout.punctures],
        "ports" => [Dict("id"=>p.id, "intersection_id"=>p.intersection_id,
                         "boundary_id"=>p.boundary_id, "beta_id"=>p.beta_id,
                         "parameter"=>p.parameter, "point"=>[p.point...]) for p in layout.ports],
        "routes" => [Dict("family"=>String(r.family), "curve_id"=>r.curve_id,
                          "segment_id"=>r.segment_id,
                          "start_port"=>r.start_port, "end_port"=>r.end_port,
                          "points"=>[[point...] for point in r.points],
                          "curve"=>_curve_json(r.curve)) for r in layout.routes],
        "crossings" => [[pair...] for pair in layout.crossings],
        "backend_trace" => deepcopy(layout.backend_trace),
        "report" => Dict("valid"=>layout.report.valid,
                          "checks"=>Dict(String(key)=>value for (key,value) in layout.report.checks)),
    )
end

function _curve_json(curve::PolylineCurve)
    Dict{String,Any}("type"=>"polyline", "points"=>[[point...] for point in curve.points])
end

function _curve_json(curve::BezierChain)
    Dict{String,Any}(
        "type"=>"bezier_chain", "flatness"=>curve.flatness, "max_depth"=>curve.max_depth,
        "segments"=>[Dict("type"=>"cubic", "points"=>[[point...] for point in
                       (segment.p0, segment.p1, segment.p2, segment.p3)])
                     for segment in curve.segments])
end

function _curve_json(curve::RouteCurve)
    throw(ArgumentError("unsupported route curve representation $(typeof(curve))"))
end

function save_layout(path::AbstractString, layout::HeegaardLayout)
    endswith(path, ".hgd.layout.json") || throw(ArgumentError("layout path must end in .hgd.layout.json"))
    open(path, "w") do io
        write(io, _canonical_json(_layout_json(layout)))
        write(io, '\n')
    end
    path
end

function _point_tuple(value)
    values = collect(value)
    length(values) == 2 || throw(ArgumentError("2D coordinate must have two entries"))
    (Float64(values[1]), Float64(values[2]))
end

function _curve_from_spec(spec)
    spec === nothing && throw(ArgumentError("route is missing its curve representation"))
    curve_type = string(_spec_field(spec, :type, ""))
    if curve_type == "polyline"
        points = [_point_tuple(point) for point in _as_vector(_spec_field(spec, :points, nothing))]
        return PolylineCurve(points)
    elseif curve_type == "bezier_chain"
        segments = CubicBezierSegment[]
        for segment_spec in _as_vector(_spec_field(spec, :segments, nothing))
            string(_spec_field(segment_spec, :type, "")) == "cubic" ||
                throw(ArgumentError("unknown Bezier segment type"))
            points = [_point_tuple(point) for point in
                      _as_vector(_spec_field(segment_spec, :points, nothing))]
            length(points) == 4 || throw(ArgumentError("cubic Bezier segment must have four points"))
            push!(segments, CubicBezierSegment((points[1], points[2], points[3], points[4])))
        end
        return BezierChain(segments;
                           flatness=Float64(_spec_field(spec, :flatness, 0.0)),
                           max_depth=Int(_spec_field(spec, :max_depth, 0)))
    end
    throw(ArgumentError("unknown route curve type: $curve_type"))
end

"""Read a layout only when its stored core hash matches `H`."""
function load_layout(path::AbstractString, H::HeegaardDiagram)
    parsed = JSON3.read(read(path, String))
    schema = string(_spec_field(parsed, :schema, ""))
    schema == "akaghef.heegaard.layout/1" || throw(ArgumentError("unknown layout schema: $schema"))
    hash = string(_spec_field(parsed, :diagram_hash, ""))
    core_hash = canonical_hash(H)
    hash == core_hash || throw(ArgumentError("layout core hash does not match the diagram"))
    punctures = [begin
        center = _point_tuple(_spec_field(item, :center, []))
        Puncture(Int(_spec_field(item, :id, 0)), Int(_spec_field(item, :pair, 0)),
                 Symbol(_spec_field(item, :side, "plus")), center,
                 Float64(_spec_field(item, :radius, 0.0)), Int(_spec_field(item, :orientation, 1)))
    end for item in _as_vector(_spec_field(parsed, :punctures, []))]
    ports = [begin
        point = _point_tuple(_spec_field(item, :point, []))
        Port(Int(_spec_field(item, :id, 0)), Int(_spec_field(item, :intersection_id, 0)),
             Int(_spec_field(item, :boundary_id, 0)), Int(_spec_field(item, :beta_id, 0)),
             Float64(_spec_field(item, :parameter, 0.0)), point)
    end for item in _as_vector(_spec_field(parsed, :ports, []))]
    routes = [begin
        encoded_points = _spec_field(item, :points, nothing)
        encoded_points === nothing && throw(ArgumentError("layout route is missing flattened points"))
        points = [_point_tuple(point) for point in _as_vector(encoded_points)]
        start = _spec_field(item, :start_port, nothing)
        finish = _spec_field(item, :end_port, nothing)
        encoded_curve = _spec_field(item, :curve, nothing)
        curve = encoded_curve === nothing ? nothing : _curve_from_spec(encoded_curve)
        Route(Symbol(_spec_field(item, :family, "beta")), Int(_spec_field(item, :curve_id, 0)),
              Int(_spec_field(item, :segment_id, 0)), start === nothing ? nothing : Int(start),
              finish === nothing ? nothing : Int(finish), points; curve)
    end for item in _as_vector(_spec_field(parsed, :routes, []))]
    crossings = [(begin
        values = collect(pair)
        length(values) == 2 || throw(ArgumentError("layout crossing must contain two route ids"))
        (Int(values[1]), Int(values[2]))
    end) for pair in _as_vector(_spec_field(parsed, :crossings, []))]
    layout = HeegaardLayout(hash, Symbol(_spec_field(parsed, :view, "alpha_cut")),
        _plain_value(_spec_field(parsed, :chart, Dict{String,Any}())), punctures, ports, routes,
        crossings, _plain_value(_spec_field(parsed, :backend_trace, Dict{String,Any}())),
        LayoutValidationReport(false, ValidationIssue[], Dict{Symbol,Bool}()))
    report = validate_layout(H, layout; expected_hash=core_hash)
    HeegaardLayout(layout.diagram_hash, layout.view, layout.chart, layout.punctures,
                   layout.ports, layout.routes, layout.crossings, layout.backend_trace, report)
end

function _svg_header(layout)
    provenance = haskey(layout.backend_trace, "provenance") ?
        _svg_escape(_canonical_json(layout.backend_trace["provenance"])) : ""
    coordinates = NTuple{2,Float64}[]
    for route in layout.routes
        append!(coordinates, [(point[1], -point[2]) for point in route.points])
    end
    for puncture in layout.punctures
        push!(coordinates, (puncture.center[1] - puncture.radius,
                            -puncture.center[2] - puncture.radius))
        push!(coordinates, (puncture.center[1] + puncture.radius,
                            -puncture.center[2] + puncture.radius))
    end
    isempty(coordinates) && push!(coordinates, (0.0, 0.0))
    xmin, xmax = extrema(first.(coordinates))
    ymin, ymax = extrema(last.(coordinates))
    span = max(xmax - xmin, ymax - ymin, 1.0)
    margin = 0.08span + 0.15
    viewbox = join(_fmt.((xmin - margin, ymin - margin,
                         xmax - xmin + 2margin, ymax - ymin + 2margin)), " ")
    "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$viewbox\" " *
    "role=\"img\" data-heegaard-view=\"$(_svg_escape(layout.view))\" " *
    "data-diagram-hash=\"$(_svg_escape(layout.diagram_hash))\">" *
    "<title>Heegaard diagram — paired punctures</title>" *
    "<desc>Alpha cut: paired puncture disks in two rows; beta arcs preserve their ordered endpoints.</desc>" *
    "<metadata data-provenance=\"$provenance\"/>"
end

function _route_path(route::Route)
    curve = route.curve
    if curve isa PolylineCurve
        points = curve.points
        result = "M $(_fmt(first(points)[1])) $(_fmt(-first(points)[2]))"
        for point in points[2:end]
            result *= " L $(_fmt(point[1])) $(_fmt(-point[2]))"
        end
        return result
    elseif curve isa BezierChain
        first_segment = first(curve.segments)
        result = "M $(_fmt(first_segment.p0[1])) $(_fmt(-first_segment.p0[2]))"
        for segment in curve.segments
            result *= " C $(_fmt(segment.p1[1])) $(_fmt(-segment.p1[2])) " *
                      "$(_fmt(segment.p2[1])) $(_fmt(-segment.p2[2])) " *
                      "$(_fmt(segment.p3[1])) $(_fmt(-segment.p3[2]))"
        end
        return result
    end
    throw(ArgumentError("unsupported route curve representation $(typeof(curve))"))
end

function _render_svg(layout::HeegaardLayout)
    parts = String[_svg_header(layout), "<g data-layer=\"routes\">"]
    for route in layout.routes
        attrs = "data-curve-family=\"$(_svg_escape(route.family))\" " *
                "data-curve-id=\"$(route.curve_id)\" data-route-segment=\"$(route.segment_id)\""
        palette = ("#2459a9", "#b75029", "#27806b", "#8650a1", "#986d15")
        color = palette[mod1(route.curve_id, length(palette))]
        push!(parts, "<path $attrs d=\"$(_route_path(route))\" fill=\"none\" stroke=\"$color\" stroke-width=\"1.6\" stroke-linecap=\"round\" stroke-linejoin=\"round\" vector-effect=\"non-scaling-stroke\"/>")
    end
    push!(parts, "</g>")
    if layout.view == :alpha_cut
        push!(parts, "<g data-layer=\"punctures\">")
        for puncture in layout.punctures
            push!(parts, "<circle data-pair=\"$(puncture.pair)\" data-side=\"$(_svg_escape(puncture.side))\" " *
                         "cx=\"$(_fmt(puncture.center[1]))\" cy=\"$(_fmt(-puncture.center[2]))\" " *
                         "r=\"$(_fmt(puncture.radius))\" fill=\"#f5f3ee\" " *
                         "stroke=\"#383c43\" stroke-width=\"1.4\" vector-effect=\"non-scaling-stroke\"/>")
            label = "α$(puncture.pair)" * (puncture.side == :plus ? "+" : "−")
            push!(parts, "<text x=\"$(_fmt(puncture.center[1]))\" y=\"$(_fmt(-puncture.center[2]))\" text-anchor=\"middle\" dominant-baseline=\"central\" font-family=\"Georgia,serif\" font-size=\"$(_fmt(0.45puncture.radius))\" fill=\"#25282d\">$label</text>")
        end
        for port in layout.ports
            puncture = only(p for p in layout.punctures if p.id == port.boundary_id)
            labelpoint = _vadd(puncture.center, _vscale(0.73, _vsub(port.point,puncture.center)))
            push!(parts, "<circle data-intersection-id=\"$(port.intersection_id)\" data-port-index=\"$(port.id)\" " *
                         "cx=\"$(_fmt(port.point[1]))\" cy=\"$(_fmt(-port.point[2]))\" r=\"$(_fmt(0.048puncture.radius))\" fill=\"#383c43\"><title>Intersection $(port.intersection_id) · β$(port.beta_id)</title></circle>")
            push!(parts, "<text data-intersection-label=\"$(port.intersection_id)\" x=\"$(_fmt(labelpoint[1]))\" y=\"$(_fmt(-labelpoint[2]))\" text-anchor=\"middle\" dominant-baseline=\"central\" font-family=\"system-ui,sans-serif\" font-size=\"$(_fmt(0.21puncture.radius))\" fill=\"#535963\">$(port.intersection_id)</text>")
        end
        push!(parts, "</g>")
    end
    push!(parts, "</svg>")
    join(parts, "\n")
end

"""Render an SVG string; pass `path=` to persist it as an artifact."""
function render(layout::HeegaardLayout; backend=:cairo, path=nothing)
    backend in (:cairo, :svg, :pure_julia) || throw(ArgumentError("unsupported renderer backend $backend"))
    svg = _render_svg(layout)
    path === nothing || open(path, "w") do io
        write(io, svg)
    end
    svg
end

"""Return a QMD-displayable layout without changing an already computed plot."""
plot(layout::HeegaardLayout) = layout

"""Compute a fixed-profile alpha-cut plot. The input layout is never mutated.

Spacing, interaction range, tautness and disk clearance are fixed at their
production setting; no browser-only deformation or tuning controls are used.
Use view=:surface to show both curve families in a punctured-sphere/annulus atlas.
"""
function plot(H::HeegaardDiagram, feasible::HeegaardLayout; view=:alpha_cut)
    view in (:alpha_cut,:surface) || throw(ArgumentError("unsupported plot view: $view"))
    feasible.view == :alpha_cut || throw(ArgumentError("plot requires an alpha-cut layout"))
    validate_layout(H,feasible).valid || throw(ArgumentError("plot requires a valid layout of this diagram"))
    result = isempty(feasible.routes) ? feasible : _joint_plot_relaxation(H,feasible)
    view == :surface ? surface_chart(H,result) : result
end

"""Plot a surface using combinatorial planar realization, or supply a valid layout.

Routing failure is reported, never replaced with a crossing-filled drawing.
"""
plot(H::HeegaardDiagram; view=:alpha_cut, router=view==:surface ? :planar : :auto) =
    plot(H, layout_heegaard(H;backend=:pure_julia,router); view)

"""Display a Heegaard layout inline in Quarto, IJulia, and other rich frontends."""
function Base.show(io::IO, ::MIME"image/svg+xml", layout::HeegaardLayout)
    print(io, _render_svg(layout))
end

function render_tikz(layout::HeegaardLayout)
    parts = ["\\begin{tikzpicture}[scale=2]", "\\tikzset{hgdroute/.style={draw=blue!60!black,line width=0.35pt}}"]
    for puncture in layout.punctures
        push!(parts, "\\draw[red!70!black] ($( _fmt(puncture.center[1]) ),$( _fmt(puncture.center[2]) )) circle ($( _fmt(puncture.radius) ));")
    end
    for route in layout.routes
        curve = route.curve
        if curve isa PolylineCurve
            points = join(["($(_fmt(point[1])),$(_fmt(point[2])))" for point in curve.points], " -- ")
            push!(parts, "\\draw[hgdroute] $points;")
        elseif curve isa BezierChain
            segment = first(curve.segments)
            path = "($(_fmt(segment.p0[1])),$(_fmt(segment.p0[2])))"
            for segment in curve.segments
                path *= " .. controls ($(_fmt(segment.p1[1])),$(_fmt(segment.p1[2]))) and " *
                        "($(_fmt(segment.p2[1])),$(_fmt(segment.p2[2]))) .. " *
                        "($(_fmt(segment.p3[1])),$(_fmt(segment.p3[2])))"
            end
            push!(parts, "\\draw[hgdroute] $path;")
        else
            throw(ArgumentError("unsupported route curve representation $(typeof(curve))"))
        end
    end
    push!(parts, "\\end{tikzpicture}")
    join(parts, "\n")
end

function inspect_html(layout::HeegaardLayout)
    "<!doctype html><meta charset=\"utf-8\"><title>Heegaard layout</title>" * _render_svg(layout)
end
