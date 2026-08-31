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
    "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"-3 -3 6 6\" " *
    "role=\"img\" data-heegaard-view=\"$(_svg_escape(layout.view))\" " *
    "data-diagram-hash=\"$(_svg_escape(layout.diagram_hash))\">" *
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
        push!(parts, "<path $attrs d=\"$(_route_path(route))\" fill=\"none\" stroke=\"#17324d\" stroke-width=\"0.025\"/>")
    end
    push!(parts, "</g>")
    if layout.view == :alpha_cut
        push!(parts, "<g data-layer=\"punctures\">")
        for puncture in layout.punctures
            push!(parts, "<circle data-pair=\"$(puncture.pair)\" data-side=\"$(_svg_escape(puncture.side))\" " *
                         "cx=\"$(_fmt(puncture.center[1]))\" cy=\"$(_fmt(-puncture.center[2]))\" " *
                         "r=\"$(_fmt(puncture.radius))\" fill=\"none\" stroke=\"#b54a35\" stroke-width=\"0.025\"/>")
        end
        for port in layout.ports
            push!(parts, "<circle data-intersection-id=\"$(port.intersection_id)\" data-port-index=\"$(port.id)\" " *
                         "cx=\"$(_fmt(port.point[1]))\" cy=\"$(_fmt(-port.point[2]))\" r=\"0.045\" fill=\"#b54a35\"/>")
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
