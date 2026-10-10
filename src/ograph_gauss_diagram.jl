"""Return an SVG Gauss diagram for an unweighted o-graph presentation.

The cyclic Gauss word is placed clockwise on an oriented circle. Each vertex
is drawn as an arrow from its positive (over-passage) occurrence to its
negative (under-passage) occurrence. The arrow label records the vertex label
and vertex sign separately from those passage signs. This is a presentation
view; it does not certify closed-normality or classical realizability.
"""
_gauss_svg_attribute(value) = replace(_esc(value), "\"" => "&quot;", "'" => "&apos;")

function gauss_diagram_svg(state::OGraphState)
    code = collect(state.gauss)
    signs = collect(state.vertex_signs)
    vertex_count = length(signs)
    occurrence_count = length(code)

    radius = max(180.0, 28.0 * vertex_count)
    margin = 70.0
    header = 46.0
    width = 2radius + 2margin
    height = 2radius + 2margin + header
    center_x = width / 2
    center_y = margin + header / 2 + radius
    point_radius = 7.0
    endpoint_label_offset = 25.0

    point_at(index) = begin
        angle = -pi / 2 + 2pi * (index - 1) / max(occurrence_count, 1)
        (center_x + radius * cos(angle), center_y + radius * sin(angle), angle)
    end
    fmt(x) = string(round(Float64(x); digits=3))

    parts = String[]
    push!(parts, "<svg xmlns=\"http://www.w3.org/2000/svg\" class=\"gauss\" width=\"$(fmt(width))\" height=\"$(fmt(height))\" viewBox=\"0 0 $(fmt(width)) $(fmt(height))\" role=\"img\" aria-labelledby=\"title desc\">")
    push!(parts, "<title id=\"title\">Gauss diagram of o-graph</title>")
    push!(parts, "<desc id=\"desc\">The cyclic Gauss word is placed clockwise on an oriented circle. Arrows go from over-passage to under-passage. Endpoint labels show passage signs; arrow labels show vertex labels and vertex signs.</desc>")
    push!(parts, """<style>
.gauss{--ground:#faf9f7;--strand:#1b3bd0;--ink:#12161d}
@media (prefers-color-scheme:dark){.gauss{--ground:#0e1116;--strand:#8ba2ff;--ink:#e6e9ef}}
.gauss .bg{fill:var(--ground)}.gauss .circle{fill:none;stroke:var(--strand)}.gauss .chord{fill:none}
.gauss .endpoint.over{fill:var(--strand)}.gauss .endpoint.under{fill:var(--ground);stroke:var(--strand)}
.gauss .label{font:16px sans-serif;fill:var(--ink);stroke:var(--ground);stroke-width:5px;paint-order:stroke}
.gauss .vertex-label{font:14px sans-serif;font-weight:600;fill:var(--ink);stroke:var(--ground);stroke-width:5px;paint-order:stroke}
.gauss .title{font:20px sans-serif;font-weight:600;fill:var(--ink)}
</style>""")
    push!(parts, "<rect class=\"bg\" width=\"100%\" height=\"100%\"/>")
    push!(parts, "<text class=\"title\" x=\"$(fmt(center_x))\" y=\"32\" text-anchor=\"middle\">Gauss diagram of o-graph</text>")
    push!(parts, "<circle class=\"circle\" cx=\"$(fmt(center_x))\" cy=\"$(fmt(center_y))\" r=\"$(fmt(radius))\" stroke-width=\"2.5\"/>")

    # Mark the increasing-index direction along the circle independently of
    # the over-to-under direction on each vertex arrow.
    direction_angle = -pi / 4
    direction_point = (center_x + radius * cos(direction_angle),
                       center_y + radius * sin(direction_angle))
    tangent = (-sin(direction_angle), cos(direction_angle))
    arrow_base = (direction_point[1] - 10tangent[1],
                  direction_point[2] - 10tangent[2])
    normal = (-tangent[2], tangent[1])
    arrow_left = (arrow_base[1] + 5normal[1], arrow_base[2] + 5normal[2])
    arrow_right = (arrow_base[1] - 5normal[1], arrow_base[2] - 5normal[2])
    push!(parts, "<path d=\"M $(fmt(arrow_left[1])) $(fmt(arrow_left[2])) L $(fmt(direction_point[1])) $(fmt(direction_point[2])) L $(fmt(arrow_right[1])) $(fmt(arrow_right[2]))\" fill=\"none\" stroke=\"var(--strand)\" stroke-width=\"2.5\" stroke-linecap=\"round\"/>")

    positions = Dict{Int,Tuple{Float64,Float64,Float64}}()
    for index in eachindex(code)
        positions[index] = point_at(index)
    end

    for vertex in 1:vertex_count
        over_index = findfirst(==(vertex), code)
        under_index = findfirst(==(-vertex), code)
        over_index === nothing && throw(ArgumentError("vertex $vertex has no over-passage"))
        under_index === nothing && throw(ArgumentError("vertex $vertex has no under-passage"))
        over_x, over_y, _ = positions[over_index]
        under_x, under_y, _ = positions[under_index]
        midpoint_x = (over_x + under_x) / 2
        midpoint_y = (over_y + under_y) / 2
        end_dx = under_x - over_x
        end_dy = under_y - over_y
        end_norm = hypot(end_dx, end_dy)
        end_norm == 0 && (end_norm = 1.0)
        ux, uy = end_dx / end_norm, end_dy / end_norm
        nx, ny = -uy, ux
        tip_x, tip_y = under_x - 10ux, under_y - 10uy
        base_x, base_y = tip_x - 9ux, tip_y - 9uy
        head_left = (base_x + 5ny, base_y + 5ny)
        head_right = (base_x - 5ny, base_y - 5ny)
        color = signs[vertex] == 1 ? "#2853b8" : "#a13e32"
        push!(parts, "<path class=\"chord\" d=\"M $(fmt(over_x)) $(fmt(over_y)) L $(fmt(tip_x)) $(fmt(tip_y))\" stroke=\"$color\" stroke-width=\"2.4\"/>")
        push!(parts, "<path d=\"M $(fmt(head_left[1])) $(fmt(head_left[2])) L $(fmt(tip_x)) $(fmt(tip_y)) L $(fmt(head_right[1])) $(fmt(head_right[2])) Z\" fill=\"$color\"/>")

        offset = 15.0 + 9.0 * mod(vertex - 1, 3)
        side = isodd(vertex) ? 1.0 : -1.0
        label_x = midpoint_x + side * nx * offset
        label_y = midpoint_y + side * ny * offset
        sign_text = signs[vertex] == 1 ? "+" : "−"
        label_text = "v$vertex · σ=$sign_text"
        push!(parts, "<text class=\"vertex-label\" x=\"$(fmt(label_x))\" y=\"$(fmt(label_y))\" text-anchor=\"middle\" dominant-baseline=\"middle\">$(_esc(label_text))</text>")
    end

    for index in eachindex(code)
        value = code[index]
        x, y, angle = positions[index]
        passage = value > 0 ? "over" : "under"
        passage_text = value > 0 ? "+$(abs(value))" : "−$(abs(value))"
        endpoint_class = value > 0 ? "over" : "under"
        label_x = center_x + (radius + endpoint_label_offset) * cos(angle)
        label_y = center_y + (radius + endpoint_label_offset) * sin(angle)
        push!(parts, "<circle class=\"endpoint $endpoint_class\" cx=\"$(fmt(x))\" cy=\"$(fmt(y))\" r=\"$(fmt(point_radius))\" stroke-width=\"2\" aria-label=\"$(_esc(passage)) passage at vertex $(abs(value))\"/>")
        push!(parts, "<text class=\"label\" x=\"$(fmt(label_x))\" y=\"$(fmt(label_y))\" text-anchor=\"middle\" dominant-baseline=\"middle\">$passage_text</text>")
    end

    convention = _gauss_svg_attribute(state.convention_id)
    gauss = join(code, ",")
    vertex_signs = join(signs, ",")
    push!(parts, "<metadata><ograph convention=\"$convention\" gauss=\"$gauss\" vertex-signs=\"$vertex_signs\"/></metadata>")
    push!(parts, "</svg>")
    join(parts, "\n")
end

gauss_diagram_svg(state::ClosedNormalOGraph) = gauss_diagram_svg(state.state)
gauss_diagram_svg(link::VirtualLink; convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) =
    gauss_diagram_svg(ograph_state(link; convention_id=convention_id))

"""Write an o-graph Gauss diagram SVG and return its absolute path."""
function plot_gauss_diagram_svg(state, path::AbstractString)
    svg = gauss_diagram_svg(state)
    open(path, "w") do io
        write(io, svg)
    end
    abspath(path)
end

Base.show(io::IO, ::MIME"image/svg+xml", state::OGraphState) = print(io, gauss_diagram_svg(state))
Base.show(io::IO, ::MIME"image/svg+xml", state::ClosedNormalOGraph) = print(io, gauss_diagram_svg(state))
