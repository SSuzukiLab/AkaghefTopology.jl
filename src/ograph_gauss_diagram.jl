"""Return an SVG Gauss diagram for an unweighted o-graph presentation.

The cyclic Gauss word is placed clockwise on an oriented circle. Each vertex
is drawn as an arrow from its positive (over-passage) occurrence to its
negative (under-passage) occurrence. The arrow label records the vertex label
and vertex sign separately from those passage signs. This is a presentation
view; it does not certify closed-normality or classical realizability.
"""
_gauss_svg_attribute(value) = replace(_esc(value), "\"" => "&quot;", "'" => "&apos;")

"""
    gauss_diagram_svg(state::OGraphState) -> String

Return an SVG Gauss diagram of the supplied one-circuit presentation.

# Arguments
`state` is an `OGraphState`; closed-normal certification is not required.
The empty MST circuit state is accepted. The drawing uses the representation
in [Gauss-code displays](../specs/ograph-terminology.md#2-standard-gauss-code-displays).

# Returns
A standalone SVG `String` containing the diagram, labels, styles and metadata
with the input convention identifier, word and vertex signs. No file is written.

# Guarantees
1. The diagram has one vertex arrow per stored vertex and one endpoint per
   stored occurrence, with the signs interpreted by the linked specification.
   [proved/agent: by construction]
2. The metadata preserves the stored word, vertex signs and convention label;
   attribute-sensitive characters in the label are escaped. [proved/agent: by construction]
3. The input is unchanged. This draws its chosen presentation, rather than
   selecting an equivalence-class representative or establishing a space.
   [proved/agent: by construction]

# Not guaranteed
The output does not certify closed-normality, classical realizability,
framing or embedding of a spine. Overlapping labels are not ruled out, and
SVG byte stability across library versions is not promised.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §2 — the geometric Gauss-diagram view and its separation from the text display.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); svg = gauss_diagram_svg(ograph_state(v));

julia> Tuple(parse.(Int, split(match(r"gauss=.([0-9,-]*).", svg).captures[1], ',')))
(1, -1)
```

See also [`plot_gauss_diagram_svg`](@ref).

Role: draw
"""
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

"""
    gauss_diagram_svg(state::ClosedNormalOGraph) -> String

Return an SVG Gauss diagram of a certified closed-normal presentation.

# Arguments
`state::ClosedNormalOGraph` supplies the presentation to draw; the rendering
does not repeat its certificate checks.

# Returns
A `String` SVG with the same fields and drawing scope as
`gauss_diagram_svg(state.state)`.

# Guarantees
1. The result equals `gauss_diagram_svg(state.state)` and leaves the input unchanged.
   [proved/agent: by construction]
2. The drawing is the [Gauss-diagram presentation view](../specs/ograph-terminology.md#2-standard-gauss-code-displays)
   of that closed-normal graph, not a geometric realization of its manifold.
   [proved/agent: by construction]

# Not guaranteed
Rendering does not strengthen the certificate or prove classical realizability.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §2 — the presentation view used for certified states too.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> svg = gauss_diagram_svg(certify_closed_ograph(vl_example("s3_abalone")));

julia> Tuple(parse.(Int, split(match(r"vertex-signs=.([0-9,-]*).", svg).captures[1], ',')))
(1,)
```

See also [`plot_gauss_diagram_svg`](@ref).

Role: draw
"""
gauss_diagram_svg(state::ClosedNormalOGraph) = gauss_diagram_svg(state.state)
"""
    gauss_diagram_svg(link::VirtualLink; convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) -> String

Return an SVG Gauss diagram of an unweighted one-circuit virtual-link presentation.

# Arguments
`link` must be accepted by [`ograph_state`](@ref), including its unweighted,
single-component and no-cut requirements. Its `ArgumentError` and conversion
failures propagate. `convention_id` is passed to that conversion unchanged.

# Returns
A `String` SVG of the converted state, with the content described by
the `OGraphState` overload.

# Guarantees
1. The result equals `gauss_diagram_svg(ograph_state(link; convention_id))`.
   [proved/agent: by construction]
2. The input is unchanged. The output is the
   [Gauss-diagram view](../specs/ograph-terminology.md#2-standard-gauss-code-displays)
   of its real presentation, with no additional claim about a represented space.
   [proved/agent: by construction]

# Not guaranteed
The call does not certify C2–C3, retain virtual crossings in the real view,
or display a weighted or cut presentation.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §2 — the geometric presentation display.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> svg = gauss_diagram_svg(vl_example("s3_abalone"));

julia> Tuple(parse.(Int, split(match(r"gauss=.([0-9,-]*).", svg).captures[1], ',')))
(1, -1)
```

See also [`ograph_state`](@ref), [`plot_gauss_diagram_svg`](@ref).

Role: draw
"""
gauss_diagram_svg(link::VirtualLink; convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) =
    gauss_diagram_svg(ograph_state(link; convention_id=convention_id))

"""
    plot_gauss_diagram_svg(state, path::AbstractString) -> String

Write the geometric Gauss-diagram view of a presentation as an SVG file.

# Arguments
- `state`: an input accepted by [`gauss_diagram_svg`](@ref); that method's
  domain checks and errors propagate.
- `path`: a writable file path whose parent exists. Filesystem errors
  propagate, and an existing file is overwritten.

# Returns
The absolute output path as a `String`.

# Guarantees
1. On successful return the file contains `gauss_diagram_svg(state)`. [proved/agent: by construction]
2. The image has the presentation-level meaning of the linked Gauss-diagram view. [proved/agent: by construction]
3. The input presentation is unchanged; the output file is written. [proved/agent: by construction]

# Not guaranteed
The view is not a planar embedding or a certificate of the represented
manifold. File replacement is not atomic and parent directories are not created.

# Sources
cf:docs/specs/ograph-terminology.md — geometric Gauss-diagram display.

# Examples
```jldoctest
julia> state = OGraphState([1,-1], [1]);

julia> mktempdir() do directory
           path = plot_gauss_diagram_svg(state, joinpath(directory, "diagram.svg"))
           Tuple(parse.(Int, split(match(r"gauss=.([0-9,-]*).", read(path, String)).captures[1], ',')))
       end
(1, -1)
```

Role: draw
"""
function plot_gauss_diagram_svg(state, path::AbstractString)
    svg = gauss_diagram_svg(state)
    open(path, "w") do io
        write(io, svg)
    end
    abspath(path)
end

Base.show(io::IO, ::MIME"image/svg+xml", state::OGraphState) = print(io, gauss_diagram_svg(state))
Base.show(io::IO, ::MIME"image/svg+xml", state::ClosedNormalOGraph) = print(io, gauss_diagram_svg(state))
