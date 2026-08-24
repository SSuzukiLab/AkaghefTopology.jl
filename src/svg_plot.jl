_esc(s) = replace(replace(replace(string(s), "&"=>"&amp;"), "<"=>"&lt;"), ">"=>"&gt;")

function _urd_weight_label(weight)
    raw=string(weight)
    length(raw)<=48 && return raw
    weight isa LaurentPolynomial && return "Laurent polynomial ($(length(weight.terms)) terms)"
    weight isa LaurentRational && return "Laurent rational ($(length(weight.numerator.terms))/$(length(weight.denominator.terms)) terms)"
    "symbolic expression ($(length(raw)) chars)"
end

function _rounded_svg_path(points; radius=0.3)
    length(points) < 3 && return "M "*join(("$x $y" for (x,y) in points)," L ")
    commands=String["M $(points[1][1]) $(points[1][2])"]
    for index in 2:length(points)-1
        previous=complex(points[index-1]...); corner=complex(points[index]...); following=complex(points[index+1]...)
        incoming=corner-previous; outgoing=following-corner
        if abs(incoming) == 0 || abs(outgoing) == 0
            push!(commands,"L $(points[index][1]) $(points[index][2])")
            continue
        end
        local_radius=min(radius,0.45abs(incoming),0.45abs(outgoing))
        before=corner-local_radius*incoming/abs(incoming)
        after=corner+local_radius*outgoing/abs(outgoing)
        push!(commands,"L $(real(before)) $(imag(before))")
        push!(commands,"Q $(real(corner)) $(imag(corner)) $(real(after)) $(imag(after))")
    end
    push!(commands,"L $(points[end][1]) $(points[end][2])")
    join(commands," ")
end

"""Drawing constants for the link renderer, in layout lattice units.

Working in lattice units rather than projected pixels is what makes these
consistent across figures: a crossing gap of `gap` is the same fraction of an
edge whether the diagram spans three units or thirty. The previous renderer
applied a fixed `0.20` before projection and a corner radius after it, so gaps
narrowed as diagrams grew.
"""
const LINK_STYLE = (strand=0.10, arrow=0.085, gap=0.18, radius=0.28,
                    vertex=0.115, font=0.32, margin=0.7, px_per_unit=72)

_num(x) = string(round(x;digits=4))
_pt(z::ComplexF64) = (real(z), -imag(z))   # SVG y grows downward

"""Render a weight without a spurious decimal tail.

o-graph weights are integers in every workflow here, but they are carried as
`Float64`, so the naive `string` gave `0.0` and `-1.0` on every label.
"""
function _fmt_weight(w)
    w isa Real || return string(w)
    isfinite(w) || return string(w)
    isinteger(w) ? string(Int(w)) : string(w)
end

"""Push overlapping labels apart, in place, in layout coordinates.

Offsetting a weight label along its arc's normal separates it from *its own*
strand but not from a parallel arc one lattice line over, which is exactly what
a bundle of arcs sharing a column produces. A few relaxation passes over the
label discs fix that without needing a real label-placement solver.
"""
function _spread_labels!(placed; passes=40)
    radius(text)=LINK_STYLE.font*(0.30+0.26*length(text))
    for _ in 1:passes
        moved=false
        for i in 1:length(placed), j in i+1:length(placed)
            want=radius(placed[i].text)+radius(placed[j].text)
            delta=placed[j].pos-placed[i].pos
            gap=abs(delta)
            gap >= want && continue
            push_by = gap < 1e-6 ? complex(0.0,want/2) : (want-gap)/2*delta/gap
            placed[i]=(pos=placed[i].pos-push_by,text=placed[i].text)
            placed[j]=(pos=placed[j].pos+push_by,text=placed[j].text)
            moved=true
        end
        moved || break
    end
    placed
end

"""A text element that stays readable where it crosses a strand.

`paint-order="stroke"` draws the halo first and the glyph over it, so labels do
not have to be nudged away from every line to remain legible.
"""
function _label(x,y,text; anchor="middle", size=LINK_STYLE.font)
    "<text class=\"label\" x=\"$(_num(x))\" y=\"$(_num(y))\" text-anchor=\"$anchor\" " *
    "dominant-baseline=\"middle\" font-family=\"sans-serif\" font-size=\"$(_num(size))\" " *
    "fill=\"#12161d\" stroke=\"#faf9f7\" stroke-width=\"$(_num(size*0.28))\" " *
    "paint-order=\"stroke\" stroke-linejoin=\"round\">$(_esc(text))</text>"
end

"""Point and unit tangent at the arc-length midpoint of a polyline.

The previous renderer placed the arrow at the midpoint *index*, so it jumped
whenever the number of corners changed.
"""
function _arclength_midpoint(points)
    length(points) < 2 && return (first(points), complex(1.0,0.0))
    spans=[abs(points[i+1]-points[i]) for i in 1:length(points)-1]
    target=sum(spans)/2
    walked=0.0
    for (i,span) in enumerate(spans)
        span == 0 && continue
        if walked+span >= target || i == length(spans)
            t=(target-walked)/span
            direction=(points[i+1]-points[i])/span
            return (points[i]+t*(points[i+1]-points[i]), direction)
        end
        walked += span
    end
    (points[end], complex(1.0,0.0))
end

"""Wrap link body markup in a themed, content-fitted SVG.

The viewBox is the content bounding box plus a margin, so the drawing scales
itself and every length above stays in lattice units. Colours are emitted twice:
as presentation attributes, which any renderer honours, and as CSS classes that
override them where CSS is supported. That keeps `rsvg-convert` output correct
while still letting a browser theme the figure.
"""
function _write_link_svg(path, body, bbox; title="o-graph", desc="", meta=nothing)
    xmin,xmax,ymin,ymax = bbox
    m = LINK_STYLE.margin
    width  = (xmax-xmin)+2m
    height = (ymax-ymin)+2m
    view = "$(_num(xmin-m)) $(_num(-(ymax+m))) $(_num(width)) $(_num(height))"
    px = LINK_STYLE.px_per_unit
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" class=\"ograph\" ",
                   "width=\"$(_num(width*px))\" height=\"$(_num(height*px))\" viewBox=\"$view\" ",
                   "role=\"img\" aria-label=\"$(_esc(title))\">")
        println(io,"<title>$(_esc(title))</title>")
        isempty(desc) || println(io,"<desc>$(_esc(desc))</desc>")
        meta === nothing || println(io,"<metadata><ograph>$(_esc(meta))</ograph></metadata>")
        println(io,"""<style>
.ograph{--ground:#faf9f7;--strand:#1b3bd0;--ink:#12161d}
@media (prefers-color-scheme:dark){.ograph{--ground:#0e1116;--strand:#7d97ff;--ink:#e6e9ef}}
.ograph .bg{fill:var(--ground)}
.ograph .arc,.ograph .arrow{stroke:var(--strand);fill:none}
.ograph .vertex{fill:var(--strand)}
.ograph .label{fill:var(--ink);stroke:var(--ground)}
</style>""")
        println(io,"<rect class=\"bg\" x=\"$(_num(xmin-m))\" y=\"$(_num(-(ymax+m)))\" ",
                   "width=\"$(_num(width))\" height=\"$(_num(height))\" fill=\"#faf9f7\"/>")
        println(io, body)
        println(io,"</svg>")
    end
    abspath(path)
end

function _write_svg(path, body; width=900, height=600)
    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#faf9f7\"/>")
        println(io, body)
        println(io, "</svg>")
    end
    abspath(path)
end

function plot_svg(d::URDiagram, path::AbstractString)
    # The MATLAB URDiagram.plot implementation is intentionally minimal: one
    # left-pointing horizontal quiver of length length(V)+1.
    body=String[]
    push!(body,"<line x1=\"780\" y1=\"300\" x2=\"120\" y2=\"300\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
    push!(body,"<path d=\"M 150 280 L 120 300 L 150 320\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
    _write_svg(path,join(body,"\n"))
end

"""Draw a readable UR-diagram diagnostic, including signed endpoints and weights.

`plot_svg(::URDiagram, ...)` deliberately preserves MATLAB's one-arrow
`URDiagram.plot`.  That rendering cannot validate a reduction: this companion
view makes the pairing represented by `V` and the current edge weights visible
without changing the source-compatible output.
"""
function plot_urd_structure_svg(d::URDiagram, path::AbstractString)
    validate_data(d)
    width,height=960,480
    margin=80
    count=max(length(d.V),2)
    x(position)=margin+(width-2margin)*(position-1)/(count-1)
    baseline=260
    body=String[
        "<text x=\"$(margin)\" y=\"42\" font-family=\"ui-monospace,monospace\" font-size=\"18\" fill=\"#252525\">UR diagram: signed endpoints, paired edges, weights</text>",
        "<line x1=\"$margin\" y1=\"$baseline\" x2=\"$(width-margin)\" y2=\"$baseline\" stroke=\"#a8a29e\" stroke-width=\"2\"/>"
    ]
    for edge in vec(d.E[1,:])
        positive=vertex_index(d,edge); negative=vertex_index(d,-edge)
        left,right=minmax(x(positive),x(negative))
        lane=42+28*((edge-1)%4)
        control=(left+right)/2
        color=positive<negative ? "#2563eb" : "#9333ea"
        push!(body,"<path d=\"M $left $baseline Q $control $(baseline-lane) $right $baseline\" fill=\"none\" stroke=\"$color\" stroke-width=\"3\"/>")
        label=_urd_weight_label(d.W[edge_index(d,edge)])
        push!(body,"<text x=\"$control\" y=\"$(baseline-lane-10)\" text-anchor=\"middle\" font-family=\"ui-monospace,monospace\" font-size=\"14\" fill=\"#292524\">e$(edge): $(_esc(label))</text>")
    end
    for (position,vertex) in enumerate(d.V)
        color=vertex>0 ? "#1d4ed8" : "#7f1d1d"
        push!(body,"<circle cx=\"$(x(position))\" cy=\"$baseline\" r=\"15\" fill=\"#fff\" stroke=\"$color\" stroke-width=\"3\"/>")
        push!(body,"<text x=\"$(x(position))\" y=\"$(baseline+5)\" text-anchor=\"middle\" font-family=\"ui-monospace,monospace\" font-size=\"14\" fill=\"$color\">$vertex</text>")
    end
    _write_svg(path,join(body,"\n");width=width,height=height)
end

_crossing_sign(slot::Integer) = slot in (1,3) ? -1 : 1
_is_outgoing(slot::Integer,orientation::Integer) =
    slot == 3 || (slot == 4 && orientation == -1) || (slot == 2 && orientation != -1)

function _directed_piece(edge::EdgeLayout,orientation)
    if _is_outgoing(edge.tail_slot,orientation[edge.tail])
        source=_crossing_sign(edge.tail_slot)*edge.tail
        target=_crossing_sign(edge.head_slot)*edge.head
        return (source=source,target=target,points=edge.points,id=edge.id)
    end
    source=_crossing_sign(edge.head_slot)*edge.head
    target=_crossing_sign(edge.tail_slot)*edge.tail
    (source=source,target=target,points=reverse(edge.points),id=edge.id)
end

"""Aggregate planar PD edges into MATLAB REdgeTable arcs between real crossings.

This stage has no counterpart inside Sage: MATLAB does the equivalent work
itself after Sage returns, in `VirtualLink.calcPos` (`connectPolylines`, the
per-arc reversal check, component offsets). The two have never been compared.
See `PLOT_PIPELINE.md` section 8.
"""
function _real_arc_layouts(v::VirtualLink,layout::LinkLayout)
    # `orthogonal_layout` planarised the same input already; recomputing keeps
    # the two call sites able to drift apart. See `PLOT_PIPELINE.md` IS6.
    code = v.real_only ? first(planarize_virtual_gauss(v.gauss,v.orientation)) : v.gauss
    pieces=Dict{Tuple{Int,Int},Vector{NamedTuple}}()
    for edge in layout.edges
        piece=_directed_piece(edge,layout.orientation)
        push!(get!(pieces,(piece.source,piece.target),NamedTuple[]),piece)
    end
    arcs=NamedTuple[]
    arc_id=0
    for (component_index,component) in enumerate(code)
        component=filter(!=(0),component)
        isempty(component) && continue
        sequence=NamedTuple[]
        for index in eachindex(component)
            pair=(component[index],component[mod1(index+1,length(component))])
            candidates=get(pieces,pair,NamedTuple[])
            isempty(candidates) && throw(ArgumentError("no directed PD edge for Gauss pair $pair"))
            push!(sequence,popfirst!(candidates))
        end
        real_positions=findall(index -> layout.orientation[abs(component[index])] != 0,eachindex(component))
        for (weight_index,target_position) in enumerate(real_positions)
            previous_position=real_positions[mod1(weight_index-1,length(real_positions))]
            indices=Int[]
            cursor=previous_position
            while cursor != target_position
                push!(indices,cursor)
                cursor=mod1(cursor+1,length(component))
            end
            points=ComplexF64[]
            for index in indices
                append!(points,isempty(points) ? sequence[index].points : sequence[index].points[2:end])
            end
            arc_id += 1
            component_weights=component_index <= length(v.weights) ? v.weights[component_index] : Float64[]
            weight=weight_index <= length(component_weights) ? component_weights[weight_index] : NaN
            push!(arcs,(id=arc_id,points=points,source=component[previous_position],
                        target=component[target_position],weight=weight))
        end
    end
    arcs
end

"""Stable MATLAB `REdgeTable` arc metadata, without invoking the layout solver."""
function _replay_arc_metadata(v::VirtualLink)
    [(id=edge.id,source=edge.crossing[1],target=edge.crossing[2],weight=edge.weight)
     for edge in real_edge_table(v)]
end

"""Render a VirtualLink using either the Julia layout or replayed MATLAB arcs.

`replay_positions`, when supplied, is one polyline per real arc in the stable
order returned by `_real_arc_layouts`.  It is the direct representation of
MATLAB `REdgeTable.Position`.  It bypasses the layout solver entirely.

As of 2026-08-20 this is an optional comparison tool, not an acceptance path:
plots are accepted by property checking rather than by matching MATLAB, and no
canonical Windows fixture can be produced in any case because the MATLAB plot
path cannot run there (`PLOT_PIPELINE.md` section 7).

The full stage sequence, and the drawing conventions applied below that have
not been compared with MATLAB's renderer -- the fixed `0.20` under-strand gap
in pre-projection units, the `0.3 * scale` corner radius, and the arrow placed
at the midpoint *index* rather than the arc-length midpoint -- are recorded in
`PLOT_PIPELINE.md` section 3.
"""
function plot_svg(v::VirtualLink, path::AbstractString; bending_numbers=nothing,
                  edge_labels=:weight,replay_positions=nothing)
    body=String[]
    if isempty(v.orientation)
        S=LINK_STYLE
        push!(body,"<circle class=\"arc\" cx=\"0\" cy=\"0\" r=\"1\" fill=\"none\" ",
                   "stroke=\"#1b3bd0\" stroke-width=\"$(_num(S.strand))\"/>")
        push!(body,"<path class=\"arrow\" d=\"M 0.88 -0.20 L 1.0 0 L 1.12 -0.20\" fill=\"none\" ",
                   "stroke=\"#1b3bd0\" stroke-width=\"$(_num(S.arrow))\" stroke-linejoin=\"round\"/>")
        return _write_link_svg(path,join(body,"\n"),(-1.0,1.0,-1.0,1.0);
                               title="unknot",desc="Unknot: one closed component, no crossings.",
                               meta="crossings=0; arcs=0; weighted=false")
    end
    layout=nothing
    arcs=if replay_positions === nothing
        layout=orthogonal_layout(v;bending_numbers=bending_numbers)
        _real_arc_layouts(v,layout)
    else
        metadata=_replay_arc_metadata(v)
        length(replay_positions)==length(metadata) || throw(ArgumentError("one MATLAB Position polyline is required for each real arc"))
        [merge(arc,(points=ComplexF64.(collect(points)),)) for (arc,points) in zip(metadata,replay_positions)]
    end
    if replay_positions !== nothing
        all(length(arc.points)>=2 for arc in arcs) || throw(ArgumentError("every MATLAB Position polyline requires at least two points"))
    end
    all_points=reduce(vcat,(arc.points for arc in arcs))
    xmin,xmax=extrema(real.(all_points)); ymin,ymax=extrema(imag.(all_points))
    S=LINK_STYLE
    placed=NamedTuple[]

    for arc in arcs
        points=copy(arc.points)
        # Under-strand ends are shortened to open the crossing gap. In lattice
        # units, so the gap reads the same on every diagram.
        if arc.source < 0 && length(points)>1
            points[1] += S.gap*(points[2]-points[1])/abs(points[2]-points[1])
        end
        if arc.target < 0 && length(points)>1
            points[end] += S.gap*(points[end-1]-points[end])/abs(points[end-1]-points[end])
        end
        data=_rounded_svg_path([_pt(point) for point in points];radius=S.radius)
        push!(body,"<path class=\"arc\" data-arc=\"$(arc.id)\" d=\"$data\" fill=\"none\" ",
                   "stroke=\"#1b3bd0\" stroke-width=\"$(_num(S.strand))\" ",
                   "stroke-linejoin=\"round\" stroke-linecap=\"round\"/>")

        midpoint,unit=_arclength_midpoint(points)
        back=midpoint-unit*0.20
        left=back+unit*im*0.11; right=back-unit*im*0.11
        arrow=[_pt(left),_pt(midpoint),_pt(right)]
        arrow_data="M $(_num(arrow[1][1])) $(_num(arrow[1][2])) L $(_num(arrow[2][1])) $(_num(arrow[2][2])) L $(_num(arrow[3][1])) $(_num(arrow[3][2]))"
        push!(body,"<path class=\"arrow\" d=\"$arrow_data\" fill=\"none\" stroke=\"#1b3bd0\" ",
                   "stroke-width=\"$(_num(S.arrow))\" stroke-linejoin=\"round\" stroke-linecap=\"round\"/>")

        if v.is_weighted || edge_labels in (:id,:all)
            parts=String[]
            edge_labels in (:id,:all) && push!(parts,"E$(arc.id)")
            v.is_weighted && edge_labels in (:weight,:all) && !isnan(arc.weight) &&
                push!(parts,_fmt_weight(arc.weight))
            label=join(parts,", ")
            # Sit the label beside the arc rather than on it: offset along the
            # normal of the local tangent, consistently to the left of travel.
            isempty(label) || push!(placed,(pos=midpoint+unit*im*0.26,text=label))
        end
    end

    vertices=ComplexF64[]
    if replay_positions === nothing
        for (index,point) in enumerate(layout.crossings)
            layout.orientation[index] == 0 && continue
            push!(vertices,point)
        end
    else
        crossings=Dict{Int,ComplexF64}()
        for arc in arcs
            for (crossing,point) in ((abs(arc.source),first(arc.points)),(abs(arc.target),last(arc.points)))
                haskey(crossings,crossing) && !isapprox(crossings[crossing],point;atol=1e-9) &&
                    throw(ArgumentError("MATLAB Position polylines disagree at crossing $crossing"))
                crossings[crossing]=point
            end
        end
        vertices=[crossings[k] for k in sort!(collect(keys(crossings)))]
    end
    # Vertex labels are pushed away from the diagram's centre, so in a dense
    # figure they fan outwards instead of piling into the middle.
    centre = isempty(vertices) ? complex(0.0,0.0) : sum(vertices)/length(vertices)
    for (real_index,point) in enumerate(vertices)
        x,y=_pt(point)
        push!(body,"<circle class=\"vertex\" cx=\"$(_num(x))\" cy=\"$(_num(y))\" ",
                   "r=\"$(_num(S.vertex))\" fill=\"#1b3bd0\"/>")
        away = point - centre
        away = abs(away) < 1e-9 ? complex(0.7071,0.7071) : away/abs(away)
        push!(placed,(pos=point + away*0.30,text="V$real_index"))
    end

    _spread_labels!(placed)
    for item in placed
        lx,ly=_pt(item.pos)
        push!(body,_label(lx,ly,item.text))
        xmin=min(xmin,real(item.pos)); xmax=max(xmax,real(item.pos))
        ymin=min(ymin,imag(item.pos)); ymax=max(ymax,imag(item.pos))
    end

    name=basename(path); name=name[1:something(findlast('.',name),length(name)+1)-1]
    gauss=join(("["*join(c,",")*"]" for c in v.gauss),", ")
    desc="$(length(vertices)) real crossing(s), $(length(arcs)) arc(s). Gauss $gauss."
    _write_link_svg(path,join(body,"\n"),(xmin,xmax,ymin,ymax);
                    title=name,desc=desc,
                    meta="crossings=$(length(vertices)); arcs=$(length(arcs)); weighted=$(v.is_weighted)")
end
