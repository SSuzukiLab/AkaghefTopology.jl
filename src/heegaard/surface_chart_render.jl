"""SVG/MIME rendering for the full-curve atlas. No browser geometry is used."""

_surface_beta_color(id) = ("#2459a9","#27806b","#8650a1","#986d15")[mod1(id,4)]
_surface_signed(x) = string(x.label,x.sign==1 ? "+" : "−")

function _surface_arrow(point,tangent,size,color;attrs="")
    u=_unit(tangent); n=(-u[2],u[1])
    tip=_vadd(point,_vscale(size,u))
    back=_vsub(point,_vscale(size,u))
    a=_vadd(back,_vscale(0.60size,n)); b=_vsub(back,_vscale(0.60size,n))
    "<polygon $attrs points=\"$(_svg_point(tip)) $(_svg_point(a)) $(_svg_point(b))\" fill=\"$color\"/>"
end

function _surface_text(x,y,value;size=14,fill="#34444e",anchor="middle",attrs="")
    "<text x=\"$(_fmt(x))\" y=\"$(_fmt(y))\" font-size=\"$(_fmt(size))\" " *
    "fill=\"$fill\" text-anchor=\"$anchor\" $attrs>$(_svg_escape(value))</text>"
end

function _surface_chart_json(chart::HeegaardSurfaceChart)
    Dict("schema"=>"akaghef.heegaard.surface-chart/1","diagram_hash"=>chart.diagram_hash,
        "genus"=>chart.genus,"model"=>"punctured sphere + g annuli",
        "base_orientation"=>"-dx wedge dy","annulus_orientation"=>"dx wedge dy",
        "annuli"=>chart.annuli,"intersections"=>chart.intersections,
        "alpha_cycles"=>chart.alpha_cycles,"beta_cycles"=>chart.beta_cycles,
        "projection_crossings"=>0,"provenance"=>chart.provenance,
        "joining_disk"=>_joining_overlay_json(chart.joining_disk),
        "base_layout"=>_layout_json(chart.base))
end

function _surface_base_svg(chart,x,y,width,height)
    layout=chart.base
    coordinates=NTuple{2,Float64}[p for r in layout.routes for p in r.points]
    for p in layout.punctures
        push!(coordinates,_vadd(p.center,(p.radius,p.radius)),_vsub(p.center,(p.radius,p.radius)))
    end
    isempty(coordinates) && push!(coordinates,(0.,0.),(1.,1.))
    xmin,xmax=extrema(first.(coordinates)); ymin,ymax=extrema(last.(coordinates))
    margin=0.20max(xmax-xmin,ymax-ymin,1.)
    vb=join(_fmt.((xmin-margin,-ymax-margin,xmax-xmin+2margin,ymax-ymin+2margin))," ")
    parts=String["<svg x=\"$x\" y=\"$y\" width=\"$width\" height=\"$height\" viewBox=\"$vb\" data-chart=\"sphere\">"]
    if chart.joining_disk!==nothing
        d=chart.joining_disk
        path=d.base_path
        if d.base_exterior
            # D contains infinity in this planar chart. Even-odd filling of the
            # viewport complement preserves it; puncture interiors are masked below.
            path="M $(_fmt(xmin-margin)) $(_fmt(-ymax-margin)) h $(_fmt(xmax-xmin+2margin)) v $(_fmt(ymax-ymin+2margin)) h $(_fmt(-(xmax-xmin+2margin))) Z "*path
        end
        push!(parts,"<path data-joining-disk=\"$(d.disk.region)\" data-disk-piece=\"base\" d=\"$path\" fill=\"#f3ca79\" fill-opacity=\"0.68\" fill-rule=\"evenodd\"/>")
        if d.base_label!==nothing
            q=d.base_label
            push!(parts,_surface_text(q[1],-q[2],"D";size=0.13,fill="#865500",attrs="data-joining-disk=\"$(d.disk.region)\""))
        end
    end
    for r in layout.routes
        color=_surface_beta_color(r.curve_id)
        attrs="data-curve-family=\"beta\" data-curve-id=\"$(r.curve_id)\" data-route-segment=\"$(r.segment_id)\""
        push!(parts,"<path $attrs d=\"$(_route_path(r))\" fill=\"none\" stroke=\"$color\" stroke-width=\"2\" stroke-linecap=\"round\" vector-effect=\"non-scaling-stroke\"/>")
        k=clamp(length(r.points)÷2,1,length(r.points)-1)
        push!(parts,_surface_arrow(_midpoint(r.points[k],r.points[k+1]),_vsub(r.points[k+1],r.points[k]),
                                  0.032,color;attrs))
    end
    for p in layout.punctures
        a=only(a for a in chart.annuli if a.pair==p.pair)
        label="B$(a.alpha_id)"*(p.side==:plus ? "+" : "−")
        push!(parts,"<circle cx=\"$(_fmt(p.center[1]))\" cy=\"$(_fmt(-p.center[2]))\" r=\"$(_fmt(p.radius))\" fill=\"#fff\" stroke=\"#80918e\" stroke-width=\"1.4\" stroke-dasharray=\"5 4\" vector-effect=\"non-scaling-stroke\" data-boundary-id=\"$(p.id)\"/>")
        push!(parts,_surface_text(p.center[1],-p.center[2]+0.030,label;size=0.36p.radius,fill="#637670"))
        for port in filter(q->q.boundary_id==p.id,layout.ports)
            q=port.point
            point=only(i for i in chart.intersections if i.id==port.intersection_id)
            attachment=point.label==string(point.id) ? "x$(point.id)" : point.label
            labelpoint=_vadd(p.center,_vscale(0.71,_vsub(q,p.center)))
            attrs="data-link-id=\"$(port.intersection_id)\" data-port-id=\"$(port.id)\" data-curve-family=\"beta\" data-curve-id=\"$(port.beta_id)\""
            push!(parts,"<circle $attrs cx=\"$(_fmt(q[1]))\" cy=\"$(_fmt(-q[2]))\" r=\"$(_fmt(0.045p.radius))\" fill=\"white\" stroke=\"$(_surface_beta_color(port.beta_id))\" stroke-width=\"1.4\" vector-effect=\"non-scaling-stroke\"><title>Attachment $(_svg_escape(attachment)) at $label; not an additional intersection</title></circle>")
            push!(parts,_surface_text(labelpoint[1],-labelpoint[2]+0.018,attachment;size=0.25p.radius,attrs))
        end
    end
    push!(parts,"</svg>")
    join(parts,"\n")
end

function _surface_annulus_svg(chart,a,x,y)
    parts=String["<g data-chart=\"annulus\" data-alpha-id=\"$(a.alpha_id)\" transform=\"translate($x $y)\">",
                 _surface_text(0,-130,"Annulus $(a.alpha_id) · α$(a.alpha_id)";size=18),
                 "<g transform=\"scale(106)\">"]
    push!(parts,"<circle r=\"1\" fill=\"#eef4ee\" stroke=\"#80918e\" stroke-width=\"1.3\" stroke-dasharray=\"5 4\" vector-effect=\"non-scaling-stroke\"/>")
    push!(parts,"<circle r=\"0.43\" fill=\"white\" stroke=\"#80918e\" stroke-width=\"1.3\" stroke-dasharray=\"5 4\" vector-effect=\"non-scaling-stroke\"/>")
    if chart.joining_disk!==nothing
        d=chart.joining_disk
        for sector in filter(s->s.alpha_id==a.alpha_id,d.sectors)
            push!(parts,"<path data-joining-disk=\"$(d.disk.region)\" data-disk-piece=\"annulus\" data-boundary-dart=\"$(sector.dart)\" d=\"$(sector.path)\" fill=\"#f3ca79\" fill-opacity=\"0.82\"/>")
            q=sector.label
            push!(parts,_surface_text(q[1],-q[2]+0.035,"D";size=0.13,fill="#865500",attrs="data-joining-disk=\"$(d.disk.region)\""))
        end
    end
    attrs="data-curve-family=\"alpha\" data-curve-id=\"$(a.alpha_id)\""
    push!(parts,"<circle $attrs r=\"0.70\" fill=\"none\" stroke=\"#be4c45\" stroke-width=\"2.8\" vector-effect=\"non-scaling-stroke\"/>")
    push!(parts,_surface_arrow((0.70,0.),(0.,-1.),0.036,"#be4c45";attrs))
    push!(parts,_surface_text(0,-0.08,"α$(a.alpha_id)";size=0.21,fill="#be4c45",attrs))
    push!(parts,_surface_text(0,0.16,"B$(a.alpha_id)−";size=0.14,fill="#637670"))
    for p in filter(p->p.alpha_id==a.alpha_id,chart.intersections)
        attrs="data-curve-family=\"beta\" data-curve-id=\"$(p.beta_id)\" data-link-id=\"$(p.id)\""
        color=_surface_beta_color(p.beta_id)
        s,f=p.beta_start,p.beta_end
        # A radial interval is an exact geodesic in this annulus chart.
        push!(parts,"<path $attrs d=\"M $(_svg_point(s)) L $(_svg_point(f))\" fill=\"none\" stroke=\"$color\" stroke-width=\"2\" vector-effect=\"non-scaling-stroke\"/>")
        arrowpoint=_surface_radial(p.sign==1 ? 0.87 : 0.53,p.parameter)
        push!(parts,_surface_arrow(arrowpoint,p.beta_tangent,0.030,color;attrs))
        for q in (s,f)
            push!(parts,"<circle $attrs cx=\"$(_fmt(q[1]))\" cy=\"$(_fmt(-q[2]))\" r=\"0.021\" fill=\"white\" stroke=\"$color\" stroke-width=\"1.2\" vector-effect=\"non-scaling-stroke\"/>")
        end
        q=p.point
        push!(parts,"<circle data-intersection-id=\"$(p.id)\" data-sign=\"$(p.sign)\" data-alpha-id=\"$(p.alpha_id)\" data-beta-id=\"$(p.beta_id)\" data-link-id=\"$(p.id)\" cx=\"$(_fmt(q[1]))\" cy=\"$(_fmt(-q[2]))\" r=\"0.030\" fill=\"#24333c\" stroke=\"white\" stroke-width=\"1\" vector-effect=\"non-scaling-stroke\"><title>Intersection $(_surface_signed(p)) · α$(p.alpha_id) visit $(p.alpha_position) · β$(p.beta_id) visit $(p.beta_position)</title></circle>")
        labelpoint=_surface_radial(1.16,p.parameter)
        push!(parts,_surface_text(labelpoint[1],-labelpoint[2]+0.035,_surface_signed(p);size=0.12,attrs="data-link-id=\"$(p.id)\""))
    end
    push!(parts,"</g>",_surface_text(0,139,"B$(a.alpha_id)+ outer ↔ sphere +";size=12,fill="#637670"),
                 _surface_text(0,158,"B$(a.alpha_id)− inner ↔ sphere −";size=12,fill="#637670"),"</g>")
    join(parts,"\n")
end

function _render_surface_chart_svg(chart::HeegaardSurfaceChart)
    height=max(510,70+340chart.genus)
    meta=Dict("model"=>"punctured sphere + annuli","genus"=>chart.genus,
        "base_orientation"=>"-dx wedge dy","annulus_orientation"=>"dx wedge dy",
        "annuli"=>chart.annuli,"alpha_cycles"=>chart.alpha_cycles,
        "beta_cycles"=>chart.beta_cycles,"provenance"=>chart.provenance,
        "joining_disk"=>_joining_overlay_json(chart.joining_disk))
    parts=String["<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 1180 $height\" role=\"img\" data-heegaard-view=\"surface\" data-genus=\"$(chart.genus)\" data-diagram-hash=\"$(chart.diagram_hash)\" font-family=\"system-ui,sans-serif\">",
        "<title>Heegaard surface atlas · genus $(chart.genus)</title>",
        "<desc>A sphere with $(2chart.genus) holes, joined by $(chart.genus) annuli. Red alpha cores cross blue/green beta curves at the signed black points. Open points mark chart attachments. No projection crossings.</desc>",
        "<metadata>$(_svg_escape(JSON3.write(meta)))</metadata>",
        "<rect x=\"8\" y=\"8\" width=\"782\" height=\"$(height-16)\" rx=\"16\" fill=\"#f6f7f3\"/>",
        _surface_text(34,41,"S² ∖ $(2chart.genus) disks";size=21,anchor="start"),
        _surface_text(34,66,"The outside continues through ∞. Dashed circles attach to the annuli.";size=13,fill="#6b787d",anchor="start"),
        _surface_base_svg(chart,15,78,768,height-145),
        _surface_text(34,height-28,"○ chart attachment xᵢ        ● signed intersection i± (in annulus)";size=13,anchor="start")]
    for (index,a) in enumerate(chart.annuli)
        push!(parts,_surface_annulus_svg(chart,a,980,182+340(index-1)))
    end
    push!(parts,"</svg>")
    join(parts,"\n")
end

function render(chart::HeegaardSurfaceChart;backend=:svg,path=nothing)
    backend in (:svg,:pure_julia,:cairo) || throw(ArgumentError("unsupported renderer backend $backend"))
    svg=_render_surface_chart_svg(chart)
    path===nothing || open(io->write(io,svg),path,"w")
    svg
end

Base.show(io::IO,::MIME"image/svg+xml",chart::HeegaardSurfaceChart)=print(io,render(chart))
inspect_html(chart::HeegaardSurfaceChart)="<!doctype html><meta charset=\"utf-8\"><title>Heegaard surface atlas</title>"*render(chart)
