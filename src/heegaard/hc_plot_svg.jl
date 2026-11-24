# cf: lib/topology.jl/docs/specs/heegaard-plot.md — polyline output and explicit role visibility
# cf: lib/topology.jl/docs/specs/lickorish-code.md §2 — display frame: placement coordinates are screen
#     coordinates (y downwards); l_i^+ is the lower row and l_i turns clockwise around it.
_hp_xml(s)=replace(string(s),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;",'\''=>"&apos;")
function _hp_svg(p::HCPlotPlacement;show=:both,title="Heegaard diagram · placement")
    show in (:both,:alpha,:beta) || throw(ArgumentError("show must be :both, :alpha or :beta"))
    pts=[q for r in p.runs for q in r.points]
    xmin=min(0.0,minimum((q[1] for q in pts);init=0.0))-60
    xmax=max(p.width,maximum((q[1] for q in pts);init=p.width))+60
    ymin=min(0.0,minimum((q[2] for q in pts);init=0.0))-60
    ymax=max(p.height,maximum((q[2] for q in pts);init=p.height))+60
    io=IOBuffer();println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$xmin $ymin $(xmax-xmin) $(ymax-ymin)\">")
    println(io,"<rect x=\"$xmin\" y=\"$ymin\" width=\"$(xmax-xmin)\" height=\"$(ymax-ymin)\" fill=\"#fcfbf7\"/>")
    println(io,"<text x=\"$(xmin+16)\" y=\"$(ymin+28)\" font-size=\"20\" font-family=\"sans-serif\">$(_hp_xml(title))</text>")
    if p.square
        println(io,"<path data-chart=\"square\" d=\"M0,0 L$(p.width),0 L$(p.width),$(p.height) L0,$(p.height) Z\" fill=\"none\" stroke=\"#777\" stroke-width=\"1.2\"/>")
        for (x,y,label,anchor) in ((p.width/2,p.height+24.,"ℓ1⁺","middle"),(p.width/2,-12.,"ℓ1⁻","middle"),(-12.,p.height/2,"m1⁻","end"),(p.width+12.,p.height/2,"m1⁺","start"))
            println(io,"<text x=\"$x\" y=\"$y\" text-anchor=\"$anchor\" font-size=\"13\">$label</text>")
        end
    end
    for r in p.runs
        show!=:both && r.role!=show && continue
        color=r.role==:alpha ? "#ba433f" : r.role==:beta ? "#255dba" : "#555555"
        data="M"*join(("$(q[1]),$(q[2])" for q in r.points)," L")*(r.closed ? " Z" : "")
        println(io,"<path d=\"$data\" fill=\"none\" stroke=\"$color\" stroke-width=\"2.5\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>")
    end
    for (name,(c,r)) in p.holes
        println(io,"<circle cx=\"$(c[1])\" cy=\"$(c[2])\" r=\"$r\" fill=\"#fcfbf7\" stroke=\"#333\" stroke-width=\"1.5\"/>")
        label=replace(name,"top"=>"l", "bot"=>"l")*(startswith(name,"top") ? "⁺" : "⁻")
        println(io,"<text x=\"$(c[1])\" y=\"$(c[2]+4)\" text-anchor=\"middle\" font-size=\"13\">$(_hp_xml(label))</text>")
    end
    println(io,"</svg>");String(take!(io))
end
