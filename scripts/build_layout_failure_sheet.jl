using Topology, JSON

const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"examples","OgraphPlotSuite.jl"))

xml_escape(value)=replace(string(value),"&"=>"&amp;","<"=>"&lt;",">"=>"&gt;","\""=>"&quot;")
num(x)=string(round(x;digits=3))
svg_point(z)="$(num(real(z))),$(num(-imag(z)))"

function failure_svg(name,link,path)
    layout=orthogonal_layout(link)
    diagnostics=check_layout(layout)
    failing=[d for d in diagnostics if d.severity===:error]
    failing_edges=Set(edge for d in failing for edge in d.edges)
    marker_keys=Set{Tuple{Symbol,Float64,Float64}}()
    markers=NamedTuple[]
    for d in failing
        d.at===nothing && continue
        key=(d.rule,round(real(d.at);digits=6),round(imag(d.at);digits=6))
        key in marker_keys && continue
        push!(marker_keys,key)
        push!(markers,(rule=d.rule,point=d.at,message=d.message,edges=d.edges))
    end
    all_points=reduce(vcat,[edge.points for edge in layout.edges])
    xmin=min(real.(all_points)...); xmax=max(real.(all_points)...)
    ymin=min(imag.(all_points)...); ymax=max(imag.(all_points)...)
    margin=1.0
    x0=xmin-margin; y0=-(ymax+margin)
    width=(xmax-xmin)+2margin; height=(ymax-ymin)+2margin
    body=String[]
    push!(body,"<rect width=\"100%\" height=\"100%\" fill=\"#fffdf8\"/>")
    for x in floor(xmin-margin):ceil(xmax+margin)
        push!(body,"<line x1=\"$(num(x))\" y1=\"$(num(y0))\" x2=\"$(num(x))\" y2=\"$(num(-(ymin-margin)))\" stroke=\"#edf0f4\" stroke-width=\"0.015\"/>")
    end
    for y in floor(ymin-margin):ceil(ymax+margin)
        push!(body,"<line x1=\"$(num(x0))\" y1=\"$(num(-y))\" x2=\"$(num(xmax+margin))\" y2=\"$(num(-y))\" stroke=\"#edf0f4\" stroke-width=\"0.015\"/>")
    end
    for edge in layout.edges
        color=edge.id in failing_edges ? "#d92d20" : "#8ca0b5"
        stroke_width=edge.id in failing_edges ? "0.105" : "0.055"
        points_text=join(svg_point.(edge.points)," ")
        push!(body,"<polyline points=\"$points_text\" fill=\"none\" stroke=\"$color\" stroke-width=\"$stroke_width\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>")
        midpoint=edge.points[cld(length(edge.points),2)]
        push!(body,"<text x=\"$(num(real(midpoint)))\" y=\"$(num(-imag(midpoint)+0.16))\" font-size=\"0.22\" fill=\"$color\" text-anchor=\"middle\">e$(edge.id)</text>")
    end
    for (index,point) in enumerate(layout.crossings)
        real_crossing=layout.orientation[index]!=0
        color=real_crossing ? "#155eef" : "#667085"
        push!(body,"<circle cx=\"$(num(real(point)))\" cy=\"$(num(-imag(point)))\" r=\"0.13\" fill=\"white\" stroke=\"$color\" stroke-width=\"0.055\"/>")
        push!(body,"<text x=\"$(num(real(point)+0.16))\" y=\"$(num(-imag(point)-0.16))\" font-size=\"0.2\" fill=\"$color\">v$index</text>")
    end
    for marker in markers
        x=real(marker.point); y=-imag(marker.point)
        color=marker.rule===:P7 ? "#f79009" : "#d92d20"
        push!(body,"<circle cx=\"$(num(x))\" cy=\"$(num(y))\" r=\"0.22\" fill=\"none\" stroke=\"$color\" stroke-width=\"0.07\"/>")
        push!(body,"<path d=\"M $(num(x-0.16)) $(num(y-0.16)) L $(num(x+0.16)) $(num(y+0.16)) M $(num(x+0.16)) $(num(y-0.16)) L $(num(x-0.16)) $(num(y+0.16))\" stroke=\"$color\" stroke-width=\"0.06\"/>")
        push!(body,"<text x=\"$(num(x+0.28))\" y=\"$(num(y-0.25))\" font-size=\"0.24\" fill=\"$color\">$(marker.rule)</text>")
    end
    description=join(unique(d.message for d in failing)," | ")
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"$x0 $y0 $width $height\" role=\"img\" aria-label=\"$(xml_escape(name)) failure diagram\">")
        println(io,"<title>$(xml_escape(name)) — layout failure</title>")
        println(io,"<desc>$(xml_escape(description))</desc>")
        println(io,join(body,"\n"))
        println(io,"</svg>")
    end
    (name=name,svg=basename(path),rules=string.(unique(d.rule for d in failing)),
     edges=sort(collect(failing_edges)),markers=[Dict("rule"=>String(m.rule),"at"=>[real(m.point),imag(m.point)],"edges"=>m.edges,"message"=>m.message) for m in markers])
end

function build_failure_sheet(;output_dir=raw"C:/Users/Akaghef/AppData/Local/Temp/layout-failure-review")
    mkpath(output_dir)
    groups=run_ograph_plot_suite()
    targets=[
        ("projects/invariant/C260125MSTinv_uqsl2/connected_sum_2",groups["projects/invariant/C260125MSTinv_uqsl2"][2]),
        ("projects/invariant/C260125MSTinv_uqsl2/connected_sum_3",groups["projects/invariant/C260125MSTinv_uqsl2"][3]),
        ("projects/invariant/C260125MSTinv_uqsl2/connected_sum_2_again",groups["projects/invariant/C260125MSTinv_uqsl2"][7]),
        ("projects/invariant/C250818MSTinv_uqsl2/connected_sum_2",groups["projects/invariant/C250818MSTinv_uqsl2"][4]),
    ]
    records=Any[]
    for (name,call) in targets
        path=joinpath(output_dir,replace(name,'/' => '_')*".svg")
        push!(records,failure_svg(name,call.link,path))
    end
    cards=join(["""
      <article class="card">
        <h2>$(xml_escape(record.name))</h2>
        <img src="$(record.svg)" alt="$(xml_escape(record.name))">
        <div class="meta"><b>Rules:</b> $(join(record.rules,", "))<br><b>Edges:</b> $(join(record.edges,", "))</div>
        <pre>$(xml_escape(join(["$(m["rule"]) at $(m["at"]) edges=$(m["edges"]) — $(m["message"])" for m in record.markers],"\n")))</pre>
      </article>
      """ for record in records],"\n")
    html="""
    <!doctype html><html lang=\"ja\"><head><meta charset=\"utf-8\"><title>Layout failure review</title>
    <style>body{font-family:system-ui,sans-serif;background:#f6f8fb;color:#172033;margin:24px}.grid{display:grid;grid-template-columns:repeat(2,minmax(360px,1fr));gap:20px}.card{background:#fff;padding:16px;border:1px solid #d7ddea;border-radius:12px}h1{margin-top:0}h2{font-size:15px}.card img{width:100%;height:420px;object-fit:contain;background:#fffdf8;border:1px solid #edf0f4}.meta{margin:10px 0;color:#536078}pre{white-space:pre-wrap;font-size:11px;color:#b42318}@media(max-width:800px){.grid{grid-template-columns:1fr}}</style></head>
    <body><h1>Julia layout failure review</h1><p>赤=対象edge、赤×=P6接触、橙×=P7通過点。Layer8経路は使用していない。</p><div class=grid>$cards</div></body></html>
    """
    html_path=joinpath(output_dir,"review.html")
    write(html_path,html)
    open(joinpath(output_dir,"failure-status.json"),"w") do io JSON.print(io,records,2) end
    println("review_html=",html_path)
    println("failure_count=",length(records))
    html_path
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    build_failure_sheet()
end
