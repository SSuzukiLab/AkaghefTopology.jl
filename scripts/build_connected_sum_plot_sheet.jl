using Topology

const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"examples","OgraphPlotSuite.jl"))

_escape(s)=replace(string(s),"&"=>"&amp;","<"=>"&lt;",">"=>"&gt;","\""=>"&quot;")

function build_connected_sum_plot_sheet(;output_dir=raw"C:/Users/Akaghef/AppData/Local/Temp/connected-sum-plots")
    mkpath(output_dir)
    groups=run_ograph_plot_suite()
    selections=[
        ("projects/invariant/C260125MSTinv_uqsl2","connected_sum_2",2),
        ("projects/invariant/C260125MSTinv_uqsl2","connected_sum_3",3),
        ("projects/invariant/C260125MSTinv_uqsl2","connected_sum_2_again",7),
        ("projects/invariant/C250818MSTinv_uqsl2","connected_sum_2",4),
    ]
    cards=String[]
    for (group,name,index) in selections
        call=groups[group][index]
        filename=replace(group*"_"*name,'/' => '_')*".svg"
        plot_svg(call.link,joinpath(output_dir,filename);bending_numbers=call.bending)
        push!(cards,"<figure><figcaption>$(_escape(group*" / "*name))</figcaption><img src=\"$filename\" alt=\"$(_escape(name))\"></figure>")
    end
    html="""
    <!doctype html><html lang="ja"><head><meta charset="utf-8"><title>Connected sum plots</title>
    <style>body{margin:24px;background:#f6f8fb;font-family:system-ui,sans-serif;color:#172033}.grid{display:grid;grid-template-columns:repeat(2,minmax(360px,1fr));gap:24px}figure{margin:0;padding:12px;background:white;border:1px solid #d7ddea;border-radius:10px}figcaption{font-weight:700;margin-bottom:10px;font-size:14px}img{display:block;width:100%;height:480px;object-fit:contain;background:#fffdf8}@media(max-width:800px){.grid{grid-template-columns:1fr}}</style></head><body><div class="grid">$(join(cards,"\n"))</div></body></html>
    """
    html_path=joinpath(output_dir,"plots.html")
    write(html_path,html)
    println("html=",html_path)
    println("plots=",length(selections))
    html_path
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && build_connected_sum_plot_sheet()
