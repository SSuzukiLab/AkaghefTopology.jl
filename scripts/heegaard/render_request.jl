using Topology
using JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
length(ARGS)==2 || error("usage: render_request.jl REQUEST.json RESPONSE.json")
request=JSON3.read(read(ARGS[1],String))
if get(request.diagram,:schema,"")=="marked-cut-view/1"
    view=JSON3.read(JSON3.write(request.diagram),Dict{String,Any})
    result=H.plot(view)
    write(ARGS[2],JSON3.write(result))
    exit()
end
if haskey(request.diagram,:patches) && haskey(request.diagram,:surface)
    doc=JSON3.read(JSON3.write(request.diagram),Dict{String,Any})
    options=get(request,:display,Dict())
    result=H.plot(doc;cut_curves=get(options,:cut_curves,nothing),
        outer_triangle=get(options,:infinity_triangle,nothing),visible_curves=get(options,:visible_curves,nothing))
    write(ARGS[2],JSON3.write(result))
    exit()
end
mktempdir() do dir
    corepath=joinpath(dir,"input.hgd.json")
    write(corepath,JSON3.write(request.diagram))
    object=if get(request.diagram,:schema,"")=="akaghef.heegaard.punctured/1"
        H.load_punctured(corepath)
    else
        H.certify_heegaard(H.load_heegaard(corepath))
    end
    diagram=H._base_diagram(object)
    view=Symbol(get(request,:view,"alpha_cut"))
    view in (:surface,:alpha_cut) || error("Unknown Heegaard view")
    before=if haskey(request,:layout) && request.layout !== nothing
        path=joinpath(dir,"input.hgd.layout.json")
        write(path,JSON3.write(request.layout)); H.load_layout(path,diagram)
    else
        H.layout_heegaard(diagram;backend=:pure_julia,router=:planar)
    end
    # Surface charts accept any valid non-overlapping polar placement. Preserve
    # the legacy alpha-cut endpoint's existing placement contract.
    if view==:alpha_cut
        expected=H.default_puncture_placement(H.cut_along_alpha(diagram))
        all(any(p.id==q.id && p.center==q.center && p.radius==q.radius for q in expected)
            for p in before.punctures) || error("Please supply the standard 2×g paired-puncture layout.")
    end
    after=get(get(before.backend_trace,"plot",Dict()),"profile",nothing)=="fixed-100" ?
        before : H.plot(diagram,before)
    result=if view==:surface || object isa H.PuncturedHeegaardDiagram || get(request,:surface_payload,false)
        surface_payload(object,after;view,title=object isa H.PuncturedHeegaardDiagram ?
            "Imported Heegaard diagram + D" : "Imported Heegaard diagram",
            seed=haskey(request,:layout) && request.layout!==nothing ?
                "User-supplied layout, validated by the common plot pipeline" : "組合せ平面埋込み → 極座標ポート → 共通固定100描画（配置fixtureなし）")
    else
      Dict("id"=>"imported","title"=>"Imported Heegaard diagram",
        "genus"=>length(diagram.code.code.alpha),"intersections"=>length(diagram.code.code.intersections),
        "beta"=>length(diagram.code.code.beta),"source"=>"User-supplied diagram",
        "seed"=>haskey(request,:layout) && request.layout !== nothing ? "User-supplied feasible layout" : "Combinatorial planar embedding and fixed-100 plot",
        "diagram"=>request.diagram,"layout"=>H._layout_json(after),"before"=>H._layout_json(before),
        "svg"=>H.render(after;backend=:svg),"before_svg"=>H.render(before;backend=:svg))
    end
    write(ARGS[2],JSON3.write(result))
end

