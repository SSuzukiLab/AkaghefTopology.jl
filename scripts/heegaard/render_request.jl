using Topology
using JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
length(ARGS)==2 || error("usage: render_request.jl REQUEST.json RESPONSE.json")
request=JSON3.read(read(ARGS[1],String))
mktempdir() do dir
    corepath=joinpath(dir,"input.hgd.json")
    write(corepath,JSON3.write(request.diagram))
    object=if get(request.diagram,:schema,"")=="akaghef.heegaard.marked/1"
        H.load_marked(corepath)
    elseif get(request.diagram,:schema,"")=="akaghef.heegaard.punctured/1"
        H.load_punctured(corepath)
    else
        H.certify_heegaard(H.load_heegaard(corepath))
    end
    diagram=H._base_diagram(object)
    unmarked=object isa H.MarkedHeegaardDiagram ? object.object : object
    view=Symbol(get(request,:view,object isa H.PuncturedHeegaardDiagram ? "surface" : "alpha_cut"))
    view in (:surface,:alpha_cut) || error("Unknown Heegaard view")
    unmarked isa H.PuncturedHeegaardDiagram && view!=:surface &&
        error("Select surface view to preserve joining disk D")
    before=if haskey(request,:layout) && request.layout !== nothing
        path=joinpath(dir,"input.hgd.layout.json")
        write(path,JSON3.write(request.layout)); H.load_layout(path,diagram)
    else
        object isa H.MarkedHeegaardDiagram && error("Marked input requires its realization layout. Import the marked project or supply the separate layout; no unmarked rerouting is allowed.")
        get(request,:mark,false) && error("Selecting a marking requires the displayed layout, not an automatic representative")
        H.layout_heegaard(diagram;backend=:pure_julia,router=view==:surface ? :planar : :auto)
    end
    if get(request,:mark,false)
        object isa H.MarkedHeegaardDiagram && error("Already marked; re-marking is not implicit")
        object=H.mark_heegaard(object,before)
    end
    # Surface charts accept any valid non-overlapping polar placement. Preserve
    # the legacy alpha-cut endpoint's existing placement contract.
    if view==:alpha_cut
        expected=H.default_puncture_placement(H.cut_along_alpha(diagram))
        all(any(p.id==q.id && p.center==q.center && p.radius==q.radius for q in expected)
            for p in before.punctures) || error("Please supply the standard 2×g paired-puncture layout.")
    end
    after=view==:surface && get(get(before.backend_trace,"plot",Dict()),"profile",nothing)=="fixed-100" ?
        before : H.plot(object isa H.MarkedHeegaardDiagram ? object : diagram,before)
    # Marked plot validates before/after solving; marked surface_payload calls
    # surface_chart, which validates even when a saved fixed-100 plot is reused.
    # Do not repeat expensive unmarked canonicalization at this outer layer.
    result=if view==:surface
        surface_payload(object,after;title=object isa H.PuncturedHeegaardDiagram ?
            "Imported Heegaard diagram + D" : "Imported Heegaard diagram",
            seed=haskey(request,:layout) && request.layout!==nothing ?
                "User-supplied layout, validated by the common plot pipeline" : "組合せ平面埋込み → 極座標ポート → 共通固定100描画（配置fixtureなし）")
    else
      Dict("id"=>"imported","title"=>"Imported Heegaard diagram",
        "genus"=>length(diagram.code.code.alpha),"intersections"=>length(diagram.code.code.intersections),
        "beta"=>length(diagram.code.code.beta),"source"=>"User-supplied diagram",
        "seed"=>haskey(request,:layout) && request.layout !== nothing ? "User-supplied feasible layout" : "Bounded fixed-endpoint router",
        "diagram"=>request.diagram,"layout"=>H._layout_json(after),"before"=>H._layout_json(before),
        "svg"=>H.render(after;backend=:svg),"before_svg"=>H.render(before;backend=:svg))
    end
    write(ARGS[2],JSON3.write(result))
end

