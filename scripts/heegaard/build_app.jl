using Topology
using JSON3
const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT = joinpath(ROOT,"artifacts","heegaard","app")
mkpath(OUTPUT)

function lens(p,q)
    code=H.heegaard_code([collect(1:p)],[[mod1(1+(i-1)*q,p) for i in 1:p]],
        [H.Intersection(i,1,1,1) for i in 1:p])
    H.certify_heegaard(H.decorate(code,[[i] for i in 1:p]))
end

examples=[("Poincare","Poincaré homology sphere",nothing),
          ("S3","S³",(1,1)),("L2_1","L(2,1)",(2,1)),
          ("L5_1","L(5,1)",(5,1)),("L5_2","L(5,2)",(5,2)),
          ("L7_1","L(7,1)",(7,1)),("L7_2","L(7,2)",(7,2))]
entries=Dict{String,Any}[]
for (id,title,pq) in examples
    diagram = pq === nothing ? H.certify_heegaard(H.load_heegaard(joinpath(ROOT,
        "test","fixtures","poincare_homology_sphere_v1.hgd.json"))) : lens(pq...)
    source = pq === nothing ? joinpath(ROOT,"test","fixtures","poincare_homology_sphere_v1.hgd.layout.json") :
        joinpath(ROOT,"artifacts","heegaard","review_alpha_cut","$(id)__alpha_cut.hgd.layout.json")
    before=H.load_layout(source,diagram)
    # Reuse only an explicitly requested current result, not an old preview.
    cache=joinpath(OUTPUT,"$id.hgd.layout.json")
    cached = if "--reuse" in ARGS && isfile(cache)
        try
            candidate=H.load_layout(cache,diagram)
            smoothing=get(candidate.backend_trace,"smoothing",Dict())
            get(smoothing,"method",nothing)=="elastic-normal-boundary" &&
                get(get(smoothing,"parameters",Dict()),"bending",nothing)==0.04 ? candidate : nothing
        catch
            nothing # An interrupted previous build may have left an incomplete file.
        end
    else
        nothing
    end
    after = cached===nothing ? H.plot(diagram,before) : cached
    H.save_layout(cache,after)
    H.save_heegaard(joinpath(OUTPUT,"$id.hgd.json"),diagram)
    H.render(before;backend=:svg,path=joinpath(OUTPUT,"$id-before.svg"))
    H.render(after;backend=:svg,path=joinpath(OUTPUT,"$id.svg"))
    entry=Dict{String,Any}("id"=>id,"title"=>title,"genus"=>length(diagram.code.code.alpha),
        "intersections"=>length(diagram.code.code.intersections),"beta"=>length(diagram.code.code.beta),
        "source"=>pq===nothing ? "Klebanov Fig. 2.15 · input code transcribed; region decoration is a candidate" :
                                 "Lens-space cyclic-word fixture",
        "seed"=>pq===nothing ? "Saved fixture-scoped feasible layout" : "Saved bounded-router feasible layout",
        "diagram"=>H._plain_value(JSON3.read(read(joinpath(OUTPUT,"$id.hgd.json"),String))),
        "layout"=>H._layout_json(after),"before"=>H._layout_json(before),
        "svg"=>H.render(after;backend=:svg),"before_svg"=>H.render(before;backend=:svg))
    open(joinpath(OUTPUT,"$id.json"),"w") do io; JSON3.write(io,entry); end
    push!(entries,Dict(k=>entry[k] for k in ("id","title","genus","intersections","beta","source","seed")))
    println("Rendered $id: ",after.report.valid," · ",after.backend_trace["smoothing"]["stop_reason"])
end
open(joinpath(OUTPUT,"manifest.json"),"w") do io; JSON3.write(io,entries); end
cp(joinpath(ROOT,"examples","HeegaardPlot.qmd"),joinpath(OUTPUT,"HeegaardPlot.qmd");force=true)
