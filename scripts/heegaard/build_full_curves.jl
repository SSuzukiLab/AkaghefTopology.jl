using Topology, JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
include("../../examples/HeegaardPoincareJoiningDisk.jl")
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const SOURCE=joinpath(ROOT,"artifacts","heegaard","app")
const OUTPUT=joinpath(ROOT,"artifacts","heegaard","full_curves")
mkpath(OUTPUT)

# The baseline's seven fixed-100 plots are local, versioned inputs. This build
# neither reads another worktree nor silently solves a different router.
entries=H._plain_value(JSON3.read(read(joinpath(SOURCE,"manifest.json"),String)))
for entry in entries
    id=entry["id"]
    diagram=H.certify_heegaard(H.load_heegaard(joinpath(SOURCE,"$id.hgd.json")))
    layout=H.load_layout(joinpath(SOURCE,"$id.hgd.layout.json"),diagram)
    get(get(layout.backend_trace,"plot",Dict()),"profile",nothing)=="fixed-100" ||
        error("$id is not the baseline's fixed-production plot")
    chart=H.surface_chart(diagram,layout)
    svg=H.render(chart;path=joinpath(OUTPUT,"$id.svg"))
    atlas=H._surface_chart_json(chart)
    open(io->JSON3.write(io,atlas),joinpath(OUTPUT,"$id.surface-chart.json"),"w")
    data=merge(entry,Dict("svg"=>svg,"base_svg"=>H.render(layout),
        "chart"=>atlas,"diagram"=>H._json_dict(diagram),
        "source_status"=>"certified by the existing software checks; source provenance retained"))
    open(io->JSON3.write(io,data),joinpath(OUTPUT,"$id.json"),"w")
    println("Rendered $id: genus=$(chart.genus), alpha=$(length(chart.alpha_cycles)), beta=$(length(chart.beta_cycles)), signed intersections=$(length(chart.intersections))")
end

# The same persisted punctured input accepted by the app's import endpoint.
# No call to the example-specific torus renderer or its coordinates.
punctured_path=joinpath(OUTPUT,"abalone","S3.punctured.hgd.json")
if isfile(punctured_path)
    object=H.load_punctured(punctured_path)
    cache=joinpath(OUTPUT,"Abalone.hgd.layout.json")
    layout=if isfile(cache)
        candidate=H.load_layout(cache,object.diagram)
        get(get(candidate.backend_trace,"plot",Dict()),"profile",nothing)=="fixed-100" ||
            error("Abalone cached layout is not a common fixed-100 plot")
        candidate
    else
        seed=H.layout_heegaard(object.diagram;backend=:pure_julia)
        H.save_layout(joinpath(OUTPUT,"Abalone-before.hgd.layout.json"),seed)
        after=H.plot(object.diagram,seed)
        H.save_layout(cache,after)
        after
    end
    data=surface_payload(object,layout;id="Abalone",title="S³ · abalone + D",
        source="Persisted S³ Heegaard diagram with joining disk",
        seed="Common bounded router and fixed-100 plot")
    data["visual_review"]="owner-confirmed"
    write(joinpath(OUTPUT,"Abalone.svg"),data["svg"])
    write(joinpath(OUTPUT,"Abalone.surface-chart.json"),JSON3.write(data["chart"]))
    write(joinpath(OUTPUT,"Abalone.json"),JSON3.write(data))
    push!(entries,Dict(k=>data[k] for k in ("id","title","genus","intersections","beta","source","seed")))
    println("Rendered Abalone through common surface_chart: joining disk region=",object.joining_disk.region)
end
let
    object=HeegaardPoincareJoiningDiskExample.fixture()
    # Same signed labelled incidences as the baseline; load_layout checks the
    # canonical hash and full geometry. This is a fixture, not a new router claim.
    layout=H.load_layout(joinpath(SOURCE,"Poincare.hgd.layout.json"),object.diagram)
    data=surface_payload(object,layout;id="PoincareD",title="Poincaré + D₀",
        source="HeegaardDiagramProgram · Poincaré例のD₀（x4, x6, y4, y1）",
        seed="既存Poincaréの検証済み固定100配置（個別fixture）")
    H.save_punctured(joinpath(OUTPUT,"PoincareD.punctured.hgd.json"),object)
    H.save_layout(joinpath(OUTPUT,"PoincareD.hgd.layout.json"),layout)
    write(joinpath(OUTPUT,"PoincareD.svg"),data["svg"])
    write(joinpath(OUTPUT,"PoincareD.surface-chart.json"),JSON3.write(data["chart"]))
    write(joinpath(OUTPUT,"PoincareD.json"),JSON3.write(data))
    push!(entries,Dict(k=>data[k] for k in ("id","title","genus","intersections","beta","source","seed")))
    println("Rendered PoincareD through common surface_chart: joining disk region=",object.joining_disk.region)
end
open(io->JSON3.write(io,entries),joinpath(OUTPUT,"manifest.json"),"w")
cp(joinpath(ROOT,"examples","HeegaardFullCurves.qmd"),joinpath(OUTPUT,"HeegaardFullCurves.qmd");force=true)

# A checkout-local kernel avoids the machine's existing Julia kernel, which
# may hard-code a different worktree. JUPYTER_PATH selects this directory.
kernel_dir=joinpath(OUTPUT,"jupyter","kernels","heegaard-full-curves")
mkpath(kernel_dir)
kernel=Dict("display_name"=>"Julia · Heegaard full curves","language"=>"julia",
    "argv"=>[first(Base.julia_cmd().exec),"-i","--color=yes","--startup-file=no",
        "--compiled-modules=existing","--pkgimages=no","-O0","--project=$ROOT",
        "-e","import IJulia; IJulia.run_kernel()","{connection_file}"],
    "env"=>Dict(),"interrupt_mode"=>"message")
open(io->JSON3.write(io,kernel),joinpath(kernel_dir,"kernel.json"),"w")
