using Topology, JSON3, LinearAlgebra
const H=Topology.Heegaard
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT=joinpath(ROOT,"artifacts","heegaard","app")
include("gallery_examples.jl")
mkpath(OUTPUT)
BLAS.set_num_threads(1)

function write_json(path,value)
    temporary=path*".tmp"
    open(temporary,"w") do io; JSON3.write(io,value); end
    mv(temporary,path;force=true)
end

function build_example(spec)
    id=spec.id
    diagram = id=="T3" ? H.certify_heegaard(H.load_heegaard(joinpath(ROOT,"test","fixtures","t3_octahedron_v1.hgd.json"))) :
        id=="Poincare" ? H.certify_heegaard(H.load_heegaard(joinpath(ROOT,"test","fixtures","poincare_homology_sphere_v1.hgd.json"))) :
        length(spec.factors)==1 ? lens(only(spec.factors)...) : lens_sum(spec.factors)
    old=joinpath(OUTPUT,"$id.json")
    saved=isfile(old) ? JSON3.read(read(old,String)) : nothing
    cache=joinpath(OUTPUT,"$id.hgd.layout.json")
    after=nothing
    if !("--fresh" in ARGS) && isfile(cache)
        try
            candidate=H.load_layout(cache,diagram)
            smoothing=get(candidate.backend_trace,"smoothing",Dict())
            if get(smoothing,"method",nothing)=="elastic-normal-boundary" &&
                get(get(smoothing,"parameters",Dict()),"bending",nothing)==0.04
                after=candidate
            end
        catch
            # Recompute an interrupted/stale result instead of publishing it.
        end
    end
    beforepath=joinpath(OUTPUT,"$id-before.hgd.layout.json")
    if isfile(beforepath)
        before=H.load_layout(beforepath,diagram)
    elseif saved!==nothing
        write_json(beforepath,saved.before)
        before=H.load_layout(beforepath,diagram)
    elseif id=="T3"
        before=H.load_layout(joinpath(ROOT,"test","fixtures","t3_octahedron_v1.hgd.layout.json"),diagram)
    elseif id=="Poincare"
        before=H.load_layout(joinpath(ROOT,"test","fixtures","poincare_homology_sphere_v1.hgd.layout.json"),diagram)
    else
        before=annulus_seed(diagram;narrow=length(spec.factors)>1)
    end
    after===nothing && (after=H.plot(diagram,before))
    homology=H.first_homology(diagram)
    genus=length(diagram.code.code.alpha)
    H.surface_genus(diagram)==genus || error("surface genus mismatch")
    homology.free_rank==(id=="T3" ? 3 : 0) || error("unexpected free homology")
    expected_order=id in ("Poincare","T3") ? 1 : prod(first,spec.factors;init=big(1))
    prod(homology.torsion;init=big(1))==expected_order || error("H1 order mismatch")
    before.report.valid && after.report.valid && isempty(after.crossings) || error("invalid geometry")
    H._exit_angle_degrees(after)<=15 || error("normal-exit constraint failed")
    H.save_layout(cache,after); H.save_layout(beforepath,before)
    H.save_heegaard(joinpath(OUTPUT,"$id.hgd.json"),diagram)
    H.render(before;backend=:svg,path=joinpath(OUTPUT,"$id-before.svg"))
    H.render(after;backend=:svg,path=joinpath(OUTPUT,"$id.svg"))
    source=spec.construction
    seed=get(before.backend_trace,"seed",nothing)=="analytic-annulus-family" ?
        "Constructive annulus-family seed; validated before relaxation" :
        id=="T3" ? "Source octahedron seed transported to 2×3 disks by an explicit annular twist" :
        id=="Poincare" ? "Saved fixture-scoped feasible layout" : "Saved bounded-router feasible layout"
    entry=Dict{String,Any}("id"=>id,"title"=>spec.title,"category"=>spec.category,"genus"=>genus,
        "intersections"=>length(diagram.code.code.intersections),"beta"=>length(diagram.code.code.beta),
        "source"=>source,"seed"=>seed,"factors"=>[collect(f) for f in spec.factors],
        "h1"=>Dict("free_rank"=>homology.free_rank,"torsion"=>homology.torsion),
        "diagram"=>H._plain_value(JSON3.read(read(joinpath(OUTPUT,"$id.hgd.json"),String))),
        "layout"=>H._layout_json(after),"before"=>H._layout_json(before),
        "svg"=>H.render(after;backend=:svg),"before_svg"=>H.render(before;backend=:svg))
    write_json(old,entry)
    Dict(k=>entry[k] for k in ("id","title","category","genus","intersections","beta","source","seed","factors","h1"))
end

function main()
    catalog=gallery_catalog(;extended="--extended" in ARGS)
    selected=findfirst(a->startswith(a,"--ids="),ARGS)
    ids=selected===nothing ? nothing : Set(split(split(ARGS[selected],"=";limit=2)[2],","))
    todo=ids===nothing ? catalog : filter(e->e.id in ids,catalog)
    manifest=joinpath(OUTPUT,"manifest.json")
    existing=isfile(manifest) ? H._plain_value(JSON3.read(read(manifest,String))) : Any[]
    entries=Dict(e["id"]=>e for e in existing)
    # Preserve existing variant files/links without rebuilding or promoting them.
    catalog_order=gallery_catalog(;extended=true)
    published_catalog="--extended" in ARGS ? catalog_order : catalog
    failures=Dict{String,String}(); guard=ReentrantLock(); completed=0
    Threads.@threads for spec in todo
        try
            entry=build_example(spec)
            lock(guard) do
                entries[spec.id]=entry; completed+=1
                write_json(manifest,[entries[s.id] for s in published_catalog if haskey(entries,s.id)])
                println("Rendered $(spec.id) ($completed/$(length(todo)))"); flush(stdout)
            end
        catch exception
            lock(guard) do
                failures[spec.id]=sprint(showerror,exception)
                println("FAILED $(spec.id): ",failures[spec.id]); flush(stdout)
            end
        end
    end
    write_json(joinpath(OUTPUT,"build-status.json"),Dict("requested"=>length(todo),"completed"=>completed,"failures"=>failures))
    cp(joinpath(ROOT,"examples","HeegaardPlot.qmd"),joinpath(OUTPUT,"HeegaardPlot.qmd");force=true)
    cp(joinpath(ROOT,"examples","HeegaardT3.qmd"),joinpath(OUTPUT,"HeegaardT3.qmd");force=true)
    isempty(failures) || error("$(length(failures)) examples failed; see build-status.json")
    println("Selected gallery ready")
end
abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()
