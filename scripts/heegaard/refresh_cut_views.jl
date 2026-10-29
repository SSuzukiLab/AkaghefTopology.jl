# Refresh view payloads from validated saved layouts. This is not a routing test;
# fixture-free generation is exercised separately through render_request.jl.
using Topology, JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
out=normpath(joinpath(@__DIR__,"../../artifacts/heegaard/full_curves"))
for entry in JSON3.read(read(joinpath(out,"manifest.json"),String))
    id=String(entry.id)
    !isempty(ARGS) && !(id in ARGS) && continue
    path=joinpath(out,id*".json")
    old=JSON3.read(read(path,String),Dict{String,Any})
    mktempdir() do dir
        core=joinpath(dir,"input.json"); write(core,JSON3.write(old["diagram"]))
        object=old["diagram"]["schema"]=="akaghef.heegaard.punctured/1" ?
            H.load_punctured(core) : H.certify_heegaard(H.load_heegaard(core))
        lp=joinpath(dir,"layout.json")
        write(lp,JSON3.write(get(old,"layout",old["chart"]["base_layout"])))
        layout=H.load_layout(lp,H._base_diagram(object))
        payload=surface_payload(object,layout;id,title=old["title"],source=old["source"],seed=old["seed"])
        # Retain existing conversion results/review flags; only refresh view fields.
        merge!(old,payload)
        write(path,JSON3.write(old))
        write(joinpath(out,id*"-alpha-cut.svg"),old["alpha_cut_svg"])
        println("REFRESHED ",id); flush(stdout)
    end
end
