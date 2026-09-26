using Topology, JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
root=normpath(joinpath(@__DIR__,"../../artifacts/heegaard/full_curves"))
manifest=H._plain_value(JSON3.read(read(joinpath(root,"manifest.json"),String)))
for base in (isempty(ARGS) ? ["Abalone","PoincareD"] : ARGS)
    base in ("Abalone","PoincareD") || error("unknown example")
    ph=H.load_punctured(joinpath(root,base=="Abalone" ? "abalone/S3.punctured.hgd.json" : base*".punctured.hgd.json"))
    layout=H.load_layout(joinpath(root,base*"Auto.hgd.layout.json"),ph.diagram)
    marked=H.mark_heegaard(ph,layout)
    id=base*"Marked"; title=(base=="Abalone" ? "S³ · abalone + D" : "Poincaré + D₀")*" · Marked"
    payload=surface_payload(marked,layout;id,title,
        source="Current automatic representative explicitly selected; not a recovered historic marking",
        seed="Stored relative seam words with independently checked realization")
    H.save_marked(joinpath(root,id*".hgd.json"),marked)
    write(joinpath(root,id*".json"),JSON3.write(payload))
    write(joinpath(root,id*".marked-project.json"),JSON3.write(Dict(
        "schema"=>"akaghef.heegaard.project/1","diagram"=>H._json_dict(marked),"layout"=>H._layout_json(layout))))
    filter!(e->e["id"]!=id,manifest); push!(manifest,Dict("id"=>id,"title"=>title))
    write(joinpath(root,"manifest.json"),JSON3.write(manifest))
    println("GENERATED ",id," ",H.marked_hash(marked))
end
