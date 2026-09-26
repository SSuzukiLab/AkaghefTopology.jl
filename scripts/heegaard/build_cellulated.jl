include("cellulated_payload.jl")
root=normpath(joinpath(@__DIR__,"../.."))
source=joinpath(root,"vendor/marked_surface_diagrams/examples")
out=joinpath(root,"artifacts/heegaard/cellulated")
mkpath(out)
entries=[]
for name in sort(filter(n->endswith(n,".toml"),readdir(source)))
    id=splitext(name)[1]
    data=cellulated_payload(MSD.TOML.parsefile(joinpath(source,name));id=id,title=id)
    MSD.write_json(joinpath(out,id*".json"),data)
    push!(entries,Dict("id"=>id,"title"=>id))
    println("built ",id); flush(stdout)
end
# Explicitly derived region selection; not an abalone identification or Koda D.
doc=MSD.TOML.parsefile(joinpath(source,"torus_three_crossings.toml"))
doc["region_marks"]=[Dict("id"=>"D0","seed"=>["P",68],"claim"=>"embedded_closed_disk")]
data=cellulated_payload(doc;id="torus_three_crossings_D",title="3 crossings + D₀ (geometric disk)")
MSD.write_json(joinpath(out,data["id"]*".json"),data)
push!(entries,Dict("id"=>data["id"],"title"=>data["title"]))
MSD.write_json(joinpath(out,"manifest.json"),entries)
