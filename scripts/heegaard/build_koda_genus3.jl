using JSON3
include("cellulated_payload.jl")
root=normpath(joinpath(@__DIR__,"../.."))
doc=JSON3.read(read(joinpath(root,"vendor/koda_genus3/data/figure6.heegaard.marked.json"),String),Dict{String,Any})
data=cellulated_payload(doc;id="KodaGenus3",title="Koda · genus 3 · joining disk D")
out=joinpath(root,"artifacts/heegaard/cellulated")
MSD.write_json(joinpath(out,"KodaGenus3.json"),data)
manifest=joinpath(out,"manifest.json")
entries=isfile(manifest) ? JSON3.read(read(manifest,String),Vector{Dict{String,Any}}) : Dict{String,Any}[]
filter!(e->e["id"]!="KodaGenus3",entries)
pushfirst!(entries,Dict("id"=>data["id"],"title"=>data["title"]))
MSD.write_json(manifest,entries)
println("built KodaGenus3: ",data["genus"]," / ",data["intersections"]," crossings")
