using JSON3, SHA
include("../../vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl")
include("../../src/heegaard/reference_surface.jl")
include("../../src/heegaard/marked_display.jl")
const M=MarkedSurfaceDiagrams
root=normpath(joinpath(@__DIR__,"../.."))
out=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(root,"artifacts/heegaard/reference")
mkpath(out)
for g in 1:3
    doc=lickorish_reference(g)
    M.write_json(joinpath(out,"LickorishG$g.marked.json"),doc)
    M.write_json(joinpath(out,"LickorishG$g.display.json"),prepare_marked_display(doc;cut_curves=["l$i" for i in 1:g]))
end
source=joinpath(root,"vendor/koda_genus3/data/figure6.fullmarked.json")
doc=JSON3.read(read(source,String),Dict{String,Any})
prepared=prepare_marked_display(doc;cut_curves=["L","M","R"])
prepared["source_file_sha256"]=bytes2hex(sha256(read(source)))
prepared["reference_status"]="supplied Koda marking retained; NOT silently re-labelled as Lickorish marking"
M.write_json(joinpath(out,"KodaPrepared.json"),prepared)
println("Koda: ",length(prepared["curves"])," curves, ",length(prepared["crossings"]),
    " crossings, ",length(prepared["boundaries"])," boundaries; D pieces ",length(prepared["region_marks"][1]["cut_pieces"]))
println("Reference and display complexes written to ",out)
