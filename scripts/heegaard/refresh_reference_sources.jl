#= Refresh literature provenance without changing accepted geometry or topology.

The managed crop remains in its library location. The script verifies its bytes,
then checks that regenerating the reference changes only metadata. Saved views,
seeds and prepared complexes keep their original coordinates and edge words.
=#
using Topology, JSON3, SHA
const H=Topology.Heegaard
length(ARGS)==2 || error("usage: refresh_reference_sources.jl LIBRARY_ROOT OUTPUT_DIRECTORY")
library,out=abspath.(ARGS)
readjson(p)=JSON3.read(read(p,String),Dict{String,Any})
for g in 1:4
    name="LickorishG$g"
    doc=H.lickorish_reference(g)
    figure=doc["metadata"]["source_figure"]
    asset=joinpath(library,figure["asset"])
    bytes2hex(sha256(read(asset)))==figure["asset_sha256"] || error("Managed figure differs: $asset")
    path=joinpath(out,"$name.marked.json")
    old=readjson(path)
    a,b=deepcopy(old),deepcopy(doc)
    delete!(a,"metadata");delete!(b,"metadata")
    a==b || error("Reference mathematical data changed; provenance-only refresh refused")
    prepared=H.prepare_marked_display(doc;cut_curves=["l$i" for i in 1:g])
    before=readjson(joinpath(out,"$name.display.json"))
    for key in keys(before)
        key=="source_representation_sha256" && continue
        before[key]==prepared[key] || error("Prepared complex changed at $key")
    end
    write(path,JSON3.write(doc))
    write(joinpath(out,"$name.display.json"),JSON3.write(prepared))
    for suffix in (".view.json","Common.json","Common.seed.json")
        viewpath=joinpath(out,name*suffix)
        view=readjson(viewpath)
        view["source_sha256"]=bytes2hex(sha256(read(path)))
        view["source_figure"]=figure
        view["source_status"]=doc["metadata"]["source_status"]
        write(viewpath,JSON3.write(view))
    end
    println(name,": managed crop hash checked; source metadata refreshed; mathematical data and geometry unchanged")
end
