using Topology, JSON3
const H=Topology.Heegaard
include("koda_payload.jl")
length(ARGS)==2 || error("usage: convert_request.jl REQUEST.json RESPONSE.json")
request=JSON3.read(read(ARGS[1],String))
if haskey(request.diagram,:patches) || get(request.diagram,:schema,"")=="akaghef.heegaard.cellulated-project/1"
    raw=JSON3.read(JSON3.write(request.diagram),Dict{String,Any})
    doc=haskey(raw,"patches") ? raw : raw["diagram"]
    projection=H.cellulated_punctured(doc;disk_id=get(request,:disk_id,nothing))
    result=koda_payload(projection.object)
    result["cellulated_projection"]=projection.provenance
    write(ARGS[2],JSON3.write(result))
    exit()
end
marked=get(request.diagram,:schema,"")=="akaghef.heegaard.marked/1" ? H._marked_from_json(request.diagram) : nothing
object=marked===nothing ? H._punctured_from_json(request.diagram) : marked.object
object isa H.PuncturedHeegaardDiagram || error("Koda conversion requires a joining disk")
result=koda_payload(object)
if marked!==nothing
    result["source_marked_hash"]=H.marked_hash(marked)
    result["marking_realization"]="not checked by this combinatorial conversion"
    result["marking_policy"]="Koda construction uses the underlying (H,D); marking identity retained as provenance"
end
write(ARGS[2],JSON3.write(result))
