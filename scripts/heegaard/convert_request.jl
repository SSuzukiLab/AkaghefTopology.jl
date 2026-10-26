using Topology, JSON3
const H=Topology.Heegaard
include("koda_payload.jl")
length(ARGS)==2 || error("usage: convert_request.jl REQUEST.json RESPONSE.json")
request=JSON3.read(read(ARGS[1],String))
object=if (haskey(request.diagram,:patches) && haskey(request.diagram,:surface)) ||
          get(request.diagram,:schema,"")=="marked-cut-view/1"
    JSON3.read(JSON3.write(request.diagram),Dict{String,Any})
else
    H._punctured_from_json(request.diagram)
end
object isa Union{H.PuncturedHeegaardDiagram,AbstractDict} || error("Koda conversion requires a joining disk")
result=koda_payload(object)
write(ARGS[2],JSON3.write(result))
