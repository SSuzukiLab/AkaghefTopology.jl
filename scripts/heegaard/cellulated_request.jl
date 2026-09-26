using JSON3
include("cellulated_payload.jl")
length(ARGS)==2 || error("usage: cellulated_request.jl REQUEST RESPONSE")
request=JSON3.read(read(ARGS[1],String),Dict{String,Any})
result=cellulated_payload(request["diagram"];layout=get(request,"layout",nothing))
MSD.write_json(ARGS[2],result)
