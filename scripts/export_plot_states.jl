using Topology

include(joinpath(@__DIR__,"..","examples","VLExample.jl"))
include(joinpath(@__DIR__,"..","examples","OgraphPlotSuite.jl"))

json(value::Bool)=value ? "true" : "false"
json(value::Number)=isfinite(value) ? string(value) : "null"
json(value::AbstractString)="\""*replace(replace(value,"\\"=>"\\\\"),"\""=>"\\\"")*"\""
json(value::AbstractVector)="["*join(json.(value),",")*"]"

function record(group,index,link)
    code,orientation=real_gauss_code(link)
    "{"*join([
        "\"group\":"*json(group),
        "\"index\":"*json(index),
        "\"rgauss\":"*json(code),
        "\"orientation\":"*json(orientation),
        "\"isWeighted\":"*json(link.is_weighted),
        "\"weights\":"*json(link.weights),
    ],",")*"}"
end

records=String[]
for (index,call) in enumerate(run_vl_plot_calls())
    push!(records,record("VLExample",index,call.link))
end
for (group,calls) in run_ograph_plot_suite()
    for (index,call) in enumerate(calls)
        push!(records,record(group,index,call.link))
    end
end

output=joinpath(@__DIR__,"..","artifacts","baseline","julia_plot_states.json")
mkpath(dirname(output))
open(output,"w") do io
    print(io,"[",join(records,","),"]")
end
println(output)
