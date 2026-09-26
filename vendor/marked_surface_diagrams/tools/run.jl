#!/usr/bin/env julia
include(joinpath(@__DIR__,"..","src","MarkedSurfaceDiagrams.jl"))
using .MarkedSurfaceDiagrams
function main(args)
    isempty(args) && error("usage: julia tools/run.jl INPUT.toml [--svg OUT.svg] [--report OUT.json]")
    st=load_diagram(args[1]);pos=realize(st);r=report(st,pos);destination=nothing;i=2
    while i<=length(args)
        i<length(args) || error("missing value for $(args[i])")
        if args[i]=="--svg";r["svg"]=write_svg(st,pos,args[i+1])
        elseif args[i]=="--report";destination=args[i+1]
        else;error("unknown option $(args[i])");end
        i+=2
    end
    if destination===nothing;println(MarkedSurfaceDiagrams.json(r))
    else;write_json(destination,r);end
end
main(ARGS)
