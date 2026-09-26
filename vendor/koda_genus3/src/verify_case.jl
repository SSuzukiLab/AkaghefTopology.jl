# Standalone Julia >= 1.10, stdlibs only. This script has NOT been run by ChatGPT.
using Test
include(joinpath(@__DIR__,"MarkedSurfaceDiagrams.jl"))
using .MarkedSurfaceDiagrams
const M=MarkedSurfaceDiagrams
const ROOT=normpath(joinpath(@__DIR__,".."))
const D_RIGHT=[("beta1","c1","c2"),("alpha3","c2","c3"),
 ("beta2","c3","c4"),("alpha2","c4","c5"),
 ("beta3","c5","c6"),("alpha1","c6","c1")]

function derive_odata(core)
    cycles=Dict(c["id"]=>c["cycle"] for family in ("alpha","beta") for c in core[family])
    removed=Set(D_RIGHT);corners=Set(x[2] for x in D_RIGHT)
    edges=Dict{Tuple{String,String},Tuple{String,String,String}}()
    outgoing=Dict{String,Vector{Tuple{String,String}}}()
    for (c,vs) in cycles
        for i in eachindex(vs)
            a,b=vs[i],vs[mod1(i+1,length(vs))]
            if (c,a,b) in removed;continue;end
            edges[(c,a)]=(c,a,b)
            push!(get!(outgoing,a,Tuple{String,String}[]),(c,a))
        end
    end
    @assert all(length(outgoing[x])==1 for x in corners)
    start=("alpha1","c1");z=start;seen=Set{Tuple{String,String}}();word=Int[]
    while !(z in seen)
        push!(seen,z);c,a,b=edges[z]
        if !(a in corners)
            push!(word,(startswith(c,"alpha") ? -1 : 1)*parse(Int,a[2:end]))
        end
        z=b in corners ? only(outgoing[b]) : (c,b)
    end
    @assert z==start && length(seen)==length(edges)
    xs=Dict(x["id"]=>x for x in core["intersections"])
    ([word],[-xs["x$i"]["sign"] for i in 1:5])
end

@testset "Image-derived genus-3 joining-disk fixture" begin
    st=load_diagram(joinpath(ROOT,"data","figure6.heegaard.marked.toml"))
    @test st.surface_report["genus"]==3
    @test st.surface_report["chi"]== -4
    @test family_cut_report(st,"alpha")["cut_system"]
    @test family_cut_report(st,"beta")["cut_system"]
    @test length(st.doc["crossings"])==11
    @test check_joining_disk(st,["F_11",0])["koda_joining"]
    rs=regions(st)
    @test length(rs)==7
    @test all(r["normalized_disk"] for r in rs)
    odata=derive_odata(legacy_projection(st))
    @test odata==([[-1,-2,3,1,4,-5,-4,5,-3,2]],[-1,1,1,1,-1])
    @test odata!=([[-1,-2,3,1,-4,5,4,-5,-3,2]],[-1,1,1,1,-1])
    p=realize(st)
    @test certify_layout(st,p)
    println("Computed o-data = ",odata)
    println("The user-proposed o-data is intentionally NOT the regression expectation.")
end
