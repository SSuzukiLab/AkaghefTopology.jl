module MarkedParallelTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const M=H.MarkedSurfaceDiagrams

# Source: native builder, torus orientation in the established meridian/
# longitude fixture. A same-oriented parallel to a has intersection (b,a')=-1.
# This checks a generated annulus frontier, not a caller-provided copy route.
@testset "Constructed oriented annular auxiliary parallel copy" begin
    doc=B.grid_document(2)
    B.add_curve!(doc,"a","alpha",true,[[(0,1),(2,1)]])
    B.add_curve!(doc,"b","beta",true,[[(1,0),(1,2)]])
    B.declare_crossing!(doc,"x",1,1,"a","b",1)
    doc["metadata"]=Dict("source"=>"native oriented torus","display"=>Dict("infinity_triangle"=>["P",0]))
    original=deepcopy(doc)
    spec=Dict("curve"=>"a","id"=>"a_parallel","side"=>"left")
    result=H.marked_parallel_copy(doc,spec)
    @test doc==original
    @test H.validate(result.document;level=:surface).valid
    @test H.surface_genus(result.document)==1
    @test H.intersection_matrix(result.document)==[1;;]
    @test length(result.document["curves"])==3
    @test length(result.document["crossings"])==2
    new=only(x for x in result.document["crossings"] if x["id"]!="x")
    @test new["curves"]==["b","a_parallel"] && new["sign"]==-1
    @test all(!("a" in x["curves"] && "a_parallel" in x["curves"]) for x in result.document["crossings"])
    st=H._marked_states(result.document).geometry
    @test M.family_cut_report(st,"aux")["cut_system"]
    words=H.marked_curve_words(result.document)
    abelian(cid)=begin
        v=zeros(Int,length(words.presentation.generators))
        for (i,s) in words.curve_words[cid];v[i]+=s;end
        v
    end
    @test abelian("a")==abelian("a_parallel")
    @test result.document["metadata"]["source"]==doc["metadata"]["source"]
    # Every inverse step collapses its actual generated curve/mesh data.
    inverse=JSON3.read(JSON3.write(result.witness["inverse_specification"]),Dict{String,Any})
    @test H.marked_remove_parallel(result.document,inverse).document==doc
    @test_throws ArgumentError H.marked_parallel_copy(doc,Dict("curve"=>"a","id"=>"a","side"=>"left"))
    @test_throws ArgumentError H.marked_parallel_copy(doc,Dict("curve"=>"a","id"=>"c","side"=>"diagonal"))
end
end
