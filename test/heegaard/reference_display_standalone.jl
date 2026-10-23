using Test, JSON3
include("../../vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl")
include("../../vendor/marked_surface_diagrams/src/builders.jl")
include("../../src/heegaard/reference_surface.jl")
include("../../src/heegaard/marked_display.jl")
const M=MarkedSurfaceDiagrams
@testset "fixed reference and geometry-free display" begin
    for g in (1,2,3,4)
        doc=lickorish_reference(g);st=M.validate(doc)
        @test st.surface_report["genus"]==g
        @test length(doc["crossings"])==3g-2
        @test length(M.regions(st))==g
        @test all(r["normalized_disk"] for r in M.regions(st))
        p=prepare_marked_display(doc;cut_curves=["l$i" for i in 1:g])
        @test length(p["boundaries"])==2g
        @test Set((b["row"],b["column"]) for b in p["boundaries"])==Set((r,c) for r in 0:1 for c in 0:g-1)
        @test isempty(p["curves"]) # no forced alpha/beta family, not even r=s=g
    end
    @test_throws M.DiagramError prepare_marked_display(Dict("genus"=>3);cut_curves=["L","M","R"])
    src=joinpath(@__DIR__,"../../vendor/koda_genus3/data/figure6.fullmarked.json")
    doc=JSON3.read(read(src,String),Dict{String,Any});before=M.representation_digest(doc)
    p=prepare_marked_display(doc;cut_curves=["L","M","R"])
    @test M.representation_digest(doc)==before
    @test length(p["curves"])==6 && length(p["crossings"])==11
    @test p["region_marks"][1]["source_region"]["embedded_closed_disk"]
    @test length(p["region_marks"][1]["cut_pieces"])==1
    @test_throws M.DiagramError prepare_marked_display(doc;cut_curves=["alpha1","alpha2","alpha3"])
    @test_throws M.DiagramError prepare_marked_display(doc;cut_curves=["L","L","R"])
    # Deleting a curve does not change the reference surface or force a cut.
    reduced=deepcopy(doc)
    filter!(c->c["id"]!="alpha3",reduced["curves"])
    filter!(x->!("alpha3" in x["curves"]),reduced["crossings"])
    delete!(reduced,"region_marks")
    q=prepare_marked_display(reduced;cut_curves=["L","M","R"])
    @test count(c->c["role"]=="alpha",q["curves"])==2
    @test count(c->c["role"]=="beta",q["curves"])==3
    @test q["cut_topology"]==p["cut_topology"]
    # A marked disk straddles the auxiliary cut: it must stay ONE semantic D,
    # represented by two cut pieces. Checking D against all aux curves is wrong.
    B=MarkedDiagramBuilders; d=B.grid_document(8)
    B.add_curve!(d,"a","alpha",true,[[(2,2),(6,2),(6,6),(2,6),(2,2)]])
    B.add_curve!(d,"cut","aux",true,[[(4,0),(4,8)]])
    B.declare_crossing!(d,"p",4,2,"a","cut",1)
    B.declare_crossing!(d,"q",4,6,"a","cut",-1)
    d["region_marks"]=[Dict("id"=>"D","seed"=>["P",54],"claim"=>"embedded_closed_disk")]
    p=prepare_marked_display(d;cut_curves=["cut"])
    @test p["region_marks"][1]["source_region"]["embedded_closed_disk"]
    @test length(p["region_marks"][1]["cut_pieces"])==2
    @test sum(length(r["triangles"]) for r in p["region_marks"][1]["cut_pieces"])==32
end
