module MarkedSubdivisionTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const B=H.MarkedDiagramBuilders
include("marked_fixtures.jl")
jsonround(x)=JSON3.read(JSON3.write(x),Dict{String,Any})
fixture(name)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})

# Source: existing independent native fixture/oracle family. Subdivision is
# checked as a full marked identity transport, not merely by equal topology.
@testset "Native barycentric subdivision and genuine generated collapse" begin
    for name in ("torus_mer_lon","torus_parallel_zero_crossings","nested_zero_crossings","sphere_empty")
        doc=fixture(name);snapshot=deepcopy(doc)
        doc["metadata"]=Dict("source"=>Dict("fixture"=>name,"decimal"=>0.125),"display"=>Dict("infinity_triangle"=>[doc["patches"][1]["face"],0]))
        original=deepcopy(doc);refined=H.marked_subdivide(doc)
        @test doc==original
        @test H.validate(refined.document;level=:surface).valid
        @test H.surface_genus(refined.document)==H.surface_genus(doc)
        @test H.intersection_matrix(refined.document)==H.intersection_matrix(doc)
        @test refined.document["crossings"]==doc["crossings"]
        @test refined.document["metadata"]["source"]==doc["metadata"]["source"]
        @test sum(length(p["triangles"]) for p in refined.document["patches"])==6sum(length(p["triangles"]) for p in doc["patches"])
        @test H.replay_marked_move(doc,jsonround(refined.witness)).document==refined.document
        inverse=H.inverse_marked_move(refined.witness)
        @test H.replay_marked_move(refined.document,jsonround(inverse)).document==doc
        @test H.marked_collapse_subdivision(refined.document,jsonround(refined.witness["inverse_specification"])).document==doc
        @test H.replay_marked_move(doc,H.inverse_marked_move(inverse)).document==refined.document
        altered=deepcopy(refined.document);altered["metadata"]["source"]["fixture"]="changed"
        @test_throws ArgumentError H.marked_collapse_subdivision(altered,refined.witness["inverse_specification"])
    end

    # D is an alpha-beta region even when the aux loop/arc cuts it visually.
    split=B.grid_document(6)
    B.add_curve!(split,"a","alpha",true,[[(1,1),(5,1),(5,5),(1,5),(1,1)]])
    B.add_curve!(split,"cut","aux",true,[[(3,0),(3,6)]])
    B.declare_crossing!(split,"p",3,1,"a","cut",1)
    B.declare_crossing!(split,"q",3,5,"a","cut",-1)
    split["region_marks"]=[Dict("id"=>"D","seed"=>["P",28],"claim"=>"embedded_closed_disk")]
    split["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",28]))
    refined=H.marked_subdivide(split)
    @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in refined.document["patches"])==241
    @test refined.document["region_marks"][1]["seed"]==["P",168]
    @test length(refined.witness["correspondence"]["regions"])==2
    @test length(refined.witness["correspondence"]["geometric_regions"])==3
    @test H.replay_marked_move(refined.document,H.inverse_marked_move(refined.witness)).document==split
    forged=deepcopy(refined.witness["inverse_specification"])
    forged["patches"][1]["triangles"][1]["corners"]=reverse(forged["patches"][1]["triangles"][1]["corners"])
    @test_throws Exception H.marked_collapse_subdivision(refined.document,forged)
    # Resource rejection preserves the source and does not relax the vendor.
    too_large=B.grid_document(8);snapshot=deepcopy(too_large)
    @test_throws ArgumentError H.marked_subdivide(too_large)
    @test too_large==snapshot
    @test_throws ArgumentError H.marked_subdivide(split,Dict("replacement"=>"not allowed"))
end

@testset "Constructed quotient vertex-star atlas and verified route inverse" begin
    for doc in (fixture("torus_mer_lon"),fixture("torus_parallel_zero_crossings"),capped_nested_sphere())
        doc["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
        original=deepcopy(doc);stars=H.marked_dual_star(doc)
        @test doc==original
        @test H.validate(stars.document;level=:surface).valid
        @test H.surface_genus(stars.document)==H.surface_genus(doc)
        @test H.intersection_matrix(stars.document)==H.intersection_matrix(doc)
        @test sum(length(p["triangles"]) for p in stars.document["patches"])==6sum(length(p["triangles"]) for p in doc["patches"])
        @test all(length(p["side_order"])>=3 for p in stars.document["patches"])
        @test length(stars.witness["correspondence"]["regions"])==length(H.complement_regions(doc))
        @test H.replay_marked_move(doc,jsonround(stars.witness)).document==stars.document
        @test H.replay_marked_move(stars.document,jsonround(H.inverse_marked_move(stars.witness))).document==doc
    end
    # Original triangle corners on this two-face sphere have valence two;
    # subdivision first gives valid disk stars without inventing side corners.
    sphere=fixture("sphere_empty")
    @test_throws ArgumentError H.marked_dual_star(sphere)
    fine=H.marked_subdivide(sphere).document
    @test H.validate(H.marked_dual_star(fine).document;level=:surface).valid
    @test_throws ArgumentError H.marked_dual_star(B.grid_document(4;surface=:disk))
    # The initial failed batch put this physical-boundary disk in the closed
    # positive loop. Preserve its unsupported evidence, not a changed meaning.
    disk=fixture("nested_zero_crossings")
    failure=try H.marked_dual_star(disk);nothing catch err;err end
    @test failure isa ArgumentError && occursin("closed surface",sprint(showerror,failure))
    sphere=capped_nested_sphere()
    @test sphere["patches"][1]==disk["patches"][1]
    @test sphere["curves"]==disk["curves"]
    sphere_state=H._marked_states(sphere).geometry
    @test sphere_state.surface_report["boundary_components"]==0
    @test length(unique(M.root!(sphere_state.vertex_dsu,v) for v in keys(sphere_state.vertex_dsu.p)))==50
    @test_throws ArgumentError H.marked_dual_star(B.grid_document(9))
end
end
