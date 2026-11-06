using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const M=H.MarkedSurfaceDiagrams
loadcase(name)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})

@testset "actual marked compression topology" begin
    # Same genus and curve count, different actual compression boundaries.
    essential=B.grid_document(6)
    B.add_curve!(essential,"a","alpha",true,[[(0,1),(6,1)]])
    B.add_curve!(essential,"b","beta",true,[[(0,5),(6,5)]])
    r=H.compression_boundary_report(essential)
    @test r.closed && r.alpha_complete && r.beta_complete
    @test r.alpha_boundary_genera==[0]
    contractible=B.grid_document(6)
    B.add_curve!(contractible,"a","alpha",true,[[(1,1),(5,1),(5,5),(1,5),(1,1)]])
    B.add_curve!(contractible,"b","beta",true,[[(2,2),(4,2),(4,4),(2,4),(2,2)]])
    bad=H.compression_boundary_report(contractible)
    @test bad.surface_genus==r.surface_genus==1
    @test length(bad.alpha.curves)==length(r.alpha.curves)==1
    @test sort(bad.alpha_boundary_genera)==[0,1]
    @test !bad.closed && !bad.alpha_complete

    # Redundant parallel meridians compress to two spheres, not one torus.
    redundant=deepcopy(essential)
    B.add_curve!(redundant,"a2","alpha",true,[[(0,3),(6,3)]])
    r=H.compression_boundary_report(redundant)
    @test r.closed && !r.alpha_complete
    @test r.alpha_boundary_genera==[0,0]
    @test sum(c["attached_caps"] for c in r.alpha.components)==4
    @test H.compression_boundary_report(loadcase("sphere_empty")).closed
    @test H.is_closed_heegaard(loadcase("torus_mer_lon"))

    # An original boundary is never capped by the compression operation.
    annulus=B.grid_document(4;surface=:annulus)
    B.add_curve!(annulus,"a","alpha",true,[[(0,2),(4,2)]])
    r=H.compression_boundary_report(annulus)
    @test !r.closed
    @test [c["boundary_components"] for c in r.alpha.components]==[1,1]

    # D is a region of alpha U beta, even if an auxiliary cut crosses it.
    d=B.grid_document(8)
    B.add_curve!(d,"a","alpha",true,[[(2,2),(6,2),(6,6),(2,6),(2,2)]])
    B.add_curve!(d,"cut","aux",true,[[(4,0),(4,8)]])
    B.declare_crossing!(d,"p",4,2,"a","cut",1)
    B.declare_crossing!(d,"q",4,6,"a","cut",-1)
    d["region_marks"]=[Dict("id"=>"D","seed"=>["P",54],"claim"=>"embedded_closed_disk")]
    before=M.representation_digest(d)
    report=H.marked_topology_report(d)
    @test only(report["joining_disks"])["embedded_closed_disk"]
    @test M.representation_digest(d)==before
    @test length(H.prepare_marked_display(d;cut_curves=["cut"])["region_marks"][1]["cut_pieces"])==2
    @test H.lickorish_reference(3)["metadata"]["source_figure"]["asset_sha256"]=="c385ea4b13a46c869d40c25e96742f761d9b7121f1413931dac28e96e031dc73"
end
