module OGraphReverseTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
include(joinpath(@__DIR__,"../../examples/HeegaardAbalone.jl"))

@testset "Independent o-data to marked Heegaard construction" begin
    # Independent registered S3 o-graph; no Heegaard source is provided.
    g=H.OData([[1,-1]],[1],Dict{String,Any}())
    d=H.to_heegaard(g)
    @test H.surface_genus(d)==2 # standard dual-graph splitting, not source genus
    c=H.compression_boundary_report(d)
    @test c.alpha_complete && c.beta_complete
    @test H.validate(d;level=:closed_heegaard).valid
    witness=d["metadata"]["correspondence"]
    @test length(witness["flag_quadrilaterals"])==24
    @test length(witness["tetrahedron_faces"])==2
    @test length(witness["overlay_face_darts"])==4
    @test all(x["chi"]==2 for x in witness["tetrahedron_vertex_links"])
    @test length(witness["crossings"])==6
    @test all(!d["metadata"][k] for k in ("retained_source_used","source_marking_recovered","joining_disk_recovered","branching_field_recovered"))
    poisoned=H.OData(g.gauss,g.orientation,Dict{String,Any}("marked_source"=>Dict("invalid"=>true),"source_hash"=>"false source"))
    @test H.to_heegaard(poisoned)==d
    @test g.metadata==Dict{String,Any}()
    # Independent finite manifold evidence: abalone is the registered S3 graph.
    h=H.first_homology(d)
    @test h.free_rank==0 && isempty(h.torsion)
    # A forward input's source can be removed before invoking the reverse.
    ex=HeegaardAbaloneExample.fixture();b=H.to_odata(ex.diagram;joining_disk=ex.disk)
    raw=H.OData(deepcopy(b.ograph.gauss),copy(b.ograph.orientation),Dict{String,Any}())
    @test H.to_heegaard(raw)==d
    @test H.to_heegaard(b.spine;strategy=:retained_source)===ex.diagram
    @test H.to_heegaard(b.spine;strategy=:dual_tetrahedra)==d
    @test_throws ArgumentError H.to_heegaard(g;strategy=:retained_source)
    @test_throws Exception H.to_heegaard(H.OData([[1,-1]],[0],Dict{String,Any}()))
    @test_throws ArgumentError H.to_heegaard(H.OData(Vector{Int}[],Int[],Dict{String,Any}()))
    # Fixed mixed-sign supplied input, not a result of this reverse converter.
    source=JSON3.read(read(joinpath(@__DIR__,"../../vendor/koda_genus3/data/figure6.odata.json"),String))
    kg=H.OData([reverse(Int.(w)) for w in source[1]],Int.(source[2]),Dict{String,Any}())
    kd=H.to_heegaard(kg);kc=H.compression_boundary_report(kd)
    @test kc.alpha_complete && kc.beta_complete
    @test H.surface_genus(kd)==6
    @test length(kd["metadata"]["correspondence"]["flag_quadrilaterals"])==120
    kh=H.first_homology(kd)
    @test kh.free_rank==0 && kh.torsion==[2,2]
    # Public marked source has no source projection metadata and survives JSON.
    @test M.validate(JSON3.read(JSON3.write(kd),Dict{String,Any})).surface_report["genus"]==6
    # Independent hand transcript already recorded in docs/heegaard/joining-disk.md.
    # Eight vertices exercise native chart compaction (initially 80 charts).
    pg=H.OData([[-1,-2,-3,-4,-5,6,4,7,-7,-8,-6,5,8,3,2,1]],
        [-1,-1,-1,1,1,1,-1,1],Dict{String,Any}())
    pd=H.to_heegaard(pg);pc=H.compression_boundary_report(pd)
    @test pc.alpha_complete && pc.beta_complete
    @test H.surface_genus(pd)==9
    @test length(pd["patches"])<=64
    @test !isempty(pd["metadata"]["correspondence"]["chart_merges"])
    @test length(pd["metadata"]["correspondence"]["overlay_face_tetrahedron_vertices"])==32
    for (input,base) in ((g,d),(kg,kd))
        display=H.to_heegaard(input;display_cuts=true)
        cuts=display["metadata"]["display_cuts"]
        @test display["metadata"]["display"]["cut_curves"]==cuts
        @test length(cuts)==H.surface_genus(base)
        @test all(c["role"]=="aux" for c in display["curves"] if c["id"] in cuts)
        @test M.intersection_cycles(H._marked_states(display).mathematics)==M.intersection_cycles(H._marked_states(base).mathematics)
        @test H.compression_boundary_report(display).closed
        prep=H.prepare_marked_display(display;cut_curves=cuts)
        @test length(prep["boundaries"])==2length(cuts)
        infinity=display["metadata"]["display"]["infinity_triangle"]
        @test any(t["source"]==infinity for t in prep["triangles"])
        extension=display["metadata"]["display"]["auxiliary_construction"]
        @test length(extension["rectangles"])==length(extension["crossing_copies"])
        @test extension["original_alpha_beta_preserved"]
    end
end
end
