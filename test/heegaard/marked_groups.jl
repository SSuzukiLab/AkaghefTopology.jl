using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
loadcase(name)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})
@testset "marked cellular groups and curve words" begin
    for g in 1:4
        d=H.lickorish_reference(g);p=H.surface_group_presentation(d)
        @test length(p.generators)==2g
        @test length(p.relators)==1
        @test all(iszero,H.abelianization_matrix(p))
        @test H.first_homology(p).free_rank==2g
        @test length(H.marked_curve_words(d).curve_words)==3g-1
    end
    for (surface,rank) in ((:disk,0),(:annulus,1))
        p=H.surface_group_presentation(B.grid_document(4;surface))
        @test isempty(p.relators)
        @test H.first_homology(p).free_rank==rank
    end
    @test H.first_homology(loadcase("sphere_empty")).free_rank==0
    @test H.first_homology(loadcase("torus_mer_lon")).free_rank==0
    parallel=loadcase("torus_parallel_zero_crossings")
    @test H.first_homology(parallel).free_rank==1
    @test isempty(H.first_homology(parallel).torsion)
    @test H.intersection_matrix(parallel)==zeros(Int,1,1)
    @test H.intersection_matrix(loadcase("torus_mer_lon"))==[1;;]
    # A repeated meridian is an extra relation, never an extra handle generator.
    B.add_curve!(parallel,"a2","alpha",true,[[(0,2),(4,2)]])
    @test H.first_homology(parallel).free_rank==1
    @test !H.compression_boundary_report(parallel).alpha_complete
end
