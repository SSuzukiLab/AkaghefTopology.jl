using Test, Topology
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
@testset "ribbon compression vs independent marked cut" begin
    code=H.heegaard_code([[1,2]],[[1,2]],
        [H.Intersection(1,1,1,1),H.Intersection(2,1,1,-1)])
    # One annular complement and two disks. The alpha circle bounds a disk;
    # the beta circle is a nonseparating meridian on the torus.
    decorated=H.decorate(code,[[1,3],[2],[4]])
    @test H.surface_genus(decorated)==1
    ribbon=H.compression_boundary_report(decorated)
    @test sort(ribbon.alpha_boundary_genera)==[0,1]
    @test ribbon.beta_boundary_genera==[0]
    @test !ribbon.closed
    @test !H.validate(decorated).valid
    @test_throws H.ValidationError H.certify_heegaard(decorated)

    # A separately authored triangulated torus with the same embedded curves.
    d=B.grid_document(8)
    B.add_curve!(d,"a","alpha",true,[[(2,2),(6,2),(6,6),(2,6),(2,2)]])
    B.add_curve!(d,"b","beta",true,[[(4,0),(4,8)]])
    B.declare_crossing!(d,"lo",4,2,"a","b",1)
    B.declare_crossing!(d,"hi",4,6,"a","b",-1)
    marked=H.compression_boundary_report(d)
    @test sort(marked.alpha_boundary_genera)==sort(ribbon.alpha_boundary_genera)
    @test marked.beta_boundary_genera==ribbon.beta_boundary_genera

    # Established S3/lens models remain complete by actual cut reconstruction.
    for p in (1,2,5)
        h=H.decorate(H.heegaard_code([collect(1:p)],[collect(1:p)],
            [H.Intersection(i,1,1,1) for i in 1:p]);cellular=true)
        r=H.compression_boundary_report(h)
        @test r.closed && r.alpha_complete && r.beta_complete
        @test r.alpha_boundary_genera==r.beta_boundary_genera==[0]
    end
    empty=H.decorate(H.heegaard_code([Int[]],[Int[]],H.Intersection[]);cellular=true)
    @test_throws r"MARKED_INPUT_REQUIRED" H.surface_genus(empty)
    @test !H.validate(empty).valid
    @test_throws H.ValidationError H.certify_heegaard(empty)
end
