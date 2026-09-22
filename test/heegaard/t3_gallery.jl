using Test, Topology
const H=Topology.Heegaard
root=normpath(joinpath(@__DIR__,"..",".."))
@testset "Source T3, not a stabilized lens diagram" begin
    diagram=H.certify_heegaard(H.load_heegaard(joinpath(root,"test","fixtures","t3_octahedron_v1.hgd.json")))
    @test H.surface_genus(diagram)==3
    @test length(diagram.code.regions.boundary_cycle_blocks)==8
    homology=H.first_homology(diagram)
    @test homology.free_rank==3 && isempty(homology.torsion)
    geometric=zeros(Int,3,3)
    for p in diagram.code.code.intersections; geometric[p.alpha,p.beta]+=1;end
    @test geometric==[2 0 2;2 2 0;0 2 2]
    @test H.canonical_hash(diagram)=="2c5b90a5e4c1fbc6673de8a0966799c5e28b30bcb0068f2fdfd18c61d5058ff7"
    layout=H.load_layout(joinpath(root,"artifacts","heegaard","app","T3.hgd.layout.json"),diagram)
    @test layout.report.valid && isempty(layout.crossings)
    @test [p.center for p in layout.punctures]==[(1.,1.),(1.,-1.),(2.,1.),(2.,-1.),(3.,1.),(3.,-1.)]
    @test H._exit_angle_degrees(layout)<1e-4
    @test repr(MIME"image/svg+xml"(),layout)==read(joinpath(root,"artifacts","heegaard","app","T3.svg"),String)
    # Existing saved fixture identity must survive the canonicalization speedup.
    poincare=H.certify_heegaard(H.load_heegaard(joinpath(root,"test","fixtures","poincare_homology_sphere_v1.hgd.json")))
    @test H.load_layout(joinpath(root,"test","fixtures","poincare_homology_sphere_v1.hgd.layout.json"),poincare).report.valid
end
