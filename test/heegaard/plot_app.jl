using Test, Topology, JSON3
const H=Topology.Heegaard

@testset "Joint plot objective and continuous constraints" begin
    paths=[[(0.,0.),(.1,.1),(.4,.3),(.8,.1),(1.,0.)],
           [(0.,1.),(.1,.9),(.4,.7),(.8,.9),(1.,1.)]]
    masses=H._plot_masses(paths)
    params=(spacing=2.2,spread=.36,tautness=4.,clearance=.14,bending=.04)
    e,f=H._plot_energy_gradient(paths,H.Puncture[],masses;params...)
    eps=1e-6
    for axis in 1:2
        a=deepcopy(paths);b=deepcopy(paths)
        p=paths[1][3]
        a[1][3]=(p[1]+(axis==1 ? eps : 0),p[2]+(axis==2 ? eps : 0))
        b[1][3]=(p[1]-(axis==1 ? eps : 0),p[2]-(axis==2 ? eps : 0))
        ea=first(H._plot_energy_gradient(a,H.Puncture[],masses;params...))
        eb=first(H._plot_energy_gradient(b,H.Puncture[],masses;params...))
        @test (ea-eb)/(2eps) ≈ -f[1][3][axis] atol=1e-6
    end
    disk=H.Puncture(1,1,:plus,(0.,0.),.2,1)
    @test !H._plot_motion_clear([[(-1.,1.),(1.,1.)]],[[(-1.,-1.),(1.,-1.)]],[disk])
    @test H._plot_motion_clear([[(-1.,1.),(1.,1.)]],[[(-1.,1.2),(1.,1.2)]],[disk])
    @test !H._plot_motion_clear([[(-1.,1.),(1.,1.)],[(-1.,0.),(1.,0.)]],
        [[(-1.,-1.),(1.,-1.)],[(-1.,0.),(1.,0.)]],H.Puncture[])
end

@testset "Actual saved plot boundary" begin
    root=normpath(joinpath(@__DIR__,"..",".."))
    output=joinpath(root,"artifacts","heegaard","app")
    entries=JSON3.read(read(joinpath(output,"manifest.json"),String))
    selection=JSON3.read(read(joinpath(root,"web","heegaard","gallery-selection.json"),String))
    selected=vcat(collect(selection.representatives),collect(selection.comparisons))
    ids=Set(String(e.id) for e in selected)
    @test length(ids)==length(selected)
    @test all(id->any(e->e.id==id,entries),ids)
    @test length(unique(e.id for e in entries))==length(entries)
    # Selection is deliberate. Increasing fixture count is not an acceptance gate.
    for entry in filter(e->e.id in ids,entries)
        diagram=H.certify_heegaard(H.load_heegaard(joinpath(output,"$(entry.id).hgd.json")))
        layout=H.load_layout(joinpath(output,"$(entry.id).hgd.layout.json"),diagram)
        @test layout.report.valid
        @test isempty(layout.crossings)
        @test length(layout.punctures)==2entry.genus
        @test length(unique(p.center[2] for p in layout.punctures))==2
        @test all(only(p.center[1] for p in layout.punctures if p.pair==i && p.side==:plus)==
                  only(p.center[1] for p in layout.punctures if p.pair==i && p.side==:minus)
                  for i in 1:entry.genus)
        @test H.render(layout;backend=:svg)==read(joinpath(output,"$(entry.id).svg"),String)
        @test repr(MIME"image/svg+xml"(),layout)==H.render(layout;backend=:svg)
        separation=layout.backend_trace["smoothing"]["minimum_flattened_separation"]
        @test length(layout.routes)==1 ? separation===nothing : separation>0
        @test H._exit_angle_degrees(layout)<1e-4
        @test layout.backend_trace["smoothing"]["minimum_normal_exit_length"]>0
        @test H.first_homology(diagram).torsion==Int.(entry.h1.torsion)
        for route in layout.routes
            guard=layout.backend_trace["smoothing"]["minimum_normal_exit_length"]
            @test H._vnorm(H._vsub(first(route.curve.segments).p3,first(route.curve.segments).p0))>=0.999guard
            @test H._vnorm(H._vsub(last(route.curve.segments).p3,last(route.curve.segments).p0))>=0.999guard
        end
    end
end
