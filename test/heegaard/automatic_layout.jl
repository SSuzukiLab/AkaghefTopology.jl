module AutomaticLayoutTests
using Test, Topology
const H=Topology.Heegaard
@testset "Fixture-free layout and disconnected genus three" begin
    code=H.heegaard_code([[1],[2],[3]],[[1],[2],[3]],
        [H.Intersection(i,i,i,1) for i in 1:3])
    h=H.certify_heegaard(H.decorate(code,[[1,2,3]]))
    source=H._labelled_hash(h)
    layout=H.layout_heegaard(h;backend=:pure_julia,router=:planar)
    @test H.surface_genus(h)==3
    @test layout.report.valid
    @test isempty(layout.crossings)
    @test layout.backend_trace["terminal_embedding"]["components"]==3
    @test layout.backend_trace["fixture_used"]===false
    @test [(p.center,p.radius) for p in layout.punctures]==
        [(p.center,p.radius) for p in H.default_puncture_placement(H.cut_along_alpha(h))]
    @test H._labelled_hash(h)==source
    @test all(hypot((a.point.-b.point)...)<1e-8 for (a,b) in
        zip(layout.ports,H.place_ports(H.cut_along_alpha(h),layout.punctures)))
    @test_throws ArgumentError H.layout_heegaard(h;router=:planar,placement=layout.punctures)
    @test_throws ArgumentError H.layout_heegaard(h;router=:planar,route_certificate="fixture.json")
    chart=H.plot(h;view=:surface)
    @test chart.genus==3
    @test length(chart.annuli)==3
    @test chart.base.backend_trace["fixture_used"]===false
end
end
