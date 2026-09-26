module JoiningDiskSurfaceTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT=joinpath(ROOT,"artifacts","heegaard","full_curves")

@testset "Joining disk on the common surface atlas" begin
    ph=H.load_punctured(joinpath(OUTPUT,"abalone","S3.punctured.hgd.json"))
    h,d=ph.diagram,ph.joining_disk
    layout=H.load_layout(joinpath(OUTPUT,"Abalone.hgd.layout.json"),h)
    source=JSON3.write(H._json_dict(ph)); geometry=JSON3.write(H._layout_json(layout))
    chart=H.surface_chart(ph,layout); overlay=chart.joining_disk
    @test H.render(chart)==H.render(H.surface_chart(h,layout;joining_disk=d))
    @test overlay.disk.boundary_darts==[1,8]
    @test [(b.family,b.from,b.to) for b in overlay.boundary]==[(:alpha,1,2),(:beta,2,1)]
    @test only(overlay.sectors).side==:plus
    @test !overlay.base_exterior
    route=only(r for r in layout.routes if r.curve_id==1 && r.segment_id==1)
    beta=only(b for b in overlay.boundary if b.family==:beta)
    @test beta.points==reverse(route.points)
    last_segment=last(route.curve.segments)
    @test startswith(beta.tail,"C $(H._svg_point(last_segment.p2)) $(H._svg_point(last_segment.p1)) $(H._svg_point(last_segment.p0))")
    for k in eachindex(overlay.boundary)
        a=last(overlay.boundary[k].points)
        b=first(overlay.boundary[mod1(k+1,length(overlay.boundary))].points)
        @test hypot(a[1]-b[1],a[2]-b[2])<1e-8
    end
    @test JSON3.write(H._json_dict(ph))==source
    @test JSON3.write(H._layout_json(layout))==geometry
    @test H.surface_chart(h,layout).joining_disk===nothing
    @test_throws ArgumentError H.surface_chart(h,layout;joining_disk=H.JoiningDisk(d.labelled_diagram_hash,d.region,[1,7]))
    @test_throws ArgumentError H.plot(ph,layout;view=:alpha_cut)
    # The other valid bigon exercises the minus attachment and inner annulus.
    alternate=H.surface_chart(h,layout;joining_disk=H.joining_disk(h,3)).joining_disk
    @test only(alternate.sectors).side==:minus
    @test !alternate.base_exterior
    @test Set(b.from for b in alternate.boundary)==Set([2,3])
    # Use the real common plot entrance, not just a hand-built chart.
    seed=H.load_layout(joinpath(OUTPUT,"Abalone-before.hgd.layout.json"),h)
    plotted=H.plot(ph,seed;view=:surface)
    @test plotted isa H.HeegaardSurfaceChart
    @test H.validate_layout(h,plotted.base).valid
    @test plotted.joining_disk.disk.boundary_darts==d.boundary_darts
    @test occursin("data-joining-disk=\"1\"",H.render(plotted))
end
end
