module PoincareJoiningDiskTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUTPUT=joinpath(ROOT,"artifacts","heegaard","full_curves")
include(joinpath(ROOT,"examples","HeegaardPoincareJoiningDisk.jl"))

@testset "Referenced Poincare diagram and D0 in the common atlas" begin
    ph=HeegaardPoincareJoiningDiskExample.fixture()
    h,d=ph.diagram,ph.joining_disk
    code=h.code.code
    labels=Dict(p.id=>p.source_label for p in code.intersections)
    @test [labels[i] for i in code.alpha[1].visits]==["x1","x2","x3","x4","x6","x5","x7"]
    @test [labels[i] for i in code.alpha[2].visits]==["y1","y5","y2","y3","y4"]
    @test [labels[i] for i in code.beta[1].visits]==["x1","x2","x3","x4","y1","x5","y2"]
    @test [labels[i] for i in code.beta[2].visits]==["x7","y3","y4","x6","y5"]
    @test H.surface_genus(h)==2
    @test length(H.boundary_cycles(h))==10
    @test sort(length.(H.boundary_cycles(h)))==[4,4,4,4,4,4,6,6,6,6]
    @test H.intersection_matrix(h)==[-3 2;2 -1]
    @test H.fundamental_group_presentation(h).relators==[
        [(1,-1),(1,-1),(1,-1),(1,-1),(2,1),(1,1),(2,1)],
        [(1,1),(2,-1),(2,-1),(1,1),(2,1)]]
    @test H.validate_joining_disk(h,d)
    @test d.region==5
    @test d.boundary_darts==[31,24,43,8]
    layout=H.load_layout(joinpath(OUTPUT,"PoincareD.hgd.layout.json"),h)
    before=H._canonical_json(H._json_dict(ph))
    chart=H.surface_chart(ph,layout)
    @test H.validate_layout(h,layout).valid
    @test length(chart.intersections)==12
    @test [p.label for p in chart.intersections]==HeegaardPoincareJoiningDiskExample.LABELS
    overlay=chart.joining_disk
    @test [(labels[b.from],b.family,labels[b.to]) for b in overlay.boundary]==[
        ("x4",:beta,"y1"),("y1",:alpha,"y4"),("y4",:beta,"x6"),("x6",:alpha,"x4")]
    @test Set((s.alpha_id,s.side) for s in overlay.sectors)==Set([(1,:minus),(2,:minus)])
    @test !overlay.base_exterior
    for k in eachindex(overlay.boundary)
        a=last(overlay.boundary[k].points)
        b=first(overlay.boundary[mod1(k+1,length(overlay.boundary))].points)
        @test hypot(a[1]-b[1],a[2]-b[2])<1e-8
    end
    @test H._canonical_json(H._json_dict(ph))==before
    payload=JSON3.read(read(joinpath(OUTPUT,"PoincareD.json"),String))
    @test payload.svg==H.render(chart)
    @test payload.base_svg==replace(read(joinpath(ROOT,"artifacts","heegaard","app","Poincare.svg"),String),"\r\n"=>"\n")
    @test count("data-intersection-id=",payload.svg)==12
    @test count("data-disk-piece=\"annulus\"",payload.svg)==2
    @test count("data-disk-piece=\"base\"",payload.svg)==1
    reloaded=H.load_punctured(joinpath(OUTPUT,"PoincareD.punctured.hgd.json"))
    @test H._canonical_json(H._json_dict(reloaded))==before
    # Existing abalone topology and D remain the accepted example, unchanged.
    abalone=JSON3.read(read(joinpath(OUTPUT,"Abalone.json"),String))
    @test abalone.diagram.joining_disk.boundary_darts==[1,8]
    @test abalone.visual_review=="owner-confirmed"
end
end
