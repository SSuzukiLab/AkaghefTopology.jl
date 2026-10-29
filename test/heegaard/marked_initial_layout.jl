module MarkedInitialLayoutTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const ROOT=normpath(joinpath(@__DIR__,"../.."))
include("marked_fixtures.jl")

# Source: test/heegaard/marked_topology.jl's auxiliary-cut D regression.
# This contractible alpha circle is a rendering boundary case, not a complete
# genus-one Heegaard splitting or evidence for arbitrary numerical inputs.
function split_disk()
    d=B.grid_document(8)
    B.add_curve!(d,"a","alpha",true,[[(2,2),(6,2),(6,6),(2,6),(2,2)]])
    B.add_curve!(d,"cut","aux",true,[[(4,0),(4,8)]])
    B.declare_crossing!(d,"p",4,2,"a","cut",1)
    B.declare_crossing!(d,"q",4,6,"a","cut",-1)
    d["region_marks"]=[Dict("id"=>"D","seed"=>["P",54],"claim"=>"embedded_closed_disk")]
    d
end

# Source: examples/HeegaardInput.jl's torus pair; the parallel shifted cut is
# independent auxiliary data and does not replace either mathematical curve.
function crossing_torus()
    d=B.grid_document(6)
    B.add_curve!(d,"a","alpha",true,[[(0,2),(6,2)]])
    B.add_curve!(d,"b","beta",true,[[(3,0),(3,6)]])
    B.add_curve!(d,"cut","aux",true,[[(0,4),(6,4)]])
    B.declare_crossing!(d,"x",3,2,"a","b",1)
    B.declare_crossing!(d,"port",3,4,"cut","b",1)
    d
end

@testset "whole marked source display preparation" begin
    d=split_disk();before=deepcopy(d)
    p=H.prepare_marked_display(d;cut_curves=["cut"])
    @test d==before
    pieces=only(p["region_marks"])["cut_pieces"]
    @test length(pieces)==2
    @test all(any(x["kind"]=="cut" for cycle in piece["display_boundary_cycles"] for x in cycle) for piece in pieces)
    @test Set(Tuple(t) for piece in pieces for t in piece["triangles"])==Set(Tuple(t) for t in only(p["region_marks"])["source_region"]["triangles"])
    @test_throws H.MarkedSurfaceDiagrams.DiagramError H.prepare_marked_display(d;cut_curves=["cut"],visible_curves=String[])
    crossed=H.prepare_marked_display(crossing_torus();cut_curves=["cut"])
    @test length(crossed["crossings"])==1 && crossed["crossings"][1]["id"]=="x" && crossed["crossings"][1]["sign"]==1
    @test_throws ArgumentError H.layout_marked(d;cut_curves=["cut"])
    @test_throws ArgumentError H.layout_marked(d;cut_curves=["cut"],outer_triangle=["P",true])
    @test_throws ArgumentError H.layout_marked(d;cut_curves=["cut"],outer_triangle=["missing",0])
    sphere=JSON3.read(read(joinpath(ROOT,"vendor/marked_surface_diagrams/examples/sphere_empty.json"),String),Dict{String,Any})
    p0=H.prepare_marked_display(sphere;cut_curves=String[])
    @test p0["genus"]==0 && isempty(p0["boundaries"]) && isempty(p0["curves"])
    nested=JSON3.read(read(joinpath(ROOT,"vendor/marked_surface_diagrams/examples/nested_zero_crossings.json"),String),Dict{String,Any})
    @test H.MarkedSurfaceDiagrams.validate(nested).surface_report["boundary_components"]==1
    rejected=try H.prepare_marked_display(nested;cut_curves=String[]) catch e; e end
    @test rejected isa H.MarkedSurfaceDiagrams.DiagramError && rejected.code=="DISPLAY_SURFACE"
    # Source: marked_fixtures.jl::capped_nested_sphere. The old disk remains a
    # negative case; the new closed surface keeps that source patch and words.
    closed_nested=capped_nested_sphere()
    report=H.MarkedSurfaceDiagrams.validate(closed_nested).surface_report
    @test report["chi"]==2 && report["genus"]==0 && report["boundary_components"]==0
    @test closed_nested["patches"][1]==nested["patches"][1] && closed_nested["curves"]==nested["curves"]
    nested_prepared=H.prepare_marked_display(closed_nested;cut_curves=String[])
    @test nested_prepared["genus"]==0 && length(nested_prepared["curves"])==2
    # Explicit rendering choice wins, while retained native source stays exact.
    sphere["metadata"]=Dict("display"=>Dict("cut_curves"=>String[],
        "infinity_triangle"=>p0["triangles"][end]["source"]))
    view=H.plot(sphere;cut_curves=String[],outer_triangle=p0["triangles"][1]["source"])
    @test view["source_document"]==sphere
    @test view["display_request"]["infinity_triangle"]==p0["triangles"][1]["source"]
    @test view["infinity_source_triangle"]!=sphere["metadata"]["display"]["infinity_triangle"]
    @test view["smooth_check"]["empty_visible_family"]
    @test H.marked_display_source(view)==sphere
    @test occursin("<svg",H.render(view))
    bad=deepcopy(view);bad["infinity_source_triangle"]=["missing",0]
    @test_throws ArgumentError H.marked_display_source(bad)
    bad=deepcopy(view);bad["disk"]=[0,1,2]
    @test_throws ArgumentError H.marked_display_source(bad)
    bad=deepcopy(view);bad["marking_id"]="another marking"
    @test_throws ArgumentError H.marked_display_source(bad)
    physical=B.grid_document(4;surface=:disk)
    @test_throws H.MarkedSurfaceDiagrams.DiagramError H.prepare_marked_display(physical;cut_curves=String[])
end

# Optional emitted inputs let the numerical backend be tested independently of
# Julia's compiler slot. No fixture selects a replacement path or exterior.
if haskey(ENV,"HEEGAARD_LAYOUT_FIXTURES")
    directory=ENV["HEEGAARD_LAYOUT_FIXTURES"];mkpath(directory)
    cases=Any[("split",split_disk(),["cut"],nothing),("crossing",crossing_torus(),["cut"],nothing)]
    push!(cases,("sphere_nested",capped_nested_sphere(),String[],nothing))
    for g in 1:4
        d=H.lickorish_reference(g)
        push!(cases,("reference_g$g",d,["l$i" for i in 1:g],[c["id"] for c in d["curves"] if !(c["id"] in ["l$i" for i in 1:g])]))
    end
    for (name,d,cuts,visible) in cases
        p=H.prepare_marked_display(d;cut_curves=cuts,visible_curves=visible)
        outer=name in ("split","crossing") ? ["P",0] : p["triangles"][1]["source"]
        p["source_document"]=d;p["outer_triangle"]=outer
        p["display_request"]=Dict("cut_curves"=>p["cut_curves"],"visible_curves"=>p["visible_curves"],"infinity_triangle"=>outer)
        write(joinpath(directory,name*".prepared.json"),JSON3.write(p))
        if name=="split"
            p["outer_triangle"]=only(p["region_marks"])["source_region"]["triangles"][1]
            p["display_request"]["infinity_triangle"]=p["outer_triangle"]
            write(joinpath(directory,"split_infinity.prepared.json"),JSON3.write(p))
        end
    end
end
end
