using Test, TOML
include(joinpath(@__DIR__,"..","src","MarkedSurfaceDiagrams.jl"))
include(joinpath(@__DIR__,"..","src","builders.jl"))
using .MarkedSurfaceDiagrams, .MarkedDiagramBuilders
const M=MarkedSurfaceDiagrams
const ROOT=normpath(joinpath(@__DIR__,".."))
const EXPECT=TOML.parsefile(joinpath(@__DIR__,"expected.toml"))["cases"]
loadcase(n)=load_diagram(joinpath(ROOT,"examples",n*".toml"))
@testset "Exact shared fixtures — Julia against executed Python reference" begin
    for name in sort(collect(keys(EXPECT)))
        @testset "$name" begin
            e=EXPECT[name];s=loadcase(name);p=realize(s);rs=regions(s)
            @test representation_digest(s.doc)==e["representation_sha256"]
            @test marking_digest(s.doc)==e["marking_sha256"]
            @test [s.surface_report[k] for k in ("genus","boundary_components","chi","faces")]==e["surface"]
            @test sort([join([r["genus"],r["boundary_components"],r["chi"],Int(r["normalized_disk"]),Int(r["embedded_closed_disk"])],"/") for r in rs])==e["regions"]
            coords=[f*"/"*v*"="*join(M.qstr.(collect(p[f][v])),",") for f in sort(collect(keys(p))) for v in sort(collect(keys(p[f])))]
            @test coords==e["coordinates"]
            cycles=intersection_cycles(s)
            @test [c*":"*join([x["id"]*(x["sign"]==1 ? "+" : "-") for x in cycles[c]],",") for c in sort(collect(keys(cycles)))]==e["cycles"]
            @test [c*":"*join([join(x["side"],"/")*"@"*string(x["gate_rank"])*(x["direction"]==1 ? "+" : "-") for x in seam_word(s,c)],",") for c in sort(collect(keys(s.curves)))]==e["seams"]
            @test family_cut_report(s,"alpha")["cut_system"]==e["alpha_cut"]
            @test family_cut_report(s,"beta")["cut_system"]==e["beta_cut"]
            @test certify_layout(s,p)
            for (_,(_,ds)) in s.curves;@test parse_word(print_word(ds))==ds;end
            mktempdir() do dir
                result=write_svg(s,p,joinpath(dir,"diagram.svg"))
                @test result["certified_centerlines"]
                @test occursin("centerlines_certified=true",read(joinpath(dir,"diagram.svg"),String))
            end
        end
    end
end
@testset "Winding and true fixed boundary anchors" begin
    cases=[loadcase("annulus_winding_"*s) for s in ("m1","p0","p1")]
    @test [sum(x["direction"] for x in seam_word(s,"gamma");init=0) for s in cases]==[-1,0,1]
    @test endpoint_anchors(cases[1],"gamma")==endpoint_anchors(cases[2],"gamma")==endpoint_anchors(cases[3],"gamma")
end
@testset "Koda convention and immutable semantic snapshot" begin
    s=loadcase("torus_three_crossings");old=deepcopy(s.doc)
    @test !any(check_joining_disk(s,r["seed"])["koda_joining"] for r in regions(s))
    @test count(r->check_joining_disk(s,r["seed"];orientations=Dict("alpha1"=>-1))["koda_joining"],regions(s))==1
    @test s.doc==old
    d=deepcopy(s.doc);d["region_marks"]=[Dict("id"=>"D0","seed"=>["P",68],"claim"=>"embedded_closed_disk")];validate(d)
    d["region_marks"][1]["claim"]="koda_joining";@test_throws DiagramError validate(d)
    s.doc["curves"][1]["word"]="P{a b}";@test_throws DiagramError realize(s)
end
@testset "Fail closed" begin
    d=deepcopy(loadcase("torus_mer_lon").doc);d["crossings"][1]["sign"]=-1;@test_throws DiagramError validate(d)
    d=deepcopy(loadcase("torus_mer_lon").doc);d["crossings"]=[];@test_throws DiagramError validate(d)
    d=deepcopy(loadcase("torus_mer_lon").doc);d["surface"]["genus"]=2;@test_throws DiagramError validate(d)
    d=deepcopy(loadcase("torus_mer_lon").doc);d["curves"][2]["role"]="alpha";@test_throws DiagramError validate(d)
    d=deepcopy(loadcase("torus_mer_lon").doc);d["isotopy_canonical"]=true;@test_throws DiagramError validate(d)
    s=loadcase("torus_mer_lon");p=realize(s)
    @test_throws DiagramError move_vertex(s,p,"P","v0_2",["0","0"])
    old=p["P"]["v2_2"];q,w=move_vertex(s,p,"P","v2_2",[M.qstr(old[1]+big(1)//big(10000)),M.qstr(old[2])])
    @test certify_layout(s,q);@test w["same_representation"]==representation_digest(s.doc)
    for k in 0:20
        t=big(k)//big(20);mid=deepcopy(p);mid["P"]["v2_2"]=Tuple((1-t)*old[j]+t*q["P"]["v2_2"][j] for j in 1:2);@test certify_layout(s,mid)
    end
end
@testset "Symbolic authoring builder" begin
    d=grid_document(4;surface=:annulus)
    add_curve!(d,"gamma","aux",false,[[(1,0),(4,3)],[(0,3),(1,4)]])
    s=validate(d);t=loadcase("annulus_winding_p1")
    @test representation_digest(s.doc)==representation_digest(t.doc)
end
