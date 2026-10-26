module MarkedStabilizationTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const B=H.MarkedDiagramBuilders
roundjson(x)=JSON3.read(JSON3.write(x),Dict{String,Any})
fixture(n)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",n*".json"),String),Dict{String,Any})

# Source: marked_fixtures.jl::capped_nested_sphere's opposite oriented cap
# gluing. A smaller authored disk retains a protected triangular D and one
# genuinely constructed old zero-intersection component, not a post-move answer.
function disk_marked_sphere()
    d=B.grid_document(3;surface=:disk)
    B.add_curve!(d,"old_alpha","alpha",true,[[(1,1),(2,2),(2,1),(1,1)]])
    p=only(d["patches"]);bd=[v for s in p["side_order"] for v in p["sides"][s][1:end-1]]
    push!(d["patches"],Dict("face"=>"Q","side_order"=>deepcopy(p["side_order"]),"sides"=>deepcopy(p["sides"]),
        "triangles"=>[["cap_center",bd[i],bd[mod1(i+1,length(bd))]] for i in eachindex(bd)]))
    d["surface"]=Dict("marking_id"=>"protected-triangle-sphere","genus"=>0,"boundary"=>Any[],
        "gluings"=>[[["P","south"],["Q","south"]],[["P","east"],["Q","west"]],
                    [["P","north"],["Q","north"]],[["P","west"],["Q","east"]]])
    d["region_marks"]=[Dict("id"=>"D","seed"=>["P",8],"claim"=>"embedded_closed_disk")]
    d["metadata"]=Dict("provenance"=>Dict("preserved"=>true,"decimal"=>0.125),"display"=>Dict("infinity_triangle"=>["P",0]))
    d
end

@testset "Actual sphere stabilization and generated cut-cap inverse" begin
    source=fixture("sphere_empty")
    source["metadata"]=Dict("source"=>"native two-triangle sphere","display"=>Dict("infinity_triangle"=>["Q",0]))
    original=deepcopy(source);result=H.marked_stabilize(source,Dict("support"=>["P",0]))
    target=result.document
    @test source==original
    @test H.validate(target;level=:closed_heegaard).valid
    @test H.surface_genus(target)==1
    compression=H.compression_boundary_report(target)
    @test compression.alpha_complete && compression.beta_complete
    @test length(target["crossings"])==1 && only(target["crossings"])["sign"]==1
    @test [c["role"] for c in target["curves"]]==["alpha","beta"]
    @test length(target["patches"])<=64
    @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in target["patches"])<=256
    neck=result.witness["correspondence"]["neck"]
    @test neck["torus"]["genus"]==1 && neck["torus"]["boundary_components"]==1
    @test neck["outside"]["genus"]==0 && neck["outside"]["boundary_components"]==1
    @test target["metadata"]["source"]==source["metadata"]["source"]
    h=H.first_homology(target)
    @test h.free_rank==0 && isempty(h.torsion) # standard meridian-longitude S3 splitting
    inverse=roundjson(result.witness["inverse_specification"])
    @test H.marked_destabilize(target,inverse).document==source
    @test H.replay_marked_move(source,roundjson(result.witness)).document==target
    @test H.replay_marked_move(target,roundjson(H.inverse_marked_move(result.witness))).document==source
    forged=deepcopy(inverse)
    forged["neck_word"]=only(c for c in target["curves"] if c["role"]=="alpha")["word"]
    @test_throws ArgumentError H.marked_destabilize(target,forged)
    altered=deepcopy(target);altered["metadata"]["source"]="changed"
    @test_throws ArgumentError H.marked_destabilize(altered,inverse)
    @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>["Q",0]))
    for malformed in (Any["P",true],Any["P",0.0],Any["P",0,"extra"],Any["P",-1],Any[1,0])
        @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>malformed))
    end
    for key in ("alpha_id","beta_id","crossing_id")
        @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>["P",0],key=>false))
    end
    numeric_curve=fixture("torus_parallel_zero_crossings");numeric_curve["curves"][1]["id"]=7
    @test_throws ArgumentError H.marked_stabilize(numeric_curve,Dict("support"=>["P",0]))
    numeric_crossing=fixture("torus_mer_lon");numeric_crossing["crossings"][1]["id"]=7
    @test_throws ArgumentError H.marked_stabilize(numeric_crossing,Dict("support"=>["P",0]))
    numeric_mark=disk_marked_sphere();numeric_mark["region_marks"][1]["id"]=7
    @test_throws ArgumentError H.marked_stabilize(numeric_mark,Dict("support"=>["Q",0]))
end

@testset "Genus-one source curves and crossings survive connected sum" begin
    source=fixture("torus_mer_lon");source["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",31]))
    result=H.marked_stabilize(source,Dict("support"=>["P",0],"alpha_id"=>"extra_alpha","beta_id"=>"extra_beta","crossing_id"=>"extra_crossing"))
    @test H.surface_genus(result.document)==2
    c=H.compression_boundary_report(result.document)
    @test c.alpha_complete && c.beta_complete
    @test length(result.document["curves"])==4 && length(result.document["crossings"])==2
    @test result.document["crossings"][1]["id"]==source["crossings"][1]["id"]
    @test result.document["crossings"][1]["sign"]==source["crossings"][1]["sign"]
    @test H.marked_destabilize(result.document,roundjson(result.witness["inverse_specification"])).document==source
    @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>["P",10])) # occupied alpha/beta vertices
end

@testset "Protected D and infinity are source bound; joining claim is never downgraded" begin
    source=disk_marked_sphere();original=deepcopy(source)
    @test H.validate(source;level=:surface).valid
    result=H.marked_stabilize(source,Dict("support"=>["Q",0]))
    @test source==original
    @test only(result.document["region_marks"])["claim"]=="embedded_closed_disk"
    @test only(H.marked_topology_report(result.document)["joining_disks"])["embedded_closed_disk"]
    @test result.document["metadata"]["provenance"]==source["metadata"]["provenance"]
    @test length(result.witness["correspondence"]["triangles"])==sum(length(p["triangles"]) for p in source["patches"])-1
    @test H.marked_destabilize(result.document,roundjson(result.witness["inverse_specification"])).document==source
    @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>["P",8]))
    joining=deepcopy(source);joining["region_marks"][1]["claim"]="koda_joining"
    @test H.validate(joining;level=:surface).valid
    failure=try H.marked_stabilize(joining,Dict("support"=>["Q",0]));nothing catch err;err end
    @test failure isa ArgumentError && occursin("unchanged Koda joining disk",sprint(showerror,failure))
    @test joining["region_marks"][1]["claim"]=="koda_joining"
    @test_throws ArgumentError H.marked_stabilize(B.grid_document(3;surface=:disk),Dict("support"=>["P",0]))
end

@testset "Actual triangle-atlas compaction retains inverse cell correspondence" begin
    source=B.grid_document(5)
    B.add_curve!(source,"old_alpha","alpha",true,[[(0,1),(5,1)]])
    B.add_curve!(source,"old_beta","beta",true,[[(1,0),(1,5)]])
    B.declare_crossing!(source,"old_crossing",1,1,"old_alpha","old_beta",1)
    source["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    # 50 source faces +16 inserted faces =66 temporary disk charts: this
    # actually exercises the accepted merger and its absorbed-cell offsets.
    result=H.marked_stabilize(source,Dict("support"=>["P",48]))
    @test length(result.witness["inverse_specification"]["merges"])==2
    @test H.surface_genus(result.document)==2
    @test H.compression_boundary_report(result.document).closed
    @test H.marked_destabilize(result.document,roundjson(result.witness["inverse_specification"])).document==source
end

@testset "Explicit independent auxiliary display-cut extension" begin
    # Independent source authoring: the existing alpha is never reclassified
    # as aux. (beta,old_cut) has the independently known sign (north,east)=-1.
    source=B.grid_document(4)
    B.add_curve!(source,"old_alpha","alpha",true,[[(0,1),(4,1)]])
    B.add_curve!(source,"old_beta","beta",true,[[(1,0),(1,4)]])
    B.add_curve!(source,"old_cut","aux",true,[[(0,2),(4,2)]])
    B.declare_crossing!(source,"old_ab",1,1,"old_alpha","old_beta",1)
    B.declare_crossing!(source,"old_bc",1,2,"old_beta","old_cut",-1)
    source["metadata"]=Dict("display"=>Dict("cut_curves"=>["old_cut"],"infinity_triangle"=>["P",0]))
    snapshot=deepcopy(source)
    spec=Dict("support"=>["P",28],"extend_display_cuts"=>true,
        "auxiliary_id"=>"extra_cut","auxiliary_crossing_id"=>"extra_bc")
    result=H.marked_stabilize(source,spec);target=result.document
    @test source==snapshot
    @test H.surface_genus(target)==2
    @test target["metadata"]["display"]["cut_curves"]==["old_cut","extra_cut"]
    @test only(c for c in target["curves"] if c["id"]=="stabilization_alpha")["role"]=="alpha"
    @test only(c for c in target["curves"] if c["id"]=="extra_cut")["role"]=="aux"
    extra=only(x for x in target["crossings"] if x["id"]=="extra_bc")
    @test extra["curves"]==["stabilization_beta","extra_cut"] && extra["sign"]==-1
    report=result.witness["correspondence"]["target_display_cut"]
    @test report["cut_genus"]==0 && report["boundary_components"]==4 && report["components"]==1
    @test length(result.witness["correspondence"]["generated_geometric_regions"])==1
    @test result.witness["correspondence"]["infinity"]["source_seed"]==source["metadata"]["display"]["infinity_triangle"]
    image=only(x for x in result.witness["correspondence"]["triangles"] if x["source"]==["P",0])
    @test target["metadata"]["display"]["infinity_triangle"]==first(image["target"])
    prepared=H.prepare_marked_display(target;cut_curves=["old_cut","extra_cut"])
    @test prepared["cut_topology"]==Dict("genus"=>0,"boundary_components"=>4,"components"=>1)
    @test H.marked_destabilize(target,roundjson(result.witness["inverse_specification"])).document==source
    # D continues to mean alpha union beta, even with the new auxiliary cut.
    marked=disk_marked_sphere();marked["metadata"]["display"]["cut_curves"]=String[]
    dresult=H.marked_stabilize(marked,Dict("support"=>["Q",0],"extend_display_cuts"=>true))
    @test only(H.marked_topology_report(dresult.document)["joining_disks"])["embedded_closed_disk"]
    @test H.marked_destabilize(dresult.document,roundjson(dresult.witness["inverse_specification"])).document==marked
    nocuts=fixture("sphere_empty");nocuts["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["Q",0]))
    @test_throws ArgumentError H.marked_stabilize(nocuts,Dict("support"=>["P",0],"extend_display_cuts"=>true))
    noinf=deepcopy(source);delete!(noinf["metadata"]["display"],"infinity_triangle")
    @test_throws ArgumentError H.marked_stabilize(noinf,spec)
    wrong=deepcopy(source);wrong["metadata"]["display"]["cut_curves"]=["old_alpha"]
    @test_throws ArgumentError H.marked_stabilize(wrong,spec)
    @test_throws ArgumentError H.marked_stabilize(source,Dict("support"=>["P",28],"extend_display_cuts"=>1))
end
end
