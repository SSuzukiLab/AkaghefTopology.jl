module MarkedTwistTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const M=H.MarkedSurfaceDiagrams

function twist_torus()
    doc=B.grid_document(2)
    B.add_curve!(doc,"a","alpha",true,[[(0,1),(2,1)]])
    B.add_curve!(doc,"b","beta",true,[[(1,0),(1,2)]])
    B.declare_crossing!(doc,"x",1,1,"a","b",1)
    doc["metadata"]=Dict("source"=>"oriented meridian/longitude test","display"=>Dict("infinity_triangle"=>["P",7]))
    doc
end

# Source: exact whole-annulus half-shear convention. This independently checks
# [b] -> [b]+[a] in fixed PRE-TWIST PREPARED fine-mesh generators, not a source hash or a
# basis reconstructed anew from the output triangulation.
@testset "Whole-diagram fixed-marking Dehn twist and same-product inverse" begin
    doc=twist_torus();before=deepcopy(doc)
    println("Positive whole-diagram twist started");flush(stdout)
    positive_started=time()
    result=H.marked_dehn_twist(doc,Dict("core"=>"a","sign"=>1))
    println("Positive twist constructor seconds: ",time()-positive_started);flush(stdout)
    @test doc==before
    @test H.validate(result.document;level=:surface).valid
    # Source: the original 991-vertex failure on THIS unchanged input. The
    # repair must really rechart and satisfy the unchanged native limits.
    correspondence=result.witness["correspondence"]
    @test correspondence["rechart"]!==nothing
    @test length(result.document["patches"])<=64
    @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in result.document["patches"])<=256
    prepared=result.witness["inverse_specification"]["fine_source_mesh"]
    old_faces=Set(p["face"] for p in prepared["patches"]);native_faces=Set(p["face"] for p in result.document["patches"])
    @test all(r["source"][1] in old_faces && r["target"][1] in old_faces &&
        r["child"][1] in native_faces for r in correspondence["child_cells"])
    @test H.is_closed_heegaard(result.document)
    @test H.intersection_matrix(result.document)==[1;;]
    @test length(result.document["crossings"])==1
    @test result.document["crossings"][1]["id"]=="x" && result.document["crossings"][1]["sign"]==1
    action=result.witness["correspondence"]["old_basis_action"]
    a=only(r for r in action if r["curve"]=="a");b=only(r for r in action if r["curve"]=="b")
    @test a["image_exponents"]==a["source_exponents"]
    @test b["intersection_with_core"]==1
    @test b["image_exponents"]==b["source_exponents"]+a["source_exponents"]
    @test b["image_exponents"]!=b["source_exponents"]
    inverse=JSON3.read(JSON3.write(result.witness["inverse_specification"]),Dict{String,Any})
    println("Positive twist generated inverse started");flush(stdout)
    inverse_started=time()
    restored=H.marked_inverse_dehn_twist(result.document,inverse)
    println("Positive twist inverse seconds: ",time()-inverse_started);flush(stdout)
    @test restored.document==doc
    @test restored.witness["correspondence"]["native_child_cells"]==correspondence["child_cells"]
    @test restored.witness["correspondence"]["raw_child_cells"]==inverse["child_cells"]
    target=deepcopy(result.document);target["metadata"]["source"]="altered target"
    @test_throws ArgumentError H.marked_inverse_dehn_twist(target,inverse)
    wrong=deepcopy(inverse);wrong["product_triangles"][1]["points"][1][1]="999"
    @test_throws ArgumentError H.marked_inverse_dehn_twist(result.document,wrong)
    @test_throws ArgumentError H.marked_dehn_twist(doc,Dict("core"=>"a","sign"=>true))
    @test_throws ArgumentError H.marked_dehn_twist(doc,Dict("core"=>17,"sign"=>1))
    malformed=deepcopy(doc);malformed["curves"][1]["id"]=17
    @test_throws ArgumentError H.marked_dehn_twist(malformed,Dict("core"=>"a","sign"=>1))
    @test_throws ArgumentError H.marked_dehn_twist(doc,Dict("core"=>"a","sign"=>1,"family"=>"beta"))
    println("Negative whole-diagram twist started");flush(stdout)
    negative_started=time()
    negative=H.marked_dehn_twist(doc,Dict("core"=>"a","sign"=>-1))
    println("Negative twist constructor seconds: ",time()-negative_started);flush(stdout)
    action=negative.witness["correspondence"]["old_basis_action"]
    a=only(r for r in action if r["curve"]=="a");b=only(r for r in action if r["curve"]=="b")
    @test b["image_exponents"]==b["source_exponents"]-a["source_exponents"]
end

@testset "Whole twist carries embedded D across an auxiliary core" begin
    doc=B.grid_document(4)
    B.add_curve!(doc,"a","alpha",true,[[(1,1),(3,1),(3,3),(1,3),(1,1)]])
    B.add_curve!(doc,"b","beta",true,[[(0,2),(4,2)]])
    B.add_curve!(doc,"cut","aux",true,[[(2,0),(2,4)]])
    B.declare_crossing!(doc,"ab_right",3,2,"a","b",-1)
    B.declare_crossing!(doc,"ab_left",1,2,"a","b",1)
    B.declare_crossing!(doc,"ac_bottom",2,1,"a","cut",1)
    B.declare_crossing!(doc,"ac_top",2,3,"a","cut",-1)
    B.declare_crossing!(doc,"bc",2,2,"b","cut",1)
    doc["region_marks"]=[Dict("id"=>"D","seed"=>["P",10],"claim"=>"embedded_closed_disk")]
    doc["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    @test H.validate(doc;level=:surface).valid
    println("Embedded-D whole-diagram twist started");flush(stdout)
    disk_started=time()
    result=H.marked_dehn_twist(doc,Dict("core"=>"cut","sign"=>1))
    println("Embedded-D twist constructor seconds: ",time()-disk_started);flush(stdout)
    @test H.validate(result.document;level=:surface).valid
    @test result.document["region_marks"][1]["claim"]=="embedded_closed_disk"
    @test [x["id"] for x in result.document["crossings"]]==[x["id"] for x in doc["crossings"]]
    @test length(result.witness["correspondence"]["regions"])==length(H.complement_regions(doc))
    println("Embedded-D twist generated inverse started");flush(stdout)
    disk_inverse_started=time()
    @test H.marked_inverse_dehn_twist(result.document,result.witness["inverse_specification"]).document==doc
    println("Embedded-D twist inverse check seconds: ",time()-disk_inverse_started);flush(stdout)
end
end
