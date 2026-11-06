module MarkedSlideTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const M=H.MarkedSurfaceDiagrams

function redundant_torus()
    doc=B.grid_document(3)
    B.add_curve!(doc,"a1","alpha",true,[[(0,1),(3,1)]])
    B.add_curve!(doc,"a2","alpha",true,[[(0,2),(3,2)]])
    B.add_curve!(doc,"b","beta",true,[[(2,0),(2,3)]])
    B.declare_crossing!(doc,"x1",2,1,"a1","b",1)
    B.declare_crossing!(doc,"x2",2,2,"a2","b",1)
    doc["region_marks"]=[Dict("id"=>"R","seed"=>["P",8],"claim"=>"region")]
    doc["metadata"]=Dict("source"=>"two parallel oriented meridians and a longitude",
        "display"=>Dict("infinity_triangle"=>["P",17]))
    doc
end

# Source: marked_slides.jl's actual band-neighborhood method. The two oriented
# torus meridians have the same nonzero class; their embedded band-sum frontier
# is contractible, and the retained second meridian still meets beta once.
# No post-slide curve, region partition or predetermined crossing list is input.
@testset "Constructed handle slide, explicit disk policy and geometric inverse" begin
    original=redundant_torus();snapshot=deepcopy(original)
    spec=Dict("source"=>"a1","target"=>"a2","band"=>"P{v1_1 v1_2}")
    @test H.validate(original;level=:surface).valid
    @test H.intersection_matrix(original)==reshape([1,1],2,1)
    @test_throws ArgumentError H.marked_handle_slide(original,spec) # preserve may not silently change R
    @test original==snapshot
    spec["disk_policy"]="follow_seed"
    result=H.marked_handle_slide(original,spec)
    @test original==snapshot
    @test H.validate(result.document;level=:surface).valid
    @test H.intersection_matrix(result.document)==reshape([0,1],2,1)
    @test length(result.document["crossings"])==3
    before_h1=H.first_homology(original);after_h1=H.first_homology(result.document)
    @test before_h1.free_rank==after_h1.free_rank==0
    @test isempty(after_h1.torsion)
    @test any(!r["region_preserved"] && r["policy"]=="follow_seed" for r in result.witness["correspondence"]["marked_regions"])
    @test length(result.witness["correspondence"]["pair_of_pants_boundary_cycles"])==3
    inverse=JSON3.read(JSON3.write(result.witness["inverse_specification"]),Dict{String,Any})
    @test H.marked_inverse_handle_slide(result.document,inverse).document==original
    altered=deepcopy(result.document);altered["metadata"]["source"]="changed target"
    @test_throws ArgumentError H.marked_inverse_handle_slide(altered,inverse)
    wrong_core=deepcopy(inverse);wrong_core["core_graph"]["source"][1][3]="nonexistent_vertex"
    @test_throws ArgumentError H.marked_inverse_handle_slide(result.document,wrong_core)
    wrong_pants=deepcopy(inverse);wrong_pants["pair_of_pants_triangles"][1]=["no_such_chart",0]
    @test_throws ArgumentError H.marked_inverse_handle_slide(result.document,wrong_pants)
    invalid=deepcopy(spec);invalid["target"]="a1"
    @test_throws ArgumentError H.marked_handle_slide(original,invalid)
    invalid=deepcopy(spec);invalid["target"]="b"
    @test_throws ArgumentError H.marked_handle_slide(original,invalid)
    invalid=deepcopy(spec);invalid["band"]="P{v1_1 v2_1}"
    @test_throws ArgumentError H.marked_handle_slide(original,invalid)
end

# Source: existing lickorish_reference(2) Fig.12.3 model, with agent test roles
# alpha=l1,l2; beta=m1,m2; n1 remains auxiliary. These roles are not attributed
# to the book. HA6's independently specified band runs right-to-right between
# the two alpha curves. Thus the positive source band sum is l1+l2, with
# intersection rows I2 -> [1 1; 0 1]. The complete cut systems must remain real
# complete systems under actual normalized cut/cap topology.
@testset "Genus-two nonseparating handle slide on a closed splitting" begin
    doc=H.lickorish_reference(2)
    for c in doc["curves"]
        startswith(c["id"],"l") && (c["role"]="alpha")
        startswith(c["id"],"m") && (c["role"]="beta")
    end
    doc=H.marked_subdivide(doc).document
    band="V_ln_1_1{subedge21 subface1 subedge3 subface12 subedge1} | F1{subedge4 subface2 z subface3 subedge6} | V_ln_2_1{subedge14 subface9 subedge16 subface8 subedge23}"
    before=H.marked_topology_report(doc)
    @test before["compression"].alpha_complete && before["compression"].beta_complete
    @test H.intersection_matrix(doc)==[1 0;0 1]
    input_state=H._marked_states(doc).geometry
    @test length(unique(M.root!(input_state.vertex_dsu,v) for v in keys(input_state.vertex_dsu.p)))==190
    original=deepcopy(doc)
    result=H.marked_handle_slide(doc,Dict("source"=>"l1","target"=>"l2","band"=>band))
    @test doc==original
    @test H.intersection_matrix(result.document)==[1 1;0 1]
    report=H.marked_topology_report(result.document)
    @test report["compression"].alpha_complete && report["compression"].beta_complete
    @test H.surface_genus(result.document)==2
    @test length(result.document["patches"])<=64
    atlas=first(result.witness["correspondence"]["construction"])["correspondence"]
    @test !isempty(atlas["chart_merges"])
    @test all(c["one_subdivision_vertices"]<=256 for c in atlas["chart_correspondence"]["cell_counts"])
    st=H._marked_states(result.document).geometry
    barriers=Set(e for (e,owners) in st.used if "l1" in owners)
    cut=M.regions(st;barriers=barriers)
    @test length(cut)==1 && only(cut)["genus"]==1 && only(cut)["boundary_components"]==2
    @test H.first_homology(result.document).free_rank==0
    @test H.marked_inverse_handle_slide(result.document,result.witness["inverse_specification"]).document==doc
end
end
