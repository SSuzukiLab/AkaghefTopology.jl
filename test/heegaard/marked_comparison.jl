module MarkedComparisonTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
loadcase(name)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})
wire(x)=JSON3.read(JSON3.write(x),Dict{String,Any})

# Sources: existing native fixtures and marked_moves.jl's independently checked
# quarter-turn mesh. The expected distinction is topological: exchanging the
# two labelled torus slopes is not identity relative to the fixed marking.
function rotation_fixture()
    n=4;doc=B.grid_document(n);triangles=Any[]
    for y in 0:n-1,x in 0:n-1
        a,b,c,d=B.vertex(x,y),B.vertex(x+1,y),B.vertex(x+1,y+1),B.vertex(x,y+1)
        append!(triangles,iseven(x+y) ? [[a,b,c],[a,c,d]] : [[a,b,d],[b,c,d]])
    end
    doc["patches"][1]["triangles"]=triangles
    B.add_curve!(doc,"a","alpha",true,[[(0,2),(4,2)]])
    B.add_curve!(doc,"b","beta",true,[[(2,0),(2,4)]])
    B.declare_crossing!(doc,"x",2,2,"a","b",1)
    doc,Dict("vertices"=>Dict("P"=>Dict(B.vertex(x,y)=>B.vertex(n-y,x) for y in 0:n for x in 0:n)))
end

# Source: test/heegaard/marked_slides.jl::redundant_torus. This independently
# authored input has two parallel same-family meridians and a transverse beta;
# only the band is supplied to the actual constructor, not a replacement word.
function slide_fixture()
    d=B.grid_document(3)
    B.add_curve!(d,"a1","alpha",true,[[(0,1),(3,1)]])
    B.add_curve!(d,"a2","alpha",true,[[(0,2),(3,2)]])
    B.add_curve!(d,"b","beta",true,[[(2,0),(2,3)]])
    B.declare_crossing!(d,"x1",2,1,"a1","b",1)
    B.declare_crossing!(d,"x2",2,2,"a2","b",1)
    d["region_marks"]=[Dict("id"=>"R","seed"=>["P",8],"claim"=>"region")]
    d["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",17]))
    d
end

# Source: test/heegaard/marked_stabilization.jl::disk_marked_sphere. Its cap is
# glued with opposite orientation; the protected selected triangle is disjoint
# from the actual connected-sum support. No post-stabilization object is input.
function stabilization_fixture()
    d=B.grid_document(3;surface=:disk)
    B.add_curve!(d,"old_alpha","alpha",true,[[(1,1),(2,2),(2,1),(1,1)]])
    p=only(d["patches"]);bd=[v for s in p["side_order"] for v in p["sides"][s][1:end-1]]
    push!(d["patches"],Dict("face"=>"Q","side_order"=>deepcopy(p["side_order"]),"sides"=>deepcopy(p["sides"]),
        "triangles"=>[["cap_center",bd[i],bd[mod1(i+1,length(bd))]] for i in eachindex(bd)]))
    d["surface"]=Dict("marking_id"=>"protected-triangle-sphere","genus"=>0,"boundary"=>Any[],
        "gluings"=>[[["P","south"],["Q","south"]],[["P","east"],["Q","west"]],
                    [["P","north"],["Q","north"]],[["P","west"],["Q","east"]]])
    d["region_marks"]=[Dict("id"=>"D","seed"=>["P",8],"claim"=>"embedded_closed_disk")]
    d["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    d
end

@testset "Data identity, fixed marking and constructive comparison" begin
    source=loadcase("torus_mer_lon");snapshot=deepcopy(source)
    same=H.marked_comparison_report(source,deepcopy(source))
    @test same["identity"]["full_data_equal"] && same["identity"]["mathematical_representation_equal"]
    @test same["marked_isotopy"]["status"]=="witnessed"
    @test same["represented_manifold"]["status"]=="witnessed"
    metadata=deepcopy(source);metadata["metadata"]=Dict{String,Any}("source"=>"another locator")
    different=H.marked_comparison_report(source,metadata)
    @test !different["identity"]["full_data_equal"] && different["identity"]["mathematical_representation_equal"]
    @test different["marked_isotopy"]["status"]=="witnessed"
    infinity=deepcopy(metadata);infinity["metadata"]["display"]=Dict("infinity_triangle"=>["P",0])
    @test H.marked_comparison_report(source,infinity)["marked_isotopy"]["status"]=="unknown"
    moved=H.marked_triangle_isotopy(source,Dict("curve"=>"beta1","face"=>"P","triangle"=>10,"path"=>["v2_1","v2_2"]))
    absent=H.marked_comparison_report(source,moved.document)
    @test !absent["identity"]["mathematical_representation_equal"]
    @test absent["marked_isotopy"]["status"]=="unknown"
    @test absent["represented_manifold"]["status"]=="consistent"
    report=H.marked_comparison_report(source,moved.document;witnesses=[wire(moved.witness)])
    @test report["witness_replay"]["status"]=="verified"
    @test report["marked_isotopy"]["status"]=="witnessed"
    @test report["move_equivalence"]["relation"]=="fixed_model_ambient_isotopy"
    @test report["represented_manifold"]["status"]=="witnessed"
    @test !report["represented_manifold"]["witness_invariant_conflict"]
    nonobject=H.marked_comparison_report(source,moved.document;witnesses=[7])
    @test nonobject["witness_replay"]["status"]=="rejected"
    @test occursin("must be an object",nonobject["witness_replay"]["reason"])
    @test nonobject["move_equivalence"]["status"]=="unknown"
    for field in ("specification","inverse_specification"),value in (7,Any[])
        malformed=deepcopy(moved.witness);malformed[field]=value
        malformed_report=H.marked_comparison_report(source,moved.document;witnesses=[malformed])
        @test malformed_report["witness_replay"]["status"]=="rejected"
        @test occursin("$field must be an object",malformed_report["witness_replay"]["reason"])
    end
    @test H.marked_comparison_report(source,source;witnesses=[moved.witness,H.inverse_marked_move(moved.witness)])["witness_replay"]["status"]=="verified"
    bad=deepcopy(moved.witness);bad["inverse_specification"]["path"]=["v2_1","v2_2"]
    rejected=H.marked_comparison_report(source,moved.document;witnesses=[bad])
    @test rejected["witness_replay"]["status"]=="rejected"
    @test rejected["move_equivalence"]["status"]=="unknown"
    bad=deepcopy(moved.witness);bad["correspondence"]["triangles"][1]["target"]=["P",3]
    @test H.marked_comparison_report(source,moved.document;witnesses=[bad])["witness_replay"]["status"]=="rejected"
    @test H.marked_comparison_report(source,source;witnesses=[moved.witness])["witness_replay"]["status"]=="rejected"
    doc,map=rotation_fixture();action=H.marked_simplicial_action(doc,map)
    rotation=H.marked_comparison_report(doc,action.document;witnesses=[action.witness])
    @test rotation["marked_isotopy"]["status"]=="unknown"
    @test rotation["move_equivalence"]["relation"]=="oriented_surface_map"
    @test rotation["represented_manifold"]["status"]=="witnessed"
    @test source==snapshot
    invalid=deepcopy(source);invalid["crossings"][1]["sign"]*=-1
    @test H.marked_comparison_report(source,invalid)["marked_isotopy"]["status"]=="inapplicable"
    # Source: marked_topology.jl::_marked_states reserved metadata validation.
    # Invalid infinity metadata is not silently treated as an absent choice.
    for metadata in (Any[],7,nothing,Dict("display"=>Any[]),Dict("display"=>7),
            Dict("display"=>nothing),Dict("display"=>Dict("infinity_triangle"=>["P",999999])),
            Dict("display"=>Dict("infinity_triangle"=>["P",0,1])),
            Dict("display"=>Dict("infinity_triangle"=>["P","0"])))
        malformed=deepcopy(source);malformed["metadata"]=metadata
        invalid_report=H.marked_comparison_report(malformed,malformed)
        @test invalid_report["identity"]["status"]=="inapplicable"
        @test invalid_report["identity"]["explicit_infinity_equal"]===nothing
        @test invalid_report["marked_isotopy"]["status"]=="inapplicable"
        @test invalid_report["move_equivalence"]["status"]=="inapplicable"
        @test invalid_report["represented_manifold"]["status"]=="inapplicable"
    end
    @test !H.marked_comparison_report(B.grid_document(4;surface=:disk),source)["represented_manifold"]["applicable"]
    obstruction=H.marked_comparison_report(source,loadcase("torus_parallel_zero_crossings"))
    @test obstruction["represented_manifold"]["status"]=="refuted"
    @test !obstruction["represented_manifold"]["h1_equal"]
end

@testset "Underlying handle-slide equivalence and explicit region continuation" begin
    source=slide_fixture()
    moved=H.marked_handle_slide(source,Dict("source"=>"a1","target"=>"a2","band"=>"P{v1_1 v1_2}","disk_policy"=>"follow_seed"))
    report=H.marked_comparison_report(source,moved.document;witnesses=[wire(moved.witness)])
    @test report["witness_replay"]["status"]=="verified"
    @test report["move_equivalence"]["relation"]=="underlying_heegaard_move"
    @test report["marked_isotopy"]["status"]=="unknown"
    @test report["marked_identification"]["status"]=="unknown"
    @test report["represented_manifold"]["status"]=="witnessed"
    @test !report["represented_manifold"]["witness_invariant_conflict"]
    @test report["region_decoration"]["status"]=="continued"
    @test report["region_decoration"]["original_regions_preserved"]===false
    @test any(!r["region_preserved"] for r in report["region_decoration"]["steps"][1]["evidence"]["regions"])
    inverse=H.marked_comparison_report(moved.document,source;witnesses=[wire(H.inverse_marked_move(moved.witness))])
    @test inverse["move_equivalence"]["relation"]=="underlying_heegaard_move"
    @test inverse["represented_manifold"]["status"]=="witnessed"
    forged=deepcopy(moved.witness);forged["correspondence"]["marked_regions"][1]["region_preserved"]=true
    @test H.marked_comparison_report(source,moved.document;witnesses=[forged])["witness_replay"]["status"]=="rejected"
end

@testset "Underlying stabilization preserves manifold and protected D, not surface genus" begin
    source=stabilization_fixture();moved=H.marked_stabilize(source,Dict("support"=>["Q",0]))
    report=H.marked_comparison_report(source,moved.document;witnesses=[wire(moved.witness)])
    @test report["witness_replay"]["status"]=="verified"
    @test report["move_equivalence"]["relation"]=="underlying_stabilization_move"
    @test report["marked_isotopy"]["status"]=="unknown"
    @test report["marked_identification"]["status"]=="unknown"
    @test report["represented_manifold"]["status"]=="witnessed"
    @test !report["represented_manifold"]["witness_invariant_conflict"]
    @test report["region_decoration"]["status"]=="preserved"
    @test report["region_decoration"]["original_regions_preserved"]===true
    @test H.surface_genus(moved.document)==H.surface_genus(source)+1
    inverse=H.marked_comparison_report(moved.document,source;witnesses=[wire(H.inverse_marked_move(moved.witness))])
    @test inverse["move_equivalence"]["relation"]=="underlying_stabilization_move"
    @test inverse["represented_manifold"]["status"]=="witnessed"
    @test inverse["region_decoration"]["status"]=="preserved"
end

@testset "Explicit marked identification differs from fixed-model isotopy" begin
    # Source: accepted marked_subdivision.jl constructors. Their triangle and
    # quotient-star maps identify the old surface, while the native mesh and
    # in the dual-star case marking_id change. Test actual replay, not hashes.
    source=loadcase("torus_mer_lon")
    source["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    for transformed in (H.marked_subdivide(source),H.marked_dual_star(source))
        report=H.marked_comparison_report(source,transformed.document;witnesses=[wire(transformed.witness)])
        @test report["witness_replay"]["status"]=="verified"
        @test report["marked_identification"]["status"]=="witnessed"
        @test report["marked_isotopy"]["status"]=="unknown"
        @test !report["marked_isotopy"]["fixed_model_equal"]
        @test report["move_equivalence"]["relation"]=="explicit_marked_identification"
        @test report["represented_manifold"]["status"]=="witnessed"
        inverse=H.inverse_marked_move(transformed.witness)
        restored=H.marked_comparison_report(transformed.document,source;witnesses=[wire(inverse)])
        @test restored["marked_identification"]["status"]=="witnessed"
        forged=deepcopy(transformed.witness)
        forged["correspondence"]["triangles"][1]["target"]=[["nonexistent",0]]
        rejected=H.marked_comparison_report(source,transformed.document;witnesses=[forged])
        @test rejected["witness_replay"]["status"]=="rejected"
        @test rejected["marked_identification"]["status"]=="unknown"
    end
    refined=H.marked_subdivide(source)
    roundtrip=H.marked_comparison_report(source,source;
        witnesses=[refined.witness,H.inverse_marked_move(refined.witness)])
    @test roundtrip["witness_replay"]["status"]=="verified"
    @test roundtrip["identity"]["full_data_equal"]

    rotated,rotation=rotation_fixture();subdivided=H.marked_subdivide(rotated)
    acted=H.marked_simplicial_action(rotated,rotation)
    mixed=H.marked_comparison_report(rotated,acted.document;
        witnesses=[subdivided.witness,H.inverse_marked_move(subdivided.witness),acted.witness])
    @test mixed["witness_replay"]["status"]=="verified"
    @test mixed["move_equivalence"]["relation"]=="oriented_surface_map_under_identification"
    @test mixed["marked_isotopy"]["status"]=="unknown"

    # Source: vendor realize and marked_moves.jl::marked_vertex_isotopy. The
    # second witness is valid on its refined mesh. The comparison must require
    # that mesh's positions instead of carrying the original realization over.
    M=H.MarkedSurfaceDiagrams
    original_positions=M.realize(H._marked_states(source).geometry)
    refined_positions=M.realize(H._marked_states(refined.document).geometry)
    point=refined_positions["P"]["v1_1"]
    moved_vertex=H.marked_vertex_isotopy(refined.document,refined_positions,
        Dict("face"=>"P","vertex"=>"v1_1","target"=>[M.qstr(point[1]+1//1000000000),M.qstr(point[2])]))
    stale=H.marked_comparison_report(source,moved_vertex.document;
        witnesses=[refined.witness,moved_vertex.witness],positions=original_positions)
    @test stale["witness_replay"]["status"]=="rejected"
    @test stale["witness_replay"]["failed_step"]==2
    @test occursin("requires source positions",stale["witness_replay"]["reason"])
end

@testset "Independent original o-data relations and truthful H1 evidence" begin
    # Source: registered abalone graph and supplied Koda figure6.odata.json,
    # already used by the independent reverse lane. No reverse output feeds
    # this extraction; the genus-3 marked original is a separate fixed input.
    graph=H.OData([[1,-1]],[1],Dict{String,Any}())
    original=deepcopy(graph.gauss);raw=H.odata_computation_report(graph)
    @test raw["applicable"] && !raw["source_metadata_used"]
    @test raw["first_homology"]["free_rank"]==0 && isempty(raw["first_homology"]["torsion"])
    @test raw["presentation"]["metadata"]["pd_row_slots"]==["g","k","l","h"]
    @test length(raw["presentation"]["relators"])==2
    poisoned=H.OData(graph.gauss,graph.orientation,Dict{String,Any}("marked_source"=>Dict("invalid"=>true),"source_hash"=>"unrelated"))
    @test H.odata_computation_report(poisoned)==raw
    report=H.marked_comparison_report(loadcase("torus_mer_lon"),poisoned)
    @test report["represented_manifold"]["status"]=="consistent"
    @test report["marked_isotopy"]["status"]=="unknown"
    @test !report["represented_manifold"]["group_isomorphism_decided"]
    @test !report["represented_manifold"]["manifold_recognition_performed"]
    @test !report["represented_manifold"]["orientation_distinguished_by_h1"]
    mismatch=H.marked_comparison_report(loadcase("torus_parallel_zero_crossings"),graph)
    @test mismatch["represented_manifold"]["status"]=="refuted"
    source=JSON3.read(read(joinpath(@__DIR__,"../../vendor/koda_genus3/data/figure6.odata.json"),String))
    kg=H.OData([reverse(Int.(w)) for w in source[1]],Int.(source[2]),Dict{String,Any}())
    kh=H.odata_computation_report(kg)
    @test kh["applicable"]
    @test kh["first_homology"]["free_rank"]==0 && kh["first_homology"]["torsion"]==[2,2]
    marked=JSON3.read(read(joinpath(@__DIR__,"../../vendor/koda_genus3/data/figure6.fullmarked.json"),String),Dict{String,Any})
    @test H.marked_comparison_report(marked,kg)["represented_manifold"]["status"]=="consistent"
    # Reverse-conversion integration: the left object comes from actual dual tetrahedra, but
    # the right invariant is still extracted from ORIGINAL raw odata. Its mixed
    # signs exercise a distinct relation path; no reverse output feeds both.
    independent=H.marked_comparison_report(H.to_heegaard(kg),kg)
    @test independent["represented_manifold"]["status"]=="consistent"
    @test independent["represented_manifold"]["left_h1"]==Dict("free_rank"=>0,"torsion"=>[2,2])
    @test independent["represented_manifold"]["right_h1"]==Dict("free_rank"=>0,"torsion"=>[2,2])
    @test !H.odata_computation_report(H.OData(Vector{Int}[],Int[],Dict{String,Any}()))["applicable"]
    @test !H.odata_computation_report(H.OData([[1,-1]],[0],Dict{String,Any}()))["applicable"]
    @test graph.gauss==original
    @test wire(report)["independence"]["retained_source_used"]==false
end
end
