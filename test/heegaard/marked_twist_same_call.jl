module MarkedTwistSameCallTests
using Test,Topology
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams

# Explicit identity affine proof for this SMALL geometry/feature recomposition
# test. Actual Dehn inverse/D acceptance remains in the unchanged main tests.
function identity_raw(doc)
    st=H._marked_states(doc).geometry;records=Any[]
    unit=[[i==j ? "1" : "0" for j in 1:3] for i in 1:3]
    children=Dict((f,i-1)=>[(f,i-1)] for (f,i) in keys(st.triangles))
    for ((f,i),t) in sort(collect(st.triangles);by=first)
        push!(records,Dict("child"=>[f,i-1],"source"=>[f,i-1],"target"=>[f,i-1],"shift"=>0,
            "vertices"=>collect(t),"source_weights"=>deepcopy(unit),"target_weights"=>deepcopy(unit)))
    end
    images=[Dict("curve"=>c["id"],"source_word_format"=>c["word"],
        "edge_images"=>[Dict("source"=>[d.face,d.u,d.v],"target"=>[[d.face,d.u,d.v]],
            "partition"=>[Dict("from"=>"0","to"=>"1","dart"=>[d.face,d.u,d.v])]) for d in M.parse_word(c["word"])]) for c in doc["curves"]]
    (;construction_only=true,complex=deepcopy(doc),records,vertex_records=Any[],
        image_records=images,source_children=children,target_children=deepcopy(children))
end

@testset "Private SAME-CALL geometry reuse retains actual map and feature checks" begin
    B=H.MarkedDiagramBuilders;doc=B.grid_document(2)
    B.add_curve!(doc,"a","alpha",true,[[(0,1),(2,1)]])
    B.add_curve!(doc,"b","beta",true,[[(1,0),(1,2)]])
    B.declare_crossing!(doc,"x",1,1,"a","b",1)
    doc["metadata"]=Dict("source"=>"fresh feature recomposition","display"=>Dict("infinity_triangle"=>["P",0]))
    raw=identity_raw(doc);normal=H._mt_rechart_overlay(doc,raw)
    recovered=H._marked_unrechart_geometry(normal.document,normal.rechart.inverse_specification)
    geometry=Dict("surface"=>raw.complex["surface"],"patches"=>raw.complex["patches"])
    pairs=H._marked_rechart_opposites(geometry)
    @test recovered.fresh_forward!==normal.rechart
    @test recovered.fresh_forward.correspondence===recovered.correspondence
    @test H._mt_same_call_rechart(geometry,pairs,Dict(),recovered)===recovered.fresh_forward
    reused=H._mt_rechart_overlay(doc,raw,recovered)
    @test reused.document==normal.document
    @test reused.records==normal.records
    @test reused.image_records==normal.image_records
    @test reused.source_children==normal.source_children && reused.target_children==normal.target_children
    @test reused.rechart.inverse_specification==normal.rechart.inverse_specification
    @test H.validate(reused.document;level=:surface).valid
    # Native-small acceptance keeps its existing no-rechart path even if
    # a local recovered result exists. The actual oversized D inverse is
    # the separate runtime gate for propagation through common acceptance.
    @test H._mt_accept_map_complex(doc,raw,recovered).document==doc

    # Geometry equality does NOT permit reusing the final feature document or
    # S/T. Metadata is attached freshly; a folded fresh source map must reject.
    changed=deepcopy(raw);changed.complex["metadata"]["source"]="changed fresh metadata"
    @test H._mt_rechart_overlay(doc,changed,recovered).document["metadata"]["source"]=="changed fresh metadata"
    @test normal.document["metadata"]["source"]=="fresh feature recomposition"
    folded=deepcopy(raw);folded.records[1]["source_weights"][1]=["0","1","0"]
    @test_throws ArgumentError H._mt_rechart_overlay(doc,folded,recovered)
    other=deepcopy(geometry);other["surface"]["marking_id"]*="-changed"
    @test_throws ArgumentError H._mt_same_call_rechart(other,pairs,Dict(),recovered)
    wrongpairs=deepcopy(pairs);o=first(keys(wrongpairs));wrongpairs[o]=o
    @test_throws ArgumentError H._mt_same_call_rechart(geometry,wrongpairs,Dict(),recovered)
    @test_throws ArgumentError H._mt_same_call_rechart(geometry,pairs,Dict("strategy"=>"proper_arc"),recovered)
    changedmaps=deepcopy(recovered);changedmaps.correspondence["cell_semantics"]="changed after verification"
    @test_throws ArgumentError H._mt_same_call_rechart(geometry,pairs,Dict(),changedmaps)
end
end
