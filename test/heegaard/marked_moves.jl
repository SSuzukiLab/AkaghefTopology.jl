module MarkedMoveTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const B=H.MarkedDiagramBuilders

# Source: native fixtures/builders and marked_topology.jl::_marked_states.
# These probes use actual changed edge routes, empty-disk support and a
# nontrivial fixed-marking torus action; expected outputs are not move inputs.
function disk_region_fixture()
    doc=B.grid_document(6)
    B.add_curve!(doc,"a","alpha",true,[[(1,2),(4,2),(4,4),(1,4),(1,2)]])
    doc["region_marks"]=[Dict("id"=>"D","seed"=>["P",28],"claim"=>"embedded_closed_disk")]
    doc["metadata"]=Dict("source"=>Dict("locator"=>"marked move test"),
        "display"=>Dict("infinity_triangle"=>["P",28]))
    doc
end

function symmetric_torus()
    n=4;doc=B.grid_document(n)
    triangles=Any[]
    for y in 0:n-1,x in 0:n-1
        a,b,c,d=B.vertex(x,y),B.vertex(x+1,y),B.vertex(x+1,y+1),B.vertex(x,y+1)
        append!(triangles,iseven(x+y) ? [[a,b,c],[a,c,d]] : [[a,b,d],[b,c,d]])
    end
    doc["patches"][1]["triangles"]=triangles
    B.add_curve!(doc,"a","alpha",true,[[(0,2),(4,2)]])
    B.add_curve!(doc,"b","beta",true,[[(2,0),(2,4)]])
    B.declare_crossing!(doc,"x",2,2,"a","b",1)
    doc["region_marks"]=[Dict("id"=>"R","seed"=>["P",0],"claim"=>"region")]
    doc["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    rotation=Dict("vertices"=>Dict("P"=>Dict(B.vertex(x,y)=>B.vertex(n-y,x) for y in 0:n for x in 0:n)))
    doc,rotation
end

jsonround(x)=JSON3.read(JSON3.write(x),Dict{String,Any})
function abelian_word(data,cid)
    out=zeros(Int,length(data.presentation.generators))
    for (i,s) in data.curve_words[cid];out[i]+=s;end
    out
end

@testset "Constructive marked triangle isotopy, selected D and replay" begin
    original=disk_region_fixture();snapshot=deepcopy(original)
    spec=Dict("curve"=>"a","face"=>"P","triangle"=>28,"path"=>["v2_2","v3_2"])
    @test H.validate(original;level=:surface).valid
    moved=H.marked_triangle_isotopy(original,spec)
    @test original==snapshot
    @test moved.document["curves"][1]["word"]!=original["curves"][1]["word"]
    @test H.validate(moved.document;level=:surface).valid
    @test moved.document["region_marks"][1]["seed"]!=["P",28]
    @test moved.document["metadata"]["display"]["infinity_triangle"]!=["P",28]
    @test moved.document["metadata"]["source"]==original["metadata"]["source"]
    @test H.surface_genus(moved.document)==H.surface_genus(original)==1
    @test H.intersection_matrix(moved.document)==H.intersection_matrix(original)
    @test H.marked_curve_words(moved.document).curve_words==H.marked_curve_words(original).curve_words
    @test length(moved.witness["correspondence"]["triangles"])==72
    @test length(moved.witness["correspondence"]["regions"])==2
    @test H.replay_marked_move(original,jsonround(moved.witness)).document==moved.document
    altered_witness=deepcopy(moved.witness)
    altered_witness["correspondence"]["triangles"][1]["target"]=["P",3]
    @test_throws ArgumentError H.replay_marked_move(original,altered_witness)
    altered_witness=deepcopy(moved.witness)
    altered_witness["inverse_specification"]["path"]=["v2_2","v3_2"]
    @test_throws ArgumentError H.replay_marked_move(original,altered_witness)
    restored=H.replay_marked_move(moved.document,jsonround(H.inverse_marked_move(moved.witness)))
    @test restored.document==original
    tampered=deepcopy(original);tampered["metadata"]["source"]["locator"]="changed source"
    @test_throws ArgumentError H.replay_marked_move(tampered,moved.witness)
    bad=deepcopy(spec);bad["triangle"]=0
    @test_throws ArgumentError H.marked_triangle_isotopy(original,bad)
    bad=deepcopy(spec);bad["path"]=["v2_2","v4_2"]
    @test_throws ArgumentError H.marked_triangle_isotopy(original,bad)

    torus=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples/torus_mer_lon.json"),String),Dict{String,Any})
    # Occupied alternate vertex on beta cannot be swept by alpha.
    @test_throws ArgumentError H.marked_triangle_isotopy(torus,
        Dict("curve"=>"alpha1","face"=>"P","triangle"=>18,"path"=>["v1_2","v2_2"]))
    # A legal triangle move adjoining an existing crossing retains its sign.
    crossing_move=H.marked_triangle_isotopy(torus,
        Dict("curve"=>"beta1","face"=>"P","triangle"=>10,"path"=>["v2_1","v2_2"]))
    @test crossing_move.document["crossings"]==torus["crossings"]
    @test H.replay_marked_move(crossing_move.document,H.inverse_marked_move(crossing_move.witness)).document==torus
    # Infinity follows the full geometric region when an aux loop splits the
    # single alpha-union-beta region; it may not jump to minimum triangle 0.
    auxiliary=disk_region_fixture();delete!(auxiliary,"region_marks")
    auxiliary["curves"][1]["role"]="aux"
    moved_aux=H.marked_triangle_isotopy(auxiliary,spec)
    infinity=moved_aux.document["metadata"]["display"]["infinity_triangle"]
    geometry=H._marked_states(moved_aux.document).geometry
    selected=only(r for r in M.regions(geometry) if infinity in r["triangles"])
    @test selected["embedded_closed_disk"]
    @test H.replay_marked_move(moved_aux.document,H.inverse_marked_move(moved_aux.witness)).document==auxiliary
end

@testset "Exact realization ambient isotopy is independently labelled" begin
    doc=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples/torus_mer_lon.json"),String),Dict{String,Any})
    positions=M.realize(H._marked_states(doc).geometry);old=positions["P"]["v1_1"]
    qstr=M.qstr
    spec=Dict("face"=>"P","vertex"=>"v1_1","target"=>[qstr(old[1]+1//10000),qstr(old[2])])
    moved=H.marked_vertex_isotopy(doc,positions,spec)
    @test moved.document==doc
    @test moved.positions!=positions
    @test M.certify_layout(H._marked_states(doc).geometry,moved.positions)
    @test H.replay_marked_move(doc,jsonround(moved.witness);positions=positions).positions==moved.positions
    @test H.replay_marked_move(doc,H.inverse_marked_move(moved.witness);positions=moved.positions).positions==positions
    @test_throws ArgumentError H.replay_marked_move(doc,moved.witness;positions=moved.positions)
    @test_throws ArgumentError H.replay_marked_move(doc,moved.witness)
    @test_throws M.DiagramError H.marked_vertex_isotopy(doc,positions,
        Dict("face"=>"P","vertex"=>"v0_1","target"=>["0","0"]))
    @test_throws M.DiagramError H.marked_vertex_isotopy(doc,positions,
        Dict("face"=>"P","vertex"=>"v1_1","target"=>["100","100"]))
end

@testset "Nontrivial oriented action relative to the fixed torus marking" begin
    doc,rotation=symmetric_torus();snapshot=deepcopy(doc)
    before=H.marked_curve_words(doc);moved=H.marked_simplicial_action(doc,rotation)
    @test doc==snapshot
    @test H.validate(moved.document;level=:surface).valid
    @test moved.document["surface"]==doc["surface"]
    @test moved.document["patches"]==doc["patches"]
    @test M.marking_digest(moved.document)==M.marking_digest(doc)
    @test moved.document["curves"][1]["word"]==doc["curves"][2]["word"]
    after=H.marked_curve_words(moved.document)
    @test abelian_word(after,"a")==abelian_word(before,"b")
    @test abelian_word(after,"b")==-abelian_word(before,"a")
    @test H.intersection_matrix(moved.document)==[1;;]
    @test moved.document["crossings"]==doc["crossings"]
    @test moved.document["metadata"]["display"]["infinity_triangle"]!=["P",0]
    @test H.replay_marked_move(doc,jsonround(moved.witness)).document==moved.document
    @test H.replay_marked_move(moved.document,jsonround(H.inverse_marked_move(moved.witness))).document==doc
    iterated=doc
    for _ in 1:4;iterated=H.marked_simplicial_action(iterated,rotation).document;end
    @test iterated==doc
    malformed=deepcopy(rotation);malformed["vertices"]["P"]["v1_1"]="v1_1"
    @test_throws ArgumentError H.marked_simplicial_action(doc,malformed)
    reflection=Dict("vertices"=>Dict("P"=>Dict(B.vertex(x,y)=>B.vertex(4-x,y) for y in 0:4 for x in 0:4)))
    @test_throws ArgumentError H.marked_simplicial_action(doc,reflection)
end
end
