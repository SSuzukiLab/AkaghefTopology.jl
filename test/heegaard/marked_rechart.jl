module MarkedRechartTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
fixture(n)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",n*".json"),String),Dict{String,Any})
jsonround(x)=JSON3.read(JSON3.write(x),Dict{String,Any})
rawgeometry(d)=Dict{String,Any}("surface"=>deepcopy(d["surface"]),"patches"=>deepcopy(d["patches"]))
function reversepath(xs)
    [[x[1],x[3],x[2]] for x in reverse(xs)]
end

@testset "Private identity rechart, exact flags and canonical seam edges" begin
    source=fixture("torus_mer_lon");raw=rawgeometry(source);original=deepcopy(raw)
    pairs=H._marked_rechart_opposites(raw)
    result=H._marked_rechart_geometry(raw,pairs);doc=result.document;corr=result.correspondence
    @test raw==original
    @test H.validate(doc;level=:surface).valid
    @test H.surface_genus(doc)==1
    @test length(doc["patches"])<=64
    @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in doc["patches"])<=256
    @test isempty(doc["curves"]) && isempty(doc["crossings"])
    cells=Dict(Tuple(r["source"])=>r["target"] for r in corr["cell_children"])
    @test length(cells)==32 && all(length(v)==6 for v in values(cells))
    @test length(corr["child_barycentric_rows"])==192
    public=M.validate(doc)
    points=Dict(Tuple(r["source"])=>r["target"] for r in corr["vertex_images"])
    # Independent identity anchors: unit barycentric corners must be the actual
    # public images of the corresponding RAW vertex, not a permuted basis.
    flags_actual=true;positive_sixths=true;anchors_actual=true
    for row in corr["child_barycentric_rows"]
        f,k=row["source"];triangle=only(p for p in raw["patches"] if p["face"]==f)["triangles"][k+1]
        target=Tuple(row["target"])
        flags_actual &= Tuple(row["vertices"])==public.triangles[(target[1],target[2]+1)]
        w=[[M.qnum(x) for x in r] for r in row["weights"]]
        positive_sixths &= H._mr_det(w)==1//6
        for (j,weights) in enumerate(w)
            corner=findfirst(==(1),weights);corner===nothing && continue
            expected=points[(f,triangle[corner])]
            anchors_actual &= M.root!(public.vertex_dsu,(target[1],row["vertices"][j]))==M.root!(public.vertex_dsu,Tuple(expected))
        end
    end
    @test flags_actual
    @test positive_sixths
    @test anchors_actual
    edges=Dict(Tuple(r["source"])=>r["target"] for r in corr["edge_images"])
    # Opposite local seam occurrences are ONE global edge with reversed
    # orientation. Distinct raw endpoints alone are not a valid gauge.
    @test all(edges[(o[1],o[3],o[4])]==reversepath(edges[(p[1],p[3],p[4])]) for (o,p) in pairs)
    @test length(corr["canonical_edges"])==length(pairs)÷2
    @test length(corr["edge_occurrences"])==length(pairs)
    # Restore actual features through the public maps, then validate them.
    # The helper never received these records or a retained source diagram.
    featured=deepcopy(doc);featured["curves"]=deepcopy(source["curves"]);featured["crossings"]=deepcopy(source["crossings"])
    for c in featured["curves"]
        ds=vcat([H._st_darts(edges[(d.face,d.u,d.v)]) for d in M.parse_word(c["word"])]...)
        c["word"]=M.print_word(ds)
    end
    for x in featured["crossings"];x["at"]=deepcopy(points[Tuple(x["at"])]);end
    @test H.validate(featured;level=:surface).valid
    @test H.compression_boundary_report(featured).alpha_complete && H.compression_boundary_report(featured).beta_complete
    restored=H._marked_unrechart_geometry(featured,jsonround(result.inverse_specification))
    @test restored.construction_only
    @test Set(keys(restored.geometry))==Set(["surface","patches"])
    @test restored.geometry==raw && restored.opposite_edges==pairs
    @test H._mm_hash(restored.correspondence)==H._mm_hash(corr)
    # The optional internal result is THIS call's freshly checked regeneration.
    # It is neither the incoming serialized inverse nor a retained source doc.
    @test isempty(restored.forward_specification)
    @test restored.fresh_forward.document==H._st_empty(featured)
    @test restored.fresh_forward.inverse_specification==result.inverse_specification
    @test restored.fresh_forward.correspondence===restored.correspondence
    @test restored.fresh_forward.correspondence!==result.correspondence
    bad=deepcopy(result.inverse_specification);bad["affine_children"][1]["weights"][1]=["0","0","0"]
    @test_throws ArgumentError H._marked_unrechart_geometry(doc,bad)
    # Keep a well-formed four-endpoint involution but pair two actual edges
    # differently. The inverse must reject CURRENT adjacency, not just a cache.
    bad_pairs=deepcopy(result.inverse_specification)
    a=deepcopy(bad_pairs["opposite_edges"][1]["target"])
    bad_pairs["opposite_edges"][1]["target"]=deepcopy(bad_pairs["opposite_edges"][2]["target"])
    bad_pairs["opposite_edges"][2]["target"]=a
    @test_throws ArgumentError H._marked_unrechart_geometry(doc,bad_pairs)
    wrong=deepcopy(pairs);o=first(keys(wrong));wrong[o]=o
    @test_throws ArgumentError H._marked_rechart_geometry(raw,wrong)
    @test_throws ArgumentError H._marked_rechart_geometry(merge(raw,Dict("curves"=>[])),pairs)
    boundary=rawgeometry(H.MarkedDiagramBuilders.grid_document(3;surface=:disk))
    @test_throws ArgumentError H._marked_rechart_opposites(boundary)
end

# This case is explicit acceptance input, not an optional silent fallback.
# The lead supplies the real positive common-refinement receipt after the
# bounded packing probe. Without it, only the small identity scope above ran.
if haskey(ENV,"MARKED_RECHART_PROOF")
    @testset "Actual over-limit twist proof is recharted and genuinely reassembled" begin
        bundle=JSON3.read(read(ENV["MARKED_RECHART_PROOF"],String),Dict{String,Any})
        @test bundle["construction_only"] && !bundle["native_validity_claimed"]
        raw=rawgeometry(bundle["complex"]);snapshot=deepcopy(raw)
        @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in raw["patches"])==511
        @test sum(length(p["triangles"]) for p in raw["patches"])==2016
        pairs=H._marked_rechart_opposites(raw)
        # Source: this actual over-limit regression. Log expensive boundaries
        # without changing fixtures, checks or construction/inverse semantics.
        println("Actual rechart forward started");flush(stdout)
        forward_started=time()
        result=H._marked_rechart_geometry(raw,pairs)
        println("Actual rechart forward seconds: ",time()-forward_started);flush(stdout)
        @test raw==snapshot
        @test H.validate(result.document;level=:surface).valid
        @test length(result.document["patches"])==64
        @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in result.document["patches"])<=256
        @test sum(length(p["triangles"]) for p in result.document["patches"])==6*2016
        @test length(result.correspondence["child_barycentric_rows"])==6*2016
        @test all(length(row["target"])==6 for row in result.correspondence["cell_children"])
        println("Actual rechart inverse started");flush(stdout)
        inverse_started=time()
        restored=H._marked_unrechart_geometry(result.document,result.inverse_specification)
        println("Actual rechart inverse seconds: ",time()-inverse_started);flush(stdout)
        @test restored.geometry==raw && restored.opposite_edges==pairs
        @test restored.construction_only && Set(keys(restored.geometry))==Set(["surface","patches"])
    end
else
    println("Actual over-limit rechart scope not executed: MARKED_RECHART_PROOF was not supplied")
end

function selective_features(source,result)
    corr=result.correspondence;doc=deepcopy(result.document)
    edges=Dict(Tuple(r["source"])=>r["target"] for r in corr["edge_images"])
    vertices=Dict(Tuple(r["source"])=>r["target"] for r in corr["vertex_images"])
    children=Dict(Tuple(r["source"])=>r["target"] for r in corr["cell_children"])
    doc["curves"]=deepcopy(source["curves"]);doc["crossings"]=deepcopy(source["crossings"])
    for c in doc["curves"]
        ds=M.Dart[]
        for d in M.parse_word(c["word"]);append!(ds,H._st_darts(edges[(d.face,d.u,d.v)]));end
        c["word"]=M.print_word(ds)
    end
    for x in doc["crossings"];x["at"]=deepcopy(vertices[Tuple(x["at"])]);end
    haskey(source,"metadata") && (doc["metadata"]=deepcopy(source["metadata"]))
    if haskey(source,"region_marks")
        doc["region_marks"]=deepcopy(source["region_marks"])
        for mark in doc["region_marks"];mark["seed"]=deepcopy(first(children[Tuple(mark["seed"])]));end
    end
    infinity=H._mm_infinity(source)
    infinity!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=deepcopy(first(children[Tuple(infinity)])))
    doc
end

# Source: the existing grid builder and directed meridian/longitude example.
# This is an independently generated oversized raw disk, not the D answer
# fixture. Existing 64/256 native bounds and all default tests remain unchanged.
@testset "Selective proper-arc identity partition and actual inverse" begin
    B=H.MarkedDiagramBuilders;source=B.grid_document(16)
    B.add_curve!(source,"a","alpha",true,[[(0,1),(16,1)]])
    B.add_curve!(source,"b","beta",true,[[(1,0),(1,16)]])
    B.declare_crossing!(source,"x",1,1,"a","b",1)
    raw=rawgeometry(source);snapshot=deepcopy(raw);pairs=H._marked_rechart_opposites(raw)
    result=H._marked_rechart_geometry(raw,pairs;specification=Dict("strategy"=>"proper_arc"))
    @test raw==snapshot
    @test H.validate(result.document;level=:surface).valid
    @test length(result.document["patches"])<=64
    @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in result.document["patches"])<=256
    @test !isempty(result.correspondence["proper_cuts"])
    @test sum(length(p["triangles"]) for p in result.document["patches"])<6*sum(length(p["triangles"]) for p in raw["patches"])
    @test result.inverse_specification["schema"]=="internal-triangle-partition/1"
    featured=selective_features(source,result)
    @test H.validate(featured;level=:surface).valid
    @test H.intersection_matrix(featured)==[1;;]
    restored=H._marked_unrechart_geometry(featured,jsonround(result.inverse_specification))
    @test restored.geometry==raw && restored.opposite_edges==pairs
    @test restored.construction_only && H._mm_hash(restored.correspondence)==H._mm_hash(result.correspondence)
    @test restored.forward_specification==Dict("strategy"=>"proper_arc")
    @test restored.fresh_forward.document==H._st_empty(featured)
    @test restored.fresh_forward.inverse_specification==result.inverse_specification
    @test restored.fresh_forward.correspondence===restored.correspondence
    @test restored.fresh_forward.correspondence!==result.correspondence
    bad=deepcopy(result.inverse_specification)
    bad["affine_children"][1]["weights"]=circshift(bad["affine_children"][1]["weights"],1)
    @test_throws ArgumentError H._marked_unrechart_geometry(featured,bad)
    bad=deepcopy(result.inverse_specification)
    a=deepcopy(bad["opposite_edges"][1]["target"])
    bad["opposite_edges"][1]["target"]=deepcopy(bad["opposite_edges"][2]["target"])
    bad["opposite_edges"][2]["target"]=a
    @test_throws ArgumentError H._marked_unrechart_geometry(featured,bad)
    @test_throws ArgumentError H._marked_rechart_geometry(raw,pairs;specification=Dict("strategy"=>"proper_arc","maximum_vertices"=>999))
    # Source: selective inverse consumes JSON, not a trusted in-memory receipt.
    # Malformed schema must fail before indexing/allocation, with field context.
    malformed=[
        ("triangle_ancestry",s->(s["triangle_ancestry"][1]["source_vertices"]=["x"])),
        ("affine_children",s->(s["affine_children"]=nothing)),
        ("fields",s->delete!(s,"affine_children")),
        ("triangle_count",s->(s["raw_charts"][1]["triangle_count"]=true)),
        ("raw_surface",s->(s["raw_surface"]=nothing)),
        ("rational",s->(s["affine_children"][1]["weights"][1][1]="not-a-rational")),
        ("opposite_edges",s->(s["opposite_edges"][1]["source"]=["x",false,"a","b"]))]
    for (field,change) in malformed
        bad=deepcopy(result.inverse_specification);change(bad)
        error=try
            H._marked_unrechart_geometry(featured,bad);nothing
        catch caught
            caught
        end
        @test error isa ArgumentError
        @test error isa Exception && occursin(field,sprint(showerror,error))
    end
end

# The lead supplies the unchanged actual raw D receipt. This bypasses only
# recomputation of its expensive accepted raw producer, never native capacity,
# features, exact child incidence, or the independently reconstructed inverse.
if haskey(ENV,"MARKED_PROPER_ARC_PROOF")
    @testset "Actual embedded-D raw proof selective native realization" begin
        bundle=JSON3.read(read(ENV["MARKED_PROPER_ARC_PROOF"],String),Dict{String,Any})
        @test bundle["construction_only"] && !bundle["native_validity_claimed"]
        source=bundle["complex"];raw=rawgeometry(source);snapshot=deepcopy(raw)
        @test sum(length(p["triangles"]) for p in raw["patches"])==8704
        @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in raw["patches"])==1055
        pairs=H._marked_rechart_opposites(raw)
        println("Actual selective rechart forward started");flush(stdout);started=time()
        result=H._marked_rechart_geometry(raw,pairs;specification=Dict("strategy"=>"proper_arc"))
        println("Actual selective rechart forward seconds: ",time()-started);flush(stdout)
        @test raw==snapshot
        @test length(result.document["patches"])<=64
        @test maximum(length(Set(v for t in p["triangles"] for v in t)) for p in result.document["patches"])<=256
        featured=selective_features(source,result)
        report=H.validate(featured;level=:surface)
        @test report.valid
        @test length(featured["region_marks"])==length(source["region_marks"])
        @test H._mm_infinity(featured)!==nothing
        println("Actual selective rechart inverse started");flush(stdout);started=time()
        restored=H._marked_unrechart_geometry(featured,jsonround(result.inverse_specification))
        println("Actual selective rechart inverse seconds: ",time()-started);flush(stdout)
        @test restored.geometry==raw && restored.opposite_edges==pairs
        @test H._mm_hash(restored.correspondence)==H._mm_hash(result.correspondence)
    end
else
    println("Actual D selective scope not executed: MARKED_PROPER_ARC_PROOF was not supplied")
end
end
