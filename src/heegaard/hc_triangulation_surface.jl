# cf: docs/specs/heegaard-snappy.md — inverse certificates and the unresolved PL-map boundary.
# cf: koda_interpreter.jl::koda_table_decode, koda_hc_from_boundary — reuse local tables.

"""Certify the boundary pattern and reconstruct Koda curves using existing table primitives.

The returned punctured diagram is cellular Koda data, not resolved LC/HC.
Boundary counts are derived from the actual common B/W edges; the decoder's
reported B boundary_count is retained separately, never overwritten.
"""
function hc_flow_certificate(input;orders=nothing,max_nodes::Int=100_000)
    encoded=hc_triangulation_to_odata(input;orders=orders,max_nodes=max_nodes)
    spine=odata_to_spine(encoded.value); data=spine.metadata["table_spine"]
    boundary=data["boundary"]
    length(boundary)==1 && only(boundary)["chi"]==2 && only(boundary)["genus"]==0 ||
        throw(HCTriangulationError(:not_sphere_boundary,string(boundary)))
    length(data["e_cycles"])==1 || throw(HCTriangulationError(:not_one_suture,string(length(data["e_cycles"]))))
    cycle=only(data["e_cycles"])
    length(unique(cycle))==length(cycle) || throw(HCTriangulationError(:suture_not_simple,"repeated local boundary vertex"))
    bw=data["BW"]
    length(bw)==2 && sort([only(r["colors"]) for r in bw])==["B","W"] && all(r->r["chi"]==1,bw) ||
        throw(HCTriangulationError(:bw_not_two_disks,string(bw)))
    common=intersect(Set(bw[1]["edges"]),Set(bw[2]["edges"]))
    models=load_koda_local_tables().edge["models"]
    edge_pairs=Any[]; degrees=Dict{Tuple{Int,Int},Int}()
    for (eid,rank) in sort(collect(common))
        edge=data["edges"][eid+1]
        r,s=models[edge["model"]]["boundary_lift_map"][rank+1]
        a=(edge["source"][1],r); b=(edge["target"][1],s)
        push!(edge_pairs,minmax(a,b))
        degrees[a]=get(degrees,a,0)+1; degrees[b]=get(degrees,b,0)+1
    end
    expected=sort([minmax(cycle[i],cycle[mod1(i+1,length(cycle))]) for i in eachindex(cycle)])
    sort(edge_pairs)==expected && all(==(2),values(degrees)) && Set(keys(degrees))==Set(cycle) ||
        throw(HCTriangulationError(:suture_incidence_mismatch,"B/W common edges differ from the simple suture cycle"))
    data["spine_chi"]==1 || throw(HCTriangulationError(:spine_euler,string(data["spine_chi"])))
    raw=koda_hc_from_boundary(data,spine.vertex_signs)
    length(raw["alpha"])==length(raw["beta"])==raw["genus"] ||
        throw(HCTriangulationError(:koda_curve_count,"curve count does not equal genus"))
    diagram=koda_heegaard_from_hc(raw)
    region=_koda_select_joining_region(diagram,raw)
    disk=joining_disk(diagram,region)
    # This supported identity option avoids combinatorial canonical-label enumeration.
    punctured=make_punctured(diagram;joining_disk=disk,source_identity=:labelled)
    # Exercise the actual forward entrypoint. The separate diagnostic
    # koda_forward_gauss_from_hc currently drops its sign map (reported upstream).
    forward=hc_to_odata(diagram;joining_disk=disk)
    back=forward.correspondence["odata"]
    comparison=compare_odata_words(encoded.value,back)
    comparison["equal_up_to_metadata"] || throw(HCTriangulationError(:koda_roundtrip_mismatch,string(comparison)))
    certificate=Dict{String,Any}("boundary_components"=>1,"boundary_chi"=>2,
        "suture_components"=>1,"suture_simple"=>true,"suture_edge_count"=>length(common),
        "bw_chi"=>[r["chi"] for r in bw],"bw_boundary_count_derived"=>[1,1],
        "bw_boundary_count_reported"=>[r["boundary_count"] for r in bw],
        "bw_count_method"=>"common edge endpoints equal the one simple suture cycle; each color chi=1",
        "vertical_boundary"=>"annulus","koda_roundtrip"=>comparison,
        "roundtrip_entrypoint"=>"hc_to_odata(CellularHeegaardDiagram)",
        "source_identity"=>"labelled","source"=>"existing Koda tables and exact oriented face correspondence")
    (;encoded,spine,data,raw,diagram,region,punctured,certificate)
end

"""Produce a checked tree-cotree polygon, retaining edge ancestry.

This is a surface cut certificate. It does not yet provide the PL map needed
to pull back the standard Lickorish curves, so it is not an HC constructor.
"""
function hc_surface_cut_system(diagram)
    graph=ribbon_graph(diagram); faces=boundary_cycles(graph)
    isempty(graph.vertices) && throw(HCTriangulationError(:empty_surface_graph,"no intersection vertices"))
    primal=_KodaDSU(); tree=Int[]
    for v in graph.vertices; _koda_dsu_find!(primal,v); end
    for e in graph.edges
        if _koda_dsu_find!(primal,e.from)!=_koda_dsu_find!(primal,e.to)
            _koda_dsu_union!(primal,e.from,e.to); push!(tree,e.id)
        end
    end
    length(tree)==length(graph.vertices)-1 || throw(HCTriangulationError(:surface_graph_disconnected,"no spanning tree"))
    face_of=Dict(d=>i for (i,f) in enumerate(faces) for d in f)
    dual=_KodaDSU(); cotree=Int[]
    for i in eachindex(faces); _koda_dsu_find!(dual,i); end
    for e in graph.edges
        e.id in tree && continue
        a,b=face_of[2*e.id-1],face_of[2*e.id]
        if _koda_dsu_find!(dual,a)!=_koda_dsu_find!(dual,b)
            _koda_dsu_union!(dual,a,b); push!(cotree,e.id)
        end
    end
    length(cotree)==length(faces)-1 || throw(HCTriangulationError(:dual_tree_failure,"dual complement disconnected"))
    residual=setdiff([e.id for e in graph.edges],union(tree,cotree))
    genus=surface_genus(diagram)
    length(residual)==2genus || throw(HCTriangulationError(:surface_cut_euler,"residual edge count differs from 2g"))
    successor=Dict(i=>graph.next_dart[i] for i in eachindex(graph.darts))
    for e in cotree
        a,b=2*e-1,2*e
        # Merge the two oriented disk boundaries. Explicit walks also handle
        # a monogon, whose predecessor is the dart being removed itself.
        merged=Int[]
        for start in (a,b)
            d=successor[start]; steps=0
            while d!=start
                steps+=1
                steps>length(successor) && throw(HCTriangulationError(:surface_boundary_cycle,"nonclosing face"))
                push!(merged,d);d=successor[d]
            end
        end
        length(unique(merged))==length(merged) ||
            throw(HCTriangulationError(:dual_tree_failure,"cotree edge has both sides on one current disk"))
        for i in eachindex(merged); successor[merged[i]]=merged[mod1(i+1,length(merged))]; end
        delete!(successor,a);delete!(successor,b)
    end
    boundary=Int[]; first_dart=minimum(keys(successor)); d=first_dart
    while !(d in boundary)
        push!(boundary,d);d=successor[d]
    end
    d==first_dart && length(boundary)==length(successor) ||
        throw(HCTriangulationError(:surface_cut_not_disk,"dual tree did not produce one boundary"))
    word=[isodd(d) ? graph.darts[d].edge : -graph.darts[d].edge for d in boundary if graph.darts[d].edge in residual]
    all(e->count(==(e),word)==count(==(-e),word)==1,residual) ||
        throw(HCTriangulationError(:surface_cut_pairing,"residual edges must occur with opposite orientations"))
    corners=_KodaDSU(); n=length(word)
    for i in 1:n; _koda_dsu_find!(corners,i); end
    for e in residual
        i,j=findfirst(==(e),word),findfirst(==(-e),word)
        _koda_dsu_union!(corners,i,mod1(j+1,n));_koda_dsu_union!(corners,mod1(i+1,n),j)
    end
    length(Set(_koda_dsu_find!(corners,i) for i in 1:n))==1 ||
        throw(HCTriangulationError(:surface_polygon_vertices,"contracted polygon must have one vertex"))
    Dict{String,Any}("genus"=>genus,"primal_tree"=>tree,"dual_tree_edges"=>cotree,
        "residual_edges"=>residual,"disk_boundary_darts"=>boundary,"polygon_word"=>word,
        "edge_ancestry"=>[Dict("edge"=>e.id,"family"=>string(e.family),"curve"=>e.curve_id,
            "from"=>e.from,"to"=>e.to) for e in graph.edges],
        "checked"=>"spanning trees, one disk boundary, opposite edge pairs, one polygon vertex, 2g residual edges",
        "pl_identification"=>nothing,"stop"=>"polygon_normalization_and_reference_overlay_required")
end

export hc_each_branching, hc_flow_certificate, hc_surface_cut_system
