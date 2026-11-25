# cf: docs/specs/heegaard-snappy.md — inverse certificates and the unresolved PL-map boundary.
# cf: koda_interpreter.jl::koda_table_decode, koda_hc_from_boundary — reuse local tables.

"""
    hc_flow_certificate(input; orders=nothing, max_nodes::Int=100_000) -> NamedTuple

Return a checked Koda boundary-flow reconstruction of an oriented triangulation.

# Arguments
`input` and optional branching `orders` satisfy the
[triangulation-to-Koda domain](../specs/heegaard-snappy.md).
Encoding and branching checks propagate. The sphere-boundary, one-suture,
B/W disk, incidence, spine-Euler, curve-count and forward-word round-trip
conditions are checked; failures raise `HCTriangulationError` naming the stage.
Joining-disk and Koda reconstruction failures propagate.

# Returns
A named tuple with `encoded` o-data, `spine`, decoded `data`, raw curve data
`raw`, cellular `diagram`, selected joining `region`, `punctured` diagram and
`certificate` of boundary, suture and forward-round-trip checks.

# Guarantees
1. Successful return records the linked boundary-flow conditions and an
   exact forward o-data comparison up to metadata.
   [proved/agent: by construction]
2. The reconstructed object is a punctured cellular Koda diagram; its
   source identity is labelled. [proved/agent: by construction]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The result does not yet place curves on the fixed Lickorish reference
surface, nor independently certify a general PL homeomorphism between the
capped reconstruction and the input manifold. Equal o-data round trips do
not supply that missing geometric proof.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — inverse flow domain and unresolved PL-map boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> t = H.hc_triangulation(Topology.lens_space_link(2,1); signature=false); f = H.hc_flow_certificate(t; orders=t.branch_orders);

julia> (f.certificate["boundary_components"], f.certificate["suture_components"])
(1, 1)
```

Role: convert
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

"""
    hc_surface_cut_system(diagram) -> Dict{String,Any}

Return a checked tree–cotree polygon for a cellular Heegaard surface.

# Arguments
`diagram` must supply a nonempty connected ribbon graph and the surface data
of [the surface-cut domain](../specs/heegaard-snappy.md).
Spanning-tree, Euler, single-boundary, opposite-pair and polygon-vertex checks
raise stage-specific `HCTriangulationError`. Ribbon and surface errors propagate.
The caller supplies a valid cellular embedding for the geometric surface meaning.

# Returns
A dictionary with `genus`, primal and dual trees, residual edges, disk-boundary
darts, polygon word and edge ancestry, plus check and unresolved-transport fields.

# Guarantees
1. On success the retained polygon has the checked edge pairing and vertex
   structure of the input surface, with `2genus` residual edges.
   [proved/agent: by construction]
2. The original diagram is not modified. [proved/agent: by construction]

# Not guaranteed
The certificate leaves `pl_identification` absent and does not itself place
Lickorish reference curves or produce a resolved HC. The input's manifold
interpretation is not certified by this surface calculation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — tree–cotree cut and reference-overlay boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]); h = H.certify_heegaard(H.decorate(c,[[1]])); r = H.hc_surface_cut_system(h);

julia> (r["genus"], length(r["residual_edges"]), r["pl_identification"])
(1, 2, nothing)
```

Role: convert
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
