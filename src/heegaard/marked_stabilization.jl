"""Generated connected-sum stabilization on the native marked cellulation.

Sources: marked_subdivision.jl::marked_dual_star/marked_collapse_dual_star;
ograph_reverse.jl::_reverse_merge_charts!; vendor MarkedSurfaceDiagrams
fine_edge_occurrences/regions/vertex_dsu and builders.jl::grid_document.
The connected-sum and inverse argument is agent-derived. No original diagram
is retrieved as an inverse. The old marking is included off the removed disk.
"""

function _st_torus(;display_cut::Bool=false)
    B=MarkedDiagramBuilders
    n=display_cut ? 4 : 3;d=B.grid_document(n)
    B.add_curve!(d,"new_alpha","alpha",true,[[(0,1),(n,1)]])
    B.add_curve!(d,"new_beta","beta",true,[[(1,0),(1,n)]])
    B.declare_crossing!(d,"new_crossing",1,1,"new_alpha","new_beta",1)
    if display_cut
        # This is a DISTINCT constructed parallel meridian, not promotion of
        # alpha to a drawing cut. Ordered (beta,aux) has sign - for (north,east).
        B.add_curve!(d,"new_auxiliary","aux",true,[[(0,2),(n,2)]])
        B.declare_crossing!(d,"new_aux_crossing",1,2,"new_beta","new_auxiliary",-1)
    end
    d
end

_st_empty(d)=Dict{String,Any}("version"=>1,"surface"=>deepcopy(d["surface"]),
    "patches"=>deepcopy(d["patches"]),"curves"=>Any[],"crossings"=>Any[])
_st_darts(xs)=[MarkedSurfaceDiagrams.Dart(String(x[1]),String(x[2]),String(x[3])) for x in xs]
_st_records(ds)=[[d.face,d.u,d.v] for d in ds]

# Source: native oriented triangles and quotient vertex links. Empty means the
# CLOSED support is disjoint from every curve; a triangle interior alone is not
# enough. Selected D is protected as an entire mathematical complementary disk.
function _st_support(states,cell)
    M=MarkedSurfaceDiagrams;st=states.geometry;f,k=cell
    _mm_need(haskey(st.triangles,(f,k+1)),"unknown stabilization support triangle")
    t=st.triangles[(f,k+1)]
    roots=Set(M.root!(st.vertex_dsu,(f,v)) for v in t)
    _mm_need(length(roots)==3,"stabilization support does not have an embedded three-vertex closure")
    occupied=Set(M.root!(st.vertex_dsu,(d.face,v)) for (_,ds) in values(st.curves) for d in ds for v in (d.u,d.v))
    _mm_need(isdisjoint(roots,occupied),"stabilization support touches an existing curve")
    rs=M.regions(states.mathematics)
    for mark in get(states.document,"region_marks",[])
        r=only(r for r in rs if mark["seed"] in r["triangles"])
        _mm_need(!(collect(cell) in r["triangles"]),"stabilization support meets a protected selected region")
    end
    inf=_mm_infinity(states.document)
    _mm_need(inf===nothing || _mm_cell(inf)!=cell,"stabilization support contains the protected infinity triangle")
    t
end

function _st_star_edge_lookup(w)
    M=MarkedSurfaceDiagrams;lookup=Dict{Tuple{String,String,String},Any}()
    for r in w["correspondence"]["edges"]
        f,u,v=String.(r["source"]);ds=_st_darts(r["target"])
        lookup[(f,u,v)]=ds
        lookup[(f,v,u)]=[M.Dart(d.face,d.v,d.u) for d in reverse(ds)]
    end
    lookup
end

# All temporary patch vertices are prefixed by the existing merger, even when
# no merge is needed. Its triangle offset is ZERO based, as are exported cells.
function _st_merge_locations(rows,merges)
    loc=Dict{String,Any}()
    for row in rows
        f=row["temporary_face"]
        loc[f]=Dict("face"=>f,"triangle"=>0,"vertices"=>Dict(v=>f*"_"*v for v in ("a","b","c")))
    end
    for m in merges, r in values(loc)
        r["face"]==m["absorbed_chart"] || continue
        r["face"]=m["kept_chart"];r["triangle"]+=m["absorbed_triangle_index_offset"]
        r["vertices"]=Dict(v=>get(m["identified_vertices"],x,x) for (v,x) in r["vertices"])
    end
    loc
end

"""Cut the actual neck and inspect its two normalized surface components."""
function _st_neck(st,word,torus_cells,old_genus)
    M=MarkedSurfaceDiagrams
    # An independent neck-only document verifies a simple closed embedded loop.
    neck=_st_empty(st.doc)
    neck["curves"]=[Dict("id"=>"neck","role"=>"aux","closed"=>true,"word"=>word)]
    M.validate(neck)
    barriers=Set((d.face,M.edge(d.u,d.v)) for d in M.parse_word(word))
    rs=M.regions(st;barriers=barriers)
    _mm_need(length(rs)==2,"stabilization neck does not separate exactly two components")
    hit=findall(r->_mm_cells(r)==torus_cells,rs)
    _mm_need(length(hit)==1,"neck does not isolate the recorded torus cells")
    tor=rs[only(hit)];outside=rs[3-only(hit)]
    _mm_need(tor["genus"]==1 && tor["boundary_components"]==1 && tor["chi"]==-1,
        "generated handle is not an actual once-punctured torus")
    _mm_need(outside["genus"]==old_genus && outside["boundary_components"]==1,
        "neck outside does not recover the punctured source surface")
    Dict("torus"=>tor,"outside"=>outside)
end

function _st_region_transport(before,after,children,support;generated=Set{Tuple{String,Int}}())
    out=Any[];used=Set{Int}()
    for (i,r) in enumerate(before)
        source=_mm_cells(r)
        image=Set(c for t in source if t!=support for c in children[t])
        _mm_need(!isempty(image),"stabilization removes an entire source region")
        candidates=findall(s->issubset(image,_mm_cells(s)),after)
        _mm_need(length(candidates)==1,"stabilization source region is split or ambiguously transported")
        j=only(candidates);push!(used,j)
        support in source || _mm_need(image==_mm_cells(after[j]),"stabilization changes a region away from its support")
        push!(out,Dict("source_region"=>i,"target_region"=>j,"support_region"=>support in source,
            "source_triangles"=>[collect(c) for c in sort(collect(source))],
            "target_triangles"=>deepcopy(after[j]["triangles"])))
    end
    _mm_need(length(used)==length(out),"stabilization identifies distinct source regions")
    # An explicitly added auxiliary parallel may create new geometric regions
    # entirely INSIDE the generated block. Old region continuation remains
    # source bound; mathematical D is still an exact alpha-union-beta image.
    extra=setdiff(Set(eachindex(after)),used)
    _mm_need(all(issubset(_mm_cells(after[j]),generated) for j in extra),"stabilization creates an unexplained region outside its generated block")
    out
end

"""Check the actual cut complement of explicitly selected auxiliary circles.

Source: MarkedSurfaceDiagrams.regions/family_cut_report, using the named subset
of aux curves rather than all aux or the mathematical alpha family.
"""
function _st_display_cut(st,cuts)
    M=MarkedSurfaceDiagrams;g=st.surface_report["genus"]
    _mm_need(cuts isa AbstractVector && all(c isa AbstractString for c in cuts),"display cut_curves must be an explicit string vector")
    _mm_need(length(cuts)==g && length(Set(cuts))==g,"display auxiliary cut count must equal genus")
    _mm_need(all(haskey(st.curves,c) && st.curves[c][1]["role"]=="aux" && st.curves[c][1]["closed"] for c in cuts),"display cuts must be distinct closed auxiliary curves")
    ids=Set(cuts)
    _mm_need(all(length(intersect(Set(x["curves"]),ids))<2 for x in st.doc["crossings"]),"display auxiliary cuts intersect")
    barriers=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),ids))
    rs=M.regions(st;barriers=barriers)
    _mm_need(st.surface_report["boundary_components"]==0 && length(rs)==1 && only(rs)["genus"]==0 && only(rs)["boundary_components"]==2g,
        "specified auxiliary curves do not cut the surface to the required punctured sphere")
    Dict("curves"=>copy(cuts),"cut_genus"=>0,"boundary_components"=>2g,"components"=>1,"source_triangles"=>only(rs)["triangles"])
end

"""Attach an explicit torus in an empty disk, returning a generated inverse.

`support=[face,zero_based_triangle]` is required. Optional `alpha_id`, `beta_id`
and `crossing_id` name the newly constructed pair. Closed surfaces only; the
public dual-star resource limits are retained. A strictly unchanged Koda
joining disk is impossible for a disjoint new pair and is rejected explicitly.
Default behavior constructs no auxiliary display cut. With explicit
`extend_display_cuts=true`, an already supplied complete `metadata.display.cut_curves`
system is extended by a distinct constructed auxiliary meridian. Optional
`auxiliary_id` and `auxiliary_crossing_id` name its actual new data. The supplied
source cut and the resulting union are checked by actual cut topology.
"""
function marked_stabilize(input::AbstractDict,specification::AbstractDict)
    _md_guard_periodic_operation(input)
    # Source (app-root relative): ../topology-heegaard-cellulated/vendor/
    # marked_surface_diagrams/schema/diagram.schema.json, items.properties.id
    # of curves/crossings/region_marks; types.jl::_plain_value key normalization.
    # Reject unsupported scalar IDs BEFORE
    # kernel incidence checks or inverse maps could mix integer/string keys.
    raw=_plain_value(input)
    for field in ("curves","crossings","region_marks")
        records=get(raw,field,Any[])
        (records isa AbstractVector || records isa Tuple) || continue
        for r in records
            r isa AbstractDict && haskey(r,"id") || continue
            _mm_need(r["id"] isa AbstractString && !isempty(r["id"]),"stabilization source $field IDs must be nonempty strings")
        end
    end
    M=MarkedSurfaceDiagrams;original=_marked_document(raw);states=_marked_states(original);st=states.geometry
    spec=_plain_value(specification)
    _mm_need(Set(keys(spec)) ⊆ Set(["support","alpha_id","beta_id","crossing_id","extend_display_cuts","auxiliary_id","auxiliary_crossing_id"]),"unknown stabilization specification field")
    _mm_need(haskey(spec,"support"),"stabilization requires an explicit empty support triangle")
    # Source: marked_topology.jl's central infinity-triangle shape contract.
    # Reject before Int/String coercion: Bool and Float64 are not cell indices.
    cell=spec["support"]
    _mm_need(cell isa AbstractVector && length(cell)==2 && cell[1] isa AbstractString &&
        cell[2] isa Integer && !(cell[2] isa Bool) && cell[2]>=0,
        "stabilization support must be [face string, nonnegative integer triangle index]")
    _mm_need(all(!haskey(spec,k) || spec[k] isa AbstractString for k in ("alpha_id","beta_id","crossing_id","auxiliary_id","auxiliary_crossing_id")),
        "stabilization curve and crossing IDs must be strings")
    _mm_need(!haskey(spec,"extend_display_cuts") || spec["extend_display_cuts"] isa Bool,"extend_display_cuts must be an explicit Boolean")
    extend=get(spec,"extend_display_cuts",false)
    _mm_need(extend || !(haskey(spec,"auxiliary_id") || haskey(spec,"auxiliary_crossing_id")),"auxiliary IDs require explicit display-cut extension")
    _mm_need(st.surface_report["boundary_components"]==0,"stabilization currently requires a closed source surface")
    _mm_need(all(m["claim"]!="koda_joining" for m in get(original,"region_marks",[])),
        "disjoint stabilization cannot preserve an unchanged Koda joining disk; no claim is downgraded")
    support=_mm_cell(spec["support"]);cut=_st_support(states,support)
    degrees=Dict()
    for ((f,_),t) in st.triangles,v in t
        r=M.root!(st.vertex_dsu,(f,v));degrees[r]=get(degrees,r,0)+1
    end
    added_vertices=extend ? 13 : 6
    _mm_need(length(degrees)+added_vertices<=64,"stabilization needs $(length(degrees)+added_vertices) vertex stars; this bounded constructor retains the 64-patch limit")
    support_roots=Set(M.root!(st.vertex_dsu,(support[1],v)) for v in cut)
    _mm_need(all(3<=n+(r in support_roots ? 4 : 0)<=127 for (r,n) in degrees),
        "connected-sum vertex stars exceed the existing 3..127-corner disk contract")
    ids=(String(get(spec,"alpha_id","stabilization_alpha")),String(get(spec,"beta_id","stabilization_beta")),String(get(spec,"crossing_id","stabilization_crossing")))
    _mm_need(!isempty(ids[1]) && !isempty(ids[2]) && ids[1]!=ids[2] && all(!haskey(st.curves,c) for c in ids[1:2]),"new stabilization curve IDs must be distinct and unused")
    _mm_need(!isempty(ids[3]) && all(x["id"]!=ids[3] for x in original["crossings"]),"new stabilization crossing ID is already used")
    auxid=String(get(spec,"auxiliary_id","stabilization_cut"));auxcross=String(get(spec,"auxiliary_crossing_id","stabilization_cut_crossing"))
    generated_curves=extend ? [ids[1],ids[2],auxid] : [ids[1],ids[2]]
    generated_crossings=extend ? [ids[3],auxcross] : [ids[3]]
    _mm_need(length(unique(generated_curves))==length(generated_curves) && all(!isempty(c) && !haskey(st.curves,c) for c in generated_curves),"generated curve IDs must be distinct and unused")
    _mm_need(length(unique(generated_crossings))==length(generated_crossings) && all(!isempty(c) && all(x["id"]!=c for x in original["crossings"]) for c in generated_crossings),"generated crossing IDs must be distinct and unused")
    sourcecuts=nothing;source_cut_report=nothing
    if extend
        display=get(get(original,"metadata",Dict()),"display",nothing)
        _mm_need(display isa AbstractDict && haskey(display,"cut_curves"),"explicit display-cut extension requires source metadata.display.cut_curves; no cut is inferred")
        _mm_need(_mm_infinity(original)!==nothing,"explicit display-cut extension requires a source infinity triangle; no exterior region is inferred")
        sourcecuts=deepcopy(display["cut_curves"]);source_cut_report=_st_display_cut(st,sourcecuts)
    end
    # builders.jl::grid_document appends two triangles at each row-major (x,y):
    # first index=2*(y*n+x). Use (2,2) for n=3 or (2,3) for n=4. Their closed
    # quotient triangles avoid every constructed curve, checked again below.
    block=_st_torus(display_cut=extend);bs=M.validate(block);bcut=("P",extend ? 28 : 16)
    _st_support((geometry=bs,mathematics=bs,document=block),bcut)
    source_states=Dict("old"=>st,"torus"=>bs);removed=Dict("old"=>support,"torus"=>bcut)
    rows=Any[];by_origin=Dict();sides=Dict();patches=Any[];gluings=Any[]
    for kind in ("old","torus"),((f,i),t) in sort(collect(source_states[kind].triangles);by=first)
        (f,i-1)==removed[kind] && continue
        name=(kind=="old" ? "OldT" : "TorusT")*string(length(rows)+1)
        row=Dict("kind"=>kind,"source"=>[f,i-1],"source_vertices"=>collect(t),"temporary_face"=>name)
        push!(rows,row);by_origin[(kind,f,i)]=row
        push!(patches,Dict("face"=>name,"side_order"=>["e1","e2","e3"],
            "sides"=>Dict("e1"=>["a","b"],"e2"=>["b","c"],"e3"=>["c","a"]),"triangles"=>[["a","b","c"]]))
        for (j,(u,v)) in enumerate(M.tedges(t));sides[(kind,f,i,u,v)]=[name,"e$j"];end
    end
    pairs=Dict(kind=>M.fine_edge_occurrences(s) for (kind,s) in source_states)
    for kind in ("old","torus"),(o,p) in sort(collect(pairs[kind]);by=first)
        isless(o,p) || continue
        haskey(by_origin,(kind,o[1],o[2])) && haskey(by_origin,(kind,p[1],p[2])) || continue
        push!(gluings,[sides[(kind,o...)],sides[(kind,p...)]])
    end
    bt=bs.triangles[(bcut[1],bcut[2]+1)];boundarymap=Dict(cut[j]=>bt[[1,3,2][j]] for j in 1:3)
    for (u,v) in M.tedges(cut)
        o=pairs["old"][(support[1],support[2]+1,u,v)]
        a,b=boundarymap[v],boundarymap[u]
        p=pairs["torus"][(bcut[1],bcut[2]+1,a,b)]
        push!(gluings,[sides[("old",o...)],sides[("torus",p...)]])
    end
    atlas=Dict{String,Any}("version"=>1,"surface"=>Dict("marking_id"=>original["surface"]["marking_id"]*"/stabilized",
        "genus"=>st.surface_report["genus"]+1,"gluings"=>gluings,"boundary"=>Any[]),"patches"=>patches,"curves"=>Any[],"crossings"=>Any[])
    # The existing bare merger never receives source marks/curves/metadata.
    merges=_reverse_merge_charts!(atlas);loc=_st_merge_locations(rows,merges);merged=M.validate(atlas)
    star=marked_dual_star(atlas);doc=deepcopy(star.document);edge_lookup=_st_star_edge_lookup(star.witness)
    function image_dart(kind,d)
        src=source_states[kind];occ=first(src.incidence[d.face][M.edge(d.u,d.v)])
        o=(d.face,occ[1],occ[2],occ[3]);forward=(d.u,d.v)==o[3:4]
        if !haskey(by_origin,(kind,o[1],o[2]));o=pairs[kind][o];forward=!forward;end
        row=by_origin[(kind,o[1],o[2])];r=loc[row["temporary_face"]];t=row["source_vertices"]
        a=r["vertices"][["a","b","c"][only(findall(==(o[3]),t))]]
        b=r["vertices"][["a","b","c"][only(findall(==(o[4]),t))]]
        forward || ((a,b)=(b,a))
        edge_lookup[(r["face"],a,b)]
    end
    old_edges=Any[]
    for f in sort(collect(keys(st.patches))),e in sort(collect(keys(st.incidence[f]))),(u,v) in (e,reverse(e))
        push!(old_edges,Dict("source"=>[f,u,v],"target"=>_st_records(image_dart("old",M.Dart(f,u,v)))))
    end
    doc["curves"]=deepcopy(original["curves"])
    for c in doc["curves"];c["word"]=M.print_word(vcat([image_dart("old",d) for d in M.parse_word(c["word"])]...));end
    curve_ids=Dict("new_alpha"=>ids[1],"new_beta"=>ids[2],"new_auxiliary"=>auxid)
    crossing_ids=Dict("new_crossing"=>ids[3],"new_aux_crossing"=>auxcross)
    for c in block["curves"]
        push!(doc["curves"],Dict("id"=>curve_ids[c["id"]],"role"=>c["role"],"closed"=>true,
            "word"=>M.print_word(vcat([image_dart("torus",d) for d in M.parse_word(c["word"])]...))))
    end
    centers=Dict(Tuple(r["source_quotient_vertex"])=>r["target"] for r in star.witness["correspondence"]["vertices"])
    vertex_locations=Dict()
    for row in rows
        kind=row["kind"];r=loc[row["temporary_face"]];f=String(row["source"][1])
        for (j,v) in enumerate(row["source_vertices"])
            root=M.root!(source_states[kind].vertex_dsu,(f,v));vertex_locations[(kind,root)]=(r["face"],r["vertices"][["a","b","c"][j]])
        end
    end
    function image_vertex(kind,f,v)
        rt=M.root!(source_states[kind].vertex_dsu,(f,v));mf,mv=vertex_locations[(kind,rt)]
        deepcopy(centers[M.root!(merged.vertex_dsu,(mf,mv))])
    end
    doc["crossings"]=deepcopy(original["crossings"])
    for x in doc["crossings"];x["at"]=image_vertex("old",x["at"]...);end
    for x in block["crossings"]
        push!(doc["crossings"],Dict("id"=>crossing_ids[x["id"]],"at"=>image_vertex("torus",x["at"]...),
            "curves"=>[curve_ids[c] for c in x["curves"]],"sign"=>x["sign"]))
    end
    flags=Dict(_mm_cell(r["source"])=>[_mm_cell(c) for c in r["target"]] for r in star.witness["correspondence"]["triangles"])
    children=Dict{Tuple{String,Int},Any}();torus_cells=Set{Tuple{String,Int}}()
    for row in rows
        r=loc[row["temporary_face"]];row["merged_cell"]=[r["face"],r["triangle"]]
        row["merged_vertices"]=[r["vertices"][v] for v in ("a","b","c")]
        image=flags[_mm_cell(row["merged_cell"])]
        row["target_triangles"]=[collect(c) for c in image]
        row["kind"]=="old" ? (children[_mm_cell(row["source"])]=image) : union!(torus_cells,image)
    end
    haskey(original,"metadata") && (doc["metadata"]=deepcopy(original["metadata"]))
    if haskey(original,"region_marks")
        doc["region_marks"]=deepcopy(original["region_marks"])
        for m in doc["region_marks"];m["seed"]=collect(first(children[_mm_cell(m["seed"])]));end
    end
    inf=_mm_infinity(original)
    inf!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=collect(first(children[_mm_cell(inf)])))
    extend && (doc["metadata"]["display"]["cut_curves"]=vcat(sourcecuts,[auxid]))
    after=_marked_states(doc)
    mathregions=_st_region_transport(M.regions(states.mathematics),M.regions(after.mathematics),children,support)
    target_geo=M.regions(after.geometry)
    georegions=_st_region_transport(M.regions(st),target_geo,children,support;generated=extend ? torus_cells : Set{Tuple{String,Int}}())
    generated_regions=[deepcopy(r) for (j,r) in enumerate(target_geo) if !(j in Set(x["target_region"] for x in georegions))]
    target_cut_report=extend ? _st_display_cut(after.geometry,doc["metadata"]["display"]["cut_curves"]) : nothing
    neckword=M.print_word(vcat([image_dart("old",M.Dart(support[1],u,v)) for (u,v) in M.tedges(cut)]...))
    neckreport=_st_neck(after.geometry,neckword,torus_cells,st.surface_report["genus"])
    oldgenus=st.surface_report["genus"]
    _mm_need(after.geometry.surface_report["genus"]==oldgenus+1,"stabilization did not increase genus by one")
    prior=_marked_compression(states.mathematics);next=_marked_compression(after.mathematics)
    prior.closed && _mm_need(next.closed,"stabilization lost spherical compression boundaries")
    prior.alpha_complete && prior.beta_complete && _mm_need(next.alpha_complete && next.beta_complete,"stabilization lost complete cut systems")
    inverse=Dict{String,Any}("source_document"=>_mm_hash(original),"output_document"=>_mm_hash(doc),
        "forward_specification"=>deepcopy(spec),"source_mesh"=>Dict("surface"=>deepcopy(original["surface"]),"patches"=>deepcopy(original["patches"])),
        "support"=>collect(support),"support_vertices"=>collect(cut),"old_genus"=>oldgenus,"new_ids"=>collect(ids),
        "generated_curve_ids"=>generated_curves,"generated_crossing_ids"=>generated_crossings,
        "new_crossings"=>[deepcopy(x) for x in doc["crossings"] if x["id"] in generated_crossings],
        "source_display_cuts"=>sourcecuts,"display_cut_extension"=>extend,
        "triangle_ancestry"=>rows,"merges"=>merges,"star_witness"=>star.witness,"old_edge_images"=>old_edges,
        "old_crossings"=>[Dict("id"=>x["id"],"source_at"=>deepcopy(x["at"]),"target_at"=>image_vertex("old",x["at"]...)) for x in original["crossings"]],
        "word_formats"=>Dict(c["id"]=>c["word"] for c in original["curves"]),
        "source_region_seeds"=>_mm_inverse_seeds(original),"source_infinity"=>deepcopy(inf),
        "neck_word"=>neckword,"torus_triangles"=>[collect(c) for c in sort(collect(torus_cells))],
        "new_curve_words"=>Dict(c["id"]=>c["word"] for c in doc["curves"] if c["id"] in generated_curves))
    correspondence=Dict{String,Any}("cell_semantics"=>"identity inclusion of the source outside the removed disk, followed by recorded star flags",
        "removed_triangle"=>collect(support),"triangles"=>[Dict("source"=>collect(c),"target"=>[collect(x) for x in children[c]]) for c in sort(collect(keys(children)))],
        "edges"=>old_edges,"vertices"=>[Dict("source"=>[f,v],"target"=>image_vertex("old",f,v)) for (f,v) in sort(collect(keys(st.vertex_dsu.p)))],
        "curves"=>[Dict("source"=>c["id"],"target"=>c["id"],"source_word"=>c["word"],
            "target_word"=>only(x for x in doc["curves"] if x["id"]==c["id"])["word"]) for c in original["curves"]],
        "crossings"=>inverse["old_crossings"],"new_crossing"=>deepcopy(only(x for x in doc["crossings"] if x["id"]==ids[3])),
        "regions"=>mathregions,"geometric_regions"=>georegions,"neck"=>neckreport,"new_ids"=>collect(ids),
        "display_cut_extension_constructed"=>extend,"source_display_cut"=>source_cut_report,"target_display_cut"=>target_cut_report,
        "generated_geometric_regions"=>generated_regions,
        "selected_regions"=>[Dict("id"=>m["id"],"claim"=>m["claim"],"source_seed"=>m["seed"],
            "target_seed"=>only(x for x in doc["region_marks"] if x["id"]==m["id"])["seed"],"entire_cell_image_preserved"=>true) for m in get(original,"region_marks",[])],
        "infinity"=>(inf===nothing ? nothing : Dict("source_seed"=>inf,"target_seed"=>_mm_infinity(doc),
            "semantics"=>"retained source triangle child in its continued full geometric region")))
    w=_mm_witness("stabilization",original,doc,spec,inverse,correspondence);w["inverse_kind"]="generated_destabilization"
    (;document=doc,witness=w)
end

"""Remove only a verified generated handle by tracing, cutting and capping.

The current curves are collapsed through actual edge images. Geometric chart
provenance supplies local names; every surviving old triangle is reconstructed
from the current collapsed mesh, and the missing face is capped on the actual
neck. A fresh forward construction must reproduce target AND inverse data.
"""
function marked_destabilize(input::AbstractDict,specification::AbstractDict)
    _md_guard_periodic_operation(input)
    M=MarkedSurfaceDiagrams;current=_marked_document(input);spec=_plain_value(specification)
    _mm_need(_mm_hash(current)==spec["output_document"],"generated destabilization target has changed")
    st=_marked_states(current).geometry;ids=String.(spec["new_ids"])
    generated_curves=String.(spec["generated_curve_ids"]);generated_crossings=String.(spec["generated_crossing_ids"])
    torus=Set(_mm_cell(c) for c in spec["torus_triangles"])
    neck=_st_neck(st,spec["neck_word"],torus,Int(spec["old_genus"]))
    cap=String.(spec["support_vertices"]);support_face=String(spec["support"][1])
    edge_images=Dict(Tuple(r["source"])=>r["target"] for r in spec["old_edge_images"])
    reconstructed_neck=vcat([_st_darts(edge_images[(support_face,u,v)]) for (u,v) in M.tedges(cap)]...)
    _mm_need(M.parse_word(spec["neck_word"])==reconstructed_neck,"cap boundary does not trace the actual neck")
    for c in current["curves"]
        inside=c["id"] in generated_curves
        for d in M.parse_word(c["word"]),(i,_,_) in st.incidence[d.face][M.edge(d.u,d.v)]
            _mm_need(((d.face,i-1) in torus)==inside,"a curve crosses the generated neck or occupies the wrong summand")
        end
        c["id"] in generated_curves || continue
        _mm_need(c["word"]==spec["new_curve_words"][c["id"]],"generated torus meridian is altered")
    end
    x=only(x for x in current["crossings"] if x["id"]==ids[3])
    _mm_need(x["curves"]==ids[1:2] && x["sign"]==1,"generated torus crossing is altered")
    actual_new=[x for x in current["crossings"] if any(c->c in generated_curves,x["curves"])]
    _mm_need(_canonical_json(actual_new)==_canonical_json(spec["new_crossings"]),"new handle has altered or extra curve intersections")
    if spec["display_cut_extension"]
        _mm_need(_mm_infinity(current)!==nothing || spec["source_infinity"]===nothing,"generated display extension lost source infinity")
        targetcuts=current["metadata"]["display"]["cut_curves"]
        _mm_need(targetcuts==vcat(spec["source_display_cuts"],[last(generated_curves)]),"generated display-cut extension has changed")
        _st_display_cut(st,targetcuts)
    end
    # Source: marked_dual_star's two-dart old-edge image. Only original route
    # edges are inverse candidates; temporary seam occurrences are not answers.
    lookup=Dict{String,Any}()
    for r in spec["old_edge_images"]
        f,u,v=String.(r["source"])
        mesh=only(p for p in spec["source_mesh"]["patches"] if p["face"]==f)
        boundary=Set(M.edge(a,b) for s in mesh["side_order"] for (a,b) in zip(mesh["sides"][s][1:end-1],mesh["sides"][s][2:end]))
        M.edge(u,v) in boundary && continue
        key=_canonical_json(r["target"])
        _mm_need(!haskey(lookup,key),"old route inverse is ambiguous")
        lookup[key]=r["source"]
    end
    curves=deepcopy([c for c in current["curves"] if !(c["id"] in generated_curves)])
    for c in curves
        ds=M.parse_word(c["word"]);out=M.Dart[]
        _mm_need(iseven(length(ds)),"old curve is not a generated doubled-dart route")
        for j in 1:2:length(ds)
            key=_canonical_json(_st_records(ds[j:j+1]));_mm_need(haskey(lookup,key),"old curve enters the inserted handle or leaves its original edge image")
            f,u,v=lookup[key];push!(out,M.Dart(f,u,v))
        end
        word=spec["word_formats"][c["id"]]
        _mm_need(M.parse_word(word)==out,"lexical source word differs from reconstructed old curve")
        c["word"]=word
    end
    # Collapse an ACTUAL geometry-only target through the accepted generated
    # inverse. Metadata/curves are held above, not misread on temporary seams.
    empty=_st_empty(current)
    raw=marked_collapse_dual_star(empty,spec["star_witness"]["inverse_specification"]).document
    rawstate=M.validate(raw);meshes=Dict(p["face"]=>p for p in raw["patches"])
    patches=deepcopy(spec["source_mesh"]["patches"]);byface=Dict(p["face"]=>p for p in patches)
    for p in patches;p["triangles"]=Any[nothing for _ in p["triangles"]];end
    for row in spec["triangle_ancestry"]
        mf,k=row["merged_cell"];actual=meshes[mf]["triangles"][k+1]
        expected=String.(row["merged_vertices"])
        _mm_need(_mm_triangle(actual)==_mm_triangle(expected),"current merged triangle differs from its origin")
        row["kind"]=="old" || continue
        back=Dict(expected[j]=>row["source_vertices"][j] for j in 1:3)
        recovered=[back[v] for v in actual]
        _mm_need(_mm_triangle(recovered)==_mm_triangle(row["source_vertices"]),"old triangle orientation is altered")
        f,i=row["source"];byface[f]["triangles"][i+1]=recovered
    end
    f,k=spec["support"]
    # The neck report above verifies this three-edge boundary separates the
    # exact generated torus; orienting the cap as the deleted face closes it.
    byface[f]["triangles"][k+1]=cap
    _mm_need(all(t!==nothing for p in patches for t in p["triangles"]),"inverse cap leaves a missing old triangle")
    doc=deepcopy(current);doc["surface"]=deepcopy(spec["source_mesh"]["surface"]);doc["patches"]=patches;doc["curves"]=curves
    filter!(x->!(x["id"] in generated_crossings),doc["crossings"])
    locations=Dict(r["id"]=>r for r in spec["old_crossings"])
    for x in doc["crossings"]
        r=locations[x["id"]];_mm_need(x["at"]==r["target_at"],"old crossing leaves its transported vertex")
        x["at"]=deepcopy(r["source_at"])
    end
    for m in get(doc,"region_marks",[]);m["seed"]=deepcopy(spec["source_region_seeds"][m["id"]]);end
    spec["source_infinity"]!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=deepcopy(spec["source_infinity"]))
    spec["display_cut_extension"] && (doc["metadata"]["display"]["cut_curves"]=deepcopy(spec["source_display_cuts"]))
    _marked_states(doc)
    _mm_need(_mm_hash(doc)==spec["source_document"],"cut/cap and route reconstruction do not restore the source")
    generated=marked_stabilize(doc,spec["forward_specification"])
    _mm_need(generated.document==current && _canonical_json(generated.witness["inverse_specification"])==_canonical_json(spec),
        "forward regeneration does not reproduce this handle and witness")
    witness=_mm_witness("generated_destabilization",current,doc,spec,spec["forward_specification"],
        Dict("actual_neck"=>neck,"capped_triangle"=>[f,k],"reconstructed_curves"=>[c["id"] for c in curves]))
    witness["inverse_kind"]="stabilization"
    (;document=doc,witness)
end
