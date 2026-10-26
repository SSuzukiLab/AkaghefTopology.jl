"""Private identity rechart of an internal oriented triangle proof complex.

This is not an input format or another diagram type. `raw_geometry` contains
only `surface` and `patches`, and may exceed native per-patch limits. All
features stay with the caller until the actual bare atlas is native-valid.

Reuse boundary: marked_stabilization.jl::_st_merge_locations and
_st_star_edge_lookup; ograph_reverse.jl::_reverse_merge_charts!;
marked_subdivision.jl::marked_dual_star/marked_collapse_dual_star and its
correspondence.child_barycentric_rows. The immutable native kernel retains
64-chart/256-local-vertex checks. No cached-document inverse is permitted.
"""

_mr_json_occ(o)=[o[1],o[2]-1,o[3],o[4]]
function _mr_occ(o)
    _mm_need((o isa Tuple || o isa AbstractVector) && length(o)==4 && o[1] isa AbstractString &&
        o[2] isa Integer && !(o[2] isa Bool) && o[2]>=1 && o[3] isa AbstractString && o[4] isa AbstractString,
        "opposite edge must be (face,one-based triangle,u,v)")
    (String(o[1]),Int(o[2]),String(o[3]),String(o[4]))
end
function _mr_pairs_json(pairs)
    [Dict("source"=>_mr_json_occ(o),"target"=>_mr_json_occ(p)) for (o,p) in sort(collect(pairs);by=first) if isless(o,p)]
end
function _mr_pairs_live(rows)
    result=Dict()
    for row in rows
        a,b=row["source"],row["target"]
        o=_mr_occ((a[1],a[2]+1,a[3],a[4]));p=_mr_occ((b[1],b[2]+1,b[3],b[4]))
        _mm_need(!haskey(result,o) && !haskey(result,p) && o!=p,"duplicate inverse edge-pair occurrence")
        result[o]=p;result[p]=o
    end
    result
end

# Source: vendor patch_info, side_parameters and fine_edge_occurrences. The
# latter is duck typed; this INTERNAL NamedTuple is deliberately no DiagramState.
# Local disks and exact old seam parameters are checked before constructing a
# small native atlas that certifies the same global triangle quotient.
function _mr_proof(raw_geometry)
    M=MarkedSurfaceDiagrams;raw=_plain_value(raw_geometry)
    _mm_need(Set(keys(raw))==Set(["surface","patches"]),"construction geometry contains only surface and patches; hold features externally")
    surface=raw["surface"];M.check_keys(surface,("marking_id","genus","gluings","boundary"))
    _mm_need(surface["marking_id"] isa AbstractString && !isempty(surface["marking_id"]),"raw construction needs its explicit source marking")
    _mm_need(surface["genus"] isa Integer && !(surface["genus"] isa Bool) && surface["genus"]>=0,"raw construction genus must be a nonnegative integer")
    _mm_need(surface["boundary"] isa AbstractVector && isempty(surface["boundary"]),"triangle rechart currently requires a closed proof complex")
    _mm_need(raw["patches"] isa AbstractVector && !isempty(raw["patches"]),"raw construction needs disk patches")
    patches=Dict();incidence=Dict();triangles=Dict();sidepaths=Dict();parameters=Dict()
    for p in raw["patches"]
        M.check_keys(p,("face","side_order","sides","triangles"),("side_parameters",))
        f=p["face"];_mm_need(f isa AbstractString && !isempty(f) && !haskey(patches,f),"raw construction chart names must be distinct nonempty strings")
        _mm_need(all(v isa AbstractString && !isempty(v) for t in p["triangles"] for v in t),"raw construction vertex names must be nonempty strings")
        _,inc,_=M.patch_info(p);patches[f]=p;incidence[f]=inc
        supplied_parameters=get(p,"side_parameters",Dict())
        _mm_need(supplied_parameters isa AbstractDict && issubset(Set(keys(supplied_parameters)),Set(p["side_order"])),"raw construction has parameters for an unknown old side")
        for (i,t) in enumerate(p["triangles"]);triangles[(f,i)]=Tuple(t);end
        for sid in p["side_order"]
            path=p["sides"][sid];ts=M.side_parameters(p,sid)
            _mm_need(length(ts)==length(path) && first(ts)==0 && last(ts)==1 && all(ts[i]<ts[i+1] for i in 1:length(ts)-1),"raw seam parameters are not a strict oriented unit interval")
            sidepaths[(f,sid)]=path;parameters[(f,sid)]=ts
        end
    end
    paired=Set();seam=Dict()
    for ab in surface["gluings"]
        _mm_need(ab isa AbstractVector && length(ab)==2 && all(x->x isa AbstractVector && length(x)==2 && all(y->y isa AbstractString,x),ab),"raw old-chart gluing has invalid shape")
        a,b=Tuple(ab[1]),Tuple(ab[2])
        _mm_need(a!=b && haskey(sidepaths,a) && haskey(sidepaths,b) && !(a in paired) && !(b in paired),"raw old-chart side is unknown or paired more than once")
        push!(paired,a,b);va,vb=sidepaths[a],sidepaths[b]
        _mm_need(length(va)==length(vb) && parameters[a]==[1-t for t in reverse(parameters[b])],"raw old-chart seam subdivisions/parameters do not match by reversal")
        for i in 1:length(va)-1
            j=length(vb)-i;ka=(a[1],M.edge(va[i],va[i+1]));kb=(b[1],M.edge(vb[j],vb[j+1]))
            seam[ka]=kb;seam[kb]=ka
        end
    end
    _mm_need(paired==Set(keys(sidepaths)),"closed raw construction leaves an unpaired old-chart side")
    opposite=M.fine_edge_occurrences((;incidence,seam_partner=seam))
    all_occ=Set((f,i,u,v) for ((f,i),t) in triangles for (u,v) in M.tedges(t))
    _mm_need(Set(keys(opposite))==all_occ && all(get(opposite,p,nothing)==o && p!=o for (o,p) in opposite),"raw opposite-edge table is incomplete or not an involution")
    (;geometry=raw,patches,incidence,triangles,opposite_edges=opposite)
end

"""Derive the exact raw opposite involution without claiming native validity."""
_marked_rechart_opposites(raw_geometry::AbstractDict)=_mr_proof(raw_geometry).opposite_edges

function _mr_det(w)
    w[1][1]*(w[2][2]*w[3][3]-w[2][3]*w[3][2])-w[1][2]*(w[2][1]*w[3][3]-w[2][3]*w[3][1])+w[1][3]*(w[2][1]*w[3][2]-w[2][2]*w[3][1])
end
function _mr_weights(rows)
    M=MarkedSurfaceDiagrams
    _mm_need(rows isa AbstractVector && length(rows)==3 && all(r->r isa AbstractVector && length(r)==3,rows),"flag weights must be three barycentric rows")
    w=[[x isa Rational ? M.Rat(x) : M.qnum(x) for x in r] for r in rows]
    _mm_need(all(r->sum(r)==1 && all(>=(0),r),w) && _mr_det(w)>0,"flag barycentric map is folded or not a convex affine identity")
    w
end

function _mr_chart_provenance(raw)
    [begin
        result=deepcopy(p);result["triangle_count"]=length(pop!(result,"triangles"));result
    end for p in raw["patches"]]
end

"""Build a bounded native bare atlas and its actual identity transport.

`opposite_edges` is the COMPLETE vendor-style oriented involution on
`(face,one_based_triangle,u,v)` occurrences. Source/target weights and periodic
lifts remain caller-owned. Returned flag weights are rational barycentric
rows in actual public-triangle order over their raw parent, for composition.
"""
function _marked_rechart_geometry(raw_geometry::AbstractDict,opposite_edges::AbstractDict;
                                  specification::AbstractDict=Dict{String,Any}())
    _mm_need(isempty(specification),"private identity rechart accepts no replacement geometry or layout choices")
    first_atlas=_mr_first_atlas(raw_geometry,opposite_edges)
    star=marked_dual_star(first_atlas.atlas,Dict("merge_charts"=>true))
    _mr_compose_rechart(first_atlas,star.document,star.witness["correspondence"],star.witness["inverse_specification"])
end

# F: rebuild the raw triangle atlas and its first disk merges deterministically.
# Kept separate from S (marked_dual_star) so an inverse can compare this fresh
# atlas with the geometry it ACTUALLY collapsed before reusing that call's S.
function _mr_first_atlas(raw_geometry,opposite_edges)
    M=MarkedSurfaceDiagrams;proof=_mr_proof(raw_geometry);raw=proof.geometry
    supplied=Dict(_mr_occ(o)=>_mr_occ(p) for (o,p) in opposite_edges)
    _mm_need(supplied==proof.opposite_edges,"supplied opposite edges disagree with actual raw disks and old seams")
    rows=Any[];origin=Dict();side=Dict();patches=Any[];gluings=Any[]
    for ((f,i),t) in sort(collect(proof.triangles);by=first)
        name="RechartT$(length(rows)+1)"
        row=Dict("source"=>[f,i-1],"source_vertices"=>collect(t),"temporary_face"=>name)
        push!(rows,row);origin[(f,i)]=row
        # Unit parameter on each fine side is only its affine local coordinate;
        # the checked OLD macro-side parameters remain in raw chart provenance.
        push!(patches,Dict("face"=>name,"side_order"=>["e1","e2","e3"],
            "sides"=>Dict("e1"=>["a","b"],"e2"=>["b","c"],"e3"=>["c","a"]),"triangles"=>[["a","b","c"]]))
        for (j,(u,v)) in enumerate(M.tedges(t));side[(f,i,u,v)]=[name,"e$j"];end
    end
    for (o,p) in sort(collect(supplied);by=first);isless(o,p) && push!(gluings,[side[o],side[p]]);end
    atlas=Dict{String,Any}("version"=>1,"surface"=>Dict("marking_id"=>raw["surface"]["marking_id"]*"/recharted",
        "genus"=>raw["surface"]["genus"],"gluings"=>gluings,"boundary"=>Any[]),"patches"=>patches,"curves"=>Any[],"crossings"=>Any[])
    merger_map=Dict{String,Any}();merges=_reverse_merge_charts!(atlas;correspondence=merger_map)
    locations=_st_merge_locations(rows,merges);first_state=M.validate(atlas)
    (;proof,raw,supplied,rows,origin,atlas,merges,merger_map,locations,first_state)
end

# C: recompute the full outer identity correspondence from F and a freshly
# constructed S. This retains the canonical-edge, actual-row, affine-cover,
# capacity and full inverse-specification checks on both entry paths.
function _mr_compose_rechart(first_atlas,doc,star_correspondence,star_inverse)
    M=MarkedSurfaceDiagrams
    proof=first_atlas.proof;raw=first_atlas.raw;supplied=first_atlas.supplied;rows=first_atlas.rows
    origin=first_atlas.origin;locations=first_atlas.locations;first_state=first_atlas.first_state
    merges=first_atlas.merges;merger_map=first_atlas.merger_map
    final_state=M.validate(doc)
    _mm_need(length(doc["patches"])<=64 && all(length(Set(v for t in p["triangles"] for v in t))<=256 for p in doc["patches"]),"rechart output violates native chart capacity")
    star_edges=_st_star_edge_lookup(Dict("correspondence"=>star_correspondence))
    function mapped_occ(o)
        row=origin[o[1:2]];loc=locations[row["temporary_face"]];t=row["source_vertices"]
        letters=["a","b","c"]
        a=loc["vertices"][letters[only(findall(==(o[3]),t))]]
        b=loc["vertices"][letters[only(findall(==(o[4]),t))]]
        (loc["face"],loc["triangle"]+1,a,b)
    end
    canonical=sort([o for (o,p) in supplied if isless(o,p)])
    edge_records=Any[];edge_paths=Dict();seen_paths=Dict()
    for o in canonical
        q=mapped_occ(o);ds=star_edges[(q[1],q[3],q[4])];edge_paths[o]=ds
        for (sign,path) in ((1,ds),(-1,[M.Dart(d.face,d.v,d.u) for d in reverse(ds)]))
            key=_canonical_json(_st_records(path))
            _mm_need(!haskey(seen_paths,key),"canonical raw edge has an ambiguous public image")
            seen_paths[key]=(o,sign)
        end
        push!(edge_records,Dict("edge"=>_mr_json_occ(o),"target"=>_st_records(ds)))
    end
    occurrence_records=[Dict("source"=>_mr_json_occ(o),"edge"=>_mr_json_occ(min(o,p)),"orientation"=>(o==min(o,p) ? 1 : -1)) for (o,p) in sort(collect(supplied);by=first)]
    local_edge_images=Any[]
    for f in sort(collect(keys(proof.patches))),e in sort(collect(keys(proof.incidence[f])))
        i,u,v=first(proof.incidence[f][e]);o=(f,i,u,v);c=min(o,supplied[o]);positive=edge_paths[c]
        same=(o==c)==((u,v)==e)
        ds=same ? positive : [M.Dart(d.face,d.v,d.u) for d in reverse(positive)]
        push!(local_edge_images,Dict("source"=>[f,e...],"target"=>_st_records(ds)))
        push!(local_edge_images,Dict("source"=>[f,reverse(e)...],"target"=>_st_records([M.Dart(d.face,d.v,d.u) for d in reverse(ds)])))
    end
    centers=Dict(Tuple(v["source_quotient_vertex"])=>v["target"] for v in star_correspondence["vertices"])
    vertices=Dict();input_cells=Dict()
    for row in rows
        loc=locations[row["temporary_face"]];f,i=row["source"]
        row["merged_cell"]=[loc["face"],loc["triangle"]]
        row["merged_vertices"]=[loc["vertices"][x] for x in ("a","b","c")]
        actual=first_state.triangles[(loc["face"],loc["triangle"]+1)]
        _mm_need(_mm_triangle(actual)==_mm_triangle(row["merged_vertices"]),"first merger changed a raw triangle's orientation")
        input_cells[_mm_cell(row["merged_cell"])]=row
        for j in 1:3
            key=(String(f),String(row["source_vertices"][j]))
            point=centers[M.root!(first_state.vertex_dsu,(loc["face"],row["merged_vertices"][j]))]
            _mm_need(!haskey(vertices,key) || vertices[key]==point,"one raw vertex has inconsistent public occurrence images")
            vertices[key]=point
        end
    end
    children=Dict(_mm_cell(row["source"])=>Any[] for row in rows);affine=Any[]
    by_source=Dict(_mm_cell(row["source"])=>Any[] for row in rows)
    for child in star_correspondence["child_barycentric_rows"]
        parent=_mm_cell(child["source"]);row=input_cells[parent];target=_mm_cell(child["target"])
        actual=final_state.triangles[(target[1],target[2]+1)]
        _mm_need(Tuple(child["vertices"])==actual,"public flag row is not in actual triangle order")
        w=_mr_weights(child["weights"]);input_triangle=first_state.triangles[(parent[1],parent[2]+1)]
        perm=[only(findall(==(v),row["merged_vertices"])) for v in input_triangle]
        raw_weights=[[sum(w[k][j] for j in 1:3 if perm[j]==r) for r in 1:3] for k in 1:3]
        _mr_weights(raw_weights)
        record=Dict("source"=>copy(row["source"]),"target"=>collect(target),"vertices"=>collect(actual),"weights"=>[[M.qstr(x) for x in r] for r in raw_weights])
        push!(affine,record);push!(by_source[_mm_cell(row["source"])],record)
        push!(children[_mm_cell(row["source"])],collect(target))
    end
    for (cell,ids) in children
        maps=by_source[cell]
        _mm_need(length(ids)==6 && length(unique(Tuple.(ids)))==6 && sum(_mr_det(_mr_weights(r["weights"])) for r in maps)==1,
            "six public flags do not positively cover one raw triangle")
    end
    _mm_need(length(affine)==length(final_state.triangles),"public flag ancestry omits a native triangle")
    vertex_records=[Dict("source"=>collect(k),"target"=>deepcopy(vertices[k])) for k in sort(collect(keys(vertices)))]
    cell_records=[Dict("source"=>collect(c),"target"=>children[c]) for c in sort(collect(keys(children)))]
    corr=Dict{String,Any}("cell_semantics"=>"identity rechart of internal raw triangles; explicit affine flags preserve the caller's fixed target",
        "edge_images"=>local_edge_images,"vertex_images"=>vertex_records,"cell_children"=>cell_records,
        "child_barycentric_rows"=>affine,"canonical_edges"=>edge_records,"edge_occurrences"=>occurrence_records,
        "first_chart_merges"=>merges,"first_chart_correspondence"=>merger_map,
        "second_chart_merges"=>star_correspondence["chart_merges"])
    inverse=Dict{String,Any}("schema"=>"internal-triangle-rechart/1","raw_geometry_hash"=>_mm_hash(raw),
        "public_geometry_hash"=>_mm_hash(_st_empty(doc)),"raw_surface"=>deepcopy(raw["surface"]),
        "raw_charts"=>_mr_chart_provenance(raw),"triangle_ancestry"=>rows,"opposite_edges"=>_mr_pairs_json(supplied),
        "star_inverse"=>star_inverse,"affine_children"=>affine,"correspondence_hash"=>_mm_hash(corr))
    (;document=doc,correspondence=corr,inverse_specification=inverse)
end

"""Reconstruct actual INTERNAL raw rows and edge pairing from a native target.

The result is construction-only `geometry={surface,patches}` plus actual
opposite edges, never a native DiagramState. The caller must perform its own
map inverse on this reconstructed geometry before returning native data.
"""
function _marked_unrechart_geometry(current_native::AbstractDict,inverse_specification::AbstractDict)
    M=MarkedSurfaceDiagrams;current=_marked_document(current_native);spec=_plain_value(inverse_specification)
    _mm_need(get(spec,"schema","")=="internal-triangle-rechart/1","unknown private rechart inverse")
    bare=_st_empty(current)
    _mm_need(_mm_hash(bare)==spec["public_geometry_hash"],"public rechart geometry has changed")
    state=M.validate(bare)
    for row in spec["affine_children"]
        f,i=row["target"]
        _mm_need(Tuple(row["vertices"])==state.triangles[(f,i+1)],"inverse public affine row has wrong triangle incidence")
        _mr_weights(row["weights"])
    end
    collapse=marked_collapse_dual_star(bare,spec["star_inverse"])
    collapsed=collapse.document
    actual=M.validate(collapsed);actual_pairs=M.fine_edge_occurrences(actual)
    charts=deepcopy(spec["raw_charts"]);byface=Dict(p["face"]=>p for p in charts)
    for p in charts
        n=pop!(p,"triangle_count");p["triangles"]=Any[nothing for _ in 1:n]
    end
    occurrence_pullback=Dict()
    for row in spec["triangle_ancestry"]
        mf,k=row["merged_cell"];t=actual.triangles[(mf,k+1)];expected=String.(row["merged_vertices"])
        _mm_need(_mm_triangle(t)==_mm_triangle(expected),"actual collapsed triangle differs from its raw origin")
        backwards=Dict(expected[j]=>row["source_vertices"][j] for j in 1:3)
        recovered=[backwards[v] for v in t];f,i=row["source"]
        _mm_need(_mm_triangle(recovered)==_mm_triangle(row["source_vertices"]),"raw triangle orientation changed during reassembly")
        _mm_need(byface[f]["triangles"][i+1]===nothing,"duplicate raw triangle ancestry")
        byface[f]["triangles"][i+1]=recovered
        for (u,v) in M.tedges(t)
            q=(mf,k+1,u,v)
            _mm_need(!haskey(occurrence_pullback,q),"duplicate actual triangle edge in raw reassembly")
            occurrence_pullback[q]=(String(f),Int(i)+1,String(backwards[u]),String(backwards[v]))
        end
    end
    _mm_need(all(t!==nothing for p in charts for t in p["triangles"]),"raw reassembly loses an actual triangle")
    _mm_need(Set(keys(occurrence_pullback))==Set(keys(actual_pairs)),"raw ancestry omits actual collapsed triangle edges")
    reconstructed_pairs=Dict(occurrence_pullback[o]=>occurrence_pullback[p] for (o,p) in actual_pairs)
    _mm_need(reconstructed_pairs==_mr_pairs_live(spec["opposite_edges"]),"actual reassembled opposite-edge involution differs from recorded raw incidence")
    geometry=Dict{String,Any}("surface"=>deepcopy(spec["raw_surface"]),"patches"=>charts)
    proof=_mr_proof(geometry)
    _mm_need(proof.opposite_edges==reconstructed_pairs,"reassembled raw old seams disagree with actual current adjacency")
    _mm_need(_mm_hash(geometry)==spec["raw_geometry_hash"],"actual raw reassembly differs from the input geometric proof")
    first_atlas=_mr_first_atlas(geometry,reconstructed_pairs)
    _mm_need(first_atlas.atlas==collapse.document,"fresh raw triangle atlas differs from the actual collapsed native atlas")
    # Source: marked_subdivision.jl::marked_collapse_dual_star freshly calls
    # marked_dual_star(doc, spec["forward_specification"]) and requires BOTH
    # generated.document == before.document and canonical equality of the
    # ENTIRE generated.witness["inverse_specification"] with its spec before
    # returning. _mm_witness copies that certified spec into "specification".
    # Reuse only THIS successful call's fresh S correspondence; the certified
    # same-call specification supplies the unchanged outer star_inverse field.
    # Rebuilding F and recomputing C retain every outer row, incidence, affine,
    # inverse-specification and correspondence check without a second S call.
    fresh_star=collapse.witness["correspondence"]["forward_dual_star"]
    certified_star_inverse=collapse.witness["specification"]
    regenerated=_mr_compose_rechart(first_atlas,bare,fresh_star,certified_star_inverse)
    _mm_need(regenerated.document==bare && _canonical_json(regenerated.inverse_specification)==_canonical_json(spec) &&
        _mm_hash(regenerated.correspondence)==spec["correspondence_hash"],"forward identity rechart does not regenerate current geometry and inverse data")
    (;geometry,opposite_edges=reconstructed_pairs,correspondence=regenerated.correspondence,construction_only=true)
end
