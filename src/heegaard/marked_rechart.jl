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
    # The opt-in selective strategy preserves already legal old-chart routes.
    # Its caller must check that raw external features avoid marking edges and
    # corners. Empty-spec uniform-star behavior and its inverse remain intact.
    _plain_value(specification)==Dict("strategy"=>"proper_arc") && return _mr_proper_arc_geometry(raw_geometry,opposite_edges)
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
    get(spec,"schema",nothing)=="internal-triangle-partition/1" && return _mr_unpartition_geometry(current,spec)
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
# Source: actual 8704-triangle D proof and the independently native-validated
# proper-arc probe (44 charts, maximum 203 vertices). The local split is the
# existing exact barycentric midpoint construction: an oriented (a,b,c), with
# m=(a+b)/2,n=(c+a)/2, becomes (a,m,n),(m,b,c),(m,c,n). It retains every old
# vertex and is an identity subdivision, not a change of the fixed target.
_mrp_row(parent,vertices,weights)=(parent=parent,vertices=String.(vertices),weights=weights)
_mrp_midweights(a,b)=[(x+y)/2 for (x,y) in zip(a,b)]
function _mrp_incidence(rows)
    M=MarkedSurfaceDiagrams;inc=Dict{Tuple{String,String},Vector{Tuple{Int,String,String}}}()
    for (i,row) in enumerate(rows),(u,v) in M.tedges(row.vertices)
        push!(get!(inc,M.edge(u,v),Tuple{Int,String,String}[]),(i,u,v))
    end
    inc
end
function _mrp_mesh(rows)
    M=MarkedSurfaceDiagrams;inc=_mrp_incidence(rows);adj=Dict(i=>Int[] for i in eachindex(rows));boundary=Dict{String,String}()
    for e in sort(collect(keys(inc)))
        os=inc[e];_mm_need(length(os) in (1,2),"proper-arc cell has nonmanifold incidence")
        if length(os)==2
            _mm_need(os[1][2:3]==reverse(os[2][2:3]),"proper-arc incidence reverses surface orientation")
            push!(adj[os[1][1]],os[2][1]);push!(adj[os[2][1]],os[1][1])
        else
            _,u,v=only(os);_mm_need(!haskey(boundary,u),"proper-arc boundary branches");boundary[u]=v
        end
    end
    for ns in values(adj);sort!(ns);end
    _mm_need(!isempty(boundary) && length(M.components(eachindex(rows),adj))==1,"proper-arc disk is disconnected or has no boundary")
    cycle=String[minimum(keys(boundary))]
    while boundary[last(cycle)]!=first(cycle)
        push!(cycle,boundary[last(cycle)]);_mm_need(length(cycle)<=length(boundary),"proper-arc boundary does not close")
    end
    vertices=Set(v for row in rows for v in row.vertices)
    _mm_need(length(cycle)==length(boundary) && length(vertices)-length(inc)+length(rows)==1,"proper-arc child is not a disk")
    (;inc,adj,cycle,vertices)
end
function _mrp_check_disk(rows)
    # Reuse the native disk/link validator on the actual triangle rows. Three
    # temporary boundary names are geometry-only; final sides are assembled
    # from maximal paired runs, never from these temporary corner choices.
    M=MarkedSurfaceDiagrams;s=_mrp_mesh(rows);v=s.cycle
    p=Dict("face"=>"proper_arc_disk","side_order"=>["a","b","c"],
        "sides"=>Dict("a"=>v[1:2],"b"=>v[2:3],"c"=>vcat(v[3:end],[v[1]])),
        "triangles"=>[r.vertices for r in rows])
    M.patch_info(p);s
end
function _mrp_paths(adj,start)
    # Source: Topology._shortest_path in virtual_planarization.jl. The same BFS
    # is adjacency-indexed and shared by all sampled endpoints of one source;
    # repeatedly scanning the full edge table would be quadratic per query.
    parent=Dict(start=>0);queue=[start];head=1
    while head<=length(queue)
        v=queue[head];head+=1
        for w in adj[v]
            haskey(parent,w) && continue
            parent[w]=v;push!(queue,w)
        end
    end
    parent
end
function _mrp_path(parent,start,finish)
    _mm_need(haskey(parent,finish),"proper-arc endpoint is disconnected")
    path=[finish]
    while last(path)!=start;push!(path,parent[last(path)]);end
    reverse(path)
end
function _mrp_midpoint!(state,f,u,v)
    M=MarkedSurfaceDiagrams;key=(f,M.edge(u,v))
    if !haskey(state.midpoints,key)
        state.serial[]+=1;state.midpoints[key]=state.prefix*string(state.serial[])
    end
    state.midpoints[key]
end
function _mrp_boundary_split!(state,f,u,v,m)
    M=MarkedSurfaceDiagrams;key=(f,M.edge(u,v));haskey(state.labels,key) || return
    lab,qs=state.labels[key];q=(qs[u]+qs[v])/2
    state.labels[(f,M.edge(u,m))]=(lab,Dict(u=>qs[u],m=>q))
    state.labels[(f,M.edge(m,v))]=(lab,Dict(m=>q,v=>qs[v]))
end
function _mrp_split_boundary!(state,p,e)
    rows=p.rows;os=_mrp_incidence(rows)[e];_mm_need(length(os)==1,"mirrored refinement is not on a boundary")
    i,u,v=only(os);row=rows[i];w=only(setdiff(row.vertices,collect(e)));bw=Dict(zip(row.vertices,row.weights))
    m=_mrp_midpoint!(state,p.origin,u,v);bm=_mrp_midweights(bw[u],bw[v])
    splice!(rows,i:i,[_mrp_row(row.parent,[u,m,w],[bw[u],bm,bw[w]]),_mrp_row(row.parent,[m,v,w],[bm,bw[v],bw[w]])])
    _mrp_boundary_split!(state,p.origin,u,v,m)
end
function _mrp_synchronize!(state)
    # Midpoints are exact dyadic refinements on the existing rational interval.
    # Mirror old seams AND earlier proper cuts; final side-run segmentation
    # subsequently promotes matching endpoints to corners on both sides.
    count=0
    while true
        endpoints=Dict();entries=Any[]
        for name in sort(collect(keys(state.patches)))
            p=state.patches[name]
            for (e,os) in sort(collect(_mrp_incidence(p.rows));by=first)
                length(os)==1 || continue
                lab,qs=state.labels[(p.origin,e)];union!(get!(endpoints,lab,Set()),values(qs))
                push!(entries,(name,e,lab,sort(collect(values(qs)))))
            end
        end
        selected=findfirst(x->any(x[4][1]<q<x[4][2] for q in endpoints[x[3]]),entries)
        selected===nothing && return count
        name,e,_,_=entries[selected];_mrp_split_boundary!(state,state.patches[name],e);count+=1
    end
end
function _mrp_candidate(state,p,path,e1,e2;temporary=true)
    M=MarkedSurfaceDiagrams;rows=p.rows;crossed=[e1]
    for (a,b) in zip(path,path[2:end])
        common=intersect(Set(rows[a].vertices),Set(rows[b].vertices))
        _mm_need(length(common)==2,"dual path does not cross one triangle edge")
        push!(crossed,Tuple(sort(collect(common))))
    end
    push!(crossed,e2);_mm_need(length(Set(crossed))==length(crossed),"proper arc repeats an edge")
    mids=Dict(e=>(temporary ? state.prefix*"candidate_"*string(i) : _mrp_midpoint!(state,p.origin,e...)) for (i,e) in enumerate(crossed))
    replacements=Dict();seam=Tuple{String,String}[]
    for (j,i) in enumerate(path)
        row=rows[i];a=only(intersect(Set(crossed[j]),Set(crossed[j+1])));k=only(findall(==(a),row.vertices))
        a,b,c=[row.vertices[mod1(k+r,3)] for r in 0:2];m=mids[M.edge(a,b)];n=mids[M.edge(c,a)]
        bw=Dict(zip(row.vertices,row.weights));bm=_mrp_midweights(bw[a],bw[b]);bn=_mrp_midweights(bw[c],bw[a])
        replacements[i]=[_mrp_row(row.parent,[a,m,n],[bw[a],bm,bn]),
            _mrp_row(row.parent,[m,b,c],[bm,bw[b],bw[c]]),_mrp_row(row.parent,[m,c,n],[bm,bw[c],bn])]
        push!(seam,(mids[crossed[j]],mids[crossed[j+1]]))
    end
    refined=Any[r for (i,row) in enumerate(rows) for r in get(replacements,i,[row])]
    inc=_mrp_incidence(refined);blocked=Set(M.edge(u,v) for (u,v) in seam);adj=Dict(i=>Int[] for i in eachindex(refined))
    for (e,os) in inc
        length(os)==2 && !(e in blocked) || continue
        push!(adj[os[1][1]],os[2][1]);push!(adj[os[2][1]],os[1][1])
    end
    components=M.components(eachindex(refined),adj);_mm_need(length(components)==2,"proper arc did not cut the disk into two disks")
    children=[refined[sort(collect(c))] for c in components]
    counts=[length(Set(v for r in ch for v in r.vertices)) for ch in children]
    (;children,counts,seam,mids)
end
function _mrp_choose(state,p)
    M=MarkedSurfaceDiagrams;s=_mrp_check_disk(p.rows);boundary=collect(zip(s.cycle,circshift(s.cycle,-1)));n=length(boundary)
    sample_count=min(24,n);samples=sort(unique(mod(round(Int,i*n/sample_count),n)+1 for i in 0:sample_count-1))
    best=nothing;bestscore=nothing
    for (position,i) in enumerate(samples)
        e1=M.edge(boundary[i]...);start=only(s.inc[e1])[1];parent=_mrp_paths(s.adj,start)
        for j in samples[position+1:end]
            e2=M.edge(boundary[j]...);finish=only(s.inc[e2])[1]
            start==finish && continue
            state.labels[(p.origin,e1)][1]==state.labels[(p.origin,e2)][1] && continue
            path=_mrp_path(parent,start,finish);candidate=_mrp_candidate(state,p,path,e1,e2)
            maximum(candidate.counts)<length(s.vertices) || continue
            score=(maximum(candidate.counts),length(path),i,j)
            if bestscore===nothing || isless(score,bestscore);bestscore=score;best=(path,e1,e2);end
        end
    end
    _mm_need(best!==nothing,"bounded proper-arc search found no reducing disk cut")
    best
end
function _mrp_partition(proof)
    M=MarkedSurfaceDiagrams;Q=M.Rat;patches=Dict{String,Any}();labels=Dict();midpoints=Dict();serial=Ref(0)
    prefix="partition_mid_";oldnames=Set(v for t in values(proof.triangles) for v in t)
    while any(startswith(v,prefix) for v in oldnames);prefix="p_"*prefix;end
    for f in sort(collect(keys(proof.patches)))
        rows=Any[_mrp_row((f,i-1),t,[[Q(j==k) for k in 1:3] for j in 1:3]) for (i,t) in enumerate(proof.patches[f]["triangles"])]
        patches[f]=(origin=f,rows=rows)
    end
    for (number,(a,b)) in enumerate(proof.geometry["surface"]["gluings"])
        for (side,sign) in ((a,1),(b,-1))
            p=proof.patches[side[1]];vs=p["sides"][side[2]];qs=M.side_parameters(p,side[2]);sign==-1 && (qs=1 .-qs)
            for j in 1:length(vs)-1
                labels[(side[1],M.edge(vs[j],vs[j+1]))]=("old$number",Dict(vs[j]=>qs[j],vs[j+1]=>qs[j+1]))
            end
        end
    end
    state=(;patches,labels,midpoints,serial,prefix);trace=Any[]
    while true
        sizes=Dict(name=>length(Set(v for r in p.rows for v in r.vertices)) for (name,p) in patches)
        bad=sort([name for name in keys(patches) if sizes[name]>224];by=n->(sizes[n],n))
        isempty(bad) && break
        name=last(bad);p=patches[name];path,e1,e2=_mrp_choose(state,p)
        result=_mrp_candidate(state,p,path,e1,e2;temporary=false)
        for (e,m) in result.mids;_mrp_boundary_split!(state,p.origin,e...,m);end
        sid="cut$(length(trace))"
        for (j,(u,v)) in enumerate(result.seam)
            labels[(p.origin,M.edge(u,v))]=(sid,Dict(u=>Q(j-1)/length(result.seam),v=>Q(j)/length(result.seam)))
        end
        delete!(patches,name)
        childnames=[name*"_"*string(k-1) for k in eachindex(result.children)]
        _mm_need(all(!haskey(patches,n) for n in childnames),"proper-arc child chart name collides")
        for (child,rows) in zip(childnames,result.children);_mrp_check_disk(rows);patches[child]=(origin=p.origin,rows=rows);end
        mirrored=_mrp_synchronize!(state)
        push!(trace,Dict("source"=>name,"path_triangles"=>path .-1,"entry_edge"=>collect(e1),"exit_edge"=>collect(e2),
            "children"=>childnames,"mirrored_triangle_splits"=>mirrored))
        _mm_need(length(patches)<=64,"selective disk partition exceeds native chart count")
    end
    (;state,trace)
end

function _mrp_assemble(partition,raw)
    M=MarkedSurfaceDiagrams;state=partition.state;patches=state.patches;boundaries=Dict();meshes=Dict()
    for name in sort(collect(keys(patches)))
        p=patches[name];s=_mrp_mesh(p.rows);meshes[name]=s
        for (e,os) in s.inc
            length(os)==1 || continue
            i,u,v=only(os);lab,qs=state.labels[(p.origin,e)];lo,hi=minmax(qs[u],qs[v])
            push!(get!(boundaries,(lab,lo,hi),Any[]),(name,i,u,v,qs[u],qs[v]))
        end
    end
    peer=Dict()
    for pair in values(boundaries)
        _mm_need(length(pair)==2,"selective boundary interval has no unique opposite")
        a,b=pair;_mm_need(a[5:6]==reverse(b[5:6]),"selective boundary interval orientation differs")
        peer[(a[1],M.edge(a[3],a[4]))]=b[1];peer[(b[1],M.edge(b[3],b[4]))]=a[1]
    end
    output=Any[];sidegroups=Dict();affine=Any[];children=Dict();origins=Dict()
    for name in sort(collect(keys(patches)))
        p=patches[name];s=meshes[name];boundary=collect(zip(s.cycle,circshift(s.cycle,-1)));n=length(boundary)
        kinds=[(state.labels[(p.origin,M.edge(u,v))][1],peer[(name,M.edge(u,v))]) for (u,v) in boundary]
        starts=[i for i in 1:n if kinds[i]!=kinds[mod1(i-1,n)]]
        _mm_need(length(starts)>=3,"selective child has fewer than three maximal named sides")
        order=String[];sides=Dict();parameters=Dict()
        for (k,start) in enumerate(starts)
            stop=starts[mod1(k+1,length(starts))];indices=Int[];i=start
            while true;push!(indices,i);i=mod1(i+1,n);i==stop && break;end
            vs=vcat([boundary[i][1] for i in indices],[boundary[last(indices)][2]])
            qs=[state.labels[(p.origin,M.edge(boundary[i]...))][2][boundary[i][1]] for i in indices]
            push!(qs,state.labels[(p.origin,M.edge(boundary[last(indices)]...))][2][last(vs)])
            _mm_need(all(kinds[i]==kinds[start] for i in indices),"selective side crosses a corner or opposite-chart change")
            ps=[(q-first(qs))/(last(qs)-first(qs)) for q in qs];sid="e$k";push!(order,sid);sides[sid]=vs;parameters[sid]=M.qstr.(ps)
            _mm_need(first(ps)==0 && last(ps)==1 && all(ps[i]<ps[i+1] for i in 1:length(ps)-1),"selective side parameters are not strict")
            lo,hi=minmax(first(qs),last(qs));push!(get!(sidegroups,(kinds[start][1],lo,hi),Any[]),(name,sid,vs,ps))
        end
        push!(output,Dict("face"=>name,"side_order"=>order,"sides"=>sides,"side_parameters"=>parameters,"triangles"=>[r.vertices for r in p.rows]))
        origins[name]=p.origin
        for (i,row) in enumerate(p.rows)
            record=Dict("source"=>collect(row.parent),"target"=>[name,i-1],"vertices"=>copy(row.vertices),"weights"=>[M.qstr.(w) for w in row.weights])
            push!(affine,record);push!(get!(children,row.parent,Any[]),[name,i-1])
        end
    end
    gluings=Any[]
    for (_,pair) in sort(collect(sidegroups);by=first)
        _mm_need(length(pair)==2,"selective named side has no unique opposite")
        a,b=pair;_mm_need(a[4]==[1-q for q in reverse(b[4])],"selective named-side subdivisions disagree")
        push!(gluings,[collect(a[1:2]),collect(b[1:2])])
    end
    surface=deepcopy(raw["surface"]);surface["marking_id"]*="/proper-arc";surface["gluings"]=gluings
    doc=Dict{String,Any}("version"=>1,"surface"=>surface,"patches"=>output,"curves"=>Any[],"crossings"=>Any[])
    native=M.validate(doc)
    (;document=doc,native,affine,children,origins,meshes)
end

function _mrp_expand(midpoints,f,u,v)
    M=MarkedSurfaceDiagrams;m=get(midpoints,(f,M.edge(u,v)),nothing)
    m===nothing && return [(u,v)]
    vcat(_mrp_expand(midpoints,f,u,m),_mrp_expand(midpoints,f,m,v))
end
function _mrp_correspondence(proof,partition,assembled)
    M=MarkedSurfaceDiagrams;state=partition.state;native=assembled.native;owners=Dict();vertices=Dict()
    oldvertices=Dict(f=>Set(v for t in p["triangles"] for v in t) for (f,p) in proof.patches)
    for (name,p) in state.patches
        mesh=assembled.meshes[name]
        for e in keys(mesh.inc);push!(get!(owners,(p.origin,e),String[]),name);end
        for v in mesh.vertices
            haskey(proof.incidence,p.origin) || continue
            v in oldvertices[p.origin] || continue
            key=(p.origin,v);point=collect(M.root!(native.vertex_dsu,(name,v)))
            _mm_need(!haskey(vertices,key) || vertices[key]==point,"old vertex has inconsistent selective images");vertices[key]=point
        end
    end
    function image(o)
        f,_,u,v=o;path=Any[]
        for (a,b) in _mrp_expand(state.midpoints,f,u,v)
            names=owners[(f,M.edge(a,b))];_mm_need(length(names)==1,"old edge has an ambiguous selective local image")
            push!(path,[only(names),a,b])
        end
        path
    end
    canonical=sort([o for (o,p) in proof.opposite_edges if isless(o,p)]);paths=Dict(o=>image(o) for o in canonical)
    reversepath(xs)=[[x[1],x[3],x[2]] for x in reverse(xs)]
    edge_records=[Dict("edge"=>_mr_json_occ(o),"target"=>paths[o]) for o in canonical]
    local_edges=Any[];seen=Set()
    for o in canonical,path in (paths[o],reversepath(paths[o]))
        key=_canonical_json(path);_mm_need(!(key in seen),"selective canonical raw edges share one directed path");push!(seen,key)
    end
    for f in sort(collect(keys(proof.patches))),e in sort(collect(keys(proof.incidence[f])))
        i,u,v=first(proof.incidence[f][e]);o=(f,i,u,v);c=min(o,proof.opposite_edges[o]);same=(o==c)==((u,v)==e)
        path=same ? paths[c] : reversepath(paths[c])
        push!(local_edges,Dict("source"=>[f,e...],"target"=>path));push!(local_edges,Dict("source"=>[f,reverse(e)...],"target"=>reversepath(path)))
    end
    occurrences=[Dict("source"=>_mr_json_occ(o),"edge"=>_mr_json_occ(min(o,p)),"orientation"=>(o==min(o,p) ? 1 : -1)) for (o,p) in sort(collect(proof.opposite_edges);by=first)]
    Dict{String,Any}("cell_semantics"=>"selective proper-arc identity partition; positive exact children retain each original raw parent",
        "route_domain"=>"old chart-interior curves with legal noncorner gates; old marking-edge images are geometric incidence only",
        "edge_images"=>local_edges,"vertex_images"=>[Dict("source"=>collect(k),"target"=>vertices[k]) for k in sort(collect(keys(vertices)))],
        "cell_children"=>[Dict("source"=>collect(k),"target"=>assembled.children[k]) for k in sort(collect(keys(assembled.children)))],
        "child_barycentric_rows"=>assembled.affine,"canonical_edges"=>edge_records,"edge_occurrences"=>occurrences,
        "proper_cuts"=>partition.trace)
end

# Positive areas alone do not prove subdivision. Read actual target incidence,
# cancel every parent-internal edge, require each complete oriented unit edge,
# and pull ACTUAL opposite occurrences back to the full original involution.
# This is used on construction AND inverse, not just as a saved certificate.
function _mrp_check_cells(native,affine,ancestry,expected_pairs)
    M=MarkedSurfaceDiagrams;bytarget=Dict();byparent=Dict();names=Dict()
    for row in ancestry
        cell=_mm_cell(row["source"]);_mm_need(!haskey(names,cell),"duplicate selective parent ancestry")
        names[cell]=String.(row["source_vertices"])
    end
    for row in affine
        target=_mm_cell(row["target"]);parent=_mm_cell(row["source"])
        _mm_need(haskey(names,parent) && haskey(native.triangles,(target[1],target[2]+1)),"selective affine row has unknown cell")
        _mm_need(!haskey(bytarget,target) && Tuple(row["vertices"])==native.triangles[(target[1],target[2]+1)],"selective affine row differs from actual target triangle order")
        w=_mr_weights(row["weights"]);record=(;parent,vertices=String.(row["vertices"]),weights=w)
        bytarget[target]=record;push!(get!(byparent,parent,Any[]),record)
    end
    _mm_need(Set(keys(bytarget))==Set((f,i-1) for (f,i) in keys(native.triangles)) && Set(keys(byparent))==Set(keys(names)),"selective affine ancestry omits an actual or raw cell")
    units=Dict();recovered=Dict()
    for (parent,rows) in byparent
        counts=Dict()
        _mm_need(sum(_mr_det(row.weights) for row in rows)==1,"selective child areas do not cover one raw parent")
        for row in rows
            for j in 1:3
                a,b=Tuple(row.weights[j]),Tuple(row.weights[mod1(j+1,3)])
                counts[(a,b)]=get(counts,(a,b),0)+1
            end
        end
        for (a,b) in collect(keys(counts))
            n=min(get(counts,(a,b),0),get(counts,(b,a),0));counts[(a,b)]-=n;counts[(b,a)]=get(counts,(b,a),0)-n
        end
        leftover=[(a,b) for ((a,b),n) in counts if n!=0]
        _mm_need(all(counts[(a,b)]==1 && count(j->a[j]==b[j]==0,1:3)==1 for (a,b) in leftover),"selective children do not cancel to the raw triangle boundary")
        for j in 1:3
            k=mod1(j+2,3);intervals=sort([(a[k],b[k]) for (a,b) in leftover if a[j]==b[j]==0])
            _mm_need(!isempty(intervals) && first(intervals)[1]==0 && last(intervals)[2]==1 && all(a<b for (a,b) in intervals) &&
                all(intervals[i][2]==intervals[i+1][1] for i in 1:length(intervals)-1),"selective parent edge interval has a gap, reversal or overlap")
        end
    end
    # Unit corners must be actual current vertices, consistently identified
    # across every old local occurrence, rather than restored from a cache.
    for (target,row) in bytarget,j in 1:3
        corner=findfirst(==(1),row.weights[j]);corner===nothing && continue
        key=(row.parent[1],names[row.parent][corner]);point=M.root!(native.vertex_dsu,(target[1],row.vertices[j]))
        _mm_need(!haskey(units,key) || units[key]==point,"selective inverse splits one old vertex");units[key]=point
    end
    for (parent,vs) in names
        _mm_need(all(haskey(units,(parent[1],v)) for v in vs),"selective inverse lacks an actual old unit corner")
        recovered[parent]=copy(vs)
    end
    actual=M.fine_edge_occurrences(native);derived=Dict()
    for (o,p) in actual
        isless(o,p) || continue
        a=bytarget[(o[1],o[2]-1)];b=bytarget[(p[1],p[2]-1)]
        wa=[a.weights[only(findall(==(v),a.vertices))] for v in o[3:4]]
        wb=[b.weights[only(findall(==(v),b.vertices))] for v in p[3:4]]
        a.parent==b.parent && wa==reverse(wb) && continue
        za=findall(j->wa[1][j]==wa[2][j]==0,1:3);zb=findall(j->wb[1][j]==wb[2][j]==0,1:3)
        _mm_need(length(za)==length(zb)==1,"selective current edge crosses unrelated raw triangle interiors")
        j,k=only(za),only(zb);ja,jb=mod1(j+2,3),mod1(k+2,3)
        _mm_need(wa[1][ja]+wb[2][jb]==1 && wa[2][ja]+wb[1][jb]==1,"selective raw edge parameter is not oppositely identified")
        va,vb=names[a.parent],names[b.parent]
        ra=(a.parent[1],a.parent[2]+1,va[mod1(j+1,3)],va[mod1(j+2,3)])
        rb=(b.parent[1],b.parent[2]+1,vb[mod1(k+1,3)],vb[mod1(k+2,3)])
        _mm_need(get(derived,ra,rb)==rb && get(derived,rb,ra)==ra,"selective raw edge interval has inconsistent opposite parents")
        derived[ra]=rb;derived[rb]=ra
    end
    _mm_need(derived==expected_pairs,"actual selective current adjacency does not reconstruct the complete raw opposite involution")
    (;recovered,opposite_edges=derived,unit_vertices=units)
end

function _mr_proper_arc_geometry(raw_geometry,opposite_edges)
    proof=_mr_proof(raw_geometry);supplied=Dict(_mr_occ(o)=>_mr_occ(p) for (o,p) in opposite_edges)
    _mm_need(supplied==proof.opposite_edges,"supplied opposite edges disagree with actual selective raw geometry")
    partition=_mrp_partition(proof);assembled=_mrp_assemble(partition,proof.geometry)
    ancestry=[Dict("source"=>[f,i-1],"source_vertices"=>collect(t)) for ((f,i),t) in sort(collect(proof.triangles);by=first)]
    _mrp_check_cells(assembled.native,assembled.affine,ancestry,supplied)
    corr=_mrp_correspondence(proof,partition,assembled)
    inverse=Dict{String,Any}("schema"=>"internal-triangle-partition/1","strategy"=>"proper_arc",
        "raw_geometry_hash"=>_mm_hash(proof.geometry),"public_geometry_hash"=>_mm_hash(_st_empty(assembled.document)),
        "raw_surface"=>deepcopy(proof.geometry["surface"]),"raw_charts"=>_mr_chart_provenance(proof.geometry),
        "triangle_ancestry"=>ancestry,"opposite_edges"=>_mr_pairs_json(supplied),
        "affine_children"=>assembled.affine,"proper_cuts"=>partition.trace,"correspondence_hash"=>_mm_hash(corr))
    (;document=assembled.document,correspondence=corr,inverse_specification=inverse)
end
function _mr_unpartition_geometry(current_native,spec)
    _mrp_inverse_schema(spec)
    M=MarkedSurfaceDiagrams;current=_marked_document(current_native);bare=_st_empty(current);native=M.validate(bare)
    _mm_need(_mm_hash(bare)==spec["public_geometry_hash"],"public selective geometry has changed")
    actual=_mrp_check_cells(native,spec["affine_children"],spec["triangle_ancestry"],_mr_pairs_live(spec["opposite_edges"]))
    charts=deepcopy(spec["raw_charts"])
    for p in charts
        n=pop!(p,"triangle_count");p["triangles"]=[actual.recovered[(p["face"],i-1)] for i in 1:n]
    end
    geometry=Dict{String,Any}("surface"=>deepcopy(spec["raw_surface"]),"patches"=>charts);proof=_mr_proof(geometry)
    _mm_need(proof.opposite_edges==actual.opposite_edges && _mm_hash(geometry)==spec["raw_geometry_hash"],"selective actual reassembly disagrees with old macroseams or raw geometry")
    regenerated=_mr_proper_arc_geometry(geometry,actual.opposite_edges)
    _mm_need(regenerated.document==bare && _canonical_json(regenerated.inverse_specification)==_canonical_json(spec) &&
        _mm_hash(regenerated.correspondence)==spec["correspondence_hash"],"fresh selective partition does not regenerate actual current geometry and full inverse data")
    (;geometry,opposite_edges=actual.opposite_edges,correspondence=regenerated.correspondence,construction_only=true)
end
# Source: the selective inverse consumes untrusted JSON before any allocation
# or indexing. Keep this separate from the unchanged uniform-star schema.
function _mrp_inverse_schema(spec)
    required=Set(["schema","strategy","raw_geometry_hash","public_geometry_hash","raw_surface","raw_charts",
        "triangle_ancestry","opposite_edges","affine_children","proper_cuts","correspondence_hash"])
    need(ok,field)=_mm_need(ok,"selective inverse has malformed "*field)
    text(x)=x isa AbstractString && !isempty(x)
    index(x)=x isa Integer && !(x isa Bool) && 0<=x<typemax(Int)
    strings(x,n)=x isa AbstractVector && length(x)==n && all(text,x) && length(Set(x))==n
    cell(x)=x isa AbstractVector && length(x)==2 && text(x[1]) && index(x[2])
    object(x,fields)=x isa AbstractDict && Set(keys(x))==Set(fields)
    rational(x,field)=begin
        need((x isa Rational) || (x isa Integer && !(x isa Bool)) ||
            (x isa AbstractString && occursin(r"^-?[0-9]+(?:/[1-9][0-9]*)?$",x)),field)
        try
            x isa Rational ? MarkedSurfaceDiagrams.Rat(x) : MarkedSurfaceDiagrams.qnum(x)
        catch err
            throw(ArgumentError("selective inverse has malformed "*field*": "*sprint(showerror,err)))
        end
    end
    need(Set(keys(spec))==required,"fields")
    need(spec["schema"]=="internal-triangle-partition/1" && spec["strategy"]=="proper_arc","schema/strategy")
    for field in ("raw_geometry_hash","public_geometry_hash","correspondence_hash")
        need(spec[field] isa AbstractString && occursin(r"^[0-9a-f]{64}$",spec[field]),field)
    end
    surface=spec["raw_surface"]
    need(object(surface,("marking_id","genus","gluings","boundary")),"raw_surface")
    need(text(surface["marking_id"]) && index(surface["genus"]),"raw_surface identity/genus")
    need(surface["boundary"] isa AbstractVector && isempty(surface["boundary"]),"raw_surface boundary")
    need(surface["gluings"] isa AbstractVector,"raw_surface gluings")
    for pair in surface["gluings"]
        need(pair isa AbstractVector && length(pair)==2 && all(x->x isa AbstractVector && length(x)==2 && all(text,x),pair),"raw_surface gluing pair")
    end
    for field in ("raw_charts","triangle_ancestry","opposite_edges","affine_children","proper_cuts")
        need(spec[field] isa AbstractVector && all(x->x isa AbstractDict,spec[field]),field)
    end
    charts=Dict();parents=Dict()
    for p in spec["raw_charts"]
        needed=Set(["face","side_order","sides","triangle_count"])
        need(issubset(needed,Set(keys(p))) && issubset(Set(keys(p)),union(needed,Set(["side_parameters"]))),"raw_charts fields")
        need(text(p["face"]) && !haskey(charts,p["face"]),"raw_charts face")
        need(index(p["triangle_count"]) && 0<p["triangle_count"]<=length(spec["triangle_ancestry"]),"raw_charts triangle_count")
        order=p["side_order"];sides=p["sides"]
        need(order isa AbstractVector && length(order)>=3 && all(text,order) && length(Set(order))==length(order),"raw_charts side_order")
        need(sides isa AbstractDict && Set(keys(sides))==Set(order),"raw_charts sides")
        for path in values(sides);need(path isa AbstractVector && length(path)>=2 && all(text,path) && length(Set(path))==length(path),"raw_charts side path");end
        parameters=get(p,"side_parameters",Dict())
        need(parameters isa AbstractDict && issubset(Set(keys(parameters)),Set(order)),"raw_charts side_parameters")
        for (side,qs) in parameters
            need(qs isa AbstractVector && length(qs)==length(sides[side]),"raw_charts side parameter count")
            ts=[rational(q,"raw_charts side parameter") for q in qs]
            need(first(ts)==0 && last(ts)==1 && all(ts[i]<ts[i+1] for i in 1:length(ts)-1),"raw_charts side parameter order")
        end
        charts[p["face"]]=p
    end
    need(!isempty(charts),"empty raw_charts")
    for row in spec["triangle_ancestry"]
        need(object(row,("source","source_vertices")) && cell(row["source"]) && strings(row["source_vertices"],3),"triangle_ancestry row")
        f,i=row["source"];need(haskey(charts,f) && i<charts[f]["triangle_count"] && !haskey(parents,(f,i)),"triangle_ancestry parent")
        parents[(f,i)]=row["source_vertices"]
    end
    need(all(haskey(parents,(f,i-1)) for (f,p) in charts for i in 1:p["triangle_count"]),"triangle_ancestry coverage")
    for row in spec["affine_children"]
        need(object(row,("source","target","vertices","weights")) && cell(row["source"]) && cell(row["target"]) && strings(row["vertices"],3),"affine_children row")
        need(haskey(parents,Tuple(row["source"])),"affine_children source")
        ws=row["weights"];need(ws isa AbstractVector && length(ws)==3 && all(w->w isa AbstractVector && length(w)==3,ws),"affine_children weights")
        _mr_weights([[rational(q,"affine_children rational") for q in w] for w in ws])
    end
    for row in spec["opposite_edges"]
        need(object(row,("source","target")),"opposite_edges row")
        for field in ("source","target")
            o=row[field]
            need(o isa AbstractVector && length(o)==4 && text(o[1]) && index(o[2]) && text(o[3]) && text(o[4]) && o[3]!=o[4],"opposite_edges occurrence")
            need(haskey(parents,(o[1],o[2])) && all(v->v in parents[(o[1],o[2])],o[3:4]),"opposite_edges parent vertices")
        end
    end
    for row in spec["proper_cuts"]
        need(object(row,("source","path_triangles","entry_edge","exit_edge","children","mirrored_triangle_splits")),"proper_cuts fields")
        need(text(row["source"]) && row["path_triangles"] isa AbstractVector && !isempty(row["path_triangles"]) &&
            all(index,row["path_triangles"]) && strings(row["entry_edge"],2) && strings(row["exit_edge"],2) &&
            strings(row["children"],2) && index(row["mirrored_triangle_splits"]),"proper_cuts row")
    end
    nothing
end
