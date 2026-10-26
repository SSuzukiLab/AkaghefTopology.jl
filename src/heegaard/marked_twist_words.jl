"""Private exact RIGHT composition of maps on one fixed prepared mesh.

Sources: marked_twists.jl::_mt_bary/_mt_weight_product and the Polyhedra-backed
_mt_convex_boundary_triangles; marked_polygon_overlay.jl exact intersections;
vendor fine_edge_occurrences and vertex quotient. These are INTERNAL affine
atoms, not native diagrams, public word constructors, or homology substitutes.
The caller supplies the SAME already validated closed source/target mesh.
"""

_mtw_q(x) = x isa Rational ? MarkedSurfaceDiagrams.Rat(x) : _mt_q(x)
_mtw_xy(w) = (w[2],w[3])

function _mtw_parent(x,parents)
    _mm_need((x isa AbstractVector || x isa Tuple) && length(x)==2 &&
        x[1] isa AbstractString && x[2] isa Integer && !(x[2] isa Bool),
        "word atom parent must be a fixed face/zero-based triangle")
    p=(String(x[1]),Int(x[2]))
    _mm_need(p in parents,"word atom refers to a triangle outside the fixed mesh")
    p
end

function _mtw_rows(rows)
    _mm_need((rows isa AbstractVector || rows isa Tuple) && length(rows)==3 &&
        all(r->(r isa AbstractVector || r isa Tuple) && length(r)==3,rows),
        "word atom needs three paired barycentric rows")
    w=[_mtw_q.(r) for r in rows]
    _mm_need(all(r->sum(r)==1 && all(>=(0),r),w) && _mt_ratio(w)>0,
        "word atom barycentric triangle is outside its parent or folded")
    w
end

function _mtw_atoms(input,mesh)
    _mm_need(input isa AbstractVector && !isempty(input),"word map needs a nonempty full partition")
    parents=Set((f,i-1) for (f,i) in keys(mesh.triangles))
    [begin
        _mm_need(row isa AbstractDict && all(haskey(row,k) for k in
            ("source","target","source_weights","target_weights")),"word map atom is missing paired data")
        (;source=_mtw_parent(row["source"],parents),target=_mtw_parent(row["target"],parents),
            source_weights=_mtw_rows(row["source_weights"]),target_weights=_mtw_rows(row["target_weights"]))
    end for row in input]
end

# Source: vendor oriented edge involution. The carrier includes the actual
# fine-edge occurrence, so distinct loop edges with equal quotient endpoints
# cannot be confused. Coordinates on opposite occurrences run oppositely.
function _mtw_carrier(parent,a,b,mesh,pairs)
    zeroes=[j for j in 1:3 if a[j]==0 && b[j]==0]
    if length(zeroes)==1
        t=mesh.triangles[(parent[1],parent[2]+1)];j=only(zeroes)
        u,v=only((u,v) for (u,v) in MarkedSurfaceDiagrams.tedges(t) if u!=t[j] && v!=t[j])
        occurrence=(parent[1],parent[2]+1,u,v);canonical=min(occurrence,pairs[occurrence])
        vi=only(findall(==(v),t));pa,pb=a[vi],b[vi]
        occurrence!=canonical && ((pa,pb)=(1-pa,1-pb))
        return (("edge",canonical,min(pa,pb),max(pa,pb)),pa<pb)
    end
    x,y=Tuple(a),Tuple(b)
    (("cell",parent,min(x,y),max(x,y)),isless(x,y))
end

function _mtw_point(parent,w,mesh,pairs)
    M=MarkedSurfaceDiagrams;t=mesh.triangles[(parent[1],parent[2]+1)]
    corner=findfirst(==(1),w)
    corner!==nothing && return ("vertex",M.root!(mesh.vertex_dsu,(parent[1],t[corner])))
    j=findfirst(==(0),w)
    if j!==nothing
        u,v=only((u,v) for (u,v) in M.tedges(t) if u!=t[j] && v!=t[j])
        o=(parent[1],parent[2]+1,u,v);canonical=min(o,pairs[o])
        p=w[only(findall(==(v),t))];o!=canonical && (p=1-p)
        return ("edge",canonical,p)
    end
    ("cell",parent,Tuple(w))
end

function _mtw_partition(atoms,side,mesh)
    Q=MarkedSurfaceDiagrams.Rat;parents=Set((f,i-1) for (f,i) in keys(mesh.triangles))
    groups=Dict(p=>Int[] for p in parents)
    for (i,a) in enumerate(atoms);push!(groups[getproperty(a,side)],i);end
    areas=Any[]
    for p in sort(collect(parents))
        ids=groups[p];_mm_need(!isempty(ids),"word map omits a fixed parent triangle")
        rows=[getproperty(atoms[i],Symbol(string(side),"_weights")) for i in ids]
        _mm_need(sum(_mt_ratio(w) for w in rows)==1,"word map does not cover a fixed parent in area")
        triangles=[(;id=i,points=_mtw_xy.(w)) for (i,w) in zip(ids,rows)]
        # This self-overlay is a PAIRWISE nonoverlap check, not merely a total
        # area comparison: every target already receives its own entire area.
        # Any positive intersection with a distinct atom adds positive area.
        overlay=_marked_triangle_overlay(triangles,triangles)
        _mm_need(all(c->c.source_id==c.target_id,overlay.cells) &&
            all(c->c.area_matches,overlay.coverage),"word map parent partition has overlapping interiors")
        push!(areas,Dict("parent"=>collect(p),"atoms"=>copy(ids),"area_ratio"=>MarkedSurfaceDiagrams.qstr(Q(1))))
    end
    (;groups,areas)
end

function _mtw_continuity(atoms,side,mesh,pairs)
    other=side==:source ? :target : :source
    edges=Dict{Any,Vector{Any}}();points=Dict();breaks=Dict{Any,Set{MarkedSurfaceDiagrams.Rat}}()
    for (i,a) in enumerate(atoms)
        p=getproperty(a,side);q=getproperty(a,other)
        w=getproperty(a,Symbol(string(side),"_weights"));z=getproperty(a,Symbol(string(other),"_weights"))
        for j in 1:3
            key=_mtw_point(p,w[j],mesh,pairs);value=_mtw_point(q,z[j],mesh,pairs)
            _mm_need(!haskey(points,key) || points[key]==value,"word map has inconsistent vertex images")
            points[key]=value
        end
        for (j,k) in ((1,2),(2,3),(3,1))
            key,direction=_mtw_carrier(p,w[j],w[k],mesh,pairs)
            image,image_direction=_mtw_carrier(q,z[j],z[k],mesh,pairs)
            push!(get!(edges,key,Any[]),(;atom=i,direction,image,image_direction))
            if first(key)=="edge"
                ps=get!(breaks,key[2],Set{MarkedSurfaceDiagrams.Rat}())
                push!(ps,key[3],key[4])
            end
        end
    end
    for occurrences in values(edges)
        _mm_need(length(occurrences)==2,"word partition has an unmatched edge, seam breakpoint or T junction")
        a,b=occurrences
        _mm_need(a.direction!=b.direction,"word partition edge orientations do not cancel")
        _mm_need(a.image==b.image && a.image_direction!=b.image_direction,
            "word affine maps disagree across an actual partition edge or fixed seam")
    end
    Dict("vertices"=>length(points),"paired_edges"=>length(edges),
        "edge_index_convention"=>"fixed triangle occurrences are zero-based; includes internal edges and marking seams",
        "fixed_edge_breakpoints"=>[Dict("edge"=>[e[1],e[2]-1,e[3],e[4]],"parameters"=>MarkedSurfaceDiagrams.qstr.(sort(collect(ps))))
            for (e,ps) in sort(collect(breaks);by=first)])
end

"""Certify full positive exact partitions and both directions of continuity.

Source: exact convex intersection supplies pairwise nonoverlap; every parent
has area one, and every actual partition edge is paired, including old seams.
Both endpoint maps and affine edge carriers agree on each pairing. Together
these establish the piecewise-affine homeomorphism, not area sums alone.
"""
function _mtw_certify(input,mesh)
    _mm_need(mesh.surface_report["boundary_components"]==0,"private RIGHT maps currently require a closed fixed mesh")
    atoms=_mtw_atoms(input,mesh);pairs=MarkedSurfaceDiagrams.fine_edge_occurrences(mesh)
    source=_mtw_partition(atoms,:source,mesh);target=_mtw_partition(atoms,:target,mesh)
    forward=_mtw_continuity(atoms,:source,mesh,pairs);backward=_mtw_continuity(atoms,:target,mesh,pairs)
    (;atoms,source,target,correspondence=Dict("scope"=>"private exact atoms on one fixed prepared mesh; no native word result",
        "source_partition"=>source.areas,"target_partition"=>target.areas,
        "forward_continuity"=>forward,"inverse_continuity"=>backward))
end

# Source: marked_twists.jl's existing Polyhedra-backed polygon triangulation.
# Select BOTH row lists with the SAME returned vertex indices; independently
# triangulating/reordering source and target would describe a different map.
function _mtw_paired_triangles(points,source_rows,target_rows)
    _mm_need(length(points)==length(source_rows)==length(target_rows),"paired polygon row lengths differ")
    [begin
        s=source_rows[indices];t=target_rows[indices]
        _mm_need(_mt_ratio(s)>0 && _mt_ratio(t)>0,"paired composition triangle folded")
        (;source_weights=s,target_weights=t,indices)
    end for indices in _mt_convex_boundary_triangles(points)]
end

"""Compose LEFT F after RIGHT G by G.target/F.source intersections.

All parents refer to the SAME fixed mesh. No angular shifts from different
annuli are added. The resulting internal atoms have explicit factor ancestry
and are re-certified before return; no native realization is attempted.
"""
function _mtw_compose_maps(left,right,mesh)
    F=_mtw_certify(left,mesh);G=_mtw_certify(right,mesh)
    output=Any[];pieces=Any[]
    for middle in sort(collect(keys(F.source.groups)))
        fs=F.source.groups[middle];gs=G.target.groups[middle]
        image=[(;id=i,points=_mtw_xy.(G.atoms[i].target_weights)) for i in gs]
        domain=[(;id=i,points=_mtw_xy.(F.atoms[i].source_weights)) for i in fs]
        overlay=_marked_triangle_overlay(image,domain)
        _mm_need(all(r->r.area_matches,overlay.coverage),"RIGHT overlay lost an intermediate target cell")
        for cell in overlay.cells
            g=G.atoms[cell.source_id];f=F.atoms[cell.target_id]
            sw=[begin a=_mt_bary(p,_mtw_xy.(g.target_weights));
                [sum(a[j]*g.source_weights[j][k] for j in 1:3) for k in 1:3] end for p in cell.vertices]
            tw=[begin b=_mt_bary(p,_mtw_xy.(f.source_weights));
                [sum(b[j]*f.target_weights[j][k] for j in 1:3) for k in 1:3] end for p in cell.vertices]
            for triangle in _mtw_paired_triangles(cell.vertices,sw,tw)
                push!(output,Dict("source"=>collect(g.source),"target"=>collect(f.target),
                    "source_weights"=>_mt_json_bary.(triangle.source_weights),"target_weights"=>_mt_json_bary.(triangle.target_weights)))
                push!(pieces,Dict("atom"=>length(output),"middle_parent"=>collect(middle),
                    "left_atom"=>cell.target_id,"right_atom"=>cell.source_id,
                    "intersection_vertices"=>_mt_json_point.(cell.vertices),"selected_vertices"=>triangle.indices))
            end
        end
    end
    checked=_mtw_certify(output,mesh)
    (;atoms=output,correspondence=Dict("composition"=>"F composed with G; RIGHT multiplication by G",
        "left"=>F.correspondence,"right"=>G.correspondence,"pieces"=>pieces,"result"=>checked.correspondence))
end

"""Simultaneous RIGHT transvection in the fixed PRE-WORD curve basis.

Source: T_(F(c)) F = F T_c and the oriented signed twist formula. Original
intersection counts are supplied separately from the actual output map.
This audits saved closed-curve vectors, not an unproved full H1 basis.
"""
function _mtw_right_update(vectors,core,sign,original_intersections)
    _mm_need(sign isa Integer && !(sign isa Bool) && sign in (-1,1),"RIGHT factor sign must be +1 or -1")
    _mm_need(haskey(vectors,core),"RIGHT core is missing from the original curve vectors")
    n=length(vectors[core])
    _mm_need(all(v->v isa AbstractVector && length(v)==n && all(x->x isa Integer && !(x isa Bool),v),values(vectors)),
        "RIGHT vectors must share one integer source basis")
    # Source: the exact integer contract above. Machine-Int inputs must not
    # silently wrap when a word multiplies their transvection coefficients.
    previous=Dict(id=>BigInt.(v) for (id,v) in vectors)
    _mm_need(all(haskey(original_intersections,(core,id)) for id in keys(previous)),"original oriented intersection count is missing")
    _mm_need(original_intersections[(core,core)]==0,"a core has nonzero self-intersection")
    for id in keys(previous)
        x=original_intersections[(core,id)]
        _mm_need(x isa Integer && !(x isa Bool),"original intersection counts must be integers")
    end
    Dict(id=>previous[id]+BigInt(sign)*BigInt(original_intersections[(core,id)])*previous[core] for id in keys(previous))
end

"""Prepare once and certify original-core products on that SAME mesh pair.

Sources: marked_twists.jl::_mt_prepare supplies the existing accepted
dual-star/one-subdivision history; marked_annulus_product.jl independently
checks exact flags, full core, minimum length, annulus, orientation and lifts.
Only candidate support membership is selected here. No parallel frontier,
neighborhood proof, extra refinement or substitute core is constructed.
"""
function _mtw_prepare_original_products(input,requested_cores)
    _mm_need(requested_cores isa AbstractVector && !isempty(requested_cores) &&
        all(c->c isa AbstractString && !isempty(c),requested_cores),
        "private word preparation needs nonempty original core IDs")
    original=_marked_document(input);cores=unique(String.(requested_cores))
    for field in ("curves","crossings","region_marks"),row in get(original,field,[])
        _mm_need(row["id"] isa AbstractString && !isempty(row["id"]),
            "private word preparation requires nonempty string source IDs")
    end
    curves=Dict(c["id"]=>c for c in original["curves"])
    for cid in cores
        _mm_need(haskey(curves,cid) && curves[cid]["closed"],
            "private word preparation core must be an existing original closed curve")
    end
    # Exactly one call constructs the shared mesh pair. Do not recover from a
    # later product rejection by silently preparing/refining that core anew.
    preparation=_mt_prepare(original,first(cores))
    coarse=_marked_states(preparation.coarse).geometry
    fine=_marked_states(preparation.fine).geometry
    M=MarkedSurfaceDiagrams;products=Dict{String,Any}();supports=Dict{String,Any}()
    for cid in cores
        ds=coarse.curves[cid][2]
        roots=Set(M.root!(coarse.vertex_dsu,(d.face,d.u)) for d in ds)
        # Defining candidate derived-star membership only. The existing product
        # checker RECONSTRUCTS all six actual flags and verifies exact equality
        # with this candidate, rather than trusting this small selector.
        selected=Set((f,i-1) for ((f,i),t) in fine.triangles if
            any(haskey(coarse.vertex_dsu.p,(f,v)) &&
                M.root!(coarse.vertex_dsu,(f,v)) in roots for v in t))
        if cid==first(cores)
            _mm_need(selected==preparation.support && preparation.product.support_cells==selected,
                "first core product does not match the single preparation support")
            products[cid]=preparation.product
        else
            products[cid]=_marked_annulus_product(coarse,fine,cid,selected)
        end
        supports[cid]=selected
    end
    basis=marked_curve_words(preparation.fine)
    correspondence=Dict("scope"=>"private original-core products on one shared prepared mesh; no native word result",
        "requested_cores"=>String.(requested_cores),"product_cores"=>cores,
        "source_hash"=>_mm_hash(original),"coarse_source_hash"=>_mm_hash(preparation.coarse),
        "fine_source_hash"=>_mm_hash(preparation.fine),
        "fine_marking_hash"=>M.marking_digest(preparation.fine),
        "basis_scope"=>"fixed PRE-WORD prepared presentation; input marking identified by preparation history",
        "construction"=>preparation.history,
        "products"=>[Dict("core"=>cid,"support_cells"=>[collect(c) for c in sort(collect(supports[cid]))],
            "certificate"=>products[cid].certificate) for cid in cores])
    (;preparation,coarse,fine,products,supports,basis,correspondence)
end


"""Private realization of a certified total map, not an ordered-word API.

Source: the shared fixed-vertex assembler/native acceptance extraction. Only
paired triangle rows are passed; no h-shifts from different factors are added.
Public word semantics, signed factors and full independent replay are pending.
"""
function _mtw_realize_total_map(prepared,atoms)
    checked=_mtw_certify(atoms,prepared.fine)
    paired=[(;source=a.source,target=a.target,shift=0,
        source_weights=a.source_weights,target_weights=a.target_weights,
        triangles=[[1,2,3]]) for a in checked.atoms]
    raw=_mt_realize_map_atoms(prepared.preparation.fine,paired)
    result=_mt_accept_map_complex(prepared.preparation.fine,raw)
    (;result,map_certificate=checked.correspondence,
        scope="private total-map realization; no public ordered-word witness")
end

"""Pull the ACTUAL target back only to its prepared source candidate.

Source: existing actual feature unchart and shared source-edge pullback.
Full current-region cell sets use THIS call's regenerated rechart maps. The
result cannot be returned as a public word inverse: the preparation must still
be genuinely collapsed and the entire signed word independently regenerated.
"""
function _mtw_restore_total_source_candidate(input,specification)
    M=MarkedSurfaceDiagrams;current=_marked_document(input);spec=_plain_value(specification)
    _mm_need(_mm_hash(current)==spec["target_hash"],"total-map target has changed")
    recovered=_mt_unrechart_features_data(current,spec);proof=recovered.proof
    geometry=_mt_map_inverse_geometry(proof,spec)
    certified=_mtw_certify(spec["child_cells"],geometry.old)
    pulled=_mt_pullback_map_source(proof,spec,geometry)
    before=_marked_states(pulled.candidate);after=_marked_states(current)
    native=recovered.rechart===nothing ? nothing : Dict(_mm_cell(r["source"])=>
        [_mm_cell(c) for c in r["target"]] for r in recovered.rechart.correspondence["cell_children"])
    children=Dict((f,i-1)=>Tuple{String,Int}[] for (f,i) in keys(geometry.old.triangles))
    for r in spec["child_cells"]
        child=_mm_cell(r["child"])
        append!(children[_mm_cell(r["source"])],native===nothing ? [child] : native[child])
    end
    source_math=M.regions(before.mathematics);target_math=M.regions(after.mathematics)
    source_geo=M.regions(before.geometry);target_geo=M.regions(after.geometry)
    regions=_ms_regions(source_math,target_math,children)
    geometric=_ms_regions(source_geo,target_geo,children)
    function same_region(source_seed,target_seed,oldregions,newregions,maps)
        source=findall(r->_mm_cell(source_seed) in _mm_cells(r),oldregions)
        target=findall(r->_mm_cell(target_seed) in _mm_cells(r),newregions)
        _mm_need(length(source)==length(target)==1,"total-map protected seed is not in one region")
        actual=only(r["target_region"] for r in maps if r["source_region"]==only(source))
        _mm_need(actual==only(target),"total-map inverse changes the selected full region")
    end
    marks=Dict(m["id"]=>m for m in get(current,"region_marks",[]))
    for m in get(pulled.candidate,"region_marks",[])
        same_region(m["seed"],marks[m["id"]]["seed"],source_math,target_math,regions)
    end
    infinity=_mm_infinity(pulled.candidate)
    infinity!==nothing && same_region(infinity,_mm_infinity(current),source_geo,target_geo,geometric)
    candidate=_mt_finish_map_source(pulled.candidate,pulled.word_formats,spec)
    (;candidate,map_certificate=certified.correspondence,regions,geometric_regions=geometric,
        forward_word_regeneration_required=true)
end
