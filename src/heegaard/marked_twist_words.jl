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
    actual=_mtw_compose_fresh(F,G,mesh)
    (;atoms=actual.atoms,correspondence=actual.correspondence)
end

# Source: _mtw_certify above. ONLY this constructor call's fresh certificates
# enter here; this is not authentication of caller-supplied NamedTuples. The
# ordinary atom-input entry point still runs both full input certifications.
# This helper never mutates its fresh input atoms/maps. Newly composed output
# still runs the complete positive partition/continuity certificate below.
function _mtw_compose_fresh(F,G,mesh)
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
    (;atoms=output,checked,correspondence=Dict("composition"=>"F composed with G; RIGHT multiplication by G",
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
    _mtw_realize_fresh(prepared,_mtw_certify(atoms,prepared.fine))
end

# Source: the SAME-call fresh factor or composed-output certificate. No map,
# feature or basis result is reused; the actual total native mesh, all roles,
# D/infinity and old-basis action are still constructed and checked normally.
# The public constructor never accepts an external certificate/skip option.
function _mtw_realize_fresh(prepared,checked)
    paired=[(;source=a.source,target=a.target,shift=0,
        source_weights=a.source_weights,target_weights=a.target_weights,
        triangles=[[1,2,3]]) for a in checked.atoms]
    raw=_mt_realize_map_atoms(prepared.preparation.fine,paired)
    result=_mtw_accept_map_complex(prepared.preparation.fine,raw)
    (;result,map_certificate=checked.correspondence,
        scope="private total-map realization; no public ordered-word witness")
end

# Source: the exact private chart-budget guard in marked_rechart.jl. The
# SAME freshly built raw map is retried only for that one failure. All other
# geometry/schema/feature errors propagate, and old successful outputs remain
# unchanged. Single-twist and other legacy callers do not use this fallback.
function _mtw_accept_map_complex(fine_doc,raw)
    try
        _mt_accept_map_complex(fine_doc,raw)
    catch err
        err isa ArgumentError && err.msg=="selective disk partition exceeds native chart count" || rethrow()
        _mt_accept_map_complex(fine_doc,raw;partition_policy="chart_demand/1")
    end
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


"""Normalize a nonempty signed word; IDs name ORIGINAL saved cores.

Sequence [a,b] means T_a*T_b with the rightmost map acting first. No family
selector, reversed-negative inverse shortcut or algebraic relation solver is
part of this whole-diagram action.
"""
function _mtw_word_specification(specification)
    spec=_plain_value(specification)
    _mm_need(spec isa AbstractDict && Set(keys(spec))==Set(["word"]),"word action accepts only word")
    word=spec["word"]
    _mm_need(word isa AbstractVector && !isempty(word),"word action needs nonempty signed factors")
    steps=Any[]
    for step in word
        _mm_need(step isa AbstractDict && issubset(Set(keys(step)),Set(["core","sign"])) &&
            haskey(step,"core") && step["core"] isa AbstractString && !isempty(step["core"]),
            "word factor needs an original nonempty core ID")
        sign=get(step,"sign",1)
        _mm_need(sign isa Integer && !(sign isa Bool) && sign in (-1,1),"word factor sign must be +1 or -1")
        push!(steps,Dict("core"=>String(step["core"]),"sign"=>Int(sign)))
    end
    Dict("word"=>steps)
end

function _mtw_source_action(prepared,steps)
    data=prepared.basis;dimension=length(data.presentation.generators)
    vectors=Dict{String,Vector{BigInt}}()
    for (cid,word) in data.curve_words
        v=zeros(BigInt,dimension)
        for (i,s) in word;v[i]+=s;end
        vectors[cid]=v
    end
    cores=unique([s["core"] for s in steps])
    intersection=Dict((c,id)=>BigInt(0) for c in cores for id in keys(vectors))
    for crossing in prepared.preparation.fine["crossings"]
        a,b=crossing["curves"];s=BigInt(crossing["sign"])
        a in cores && (intersection[(a,b)]+=s)
        b in cores && (intersection[(b,a)]-=s)
    end
    (;vectors,intersection)
end

function _mtw_encoded_atoms(atoms)
    [Dict("source"=>collect(a.source),"target"=>collect(a.target),
        "source_weights"=>_mt_json_bary.(a.source_weights),
        "target_weights"=>_mt_json_bary.(a.target_weights)) for a in atoms]
end

function _mtw_word_inverse_data(original,prepared,result,spec,correspondence)
    fine=prepared.preparation.fine;raw=result.raw
    inverse=Dict{String,Any}("schema"=>"marked-dehn-word-inverse/1",
        "forward_specification"=>deepcopy(spec),"word"=>deepcopy(spec["word"]),
        "construction"=>prepared.preparation.history,
        "fine_source_mesh"=>Dict("surface"=>deepcopy(fine["surface"]),"patches"=>deepcopy(fine["patches"])),
        "fine_source_hash"=>_mm_hash(fine),"source_hash"=>_mm_hash(original),"target_hash"=>_mm_hash(result.document),
        "child_cells"=>raw.records,"new_vertices"=>raw.vertex_records,"curve_images"=>raw.image_records,
        "proof_complex_scope"=>"internal forward pre-rechart total-map proof; not current native incidence",
        "native_child_cells"=>result.records,"native_curve_images"=>result.image_records,
        "rechart_inverse_specification"=>result.rechart===nothing ? nothing : result.rechart.inverse_specification,
        "proof_crossing_locations"=>Dict(x["id"]=>deepcopy(x["at"]) for x in raw.complex["crossings"]),
        "proof_region_seeds"=>_mm_inverse_seeds(raw.complex),"region_seeds"=>_mm_inverse_seeds(fine),
        "forward_correspondence"=>deepcopy(correspondence))
    _mm_infinity(fine)!==nothing && (inverse["infinity_triangle"]=deepcopy(_mm_infinity(fine)))
    _mm_infinity(raw.complex)!==nothing && (inverse["proof_infinity_triangle"]=deepcopy(_mm_infinity(raw.complex)))
    inverse
end

"""Apply an actual whole-diagram RIGHT word on one original prepared mesh.

Source: accepted prepare-once/core composition, exact signed annular maps and
the shared total-map assembler. T_(F(c)) composed with F equals F composed
with T_c on transported support F(N_c). ALL roles, D and infinity move.
Only the FINAL total map is natively realized/recharted. Resource failure is
explicit; no homology result or retained source can replace that realization.
"""
function marked_dehn_twist_word(input::AbstractDict,specification::AbstractDict)
    original=_marked_document(input);spec=_mtw_word_specification(specification);steps=spec["word"]
    prepared=_mtw_prepare_original_products(original,[s["core"] for s in steps])
    action=_mtw_source_action(prepared,steps);expected=deepcopy(action.vectors)
    factor_table=Dict{Tuple{String,Int},Int}();factors=Any[];stages=Any[];total=nothing
    # Local proof values only; never serialized, returned, cached or accepted
    # from the API. Their atom arrays remain read-only for this whole call.
    factor_checks=Any[];total_checked=nothing
    for (i,step) in enumerate(steps)
        cid,sign=step["core"],step["sign"];key=(cid,sign)
        if !haskey(factor_table,key)
            product=prepared.products[cid]
            overlay=_mt_periodic_overlay(product,sign)
            raw=_mt_overlay_complex(prepared.preparation.fine,product,overlay,sign)
            certified=_mtw_certify(raw.records,prepared.fine)
            atoms=_mtw_encoded_atoms(certified.atoms)
            record=Dict("core"=>cid,"sign"=>sign,"atoms"=>atoms,
                "product_triangles"=>[Dict("cell"=>deepcopy(t.id),"points"=>_mt_json_point.(t.points)) for t in product.triangles],
                "product_certificate"=>product.certificate,"map_certificate"=>certified.correspondence)
            push!(factors,record);push!(factor_checks,certified);factor_table[key]=length(factors)
        end
        factor=factor_table[key];atoms=factors[factor]["atoms"]
        if total===nothing
            total=deepcopy(atoms);total_checked=factor_checks[factor]
            composition=Dict("kind"=>"first actual elementary map","factor"=>factor)
        else
            composed=_mtw_compose_fresh(total_checked,factor_checks[factor],prepared.fine)
            total=composed.atoms;total_checked=composed.checked;composition=composed.correspondence
        end
        expected=_mtw_right_update(expected,cid,sign,action.intersection)
        push!(stages,Dict("step"=>i,"core"=>cid,"sign"=>sign,"factor"=>factor,
            "composition"=>composition,"output_atoms"=>length(total),
            "expected_vectors"=>[Dict("curve"=>id,"exponents"=>expected[id]) for id in sort(collect(keys(expected)))]))
    end
    realized=_mtw_realize_fresh(prepared,total_checked);result=realized.result
    annotations=Dict(id=>Dict("word_semantics"=>"fixed PRE-WORD prepared surface presentation; source marking identified by construction; exponent equality checked",
        "original_factor_intersections"=>[action.intersection[(s["core"],id)] for s in steps]) for id in keys(expected))
    actual_action=_mt_map_old_basis_images(prepared.preparation.fine,result,prepared.basis,expected,annotations)
    M=MarkedSurfaceDiagrams
    correspondence=Dict("word"=>deepcopy(steps),
        "word_convention"=>"RIGHT multiplication: sequence [a,b] is T_a*T_b; rightmost map acts first; later cores/supports are transported",
        "action_scope"=>"whole alpha/beta/aux diagram and D/infinity; not beta-only",
        "construction"=>prepared.preparation.history,"elementary_maps"=>factors,"composition_stages"=>stages,
        "background_marking"=>"fixed original prepared target cells via total S/T and explicit identity rechart; public chart names may differ",
        "child_cells"=>result.records,"raw_proof_vertices"=>result.raw.vertex_records,
        "vertex_scope"=>"actual native children carry old fixed parents; raw_proof_vertices are pre-rechart only",
        "rechart"=>result.rechart===nothing ? nothing : result.rechart.correspondence,
        "curve_images"=>result.image_records,"regions"=>result.region_map,"geometric_regions"=>result.geometry_regions,
        "total_map_certificate"=>realized.map_certificate,"old_basis_action"=>actual_action,
        "basis_scope"=>"ONE fixed PRE-WORD prepared presentation; not output basis or numeric input generator indices",
        "basis_source_hash"=>M.representation_digest(prepared.preparation.fine),
        "basis_geometry_hash"=>M.representation_digest(_mm_bare(prepared.preparation.fine)),
        "basis_marking_hash"=>M.marking_digest(prepared.preparation.fine),
        "crossings"=>[Dict("source"=>x["id"],"target"=>x["id"],"at"=>x["at"],"sign"=>x["sign"]) for x in result.document["crossings"]])
    inverse=_mtw_word_inverse_data(original,prepared,result,spec,correspondence)
    witness=_mm_witness("dehn_twist_word",original,result.document,spec,inverse,correspondence)
    witness["inverse_kind"]="generated_twist_word_inverse"
    (;document=result.document,witness)
end

# Source: schema.jl native v1 contract; marked_subdivision.jl generated
# preparation inverses; marked_rechart.jl::_mrp_inverse_schema; the exact
# word-inverse/1 fields emitted above. This guard is NEW-WORD-ONLY. It never
# replaces geometric reconstruction or full factor/word regeneration.
function _mtwi_json_shape(input)
    stack=Any[(input,0,false)];active=Base.IdSet{Any}()
    while !isempty(stack)
        x,depth,leaving=pop!(stack)
        leaving && (delete!(active,x);continue)
        _mm_need(depth<=64,"word inverse JSON exceeds the fixed schema nesting depth")
        if x isa AbstractDict || x isa AbstractVector || x isa Tuple
            _mm_need(!(x in active),"word inverse JSON contains a container cycle")
            push!(active,x);push!(stack,(x,depth,true))
            if x isa AbstractDict
                _mm_need(all(k->k isa AbstractString || k isa Symbol,keys(x)),"word inverse JSON has a nontext key")
                for value in values(x);push!(stack,(value,depth+1,false));end
            else
                for value in x;push!(stack,(value,depth+1,false));end
            end
        else
            _mm_need(x===nothing || x isa AbstractString || x isa Integer,
                "word inverse JSON has a nonexact or unsupported scalar")
        end
    end
    nothing
end

function _mtwi_inverse_schema(input,current)
    _mtwi_json_shape(input);spec=_plain_value(input);M=MarkedSurfaceDiagrams
    fail(path,message)=throw(ArgumentError("word inverse "*path*": "*message))
    need(ok,path,message="malformed value")=ok ? nothing : fail(path,message)
    text(x,path)=need(x isa AbstractString && !isempty(x),path,"expected nonempty text")
    integer(x,path;lo=0,hi=typemax(Int)-1)=need(x isa Integer && !(x isa Bool) && lo<=x<=hi,path,"integer outside schema bounds")
    function object(x,fields,path;optional=())
        need(x isa AbstractDict,path,"expected object")
        required=Set(fields);actual=Set(keys(x))
        need(issubset(required,actual) && issubset(actual,union(required,Set(optional))),path,"unexpected or missing keys")
        x
    end
    function vector(x,path;lo=0,hi=typemax(Int)-1)
        need(x isa AbstractVector && lo<=length(x)<=hi,path,"array size/type outside schema bounds");x
    end
    function hash(x,path)
        need(x isa AbstractString && occursin(r"^[0-9a-f]{64}$",x),path,"expected SHA256 text")
    end
    function rational(x,path)
        need((x isa Integer && !(x isa Bool)) || x isa Rational ||
            (x isa AbstractString && occursin(r"^-?[0-9]+(?:/[1-9][0-9]*)?$",x)),path,"invalid rational literal")
        try _mtw_q(x) catch err;fail(path,sprint(showerror,err));end
    end
    function strings(x,n,path;distinct=false)
        vector(x,path;lo=n,hi=n);for v in x;text(v,path);end
        need(!distinct || length(Set(x))==n,path,"duplicate names");x
    end
    function cell(x,path,counts=nothing)
        vector(x,path;lo=2,hi=2);text(x[1],path);integer(x[2],path)
        counts!==nothing && need(haskey(counts,x[1]) && x[2]<counts[x[1]],path,"unknown face/cell")
        (String(x[1]),Int(x[2]))
    end
    function point(x,path,vertices=nothing)
        strings(x,2,path)
        vertices!==nothing && need(haskey(vertices,x[1]) && x[2] in vertices[x[1]],path,"unknown vertex")
        x
    end
    function dart(x,path,shape=nothing)
        strings(x,3,path);need(x[2]!=x[3],path,"collapsed dart")
        if shape!==nothing
            need(haskey(shape.edges,x[1]) && M.edge(x[2],x[3]) in shape.edges[x[1]],path,"unknown edge")
        end
        x
    end
    function rows(x,path)
        vector(x,path;lo=3,hi=3);ws=Any[]
        for row in x
            vector(row,path;lo=3,hi=3);w=[rational(q,path) for q in row]
            need(sum(w)==1 && all(>=(0),w),path,"invalid barycentric row");push!(ws,w)
        end
        need(_mt_ratio(ws)>0,path,"folded barycentric triangle");ws
    end
    function mesh(x,path)
        object(x,("surface","patches"),path)
        surface=object(x["surface"],("marking_id","genus","gluings","boundary"),path*".surface")
        text(surface["marking_id"],path);integer(surface["genus"],path)
        vector(surface["boundary"],path;hi=0)
        patches=vector(x["patches"],path*".patches";lo=1,hi=64)
        counts=Dict{String,Int}();vertices=Dict();triangles=Dict();edges=Dict();sides=Dict()
        for p in patches
            object(p,("face","side_order","sides","triangles"),path*".patch";optional=("side_parameters",))
            f=p["face"];text(f,path);need(!haskey(counts,f),path,"duplicate face")
            order=vector(p["side_order"],path;lo=1,hi=256)
            for s in order;text(s,path);end
            need(length(Set(order))==length(order),path,"duplicate side")
            object(p["sides"],order,path*".sides")
            ts=vector(p["triangles"],path*".triangles";lo=1,hi=507)
            counts[f]=length(ts);vs=Set{String}();es=Set()
            for (i,t) in enumerate(ts)
                strings(t,3,path;distinct=true);union!(vs,String.(t));triangles[(f,i-1)]=t
                for (u,v) in M.tedges(t);push!(es,M.edge(u,v));end
            end
            need(length(vs)<=256,path,"native source vertex bound exceeded")
            vertices[f]=vs;edges[f]=es
            for s in order
                chain=vector(p["sides"][s],path;lo=2,hi=257)
                need(all(v->v isa AbstractString && v in vs,chain),path,"unknown side vertex")
                sides[(f,s)]=chain
            end
            ps=get(p,"side_parameters",Dict())
            need(ps isa AbstractDict && issubset(Set(keys(ps)),Set(order)),path,"unknown side parameter key")
            for (s,qs) in ps
                vector(qs,path;lo=length(p["sides"][s]),hi=length(p["sides"][s]))
                q=[rational(v,path) for v in qs]
                need(first(q)==0 && last(q)==1 && all(q[i]<q[i+1] for i in 1:length(q)-1),path,"unordered side parameters")
            end
        end
        need(surface["genus"]<=sum(values(counts)),path,"genus exceeds native face bound")
        for pair in vector(surface["gluings"],path)
            vector(pair,path;lo=2,hi=2)
            for side in pair;strings(side,2,path);need(haskey(sides,Tuple(side)),path,"unknown glued side");end
        end
        (;counts,vertices,triangles,edges,sides,total=sum(values(counts)))
    end
    function dictionary(x,ids,path,validator)
        object(x,ids,path)
        for (id,value) in x;validator(value,path*"["*id*"]");end
    end
    function formats(x,ids,path,shape)
        dictionary(x,ids,path,(word,p)->begin
            text(word,p)
            ds=try M.parse_word(word) catch err;fail(p,sprint(showerror,err));end
            need(length(ds)<=3shape.total,p,"source word exceeds source edge count")
            for d in ds;dart([d.face,d.u,d.v],p,shape);end
        end)
    end
    function seedmap(x,ids,path,counts)
        dictionary(x,ids,path,(seed,p)->cell(seed,p,counts))
    end
    function star_inverse(x,path,curveids,crossids,markids,hasinfinity)
        object(x,("source_mesh","route_edges","crossing_locations","region_seeds","word_formats",
            "source_document","star_document","forward_specification"),path;optional=("infinity_triangle",))
        old=mesh(x["source_mesh"],path*".source_mesh")
        for key in ("source_document","star_document");hash(x[key],path*"."*key);end
        opts=object(x["forward_specification"],("merge_charts",),path*".forward_specification";optional=("subdivision_budget",))
        need(opts["merge_charts"]===true,path,"unexpected star merge policy")
        haskey(opts,"subdivision_budget") && integer(opts["subdivision_budget"],path;lo=3,hi=256)
        for r in vector(x["route_edges"],path;hi=6old.total)
            object(r,("source","target"),path);dart(r["source"],path,old)
            for d in vector(r["target"],path;lo=2,hi=2);dart(d,path);end
        end
        dictionary(x["crossing_locations"],crossids,path,(p,q)->point(p,q,old.vertices))
        seedmap(x["region_seeds"],markids,path,old.counts);formats(x["word_formats"],curveids,path,old)
        need(haskey(x,"infinity_triangle")==hasinfinity,path,"infinity presence differs")
        hasinfinity && cell(x["infinity_triangle"],path,old.counts)
        old
    end

    required=("schema","forward_specification","word","construction","fine_source_mesh","fine_source_hash",
        "source_hash","target_hash","child_cells","new_vertices","curve_images","proof_complex_scope",
        "native_child_cells","native_curve_images","rechart_inverse_specification","proof_crossing_locations",
        "proof_region_seeds","region_seeds","forward_correspondence")
    object(spec,required,"fields";optional=("infinity_triangle","proof_infinity_triangle"))
    need(spec["schema"]=="marked-dehn-word-inverse/1","schema")
    for key in ("fine_source_hash","source_hash","target_hash");hash(spec[key],key);end
    text(spec["proof_complex_scope"],"proof_complex_scope")
    function word(x,path)
        for r in vector(x,path;lo=1)
            object(r,("core","sign"),path);text(r["core"],path)
            need(r["sign"] isa Integer && !(r["sign"] isa Bool) && r["sign"] in (-1,1),path,"invalid factor sign")
        end
    end
    object(spec["forward_specification"],("word",),"forward_specification")
    word(spec["word"],"word");word(spec["forward_specification"]["word"],"forward_specification.word")
    fine=mesh(spec["fine_source_mesh"],"fine_source_mesh")
    native=mesh(Dict("surface"=>current["surface"],"patches"=>current["patches"]),"current")
    curveids=Set(c["id"] for c in current["curves"]);crossids=Set(c["id"] for c in current["crossings"])
    markids=Set(m["id"] for m in get(current,"region_marks",[]))
    for id in union(curveids,crossids,markids);text(id,"current IDs");end
    hasinfinity=_mm_infinity(current)!==nothing
    need(haskey(spec,"infinity_triangle")==hasinfinity && haskey(spec,"proof_infinity_triangle")==hasinfinity,
        "infinity","optional seed presence differs from target")
    hasinfinity && cell(spec["infinity_triangle"],"infinity_triangle",fine.counts)

    rawrows=vector(spec["child_cells"],"child_cells";lo=1,hi=native.total)
    rawcounts=Dict{String,Int}();rawtriangles=Dict();rawvertices=Dict(f=>Set{String}() for f in keys(fine.counts));rawedges=Dict(f=>Set() for f in keys(fine.counts))
    function maprow(r,path,counts=nothing)
        object(r,("child","source","target","shift","vertices","source_weights","target_weights"),path)
        c=cell(r["child"],path,counts);cell(r["source"],path,fine.counts);cell(r["target"],path,fine.counts)
        integer(r["shift"],path;lo=0,hi=0);strings(r["vertices"],3,path;distinct=true)
        rows(r["source_weights"],path);rows(r["target_weights"],path);c
    end
    for r in rawrows
        c=maprow(r,"child_cells")
        need(haskey(fine.counts,c[1]) && c[2]<length(rawrows) && !haskey(rawtriangles,c),"child_cells","invalid/duplicate raw child")
        rawtriangles[c]=r["vertices"];rawcounts[c[1]]=get(rawcounts,c[1],0)+1
        union!(rawvertices[c[1]],String.(r["vertices"]))
        for (a,b) in M.tedges(r["vertices"]);push!(rawedges[c[1]],M.edge(a,b));end
    end
    need(Set(keys(rawcounts))==Set(keys(fine.counts)) &&
        all(haskey(rawtriangles,(f,i-1)) for (f,n) in rawcounts for i in 1:n),"child_cells","noncontiguous raw cell indices")
    raw=(;counts=rawcounts,vertices=rawvertices,triangles=rawtriangles,edges=rawedges,total=length(rawrows))
    function occurrence(x,path,shape)
        vector(x,path;lo=4,hi=4);c=cell(x[1:2],path,shape.counts)
        strings(x[3:4],2,path;distinct=true)
        need(all(v->v in shape.triangles[c],x[3:4]),path,"unknown occurrence vertex");x
    end
    function mapcertificate(x,path,n)
        object(x,("scope","source_partition","target_partition","forward_continuity","inverse_continuity"),path)
        text(x["scope"],path)
        for key in ("source_partition","target_partition")
            seen=Set();ids=Set()
            for r in vector(x[key],path*"."*key;lo=fine.total,hi=fine.total)
                object(r,("parent","atoms","area_ratio"),path);p=cell(r["parent"],path,fine.counts)
                need(!(p in seen),path,"duplicate partition parent");push!(seen,p)
                for i in vector(r["atoms"],path;lo=1,hi=n)
                    integer(i,path;lo=1,hi=n);need(!(i in ids),path,"duplicate partition atom");push!(ids,i)
                end
                need(rational(r["area_ratio"],path)==1,path,"invalid parent area ratio")
            end
            need(length(ids)==n,path,"missing partition atom")
        end
        for key in ("forward_continuity","inverse_continuity")
            c=object(x[key],("vertices","paired_edges","edge_index_convention","fixed_edge_breakpoints"),path)
            integer(c["vertices"],path;lo=1,hi=3n);integer(c["paired_edges"],path;lo=1,hi=3n)
            text(c["edge_index_convention"],path)
            for r in vector(c["fixed_edge_breakpoints"],path;hi=3fine.total)
                object(r,("edge","parameters"),path);occurrence(r["edge"],path,fine)
                qs=[rational(q,path) for q in vector(r["parameters"],path;lo=2,hi=3n)]
                need(first(qs)==0 && last(qs)==1 && all(qs[i]<qs[i+1] for i in 1:length(qs)-1),path,"invalid seam breakpoints")
            end
        end
    end
    function maprows(x,path,shape)
        vector(x,path;lo=shape.total,hi=shape.total);seen=Set()
        for r in x
            c=maprow(r,path,shape.counts);need(!(c in seen),path,"duplicate child");push!(seen,c)
            need(r["vertices"]==shape.triangles[c],path,"child vertices disagree with actual target")
        end
    end
    maprows(spec["native_child_cells"],"native_child_cells",native)
    seenvertices=Set()
    for r in vector(spec["new_vertices"],"new_vertices";hi=3raw.total)
        object(r,("vertex","ancestry"),"new_vertices");point(r["vertex"],"new_vertices",raw.vertices)
        key=Tuple(r["vertex"]);need(!(key in seenvertices),"new_vertices","duplicate vertex");push!(seenvertices,key)
        a=r["ancestry"];need(a isa AbstractDict && haskey(a,"kind"),"new_vertices.ancestry")
        if a["kind"]=="edge"
            object(a,("kind","edge","parameter"),"new_vertices.edge");strings(a["edge"],2,"new_vertices.edge";distinct=true)
            need(M.edge(a["edge"]...) in fine.edges[key[1]],"new_vertices.edge","unknown old edge")
            q=rational(a["parameter"],"new_vertices.parameter");need(0<q<1,"new_vertices.parameter")
        else
            object(a,("kind","cell","weights"),"new_vertices.triangle");need(a["kind"]=="triangle","new_vertices.kind")
            cell(a["cell"],"new_vertices.cell",fine.counts);vector(a["weights"],"new_vertices.weights";lo=3,hi=3)
            w=[rational(v,"new_vertices.weights") for v in a["weights"]]
            need(sum(w)==1 && all(>(0),w),"new_vertices.weights")
        end
    end
    need(all(v in fine.vertices[f] || (f,v) in seenvertices for (f,vs) in raw.vertices for v in vs),
        "new_vertices","missing actual raw vertex ancestry")
    function intervals(x,path,maximum)
        for r in vector(x,path;hi=maximum)
            object(r,("from","to","dart"),path);a=rational(r["from"],path);b=rational(r["to"],path)
            need(0<=a<b<=1,path,"invalid edge interval");dart(r["dart"],path)
        end
    end
    function images(x,path,shape;expanded=false)
        vector(x,path;lo=length(curveids),hi=length(curveids));seen=Set()
        for r in x
            object(r,("curve","source_word_format","edge_images"),path);text(r["curve"],path)
            need(r["curve"] in curveids && !(r["curve"] in seen),path,"unknown/duplicate curve");push!(seen,r["curve"])
            text(r["source_word_format"],path)
            ds=try M.parse_word(r["source_word_format"]) catch err;fail(path,sprint(showerror,err));end
            pieces=vector(r["edge_images"],path;lo=1,hi=3fine.total)
            need(length(pieces)==length(ds),path,"source edge count differs from lexical word")
            for (p,d) in zip(pieces,ds)
                fields=expanded ? ("source","target","partition","proof_target","proof_partition") : ("source","target","partition")
                object(p,fields,path);dart(p["source"],path,fine)
                need(p["source"]==[d.face,d.u,d.v],path,"source word/dart mismatch")
                for t in vector(p["target"],path;lo=1,hi=3shape.total);dart(t,path,shape);end
                intervals(p["partition"],path,length(p["target"]))
                if expanded
                    for t in vector(p["proof_target"],path;lo=1,hi=3raw.total);dart(t,path,raw);end
                    intervals(p["proof_partition"],path,length(p["proof_target"]))
                end
            end
        end
    end
    images(spec["curve_images"],"curve_images",raw)
    images(spec["native_curve_images"],"native_curve_images",native;expanded=spec["rechart_inverse_specification"]!==nothing)
    dictionary(spec["proof_crossing_locations"],crossids,"proof_crossing_locations",(p,q)->point(p,q,raw.vertices))
    seedmap(spec["proof_region_seeds"],markids,"proof_region_seeds",raw.counts)
    seedmap(spec["region_seeds"],markids,"region_seeds",fine.counts)
    hasinfinity && cell(spec["proof_infinity_triangle"],"proof_infinity_triangle",raw.counts)

    history=vector(spec["construction"],"construction";lo=2,hi=2)
    witnessfields=("schema","kind","source_document","target_document","source_marking","target_marking",
        "specification","inverse_specification","correspondence","inverse_kind")
    for (i,w) in enumerate(history)
        object(w,witnessfields,"construction[$i]")
        need(w["schema"]=="marked-move/1" && w["kind"]==(i==1 ? "dual_star" : "barycentric_subdivision") &&
            w["inverse_kind"]==(i==1 ? "dual_star_collapse" : "barycentric_collapse"),"construction[$i].kind")
        for key in ("source_document","target_document","source_marking","target_marking");hash(w[key],"construction[$i]."*key);end
        need(w["specification"] isa AbstractDict && w["inverse_specification"] isa AbstractDict && w["correspondence"] isa AbstractDict,
            "construction[$i]","witness payloads must be objects")
        common=("cell_semantics","vertices","introduced_vertices","edges","triangles","regions","geometric_regions","curves","crossings")
        additional=i==1 ? ("chart_merges","chart_correspondence","child_barycentric_rows") : ()
        wc=object(w["correspondence"],(common...,additional...),"construction[$i].correspondence")
        text(wc["cell_semantics"],"construction[$i].cell_semantics")
        for key in common[2:end]
            for r in vector(wc[key],"construction[$i]."*key;hi=6fine.total)
                need(r isa AbstractDict,"construction[$i]."*key,"expected generated proof row")
            end
        end
        if i==1
            for key in ("chart_merges","child_barycentric_rows")
                for r in vector(wc[key],"construction[$i]."*key;hi=6fine.total)
                    need(r isa AbstractDict,"construction[$i]."*key,"expected generated proof row")
                end
            end
            need(wc["chart_correspondence"] isa AbstractDict,"construction[$i].chart_correspondence")
        end
    end
    object(history[1]["specification"],("merge_charts","subdivision_budget"),"construction[1].specification")
    need(history[1]["specification"]["merge_charts"]===true &&
        history[1]["specification"]["subdivision_budget"] isa Integer && !(history[1]["specification"]["subdivision_budget"] isa Bool) &&
        history[1]["specification"]["subdivision_budget"]==256,
        "construction[1].specification","unexpected preparation options")
    object(history[2]["specification"],(),"construction[2].specification")
    star_inverse(history[1]["inverse_specification"],"construction[1].inverse",curveids,crossids,markids,hasinfinity)
    need(history[1]["inverse_specification"]["forward_specification"]==history[1]["specification"],
        "construction[1].inverse.forward_specification","copied star options differ")
    b=object(history[2]["inverse_specification"],("patches","coarse_document","refined_document","region_seeds","word_formats"),
        "construction[2].inverse";optional=("infinity_triangle",))
    hash(b["coarse_document"],"construction[2].coarse_document");hash(b["refined_document"],"construction[2].refined_document")
    patches=vector(b["patches"],"construction[2].patches";lo=length(fine.counts),hi=length(fine.counts));coarsecounts=Dict();seen=Set()
    for p in patches
        object(p,("face","edges","triangles","side_parameters"),"construction[2].patch")
        f=p["face"];text(f,"construction[2].face");need(haskey(fine.counts,f) && !(f in seen),"construction[2].face");push!(seen,f)
        need(fine.counts[f]%6==0,"construction[2].triangles","fine flags are not divisible by six")
        n=div(fine.counts[f],6);coarsecounts[f]=n
        for r in vector(p["edges"],"construction[2].edges";hi=3n)
            object(r,("edge","midpoint"),"construction[2].edge");strings(r["edge"],2,"construction[2].edge";distinct=true);text(r["midpoint"],"construction[2].midpoint")
        end
        for r in vector(p["triangles"],"construction[2].triangles";lo=n,hi=n)
            object(r,("corners","center"),"construction[2].triangle");strings(r["corners"],3,"construction[2].corners";distinct=true);text(r["center"],"construction[2].center")
        end
        ps=p["side_parameters"]
        need(ps===nothing || ps isa AbstractDict,"construction[2].side_parameters")
        if ps!==nothing
            for (s,qs) in ps
                text(s,"construction[2].side_parameters");vector(qs,"construction[2].side_parameters";lo=2,hi=256)
                for q in qs;rational(q,"construction[2].side_parameters");end
            end
        end
    end
    seedmap(b["region_seeds"],markids,"construction[2].region_seeds",coarsecounts)
    dictionary(b["word_formats"],curveids,"construction[2].word_formats",(w,p)->begin
        text(w,p);try M.parse_word(w) catch err;fail(p,sprint(showerror,err));end
    end)
    need(haskey(b,"infinity_triangle")==hasinfinity,"construction[2].infinity")
    hasinfinity && cell(b["infinity_triangle"],"construction[2].infinity",coarsecounts)

    # Non-consumed proof annotations remain JSON objects/row arrays here;
    # their EXACT contents are compared with independent full regeneration.
    # Do not bind stage core order to word labels in a SHAPE guard: coherent
    # label tampering must reach the actual full-word regeneration test.
    corrfields=("word","word_convention","action_scope","construction","elementary_maps","composition_stages",
        "background_marking","child_cells","raw_proof_vertices","vertex_scope","rechart","curve_images","regions",
        "geometric_regions","total_map_certificate","old_basis_action","basis_scope","basis_source_hash",
        "basis_geometry_hash","basis_marking_hash","crossings")
    corr=object(spec["forward_correspondence"],corrfields,"forward_correspondence")
    word(corr["word"],"forward_correspondence.word")
    for key in ("word_convention","action_scope","background_marking","vertex_scope","basis_scope");text(corr[key],"forward_correspondence."*key);end
    for key in ("basis_source_hash","basis_geometry_hash","basis_marking_hash");hash(corr[key],"forward_correspondence."*key);end
    vector(corr["construction"],"forward_correspondence.construction";lo=2,hi=2)
    need(corr["construction"]==history,"forward_correspondence.construction","copied history differs")
    maprows(corr["child_cells"],"forward_correspondence.child_cells",native)
    images(corr["curve_images"],"forward_correspondence.curve_images",native;expanded=spec["rechart_inverse_specification"]!==nothing)
    vector(corr["raw_proof_vertices"],"forward_correspondence.raw_proof_vertices";hi=3raw.total)
    need(corr["raw_proof_vertices"]==spec["new_vertices"],"forward_correspondence.raw_proof_vertices","copied vertex rows differ")
    for key in ("regions","geometric_regions")
        rr=vector(corr[key],"forward_correspondence."*key;lo=1,hi=fine.total)
        for row in rr
            object(row,("source_region","target_region","source_triangles","target_triangles"),"forward_correspondence."*key)
            integer(row["source_region"],key;lo=1,hi=length(rr));integer(row["target_region"],key;lo=1,hi=length(rr))
            for c in vector(row["source_triangles"],key;lo=1,hi=fine.total);cell(c,key,fine.counts);end
            for c in vector(row["target_triangles"],key;lo=1,hi=native.total);cell(c,key,native.counts);end
        end
    end
    mapcertificate(corr["total_map_certificate"],"forward_correspondence.total_map_certificate",raw.total)
    dimension=2Int(current["surface"]["genus"])
    function integers(x,n,path)
        for v in vector(x,path;lo=n,hi=n);need(v isa Integer && !(v isa Bool),path,"noninteger coordinate");end
    end
    function generatorword(x,path)
        for pair in vector(x,path)
            vector(pair,path;lo=2,hi=2);integer(pair[1],path;lo=1,hi=dimension)
            need(pair[2] isa Integer && !(pair[2] isa Bool) && pair[2] in (-1,1),path,"invalid generator sign")
        end
    end
    for r in vector(corr["old_basis_action"],"old_basis_action";lo=length(curveids),hi=length(curveids))
        object(r,("curve","source_word","image_word_old_basis","source_exponents","image_exponents","expected_exponents",
            "edge_retractions","word_semantics","original_factor_intersections"),"old_basis_action.row")
        need(r["curve"] in curveids,"old_basis_action.curve");text(r["word_semantics"],"old_basis_action.word_semantics")
        generatorword(r["source_word"],"old_basis_action.source_word");generatorword(r["image_word_old_basis"],"old_basis_action.image_word")
        for k in ("source_exponents","image_exponents","expected_exponents");integers(r[k],dimension,"old_basis_action."*k);end
        integers(r["original_factor_intersections"],length(spec["word"]),"old_basis_action.intersections")
        for e in vector(r["edge_retractions"],"old_basis_action.edge_retractions";hi=3native.total)
            object(e,("image_dart","old_triangle","anchors","old_edge_word"),"old_basis_action.edge_retraction")
            dart(e["image_dart"],"old_basis_action.image_dart",native)
            c=cell(e["old_triangle"],"old_basis_action.old_triangle",fine.counts);strings(e["anchors"],2,"old_basis_action.anchors")
            need(all(v->v in fine.triangles[c],e["anchors"]),"old_basis_action.anchors")
            generatorword(e["old_edge_word"],"old_basis_action.old_edge_word")
        end
    end
    for r in vector(corr["crossings"],"forward_correspondence.crossings";lo=length(crossids),hi=length(crossids))
        object(r,("source","target","at","sign"),"forward_correspondence.crossing")
        need(r["source"] in crossids && r["target"] in crossids,"forward_correspondence.crossing_id")
        point(r["at"],"forward_correspondence.crossing_at",native.vertices)
        need(r["sign"] isa Integer && !(r["sign"] isa Bool) && r["sign"] in (-1,1),"forward_correspondence.crossing_sign")
    end
    need(corr["rechart"]===nothing || corr["rechart"] isa AbstractDict,"forward_correspondence.rechart")
    need((corr["rechart"]===nothing)==(spec["rechart_inverse_specification"]===nothing),"forward_correspondence.rechart","inverse/correspondence presence differs")
    if corr["rechart"]!==nothing
        need(spec["rechart_inverse_specification"] isa AbstractDict,"rechart_inverse_specification")
        native_schema=get(spec["rechart_inverse_specification"],"schema",nothing)
        proper=native_schema in ("internal-triangle-partition/1","internal-triangle-partition/2")
        policy2=native_schema=="internal-triangle-partition/2"
        common=("cell_semantics","edge_images","vertex_images","cell_children","child_barycentric_rows","canonical_edges","edge_occurrences")
        extras=proper ? ("route_domain","proper_cuts") : ("first_chart_merges","first_chart_correspondence","second_chart_merges")
        policy2 && (extras=(extras...,"partition_policy"))
        rc=object(corr["rechart"],(common...,extras...),"forward_correspondence.rechart")
        policy2 && need(rc["partition_policy"]=="chart_demand/1" &&
            get(spec["rechart_inverse_specification"],"partition_policy",nothing)==rc["partition_policy"],
            "rechart.partition_policy","correspondence/inverse policy differs")
        text(rc["cell_semantics"],"rechart.cell_semantics")
        for r in vector(rc["edge_images"],"rechart.edge_images";hi=6raw.total)
            object(r,("source","target"),"rechart.edge_image");dart(r["source"],"rechart.edge_source",raw)
            for d in vector(r["target"],"rechart.edge_target";lo=1,hi=3native.total);dart(d,"rechart.edge_target",native);end
        end
        for r in vector(rc["vertex_images"],"rechart.vertex_images";hi=3raw.total)
            object(r,("source","target"),"rechart.vertex_image");point(r["source"],"rechart.vertex_source",raw.vertices);point(r["target"],"rechart.vertex_target",native.vertices)
        end
        for r in vector(rc["cell_children"],"rechart.cell_children";lo=raw.total,hi=raw.total)
            object(r,("source","target"),"rechart.cell_children.row");cell(r["source"],"rechart.cell_source",raw.counts)
            for c in vector(r["target"],"rechart.cell_targets";lo=1,hi=native.total);cell(c,"rechart.cell_target",native.counts);end
        end
        for r in vector(rc["child_barycentric_rows"],"rechart.child_barycentric_rows";lo=native.total,hi=native.total)
            object(r,("source","target","vertices","weights"),"rechart.affine_row")
            cell(r["source"],"rechart.affine_source",raw.counts);cell(r["target"],"rechart.affine_target",native.counts)
            strings(r["vertices"],3,"rechart.affine_vertices";distinct=true);rows(r["weights"],"rechart.affine_weights")
        end
        for r in vector(rc["canonical_edges"],"rechart.canonical_edges";hi=3raw.total)
            object(r,("edge","target"),"rechart.canonical_edge");occurrence(r["edge"],"rechart.canonical_occurrence",raw)
            for d in vector(r["target"],"rechart.canonical_target";lo=1,hi=3native.total);dart(d,"rechart.canonical_dart",native);end
        end
        for r in vector(rc["edge_occurrences"],"rechart.edge_occurrences";lo=3raw.total,hi=3raw.total)
            object(r,("source","edge","orientation"),"rechart.edge_occurrence")
            occurrence(r["source"],"rechart.source_occurrence",raw);occurrence(r["edge"],"rechart.edge_occurrence",raw)
            need(r["orientation"] isa Integer && !(r["orientation"] isa Bool) && r["orientation"] in (-1,1),"rechart.edge_orientation")
        end
        if proper
            text(rc["route_domain"],"rechart.route_domain")
            for r in vector(rc["proper_cuts"],"rechart.proper_cuts";hi=native.total)
                need(r isa AbstractDict,"rechart.proper_cuts","expected generated cut row")
            end
        else
            for key in ("first_chart_merges","second_chart_merges")
                for r in vector(rc[key],"rechart."*key;hi=3raw.total)
                    need(r isa AbstractDict,"rechart."*key,"expected generated merge row")
                end
            end
            need(rc["first_chart_correspondence"] isa AbstractDict,"rechart.first_chart_correspondence")
        end
    end
    factors=vector(corr["elementary_maps"],"elementary_maps";lo=1,hi=length(spec["word"]))
    for factor in factors
        object(factor,("core","sign","atoms","product_triangles","product_certificate","map_certificate"),"elementary_maps.factor")
        text(factor["core"],"elementary_maps.core");need(factor["sign"] isa Integer && !(factor["sign"] isa Bool) && factor["sign"] in (-1,1),"elementary_maps.sign")
        for atom in vector(factor["atoms"],"elementary_maps.atoms";lo=1,hi=raw.total)
            object(atom,("source","target","source_weights","target_weights"),"elementary_maps.atom")
            cell(atom["source"],"elementary_maps.source",fine.counts);cell(atom["target"],"elementary_maps.target",fine.counts)
            rows(atom["source_weights"],"elementary_maps.source_weights");rows(atom["target_weights"],"elementary_maps.target_weights")
        end
        for triangle in vector(factor["product_triangles"],"elementary_maps.product_triangles";lo=1,hi=fine.total)
            object(triangle,("cell","points"),"elementary_maps.product_triangle");cell(triangle["cell"],"elementary_maps.product_cell",fine.counts)
            for p in vector(triangle["points"],"elementary_maps.product_points";lo=3,hi=3)
                vector(p,"elementary_maps.point";lo=2,hi=2);for q in p;rational(q,"elementary_maps.rational");end
            end
        end
        pc=object(factor["product_certificate"],("construction","core","core_order","period","left_height","right_height",
            "fans","edge_transitions","frontiers","area","core_winding","triangle_coordinates","inverse_support"),"elementary_maps.product_certificate")
        for k in ("construction","core","inverse_support");text(pc[k],"elementary_maps.product_certificate."*k);end
        need(pc["core"]==factor["core"] && pc["period"] isa Integer && pc["period"]==1 && !(pc["period"] isa Bool) &&
            pc["left_height"] isa Integer && pc["left_height"]==1 && !(pc["left_height"] isa Bool) &&
            pc["right_height"] isa Integer && pc["right_height"]==-1 && !(pc["right_height"] isa Bool),"elementary_maps.product_heights")
        core=vector(pc["core_order"],"elementary_maps.core_order";lo=3,hi=3fine.total)
        for d in core;dart(d,"elementary_maps.core_dart");point(d[1:2],"elementary_maps.core_vertex",fine.vertices);end
        need(rational(pc["area"],"elementary_maps.area")==2 && rational(pc["core_winding"],"elementary_maps.winding")==1,"elementary_maps.product_normalization")
        for r in vector(pc["fans"],"elementary_maps.fans";lo=2length(core),hi=2length(core))
            object(r,("core_vertex","index","side","frontier_edges","link_path","triangles"),"elementary_maps.fan")
            point(r["core_vertex"],"elementary_maps.fan_vertex",fine.vertices);integer(r["index"],"elementary_maps.fan_index";hi=length(core)-1)
            need(r["side"] isa Integer && !(r["side"] isa Bool) && r["side"] in (-1,1),"elementary_maps.fan_side")
            integer(r["frontier_edges"],"elementary_maps.frontier_edges";lo=1,hi=3fine.total)
            for v in vector(r["link_path"],"elementary_maps.link_path";lo=4,hi=3fine.total);point(v,"elementary_maps.link_vertex",fine.vertices);end
            for c in vector(r["triangles"],"elementary_maps.fan_triangles";lo=1,hi=fine.total);cell(c,"elementary_maps.fan_cell",fine.counts);end
        end
        for r in vector(pc["edge_transitions"],"elementary_maps.edge_transitions";hi=3fine.total)
            object(r,("source","target","u_shift"),"elementary_maps.edge_transition")
            occurrence(r["source"],"elementary_maps.transition_source",fine);occurrence(r["target"],"elementary_maps.transition_target",fine)
            need(r["u_shift"] isa AbstractString && occursin(r"^-?[0-9]+$",r["u_shift"]),"elementary_maps.u_shift")
        end
        for r in vector(pc["frontiers"],"elementary_maps.frontiers";lo=2,hi=2)
            object(r,("side","cycle"),"elementary_maps.frontier");need(r["side"] isa Integer && !(r["side"] isa Bool) && r["side"] in (-1,1),"elementary_maps.frontier_side")
            for o in vector(r["cycle"],"elementary_maps.frontier_cycle";lo=1,hi=3fine.total);occurrence(o,"elementary_maps.frontier_occurrence",fine);end
        end
        coordinates=vector(pc["triangle_coordinates"],"elementary_maps.triangle_coordinates";lo=length(factor["product_triangles"]),hi=length(factor["product_triangles"]))
        for (a,b) in zip(coordinates,factor["product_triangles"])
            object(a,("id","points"),"elementary_maps.triangle_coordinate")
            need(a["id"]==b["cell"] && a["points"]==b["points"],"elementary_maps.triangle_coordinates","copied product rows differ")
        end
        mapcertificate(factor["map_certificate"],"elementary_maps.map_certificate",length(factor["atoms"]))
    end
    stages=vector(corr["composition_stages"],"composition_stages";lo=length(spec["word"]),hi=length(spec["word"]))
    previous=0
    for (i,stage) in enumerate(stages)
        object(stage,("step","core","sign","factor","composition","output_atoms","expected_vectors"),"composition_stages.row")
        integer(stage["step"],"composition_stages.step";lo=1,hi=length(stages));text(stage["core"],"composition_stages.core")
        need(stage["sign"] isa Integer && !(stage["sign"] isa Bool) && stage["sign"] in (-1,1),"composition_stages.sign")
        integer(stage["factor"],"composition_stages.factor";lo=1,hi=length(factors));integer(stage["output_atoms"],"composition_stages.output_atoms";lo=1,hi=raw.total)
        comp=stage["composition"]
        if i==1
            object(comp,("kind","factor"),"composition_stages.first")
            text(comp["kind"],"composition_stages.first.kind");integer(comp["factor"],"composition_stages.first.factor";lo=1,hi=length(factors))
        else
            object(comp,("composition","left","right","pieces","result"),"composition_stages.composition")
            text(comp["composition"],"composition_stages.composition")
            nright=length(factors[stage["factor"]]["atoms"])
            mapcertificate(comp["left"],"composition_stages.left",previous)
            mapcertificate(comp["right"],"composition_stages.right",nright)
            mapcertificate(comp["result"],"composition_stages.result",stage["output_atoms"])
            for r in vector(comp["pieces"],"composition_stages.pieces";lo=stage["output_atoms"],hi=stage["output_atoms"])
                object(r,("atom","middle_parent","left_atom","right_atom","intersection_vertices","selected_vertices"),"composition_stages.piece")
                integer(r["atom"],"composition_stages.atom";lo=1,hi=stage["output_atoms"])
                cell(r["middle_parent"],"composition_stages.middle_parent",fine.counts)
                integer(r["left_atom"],"composition_stages.left_atom";lo=1,hi=previous);integer(r["right_atom"],"composition_stages.right_atom";lo=1,hi=nright)
                points=vector(r["intersection_vertices"],"composition_stages.intersection_vertices";lo=3,hi=stage["output_atoms"]+2)
                for p in points;vector(p,"composition_stages.point";lo=2,hi=2);for q in p;rational(q,"composition_stages.rational");end;end
                for j in vector(r["selected_vertices"],"composition_stages.selected_vertices";lo=3,hi=3);integer(j,"composition_stages.selected_vertex";lo=1,hi=length(points));end
            end
        end
        for row in vector(stage["expected_vectors"],"composition_stages.expected_vectors";lo=length(curveids),hi=length(curveids))
            object(row,("curve","exponents"),"composition_stages.vector");need(row["curve"] in curveids,"composition_stages.curve")
            for x in vector(row["exponents"],"composition_stages.exponents";lo=2current["surface"]["genus"],hi=2current["surface"]["genus"])
                need(x isa Integer && !(x isa Bool),"composition_stages.exponent")
            end
        end
        previous=Int(stage["output_atoms"])
    end
    # Existing selective schema owns its detailed raw chart/template guards.
    rechart=spec["rechart_inverse_specification"]
    if rechart!==nothing
        need(rechart isa AbstractDict && haskey(rechart,"schema"),"rechart_inverse_specification")
        if rechart["schema"] in ("internal-triangle-partition/1","internal-triangle-partition/2")
            fields=("schema","strategy","raw_geometry_hash","public_geometry_hash","raw_surface","raw_charts",
                "triangle_ancestry","opposite_edges","affine_children","proper_cuts","correspondence_hash")
            if rechart["schema"]=="internal-triangle-partition/2"
                fields=(fields...,"partition_policy")
            end
            object(rechart,fields,"rechart_inverse_specification")
            rechart["schema"]=="internal-triangle-partition/2" &&
                need(rechart["partition_policy"]=="chart_demand/1","rechart.partition_policy")
            vector(rechart["raw_charts"],"rechart.raw_charts";lo=length(raw.counts),hi=length(raw.counts))
            vector(rechart["triangle_ancestry"],"rechart.triangle_ancestry";lo=raw.total,hi=raw.total)
            vector(rechart["affine_children"],"rechart.affine_children";lo=native.total,hi=native.total)
            vector(rechart["opposite_edges"],"rechart.opposite_edges";hi=3raw.total)
            vector(rechart["proper_cuts"],"rechart.proper_cuts";hi=native.total)
            faces=Set()
            for p in rechart["raw_charts"]
                object(p,("face","side_order","sides","triangle_count"),"rechart.raw_chart";optional=("side_parameters",))
                text(p["face"],"rechart.raw_chart.face")
                need(haskey(raw.counts,p["face"]) && !(p["face"] in faces),"rechart.raw_chart.face");push!(faces,p["face"])
                integer(p["triangle_count"],"rechart.raw_chart.triangle_count";lo=1,hi=raw.total)
                need(p["triangle_count"]==raw.counts[p["face"]],"rechart.raw_chart.triangle_count")
                need(p["sides"] isa AbstractDict,"rechart.raw_chart.sides")
                for chain in values(p["sides"]);vector(chain,"rechart.raw_chart.side";lo=2,hi=3raw.total);end
            end
            try _mrp_inverse_schema(rechart) catch err;fail("rechart_inverse_specification",sprint(showerror,err));end
        else
            object(rechart,("schema","raw_geometry_hash","public_geometry_hash","raw_surface","raw_charts",
                "triangle_ancestry","opposite_edges","star_inverse","affine_children","correspondence_hash"),"rechart_inverse_specification")
            need(rechart["schema"]=="internal-triangle-rechart/1","rechart_inverse_specification.schema")
            for key in ("raw_geometry_hash","public_geometry_hash","correspondence_hash");hash(rechart[key],"rechart."*key);end
            object(rechart["raw_surface"],("marking_id","genus","gluings","boundary"),"rechart.raw_surface")
            rs=rechart["raw_surface"];text(rs["marking_id"],"rechart.raw_surface.marking_id")
            integer(rs["genus"],"rechart.raw_surface.genus");vector(rs["boundary"],"rechart.raw_surface.boundary";hi=0)
            vector(rechart["raw_charts"],"rechart.raw_charts";lo=length(raw.counts),hi=length(raw.counts))
            chartfaces=Set();chartsides=Set()
            for p in rechart["raw_charts"]
                object(p,("face","side_order","sides","triangle_count"),"rechart.raw_chart";optional=("side_parameters",))
                text(p["face"],"rechart.raw_chart.face");need(haskey(raw.counts,p["face"]) && !(p["face"] in chartfaces),"rechart.raw_chart.face")
                push!(chartfaces,p["face"])
                integer(p["triangle_count"],"rechart.raw_chart.triangle_count";lo=1,hi=raw.total)
                need(p["triangle_count"]==raw.counts[p["face"]],"rechart.raw_chart.triangle_count")
                order=vector(p["side_order"],"rechart.raw_chart.side_order";lo=1,hi=3raw.total)
                for s in order;text(s,"rechart.raw_chart.side_order");end
                need(length(Set(order))==length(order),"rechart.raw_chart.side_order","duplicate side")
                object(p["sides"],order,"rechart.raw_chart.sides")
                for (side,chain) in p["sides"]
                    push!(chartsides,(p["face"],side));vector(chain,"rechart.raw_chart.side";lo=2,hi=3raw.total)
                    need(all(v->v isa AbstractString && v in raw.vertices[p["face"]],chain),"rechart.raw_chart.vertex")
                end
                need(!haskey(p,"side_parameters") || p["side_parameters"] isa AbstractDict,"rechart.raw_chart.parameters")
                for (side,qs) in get(p,"side_parameters",Dict())
                    need(side in order,"rechart.raw_chart.parameters","unknown side")
                    vector(qs,"rechart.raw_chart.parameters";lo=length(p["sides"][side]),hi=length(p["sides"][side]))
                    q=[rational(v,"rechart.raw_chart.parameter") for v in qs]
                    need(first(q)==0 && last(q)==1 && all(q[i]<q[i+1] for i in 1:length(q)-1),"rechart.raw_chart.parameters")
                end
            end
            for pair in vector(rs["gluings"],"rechart.raw_surface.gluings")
                vector(pair,"rechart.raw_surface.gluing";lo=2,hi=2)
                for side in pair;strings(side,2,"rechart.raw_surface.side");need(Tuple(side) in chartsides,"rechart.raw_surface.side");end
            end
            atlas=star_inverse(rechart["star_inverse"],"rechart.star_inverse",Set{String}(),Set{String}(),Set{String}(),false)
            vector(rechart["triangle_ancestry"],"rechart.triangle_ancestry";lo=raw.total,hi=raw.total)
            for r in rechart["triangle_ancestry"]
                object(r,("source","source_vertices","temporary_face","merged_cell","merged_vertices"),"rechart.triangle_ancestry.row")
                cell(r["source"],"rechart.source",raw.counts);strings(r["source_vertices"],3,"rechart.source_vertices";distinct=true)
                text(r["temporary_face"],"rechart.temporary_face");cell(r["merged_cell"],"rechart.merged_cell",atlas.counts)
                strings(r["merged_vertices"],3,"rechart.merged_vertices";distinct=true)
            end
            for r in vector(rechart["opposite_edges"],"rechart.opposite_edges";lo=1,hi=3raw.total)
                object(r,("source","target"),"rechart.opposite_edges.row")
                for key in ("source","target")
                    o=r[key];vector(o,"rechart.occurrence";lo=4,hi=4);c=cell(o[1:2],"rechart.occurrence",raw.counts)
                    strings(o[3:4],2,"rechart.occurrence";distinct=true)
                    need(all(v->v in raw.triangles[c],o[3:4]),"rechart.occurrence.vertices")
                end
            end
            for r in vector(rechart["affine_children"],"rechart.affine_children";lo=native.total,hi=native.total)
                object(r,("source","target","vertices","weights"),"rechart.affine_child")
                cell(r["source"],"rechart.affine_source",raw.counts);cell(r["target"],"rechart.affine_target",native.counts)
                strings(r["vertices"],3,"rechart.affine_vertices";distinct=true);rows(r["weights"],"rechart.affine_weights")
            end
        end
    end
    spec
end

"""Invert the ACTUAL generated total map, then regenerate the entire word.

No reversed negative-word representative is applied to the current target.
Actual current geometry/features are unrecharted, SAME TOTAL S/T is pulled
back, protected full regions are checked, actual preparation is collapsed,
and every original annulus/factor/composition/native target is regenerated.
"""
function marked_inverse_dehn_twist_word(input::AbstractDict,specification::AbstractDict)
    current=_marked_document(input);spec=_mtwi_inverse_schema(specification,current)
    try
        _mtw_checked_inverse_word(current,spec)
    catch err
        # This NEW public boundary reports malformed references as contextual
        # ArgumentError; it never converts failure into an accepted document.
        # Positive inverse tests still require every geometric/replay check.
        if err isa KeyError || err isa BoundsError || err isa MethodError || err isa MarkedSurfaceDiagrams.DiagramError
            throw(ArgumentError("word inverse reconstruction: "*sprint(showerror,err)))
        end
        rethrow()
    end
end

function _mtw_checked_inverse_word(current,spec)
    forward=_mtw_word_specification(spec["forward_specification"])
    _mm_need(spec["word"]==forward["word"] && spec["forward_correspondence"]["word"]==forward["word"],
        "recorded RIGHT word order/signs disagree")
    pulled=_mtw_restore_total_source_candidate(current,spec)
    candidate=pulled.candidate;history=spec["construction"]
    _mm_need(!isempty(history),"generated word has no actual preparation history")
    original=marked_collapse_subdivision(candidate,last(history)["inverse_specification"]).document
    for w in reverse(history[1:end-1])
        original=replay_marked_move(original,inverse_marked_move(w)).document
    end
    _mm_need(_mm_hash(original)==spec["source_hash"],"actual word preparation collapse does not restore the source")
    regenerated=marked_dehn_twist_word(original,forward)
    _mm_need(regenerated.document==current &&
        _canonical_json(regenerated.witness["inverse_specification"])==_canonical_json(spec) &&
        _canonical_json(regenerated.witness["correspondence"])==_canonical_json(spec["forward_correspondence"]),
        "full factor/RIGHT-word regeneration does not reproduce actual target and inverse/correspondence")
    correspondence=Dict("word_convention"=>spec["forward_correspondence"]["word_convention"],
        "inverse_map"=>"actual native to raw to SAME TOTAL S/T to prepared source; real preparation collapse; full factor/word regeneration",
        "raw_child_cells"=>spec["child_cells"],"raw_cell_scope"=>spec["proof_complex_scope"],
        "native_child_cells"=>spec["native_child_cells"],
        "native_cell_scope"=>"FORWARD target native children; old prepared source/target parents retained",
        "regions"=>pulled.regions,"geometric_regions"=>pulled.geometric_regions,
        "region_direction"=>"prepared source to actual forward target, independently checked during inverse",
        "inverse_preparation"=>history,"regenerated_forward_correspondence"=>regenerated.witness["correspondence"])
    witness=_mm_witness("generated_twist_word_inverse",current,original,spec,forward,correspondence)
    witness["inverse_kind"]="dehn_twist_word"
    (;document=original,witness)
end
