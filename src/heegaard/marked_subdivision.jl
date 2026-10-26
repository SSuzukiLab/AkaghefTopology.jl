"""Identity marked transport by barycentric subdivision of existing disk patches.

Source: vendor MarkedSurfaceDiagrams `patch_info`, `side_parameters`, `Dart`,
`regions` and the fixed gluing convention t <-> 1-t. The six oriented flags
per triangle and their barycentric identity map are constructed here. This
does not choose a new marking, draw a replacement route, or implement a band
neighborhood. Vendor limits (64 patches, 256 vertices per patch) are retained.
"""
function _ms_name!(used,prefix)
    n=1
    while "$prefix$n" in used;n+=1;end
    value="$prefix$n";push!(used,value);value
end

function _ms_regions(before,after,children)
    target=[_mm_cells(r) for r in after];result=Any[]
    for (i,region) in enumerate(before)
        images=Set(c for old in _mm_cells(region) for c in children[old])
        matches=findall(==(images),target)
        _mm_need(length(matches)==1,"subdivision does not transport an entire region exactly")
        push!(result,Dict("source_region"=>i,"target_region"=>only(matches),
            "source_triangles"=>[collect(c) for c in sort(collect(_mm_cells(region)))],
            "target_triangles"=>[collect(c) for c in sort(collect(images))]))
    end
    _mm_need(length(result)==length(target) && length(unique(r["target_region"] for r in result))==length(target),"subdivision changes region incidence")
    result
end

"""Subdivide each whole patch into six barycentric flags per old triangle.

Returns the actual native document and replayable identity-transport witness.
Old directed darts become two darts through their edge midpoint; old crossing
vertices stay fixed. Region marks and infinity map to the first child of their
saved triangle; their entire region correspondence is independently checked.
"""
function marked_subdivide(input::AbstractDict,specification::AbstractDict=Dict{String,Any}())
    _mm_need(isempty(specification),"barycentric subdivision takes no replacement or layout specification")
    before=_marked_states(_marked_document(input));st=before.geometry;M=MarkedSurfaceDiagrams
    doc=deepcopy(before.document);midpoints=Dict();children=Dict{Tuple{String,Int},Vector{Tuple{String,Int}}}()
    provenance=Any[];introduced=Any[];edges=Any[];triangles=Any[]
    for patch in doc["patches"]
        f=patch["face"];original=st.patches[f]
        used=Set(v for t in original["triangles"] for v in t)
        count_new=length(used)+length(st.incidence[f])+length(original["triangles"])
        _mm_need(count_new<=256,"barycentric subdivision requires $count_new vertices in patch $f; existing kernel limit is 256")
        edge_records=Any[];triangle_records=Any[]
        for e in sort(collect(keys(st.incidence[f])))
            v=_ms_name!(used,"subedge");midpoints[(f,e)]=v
            push!(edge_records,Dict("edge"=>collect(e),"midpoint"=>v))
            push!(introduced,Dict("target"=>[f,v],"source_edge"=>[f,e...],"barycentric_weights"=>["1/2","1/2"]))
            push!(edges,Dict("source"=>[f,e...],"target"=>[[f,e[1],v],[f,v,e[2]]]))
        end
        fine=Any[]
        for (i,t) in enumerate(original["triangles"])
            a,b,c=t;ab=midpoints[(f,M.edge(a,b))];bc=midpoints[(f,M.edge(b,c))];ca=midpoints[(f,M.edge(c,a))]
            z=_ms_name!(used,"subface")
            flags=[[a,ab,z],[ab,b,z],[b,bc,z],[bc,c,z],[c,ca,z],[ca,a,z]]
            first_index=length(fine);append!(fine,flags)
            children[(f,i-1)]=[(f,first_index+j) for j in 0:5]
            push!(triangle_records,Dict("corners"=>copy(t),"center"=>z))
            push!(introduced,Dict("target"=>[f,z],"source_triangle"=>[f,i-1],"source_vertices"=>copy(t),"barycentric_weights"=>["1/3","1/3","1/3"]))
            push!(triangles,Dict("source"=>[f,i-1],"target"=>[collect(x) for x in children[(f,i-1)]]))
        end
        old_parameters=deepcopy(get(original,"side_parameters",nothing))
        parameters=Dict{String,Any}()
        for side in patch["side_order"]
            path=original["sides"][side];ts=M.side_parameters(original,side)
            refined=String[];ps=String[]
            for j in 1:length(path)-1
                append!(refined,[path[j],midpoints[(f,M.edge(path[j],path[j+1]))]])
                append!(ps,[M.qstr(ts[j]),M.qstr((ts[j]+ts[j+1])/2)])
            end
            push!(refined,path[end]);push!(ps,M.qstr(ts[end]))
            patch["sides"][side]=refined;parameters[side]=ps
        end
        patch["side_parameters"]=parameters;patch["triangles"]=fine
        push!(provenance,Dict("face"=>f,"edges"=>edge_records,"triangles"=>triangle_records,"side_parameters"=>old_parameters))
    end
    for c in doc["curves"]
        darts=M.Dart[]
        for d in M.parse_word(c["word"])
            v=midpoints[(d.face,M.edge(d.u,d.v))]
            append!(darts,[M.Dart(d.face,d.u,v),M.Dart(d.face,v,d.v)])
        end
        c["word"]=M.print_word(darts)
    end
    for mark in get(doc,"region_marks",[])
        mark["seed"]=collect(first(children[_mm_cell(mark["seed"])]))
    end
    infinity=_mm_infinity(before.document)
    infinity!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=collect(first(children[_mm_cell(infinity)])))
    after=_marked_states(doc)
    correspondence=Dict{String,Any}("cell_semantics"=>"barycentric identity map: each new flag maps into its stated original triangle",
        "vertices"=>[Dict("source"=>[f,v],"target"=>[f,v]) for (f,v) in sort(collect(keys(st.vertex_dsu.p)))],
        "introduced_vertices"=>introduced,"edges"=>edges,"triangles"=>triangles,
        "regions"=>_ms_regions(M.regions(before.mathematics),M.regions(after.mathematics),children),
        "geometric_regions"=>_ms_regions(M.regions(before.geometry),M.regions(after.geometry),children),
        "curves"=>[Dict("source"=>c["id"],"target"=>c["id"],"source_word"=>c["word"],"target_word"=>only(x for x in doc["curves"] if x["id"]==c["id"])["word"]) for c in before.document["curves"]],
        "crossings"=>[Dict("source"=>x["id"],"target"=>x["id"],"at"=>deepcopy(x["at"])) for x in doc["crossings"]])
    inverse=Dict{String,Any}("patches"=>provenance,"coarse_document"=>_mm_hash(before.document),"refined_document"=>_mm_hash(doc),
        "region_seeds"=>_mm_inverse_seeds(before.document),"word_formats"=>Dict(c["id"]=>c["word"] for c in before.document["curves"]))
    infinity!==nothing && (inverse["infinity_triangle"]=deepcopy(infinity))
    witness=_mm_witness("barycentric_subdivision",before.document,doc,Dict{String,Any}(),inverse,correspondence)
    witness["inverse_kind"]="barycentric_collapse"
    (;document=doc,witness)
end

"""Collapse only a verified generated barycentric subdivision.

The saved local provenance names each original triangle's corners and its
edge/face subdivision vertices. Curves are actually collapsed in doubled-dart
pairs. The reconstructed candidate is then independently refined and compared
to the entire current document, proving applicability before returning it.
This is a generated-collapse operation, not arbitrary mesh simplification.
"""
function marked_collapse_subdivision(input::AbstractDict,specification::AbstractDict)
    before=_marked_states(_marked_document(input));doc=deepcopy(before.document);M=MarkedSurfaceDiagrams
    spec=_plain_value(specification)
    _mm_need(_mm_hash(doc)==spec["refined_document"],"subdivision collapse source has changed")
    _mm_need(Set(p["face"] for p in spec["patches"])==Set(p["face"] for p in doc["patches"]),"collapse provenance patch set differs")
    mids=Dict();parent=Dict{Tuple{String,Int},Tuple{String,Int}}()
    for p in spec["patches"]
        f=p["face"];patch=only(x for x in doc["patches"] if x["face"]==f)
        originals=Set(v for t in p["triangles"] for v in t["corners"])
        for e in p["edges"]
            _mm_need(!haskey(mids,(f,e["midpoint"])),"duplicate collapse midpoint")
            mids[(f,e["midpoint"])]=Set(String.(e["edge"]))
        end
        patch["triangles"]=[copy(t["corners"]) for t in p["triangles"]]
        for j in eachindex(p["triangles"]),k in 0:5;parent[(f,6(j-1)+k)]=(f,j-1);end
        for side in patch["side_order"]
            patch["sides"][side]=filter(v->v in originals,patch["sides"][side])
        end
        if p["side_parameters"]===nothing
            delete!(patch,"side_parameters")
        else
            patch["side_parameters"]=deepcopy(p["side_parameters"])
        end
    end
    for c in doc["curves"]
        ds=M.parse_word(c["word"]);out=M.Dart[]
        _mm_need(iseven(length(ds)),"collapse curve does not consist of doubled darts")
        for j in 1:2:length(ds)
            a,b=ds[j:j+1]
            _mm_need(a.face==b.face && a.v==b.u && get(mids,(a.face,a.v),nothing)==Set([a.u,b.v]),"curve does not follow an original edge through its generated midpoint")
            push!(out,M.Dart(a.face,a.u,b.v))
        end
        word=spec["word_formats"][c["id"]]
        _mm_need(M.parse_word(word)==out,"saved lexical format differs from the constructed coarse curve")
        c["word"]=word
    end
    for mark in get(doc,"region_marks",[])
        expected=parent[_mm_cell(mark["seed"])];preferred=_mm_cell(spec["region_seeds"][mark["id"]])
        _mm_need(expected==preferred,"collapse selected-region seed is not a child of its original triangle")
        mark["seed"]=collect(preferred)
    end
    infinity=_mm_infinity(before.document)
    if infinity!==nothing
        target=_mm_cell(spec["infinity_triangle"])
        _mm_need(parent[_mm_cell(infinity)]==target,"collapse infinity is not a child of its original triangle")
        doc["metadata"]["display"]["infinity_triangle"]=collect(target)
    end
    _marked_states(doc)
    _mm_need(_mm_hash(doc)==spec["coarse_document"],"constructed coarse document differs from recorded identity")
    generated=marked_subdivide(doc)
    _mm_need(generated.document==before.document && _canonical_json(generated.witness["inverse_specification"])==_canonical_json(spec),"current mesh/provenance is not the recorded barycentric subdivision")
    witness=_mm_witness("barycentric_collapse",before.document,doc,spec,Dict{String,Any}(),
        Dict("cell_semantics"=>"inverse of the verified barycentric identity transport", "forward_subdivision"=>generated.witness["correspondence"]))
    witness["inverse_kind"]="barycentric_subdivision"
    (;document=doc,witness)
end

"""Repatch a closed quotient mesh by oriented disks around its vertices.

Source: vendor `fine_edge_occurrences`, `vertex_dsu`, triangle corner links,
and `reference_surface.jl::lickorish_reference`'s crossing-disk/gate pattern.
Every source triangle contributes its six barycentric flags to three stars.
Crossings become star centers; each original curve edge crosses a midpoint
gate. This is an explicit marked identity transport, not a new inferred marking.
Requires 3..127 corners per star. Default requires at most 64 stars. Explicit
`merge_charts=true` reuses the existing verified disk merger on the bare atlas;
optional `subdivision_budget` conservatively bounds its next refinement.
True-boundary half-disks remain outside this operation's current scope.
"""
function marked_dual_star(input::AbstractDict,specification::AbstractDict=Dict{String,Any}())
    spec=_plain_value(specification)
    _mm_need(issubset(Set(keys(spec)),Set(["merge_charts","subdivision_budget"])),"unknown dual-star option; replacement data is not accepted")
    merge_charts=get(spec,"merge_charts",false);budget=get(spec,"subdivision_budget",nothing)
    _mm_need(merge_charts isa Bool,"merge_charts must be boolean")
    _mm_need(budget===nothing || (budget isa Integer && !(budget isa Bool) && 3<=budget<=256),"subdivision_budget must lie in 3:256")
    _mm_need(budget===nothing || merge_charts,"subdivision_budget requires explicit merge_charts")
    before=_marked_states(_marked_document(input));st=before.geometry;M=MarkedSurfaceDiagrams
    _mm_need(st.surface_report["boundary_components"]==0,"dual-star repatching currently requires a closed surface; boundary half-disks are not implemented")
    pairs=M.fine_edge_occurrences(st);groups=Dict()
    for ((f,i),t) in st.triangles,j in 1:3
        rt=M.root!(st.vertex_dsu,(f,t[j]))
        push!(get!(groups,rt,Tuple{String,Int,Int}[]),(f,i,j))
    end
    _mm_need(length(groups)<=64 || merge_charts,"dual-star atlas needs $(length(groups)) disk patches; explicit merge_charts is required above the existing limit 64")
    star=Dict(rt=>"Star$k" for (k,rt) in enumerate(sort(collect(keys(groups)))))
    outgoing=Dict();corners=Dict();patches=Any[];gluings=Any[];flag_map=Dict{Tuple{String,Int},Vector{Tuple{String,Int}}}()
    introduced=Any[];source_vertices=Any[]
    for rt in sort(collect(keys(groups)))
        first_corner=minimum(groups[rt]);cycle=Tuple{String,Int,Int}[];c=first_corner
        while isempty(cycle) || c!=first_corner
            _mm_need(!(c in cycle),"quotient vertex link does not form one disk-star cycle")
            push!(cycle,c);f,i,j=c;t=st.triangles[(f,i)]
            incoming=(f,i,t[mod1(j-1,3)],t[j])
            _mm_need(haskey(pairs,incoming),"closed vertex-star link has an unpaired edge")
            q=pairs[incoming];next_triangle=st.triangles[q[1:2]]
            c=(q[1],q[2],only(findall(==(q[3]),next_triangle)))
        end
        _mm_need(Set(cycle)==Set(groups[rt]),"quotient vertex link has multiple star cycles")
        n=length(cycle)
        _mm_need(n>=3,"dual-star disk has fewer than three sides; refine the source before repatching")
        _mm_need(2n+1<=256,"dual-star disk requires $(2n+1) vertices; existing patch limit is 256")
        name=star[rt];sides=Dict();order=String[];triangles=Any[]
        for (k,(f,i,j)) in enumerate(cycle)
            t=st.triangles[(f,i)];side="e$k";gate="p$k";a="t$(mod1(k-1,n))";b="t$k"
            sides[side]=[a,gate,b];push!(order,side)
            outgoing[(f,i,t[j],t[mod1(j+1,3)])]=(name,side,gate)
            corners[(f,i,j)]=(name,k)
            append!(triangles,[["z",gate,b],["z",b,"p$(mod1(k+1,n))"]])
            append!(get!(flag_map,(f,i-1),Tuple{String,Int}[]),[(name,2(k-1)),(name,2(k-1)+1)])
            push!(introduced,Dict("target"=>[name,gate],"source_edge_occurrence"=>[f,i-1,t[j],t[mod1(j+1,3)]],"barycentric_weights"=>["1/2","1/2"]))
            push!(introduced,Dict("target"=>[name,b],"source_triangle"=>[f,i-1],"barycentric_weights"=>["1/3","1/3","1/3"]))
        end
        push!(patches,Dict("face"=>name,"side_order"=>order,"sides"=>sides,"triangles"=>triangles))
        push!(source_vertices,Dict("source_quotient_vertex"=>collect(rt),"target"=>[name,"z"]))
    end
    for (o,p) in sort(collect(pairs);by=first)
        isless(o,p) || continue
        a,b=outgoing[o],outgoing[p];push!(gluings,[collect(a[1:2]),collect(b[1:2])])
    end
    # An old oriented fine edge is represented by two chart-interior darts.
    function image_dart(d)
        i,u,v=first(st.incidence[d.face][M.edge(d.u,d.v)])
        o=(d.face,i,u,v);p=pairs[o]
        a,b=(d.u==u ? (outgoing[o],outgoing[p]) : (outgoing[p],outgoing[o]))
        [M.Dart(a[1],"z",a[3]),M.Dart(b[1],b[3],"z")]
    end
    route_provenance=Any[];edge_map=Any[]
    for f in sort(collect(keys(st.patches))),e in sort(collect(keys(st.incidence[f])))
        # Routes on old marking sides are prohibited, so inverse curve recovery
        # only needs interior source edges; seam edges still have cell images.
        a,b=e;image=image_dart(M.Dart(f,a,b))
        push!(edge_map,Dict("source"=>[f,a,b],"target"=>[[d.face,d.u,d.v] for d in image]))
        if length(st.incidence[f][e])==2
            for (u,v) in ((a,b),(b,a))
                ds=image_dart(M.Dart(f,u,v))
                push!(route_provenance,Dict("source"=>[f,u,v],"target"=>[[d.face,d.u,d.v] for d in ds]))
            end
        end
    end
    curves=deepcopy(before.document["curves"])
    for c in curves
        ds=M.Dart[]
        for d in M.parse_word(c["word"]);append!(ds,image_dart(d));end
        c["word"]=M.print_word(ds)
    end
    crossings=deepcopy(before.document["crossings"])
    for x in crossings
        f,v=x["at"];rt=M.root!(st.vertex_dsu,(f,v));x["at"]=[star[rt],"z"]
    end
    doc=Dict{String,Any}("version"=>1,"surface"=>Dict("marking_id"=>before.document["surface"]["marking_id"]*"/dual-star",
        "genus"=>st.surface_report["genus"],"gluings"=>gluings,"boundary"=>[]),"patches"=>patches,"curves"=>curves,"crossings"=>crossings)
    atlas_map=Dict{String,Any}();merges=Any[]
    if merge_charts
        # Source: ograph_reverse.jl::_reverse_merge_charts!, independently
        # exercised on >64 bare charts. Reuse its exact prefix, side gluing and
        # triangle-offset maps; do not implement another chart merger here.
        refinement_limit=budget===nothing ? nothing : Int(budget)
        merges=_reverse_merge_charts!(doc;maximum_refined_vertices=refinement_limit,correspondence=atlas_map)
        vm=Dict(Tuple(e["source"])=>Tuple(e["target"]) for e in atlas_map["vertices"])
        tm=Dict(_mm_cell(e["source"])=>_mm_cell(e["target"]) for e in atlas_map["triangles"])
        function mapped_dart(x)
            a=vm[(x[1],x[2])];b=vm[(x[1],x[3])]
            _mm_need(a[1]==b[1],"chart merger split one local edge across two charts")
            [a[1],a[2],b[2]]
        end
        for c in keys(flag_map);flag_map[c]=[tm[x] for x in flag_map[c]];end
        for row in source_vertices;row["target"]=collect(vm[Tuple(row["target"])]);end
        for row in introduced;row["target"]=collect(vm[Tuple(row["target"])]);end
        for row in edge_map;row["target"]=[mapped_dart(x) for x in row["target"]];end
        for row in route_provenance;row["target"]=[mapped_dart(x) for x in row["target"]];end
    end
    haskey(before.document,"metadata") && (doc["metadata"]=deepcopy(before.document["metadata"]))
    if haskey(before.document,"region_marks")
        doc["region_marks"]=deepcopy(before.document["region_marks"])
        for mark in doc["region_marks"];mark["seed"]=collect(first(flag_map[_mm_cell(mark["seed"])]));end
    end
    infinity=_mm_infinity(before.document)
    infinity!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=collect(first(flag_map[_mm_cell(infinity)])))
    after=_marked_states(doc)
    correspondence=Dict{String,Any}("cell_semantics"=>"explicit quotient-vertex star atlas; six barycentric flags map into each source triangle",
        "vertices"=>source_vertices,"introduced_vertices"=>introduced,"edges"=>edge_map,"chart_merges"=>merges,"chart_correspondence"=>atlas_map,
        "triangles"=>[Dict("source"=>collect(c),"target"=>[collect(x) for x in flag_map[c]]) for c in sort(collect(keys(flag_map)))],
        "regions"=>_ms_regions(M.regions(before.mathematics),M.regions(after.mathematics),flag_map),
        "geometric_regions"=>_ms_regions(M.regions(before.geometry),M.regions(after.geometry),flag_map),
        "curves"=>[Dict("source"=>c["id"],"target"=>c["id"],"source_word"=>c["word"],"target_word"=>only(x for x in curves if x["id"]==c["id"])["word"]) for c in before.document["curves"]],
        "crossings"=>[Dict("source"=>x["id"],"target"=>x["id"],"source_at"=>x["at"],"target_at"=>only(y for y in crossings if y["id"]==x["id"])["at"]) for x in before.document["crossings"]])
    inverse=Dict{String,Any}("source_mesh"=>Dict("surface"=>deepcopy(before.document["surface"]),"patches"=>deepcopy(before.document["patches"])),
        "route_edges"=>route_provenance,"crossing_locations"=>Dict(x["id"]=>deepcopy(x["at"]) for x in before.document["crossings"]),
        "region_seeds"=>_mm_inverse_seeds(before.document),"word_formats"=>Dict(c["id"]=>c["word"] for c in before.document["curves"]),
        "source_document"=>_mm_hash(before.document),"star_document"=>_mm_hash(doc),"forward_specification"=>spec)
    infinity!==nothing && (inverse["infinity_triangle"]=deepcopy(infinity))
    witness=_mm_witness("dual_star",before.document,doc,spec,inverse,correspondence)
    witness["inverse_kind"]="dual_star_collapse"
    (;document=doc,witness)
end

"""Undo a verified generated dual-star atlas by collapsing paired route darts.

Original mesh incidence is local geometric provenance; curve records are
reconstructed from actual paired darts, not a saved source diagram. A fresh
forward repatching must reproduce the entire current document and provenance.
"""
function marked_collapse_dual_star(input::AbstractDict,specification::AbstractDict)
    before=_marked_states(_marked_document(input));spec=_plain_value(specification);M=MarkedSurfaceDiagrams
    _mm_need(_mm_hash(before.document)==spec["star_document"],"dual-star collapse source has changed")
    lookup=Dict{String,Any}()
    for r in spec["route_edges"]
        key=_canonical_json(r["target"])
        _mm_need(!haskey(lookup,key),"ambiguous original route edge in dual-star collapse")
        lookup[key]=r["source"]
    end
    doc=deepcopy(before.document);doc["surface"]=deepcopy(spec["source_mesh"]["surface"]);doc["patches"]=deepcopy(spec["source_mesh"]["patches"])
    for c in doc["curves"]
        ds=M.parse_word(c["word"]);out=M.Dart[]
        _mm_need(iseven(length(ds)),"dual-star route is not a sequence of source edge images")
        for j in 1:2:length(ds)
            key=_canonical_json([[d.face,d.u,d.v] for d in ds[j:j+1]])
            _mm_need(haskey(lookup,key),"dual-star curve pair is not an original source edge")
            f,u,v=lookup[key];push!(out,M.Dart(f,u,v))
        end
        word=spec["word_formats"][c["id"]]
        _mm_need(M.parse_word(word)==out,"saved lexical format changes the reconstructed original route")
        c["word"]=word
    end
    for x in doc["crossings"];x["at"]=deepcopy(spec["crossing_locations"][x["id"]]);end
    for mark in get(doc,"region_marks",[]);mark["seed"]=deepcopy(spec["region_seeds"][mark["id"]]);end
    _mm_infinity(before.document)!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=deepcopy(spec["infinity_triangle"]))
    _marked_states(doc)
    _mm_need(_mm_hash(doc)==spec["source_document"],"constructed source chart document differs from its identity")
    generated=marked_dual_star(doc,spec["forward_specification"])
    _mm_need(generated.document==before.document && _canonical_json(generated.witness["inverse_specification"])==_canonical_json(spec),"current chart document/provenance is not the recorded dual-star transport")
    witness=_mm_witness("dual_star_collapse",before.document,doc,spec,spec["forward_specification"],
        Dict("cell_semantics"=>"inverse of verified dual-star identity atlas", "forward_dual_star"=>generated.witness["correspondence"]))
    witness["inverse_kind"]="dual_star"
    (;document=doc,witness)
end

# Source: MarkedSurfaceDiagrams._validate's oriented link crossing test.
# Reused for newly constructed intersections, with the existing curve first
# and the new parallel curve second; its ordered pair determines the sign.
function _ms_add_crossings!(doc,new_id)
    M=MarkedSurfaceDiagrams;bare=deepcopy(doc)
    filter!(c->c["id"]!=new_id,bare["curves"])
    st=_marked_states(bare).geometry
    all_curves=Dict(c["id"]=>M.parse_word(c["word"]) for c in doc["curves"])
    function passes(ds)
        Dict(M.root!(st.vertex_dsu,(b.face,b.u))=>(a,b) for (a,b) in zip(circshift(ds,1),ds))
    end
    new_pass=passes(all_curves[new_id]);ids=Set(x["id"] for x in doc["crossings"]);added=String[]
    for c in doc["curves"]
        cid=c["id"];cid==new_id && continue;old_pass=passes(all_curves[cid])
        for rt in sort(collect(intersect(Set(keys(new_pass)),Set(keys(old_pass)))))
            f,v=rt
            _mm_need(!(v in st.boundaries[f]),"constructed crossing lies on a marking seam; local refinement is required")
            a0,a1=old_pass[rt];b0,b1=new_pass[rt]
            _mm_need(all(d->d.face==f,(a0,a1,b0,b1)),"constructed crossing germs use inconsistent charts")
            links=Dict()
            for t in st.patches[f]["triangles"]
                j=findfirst(==(v),t);j===nothing && continue
                links[t[mod1(j+1,3)]]=t[mod1(j+2,3)]
            end
            ao,ai,bo,bi=a1.v,a0.u,b1.v,b0.u;order=String[];w=ao
            while isempty(order) || w!=ao
                _mm_need(!(w in order),"constructed crossing link is not cyclic")
                push!(order,w);w=links[w]
            end
            reduced=filter(x->x in (ao,ai,bo,bi),order)
            _mm_need(reduced in ([ao,bo,ai,bi],[ao,bi,ai,bo]),"constructed parallel curve is not transverse")
            n=1;id="parallel_$n"
            while id in ids;n+=1;id="parallel_$n";end
            push!(ids,id);push!(added,id)
            push!(doc["crossings"],Dict("id"=>id,"at"=>[f,v],"curves"=>[cid,new_id],"sign"=>(reduced==[ao,bo,ai,bi] ? 1 : -1)))
        end
    end
    added
end

"""Construct a distinct oriented auxiliary parallel copy through an annulus.

The curve and requested left/right side are operation parameters; no output
route is accepted. Dual-star repatching makes the curve a full subcomplex;
one further derived subdivision constructs its closed regular star.
The selected star must independently normalize to an annulus. Its frontier
cycle and all new crossing signs are computed from the same marked mesh.
Resource/atlas restrictions of the preceding constructors remain explicit.
"""
function marked_parallel_copy(input::AbstractDict,specification::AbstractDict)
    original=_marked_document(input);spec=_plain_value(specification);M=MarkedSurfaceDiagrams
    cid=String(spec["curve"]);new_id=String(spec["id"]);side=String(get(spec,"side","left"))
    _mm_need(side in ("left","right"),"parallel side must be left or right")
    _mm_need(any(c->c["id"]==cid,original["curves"]),"unknown parallel source curve")
    _mm_need(!any(c->c["id"]==new_id,original["curves"]),"parallel curve id already exists")
    _mm_need(!isempty(new_id),"parallel curve id must be nonempty")
    # Agent proof from marked_dual_star's actual flags: mapped edge graphs use
    # z/p vertices, never t; a p gate is in the graph iff both corresponding
    # z-p halfedges are owned. Gluings preserve these types. Thus the graph is
    # already full; retain the independent fullness checks below.
    star=marked_dual_star(original,Dict("merge_charts"=>true,"subdivision_budget"=>256))
    coarse=_marked_states(star.document).geometry
    _mm_need(coarse.curves[cid][1]["closed"],"parallel source must be closed")
    core=coarse.curves[cid][2]
    vertices=Set(M.root!(coarse.vertex_dsu,(d.face,d.u)) for d in core)
    owned=Set((d.face,M.edge(d.u,d.v)) for d in core)
    for (f,inc) in coarse.incidence,(e,_) in inc
        a,b=e
        if M.root!(coarse.vertex_dsu,(f,a)) in vertices && M.root!(coarse.vertex_dsu,(f,b)) in vertices
            _mm_need((f,e) in owned,"chosen curve is not a full subcomplex: unused chord edge")
        end
    end
    for ((f,_),t) in coarse.triangles
        _mm_need(!all(v->M.root!(coarse.vertex_dsu,(f,v)) in vertices,t),"chosen curve subcomplex contains a whole triangle")
    end
    second_refinement=marked_subdivide(star.document);doc=deepcopy(second_refinement.document)
    fine=_marked_states(doc).geometry;selected=Set{Tuple{String,Int}}()
    for ((f,i),t) in fine.triangles
        old=[v for v in t if haskey(coarse.vertex_dsu.p,(f,v))]
        _mm_need(length(old)==1,"derived flag does not have one original vertex")
        M.root!(coarse.vertex_dsu,(f,only(old))) in vertices && push!(selected,(f,i-1))
    end
    pairing=M.fine_edge_occurrences(fine);frontier=Set()
    for (a,b) in pairing
        (a[1],a[2]-1) in selected && !((b[1],b[2]-1) in selected) && push!(frontier,(a[1],M.edge(a[3],a[4])))
    end
    regions=M.regions(fine;barriers=frontier)
    candidates=[r for r in regions if _mm_cells(r)==selected]
    _mm_need(length(candidates)==1,"derived regular star is not one normalized component")
    annulus=only(candidates)
    _mm_need(annulus["genus"]==0 && annulus["boundary_components"]==2,"derived regular star is not an annulus")
    cycles=annulus["boundary_cycles"]
    d=first(core);left=only(o for o in coarse.incidence[d.face][M.edge(d.u,d.v)] if o[2]==d.u && o[3]==d.v)
    # Source: vendor `regions` exports boundary-cycle triangle indices as
    # zero-based, whereas incidence entries (`left[1]`) remain one-based.
    left_cycle=findall(c->any(o->o[1]==d.face && div(o[2],6)==left[1]-1,c),cycles)
    _mm_need(length(left_cycle)==1,"annulus frontier does not determine one oriented left side")
    chosen=side=="left" ? only(left_cycle) : only(setdiff(1:2,left_cycle))
    ds=[M.Dart(f,u,v) for (f,_,u,v) in cycles[chosen]]
    side=="left" && (ds=[M.Dart(d.face,d.v,d.u) for d in reverse(ds)])
    push!(doc["curves"],Dict("id"=>new_id,"role"=>"aux","closed"=>true,"word"=>M.print_word(ds)))
    added=_ms_add_crossings!(doc,new_id)
    _marked_states(doc)
    history=[star.witness,second_refinement.witness]
    inverse=Dict("curve"=>cid,"side"=>side,"new_curve"=>new_id,"new_crossings"=>added,"construction"=>history,
        "refined_document"=>_mm_hash(second_refinement.document),"output_document"=>_mm_hash(doc),"source_document"=>_mm_hash(original))
    correspondence=Dict("construction"=>history,"core_curve"=>cid,"new_curve"=>new_id,"side"=>side,
        "annulus_triangles"=>[collect(c) for c in sort(collect(selected))],"annulus_boundary_words"=>[M.print_word([M.Dart(f,u,v) for (f,_,u,v) in c]) for c in cycles],
        "annulus_boundary_cycles"=>cycles,"selected_boundary"=>chosen,"source_left_triangle"=>[d.face,left[1]-1],
        "new_crossings"=>added,"infinity_semantics"=>"retained transported source point; added auxiliary curve may split its geometric region")
    witness=_mm_witness("parallel_copy",original,doc,spec,inverse,correspondence);witness["inverse_kind"]="parallel_remove"
    (;document=doc,witness)
end

"""Remove a generated auxiliary parallel, then undo its verified mesh history."""
function marked_remove_parallel(input::AbstractDict,specification::AbstractDict)
    doc=_marked_document(input);spec=_plain_value(specification)
    _mm_need(_mm_hash(doc)==spec["output_document"],"generated parallel output has changed")
    before=deepcopy(doc);cid=spec["new_curve"]
    _mm_need(count(c->c["id"]==cid && c["role"]=="aux",doc["curves"])==1,"generated auxiliary curve is missing")
    filter!(c->c["id"]!=cid,doc["curves"])
    removed=Set(spec["new_crossings"])
    filter!(x->!(x["id"] in removed),doc["crossings"])
    _mm_need(_mm_hash(doc)==spec["refined_document"],"removing generated data does not restore the recorded refined object")
    for w in reverse(spec["construction"])
        doc=replay_marked_move(doc,inverse_marked_move(w)).document
    end
    _mm_need(_mm_hash(doc)==spec["source_document"],"inverse mesh history does not restore the parallel source")
    # Forward parameters are recovered from the actual construction's core and
    # annulus record by the top-level replay witness, not a retained diagram.
    inverse=Dict("curve"=>spec["curve"],"id"=>cid,"side"=>spec["side"])
    generated=marked_parallel_copy(doc,inverse)
    _mm_need(generated.document==before && _canonical_json(generated.witness["inverse_specification"])==_canonical_json(spec),"parallel provenance/forward parameters do not reconstruct the current output")
    witness=_mm_witness("parallel_remove",before,doc,spec,inverse,Dict("inverse_construction"=>spec["construction"]))
    witness["inverse_kind"]="parallel_copy"
    (;document=doc,witness)
end
