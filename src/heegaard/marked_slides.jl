"""Constructive marked handle slides through an explicit band neighborhood.

Sources: marked_subdivision.jl::marked_dual_star/marked_subdivide and its
verified derived regular-neighborhood construction; vendor
MarkedSurfaceDiagrams._validate/regions for word, link and cut conventions.
The new curve is computed as the pair-of-pants frontier carrying both cores.
No replacement curve or post-slide region partition is an operation input.
"""
function _mh_transport_darts(ds,history)
    M=MarkedSurfaceDiagrams;out=copy(ds)
    for w in history
        edges=Dict((e["source"][1],M.edge(e["source"][2],e["source"][3]))=>e for e in w["correspondence"]["edges"])
        next=M.Dart[]
        for d in out
            e=edges[(d.face,M.edge(d.u,d.v))]
            image=[M.Dart(x[1],x[2],x[3]) for x in e["target"]]
            d.u!=e["source"][2] && (image=[M.Dart(x.face,x.v,x.u) for x in reverse(image)])
            append!(next,image)
        end
        out=next
    end
    out
end

# Source: vendor crossing-link test in `_validate`. Here the second path is
# an embedded band core with endpoints on distinct same-family curves rather
# than a saved proper arc on physical boundary. Its interior must be transverse.
function _mh_transverse(st,a0,a1,b0,b1)
    M=MarkedSurfaceDiagrams;f=a1.face;v=a1.u
    _mm_need(all(d->d.face==f,(a0,a1,b0,b1)),"band crossing germs must lie in one interior chart")
    _mm_need(!(v in st.boundaries[f]),"band crosses another curve on a marking seam")
    links=Dict()
    for t in st.patches[f]["triangles"]
        j=findfirst(==(v),t);j===nothing && continue
        links[t[mod1(j+1,3)]]=t[mod1(j+2,3)]
    end
    ao,ai,bo,bi=a1.v,a0.u,b1.v,b0.u;order=String[];w=ao
    while isempty(order) || w!=ao
        _mm_need(!(w in order),"band crossing link is not cyclic")
        push!(order,w);w=links[w]
    end
    reduced=filter(x->x in (ao,ai,bo,bi),order)
    _mm_need(reduced in ([ao,bo,ai,bi],[ao,bi,ai,bo]),"band is tangent rather than transverse to an existing curve")
    true
end

function _mh_band(st,source,target,word)
    M=MarkedSurfaceDiagrams
    _mm_need(source!=target && haskey(st.curves,source) && haskey(st.curves,target),"slide needs two distinct existing curves")
    s,sd=st.curves[source];t,td=st.curves[target]
    _mm_need(s["role"]==t["role"] && s["role"] in ("alpha","beta"),"slide curves must be in the same alpha/beta family")
    _mm_need(s["closed"] && t["closed"],"slide curves must be closed")
    band=M.parse_word(word)
    global_vertex(f,v)=M.root!(st.vertex_dsu,(f,v))
    for d in band
        _mm_need(haskey(st.incidence,d.face) && haskey(st.incidence[d.face],M.edge(d.u,d.v)),"band references a missing mesh edge")
        _mm_need(length(st.incidence[d.face][M.edge(d.u,d.v)])==2,"band may cross a marking seam but cannot run along it")
        _mm_need(!haskey(st.used,(d.face,M.edge(d.u,d.v))),"band overlaps an existing curve edge")
        corners=Set(st.patches[d.face]["sides"][x][1] for x in st.patches[d.face]["side_order"])
        _mm_need(!(d.u in corners) && !(d.v in corners),"band passes through a marked patch corner")
    end
    for (a,b) in zip(band[1:end-1],band[2:end])
        _mm_need(global_vertex(a.face,a.v)==global_vertex(b.face,b.u),"band edge word is discontinuous")
        _mm_need((a.face,a.v)!=(b.face,b.u) || !(a.v in st.boundaries[a.face]),"band turns back at a marking gate")
    end
    vertices=[global_vertex(d.face,d.u) for d in band];push!(vertices,global_vertex(band[end].face,band[end].v))
    _mm_need(length(vertices)==length(Set(vertices)),"band is not embedded")
    source_vertices=Set(global_vertex(d.face,d.u) for d in sd)
    target_vertices=Set(global_vertex(d.face,d.u) for d in td)
    _mm_need(first(vertices) in source_vertices && last(vertices) in target_vertices,"oriented band must attach from source to target")
    visits=Dict()
    for (cid,(c,ds)) in st.curves
        for (a,b) in zip(circshift(ds,1),ds)
            push!(get!(visits,global_vertex(b.face,b.u),Any[]),(cid,a,b))
        end
    end
    _mm_need(length(get(visits,first(vertices),[]))==1 && length(get(visits,last(vertices),[]))==1,"band attachment points must avoid existing crossings")
    for j in 2:length(vertices)-1
        at=get(visits,vertices[j],[])
        _mm_need(length(at)<=1,"band would create a triple point")
        for (cid,a,b) in at
            _mm_need(st.curves[cid][1]["role"]!=s["role"],"band interior meets its own curve family")
            _mh_transverse(st,a,b,band[j-1],band[j])
        end
    end
    band
end

function _mh_neighborhood(coarse,fine,paths)
    M=MarkedSurfaceDiagrams
    vertices=Set(M.root!(coarse.vertex_dsu,(d.face,d.u)) for ds in values(paths) for d in ds)
    # Include the terminal band endpoint as well as all initial edge vertices.
    for ds in values(paths),d in ds;push!(vertices,M.root!(coarse.vertex_dsu,(d.face,d.v)));end
    owned=Dict{Any,String}()
    for (label,ds) in paths,d in ds
        key=(d.face,M.edge(d.u,d.v))
        _mm_need(!haskey(owned,key),"neighborhood graph contains overlapping labelled edges")
        owned[key]=label
    end
    for (f,inc) in coarse.incidence,(e,_) in inc
        if all(v->M.root!(coarse.vertex_dsu,(f,v)) in vertices,e)
            _mm_need(haskey(owned,(f,e)),"slide graph is not a full subcomplex: unused chord edge")
        end
    end
    labels=Dict{Tuple{String,Int},String}()
    for ((f,i),t) in coarse.triangles
        _mm_need(!all(v->M.root!(coarse.vertex_dsu,(f,v)) in vertices,t),"slide graph contains a whole triangle")
        kinds=[owned[(f,M.edge(a,b))] for (a,b) in M.tedges(t) if haskey(owned,(f,M.edge(a,b)))]
        _mm_need(length(kinds)<=1,"full graph has multiple edges of a coarse triangle")
        !isempty(kinds) && (labels[(f,i-1)]=only(kinds))
    end
    selected=Set{Tuple{String,Int}}()
    for ((f,i),t) in fine.triangles
        old=[v for v in t if haskey(coarse.vertex_dsu.p,(f,v))]
        _mm_need(length(old)==1,"derived slide flag must have one coarse vertex")
        M.root!(coarse.vertex_dsu,(f,only(old))) in vertices && push!(selected,(f,i-1))
    end
    frontier=Set()
    for (a,b) in M.fine_edge_occurrences(fine)
        (a[1],a[2]-1) in selected && !((b[1],b[2]-1) in selected) && push!(frontier,(a[1],M.edge(a[3],a[4])))
    end
    candidates=[r for r in M.regions(fine;barriers=frontier) if _mm_cells(r)==selected]
    _mm_need(length(candidates)==1,"slide graph regular neighborhood is not one normalized component")
    pants=only(candidates)
    _mm_need(pants["genus"]==0 && pants["boundary_components"]==3,"slide graph regular neighborhood is not a pair of pants")
    cycles=pants["boundary_cycles"]
    tags=[Set(get(labels,(o[1],div(o[2],6)),"") for o in cycle) for cycle in cycles]
    chosen=findall(t->"source" in t && "target" in t,tags)
    _mm_need(length(chosen)==1,"pair-of-pants frontier does not select one source-target band sum")
    index=only(chosen);cycle=cycles[index];d=first(paths["source"])
    left=only(o for o in coarse.incidence[d.face][M.edge(d.u,d.v)] if o[2]==d.u && o[3]==d.v)
    right=only(o for o in coarse.incidence[d.face][M.edge(d.u,d.v)] if o[2]==d.v && o[3]==d.u)
    on_left=any(o->o[1]==d.face && div(o[2],6)==left[1]-1,cycle)
    on_right=any(o->o[1]==d.face && div(o[2],6)==right[1]-1,cycle)
    _mm_need(xor(on_left,on_right),"band-sum orientation is not determined at the source core")
    ds=[M.Dart(f,u,v) for (f,_,u,v) in cycle]
    on_left && (ds=[M.Dart(d.face,d.v,d.u) for d in reverse(ds)])
    (;pants,selected,cycles,selected_boundary=index,boundary_tags=[sort(collect(t)) for t in tags],darts=ds)
end

function _mh_disk_policy!(doc,before,policy)
    M=MarkedSurfaceDiagrams
    _mm_need(policy in ("preserve","follow_seed"),"disk_policy must be preserve or follow_seed")
    old=M.regions(before.mathematics);next=_marked_states(_mm_bare(doc));new=M.regions(next.mathematics);records=Any[]
    for mark in get(doc,"region_marks",[])
        seed=_mm_cell(mark["seed"])
        a=only(r for r in old if seed in _mm_cells(r));b=only(r for r in new if seed in _mm_cells(r))
        unchanged=_mm_cells(a)==_mm_cells(b)
        _mm_need(policy=="follow_seed" || unchanged,"slide changes marked region $(mark["id"]); default preserve rejects, explicit follow_seed is required")
        push!(records,Dict("id"=>mark["id"],"policy"=>policy,"region_preserved"=>unchanged,"claim"=>mark["claim"],
            "protected_seed"=>copy(mark["seed"]),"old_triangles"=>a["triangles"],"new_triangles"=>b["triangles"],
            "old_boundary_cycles"=>a["boundary_cycles"],"new_boundary_cycles"=>b["boundary_cycles"]))
    end
    _marked_states(doc) # Revalidate original D disk/joining claims after policy.
    records
end

"""Construct an oriented handle slide from two curves and an embedded band.

Specification: `source`, `target`, `band` (explicit original edge word), and
optional `disk_policy`. Default preserve demands exact complete transported
D cell membership. Explicit follow_seed continues its protected point and
rechecks the same disk/joining claim, recording the changed region boundary.
"""
function marked_handle_slide(input::AbstractDict,specification::AbstractDict)
    _md_guard_periodic_operation(input)
    original=_marked_document(input);spec=_plain_value(specification);M=MarkedSurfaceDiagrams
    source=String(spec["source"]);target=String(spec["target"])
    initial=_marked_states(original);band=_mh_band(initial.geometry,source,target,spec["band"])
    policy=String(get(spec,"disk_policy","preserve"))
    _mm_need(policy in ("preserve","follow_seed"),"disk_policy must be preserve or follow_seed")
    # Source: marked_dual_star flags [z,p,t]. Any transported source edge graph
    # (including the band) uses z/p, never t, and is already a full subcomplex.
    # Keep `_mh_neighborhood`'s actual fullness and one-pants checks unchanged.
    star=marked_dual_star(original,Dict("merge_charts"=>true,"subdivision_budget"=>256))
    first_history=[star.witness]
    coarse=_marked_states(star.document).geometry
    paths=Dict("source"=>coarse.curves[source][2],"target"=>coarse.curves[target][2],"band"=>_mh_transport_darts(band,first_history))
    second_refinement=marked_subdivide(star.document);before=_marked_states(second_refinement.document)
    neighborhood=_mh_neighborhood(coarse,before.geometry,paths)
    doc=deepcopy(before.document);old_source=only(c for c in doc["curves"] if c["id"]==source)
    old_word=old_source["word"];old_source["word"]=M.print_word(neighborhood.darts)
    removed=deepcopy([x for x in doc["crossings"] if source in x["curves"]])
    removed_positions=[Dict("index"=>i,"crossing"=>deepcopy(x)) for (i,x) in enumerate(doc["crossings"]) if source in x["curves"]]
    filter!(x->!(source in x["curves"]),doc["crossings"])
    marks=pop!(doc,"region_marks",nothing)
    added=_ms_add_crossings!(doc,source)
    marks!==nothing && (doc["region_marks"]=marks)
    disks=_mh_disk_policy!(doc,before,policy)
    history=vcat(first_history,[second_refinement.witness])
    inverse=Dict{String,Any}("source"=>source,"target"=>target,"forward_specification"=>spec,"construction"=>history,
        "core_graph"=>Dict(label=>[[d.face,d.u,d.v] for d in ds] for (label,ds) in paths),
        "old_source_refined_word"=>old_word,"removed_crossings"=>removed,"added_crossings"=>added,
        "removed_crossing_positions"=>removed_positions,
        "pair_of_pants_triangles"=>[collect(c) for c in sort(collect(neighborhood.selected))],
        "pair_of_pants_boundary_cycles"=>neighborhood.cycles,"selected_boundary"=>neighborhood.selected_boundary,
        "refined_document"=>_mm_hash(before.document),"source_document"=>_mm_hash(original),"output_document"=>_mm_hash(doc))
    correspondence=Dict("construction"=>history,"source_curve"=>source,"target_curve"=>target,
        "pair_of_pants_triangles"=>[collect(c) for c in sort(collect(neighborhood.selected))],
        "pair_of_pants_boundary_cycles"=>neighborhood.cycles,"boundary_tags"=>neighborhood.boundary_tags,
        "selected_boundary"=>neighborhood.selected_boundary,"core_graph"=>inverse["core_graph"],
        "old_source_word"=>old_word,"new_source_word"=>M.print_word(neighborhood.darts),
        "removed_crossings"=>removed,"added_crossings"=>added,"marked_regions"=>disks,
        "infinity_semantics"=>"transported protected source point; slide may change its complementary region")
    witness=_mm_witness("handle_slide",original,doc,spec,inverse,correspondence);witness["inverse_kind"]="handle_slide_inverse"
    (;document=doc,witness)
end

function _mh_core_cycle(records,st)
    M=MarkedSurfaceDiagrams
    ds=[M.Dart(x[1],x[2],x[3]) for x in records]
    _mm_need(!isempty(ds),"inverse source core graph is empty")
    outgoing=Dict();incoming=Dict()
    for d in ds
        _mm_need(haskey(st.incidence,d.face) && haskey(st.incidence[d.face],M.edge(d.u,d.v)),"inverse source core graph references a missing edge")
        u=M.root!(st.vertex_dsu,(d.face,d.u));v=M.root!(st.vertex_dsu,(d.face,d.v))
        _mm_need(!haskey(outgoing,u) && !haskey(incoming,v),"inverse source core graph is not an oriented embedded cycle")
        outgoing[u]=d;incoming[v]=d
    end
    _mm_need(Set(keys(outgoing))==Set(keys(incoming)),"inverse source core graph is not closed")
    result=M.Dart[];start=first(ds);d=start
    while isempty(result) || d!=start
        _mm_need(!(d in result),"inverse source core graph has multiple cycles")
        push!(result,d);v=M.root!(st.vertex_dsu,(d.face,d.v));d=outgoing[v]
    end
    _mm_need(length(result)==length(ds),"inverse source core graph is disconnected")
    result
end

"""Invert a generated slide through its verified labelled pair-of-pants core.

The old source is reconstructed as the graph's oriented cycle, refined through
the recorded edge map. The saved word is lexical formatting only and must
parse to that constructed route. The current frontier, restored intersections,
and entire neighborhood are checked before genuine mesh collapses. Independent
forward regeneration must reproduce the current target and all inverse data.
This is a generated geometric inverse, not an arbitrary inverse-band search.
"""
function marked_inverse_handle_slide(input::AbstractDict,specification::AbstractDict)
    _md_guard_periodic_operation(input)
    current=_marked_document(input);spec=_plain_value(specification);M=MarkedSurfaceDiagrams
    _mm_need(_mm_hash(current)==spec["output_document"],"inverse slide target has changed")
    cells=Set((p["face"],i-1) for p in current["patches"] for i in eachindex(p["triangles"]))
    _mm_need(all(c->_mm_cell(c) in cells,spec["pair_of_pants_triangles"]),"inverse neighborhood references a missing target triangle")
    source=spec["source"];target=spec["target"];history=spec["construction"]
    _mm_need(length(history)==2,"inverse slide needs its exact atlas and one-refinement construction")
    # Reconstruct coarse mesh incidence from subdivision provenance, without
    # supplying an original curve. Its empty curves are legitimate surface data.
    coarse_mesh=deepcopy(current);empty!(coarse_mesh["curves"]);empty!(coarse_mesh["crossings"]);delete!(coarse_mesh,"region_marks")
    refinement=last(history);atlas_history=history[1:end-1]
    provenance=refinement["inverse_specification"]["patches"]
    for p in provenance
        f=p["face"];patch=only(x for x in coarse_mesh["patches"] if x["face"]==f)
        original_vertices=Set(v for t in p["triangles"] for v in t["corners"])
        patch["triangles"]=[copy(t["corners"]) for t in p["triangles"]]
        for side in patch["side_order"];patch["sides"][side]=filter(v->v in original_vertices,patch["sides"][side]);end
        p["side_parameters"]===nothing ? delete!(patch,"side_parameters") : (patch["side_parameters"]=deepcopy(p["side_parameters"]))
    end
    # Source: marked_topology.jl central infinity-triangle guard. This temporary
    # coarse mesh uses its recorded coarse point, not the fine child index.
    # The actual restored document's metadata still travels through real collapse.
    if _mm_infinity(coarse_mesh)!==nothing
        coarse_mesh["metadata"]["display"]["infinity_triangle"]=deepcopy(refinement["inverse_specification"]["infinity_triangle"])
    end
    coarse_state=_marked_states(coarse_mesh).geometry
    core=_mh_core_cycle(spec["core_graph"]["source"],coarse_state)
    refined_core=_mh_transport_darts(core,[refinement])
    old_word=spec["old_source_refined_word"]
    _mm_need(M.parse_word(old_word)==refined_core,"inverse lexical word changes the reconstructed graph core")
    restored=deepcopy(current);curve=only(c for c in restored["curves"] if c["id"]==source);curve["word"]=old_word
    added=Set(spec["added_crossings"])
    _mm_need(all(x->x["id"] in added || !(source in x["curves"]),restored["crossings"]),"target has unexpected source-frontier intersections")
    filter!(x->!(x["id"] in added),restored["crossings"])
    marks=pop!(restored,"region_marks",nothing)
    generated_ids=_ms_add_crossings!(restored,source)
    actual=Dict(Tuple(x["at"])=>x for x in restored["crossings"] if x["id"] in generated_ids)
    _mm_need(length(actual)==length(spec["removed_crossings"]),"inverse core has a different number of intersections")
    for expected in spec["removed_crossings"]
        key=Tuple(expected["at"])
        _mm_need(haskey(actual,key),"inverse core does not recover an original crossing point")
        x=actual[key]
        _mm_need(Set(x["curves"])==Set(expected["curves"]),"inverse core crossing has different curve incidence")
        sign=x["curves"]==expected["curves"] ? x["sign"] : -x["sign"]
        _mm_need(sign==expected["sign"],"inverse core crossing orientation differs")
    end
    filter!(x->!(x["id"] in generated_ids),restored["crossings"])
    for record in spec["removed_crossing_positions"]
        insert!(restored["crossings"],record["index"],deepcopy(record["crossing"]))
    end
    marks!==nothing && (restored["region_marks"]=marks)
    _marked_states(restored)
    _mm_need(_mm_hash(restored)==spec["refined_document"],"constructed inverse core does not restore the recorded refined source")
    coarse_result=marked_collapse_subdivision(restored,refinement["inverse_specification"])
    coarse=_marked_states(coarse_result.document).geometry
    paths=Dict(label=>[M.Dart(x[1],x[2],x[3]) for x in records] for (label,records) in spec["core_graph"])
    _mm_need(paths["source"]==core && paths["target"]==coarse.curves[target][2],"inverse labelled cores disagree with actual coarse curves")
    _mh_band(coarse,source,target,M.print_word(paths["band"]))
    neighborhood=_mh_neighborhood(coarse,_marked_states(restored).geometry,paths)
    _mm_need([collect(c) for c in sort(collect(neighborhood.selected))]==spec["pair_of_pants_triangles"] &&
        _canonical_json(neighborhood.cycles)==_canonical_json(spec["pair_of_pants_boundary_cycles"]) &&
        neighborhood.selected_boundary==spec["selected_boundary"],"inverse pair-of-pants neighborhood differs from the witness")
    _mm_need(M.print_word(neighborhood.darts)==only(c for c in current["curves"] if c["id"]==source)["word"],"current source is not the generated band-sum frontier")
    original=coarse_result.document
    for w in reverse(atlas_history);original=replay_marked_move(original,inverse_marked_move(w)).document;end
    _mm_need(_mm_hash(original)==spec["source_document"],"inverse atlas/refinement history does not restore the source")
    generated=marked_handle_slide(original,spec["forward_specification"])
    _mm_need(generated.document==current && _canonical_json(generated.witness["inverse_specification"])==_canonical_json(spec),"forward regeneration does not reproduce the current slide/witness")
    # Source: the independently regenerated forward witness above. Its complete
    # region-cell comparison lives in the forward prepared common refinement;
    # it is NOT an inverse-native cell map. Region equality/non-equality is
    # symmetric, so this scalar decoration relation is valid in reverse too.
    forward_regions=deepcopy(generated.witness["correspondence"]["marked_regions"])
    decoration=Dict("direction"=>"inverse","cell_scope"=>"forward prepared common refinement; not native inverse endpoint cell indices",
        "forward_regions"=>forward_regions,"region_count"=>length(forward_regions),
        "all_regions_preserved"=>isempty(forward_regions) ? nothing : all(r["region_preserved"] for r in forward_regions),
        "policy"=>get(spec["forward_specification"],"disk_policy","preserve"))
    witness=_mm_witness("handle_slide_inverse",current,original,spec,spec["forward_specification"],
        Dict("pair_of_pants_triangles"=>spec["pair_of_pants_triangles"],"reconstructed_core"=>[[d.face,d.u,d.v] for d in core],"inverse_construction"=>history,"decoration_evidence"=>decoration))
    witness["inverse_kind"]="handle_slide"
    (;document=original,witness)
end
