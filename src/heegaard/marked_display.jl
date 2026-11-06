"""
    prepare_marked_display(doc; cut_curves)

Geometry-free, explicit triangle-corner quotient of a marked diagram. Only named
auxiliary cuts are left open; all other patch seams are glued back. The caller
must supply the marking: signed alpha/beta cyclic words alone are rejected.
Triangle indices in the returned interchange record are zero based.
"""
function prepare_marked_display(doc::Dict;cut_curves,visible_curves=nothing)
    M=MarkedSurfaceDiagrams
    M.need(haskey(doc,"patches") && haskey(doc,"surface") && haskey(doc["surface"],"marking_id"),
        "MARKING_REQUIRED","full disk cellulation required; routing will not be inferred")
    # The legacy kernel interprets region marks relative to ALL declared curves.
    # Reference curves must not split the semantic D: validate geometry without
    # marks, then validate the SAME marks on alpha U beta only.
    geometry=deepcopy(doc);delete!(geometry,"region_marks");st=M.validate(geometry)
    mathdoc=deepcopy(doc)
    filter!(c->c["role"] in ("alpha","beta"),mathdoc["curves"])
    mathematical_ids=Set(c["id"] for c in mathdoc["curves"])
    filter!(x->issubset(Set(x["curves"]),mathematical_ids),mathdoc["crossings"])
    M.validate(mathdoc)
    g=st.surface_report["genus"];cuts=String.(collect(cut_curves));cs=Set(cuts)
    M.need(g>=0 && st.surface_report["boundary_components"]==0,"DISPLAY_SURFACE","standard display requires a closed surface")
    M.need(length(cuts)==g && length(cs)==g,"DISPLAY_CUT_COUNT",g)
    M.need(all(haskey(st.curves,c) && st.curves[c][1]["role"]=="aux" && st.curves[c][1]["closed"] for c in cuts),
        "DISPLAY_CUT_ROLE","cuts must be explicit auxiliary curves, independent of alpha/beta")
    M.need(all(length(intersect(Set(x["curves"]),cs))<2 for x in doc["crossings"]),"DISPLAY_CUT_INTERSECTION",cuts)
    blocked=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),cs))
    rs=M.regions(st;barriers=blocked)
    M.need(length(rs)==1 && rs[1]["genus"]==0 && rs[1]["boundary_components"]==2g,"DISPLAY_NOT_PUNCTURED_SPHERE",cuts)
    pairs=M.fine_edge_occurrences(st)
    corners=[(f,i,j) for (f,i) in keys(st.triangles) for j in 1:3];vd=M.UnionFind(corners)
    ed=M.UnionFind(collect(keys(pairs)))
    corner(a,v)=(a[1],a[2],findfirst(==(v),st.triangles[a[1:2]]))
    blocked_at(a)=(a[1],M.edge(a[3],a[4])) in blocked
    for (a,b) in pairs
        if !blocked_at(a) && !blocked_at(b)
            M.unite!(ed,a,b)
            M.unite!(vd,corner(a,a[3]),corner(b,b[4]));M.unite!(vd,corner(a,a[4]),corner(b,b[3]))
        end
    end
    roots=sort(unique(M.root!(vd,c) for c in corners));vi=Dict(r=>i-1 for (i,r) in enumerate(roots))
    edge_roots=sort(unique(M.root!(ed,a) for a in keys(pairs)))
    ei=Dict(r=>i-1 for (i,r) in enumerate(edge_roots))
    edge_id(a)=ei[M.root!(ed,a)]
    edge_sign(a)=a==M.root!(ed,a) ? 1 : -1
    vertex(a,v)=vi[M.root!(vd,corner(a,v))]
    triangles=[Dict("source"=>[f,i-1],"vertices"=>[vi[M.root!(vd,(f,i,j))] for j in 1:3],
        "edges"=>[edge_id((f,i,u,v)) for (u,v) in M.tedges(st.triangles[(f,i)])],
        "edge_signs"=>[edge_sign((f,i,u,v)) for (u,v) in M.tedges(st.triangles[(f,i)])])
        for (f,i) in sort(collect(keys(st.triangles)))]
    boundaries=[];vertex_pairs=Dict();seen_sides=Set()
    for cycle in rs[1]["boundary_cycles"]
        labels=Set(c for (f,_,u,v) in cycle for c in st.used[(f,M.edge(u,v))] if c in cs)
        M.need(length(labels)==1,"DISPLAY_BOUNDARY_LABEL",labels);cid=only(labels)
        f,iz,u,v=cycle[1];d=only(filter(d->d.face==f && M.edge(d.u,d.v)==M.edge(u,v),st.curves[cid][2]))
        side=(d.u,d.v)==(u,v) ? "left" : "right"
        M.need(!((cid,side) in seen_sides),"DISPLAY_BOUNDARY_SIDE",(cid,side));push!(seen_sides,(cid,side))
        vs=[];es=Int[];ss=Int[]
        for (f,iz,u,v) in cycle
            a=(f,iz+1,u,v);w=vertex(a,u);push!(vs,w)
            push!(es,edge_id(a))
            push!(ss,edge_sign(a))
            key=M.root!(st.vertex_dsu,(f,u));push!(get!(vertex_pairs,key,Set{Int}()),w)
        end
        push!(boundaries,Dict("cut"=>cid,"side"=>side,"row"=>(side=="right" ? 0 : 1),
            "column"=>findfirst(==(cid),cuts)-1,"vertices"=>vs,"edges"=>es,"edge_signs"=>ss,"source"=>cycle))
    end
    # Each cut boundary has an exact vertex-by-vertex involution. Endpoint labels
    # are NOT a replacement for this pairing or for the retained disk embedding.
    pairing=[sort(collect(vs)) for (_,vs) in sort(collect(vertex_pairs);by=first)]
    M.need(all(length(p)==2 for p in pairing),"DISPLAY_PORT_PAIRING","two copies per boundary vertex")
    selected=visible_curves===nothing ? Set(c["id"] for c in doc["curves"] if c["role"] in ("alpha","beta")) : Set(String.(visible_curves))
    M.need(isempty(intersect(selected,cs)) && all(haskey(st.curves,c) for c in selected),"DISPLAY_CURVE_SELECTION",selected)
    curves=[]
    for c in doc["curves"]
        c["id"] in selected || continue
        segments=[]
        for d in st.curves[c["id"]][2]
            oc=first(st.incidence[d.face][M.edge(d.u,d.v)]);a=(d.face,oc[1],d.u,d.v)
            push!(segments,Dict("source"=>[d.face,d.u,d.v],"vertices"=>[vertex(a,d.u),vertex(a,d.v)],
                "edge"=>edge_id((d.face,oc[1],oc[2],oc[3])),
                "edge_sign"=>edge_sign((d.face,oc[1],oc[2],oc[3]))*(d.u==oc[2] ? 1 : -1)))
        end
        push!(curves,Dict("id"=>c["id"],"role"=>c["role"],"closed"=>c["closed"],"word"=>c["word"],
            "segments"=>segments,"seam_word"=>M.seam_word(st,c["id"])))
    end
    target_ids=Set(c["id"] for c in curves)
    xs=[x for x in doc["crossings"] if issubset(Set(x["curves"]),target_ids)]
    # D is a region of alpha U beta, NOT a region after auxiliary cuts. Retain
    # every source triangle and its cut pieces, even if a display cut splits D.
    math_edges=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),mathematical_ids))
    math_regions=M.regions(st;barriers=math_edges)
    split_regions=M.regions(st;barriers=union(math_edges,blocked));marks=[]
    for mark in get(doc,"region_marks",[])
        r=only(filter(r->mark["seed"] in r["triangles"],math_regions))
        cells=Set(Tuple(t) for t in r["triangles"])
        pieces=[p for p in split_regions if Tuple(p["seed"]) in cells]
        M.need(Set(Tuple(t) for p in pieces for t in p["triangles"])==cells,"DISPLAY_REGION_COVER",mark["id"])
        push!(marks,Dict("id"=>mark["id"],"claim"=>mark["claim"],"seed"=>mark["seed"],
            "source_region"=>r,"cut_pieces"=>pieces))
    end
    # Source: MarkedSurfaceDiagrams.regions boundary_cycles (interior on the
    # left), and the cut-side cycles above. A fill follows existing curve words
    # and cut circles; it never chooses an independent path between endpoints.
    # Keep the triangle partition as well as these ordered darts for read-back.
    curve_darts=Dict()
    for c in curves, s in c["segments"]
        f,u,v=s["source"]
        curve_darts[(f,u,v)]=Dict("kind"=>"curve","curve"=>c["id"],"source"=>s["source"],"direction"=>1)
        curve_darts[(f,v,u)]=Dict("kind"=>"curve","curve"=>c["id"],"source"=>s["source"],"direction"=>-1)
    end
    cut_darts=Dict(Tuple(d)=>(b,k) for b in boundaries for (k,d) in enumerate(b["source"]))
    for mark in marks, piece in mark["cut_pieces"]
        cycles=[]
        for cycle in piece["boundary_cycles"]
            first_curve=findfirst(d->!haskey(cut_darts,Tuple(d)),cycle)
            M.need(first_curve!==nothing,"DISPLAY_D_BOUNDARY","a disk piece requires a mathematical boundary interval")
            ordered=vcat(cycle[first_curve:end],cycle[1:first_curve-1]);parts=[];i=1
            while i<=length(ordered)
                d=ordered[i];f,_,u,v=d
                if haskey(cut_darts,Tuple(d))
                    b,k=cut_darts[Tuple(d)];run=[];j=i
                    while j<=length(ordered) && haskey(cut_darts,Tuple(ordered[j]))
                        bb,kk=cut_darts[Tuple(ordered[j])]
                        bb===b || break
                        M.need(kk==mod1(k+j-i,length(b["source"])),"DISPLAY_D_CUT_ORDER",mark["id"])
                        push!(run,ordered[j]);j+=1
                    end
                    push!(parts,Dict("kind"=>"cut","cut"=>b["cut"],"side"=>b["side"],
                        "from_source_vertex"=>b["vertices"][k],
                        "to_source_vertex"=>b["vertices"][mod1(k+length(run),length(b["vertices"]))],
                        "source_darts"=>run));i=j
                else
                    M.need(haskey(curve_darts,(f,u,v)),"DISPLAY_D_VISIBILITY","every curve bounding D must be visible")
                    push!(parts,curve_darts[(f,u,v)]);i+=1
                end
            end
            push!(cycles,parts)
        end
        piece["display_boundary_cycles"]=cycles
    end
    Dict("schema"=>"marked-display-complex/1","source_representation_sha256"=>M.representation_digest(doc),
        "marking_id"=>doc["surface"]["marking_id"],"genus"=>g,"cut_curves"=>cuts,
        "visible_curves"=>[c["id"] for c in curves],
        "grid"=>Dict("rows"=>2,"columns"=>g,"equal_radii"=>true),
        "cut_topology"=>Dict("genus"=>0,"boundary_components"=>2g,"components"=>1),
        "reference_disks"=>doc["patches"],"reference_gluings"=>doc["surface"]["gluings"],
        "triangles"=>triangles,"boundaries"=>boundaries,"boundary_vertex_pairs"=>pairing,
        "source_cut_corner_classes"=>Dict(join([f,i-1,j-1],",")=>vi[M.root!(vd,(f,i,j))] for (f,i,j) in corners),
        "curves"=>curves,"crossings"=>xs,"region_marks"=>marks,
        "layout_policy"=>"coordinates may realize this embedding; they may not choose or replace routes")
end

"""Recover a retained marked source only after checking its displayed incidence.

This checks source edge words, cuts, crossing signs and the entire D boundary.
It does not certify arbitrary SVG coordinates or infer a source from geometry.
"""
function marked_display_source(input::AbstractDict)
    d=_plain_value(input)
    get(d,"schema","")=="marked-cut-view/1" || throw(ArgumentError("Expected marked-cut-view/1"))
    haskey(d,"source_document") || throw(ArgumentError("The full marked source document is missing"))
    doc=d["source_document"]
    p=prepare_marked_display(doc;cut_curves=d["cuts"],visible_curves=get(d,"visible_curves",nothing))
    check(ok,message)=ok || throw(ArgumentError("Display/source mismatch: "*message))
    check(d["genus"]==p["genus"],"genus")
    check(get(d,"marking_id",nothing)==p["marking_id"],"marking identifier")
    check(d["preparation_representation_sha256"]==p["source_representation_sha256"],"source digest")
    check(d["source_cut_corner_classes"]==p["source_cut_corner_classes"],"cut quotient")
    check(d["boundary_vertex_pairs"]==p["boundary_vertex_pairs"],"boundary pairing")
    curve_edges(x)=Dict(c["id"]=>(c["role"],[s["source"] for s in c["segments"]]) for c in x["curves"])
    check(curve_edges(d)==curve_edges(p),"curve edge words")
    crossings(x)=[Dict(k=>c[k] for k in ("id","at","curves","sign")) for c in x["crossings"]]
    check(crossings(d)==crossings(p),"crossings")
    boundary_key(b)=(b["cut"],b["side"],b["row"],b["column"],b["source"])
    check([boundary_key(b) for b in d["boundaries"]]==[boundary_key(b) for b in p["boundaries"]],"cut boundaries")
    if haskey(d,"display_request")
        request=d["display_request"]
        semantic(c)=([c[k] for k in ("id","role","closed","word","seam_word")],
            [[s[k] for k in ("source","edge","edge_sign")] for s in c["segments"]])
        check(semantic.(d["curves"])==semantic.(p["curves"]),"directed source curve intervals")
        check(request["cut_curves"]==d["cuts"] && request["visible_curves"]==d["visible_curves"],"display request")
        outer=request["infinity_triangle"]
        check(outer==d["infinity_source_triangle"] && any(t["source"]==outer for t in p["triangles"]),"source infinity triangle")
        check(get(d,"disk_region_marks",nothing)==p["region_marks"],"D source partition and directed boundaries")
        expected=[Dict("mark"=>mark["id"],"source_triangles"=>piece["triangles"],
            "boundary_cycles"=>piece["display_boundary_cycles"],"contains_infinity"=>(outer in piece["triangles"]))
            for mark in p["region_marks"] for piece in mark["cut_pieces"]]
        if haskey(d,"disk_pieces")
            check(d["disk_pieces"]==expected,"D fill source intervals")
            check(isempty(get(d,"disk",[])),"unexpected additional D boundary")
            return doc
        end
    end
    if isempty(p["region_marks"])
        check(isempty(get(d,"disk",[])) && isempty(get(d,"disk_pieces",[])),"unexpected D boundary")
        return doc
    end
    length(p["region_marks"])==1 || throw(ArgumentError("Select one stored D for marked display conversion"))
    mark=only(p["region_marks"])
    length(mark["cut_pieces"])==1 || throw(ArgumentError("Split-D display conversion is not supported yet"))
    paths=Dict(Tuple(s["source"])=>s["vertices"] for c in d["curves"] for s in c["segments"])
    boundary=Int[];last_vertex=nothing
    for (f,_,u,v) in only(mark["source_region"]["boundary_cycles"])
        vs=haskey(paths,(f,u,v)) ? paths[(f,u,v)] : reverse(paths[(f,v,u)])
        check(last_vertex===nothing || last_vertex==first(vs),"D has a disconnected source path")
        append!(boundary,vs[1:end-1]);last_vertex=last(vs)
    end
    check(!isempty(boundary) && last_vertex==first(boundary),"D is not closed")
    check(_cyclic_equal(boundary,Int.(d["disk"])),"oriented D boundary")
    doc
end
