"""
    prepare_marked_display(doc; cut_curves)

Geometry-free, explicit triangle-corner quotient of a marked diagram. Only named
auxiliary cuts are left open; all other patch seams are glued back. The caller
must supply the marking: signed alpha/beta cyclic words alone are rejected.
Triangle indices in the returned interchange record are zero based.
"""
function prepare_marked_display(doc::Dict;cut_curves)
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
    M.need(g>=1 && st.surface_report["boundary_components"]==0,"DISPLAY_SURFACE","closed positive-genus surface")
    M.need(length(cuts)==g && length(cs)==g,"DISPLAY_CUT_COUNT",g)
    M.need(all(haskey(st.curves,c) && st.curves[c][1]["role"]=="aux" && st.curves[c][1]["closed"] for c in cuts),
        "DISPLAY_CUT_ROLE","cuts must be explicit auxiliary curves, independent of alpha/beta")
    M.need(all(length(intersect(Set(x["curves"]),cs))<2 for x in doc["crossings"]),"DISPLAY_CUT_INTERSECTION",cuts)
    blocked=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),cs))
    rs=M.regions(st;barriers=blocked)
    M.need(length(rs)==1 && rs[1]["genus"]==0 && rs[1]["boundary_components"]==2g,"DISPLAY_NOT_PUNCTURED_SPHERE",cuts)
    pairs=M.fine_edge_occurrences(st)
    corners=[(f,i,j) for (f,i) in keys(st.triangles) for j in 1:3];vd=M.UnionFind(corners)
    corner(a,v)=(a[1],a[2],findfirst(==(v),st.triangles[a[1:2]]))
    blocked_at(a)=(a[1],M.edge(a[3],a[4])) in blocked
    for (a,b) in pairs
        if !blocked_at(a) && !blocked_at(b)
            M.unite!(vd,corner(a,a[3]),corner(b,b[4]));M.unite!(vd,corner(a,a[4]),corner(b,b[3]))
        end
    end
    roots=sort(unique(M.root!(vd,c) for c in corners));vi=Dict(r=>i-1 for (i,r) in enumerate(roots))
    vertex(a,v)=vi[M.root!(vd,corner(a,v))]
    triangles=[Dict("source"=>[f,i-1],"vertices"=>[vi[M.root!(vd,(f,i,j))] for j in 1:3]) for (f,i) in sort(collect(keys(st.triangles)))]
    boundaries=[];vertex_pairs=Dict();seen_sides=Set()
    for cycle in rs[1]["boundary_cycles"]
        labels=Set(c for (f,_,u,v) in cycle for c in st.used[(f,M.edge(u,v))] if c in cs)
        M.need(length(labels)==1,"DISPLAY_BOUNDARY_LABEL",labels);cid=only(labels)
        f,iz,u,v=cycle[1];d=only(filter(d->d.face==f && M.edge(d.u,d.v)==M.edge(u,v),st.curves[cid][2]))
        side=(d.u,d.v)==(u,v) ? "left" : "right"
        M.need(!((cid,side) in seen_sides),"DISPLAY_BOUNDARY_SIDE",(cid,side));push!(seen_sides,(cid,side))
        vs=[]
        for (f,iz,u,v) in cycle
            a=(f,iz+1,u,v);w=vertex(a,u);push!(vs,w)
            key=M.root!(st.vertex_dsu,(f,u));push!(get!(vertex_pairs,key,Set{Int}()),w)
        end
        push!(boundaries,Dict("cut"=>cid,"side"=>side,"row"=>(side=="right" ? 0 : 1),
            "column"=>findfirst(==(cid),cuts)-1,"vertices"=>vs,"source"=>cycle))
    end
    # Each cut boundary has an exact vertex-by-vertex involution. Endpoint labels
    # are NOT a replacement for this pairing or for the retained disk embedding.
    pairing=[sort(collect(vs)) for (_,vs) in sort(collect(vertex_pairs);by=first)]
    M.need(all(length(p)==2 for p in pairing),"DISPLAY_PORT_PAIRING","two copies per boundary vertex")
    curves=[]
    for c in doc["curves"]
        c["role"] in ("alpha","beta") || continue
        segments=[]
        for d in st.curves[c["id"]][2]
            oc=first(st.incidence[d.face][M.edge(d.u,d.v)]);a=(d.face,oc[1],d.u,d.v)
            push!(segments,Dict("source"=>[d.face,d.u,d.v],"vertices"=>[vertex(a,d.u),vertex(a,d.v)]))
        end
        push!(curves,Dict("id"=>c["id"],"role"=>c["role"],"closed"=>c["closed"],"word"=>c["word"],
            "segments"=>segments,"seam_word"=>M.seam_word(st,c["id"])))
    end
    target_ids=Set(c["id"] for c in curves)
    xs=[x for x in doc["crossings"] if issubset(Set(x["curves"]),target_ids)]
    # D is a region of alpha U beta, NOT a region after auxiliary cuts. Retain
    # every source triangle and its cut pieces, even if a display cut splits D.
    math_edges=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),target_ids))
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
    Dict("schema"=>"marked-display-complex/1","source_representation_sha256"=>M.representation_digest(doc),
        "marking_id"=>doc["surface"]["marking_id"],"genus"=>g,"cut_curves"=>cuts,
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
    p=prepare_marked_display(doc;cut_curves=d["cuts"])
    check(ok,message)=ok || throw(ArgumentError("Display/source mismatch: "*message))
    check(d["genus"]==p["genus"],"genus")
    check(d["preparation_representation_sha256"]==p["source_representation_sha256"],"source digest")
    check(d["source_cut_corner_classes"]==p["source_cut_corner_classes"],"cut quotient")
    check(d["boundary_vertex_pairs"]==p["boundary_vertex_pairs"],"boundary pairing")
    curve_edges(x)=Dict(c["id"]=>(c["role"],[s["source"] for s in c["segments"]]) for c in x["curves"])
    check(curve_edges(d)==curve_edges(p),"curve edge words")
    crossings(x)=[Dict(k=>c[k] for k in ("id","at","curves","sign")) for c in x["crossings"]]
    check(crossings(d)==crossings(p),"crossings")
    boundary_key(b)=(b["cut"],b["side"],b["row"],b["column"],b["source"])
    check([boundary_key(b) for b in d["boundaries"]]==[boundary_key(b) for b in p["boundaries"]],"cut boundaries")
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
