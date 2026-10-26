"""Whole marked-diagram Dehn twists on a fixed background surface.

Sources: marked_annulus_product.jl constructs the exact periodic product h;
marked_polygon_overlay.jl delegates exact convex intersection to Polyhedra;
marked_groups.jl supplies old-basis edge words. The half-supported shear is
F(u,t)=(u+s*max(t,0),t), fixing both frontiers modulo period and the core.
All alpha/beta/aux curves and selected regions move together. A beta-only
twist with alpha fixed is a different operation and is not implemented here.
"""
_mt_q(x)=MarkedSurfaceDiagrams.qnum(x)
_mt_json_point(p)=MarkedSurfaceDiagrams.qstr.(collect(p))
_mt_json_bary(p)=MarkedSurfaceDiagrams.qstr.(collect(p))
_mt_live_bary(p)=_mt_q.(p)
function _mt_bary(p,triangle)
    M=MarkedSurfaceDiagrams;den=M.orient(triangle...)
    _mm_need(den>0,"twist product triangle is not positively oriented")
    a,b,c=triangle
    weights=[M.orient(p,b,c)/den,M.orient(a,p,c)/den,M.orient(a,b,p)/den]
    _mm_need(sum(weights)==1 && all(>=(0),weights),"overlay point is outside its recorded triangle")
    weights
end
_mt_point(weights,points)=Tuple(sum(weights[j]*points[j][k] for j in 1:3) for k in 1:2)
_mt_shear(p,sign)=(p[1]+sign*max(p[2],0),p[2])
_mt_ratio(weights)=MarkedSurfaceDiagrams.orient([(w[2],w[3]) for w in weights]...)

"""Triangulate a convex overlay boundary without adding polygon centers.

Source: installed Polyhedra 0.8.3 triangulation.jl::triangulation and
redundancy.jl::removevredundancy! (exact planar hull). The library chooses
extreme-vertex triangles. The thin assembly below restores every labelled
collinear boundary point by subdividing its boundary-edge-owning triangle.
No clipping, hull or general polygon triangulator is reimplemented here.
"""
function _mt_convex_boundary_triangles(points)
    M=MarkedSurfaceDiagrams;Q=M.Rat;n=length(points)
    _mm_need(n>=3 && length(Set(points))==n && _mpo_area(points)>0,"overlay polygon boundary is invalid")
    _mm_need(all(M.orient(points[mod1(i-1,n)],points[i],points[mod1(i+1,n)])>=0 for i in 1:n),"overlay polygon is not convex counterclockwise")
    poly=Polyhedra.polyhedron(Polyhedra.vrep([Q[p[1],p[2]] for p in points]),Polyhedra.DefaultLibrary{Q}())
    Polyhedra.removevredundancy!(poly)
    index=Dict(p=>i for (i,p) in enumerate(points));triangles=Vector{Int}[]
    for piece in Polyhedra.triangulation(poly;tol=0)
        corners=[Tuple(Q.(p)) for p in Polyhedra.points(piece)]
        _mm_need(length(corners)==3 && all(haskey(index,p) for p in corners),"library convex triangulation added or lost a vertex")
        triangle=[index[p] for p in corners]
        M.orient((points[i] for i in triangle)...)<0 && (triangle=triangle[[1,3,2]])
        _mm_need(M.orient((points[i] for i in triangle)...)>0,"library convex triangulation produced a degenerate cell")
        push!(triangles,triangle)
    end
    extremes=[i for i in 1:n if M.orient(points[mod1(i-1,n)],points[i],points[mod1(i+1,n)])>0]
    for (a,b) in zip(extremes,circshift(extremes,-1))
        chain=[a];j=a
        while j!=b;j=mod1(j+1,n);push!(chain,j);end
        length(chain)<=2 && continue
        owners=findall(t->a in t && b in t,triangles)
        _mm_need(length(owners)==1,"convex boundary edge does not have a unique triangle owner")
        owner=only(owners);third=only(setdiff(triangles[owner],[a,b]));deleteat!(triangles,owner)
        for (u,v) in zip(chain[1:end-1],chain[2:end]);push!(triangles,[third,u,v]);end
    end
    edge_counts=Dict();area=Q(0)
    for triangle in triangles
        determinant=M.orient((points[i] for i in triangle)...)
        _mm_need(determinant>0,"collinear label reinsertion produced a degenerate triangle")
        area+=determinant/2
        for (a,b) in M.tedges(triangle);e=minmax(a,b);edge_counts[e]=get(edge_counts,e,0)+1;end
    end
    boundary=Set(minmax(i,mod1(i+1,n)) for i in 1:n)
    _mm_need(area==_mpo_area(points) && all(get(edge_counts,e,0)==1 for e in boundary) &&
        all(count==(e in boundary ? 1 : 2) for (e,count) in edge_counts),"convex cell triangulation lost boundary labels or coverage")
    triangles
end

function _mt_vertex_bary(f,v,cell,st,new_vertices)
    Q=MarkedSurfaceDiagrams.Rat;t=st.triangles[(f,cell[2]+1)];where=findfirst(==(v),t)
    where!==nothing && return Q[where==j ? 1 : 0 for j in 1:3]
    _mm_need(haskey(new_vertices,(f,v)),"common-refinement vertex has no old-cell ancestry")
    record=new_vertices[(f,v)]
    if record["kind"]=="edge"
        a,b=record["edge"];r=_mt_q(record["parameter"])
        _mm_need(0<r<1 && a in t && b in t,"new edge vertex is outside its target cell")
        return Q[t[j]==a ? 1-r : t[j]==b ? r : 0 for j in 1:3]
    end
    _mm_need(record["kind"]=="triangle" && _mm_cell(record["cell"])==cell,"new interior point has inconsistent target ancestry")
    weights=_mt_live_bary(record["weights"])
    _mm_need(length(weights)==3 && sum(weights)==1 && all(>(0),weights),"new interior vertex has invalid barycentric coordinates")
    weights
end

"""Reuse the verified annular constructor without retaining its probe curve."""
function _mt_prepare(input,core)
    doc=_marked_document(input);ids=Set(c["id"] for c in doc["curves"])
    n=1;probe="twist_support_$n"
    while probe in ids;n+=1;probe="twist_support_$n";end
    annulus=marked_parallel_copy(doc,Dict("curve"=>core,"id"=>probe,"side"=>"left"))
    fine=deepcopy(annulus.document);filter!(c->c["id"]!=probe,fine["curves"])
    removed=Set(annulus.witness["inverse_specification"]["new_crossings"])
    filter!(x->!(x["id"] in removed),fine["crossings"])
    _mm_need(_mm_hash(fine)==annulus.witness["inverse_specification"]["refined_document"],"annular probe removal changed source data")
    history=annulus.witness["correspondence"]["construction"]
    coarse=marked_collapse_subdivision(fine,last(history)["inverse_specification"]).document
    support=Set(_mm_cell(c) for c in annulus.witness["correspondence"]["annulus_triangles"])
    product=_marked_annulus_product(_marked_states(coarse).geometry,_marked_states(fine).geometry,core,support)
    (;fine,coarse,history,support,product)
end

function _mt_periodic_overlay(product,sign)
    Q=MarkedSurfaceDiagrams.Rat
    _mm_need(product.period==1,"twist product must use the certified unit angular period")
    targets=product.triangles;lower=minimum(p[1] for t in targets for p in t.points)
    upper=maximum(p[1] for t in targets for p in t.points);sources=Any[]
    for triangle in targets
        sheared=[_mt_shear(p,sign) for p in triangle.points]
        _mm_need(MarkedSurfaceDiagrams.orient(sheared...)>0,"annular shear folded a source cell")
        first_shift=ceil(Int,lower-maximum(p[1] for p in sheared))
        last_shift=floor(Int,upper-minimum(p[1] for p in sheared))
        for k in first_shift:last_shift
            push!(sources,(id=Dict("cell"=>deepcopy(triangle.id),"shift"=>k),
                points=[(p[1]+Q(k),p[2]) for p in sheared]))
        end
    end
    result=_marked_triangle_overlay(sources,targets)
    _mm_need(all(r.area_matches for r in result.coverage),"periodic image cells do not cover every fixed target triangle exactly in area")
    # No-overlap is supplied by the independently certified product chart and
    # the explicit injective periodic shear (u,t)->(u+phi(t),t), not area alone.
    result
end

function _mt_overlay_complex(fine_doc,product,overlay,sign)
    M=MarkedSurfaceDiagrams;Q=M.Rat;before=_marked_states(fine_doc);st=before.geometry
    doc=deepcopy(fine_doc);h=Dict(_mm_cell(t.id)=>t.points for t in product.triangles)
    support=Set(keys(h));used=Dict(f=>Set(v for t in p["triangles"] for v in t) for (f,p) in st.patches)
    edge_points=Dict((f,e)=>Dict(Q(0)=>e[1],Q(1)=>e[2]) for (f,inc) in st.incidence for e in keys(inc))
    interiors=Dict();new_vertices=Dict{Tuple{String,String},Any}();polygons=Dict{Tuple{String,Int},Vector{Any}}()
    function register(weights,cell)
        f,i=cell;t=st.triangles[(f,i+1)]
        index=findfirst(==(1),weights)
        index!==nothing && return t[index]
        zero_index=findfirst(==(0),weights)
        if zero_index!==nothing
            ends=[j for j in 1:3 if j!=zero_index];e=M.edge(t[ends[1]],t[ends[2]])
            r=weights[only(findall(==(e[2]),t))];points=edge_points[(f,e)]
            haskey(points,r) && return points[r]
            name=_ms_name!(used[f],"twist_vertex");points[r]=name
            new_vertices[(f,name)]=Dict("kind"=>"edge","edge"=>collect(e),"parameter"=>M.qstr(r))
            return name
        end
        key=(f,i,Tuple(weights));haskey(interiors,key) && return interiors[key]
        name=_ms_name!(used[f],"twist_vertex");interiors[key]=name
        new_vertices[(f,name)]=Dict("kind"=>"triangle","cell"=>collect(cell),"weights"=>_mt_json_bary(weights))
        name
    end
    bary_vertex(f,v,cell)=_mt_vertex_bary(f,v,cell,st,new_vertices)
    for cell in overlay.cells
        target=_mm_cell(cell.target_id);source=_mm_cell(cell.source_id["cell"]);shift=Int(cell.source_id["shift"])
        names=[register(_mt_bary(p,h[target]),target) for p in cell.vertices]
        triangles=_mt_convex_boundary_triangles(cell.vertices)
        push!(get!(polygons,target,Any[]),(;source,shift,names,triangles))
    end
    source_children=Dict(c=>Tuple{String,Int}[] for c in [(f,i-1) for (f,i) in keys(st.triangles)])
    target_children=Dict(c=>Tuple{String,Int}[] for c in keys(source_children))
    records=Any[];source_areas=Dict(c=>Q(0) for c in keys(source_children));target_areas=copy(source_areas)
    source_segments=Dict{Any,Vector{Any}}()
    function source_bary(weights,target,source,shift)
        target in support || return weights
        p=_mt_point(weights,h[target]);inverse=_mt_shear(p,-sign)
        _mt_bary((inverse[1]-shift,inverse[2]),h[source])
    end
    function record_triangle!(triangles,names,target,source,shift)
        f=target[1];tb=[bary_vertex(f,v,target) for v in names]
        sb=[source_bary(w,target,source,shift) for w in tb]
        ta=_mt_ratio(tb);sa=_mt_ratio(sb)
        _mm_need(ta>0 && sa>0,"common-refinement child is degenerate or folded")
        k=length(triangles);push!(triangles,names);child=(f,k)
        push!(source_children[source],child);push!(target_children[target],child)
        source_areas[source]+=sa;target_areas[target]+=ta
        push!(records,Dict("child"=>collect(child),"source"=>collect(source),"target"=>collect(target),"shift"=>shift,
            "vertices"=>copy(names),"source_weights"=>_mt_json_bary.(sb),"target_weights"=>_mt_json_bary.(tb)))
    end
    for patch in doc["patches"]
        f=patch["face"];triangles=Any[]
        for (i,t) in enumerate(st.patches[f]["triangles"])
            target=(f,i-1)
            if target in support
                _mm_need(haskey(polygons,target),"overlay lost a supported target cell")
                for polygon in polygons[target]
                    source,shift=polygon.source,polygon.shift
                    for indices in polygon.triangles
                        record_triangle!(triangles,polygon.names[indices],target,source,shift)
                    end
                    for (u,v) in zip(polygon.names,circshift(polygon.names,-1))
                        su=source_bary(bary_vertex(f,u,target),target,source,shift)
                        sv=source_bary(bary_vertex(f,v,target),target,source,shift)
                        common=[j for j in 1:3 if su[j]==0 && sv[j]==0]
                        if length(common)==1
                            old=st.triangles[(source[1],source[2]+1)];ends=[j for j in 1:3 if j!=only(common)]
                            edge=M.edge(old[ends[1]],old[ends[2]]);rindex=only(findall(==(edge[2]),old));a,b=su[rindex],sv[rindex]
                            a==b && continue
                            lo,hi=minmax(a,b);dart=a<b ? [f,u,v] : [f,v,u]
                            push!(get!(source_segments,(source,edge),Any[]),(;lo,hi,dart))
                        end
                    end
                end
            else
                boundary=String[]
                for (a,b) in M.tedges(t)
                    e=M.edge(a,b);points=edge_points[(f,e)];parameters=sort(collect(keys(points)))
                    a!=e[1] && reverse!(parameters)
                    append!(boundary,[points[r] for r in parameters[1:end-1]])
                end
                if length(boundary)==3
                    record_triangle!(triangles,copy(t),target,target,0)
                else
                    coordinates=[begin w=bary_vertex(f,v,target);(w[2],w[3]) end for v in boundary]
                    for indices in _mt_convex_boundary_triangles(coordinates)
                        record_triangle!(triangles,boundary[indices],target,target,0)
                    end
                end
            end
        end
        patch["triangles"]=triangles
        parameters=Dict{String,Any}()
        for side in patch["side_order"]
            old=st.patches[f]["sides"][side];ts=M.side_parameters(st.patches[f],side);path=String[];ps=String[]
            for j in 1:length(old)-1
                a,b=old[j:j+1];e=M.edge(a,b);points=edge_points[(f,e)];rs=sort(collect(keys(points)))
                a!=e[1] && reverse!(rs)
                for r in rs[1:end-1]
                    q=a==e[1] ? r : 1-r
                    push!(path,points[r]);push!(ps,M.qstr((1-q)*ts[j]+q*ts[j+1]))
                end
            end
            push!(path,last(old));push!(ps,M.qstr(last(ts)));patch["sides"][side]=path;parameters[side]=ps
        end
        patch["side_parameters"]=parameters
    end
    _mm_need(all(==(1),values(source_areas)) && all(==(1),values(target_areas)),"twist common-refinement cells do not cover source and fixed target triangles exactly")
    for mark in get(doc,"region_marks",[]);mark["seed"]=collect(first(source_children[_mm_cell(mark["seed"])]));end
    infinity=_mm_infinity(fine_doc)
    infinity!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=collect(first(source_children[_mm_cell(infinity)])))
    # This is an internal construction complex, not an accepted DiagramState.
    # A common refinement may exceed native resource limits. Existing local
    # disk incidence suffices to construct external image words; native bounds
    # are enforced only by the document/rechart acceptance boundary below.
    incidence=Dict(p["face"]=>M.patch_info(p)[2] for p in doc["patches"])
    image_records=Any[]
    for curve in doc["curves"]
        old=M.parse_word(curve["word"]);out=M.Dart[];pieces=Any[]
        for d in old
            edge=M.edge(d.u,d.v);inside=[(d.face,o[1]-1) for o in st.incidence[d.face][edge] if (d.face,o[1]-1) in support]
            image=M.Dart[];intervals=Any[]
            if isempty(inside)
                points=edge_points[(d.face,edge)];rs=sort(collect(keys(points)));d.u!=edge[1] && reverse!(rs)
                image=[M.Dart(d.face,points[a],points[b]) for (a,b) in zip(rs[1:end-1],rs[2:end])]
            else
                source=first(inside);unique_segments=Dict()
                for s in get(source_segments,(source,edge),Any[])
                    f,u,v=s.dart
                    length(incidence[f][M.edge(u,v)])==2 || continue
                    unique_segments[(s.lo,s.hi,f,u,v)]=s
                end
                segments=sort(collect(values(unique_segments));by=s->(s.lo,s.hi,Tuple(s.dart)))
                _mm_need(!isempty(segments) && first(segments).lo==0 && last(segments).hi==1,"twist curve image cannot be realized off the marking sides")
                _mm_need(all(a.hi==b.lo for (a,b) in zip(segments[1:end-1],segments[2:end])),"twist curve image has an ambiguous or incomplete edge partition")
                image=[M.Dart(s.dart...) for s in segments]
                intervals=[Dict("from"=>M.qstr(s.lo),"to"=>M.qstr(s.hi),"dart"=>s.dart) for s in segments]
                d.u!=edge[1] && (image=[M.Dart(x.face,x.v,x.u) for x in reverse(image)])
            end
            append!(out,image);push!(pieces,Dict("source"=>[d.face,d.u,d.v],"target"=>[[x.face,x.u,x.v] for x in image],"partition"=>intervals))
        end
        push!(image_records,Dict("curve"=>curve["id"],"source_word_format"=>curve["word"],"edge_images"=>pieces))
        curve["word"]=M.print_word(out)
    end
    vertex_records=[Dict("vertex"=>[f,v],"ancestry"=>r) for ((f,v),r) in sort(collect(new_vertices);by=first)]
    (;construction_only=true,complex=doc,records,vertex_records,image_records,source_children,target_children)
end

function _mt_weight_product(w,rows)
    [[sum(w[i][j]*rows[j][k] for j in 1:3) for k in 1:3] for i in 1:3]
end

function _mt_pullback_interval(d,source_dart,new,records,old)
    M=MarkedSurfaceDiagrams;f,u,v=source_dart.face,source_dart.u,source_dart.v;intervals=Set()
    for (i,_,_) in new.incidence[d.face][M.edge(d.u,d.v)]
        child=records[(d.face,i-1)];parent=_mm_cell(child["source"])
        parent[1]==f || continue;t=old.triangles[(f,parent[2]+1)]
        u in t && v in t || continue
        a=only(findall(==(d.u),child["vertices"]));b=only(findall(==(d.v),child["vertices"]))
        wa=_mt_live_bary(child["source_weights"][a]);wb=_mt_live_bary(child["source_weights"][b])
        outside=only(j for j in 1:3 if t[j]!=u && t[j]!=v)
        wa[outside]==0 && wb[outside]==0 || continue
        vi=only(findall(==(v),t));push!(intervals,(wa[vi],wb[vi]))
    end
    _mm_need(length(intervals)==1,"curve edge does not map to one original directed source edge")
    only(intervals)
end

function _mt_rechart_overlay(fine_doc,raw)
    M=MarkedSurfaceDiagrams;Q=M.Rat
    geometry=Dict("surface"=>deepcopy(raw.complex["surface"]),"patches"=>deepcopy(raw.complex["patches"]))
    # Source: marked_rechart.jl's default uniform construction and explicit
    # proper_arc identity partition. A disk with V <= 256 has F <= 507;
    # 64 native disks therefore cannot contain more than 32448 triangles.
    # Preserve the accepted empty-spec path unless six flags already exceed
    # that absolute bound. Features stay external until geometry is valid.
    if big(6)*sum(length(p["triangles"]) for p in geometry["patches"]) > 64*507
        # Source: the proper-arc helper retains old macro sides and avoids all
        # old vertices. Its interior-edge images are legal curve routes only
        # when the ACTUAL raw features already avoid old seams and corners.
        charts=Dict()
        for p in geometry["patches"]
            _,inc,vertices=M.patch_info(p)
            charts[p["face"]]=(;incidence=inc,vertices=Set(vertices),
                corners=Set(first(p["sides"][s]) for s in p["side_order"]),
                boundary=Set(v for path in values(p["sides"]) for v in path))
        end
        for c in raw.complex["curves"],d in M.parse_word(c["word"])
            _mm_need(haskey(charts,d.face) && haskey(charts[d.face].incidence,M.edge(d.u,d.v)),
                "proper-arc raw curve references a missing old edge")
            chart=charts[d.face]
            _mm_need(length(chart.incidence[M.edge(d.u,d.v)])==2,
                "proper-arc raw curves may not run along an old macro seam")
            _mm_need(!(d.u in chart.corners) && !(d.v in chart.corners),
                "proper-arc raw curves may not pass through an old macro corner")
        end
        for x in raw.complex["crossings"]
            f,v=x["at"]
            _mm_need(haskey(charts,f) && v in charts[f].vertices && !(v in charts[f].boundary),
                "proper-arc raw crossing must lie inside its old chart")
        end
        charted=_marked_rechart_geometry(geometry,_marked_rechart_opposites(geometry);
            specification=Dict("strategy"=>"proper_arc"))
    else
        charted=_marked_rechart_geometry(geometry,_marked_rechart_opposites(geometry))
    end
    corr=charted.correspondence;edge_map=Dict(Tuple(r["source"])=>r["target"] for r in corr["edge_images"])
    vertex_map=Dict(Tuple(r["source"])=>r["target"] for r in corr["vertex_images"])
    children=Dict(_mm_cell(r["source"])=>[_mm_cell(c) for c in r["target"]] for r in corr["cell_children"])
    raw_records=Dict(_mm_cell(r["child"])=>r for r in raw.records)
    records=Any[]
    for child in corr["child_barycentric_rows"]
        source=_mm_cell(child["source"]);parent=raw_records[source]
        weights=_mt_live_bary.(child["weights"])
        sw=_mt_weight_product(weights,_mt_live_bary.(parent["source_weights"]))
        tw=_mt_weight_product(weights,_mt_live_bary.(parent["target_weights"]))
        _mm_need(_mt_ratio(sw)>0 && _mt_ratio(tw)>0,"identity rechart folded the fixed source/target cell map")
        push!(records,Dict("child"=>deepcopy(child["target"]),"vertices"=>deepcopy(child["vertices"]),
            "source"=>deepcopy(parent["source"]),"target"=>deepcopy(parent["target"]),"shift"=>parent["shift"],
            "source_weights"=>_mt_json_bary.(sw),"target_weights"=>_mt_json_bary.(tw)))
    end
    doc=deepcopy(charted.document);doc["curves"]=deepcopy(raw.complex["curves"])
    image_records=deepcopy(raw.image_records)
    for (curve,image) in zip(doc["curves"],image_records)
        _mm_need(curve["id"]==image["curve"],"rechart curve correspondence changed order")
        out=M.Dart[]
        for piece in image["edge_images"]
            expanded=M.Dart[]
            for d in piece["target"]
                _mm_need(haskey(edge_map,Tuple(d)),"rechart omitted an external curve edge")
                append!(expanded,[M.Dart(x...) for x in edge_map[Tuple(d)]])
            end
            piece["proof_target"]=piece["target"]
            piece["proof_partition"]=piece["partition"]
            piece["target"]=[[d.face,d.u,d.v] for d in expanded]
            append!(out,expanded)
        end
        curve["word"]=M.print_word(out)
    end
    doc["crossings"]=deepcopy(raw.complex["crossings"])
    for x in doc["crossings"]
        _mm_need(haskey(vertex_map,Tuple(x["at"])),"rechart omitted an actual crossing vertex")
        x["at"]=deepcopy(vertex_map[Tuple(x["at"])])
    end
    if haskey(raw.complex,"metadata");doc["metadata"]=deepcopy(raw.complex["metadata"]);end
    if haskey(raw.complex,"region_marks")
        doc["region_marks"]=deepcopy(raw.complex["region_marks"])
        for mark in doc["region_marks"];mark["seed"]=collect(first(children[_mm_cell(mark["seed"])]));end
    end
    infinity=_mm_infinity(raw.complex)
    infinity!==nothing && (doc["metadata"]["display"]["infinity_triangle"]=collect(first(children[_mm_cell(infinity)])))
    after=_marked_states(doc);old=_marked_states(fine_doc).geometry
    by_child=Dict(_mm_cell(r["child"])=>r for r in records)
    for image in image_records,piece in image["edge_images"]
        source=M.Dart(piece["source"]...);partition=Any[];previous=Q(0)
        for x in piece["target"]
            d=M.Dart(x...);a,b=_mt_pullback_interval(d,source,after.geometry,by_child,old)
            _mm_need(a==previous && a<b,"recharted curve did not retain its oriented source-edge partition")
            push!(partition,Dict("from"=>M.qstr(a),"to"=>M.qstr(b),"dart"=>x));previous=b
        end
        _mm_need(previous==1,"rechart lost the end of a source curve edge")
        piece["partition"]=partition
    end
    source_children=Dict(c=>[x for r in raw.source_children[c] for x in children[r]] for c in keys(raw.source_children))
    target_children=Dict(c=>[x for r in raw.target_children[c] for x in children[r]] for c in keys(raw.target_children))
    source_area=Dict(c=>Q(0) for c in keys(source_children));target_area=Dict(c=>Q(0) for c in keys(target_children))
    for r in records
        source_area[_mm_cell(r["source"])]+=_mt_ratio(_mt_live_bary.(r["source_weights"]))
        target_area[_mm_cell(r["target"])]+=_mt_ratio(_mt_live_bary.(r["target_weights"]))
    end
    _mm_need(all(==(1),values(source_area)) && all(==(1),values(target_area)),"identity rechart lost complete fixed source/target cell coverage")
    (;document=doc,records,image_records,source_children,target_children,rechart=charted)
end

function _mt_overlay_document(fine_doc,product,overlay,sign)
    M=MarkedSurfaceDiagrams;raw=_mt_overlay_complex(fine_doc,product,overlay,sign)
    oversized=length(raw.complex["patches"])>64 || any(length(Set(v for t in p["triangles"] for v in t))>256 for p in raw.complex["patches"])
    mapped=oversized ? _mt_rechart_overlay(fine_doc,raw) :
        (;document=raw.complex,records=raw.records,image_records=raw.image_records,source_children=raw.source_children,target_children=raw.target_children,rechart=nothing)
    doc=mapped.document;before=_marked_states(fine_doc);after=_marked_states(doc)
    region_map=_ms_regions(M.regions(before.mathematics),M.regions(after.mathematics),mapped.source_children)
    geometry_regions=_ms_regions(M.regions(before.geometry),M.regions(after.geometry),mapped.source_children)
    (;document=doc,records=mapped.records,image_records=mapped.image_records,
        source_children=mapped.source_children,target_children=mapped.target_children,region_map,geometry_regions,raw,rechart=mapped.rechart)
end

function _mt_inverse_geometry(current,spec)
    M=MarkedSurfaceDiagrams;Q=M.Rat
    mesh=spec["fine_source_mesh"]
    old_doc=Dict("version"=>1,"surface"=>deepcopy(mesh["surface"]),"patches"=>deepcopy(mesh["patches"]),"curves"=>[],"crossings"=>[])
    old=_marked_states(old_doc).geometry
    # `current` here is the ACTUALLY unrecharted internal proof complex, not
    # an over-limit accepted native document. The public target was validated
    # by _marked_unrechart_geometry before features were pulled back. Reuse
    # existing local disk/old-seam checks, not an artificial DiagramState.
    proof=_mr_proof(Dict("surface"=>current["surface"],"patches"=>current["patches"]))
    vertex_keys=Set((f,v) for (f,p) in proof.patches for t in p["triangles"] for v in t)
    new=(;patches=proof.patches,incidence=proof.incidence,triangles=proof.triangles,vertex_keys)
    _mm_need(current["surface"]==mesh["surface"] && Set(keys(old.patches))==Set(keys(new.patches)),"twist inverse target changed its fixed background model")
    h=Dict(_mm_cell(r["cell"])=>[Tuple(_mt_q.(p)) for p in r["points"]] for r in spec["product_triangles"])
    _mm_need(Set(keys(h))==Set(_mm_cell(c) for c in spec["support_cells"]),"twist inverse product support differs")
    vertices=Dict{Tuple{String,String},Any}()
    for record in spec["new_vertices"]
        key=Tuple(record["vertex"])
        _mm_need(!haskey(vertices,key) && key in new.vertex_keys && !haskey(old.vertex_dsu.p,key),"twist inverse has duplicate or invalid new-vertex ancestry")
        vertices[key]=record["ancestry"]
    end
    _mm_need(all(v->haskey(old.vertex_dsu.p,v) || haskey(vertices,v),new.vertex_keys),"twist inverse omitted a new vertex")
    _mm_need(all(v->v in new.vertex_keys,keys(old.vertex_dsu.p)),"twist inverse removed a fixed old vertex")
    records=Dict{Tuple{String,Int},Any}();source_area=Dict((f,i-1)=>Q(0) for (f,i) in keys(old.triangles));target_area=copy(source_area)
    for r in spec["child_cells"]
        cell=_mm_cell(r["child"]);source=_mm_cell(r["source"]);target=_mm_cell(r["target"])
        _mm_need(!haskey(records,cell) && haskey(new.triangles,(cell[1],cell[2]+1)),"twist inverse has duplicate or missing child triangle")
        _mm_need(Tuple(r["vertices"])==new.triangles[(cell[1],cell[2]+1)],"twist inverse child vertices do not match actual incidence")
        _mm_need(haskey(source_area,source) && haskey(target_area,target) && target[1]==cell[1],"twist inverse triangle parent is invalid")
        tb=[_mt_vertex_bary(cell[1],v,target,old,vertices) for v in r["vertices"]]
        _mm_need(tb==_mt_live_bary.(r["target_weights"]),"twist inverse fixed-target coordinates disagree with vertex ancestry")
        sb=if haskey(h,target)
            _mm_need(haskey(h,source),"supported twist child has an unsupported source parent")
            [begin
                p=_mt_shear(_mt_point(w,h[target]),-spec["sign"])
                _mt_bary((p[1]-r["shift"],p[2]),h[source])
            end for w in tb]
        else
            _mm_need(source==target && r["shift"]==0,"twist inverse moved a cell outside its support")
            tb
        end
        _mm_need(sb==_mt_live_bary.(r["source_weights"]),"recorded inverse shear does not recover source child coordinates")
        ta,sa=_mt_ratio(tb),_mt_ratio(sb)
        _mm_need(ta>0 && sa>0,"inverse twist folds a common-refinement triangle")
        target_area[target]+=ta;source_area[source]+=sa;records[cell]=r
    end
    _mm_need(length(records)==length(new.triangles) && all(==(1),values(source_area)) && all(==(1),values(target_area)),"inverse twist children do not cover fixed source/target triangles")
    # The actual side paths must refine the OLD marked sides with their exact
    # old parameters, making the fixed-background map continuous across seams.
    for (f,p) in new.patches
        original=old.patches[f]
        _mm_need(p["side_order"]==original["side_order"],"twist inverse changed marked side order")
        for side in p["side_order"]
            path=original["sides"][side];ts=M.side_parameters(original,side)
            old_parameters=Dict(v=>ts[i] for (i,v) in enumerate(path));actual=M.side_parameters(p,side)
            _mm_need(issubset(Set(path),Set(p["sides"][side])),"twist inverse lost an old gate")
            for (i,v) in enumerate(p["sides"][side])
                value=if haskey(old_parameters,v)
                    old_parameters[v]
                else
                    record=vertices[(f,v)];_mm_need(record["kind"]=="edge","new gate is not on an old marked edge")
                    e=Tuple(record["edge"]);j=only(k for k in 1:length(path)-1 if M.edge(path[k],path[k+1])==e)
                    r=_mt_q(record["parameter"]);q=path[j]==e[1] ? r : 1-r
                    (1-q)*ts[j]+q*ts[j+1]
                end
                _mm_need(actual[i]==value,"twist inverse changed the fixed marking's gate parameters")
            end
        end
    end
    (;old,new,h,vertices,records)
end

function _mt_unrechart_features(current,spec)
    M=MarkedSurfaceDiagrams;inverse=get(spec,"rechart_inverse_specification",nothing)
    inverse===nothing && return deepcopy(current)
    # This helper reconstructs geometry from ACTUAL current native triangles,
    # then independently regenerates its maps. It never returns cached input.
    recovered=_marked_unrechart_geometry(current,inverse);corr=recovered.correspondence
    raw=Dict{String,Any}("version"=>1,"surface"=>deepcopy(recovered.geometry["surface"]),
        "patches"=>deepcopy(recovered.geometry["patches"]),"curves"=>deepcopy(current["curves"]),"crossings"=>deepcopy(current["crossings"]))
    edges=Dict(Tuple(r["source"])=>r["target"] for r in corr["edge_images"])
    vertices=Dict(Tuple(r["source"])=>r["target"] for r in corr["vertex_images"])
    children=Dict(_mm_cell(r["source"])=>Set(_mm_cell(c) for c in r["target"]) for r in corr["cell_children"])
    images=Dict(r["curve"]=>r for r in spec["curve_images"])
    _mm_need(Set(keys(images))==Set(c["id"] for c in current["curves"]),"identity rechart external curve set changed")
    for curve in raw["curves"]
        old=M.Dart[];expected=M.Dart[]
        for piece in images[curve["id"]]["edge_images"],d in piece["target"]
            _mm_need(haskey(edges,Tuple(d)),"actual reassembly omitted an external proof edge")
            push!(old,M.Dart(d...));append!(expected,[M.Dart(x...) for x in edges[Tuple(d)]])
        end
        _mm_need(M.parse_word(curve["word"])==expected,"actual curve is not the regenerated rechart edge image")
        curve["word"]=M.print_word(old) # old occurrence spelling was verified above.
    end
    locations=spec["proof_crossing_locations"]
    _mm_need(Set(keys(locations))==Set(x["id"] for x in current["crossings"]),"identity rechart crossing set changed")
    for x in raw["crossings"]
        point=locations[x["id"]]
        _mm_need(get(vertices,Tuple(point),nothing)==x["at"],"actual crossing vertex does not pull back through the regenerated identity map")
        x["at"]=deepcopy(point)
    end
    if haskey(current,"metadata");raw["metadata"]=deepcopy(current["metadata"]);end
    if haskey(current,"region_marks")
        raw["region_marks"]=deepcopy(current["region_marks"])
        for mark in raw["region_marks"]
            seed=spec["proof_region_seeds"][mark["id"]]
            _mm_need(_mm_cell(mark["seed"]) in children[_mm_cell(seed)],"actual marked seed is outside its regenerated proof-cell image")
            mark["seed"]=deepcopy(seed)
        end
    end
    infinity=_mm_infinity(current)
    if infinity!==nothing
        seed=spec["proof_infinity_triangle"]
        _mm_need(_mm_cell(infinity) in children[_mm_cell(seed)],"actual infinity is outside its regenerated proof-cell image")
        raw["metadata"]["display"]["infinity_triangle"]=deepcopy(seed)
    end
    raw # Internal construction complex only; native source is returned later.
end

"""Invert a generated whole-diagram twist using the SAME recorded h and -s.

The inverse shear is applied to actual common-refinement vertices/triangles.
Its positive source-child maps and source-edge curve partitions are checked
before collapsing to the old mesh. Rebuilding h occurs only as an independent
verification of the recorded support, never as a replacement inverse action.
"""
function marked_inverse_dehn_twist(input::AbstractDict,specification::AbstractDict)
    current=_marked_document(input);spec=_plain_value(specification);M=MarkedSurfaceDiagrams
    _mm_need(_mm_hash(current)==spec["target_hash"],"generated twist target has changed")
    _mm_need(spec["sign"] isa Integer && !(spec["sign"] isa Bool) && spec["sign"] in (-1,1),"invalid recorded twist sign")
    proof=_mt_unrechart_features(current,spec)
    geometry=_mt_inverse_geometry(proof,spec);old,new=geometry.old,geometry.new
    candidate=deepcopy(proof);candidate["surface"]=deepcopy(spec["fine_source_mesh"]["surface"]);candidate["patches"]=deepcopy(spec["fine_source_mesh"]["patches"])
    records_by_curve=Dict(r["curve"]=>r for r in spec["curve_images"])
    _mm_need(Set(keys(records_by_curve))==Set(c["id"] for c in proof["curves"]),"inverse twist curve correspondence differs")
    for curve in candidate["curves"]
        record=records_by_curve[curve["id"]];actual=M.parse_word(curve["word"]);expected=M.Dart[];source=M.Dart[]
        for piece in record["edge_images"]
            f,u,v=piece["source"];push!(source,M.Dart(f,u,v));previous=MarkedSurfaceDiagrams.Rat(0)
            images=[M.Dart(x...) for x in piece["target"]]
            _mm_need(!isempty(images),"inverse twist omitted a source curve edge")
            for d in images
                a,b=_mt_pullback_interval(d,M.Dart(f,u,v),new,geometry.records,old)
                _mm_need(a==previous && a<b,"inverse curve edge partition is incomplete or reverses orientation")
                previous=b
            end
            _mm_need(previous==1,"inverse curve did not traverse the entire original edge")
            append!(expected,images)
        end
        _mm_need(actual==expected,"current twist curve is not its recorded geometric edge image")
        _mm_need(M.parse_word(record["source_word_format"])==source,"inverse lexical word differs from geometrically reconstructed edges")
        curve["word"]=record["source_word_format"]
    end
    for mark in get(candidate,"region_marks",[]);mark["seed"]=deepcopy(spec["region_seeds"][mark["id"]]);end
    _mm_infinity(proof)!==nothing && (candidate["metadata"]["display"]["infinity_triangle"]=deepcopy(spec["infinity_triangle"]))
    _marked_states(candidate)
    _mm_need(_mm_hash(candidate)==spec["fine_source_hash"],"inverse shear/common-refinement collapse does not recover prepared source")
    history=spec["construction"];_mm_need(!isempty(history),"twist inverse has no preparation history")
    coarse=marked_collapse_subdivision(candidate,last(history)["inverse_specification"]).document
    product=_marked_annulus_product(_marked_states(coarse).geometry,_marked_states(candidate).geometry,spec["core"],Set(_mm_cell(c) for c in spec["support_cells"]))
    product_rows=[Dict("cell"=>deepcopy(t.id),"points"=>_mt_json_point.(t.points)) for t in product.triangles]
    _mm_need(_canonical_json(product_rows)==_canonical_json(spec["product_triangles"]) && _canonical_json(product.certificate)==_canonical_json(spec["product_certificate"]),"recorded annular product/support differs from actual source construction")
    original=coarse
    for w in reverse(history[1:end-1]);original=replay_marked_move(original,inverse_marked_move(w)).document;end
    _mm_need(_mm_hash(original)==spec["source_hash"],"inverse preparation history does not restore source")
    regenerated=marked_dehn_twist(original,spec["forward_specification"])
    _mm_need(regenerated.document==current && _canonical_json(regenerated.witness["inverse_specification"])==_canonical_json(spec),"independent forward twist does not regenerate target and recorded inverse")
    # Source: the forward common-refinement/rechart certificate above. Raw
    # proof cells and current native cells use different charts; retain their
    # contexts explicitly instead of presenting raw IDs as native incidence.
    witness=_mm_witness("generated_twist_inverse",current,original,spec,spec["forward_specification"],
        Dict("inverse_map"=>"same recorded annular h, shear with opposite sign, then verified common-refinement collapse",
            "raw_child_cells"=>spec["child_cells"],"raw_cell_scope"=>spec["proof_complex_scope"],
            "native_child_cells"=>spec["native_child_cells"],
            "native_cell_scope"=>"forward target native triangles; source/target parents remain fixed pre-twist prepared cells",
            "rechart_inverse_specification"=>spec["rechart_inverse_specification"],"inverse_construction"=>history))
    witness["inverse_kind"]="dehn_twist"
    (;document=original,witness)
end

"""Read images in the fixed PRE-TWIST PREPARED fine-mesh presentation basis.

Source: marked_groups.jl per-oriented-cell edge words. Canonical anchors on an
old edge use a chosen edge occurrence endpoint, not merely minimum quotient
vertex (which is ambiguous on loop edges). Interior anchors are old triangle
corners. Consecutive anchor legs cancel; triangle choices agree modulo the
prepared surface relator. Preparation's explicit identity transports link this
basis to the input marking; its generator numbers need not match input numbers.
Only abelianized equality is asserted here.
"""
function _mt_old_basis_images(fine_doc,result,core,sign)
    M=MarkedSurfaceDiagrams;old=_marked_states(fine_doc).geometry;new=_marked_states(result.document).geometry
    data=marked_curve_words(fine_doc);pairing=M.fine_edge_occurrences(old)
    records=Dict(_mm_cell(r["child"])=>r for r in result.records)
    function anchor(v,record)
        j=only(findall(==(v),record["vertices"]));weights=_mt_live_bary(record["target_weights"][j])
        cell=_mm_cell(record["target"]);f=cell[1];t=old.triangles[(f,cell[2]+1)]
        original=findfirst(==(1),weights);original!==nothing && return t[original]
        z=findfirst(==(0),weights)
        z===nothing && return t[1]
        ends=[k for k in 1:3 if k!=z];edge=M.edge(t[ends[1]],t[ends[2]])
        a,b=only((a,b) for (a,b) in M.tedges(t) if M.edge(a,b)==edge)
        occurrence=(f,cell[2]+1,a,b);partner=pairing[occurrence]
        isless(occurrence,partner) ? a : b
    end
    dimension=length(data.presentation.generators)
    exponents(word)=begin
        out=zeros(Int,dimension)
        for (i,s) in word;out[i]+=s;end
        out
    end
    reports=Any[];core_class=exponents(data.curve_words[core])
    for c in result.document["curves"]
        word=Int[];retractions=Any[]
        for d in M.parse_word(c["word"])
            i,_,_=first(new.incidence[d.face][M.edge(d.u,d.v)])
            record=records[(d.face,i-1)];cell=_mm_cell(record["target"])
            u=anchor(d.u,record);v=anchor(d.v,record)
            edge_word=u==v ? Tuple{Int,Int}[] : data.edge_words[(cell[1],cell[2],u,v)]
            append!(word,[i*s for (i,s) in edge_word])
            push!(retractions,Dict("image_dart"=>[d.face,d.u,d.v],"old_triangle"=>collect(cell),"anchors"=>[u,v],"old_edge_word"=>[[i,s] for (i,s) in edge_word]))
        end
        reduced=_marked_reduce(word);image=[(abs(x),Base.sign(x)) for x in reduced]
        intersection=0
        if c["id"]!=core
            for x in fine_doc["crossings"]
                Set(x["curves"])==Set([core,c["id"]]) || continue
                intersection+=x["curves"][1]==core ? x["sign"] : -x["sign"]
            end
        end
        expected=exponents(data.curve_words[c["id"]])+sign*intersection*core_class
        actual=exponents(image)
        _mm_need(actual==expected,"constructed twist image has the wrong homology action in the fixed old basis")
        push!(reports,Dict("curve"=>c["id"],"source_word"=>[[i,s] for (i,s) in data.curve_words[c["id"]]],
            "image_word_old_basis"=>[[i,s] for (i,s) in image],"intersection_with_core"=>intersection,
            "source_exponents"=>exponents(data.curve_words[c["id"]]),"image_exponents"=>actual,"expected_exponents"=>expected,
            "edge_retractions"=>retractions,"word_semantics"=>"fixed pre-twist prepared fine-mesh surface presentation; modulo its relator; exponent equality checked"))
    end
    reports
end

"""Apply a signed whole-diagram twist using one recorded annular product.

Positive convention: if the oriented core a crosses b positively, [b] maps
to [b]+[a] in the fixed pre-twist prepared basis. The fixed background is the same prepared marked mesh,
with a common refinement; every curve, crossing, D and infinity moves by F.
"""
function marked_dehn_twist(input::AbstractDict,specification::AbstractDict)
    # Source: marked_stabilization.jl's narrow ID preflight. Native topology
    # allows scalar IDs, but JSON inverse lookup objects require string keys.
    # Reject unsupported IDs rather than silently coercing mathematical labels.
    raw=_plain_value(input)
    for field in ("curves","crossings","region_marks")
        rows=get(raw,field,[])
        rows isa AbstractVector || rows isa Tuple || continue
        for row in rows
            row isa AbstractDict && haskey(row,"id") || continue
            _mm_need(row["id"] isa AbstractString && !isempty(row["id"]),"twist source $field IDs must be nonempty strings")
        end
    end
    spec=_plain_value(specification)
    _mm_need(issubset(Set(keys(spec)),Set(["core","sign"])),"twist accepts only core/sign: selected-family or beta-only action is a different operation")
    _mm_need(haskey(spec,"core") && spec["core"] isa AbstractString && !isempty(spec["core"]),"twist core must be a nonempty string ID")
    original=_marked_document(raw)
    core=String(spec["core"]);sign=get(spec,"sign",1)
    _mm_need(sign isa Integer && !(sign isa Bool) && sign in (-1,1),"twist sign must be +1 or -1")
    prepared=_mt_prepare(original,core);overlay=_mt_periodic_overlay(prepared.product,sign)
    result=_mt_overlay_document(prepared.fine,prepared.product,overlay,sign)
    basis=_mt_old_basis_images(prepared.fine,result,core,sign)
    h_records=[Dict("cell"=>deepcopy(t.id),"points"=>_mt_json_point.(t.points)) for t in prepared.product.triangles]
    inverse=Dict{String,Any}("forward_specification"=>spec,"core"=>core,"sign"=>sign,"construction"=>prepared.history,
        "fine_source_mesh"=>Dict("surface"=>deepcopy(prepared.fine["surface"]),"patches"=>deepcopy(prepared.fine["patches"])),
        "fine_source_hash"=>_mm_hash(prepared.fine),"source_hash"=>_mm_hash(original),"target_hash"=>_mm_hash(result.document),
        "product_triangles"=>h_records,"product_certificate"=>prepared.product.certificate,
        "support_cells"=>[collect(c) for c in sort(collect(prepared.support))],
        "child_cells"=>result.raw.records,"new_vertices"=>result.raw.vertex_records,"curve_images"=>result.raw.image_records,
        "proof_complex_scope"=>"internal pre-rechart fixed-target common refinement; not an accepted native diagram",
        "native_child_cells"=>result.records,"native_curve_images"=>result.image_records,
        "rechart_inverse_specification"=>result.rechart===nothing ? nothing : result.rechart.inverse_specification,
        "proof_crossing_locations"=>Dict(x["id"]=>deepcopy(x["at"]) for x in result.raw.complex["crossings"]),
        "proof_region_seeds"=>_mm_inverse_seeds(result.raw.complex),
        "region_seeds"=>_mm_inverse_seeds(prepared.fine))
    _mm_infinity(prepared.fine)!==nothing && (inverse["infinity_triangle"]=deepcopy(_mm_infinity(prepared.fine)))
    _mm_infinity(result.raw.complex)!==nothing && (inverse["proof_infinity_triangle"]=deepcopy(_mm_infinity(result.raw.complex)))
    correspondence=Dict("construction"=>prepared.history,"background_marking"=>"fixed prepared target surface via common refinement and explicit marked identity rechart; public chart names may differ",
        "annular_product"=>prepared.product.certificate,"child_cells"=>result.records,
        "raw_proof_vertices"=>result.raw.vertex_records,
        "rechart"=>result.rechart===nothing ? nothing : result.rechart.correspondence,
        "vertex_scope"=>"all actual native vertex occurrences are described by child_cells; raw_proof_vertices are pre-rechart context only",
        "curve_images"=>result.image_records,"regions"=>result.region_map,"geometric_regions"=>result.geometry_regions,
        "old_basis_action"=>basis,"basis_scope"=>"fixed pre-twist prepared fine-mesh generators; linked to input marking by construction, not original input generator numbers",
        "basis_source_hash"=>MarkedSurfaceDiagrams.representation_digest(prepared.fine),
        "basis_geometry_hash"=>MarkedSurfaceDiagrams.representation_digest(_mm_bare(prepared.fine)),
        "basis_marking_hash"=>MarkedSurfaceDiagrams.marking_digest(prepared.fine),
        "crossings"=>[Dict("source"=>x["id"],"target"=>x["id"],"at"=>x["at"],"sign"=>x["sign"]) for x in result.document["crossings"]],
        "action_scope"=>"whole alpha/beta/aux diagram and selected regions; not beta-only action")
    witness=_mm_witness("dehn_twist",original,result.document,spec,inverse,correspondence);witness["inverse_kind"]="generated_twist_inverse"
    (;document=result.document,witness)
end
