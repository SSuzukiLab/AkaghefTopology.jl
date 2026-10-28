# Exact full-cell maps used by the marked-display reconstruction recipe.
# Source: marked_polygon_overlay.jl (Polyhedra), MarkedSurfaceDiagrams exact
# orientation/segment predicates, and the retained source/collar triangulations.
# A saved diagnostic receipt is never a substitute for executing these checks.
# Heegaard.jl includes this helper for the common native and saved-view path.

const _MDSMQ = Rational{BigInt}
const _MDSMP = Tuple{_MDSMQ,_MDSMQ}
_mdsm_q(x::Integer)=_MDSMQ(x)
_mdsm_q(x::Rational)=BigInt(numerator(x))//BigInt(denominator(x))
# Source: Julia 1.12 base/mpfr.jl::Rational{BigInt}(::AbstractFloat) uses
# precision/exponent/ldexp to retain the exact binary value. Continued-fraction
# rationalize(tol=0) can differ even at zero tolerance for tiny coordinates.
_mdsm_q(x::AbstractFloat)=isfinite(x) ? _MDSMQ(x) : throw(ArgumentError("Nonfinite source-map coordinate"))
function _mdsm_q(x::AbstractString)
    parts=split(replace(x,"//"=>"/"),'/')
    length(parts)==1 && return _MDSMQ(parse(BigInt,parts[1]))
    length(parts)==2 || throw(ArgumentError("Invalid rational source-map coordinate"))
    parse(BigInt,parts[1])//parse(BigInt,parts[2])
end

# Source: MarkedSurfaceDiagrams.orient/segments_meet and the directed boundary
# cycles validated by _mdsm_mesh. Winding is used only after simplicity; a
# boundary contact is rejected instead of assigned an arbitrary inside value.
function _mdsm_strictly_inside(p,polygon)
    M=MarkedSurfaceDiagrams;winding=0
    for (a,b) in zip(polygon,circshift(polygon,-1))
        turn=M.orient(a,b,p)
        iszero(turn) && min(a[1],b[1])<=p[1]<=max(a[1],b[1]) &&
            min(a[2],b[2])<=p[2]<=max(a[2],b[2]) && return false
        a[2]<=p[2]<b[2] && turn>0 && (winding+=1)
        b[2]<=p[2]<a[2] && turn<0 && (winding-=1)
    end
    abs(winding)==1
end

"""Bind disjoint collar maps to both complete ordered hole boundaries.

Source: _mdsm_collar_layer's full annuli. Each outer support is fixed, simple,
contained in the declared exterior, and disjoint from every other support.
Thus the annular maps extend by identity and replace exactly these holes.
"""
function _mdsm_bind_collars(layers,source_holes,target_holes,exterior)
    word(q)=_mdsm_boundary_word(Dict(enumerate(q)),collect(eachindex(q)))
    _mdsm_need(Set(word.(source_holes))==Set(word([a.source[v] for v in a.inner]) for a in layers),"Collar source holes differ from full domain")
    _mdsm_need(Set(word.(target_holes))==Set(word([a.target[v] for v in a.inner]) for a in layers),"Collar target holes differ from full domain")
    supports=[[a.source[v] for v in a.outer] for a in layers]
    meet(a,b)=any(MarkedSurfaceDiagrams.segments_meet(u,v,x,y) for
        (u,v) in zip(a,circshift(a,-1)) for (x,y) in zip(b,circshift(b,-1)))
    for (i,p) in enumerate(supports)
        _mdsm_need(!meet(p,exterior) && _mdsm_strictly_inside(first(p),exterior),"Collar support leaves its fixed exterior chart")
        for q in supports[i+1:end]
            _mdsm_need(!meet(p,q) && !_mdsm_strictly_inside(first(p),q) && !_mdsm_strictly_inside(first(q),p),"Collar supports overlap or nest")
        end
    end
    length(layers)
end

"""Check the actual infinity-cap fan and both exterior attachment joins.

Source: build_infinity_extension's reciprocal-radius construction and the
already validated fixed PL cells. Positive native fan wedges, reversed planar
fan orientation, and exact exterior homothety give the Jacobian factor / t^3.
The t=1 join uses the original frame; t=1/2 uses the actual fixed PL boundary.
The native barycentric center is the distinguished one-point compactification.
"""
function _mdsm_verify_infinity(recipe,native,fixed_source,fixed_target)
    M=MarkedSurfaceDiagrams;chart=recipe["source_chart"];inf=chart["infinity_chart"]
    extension=recipe["infinity_extension"];P=_mdsm_points(chart["positions"])
    S=fixed_source.points;R=fixed_target.points;frame=Int.(chart["chart_outer_boundary"])
    ext=Int.(recipe["fixed_map"]["source"]["outer"])
    _mdsm_need(extension["source_triangle"]==inf["source_triangle"]==recipe["display_request"]["infinity_triangle"],"Infinity source triangle changed")
    _mdsm_need(extension["infinity_point_native_barycentric"]==inf["point_barycentric"] &&
        extension["native_frame_barycentric"]==inf["frame_barycentric"],"Native infinity fan changed")
    _mdsm_need(extension["source_frame_vertices"]==frame && extension["source_exterior_vertices"]==ext,"Infinity frame correspondence changed")
    rows=[_mdsm_q.(q) for q in vcat([inf["point_barycentric"]],inf["frame_barycentric"])]
    _mdsm_need(all(length(q)==3 && sum(q)==1 && minimum(q)>0 for q in rows),"Invalid native infinity barycentric point")
    U=(rows[1][2],rows[1][3]);B=[(q[2],q[3]) for q in rows[2:end]]
    _mdsm_need(length(B)==length(frame)==length(ext) && length(frame)>=3,"Infinity wedge count mismatch")
    cap_points=Dict{Int,_MDSMP}(0=>U)
    for (i,p) in enumerate(B);cap_points[i]=p;end
    _mdsm_mesh(cap_points,[[0,i,mod1(i+1,length(B))] for i in eachindex(B)];expected_boundary=[collect(eachindex(B))])
    ordered=_mpo_area([P[v] for v in frame])>0 ? frame : reverse(frame)
    _mdsm_need(_mdsm_boundary_word(P,ordered)==_mdsm_boundary_word(native.points,native.outer),"Native infinity frame differs from physical source boundary")
    _mdsm_need(_mdsm_boundary_word(S,ext)==_mdsm_boundary_word(S,fixed_source.outer),"Exterior frame differs from fixed source boundary")
    origin=_mdsm_point(extension["source_plane_origin"]);target_origin=_mdsm_point(extension["target_plane_origin"])
    _mdsm_need(all(S[v]==P[v] for v in frame),"Fixed source moved the native frame")
    direct=isempty(chart["boundaries"])
    _mdsm_need(get(extension,"attachment_kind","homothetic_collar")== (direct ? "direct_no_cut" : "homothetic_collar"),"Infinity attachment kind differs from the native cut surface")
    if direct
        _mdsm_need(chart["genus"]==0 && isempty(chart["cuts"]) && ordered==ext,"No-cut infinity frame correspondence changed")
        _mdsm_same_domain(native,fixed_source)
        _mdsm_need(all(S[v]==P[v]==R[v] for v in frame),"Direct infinity join is not pointwise identity")
        _mdsm_need(origin==target_origin==(_mdsm_q(0),_mdsm_q(0)),"Direct infinity origin differs from the verified exterior family")
    else
        _mdsm_need(all(ntuple(k->2*P[v][k]-S[w][k],2)==origin for (v,w) in zip(ordered,ext)),"Exterior source homothety or wedge correspondence changed")
    end
    _mdsm_need(extension["native_frame_exterior_vertices"]==[ext[only(findall(==(v),ordered))] for v in frame],"Infinity wedge-to-exterior correspondence changed")
    _mdsm_need(_mdsm_point.(extension["source_exterior_positions"])==[S[v] for v in ext] &&
        _mdsm_point.(extension["target_exterior_positions"])==[R[v] for v in ext],"Infinity attachment uses unrelated coordinates")
    factors=_MDSMQ[]
    for j in eachindex(frame)
        k=mod1(j+1,length(frame));a=M.orient(U,B[j],B[k]);b=M.orient(origin,P[frame[j]],P[frame[k]])
        _mdsm_need(a>0 && b<0,"Infinity reciprocal fan loses orientation")
        push!(factors,-b/a)
    end
    _mdsm_need(factors==_mdsm_q.(extension["positive_jacobian_factors"]),"Infinity Jacobian factors changed")
    _mdsm_need(all(M.orient(origin,S[u],S[v])>0 && M.orient(target_origin,R[u],R[v])>0 for (u,v) in zip(ext,circshift(ext,-1))),"Infinity exterior fans are not strictly star shaped")
    family=recipe["fixed_quotient"]["exterior_extension"]
    start=_mdsm_point.(first(family["outer_stages"]))
    _mdsm_need(_mdsm_boundary_word(Dict(enumerate(start)),collect(eachindex(start)))==_mdsm_boundary_word(R,ext),"Infinity extension does not join the fixed-quotient family")
    _mdsm_need(_mdsm_point(family["origin"])==target_origin,"Infinity target ray origins differ")
    _mdsm_need(extension["finite_frame_join_parameter"]=="1" && extension["outer_collar_join_parameter"]==(direct ? "1" : "1/2"),"Infinity radial join parameters changed")
    length(factors)
end

"""Validate every full-cell component supplied by the common Python consumer.

This is the exact-cell half of the SAME-process pipeline, not a standalone
acceptance of arbitrary saved checks. The Python consumer freshly binds native
partitions/marked chains and verifies paired-lambda, boundary and Bezier steps;
this function independently validates complete domains and exact collar maps.
The compact recipe retains all cells and deterministic point-location maps.
"""
function _mdsm_verify_recipe_cells(recipe;progress=identity)
    _mdsm_need(recipe["schema"]=="marked-display-source-map/1","Unknown source-map recipe")
    chart=recipe["source_chart"];P=_mdsm_points(chart["positions"])
    old=_mdsm_mesh(P,chart["reduction_source_triangles"])
    reduced=_mdsm_mesh(P,chart["triangles"])
    _mdsm_same_domain(old,reduced)
    progress("source/reduced domains")
    physical=[t for (t,label) in zip(chart["reduction_source_triangles"],chart["reduction_source_triangle_labels"]) if label!==nothing]
    native=_mdsm_mesh(P,physical)
    source=recipe["fixed_map"]["source"];S=_mdsm_points(source["source_positions"])
    Q=_mdsm_points(recipe["fixed_map"]["target_positions"])
    fixed_source=_mdsm_mesh(S,source["triangles"])
    fixed_target=_mdsm_mesh(Q,source["triangles"];expected_boundary=fixed_source.cycles)
    # Source: fixed_mesh's unchanged final target acceptance threshold. Exact
    # intermediate refinements may be smaller; the actual fixed target may not.
    _mdsm_need(all(MarkedSurfaceDiagrams.orient((Q[v] for v in t)...)>_mdsm_q(1e-8) for t in fixed_target.triangles),"Fixed target misses the existing mesh-area threshold")
    progress("fixed source/target meshes")
    infinity_wedges=_mdsm_verify_infinity(recipe,native,fixed_source,fixed_target)
    progress("native infinity cap, exterior collar and quotient-family joins")
    layers=recipe["collars"]
    _mdsm_need(length(layers)==2 && layers[1]["layer"]=="native_to_clipped_source" &&
        layers[2]["layer"]=="final_mesh_polygon_to_circle_polygon","Source-map collar order changed")
    expected=Set((b["cut"],b["side"]) for b in chart["boundaries"])
    for layer in layers
        keys=[(b["cut"],b["side"]) for b in layer["boundaries"]]
        _mdsm_need(length(keys)==length(expected) && Set(keys)==expected,"Collar marked side binding changed")
    end
    native_layers=_mdsm_collar_layer(layers[1];progress=progress)
    circle_layers=_mdsm_collar_layer(layers[2];progress=progress)
    holes(m)=[[m.points[v] for v in c] for c in m.cycles if c!=m.outer]
    outside(m)=[m.points[v] for v in m.outer]
    _mdsm_bind_collars(native_layers,holes(native),holes(fixed_source),outside(native))
    # fixed_mesh adds the explicitly retained exterior collar before moving the
    # grid. Its source frame contains the old native frame with no crossings.
    outer=outside(fixed_source);inner=outside(native)
    if isempty(chart["boundaries"])
        _mdsm_same_domain(native,fixed_source)
    else
        _mdsm_need(all(_mdsm_strictly_inside(p,outer) for p in inner),"Native infinity frame leaves the retained exterior collar")
        _mdsm_need(!any(MarkedSurfaceDiagrams.segments_meet(a,b,c,d) for
            (a,b) in zip(inner,circshift(inner,-1)) for (c,d) in zip(outer,circshift(outer,-1))),"Infinity collar boundaries cross")
    end
    domain=recipe["port_domain"];D=_mdsm_points(domain["source_positions"])
    port_domain=_mdsm_mesh(D,domain["triangles"])
    _mdsm_need(_mdsm_boundary_word(Q,fixed_target.outer)==_mdsm_boundary_word(D,port_domain.outer),"Port map changes the exterior frame")
    _mdsm_bind_collars(circle_layers,holes(fixed_target),holes(port_domain),outside(fixed_target))
    map=recipe["port_map"];U=_mdsm_points(map["source_positions"]);V=_mdsm_points(map["target_positions"])
    cdt_source=_mdsm_mesh(U,map["triangles"]);cdt_target=_mdsm_mesh(V,map["triangles"];expected_boundary=cdt_source.cycles)
    _mdsm_same_domain(port_domain,cdt_source)
    bound=_mdsm_linear_face_bound(U,V,map["triangles"])
    progress("complete port domains and whole-step determinants")
    Dict("full_cell_maps_verified"=>true,"native_faces"=>length(native.triangles),
        "fixed_faces"=>length(fixed_source.triangles),"port_faces"=>length(cdt_source.triangles),
        "specified_infinity_wedges"=>infinity_wedges,
        "collar_sectors"=>sum((length(b["sectors"]) for layer in layers for b in layer["boundaries"]);init=0),
        "collar_common_cells"=>sum((a.common_cells for a in vcat(native_layers,circle_layers));init=0),
        "port_linear_minimum"=>string(bound),"expanded_global_overlay_required"=>false)
end
_mdsm_point(x)=(_mdsm_q(x[1]),_mdsm_q(x[2]))
_mdsm_need(ok,message)=ok ? nothing : throw(ArgumentError(message))
_mdsm_points(raw)=Dict(parse(Int,string(v))=>_mdsm_point(p) for (v,p) in raw)

"""Check a complete oriented planar manifold mesh and its exact boundary.

Positive faces, connected manifold links and a simple correctly oriented full
boundary give degree one on the bounded planar domain. This certifies global
coverage and injectivity, rather than relying only on positive local Jacobians.
The complete cells stay available for exact composition and reconstruction.
"""
function _mdsm_mesh(points,triangles; expected_boundary=nothing)
    M=MarkedSurfaceDiagrams
    T=[Int.(t) for t in triangles]
    edges=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    links=Dict{Int,Dict{Int,Vector{Int}}}()
    for t in T
        _mdsm_need(length(t)==3 && length(Set(t))==3,"Source-map face is not a triangle")
        _mdsm_need(M.orient((points[v] for v in t)...)>0,"Source-map face is not strictly positive")
        for j in 1:3
            u=t[j];v=t[mod1(j+1,3)];w=t[mod1(j+2,3)]
            push!(get!(edges,minmax(u,v),Tuple{Int,Int}[]),(u,v))
            link=get!(links,u,Dict{Int,Vector{Int}}())
            push!(get!(link,v,Int[]),w);push!(get!(link,w,Int[]),v)
        end
    end
    _mdsm_need(!isempty(T),"Empty source-map mesh")
    _mdsm_need(all(length(e)==1 || (length(e)==2 && e[1]==reverse(e[2])) for e in values(edges)),"Source-map edge incidence is inconsistent")
    boundary=[only(e) for e in values(edges) if length(e)==1]
    next=Dict(u=>v for (u,v) in boundary)
    _mdsm_need(length(next)==length(boundary) && Set(keys(next))==Set(values(next)),"Source-map boundary branches")
    boundary_nodes=Set(keys(next));cycles=Vector{Int}[];unseen=copy(boundary_nodes)
    while !isempty(unseen)
        firstnode=minimum(unseen);cycle=[firstnode];v=next[firstnode]
        while v!=firstnode
            _mdsm_need(v in unseen && !(v in cycle),"Source-map boundary is not a disjoint cycle")
            push!(cycle,v);v=next[v]
        end
        setdiff!(unseen,cycle);push!(cycles,cycle)
    end
    _mdsm_need(length(links)-length(edges)+length(T)==2-length(cycles),"Source-map Euler characteristic disagrees with its boundary")
    for (v,link) in links
        degrees=length.(values(link))
        _mdsm_need(all(x in (1,2) for x in degrees) && count(==(1),degrees)==(v in boundary_nodes ? 2 : 0),"Source-map vertex link is not a manifold")
        seen=Set([first(keys(link))]);todo=collect(seen)
        while !isempty(todo)
            for w in link[pop!(todo)]
                w in seen && continue
                push!(seen,w);push!(todo,w)
            end
        end
        _mdsm_need(seen==Set(keys(link)),"Source-map vertex link is disconnected")
    end
    seen=Set([first(keys(links))]);todo=collect(seen)
    while !isempty(todo)
        for w in keys(links[pop!(todo)])
            w in seen && continue
            push!(seen,w);push!(todo,w)
        end
    end
    _mdsm_need(seen==Set(keys(links)),"Source-map mesh is disconnected")
    areas=[_mpo_area([points[v] for v in c]) for c in cycles]
    _mdsm_need(count(>(0),areas)==1 && all(x->!iszero(x),areas),"Source-map outer/hole orientations disagree")
    for i in eachindex(boundary), j in i+1:length(boundary)
        a,b=boundary[i];c,d=boundary[j]
        if !isempty(intersect(Set((a,b)),Set((c,d))))
            # Adjacent edges must not overlap along a ray at their common node.
            common=only(intersect(Set((a,b)),Set((c,d))))
            u=a==common ? b : a;v=c==common ? d : c
            q=points[common];x=points[u];y=points[v]
            _mdsm_need(M.orient(q,x,y)!=0 || (x[1]-q[1])*(y[1]-q[1])+(x[2]-q[2])*(y[2]-q[2])<0,"Source-map adjacent boundary edges overlap")
        else
            _mdsm_need(!M.segments_meet(points[a],points[b],points[c],points[d]),"Source-map boundary intersects")
        end
    end
    if expected_boundary!==nothing
        _mdsm_need(Set(boundary)==Set((Int(u),Int(v)) for cycle in expected_boundary for (u,v) in zip(cycle,circshift(cycle,-1))),"Source-map boundary differs from retained source")
    end
    (;points,triangles=T,cycles,edges,outer=cycles[only(findall(>(0),areas))])
end

"""Exact affine point image, with deterministic and consistent boundary ties.

Both coordinate meshes must first pass `_mdsm_mesh` with the same cells. The
lowest face index is the canonical tie; every other containing face must give
the identical image. No numerical inverse or rounded barycentric guess enters.
"""
function _mdsm_map_point(point,source,target,triangles)
    M=MarkedSurfaceDiagrams;p=_mdsm_point(point);image=nothing;chosen=nothing
    for (index,t) in enumerate(triangles)
        a,b,c=(source[v] for v in t);den=M.orient(a,b,c)
        _mdsm_need(den>0,"Unvalidated source-map triangle")
        weights=(M.orient(p,b,c)/den,M.orient(a,p,c)/den,M.orient(a,b,p)/den)
        all(>=(0),weights) || continue
        q=ntuple(k->sum(weights[j]*target[t[j]][k] for j in 1:3),2)
        _mdsm_need(image===nothing || image==q,"Source-map boundary tie gives different images")
        if image===nothing;image=q;chosen=index;end
    end
    _mdsm_need(image!==nothing,"Point is outside the retained source-map domain")
    (;point=image,triangle=chosen)
end

"""Canonical directed boundary word, omitting only forward collinear nodes.

Source: the complete-domain composition check in marked-display reconstruction.
This compares literal exact domains when two triangulations have different
boundary subdivisions; it never replaces a marked edge or its parameter map.
"""
function _mdsm_boundary_word(points,cycle)
    M=MarkedSurfaceDiagrams;q=[points[v] for v in cycle]
    _mdsm_need(length(Set(q))==length(q),"Repeated source-map boundary point")
    changed=true
    while changed
        changed=false;out=_MDSMP[]
        for i in eachindex(q)
            u=q[mod1(i-1,length(q))];v=q[i];w=q[mod1(i+1,length(q))]
            if iszero(M.orient(u,v,w))
                _mdsm_need((v[1]-u[1])*(w[1]-v[1])+(v[2]-u[2])*(w[2]-v[2])>0,"Boundary subdivision reverses direction")
                changed=true
            else
                push!(out,v)
            end
        end
        q=out;_mdsm_need(length(q)>=3,"Collapsed source-map boundary")
    end
    firstnode=argmin(q)
    Tuple(vcat(q[firstnode:end],q[1:firstnode-1]))
end

"""Bind two complete validated meshes to exactly the same planar domain.

The maps are composed as functions using exact point location. No expanded
intersection table is required: degree-one mesh validation, identical complete
directed boundary words and consistent affine edge ties provide full coverage.
Marked curve/D words are an additional binding and are not inferred here.
"""
function _mdsm_same_domain(first_mesh,second_mesh)
    a=[_mdsm_boundary_word(first_mesh.points,c) for c in first_mesh.cycles]
    b=[_mdsm_boundary_word(second_mesh.points,c) for c in second_mesh.cycles]
    _mdsm_need(length(Set(a))==length(a)==length(b) && Set(a)==Set(b),"Complete source-map domains differ")
    length(a)
end

"""Minimum of an exact determinant quadratic over the entire unit interval.

Source: linear interpolation of vertex coordinates, as in the existing
check_fairing_isotopy polynomial checks. Endpoint positivity alone is not used.
This is a local orientation check; all-time embedded boundary/paired gluing
conditions remain separate requirements for a whole-surface isotopy.
"""
function _mdsm_quadratic_minimum(c,b,a)
    result=min(c,c+b+a)
    if a>0 && 0 < -b < 2a
        result=min(result,c-b*b/(4a))
    end
    result
end
function _mdsm_linear_face_bound(source,target,triangles)
    cross(u,v)=u[1]*v[2]-u[2]*v[1]
    sub(u,v)=(u[1]-v[1],u[2]-v[2])
    bound=nothing
    for t in triangles
        a,b,c=(source[v] for v in t);x,y,z=(target[v] for v in t)
        e=sub(b,a);f=sub(c,a);de=sub(sub(y,x),e);df=sub(sub(z,x),f)
        value=_mdsm_quadratic_minimum(cross(e,f),cross(de,f)+cross(e,df),cross(de,df))
        _mdsm_need(value>0,"Source-map linear step loses orientation")
        bound=bound===nothing ? value : min(bound,value)
    end
    _mdsm_need(bound!==nothing,"Empty source-map linear step")
    bound
end

"""Reconstruct one collar sector with the existing exact Polyhedra overlay.

Source/target have the same ordered boundary points, possibly different disk
triangulations. Both triangulations are pulled through the same rational convex
polygon; every common cell is transported by its exact affine barycentric map.
"""
function _mdsm_collar_sector(sector)
    M=MarkedSurfaceDiagrams
    A=_mdsm_point.(sector["source"]);B=_mdsm_point.(sector["target"]);n=length(A)
    _mdsm_need(n==length(B),"Collar boundary correspondence length mismatch")
    TA=[Int.(t).+1 for t in sector["source_triangles"]];TB=[Int.(t).+1 for t in sector["target_triangles"]]
    boundary=[collect(1:n)]
    C=M.circle_boundary(n)
    for (points,T) in ((A,TA),(B,TB),(C,TA),(C,TB))
        _mdsm_mesh(Dict(enumerate(points)),T;expected_boundary=boundary)
    end
    source=[(;id=i,points=C[t]) for (i,t) in enumerate(TA)]
    target=[(;id=i,points=C[t]) for (i,t) in enumerate(TB)]
    overlay=_marked_triangle_overlay(source,target)
    _mdsm_need(all(c.area_matches for c in overlay.coverage),"Collar common refinement does not cover its source")
    maps=Any[]
    for cell in overlay.cells
        points=cell.vertices;source_points=_MDSMP[];target_points=_MDSMP[]
        for p in points
            push!(source_points,_mdsm_map_point(p,Dict(enumerate(C)),Dict(enumerate(A)),[TA[cell.source_id]]).point)
            push!(target_points,_mdsm_map_point(p,Dict(enumerate(C)),Dict(enumerate(B)),[TB[cell.target_id]]).point)
        end
        _mdsm_need(_mpo_area(source_points)>0 && _mpo_area(target_points)>0,"Collar common cell loses orientation")
        push!(maps,(;source=source_points,target=target_points))
    end
    maps
end

"""Reconstruct a whole collar layer and check both complete annuli.

Source: `_mdsm_collar_sector` and `_mdsm_mesh`. Each exact common polygon is
fanned about the mean of its vertices; both affine images use the SAME mean.
Coordinate-keyed vertices force all shared radial bridges to agree pointwise.
The complete directed outer/inner words are checked, not only an area sum.
"""
function _mdsm_collar_layer(input;progress=identity)
    result=Any[]
    geometric_word(q)=_mdsm_boundary_word(Dict(enumerate(q)),collect(eachindex(q)))
    for boundary in input["boundaries"]
        P=Dict{Int,_MDSMP}();Q=Dict{Int,_MDSMP}()
        ids=Dict{_MDSMP,Int}();inverse=Dict{_MDSMP,_MDSMP}();triangles=Vector{Int}[]
        function vertex(p,q)
            _mdsm_need(!haskey(inverse,q) || inverse[q]==p,"Collar target identifies different source points")
            inverse[q]=p
            if haskey(ids,p)
                i=ids[p];_mdsm_need(Q[i]==q,"Shared collar bridge has inconsistent images")
                return i
            end
            i=length(ids)+1;ids[p]=i;P[i]=p;Q[i]=q;i
        end
        cells=0
        for sector in boundary["sectors"]
            for cell in _mdsm_collar_sector(sector)
                n=length(cell.source)
                a=ntuple(k->sum(p[k] for p in cell.source)/n,2)
                b=ntuple(k->sum(p[k] for p in cell.target)/n,2)
                center=vertex(a,b)
                ring=[vertex(p,q) for (p,q) in zip(cell.source,cell.target)]
                for (u,v) in zip(ring,circshift(ring,-1));push!(triangles,[center,u,v]);end
                cells+=1
            end
            for (p,q) in zip(_mdsm_point.(sector["source"]),_mdsm_point.(sector["target"]))
                _mdsm_need(haskey(ids,p) && Q[ids[p]]==q,"Collar lost a prescribed boundary parameter point")
            end
        end
        source=_mdsm_mesh(P,triangles);target=_mdsm_mesh(Q,triangles;expected_boundary=source.cycles)
        _mdsm_need(length(source.cycles)==2,"Collar is not one annulus")
        outer=source.outer;inner=only(filter(c->c!=outer,source.cycles))
        _mdsm_need(all(P[v]==Q[v] for v in outer),"Collar exterior join is not pointwise fixed")
        outer_word=geometric_word(_mdsm_point.(boundary["source_outer"]))
        _mdsm_need(_mdsm_boundary_word(P,outer)==outer_word && _mdsm_boundary_word(Q,outer)==outer_word,"Collar exterior differs from retained chart")
        _mdsm_need(_mdsm_boundary_word(P,inner)==geometric_word(_mdsm_point.(boundary["source_inner"])),"Collar source boundary differs from retained chart")
        _mdsm_need(_mdsm_boundary_word(Q,inner)==geometric_word(_mdsm_point.(boundary["target_inner"])),"Collar target boundary differs from retained chart")
        push!(result,(;cut=boundary["cut"],side=boundary["side"],source=P,target=Q,triangles,outer,inner,common_cells=cells))
        progress(string(input["layer"]," ",boundary["cut"]," ",boundary["side"]," ",cells," common cells"))
    end
    result
end
