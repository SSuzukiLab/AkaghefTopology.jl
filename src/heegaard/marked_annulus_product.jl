# Source: marked_subdivision.jl::marked_subdivide's six oriented flags and
# ::marked_parallel_copy's full-core derived closed star; vendor
# MarkedSurfaceDiagrams vertex links, fine_edge_occurrences and regions.
# Construction/proof: one vertex fan per side maps onto a rational rectangle.
# The source cells, their orientation and integer transition lifts are checked
# below. This chart does not itself apply a twist or change a native document.

_ann_q(a,b=1)=BigInt(a)//BigInt(b)
_ann_det(p)=((p[2][1]-p[1][1])*(p[3][2]-p[1][2])-
             (p[2][2]-p[1][2])*(p[3][1]-p[1][1]))

"""Construct an exact oriented product chart on a full-core derived annulus.

`fine` must be exactly one `marked_subdivide` of `coarse`; selected cell
indices are zero-based. Returned triangle points are Rational{BigInt} pairs
in the SAME vertex order as each positively oriented source triangle. Their
u coordinates are occurrence-level lifts, periodic with period one; t=+1
is the left frontier and t=-1 the right. The first saved core dart starts
at u=0. Across any interior edge the two lifts differ by one constant integer.

The chart is geometric provenance for a fixed support. Inverse twists must
use that same chart with negative shear, not reparameterize the output core.
"""
function _marked_annulus_product(coarse,fine,core_id,selected_cells)
    M=MarkedSurfaceDiagrams;cid=String(core_id)
    _mm_need(haskey(coarse.curves,cid) && coarse.curves[cid][1]["closed"],"annular chart requires an existing closed core")
    _mm_need(coarse.surface_report["boundary_components"]==0,"annular chart currently requires a closed source surface")
    _mm_need(coarse.doc["surface"]==fine.doc["surface"],"annular chart subdivision changed the surface marking")
    _mm_need(Set(keys(coarse.patches))==Set(keys(fine.patches)),"annular chart subdivision changed patch identities")
    core=coarse.curves[cid][2];n=length(core)
    _mm_need(n>=3,"annular core needs at least three coarse edge occurrences")
    cr(f,v)=M.root!(coarse.vertex_dsu,(f,v))
    fr(f,v)=M.root!(fine.vertex_dsu,(f,v))
    roots=[cr(d.face,d.u) for d in core];rootset=Set(roots)
    _mm_need(length(rootset)==n,"annular core repeats a quotient vertex before closing")
    _mm_need(all(cr(d.face,d.v)==roots[mod1(i+1,n)] for (i,d) in enumerate(core)),"annular core edge order is not a closed cycle")
    owned=Set((d.face,M.edge(d.u,d.v)) for d in core)
    for (f,inc) in coarse.incidence,(e,_) in inc
        if cr(f,e[1]) in rootset && cr(f,e[2]) in rootset
            _mm_need((f,e) in owned,"annular core has an unused chord edge")
        end
    end
    _mm_need(all(!all(cr(f,v) in rootset for v in t) for ((f,_),t) in coarse.triangles),"annular core contains a whole coarse triangle")

    # Recover the actual six flags, including their orientation, rather than
    # inferring barycenters from generated name prefixes.
    mids=Dict();centers=Dict();apex=Dict();expected=Set{Tuple{String,Int}}()
    for (f,p) in coarse.patches
        fp=fine.patches[f]
        _mm_need(length(fp["triangles"])==6length(p["triangles"]),"annular chart fine triangle count is not one subdivision")
        old=Set(v for t in p["triangles"] for v in t);newcenters=Set()
        for (i,t) in enumerate(p["triangles"])
            a,b,c=t;flags=[fine.triangles[(f,6*(i-1)+j)] for j in 1:6]
            ab,z=flags[1][2:3];bc=flags[3][2];ca=flags[5][2]
            # DiagramState stores oriented triangles as tuples (vendor
            # _validate), although the input document stores JSON arrays.
            _mm_need(flags==[(a,ab,z),(ab,b,z),(b,bc,z),(bc,c,z),(c,ca,z),(ca,a,z)],"annular chart fine cells do not match oriented barycentric flags")
            _mm_need(length(Set([ab,bc,ca,z]))==4 && all(!(v in old) for v in (ab,bc,ca,z)) && !(z in newcenters),"annular chart subdivision vertices are not distinct")
            push!(newcenters,z);centers[(f,i)]=z
            for (e,m) in ((M.edge(a,b),ab),(M.edge(b,c),bc),(M.edge(c,a),ca))
                _mm_need(!haskey(mids,(f,e)) || mids[(f,e)]==m,"annular chart has inconsistent edge midpoints")
                mids[(f,e)]=m
            end
            for (j,v) in enumerate((a,b,b,c,c,a))
                cell=(f,6*(i-1)+j-1);apex[cell]=cr(f,v)
                cr(f,v) in rootset && push!(expected,cell)
            end
        end
        ms=[m for ((g,_),m) in mids if g==f]
        _mm_need(length(unique(ms))==length(ms) && isempty(intersect(Set(ms),newcenters)),"annular chart midpoint or barycenter identities collide")
        _mm_need(p["side_order"]==fp["side_order"],"annular chart subdivision changed marking sides")
        for side in p["side_order"]
            path=p["sides"][side];refined=String[];params=M.side_parameters(p,side);ps=typeof(_ann_q(0))[]
            for j in 1:length(path)-1
                append!(refined,[path[j],mids[(f,M.edge(path[j],path[j+1]))]])
                append!(ps,[params[j],(params[j]+params[j+1])/2])
            end
            push!(refined,path[end]);push!(ps,params[end])
            _mm_need(fp["sides"][side]==refined && M.side_parameters(fp,side)==ps,"annular chart subdivision changed seam parameterization")
        end
    end
    selected=Set((String(c[1]),Int(c[2])) for c in selected_cells)
    _mm_need(selected==expected,"annular chart support is not the exact full-core derived star")
    _mm_need(haskey(fine.curves,cid),"fine annular chart has no core curve")
    doubled=M.Dart[]
    for d in core
        m=mids[(d.face,M.edge(d.u,d.v))]
        append!(doubled,[M.Dart(d.face,d.u,m),M.Dart(d.face,m,d.v)])
    end
    _mm_need(fine.curves[cid][2]==doubled,"fine annular core does not follow its actual subdivision")

    pairing=M.fine_edge_occurrences(fine);barriers=Set()
    for (a,b) in pairing
        (a[1],a[2]-1) in selected && !((b[1],b[2]-1) in selected) && push!(barriers,(a[1],M.edge(a[3],a[4])))
    end
    regions=M.regions(fine;barriers=barriers)
    matches=[r for r in regions if _mm_cells(r)==selected]
    _mm_need(length(matches)==1,"annular chart support is not one normalized region")
    annulus=only(matches)
    _mm_need(annulus["genus"]==0 && annulus["boundary_components"]==2,"annular chart support is not an annulus")

    # The two paths in the fine vertex link are the left/right fans. Every
    # flag has one old apex; its opposite edge is an actual link edge.
    cells_by_root=Dict(r=>Tuple{String,Int}[] for r in roots)
    for c in selected;push!(cells_by_root[apex[c]],c);end
    coordinates=Dict();fan_records=Any[]
    for (k,r) in enumerate(roots)
        i=k-1;adj=Dict();linkcell=Dict()
        for cell in cells_by_root[r]
            f,j=cell;t=fine.triangles[(f,j+1)]
            other=[fr(f,v) for v in t if !(haskey(coarse.vertex_dsu.p,(f,v)) && cr(f,v)==r)]
            _mm_need(length(other)==2 && other[1]!=other[2],"annular flag does not have exactly one old apex")
            a,b=other;e=isless(a,b) ? (a,b) : (b,a)
            _mm_need(!haskey(linkcell,e),"annular vertex fan repeats a link edge")
            linkcell[e]=cell;push!(get!(adj,a,Set()),b);push!(get!(adj,b,Set()),a)
        end
        _mm_need(all(length(ns)==2 for ns in values(adj)),"annular apex link is not a circular fan")
        previous=core[mod1(k-1,n)];following=core[k]
        min=fr(previous.face,mids[(previous.face,M.edge(previous.u,previous.v))])
        mout=fr(following.face,mids[(following.face,M.edge(following.u,following.v))])
        _mm_need(min!=mout,"annular fan has coincident core midpoint endpoints")
        for side in (1,-1)
            incident=coarse.incidence[previous.face][M.edge(previous.u,previous.v)]
            onside=[o for o in incident if (o[2]==previous.u && o[3]==previous.v)==(side==1)]
            _mm_need(length(onside)==1,"annular core edge has no unique oriented side triangle")
            firstcenter=fr(previous.face,centers[(previous.face,only(onside)[1])])
            _mm_need(firstcenter in get(adj,min,Set()),"annular fan seam is not an actual flag edge")
            path=[min,firstcenter]
            while path[end]!=mout
                next=only(setdiff(adj[path[end]],Set([path[end-1]])))
                _mm_need(!(next in path),"annular side link closes before reaching the next core edge")
                push!(path,next)
            end
            _mm_need(length(path)>=4,"annular fan has no nondegenerate frontier interval")
            q=length(path)-3
            points=Dict(min=>(_ann_q(2i-1,2n),_ann_q(0)),mout=>(_ann_q(2i+1,2n),_ann_q(0)))
            for j in 0:q;points[path[j+2]]=(_ann_q(2i-1,2n)+_ann_q(j,n*q),_ann_q(side));end
            localcells=Any[]
            for j in 1:length(path)-1
                a,b=path[j:j+1];e=isless(a,b) ? (a,b) : (b,a);cell=linkcell[e]
                _mm_need(!haskey(coordinates,cell),"annular chart covers a fine triangle twice")
                f,idx=cell;t=fine.triangles[(f,idx+1)]
                p=[haskey(coarse.vertex_dsu.p,(f,v)) && cr(f,v)==r ? (_ann_q(i,n),_ann_q(0)) : points[fr(f,v)] for v in t]
                _mm_need(_ann_det(p)>0,"annular fan chart has a non-positive oriented triangle")
                coordinates[cell]=p;push!(localcells,collect(cell))
            end
            push!(fan_records,Dict("core_vertex"=>collect(r),"index"=>i,"side"=>side,"frontier_edges"=>q,
                "link_path"=>[collect(v) for v in path],"triangles"=>localcells))
        end
    end
    _mm_need(Set(keys(coordinates))==selected,"annular chart does not cover exactly its support")
    point(cell,v)=coordinates[cell][only(findall(==(v),fine.triangles[(cell[1],cell[2]+1)]))]
    transitions=Any[]
    for (a,b) in sort(collect(pairing);by=first)
        ca=(a[1],a[2]-1);cb=(b[1],b[2]-1)
        ca in selected && cb in selected && isless(a,b) || continue
        p,q=point(ca,a[3]),point(ca,a[4]);r,s=point(cb,b[4]),point(cb,b[3])
        shift=r[1]-p[1]
        _mm_need(denominator(shift)==1 && s[1]-q[1]==shift && r[2]==p[2] && s[2]==q[2],"annular chart edge coordinates do not have one integer transition lift")
        push!(transitions,Dict("source"=>[a[1],a[2]-1,a[3],a[4]],"target"=>[b[1],b[2]-1,b[4],b[3]],"u_shift"=>string(numerator(shift))))
    end
    frontier_records=Any[]
    for cycle in annulus["boundary_cycles"]
        levels=Set()
        for (f,j,u,v) in cycle
            p,q=point((f,j),u),point((f,j),v)
            _mm_need(p[2]==q[2] && abs(p[2])==1,"annular frontier is not at constant unit height")
            push!(levels,p[2])
        end
        _mm_need(length(levels)==1,"annular frontier changes sides")
        push!(frontier_records,Dict("side"=>Int(only(levels)),"cycle"=>cycle))
    end
    _mm_need(Set(x["side"] for x in frontier_records)==Set([-1,1]),"annular frontiers do not have opposite sides")
    area=sum(_ann_det(p)/2 for p in values(coordinates))
    _mm_need(area==2,"annular product chart does not have unit period and total area two")
    winding=_ann_q(0)
    for d in doubled
        candidates=[(d.face,o[1]-1) for o in fine.incidence[d.face][M.edge(d.u,d.v)] if (d.face,o[1]-1) in selected]
        _mm_need(!isempty(candidates),"annular core edge is outside the support")
        c=first(candidates);a,b=point(c,d.u),point(c,d.v)
        _mm_need(a[2]==b[2]==0 && b[1]-a[1]==_ann_q(1,2n),"annular core coordinate does not follow its positive unit winding")
        winding+=b[1]-a[1]
    end
    _mm_need(winding==1,"annular chart core winding is not one")
    triangles=[(id=collect(c),points=coordinates[c]) for c in sort(collect(selected))]
    certificate=Dict("construction"=>"oriented full-core barycentric vertex fans","core"=>cid,
        "core_order"=>[[d.face,d.u,d.v] for d in core],"period"=>1,"left_height"=>1,"right_height"=>-1,
        "fans"=>fan_records,"edge_transitions"=>transitions,"frontiers"=>frontier_records,
        "area"=>M.qstr(area),"core_winding"=>M.qstr(winding),
        "triangle_coordinates"=>[Dict("id"=>row.id,"points"=>[[M.qstr(v) for v in p] for p in row.points]) for row in triangles],
        "inverse_support"=>"reuse these fixed source triangles and integer lifts; apply the opposite shear")
    (;triangles,period=_ann_q(1),support_cells=selected,
      core_order=certificate["core_order"],seam_lifts=transitions,certificate)
end
