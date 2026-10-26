"""Exact triangle intersections for the marked annulus action's common mesh.

This is a thin labelled-cell adapter, not a clipping implementation. Polyhedra
owns intersection/H-to-V conversion. The caller owns periodic lifts, the annular
map, source/target tiling checks and native chart assembly. Dependency versions
belong to the package Project.toml/Manifest.toml.
"""

function _mpo_point(point)
    _mm_need((point isa Tuple || point isa AbstractVector) && length(point)==2,
        "overlay points must have two rational coordinates")
    _mm_need(all(x->x isa Rational,point),"overlay coordinates must be exact rationals; floats are not accepted")
    (MarkedSurfaceDiagrams.Rat(point[1]),MarkedSurfaceDiagrams.Rat(point[2]))
end

function _mpo_triangle(item)
    _mm_need(hasproperty(item,:id) && hasproperty(item,:points),"overlay triangles need id and points")
    points=_mpo_point.(collect(item.points))
    _mm_need(length(points)==3,"overlay input must be a triangle")
    _mm_need(MarkedSurfaceDiagrams.orient(points...)>0,"overlay triangles must be nondegenerate and counterclockwise")
    (;id=deepcopy(item.id),points)
end

_mpo_area(points)=sum(a[1]*b[2]-a[2]*b[1] for (a,b) in zip(points,circshift(points,-1)))/2

# Source: Polyhedra.HalfSpace/hrep and Base.intersect documented at
# https://juliapolyhedra.github.io/Polyhedra.jl/stable/representation/
# and /stable/utilities/. For each CCW edge, its right normal gives n*x <= n*p.
# All coefficients are Rational{BigInt}; do not route through a Float64 hull.
function _mpo_halfspaces(points)
    Q=MarkedSurfaceDiagrams.Rat
    constraints=map(zip(points,circshift(points,-1))) do (a,b)
        nx,ny=b[2]-a[2],a[1]-b[1]
        Polyhedra.HalfSpace(Q[nx,ny],nx*a[1]+ny*a[2])
    end
    Polyhedra.hrep(constraints)
end

# Only order the vertices returned by the library. This exact polar comparator
# performs no intersection or hull computation and retains all supplied points.
function _mpo_order(points)
    Q=MarkedSurfaceDiagrams.Rat
    center=(sum(p[1] for p in points)/length(points),sum(p[2] for p in points)/length(points))
    upper(p)=p[2]>center[2] || (p[2]==center[2] && p[1]>=center[1])
    function earlier(a,b)
        upper(a)!=upper(b) && return upper(a)
        cross=(a[1]-center[1])*(b[2]-center[2])-(a[2]-center[2])*(b[1]-center[1])
        cross!=0 && return cross>0
        da=(a[1]-center[1])^2+(a[2]-center[2])^2
        db=(b[1]-center[1])^2+(b[2]-center[2])^2
        da<db
    end
    ordered=sort(copy(points);lt=earlier)
    _mm_need(_mpo_area(ordered)>zero(Q),"Polyhedra intersection boundary has no positive oriented area")
    ordered
end

function _mpo_segment_points(a,b,registry)
    _mm_need(a!=b,"overlay boundary segment has coincident endpoints")
    axis=a[1]!=b[1] ? 1 : 2
    points=[p for p in registry if MarkedSurfaceDiagrams.onsegment(a,b,p)]
    sort!(points;by=p->(p[axis]-a[axis])/(b[axis]-a[axis]))
    _mm_need(first(points)==a && last(points)==b,"overlay registry lost a boundary endpoint")
    points
end

# Source: Polyhedra.doubledescription and points. The official implementation
# src/doubledescription.jl at 3df3f01b048ab4cd4b1f48229eae907efb0313b3
# returns a nonredundant V-representation; src/comp.jl uses exact comparisons
# for rational coefficients. Hence every omitted collinear marked subdivision
# point must be reinserted from our shared registry after ALL pairs are computed.
function _mpo_intersection(source,target)
    # Source: coordinate extrema bound every convex combination of a triangle's
    # vertices. Strict separation on either axis therefore proves emptiness.
    # Equality is NOT rejection: all point/edge contacts still reach Polyhedra.
    # Inputs have already passed _mpo_triangle; retain the full pair iteration.
    for axis in 1:2
        slo,shi=extrema(p[axis] for p in source.points)
        tlo,thi=extrema(p[axis] for p in target.points)
        (shi<tlo || thi<slo) && return (;dimension=-1,points=MarkedSurfaceDiagrams.Point[])
    end
    intersection=intersect(_mpo_halfspaces(source.points),_mpo_halfspaces(target.points))
    vertices=Polyhedra.doubledescription(intersection;tol=0)
    _mm_need(isempty(collect(Polyhedra.rays(vertices))) && isempty(collect(Polyhedra.lines(vertices))),
        "bounded triangle intersection unexpectedly contains a ray or line")
    points=sort!(unique(_mpo_point.(collect(Polyhedra.points(vertices)))))
    isempty(points) && return (;dimension=-1,points)
    length(points)==1 && return (;dimension=0,points)
    collinear=all(p->MarkedSurfaceDiagrams.orient(first(points),last(points),p)==0,points)
    collinear && return (;dimension=1,points=[first(points),last(points)])
    (;dimension=2,points=_mpo_order(points))
end

"""
    _marked_triangle_overlay(source_triangles, target_triangles; registry_points=[])

Each input item has `id` (an opaque JSON label) and `points` (three CCW rational
2D points). Source copies must already have their periodic shifts applied and
encoded in their labels. No angular or marking semantics are inferred here.

Returns `(cells, contacts, registry, coverage)`. Positive-area cells contain
`source_id`, `target_id`, a CCW `vertices` boundary retaining every registered
collinear point, exact `centroid`, and exact `area`. Contacts contain source and
target labels, `dimension` 0 or 1 and ordered `vertices`. Every input endpoint,
extra registry point and intersection/contact vertex is in the sorted registry.
`coverage` reports target area versus the sum of intersection-cell areas, with
`area_matches`; this is not a no-overlap/tiling proof for arbitrary caller input.
"""
function _marked_triangle_overlay(source_triangles,target_triangles;registry_points=[])
    sources=_mpo_triangle.(collect(source_triangles))
    targets=_mpo_triangle.(collect(target_triangles))
    for (name,items) in (("source",sources),("target",targets))
        keys=[_canonical_json(t.id) for t in items]
        _mm_need(length(keys)==length(Set(keys)),"overlay $name triangle labels must be unique")
    end
    Point=MarkedSurfaceDiagrams.Point;Q=MarkedSurfaceDiagrams.Rat
    registry=Set{Point}(_mpo_point.(collect(registry_points)))
    for triangle in [sources;targets],point in triangle.points;push!(registry,point);end
    raw=Any[]
    for (i,source) in enumerate(sources),(j,target) in enumerate(targets)
        result=_mpo_intersection(source,target)
        result.dimension<0 && continue
        union!(registry,result.points)
        push!(raw,(;source_index=i,target_index=j,dimension=result.dimension,points=result.points))
    end
    shared=sort!(collect(registry));cells=Any[];contacts=Any[];covered=zeros(Q,length(targets))
    for piece in raw
        source_id=deepcopy(sources[piece.source_index].id);target_id=deepcopy(targets[piece.target_index].id)
        if piece.dimension==2
            boundary=Point[]
            for (a,b) in zip(piece.points,circshift(piece.points,-1))
                append!(boundary,_mpo_segment_points(a,b,shared)[1:end-1])
            end
            area=_mpo_area(boundary)
            _mm_need(area>0,"overlay output contains a nonpositive cell")
            # The average of convex extreme vertices is strictly interior.
            # Compute it before collinear subdivision insertion so a caller's
            # extra boundary points do not move the chosen triangulation center.
            centroid=(sum(p[1] for p in piece.points)/length(piece.points),sum(p[2] for p in piece.points)/length(piece.points))
            _mm_need(all(MarkedSurfaceDiagrams.orient(a,b,centroid)>0 for (a,b) in zip(boundary,circshift(boundary,-1))),
                "overlay centroid is not strictly inside its convex cell")
            covered[piece.target_index]+=area
            push!(cells,(;source_id,target_id,vertices=boundary,centroid,area))
        else
            vertices=piece.dimension==0 ? copy(piece.points) : _mpo_segment_points(first(piece.points),last(piece.points),shared)
            push!(contacts,(;source_id,target_id,dimension=piece.dimension,vertices))
        end
    end
    coverage=[(;target_id=deepcopy(t.id),target_area=_mpo_area(t.points),covered_area=covered[j],
        area_matches=covered[j]==_mpo_area(t.points)) for (j,t) in enumerate(targets)]
    (;cells,contacts,registry=shared,coverage)
end
