"""Joining-disk geometry derived from boundary darts and the common surface atlas.

No example-specific coordinates. Alpha boundary intervals become a boundary
circle arc in the base and an annular sector; beta intervals reuse the exact
layout Beziers, reversing their controls when traversed backwards.
"""

function _disk_circle_tail(center,radius,theta,delta)
    n=max(1,ceil(Int,abs(delta)/π))
    join([begin
        a=theta+delta*k/n
        p=_vadd(center,(radius*cos(a),radius*sin(a)))
        "A $(_fmt(radius)) $(_fmt(radius)) 0 0 $(delta>0 ? 0 : 1) $(_svg_point(p))"
    end for k in 1:n]," ")
end

function _disk_route_geometry(route,forward)
    points=forward ? copy(route.points) : reverse(route.points)
    if route.curve isa BezierChain
        segments=forward ? route.curve.segments : reverse(route.curve.segments)
        tail=join([begin
            a,b,c=forward ? (s.p1,s.p2,s.p3) : (s.p2,s.p1,s.p0)
            "C $(_svg_point(a)) $(_svg_point(b)) $(_svg_point(c))"
        end for s in segments]," ")
    else
        tail=join(("L $(_svg_point(p))" for p in points[2:end])," ")
    end
    (;points,tail)
end

function _disk_contains(points,p)
    inside=false
    for i in eachindex(points)
        a,b=points[i],points[mod1(i+1,length(points))]
        if (a[2]>p[2]) != (b[2]>p[2]) &&
           p[1] < (b[1]-a[1])*(p[2]-a[2])/(b[2]-a[2])+a[1]
            inside=!inside
        end
    end
    inside
end

function _disk_base_label(points,exterior)
    xmin,xmax=extrema(first.(points)); ymin,ymax=extrema(last.(points))
    exterior && return (xmin-0.06max(xmax-xmin,1.),ymax+0.06max(ymax-ymin,1.))
    best=nothing; distance=-Inf
    for x in range(xmin,xmax;length=21),y in range(ymin,ymax;length=21)
        p=(x,y)
        _disk_contains(points,p) || continue
        d=minimum(_segment_distance_squared(p,points[i],points[mod1(i+1,length(points))])
                  for i in eachindex(points))
        if d>distance
            best=p; distance=d
        end
    end
    best # Very narrow regions may have no sampled label position; don't invent one.
end

function _joining_disk_layout(H::HeegaardDiagram,D::JoiningDisk,chart::HeegaardSurfaceChart)
    validate_joining_disk(H,D)
    canonical_hash(H)==chart.diagram_hash || throw(ArgumentError("joining disk chart hash mismatch"))
    graph=ribbon_graph(H)
    curves=Dict((c.family,c.id)=>c for c in _curve_vector(H.code.code))
    intersections=Dict(p.id=>p for p in chart.intersections)
    punctures=Dict(p.id=>p for p in chart.base.punctures)
    ports=Dict(p.id=>p for p in chart.base.ports)
    boundary=NamedTuple[]; sectors=NamedTuple[]
    for id in D.boundary_darts
        dart=graph.darts[id]; edge=graph.edges[dart.edge]
        curve=curves[(edge.family,edge.curve_id)]
        forward=isodd(id)
        if edge.family==:beta
            k=findfirst(==(edge.from),curve.visits)
            route=only(r for r in chart.base.routes if r.curve_id==edge.curve_id && r.segment_id==k)
            data=_disk_route_geometry(route,forward)
            push!(boundary,(dart=id,family=:beta,curve_id=edge.curve_id,
                from=dart.origin,to=dart.target,points=data.points,tail=data.tail))
        else
            start=intersections[dart.origin]; finish=intersections[dart.target]
            a=only(a for a in chart.annuli if a.alpha_id==edge.curve_id)
            side=forward ? :plus : :minus
            port1=ports[forward ? start.plus_port : start.minus_port]
            port2=ports[forward ? finish.plus_port : finish.minus_port]
            p=punctures[port1.boundary_id]
            port2.boundary_id==p.id || error("alpha interval changes cut boundary")
            theta=2π*port1.parameter; delta=2π/length(curve.visits)
            count=max(8,ceil(Int,delta/(π/64)))
            points=[_vadd(p.center,(p.radius*cos(t),p.radius*sin(t)))
                for t in range(theta,theta+delta;length=count+1)]
            _vnorm(_vsub(last(points),port2.point))<1e-8 || error("alpha boundary cyclic order mismatch")
            points[1]=port1.point; points[end]=port2.point
            tail=_disk_circle_tail(p.center,p.radius,theta,delta)
            push!(boundary,(dart=id,family=:alpha,curve_id=edge.curve_id,
                from=dart.origin,to=dart.target,points=points,tail=tail))
            # alpha+ has D on the outer side of the core; alpha- on the inner.
            radius=forward ? a.outer_radius : a.inner_radius
            t0=-2π*start.parameter; dt=(forward ? -1 : 1)*2π/length(curve.visits)
            outer_end=(radius*cos(t0+dt),radius*sin(t0+dt))
            path="M $(_svg_point(start.point)) " *
                _disk_circle_tail((0.,0.),a.core_radius,t0,dt) *
                " L $(_svg_point(outer_end)) " *
                _disk_circle_tail((0.,0.),radius,t0+dt,-dt) * " Z"
            label_radius=(radius+a.core_radius)/2
            push!(sectors,(dart=id,alpha_id=edge.curve_id,side=side,
                from=dart.origin,to=dart.target,path=path,
                label=(label_radius*cos(t0+dt/2),label_radius*sin(t0+dt/2))))
        end
    end
    for k in eachindex(boundary)
        _vnorm(_vsub(last(boundary[k].points),first(boundary[mod1(k+1,length(boundary))].points)))<1e-8 ||
            throw(ArgumentError("joining boundary is discontinuous in the alpha-cut chart"))
    end
    points=reduce(vcat,(b.points[1:end-1] for b in boundary))
    area=sum(points[i][1]*points[mod1(i+1,length(points))][2]-
             points[i][2]*points[mod1(i+1,length(points))][1] for i in eachindex(points))/2
    abs(area)>1e-10 || throw(ArgumentError("joining boundary has zero planar area"))
    # The base chart orientation is -dx wedge dy: bounded disk boundary is clockwise.
    exterior=area>0
    path="M $(_svg_point(first(boundary).points[1])) " * join((b.tail for b in boundary)," ") * " Z"
    (disk=deepcopy(D),boundary=boundary,sectors=sectors,base_path=path,
        base_points=points,base_exterior=exterior,base_label=_disk_base_label(points,exterior))
end

function _joining_overlay_json(overlay)
    overlay===nothing && return nothing
    Dict("disk"=>_disk_json(overlay.disk),"boundary"=>[
        Dict("dart"=>b.dart,"family"=>String(b.family),"curve_id"=>b.curve_id,
             "from"=>b.from,"to"=>b.to) for b in overlay.boundary],
        "annulus_sectors"=>overlay.sectors,"base_path"=>overlay.base_path,
        "base_exterior"=>overlay.base_exterior,"base_label"=>overlay.base_label)
end
