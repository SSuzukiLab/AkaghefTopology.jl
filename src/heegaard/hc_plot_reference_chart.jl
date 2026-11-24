# cf: lib/topology.jl/docs/specs/heegaard-plot.md — common chart transport for fixed reference geometry.
# cf: .research-scratch/PJ260922/astra-lanes/hc-plot/answer-5c.md — persistent center anchors and arrangement preflight; reconstructed.

function _hp_reference_circle_events!(ys,a,b,c,r)
    u=_hp_sub(b,a);v=_hp_sub(a,c);A=_hp_dot(u,u)
    A>1e-20 || return
    B=2*_hp_dot(v,u);C=_hp_dot(v,v)-r*r;D=B*B-4A*C
    # Recompute a near-tangent discriminant at higher precision, before
    # deciding whether two arrangement events exist.
    if abs(D)<64eps(Float64)*max(B*B,abs(4A*C),1.)
        D=Float64(setprecision(BigFloat,256) do
            ux,uy=BigFloat(u[1]),BigFloat(u[2]);vx,vy=BigFloat(v[1]),BigFloat(v[2])
            aa=ux*ux+uy*uy;bb=2(vx*ux+vy*uy);cc=vx*vx+vy*vy-BigFloat(r)^2
            bb*bb-4aa*cc
        end)
    end
    D>=0 || return
    z=sqrt(D)
    for t in ((-B-z)/(2A),(-B+z)/(2A))
        -1e-10<=t<=1+1e-10 && push!(ys,a[2]+clamp(t,0.,1.)*u[2])
    end
end

function _hp_reference_straighten(ctx)
    p=ctx.placement
    p.square && return ctx
    spines=[r for r in p.runs if startswith(r.name,"@K:m")]
    isempty(spines) && return ctx
    for i in eachindex(p.holes),j in i+1:length(p.holes)
        a,ra=p.holes[i].second;b,rb=p.holes[j].second
        _hp_dist(a,b)>ra+rb+1e-8 || throw(ArgumentError("reference chart needs disjoint hole disks"))
    end
    all(all(r.points[j][2]>r.points[j+1][2] for j in 1:length(r.points)-1) for r in spines) ||
        throw(ArgumentError("reference meridian is not a vertical graph"))
    xmin=minimum((q[1] for r in p.runs for q in r.points);init=0.)-2p.width-1000.
    xmax=maximum((q[1] for r in p.runs for q in r.points);init=p.width)+2p.width+1000.
    function spine_x(r,y)
        P=r.points
        y>=P[1][2] && return P[1][1]
        y<=P[end][2] && return P[end][1]
        for j in 1:length(P)-1
            a,b=P[j],P[j+1]
            b[2]<=y<=a[2] || continue
            return a[1]+(b[1]-a[1])*(y-a[2])/(b[2]-a[2])
        end
        error("meridian graph interval missing")
    end
    function chart_anchors(y)
        anchors=Tuple{Float64,Float64}[(xmin,xmin),(xmax,xmax)]
        for r in spines
            target=r.points[1][1]
            abs(target-r.points[end][1])<1e-7 || throw(ArgumentError("meridian ports do not have the same x coordinate"))
            push!(anchors,(spine_x(r,y),target))
        end
        for (_,(c,r)) in p.holes
            # This anchor persists when the two circle intersections disappear
            # at an extremum; otherwise the fiber map can jump there.
            push!(anchors,(c[1],c[1]))
            abs(y-c[2])<=r || continue
            dx=sqrt(max(0.,r*r-(y-c[2])^2))
            push!(anchors,(c[1]-dx,c[1]-dx),(c[1]+dx,c[1]+dx))
        end
        sort!(anchors;by=first);unique_anchors=Tuple{Float64,Float64}[]
        for a in anchors
            if !isempty(unique_anchors) && abs(a[1]-unique_anchors[end][1])<1e-8
                abs(a[2]-unique_anchors[end][2])<1e-7 || throw(ArgumentError("incompatible chart anchors"))
                continue
            end
            push!(unique_anchors,a)
        end
        all(unique_anchors[i][2]<unique_anchors[i+1][2] for i in 1:length(unique_anchors)-1) ||
            throw(ArgumentError("reference chart transport is not increasing"))
        unique_anchors
    end
    ys=Float64[q[2] for r in spines for q in r.points]
    for (_,(c,rad)) in p.holes
        append!(ys,(c[2]-rad,c[2]+rad))
        for r in spines
            for j in 1:length(r.points)-1
                _hp_reference_circle_events!(ys,r.points[j],r.points[j+1],c,rad)
            end
            t=r.points[1][1];z=rad*rad-(t-c[1])^2
            z>=0 && append!(ys,(c[2]-sqrt(z),c[2]+sqrt(z)))
        end
    end
    sort!(ys);events=Float64[]
    for y in ys;(isempty(events) || abs(y-events[end])>1e-9) && push!(events,y);end
    probes=vcat(events,[(events[i]+events[i+1])/2 for i in 1:length(events)-1])
    foreach(chart_anchors,probes)
    function transport(q)
        x,y=q;unique_anchors=chart_anchors(y)
        for i in 1:length(unique_anchors)-1
            a,b=unique_anchors[i],unique_anchors[i+1]
            a[1]-1e-9<=x<=b[1]+1e-9 || continue
            return (a[2]+(b[2]-a[2])*(x-a[1])/(b[1]-a[1]),y)
        end
        q
    end
    function mapped(r)
        if startswith(r.name,"@K:m")
            return _hp_run(r,[r.points[1],r.points[end]])
        end
        P=r.points;points=HCPlotPoint[]
        for i in 1:length(P)-1
            a,b=P[i],P[i+1];count=max(1,ceil(Int,_hp_dist(a,b)/4.))
            for j in 0:count-1;push!(points,transport(_hp_lerp(a,b,j/count)));end
        end
        push!(points,transport(P[end]))
        # Horizontal K edges stay horizontal under this map; retain only their
        # ends. Curved collars and all beta vertices keep their mapped samples.
        if r.name in ctx.fixed
            clean=HCPlotPoint[]
            for q in points
                while length(clean)>=2
                    a,b=clean[end-1],clean[end];u=_hp_sub(b,a);v=_hp_sub(q,b)
                    abs(_hp_cross(u,v))<=1e-10 && _hp_dot(u,v)>=0 || break
                    pop!(clean)
                end
                push!(clean,q)
            end
            points=clean
        end
        _hp_run(r,points)
    end
    candidate=_hp_rebuild(p,[mapped(r) for r in p.runs])
    # The analytic fiber map fixes every circle boundary and is increasing in x.
    # Its polygonal sampling is independently checked against the source LC by
    # the caller; a failed approximation is an error, never a fallback drawing.
    (placement=candidate,contacts=ctx.contacts,fixed=ctx.fixed)
end
