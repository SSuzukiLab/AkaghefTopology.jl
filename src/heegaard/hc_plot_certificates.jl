# cf: .research-scratch/PJ260922/astra-lanes/hc-plot/answer-3.md §§3–7,11 — pre-motion sufficient predicates; exact rational fallback for binary input coordinates
# These predicates certify local claims only. Integration into a full certified
# smoothing loop is tracked separately; their presence does not enable production.
const HPCertRational = Rational{BigInt}
_hp_q(x::Real)=HPCertRational(x)
_hp_qpoint(p)=(_hp_q(p[1]),_hp_q(p[2]))

"""Outward bounds for basic IEEE Float64 operations; no approximate epsilon decisions."""
struct HPFloatBounds
    lo::Float64
    hi::Float64
end
HPFloatBounds(x::Float64)=HPFloatBounds(x,x)
function HPFloatBounds(x::Real)
    y=Float64(x)
    isfinite(y) || return HPFloatBounds(-Inf,Inf)
    q=_hp_q(x);qy=_hp_q(y)
    HPFloatBounds(qy>q ? prevfloat(y) : y,qy<q ? nextfloat(y) : y)
end
Base.:+(a::HPFloatBounds,b::HPFloatBounds)=HPFloatBounds(prevfloat(a.lo+b.lo),nextfloat(a.hi+b.hi))
Base.:-(a::HPFloatBounds,b::HPFloatBounds)=HPFloatBounds(prevfloat(a.lo-b.hi),nextfloat(a.hi-b.lo))
function Base.:*(a::HPFloatBounds,b::HPFloatBounds)
    p=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi)
    any(isnan,p) && return HPFloatBounds(-Inf,Inf)
    HPFloatBounds(prevfloat(minimum(p)),nextfloat(maximum(p)))
end
function Base.:/(a::HPFloatBounds,b::HPFloatBounds)
    b.lo<=0<=b.hi && return HPFloatBounds(-Inf,Inf)
    a*HPFloatBounds(prevfloat(1/b.hi),nextfloat(1/b.lo))
end
_hp_bpoint(p)=(HPFloatBounds(p[1]),HPFloatBounds(p[2]))
function _hp_orientation_sign(a,b,c)
    d=_hp_cross(_hp_sub(_hp_bpoint(b),_hp_bpoint(a)),_hp_sub(_hp_bpoint(c),_hp_bpoint(a)))
    d.lo>0 && return 1
    d.hi<0 && return -1
    Int(sign(_hp_cross(_hp_sub(_hp_qpoint(b),_hp_qpoint(a)),_hp_sub(_hp_qpoint(c),_hp_qpoint(a)))))
end
function _hp_closed_segment_contact(a,b,c,d)
    max(a[1],b[1])<min(c[1],d[1]) && return false
    max(c[1],d[1])<min(a[1],b[1]) && return false
    max(a[2],b[2])<min(c[2],d[2]) && return false
    max(c[2],d[2])<min(a[2],b[2]) && return false
    s1=_hp_orientation_sign(a,b,c);s2=_hp_orientation_sign(a,b,d)
    s3=_hp_orientation_sign(c,d,a);s4=_hp_orientation_sign(c,d,b)
    s1*s2<=0 && s3*s4<=0
end
function _hp_point_segment_d2q(p,a,b)
    u=_hp_sub(b,a);v=_hp_sub(p,a);l2=_hp_dot(u,u)
    l2==0 && return _hp_dot(v,v)
    t=clamp(_hp_dot(v,u)/l2,zero(l2),one(l2));w=_hp_sub(v,_hp_mul(u,t))
    _hp_dot(w,w)
end
function _hp_segment_d2q(a,b,c,d)
    _hp_closed_segment_contact(a,b,c,d) && return _hp_q(0)
    a,b,c,d=_hp_qpoint.((a,b,c,d))
    minimum((_hp_point_segment_d2q(a,c,d),_hp_point_segment_d2q(b,c,d),
             _hp_point_segment_d2q(c,a,b),_hp_point_segment_d2q(d,a,b)))
end
function _hp_point_segment_d2bound(p,a,b)
    p,a,b=_hp_bpoint.((p,a,b));u=_hp_sub(b,a);v=_hp_sub(p,a)
    l2=_hp_dot(u,u);t=_hp_dot(v,u)
    t.hi<0 && return _hp_dot(v,v)
    z=_hp_sub(p,b);_hp_dot(z,u).lo>0 && return _hp_dot(z,z)
    area=_hp_cross(u,v);area*area/l2
end
function _hp_segment_distance_lower(a,b,c,d)
    _hp_closed_segment_contact(a,b,c,d) && return 0.0
    d2=minimum((_hp_point_segment_d2bound(a,c,d).lo,_hp_point_segment_d2bound(b,c,d).lo,
                _hp_point_segment_d2bound(c,a,b).lo,_hp_point_segment_d2bound(d,a,b).lo))
    d2>0 ? max(0.0,prevfloat(sqrt(d2))) : 0.0
end
"""Certify separated swept tubes using supplied upper bounds on endpoint path lengths."""
function _hp_cert_swept_pair(a,b,c,d,ma,mb,mc,md)
    all(x->isfinite(x) && x>=0,(ma,mb,mc,md)) || return false
    mu=max(ma,mb,mc,md);lower=_hp_segment_distance_lower(a,b,c,d)
    lower>0 && mu<=prevfloat(.3lower) && return true
    d2=_hp_segment_d2q(a,b,c,d)
    d2>0 && 100*_hp_q(mu)^2<=9d2
end

struct HPCrossingChart
    origin::HCPlotPoint
    transform::NTuple{4,Float64}
    ka::Float64
    kb::Float64
    x_inner::Float64
    y_inner::Float64
    sign::Int
end
function _hp_chart_coordinate(ch,p,cast)
    a,b,c,d=cast.(ch.transform);x,y=_hp_sub(cast.(p),cast.(ch.origin))
    (a*x+b*y,c*x+d*y)
end
function _hp_chart_graph_q(ch,A,B)
    ka,kb,xi,yi=_hp_q.((ch.ka,ch.kb,ch.x_inner,ch.y_inner))
    0<=ka && 0<=kb && ka*kb<1 && xi>0 && yi>0 || return false
    length(A)>=2 && length(B)>=2 || return false
    x=[_hp_chart_coordinate(ch,p,_hp_q) for p in A]
    y=[_hp_chart_coordinate(ch,p,_hp_q) for p in B]
    first(x)[1]<-xi && last(x)[1]>xi && first(y)[2]<-yi && last(y)[2]>yi || return false
    all(p->abs(p[2])<=yi,x) && all(p->abs(p[1])<=xi,y) || return false
    all(begin v=_hp_sub(x[i+1],x[i]);v[1]>0 && abs(v[2])<=ka*v[1] end for i in 1:length(x)-1) || return false
    all(begin v=_hp_sub(y[i+1],y[i]);v[2]>0 && abs(v[1])<=kb*v[2] end for i in 1:length(y)-1)
end
"""Exact cone-ray angle test. min_sine2=1/4 represents the unchanged 30-degree criterion."""
function _hp_chart_angle_certificate(ch;min_sine2=1//4)
    0<min_sine2<=1 || return false
    a,b,c,d=_hp_q.(ch.transform);det=a*d-b*c
    det!=0 && sign(det)==ch.sign || return false
    ka,kb=_hp_q(ch.ka),_hp_q(ch.kb)
    ka>=0 && kb>=0 && ka*kb<1 || return false
    world(x,y)=((d*x-b*y)/det,(-c*x+a*y)/det)
    for s in (-ka,ka),t in (-kb,kb)
        u=world(1,s);v=world(t,1)
        _hp_cross(u,v)^2>=_hp_q(min_sine2)*_hp_dot(u,u)*_hp_dot(v,v) || return false
    end
    true
end
function _hp_cert_chart_linear_motion(ch,A,B,newA,newB;min_sine2=1//4)
    length(A)==length(newA) && length(B)==length(newB) || return false
    _hp_chart_angle_certificate(ch;min_sine2) && _hp_chart_graph_q(ch,A,B) && _hp_chart_graph_q(ch,newA,newB)
end

# Graph-ribbon triangle clipping is exact for the supplied binary coordinates.
# Allowed endpoints mean isolated contacts only, never overlapping intervals.
function _hp_segment_meets_triangle(a,b,c,p,q;allowed=())
    a,b,c,p,q=_hp_qpoint.((a,b,c,p,q));s=sign(_hp_cross(_hp_sub(b,a),_hp_sub(c,a)))
    s==0 && throw(ArgumentError("triangle must have nonzero area"))
    low,high=_hp_q(0),_hp_q(1)
    for (u,v) in ((a,b),(b,c),(c,a))
        f0=s*_hp_cross(_hp_sub(v,u),_hp_sub(p,u));f1=s*_hp_cross(_hp_sub(v,u),_hp_sub(q,u))
        slope=f1-f0
        if slope==0
            f0<0 && return false
        elseif slope>0
            low=max(low,-f0/slope)
        else
            high=min(high,-f0/slope)
        end
        low>high && return false
    end
    if low==high
        hit=_hp_lerp(p,q,low)
        any(x->_hp_qpoint(x)==hit,allowed) && return false
    end
    true
end
"""Certify deletion of one graph vertex; crossing partners require a separate chart certificate."""
function _hp_cert_graph_delete(a,b,c,forbidden;holes=())
    qa,qb,qc=_hp_qpoint.((a,b,c));u=_hp_sub(qc,qa);l2=_hp_dot(u,u)
    l2>0 || return false
    t=_hp_dot(_hp_sub(qb,qa),u)
    0<t<l2 || return false
    area=_hp_cross(u,_hp_sub(qb,qa))
    area==0 && return true # exact same oriented geometry; assumes valid input state
    for (p,q) in forbidden
        _hp_segment_meets_triangle(a,b,c,p,q;allowed=(a,c)) && return false
    end
    for (center,radius) in holes
        o=_hp_qpoint(center);r=_hp_q(radius)
        _hp_segment_meets_triangle(a,b,c,center,center) && return false
        minimum((_hp_point_segment_d2q(o,qa,qb),_hp_point_segment_d2q(o,qb,qc),_hp_point_segment_d2q(o,qc,qa)))>r*r || return false
    end
    true
end

function _hp_port_collar_certificate(p,q,c,r,dtheta,collar_path_bound)
    all(isfinite,(p...,q...,c...,r,dtheta,collar_path_bound)) || return false
    r>0 && collar_path_bound>=0 || return false
    # 2sin(|theta|/2) <= |theta|. The rational weaker bound avoids unverified libm rounding.
    radius=_hp_q(r);radial=_hp_sub(_hp_qpoint(p),_hp_qpoint(c))
    _hp_dot(radial,radial)==radius^2 || return false
    n=_hp_mul(radial,1/radius)
    e=_hp_sub(_hp_qpoint(q),_hp_qpoint(p));theta=abs(_hp_q(dtheta));dq=_hp_q(collar_path_bound)
    # |e| <= |e.x|+|e.y|, a conservative exact upper bound.
    _hp_dot(n,e)-theta*(abs(e[1])+abs(e[2]))-dq-_hp_q(r)*theta>0
end
