# cf: .research-scratch/PJ260922/261109-lickorish-drawing-mock/LIB-REQUIREMENTS.md §3.1–3.2, §3.5, §5 — signed crossing identities and finite geometry checks
struct HCPlotCrossing
    run_a::Int
    run_b::Int
    station_a::Float64
    station_b::Float64
    sign::Int
    angle::Float64
    point::HCPlotPoint
end
function _hp_segments(runs)
    segments=Tuple{Int,Int,HCPlotPoint,HCPlotPoint}[]
    for (ri,r) in enumerate(runs)
        n=length(r.points)
        for i in 1:n-1
            _hp_dist(r.points[i],r.points[i+1])>1e-12 && push!(segments,(ri,i,r.points[i],r.points[i+1]))
        end
        r.closed && _hp_dist(r.points[end],r.points[1])>1e-12 && push!(segments,(ri,n,r.points[end],r.points[1]))
    end
    segments
end
function _hp_grid_pairs(segments;cell=40.0)
    grid=Dict{Tuple{Int,Int},Vector{Int}}();pairs=Set{Tuple{Int,Int}}()
    for (i,(_,_,a,b)) in enumerate(segments)
        for x in floor(Int,min(a[1],b[1])/cell):floor(Int,max(a[1],b[1])/cell),
            y in floor(Int,min(a[2],b[2])/cell):floor(Int,max(a[2],b[2])/cell)
            bucket=get!(grid,(x,y),Int[])
            for j in bucket;push!(pairs,(j,i));end
            push!(bucket,i)
        end
    end
    sort!(collect(pairs))
end
function _hp_crossings(runs)
    segs=_hp_segments(runs);crossings=HCPlotCrossing[];overlaps=Tuple{Int,Int}[]
    for (i,j) in _hp_grid_pairs(segs)
        ra,ia,a,b=segs[i];rb,ib,c,d=segs[j]
        adjacent=ra==rb && (abs(ia-ib)==1 || (runs[ra].closed && Set((ia,ib))==Set((1,length(runs[ra].points)))))
        u=_hp_sub(b,a);v=_hp_sub(d,c);nu,nv=_hp_norm(u),_hp_norm(v)
        if abs(_hp_cross(u,v))<1e-10*nu*nv && abs(_hp_cross(u,_hp_sub(c,a)))<1e-9*nu
            ts=[_hp_dot(_hp_sub(p,a),u)/(nu*nu) for p in (c,d)]
            min(1,maximum(ts))-max(0,minimum(ts))>1e-8 && push!(overlaps,(ra,rb))
        end
        adjacent && continue
        hit=_hp_intersection(a,b,c,d;halfopen=true);isnothing(hit) && continue
        t,s,sign=hit
        angle=acosd(clamp(abs(_hp_dot(u,v))/(nu*nv),0,1))
        push!(crossings,HCPlotCrossing(ra,rb,ia-1+t,ib-1+s,sign,angle,_hp_lerp(a,b,t)))
    end
    crossings,overlaps
end

"""Match geometric crossings to input labels, allowing only cyclic starts on closed runs.

The bijection is shared by all runs. Equal partner/sign strings alone are insufficient.
"""
function _hp_match_words(runs,crossings)
    visits=[Tuple{Float64,Int,Int}[] for _ in runs]
    for (id,c) in enumerate(crossings)
        push!(visits[c.run_a],(c.station_a,id,c.sign));push!(visits[c.run_b],(c.station_b,id,-c.sign))
    end
    for v in visits;sort!(v;by=first);end
    any(length(visits[i])!=length(runs[i].word) for i in eachindex(runs)) && return false,Dict{Int,Int}()
    order=sort!(collect(eachindex(runs));by=i->(runs[i].closed,length(visits[i])))
    function match(at,forward,backward)
        at>length(order) && return forward
        ri=order[at];r=runs[ri];events=visits[ri];n=length(events)
        shifts=r.closed && n>0 ? (0:n-1) : (0:0)
        for shift in shifts
            f=copy(forward);b=copy(backward);ok=true
            for j in 1:n
                _,id,sign=events[j];label=r.word[mod1(j+shift,n)]
                if sign!=Base.sign(label) || get(f,id,abs(label))!=abs(label) || get(b,abs(label),id)!=id
                    ok=false;break
                end
                f[id]=abs(label);b[abs(label)]=id
            end
            ok || continue
            result=match(at+1,f,b);isnothing(result) || return result
        end
        nothing
    end
    mapping=match(1,Dict{Int,Int}(),Dict{Int,Int}())
    isnothing(mapping) ? (false,Dict{Int,Int}()) : (true,mapping)
end
# cf: three-examples-261111/angles.py — three-segment tangent used in the owner's prototype comparison
function _hp_stencil_angle(runs,c)
    function tangent(ri,station)
        r=runs[ri];n=length(r.points);i=floor(Int,station)+1
        a=r.closed ? mod1(i-1,n) : max(i-1,1)
        b=r.closed ? mod1(i+2,n) : min(i+2,n)
        v=_hp_sub(r.points[b],r.points[a]);_hp_mul(v,1/max(_hp_norm(v),1e-12))
    end
    a,b=tangent(c.run_a,c.station_a),tangent(c.run_b,c.station_b)
    # A zero stencil has no tangent direction; never report it as a right angle.
    min(_hp_norm(a),_hp_norm(b))<.5 && return 0.0
    acosd(clamp(abs(_hp_dot(a,b)),0,1))
end

# cf: LIB-REQUIREMENTS.md §5 — independent before/after check of port cyclic order
function _hp_port_snapshot(p)
    byhole=Dict(name=>Tuple{Tuple{String,Symbol},Float64}[] for (name,_) in p.holes)
    offcircle=0
    for r in p.runs
        r.closed && continue
        for (side,pt) in ((:start,first(r.points)),(:end,last(r.points)))
            isempty(p.holes) && (offcircle+=1;continue)
            gaps=[abs(_hp_dist(pt,c)-radius) for (_,(c,radius)) in p.holes]
            hi=argmin(gaps);gaps[hi]>1e-5 && (offcircle+=1)
            name,(c,_)=p.holes[hi];v=_hp_sub(pt,c)
            push!(byhole[name],((r.name,side),mod(atan(v[2],v[1]),2pi)))
        end
    end
    for rows in values(byhole);sort!(rows;by=last);end
    byhole,offcircle
end
function _hp_port_audit(before,after;pairs=_hp_ports(before))
    before.square && return _hp_square_endpoint_audit(before,after;pairs)
    old,oldoff=_hp_port_snapshot(before);new,newoff=_hp_port_snapshot(after)
    cyclic(a,b)=length(a)==length(b) && (isempty(a) || any(circshift(a,k)==b for k in 0:length(a)-1))
    order=Set(keys(old))==Set(keys(new)) && all(cyclic(first.(old[k]),first.(new[k])) for k in keys(old))
    gapmap(rows)=Dict(rows[i][1]=>mod(rows[mod1(i+1,length(rows))][2]-rows[i][2],2pi) for i in eachindex(rows))
    gaps=Float64[];nondecrease=true
    for k in keys(old)
        length(old[k])<2 && continue
        a,b=gapmap(old[k]),gapmap(new[k]);append!(gaps,values(b))
        nondecrease &= all(get(b,id,-Inf)>=min(gap,deg2rad(18.0))-1e-10 for (id,gap) in a)
    end
    mirror=0.0;final=Dict(r.name=>r for r in after.runs)
    for pair in pairs
        coordinates=HCPlotPoint[]
        for (ri,i,_,_) in pair
            r=before.runs[ri];fr=final[r.name]
            push!(coordinates,i==1 ? first(fr.points) : last(fr.points))
        end
        a,b=coordinates;mirror=max(mirror,_hp_dist((a[1],before.height-a[2]),b))
    end
    (cyclic_order=order,offcircle=newoff,mirror_error=mirror,
     minimum_gap_degrees=isempty(gaps) ? nothing : rad2deg(minimum(gaps)),
     minimum_spacing=isempty(gaps) || minimum(gaps)>=deg2rad(18.0)-1e-10,
     existing_small_gaps_not_worsened=nondecrease,
     ok=order && oldoff==0 && newoff==0 && mirror<1e-5 && nondecrease)
end
function _hp_audit(placement::HCPlotPlacement)
    runs=placement.runs;crossings,overlaps=_hp_crossings(runs)
    structure,mapping=_hp_match_words(runs,crossings)
    same_family=count(c->runs[c.run_a].role==runs[c.run_b].role && runs[c.run_a].role!=:none,crossings)
    hole_entries=Tuple{String,String,Int}[]
    for r in runs
        for i in 1:(r.closed ? length(r.points) : length(r.points)-1)
            a=r.points[i];b=r.points[mod1(i+1,length(r.points))]
            for (name,(center,radius)) in placement.holes
                _hp_nearest(center,a,b)[1]<radius-1e-6 && push!(hole_entries,(r.name,name,i))
            end
        end
    end
    oc,radius=_hp_frame(placement)
    frame_entries=[(r.name,i) for r in runs for i in eachindex(r.points) if
        placement.square ? !( -1e-6<=r.points[i][1]<=placement.width+1e-6 && -1e-6<=r.points[i][2]<=placement.height+1e-6) : _hp_dist(r.points[i],oc)>radius+1e-6]
    stencil_angles=[_hp_stencil_angle(runs,c) for c in crossings]
    (structure=structure,crossings=crossings,labels=mapping,stencil_angles=stencil_angles,
     minimum_stencil_angle=isempty(stencil_angles) ? nothing : minimum(stencil_angles),
     minimum_angle=isempty(crossings) ? nothing : minimum(c.angle for c in crossings),
     same_family=same_family,overlaps=overlaps,hole_entries=hole_entries,frame_entries=frame_entries,
     ok=structure && same_family==0 && isempty(overlaps) && isempty(hole_entries) && isempty(frame_entries))
end
