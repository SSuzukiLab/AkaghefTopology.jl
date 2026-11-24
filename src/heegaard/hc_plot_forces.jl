# cf: .research-scratch/PJ260922/astra-lanes/hc-plot/answer-3.md §9 — reusable grid and force workspaces; force law remains LIB-REQUIREMENTS.md §3
# cf: lib/topology.jl/docs/specs/heegaard-plot.md 「二角形」 — bigon corners: narrow crossing window and local dilation
const _HP_BIGON_WINDOW=10.  # arc-length window of a bigon corner; other crossings use 30
const _HP_BIGON_REACH=130.  # a bigon is enlarged while its corners are closer than this
const _HP_BIGON_SPREAD=15.  # arc-length scale over which the push fades beyond the corners
mutable struct HCPlotForceWorkspace
    cell::Float64
    origin::HCPlotPoint
    nx::Int
    ny::Int
    counts::Vector{Int}
    offsets::Vector{Int}
    cursor::Vector{Int}
    entries::Vector{Int}
    neighbors::Vector{Vector{Int}}
    valid::BitVector
    seen::Vector{Int}
    segments::Vector{Tuple{Int,Int,HCPlotPoint,HCPlotPoint}}
    vectors::Vector{HCPlotPoint}
    length2::Vector{Float64}
    lengths::Vector{Float64}
    stations::Vector{Vector{Float64}}
    totals::Vector{Float64}
    events::Matrix{Vector{Tuple{Float64,Float64}}}
    window2::Matrix{Vector{Float64}} # parallel to events: squared arc-length window of the crossing
    bigons::Vector{Tuple{Int,Int}}
    ghost::BitVector # parallel to segments: a translated copy across the identified edges of the square
    component_id::Vector{Int}
    component_offset::Vector{Float64}
    component_total::Vector{Float64}
    event_source::Vector{Vector{Float64}}
    force::Vector{Vector{HCPlotPoint}}
    cap::Vector{Vector{Float64}}
    prox::Vector{Vector{Float64}}
end
function HCPlotForceWorkspace(p;cell=55.)
    center,r=_hp_frame(p);nx=ceil(Int,2r/cell)+3;ny=nx
    n=length(p.runs)
    HCPlotForceWorkspace(cell,(center[1]-r-cell,center[2]-r-cell),nx,ny,
        zeros(Int,nx*ny),zeros(Int,nx*ny+1),zeros(Int,nx*ny),Int[],[Int[] for _ in 1:nx*ny],falses(nx*ny),Int[],
        Tuple{Int,Int,HCPlotPoint,HCPlotPoint}[],HCPlotPoint[],Float64[],Float64[],[Float64[] for _ in 1:n],zeros(n),
        [Tuple{Float64,Float64}[] for _ in 1:n,_ in 1:n],[Float64[] for _ in 1:n,_ in 1:n],Tuple{Int,Int}[],falses(0),zeros(Int,n),zeros(n),zeros(n),[Float64[] for _ in 1:n],
        [HCPlotPoint[] for _ in 1:n],[Float64[] for _ in 1:n],[Float64[] for _ in 1:n])
end
_hp_fastnorm(v)=sqrt(v[1]*v[1]+v[2]*v[2])
_hp_force_cell(ws,p)=(clamp(floor(Int,(p[1]-ws.origin[1])/ws.cell)+1,1,ws.nx),clamp(floor(Int,(p[2]-ws.origin[2])/ws.cell)+1,1,ws.ny))
_hp_force_index(ws,x,y)=x+(y-1)*ws.nx
function _hp_force_prepare!(ws,p,crossings)
    empty!(ws.segments);empty!(ws.vectors);empty!(ws.length2);empty!(ws.lengths)
    @inbounds for (ri,r) in enumerate(p.runs)
        n=length(r.points);resize!(ws.stations[ri],n);ws.stations[ri][1]=0.
        resize!(ws.force[ri],n);fill!(ws.force[ri],(0.,0.))
        resize!(ws.cap[ri],n);fill!(ws.cap[ri],Inf)
        resize!(ws.prox[ri],n);fill!(ws.prox[ri],Inf)
        for i in 1:(r.closed ? n : n-1)
            a,b=r.points[i],r.points[mod1(i+1,n)];v=_hp_sub(b,a);l2=_hp_dot(v,v);len=sqrt(l2)
            len>1e-12 || throw(ArgumentError("zero length force segment"))
            push!(ws.segments,(ri,i,a,b));push!(ws.vectors,v);push!(ws.length2,l2);push!(ws.lengths,len)
            i<n && (ws.stations[ri][i+1]=ws.stations[ri][i]+len)
        end
        ws.totals[ri]=ws.stations[ri][end]+(r.closed ? ws.lengths[end] : 0.)
    end
    resize!(ws.ghost,length(ws.segments));fill!(ws.ghost,false)
    p.square && _hp_square_ghosts!(ws,p)
    for v in ws.events;empty!(v);end
    for v in ws.window2;empty!(v);end
    _hp_bigon_pairs!(ws.bigons,p.runs,crossings)
    @inbounds for (id,c) in enumerate(crossings)
        function at(ri,t)
            r=p.runs[ri]
            # A half-open crossing on the closing segment can round to n.
            # Its station is the same point as zero on a closed run.
            r.closed && (t=mod(t,length(r.points)))
            i=floor(Int,t)+1;u=t-(i-1)
            ws.stations[ri][i]+u*_hp_fastnorm(_hp_sub(r.points[mod1(i+1,length(r.points))],r.points[i]))
        end
        a,b=at(c.run_a,c.station_a),at(c.run_b,c.station_b)
        push!(ws.events[c.run_a,c.run_b],(a,b));push!(ws.events[c.run_b,c.run_a],(b,a))
        window2=any(id in pair for pair in ws.bigons) ? _HP_BIGON_WINDOW^2 : 900.
        push!(ws.window2[c.run_a,c.run_b],window2);push!(ws.window2[c.run_b,c.run_a],window2)
    end
    _hp_grid_rebuild!(ws)
end
# On the genus 1 square the opposite edges are identified, so a strand near an edge has neighbours
# across it. Without them a family of parallel strands is pushed apart towards the corners and bows.
# The translated copies of all segments act as repellers; they are never moved themselves. The runs
# of one component follow each other across the edges, so arc length is measured along the component:
# the copy of the preceding run ends where a run starts and is its continuation, not a neighbour.
function _hp_square_ghosts!(ws,p)
    ids=Dict{String,Int}();offset=Dict{String,Float64}()
    for (ri,r) in enumerate(p.runs)
        ws.component_id[ri]=get!(ids,r.component,length(ids)+1)
        ws.component_offset[ri]=get(offset,r.component,0.)
        offset[r.component]=ws.component_offset[ri]+ws.totals[ri]
    end
    for (ri,r) in enumerate(p.runs);ws.component_total[ri]=offset[r.component];end
    real=length(ws.segments)
    for dx in (-p.width,0.,p.width),dy in (-p.height,0.,p.height)
        dx==0 && dy==0 && continue
        for j in 1:real
            ri,i,a,b=ws.segments[j];a=(a[1]+dx,a[2]+dy);b=(b[1]+dx,b[2]+dy)
            (max(a[1],b[1])<-220 || min(a[1],b[1])>p.width+220 ||
             max(a[2],b[2])<-220 || min(a[2],b[2])>p.height+220) && continue
            push!(ws.segments,(ri,i,a,b));push!(ws.vectors,ws.vectors[j])
            push!(ws.length2,ws.length2[j]);push!(ws.lengths,ws.lengths[j]);push!(ws.ghost,true)
        end
    end
end
# Two crossings of the same two runs that are neighbours on both runs are the
# corners of a bigon. Faces with more corners have two sides of one family,
# which repel at every distance; a bigon has none, so tension alone closes it.
function _hp_bigon_pairs!(pairs,runs,crossings)
    empty!(pairs);length(crossings)<2 && return pairs
    orders=_hp_gate_orders(runs,crossings)
    function neighbours(ri,a,b)
        o=orders[ri];n=length(o);i=findfirst(==(a),o);j=findfirst(==(b),o)
        abs(i-j)==1 || (runs[ri].closed && abs(i-j)==n-1)
    end
    for a in eachindex(crossings),b in a+1:length(crossings)
        ca,cb=crossings[a],crossings[b]
        minmax(ca.run_a,ca.run_b)==minmax(cb.run_a,cb.run_b) || continue
        neighbours(ca.run_a,a,b) && neighbours(ca.run_b,a,b) && push!(pairs,(a,b))
    end
    pairs
end
_hp_bigon_pairs(runs,crossings)=_hp_bigon_pairs!(Tuple{Int,Int}[],runs,crossings)
# A small bigon is enlarged by a local dilation about its midpoint. A dilation
# maps the two arcs and their two crossings together, so it separates the
# corners without changing the crossing structure, also when the bigon has
# collapsed to less than one sample spacing. On the two arcs of the bigon the
# push is the dilation; beyond the corners it points away from the midpoint and
# fades with arc length. Only the two runs of the bigon are moved.
function _hp_bigon_push!(ws,p,crossings)
    isempty(ws.bigons) && return
    function station(ri,c)
        r=p.runs[ri];n=length(r.points);t=c.run_a==ri ? c.station_a : c.station_b
        r.closed && (t=mod(t,n))
        i=floor(Int,t)+1;u=t-(i-1)
        ws.stations[ri][i]+u*_hp_fastnorm(_hp_sub(r.points[mod1(i+1,n)],r.points[i]))
    end
    for (a,b) in ws.bigons
        ca,cb=crossings[a],crossings[b];d=_hp_fastnorm(_hp_sub(ca.point,cb.point))
        d<_HP_BIGON_REACH || continue
        mid=_hp_mul(_hp_add(ca.point,cb.point),.5);half=max(.5d,1e-9)
        strength=(_HP_BIGON_REACH/max(d,3.))^2-1
        for ri in (ca.run_a,ca.run_b)
            r=p.runs[ri];n=length(r.points);total=ws.totals[ri]
            lo,hi=minmax(station(ri,ca),station(ri,cb))
            # The arc of the bigon is the one between the corners that carries no other crossing.
            inner=true
            if r.closed
                count=0;between=false
                for rj in eachindex(p.runs),(s,_) in ws.events[ri,rj]
                    count+=1;between |= lo<s<hi
                end
                inner=count==2 ? 2(hi-lo)<=total : !between
            end
            for k in 1:n
                !r.closed && k in (1,n) && continue
                s=ws.stations[ri][k]
                outside=if inner
                    lo<=s<=hi ? 0. : r.closed ? min(min(abs(s-lo),total-abs(s-lo)),min(abs(s-hi),total-abs(s-hi))) : min(abs(s-lo),abs(s-hi))
                else
                    (s<=lo || s>=hi) ? 0. : min(s-lo,hi-s)
                end
                outside<2.5_HP_BIGON_SPREAD || continue
                rel=_hp_sub(r.points[k],mid)
                weight=strength*exp(-(outside/_HP_BIGON_SPREAD)^2)/max(_hp_fastnorm(rel),half)
                ws.force[ri][k]=_hp_add(ws.force[ri][k],_hp_mul(rel,weight))
            end
        end
    end
end
function _hp_grid_rebuild!(ws)
    fill!(ws.counts,0);fill!(ws.valid,false)
    @inbounds for (_,_,a,b) in ws.segments
        ax,ay=_hp_force_cell(ws,a);bx,by=_hp_force_cell(ws,b)
        for y in min(ay,by):max(ay,by),x in min(ax,bx):max(ax,bx);ws.counts[_hp_force_index(ws,x,y)]+=1;end
    end
    ws.offsets[1]=1
    @inbounds for i in eachindex(ws.counts);ws.offsets[i+1]=ws.offsets[i]+ws.counts[i];ws.cursor[i]=ws.offsets[i];end
    resize!(ws.entries,ws.offsets[end]-1)
    @inbounds for (j,(_,_,a,b)) in enumerate(ws.segments)
        ax,ay=_hp_force_cell(ws,a);bx,by=_hp_force_cell(ws,b)
        for y in min(ay,by):max(ay,by),x in min(ax,bx):max(ax,bx)
            k=_hp_force_index(ws,x,y);ws.entries[ws.cursor[k]]=j;ws.cursor[k]+=1
        end
    end
    resize!(ws.seen,length(ws.segments));fill!(ws.seen,0)
end
function _hp_force_candidates!(ws,pt)
    gx,gy=_hp_force_cell(ws,pt);key=_hp_force_index(ws,gx,gy)
    ws.valid[key] && return ws.neighbors[key]
    dest=ws.neighbors[key];empty!(dest)
    reach=ceil(Int,220/ws.cell)
    @inbounds for y in max(1,gy-reach):min(ws.ny,gy+reach),x in max(1,gx-reach):min(ws.nx,gx+reach)
        dx=max(abs(x-gx)-1,0)*ws.cell;dy=max(abs(y-gy)-1,0)*ws.cell
        dx*dx+dy*dy>220.0^2 && continue
        k=_hp_force_index(ws,x,y)
        for at in ws.offsets[k]:ws.offsets[k+1]-1
            j=ws.entries[at]
            ws.seen[j]=key
        end
    end
    # The stamp scan retains the prototype's segment order without sorting
    # and allocating radix-sort scratch for every occupied query cell.
    for j in eachindex(ws.seen);ws.seen[j]==key && push!(dest,j);end
    ws.valid[key]=true;dest
end
function _hp_force_state(p,crossings,partners,ws::HCPlotForceWorkspace=HCPlotForceWorkspace(p))
    _hp_force_prepare!(ws,p,crossings)
    runs=p.runs;stations=ws.stations;totals=ws.totals
    force,cap,prox=ws.force,ws.cap,ws.prox;oc,radius=_hp_frame(p)
    cyclic(d,ri)=runs[ri].closed ? min(abs(d),totals[ri]-abs(d)) : abs(d)
    function hole_of(r,k)
        best=Inf;result=0
        for (i,(_, (c,rad))) in enumerate(p.holes)
            gap=abs(_hp_fastnorm(_hp_sub(r.points[k],c))-rad)
            if gap<best;best=gap;result=i;end
        end
        result
    end
    @inbounds for (ri,r) in enumerate(runs)
        n=length(r.points);ports=r.closed ? (0,0) : (hole_of(r,1),hole_of(r,n))
        for i in 1:n
            pt=r.points[i];prev=r.points[mod1(i-1,n)];nxt=r.points[mod1(i+1,n)];fixed=!r.closed && i in (1,n)
            # The source station depends on this point and crossing, not on
            # the candidate segment. Retain the exact arithmetic and event
            # order while reusing this term throughout the neighbor query.
            for rj in eachindex(runs)
                events=ws.events[ri,rj];buffer=ws.event_source[rj];resize!(buffer,length(events))
                for k in eachindex(events);buffer[k]=cyclic(stations[ri][i]-events[k][1],ri)^2;end
            end
            tangent=fixed ? (i==1 ? _hp_sub(nxt,pt) : _hp_sub(pt,prev)) : _hp_sub(nxt,prev)
            tangent=_hp_mul(tangent,1/max(_hp_fastnorm(tangent),1e-12))
            f=fixed ? (0.,0.) : _hp_mul(_hp_sub(_hp_add(prev,nxt),_hp_mul(pt,2)),.4)
            for j in _hp_force_candidates!(ws,pt)
                rj,ij,a,b=ws.segments[j];seglen=ws.lengths[j]
                idxdist=abs(i-ij-.5);r.closed && ri==rj && (idxdist=min(idxdist,n-idxdist))
                ghost=ws.ghost[j]
                !ghost && ri==rj && idxdist<=2.5 && continue
                u=ws.vectors[j];v=_hp_sub(pt,a);t=clamp(_hp_dot(v,u)/ws.length2[j],0.,1.)
                rel=_hp_sub(v,_hp_mul(u,t));d2=_hp_dot(rel,rel);d2>220.0^2 && continue
                dist=sqrt(d2);smid=stations[rj][ij]+.5seglen
                # Prototype adaptive_spacing excludes index neighbours only;
                # arc-length continuation is an exclusion for force/caps only.
                if ghost && ws.component_id[ri]==ws.component_id[rj]
                    along=abs(ws.component_offset[ri]+stations[ri][i]-ws.component_offset[rj]-smid)
                    min(along,ws.component_total[ri]-along)<3dist && continue
                end
                prox[ri][i]=min(prox[ri][i],dist)
                !ghost && ri==rj && cyclic(stations[ri][i]-smid,ri)<3dist && continue
                partner=partners[ri,rj];active=!partner || dist<=28
                fade=active ? 1. : 0.;near=false
                window2=ws.window2[ri,rj]
                for (k,(_,s2)) in enumerate(ws.events[ri,rj])
                    source2=ws.event_source[rj][k]
                    # With zero force weight only the window predicate is
                    # needed. A nonnegative target square cannot lower it.
                    !active && source2>=900 && continue
                    eventd2=source2+cyclic(smid-s2,rj)^2
                    active && (fade*=1-exp(-eventd2/window2[k]));near |= eventd2<900
                    !active && near && break
                end
                if partner && active
                    align=_hp_dot(tangent,_hp_mul(u,1/seglen))^2
                    fade*=align
                end
                if dist>1e-12 && fade>0
                    weight=3*(14/max(dist,3))^2/dist*fade*(seglen/40)
                    f=_hp_add(f,_hp_mul(rel,weight))
                end
                if !near
                    cap[ri][i]=min(cap[ri][i],.3dist);cap[rj][ij]=min(cap[rj][ij],.3dist)
                    jj=mod1(ij+1,length(runs[rj].points));cap[rj][jj]=min(cap[rj][jj],.3dist)
                end
            end
            for (hi,(_, (c,rad))) in enumerate(p.holes)
                rel=_hp_sub(pt,c);rr=_hp_fastnorm(rel);gap=rr-rad;fade=1.
                if !r.closed
                    hi==ports[1] && (fade*=1-exp(-(stations[ri][i]/45)^2))
                    hi==ports[2] && (fade*=1-exp(-((totals[ri]-stations[ri][i])/45)^2))
                end
                gap<220 && (f=_hp_add(f,_hp_mul(rel,6*(14/max(gap,3))^2/max(rr,1e-12)*fade)))
                cap[ri][i]=min(cap[ri][i],.3max(gap,1.));prox[ri][i]=min(prox[ri][i],gap)
            end
            rr=_hp_fastnorm(_hp_sub(pt,oc));gap=radius-rr
            gap<220 && (f=_hp_sub(f,_hp_mul(_hp_sub(pt,oc),6*(14/max(gap,3))^2/max(rr,1e-12))))
            cap[ri][i]=min(cap[ri][i],.3max(gap,1.));force[ri][i]=fixed ? (0.,0.) : f
        end
    end
    _hp_bigon_push!(ws,p,crossings)
    force,cap,prox
end
