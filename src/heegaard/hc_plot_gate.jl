# cf: .research-scratch/PJ260922/261109-lickorish-drawing-mock/three-examples-261111/relax23.py — intersections, valid, offenders, remesh
# lib/topology.jl/docs/specs/heegaard-plot.md — local proposal checks and one final input Gauss comparison
mutable struct HCPlotGateWorkspace
    grid::HCPlotForceWorkspace
    seen::Vector{UInt32}
    epoch::UInt32
end
HCPlotGateWorkspace(p)=HCPlotGateWorkspace(HCPlotForceWorkspace(p;cell=40.),UInt32[],0)
struct HCPlotGateState
    crossings::Vector{HCPlotCrossing} # persistent local IDs index this vector
    orders::Vector{Vector{Int}}
    pairs::Matrix{Vector{Int}}
end
function _hp_gate_orders(runs,crossings,ids=eachindex(crossings))
    visits=[Tuple{Float64,Int}[] for _ in runs]
    for (c,id) in zip(crossings,ids)
        push!(visits[c.run_a],(c.station_a,id));push!(visits[c.run_b],(c.station_b,id))
    end
    [last.(sort!(v;by=first)) for v in visits]
end
function HCPlotGateState(runs,crossings)
    n=length(runs);pairs=[Int[] for _ in 1:n,_ in 1:n]
    for (id,c) in enumerate(crossings);push!(pairs[c.run_a,c.run_b],id);end
    HCPlotGateState(crossings,_hp_gate_orders(runs,crossings),pairs)
end
function _hp_same_event_order(a,b,closed)
    length(a)==length(b) || return false
    a==b && return true
    closed && !isempty(a) || return false
    start=findfirst(==(a[1]),b);isnothing(start) && return false
    all(a[i]==b[mod1(start+i-1,length(b))] for i in eachindex(a))
end
function _hp_gate_geometry(p,workspace::HCPlotGateWorkspace)
    ws=workspace.grid;empty!(ws.segments);held=Tuple{Int,Int}[]
    for (ri,r) in enumerate(p.runs),i in 1:(r.closed ? length(r.points) : length(r.points)-1)
        a,b=r.points[i],r.points[mod1(i+1,length(r.points))]
        if _hp_dot(_hp_sub(b,a),_hp_sub(b,a))<=1e-24
            push!(held,(ri,i),(ri,mod1(i+1,length(r.points))));continue
        end
        push!(ws.segments,(ri,i,a,b))
    end
    _hp_grid_rebuild!(ws);m=length(ws.segments);nseen=m*(m-1)÷2
    if length(workspace.seen)<nseen
        old=length(workspace.seen);resize!(workspace.seen,nseen);fill!(view(workspace.seen,old+1:nseen),0)
    end
    if workspace.epoch==typemax(UInt32);fill!(workspace.seen,0);workspace.epoch=0;end
    workspace.epoch+=1;epoch=workspace.epoch
    crossings=HCPlotCrossing[];overlaps=0;same_family=0
    function hold(ri,i)
        push!(held,(ri,i),(ri,mod1(i+1,length(p.runs[ri].points))))
    end
    @inbounds for cell in eachindex(ws.counts)
        lo,hi=ws.offsets[cell],ws.offsets[cell+1]-1
        for bpos in lo+1:hi,apos in lo:bpos-1
            x,y=ws.entries[apos],ws.entries[bpos];i,j=minmax(x,y)
            key=(j-1)*(j-2)÷2+i
            workspace.seen[key]==epoch && continue
            workspace.seen[key]=epoch
            ra,ia,a,b=ws.segments[i];rb,ib,c,d=ws.segments[j]
            max(a[1],b[1])<min(c[1],d[1]) && continue
            max(c[1],d[1])<min(a[1],b[1]) && continue
            max(a[2],b[2])<min(c[2],d[2]) && continue
            max(c[2],d[2])<min(a[2],b[2]) && continue
            u=_hp_sub(b,a);v=_hp_sub(d,c);l2=_hp_dot(u,u);k2=_hp_dot(v,v);den=_hp_cross(u,v)
            if den*den<1e-20*l2*k2 && _hp_cross(u,_hp_sub(c,a))^2<1e-18*l2
                t1=_hp_dot(_hp_sub(c,a),u)/l2;t2=_hp_dot(_hp_sub(d,a),u)/l2
                if min(1.,max(t1,t2))-max(0.,min(t1,t2))>1e-8
                    overlaps+=1;hold(ra,ia);hold(rb,ib)
                end
            end
            adjacent=ra==rb && (abs(ia-ib)==1 || (p.runs[ra].closed && min(ia,ib)==1 && max(ia,ib)==length(p.runs[ra].points)))
            adjacent && continue
            abs(den)<1e-12 && continue
            q=_hp_sub(c,a);t=_hp_cross(q,v)/den;s=_hp_cross(q,u)/den
            0<=t<1 && 0<=s<1 || continue
            push!(crossings,HCPlotCrossing(ra,rb,ia-1+t,ib-1+s,Int(sign(den)),0.,_hp_lerp(a,b,t)))
            if ra==rb || (p.runs[ra].role==p.runs[rb].role && p.runs[ra].role!=:none)
                same_family+=1;hold(ra,ia);hold(rb,ib)
            end
        end
    end
    holes=0;frame=0;center,radius=_hp_frame(p);frame2=(radius+1e-6)^2
    @inbounds for (ri,i,a,b) in ws.segments
        v=_hp_sub(b,a);l2=_hp_dot(v,v)
        for (_,(c,r)) in p.holes
            min(a[1],b[1])>c[1]+r && continue
            max(a[1],b[1])<c[1]-r && continue
            min(a[2],b[2])>c[2]+r && continue
            max(a[2],b[2])<c[2]-r && continue
            t=clamp(_hp_dot(_hp_sub(c,a),v)/l2,0.,1.);q=_hp_sub(c,_hp_lerp(a,b,t))
            if _hp_dot(q,q)<(r-1e-6)^2;holes+=1;hold(ri,i);end
        end
    end
    @inbounds for (ri,r) in enumerate(p.runs),i in eachindex(r.points)
        q=r.points[i];d=_hp_sub(q,center)
        outside=p.square ? !( -1e-6<=q[1]<=p.width+1e-6 && -1e-6<=q[2]<=p.height+1e-6) : _hp_dot(d,d)>frame2
        if outside;frame+=1;push!(held,(ri,i));end
    end
    (crossings=crossings,held=held,overlaps=overlaps,same_family=same_family,holes=holes,frame=frame,
     ok=isempty(held) && overlaps==0 && same_family==0 && holes==0 && frame==0)
end
function _hp_gate_trial(p,geometry,old::HCPlotGateState)
    n=length(p.runs);new=geometry.crossings;newpairs=[Int[] for _ in 1:n,_ in 1:n]
    for (id,c) in enumerate(new);push!(newpairs[c.run_a,c.run_b],id);end
    ids=zeros(Int,length(new));used=falses(length(old.crossings));held=copy(geometry.held);structure=true
    function hold(c)
        for (ri,t) in ((c.run_a,c.station_a),(c.run_b,c.station_b))
            i=floor(Int,t)+1;push!(held,(ri,i),(ri,mod1(i+1,length(p.runs[ri].points))))
        end
    end
    for ra in 1:n,rb in ra:n
        before,after=old.pairs[ra,rb],newpairs[ra,rb]
        candidates=Tuple{Float64,Int,Int}[]
        for a in before,b in after
            d=_hp_sub(old.crossings[a].point,new[b].point);push!(candidates,(_hp_dot(d,d),a,b))
        end
        sort!(candidates)
        for (d,a,b) in candidates
            !used[a] && ids[b]==0 && d<60.0^2 || continue
            ids[b]=a;used[a]=true
            if old.crossings[a].sign!=new[b].sign;structure=false;hold(new[b]);end
        end
        if length(before)!=length(after);structure=false;end
        for a in before
            if !used[a]
                structure=false
                # Old and candidate have the same segment indices for a motion.
                # Remeshing is rejected as a whole and never consumes these indices.
                hold(old.crossings[a])
            end
        end
        for b in after;if ids[b]==0;structure=false;hold(new[b]);end;end
    end
    state=old
    if structure
        ordered=similar(old.crossings)
        for (c,id) in zip(new,ids);ordered[id]=c;end
        # Each retained pair of neighbouring events must remain in order.
        # No candidate Gauss word or input-label bijection is built here.
        for ri in 1:n
            neighbours=old.orders[ri];length(neighbours)<2 && continue
            station(id)=ordered[id].run_a==ri ? ordered[id].station_a : ordered[id].station_b
            anchor=station(neighbours[1]);previous=0.
            for id in neighbours[2:end]
                current=station(id)-anchor
                p.runs[ri].closed && (current=mod(current,length(p.runs[ri].points)))
                if current<=previous;structure=false;break;end
                previous=current
            end
        end
        structure && (state=HCPlotGateState(ordered,old.orders,old.pairs))
    end
    (ok=geometry.ok && structure,structure=structure,state=state,held=held,
     overlaps=geometry.overlaps,same_family=geometry.same_family,holes=geometry.holes,frame=geometry.frame)
end
