# cf: .research-scratch/PJ260922/261109-lickorish-drawing-mock/LIB-REQUIREMENTS.md §3 — independent runs, local partner repulsion, shared sliding ports and adaptive sampling
# Production and DEBUG_MODE share locally checked proposals. DEBUG_MODE also
# runs the independent input-word audit after accepted motions and remeshes.
function _hp_rebuild(p,runs)
    HCPlotPlacement(p.width,p.height,p.holes,runs,p.methods,p.center,p.delta,p.square)
end
function _hp_sample(points,spacing,closed)
    P=copy(points);sp=copy(spacing)
    if closed
        if _hp_dist(P[1],P[end])<1e-8;pop!(P);pop!(sp);end
        push!(P,P[1]);push!(sp,sp[1])
    end
    len=[_hp_dist(P[i],P[i+1]) for i in 1:length(P)-1];s=vcat(0.0,cumsum(len));total=s[end]
    total>1e-8 || throw(ArgumentError("zero-length run"))
    out=HCPlotPoint[P[1]];t=0.0;i=1
    while true
        while i<length(len) && s[i+1]<=t;i+=1;end
        u=len[i]>0 ? (t-s[i])/len[i] : 0.0
        step=min(total/3,(1-u)*sp[i]+u*sp[i+1]);t+=step
        t>=total-.5step && break
        while i<length(len) && s[i+1]<t;i+=1;end
        push!(out,_hp_lerp(P[i],P[i+1],(t-s[i])/len[i]))
    end
    closed || push!(out,P[end])
    out
end
function _hp_run(r,points)
    HCPlotRun(r.name,r.component,r.role,r.closed,points,r.word)
end
function _hp_stations(r)
    vcat(0.0,cumsum([_hp_dist(r.points[i],r.points[i+1]) for i in 1:length(r.points)-1]))
end
function _hp_ports(p)
    ends=Tuple{Int,Int,Int,Int}[]
    isempty(p.holes) && return Tuple{NTuple{4,Int},NTuple{4,Int}}[]
    for (ri,r) in enumerate(p.runs)
        r.closed && continue;n=length(r.points)
        for (i,nb) in ((1,2),(n,n-1))
            ds=[abs(_hp_dist(r.points[i],c)-rad) for (_, (c,rad)) in p.holes]
            minimum(ds)<1e-5 || throw(ArgumentError("open endpoint is not on a hole"))
            push!(ends,(ri,i,nb,argmin(ds)))
        end
    end
    used=Set{Int}();pairs=Tuple{NTuple{4,Int},NTuple{4,Int}}[]
    for a in eachindex(ends)
        a in used && continue
        ri,i,_,_=ends[a];pt=p.runs[ri].points[i];mirror=(pt[1],p.height-pt[2])
        candidates=[b for b in eachindex(ends) if b!=a && !(b in used)]
        ds=[_hp_dist(mirror,p.runs[ends[b][1]].points[ends[b][2]]) for b in candidates]
        isempty(ds) && throw(ArgumentError("unpaired port"))
        minimum(ds)<1e-5 || throw(ArgumentError("unpaired mirrored port"))
        b=candidates[argmin(ds)];push!(pairs,(ends[a],ends[b]));push!(used,a,b)
    end
    pairs
end
# cf: answer-3.md §8 — negative gradient of squared cyclic gaps on the whole circle
_hp_port_spread(previous_gap,next_gap)=(next_gap-previous_gap)/pi
function _hp_port_motion(p,pairs)
    motions=[fill((0.0,0.0),length(r.points)) for r in p.runs]
    ends=[e for pair in pairs for e in pair];angles=Dict{Tuple{Int,Int},Float64}()
    for (ri,i,_,hi) in ends
        c,_=p.holes[hi].second;rel=_hp_sub(p.runs[ri].points[i],c);angles[(ri,i)]=atan(rel[2],rel[1])
    end
    function gaps(ri,i,hi)
        a=angles[(ri,i)];others=[mod(angles[(r,j)]-a,2pi) for (r,j,_,h) in ends if h==hi && (r,j)!=(ri,i)]
        isempty(others) ? (2pi,2pi) : (2pi-maximum(others),minimum(others))
    end
    function pull(e)
        ri,i,nb,_=e;a=angles[(ri,i)];rel=_hp_sub(p.runs[ri].points[nb],p.runs[ri].points[i])
        _hp_dot(rel,(-sin(a),cos(a)))/max(_hp_norm(rel),1e-12)
    end
    function spread(e)
        ri,i,_,hi=e;prev,nxt=gaps(ri,i,hi)
        _hp_port_spread(prev,nxt)
    end
    for (a,b) in pairs
        dth=clamp(.3*(pull(a)-pull(b))+.05*(spread(a)-spread(b)),-.02,.02)
        oldg=[gaps(r,i,h) for (r,i,_,h) in ends]
        ra,ia,_,ha=a;rb,ib,_,hb=b
        angles[(ra,ia)]+=dth;angles[(rb,ib)]-=dth
        newg=[gaps(r,i,h) for (r,i,_,h) in ends]
        if any(minimum(newg[k])<deg2rad(18.0) && minimum(newg[k])<minimum(oldg[k])-1e-12 for k in eachindex(ends))
            angles[(ra,ia)]-=dth;angles[(rb,ib)]+=dth;continue
        end
        c,rad=p.holes[ha].second;t=angles[(ra,ia)]
        target=_hp_add(c,(rad*cos(t),rad*sin(t)))
        motions[ra][ia]=_hp_sub(target,p.runs[ra].points[ia])
        motions[rb][ib]=_hp_sub((target[1],p.height-target[2]),p.runs[rb].points[ib])
    end
    motions
end

function _hp_adaptive(p,prox;fixed=Set{String}())
    runs=HCPlotRun[]
    for (ri,r) in enumerate(p.runs)
        if r.name in fixed;push!(runs,r);continue;end
        P=r.points;n=length(P);spacing=fill(40.0,n)
        for i in 1:n
            if r.closed || !(i in (1,n))
                a=_hp_sub(P[i],P[mod1(i-1,n)]);b=_hp_sub(P[mod1(i+1,n)],P[i])
                angle=abs(atan(_hp_cross(a,b),_hp_dot(a,b)));len=.5*(_hp_norm(a)+_hp_norm(b))
                spacing[i]=clamp(.35/max(angle/max(len,1e-9),1e-9),8,40)
            end
            spacing[i]=min(spacing[i],clamp(.8prox[ri][i],8,40))
        end
        for _ in 1:3
            old=copy(spacing)
            for i in (r.closed ? (1:n) : (2:n-1));spacing[i]=min(old[i],.5*(old[mod1(i-1,n)]+old[mod1(i+1,n)]));end
        end
        push!(runs,_hp_run(r,_hp_sample(P,spacing,r.closed)))
    end
    _hp_rebuild(p,runs)
end

"""Diagnostic smoother: every proposal is audited; never call this the production certificate."""
# cf: LIB-REQUIREMENTS.md §3.2, §7 — retain the hole constraint while localizing a failed move
function _hp_hold_obstacles!(move,ports,p,audit,pairs)
    indices=Dict(r.name=>i for (i,r) in enumerate(p.runs));held=Set{Tuple{Int,Int}}()
    for (name,_,i) in audit.hole_entries
        ri=indices[name];push!(held,(ri,i),(ri,mod1(i+1,length(p.runs[ri].points))))
    end
    for (name,i) in audit.frame_entries;push!(held,(indices[name],i));end
    for pair in pairs
        if any((e[1],e[2]) in held for e in pair)
            for e in pair;push!(held,(e[1],e[2]));end
        end
    end
    changed=false
    for (ri,i) in held
        changed |= _hp_norm(move[ri][i])+_hp_norm(ports[ri][i])>1e-12
        move[ri][i]=(0.0,0.0);ports[ri][i]=(0.0,0.0)
    end
    changed
end
function _hp_displacement(before,after)
    maximum((_hp_dist(a,b) for (r,s) in zip(before.runs,after.runs)
        for (a,b) in zip(r.points,s.points));init=0.0)
end

function _hp_gate_hold!(move,velocity,ports,p,held,pairs)
    frozen=[falses(length(r.points)) for r in p.runs]
    for (ri,i) in held
        1<=i<=length(frozen[ri]) && (frozen[ri][i]=true)
    end
    for (a,b) in pairs
        if frozen[a[1]][a[2]] || frozen[b[1]][b[2]]
            frozen[a[1]][a[2]]=true;frozen[b[1]][b[2]]=true
        end
    end
    count=0
    for ri in eachindex(frozen),i in eachindex(frozen[ri])
        if frozen[ri][i]
            move[ri][i]=(0.,0.);velocity[ri][i]=(0.,0.);ports[ri][i]=(0.,0.);count+=1
        end
    end
    count
end

# cf: relax23.py local_steps/remesh; owner turn 4 — no intermediate 30-degree gate
function _hp_relax(placement;rounds=24,local_steps=20,DEBUG_MODE=false,log=devnull,
                   fixed=Set{String}(),port_check=_hp_port_audit,check_final_angle=true,
                   placement_check=nothing)
    started=time_ns();force_ns=UInt64(0);gate_ns=UInt64(0);remesh_ns=UInt64(0);debug_ns=UInt64(0)
    local_checks=0;global_checks=0;placement_checks=0;accepted=0;rejected=0;remeshed=0;frozen=0;history=NamedTuple[]
    p=_hp_rebuild(placement,[r.name in fixed ? r : _hp_run(r,_hp_sample(r.points,fill(10.,length(r.points)),r.closed)) for r in placement.runs])
    gate_workspace=HCPlotGateWorkspace(p);force_workspace=HCPlotForceWorkspace(p)
    t=time_ns();geometry=_hp_gate_geometry(p,gate_workspace);gate_ns+=time_ns()-t;local_checks+=1
    geometry.ok || throw(ArgumentError("fine sampling failed the local hole/frame/intersection check"))
    if !isnothing(placement_check)
        placement_checks+=1
        placement_check(p) || throw(ArgumentError("fine sampling changed an empty-curve placement"))
    end
    state=HCPlotGateState(p.runs,geometry.crossings)
    n=length(p.runs);partners=falses(n,n)
    for i in 1:n,j in 1:n
        partners[i,j]=i!=j && !isempty(intersect(abs.(p.runs[i].word),abs.(p.runs[j].word)))
    end
    # Prototype: initial adaptive pass, then twenty local steps and one remesh per round.
    for round in 0:rounds
        if round>0
            pairs=_hp_ports(p);velocity=[fill((0.,0.),length(r.points)) for r in p.runs]
            for step in 1:local_steps
                t=time_ns();force,cap,_=_hp_force_state(p,state.crossings,partners,force_workspace);force_ns+=time_ns()-t
                move=[[begin
                    v=_hp_add(_hp_mul(velocity[ri][i],.8),_hp_mul(force[ri][i],.5));len=_hp_fastnorm(v)
                    _hp_mul(v,min(1.,cap[ri][i]/max(len,1e-12)))
                end for i in eachindex(r.points)] for (ri,r) in enumerate(p.runs)]
                for (ri,r) in enumerate(p.runs)
                    r.name in fixed && fill!(move[ri],(0.,0.))
                end
                maximum((_hp_fastnorm(v) for row in move for v in row);init=0.)<1e-3 && break
                dp=_hp_port_motion(p,pairs);success=false
                for (ri,r) in enumerate(p.runs)
                    r.name in fixed && fill!(dp[ri],(0.,0.))
                end
                for attempt in 0:12
                    scale=2.0^-attempt
                    # LIB-REQUIREMENTS.md §7: retain port motion on retries.
                    trial=[_hp_run(r,[_hp_add(r.points[i],_hp_mul(_hp_add(move[ri][i],dp[ri][i]),scale))
                        for i in eachindex(r.points)]) for (ri,r) in enumerate(p.runs)]
                    for (a,b) in pairs
                        ri,i,_,hi=a;c,rad=p.holes[hi].second;v=_hp_sub(trial[ri].points[i],c)
                        p.runs[ri].name in fixed && continue
                        target=_hp_add(c,_hp_mul(v,rad/_hp_fastnorm(v)))
                        trial[ri].points[i]=target
                        rj,j,_,_=b;trial[rj].points[j]=(target[1],p.height-target[2])
                    end
                    candidate=_hp_rebuild(p,trial)
                    _hp_displacement(p,candidate)<1e-12 && break
                    if !isnothing(placement_check)
                        placement_checks+=1
                        placement_check(candidate) || continue
                    end
                    t=time_ns();geometry=_hp_gate_geometry(candidate,gate_workspace)
                    check=_hp_gate_trial(candidate,geometry,state);gate_ns+=time_ns()-t;local_checks+=1
                    if check.ok && !port_check(p,candidate;pairs).ok
                        for row in dp;fill!(row,(0.,0.));end
                        continue
                    end
                    if check.ok
                        if DEBUG_MODE
                            t=time_ns();a=_hp_audit(candidate);debug_ns+=time_ns()-t;global_checks+=1
                            a.ok || error("local plot gate disagrees with independent DEBUG_MODE audit")
                        end
                        p=candidate;state=check.state;velocity=[[_hp_mul(v,scale) for v in row] for row in move]
                        accepted+=1;success=true;break
                    end
                    if !isempty(check.held)
                        frozen+=_hp_gate_hold!(move,velocity,dp,p,check.held,pairs)
                    end
                end
                if !success
                    rejected+=1
                    for row in velocity;fill!(row,(0.,0.));end
                end
                row=(round=round,step=step,accepted=success,points=sum((length(r.points) for r in p.runs);init=0))
                push!(history,row)
                if DEBUG_MODE;println(log,row);flush(log);end
            end
        end
        t=time_ns();_,_,prox=_hp_force_state(p,state.crossings,partners,force_workspace);force_ns+=time_ns()-t
        t=time_ns();candidate=_hp_adaptive(p,prox;fixed);remesh_ns+=time_ns()-t
        t=time_ns();geometry=_hp_gate_geometry(candidate,gate_workspace)
        check=_hp_gate_trial(candidate,geometry,state);gate_ns+=time_ns()-t;local_checks+=1
        place_ok=true
        if !isnothing(placement_check)
            placement_checks+=1;place_ok=placement_check(candidate)
        end
        if check.ok && place_ok
            if DEBUG_MODE
                t=time_ns();a=_hp_audit(candidate);debug_ns+=time_ns()-t;global_checks+=1
                a.ok || error("local remesh gate disagrees with independent DEBUG_MODE audit")
            end
            p=candidate;state=check.state;remeshed+=1
        else
            # All incident hole segments are checked, including port collars
            # which the prototype exempts. Preserve only the runs whose
            # resampling is blocked instead of retaining every fine run.
            println(log,(round=round,batch_remesh=false,structure=check.structure,holes=check.holes,
                overlaps=check.overlaps,held_runs=unique(first.(check.held))))
            for ri in 1:n
                runs=copy(p.runs);runs[ri]=candidate.runs[ri];localp=_hp_rebuild(p,runs)
                t=time_ns();lg=_hp_gate_geometry(localp,gate_workspace)
                lc=_hp_gate_trial(localp,lg,state);gate_ns+=time_ns()-t;local_checks+=1
                local_place_ok=true
                if !isnothing(placement_check)
                    placement_checks+=1;local_place_ok=placement_check(localp)
                end
                if lc.ok && local_place_ok
                    if DEBUG_MODE
                        t=time_ns();a=_hp_audit(localp);debug_ns+=time_ns()-t;global_checks+=1
                        a.ok || error("local per-run remesh disagrees with independent DEBUG_MODE audit")
                    end
                    p=localp;state=lc.state;remeshed+=1
                end
            end
        end
        println(log,(round=round,remesh=check.ok && place_ok,accepted=accepted,rejected=rejected,points=sum((length(r.points) for r in p.runs);init=0)))
    end
    t=time_ns();audit=_hp_audit(p);global_checks+=1;ports=port_check(placement,p)
    !isnothing(placement_check) && (placement_checks+=1)
    place_ok=isnothing(placement_check) || placement_check(p)
    angle_ok=!check_final_angle || isnothing(audit.minimum_stencil_angle) || audit.minimum_stencil_angle >= _HP_MIN_CROSSING_ANGLE
    audit=merge(audit,(ports=ports,angle_ok=angle_ok,empty_places=place_ok,ok=audit.ok && ports.ok && angle_ok && place_ok))
    final_ns=time_ns()-t
    (placement=p,audit=audit,steps=accepted,rejected=rejected,remeshed=remeshed,frozen=frozen,
     seconds=(time_ns()-started)/1e9,history=history,mode=DEBUG_MODE ? :debug : :production,production_ready=true,
     checks=(local_checks=local_checks,global_gauss=global_checks,placement_checks=placement_checks),
     timings=(force=force_ns/1e9,gate=gate_ns/1e9,remesh=remesh_ns/1e9,debug=debug_ns/1e9,final=final_ns/1e9))
end
