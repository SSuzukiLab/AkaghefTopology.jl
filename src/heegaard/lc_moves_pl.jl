# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-3-full.md §2 — certified normal boundaries and periodic PL extraction; reconstructed.
# Boundary events follow increasing s on BOTH circles; the annulus is left of t=0 and right of t=1.

_move_plpoint(p)=(_move_rational(p[1]),_move_rational(p[2]))
_move_portkey(p)=(p.component,p.edge,p.u)

"""All periodic segment events, including returning, winding and closed runs.
The first tangent receives -sign(det); deck ranges are exact finite bounds.
"""
function _move_pl_events(runs)
    events=[Any[] for _ in runs];hits=Any[]
    segs=[(run=i,seg=j,p0=r.points[j],p1=r.points[j+1]) for(i,r) in enumerate(runs) for j in 1:length(r.points)-1]
    function joins(x,y,deck)
        x.run==y.run || return false
        shift(p)=(p[1]+deck,p[2])
        x.seg==y.seg+1 && x.p0==shift(y.p1) && return true
        y.seg==x.seg+1 && x.p1==shift(y.p0) && return true
        if runs[x.run].closed
            lastseg=length(runs[x.run].points)-1
            x.seg==1 && y.seg==lastseg && x.p0==shift(y.p1) && return true
            y.seg==1 && x.seg==lastseg && x.p1==shift(y.p0) && return true
        end
        false
    end
    for a in eachindex(segs),b in 1:a
        x,y=segs[a],segs[b];u=_move_sub2(x.p1,x.p0);v=_move_sub2(y.p1,y.p0)
        low=ceil(Int,min(x.p0[1],x.p1[1])-max(y.p0[1],y.p1[1]))
        high=floor(Int,max(x.p0[1],x.p1[1])-min(y.p0[1],y.p1[1]))
        for deck in low:high
            a==b && deck<=0 && continue
            offset=(y.p0[1]+deck-x.p0[1],y.p0[2]-x.p0[2]);det=_move_cross2(u,v)
            if det==0
                _move_cross2(offset,u)!=0 && continue
                axis=u[1]!=0 ? 1 : 2
                lo=max(min(x.p0[axis],x.p1[axis]),min(y.p0[axis],y.p1[axis])+(axis==1 ? deck : 0))
                hi=min(max(x.p0[axis],x.p1[axis]),max(y.p0[axis],y.p1[axis])+(axis==1 ? deck : 0))
                lo<hi && throw(LCMoveError(:affine_overlay_degenerate,"overlapping PL segments",Dict()))
                lo==hi && !joins(x,y,deck) && throw(LCMoveError(:affine_overlay_degenerate,"tangential PL contact",Dict()))
                continue
            end
            t=_move_cross2(offset,v)/det;z=_move_cross2(offset,u)/det
            0<=t<=1 && 0<=z<=1 || continue
            if t in (0,1) || z in (0,1)
                joins(x,y,deck) && t in (0,1) && z in (0,1) && continue
                throw(LCMoveError(:affine_overlay_degenerate,"PL vertex or boundary-port collision",Dict()))
            end
            key=(x.run,x.seg,y.run,y.seg,deck);sgn=-Int(sign(det));h=length(hits)+1
            push!(hits,(key=key,a=x.run,b=y.run,sign=sgn))
            push!(events[x.run],(x.seg-1+t,sgn,h));push!(events[y.run],(y.seg-1+z,-sgn,h))
        end
    end
    for ev in events
        sort!(ev;by=first)
        length(unique(first.(ev)))==length(ev) || throw(LCMoveError(:affine_overlay_degenerate,"multiple crossing",Dict()))
    end
    (;events,hits)
end

function _move_annulus_certificate(lc,raw)
    isempty(lc.places) || throw(LCMoveError(:annulus_place_embedding_required,"placed circles need a normal-cell embedding",Dict()))
    get(raw,"kind",nothing)=="embedded_annulus" || throw(LCMoveError(:annulus_embedding_required,"expected embedded_annulus",Dict()))
    faces,right=_lc_faces(_lc_nonempty_words(lc.words));occ=_lc_rotation_occ(lc.words)
    rank(d,u)=findfirst(==(d),only(f for f in faces if f.id==right[d]).corners)-1+(_move_dir(d)==1 ? u : 1-u)
    points=Any[];cycles=Vector{Int}[];arcs=Dict{String,Vector{Any}}();ids=Set{String}()
    for level in 0:1
        key="boundary$level";haskey(raw,key) || throw(LCMoveError(:annulus_embedding_required,"missing boundary",Dict("boundary"=>key)))
        cycle=Int[]
        for event in raw[key]["events"]
            id=String(event["id"]);id in ids && throw(LCMoveError(:annulus_port_collision,"duplicate port id",Dict()));push!(ids,id)
            p=_move_route_point(lc,occ,right,event);s=_move_rational(event["s"])
            0<=s<1 || throw(LCMoveError(:annulus_embedding_required,"boundary s must be in [0,1)",Dict()))
            push!(points,merge(p,(id=id,s=s,level=level)));push!(cycle,length(points))
        end
        isempty(cycle) && throw(LCMoveError(:annulus_boundary_embedding_required,"a boundary inside one face needs an explicit cell patch",Dict()))
        ss=[points[i].s for i in cycle]
        length(unique(ss))==length(ss) && (length(ss)==1 || count(i->ss[mod1(i+1,length(ss))]<ss[i],eachindex(ss))==1) ||
            throw(LCMoveError(:annulus_boundary_not_embedded,"boundary parameter order disagrees",Dict()))
        for i in eachindex(cycle)
            p=points[cycle[i]];q=points[cycle[mod1(i+1,length(cycle))]];d=_lc_tau(lc.words,occ,p.dart)
            right[d]==q.face || throw(LCMoveError(:annulus_boundary_not_embedded,"boundary leaves the current face",Dict()))
            push!(get!(arcs,q.face,Any[]),(rank(d,p.u),rank(q.dart,q.u)))
        end
        push!(cycles,cycle)
    end
    length(unique(_move_portkey.(points)))==length(points) || throw(LCMoveError(:annulus_port_collision,"coincident edge ports",Dict()))
    for (f,cs) in arcs,i in eachindex(cs),j in 1:i-1
        a,b=sort(collect(cs[i]));c,d=cs[j]
        xor(a<c<b,a<d<b) && throw(LCMoveError(:annulus_boundary_not_embedded,"boundary chords alternate",Dict("face"=>f)))
    end
    names=_component_names(lc);bn=["__move_boundary0","__move_boundary1"]
    any(haskey(lc.words,n) for n in bn) && throw(LCMoveError(:component_name,"reserved support component is occupied",Dict()))
    append!(names,bn);index=Dict(n=>i for(i,n) in enumerate(names));fresh=maximum(vcat([0],abs.(reduce(vcat,values(lc.words);init=Int[]))))
    portlabels=Dict(fresh+i=>p.id for(i,p) in enumerate(points))
    before=[Int[],Int[]];after=[[ _move_dir(points[i].dart)*(fresh+i) for i in cycle] for cycle in cycles]
    locations=[(index[n],1) for n in bn];ops=sort([(p.component,p.edge,p.u,-_move_dir(p.dart)*(fresh+i)) for(i,p) in enumerate(points)];by=x->x[1:3])
    _emit_gap_inserts!(locations,before,after,index,[(n,e+1,-x) for(n,e,u,x) in ops])
    labels=Set(abs.(reduce(vcat,values(lc.words);init=Int[])))
    result=_kernel_apply(names,lc.words,before,after,locations,Int[],collect(1:length(before));label_map=Dict(k=>k for k in labels))[1]
    words=Dict(n=>result[i] for(i,n) in enumerate(names));fs,rf=_lc_faces(words);oo=_lc_rotation_occ(words)
    length(fs)-length(oo)==2-2lc.genus || throw(LCMoveError(:annulus_boundary_not_embedded,"normal overlay changes the surface",Dict()))
    tau(d)=_lc_tau(words,oo,d)
    firsttoken=first(words[bn[1]]);d=(abs(firsttoken),firsttoken>0 ? "+" : "-","out")
    selected=Set([rf[tau(d)]]);queue=collect(selected)
    for f in queue
        for c in only(x for x in fs if x.id==f).corners
            oo[c[1]][c[2]][1] in bn && continue
            nxt=rf[tau(c)]
            nxt in selected || (push!(selected,nxt);push!(queue,nxt))
        end
    end
    for (level,n) in enumerate(bn),x in words[n]
        c=(abs(x),x>0 ? "+" : "-","out")
        inside=level==1 ? rf[tau(c)] : rf[c];outside=level==1 ? rf[c] : rf[tau(c)]
        inside in selected && !(outside in selected) || throw(LCMoveError(:annulus_not_annular,"selected side does not have the two oriented boundaries",Dict()))
    end
    corners=[c for f in fs if f.id in selected for c in f.corners]
    vertices=Set(c[1] for c in corners);edges=Set(min(c,tau(c)) for c in corners)
    length(vertices)-length(edges)+length(selected)==0 || throw(LCMoveError(:annulus_not_annular,"selected surface has nonzero Euler characteristic",Dict()))
    intervals=Dict{Tuple{String,String,String},Any}();streams=Dict{String,Any}()
    for n in _component_names(lc)
        w=words[n];at=findall(x->haskey(portlabels,abs(x)),w)
        if isempty(at)
            isempty(w) && continue
            x=first(w);inside=rf[(abs(x),x>0 ? "+" : "-","out")] in selected
            if inside;intervals[(n,"","")]=(word=copy(w),closed=true);end
            streams[n]=(word=w,ports=at)
            continue
        end
        for (j,pos) in enumerate(at)
            stop=at[mod1(j+1,length(at))];tokens=Int[];k=mod1(pos+1,length(w))
            while k!=stop;push!(tokens,w[k]);k=mod1(k+1,length(w));end
            x=w[pos];inside=rf[(abs(x),x>0 ? "+" : "-","out")] in selected
            inside && (intervals[(n,portlabels[abs(x)],portlabels[abs(w[stop])])]=(word=tokens,closed=false))
        end
        streams[n]=(word=w,ports=at)
    end
    lookup=Dict(p.id=>p for p in points);runs=Any[];keysseen=Set{Tuple{String,String,String}}()
    for r in raw["runs"]
        get(r,"direction",1)==1 || throw(LCMoveError(:annulus_source_mismatch,"runs must follow the component orientation",Dict()))
        n=String(r["component"]);closed=Bool(get(r,"closed",false));a=closed ? "" : String(r["start_port"]);b=closed ? "" : String(r["finish_port"])
        key=(n,a,b);haskey(intervals,key) && !(key in keysseen) || throw(LCMoveError(:annulus_source_mismatch,"run is not a unique interior source interval",Dict("run"=>key)))
        push!(keysseen,key);ps=_move_plpoint.(r["points"])
        length(ps)>=2 && all(ps[i]!=ps[i+1] for i in 1:length(ps)-1) || throw(LCMoveError(:affine_overlay_degenerate,"empty PL segment",Dict()))
        all(p->0<p[2]<1,ps[2:end-1]) || throw(LCMoveError(:annulus_source_mismatch,"interior run vertex is on/outside boundary",Dict()))
        if closed
            first(ps)[2]==last(ps)[2] && 0<first(ps)[2]<1 && denominator(last(ps)[1]-first(ps)[1])==1 || throw(LCMoveError(:annulus_source_mismatch,"closed lift does not close",Dict()))
        else
            for (p,id) in ((first(ps),a),(last(ps),b))
                ref=lookup[id];p[2]==ref.level && mod(p[1],1)==ref.s || throw(LCMoveError(:annulus_source_mismatch,"PL endpoint differs from boundary port",Dict("port"=>id)))
            end
        end
        push!(runs,(component=n,start=a,finish=b,closed=closed,points=ps,source_word=intervals[key].word))
    end
    keysseen==Set(keys(intervals)) || throw(LCMoveError(:annulus_source_mismatch,"some interior source intervals are omitted",Dict()))
    ev=_move_pl_events(runs);eventlabels=Dict{Any,Int}()
    for (i,r) in enumerate(runs)
        length(ev.events[i])==length(r.source_word) || throw(LCMoveError(:annulus_source_mismatch,"PL crossing count differs",Dict("curve"=>r.component)))
        for ((_,sgn,h),x) in zip(ev.events[i],r.source_word)
            sgn==sign(x) || throw(LCMoveError(:annulus_source_mismatch,"PL passage sign differs",Dict()))
            key=ev.hits[h].key
            get(eventlabels,key,abs(x))==abs(x) || throw(LCMoveError(:annulus_source_mismatch,"paired PL crossings differ",Dict()))
            eventlabels[key]=abs(x)
        end
    end
    length(unique(values(eventlabels)))==length(eventlabels) || throw(LCMoveError(:annulus_source_mismatch,"crossing occurrence reused",Dict()))
    (;raw,runs,eventlabels,streams,portlabels,points,source_faces=Set(keys(arcs)),selected,chi=0,boundaries=2)
end

"""Compile ordered port intervals into full fragments for the existing common kernel.
Port identities survive the assembly and provide the inverse edge occurrences.
"""
function _move_compile_port_intervals(lc,A,targets)
    bykey=Dict((r.component,r.start,r.finish)=>targets[i] for(i,r) in enumerate(A.runs))
    words=Dict{String,Vector{Int}}();portanchors=Dict{String,Any}()
    for n in _component_names(lc)
        w=A.streams[n].word;at=A.streams[n].ports
        if isempty(at)
            words[n]=get(bykey,(n,"",""),copy(lc.words[n]));continue
        end
        stream=Any[]
        for (j,pos) in enumerate(at)
            stop=at[mod1(j+1,length(at))];a=A.portlabels[abs(w[pos])];b=A.portlabels[abs(w[stop])]
            push!(stream,a)
            if haskey(bykey,(n,a,b));append!(stream,bykey[(n,a,b)])
            else
                k=mod1(pos+1,length(w));while k!=stop;push!(stream,w[k]);k=mod1(k+1,length(w));end
            end
        end
        ts=findall(x->x isa Integer,stream)
        isempty(ts) && throw(LCMoveError(:missing_place_transport,"target component becomes empty inside support",Dict("component"=>n)))
        words[n]=Int[stream[k] for k in ts]
        for (i,pos) in enumerate(ts)
            stop=ts[mod1(i+1,length(ts))];ids=String[];k=mod1(pos+1,length(stream))
            while k!=stop;push!(ids,stream[k]);k=mod1(k+1,length(stream));end
            for (j,id) in enumerate(ids)
                p=only(p for p in A.points if p.id==id)
                x=_move_dir(p.dart)==1 ? stream[pos] : stream[stop]
                portanchors[id]=Dict("dart"=>[abs(x),x>0 ? "+" : "-",_move_dir(p.dart)==1 ? "out" : "in"],"u"=>string(j//(length(ids)+1)))
            end
        end
    end
    words,portanchors
end

function _move_embedded_affine(input,row,raw,acted,q)
    lc=_lc_of(input);_move_require_separated(input,acted)
    all(haskey(lc.words,n) && !_is_k(lc,n) for n in acted) || throw(LCMoveError(:fixed_reference,"acted components must be non-K input curves",Dict()))
    A=_move_annulus_certificate(lc,raw)
    mapped=[r.component in acted ? merge(r,(points=[(s+q*t,t) for(s,t) in r.points],)) : r for r in A.runs]
    target=_move_pl_events(mapped);labels=Dict{Any,Int}();fresh=maximum(vcat([0],abs.(reduce(vcat,values(lc.words);init=Int[]))))
    for h in target.hits
        a=mapped[h.a].component;b=mapped[h.b].component
        if input isa HeegaardCode && get(input.roles,a,nothing)!==nothing && get(input.roles,a,nothing)==get(input.roles,b,nothing)
            throw(LCMoveError(:same_family_crossing,"twist creates same-family intersection",Dict("curves"=>[a,b])))
        end
        old=get(A.eventlabels,h.key,nothing)
        if old!==nothing && ((a in acted)==(b in acted));labels[h.key]=old
        else;fresh+=1;labels[h.key]=fresh;end
    end
    targets=[[sgn*labels[target.hits[h].key] for(_,sgn,h) in ev] for ev in target.events]
    words,anchors=_move_compile_port_intervals(lc,A,targets)
    inverse_support=deepcopy(raw)
    for level in 0:1,event in inverse_support["boundary$level"]["events"]
        merge!(event,anchors[String(event["id"])])
    end
    for (i,r) in enumerate(inverse_support["runs"])
        r["points"]=[[string(s),string(t)] for(s,t) in mapped[i].points]
    end
    names=_component_names(lc);oldlabels=Set(abs.(reduce(vcat,values(lc.words);init=Int[])));newlabels=Set(abs.(reduce(vcat,values(words);init=Int[])))
    _witness_from(input,row,names,[(i,1) for i in eachindex(names)],Int[],[lc.words[n] for n in names],[words[n] for n in names],collect(eachindex(names));
        extra=Dict{String,Any}("label_map"=>Dict(k=>k for k in intersect(oldlabels,newlabels)),"global_inverse"=>true,
            "embedded_annulus_after"=>inverse_support,"acted"=>acted,"q"=>q,
            "support"=>Dict("kind"=>"embedded_annulus","certificate"=>raw,"chi"=>0,"boundary_components"=>2)))
end
