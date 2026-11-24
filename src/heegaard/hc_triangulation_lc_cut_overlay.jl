# cf: .research-scratch/PJ260922/astra-lanes/snappy-bridge/answer-6.md — oriented standard octagon and ordered crossing events.
# Adjustments: preserve actual blown-up endpoints, reject shared support, and bound genus-two search; docs/specs/heegaard-snappy.md.
# Geometry uses hc_triangulation_lc_cut.jl::_hc_cut_intersection and its exact determinant helpers.
function _hc_standard_polygon_overlay(polygon_word::AbstractVector{<:Integer}, arcs::AbstractVector;
                                      first_label::Integer=0, endpoint_point=nothing, arc_points=nothing)
    # arc = (id=..., from=(side=..., rank=...), to=(side=..., rank=...));
    # rank is 1:N in the boundary orientation of that signed side occurrence.
    W=Int.(polygon_word); length(W)%4==0 || throw(ArgumentError("polygon_word length must be 4g"))
    g=length(W)÷4; g in (1,2) || throw(ArgumentError("this overlay implements only g=1,2"))
    A=Int[]; B=Int[]
    for i in 1:g
        a,b,c,d=W[4i-3:4i]; c==-a && d==-b || throw(ArgumentError("block $i is not a b -a -b"))
        push!(A,a); push!(B,b)
    end
    length(unique(abs.(vcat(A,B))))==2g || throw(ArgumentError("edge pairs are not distinct"))
    side_index=Dict(s=>i for (i,s) in enumerate(W)); length(side_index)==4g || throw(ArgumentError("signed sides repeat"))
    ids=[a.id for a in arcs]; length(unique(ids))==length(ids) || throw(ArgumentError("arc ids repeat"))
    Rat=Rational{BigInt}; Z=BigInt(0)//BigInt(1); O=BigInt(1)//BigInt(1)
    Q(n::Integer)=BigInt(n)//BigInt(1)
    byside=Dict{Int,Vector{Tuple{Int,Symbol,Int}}}()
    for (j,a) in enumerate(arcs), w in (:from,:to)
        e=getproperty(a,w); s=Int(e.side); r=Int(e.rank)
        haskey(side_index,s) || throw(ArgumentError("arc $(a.id): side $s is absent")); r>0 || throw(ArgumentError("rank must be positive"))
        push!(get!(byside,s,Tuple{Int,Symbol,Int}[]),(j,w,r))
    end
    endq=Dict{Tuple{Int,Symbol},Rat}()
    for (s,v) in byside
        sort([x[3] for x in v])==collect(1:length(v)) || throw(ArgumentError("side $s ranks must be exactly 1:N"))
        for (j,w,r) in v; endq[(j,w)]=BigInt(r)//BigInt(length(v)+1); end
    end
    xof(s,u)=Q(side_index[s]-1)+u
    xy(p)=(x=p[1],y=p[2])
    point(s,u)=endpoint_point===nothing ? (x=xof(s,u),y=xof(s,u)^2) : xy(endpoint_point(s,u))
    sub(p,q)=_hc_cut_sub((p.x,p.y),(q.x,q.y)); det=_hc_cut_det
    crossdata(S,T)=_hc_cut_intersection((S.p0.x,S.p0.y),(S.p1.x,S.p1.y),(T.p0.x,T.p0.y),(T.p1.x,T.p1.y))
    isegs=Any[]
    for (j,a) in enumerate(arcs)
        f=(side=Int(a.from.side),u=endq[(j,:from)]); t=(side=Int(a.to.side),u=endq[(j,:to)])
        f==t && throw(ArgumentError("arc $(a.id) has equal endpoints"))
        push!(isegs,(id=a.id,index=j,from=f,to=t,p0=arc_points===nothing ? point(f.side,f.u) : xy(arc_points[a.id][1]),p1=arc_points===nothing ? point(t.side,t.u) : xy(arc_points[a.id][2])))
    end
    for i in eachindex(isegs), j in i+1:length(isegs)
        x=crossdata(isegs[i],isegs[j]); x===nothing && continue
        0<x.t<1 && 0<x.s<1 && throw(ArgumentError("input arcs $(isegs[i].id), $(isegs[j].id) alternate on the boundary"))
    end
    kmax=3g-2; first_label==0 && (first_label=kmax+1); first_label>kmax || throw(ArgumentError("first_label must exceed $kmax"))
    order=g==1 ? ["l1","m1"] : ["l1","l2","m1","m2","n1"]
    constants=g==1 ? Dict("l1"=>[-1],"m1"=>[1]) :
        Dict("l1"=>[-1,-3],"l2"=>[-2,4],"m1"=>[1],"m2"=>[2],"n1"=>[3,-4])
    key(a,b)=isless(a,b) ? (a,b) : (b,a)
    expected=Dict{Tuple{String,String},Any}()
    for i in 1:g
        expected[key("l$i","m$i")]=(label=i,tokens=Dict("l$i"=>-i,"m$i"=>i))
    end
    g==2 && (expected[key("l1","n1")]=(label=3,tokens=Dict("l1"=>-3,"n1"=>3));
             expected[key("l2","n1")]=(label=4,tokens=Dict("l2"=>4,"n1"=>-4)))
    EP(s,u)=(side=Int(s),u=u,point=point(Int(s),u))
    function makepieces(P)
        out=Any[]
        for i in 1:g
            s=A[i]; u=P.ul[i]; f=isodd(i) ? EP(s,u) : EP(-s,O-u); t=isodd(i) ? EP(-s,O-u) : EP(s,u)
            push!(out,(curve="l$i",piece=1,from=f,to=t,p0=f.point,p1=t.point))
        end
        for i in 1:g
            s=B[i]; u=P.um[i]; f=isodd(i) ? EP(s,u) : EP(-s,O-u); t=isodd(i) ? EP(-s,O-u) : EP(s,u)
            push!(out,(curve="m$i",piece=1,from=f,to=t,p0=f.point,p1=t.point))
        end
        if g==2
            # Two nested arcs; the even handle was reversed above to obtain n1 signs +3,-4.
            f1=EP(B[1],P.un[1]); t1=EP(-B[2],O-P.un[2]); f2=EP(B[2],P.un[2]); t2=EP(-B[1],O-P.un[1])
            push!(out,(curve="n1",piece=1,from=f1,to=t1,p0=f1.point,p1=t1.point))
            push!(out,(curve="n1",piece=2,from=f2,to=t2,p0=f2.point,p1=t2.point))
        end
        out
    end
    inputloc=Set((p.x,p.y) for S in isegs for p in (S.p0,S.p1))
    function admissible(pieces)
        seen=Set{Tuple{Int,Rat}}()
        for p in pieces, e in (p.from,p.to)
            (e.point.x,e.point.y) in inputloc && return nothing
            (e.side,e.u) in seen && return nothing
            push!(seen,(e.side,e.u))
        end
        kk=Any[]; count=Dict(k=>0 for k in keys(expected))
        for i in eachindex(pieces), j in i+1:length(pieces)
            x=crossdata(pieces[i],pieces[j]); if x===nothing
                det(sub(pieces[j].p1,pieces[j].p0),sub(pieces[i].p0,pieces[j].p0))==0 && return nothing; continue
            end
            if 0<x.t<1 && 0<x.s<1
                pieces[i].curve==pieces[j].curve && return nothing
                k=key(pieces[i].curve,pieces[j].curve); haskey(expected,k) || return nothing
                E=expected[k]; ti=-x.sign*E.label; tj=x.sign*E.label
                ti==E.tokens[pieces[i].curve] && tj==E.tokens[pieces[j].curve] || return nothing
                count[k]+=1; count[k]==1 || return nothing
                push!(kk,(kpiece=i,t=x.t,curve=pieces[i].curve,token=ti,label=E.label,other=pieces[j].curve))
                push!(kk,(kpiece=j,t=x.s,curve=pieces[j].curve,token=tj,label=E.label,other=pieces[i].curve))
            end
        end
        all(v==1 for v in values(count)) || return nothing
        hits=Any[]
        for (j,S) in enumerate(isegs), (k,T) in enumerate(pieces)
            x=crossdata(S,T); if x===nothing
                (g==1 || det(sub(T.p1,T.p0),sub(S.p0,T.p0))==0) && return nothing; continue
            end
            if ((Z<=x.t<=O) && (x.s==Z || x.s==O)) || ((Z<=x.s<=O) && (x.t==Z || x.t==O)); return nothing; end
            0<x.t<1 && 0<x.s<1 && push!(hits,(arc=j,kpiece=k,old_t=x.t,k_t=x.s,old_sign=-x.sign,k_sign=x.sign))
        end
        for k in eachindex(pieces)
            ts=vcat([e.t for e in kk if e.kpiece==k],[h.k_t for h in hits if h.kpiece==k])
            length(ts)==length(unique(ts)) || return nothing
        end
        for j in eachindex(isegs)
            ts=[h.old_t for h in hits if h.arc==j]; length(ts)==length(unique(ts)) || return nothing
        end
        (pieces=pieces,kk=kk,hits=hits)
    end
    n=length(arcs); bound=8(n+2)^2+17
    chosen=nothing; params=nothing; trials=0
    if g==1
        for total in 2:2bound
            chosen===nothing || break
            for il in max(1,total-bound):min(bound,total-1)
                im=total-il; trials+=1
                P=(ul=Rat[BigInt(il)//BigInt(bound+1)],um=Rat[BigInt(im)//BigInt(bound+1)],un=Rat[])
                C=admissible(makepieces(P))
                if C!==nothing; chosen=C; params=P; break; end
            end
        end
    else
        r=1
        while chosen===nothing && r<=bound
            trials+=1; x=BigInt(1)//BigInt(r+1)
            P=(ul=Rat[(8+3x)/20,(9+4x)/20],um=Rat[(14+5x)/20,(15+3x)/20],un=Rat[(2+x)/20,(3+2x)/20])
            C=admissible(makepieces(P))
            if C!==nothing; chosen=C; params=P; end
            r+=1
        end
    end
    chosen===nothing && error("fixed rational grid exhausted (genus=$g, bound=$bound)")
    pieces=chosen.pieces; kevents=Dict(c=>Any[] for c in order)
    for e in chosen.kk
        p=pieces[e.kpiece]; push!(kevents[e.curve],(t=Q(p.piece-1)+e.t,token=e.token,passage_sign=sign(e.token),kind=:K_K,label=e.label,other=e.other))
    end
    edge_events=Dict(e=>Any[] for e in unique(abs.(W))); arc_events=Dict(a.id=>Any[] for a in arcs)
    endpoint_pairing=Any[]; next=Int(first_label)
    for c in order
        inds=sort([i for i in eachindex(pieces) if pieces[i].curve==c];by=i->pieces[i].piece)
        for (z,k) in enumerate(inds)
            p=pieces[k]; prev=pieces[inds[mod1(z-1,length(inds))]].to
            prev.side==-p.from.side && prev.u==O-p.from.u || error("internal K endpoint pairing failure")
            label=next; next+=1; sg=sign(p.from.side); ktok=sg*label; etok=-sg*label
            phase=Q(p.piece-1); et=p.from.side>0 ? p.from.u : O-p.from.u
            push!(kevents[c],(t=phase,token=ktok,passage_sign=sign(ktok),kind=:side,label=label,edge=abs(p.from.side)))
            push!(edge_events[abs(p.from.side)],(t=et,token=etok,curve=c,label=label,passage_sign=sign(etok)))
            push!(endpoint_pairing,(curve=c,piece=p.piece,outgoing=(side=p.from.side,u=p.from.u),incoming=(side=prev.side,u=prev.u),label=label))
            for j in eachindex(arcs)
                h=findfirst(x->x.arc==j && x.kpiece==k,chosen.hits); h===nothing && continue
                H=chosen.hits[h]; label=next; next+=1
                push!(kevents[c],(t=phase+H.k_t,token=H.k_sign*label,passage_sign=H.k_sign,kind=:input,label=label,arc=arcs[j].id))
                push!(arc_events[arcs[j].id],(t=H.old_t,token=H.old_sign*label,curve=c,label=label,passage_sign=sign(H.old_sign)))
            end
        end
    end
    words=Dict{String,Vector{Int}}(); rotations=Dict{String,Int}(); basepoints=Dict{String,Rat}(); ordered_kevents=Dict{String,Any}()
    for c in order
        raw=sort(kevents[c];by=e->e.t); shift=nothing
        for s in 0:length(raw)-1
            rr=vcat(raw[s+1:end],raw[1:s])
            if [e.token for e in rr if abs(e.token)<=kmax]==constants[c]; shift=s; break; end
        end
        shift===nothing && error("K constant word mismatch on $c")
        rr=vcat(raw[shift+1:end],raw[1:shift]); t0=rr[1].t
        period=Q(maximum(p.piece for p in pieces if p.curve==c))
        rr=[merge(e,(t=(e.t>=t0 ? e.t-t0 : e.t-t0+period),)) for e in rr]
        ordered_kevents[c]=rr; words[c]=[e.token for e in rr]; rotations[c]=shift; basepoints[c]=t0
    end
    foreach(v->sort!(v;by=e->e.t),values(edge_events)); foreach(v->sort!(v;by=e->e.t),values(arc_events))
    reference_arcs=Dict(c=>[(piece=pieces[i].piece,from=(side=pieces[i].from.side,u=pieces[i].from.u),to=(side=pieces[i].to.side,u=pieces[i].to.u))
                             for i in eachindex(pieces) if pieces[i].curve==c] for c in order)
    input_realization=Dict(S.id=>(from=S.from,to=S.to,from_point=S.p0,to_point=S.p1) for S in isegs)
    (;genus=g,input_realization,reference_arcs,endpoint_pairing,edge_events,arc_events,
      arc_words=Dict(id=>[e.token for e in arc_events[id]] for id in ids),kevents=ordered_kevents,words,next_label=next,
      correspondence=(input_arcs=Dict(id=>id for id in ids),boundary_endpoints=:fixed,regions=:identity_before_K_subdivision),
      certificate=(polygon_word=W,K_constants=constants,word_rotations=rotations,original_basepoint_t=basepoints,grid_bound=bound,grid_trials=trials,
                   parameters=params,grid_rule=(g==1 ? :legacy_diagonal : :x_equals_1_over_r_plus_1),coordinate_model="strictly convex boundary x -> (x,x^2); isotopy rel boundary"))
end

function _hc_cut_overlay_arc_specs(graph,cut,side,port,n)
    word=Int.(cut["polygon_word"]);cotree=Int.(cut["dual_tree_edges"])
    starts=Dict(s=>side[s>0 ? 2s-1 : -2s][1] for s in word)
    groups=Dict(s=>Int[] for s in word)
    for half in keys(port)
        s=argmin(s->mod(starts[s]-port[half],n),word)
        push!(groups[s],half)
    end
    which=Dict{Int,Any}()
    for s in word
        group=sort(groups[s];by=h->mod(port[h]-starts[s],n))
        for (rank,h) in enumerate(group);which[h]=(side=s,rank=rank);end
    end
    [(id=e,from=which[2*e-1],to=which[2*e]) for e in cotree]
end

function _hc_cut_overlay_layout(graph,cut,side,port,points,chords,endpoint;actual_geometry::Bool)
    g=Int(cut["genus"]);word=Int.(cut["polygon_word"]);cotree=Int.(cut["dual_tree_edges"])
    arcs=_hc_cut_overlay_arc_specs(graph,cut,side,port,length(points))
    result=try
        actual_geometry ? _hc_standard_polygon_overlay(word,arcs;first_label=length(graph.vertices)+3g-1,
        endpoint_point=endpoint,arc_points=chords) : _hc_standard_polygon_overlay(word,arcs;first_label=length(graph.vertices)+3g-1)
    catch err
        if err isa ErrorException && startswith(err.msg,"fixed rational grid exhausted")
            throw(HCTriangulationError(:overlay_grid_exhausted,err.msg))
        end
        rethrow()
    end
    events=Dict(e.id=>Any[] for e in graph.edges)
    for (e,v) in result.edge_events;append!(events[e],[(t=x.t,token=x.token) for x in v]);end
    for (e,v) in result.arc_events;append!(events[e],[(t=x.t,token=x.token) for x in v]);end
    labels=Dict(v=>i+3g-2 for (i,v) in enumerate(graph.vertices))
    (;labels,events,words=result.words,chosen=(us=vcat(result.certificate.parameters.ul,result.certificate.parameters.um),),
        trials=result.certificate.grid_trials,overlay=result)
end
