# cf: lib/topology.jl/scripts/heegaard/lc_placement.py place — accepted cut geometry and deterministic sampling
# cf: lib/topology.jl/docs/specs/lickorish-code.md §6.2–6.3 — face construction, gluing and incircles
# `with_contacts=true` also returns the LC label at each end of every open run, keyed by
# (run name, :start or :end). Reference placements always return it.
function _hp_place(lc::LickorishCode, roles=Dict{String,Symbol}(); names=hc_component_names(lc), reference=false, with_contacts=false,infinity_face=nothing)
    g=lc.genus;words=lc.words;knames=_lc_k_names(g)
    curves=[c for c in names if !(c in knames) && !isempty(words[c])]
    intervals=Dict{Tuple{String,Int},Vector{Vector{HCPlotPoint}}}()
    if g<=1
        isnothing(infinity_face) || throw(ArgumentError("the exterior disk chart requires genus at least two"))
        return _hp_square_sphere(lc,roles,curves,intervals,names;reference,with_contacts)
    end
    W=520.0*(g-1);H=520.0;mirror(p)=(p[1],H-p[2]);O=(W/2,H/2)
    kv(i,e)=first(abs(x) for x in words["l$i"] if abs(x) in abs.(words[e]))
    kvertices=Set(abs(x) for k in knames for x in words[k] if abs(x)<=3g-2)
    count_l(i)=count(x->!(abs(x) in kvertices),words["l$i"])
    A=max(130.0,min(220.0,40.0*max(1,maximum(count_l(i) for i in 1:g))+120))
    triangles=Vector{HCPlotPoint}[]
    pieces=Vector{Tuple{Int,Int,Vector{HCPlotPoint},String}}[]
    for i in 1:g
        x=520.0*(i-1)
        if i==1
            a=(x,H-A);b=(x+A,H);c=(x,H);push!(triangles,[c,a,b])
            vm,vn=kv(i,"m$i"),kv(i,"n$i")
            push!(pieces,[(vm,vn,[a,b],"F1"),(vn,vm,[b,c,a],"Finf")])
        elseif i==g
            a=(x,H-A);b=(x-A,H);c=(x,H);push!(triangles,[c,a,b])
            vm,vn=kv(i,"m$i"),kv(i,"n$(i-1)")
            push!(pieces,[(vn,vm,[b,a],"F$(g-1)"),(vm,vn,[a,c,b],"Finf")])
        else
            bl=(x-A,H);br=(x+A,H);a=(x,H-A*0.9);push!(triangles,[bl,br,a])
            vm,vl,vr=kv(i,"m$i"),kv(i,"n$(i-1)"),kv(i,"n$i")
            push!(pieces,[(vr,vl,[br,bl],"Finf"),(vl,vm,[bl,a],"F$(i-1)"),(vm,vr,[a,br],"F$i")])
        end
    end
    edges=Dict{String,Vector{Tuple{HCPlotPoint,HCPlotPoint}}}();sides=Dict{String,Tuple{String,String}}()
    enames=vcat(["m$i" for i in 1:g],["n$i" for i in 1:g-1])
    for i in 1:g
        p=first(p for p in triangles[i] if p[2]<H)
        edges["m$i"]=[(p,mirror(p))];sides["m$i"]=(i<g ? "F$i" : "Finf", i>1 ? "F$(i-1)" : "Finf")
    end
    for i in 1:g-1
        a=last(sort!([p for p in triangles[i] if p[2]==H];by=first))
        b=first(sort!([p for p in triangles[i+1] if p[2]==H];by=first))
        edges["n$i"]=[(a,b),(mirror(b),mirror(a))];sides["n$i"]=("Finf","F$i")
    end
    pos=Dict{Int,HCPlotPoint}();lface=Dict{Int,String}()
    reference_edges=NamedTuple[]
    signed(c,l)=only(x for x in words[c] if abs(x)==l)
    for e in enames
        w=words[e];i0=findfirst(x->x>0 && x in kvertices,w);labs=abs.(vcat(w[i0:end],w[1:i0-1]))
        starts=findall(in(kvertices),labs)
        length(starts)==length(edges[e]) || throw(ArgumentError("K edge partition: $e"))
        for (vi,(p,q)) in enumerate(edges[e])
            firstidx=starts[vi];lastidx=vi==length(starts) ? length(labs)+1 : starts[vi+1]
            labels=labs[firstidx+1:lastidx-1]
            reference && push!(reference_edges,(component=e,index=vi,p=p,q=q,
                first=labs[firstidx],last=labs[mod1(lastidx,length(labs))],word=[signed(e,l) for l in labels]))
            _hp_store!(intervals,e,[signed(e,l) for l in vcat(labs[firstidx],labels)],[p,q])
            for (j,l) in enumerate(labels);pos[l]=_hp_lerp(p,q,j/(length(labels)+1));end
        end
    end
    for i in 1:g
        w=abs.(words["l$i"])
        for (v0,v1,path,face) in pieces[i]
            a=findfirst(==(v0),w);ww=vcat(w[a:end],w[1:a-1]);b=findfirst(==(v1),ww);labels=ww[2:b-1]
            sig=[signed("l$i",l) for l in vcat(v0,labels)]
            _hp_store!(intervals,"l$i",sig,path);_hp_store!(intervals,"l$i",sig,mirror.(path))
            total=sum(_hp_dist(path[k],path[k+1]) for k in 1:length(path)-1)
            for (j,l) in enumerate(labels)
                t=j/(length(labels)+1)*total
                pos[l]=first(_hp_path_slice(path,t,t));lface[l]=face
            end
        end
    end
    labedge=Dict(abs(x)=>e for e in knames for x in words[e] if !(abs(x) in kvertices))
    perimeter=2(W+H)
    function fraction(p)
        x,y=p
        abs(y)<1e-6 && return x/perimeter
        abs(x-W)<1e-6 && return (W+y)/perimeter
        abs(y-H)<1e-6 && return (W+H+(W-x))/perimeter
        (2W+H+(H-y))/perimeter
    end
    outline(f)=_hp_outline(f,W,H)
    arcs=HCPlotArc[]
    for c in curves
        w=words[c];n=length(w);indices=findall(x->haskey(labedge,abs(x)),w)
        isempty(indices) && throw(ArgumentError("curve without K crossings: $c"))
        for k in eachindex(indices)
            a=indices[k];b=k==length(indices) ? indices[1]+n : indices[k+1]
            s,e=abs(w[a]),abs(w[mod1(b,n)]);es,ee=labedge[s],labedge[e]
            face=startswith(es,"l") ? lface[s] : sides[es][signed(es,s)>0 ? 1 : 2]
            before=startswith(ee,"l") ? lface[e] : sides[ee][signed(ee,e)>0 ? 2 : 1]
            face==before || throw(ArgumentError("face mismatch $c $s $e: $face / $before"))
            p=startswith(es,"l") && signed(c,s)<0 ? mirror(pos[s]) : pos[s]
            q=startswith(ee,"l") && signed(c,e)>0 ? mirror(pos[e]) : pos[e]
            model(p)=face=="Finf" ? (cos(2pi*fraction(p)),sin(2pi*fraction(p))) : p
            push!(arcs,HCPlotArc(c,s,e,[abs(w[mod1(j,n)]) for j in a+1:b-1],face,p,q,model(p),model(q)))
        end
    end
    faces=vcat(["F$i" for i in 1:g-1],"Finf")
    boundaries=Dict("Finf"=>HCPlotPoint[(1,0),(0,1),(-1,0),(0,-1)])
    for i in 1:g-1
        p=edges["m$i"][1][1];q=edges["m$(i+1)"][1][1];a,b=edges["n$i"][1]
        boundaries["F$i"]=[p,mirror(p),mirror(a),mirror(b),mirror(q),q,b,a]
    end
    methods,meshes=_hp_faces!(arcs,words,faces,boundaries)
    fsegs=[(a.model[k],a.model[k+1]) for a in arcs if a.face=="Finf" for k in 1:length(a.model)-1]
    best=-Inf;CC=(0.0,0.0)
    if methods["Finf"]==:tutte
        mesh=meshes["Finf"]
        candidates=sort!([(sum(mesh.points[v][1] for v in t)/3,sum(mesh.points[v][2] for v in t)/3) for t in mesh.triangles])
    else
        candidates=HCPlotPoint[(-1+2*(i+.5)/160,-1+2*(j+.5)/160) for i in 0:159 for j in 0:159 if hypot(-1+2*(i+.5)/160,-1+2*(j+.5)/160)<.97]
    end
    if !isnothing(infinity_face)
        candidates = _hp_infinity_candidates(arcs, words, candidates, infinity_face)
    end
    for p in candidates
        d=1-_hp_norm(p)
        for (a,b) in fsegs;d=min(d,_hp_nearest(p,a,b)[1]);end
        if d>best;best=d;CC=p;end
    end
    delta=best
    function phi(p)
        v=_hp_sub(p,CC);nv=_hp_norm(v);u=(v[1]/nv,v[2]/nv)
        bq=_hp_dot(CC,u);cq=_hp_dot(CC,CC)-1;t=-bq+sqrt(bq*bq-cq)
        q=_hp_add(CC,_hp_mul(u,t));s=nv/t
        P=outline(mod(atan(q[2],q[1]),2pi)/(2pi));d=_hp_dist(P,O);uu=((P[1]-O[1])/d,(P[2]-O[2])/d)
        f=s<1 ? 1+0.5*(1-s)/(1-delta) : 1.0
        _hp_add(O,_hp_mul(uu,d*f))
    end
    map_point(face,p)=face=="Finf" ? phi(p) : p
    for a in arcs
        if a.face=="Finf"
            a.poly=[phi(_hp_lerp(a.model[i],a.model[i+1],k/300)) for i in 1:length(a.model)-1 for k in 0:299]
            push!(a.poly,phi(a.model[end]));a.poly[1]=a.p;a.poly[end]=a.q
        else;a.poly=copy(a.model);end
    end
    holetri=Pair{String,Vector{HCPlotPoint}}[];holes=Pair{String,Tuple{HCPlotPoint,Float64}}[]
    for i in 1:g,(name,t) in (("top$i",triangles[i]),("bot$i",mirror.(triangles[i])))
        push!(holetri,name=>t)
        a=_hp_dist(t[2],t[3]);b=_hp_dist(t[1],t[3]);c=_hp_dist(t[1],t[2]);per=a+b+c
        center=((a*t[1][1]+b*t[2][1]+c*t[3][1])/per,(a*t[1][2]+b*t[2][2]+c*t[3][2])/per)
        s=per/2;radius=sqrt(s*(s-a)*(s-b)*(s-c))/s;push!(holes,name=>(center,radius))
    end
    function extend(p)
        for ((_,t),(_, (center,r))) in zip(holetri,holes)
            minimum(_hp_nearest(p,t[i],t[mod1(i+1,3)])[1] for i in 1:3)<1e-6 || continue
            v=_hp_sub(p,center);n=_hp_dist(p,center)
            return (center[1]+r*v[1]/n,center[2]+r*v[2]/n)
        end
        throw(ArgumentError("endpoint is not on a cut triangle: $p"))
    end
    runs=HCPlotRun[];contacts=Dict{Tuple{String,Symbol},Int}()
    for c in curves
        ca=[a for a in arcs if a.curve==c];groups=[[ca[1]]]
        for a in ca[2:end]
            if _hp_dist(groups[end][end].poly[end],a.poly[1])<1e-6;push!(groups[end],a)
            else;push!(groups,[a]);end
        end
        if length(groups)>1 && _hp_dist(groups[end][end].poly[end],groups[1][1].poly[1])<1e-6
            groups[1]=vcat(pop!(groups),groups[1])
        end
        closed=length(groups)==1 && _hp_dist(groups[1][1].poly[1],groups[1][end].poly[end])<1e-6
        for (ri,group) in enumerate(groups)
            points=HCPlotPoint[]
            paths=[a.poly for a in group]
            if !closed;paths=vcat([[extend(group[1].poly[1])]],paths,[[extend(group[end].poly[end])]]);end
            for path in paths,p in path
                (isempty(points) || _hp_dist(points[end],p)>1e-6) && push!(points,p)
            end
            word=reference ? [signed(c,m) for a in group for m in vcat(startswith(labedge[a.first],"l") ? Int[] : [a.first],a.mids)] :
                [signed(c,m) for a in group for m in a.mids]
            runname=closed ? c : "$c#$(ri-1)"
            push!(runs,HCPlotRun(runname,c,get(roles,c,:none),closed,points,word))
            if !closed
                contacts[(runname,:start)]=group[1].first;contacts[(runname,:end)]=group[end].last
            end
        end
    end
    _hp_empty!(runs,lc,roles,arcs,intervals,map_point,[t for (_,t) in holetri],holes,W,H,false,names)
    placed=HCPlotPlacement(W,H,holes,runs,methods,CC,delta,false)
    reference || return with_contacts ? (placement=placed,contacts=contacts) : placed
    guides=HCPlotRun[]
    for e in reference_edges
        name="@K:"*e.component*":"*string(e.index)
        points=HCPlotPoint[]
        for q in (extend(e.p),e.p,e.q,extend(e.q))
            (isempty(points) || _hp_dist(points[end],q)>1e-9) && push!(points,q)
        end
        push!(guides,HCPlotRun(name,e.component,get(roles,e.component,:none),false,points,e.word))
        contacts[(name,:start)]=e.first;contacts[(name,:end)]=e.last
    end
    (placement=_hp_rebuild(placed,vcat(runs,guides)),contacts=contacts,fixed=Set(r.name for r in guides))
end

# cf: lib/topology.jl/scripts/heegaard/lc_square.py place_square — paired square edges
# Position of the j-th of n-1 ports on a square edge, as a fraction of the edge (j=0 and j=n are the corners).
# The ports sit at the midpoints of n-1 equal cells. With this spacing a curve of constant direction on the
# torus is drawn as parallel straight segments; the earlier spacing j/n gave segments of different slopes.
_hp_square_port_fraction(j,n)=j<=0 ? 0. : j>=n ? 1. : (j-.5)/(n-1)
# Side lengths of the genus 1 chart. A curve that crosses m1 p times and l1 q times per turn runs p widths
# across and q heights up, so with height/width = p/q a straight curve meets every edge at 45 degrees.
# The ratio is kept within 1/3..3 and the longer side is 520. The chart is still a fundamental domain
# of the torus; only its shape changes.
function _hp_square_sides(words)
    p=length(get(words,"m1",Int[]))-1;q=length(get(words,"l1",Int[]))-1
    ratio=p>0 && q>0 ? clamp(p/q,1/3,3.) : 1.
    ratio>=1 ? (520.0/ratio,520.0) : (520.0,520.0*ratio)
end
function _hp_square_sphere(lc,roles,curves,intervals,names;reference=false,with_contacts=false)
    g=lc.genus;words=lc.words;W,H=g==1 ? _hp_square_sides(words) : (520.0,520.0);arcs=HCPlotArc[];runs=HCPlotRun[]
    holes=Pair{String,Tuple{HCPlotPoint,Float64}}[];contacts=Dict{Tuple{String,Symbol},Int}()
    if g==0
        isempty(curves) || throw(ArgumentError("sphere crossing placement is not specified by the accepted placement"))
        methods=Dict("outer"=>:empty_circles)
    else
        ports=Dict{Int,Tuple{HCPlotPoint,HCPlotPoint}}()
        for name in ("l1","m1")
            w=words[name];k=findfirst(x->abs(x)==1,w);w=vcat(w[k:end],w[1:k-1]);n=length(w)
            copies=[[name=="l1" ? (W*_hp_square_port_fraction(j,n),side) : (side,H-H*_hp_square_port_fraction(j,n)) for j in 0:n] for side in (0.0,name=="l1" ? H : W)]
            for pts in copies,j in 1:n;push!(get!(intervals,(name,w[j]),Vector{HCPlotPoint}[]),pts[j:j+1]);end
            for j in 2:n;ports[abs(w[j])]=(copies[1][j],copies[2][j]);end
        end
        endpoint(x,outgoing)=ports[abs(x)][(x>0)==outgoing ? 2 : 1]
        for c in curves
            w=words[c];n=length(w);indices=findall(x->haskey(ports,abs(x)),w)
            isempty(indices) && throw(ArgumentError("curve without K crossings: $c"))
            for k in eachindex(indices)
                a=indices[k];b=k==length(indices) ? indices[1]+n : indices[k+1]
                p=endpoint(w[a],true);q=endpoint(w[mod1(b,n)],false)
                push!(arcs,HCPlotArc(c,abs(w[a]),abs(w[mod1(b,n)]),[abs(w[mod1(j,n)]) for j in a+1:b-1],"F1",p,q))
            end
        end
        methods,_=_hp_faces!(arcs,words,["F1"],Dict("F1"=>HCPlotPoint[(0,0),(W,0),(W,H),(0,H)]))
        for c in curves,(j,a) in enumerate([a for a in arcs if a.curve==c])
            a.poly=copy(a.model)
            word=[only(x for x in words[c] if abs(x)==m) for m in a.mids]
            name="$c#$(j-1)"
            push!(runs,HCPlotRun(name,c,get(roles,c,:none),false,copy(a.poly),word))
            contacts[(name,:start)]=a.first;contacts[(name,:end)]=a.last
        end
    end
    _hp_empty!(runs,lc,roles,arcs,intervals,(f,p)->p,Vector{HCPlotPoint}[],holes,W,H,g==1,names)
    placed=HCPlotPlacement(W,H,holes,runs,methods,nothing,0.0,g==1)
    reference ? (placement=placed,contacts=contacts,fixed=Set{String}()) : with_contacts ? (placement=placed,contacts=contacts) : placed
end
