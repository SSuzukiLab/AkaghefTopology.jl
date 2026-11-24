# cf: lib/topology.jl/docs/specs/lickorish-code.md §6 — cut-disk geometry in screen coordinates
const HCPlotPoint = NTuple{2,Float64}
_hp_add(a,b) = (a[1]+b[1],a[2]+b[2])
_hp_sub(a,b) = (a[1]-b[1],a[2]-b[2])
_hp_mul(a,t) = (a[1]*t,a[2]*t)
_hp_dot(a,b) = a[1]*b[1]+a[2]*b[2]
_hp_cross(a,b) = a[1]*b[2]-a[2]*b[1]
_hp_norm(a) = hypot(a[1],a[2])
_hp_dist(a,b) = _hp_norm(_hp_sub(a,b))
_hp_lerp(a,b,t) = _hp_add(a,_hp_mul(_hp_sub(b,a),t))
# cf: scripts/heegaard/lc_placement.py outline_pt — all four rectangle sides
function _hp_outline(f,W,H)
    t=mod(f,1)*2(W+H)
    t<W && return (t,0.0);t-=W
    t<H && return (W,t);t-=H
    t<W && return (W-t,H);t-=W
    (0.0,H-t)
end
function _hp_nearest(p,a,b)
    d=_hp_sub(b,a); den=_hp_dot(d,d)
    t=den==0 ? 0.0 : clamp(_hp_dot(_hp_sub(p,a),d)/den,0.0,1.0)
    q=_hp_lerp(a,b,t)
    (_hp_dist(p,q),q,t)
end
function _hp_intersection(a,b,c,d; halfopen=false)
    r=_hp_sub(b,a); s=_hp_sub(d,c); den=_hp_cross(r,s)
    abs(den)<1e-12 && return nothing
    q=_hp_sub(c,a); u=_hp_cross(q,s)/den; v=_hp_cross(q,r)/den
    good=halfopen ? (0<=u<1 && 0<=v<1) : (0<u<1 && 0<v<1)
    good ? (u,v,Int(sign(den))) : nothing
end
mutable struct HCPlotArc
    curve::String
    first::Int
    last::Int
    mids::Vector{Int}
    face::String
    p::HCPlotPoint
    q::HCPlotPoint
    mp::HCPlotPoint
    mq::HCPlotPoint
    model::Vector{HCPlotPoint}
    events::Vector{HCPlotPoint}
    eventpaths::Vector{Vector{HCPlotPoint}}
    poly::Vector{HCPlotPoint}
end
HCPlotArc(c,s,e,m,f,p,q,mp=p,mq=q)=HCPlotArc(c,s,e,m,f,p,q,mp,mq,HCPlotPoint[],HCPlotPoint[],Vector{HCPlotPoint}[],HCPlotPoint[])
struct HCPlotRun
    name::String
    component::String
    role::Symbol
    closed::Bool
    points::Vector{HCPlotPoint}
    word::Vector{Int}
end
struct HCPlotPlacement
    width::Float64
    height::Float64
    holes::Vector{Pair{String,Tuple{HCPlotPoint,Float64}}}
    runs::Vector{HCPlotRun}
    methods::Dict{String,Symbol}
    center::Union{Nothing,HCPlotPoint}
    delta::Float64
    square::Bool
end

# A display parameter derived only from the cut-hole geometry, fixed throughout
# a run. Twice the enclosing hole extent leaves an exterior annulus of that
# extent. This is a canvas choice, not a per-iterate expansion to hide intrusion.
function _hp_frame(p::HCPlotPlacement)
    center=(p.width/2,p.height/2)
    extent=isempty(p.holes) ? hypot(p.width/2,p.height/2) : maximum(_hp_dist(c,center)+r for (_, (c,r)) in p.holes)
    center,2extent+60
end

function _hp_path_slice(path,a,b)
    lengths=[_hp_dist(path[i],path[i+1]) for i in 1:length(path)-1]
    stations=vcat(0.0,cumsum(lengths))
    function at(t)
        for i in eachindex(lengths)
            t<=stations[i]+lengths[i]+1e-9 && return _hp_lerp(path[i],path[i+1],(t-stations[i])/lengths[i])
        end
        path[end]
    end
    vcat([at(a)],[path[i] for i in 2:length(path)-1 if a<stations[i]<b],[at(b)])
end
function _hp_store!(intervals,name,labels,path)
    total=sum(_hp_dist(path[i],path[i+1]) for i in 1:length(path)-1)
    n=length(labels)
    for (j,label) in enumerate(labels)
        push!(get!(intervals,(name,label),Vector{HCPlotPoint}[]),_hp_path_slice(path,total*(j-1)/n,total*j/n))
    end
end

# cf: lib/topology.jl/docs/heegaard/261109_tutte_construction.md — occurrence subdivision; same uniform harmonic system as lc_tutte.py
function _hp_tutte!(arcs,words,corners; exterior=false)
    fixed=Union{Nothing,HCPlotPoint}[]
    vertex(p=nothing)=(push!(fixed,p);length(fixed))
    endpoints=[(vertex(a.mp),vertex(a.mq)) for a in arcs]
    boundary=Int[v for pair in endpoints for v in pair]
    for p in corners
        any(v->_hp_dist(p,fixed[v])<1e-12,boundary) || push!(boundary,vertex(p))
    end
    center=(sum(p[1] for p in corners)/length(corners),sum(p[2] for p in corners)/length(corners))
    sort!(boundary;by=v->atan(fixed[v][2]-center[2],fixed[v][1]-center[1]),alg=MergeSort)
    crossings=Dict(m=>vertex() for m in sort!(unique([m for a in arcs for m in a.mids])))
    edges=Tuple{Int,Int}[]; rotation=Dict{Int,Vector{Int}}()
    function edge(u,v)
        d=2length(edges)+1;push!(edges,(u,v))
        push!(get!(rotation,u,Int[]),d);push!(get!(rotation,v,Int[]),d+1);d
    end
    flip(d)=isodd(d) ? d+1 : d-1
    target(d)=edges[(d+1)÷2][isodd(d) ? 2 : 1]
    origin(d)=edges[(d+1)÷2][isodd(d) ? 1 : 2]
    bd=[edge(boundary[i],boundary[mod1(i+1,length(boundary))]) for i in eachindex(boundary)]
    visits=Dict{Tuple{Int,Bool},Tuple{Int,Int}}(); chains=Vector{Int}[]
    for (a,(s,t)) in zip(arcs,endpoints)
        original=vcat(s,[crossings[m] for m in a.mids],t); chain=[s]
        start=findfirst(x->abs(x)==a.first,words[a.curve])
        for i in 1:length(original)-1
            u,v=original[i:i+1];x,y=vertex(),vertex()
            edge(u,x);edge(x,y);edge(y,v);append!(chain,(x,y,v))
        end
        push!(chains,chain)
        for (k,m) in enumerate(a.mids)
            v=crossings[m];incoming=only(d for d in rotation[v] if target(d)==chain[3k])
            outgoing=only(d for d in rotation[v] if target(d)==chain[3k+2])
            signed=words[a.curve][mod1(start+k,length(words[a.curve]))]
            visits[(m,signed>0)]=(outgoing,incoming)
        end
    end
    for (m,v) in crossings
        po,pi=visits[(m,true)];no,ni=visits[(m,false)]
        rotation[v]=exterior ? [po,ni,pi,no] : [po,no,pi,ni]
    end
    for (i,v) in enumerate(boundary)
        forward=bd[i];backward=flip(bd[mod1(i-1,length(bd))])
        rotation[v]=vcat(forward,[d for d in rotation[v] if d!=forward && d!=backward],backward)
    end
    successor=zeros(Int,2length(edges))
    for darts in values(rotation),i in eachindex(darts)
        successor[flip(darts[i])]=darts[mod1(i-1,length(darts))]
    end
    cells=Vector{Int}[];seen=falses(length(successor));outer=flip(bd[1])
    for start in eachindex(seen)
        seen[start] && continue
        walk=Int[];d=start
        while !seen[d]
            seen[d]=true;push!(walk,d);d=successor[d]
        end
        outer in walk || push!(cells,walk)
    end
    length(fixed)-length(edges)+length(cells)==1 || throw(ArgumentError("cut face is not a disk"))
    triangles=NTuple{3,Int}[]
    for walk in cells
        ring=origin.(walk);inner=[vertex() for _ in ring];hub=vertex()
        for i in eachindex(ring)
            j=mod1(i+1,length(ring));a,b=ring[i],ring[j];p,q=inner[i],inner[j]
            append!(triangles,[(a,b,q),(a,q,p),(p,q,hub)])
        end
    end
    neighbors=[Set{Int}() for _ in fixed]
    for t in triangles,i in 1:3
        a,b=t[i],t[mod1(i+1,3)];push!(neighbors[a],b);push!(neighbors[b],a)
    end
    free=findall(isnothing,fixed);indices=Dict(v=>i for (i,v) in enumerate(free))
    matrix=zeros(length(free),length(free));rhs=zeros(length(free),2)
    for v in free
        i=indices[v];matrix[i,i]=length(neighbors[v])
        for w in sort!(collect(neighbors[v]))
            if isnothing(fixed[w]);matrix[i,indices[w]]-=1
            else;rhs[i,1]+=fixed[w][1];rhs[i,2]+=fixed[w][2];end
        end
    end
    solution=matrix\rhs
    points=HCPlotPoint[isnothing(p) ? (solution[indices[v],1],solution[indices[v],2]) : p for (v,p) in enumerate(fixed)]
    for (a,chain) in zip(arcs,chains)
        a.model=points[chain];a.events=a.model[1:3:end]
        a.eventpaths=[a.model[i:i+3] for i in 1:3:length(chain)-3]
    end
    (points=points,triangles=triangles)
end

function _hp_faces!(arcs,words,faces,boundaries)
    methods=Dict{String,Symbol}(); meshes=Dict{String,Any}()
    for face in faces
        fa=[a for a in arcs if a.face==face];found=[Tuple{Float64,Int}[] for _ in fa]
        for i in 1:length(fa),j in i+1:length(fa)
            hit=_hp_intersection(fa[i].mp,fa[i].mq,fa[j].mp,fa[j].mq)
            isnothing(hit) && continue
            push!(found[i],(hit[1],j));push!(found[j],(hit[2],i))
        end
        ok=true
        # A chord between two ports on the same straight side of the face lies
        # on that side: it would be drawn on top of the K edge and of every
        # other such chord. This happens when a curve meets the same K edge
        # twice in a row. Such a face is laid out by the Tutte embedding.
        # The outer face is drawn in a disk model whose chords never lie on
        # the boundary circle, so the test applies to the polygonal faces only.
        corners=face=="Finf" ? HCPlotPoint[] : boundaries[face]
        for a in fa,k in eachindex(corners)
            p,q=corners[k],corners[mod1(k+1,length(corners))]
            ok &= !(_hp_nearest(a.mp,p,q)[1]<1e-9 && _hp_nearest(a.mq,p,q)[1]<1e-9)
        end
        for (i,a) in enumerate(fa)
            sort!(found[i];by=first,alg=MergeSort)
            expected=[only(c for c in keys(words) if c!=a.curve && any(x->abs(x)==m,words[c])) for m in a.mids]
            ok &= [fa[j].curve for (_,j) in found[i]]==expected
        end
        if ok
            methods[face]=:chords
            for (i,a) in enumerate(fa)
                a.model=[a.mp,a.mq];a.events=vcat([a.mp],[_hp_lerp(a.mp,a.mq,t) for (t,_) in found[i]],[a.mq])
                a.eventpaths=[[a.events[j],a.events[j+1]] for j in 1:length(a.events)-1]
            end
        else
            methods[face]=:tutte;meshes[face]=_hp_tutte!(fa,words,boundaries[face];exterior=face=="Finf")
        end
    end
    methods,meshes
end


# cf: lickorish-code.md (face orbits of ρ ∘ τ) and heegaard-solid-view.md §2.2.
# In the exterior disk model, ρ is counterclockwise. The face traversing an
# edge forward is therefore on its right. Classify a pole by a nearest
# interior edge point: the open segment to this foot crosses no other edge.
function _hp_infinity_candidates(arcs, words, candidates, face)
    _, corner_face = _lc_faces(words)
    faces = Set(values(corner_face))
    face in faces || throw(ArgumentError("unknown LC face $face"))
    sides = [(a=path[i], b=path[i+1], right=corner_face[(u, _lc_sign_char(x), "out")],
              left=corner_face[(v, _lc_sign_char(y), "in")])
             for arc in arcs if arc.face == "Finf"
             for (u,v,path) in zip([arc.first; arc.mids], [arc.mids; arc.last], arc.eventpaths)
             for x in (only(z for z in words[arc.curve] if abs(z)==u),)
             for y in (only(z for z in words[arc.curve] if abs(z)==v),)
             for i in 1:length(path)-1]
    # Seed both sides of every incident edge as well as the usual clearance
    # samples. Narrow faces can fall between all points of the regular grid.
    samples = copy(candidates)
    for edge in sides
        face in (edge.left,edge.right) || continue
        m = _hp_lerp(edge.a,edge.b,0.5)
        v = _hp_sub(edge.b,edge.a)
        sz = _hp_norm(v)
        sz > 1e-12 || continue
        gap = min(sz/2,1-_hp_norm(m))
        for other in sides
            other == edge && continue
            gap = min(gap,_hp_nearest(m,other.a,other.b)[1])
        end
        normal = (-v[2]/sz,v[1]/sz)
        for (id,sgn) in ((edge.left,1),(edge.right,-1))
            id == face && push!(samples,_hp_add(m,_hp_mul(normal,sgn*gap/2)))
        end
    end
    selected = filter(samples) do p
        isempty(sides) && return false
        _,i = findmin(e->_hp_nearest(p,e.a,e.b)[1],sides)
        edge = sides[i]
        d,_,t = _hp_nearest(p,edge.a,edge.b)
        d>1e-12 && 0<t<1 || return false
        cross = _hp_cross(_hp_sub(edge.b,edge.a),_hp_sub(p,edge.a))
        (cross>0 ? edge.left : edge.right) == face
    end
    isempty(selected) && throw(ArgumentError("LC face $face has no pole in the exterior chart"))
    selected
end
