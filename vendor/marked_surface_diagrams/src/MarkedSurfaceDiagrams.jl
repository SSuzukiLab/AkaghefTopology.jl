module MarkedSurfaceDiagrams

#= Marked oriented PL diagrams. Standalone Julia >=1.10; TOML and SHA stdlibs.
The reference implementation and shared fixtures are in tools/reference.py.
This source has NOT been executed in the delivery environment: run test/runtests.jl.
=#
using TOML, SHA
export DiagramError, Dart, DiagramState, load_diagram, validate, parse_word, print_word,
       realize, certify_layout, seam_word, endpoint_anchors, intersection_cycles,
       regions, check_joining_disk, family_cut_report, move_vertex, write_svg,
       representation_digest, marking_digest, report, write_json, legacy_projection

const Rat = Rational{BigInt}
const Point = NTuple{2,Rat}
struct DiagramError <: Exception
    code::String
    detail::String
end
Base.showerror(io::IO,e::DiagramError)=print(io,e.code,": ",e.detail)
need(ok,code,detail) = ok ? nothing : throw(DiagramError(code,string(detail)))
edge(u,v)=isless(u,v) ? (u,v) : (v,u)
tedges(t)=[(t[i],t[mod1(i+1,3)]) for i in 1:3]
function components(nodes,adj)
    unseen=Set(nodes);out=[]
    while !isempty(unseen)
        seed=minimum(unseen);delete!(unseen,seed);comp=Set([seed]);stack=[seed]
        while !isempty(stack)
            v=pop!(stack)
            for w in get(adj,v,[])
                if w in unseen
                    delete!(unseen,w);push!(comp,w);push!(stack,w)
                end
            end
        end
        push!(out,comp)
    end
    out
end
mutable struct UnionFind{T}
    p::Dict{T,T}
end
UnionFind(nodes)=UnionFind(Dict(v=>v for v in nodes))
function root!(d::UnionFind,x)
    p=d.p[x]
    if p!=x;d.p[x]=root!(d,p);end
    d.p[x]
end
function unite!(d::UnionFind,a,b)
    a=root!(d,a);b=root!(d,b)
    if a!=b;d.p[max(a,b)]=min(a,b);end
    nothing
end
struct Dart
    face::String
    u::String
    v::String
end
Base.:(==)(a::Dart,b::Dart)=(a.face,a.u,a.v)==(b.face,b.u,b.v)
struct DiagramState
    doc::Dict
    patches::Dict
    boundaries::Dict
    incidence::Dict
    triangles::Dict
    seam_partner::Dict
    vertex_dsu::UnionFind
    curves::Dict
    used::Dict
    crossing_at::Dict
    surface_report::Dict
    fingerprint::String
end
# Unicode word identifiers start with a non-digit word character; no implicit routing.
const BLOCK=r"^([^\W\d][\w.-]*)\{\s*([^\W\d][\w.-]*(?:\s+[^\W\d][\w.-]*)+)\s*\}$"
function parse_word(s::AbstractString)
    need(!isempty(strip(s)),"WORD_EMPTY",s);out=Dart[]
    for block in split(s,'|')
        m=match(BLOCK,strip(block));need(m!==nothing,"WORD_SYNTAX",block)
        f=String(m.captures[1]);vs=String.(split(m.captures[2]))
        append!(out,[Dart(f,u,v) for (u,v) in zip(vs[1:end-1],vs[2:end])])
    end
    out
end
function print_word(ds::Vector{Dart})
    blocks=Tuple{String,Vector{String}}[]
    for d in ds
        if !isempty(blocks) && blocks[end][1]==d.face && blocks[end][2][end]==d.u
            push!(blocks[end][2],d.v)
        else;push!(blocks,(d.face,[d.u,d.v]));end
    end
    join([f*"{"*join(vs," ")*"}" for (f,vs) in blocks]," | ")
end
function qnum(x)
    need((x isa Integer && !(x isa Bool)) || x isa AbstractString,"RATIONAL_LITERAL",x)
    x isa Integer && return BigInt(x)//BigInt(1)
    need(occursin(r"^-?[0-9]+(?:/[1-9][0-9]*)?$",x),"RATIONAL_LITERAL",x)
    z=split(x,'/');need(length(z) in (1,2),"RATIONAL_LITERAL",x)
    try
        length(z)==1 ? parse(BigInt,z[1])//BigInt(1) : parse(BigInt,z[1])//parse(BigInt,z[2])
    catch;throw(DiagramError("RATIONAL_LITERAL",string(x)));end
end
qstr(x::Rat)=denominator(x)==1 ? string(numerator(x)) : string(numerator(x),"/",denominator(x))
function side_parameters(p,sid)
    supplied=get(p,"side_parameters",Dict());n=length(p["sides"][sid])-1
    haskey(supplied,sid) ? qnum.(supplied[sid]) : [BigInt(k)//BigInt(n) for k in 0:n]
end
function patch_info(p)
    f=p["face"];order=p["side_order"];sides=p["sides"]
    need(length(order)>=3 && length(Set(order))==length(order),"PATCH_SIDES",f)
    need(Set(order)==Set(keys(sides)),"PATCH_SIDES",f);bd=String[]
    for (i,s) in enumerate(order)
        vs=sides[s]
        need(length(vs)>=2 && length(Set(vs))==length(vs),"SIDE_PATH","$f/$s")
        need(vs[end]==sides[order[mod1(i+1,length(order))]][1],"SIDE_CORNERS","$f/$s")
        append!(bd,vs[1:end-1])
    end
    need(length(Set(bd))==length(bd),"BOUNDARY_NOT_SIMPLE",f)
    inc=Dict{Tuple{String,String},Vector{Tuple{Int,String,String}}}();vertices=Set{String}();seen=Set()
    for (i,t) in enumerate(p["triangles"])
        need(length(t)==3 && length(Set(t))==3,"TRIANGLE_INVALID","$f/$i")
        key=Tuple(sort(t));need(!(key in seen),"TRIANGLE_DUPLICATE",f);push!(seen,key);union!(vertices,t)
        for (u,v) in tedges(t);push!(get!(inc,edge(u,v),Tuple{Int,String,String}[]),(i,u,v));end
    end
    need(issubset(Set(bd),vertices),"MISSING_BOUNDARY_VERTEX",f)
    bedges=Dict(edge(u,v)=>(u,v) for (u,v) in zip(bd,circshift(bd,-1)))
    for (e,occ) in inc
        if haskey(bedges,e)
            need(length(occ)==1 && occ[1][2:3]==bedges[e],"BOUNDARY_ORIENTATION","$f/$e")
        else
            need(length(occ)==2 && occ[1][2:3]==reverse(occ[2][2:3]),"EDGE_MANIFOLD","$f/$e")
        end
    end
    need(issubset(Set(keys(bedges)),Set(keys(inc))),"BOUNDARY_EDGE_MISSING",f)
    adj=Dict()
    for occ in values(inc)
        if length(occ)==2
            a,b=occ[1][1],occ[2][1];push!(get!(adj,a,Set()),b);push!(get!(adj,b,Set()),a)
        end
    end
    need(length(components(1:length(p["triangles"]),adj))==1,"PATCH_DISCONNECTED",f)
    need(length(vertices)-length(inc)+length(p["triangles"])==1,"PATCH_NOT_DISK",f)
    for v in vertices
        link=Dict()
        for t in p["triangles"]
            if v in t
                a,b=filter(!=(v),t);push!(get!(link,a,[]),b);push!(get!(link,b,[]),a)
            end
        end
        need(length(components(keys(link),link))==1,"VERTEX_LINK_DISCONNECTED","$f/$v")
        degrees=length.(collect(values(link)))
        if v in bd
            need(count(==(1),degrees)==2 && all(x->x in (1,2),degrees),"BOUNDARY_LINK","$f/$v")
        else;need(all(==(2),degrees),"INTERIOR_LINK","$f/$v");end
    end
    bd,inc,vertices
end

function check_keys(x,required,optional=())
    need(x isa AbstractDict,"SCHEMA_OBJECT",x)
    need(issubset(Set(required),Set(keys(x))) && issubset(Set(keys(x)),union(Set(required),Set(optional))),"SCHEMA_KEYS",collect(keys(x)))
end
function check_shape(doc)
    check_keys(doc,("version","surface","patches","curves","crossings"),("metadata","region_marks"))
    check_keys(doc["surface"],("marking_id","genus","gluings","boundary"))
    for p in doc["patches"];check_keys(p,("face","side_order","sides","triangles"),("side_parameters",));end
    for c in doc["curves"];check_keys(c,("id","role","closed","word"));end
    for x in doc["crossings"];check_keys(x,("id","at","curves","sign"));end
    for b in doc["surface"]["boundary"];check_keys(b,("side","label"));end
    for r in get(doc,"region_marks",[])
        check_keys(r,("id","seed","claim"));z=r["seed"]
        need(z isa AbstractVector && length(z)==2 && z[1] isa AbstractString && z[2] isa Integer && !(z[2] isa Bool) && z[2]>=0,"REGION_SEED",z)
    end
end
check_state(st)=need(representation_digest(st.doc)==st.fingerprint,"MATHEMATICAL_STATE_MUTATED","revalidate edited document")
function _validate(doc::Dict)
    check_shape(doc);doc=deepcopy(doc)
    need(doc["version"] isa Integer && !(doc["version"] isa Bool) && doc["version"]==1,"VERSION","expected version=1")
    surf=doc["surface"];need(!isempty(get(surf,"marking_id","")),"MARKING_REQUIRED","surface.marking_id")
    patches=Dict(p["face"]=>p for p in doc["patches"])
    need(length(patches)==length(doc["patches"]),"DUPLICATE_FACE","patches")
    need(!isempty(patches),"SURFACE_EMPTY","patches");need(length(patches)<=64,"RESOURCE_LIMIT",">64 patches")
    bds=Dict();incs=Dict();nodes=Tuple{String,String}[];tri=Dict();sidepaths=Dict();corners=Set{Tuple{String,String}}()
    for f in sort(collect(keys(patches)))
        p=patches[f];bd,inc,vs=patch_info(p);bds[f]=bd;incs[f]=inc
        need(length(vs)<=256,"RESOURCE_LIMIT","$f: >256 vertices in reference exact kernel")
        supplied=get(p,"side_parameters",Dict());need(issubset(Set(keys(supplied)),Set(keys(p["sides"]))),"SIDE_PARAMETERS",f)
        for (sid,path) in p["sides"]
            ts=side_parameters(p,sid)
            need(length(ts)==length(path) && ts[1]==0 && ts[end]==1 && all(ts[i]<ts[i+1] for i in 1:length(ts)-1),"SIDE_PARAMETERS","$f/$sid")
        end
        append!(nodes,[(f,v) for v in vs])
        for (i,t) in enumerate(p["triangles"]);tri[(f,i)]=Tuple(t);end
        for sid in p["side_order"]
            sidepaths[(f,sid)]=p["sides"][sid];push!(corners,(f,p["sides"][sid][1]))
        end
    end
    dsu=UnionFind(nodes);faceadj=Dict();paired=Set();seam=Dict();corneradj=Dict()
    for ab in get(surf,"gluings",[])
        a,b=Tuple(ab[1]),Tuple(ab[2])
        need(haskey(sidepaths,a) && haskey(sidepaths,b) && a!=b,"GLUING_SIDE","$a,$b")
        need(!(a in paired) && !(b in paired),"GLUING_DUPLICATE","$a,$b");push!(paired,a,b)
        va,vb=sidepaths[a],sidepaths[b];need(length(va)==length(vb),"GATE_COUNT","$a,$b")
        pa,pb=side_parameters(patches[a[1]],a[2]),side_parameters(patches[b[1]],b[2])
        need(pa==[1-t for t in reverse(pb)],"GATE_PARAMETERS","$a,$b")
        for (u,v) in zip(va,reverse(vb));unite!(dsu,(a[1],u),(b[1],v));end
        for i in 1:length(va)-1
            j=length(vb)-i;ka=(a[1],edge(va[i],va[i+1]));kb=(b[1],edge(vb[j],vb[j+1]));seam[ka]=kb;seam[kb]=ka
        end
        for (u,v) in [(va[1],vb[end]),(va[end],vb[1])]
            ca,cb=(a[1],u),(b[1],v);push!(get!(corneradj,ca,[]),cb);push!(get!(corneradj,cb,[]),ca)
        end
        push!(get!(faceadj,a[1],Set()),b[1]);push!(get!(faceadj,b[1],Set()),a[1])
    end
    need(length(components(keys(patches),faceadj))==1,"SURFACE_DISCONNECTED","base cells")
    unpaired=setdiff(Set(keys(sidepaths)),paired)
    labels=Dict(Tuple(x["side"])=>x["label"] for x in get(surf,"boundary",[]))
    need(Set(keys(labels))==unpaired && length(labels)==length(get(surf,"boundary",[])),"BOUNDARY_LABELS","label each unglued side exactly once")
    groups=Dict()
    for c in corners;push!(get!(groups,root!(dsu,c),Set()),c);end
    for (rt,cs) in groups
        need(length(components(cs,corneradj))==1,"SURFACE_VERTEX_LINK",rt)
        degrees=[length(get(corneradj,c,[])) for c in cs];defect=sum(2-x for x in degrees)
        need(all(x<=2 for x in degrees) && defect in (0,2),"SURFACE_VERTEX_LINK",rt)
    end
    badj=Dict();bv=Set();sideends=Dict()
    for s in unpaired
        a,b=root!(dsu,(s[1],sidepaths[s][1])),root!(dsu,(s[1],sidepaths[s][end]))
        push!(get!(badj,a,[]),b);push!(get!(badj,b,[]),a);push!(bv,a,b);sideends[s]=(a,b)
    end
    need(all(length(badj[v])==2 for v in bv),"SURFACE_BOUNDARY_LINK","boundary")
    bc=components(bv,badj);usedlabels=[]
    for c in bc
        names=Set(labels[s] for (s,ends) in sideends if ends[1] in c)
        need(length(names)==1,"BOUNDARY_COMPONENT_LABEL","one label per component");append!(usedlabels,names)
    end
    need(length(Set(usedlabels))==length(usedlabels),"BOUNDARY_LABEL_REUSED","distinct components")
    V=length(groups);E=length(sidepaths)-length(paired)÷2;F=length(patches);chi=V-E+F;g=surf["genus"]
    need(g isa Integer && !(g isa Bool) && g>=0 && chi==2-2*g-length(bc),"SURFACE_GENUS","chi=$chi,g=$g,b=$(length(bc))")
    sr=Dict("genus"=>g,"boundary_components"=>length(bc),"chi"=>chi,"faces"=>F)
    curves=Dict();used=Dict();visits=Dict();gv(f,v)=root!(dsu,(f,v))
    physical=Set(root!(dsu,(s[1],v)) for s in unpaired for v in sidepaths[s][2:end-1])
    for c in get(doc,"curves",[])
        cid=c["id"];need(c["closed"] isa Bool,"CURVE_CLOSED","closed must be boolean")
        need(!haskey(curves,cid),"CURVE_DUPLICATE",cid);need(c["role"] in ("alpha","beta","aux"),"CURVE_ROLE",cid)
        ds=parse_word(c["word"]);curves[cid]=(c,ds)
        for d in ds
            need(haskey(incs,d.face) && haskey(incs[d.face],edge(d.u,d.v)),"ROUTE_EDGE_MISSING",d)
            need(length(incs[d.face][edge(d.u,d.v)])==2,"ROUTE_ON_MARKING",d)
            need(!((d.face,d.u) in corners) && !((d.face,d.v) in corners),"ROUTE_AT_CORNER",d)
            push!(get!(used,(d.face,edge(d.u,d.v)),[]),cid)
        end
        pairs=collect(zip(ds[1:end-1],ds[2:end]));if c["closed"];push!(pairs,(ds[end],ds[1]));end
        for (a,b) in pairs
            need(gv(a.face,a.v)==gv(b.face,b.u),"ROUTE_DISCONTINUITY",cid)
            if (a.face,a.v)!=(b.face,b.u)
                need(a.v in bds[a.face] && b.u in bds[b.face],"ILLEGAL_CHART_JUMP",cid)
            else;need(!(a.v in bds[a.face]),"GATE_TURNBACK","$cid/$(a.v)");end
        end
        if !c["closed"]
            need(gv(ds[1].face,ds[1].u) in physical && gv(ds[end].face,ds[end].v) in physical,"ARC_ENDPOINT",cid)
        end
        ordered=[gv(d.face,d.u) for d in ds];if !c["closed"];push!(ordered,gv(ds[end].face,ds[end].v));end
        need(length(ordered)==length(Set(ordered)),"SELF_INTERSECTION",cid)
        for (a,b) in pairs
            rt=gv(a.face,a.v);at=get!(visits,rt,Dict());need(!haskey(at,cid),"CURVE_REPEATED_VERTEX",cid);at[cid]=(a,b)
        end
        if !c["closed"]
            get!(visits,ordered[1],Dict())[cid]=(nothing,ds[1]);get!(visits,ordered[end],Dict())[cid]=(ds[end],nothing)
        end
    end
    need(all(length(v)==1 for v in values(used)),"SHARED_ROUTE_EDGE","curves overlap along an edge")
    declared=Dict();crossingids=Set()
    for x in get(doc,"crossings",[])
        f,v=x["at"];key=(f,v);need(haskey(dsu.p,key),"CROSSING_VERTEX",key)
        need(!(v in bds[f]),"CROSSING_ON_MARKING",key)
        need(!(x["id"] in crossingids) && !haskey(declared,key),"CROSSING_DUPLICATE",key)
        push!(crossingids,x["id"]);declared[key]=x
    end
    seen=Set()
    for (rt,at) in visits
        if length(at)==1
            need(!haskey(declared,rt),"FALSE_CROSSING",rt);continue
        end
        need(length(at)==2,"TRIPLE_POINT",rt);need(haskey(declared,rt),"UNDECLARED_INTERSECTION",rt)
        x=declared[rt];push!(seen,rt);first,second=x["curves"]
        need(Set([first,second])==Set(keys(at)),"CROSSING_CURVES",x["id"])
        r1,r2=curves[first][1]["role"],curves[second][1]["role"]
        need(!(r1==r2 && r1 in ("alpha","beta")),"SAME_FAMILY_INTERSECTION",x["id"])
        a0,a1=at[first];b0,b1=at[second];need(all(d!==nothing for d in (a0,a1,b0,b1)),"ENDPOINT_INTERSECTION",x["id"])
        links=Dict()
        for t in patches[rt[1]]["triangles"]
            if rt[2] in t;i=findfirst(==(rt[2]),t);links[t[mod1(i+1,3)]]=t[mod1(i+2,3)];end
        end
        ao,ai,bo,bi=a1.v,a0.u,b1.v,b0.u;order=String[];v=ao
        while isempty(order) || v!=ao
            push!(order,v);v=links[v];need(length(order)<=length(links),"CROSSING_LINK",x["id"])
        end
        reduced=filter(v->v in (ao,ai,bo,bi),order)
        need(reduced in ([ao,bo,ai,bi],[ao,bi,ai,bo]),"NONTRANSVERSE_INTERSECTION",x["id"])
        sign=reduced==[ao,bo,ai,bi] ? 1 : -1
        need(x["sign"] isa Integer && !(x["sign"] isa Bool) && x["sign"]==sign,"CROSSING_SIGN","$(x["id"]): expected $sign")
    end
    need(seen==Set(keys(declared)),"FALSE_CROSSING","unused crossing declaration")
    st=DiagramState(doc,patches,bds,incs,tri,seam,dsu,curves,used,declared,sr,representation_digest(doc))
    ids=Set()
    for mark in get(doc,"region_marks",[])
        need(!(mark["id"] in ids),"REGION_MARK_DUPLICATE",mark["id"]);push!(ids,mark["id"])
        need(mark["claim"] in ("region","embedded_closed_disk","koda_joining"),"REGION_CLAIM",mark["claim"])
        audit=check_joining_disk(st,mark["seed"])
        need(mark["claim"]=="region" || audit[mark["claim"]],"REGION_CLAIM_FAILED",mark["id"])
    end
    st
end
function validate(doc::Dict)
    try
        _validate(doc)
    catch e
        e isa DiagramError && rethrow()
        if e isa KeyError || e isa MethodError || e isa BoundsError || e isa ArgumentError || e isa TypeError
            throw(DiagramError("SCHEMA_SHAPE",sprint(showerror,e)))
        end
        rethrow()
    end
end
load_diagram(path::AbstractString)=validate(TOML.parsefile(path))

orient(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
onsegment(a,b,p)=orient(a,b,p)==0 && min(a[1],b[1])<=p[1]<=max(a[1],b[1]) && min(a[2],b[2])<=p[2]<=max(a[2],b[2])
function segments_meet(a,b,c,d)
    s,t,u,v=orient(a,b,c),orient(a,b,d),orient(c,d,a),orient(c,d,b)
    (s*t<0 && u*v<0) || any(z==0 && onsegment(p,q,r) for (z,p,q,r) in [(s,a,b,c),(t,a,b,d),(u,c,d,a),(v,c,d,b)])
end
function circle_boundary(n)
    out=Point[]
    for i in 0:n-1
        quadrant=4*i÷n;t=BigInt(mod(4*i,n))//BigInt(n);x,y=(1-t*t)/(1+t*t),2*t/(1+t*t)
        for _ in 1:quadrant;x,y=-y,x;end
        push!(out,(x,y))
    end
    out
end
function solve_exact(A::Matrix{Rat},B::Matrix{Rat})
    n=size(A,1);n==0 && return zeros(Rat,0,size(B,2));M=hcat(A,B)
    for i in 1:n
        k=findfirst(j->M[j,i]!=0,i:n);need(k!==nothing,"LAYOUT_SINGULAR","Dirichlet matrix");k=i+k-1
        M[[i,k],:]=M[[k,i],:];pivot=M[i,i];M[i,:]./=pivot
        for j in 1:n
            if j!=i && M[j,i]!=0;a=M[j,i];M[j,:].-=a.*M[i,:];end
        end
    end
    M[:,n+1:end]
end
function certify_layout(st::DiagramState,positions)
    check_state(st)
    need(Set(keys(positions))==Set(keys(st.patches)),"LAYOUT_FACES","exact face set required")
    for (f,p) in st.patches
        pts=positions[f];allv=Set(v for t in p["triangles"] for v in t)
        need(Set(keys(pts))==allv,"LAYOUT_VERTICES",f);need(length(Set(values(pts)))==length(pts),"COINCIDENT_VERTICES",f)
        need(all(length(pt)==2 && all(x isa Rat for x in pt) for pt in values(pts)),"LAYOUT_RATIONAL",f)
        for (i,t) in enumerate(p["triangles"]);need(orient((pts[v] for v in t)...)>0,"FLIPPED_TRIANGLE","$f/$i");end
        es=sort(collect(keys(st.incidence[f])))
        for i in 1:length(es),j in i+1:length(es)
            a,b=es[i];c,d=es[j];common=intersect(Set([a,b]),Set([c,d]))
            if isempty(common)
                need(!segments_meet(pts[a],pts[b],pts[c],pts[d]),"MESH_EDGE_CROSSING","$f/$a,$b/$c,$d")
            else
                z=first(common);u=a==z ? b : a;v=c==z ? d : c
                need(!onsegment(pts[z],pts[u],pts[v]) && !onsegment(pts[z],pts[v],pts[u]),"OVERLAPPING_EDGES","$f/$z")
            end
        end
    end
    true
end
function realize(st::DiagramState)
    check_state(st)
    out=Dict{String,Dict{String,Point}}()
    for f in sort(collect(keys(st.patches)))
        bd=st.boundaries[f];pts=Dict(zip(bd,circle_boundary(length(bd))));adj=Dict()
        for (u,v) in keys(st.incidence[f]);push!(get!(adj,u,Set()),v);push!(get!(adj,v,Set()),u);end
        interior=sort(collect(setdiff(Set(keys(adj)),Set(bd))));index=Dict(v=>i for (i,v) in enumerate(interior));n=length(interior)
        A=zeros(Rat,n,n);B=zeros(Rat,n,2)
        for (v,i) in index
            A[i,i]=length(adj[v])
            for w in adj[v]
                if haskey(index,w);A[i,index[w]]-=1
                else;for j in 1:2;B[i,j]+=pts[w][j];end;end
            end
        end
        X=solve_exact(A,B);for (i,v) in enumerate(interior);pts[v]=(X[i,1],X[i,2]);end
        out[f]=pts
    end
    certify_layout(st,out);out
end

# Canonical UTF-8 JSON is for deterministic interchange, not class equivalence.
function json_string(s::AbstractString)
    io=IOBuffer();write(io,'"')
    for c in s
        if c=='"';write(io,"\\\"")
        elseif c=='\\';write(io,"\\\\")
        elseif c=='\b';write(io,"\\b")
        elseif c=='\f';write(io,"\\f")
        elseif c=='\n';write(io,"\\n")
        elseif c=='\r';write(io,"\\r")
        elseif c=='\t';write(io,"\\t")
        elseif Int(c)<0x20;write(io,"\\u"*lpad(string(Int(c),base=16),4,'0'))
        else;write(io,c);end
    end
    write(io,'"');String(take!(io))
end
function json(x)
    if x===nothing;"null"
    elseif x isa Bool;x ? "true" : "false"
    elseif x isa Integer;string(x)
    elseif x isa AbstractString;json_string(x)
    elseif x isa Rational;json_string(qstr(Rat(x)))
    elseif x isa AbstractDict;"{"*join([json_string(string(k))*":"*json(x[k]) for k in sort(collect(keys(x)),by=string)],",")*"}"
    elseif x isa AbstractVector || x isa Tuple;"["*join(json.(collect(x)),",")*"]"
    else;throw(DiagramError("JSON_TYPE",string(typeof(x))));end
end
write_json(path,x)=write(path,json(x)*"\n")
representation_digest(doc)=bytes2hex(sha256(json(Dict(k=>v for (k,v) in doc if k!="metadata"))))
function marking_digest(doc)
    o=Dict("surface"=>doc["surface"],"faces"=>[Dict("id"=>p["face"],"side_order"=>p["side_order"]) for p in doc["patches"]])
    bytes2hex(sha256(json(o)))
end
function endpoint_anchors(st,cid)
    check_state(st)
    c,ds=st.curves[cid];c["closed"] && return [];out=[]
    for (f,v) in [(ds[1].face,ds[1].u),(ds[end].face,ds[end].v)]
        p=st.patches[f];sid=only([s for s in p["side_order"] if v in p["sides"][s][2:end-1]])
        i=findfirst(==(v),p["sides"][sid]);push!(out,Dict("side"=>[f,sid],"parameter"=>qstr(side_parameters(p,sid)[i])))
    end
    out
end
function seam_word(st,cid)
    check_state(st)
    c,ds=st.curves[cid];pair=Dict();out=[]
    for ab in get(st.doc["surface"],"gluings",[]);a,b=Tuple(ab[1]),Tuple(ab[2]);pair[a]=b;pair[b]=a;end
    pairs=collect(zip(ds[1:end-1],ds[2:end]));if c["closed"];push!(pairs,(ds[end],ds[1]));end
    for (a,b) in pairs
        (a.face,a.v)==(b.face,b.u) && continue
        p=st.patches[a.face];sid=only([s for s in p["side_order"] if a.v in p["sides"][s][2:end-1]])
        s=(a.face,sid);t=pair[s];base=min(s,t);sign=s==base ? 1 : -1;vs=st.patches[base[1]]["sides"][base[2]]
        k=findfirst(==(a.v),p["sides"][sid])-1;if s!=base;k=length(vs)-1-k;end
        push!(out,Dict("side"=>collect(base),"direction"=>sign,"gate_rank"=>k))
    end
    out
end
function intersection_cycles(st)
    check_state(st)
    out=Dict()
    for (cid,(_,ds)) in st.curves
        out[cid]=[Dict("id"=>x["id"],"sign"=>x["sign"],"curves"=>x["curves"]) for d in ds for x in [get(st.crossing_at,(d.face,d.u),nothing)] if x!==nothing]
    end
    out
end
function legacy_projection(st)
    check_state(st)
    need(all(c["role"] in ("alpha","beta") && c["closed"] for (c,_) in values(st.curves)),"LEGACY_NOT_HEEGAARD","closed alpha/beta curves only")
    events=[]
    for x in st.doc["crossings"]
        a,b=x["curves"];sign=x["sign"]
        if st.curves[a][1]["role"]=="beta";a,b=b,a;sign=-sign;end
        need(st.curves[a][1]["role"]=="alpha" && st.curves[b][1]["role"]=="beta","LEGACY_NOT_HEEGAARD",x["id"])
        push!(events,Dict("id"=>x["id"],"alpha"=>a,"beta"=>b,"sign"=>sign))
    end
    cyc=intersection_cycles(st)
    Dict("genus"=>st.surface_report["genus"],
     "alpha"=>[Dict("id"=>c["id"],"cycle"=>[x["id"] for x in cyc[c["id"]]]) for c in st.doc["curves"] if c["role"]=="alpha"],
     "beta"=>[Dict("id"=>c["id"],"cycle"=>[x["id"] for x in cyc[c["id"]]]) for c in st.doc["curves"] if c["role"]=="beta"],
     "intersections"=>events,"source_representation_sha256"=>st.fingerprint,"marking_sha256"=>marking_digest(st.doc))
end
function fine_edge_occurrences(st)
    pairs=Dict();one=Dict()
    for (f,inc) in st.incidence,(e,occ) in inc
        ids=[(f,o[1],o[2],o[3]) for o in occ]
        if length(ids)==2;pairs[ids[1]]=ids[2];pairs[ids[2]]=ids[1]
        else;one[(f,e)]=ids[1];end
    end
    for (a,b) in st.seam_partner;pairs[one[a]]=one[b];end
    pairs
end
"""Normalized cut complement. Exported triangle indices are ZERO-based for interchange."""
function regions(st::DiagramState;barriers=nothing)
    check_state(st)
    blocked=barriers===nothing ? Set(keys(st.used)) : Set(barriers);epairs=fine_edge_occurrences(st);adj=Dict()
    cut(o)=(o[1],edge(o[3],o[4])) in blocked
    for (a,b) in epairs
        if !cut(a) && !cut(b);push!(get!(adj,a[1:2],Set()),b[1:2]);end
    end
    out=[]
    for comp in components(keys(st.triangles),adj)
        cn=[(f,i,j) for (f,i) in comp for j in 1:3];vd=UnionFind(cn)
        eo=[(f,i,u,v) for (f,i) in comp for (u,v) in tedges(st.triangles[(f,i)])];ed=UnionFind(eo)
        for a in eo
            b=get(epairs,a,nothing)
            if b!==nothing && b[1:2] in comp && !cut(a) && !cut(b)
                unite!(ed,a,b);ta,tb=st.triangles[a[1:2]],st.triangles[b[1:2]]
                unite!(vd,(a[1],a[2],findfirst(==(a[3]),ta)),(b[1],b[2],findfirst(==(b[4]),tb)))
                unite!(vd,(a[1],a[2],findfirst(==(a[4]),ta)),(b[1],b[2],findfirst(==(b[3]),tb)))
            end
        end
        bd=[a for a in eo if !haskey(epairs,a) || cut(a) || cut(epairs[a])];bnext=Dict()
        for a in bd
            t=st.triangles[a[1:2]];u=root!(vd,(a[1],a[2],findfirst(==(a[3]),t)));v=root!(vd,(a[1],a[2],findfirst(==(a[4]),t)))
            need(!haskey(bnext,u),"REGION_BOUNDARY_BRANCH","normalized cut complex");bnext[u]=(v,a)
        end
        cycles=[];unseen=Set(keys(bnext))
        while !isempty(unseen)
            seed=minimum(unseen);z=seed;cycle=[]
            while true
                need(z in unseen,"REGION_BOUNDARY_CYCLE","normalized cut complex")
                delete!(unseen,z);z,a=bnext[z];push!(cycle,a);z==seed && break
            end
            push!(cycles,cycle)
        end
        chi=length(Set(root!(vd,x) for x in cn))-length(Set(root!(ed,x) for x in eo))+length(comp);b=length(cycles);g=(2-b-chi)÷2
        need(2-2*g-b==chi && g>=0,"REGION_TOPOLOGY","normalized cut complex")
        simple=false
        if b==1
            ids=[root!(st.vertex_dsu,(a[1],a[3])) for a in cycles[1]];simple=length(ids)==length(Set(ids))
        end
        seed=minimum(comp)
        push!(out,Dict("seed"=>[seed[1],seed[2]-1],"triangles"=>[[f,i-1] for (f,i) in sort(collect(comp))],
            "chi"=>chi,"genus"=>g,"boundary_components"=>b,"normalized_disk"=>chi==1 && b==1,
            "embedded_closed_disk"=>chi==1 && b==1 && simple,
            "boundary_cycles"=>[[[a[1],a[2]-1,a[3],a[4]] for a in cyc] for cyc in cycles]))
    end
    out
end
function family_cut_report(st,role::AbstractString)
    ids=Set(cid for (cid,(c,_)) in st.curves if c["role"]==role)
    rs=regions(st;barriers=Set(e for (e,cs) in st.used if cs[1] in ids));g=st.surface_report["genus"];closed=st.surface_report["boundary_components"]==0
    Dict("role"=>role,"curves"=>length(ids),"components"=>length(rs),
      "cut_system"=>closed && length(ids)==g && length(rs)==1 && rs[1]["genus"]==0 && rs[1]["boundary_components"]==2*g,
      "topology"=>[Dict("genus"=>r["genus"],"boundary_components"=>r["boundary_components"]) for r in rs])
end
function check_joining_disk(st,seed;orientations=Dict())
    need(issubset(Set(keys(orientations)),Set(keys(st.curves))) && all(x isa Integer && !(x isa Bool) && x in (-1,1) for x in values(orientations)),"ORIENTATION_AUDIT","known curve -> +/-1")
    candidates=[r for r in regions(st) if seed in r["triangles"]];need(length(candidates)==1,"REGION_SEED",seed);r=only(candidates)
    result=Dict("seed"=>seed,"embedded_closed_disk"=>r["embedded_closed_disk"],"arcs_all_nonempty_connected"=>false,"boundary_outward_positive"=>false,"koda_joining"=>false)
    !r["embedded_closed_disk"] && return result
    labels=String[];against=Bool[]
    for (f,i,u,v) in r["boundary_cycles"][1]
        owners=get(st.used,(f,edge(u,v)),[]);length(owners)!=1 && return result;cid=only(owners);push!(labels,cid)
        d=only([d for d in st.curves[cid][2] if d.face==f && edge(d.u,d.v)==edge(u,v)])
        push!(against,((d.u,d.v)==(u,v) ? 1 : -1)*get(orientations,cid,1)==-1)
    end
    required=Set(cid for (cid,(c,_)) in st.curves if c["role"] in ("alpha","beta"));runs=String[]
    for x in labels;if isempty(runs) || x!=runs[end];push!(runs,x);end;end
    if length(runs)>1 && runs[1]==runs[end];pop!(runs);end
    result["arcs_all_nonempty_connected"]=Set(runs)==required && length(runs)==length(required)
    result["boundary_outward_positive"]=all(against);result["koda_joining"]=result["arcs_all_nonempty_connected"] && all(against)
    result
end
function move_vertex(st,pos,face,vertex,target)
    certify_layout(st,pos);need(haskey(st.patches,face),"MOVE_FACE",face)
    need(!(vertex in st.boundaries[face]),"MOVE_BOUNDARY","boundary motion requires a separate marked transport")
    need(haskey(pos[face],vertex),"MOVE_VERTEX",vertex);need(length(target)==2,"MOVE_TARGET","two exact coordinates")
    out=deepcopy(pos);out[face][vertex]=Tuple(qnum.(target));certify_layout(st,out)
    witness=Dict("kind"=>"interior_vertex_linear","face"=>face,"vertex"=>vertex,"old"=>qstr.(collect(pos[face][vertex])),
      "new"=>qstr.(collect(out[face][vertex])),"same_representation"=>representation_digest(st.doc))
    out,witness
end
roundq(q,d)=begin n=big(10)^d;fld(numerator(q)*n*2+denominator(q),2*denominator(q))//n;end
function decimal(q,d)
    z=roundq(q,d);n=abs(numerator(z))*big(10)^d÷denominator(z)
    (z<0 ? "-" : "")*string(n÷big(10)^d)*"."*lpad(string(n%big(10)^d),d,'0')
end
esc(s)=replace(string(s),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;",'\''=>"&#x27;")
function write_svg(st,pos,path;mesh=true)
    certify_layout(st,pos);digits=nothing;rpos=nothing
    for d in (6,12,24,48,96,192,384,768)
        candidate=Dict(f=>Dict(v=>Tuple(roundq(x,d) for x in p) for (v,p) in ps) for (f,ps) in pos)
        try
            certify_layout(st,candidate);digits=d;rpos=candidate;break
        catch e;e isa DiagramError || rethrow();end
    end
    need(digits!==nothing,"SVG_PRECISION_LIMIT","no certified decimal realization")
    fs=[p["face"] for p in st.doc["patches"]];n=length(fs);cols=min(n,2);rows=cld(n,cols);w=520*cols;h=500*rows;out=String[]
    push!(out,"<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $w $h\" role=\"img\">",
      "<style>.mesh{stroke:#bcc5d3;stroke-width:0.004;fill:none}.curve{fill:none;stroke-width:0.014;stroke-linejoin:round}.alpha{stroke:#b63837}.beta{stroke:#275db2}.aux{stroke:#7255a3}.label{font:0.065px monospace;fill:#152435}.gate{fill:white;stroke:#566479;stroke-width:0.006}</style>",
      "<metadata>representation_sha256=$(representation_digest(st.doc));decimal_digits=$digits;centerlines_certified=true</metadata>","<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
    for (ix,f) in enumerate(fs)
        ox=260+520*mod(ix-1,cols);oy=250+500*((ix-1)÷cols);p=st.patches[f];pts=rpos[f]
        point(v)=join([decimal(x,digits) for x in pts[v]],",")
        push!(out,"<text x=\"$(ox-225)\" y=\"$(oy-222)\" font-size=\"18\" font-family=\"sans-serif\">$(esc(f)) · marked disk</text>",
          "<g transform=\"translate($ox $oy) scale(195 -195)\">",
          "<polygon points=\"$(join(point.(st.boundaries[f])," "))\" fill=\"#f6f8fc\" stroke=\"#53647d\" stroke-width=\"0.01\"/>")
        for mark in get(st.doc,"region_marks",[])
            region=only([r for r in regions(st) if mark["seed"] in r["triangles"]])
            for (ff,i) in region["triangles"]
                if ff==f
                    poly=join(point.(p["triangles"][i+1])," ")
                    push!(out,"<polygon class=\"region\" data-region=\"$(esc(mark["id"]))\" points=\"$poly\" fill=\"#f1c94c\" fill-opacity=\"0.55\"><title>$(esc(mark["id"])): $(mark["claim"])</title></polygon>")
                end
            end
        end
        if mesh
            for (u,v) in sort(collect(keys(st.incidence[f])));push!(out,"<path class=\"mesh\" d=\"M $(point(u)) L $(point(v))\"/>");end
        end
        for c in st.doc["curves"]
            cid=c["id"];ds=st.curves[cid][2]
            for d in ds
                if d.face==f;push!(out,"<path class=\"curve $(c["role"])\" data-curve=\"$(esc(cid))\" d=\"M $(point(d.u)) L $(point(d.v))\"><title>$(esc(cid))</title></path>");end
            end
        end
        for v in st.boundaries[f]
            x,y=pts[v];push!(out,"<circle class=\"gate\" cx=\"$(decimal(x,digits))\" cy=\"$(decimal(y,digits))\" r=\"0.012\"/>")
        end
        labelv=union(Set(st.boundaries[f]),Set(v for (ff,v) in keys(st.crossing_at) if ff==f))
        for v in sort(collect(labelv))
            x,y=pts[v];lab=get(get(st.crossing_at,(f,v),Dict()),"id",v)
            push!(out,"<g transform=\"translate($(decimal(x,digits)) $(decimal(y,digits))) scale(1 -1)\"><text class=\"label\" x=\"0.025\" y=\"-0.025\">$(esc(lab))</text></g>")
        end
        for ((ff,v),cr) in st.crossing_at
            if ff==f
                x,y=pts[v];push!(out,"<circle cx=\"$(decimal(x,digits))\" cy=\"$(decimal(y,digits))\" r=\"0.018\" fill=\"#172333\"/>")
            end
        end
        push!(out,"</g>")
    end
    push!(out,"</svg>");write(path,join(out,"\n"))
    Dict("digits"=>digits,"certified_centerlines"=>true,"representation_sha256"=>representation_digest(st.doc))
end
function report(st,pos)
    cyc=intersection_cycles(st)
    Dict("surface"=>st.surface_report,"representation_sha256"=>representation_digest(st.doc),"marking_sha256"=>marking_digest(st.doc),
      "curves"=>Dict(cid=>Dict("seam_word"=>seam_word(st,cid),"endpoint_anchors"=>endpoint_anchors(st,cid),"crossings"=>cyc[cid]) for cid in keys(st.curves)),
      "regions"=>regions(st),"alpha_cut"=>family_cut_report(st,"alpha"),"beta_cut"=>family_cut_report(st,"beta"),
      "rational_positions"=>Dict(f=>Dict(v=>qstr.(collect(pt)) for (v,pt) in ps) for (f,ps) in pos))
end
end # module
