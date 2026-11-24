# cf: docs/specs/heegaard-moves.md §6.4 — disk replacement, exact path extraction and standard K map.
# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-2-recovered.md §§3-5 — handle chart and boundary gluing; reconstructed.

# A simple dual cycle crosses the terminal l edge once. A tree to the chosen
# empty face supplies two ordered banks when a stem is necessary. The chosen K
# identification is the resulting ribbon-graph isomorphism, extended over disks.
function _move_terminal_cycle(lc,disk)
    active=_lc_nonempty_words(lc.words);faces,right=_lc_faces(active);occ=_lc_rotation_occ(active)
    tau(c)=_lc_tau(active,occ,c)
    any(f->f.id==disk,faces) || throw(LCMoveError(:missing_selection,"disk is not an input face",Dict("disk"=>disk)))
    parent=Dict{String,Any}(disk=>nothing);queue=[disk]
    for face in queue
        f=only(x for x in faces if x.id==face)
        for c in f.corners
            name=occ[c[1]][c[2]][1];_is_k(lc,name) && continue
            next=right[tau(c)];haskey(parent,next) && continue
            parent[next]=(face,tau(c));push!(queue,next)
        end
    end
    c=(lc.genus,"-","out");a=right[c];b=right[tau(c)]
    all(haskey(parent,f) for f in (a,b)) || throw(LCMoveError(:disk_transport_required,
        "ladder_end needs a disk in the terminal K face; another face needs a specified disk-push map",Dict("disk"=>disk)))
    function to_root(f)
        ds=Any[]
        while f!=disk
            f,d=parent[f];push!(ds,d)
        end
        ds
    end
    ds=vcat([c],to_root(b),[tau(d) for d in reverse(to_root(a))])
    boundary=Dict(f.id=>Dict(c=>i for(i,c) in enumerate(f.corners)) for f in faces)
    mid=[_move_route_point(lc,occ,right,Dict("dart"=>collect(d),"u"=>"1/2")) for d in ds]
    key(p)=(p.component,p.edge)
    for side in (-1,1)
        points=[merge(p,(u=count(x->key(x)==key(p),mid)>1 ? BigInt(3+side*_move_dir(p.dart))//BigInt(6) : BigInt(1)//BigInt(2),)) for p in mid]
        arcs=Dict{String,Vector{Tuple{_MoveRational,_MoveRational}}}()
        rank(d,u)=boundary[right[d]][d]-1+(_move_dir(d)==1 ? u : 1-u)
        seams=Int[]
        for i in eachindex(points)
            x=points[i];y=points[mod1(i+1,length(points))];d=tau(x.dart)
            right[d]==right[y.dart] || throw(LCMoveError(:disconnected_path,"terminal path has disconnected faces",Dict()))
            push!(get!(arcs,right[d],Tuple{_MoveRational,_MoveRational}[]),(rank(d,x.u),rank(y.dart,y.u)))
            right[d]==disk && push!(seams,i)
        end
        simple=all(!xor(min(a,b)<c<max(a,b),min(a,b)<d<max(a,b)) for cs in values(arcs) for (i,(a,b)) in enumerate(cs) for (c,d) in cs[1:i-1])
        simple || continue
        i=first(seams);points=vcat(points[i+1:end],points[1:i])
        face=only(f for f in faces if f.id==disk)
        ports=(rank(first(points).dart,first(points).u),rank(tau(last(points).dart),last(points).u))
        return (points=points,face=face,ports=ports,arcs=arcs,tree=parent)
    end
    throw(LCMoveError(:path_self_intersection,"terminal banks have alternating endpoints",Dict()))
end

# Exact oriented local surface. The old face boundary is parametrized clockwise;
# the two route ports have phases 0 and 1/2. The product strip has ds∧dt > 0.
# Its inner circle is glued to the handle by the explicit orientation-reversing table.
function _move_puncture_patch(template,route)
    Q=_move_rational;point(v)=(Q(v[1]),Q(v[2]))
    x=Q.(template["grid"]);lo,hi=Q.(template["hole"])
    tris=Any[]
    for i in 1:length(x)-1,j in 1:length(x)-1
        lo<=x[i]<x[i+1]<=hi && lo<=x[j]<x[j+1]<=hi && continue
        a=(x[i],x[j]);b=(x[i+1],x[j]);c=(x[i+1],x[j+1]);d=(x[i],x[j+1])
        push!(tris,(a,b,c),(a,c,d))
    end
    m=length(route.face.corners);p0,p1=route.ports;gap=mod(p1-p0,m)
    0<gap<m || throw(LCMoveError(:support,"puncture ports coincide",Dict()))
    phase(r)=begin z=mod(r-p0,m);z<=gap ? z/(2gap) : 1//2+(z-gap)/(2(m-gap)) end
    outer=[phase(BigInt(i)//BigInt(1)) for i in 0:m-1]
    knots=sort(unique(vcat(outer,[BigInt(i)//BigInt(8) for i in 0:8])))
    strip=[t for i in 1:length(knots)-1 for t in (((knots[i],0//1),(knots[i+1],0//1),(knots[i+1],1//1)),
        ((knots[i],0//1),(knots[i+1],1//1),(knots[i],1//1)))]
    boundary=point.(template["boundary"]);glue=Int.(template["glue_dst"])
    all(mod(glue[mod1(i+1,8)]-glue[i],8)==7 for i in 1:8) || throw(LCMoveError(:catalog,"handle boundary gluing reverses no orientation",Dict()))
    hpoint(s)=begin i=min(7,floor(Int,8s));u=8s-i;a=boundary[glue[i+1]];b=boundary[glue[mod1(i+2,8)]];((1-u)*a[1]+u*b[1],(1-u)*a[2]+u*b[2]) end
    # Refine each boundary triangle at every gluing knot; coordinates and cells
    # remain exact, so the two sides have an explicit common subdivision.
    insertions=[hpoint(s) for s in knots if !(hpoint(s) in boundary)]
    for p in unique(insertions)
        for ti in eachindex(tris)
            t=tris[ti];done=false
            for e in 1:3
                a=t[e];b=t[mod1(e+1,3)];c=t[mod1(e+2,3)]
                _move_cross2(_move_sub2(p,a),_move_sub2(b,a))==0 || continue
                u= b[1]!=a[1] ? (p[1]-a[1])/(b[1]-a[1]) : (p[2]-a[2])/(b[2]-a[2])
                0<u<1 || continue
                onhole=(a[1]==b[1] && a[1] in (lo,hi) && lo<=a[2]<=hi && lo<=b[2]<=hi) ||
                    (a[2]==b[2] && a[2] in (lo,hi) && lo<=a[1]<=hi && lo<=b[1]<=hi)
                onhole || continue
                tris[ti]=(a,p,c);push!(tris,(p,b,c));done=true;break
            end
            done && break
        end
    end
    Dict{String,Any}("kind"=>"punctured_face_and_torus","source_face"=>route.face.id,
        "external_boundary"=>route.face.corners,"external_phases"=>outer,"annulus_triangles"=>strip,
        "annulus_period"=>1,"handle_triangles"=>tris,"handle_periods"=>template["periodic_pairs"],
        "boundary_gluing"=>[(s,hpoint(s)) for s in knots],"route_ports"=>[0//1,1//2],
        "inverse"=>"remove the attached handle and cap the same boundary; retain every exterior port")
end

# Shared exact segment overlay: paths are ordered cell-local PL pieces. Passage
# signs come only from oriented determinants; semantic labels may pin known pairs.
function _move_planar_events(paths,pinned,fresh)
    names=sort(collect(keys(paths)));events=Dict(n=>Tuple{_MoveRational,Int}[] for n in names)
    vertices=Any[]
    for (i,a) in enumerate(names),b in names[1:i-1],(j,(a0,a1)) in enumerate(paths[a]),(k,(b0,b1)) in enumerate(paths[b])
        u=_move_sub2(a1,a0);v=_move_sub2(b1,b0);det=_move_cross2(u,v)
        det==0 && continue
        off=_move_sub2(b0,a0);t=_move_cross2(off,v)/det;s=_move_cross2(off,u)/det
        (0<=t<=1 && 0<=s<=1) || continue
        (0<t<1 && 0<s<1) || throw(LCMoveError(:affine_overlay_degenerate,"planar crossing lies on a port",Dict()))
        pair=Tuple(sort([a,b]));lab=haskey(pinned,pair) ? pinned[pair] : (fresh+=1)
        sig=-Int(sign(det));push!(events[a],(j+t,sig*lab));push!(events[b],(k+s,-sig*lab))
        push!(vertices,(label=lab,curves=(a,b),pieces=(j,k),parameters=(t,s),sign=sig))
    end
    for ev in values(events)
        sort!(ev;by=first)
        length(unique(first.(ev)))==length(ev) || throw(LCMoveError(:affine_overlay_degenerate,"planar overlay has a multiple crossing",Dict()))
    end
    (events=events,vertices=vertices,fresh=fresh)
end

function _move_surface_inclusion(input,row,disk)
    lc=input.lc;g=lc.genus
    g>=1 || throw(LCMoveError(:applicability,"terminal inclusion starts at genus one",Dict()))
    isempty(lc.places) || throw(LCMoveError(:missing_place_transport,"puncture and K identification need placed-circle transport",Dict()))
    route=_move_terminal_cycle(lc,disk);template=row["expansion"]["surface"]
    patch=_move_puncture_patch(template,route)
    forward=Dict(abs(x)=>_shift_label(abs(x),g) for w in values(lc.words) for x in w)
    fresh=max(_k_limit(g+1),maximum(values(forward)));ops=Any[];nword=Int[]
    for p in route.points
        label=p.component=="l$g" ? 3g : (fresh+=1)
        token=_move_dir(p.dart)*label;push!(nword,token);push!(ops,(p.component,p.edge,p.u,-token))
    end
    function fresh_role_name(prefix)
        i=g+1
        while haskey(lc.words,string(prefix,i));i+=1;end
        string(prefix,i)
    end
    rename=Dict("l"=>"l$(g+1)","m"=>"m$(g+1)","n"=>"n$g",
        "alpha"=>fresh_role_name("alpha"),"beta"=>fresh_role_name("beta"))
    created=[rename[k] for k in ("l","m","n","alpha","beta")]
    any(haskey(lc.words,n) for n in created) && throw(LCMoveError(:component_name,"new handle component name exists",Dict("names"=>created)))
    point(v)=(_move_rational(v[1]),_move_rational(v[2]))
    paths=Dict(rename[n]=>[(point(p[1]),point(p[2])) for p in ps] for(n,ps) in template["paths"])
    pinned=Dict(Tuple(sort([rename["l"],rename["m"]]))=>g+1,Tuple(sort([rename["l"],rename["n"]]))=>3g+1)
    overlay=_move_planar_events(paths,pinned,fresh)
    names=vcat(_component_names(lc),created);index=Dict(n=>i for(i,n) in enumerate(names))
    before=Vector{Int}[];after=Vector{Int}[];locations=Tuple{Int,Int}[]
    sort!(ops;by=x->(x[1],x[2],x[3]))
    _emit_gap_inserts!(locations,before,after,index,[(n,e+1,-x) for(n,e,u,x) in ops])
    for n in created
        word=last.(overlay.events[n]);n==rename["n"] && (word=vcat(nword,word))
        push!(before,Int[]);push!(after,word);push!(locations,(index[n],1))
    end
    # The geometric traces determine a marked K graph. Identify its vertices,
    # oriented edges and right-face boundary cycles with the canonical graph.
    # Both graphs fill the oriented genus-(g+1) surface, so these maps extend
    # over the complementary disks; no whole-input backup defines this map.
    projected,_,_,_,_,_,_=_kernel_apply(names,lc.words,before,after,locations,Int[],collect(1:length(before));label_map=forward)
    target=Dict(names[i]=>projected[i] for i in eachindex(names))
    ref=_lc_k_constants(g+1);limit=3g+1
    kg=Dict(n=>filter(x->abs(x)<=limit,target[n]) for n in keys(ref))
    all(_lc_canon_cyclic(kg[n])==_lc_canon_cyclic(ref[n]) for n in keys(ref)) ||
        throw(LCMoveError(:k_identification,"extracted K graph differs from the catalogue marking",Dict()))
    sourcefaces,_=_lc_faces(kg);targetfaces,_=_lc_faces(ref)
    face_map=Dict(f.id=>only(t.id for t in targetfaces if Set(t.corners)==Set(f.corners)) for f in sourcefaces)
    edge_map=Any[]
    for n in keys(ref)
        word=target[n];starts=findall(x->abs(x)<=limit,word)
        for (i,k) in enumerate(starts)
            stop=starts[mod1(i+1,length(starts))];len=mod(stop-k,length(word))
            len==0 && (len=length(word))
            push!(edge_map,(component=n,start=word[k],finish=word[stop],path=[word[mod1(k+j,length(word))] for j in 0:len]))
        end
    end
    kmap=Dict("vertices"=>collect(1:limit),"edges"=>edge_map,"faces"=>face_map,
        "extension"=>"orientation-preserving maps of complementary disks")
    _witness_from(input,row,names,locations,Int[],before,after,collect(1:length(before));
        extra=Dict{String,Any}("label_map"=>forward,"inverse_label_map"=>Dict(v=>k for(k,v) in forward),
            "created_components"=>created,"new_roles"=>Dict(rename["alpha"]=>:alpha,rename["beta"]=>:beta),
            "support"=>patch,"terminal_route"=>route,"handle_paths"=>paths,"handle_crossings"=>overlay.vertices,
            "k_identification"=>kmap))
end
