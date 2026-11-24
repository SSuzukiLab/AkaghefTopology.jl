# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-1-full.md §0–§1 — right-side route, ribbon boundary and bank order; reconstructed
# cf: lib/topology.jl/docs/specs/heegaard-moves.md §3–§6 — finite support expanded through the shared reconnect kernel

_move_dir(c) = c[3] == "out" ? 1 : -1
function _move_rational(x)
    x isa Rational && return BigInt(numerator(x))//BigInt(denominator(x))
    x isa Integer && return BigInt(x)//BigInt(1)
    parts=split(replace(String(x),"//"=>"/"),'/')
    length(parts)==2 || throw(LCMoveError(:support,"edge position must be a rational fraction",Dict()))
    parse(BigInt,parts[1])//parse(BigInt,parts[2])
end

function _move_route_point(lc,occ,right,raw)
    raw isa AbstractDict || throw(LCMoveError(:support,"route point needs dart and u",Dict()))
    c=_move_corner(raw["dart"])
    haskey(right,c) || throw(LCMoveError(:support,"route dart is absent",Dict("dart"=>collect(c))))
    u=_move_rational(get(raw,"u","1/2"))
    0<u<1 || throw(LCMoveError(:support,"edge position must be interior",Dict()))
    name,pos0=occ[c[1]][c[2]]
    edge=_move_dir(c)==1 ? pos0+1 : mod1(pos0,length(lc.words[name]))
    (dart=c,component=name,edge=edge,u=u,face=right[c])
end

"""Certify an ordered route through disk faces; right-face darts and rational edge positions."""
function _move_certify_route(lc,band;avoid=String[])
    active=_lc_nonempty_words(lc.words); faces,right=_lc_faces(active);occ=_lc_rotation_occ(active)
    points=[_move_route_point(lc,occ,right,p) for p in vcat([band["start"]],get(band,"crossings",Any[]),[band["finish"]])]
    keys=[(p.component,p.edge,p.u) for p in points]
    length(unique(keys))==length(keys) || throw(LCMoveError(:path_order_required,"route events need distinct edge positions",Dict()))
    arcs=Dict{String,Vector{Tuple{Rational{BigInt},Rational{BigInt}}}}()
    boundary=Dict(f.id=>Dict(c=>i for (i,c) in enumerate(f.corners)) for f in faces)
    function rank(c,u)
        BigInt(boundary[right[c]][c]-1)//BigInt(1)+(_move_dir(c)==1 ? u : 1-u)
    end
    c=points[1].dart;u=points[1].u;visited=String[right[c]]
    for i in 2:length(points)
        p=points[i]
        right[c]==p.face || throw(LCMoveError(:disconnected_path,"route endpoint is not on the current face",Dict("step"=>i)))
        p.face in values(lc.places) && throw(LCMoveError(:path_order_required,"punctured face needs a normal-arc witness",Dict("face"=>p.face)))
        push!(get!(arcs,p.face,Tuple{Rational{BigInt},Rational{BigInt}}[]),(rank(c,u),rank(p.dart,p.u)))
        if i<length(points)
            p.component in avoid && throw(LCMoveError(:path_meets_family,"route crosses an excluded component",Dict("component"=>p.component)))
            c=_lc_tau(active,occ,p.dart);u=p.u;push!(visited,right[c])
        end
    end
    for (face,chords) in arcs, i in eachindex(chords), j in 1:i-1
        a,b=sort(collect(chords[i]));c,d=chords[j]
        xor(a<c<b,a<d<b) && throw(LCMoveError(:path_self_intersection,"face arcs have alternating endpoints",Dict("face"=>face)))
    end
    (start=first(points),finish=last(points),crossings=points[2:end-1],faces=visited,arcs=arcs)
end

"""Boundary orbits of an oriented ribbon graph; signed edge darts are paired by negation."""
function _move_ribbon_boundaries(rotations)
    rho=Dict{Int,Int}()
    for rotation in rotations, i in eachindex(rotation)
        d=Int(rotation[i]);haskey(rho,d) && throw(LCMoveError(:support,"repeated ribbon dart",Dict()))
        rho[d]=Int(rotation[mod1(i+1,length(rotation))])
    end
    all(haskey(rho,-d) for d in keys(rho)) || throw(LCMoveError(:support,"unpaired ribbon dart",Dict()))
    seen=Set{Int}();cycles=Vector{Int}[]
    for first in sort(collect(keys(rho)))
        first in seen && continue
        cycle=Int[];d=first
        while !(d in seen)
            push!(seen,d);push!(cycle,d);d=rho[-d]
        end
        d==first || throw(LCMoveError(:support,"ribbon boundary is not a cycle",Dict()))
        push!(cycles,cycle)
    end
    cycles
end

function _move_barbell_boundary(si,sj)
    t=load_heegaard_moves()["shared"]["barbell"]
    side(s)=s==1 ? "left" : "right"
    cycles=_move_ribbon_boundaries([t["start"][side(si)],t["finish"][side(sj)]])
    length(cycles)==3 || throw(LCMoveError(:catalog,"barbell boundary must have three components",Dict()))
    c=only(filter(c->3 in c && -3 in c,cycles))
    -1 in c && (c=reverse(-c))
    i=findfirst(==(1),c);c=vcat(c[i:end],c[1:i-1])
    c
end

function _move_face_bands(lc,a,b)
    faces,_=_lc_faces(_lc_nonempty_words(lc.words));occ=_lc_rotation_occ(lc.words)
    bands=Dict{String,Any}[]
    for face in faces
        face.id in values(lc.places) && continue
        aa=[c for c in face.corners if occ[c[1]][c[2]][1]==a]
        bb=[c for c in face.corners if occ[c[1]][c[2]][1]==b]
        for x in aa,y in bb
            push!(bands,Dict("start"=>Dict("dart"=>collect(x),"u"=>"1/2"),
                "finish"=>Dict("dart"=>collect(y),"u"=>"1/2"),"crossings"=>Any[]))
        end
    end
    bands
end

# Runs are geometry from a boundary orbit. The emitter knows parallel and bank pieces, not move kinds.
function _move_compile_boundary(input,row,target,runs,support)
    lc=_lc_of(input);names=_component_names(lc)
    created=!haskey(lc.words,target);created && push!(names,target)
    index=Dict(n=>i for(i,n) in enumerate(names))
    labels=Set(abs.(reduce(vcat,values(lc.words);init=Int[])))
    fresh=maximum(vcat([0],collect(labels)))
    out=Int[];ops=NamedTuple[]
    for run in runs
        if run.kind==:curve
            word=lc.words[run.component];n=length(word)
            inds=run.direction==1 ? [mod1(run.edge+k,n) for k in 1:n] : [mod1(run.edge-k,n) for k in 0:n-1]
            for pos in inds
                x=word[pos]
                if run.survives
                    push!(out,run.direction*x);continue
                end
                fresh+=1;token=run.direction*sign(x)*fresh;push!(out,token)
                partner=_owner_of(lc.words,abs(x),x>0 ? "-" : "+")
                pp=findfirst(==(-x),lc.words[partner]);before=run.side*sign(x)==1
                edge=before ? mod1(pp-1,length(lc.words[partner])) : pp
                push!(ops,(component=partner,edge=edge,u=BigInt(before ? 1 : 0)//BigInt(1),priority=0,token=-token))
            end
        elseif run.kind==:bank
            events=run.direction==1 ? run.route.crossings : reverse(run.route.crossings)
            for p in events
                fresh+=1;k=_move_dir(p.dart);token=run.direction*k*fresh;push!(out,token)
                push!(ops,(component=p.component,edge=p.edge,u=p.u,priority=-run.side*k,token=-token))
            end
        else
            throw(LCMoveError(:catalog,"unknown boundary run",Dict("kind"=>run.kind)))
        end
    end
    sort!(ops;by=p->(p.component,p.edge,p.u,p.priority))
    before=[copy(get(lc.words,target,Int[]))];after=[out];locations=[(index[target],1)]
    _emit_gap_inserts!(locations,before,after,index,[(p.component,p.edge+1,-p.token) for p in ops])
    _witness_from(input,row,names,locations,Int[],before,after,collect(1:length(before));
        extra=Dict{String,Any}("label_map"=>Dict(x=>x for x in labels),"global_inverse"=>true,
            "support"=>support,"runs"=>runs,"boundary_target"=>target,
            "created_components"=>(created ? [target] : String[])))
end

function _move_pants_expand(input,row,curves,band)
    lc=input.lc;a,b=String.(curves)
    (isempty(lc.words[a]) || isempty(lc.words[b])) && throw(LCMoveError(:path_order_required,"empty attaching curves need a collar witness",Dict()))
    !isempty(intersect(abs.(lc.words[a]),abs.(lc.words[b]))) && throw(LCMoveError(:applicability,"attaching curves intersect",Dict()))
    for n in (a,b)
        length(unique(abs.(lc.words[n])))==length(lc.words[n]) || throw(LCMoveError(:applicability,"attaching curve is not simple",Dict("curve"=>n)))
    end
    family=[n for(n,r) in input.roles if r==input.roles[a]]
    route=_move_certify_route(lc,band;avoid=family)
    (route.start.component,route.finish.component)==(a,b) || throw(LCMoveError(:support,"band endpoints do not match attaching curves",Dict()))
    si=-_move_dir(route.start.dart);sj=-_move_dir(route.finish.dart)
    cycle=_move_barbell_boundary(si,sj)
    runs=Any[]
    for d in cycle
        if abs(d)==1
            push!(runs,(kind=:curve,component=a,edge=route.start.edge,side=si,direction=sign(d),survives=true))
        elseif abs(d)==2
            push!(runs,(kind=:curve,component=b,edge=route.finish.edge,side=sj,direction=sign(d),survives=false))
        else
            push!(runs,(kind=:bank,route=route,side=si*sign(d),direction=sign(d)))
        end
    end
    _move_compile_boundary(input,row,a,runs,Dict("kind"=>"pants","route"=>route,
        "boundary_cycle"=>cycle,"coefficient"=>si*sj,"stationary"=>b))
end
