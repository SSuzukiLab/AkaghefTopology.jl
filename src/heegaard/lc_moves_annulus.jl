# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-1-full.md §2 — common affine annulus and moved/fixed overlay; reconstructed
# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.2 — F_q and boundary port correspondence

const _MoveRational=Rational{BigInt}
_move_cross2(a,b)=a[1]*b[2]-a[2]*b[1]
_move_sub2(a,b)=(a[1]-b[1],a[2]-b[2])

"""A sufficiently thin regular collar excludes every crossing off its simple core.
Every transverse visit is a separate run, connected to its own original exterior ports.
The rational event positions avoid core/overlay triple points for integer shears.
"""
function _move_thin_annulus(lc,core)
    word=lc.words[core];n=length(word)
    n>0 || throw(LCMoveError(:annulus_support_required,"empty core needs an explicit collar",Dict()))
    length(unique(abs.(word)))==n || throw(LCMoveError(:support,"annulus core has a self crossing",Dict()))
    isempty(lc.places) || throw(LCMoveError(:missing_place_transport,"collar needs the position of place-only curves",Dict()))
    runs=Any[]
    for (i,x) in enumerate(word)
        name=_owner_of(lc.words,abs(x),x>0 ? "-" : "+")
        eps=-sign(x);s=BigInt(i)//BigInt(2n+1)
        t0=BigInt(eps==1 ? 0 : 1)//BigInt(1);t1=1-t0
        push!(runs,(component=name,core_label=abs(x),p0=(s,t0),p1=(s,t1)))
    end
    Dict{String,Any}("core"=>core,"runs"=>runs,"labels"=>Dict{Tuple{Int,Int,Int},Int}(),
        "range"=>"thin regular collar of a simple existing component")
end

function _move_strip_intersections(a,b)
    u=_move_sub2(a.p1,a.p0);v=_move_sub2(b.p1,b.p0);det=_move_cross2(u,v)
    hits=Tuple{Int,_MoveRational,_MoveRational,Int}[]
    low=floor(Int,min(a.p0[1],a.p1[1])-max(b.p0[1],b.p1[1]))-1
    high=ceil(Int,max(a.p0[1],a.p1[1])-min(b.p0[1],b.p1[1]))+1
    for deck in low:high
        offset=(b.p0[1]+deck-a.p0[1],b.p0[2]-a.p0[2])
        if det==0
            _move_cross2(offset,u)==0 && throw(LCMoveError(:affine_overlay_degenerate,"coincident strip arcs",Dict()))
            continue
        end
        ta=_move_cross2(offset,v)/det;tb=_move_cross2(offset,u)/det
        if 0<=ta<=1 && 0<=tb<=1
            (ta==0 || ta==1 || tb==0 || tb==1) && throw(LCMoveError(:affine_overlay_degenerate,"crossing on a boundary port",Dict()))
            push!(hits,(deck,ta,tb,-Int(sign(det))))
        end
    end
    hits
end

"""Extract all crossing events from one common periodic strip, with original port identities."""
function _move_annulus_extract(atlas;fresh=0,allocate=false)
    runs=atlas["runs"];labels=copy(atlas["labels"])
    events=[Tuple{_MoveRational,Int}[] for _ in runs]
    core_events=Tuple{_MoveRational,Int}[]
    for (i,r) in enumerate(runs)
        eps=Int(sign(r.p1[2]-r.p0[2]));tc=(BigInt(1)//BigInt(2)-r.p0[2])/(r.p1[2]-r.p0[2])
        s=r.p0[1]+tc*(r.p1[1]-r.p0[1]);s=mod(s,1)
        push!(events[i],(tc,eps*r.core_label));push!(core_events,(s,-eps*r.core_label))
    end
    used=Set{Tuple{Int,Int,Int}}()
    for i in eachindex(runs),j in 1:i-1
        for (deck,ti,tj,sgn) in _move_strip_intersections(runs[i],runs[j])
            key=(i,j,deck);push!(used,key)
            if !haskey(labels,key)
                allocate || throw(LCMoveError(:witness_mismatch,"annulus intersection lacks a source label",Dict()))
                fresh+=1;labels[key]=fresh
            end
            k=labels[key];push!(events[i],(ti,sgn*k));push!(events[j],(tj,-sgn*k))
        end
    end
    for ev in vcat(events,[core_events])
        sort!(ev;by=first)
        length(unique(first.(ev)))==length(ev) || throw(LCMoveError(:affine_overlay_degenerate,"annulus has a multiple crossing",Dict()))
    end
    atlas_out=copy(atlas);atlas_out["labels"]=Dict(k=>labels[k] for k in used)
    (core=[x[2] for x in core_events],runs=[[x[2] for x in ev] for ev in events],atlas=atlas_out)
end

function _move_annulus_map(atlas,acted,q)
    out=copy(atlas)
    out["runs"]=[r.component in acted ? merge(r,(p0=(r.p0[1]+q*r.p0[2],r.p0[2]),
        p1=(r.p1[1]+q*r.p1[2],r.p1[2]))) : r for r in atlas["runs"]]
    out
end

function _move_affine_witness(input,row,atlas,acted,q)
    lc=_lc_of(input);core=atlas["core"]
    all(haskey(lc.words,c) for c in acted) || throw(LCMoveError(:missing_selection,"acted component is absent",Dict()))
    any(_is_k(lc,c) && c!=core for c in acted) && throw(LCMoveError(:fixed_reference,"acted component is a reference curve",Dict()))
    _move_require_separated(input,acted)
    source=_move_annulus_extract(atlas)
    oldlabels=Set(abs.(reduce(vcat,values(lc.words);init=Int[])))
    target=_move_annulus_extract(_move_annulus_map(atlas,Set(acted),q);fresh=maximum(oldlabels),allocate=true)
    if input isa HeegaardCode
        for (i,j,_) in keys(target.atlas["labels"])
            a=atlas["runs"][i].component;b=atlas["runs"][j].component
            get(input.roles,a,nothing)!==nothing && get(input.roles,a,nothing)==get(input.roles,b,nothing) &&
                throw(LCMoveError(:same_family_crossing,"partial twist intersects curves of the same HC family",Dict("curves"=>[a,b])))
        end
    end
    names=_component_names(lc);index=Dict(n=>i for(i,n) in enumerate(names))
    _lc_canon_cyclic(source.core)==_lc_canon_cyclic(lc.words[core]) || throw(LCMoveError(:witness_mismatch,"core order disagrees with the annulus",Dict()))
    before=[copy(lc.words[core])];after=[target.core];locations=[(index[core],1)]
    for (i,r) in enumerate(atlas["runs"])
        hits=_LC_MV._cyclic_occurrences(lc.words[r.component],source.runs[i])
        length(hits)==1 || throw(LCMoveError(:witness_mismatch,"annulus run has no unique source interval",Dict("curve"=>r.component,"run"=>i)))
        push!(locations,(index[r.component],only(hits)));push!(before,source.runs[i]);push!(after,target.runs[i])
    end
    removed=setdiff(Set(abs.(reduce(vcat,before;init=Int[]))),Set(abs.(reduce(vcat,after;init=Int[]))))
    label_map=Dict(k=>k for k in setdiff(oldlabels,removed))
    _witness_from(input,row,names,locations,collect(oldlabels),before,after,collect(1:length(before));
        extra=Dict{String,Any}("label_map"=>label_map,"global_inverse"=>true,
            "annulus_after"=>target.atlas,"acted"=>acted,"q"=>q,
            "support"=>Dict("kind"=>"affine_annulus","core"=>core,"atlas"=>atlas)))
end
