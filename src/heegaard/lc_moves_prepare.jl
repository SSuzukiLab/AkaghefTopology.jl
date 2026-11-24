# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-3-full.md §1 — separate a moving standard alpha from K before acting; reconstructed.
# cf: threads/PJ260922/ACTIVE.md Owner decisions 2026-11-12 — preparations are ordinary executable moves.

function _move_auto_separation(input,curve,side)
    stem=curve*"_"*String(input.roles[curve])*"_"*side;name=stem;i=1
    while haskey(input.lc.words,name);i+=1;name=stem*"_"*string(i);end
    only(_move_separate_role(input,heegaard_move_row("hc.role-support.separate"),Dict("curve"=>curve,"side"=>side,"target_name"=>name)))
end

function _move_prepare_slide(input,row,selections)
    curves=String.(selections["curves"]);a,b=curves
    band=selections["band"]
    band isa AbstractVector && isempty(band) || throw(LCMoveError(:shared_support_explicit_band_requires_transport,
        "an explicit band needs its ports carried onto the separated copy",Dict("curve"=>a)))
    out=LCMoveWitness[]
    for side in ("left","right")
        prep=_move_auto_separation(input,a,side);state=move(input,prep).output
        copyname=only(prep.extra["created_components"])
        next=copy(selections);next["curves"]=[copyname,b]
        for action in _expand_glue(state,row,next)
            push!(out,_move_sequence_witness(input,row,[prep,action];detail=Dict("kind"=>"separate_moving_role",
                "curves_before"=>curves,"curves_after"=>[copyname,b],"side"=>side,"role_component"=>Dict(a=>copyname))))
        end
    end
    out
end

function _move_prepare_affine(input,row,selections,acted)
    selections["support"] isa AbstractString || throw(LCMoveError(:annulus_reembedding_after_separation_required,
        "an explicit normal annulus needs the additional parallel curve in its source chart",Dict()))
    state=input;steps=LCMoveWitness[];names=copy(acted);renamed=Dict{String,String}()
    for (i,n) in enumerate(names)
        _is_k(state.lc,n) && haskey(state.roles,n) || continue
        prep=_move_auto_separation(state,n,"left");state=move(state,prep).output;push!(steps,prep)
        names[i]=only(prep.extra["created_components"]);renamed[n]=names[i]
    end
    next=copy(selections);next["acted_components"]=names
    out=LCMoveWitness[]
    for action in _expand_affine(state,row,next)
        push!(out,_move_sequence_witness(input,row,vcat(steps,[action]);detail=Dict("kind"=>"separate_moving_roles",
            "acted_before"=>acted,"acted_after"=>names,"role_component"=>renamed)))
    end
    out
end

# The attachment disk is a free selection. Move that annotation first, then use
# the existing certified terminal K chart. The requested and used cells are returned.
function _move_prepare_stabilization(input,row,selections)
    lc=input.lc;g=lc.genus
    g>=1 || throw(LCMoveError(:applicability,"terminal inclusion starts at genus one",Dict()))
    faces,right=_lc_faces(_lc_nonempty_words(lc.words));occ=_lc_rotation_occ(lc.words)
    c=(g,"-","out");canonical=right[c];component=Set([canonical]);queue=[canonical]
    for id in queue
        for d in only(f for f in faces if f.id==id).corners
            _is_k(lc,occ[d[1]][d[2]][1]) && continue
            nxt=right[_lc_tau(lc.words,occ,d)]
            nxt in component || (push!(component,nxt);push!(queue,nxt))
        end
    end
    request=String(selections["disk"]);named=haskey(input.disks,request)
    refs=named ? input.disks[request] : [request]
    valid=sort([x for x in refs if any(f.id==x for f in faces)];by=_move_ref_key)
    isempty(valid) && !named && throw(LCMoveError(:missing_selection,"attachment disk face does not exist",Dict("disk"=>request)))
    available=[x for x in valid if x in component]
    actual=isempty(available) ? canonical : first(available)
    detail=Dict{String,Any}("requested_disk"=>request,"source_refs"=>refs,"actual_face"=>actual,
        "reassigned"=>isempty(available),"policy"=>"free attachment disk followed by terminal K chart")
    if !named && actual==request
        action=_move_surface_inclusion(input,row,actual);action.extra["attachment"]=detail
        return [action]
    end
    name=named ? request : "stabilization_disk";i=1
    while !named && haskey(input.disks,name);i+=1;name="stabilization_disk_"*string(i);end
    prep=_move_mark_edit(input,heegaard_move_row("hc.mark.disk.move"),Dict("name"=>name,"target"=>actual))
    prep.extra["free_preparation"]=true
    state=move(input,prep).output;action=_move_surface_inclusion(state,row,actual)
    [_move_sequence_witness(input,row,[prep,action];detail=detail)]
end
