# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-3-full.md §1 — role-bearing reference meridian and its oriented parallel; reconstructed.

function _move_separate_role(input,row,selections)
    input isa HeegaardCode || throw(LCMoveError(:layer,"role separation requires HC",Dict()))
    lc=input.lc;curve=String(get(selections,"curve",""))
    _is_k(lc,curve) || throw(LCMoveError(:shared_support_not_reference,"select a role-bearing reference curve",Dict("curve"=>curve)))
    role=get(input.roles,curve,nothing)
    role in (:alpha,:beta) || throw(LCMoveError(:shared_support_without_role,"reference curve has no attaching role",Dict("curve"=>curve)))
    word=lc.words[curve]
    isempty(word) && throw(LCMoveError(:shared_support_empty_core,"a valid reference meridian has a K crossing",Dict()))
    length(unique(abs.(word)))==length(word) || throw(LCMoveError(:shared_support_core_not_simple,"core self-intersects",Dict()))
    isempty(lc.places) || throw(LCMoveError(:missing_place_transport,"parallel collar needs the placed-circle support",Dict()))
    sides=haskey(selections,"side") ? [String(selections["side"])] : ["left","right"]
    result=LCMoveWitness[]
    for side in sides
        side in ("left","right") || throw(LCMoveError(:support,"side must be left or right",Dict()))
        target=String(get(selections,"target_name",curve*"_"*String(role)*"_"*side))
        (isempty(target) || occursin(r"[\s:#/]",target) || target in string.(LC_META) || haskey(lc.words,target)) &&
            throw(LCMoveError(:component_name,"target must be an unused HC-safe component name",Dict("name"=>target)))
        run=(kind=:curve,component=curve,edge=length(word),side=side=="left" ? 1 : -1,direction=1,survives=false)
        w=_move_compile_boundary(input,row,target,[run],Dict("kind"=>"role_parallel","curve"=>curve,"side"=>side))
        w.extra["role_transfers"]=[(from=curve,to=target,role=role)]
        w.extra["curve"]=curve;w.extra["side"]=side
        push!(result,w)
    end
    result
end
