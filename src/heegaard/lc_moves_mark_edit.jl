# cf: threads/PJ260922/ACTIVE.md Owner decisions 2026-11-12 — separate mark-moving operations.
# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-4-1-full.md — actual cell references and additive weight push-forward; reconstructed.
# Every catalogue row compiles to the same dictionary patch format and executor.

function _move_mark_target(input,field,target)
    cells=_move_mark_cells(_lc_of(input))
    if field=="points"
        target isa AbstractString && target in vcat(cells.faces,cells.edges,cells.empties) ||
            throw(LCMoveError(:invalid_mark_target,"point target is not a current cell",Dict("target"=>target)))
        return String(target)
    elseif field=="weights"
        target isa AbstractString && target in vcat(cells.faces,cells.weight_sides) ||
            throw(LCMoveError(:invalid_mark_target,"weight target is not a current region",Dict("target"=>target)))
        return String(target)
    elseif field=="disks"
        if target isa AbstractDict && haskey(target,"joining_region")
            input isa HeegaardCode || throw(LCMoveError(:layer,"joining-region selection requires HC",Dict()))
            id=Int(target["joining_region"]);rs=[r for r in forget_reference(realize_hc(input)).regions if r.id==id]
            length(rs)==1 && only(rs).joining_disk || throw(LCMoveError(:invalid_joining_region,"selected region is not a joining disk",Dict("region"=>id)))
            return copy(only(rs).face_roots)
        end
        refs=target isa AbstractString ? [String(target)] : target isa AbstractVector ? String.(target) : String[]
        length(refs)==1 && only(refs) in cells.faces && !(only(refs) in values(_lc_of(input).places)) ||
            throw(LCMoveError(:invalid_mark_target,"use one unpunctured face or an explicit joining_region",Dict("target"=>target)))
        return refs
    end
    throw(LCMoveError(:catalog,"unknown annotation field",Dict("field"=>field)))
end

function _move_mark_edit(input,row,selections)
    _require_selection(row,selections)
    field=String(row["expansion"]["field"]);name=String(selections["name"])
    (isempty(name) || occursin(r"[\s:#/]",name)) && throw(LCMoveError(:mark_name,"mark name is not HC-safe",Dict()))
    m=_marks(input);dict=getproperty(m,Symbol(field));target=_move_mark_target(input,field,selections["target"])
    edits=Any[]
    function edit(key,present,value)
        push!(edits,(field=field,key=key,before_present=haskey(dict,key),before=deepcopy(get(dict,key,nothing)),after_present=present,after=deepcopy(value)))
    end
    if field=="weights"
        haskey(dict,name) || throw(LCMoveError(:missing_mark,"source weight does not exist",Dict("name"=>name)))
        if name!=target
            edit(name,false,nothing);edit(target,true,Base.Checked.checked_add(get(dict,target,0),dict[name]))
        end
    else
        edit(name,true,target)
    end
    _witness_from(input,row,_component_names(_lc_of(input)),Tuple{Int,Int}[],Int[],Vector{Vector{Int}}(),Vector{Vector{Int}}(),Int[];
        extra=Dict{String,Any}("mark_patch"=>edits,"mark_selection"=>copy(selections),"support"=>Dict("kind"=>"annotation_patch")))
end

function _move_apply_mark_patch(input,w)
    row=heegaard_move_row(w.row_id)
    get(row["expansion"],"op",nothing)=="mark_edit" || throw(LCMoveError(:witness_plan_mismatch,"annotation patch on another row",Dict()))
    field=String(row["expansion"]["field"]);old=_marks(input)
    maps=Dict{String,Any}("disks"=>deepcopy(old.disks),"points"=>deepcopy(old.points),"weights"=>copy(old.weights))
    inverse_patch=Any[]
    for p in w.extra["mark_patch"]
        p.field==field || throw(LCMoveError(:witness_plan_mismatch,"patch field differs from catalog",Dict()))
        dict=maps[p.field]
        (get(w.extra,"free_preparation",false) || (haskey(dict,p.key)==p.before_present && (!p.before_present || dict[p.key]==p.before))) ||
            throw(LCMoveError(:mark_patch_mismatch,"selected annotation changed",Dict("field"=>field,"key"=>p.key)))
        if p.after_present
            if get(w.extra,"free_preparation",false)
                # A preparation inverse relocates current free annotations; it need not restore their old cells.
            elseif field=="points";_move_mark_target(input,field,p.after)
            elseif field=="disks"
                validfaces=Set(_move_mark_cells(_lc_of(input)).faces)
                !isempty(p.after) && all(x in validfaces for x in p.after) || throw(LCMoveError(:invalid_mark_target,"disk image has absent cells",Dict()))
            else
                _move_mark_target(input,field,p.key);p.after isa Int || throw(LCMoveError(:invalid_mark_target,"weight is not an integer",Dict()))
            end
            dict[p.key]=deepcopy(p.after)
        else;delete!(dict,p.key)
        end
        pushfirst!(inverse_patch,(field=p.field,key=p.key,before_present=p.after_present,before=p.after,after_present=p.before_present,after=p.before))
    end
    output=_move_rebuild_marks(input,maps["disks"],maps["points"],maps["weights"])
    labels=Set(abs.(reduce(vcat,values(_lc_of(input).words);init=Int[])))
    repair=nothing
    if get(w.extra,"free_preparation",false)
        output,repair=_move_transport_marks(output,output,w,Dict(k=>k for k in labels))
    end
    extra=copy(w.extra);extra["mark_patch"]=inverse_patch
    inverse=LCMoveWitness(isempty(w.inverse_id) ? w.row_id : w.inverse_id,w.row_id,w.digest,w.relation,w.layer,w.names,
        w.locations,w.vertices,w.after_frags,w.before_frags,w.capacity,_input_fingerprint(output),w.label_policy,0,extra)
    LCMoveResult(output,inverse,Dict{String,Any}("row"=>w.row_id,"relation"=>w.relation,"mark_patch"=>w.extra["mark_patch"],
        "label_map"=>Dict(k=>k for k in labels),"created_labels"=>Int[],"removed_labels"=>Int[],"mark_reassignment"=>repair))
end
