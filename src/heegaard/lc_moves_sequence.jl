# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-3-full.md §§1,5 and threads/PJ260922/ACTIVE.md Owner decisions 2026-11-12 — compose certified preparation operations; reconstructed.
# A recipe contains executable witnesses, not saved input/output diagrams.

function _move_sequence_validate(w)
    steps=w.extra["sequence"]
    isempty(steps) && throw(LCMoveError(:invalid_recipe,"empty preparation recipe",Dict()))
    row=heegaard_move_row(get(w.extra,"origin_row_id",w.row_id))
    permitted=String.(get(row["expansion"],"preparation_rows",String[]))
    forward=get(w.extra,"support_direction",1)==1
    action=forward ? last(steps) : first(steps)
    target=forward ? row["id"] : something(row["inverse_id"],row["id"])
    action.row_id==target || throw(LCMoveError(:invalid_recipe,"recipe ends at another operation",Dict()))
    prep=forward ? steps[1:end-1] : steps[2:end]
    all(s->s.row_id in permitted && !haskey(s.extra,"sequence"),prep) ||
        throw(LCMoveError(:invalid_recipe,"preparation is not in this catalog row",Dict()))
    any(haskey(s.extra,"sequence") for s in steps) && throw(LCMoveError(:invalid_recipe,"nested recipe is not certified",Dict()))
    sum(s.genus_delta for s in steps)==w.genus_delta || throw(LCMoveError(:invalid_recipe,"genus change differs",Dict()))
    nothing
end

function _move_sequence_witness(input,row,steps;detail=Dict{String,Any}())
    extra=Dict{String,Any}("sequence"=>steps,"support"=>Dict("kind"=>"prepared_operation"),"preparation"=>detail)
    w=_witness_from(input,row,_component_names(_lc_of(input)),Tuple{Int,Int}[],Int[],Vector{Vector{Int}}(),Vector{Vector{Int}}(),Int[];extra=extra)
    _move_sequence_validate(w)
    w
end

function _move_apply_sequence(input,w)
    _input_fingerprint(input)==w.fingerprint || throw(LCMoveError(:witness_mismatch,"preparation recipe belongs to different geometry",Dict()))
    _move_sequence_validate(w)
    state=input;inverse_steps=LCMoveWitness[];receipts=Any[]
    source_labels=Set(abs.(reduce(vcat,values(_lc_of(input).words);init=Int[])))
    mapping=Dict(k=>k for k in source_labels)
    for s in w.extra["sequence"]
        r=move(state,s);state=r.output;pushfirst!(inverse_steps,r.inverse);push!(receipts,r.correspondence)
        labels=get(r.correspondence,"label_map",Dict(k=>k for k in values(mapping)))
        mapping=Dict(k=>labels[v] for(k,v) in mapping if haskey(labels,v))
    end
    inverse_extra=Dict{String,Any}("sequence"=>inverse_steps,"origin_row_id"=>w.extra["origin_row_id"],
        "support_direction"=>-get(w.extra,"support_direction",1),"support"=>w.extra["support"],"preparation"=>w.extra["preparation"])
    inverse=LCMoveWitness(isempty(w.inverse_id) ? w.row_id : w.inverse_id,w.row_id,w.digest,w.relation,w.layer,
        _component_names(_lc_of(state)),Tuple{Int,Int}[],Int[],Vector{Vector{Int}}(),Vector{Vector{Int}}(),Int[],
        _input_fingerprint(state),w.label_policy,-w.genus_delta,inverse_extra)
    target_labels=Set(abs.(reduce(vcat,values(_lc_of(state).words);init=Int[])))
    corr=Dict{String,Any}("relation"=>w.relation,"row"=>w.row_id,"support"=>w.extra["support"],
        "preparation"=>w.extra["preparation"],"steps"=>receipts,"label_map"=>mapping,
        "removed_labels"=>sort(collect(setdiff(source_labels,Set(keys(mapping))))),
        "created_labels"=>sort(collect(setdiff(target_labels,Set(values(mapping))))))
    LCMoveResult(state,inverse,corr)
end
