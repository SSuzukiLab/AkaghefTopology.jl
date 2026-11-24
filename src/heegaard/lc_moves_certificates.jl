# cf: answer-3-full.md §5 — recompile the specified geometric support, not a site search; reconstructed.
# This protects the upper operations' executable fragments. Literal lower-level
# witnesses still have the legacy trusted-constructor contract.

function _move_recompile_upper(input,w)
    support=get(w.extra,"support",nothing);support isa AbstractDict || return nothing
    kind=get(support,"kind",nothing)
    row=heegaard_move_row(get(w.extra,"origin_row_id",w.row_id))
    if kind=="embedded_annulus"
        return _move_embedded_affine(input,row,support["certificate"],w.extra["acted"],w.extra["q"])
    elseif kind=="affine_annulus"
        return _move_affine_witness(input,row,support["atlas"],w.extra["acted"],w.extra["q"])
    elseif kind=="role_parallel"
        return only(_move_separate_role(input,row,Dict("curve"=>support["curve"],"side"=>support["side"],"target_name"=>w.extra["boundary_target"])))
    elseif kind=="pants"
        route=support["route"];point(p)=Dict("dart"=>collect(p.dart),"u"=>p.u)
        band=Dict("start"=>point(route.start),"finish"=>point(route.finish),"crossings"=>point.(route.crossings))
        return _move_pants_expand(input,row,[w.extra["boundary_target"],support["stationary"]],band)
    elseif kind=="punctured_face_and_torus"
        return _move_surface_inclusion(input,row,support["source_face"])
    end
    nothing
end

function _move_execution_equal(a,b)
    all(getfield(a,k)==getfield(b,k) for k in (:row_id,:inverse_id,:relation,:layer,:names,:locations,:vertices,:before_frags,:after_frags,:capacity,:fingerprint,:label_policy,:genus_delta)) &&
        all(get(a.extra,k,nothing)==get(b.extra,k,nothing) for k in ("identity","label_map","inverse_label_map","created_components","new_roles","role_transfers"))
end

function _move_same_named_cells(a,b)
    x,y=_lc_of(a),_lc_of(b)
    x.genus==y.genus && Set(keys(x.words))==Set(keys(y.words)) && x.places==y.places &&
        all(_lc_canon_cyclic(x.words[n])==_lc_canon_cyclic(y.words[n]) for n in keys(x.words)) &&
        _roles_of(a)==_roles_of(b)
end

function _move_checked_apply(input,w)
    get(w.extra,"destabilization_apply",false) && return _apply_recorded_destabilization(input,w)
    haskey(w.extra,"sequence") && return _move_apply_sequence(input,w)
    origin=heegaard_move_row(get(w.extra,"origin_row_id",w.row_id))
    (w.row_id==origin["id"] || w.row_id==origin["inverse_id"]) || throw(LCMoveError(:witness_plan_mismatch,"origin row is not paired with this row",Dict()))
    op=get(get(origin,"expansion",Dict()),"op","literal")
    flags=(get(w.extra,"annulus_recompute",false),get(w.extra,"embedded_annulus_recompute",false))
    any(flags) && op!="affine" && throw(LCMoveError(:witness_plan_mismatch,"affine inverse flag on another operation",Dict()))
    if op in ("affine","glue","k_inclusion","role_parallel")
        (w.row_id==origin["id"] || w.row_id==origin["inverse_id"]) && w.relation==origin["relation"] && w.layer==origin["layer"] ||
            throw(LCMoveError(:witness_plan_mismatch,"row metadata differs",Dict()))
        kind=get(get(w.extra,"support",Dict()),"kind",nothing)
        expected=Dict("affine"=>("affine_annulus","embedded_annulus"),"glue"=>("pants",),"k_inclusion"=>("punctured_face_and_torus",),"role_parallel"=>("role_parallel",))
        kind in expected[op] || (op=="affine" && get(w.extra,"identity",false)) ||
            throw(LCMoveError(:witness_plan_mismatch,"geometric certificate is missing",Dict()))
        (flags[1] && kind!="affine_annulus" || flags[2] && kind!="embedded_annulus") &&
            throw(LCMoveError(:witness_plan_mismatch,"inverse chart kind differs",Dict()))
    end
    # Affine inverses reconstruct F_-q on the certified current source chart.
    if get(w.extra,"annulus_recompute",false) || get(w.extra,"embedded_annulus_recompute",false)
        return _apply_witness(input,w)
    end
    if get(w.extra,"support_direction",1)==1
        rebuilt=_move_recompile_upper(input,w)
        if rebuilt!==nothing
            _move_execution_equal(w,rebuilt) || throw(LCMoveError(:witness_plan_mismatch,"fragments differ from the specified support",Dict("row"=>w.row_id)))
        end
        return _apply_witness(input,w)
    end
    result=_apply_witness(input,w)
    rebuilt=_move_recompile_upper(result.output,w)
    if rebuilt!==nothing
        haskey(w.extra,"mark_transport") && (rebuilt.extra["mark_transport"]=w.extra["mark_transport"])
        replay=_apply_witness(result.output,rebuilt)
        _move_same_named_cells(replay.output,input) || throw(LCMoveError(:witness_plan_mismatch,"inverse candidate does not reproduce the certified forward surgery",Dict("row"=>w.row_id)))
    end
    result
end
