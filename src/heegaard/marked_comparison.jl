"""Evidence-separated comparisons of native marked diagrams.

Source: research root threads/PJ260922/ACTIVE.md, acceptance contract
(comparison); vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl::
representation_digest excludes only metadata. Compare
actual values below, rather than treating digest equality as a semantic test.
Operation claims reuse marked_moves.jl::replay_marked_move, not retained targets.
"""
_mc_math(doc)=Dict(k=>v for (k,v) in doc if k!="metadata")
_mc_model(doc)=Dict(k=>get(doc,k,nothing) for k in ("surface","patches"))
_mc_h1(h)=Dict("free_rank"=>h["free_rank"],"torsion"=>sort(collect(h["torsion"])))
_mc_claim(status,reason)=Dict{String,Any}("status"=>status,"reason"=>reason)

function _mc_summary(report)
    Dict("diagram_valid"=>report["diagram_valid"],
        "issues"=>get(report,"issues",Any[]),
        "source"=>get(report,"source",nothing),
        "manifold_applicable"=>report["manifold"]["applicable"],
        "manifold_reasons"=>report["manifold"]["reasons"])
end

# Source: marked_slides.jl::_mh_disk_policy!/marked_inverse_handle_slide.
# Slide preservation concerns the complete selected region in a common
# refinement. follow_seed validates a continuation, not necessarily that region.
# Inverse forward-region cells retain their original scope and direction.
function _mc_decoration_step(kind,spec,correspondence)
    if kind=="handle_slide_inverse"
        evidence=get(correspondence,"decoration_evidence",nothing)
        preserved=evidence===nothing ? (get(spec["forward_specification"],"disk_policy","preserve")=="preserve" ? true : nothing) :
            evidence["all_regions_preserved"]
        return Dict("direction"=>"inverse","status"=>preserved===true ? "preserved" : preserved===false ? "continued" : "unknown",
            "original_regions_preserved"=>preserved,"region_claims_revalidated"=>true,
            "evidence"=>deepcopy(evidence),"policy"=>get(spec["forward_specification"],"disk_policy","preserve"))
    elseif kind=="handle_slide"
        records=correspondence["marked_regions"];preserved=all(r["region_preserved"] for r in records)
        return Dict("direction"=>"forward","status"=>preserved ? "preserved" : "continued",
            "original_regions_preserved"=>preserved,"region_claims_revalidated"=>true,
            "cell_scope"=>"prepared common refinement before and after the slide",
            "regions"=>deepcopy(records),"policy"=>get(spec,"disk_policy","preserve"))
    end
    # Source: marked_stabilization.jl::marked_stabilize selected_regions and
    # marked_destabilize's actual generated cut/cap inverse. Accepted marked
    # claims retain their entire cell image; unsupported joining-D extension
    # is a constructor rejection, never a weakened decoration claim here.
    inverse=kind=="generated_destabilization"
    Dict("direction"=>inverse ? "inverse" : "forward","status"=>"preserved",
        "original_regions_preserved"=>true,"region_claims_revalidated"=>true,
        "cell_scope"=>inverse ? "inverse of verified outside-neck cell inclusion" : "entire source-region cell image outside the removed disk",
        "regions"=>deepcopy(get(correspondence,"selected_regions",Any[])))
end

function _mc_decorations(left,right,replay,valid,identity)
    result=Dict{String,Any}("scope"=>"supplied operation chain; endpoint data equality is a separate descriptive check",
        "endpoint_region_marks_data_equal"=>get(left,"region_marks",Any[])==get(right,"region_marks",Any[]),
        "steps"=>[Dict("index"=>s["index"],"kind"=>s["kind"],"evidence"=>s["region_decoration"])
            for s in replay["steps"] if haskey(s,"region_decoration")])
    if !valid || (isempty(get(left,"region_marks",Any[])) && isempty(get(right,"region_marks",Any[])))
        merge!(result,_mc_claim("inapplicable",!valid ? "invalid marked input" : "neither endpoint has selected regions"))
        result["original_regions_preserved"]=nothing
    elseif replay["status"]=="verified"
        values=[s["evidence"]["original_regions_preserved"] for s in result["steps"]]
        preserved=any(v->v===false,values) ? false : any(isnothing,values) ? nothing : true
        merge!(result,preserved===false ? _mc_claim("continued","the verified chain explicitly changes at least one selected region; its claim is revalidated") :
            preserved===true ? _mc_claim("preserved","the verified chain preserves the entire selected regions under its recorded transports") :
            _mc_claim("unknown","the verified underlying move has no complete original-region preservation summary"))
        result["original_regions_preserved"]=preserved
    else
        merge!(result,identity ? _mc_claim("preserved","endpoint identity independently preserves the selected regions") :
            _mc_claim("unknown","no region-preserving or explicitly continued operation chain is established"))
        result["original_regions_preserved"]=identity ? true : nothing
    end
    result
end

# Source: marked_moves.jl::_mm_execute and replay_marked_move. The whitelist
# grants only accepted constructors their established preservation claim.
# A simplicial action can be a nontrivial mapping class even on the same mesh.
# Subdivision identification must not be called fixed-model ambient isotopy.
function _mc_replay(left,right,witnesses,positions)
    isempty(witnesses) && return Dict{String,Any}("status"=>"absent","steps"=>Any[])
    current=left;realization=positions;steps=Any[]
    identifications=Set(["barycentric_subdivision","barycentric_collapse","dual_star","dual_star_collapse"])
    slides=Set(["handle_slide","handle_slide_inverse"])
    stabilizations=Set(["stabilization","generated_destabilization"])
    underlying_moves=union(slides,stabilizations)
    supported=union(Set(["triangle_isotopy","vertex_isotopy","simplicial_action"]),identifications,underlying_moves)
    try
        for (i,w) in enumerate(witnesses)
            w isa AbstractDict || throw(ArgumentError("move witness step $i must be an object"))
            kind=get(w,"kind","")
            kind in supported || return Dict{String,Any}("status"=>"unsupported",
                "reason"=>"comparison has no accepted preservation contract for move kind $kind",
                "failed_step"=>i,"steps"=>steps)
            for field in ("specification","inverse_specification")
                get(w,field,nothing) isa AbstractDict ||
                    throw(ArgumentError("move witness step $i $field must be an object"))
            end
            result=replay_marked_move(current,w;positions=realization)
            push!(steps,Dict("index"=>i,"kind"=>kind,
                "relation"=>kind in slides ? "underlying_heegaard_move" : kind in stabilizations ? "underlying_stabilization_move" :
                    kind in identifications ? "explicit_marked_identification" :
                    kind=="simplicial_action" ? "oriented_surface_map" : "fixed_model_ambient_isotopy",
                "correspondence"=>deepcopy(result.witness["correspondence"])))
            kind in underlying_moves && (steps[end]["region_decoration"]=_mc_decoration_step(kind,w["specification"],result.witness["correspondence"]))
            current=result.document
            hasproperty(result,:positions) && (realization=result.positions)
            # Source: marked_subdivision.jl::marked_subdivide/marked_dual_star
            # certify underlying surface identity through new cells. They do
            # not return a positional realization on the changed mesh.
            (kind in identifications || kind in underlying_moves) && (realization=nothing)
        end
        endpoint=current==right
        Dict{String,Any}("status"=>endpoint ? "verified" : "rejected",
            "reason"=>endpoint ? "constructive replay and inverse checks reached the supplied full target" :
                "constructed endpoint differs from the supplied full target",
            "endpoint_full_data_equal"=>endpoint,"steps"=>steps)
    catch err
        (err isa ArgumentError || err isa KeyError || err isa MarkedSurfaceDiagrams.DiagramError) || rethrow()
        Dict{String,Any}("status"=>"rejected","reason"=>sprint(showerror,err),
            "failed_step"=>length(steps)+1,"steps"=>steps)
    end
end

# Source: marked_groups.jl::marked_computation_report. H1 is a necessary
# invariant of homeomorphism, not a complete invariant or an orientation test.
# A contradiction to a reconstructed preservation witness remains visible.
function _mc_manifold(left,right;preserving_witness=false)
    applicable=left["applicable"] && right["applicable"]
    result=Dict{String,Any}("applicable"=>applicable,
        "group_isomorphism_decided"=>false,"manifold_recognition_performed"=>false,
        "orientation_distinguished_by_h1"=>false)
    if !applicable
        merge!(result,_mc_claim("inapplicable","closed represented-manifold evidence requires both applicability checks"))
        result["left_reasons"]=left["reasons"];result["right_reasons"]=right["reasons"]
        return result
    end
    a,b=_mc_h1(left["first_homology"]),_mc_h1(right["first_homology"])
    same=a==b
    result["left_h1"]=a;result["right_h1"]=b;result["h1_equal"]=same
    result["witness_invariant_conflict"]=preserving_witness && !same
    merge!(result,!same ? _mc_claim("refuted","different H1 groups obstruct homeomorphism; any preserving witness is contradictory") :
        preserving_witness ? _mc_claim("witnessed","the verified oriented-manifold-preserving construction reaches the target; H1 also agrees") :
        _mc_claim("consistent","H1 agrees; no homeomorphism or orientation decision follows"))
    result
end

"""Compare native data, marked isotopy witnesses, operations and manifold evidence.

`witnesses` is an ordered JSON-compatible move chain; each step is reconstructed
with inverse/applicability checks. `positions` is needed for realization motion.
Cellular subdivision/repatching supplies an explicit marked identification,
reported separately from ambient isotopy on one fixed triangulated model.
Handle slides and generated stabilization/inverse preserve the underlying
Heegaard manifold; they are not simultaneous surface maps. Selected-region
preservation or explicit continuation is reported separately for those moves.
Rejected or absent witnesses do not refute existence of another equivalence.
Metadata-only changes leave the mathematical representation equal, but explicit
infinity selection is additionally protected for an identity-isotopy claim.
No exhaustive isotopy search, group-isomorphism or manifold recognition is run.
"""
function marked_comparison_report(left::AbstractDict,right::AbstractDict;witnesses=Any[],positions=nothing)
    a,b=_plain_value(left),_plain_value(right)
    ar,br=marked_computation_report(a),marked_computation_report(b)
    valid=ar["diagram_valid"] && br["diagram_valid"]
    mathematical=_mc_math(a)==_mc_math(b)
    # Source: marked_topology.jl::_marked_states validates reserved display
    # metadata and supplied infinity cells. A malformed value is an invalid
    # input, never an absent infinity that could certify an identity map.
    infinity_equal=valid ? _mm_infinity(a)==_mm_infinity(b) : nothing
    identity=valid && mathematical && infinity_equal===true
    replay=valid ? _mc_replay(a,b,witnesses,positions) :
        Dict{String,Any}("status"=>"inapplicable","reason"=>"invalid marked input","steps"=>Any[])
    witnessed=replay["status"]=="verified"
    isotopies=witnessed && all(s["relation"]=="fixed_model_ambient_isotopy" for s in replay["steps"])
    identified=witnessed && any(s["relation"]=="explicit_marked_identification" for s in replay["steps"])
    surface_action=witnessed && any(s["relation"]=="oriented_surface_map" for s in replay["steps"])
    slide=witnessed && any(s["relation"]=="underlying_heegaard_move" for s in replay["steps"])
    stabilization=witnessed && any(s["relation"]=="underlying_stabilization_move" for s in replay["steps"])
    underlying=slide || stabilization
    fixed_model=_mc_model(a)==_mc_model(b)
    isotopy=!valid ? _mc_claim("inapplicable","both inputs must be valid marked embeddings") :
        identity ? _mc_claim("witnessed","identical mathematical data and explicit infinity selection give the identity isotopy") :
        isotopies && fixed_model ? _mc_claim("witnessed","all reconstructed steps are ambient isotopies on the fixed marked model") :
        _mc_claim("unknown","no fixed-model ambient isotopy witness has been established")
    isotopy["fixed_model_equal"]=fixed_model
    identification=!valid ? _mc_claim("inapplicable","both inputs must be valid marked embeddings") :
        identified && !surface_action && !underlying ? _mc_claim("witnessed","reconstructed cellular identity identifications and local isotopies identify the marked diagrams") :
        _mc_claim("unknown","no identity identification across the marked cell models has been established")
    moves=!valid ? _mc_claim("inapplicable","both inputs must be valid marked embeddings") :
        witnessed ? _mc_claim("witnessed","the supplied constructive operation chain reaches the target") :
        identity ? _mc_claim("witnessed","the identity operation suffices") :
        _mc_claim("unknown","no accepted operation chain reaches the supplied target")
    moves["relation"]=witnessed ? (underlying ? (slide && stabilization ? "underlying_heegaard_moves" : stabilization ? "underlying_stabilization_move" : "underlying_heegaard_move") :
        surface_action ? (identified ? "oriented_surface_map_under_identification" : "oriented_surface_map") :
        identified ? "explicit_marked_identification" : "fixed_model_ambient_isotopy") : identity ? "identity" : "unknown"
    moves["scope"]=underlying ? "underlying Heegaard move equivalence; selected-region continuation is separate" : "whole marked diagram under the stated relation"
    preserving=valid && (identity || witnessed)
    Dict("schema"=>"marked-comparison-report/1","left"=>_mc_summary(ar),"right"=>_mc_summary(br),
        "identity"=>Dict("status"=>valid ? "applicable" : "inapplicable",
            "full_data_equal"=>a==b,"mathematical_representation_equal"=>mathematical,
            "mathematical_convention"=>"all native fields except metadata, matching the vendor representation convention; direct values compared",
            "explicit_infinity_equal"=>infinity_equal),
        "marked_isotopy"=>isotopy,"marked_identification"=>identification,
        "move_equivalence"=>moves,"witness_replay"=>replay,
        "region_decoration"=>_mc_decorations(a,b,replay,valid,identity),
        "represented_manifold"=>_mc_manifold(ar["manifold"],br["manifold"];preserving_witness=preserving))
end

"""Compute raw o-data evidence without using any retained Heegaard source.

Adapted pure relation construction from research-root Topology.jl commit
5a372b64f3a11519d9c7787b67eefcb4db6613df,
src/ograph_invariants.jl::ograph_pi1_presentation (PD slots and relation order).
The existing app's Topology.pd_code supplies the same Sage/MATLAB slot order.
Only closed ordinary nonzero-vertex normal o-data is supported here. Database,
canonical-state and cache machinery are unnecessary for this pure consumer.
"""
function odata_computation_report(graph::OData)
    result=Dict{String,Any}("schema"=>"odata-computation-report/1","applicable"=>false,
        "reasons"=>String[],"source_metadata_used"=>false,
        "algorithm"=>"koda-edge-presentation",
        "implementation_source"=>"Topology.jl@5a372b64f3a11519d9c7787b67eefcb4db6613df/src/ograph_invariants.jl::ograph_pi1_presentation")
    isempty(graph.orientation) && (push!(result["reasons"],"zero-vertex o-data is outside this adapter");return result)
    if !all(s->s in (-1,1),graph.orientation) || any(w->isempty(w) || any(iszero,w),graph.gauss)
        push!(result["reasons"],"requires ordinary signed vertices and closed nonempty Gauss components")
        return result
    end
    T=parentmodule(@__MODULE__);v=T.VirtualLink()
    try
        T.set_data!(v;rgauss=deepcopy(graph.gauss),orientation=copy(graph.orientation))
        closed,conditions=T.is_closed(v);result["closed_conditions"]=collect(conditions)
        if !closed
            push!(result["reasons"],"normal o-data fails one or more closedness conditions C1,C2,C3")
            return result
        end
        pd=T.pd_code(v);count=maximum(pd)
        Set(vec(pd))==Set(1:count) || throw(ArgumentError("PD edge labels are not contiguous from one"))
        words=vcat([[(Int(r[1]),1),(Int(r[3]),-1)] for r in eachrow(pd)],
            [[(Int(r[4]),1),(Int(r[1]),1),(Int(r[2]),-1)] for r in eachrow(pd)])
        p=GroupPresentation(["e$i" for i in 1:count],words,Int[],
            Dict{String,Any}("pd_rows"=>[collect(r) for r in eachrow(pd)],
                "pd_row_slots"=>["g","k","l","h"],"relation_order"=>"all g*l^-1, then all h*g*k^-1"))
        h=first_homology(p)
        result["applicable"]=true
        result["presentation"]=_marked_group_report(p)
        result["first_homology"]=Dict("free_rank"=>h.free_rank,"torsion"=>h.torsion,
            "exponent_matrix_shape"=>collect(size(h.exponent_matrix)),
            "exponent_matrix_rows"=>[collect(h.exponent_matrix[i,:]) for i in axes(h.exponent_matrix,1)])
        result["input"]=Dict("gauss"=>deepcopy(graph.gauss),"orientation"=>copy(graph.orientation))
        result
    catch err
        err isa ArgumentError || rethrow()
        push!(result["reasons"],sprint(showerror,err));result
    end
end

"""Compare marked-manifold H1 with independently extracted raw o-data relations.

Source: marked_groups.jl::marked_computation_report and the pure odata adapter
above. Relation construction is independent (surface cells versus PD edges);
both reductions share groups.jl::first_homology/smith_normal_form. A match is
consistency evidence only. Graph source hashes/metadata cannot certify isotopy,
recover a source marking, or prove equality of represented oriented manifolds.
"""
function marked_comparison_report(left::AbstractDict,right::OData)
    a=marked_computation_report(left);b=odata_computation_report(right)
    unknown=_mc_claim("unknown","raw o-data provides no source marked-surface or operation witness")
    Dict("schema"=>"marked-odata-comparison-report/1","left"=>_mc_summary(a),"right"=>b,
        "identity"=>Dict("status"=>"inapplicable","reason"=>"different representation types"),
        "marked_isotopy"=>deepcopy(unknown),"marked_identification"=>deepcopy(unknown),"move_equivalence"=>deepcopy(unknown),
        "region_decoration"=>deepcopy(unknown),
        "represented_manifold"=>_mc_manifold(a["manifold"],b),
        "independence"=>Dict("relation_construction"=>"marked cellular attaching curves versus original o-data PD edge relations",
            "shared_algebra_backend"=>"groups.jl::first_homology/smith_normal_form",
            "retained_source_used"=>false,"reverse_conversion_used"=>false))
end
