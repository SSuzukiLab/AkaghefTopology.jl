"""Separate mathematical curves from the auxiliary reference before testing D.

Reuse the whole-cellulation validator, not the lossy intersection-only code.
See marked_display.jl::prepare_marked_display and the immutable vendor kernel's
regions/check_joining_disk. Neither this adapter nor validation chooses routes.
"""
function _marked_states(input::AbstractDict)
    M=MarkedSurfaceDiagrams
    doc=_plain_value(input)
    # Source: schema.jl::_marked_document (native metadata object contract),
    # marked_moves.jl::_mm_infinity and marked_display.jl::prepare_marked_display.
    # Other metadata remains opaque; these reserved fields are consumed as
    # typed geometric selections and must not silently become absent on error.
    infinity=nothing
    if haskey(doc,"metadata")
        metadata=doc["metadata"]
        M.need(metadata isa AbstractDict,"METADATA_OBJECT","metadata must be a JSON object")
        if haskey(metadata,"display")
            display=metadata["display"]
            M.need(display isa AbstractDict,"DISPLAY_OBJECT","metadata.display must be a JSON object")
            if haskey(display,"infinity_triangle")
                infinity=display["infinity_triangle"]
                M.need(infinity isa AbstractVector && length(infinity)==2 &&
                    infinity[1] isa AbstractString && infinity[2] isa Integer &&
                    !(infinity[2] isa Bool) && infinity[2]>=0,
                    "INFINITY_TRIANGLE","metadata.display.infinity_triangle must be [face, zero-based nonnegative triangle index]")
            end
        end
    end
    geometry=deepcopy(doc);delete!(geometry,"region_marks")
    st=M.validate(geometry)
    infinity!==nothing && M.need(haskey(st.triangles,(infinity[1],infinity[2]+1)),
        "INFINITY_TRIANGLE_MISSING","metadata.display.infinity_triangle does not select an existing source triangle")
    mathematical=deepcopy(doc)
    filter!(c->c["role"] in ("alpha","beta"),mathematical["curves"])
    ids=Set(c["id"] for c in mathematical["curves"])
    filter!(x->issubset(Set(x["curves"]),ids),mathematical["crossings"])
    (geometry=st,mathematics=M.validate(mathematical),document=doc)
end

"""Cut along one closed disjoint family, then cap only the new boundary circles.

The normalization is computed from oriented triangle-corner identifications by
MarkedSurfaceDiagrams.regions, with selected curve edges as barriers. A new cap
raises Euler characteristic by one; it leaves the component genus unchanged.
Original surface boundaries remain unfilled. Every output component retains its
source triangles and cut boundary cycles so the conclusion can be inspected.
"""
function _marked_family_compression(st,role::AbstractString)
    M=MarkedSurfaceDiagrams; M.check_state(st)
    ids=Set(cid for (cid,(c,_)) in st.curves if c["role"]==role)
    all(st.curves[c][1]["closed"] for c in ids) ||
        throw(ArgumentError("Compression along $role requires closed curves, not arcs"))
    barriers=Set(e for (e,owners) in st.used if !isdisjoint(Set(owners),ids))
    components=map(M.regions(st;barriers=barriers)) do region
        new_boundary=0
        for cycle in region["boundary_cycles"]
            flags=[(f,M.edge(u,v)) in barriers for (f,_,u,v) in cycle]
            all(flags) || all(!,flags) || error("Closed compression curve meets an original boundary")
            new_boundary+=all(flags)
        end
        Dict("source_triangles"=>region["triangles"],"cut_boundary_cycles"=>region["boundary_cycles"],
            "cut_genus"=>region["genus"],"cut_boundary_components"=>region["boundary_components"],
            "cut_euler_characteristic"=>region["chi"],"attached_caps"=>new_boundary,
            "genus"=>region["genus"],"boundary_components"=>region["boundary_components"]-new_boundary,
            "euler_characteristic"=>region["chi"]+new_boundary)
    end
    @assert sum(c["attached_caps"] for c in components)==2length(ids)
    @assert sum(c["euler_characteristic"] for c in components)==st.surface_report["chi"]+2length(ids)
    spherical=all(c["genus"]==0 && c["boundary_components"]==0 for c in components)
    complete=st.surface_report["boundary_components"]==0 && length(ids)==st.surface_report["genus"] &&
        length(components)==1 && spherical
    (;role=String(role),curves=sort(collect(ids)),components,
      spherical_compression=spherical,complete_cut_system=complete)
end

function _marked_compression(st)
    a=_marked_family_compression(st,"alpha");b=_marked_family_compression(st,"beta")
    (;surface_genus=st.surface_report["genus"],
      alpha_components=length(a.components),beta_components=length(b.components),
      alpha_boundary_genera=Int[c["genus"] for c in a.components],
      beta_boundary_genera=Int[c["genus"] for c in b.components],
      alpha_complete=a.complete_cut_system,beta_complete=b.complete_cut_system,
      closed=a.spherical_compression && b.spherical_compression,
      alpha=a,beta=b,method="normalized triangle-corner cut followed by capping only new boundaries")
end

"""Compute compression topology from a full marked document, without a count proxy."""
compression_boundary_report(doc::AbstractDict)=_marked_compression(_marked_states(doc).mathematics)
surface_genus(doc::AbstractDict)=_marked_states(doc).geometry.surface_report["genus"]
euler_characteristic(doc::AbstractDict)=_marked_states(doc).geometry.surface_report["chi"]
complement_regions(doc::AbstractDict)=MarkedSurfaceDiagrams.regions(_marked_states(doc).mathematics)
is_cellular(doc::AbstractDict)=all(r["normalized_disk"] for r in complement_regions(doc))
is_closed_heegaard(doc::AbstractDict)=validate(doc;level=:full).valid

# Source: validation.jl::validate/ValidationReport and the unchanged vendor
# MarkedSurfaceDiagrams._validate/check_joining_disk. The native kernel checks
# the whole embedding together; legacy intermediate levels therefore also run
# that check. :surface accepts proper arcs and incomplete compression systems;
# :full additionally asks whether both closed families compress to spheres.
function _marked_validation_issue(err,doc)
    M=MarkedSurfaceDiagrams
    err isa M.DiagramError || err isa ArgumentError || rethrow(err)
    code=err isa M.DiagramError ? Symbol(err.code) : :compression_not_applicable
    message=sprint(showerror,err)
    if code==:REGION_CLAIM_FAILED
        # Reuse the exact same disk audit for diagnostics, without accepting a
        # failed stored claim or changing its orientation/seed.
        bare=_plain_value(doc);marks=pop!(bare,"region_marks",Any[])
        st=_marked_states(bare).mathematics
        for mark in marks
            mark["id"]==err.detail || continue
            audit=M.check_joining_disk(st,mark["seed"])
            failed=[key for key in ("embedded_closed_disk","arcs_all_nonempty_connected","boundary_outward_positive")
                if !audit[key] && (mark["claim"]=="koda_joining" || key=="embedded_closed_disk")]
            message*="; failed conditions: "*join(failed,", ")
        end
    end
    ValidationIssue(:surface,code,message,"marked_document")
end

"""Validate a native marked embedding; :full also requires closed compression.

This uses the same ValidationReport as the ribbon API. A valid embedded diagram
may have arcs or positive-genus compression boundary: use level=:surface to
ask only for embedding and stored joining-disk validity.
"""
function validate(doc::AbstractDict;level=:full)
    requested=_validation_level(level)
    issues=ValidationIssue[];checks=Dict{Symbol,Bool}()
    states=try
        _marked_states(doc)
    catch err
        push!(issues,_marked_validation_issue(err,doc))
        checks[:marked_embedding]=false
        return ValidationReport(level,false,issues,checks)
    end
    checks[:marked_embedding]=true;checks[:joining_disks]=true
    if requested>=6
        try
            report=_marked_compression(states.mathematics)
            checks[:closed_heegaard]=report.closed
            report.closed || push!(issues,ValidationIssue(:closed_heegaard,:not_closed,
                "compression retains boundary: alpha genera $(report.alpha_boundary_genera), beta genera $(report.beta_boundary_genera); original physical boundaries are not capped","compression"))
        catch err
            push!(issues,_marked_validation_issue(err,doc))
            checks[:closed_heegaard]=false
        end
    end
    ValidationReport(level,isempty(issues),issues,checks)
end

# Source: cellulated_conversion.jl::cellulated_punctured, applicability gates.
# This is a preflight of that existing adapter, not a new conversion or a move
# certificate. Fine-region/ribbon bijection remains checked by the converter.
function _marked_conversion_applicability(st,compression,joining)
    conditions=Dict{String,Bool}(
        "closed_curves"=>all(c["closed"] for (c,_) in values(st.curves)),
        "complete_cut_systems"=>compression!==nothing && compression.alpha_complete && compression.beta_complete,
        "cellular_complement"=>all(r["normalized_disk"] for r in MarkedSurfaceDiagrams.regions(st)),
        "one_stored_disk"=>length(joining)==1,
        "stored_orientation_koda_disk"=>length(joining)==1 && only(joining)["koda_joining"])
    failed=[key for key in ("closed_curves","complete_cut_systems","cellular_complement","one_stored_disk","stored_orientation_koda_disk") if !conditions[key]]
    Dict("applicable"=>isempty(failed),"conditions"=>conditions,"reasons"=>failed,
        "scope"=>"cellular Koda adapter preconditions; conversion still checks fine-region/ribbon correspondence")
end

"""Public mathematical report shared by consumers of a marked diagram.

Hashes identify the saved representation/marking, not an isotopy normal form.
Invalid embeddings and false stored disk claims retain explanatory kernel issue
codes. Diagram validity, closed compression and conversion preconditions are
separate; a preflight does not assert that conversion was executed.
"""
function marked_topology_report(doc::AbstractDict)
    M=MarkedSurfaceDiagrams
    validity=validate(doc;level=:surface)
    issue_json(i)=Dict("level"=>string(i.level),"code"=>string(i.code),"message"=>i.message,"path"=>i.path)
    validity_json=Dict("valid"=>validity.valid,"issues"=>issue_json.(validity.issues))
    validity.valid || return Dict("schema"=>"marked-topology-report/1","diagram_valid"=>false,"validity"=>validity_json,
        "conversion"=>Dict("applicable"=>false,"reasons"=>["invalid_marked_diagram"]))
    s=_marked_states(doc);st=s.mathematics
    closed_curves=all(c["closed"] for (c,_) in values(st.curves))
    compression=closed_curves ? _marked_compression(st) : nothing
    joining=[merge(Dict("id"=>m["id"],"claim"=>m["claim"]),M.check_joining_disk(st,m["seed"]))
        for m in get(s.document,"region_marks",[])]
    Dict("schema"=>"marked-topology-report/1","diagram_valid"=>true,"validity"=>validity_json,
        "surface"=>s.geometry.surface_report,
        "source_representation_sha256"=>M.representation_digest(s.document),
        "source_marking_sha256"=>M.marking_digest(s.document),
        "regions"=>M.regions(st),"intersection_cycles"=>M.intersection_cycles(st),
        "compression"=>compression,"compression_applicable"=>closed_curves,
        "compression_reasons"=>(closed_curves ? String[] : ["compression requires closed alpha/beta curves, not proper arcs"]),
        "closed_heegaard"=>compression!==nothing && compression.closed,
        "conversion"=>_marked_conversion_applicability(st,compression,joining),
        "joining_disks"=>joining)
end
