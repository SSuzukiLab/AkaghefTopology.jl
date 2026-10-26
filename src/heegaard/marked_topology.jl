"""Separate mathematical curves from the auxiliary reference before testing D.

Reuse the whole-cellulation validator, not the lossy intersection-only code.
See marked_display.jl::prepare_marked_display and the immutable vendor kernel's
regions/check_joining_disk. Neither this adapter nor validation chooses routes.
"""
function _marked_states(input::AbstractDict)
    M=MarkedSurfaceDiagrams
    doc=_plain_value(input)
    geometry=deepcopy(doc);delete!(geometry,"region_marks")
    st=M.validate(geometry)
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
is_closed_heegaard(doc::AbstractDict)=compression_boundary_report(doc).closed

"""Public mathematical report shared by consumers of a marked diagram.

Hashes identify the saved representation/marking, not an isotopy normal form.
An invalid embedding raises the kernel's structured DiagramError. A well-formed
closed-curve diagram whose compression has positive genus is reported as such.
"""
function marked_topology_report(doc::AbstractDict)
    M=MarkedSurfaceDiagrams;s=_marked_states(doc);st=s.mathematics
    Dict("schema"=>"marked-topology-report/1","surface"=>s.geometry.surface_report,
        "source_representation_sha256"=>M.representation_digest(s.document),
        "source_marking_sha256"=>M.marking_digest(s.document),
        "regions"=>M.regions(st),"intersection_cycles"=>M.intersection_cycles(st),
        "compression"=>_marked_compression(st),
        "joining_disks"=>[merge(Dict("id"=>m["id"],"claim"=>m["claim"]),M.check_joining_disk(st,m["seed"]))
            for m in get(s.document,"region_marks",[])])
end
