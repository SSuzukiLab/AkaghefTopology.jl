# Source: docs/specs/heegaard-koda-conversion.md §2.2, §5–7;
# conversion_spine.jl to_flow_spine / to_ograph (punctured Heegaard only).
"""Forward pipeline: LC/HC → forgotten Heegaard → punctured → spine → o-data."""

function heegaard_to_spine(ph::PuncturedHeegaardDiagram)
    _koda_combinatorial_spine_from_punctured(ph)
end

"""§6 from α/β crossing words and §5.2 joining boundary ports (no K in the diagram)."""
struct CombinatorialPuncturedInput
    genus::Int
    crossings::Vector{Dict{String, Any}}
    alpha_words::Dict{String, Vector{Int}}
    beta_words::Dict{String, Vector{Int}}
    joining_boundary_ports::Vector{RegionPort}
end

function heegaard_to_spine(inp::CombinatorialPuncturedInput)
  isempty(inp.joining_boundary_ports) &&
    throw(KodaInputError(:heegaard_to_spine, "joining_boundary_ports required", Dict()))
  α, β = _koda_word_vectors_from_dict(inp.alpha_words, inp.beta_words)
  dγ = _koda_d_gamma_from_ports(inp.joining_boundary_ports, inp.alpha_words, inp.beta_words)
  joining = _koda_joining_corners_from_ports(inp.joining_boundary_ports)
  _koda_combinatorial_spine_core(α, β, inp.crossings, joining, dγ, inp, inp.genus)
end

function heegaard_to_spine(forgot::ForgottenHeegaard; joining_region=nothing)
  sel = joining_region
  if sel === nothing
    joining = [r.id for r in forgot.regions if r.joining_disk]
    length(joining) == 1 && (sel = only(joining))
  end
  sel === nothing && throw(KodaInputError(:heegaard_to_spine, "joining region required", Dict()))
  region = _resolve_joining_region(forgot, sel)
  region.joining_disk || throw(KodaInputError(
    :heegaard_to_spine,
    "§5.2 joining disk conditions fail",
    Dict("region" => region.id, "spec" => "heegaard-koda-conversion.md §5.2"),
  ))
  result = if forgot.diagram !== nothing
    heegaard_to_spine(as_punctured_heegaard(forgot, sel).value)
  else
    _koda_combinatorial_spine_from_forgotten_region(forgot, region)
  end
  if !haskey(result.correspondence, "spine_path")
    result.correspondence["spine_path"] = forgot.diagram === nothing ? "ports" : "diagram"
  end
  result.value.complex.checks["spine_path"] = result.correspondence["spine_path"]
  if haskey(result.correspondence, "spine_switch_error")
    result.value.complex.checks["spine_switch_error"] = result.correspondence["spine_switch_error"]
  end
  _koda_attach_sectors!(result.value, forgot, region)
  result
end

function _selected_joining_disk(H::HeegaardDiagram, selection)
    selection isa JoiningDisk && return selection
    selection isa Integer || throw(KodaInputError(
        :hc_to_odata,
        "joining disk required; pass a JoiningDisk or a region index",
        Dict("joining_disk" => selection),
    ))
    joining_disk(H, selection)
end

"""Sector list that must not be encoded as o-data.

A complementary sector is a disk only when its own Euler relation holds and the
2-complex has Euler characteristic 1. Projecting a non-disk onto genus 0 and
still emitting a VirtualLink is not that value.
"""
function _koda_nondisk_odata(spine::FlowSpine)
    bad = _non_disk_sector(spine)
    topo = spine_topology(spine)
    sectors = topo["sectors"]
    euler = nothing
    inconsistent = false
    if !isempty(sectors)
        euler = length(topo["vertices"]) - length(topo["edges"]) + sum(s["chi"] for s in sectors)
        for sec in sectors
            g = get(sec, "genus", nothing)
            b = get(sec, "boundary_count", nothing)
            chi = get(sec, "chi", nothing)
            disk = g isa Integer && b isa Integer && chi isa Integer && chi == 2 - 2 * g - b && g == 0 && b == 1
            disk && continue
            inconsistent = true
            bad === nothing && (bad = sec)
        end
        euler == 1 || (inconsistent = true)
    end
    bad === nothing && !inconsistent && return nothing
    witness = bad === nothing ? first(sectors) : bad
    Dict{String, Any}(
        "spec" => "heegaard-koda-conversion.md §6.3",
        "kind" => "not_computed",
        "sector" => get(witness, "id", nothing),
        "genus" => get(witness, "genus", nothing),
        "boundary_count" => get(witness, "boundary_count", nothing),
        "euler" => euler,
    )
end

function _koda_virtual_image(spine::FlowSpine; stage::Symbol=:hc_to_odata)
    odata = try
        spine_to_odata(spine)
    catch err
        err isa KodaNotComputed && koda_not_computed(stage, err.reason, err.detail)
        rethrow()
    end
    odata_to_virtual_link(odata), odata
end

"""Koda map from a Heegaard diagram to the VirtualLink that holds its o-graph.

cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §2.2 — the value is that VirtualLink; its Gauss code is the smoothed o-graph
"""
function hc_to_odata(H::HeegaardDiagram; joining_disk=nothing, choices=Dict{String, Any}())
    disk = _selected_joining_disk(H, joining_disk)
    ph = make_punctured_combinatorial(H, disk)
    spine = heegaard_to_spine(ph)
    nondisk = _koda_nondisk_odata(spine.value)
    nondisk !== nothing && koda_not_computed(:hc_to_odata, :non_disk_sector, nondisk)
    link, odata = _koda_virtual_image(spine.value; stage=:hc_to_odata)
    corr = merge(copy(spine.correspondence), Dict{String, Any}(
        "odata" => odata,
        "image" => "VirtualLink",
        "source_object" => "HeegaardDiagram",
        "odata_status" => "constructed",
        "joining_disk_region" => disk.region,
    ))
    KodaConversionResult(
        link,
        merge(copy(choices), Dict{String, Any}("joining_disk_region" => disk.region)),
        corr,
        :constructed,
    )
end

function _hc_disk_selection(hc::HeegaardLayer; joining_disk=nothing, disk_name=nothing, disk_id=nothing)
  joining_disk !== nothing && return joining_disk
  disk_name !== nothing && return disk_name
  disk_id !== nothing && return disk_id
  if !isempty(hc.disks)
    length(hc.disks) == 1 && return first(keys(hc.disks))
    throw(KodaInputError(
      :hc_to_odata,
      "multiple HC disk lines; pass joining_disk explicitly",
      Dict("disks" => collect(keys(hc.disks))),
    ))
  end
  nothing
end

function hc_to_odata(hc::HeegaardLayer; joining_disk=nothing, disk_name=nothing, disk_id=nothing, choices=Dict{String, Any}())
  validate_hc_input(hc)
  realized = realize_hc(hc)
  forgot = forget_reference(realized)
  sel = _hc_disk_selection(hc; joining_disk=joining_disk, disk_name=disk_name, disk_id=disk_id)
  if sel === nothing
    cands = enumerate_joining_disk_candidates(forgot)
    length(cands) == 1 && (sel = cands[1]["region_id"])
  end
  if sel === nothing
    throw(KodaInputError(
      :hc_to_odata,
      "joining disk not selected; pass joining_disk or list a candidate seed",
      Dict("candidates" => enumerate_joining_disk_candidates(forgot)),
    ))
  end
  merge!(choices, Dict("joining_disk_selection" => sel))
  region = _resolve_joining_region(forgot, sel)
  if region !== nothing && !region.joining_disk
    reason = joining_disk_failure(region, merge(forgot.alpha_words, forgot.beta_words), hc.roles)
    throw(KodaInputError(
      :hc_to_odata,
      "§5.2 joining disk conditions fail",
      Dict(
        "region" => region.id,
        "condition" => reason === nothing ? "not_a_joining_disk" : string(reason),
        "spec" => "heegaard-koda-conversion.md §5.2",
      ),
    ))
  end
  punctured_corr = Dict{String, Any}()
  spine_path = "unset"
  diagram_path_error = nothing
  spine = if forgot.diagram !== nothing
    punctured = as_punctured_heegaard(forgot, sel; choices=choices)
    punctured_corr = punctured.correspondence
    spine_path = "ribbon"
    heegaard_to_spine(punctured.value)
  else
    region === nothing && throw(KodaInputError(:hc_to_odata, "joining region missing", Dict("selection" => sel)))
    spine_path = "ports"
    diagram_path_error = get(forgot.correspondence, "diagram_build_error", "diagram not built")
    _koda_combinatorial_spine_from_forgotten_region(forgot, region)
  end
  region !== nothing && _koda_attach_sectors!(spine.value, forgot, region)
  nondisk = _koda_nondisk_odata(spine.value)
  if nondisk !== nothing
    koda_not_computed(:hc_to_odata, :non_disk_sector, nondisk)
  end
  if get(spine.correspondence, "vertex_free", false)
    koda_not_computed(:hc_to_odata, :vertex_free_odata_not_encoded, Dict{String, Any}(
      "spec" => "heegaard-koda-conversion.md §7.0 / §13.3",
      "kind" => "not_computed",
      "region" => region === nothing ? nothing : region.id,
    ))
  end
  cγ_agreement = nothing
  if spine_path == "ribbon" && region !== nothing && haskey(spine.correspondence, "d_gamma")
    port_d = try
      _koda_d_gamma_from_ports(region.boundary_ports, forgot.alpha_words, forgot.beta_words)
    catch
      nothing
    end
    if port_d isa Dict
      cγ_agreement = port_d == spine.correspondence["d_gamma"]
    end
  end
  merge!(choices, Dict(
    "spine_path" => spine_path,
    "diagram_path_error" => diagram_path_error,
    "cγ_agreement" => cγ_agreement,
  ))
  link, odata = _koda_virtual_image(spine.value; stage=:hc_to_odata)
  corr = merge(
    forgot.correspondence,
    punctured_corr,
    spine.correspondence,
    Dict(
      "spine_edges" => length(spine.value.edges),
      "odata_status" => "constructed",
      "odata" => odata,
      "image" => "VirtualLink",
      "source_object" => forgot.diagram === nothing ? "ForgottenHeegaard" : "HeegaardDiagram",
    ),
  )
  KodaConversionResult(link, choices, corr, :constructed)
end

function hc_to_odata(lc_path::AbstractString, hc_path::AbstractString; kwargs...)
  hc_to_odata(read_hc(hc_path); kwargs...)
end
