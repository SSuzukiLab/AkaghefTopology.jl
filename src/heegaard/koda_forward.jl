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
  if _non_disk_sector(spine.value) !== nothing
    bad = _non_disk_sector(spine.value)
    koda_not_computed(:hc_to_odata, :non_disk_sector, Dict{String, Any}(
      "spec" => "heegaard-koda-conversion.md §6.3",
      "kind" => "not_computed",
      "sector" => get(bad, "id", nothing),
      "genus" => get(bad, "genus", nothing),
      "boundary_count" => get(bad, "boundary_count", nothing),
    ))
  end
  odata = spine_to_odata(spine.value)
  corr = merge(
    forgot.correspondence,
    punctured_corr,
    spine.correspondence,
    Dict("spine_edges" => length(spine.value.edges), "odata_status" => "constructed"),
  )
  KodaConversionResult(odata, choices, corr, :constructed)
end

function hc_to_odata(lc_path::AbstractString, hc_path::AbstractString; kwargs...)
  hc_to_odata(read_hc(hc_path); kwargs...)
end
