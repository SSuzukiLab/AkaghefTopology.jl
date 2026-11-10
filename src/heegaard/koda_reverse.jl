# Source: docs/specs/heegaard-koda-conversion.md §8–9 (reverse boundary).
"""Koda reverse construction via finite local tables (§8; §13.1–13.3)."""

function spine_to_heegaard_koda(spine::DecodedSpine; curve_in_annulus=nothing, choices=Dict{String, Any}(), tables=load_koda_local_tables())
    isempty(spine.vertex_signs) && koda_not_computed(
        :spine_to_heegaard_koda,
        :vertex_free_odata_not_encoded,
        Dict{String, Any}("spec" => "heegaard-koda-conversion.md §7.0 / §13.3"),
    )
    data = get(spine.metadata, "table_spine", nothing)
    data === nothing && (data = koda_table_decode(spine.gauss, spine.vertex_signs; tables=tables))
    hc = koda_hc_from_boundary(data, spine.vertex_signs; tables=tables)
  get(hc, "true_intersections", 0) == 0 && koda_not_computed(
      :spine_to_heegaard_koda,
      :vertex_free_odata_not_encoded,
      Dict{String, Any}("witness_intersections" => get(hc, "intersections", nothing)),
  )
    h = koda_heegaard_from_hc(hc)
    region = _koda_select_joining_region(h, hc)
    ph = make_punctured(h; joining_disk=joining_disk(h, region))
    corr = Dict(
        "genus" => hc["genus"],
        "chi" => hc["chi"],
        "intersections" => hc["intersections"],
        "joining_region" => region,
        "table_decode" => Dict(
            "sector_count" => data["sector_count"],
            "spine_chi" => data["spine_chi"],
        ),
    )
    KodaConversionResult(ph, merge(choices, Dict("curve_in_annulus" => curve_in_annulus)), corr, :constructed)
end

function spine_to_heegaard_koda(spine::FlowSpine; kwargs...)
    G = to_ograph(spine)
    odata_to_heegaard_koda(G; kwargs...)
end

function odata_to_heegaard_koda(G::OData; curve_in_annulus=nothing, choices=Dict{String, Any}())
  spine = odata_to_spine(G)
  spine_to_heegaard_koda(spine; curve_in_annulus=curve_in_annulus, choices=choices)
end

function mark_heegaard(h::HeegaardDiagram; h_map=nothing, choices=Dict{String, Any}())
  g = surface_genus(h)
  h_map === nothing &&
    throw(KodaSpecificationPending(
      :mark_heegaard,
      :surface_identification,
      "explicit orientation-preserving PL identification h: S → Σ_g is required (§9.1; §13). Automatic selection is unresolved",
      Dict("genus" => g),
    ))
  preserving = get(h_map, "orientation_preserving", get(h_map, :orientation_preserving, nothing))
  preserving == true || throw(KodaInputError(
    :mark_heegaard,
    "h must be orientation-preserving (§9.1)",
    Dict("orientation_preserving" => preserving),
  ))
  stated = get(h_map, "genus", get(h_map, :genus, nothing))
  stated === nothing || stated == g || throw(KodaInputError(
    :mark_heegaard,
    "h genus does not match the diagram",
    Dict("h_genus" => stated, "diagram_genus" => g),
  ))
  # g = 0 ⇒ K empty is determined. g ≥ 1 needs K = h⁻¹(K_std) and the
  # transversalization / word extraction of §9.2, which still lack a cell-map
  # algorithm (§13). Do not choose h.
  throw(KodaSpecificationPending(
    :mark_heegaard,
    :lc_extraction,
    "explicit h was accepted as orientation-preserving; LC word extraction after transporting K is not implemented (§9.2, §13)",
    Dict("genus" => g, "k_empty" => g == 0),
  ))
end

function odata_to_hc(G::OData; curve_in_annulus=nothing, h_map=nothing, choices=Dict{String, Any}())
  heegaard = odata_to_heegaard_koda(G; curve_in_annulus=curve_in_annulus, choices=choices)
  diagram = heegaard.value isa PuncturedHeegaardDiagram ? heegaard.value.diagram : heegaard.value
  mark_heegaard(diagram; h_map=h_map, choices=choices)
end
