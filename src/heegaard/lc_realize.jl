# Source: docs/specs/heegaard-koda-conversion.md §3.1–3.4.
"""LC/HC realization from words via combinatorial rotation system only."""

function _lc_has_place_only_words(lc::LickorishConfiguration)
    any(isempty(w) for w in values(lc.words)) && !isempty(lc.places)
end

function realize_lc(lc::LickorishConfiguration)
    validate_lc_input(lc)
    g = lc.genus
    !isempty(lc.places) && _lc_has_place_only_words(lc) &&
        throw(KodaSpecificationPending(
            :realize_lc,
            :sphere_place_embedding,
            "empty-word placement on S^2 requires §3.2 / §14 side assignment",
            Dict("places" => copy(lc.places)),
        ))
    g == 0 && all(isempty(w) for w in values(lc.words)) &&
        throw(KodaSpecificationPending(:realize_lc, :sphere_empty, "genus 0 with no curve words", Dict()))
    realization = _lc_realize_combinatorial(lc)
    KodaConversionResult(
        Dict("lc" => lc, "genus" => g, "realization" => realization),
        Dict{String, Any}(),
        Dict("stage" => "realize_lc", "surface_chi" => realization.surface_chi, "faces" => length(realization.faces)),
        :constructed,
    )
end

function realize_hc(hc::HeegaardLayer)
    validate_hc_input(hc)
    base = realize_lc(hc.lc)
    merge!(
        base.value,
        Dict("hc" => hc, "roles_applied" => copy(hc.roles), "realization" => base.value["realization"]),
    )
    base.correspondence["hc_disks"] = copy(hc.disks)
    base.correspondence["hc_points"] = copy(hc.points)
    base.correspondence["hc_weights"] = copy(hc.weights)
    base.correspondence["roles_applied"] = copy(hc.roles)
    base
end
