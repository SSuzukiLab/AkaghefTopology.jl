# Source: docs/specs/heegaard-koda-conversion.md §3.1–3.4.
"""LC/HC realization from words via combinatorial rotation system only."""

function _lc_has_place_only_words(lc::LickorishCode)
    any(isempty(w) for w in values(lc.words)) && !isempty(lc.places)
end

"""
    realize_lc(lc::LickorishCode) -> KodaConversionResult

Return the combinatorial surface realization of a supported Lickorish code.

# Arguments
`lc` must satisfy the checked [LC input conditions](../specs/lickorish-code.md)
and the [realization domain](../specs/heegaard-koda-conversion.md#31-有限な回転系).
Validation and realization errors propagate. A nonempty placement dictionary
together with an empty component word raises `KodaSpecificationPending`
with condition `:sphere_place_embedding`; genus zero with all words empty
raises the same type with condition `:sphere_empty`.

# Returns
A constructed `KodaConversionResult` whose `value` has `"lc"`, `"genus"`
and `"realization"`, whose choices are empty, and whose correspondence
record the stage, surface Euler characteristic and face count.

# Guarantees
1. The realization is the finite surface arrangement of the supported input
   described in [Koda conversion §3](../specs/heegaard-koda-conversion.md#3-lc-の実現と-hc-の情報の付加).
   [proved/agent: by construction]
2. The input is not modified; the returned `"lc"` is that same object.
   [proved/agent: by construction]

# Not guaranteed
This result does not certify Heegaard cut systems, a joining disk or a
3-manifold. Constructed status does not resolve the excluded empty-word cases.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§3.1–3.4 — finite surface realization and its domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> r = realize_lc(lc); (r.source_status, r.value["genus"], r.correspondence["surface_chi"])
(:constructed, 1, 0)
```

See also [`realize_hc`](@ref).

Role: convert
"""
function realize_lc(lc::LickorishCode)
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

"""
    realize_hc(hc::HeegaardCode) -> KodaConversionResult

Return an LC surface realization with the Heegaard roles and marks attached.

# Arguments
`hc` must satisfy [HC input requirements](../specs/heegaard-code.md)
and the supported domain of `realize_lc`. `validate_hc_input` checks its
structural requirements; `KodaInputError`, `KodaSpecificationPending` and
realization errors propagate. Complete cut systems and geometric joining
conditions are not established by this operation.

# Returns
The constructed result has the LC realization fields plus `"hc"` and
`"roles_applied"` in `value`. Correspondence records `"hc_disks"`,
`"hc_points"`, `"hc_weights"` and `"roles_applied"`, together with the
LC realization stage, Euler characteristic and face count. Choices are empty.

# Guarantees
1. The result attaches the specified roles and marks to the realized surface
   arrangement of [Koda conversion §3.4](../specs/heegaard-koda-conversion.md#34-hc-の役割と追加行の実現).
   [proved/agent: by construction]
2. The input is not modified. The returned HC and LC references alias the
   input; annotation dictionaries are shallow copies.
   [proved/agent: by construction]

# Not guaranteed
No Heegaard certification or identification of a resulting 3-manifold is
returned. Nested disk lists can remain shared.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §3.4 — roles and marks on the LC realization.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc, Dict("m1"=>:alpha, "l1"=>:beta), Dict(), Dict(), Dict(), nothing);

julia> r = realize_hc(h); (r.source_status, r.correspondence["roles_applied"] == h.roles)
(:constructed, true)
```

See also [`forget_reference`](@ref).

Role: convert
"""
function realize_hc(hc::HeegaardCode)
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
