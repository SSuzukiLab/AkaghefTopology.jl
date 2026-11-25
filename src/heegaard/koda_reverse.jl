# Source: docs/specs/heegaard-koda-conversion.md §8–9 (reverse boundary).
"""Koda reverse construction via finite local tables (§8; §13.1–13.3)."""

"""
    spine_to_heegaard_koda(spine::DecodedSpine; curve_in_annulus=nothing, choices=Dict{String,Any}(), tables=load_koda_local_tables()) -> KodaConversionResult

Construct a punctured legacy Heegaard diagram from a decoded Koda spine.

# Arguments
`spine` must be a consistent result of `odata_to_spine` in
[the reverse-construction domain](../specs/heegaard-koda-conversion.md#8-入力-spine-そのものから-heegaard-diagramを作る).
Empty vertex data or a reconstruction without true intersections raises
`KodaNotComputed(:vertex_free_odata_not_encoded)`. Table, joining-region and
legacy certification failures propagate as Koda errors, `ArgumentError`
or `ValidationError`. The caller supplies compatible `tables` and cached
metadata. `curve_in_annulus` and `choices` are annotations, not overrides of
the geometric construction.

# Returns
A constructed `KodaConversionResult` with a `PuncturedHeegaardDiagram`,
merged choices, and correspondence recording genus, Euler characteristic,
intersection count, joining region and decoded sector/spine counts.

# Guarantees
1. The returned diagram passes legacy certification and its joining disk
   passes the current joining checks. [proved/agent: by construction]
2. The intended space is the thickening of the decoded branched spine.
   [unverified: general geometric identification and source alignment remain limited by Koda conversion §13.3]
3. The input and choices are not modified; nested choice values may remain
   shared. [proved/agent: by construction]

# Not guaranteed
This returns a legacy diagram without a Lickorish marking. The general
manifold identification in guarantee 2 is not certified by finite validation;
see [the source-alignment limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合).

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§8, 13.3 — reverse construction and limits of its identification.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> s = odata_to_spine(OData([[1,-1]], [1], Dict{String,Any}())); spine_to_heegaard_koda(s).value isa PuncturedHeegaardDiagram
true
```

See also [`mark_heegaard`](@ref).

Role: convert
"""
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

"""
    spine_to_heegaard_koda(spine::FlowSpine; kwargs...) -> KodaConversionResult

Reconstruct a punctured legacy Heegaard diagram from a flow spine's o-data.

# Arguments
`spine` must admit `to_ograph` and the checked reverse domain of
`odata_to_heegaard_koda`. Encoding, closedness, table and certification
errors propagate. Keywords are those of `odata_to_heegaard_koda`; unsupported
keywords raise `MethodError`. No regularization is requested by this method.

# Returns
The reverse conversion's constructed `KodaConversionResult`, with a punctured
diagram, choices and correspondence.

# Guarantees
1. On success the returned diagram and disk pass the reverse conversion's
   legacy checks. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The original diagram, labels and certificate are not reconstructed verbatim.
Identification with the spine's thickened 3-manifold remains within
[the reverse-construction evidence limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合).

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§7–8, 13.3 — o-data encoding and reverse geometric domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> ph = odata_to_heegaard_koda(OData([[1,-1]], [1], Dict{String,Any}())).value; s = heegaard_to_spine(ph).value;

julia> spine_to_heegaard_koda(s).value isa PuncturedHeegaardDiagram
true
```

See also [`spine_to_odata`](@ref).

Role: convert
"""
function spine_to_heegaard_koda(spine::FlowSpine; kwargs...)
    G = to_ograph(spine)
    odata_to_heegaard_koda(G; kwargs...)
end

"""
    odata_to_heegaard_koda(G::OData; curve_in_annulus=nothing, choices=Dict{String,Any}()) -> KodaConversionResult

Construct a punctured cellular Koda Heegaard diagram from closed normal o-data.

# Arguments
`G` must satisfy [the inverse Koda domain](../specs/heegaard-koda-conversion.md#8-入力-spine-そのものから-heegaard-diagramを作る).
`odata_to_spine` checks closedness C1–C3 and a nonempty orientation vector,
with `KodaInputError` on those failures. Table, legacy-certification and
joining-disk selection failures propagate as Koda errors, `ValidationError`
or `ArgumentError`. Unsupported vertex-free data raise `KodaNotComputed`.
`curve_in_annulus` is retained as a choice annotation in this implementation;
it does not override the local table's construction. `choices` supplies
additional annotations and is not a proof of their geometric validity.

# Returns
`KodaConversionResult` with `value::PuncturedHeegaardDiagram`,
`source_status == :constructed`, copied/merged choices, and correspondence
for genus, Euler characteristic, intersections, selected joining region and
table counts. Its diagram uses the legacy cellular representation.

# Guarantees
1. The returned legacy diagram passes full certification and its selected
   joining disk passes `validate_joining_disk`.
   [proved/agent: by construction]
2. It is intended to be the punctured Heegaard diagram of the thickened
   decoded spine. [unverified: general source alignment and the composed geometric identification remain limited by §13.3]
3. The o-data and supplied choices are not modified; nested choice values
   can remain shared. [proved/agent: by construction]

# Not guaranteed
This does not insert a fixed Lickorish reference system or return a resolved
LC-bound HC. Finite certification does not resolve the P1 source discrepancy
or independently prove manifold identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§8, 13.3 — inverse Koda construction and its source-alignment limits.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}());

julia> odata_to_heegaard_koda(g).value isa PuncturedHeegaardDiagram
true
```

See also [`odata_to_spine`](@ref), [`mark_heegaard`](@ref).

Role: convert
"""
function odata_to_heegaard_koda(G::OData; curve_in_annulus=nothing, choices=Dict{String, Any}())
  spine = odata_to_spine(G)
  spine_to_heegaard_koda(spine; curve_in_annulus=curve_in_annulus, choices=choices)
end

"""
    mark_heegaard(h::CellularHeegaardDiagram; h_map=nothing, choices=Dict{String,Any}()) -> Union{}

Report the unresolved reference-surface marking required to obtain an LC-bound HC.

# Arguments
`h` is a legacy certified diagram with computable surface genus.
`h_map` is the proposed identification of
[§9.1](../specs/heegaard-koda-conversion.md#91-基準曲面への同定).
Omission raises `KodaSpecificationPending(:surface_identification)`.
When supplied, its string- or symbol-keyed `orientation_preserving` field
must compare equal to `true`, and an optional `genus` must match the diagram;
violations raise `KodaInputError` with the mismatch. Passing these recorded
checks raises `KodaSpecificationPending(:lc_extraction)` because extraction
is not implemented. `choices` is accepted but currently unused.

# Returns
No normal value: the current function always throws after genus evaluation.

# Guarantees
1. No LC-bound HC is fabricated from a declaration of an identification.
   The exceptions identify the missing marking or extraction step of §9.
   [proved/agent: by construction]
2. The diagram, proposed map and choices are not modified.
   [proved/agent: by construction]

# Not guaranteed
Accepting the orientation flag is not verification of a PL map. No genus,
including genus zero, currently has a successful extraction path here.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§9, 13 — explicit marking and the unresolved extraction boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> try H.mark_heegaard(d); :returned catch e; e.topic end
:surface_identification
```

See also [`odata_to_hc`](@ref), [`hc_koda_to_standard`](@ref).

Role: convert
"""
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

"""
    odata_to_hc(G::OData; curve_in_annulus=nothing, h_map=nothing, choices=Dict{String,Any}()) -> Union{}

Attempt the inverse Koda construction and report its unresolved LC marking step.

# Arguments
`G`, `curve_in_annulus` and `choices` have the domain of
`odata_to_heegaard_koda`; its conversion and validation errors propagate.
`h_map` is passed to `mark_heegaard` under
[§9.1](../specs/heegaard-koda-conversion.md#91-基準曲面への同定).
Missing identification raises `KodaSpecificationPending(:surface_identification)`;
invalid orientation/genus declarations raise `KodaInputError`; otherwise
`KodaSpecificationPending(:lc_extraction)` reports the unimplemented step.

# Returns
No normal value in the current implementation. The constructed intermediate
punctured diagram can be obtained separately from `odata_to_heegaard_koda`.

# Guarantees
1. This entrance does not fabricate a resolved HC when the required
   reference marking or extraction is unavailable.
   [proved/agent: by construction]
2. Input o-data, map declarations and choices are not modified.
   [proved/agent: by construction]

# Not guaranteed
There is currently no successful LC-bound output path through this function,
even with an orientation-preserving declaration. Such a declaration is not
a verified PL identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§8–9, 13 — inverse construction and unresolved reference marking.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}());

julia> try odata_to_hc(g); :returned catch e; e.topic end
:surface_identification
```

See also [`odata_to_heegaard_koda`](@ref), [`mark_heegaard`](@ref).

Role: convert
"""
function odata_to_hc(G::OData; curve_in_annulus=nothing, h_map=nothing, choices=Dict{String, Any}())
  heegaard = odata_to_heegaard_koda(G; curve_in_annulus=curve_in_annulus, choices=choices)
  diagram = heegaard.value isa PuncturedHeegaardDiagram ? heegaard.value.diagram : heegaard.value
  mark_heegaard(diagram; h_map=h_map, choices=choices)
end
