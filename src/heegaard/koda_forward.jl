# Source: docs/specs/heegaard-koda-conversion.md §2.2, §5–7;
# conversion_spine.jl to_flow_spine / to_ograph (punctured Heegaard only).
"""Forward pipeline: LC/HC → forgotten Heegaard → punctured → spine → o-data."""

"""
    heegaard_to_spine(ph::PuncturedHeegaardDiagram) -> KodaConversionResult

Construct the Koda spine record of a punctured legacy Heegaard diagram.

# Arguments
`ph` must satisfy [the punctured-diagram domain, §5](../specs/heegaard-koda-conversion.md#5-heegaard-条件と-joining-disk-の入力検査).
Its joining disk is rechecked by `validate_joining_disk`; failures raise
`ArgumentError` naming the invalid disk, diagram or boundary incidence.
Inconsistent strand directions or crossing data can raise
`KodaInputError` or `KodaImplementationError` with the conversion stage.
The caller supplies the complete geometric hypotheses of §5; the record's
existence alone does not certify them.

# Returns
`KodaConversionResult` with `value::FlowSpine`, `source_status == :constructed`,
conversion choices and correspondence entries `true_vertices`, `d_gamma`
and `vertex_free`. The spine retains its source and an unverified conversion
certificate.

# Guarantees
1. The true-vertex correspondence consists of input intersections away
   from the selected joining boundary. [proved/agent: by construction]
2. The intended geometric object is the branched spine of the punctured
   diagram's 3-manifold specified in §6.3.
   [unverified: the general identification of the emitted finite record with the geometric Koda construction is not certified]
3. The punctured input is not modified; it can be retained by reference in
   the output. [proved/agent: by construction]

# Not guaranteed
Constructed status is not verified status. Disk sectors, a nonempty o-data
encoding and a general manifold-identification certificate are not guaranteed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5–6.4 — punctured domain, the Koda spine and source correspondence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); ph = odata_to_heegaard_koda(g).value;

julia> back = spine_to_odata(heegaard_to_spine(ph).value);

julia> compare_odata_words(g, back)["equal_up_to_metadata"]
true
```

See also [`spine_to_odata`](@ref), [`spine_topology`](@ref).

Role: convert
"""
function heegaard_to_spine(ph::PuncturedHeegaardDiagram)
    _koda_combinatorial_spine_from_punctured(ph)
end

"""
    CombinatorialPuncturedInput(genus, crossings, alpha_words, beta_words, joining_boundary_ports) -> CombinatorialPuncturedInput

Store punctured Heegaard crossing words and joining-boundary ports without a K marking.

# Arguments
`genus::Int`, `crossings::Vector{Dict{String,Any}}`, alpha and beta
`Dict{String,Vector{Int}}` word maps, and a vector of `RegionPort` records,
or convertible values. Conversion errors propagate. The caller supplies
consistent data in [Koda conversion §§5.2–6](../specs/heegaard-koda-conversion.md);
no embedded joining-disk or incidence validation occurs here.

# Returns
A `CombinatorialPuncturedInput` retaining the supplied genus, crossing
records, curve words and joining-boundary ports.

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers are
   shared and not modified. [proved/agent: by construction]
2. The result is an input record for the specified punctured-diagram spine
   construction; construction of this record performs no conversion.
   [proved/agent: by construction]

# Not guaranteed
The record need not describe a valid punctured surface or joining disk.
It carries no independent certificate of the associated 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5.2–6 — joining-boundary data for the forward construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> x = CombinatorialPuncturedInput(0, Dict{String,Any}[], Dict{String,Vector{Int}}(), Dict{String,Vector{Int}}(), H.RegionPort[]); (x.genus, isempty(x.joining_boundary_ports))
(0, true)
```

See also [`heegaard_to_spine`](@ref).

Role: construct
"""
struct CombinatorialPuncturedInput
    genus::Int
    crossings::Vector{Dict{String, Any}}
    alpha_words::Dict{String, Vector{Int}}
    beta_words::Dict{String, Vector{Int}}
    joining_boundary_ports::Vector{RegionPort}
end

"""
    heegaard_to_spine(inp::CombinatorialPuncturedInput) -> KodaConversionResult

Construct a Koda spine record from crossing words and a joining-boundary record.

# Arguments
`inp` supplies the crossing and boundary-port data of
[§§5.2–6](../specs/heegaard-koda-conversion.md).
The caller guarantees their geometric consistency and joining-disk conditions;
this overload does not validate an embedded joining disk on a ribbon carrier.
Empty `joining_boundary_ports` raises `KodaInputError` requesting that field.
Missing or inconsistent curve directions and crossings raise Koda input or
implementation errors; malformed dictionaries can also raise `KeyError`.

# Returns
`KodaConversionResult` with `value::FlowSpine`, `source_status == :constructed`,
conversion choices and the `true_vertices`, `d_gamma`, `vertex_free`
correspondence. The output retains its input as source data.

# Guarantees
1. The correspondence identifies the supplied crossings that remain away
   from the joining corners. [proved/agent: by construction]
2. Under the geometric domain of §§5–6, the intended result is the Koda
   branched spine of the punctured diagram's 3-manifold.
   [unverified: this overload does not certify a geometric realization of the supplied port record]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This low-level entrance does not certify a complete embedded Heegaard diagram
or a closed normal o-graph encoding. Missing geometric input cannot be inferred
from a successful word construction.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5.2–6.4 — joining boundary ports and the Koda spine construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> inp = CombinatorialPuncturedInput(1, Dict{String,Any}[], Dict{String,Vector{Int}}(), Dict{String,Vector{Int}}(), Topology.HeegaardCodes.RegionPort[]);

julia> try heegaard_to_spine(inp); false catch e; e isa KodaInputError end
true
```

See also [`as_punctured_heegaard`](@ref), [`spine_to_odata`](@ref).

Role: convert
"""
function heegaard_to_spine(inp::CombinatorialPuncturedInput)
  isempty(inp.joining_boundary_ports) &&
    throw(KodaInputError(:heegaard_to_spine, "joining_boundary_ports required", Dict()))
  α, β = _koda_word_vectors_from_dict(inp.alpha_words, inp.beta_words)
  dγ = _koda_d_gamma_from_ports(inp.joining_boundary_ports, inp.alpha_words, inp.beta_words)
  joining = _koda_joining_corners_from_ports(inp.joining_boundary_ports)
  _koda_combinatorial_spine_core(α, β, inp.crossings, joining, dγ, inp, inp.genus)
end

"""Use the original LC side owners to locate the same disk on the ribbon carrier."""
function _koda_original_punctured(forgot::ForgottenHeegaard, region::CombinatorialRegion; choices=Dict{String,Any}())
  h = forgot.diagram
  source = _koda_original_complement(forgot, h, region)
  ci = only(i for (i, cycle) in enumerate(source.cycles)
            if all(source.owner[d] == region.id for d in cycle))
  ri = only(i for (i, block) in enumerate(h.code.regions.boundary_cycle_blocks) if block == [ci])
  jd = joining_disk(h, ri)
  ph = make_punctured_combinatorial(h, jd)
  merge!(choices, Dict("joining_region" => region.id))
  KodaConversionResult(ph, choices, Dict{String,Any}(
      "joining_region" => region.id, "ribbon_joining_region" => ri,
      "joining_region_faces" => copy(region.face_roots)), :constructed)
end

"""
    heegaard_to_spine(forgot::ForgottenHeegaard; joining_region=nothing) -> KodaConversionResult

Construct the Koda spine of a forgotten HC with a selected joining region.

# Arguments
`forgot` is a result of `forget_reference`. `joining_region` selects a region
under [§5.2](../specs/heegaard-koda-conversion.md); when omitted it is selected
only if exactly one joining-disk candidate exists. Otherwise `KodaInputError`
requires a selection. A selected nonjoining region raises `KodaInputError`
identifying §5.2. Selection, source-carrier and Koda consistency errors propagate;
an invalid integer selection can raise `ArgumentError`.
The caller supplies a consistent unmodified forgotten record.

# Returns
`KodaConversionResult` with a `FlowSpine`, `source_status == :constructed`, choices,
source correspondence, and sector topology attached from the forgotten
complement. `spine_path` records which available carrier was used.

# Guarantees
1. The selected region and remaining crossings determine the returned
   source correspondence and sector records. [proved/agent: by construction]
2. The intended result is the branched spine of the punctured 3-manifold
   described by §§5–6, retaining nondisk sectors where present.
   [unverified: general geometric realization and manifold identification of the emitted record remain uncertified]
3. The forgotten input is not modified; the result can retain it as source
   data. [proved/agent: by construction]

# Not guaranteed
Vertex-free spines and nondisk sectors need not admit `spine_to_odata`.
Constructed status and finite correspondence checks do not independently
establish the manifold identification in guarantee 2.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5.2–6.4 — selected joining disk, singular circles and sector topology.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> root = normpath(joinpath(dirname(pathof(Topology)), "..", "..", ".."));

julia> h = read_hc(joinpath(root, "threads", "PJ260922", "data", "261107-lickorish-data", "group-12", "case-01.hc"));

julia> f = forget_reference(realize_hc(h));

julia> spine_topology(heegaard_to_spine(f; joining_region=3).value)["vertex_free"]
true
```

See also [`forget_reference`](@ref), [`spine_topology`](@ref), [`spine_to_odata`](@ref).

Role: convert
"""
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
    _koda_combinatorial_spine_from_punctured(_koda_original_punctured(forgot, region).value;
                                           genus=forgot.genus)
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
        "complex_euler" => euler,
        "sector_chi" => get(witness, "chi", nothing),
        "vertex_count" => length(topo["vertices"]),
        "edge_count" => length(topo["edges"]),
        "sector_chi_sum" => sum(s["chi"] for s in sectors),
    )
end

function _koda_virtual_image(spine::FlowSpine; stage::Symbol=:hc_to_odata, regularize::Bool=false)
    odata = try
        spine_to_odata(spine; regularize=regularize)
    catch err
        err isa KodaNotComputed && koda_not_computed(stage, err.reason, err.detail)
        rethrow()
    end
    odata_to_virtual_link(odata), odata
end

"""
    hc_to_odata(H::CellularHeegaardDiagram; joining_disk=nothing, choices=Dict(), regularize=false) -> KodaConversionResult

Return the VirtualLink encoding of the Koda o-graph of a punctured cellular Heegaard diagram.

# Arguments
`H` and the explicit `joining_disk` satisfy
[the punctured Koda domain](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).
A missing selection raises `KodaInputError`; an integer selects a checked
region and a `JoiningDisk` must belong to the caller's diagram.
Punctured, spine and encoding checks propagate. Unsupported non-disk sectors
raise `KodaNotComputed`; `regularize=true` requests the supported regularization,
not arbitrary repair. `choices` must be a compatible string-keyed dictionary.

# Returns
`KodaConversionResult` with `value::VirtualLink`, copied/extended `choices`,
source-to-spine correspondence including `odata` and joining region, and
`source_status=:constructed`.

# Guarantees
1. The value is the o-graph encoding of the constructed punctured Koda spine
   in the linked conversion scope. [proved/agent: by construction]
2. The explicit disk selection is retained in the returned choices and
   correspondence. [proved/agent: by construction]
3. The diagram and supplied choices are not modified.
   [proved/agent: by construction]

# Not guaranteed
Construction does not independently prove a general PL identification or
transport a combing/framing. The selected disk and regularization have their
own geometric hypotheses; unsupported sectors are not encoded as disks.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§2.2, 5–7 — punctured source and smoothed o-graph image.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]); h = H.certify_heegaard(H.decorate(c,[[1]]));

julia> try H.hc_to_odata(h) catch e; e isa H.KodaInputError && e.stage == :hc_to_odata end
true
```

Role: convert
"""
function hc_to_odata(H::HeegaardDiagram; joining_disk=nothing, choices=Dict{String, Any}(), regularize::Bool=false)
    disk = _selected_joining_disk(H, joining_disk)
    ph = make_punctured_combinatorial(H, disk)
    spine = heegaard_to_spine(ph)
    nondisk = _koda_nondisk_odata(spine.value)
    nondisk !== nothing && (!regularize || nondisk["complex_euler"] != 1) &&
        koda_not_computed(:hc_to_odata, :non_disk_sector, nondisk)
    link, odata = _koda_virtual_image(spine.value; stage=:hc_to_odata, regularize=regularize)
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

function _hc_disk_selection(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing)
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

"""cf: docs/specs/heegaard-koda-conversion.md §4.2 — preserve the stored LC face support.

Resolve a disk name only when its complete support is one forgotten region.
Integer IDs and numeric strings retain the existing region-selection precedence.
"""
function _resolve_hc_named_disk(hc::HeegaardCode, forgot::ForgottenHeegaard, selection)
  selection isa Integer && return selection
  key = string(selection)
  any(r -> string(r.id) == key, forgot.regions) && return selection
  haskey(hc.disks, key) || return selection
  wanted = Set(string.(hc.disks[key]))
  known = Set{String}()
  for r in forgot.regions
    union!(known, string.(r.face_roots))
  end
  missing = sort!(collect(setdiff(wanted, known)))
  isempty(missing) || throw(KodaInputError(
    :hc_to_odata, "named HC disk contains unknown LC faces",
    Dict("disk" => key, "missing_faces" => missing),
  ))
  ids = [r.id for r in forgot.regions if Set(string.(r.face_roots)) == wanted]
  length(ids) == 1 && return only(ids)
  msg = isempty(ids) ? "named HC disk is not exactly one forgotten region" :
    "named HC disk is ambiguous after forgetting reference"
  throw(KodaInputError(:hc_to_odata, msg,
    Dict("disk" => key, "faces" => sort!(collect(wanted)), "regions" => ids)))
end

"""
    hc_to_odata(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, choices=Dict(), regularize=false) -> KodaConversionResult

Return the VirtualLink encoding of the Koda spine associated with a resolved HC and joining disk.

# Arguments
`hc` satisfies the [forward Koda domain](../specs/heegaard-koda-conversion.md).
Disk-selection keywords follow the joining-disk selection convention on that
page. Ambiguous disk lines, unknown faces, nonregional disk support, absent
candidates or failed joining conditions raise `KodaInputError`. Other
realization and spine errors propagate. Unsupported sectors or vertex-free
encoding raise `KodaNotComputed`; `regularize=true` requests only supported
regularization. `choices` must be a mutable compatible string-keyed dictionary.

# Returns
`KodaConversionResult` with `value::VirtualLink`, updated `choices`, o-data
and source/spine correspondences, and `source_status=:constructed`.

# Guarantees
1. Successful return contains the o-graph encoding of the constructed
   punctured spine and records the selected joining region.
   [proved/agent: by construction]
2. Named disk support is matched to the whole forgotten region.
   [proved/agent: by construction]
3. The HC is not modified; the supplied `choices` dictionary is updated and
   retained in the result, including on some failing paths.
   [proved/agent: by construction]

# Not guaranteed
Successful construction is not an independent general PL homeomorphism proof.
Automatic disk selection does not preserve additional combing or Spinᶜ data.
Regularization is not a general repair procedure.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§4.2, 5.2, 6–7 — disk support, selected puncture and o-graph output.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(-1,Dict(),Dict(),Dict(),Dict(),Dict(),nothing); hc = H.HeegaardCode(lc,Dict(),Dict(),Dict(),Dict(),nothing);

julia> try H.hc_to_odata(hc) catch e; e isa H.KodaInputError end
true
```

Role: convert
"""
function hc_to_odata(hc::HeegaardCode; joining_disk=nothing, disk_name=nothing, disk_id=nothing, choices=Dict{String, Any}(), regularize::Bool=false)
  validate_hc_input(hc)
  realized = realize_hc(hc)
  forgot = forget_reference(realized)
  sel = _hc_disk_selection(hc; joining_disk=joining_disk, disk_name=disk_name, disk_id=disk_id)
  if sel === nothing
    # cf: docs/specs/heegaard-koda-conversion.md §5.2 — minimum region_id; empty stays unselected.
    cands = enumerate_joining_disk_candidates(forgot)
    if !isempty(cands)
      sel = minimum(c["region_id"] for c in cands)
    end
  end
  if sel === nothing
    throw(KodaInputError(
      :hc_to_odata,
      "joining disk not selected; pass joining_disk or list a candidate seed",
      Dict("candidates" => enumerate_joining_disk_candidates(forgot)),
    ))
  end
  requested_sel = sel
  sel = _resolve_hc_named_disk(hc, forgot, sel)
  region = _resolve_joining_region(forgot, sel)
  merge!(choices, Dict(
    "joining_disk_selection" => requested_sel,
    "joining_disk_region" => region.id,
  ))
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
    punctured = _koda_original_punctured(forgot, region; choices=choices)
    punctured_corr = punctured.correspondence
    spine_path = "ribbon"
    _koda_combinatorial_spine_from_punctured(punctured.value; genus=forgot.genus)
  else
    region === nothing && throw(KodaInputError(:hc_to_odata, "joining region missing", Dict("selection" => sel)))
    spine_path = "ports"
    diagram_path_error = get(forgot.correspondence, "diagram_build_error", "diagram not built")
    _koda_combinatorial_spine_from_forgotten_region(forgot, region)
  end
  region !== nothing && _koda_attach_sectors!(spine.value, forgot, region)
  nondisk = _koda_nondisk_odata(spine.value)
  if nondisk !== nothing && (!regularize || nondisk["complex_euler"] != 1)
    koda_not_computed(:hc_to_odata, :non_disk_sector, nondisk)
  end
  if get(spine.correspondence, "vertex_free", false) && !regularize
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
  link, odata = _koda_virtual_image(spine.value; stage=:hc_to_odata, regularize=regularize)
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
  if haskey(odata.metadata, "regularization")
    corr["regularization"] = odata.metadata["regularization"]
    choices["regularization"] = "annulus_beta_finger"
  end
  KodaConversionResult(link, choices, corr, :constructed)
end

"""
    hc_to_odata(lc_path::AbstractString, hc_path::AbstractString; kwargs...) -> KodaConversionResult

Read an HC file and return its forward Koda o-graph conversion.

# Arguments
`hc_path` must be readable and accepted by `read_hc`; the LC named by that
HC file must also be readable. `lc_path` is currently ignored.
File, parse and [forward conversion](../specs/heegaard-koda-conversion.md)
errors propagate. Keywords have the domain of `hc_to_odata(::HeegaardCode)`.

# Returns
The `KodaConversionResult` of the HC loaded from `hc_path`, whose value is
the constructed `VirtualLink`.

# Guarantees
1. The result is the forward conversion of the HC file and its declared LC
   dependency, independent of `lc_path`. [proved/agent: by construction]
2. Input files are read without modification; a supplied mutable `choices`
   keyword can be updated by the conversion. [proved/agent: by construction]

# Not guaranteed
This overload does not override or validate the HC file's LC reference
against `lc_path`; see [the file-argument discrepancy](../specs/heegaard-koda-conversion.md#forward-file-lc-argument). The underlying conversion's geometric limits apply.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md — file-to-HC forward conversion domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> mktempdir() do dir; try H.hc_to_odata("ignored.lc",joinpath(dir,"absent.hc")) catch e; e isa SystemError end; end
true
```

Role: convert
"""
function hc_to_odata(lc_path::AbstractString, hc_path::AbstractString; kwargs...)
  hc_to_odata(read_hc(hc_path); kwargs...)
end
