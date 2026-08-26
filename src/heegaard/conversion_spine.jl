"""Typed, certificate-carrying Heegaard to flow-spine boundary."""

function _conversion_snapshot(H)
    h = first_homology(H)
    Dict{String,Any}("h1_free_rank"=>h.free_rank, "h1_torsion"=>copy(h.torsion),
                     "surface_genus"=>surface_genus(H))
end

function make_punctured(H::HeegaardDiagram; strategy=:koda)
    strategy == :koda || throw(ArgumentError("only the :koda puncturing strategy is implemented"))
    cut = cut_along_alpha(H)
    joining_disk = Dict{String,Any}(
        "strategy" => "koda",
        "chart" => "sphere-at-infinity",
        "selected_region" => 1,
        "selection_is_canonical" => false,
    )
    certificate = ConversionCertificate(:koda, canonical_hash(H), joining_disk, Any[],
        "oriented alpha/beta ribbon convention", Dict{String,Any}(), nothing,
        _conversion_snapshot(H), :unverified)
    PuncturedHeegaardDiagram(H, joining_disk, cut, certificate)
end

function to_flow_spine(PH::PuncturedHeegaardDiagram)
    H = PH.diagram
    branches = NamedTuple{(:curve_family, :curve_id, :intersection_ids),Tuple{Symbol,Int,Vector{Int}}}[]
    for curve in vcat(H.code.code.alpha, H.code.code.beta)
        push!(branches, (curve_family=curve.family, curve_id=curve.id,
                         intersection_ids=copy(curve.visits)))
    end
    certificate = ConversionCertificate(:koda, canonical_hash(H), deepcopy(PH.joining_disk),
        Any[], "oriented branch follows curve orientation", Dict{String,Any}(
            "alpha_curves"=>[curve.id for curve in H.code.code.alpha],
            "beta_curves"=>[curve.id for curve in H.code.code.beta]), nothing,
        _conversion_snapshot(H), :unverified)
    FlowSpine(surface_genus(H), branches, PH, certificate)
end

function to_ograph(spine::FlowSpine)
    metadata = Dict{String,Any}(
        "representation" => "closed-normal-o-graph boundary object",
        "canonical" => false,
        "source_hash" => spine.certificate.source_hash,
        "branch_count" => length(spine.branches),
    )
    # OData is intentionally not VirtualLink: each Heegaard curve is a
    # branch word, not a knot crossing word.  Downstream conversion into a
    # VirtualLink needs a future explicit standardization choice.
    OData([copy(branch.intersection_ids) for branch in spine.branches],
          [branch.curve_family == :alpha ? 1 : -1 for branch in spine.branches], metadata)
end

function to_heegaard(PH::PuncturedHeegaardDiagram; strategy=:block)
    strategy == :block || throw(ArgumentError("only the :block reverse strategy is implemented"))
    PH.diagram
end

function to_heegaard(spine::FlowSpine; strategy=:block)
    spine.source isa PuncturedHeegaardDiagram ||
        throw(ArgumentError("flow spine has no retained punctured Heegaard source"))
    to_heegaard(spine.source; strategy=strategy)
end

function _ograph_hash(ograph::OData)
    bytes2hex(SHA.sha256(_canonical_json(Dict("gauss"=>ograph.gauss,
                                             "orientation"=>ograph.orientation,
                                             "metadata"=>ograph.metadata))))
end

function compare_representations(H::HeegaardDiagram, ograph::OData)
    expected = _conversion_snapshot(H)
    actual = get(ograph.metadata, "h1", nothing)
    result = Dict{String,Any}(
        "source_hash" => canonical_hash(H),
        "ograph_hash" => _ograph_hash(ograph),
        "h1" => expected,
        "comparison" => actual === nothing ? :not_available : :available,
        "canonical_conversion" => false,
    )
    result
end

function to_odata(H::HeegaardDiagram; strategy=:koda)
    punctured = make_punctured(H; strategy=strategy)
    spine = to_flow_spine(punctured)
    ograph = to_ograph(spine)
    ograph.metadata["h1"] = _conversion_snapshot(H)
    comparison = compare_representations(H, ograph)
    certificate = ConversionCertificate(strategy, canonical_hash(H),
        deepcopy(punctured.joining_disk), Any[],
        "oriented branch follows curve orientation",
        Dict{String,Any}("flow_spine_branches"=>length(spine.branches)),
        _ograph_hash(ograph), comparison, :unverified)
    ODataBundle(punctured, spine, ograph, certificate, comparison)
end

function _conversion_json(certificate::ConversionCertificate)
    Dict{String,Any}(
        "schema" => "akaghef.heegaard.hgd2ograph/1",
        "strategy" => String(certificate.strategy),
        "source_diagram_hash" => certificate.source_hash,
        "joining_disk" => deepcopy(certificate.joining_disk),
        "isotopy_moves" => deepcopy(certificate.isotopy_moves),
        "branch_convention" => certificate.branch_convention,
        "label_map" => deepcopy(certificate.label_map),
        "generated_hash" => certificate.generated_hash,
        "comparison_invariants" => deepcopy(certificate.comparison_invariants),
        "evidence_status" => String(certificate.evidence_status),
    )
end

function save_conversion(path::AbstractString, value)
    certificate = value isa ODataBundle ? value.certificate : value
    certificate isa ConversionCertificate || throw(ArgumentError("save_conversion expects a conversion certificate or bundle"))
    endswith(path, ".hgd2ograph.json") || throw(ArgumentError("conversion path must end in .hgd2ograph.json"))
    open(path, "w") do io
        write(io, _canonical_json(_conversion_json(certificate)))
        write(io, '\n')
    end
    path
end
