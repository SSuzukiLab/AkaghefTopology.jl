"""Construction and schema-versioned JSON interchange."""

function _spec_field(spec, key::Symbol, default=nothing)
    if spec isa NamedTuple
        return hasproperty(spec, key) ? getproperty(spec, key) : default
    elseif spec isa AbstractDict
        return haskey(spec, key) ? spec[key] : haskey(spec, string(key)) ? spec[string(key)] : default
    else
        return hasproperty(spec, key) ? getproperty(spec, key) : default
    end
end

function _curve_from_spec(spec, family::Symbol, position::Integer)
    if spec isa Curve
        spec.family == family || throw(ArgumentError("curve family does not match its container"))
        return Curve(family, spec.id, spec.visits; label=spec.label, tags=spec.tags)
    end
    id = _spec_field(spec, :id, position)
    label = _spec_field(spec, :label, nothing)
    visits = spec isa AbstractVector ? spec : _spec_field(spec, :visits, Int[])
    tags = _spec_field(spec, :tags, nothing)
    Curve(family, id, visits; label=label, tags=tags)
end

function _intersection_from_spec(spec, position::Integer)
    spec isa Intersection && return Intersection(spec.id, spec.alpha, spec.beta, spec.sign;
                                                 source_label=spec.source_label)
    Intersection(_spec_field(spec, :id, position),
                 _spec_field(spec, :alpha, _spec_field(spec, :alpha_curve, 0)),
                 _spec_field(spec, :beta, _spec_field(spec, :beta_curve, 0)),
                 _spec_field(spec, :sign, 0);
                 source_label=_spec_field(spec, :source_label, nothing))
end

"""Build an undecorated signed cyclic-word code.

`alpha` and `beta` are vectors of visit vectors, `Curve`s, or curve-like
named tuples/dictionaries.  `intersections` is a vector of `Intersection`s or
objects with `id`, `alpha`, `beta`, and `sign` fields.
"""
function heegaard_code(alpha::AbstractVector, beta::AbstractVector,
                       intersections::AbstractVector; metadata=nothing)
    a = [_curve_from_spec(spec, :alpha, i) for (i, spec) in enumerate(alpha)]
    b = [_curve_from_spec(spec, :beta, i) for (i, spec) in enumerate(beta)]
    p = [_intersection_from_spec(spec, i) for (i, spec) in enumerate(intersections)]
    HeegaardCode(a, b, p, _copy_meta(metadata))
end

function heegaard_code(; alpha=Vector{Any}(), beta=Vector{Any}(), intersections=Vector{Any}(), metadata=nothing)
    heegaard_code(alpha, beta, intersections; metadata=metadata)
end

function _code_of(object)
    object isa HeegaardDiagram && return object.code
    object isa DecoratedHeegaardCode && return object
    object isa HeegaardCode && return object
    throw(ArgumentError("expected HeegaardCode, DecoratedHeegaardCode, or HeegaardDiagram"))
end

function _undecorated(object)
    object isa HeegaardDiagram && return object.code.code
    object isa DecoratedHeegaardCode && return object.code
    object isa HeegaardCode && return object
    throw(ArgumentError("expected Heegaard code"))
end

"""Attach a complete boundary-cycle partition to a code."""
function decorate(code::HeegaardCode, regions; metadata=nothing, cellular::Bool=false)
    if regions === nothing
        cellular || throw(ArgumentError("region decoration is required; pass cellular=true explicitly"))
        regions = [[i] for i in eachindex(boundary_cycles(code))]
    end
    partition = RegionPartition(regions)
    DecoratedHeegaardCode(code, partition,
        _copy_meta(metadata === nothing ? code.metadata : metadata))
end

function decorate(code::HeegaardCode; regions=nothing, metadata=nothing, cellular::Bool=false)
    decorate(code, regions; metadata=metadata, cellular=cellular)
end

function decorate(code::DecoratedHeegaardCode, regions; metadata=nothing, cellular::Bool=false)
    decorate(code.code, regions; metadata=metadata === nothing ? code.metadata : metadata, cellular=cellular)
end

function _curve_json(curve::Curve)
    Dict{String,Any}(
        "id" => curve.id,
        "label" => curve.label,
        "visits" => copy(curve.visits),
        "tags" => deepcopy(curve.tags),
    )
end

function _intersection_json(point::Intersection)
    result = Dict{String,Any}(
        "id" => point.id,
        "alpha" => point.alpha,
        "beta" => point.beta,
        "sign" => point.sign,
    )
    point.source_label === nothing || (result["source_label"] = point.source_label)
    result
end

function _json_dict(object::HeegaardCode)
    Dict{String,Any}(
        "schema" => HGD_SCHEMA,
        "curves" => Dict{String,Any}(
            "alpha" => [_curve_json(curve) for curve in object.alpha],
            "beta" => [_curve_json(curve) for curve in object.beta],
        ),
        "intersections" => [_intersection_json(point) for point in object.intersections],
        "metadata" => deepcopy(object.metadata),
    )
end

function _json_dict(object::DecoratedHeegaardCode)
    result = _json_dict(object.code)
    result["regions"] = Dict{String,Any}(
        "boundary_cycle_blocks" => [copy(block) for block in object.regions.boundary_cycle_blocks],
    )
    result["metadata"] = deepcopy(object.metadata)
    result
end

_json_dict(object::HeegaardDiagram) = _json_dict(object.code)

function _as_vector(x)
    x === nothing && return Any[]
    collect(x)
end

"""Read a `.hgd.json` file and return an undecorated or decorated code."""
function load_heegaard(path::AbstractString)
    parsed = JSON3.read(read(path, String))
    schema = string(_spec_field(parsed, :schema, ""))
    startswith(schema, "akaghef.heegaard/") || throw(ArgumentError("unknown Heegaard schema: $schema"))
    split(schema, "/")[end] == "1" || throw(ArgumentError("unsupported Heegaard schema major version: $schema"))
    curves = _spec_field(parsed, :curves, nothing)
    curves === nothing && throw(ArgumentError("missing curves object"))
    alpha = _as_vector(_spec_field(curves, :alpha, nothing))
    beta = _as_vector(_spec_field(curves, :beta, nothing))
    intersections = _as_vector(_spec_field(parsed, :intersections, nothing))
    code = heegaard_code(alpha, beta, intersections; metadata=_spec_field(parsed, :metadata, nothing))
    regions = _spec_field(parsed, :regions, nothing)
    regions === nothing && return code
    blocks = _spec_field(regions, :boundary_cycle_blocks, nothing)
    blocks === nothing && throw(ArgumentError("regions must contain boundary_cycle_blocks"))
    decorate(code, _as_vector(blocks); metadata=code.metadata)
end

"""Write canonical JSON for a Heegaard core object."""
function save_heegaard(path::AbstractString, object)
    object isa Union{HeegaardCode,DecoratedHeegaardCode,HeegaardDiagram} ||
        throw(ArgumentError("save_heegaard expects a Heegaard core object"))
    open(path, "w") do io
        write(io, _canonical_json(_json_dict(object)))
        write(io, '\n')
    end
    path
end

source_manifest(object) = deepcopy(_undecorated(object).metadata)
