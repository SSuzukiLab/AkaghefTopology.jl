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

"""
    heegaard_code(alpha::AbstractVector, beta::AbstractVector, intersections::AbstractVector; metadata=nothing) -> CombinatorialHeegaardCode

Build an undecorated legacy combinatorial Heegaard code from curve and intersection data.

# Arguments
`alpha` and `beta` supply vectors of visits, `Curve` values, or curve-like records; `intersections` supplies `Intersection` values or intersection-like records in the [legacy representation](../specs/heegaard-terminology.md#legacy-cellular-records). Existing `Curve` values must match their containing family or raise `ArgumentError`. Record conversions and the `Curve`/`Intersection` constructor checks propagate. Optional `metadata` is copied. Global incidence and surface conditions are checked separately by `validate`.

# Returns
A `CombinatorialHeegaardCode` containing newly constructed curve and intersection vectors and copied metadata.

# Guarantees
1. The supplied visit order and intersection data are retained in the legacy signed cyclic-word representation. [proved/agent: by construction]
2. Inputs are not modified; the result records a proposed curve arrangement, without complementary-region decoration. [proved/agent: by construction]

# Not guaranteed
This construction does not certify global consistency, a closed surface, Heegaard cut systems or a represented 3-manifold. Use `decorate` and `certify_heegaard` for the legacy validation path.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cellular-records — legacy curve, intersection and decoration records.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); (length(c.alpha), length(c.beta), length(c.intersections))
(1, 1, 1)
```

Role: construct
"""
function heegaard_code(alpha::AbstractVector, beta::AbstractVector,
                       intersections::AbstractVector; metadata=nothing)
    a = [_curve_from_spec(spec, :alpha, i) for (i, spec) in enumerate(alpha)]
    b = [_curve_from_spec(spec, :beta, i) for (i, spec) in enumerate(beta)]
    p = [_intersection_from_spec(spec, i) for (i, spec) in enumerate(intersections)]
    CombinatorialHeegaardCode(a, b, p, _copy_meta(metadata))
end

"""
    heegaard_code(; alpha=Vector{Any}(), beta=Vector{Any}(), intersections=Vector{Any}(), metadata=nothing) -> CombinatorialHeegaardCode

Build an undecorated legacy combinatorial Heegaard code from keyword data.

# Arguments
`alpha` and `beta` supply vectors of visits, `Curve` values, or curve-like records; `intersections` supplies `Intersection` values or intersection-like records in the [legacy representation](../specs/heegaard-terminology.md#legacy-cellular-records). Existing `Curve` values must match their containing family or raise `ArgumentError`. Record conversions and the `Curve`/`Intersection` constructor checks propagate. Optional `metadata` is copied. Global incidence and surface conditions are checked separately by `validate`. Omitted collections are empty.

# Returns
The same record as the positional `heegaard_code(alpha,beta,intersections; metadata)` call.

# Guarantees
1. The positional overload's representation and input-preservation guarantees apply. [proved/agent: by construction]

# Not guaranteed
Empty defaults and successfully constructed records are not certificates of a Heegaard diagram or of a 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cellular-records — legacy combinatorial input.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code(); (length(c.alpha), length(c.beta), length(c.intersections))
(0, 0, 0)
```

Role: construct
"""
function heegaard_code(; alpha=Vector{Any}(), beta=Vector{Any}(), intersections=Vector{Any}(), metadata=nothing)
    heegaard_code(alpha, beta, intersections; metadata=metadata)
end

function _code_of(object)
    object isa HeegaardDiagram && return object.code
    object isa DecoratedHeegaardCode && return object
    object isa CombinatorialHeegaardCode && return object
    throw(ArgumentError("expected CombinatorialHeegaardCode, DecoratedHeegaardCode, or HeegaardDiagram"))
end

function _undecorated(object)
    object isa HeegaardDiagram && return object.code.code
    object isa DecoratedHeegaardCode && return object.code
    object isa CombinatorialHeegaardCode && return object
    throw(ArgumentError("expected Heegaard code"))
end

"""
    decorate(code::CombinatorialHeegaardCode, regions; metadata=nothing, cellular=false) -> DecoratedHeegaardCode

Attach complementary-region boundary data to a legacy combinatorial code.

# Arguments
`regions` specifies the partition in [legacy coordinates](../specs/heegaard-terminology.md#legacy-cellular-records).
If it is `nothing`, `cellular=true` is required, otherwise `ArgumentError`
requests explicit decoration. Region conversion and `RegionPartition`
checks propagate. The caller checks coverage and geometric consistency
through `validate`; this operation does not fully certify the partition.
Optional metadata replaces the decoration metadata; the code's metadata is
used when omitted.

# Returns
A `DecoratedHeegaardCode` retaining the code, with a `RegionPartition` and
copied metadata. Explicit cellular mode with no partition assigns disk regions.

# Guarantees
1. The returned decoration contains the supplied or explicitly requested
   cellular partition in the linked representation.
   [proved/agent: by construction]
2. The input code is not modified and is shared with the result.
   [proved/agent: by construction]

# Not guaranteed
Constructing a decoration does not certify its surface or Heegaard conditions.
Using cellular mode asserts the intended disk-complement interpretation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy-cellular-records — complementary regions and boundary partitions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]); d = H.decorate(c, [[1]]); Tuple(only(d.regions.boundary_cycle_blocks))
(1,)
```

See also [`certify_heegaard`](@ref).

Role: convert
"""
function decorate(code::CombinatorialHeegaardCode, regions; metadata=nothing, cellular::Bool=false)
    if regions === nothing
        cellular || throw(ArgumentError("region decoration is required; pass cellular=true explicitly"))
        regions = [[i] for i in eachindex(boundary_cycles(code))]
    end
    partition = RegionPartition(regions)
    DecoratedHeegaardCode(code, partition,
        _copy_meta(metadata === nothing ? code.metadata : metadata))
end

"""
    decorate(code::CombinatorialHeegaardCode; regions=nothing, metadata=nothing, cellular=false) -> DecoratedHeegaardCode

Attach an explicit or caller-asserted cellular region partition using keyword arguments.

# Arguments
`code`, `regions`, `metadata` and `cellular` have the positional `decorate`
overload's requirements. Omitting `regions` requires `cellular=true`, checked
with `ArgumentError`; partition and ribbon errors propagate.

# Returns
The same decorated record as `decorate(code, regions; metadata, cellular)`.

# Guarantees
1. The positional overload's partition and metadata guarantees apply.
   [proved/agent: by construction]
2. The code is retained by reference without mutation.
   [proved/agent: by construction]

# Not guaranteed
The caller's cellular assertion does not replace full validation of
[the legacy surface](../specs/heegaard-terminology.md#legacy-cellular-records).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]); length(H.decorate(c; cellular=true).regions.boundary_cycle_blocks)
1
```

See also [`certify_heegaard`](@ref).

Role: convert
"""
function decorate(code::CombinatorialHeegaardCode; regions=nothing, metadata=nothing, cellular::Bool=false)
    decorate(code, regions; metadata=metadata, cellular=cellular)
end

"""
    decorate(code::DecoratedHeegaardCode, regions; metadata=nothing, cellular=false) -> DecoratedHeegaardCode

Replace the region decoration of a legacy decorated code.

# Arguments
`regions` and `cellular` follow the combinatorial `decorate` overload's
requirements and errors. Optional metadata replaces the decoration metadata;
otherwise the existing decoration's metadata is retained by copying.
The caller validates the new region interpretation separately.

# Returns
A `DecoratedHeegaardCode` with the same underlying combinatorial code,
the newly specified partition and copied metadata.

# Guarantees
1. The new record carries the requested region decoration of the original
   curve/intersection data. [proved/agent: by construction]
2. The original decoration and code are not modified; the code is shared.
   [proved/agent: by construction]

# Not guaranteed
Changing the partition can change the represented surface and manifold;
this operation does not assert equivalence to the previous decoration.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]); d = H.decorate(c,[[1]]); H.decorate(d,[[1]]).code === d.code
true
```

See also [`validate`](@ref).

Role: move
"""
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

function _json_dict(object::CombinatorialHeegaardCode)
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

"""
    load_heegaard(path::AbstractString) -> Union{CombinatorialHeegaardCode,DecoratedHeegaardCode}

Read a legacy combinatorial Heegaard record from its JSON file.

# Arguments
`path` names readable JSON in the legacy schema of [the cellular records](../specs/heegaard-terminology.md#legacy-cellular-records). Unknown schema, unsupported major version, missing curve data or missing region blocks raises `ArgumentError`. I/O, JSON parsing and record-conversion errors propagate.

# Returns
An undecorated code when no region object is present, otherwise a decorated code. Curve order, visits, intersection data and supported metadata are retained.

# Guarantees
1. The result reconstructs the stored legacy combinatorial record and its optional region partition. [proved/agent: by construction]
2. The input file is read without modification. [proved/agent: by construction]

# Not guaranteed
Loading does not certify global diagram validity, complete cut systems or a represented 3-manifold. It does not load the modern two-file LC/HC format or arbitrary marked dictionaries; use `read_hc` for that format.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cellular-records — legacy record semantics.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> mktemp() do path, io
           close(io)
           c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)])
           H.save_heegaard(path,c)
           r = H.load_heegaard(path)
           (length(r.alpha), length(r.intersections))
       end
(1, 1)
```

Role: construct
"""
function load_heegaard(path::AbstractString)
    parsed = JSON3.read(read(path, String))
    _heegaard_from_json(parsed)
end

function _heegaard_from_json(parsed)
    parsed isa AbstractDict || throw(ArgumentError("Heegaard JSON must contain an object"))
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

"""
    save_heegaard(path::AbstractString, object) -> AbstractString

Write a legacy Heegaard record as schema-versioned JSON.

# Arguments
`object` is a `CombinatorialHeegaardCode`, `DecoratedHeegaardCode` or `CellularHeegaardDiagram`; other values raise `ArgumentError`. The caller provides serializable metadata and a writable destination with an existing parent directory. Serialization and I/O failures propagate.

# Returns
The supplied path after writing the record and a trailing newline.

# Guarantees
1. Supported curve, intersection and region data are serialized in their stored order for `load_heegaard`; a cellular wrapper is saved as its decorated record. [proved/agent: by construction]
2. The object is not modified. The destination file is created or overwritten. [proved/agent: by construction]

# Not guaranteed
No new mathematical validation or atomic file replacement is performed. An I/O failure can leave a partial destination. This is the legacy JSON format, not the two-file LC/HC format.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cellular-records — saved legacy representation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> mktemp() do path, io
           close(io)
           c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)])
           H.save_heegaard(path,c) == path && H.load_heegaard(path).alpha[1].visits == [1]
       end
true
```

Role: convert
"""
function save_heegaard(path::AbstractString, object)
    object isa Union{CombinatorialHeegaardCode,DecoratedHeegaardCode,HeegaardDiagram} ||
        throw(ArgumentError("save_heegaard expects a Heegaard core object"))
    data = _json_dict(object)
    # canonicalization.jl::_canonical_json preserves array/word order and opaque
    # JSON metadata; the vendor exact-math writer does not support float metadata.
    encoded = _canonical_json(data)
    open(path, "w") do io
        write(io, encoded)
        write(io, '\n')
    end
    path
end

"""
    source_manifest(object) -> Dict{String,Any}

Return a copy of the source metadata attached to a legacy combinatorial code.

# Arguments
`object` is a `CombinatorialHeegaardCode`, `DecoratedHeegaardCode` or
`CellularHeegaardDiagram`. Other inputs raise `ArgumentError` requiring a
Heegaard core object. Metadata contents are supplied annotations and are
not verified against their cited sources.

# Returns
A deep copy of the underlying combinatorial code's metadata dictionary.

# Guarantees
1. The result contains the underlying code's source annotations, rather than
   the enclosing decoration's separate metadata.
   [proved/agent: by construction]
2. The input is not modified; ordinary nested metadata containers are copied.
   [proved/agent: by construction]

# Not guaranteed
Recorded provenance or an evidence-status string is not a proof of its claim.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); H.source_manifest(c) == c.metadata
true
```

Role: access
"""
source_manifest(object) = deepcopy(_undecorated(object).metadata)
