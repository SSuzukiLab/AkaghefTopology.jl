"""Immutable-by-convention value types for the Heegaard subsystem.

The vectors stored in these structs are copied at construction time.  The
public API does not mutate them; transformations return a fresh value and a
certificate.  Julia cannot make a `Vector` field deeply immutable, so callers
that need a hard boundary should treat these values as persistent records.
"""

const HGD_SCHEMA = "akaghef.heegaard/1"

function _plain_value(value)
    value isa AbstractDict && return Dict{String,Any}(string(k) => _plain_value(v) for (k, v) in value)
    value isa AbstractVector && return [_plain_value(item) for item in value]
    value isa Tuple && return tuple((_plain_value(item) for item in value)...)
    value
end

_copy_meta(x) = x === nothing ? Dict{String,Any}() :
    Dict{String,Any}(string(k) => _plain_value(v) for (k, v) in x)

function _point2(value)
    values = collect(value)
    length(values) == 2 || throw(ArgumentError("2D coordinate must have two entries"))
    (Float64(values[1]), Float64(values[2]))
end

struct Curve
    family::Symbol
    id::Int
    label::String
    visits::Vector{Int}
    tags::Dict{String,Any}

    function Curve(family::Symbol, id::Integer, visits::AbstractVector{<:Integer};
                   label=nothing, tags=nothing)
        family in (:alpha, :beta) || throw(ArgumentError("curve family must be :alpha or :beta"))
        id > 0 || throw(ArgumentError("curve id must be positive"))
        new(family, Int(id), label === nothing ? string(family == :alpha ? "a" : "b", id) : string(label),
            Int.(collect(visits)), _copy_meta(tags))
    end
end

struct Intersection
    id::Int
    alpha::Int
    beta::Int
    sign::Int
    source_label::Union{Nothing,String}

    function Intersection(id::Integer, alpha::Integer, beta::Integer, sign::Integer;
                          source_label=nothing)
        id > 0 || throw(ArgumentError("intersection id must be positive"))
        alpha > 0 || throw(ArgumentError("alpha curve id must be positive"))
        beta > 0 || throw(ArgumentError("beta curve id must be positive"))
        sign in (-1, 1) || throw(ArgumentError("intersection sign must be +1 or -1"))
        new(Int(id), Int(alpha), Int(beta), Int(sign),
            source_label === nothing ? nothing : string(source_label))
    end
end

struct HeegaardCode
    alpha::Vector{Curve}
    beta::Vector{Curve}
    intersections::Vector{Intersection}
    metadata::Dict{String,Any}
end

struct RegionPartition
    boundary_cycle_blocks::Vector{Vector{Int}}
    function RegionPartition(blocks::AbstractVector)
        new([Int.(collect(block)) for block in blocks])
    end
end

struct DecoratedHeegaardCode
    code::HeegaardCode
    regions::RegionPartition
    metadata::Dict{String,Any}
end

struct ValidationIssue
    level::Symbol
    code::Symbol
    message::String
    path::String
end

struct ValidationReport
    level::Symbol
    valid::Bool
    issues::Vector{ValidationIssue}
    checks::Dict{Symbol,Bool}
end

struct ValidationError <: Exception
    report::ValidationReport
end

function Base.showerror(io::IO, err::ValidationError)
    print(io, "invalid Heegaard diagram")
    for issue in err.report.issues
        print(io, "; ", issue.level, "/", issue.code, ": ", issue.message)
    end
end

struct HeegaardDiagram
    code::DecoratedHeegaardCode
    validation::ValidationReport
    function HeegaardDiagram(code::DecoratedHeegaardCode, report::ValidationReport, ::Val{:certified})
        report.valid || throw(ArgumentError("only a valid report can certify a Heegaard diagram"))
        new(code, report)
    end
end

struct RibbonEdge
    id::Int
    family::Symbol
    curve_id::Int
    from::Int
    to::Int
    from_port::Symbol
    to_port::Symbol
end

struct RibbonDart
    id::Int
    edge::Int
    origin::Int
    target::Int
    origin_port::Symbol
    target_port::Symbol
end

struct RibbonGraph
    vertices::Vector{Int}
    edges::Vector{RibbonEdge}
    darts::Vector{RibbonDart}
    rotation::Dict{Int,Vector{Symbol}}
    next_dart::Vector{Int}
end

struct ComplementRegion
    id::Int
    boundary_cycles::Vector{Int}
    genus::Int
end

struct GroupPresentation
    generators::Vector{String}
    relators::Vector{Vector{Tuple{Int,Int}}}
    beta_ids::Vector{Int}
    metadata::Dict{String,Any}
end

struct Puncture
    id::Int
    pair::Int
    side::Symbol
    center::NTuple{2,Float64}
    radius::Float64
    orientation::Int
end

struct CutBoundary
    id::Int
    pair::Int
    side::Symbol
    curve_id::Int
    orientation::Int
end

struct CutSurface
    diagram_hash::String
    genus::Int
    boundaries::Vector{CutBoundary}
    alpha_orders::Dict{Int,Vector{Int}}
    beta_orders::Dict{Int,Vector{Int}}
    intersection_signs::Dict{Int,Int}
    intersection_betas::Dict{Int,Int}
end

struct Port
    id::Int
    intersection_id::Int
    boundary_id::Int
    beta_id::Int
    parameter::Float64
    point::NTuple{2,Float64}
end

abstract type RouteCurve end

struct PolylineCurve <: RouteCurve
    points::Vector{NTuple{2,Float64}}
    function PolylineCurve(points::AbstractVector)
        normalized = [(Float64(point[1]), Float64(point[2])) for point in points]
        length(normalized) >= 2 || throw(ArgumentError("a route curve needs at least two points"))
        all(isfinite(point[1]) && isfinite(point[2]) for point in normalized) ||
            throw(ArgumentError("route curve contains a non-finite point"))
        new(normalized)
    end
end

struct CubicBezierSegment
    p0::NTuple{2,Float64}
    p1::NTuple{2,Float64}
    p2::NTuple{2,Float64}
    p3::NTuple{2,Float64}
    function CubicBezierSegment(points::NTuple{4,NTuple{2,Float64}})
        all(isfinite(point[1]) && isfinite(point[2]) for point in points) ||
            throw(ArgumentError("Bezier segment contains a non-finite point"))
        new(points...)
    end
end

function CubicBezierSegment(p0::AbstractVector, p1::AbstractVector,
                            p2::AbstractVector, p3::AbstractVector)
    CubicBezierSegment((_point2(p0), _point2(p1), _point2(p2), _point2(p3)))
end

struct BezierChain <: RouteCurve
    segments::Vector{CubicBezierSegment}
    flatness::Float64
    max_depth::Int
    function BezierChain(segments::AbstractVector{<:CubicBezierSegment};
                         flatness=1e-4, max_depth=20)
        normalized = collect(segments)
        isempty(normalized) && throw(ArgumentError("Bezier chain needs at least one segment"))
        tolerance = Float64(flatness)
        tolerance > 0 && isfinite(tolerance) || throw(ArgumentError("Bezier flatness must be positive and finite"))
        depth = Int(max_depth)
        depth >= 1 || throw(ArgumentError("Bezier max_depth must be positive"))
        for index in 1:length(normalized)-1
            normalized[index].p3 == normalized[index + 1].p0 ||
                throw(ArgumentError("Bezier chain segments are not endpoint-continuous"))
        end
        new(normalized, tolerance, depth)
    end
end

struct WiringProblem
    diagram_hash::String
    boxes::Vector{Int}
    ports::Vector{Port}
    wires::Vector{NamedTuple{(:beta_id, :from_port, :to_port),Tuple{Int,Int,Int}}}
end

struct Route
    family::Symbol
    curve_id::Int
    segment_id::Int
    start_port::Union{Nothing,Int}
    end_port::Union{Nothing,Int}
    points::Vector{NTuple{2,Float64}}
    curve::RouteCurve

    function Route(family::Symbol, curve_id::Integer, segment_id::Integer,
                   start_port::Union{Nothing,Integer}, end_port::Union{Nothing,Integer},
                   points::Vector{NTuple{2,Float64}}, curve::RouteCurve)
        new(family, Int(curve_id), Int(segment_id),
            start_port === nothing ? nothing : Int(start_port),
            end_port === nothing ? nothing : Int(end_port), points, curve)
    end
end

function Route(family::Symbol, curve_id::Integer, segment_id::Integer,
               start_port, end_port, points::AbstractVector; curve=nothing)
    normalized = [(Float64(point[1]), Float64(point[2])) for point in points]
    actual = curve === nothing ? PolylineCurve(normalized) : curve
    flattened = curve === nothing ? normalized : _flatten_curve(actual)
    normalized == flattened || throw(ArgumentError("route points do not match its curve representation"))
    Route(family, curve_id, segment_id, start_port, end_port, flattened, actual)
end

struct LayoutValidationReport
    valid::Bool
    issues::Vector{ValidationIssue}
    checks::Dict{Symbol,Bool}
end

struct HeegaardLayout
    diagram_hash::String
    view::Symbol
    chart::Dict{String,Any}
    punctures::Vector{Puncture}
    ports::Vector{Port}
    routes::Vector{Route}
    crossings::Vector{Tuple{Int,Int}}
    backend_trace::Dict{String,Any}
    report::LayoutValidationReport
end

struct MoveCertificate
    move_type::Symbol
    attachment_data::Dict{String,Any}
    label_correspondence::Dict{String,Any}
    before_hash::String
    after_hash::String
    inverse_data::Dict{String,Any}
    preserved_invariants::Dict{String,Any}
end

struct ConversionCertificate
    strategy::Symbol
    source_hash::String
    joining_disk::Dict{String,Any}
    isotopy_moves::Vector{Any}
    branch_convention::String
    label_map::Dict{String,Any}
    generated_hash::Union{Nothing,String}
    comparison_invariants::Dict{String,Any}
    evidence_status::Symbol
end

"""Embedded disk face bound to an exact labelled diagram, with disk-left boundary darts."""
struct JoiningDisk
    labelled_diagram_hash::String
    region::Int
    boundary_darts::Vector{Int}
end

struct PuncturedHeegaardDiagram
    diagram::HeegaardDiagram
    joining_disk::JoiningDisk
    cut::CutSurface
    certificate::ConversionCertificate
end

"""Directed singular edge after removing D and smoothing its corners."""
struct SpineEdge
    source_vertex::Int
    source_family::Symbol
    target_vertex::Int
    target_family::Symbol
    source_darts::Vector{Int}
end

"""Oriented disk sector, with signed singular-edge attaching word and source cells."""
struct SpineRegion
    source_cells::Vector{String}
    source_darts::Vector{Int}
    boundary::Vector{Int}
end

"""K08 B minus D, retaining the face-gluing witness and repeated incidences."""
struct BranchedSpineComplex
    regions::Vector{SpineRegion}
    gluings::Vector{Dict{String,Any}}
    edge_incidence::Vector{Vector{Tuple{Int,Int}}}
    checks::Dict{String,Any}
end

struct FlowSpine
    genus::Int
    edges::Vector{SpineEdge}
    vertex_signs::Dict{Int,Int}
    source::Any
    certificate::ConversionCertificate
    complex::BranchedSpineComplex
end

struct OData
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    metadata::Dict{String,Any}
end

struct ODataBundle
    punctured::PuncturedHeegaardDiagram
    spine::FlowSpine
    ograph::OData
    certificate::ConversionCertificate
    comparison::Dict{String,Any}
end

StructTypes.StructType(::Type{Curve}) = StructTypes.Struct()
StructTypes.StructType(::Type{Intersection}) = StructTypes.Struct()
StructTypes.StructType(::Type{HeegaardCode}) = StructTypes.Struct()
