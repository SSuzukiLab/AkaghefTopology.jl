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

struct PuncturedHeegaardDiagram
    diagram::HeegaardDiagram
    joining_disk::Dict{String,Any}
    cut::CutSurface
    certificate::ConversionCertificate
end

struct FlowSpine
    genus::Int
    branches::Vector{NamedTuple{(:curve_family, :curve_id, :intersection_ids),Tuple{Symbol,Int,Vector{Int}}}}
    source::Any
    certificate::ConversionCertificate
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
