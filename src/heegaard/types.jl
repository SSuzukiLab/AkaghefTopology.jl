"""Immutable-by-convention value types for the Heegaard subsystem.

Copying is constructor-specific: `Curve` copies its visits and `RegionPartition`
copies its blocks, while default record constructors can retain typed vectors
and dictionaries. Consult each constructor's contract for aliasing guarantees;
an immutable struct does not make its container fields deeply immutable.
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

"""
    Curve(family::Symbol, id::Integer, visits::AbstractVector{<:Integer}; label=nothing, tags=nothing) -> Curve

Store an oriented curve's cyclic intersection visits in a combinatorial Heegaard code.

# Arguments
`family` must be `:alpha` or `:beta` and `id` positive; violations throw
`ArgumentError` naming the family or identifier requirement. `visits` is
converted to `Vector{Int}` without incidence checks. `label` is converted
with `string`; its default is the family prefix followed by `id`.
`tags` supplies metadata pairs or `nothing`. Unsupported conversions propagate
their Julia exceptions.

# Returns
`Curve` with `family`, `id`, display `label`, the copied `visits`, and
`tags::Dict{String,Any}`. Its use in a diagram is described by
[legacy curve records](../specs/heegaard-terminology.md#legacy-cellular-records).

# Guarantees
1. The result retains the supplied visit order and has a positive identifier
   in the requested family. Mutating the input visit vector does not change it. Already-typed arguments are not modified. [proved/agent: by construction]
2. The record carries curve data; it does not select a complement decoration
   or a placement in the reference surface. [proved/agent: by construction]

# Not guaranteed
Visit identifiers need not be positive, distinct or present in a diagram.
Empty visits are accepted. Construction does not establish an embedded curve,
a complete cut system or a Heegaard diagram. Arbitrary metadata objects are
not guaranteed to be deeply copied.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §4 — cyclic intersection words and their distinction from geometric placement.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([H.Curve(:alpha, 1, [1])], [H.Curve(:beta, 1, [1])], [H.Intersection(1,1,1,1)]);

julia> d = H.certify_heegaard(H.decorate(c, [[1]]));

julia> (H.is_valid(d), H.surface_genus(d))
(true, 1)
```

See also [`heegaard_code`](@ref), [`Intersection`](@ref).

Role: construct
"""
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

"""
    Intersection(id::Integer, alpha::Integer, beta::Integer, sign::Integer; source_label=nothing) -> Intersection

Store one signed alpha–beta intersection in a combinatorial Heegaard code.

# Arguments
`id`, `alpha` and `beta` must be positive, and `sign` must be `-1` or `1`.
These checks throw `ArgumentError` naming the offending identifier or sign.
Integers are converted to `Int`; conversion failures propagate.
`source_label` is `nothing` or a value accepted by `string`.

# Returns
`Intersection` with its identifier, alpha and beta curve identifiers, sign,
and optional source label. The sign belongs to
[the oriented intersection convention](../specs/heegaard-code.md#intersection-sign-translation).

# Guarantees
1. Returned identifiers are positive and the sign is a unit integer.
   Arguments are not modified. [proved/agent: by construction]
2. The record specifies an incidence and its sign, not a position on a
   realized surface. Its cyclic location is supplied separately by the
   [curve words](../specs/heegaard-terminology.md#4-交点語符号).
   [proved/agent: by construction]

# Not guaranteed
The referenced curves need not exist, and the identifier need not be unique
in a larger code. Construction does not check realizability or consistency
with the curves' visits.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §§4, 9.1 — intersection records and the distinction from cyclic order and algebraic intersection.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); d = H.certify_heegaard(H.decorate(c, [[1]]));

julia> (H.intersection_matrix(d)[1,1], H.first_homology(d).free_rank)
(1, 0)
```

See also [`Curve`](@ref), [`heegaard_code`](@ref), [`intersection_matrix`](@ref).

Role: construct
"""
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

"""
    CombinatorialHeegaardCode(alpha, beta, intersections, metadata) -> CombinatorialHeegaardCode

Store the signed curve-incidence data of an undecorated Heegaard code.

# Arguments
`alpha::Vector{Curve}`, `beta::Vector{Curve}`,
`intersections::Vector{Intersection}`, `metadata::Dict{String,Any}`, or
values convertible to these types. Conversion errors propagate. The
constructor does not check curve families, incidence or surface validity.

# Returns
`CombinatorialHeegaardCode` with the supplied two curve families,
intersection records and metadata. The scope of this legacy code is given
in [the module structure](../specs/module-structure.md#1-構成).

# Guarantees
1. Fields preserve the supplied values after conversion; typed vectors and
   dictionaries are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. The value stores combinatorial incidence data without a complement-region
   decoration. [`decorate`](@ref) supplies the additional region data needed
   for the decorated-diagram validation path. [proved/agent: by construction]

# Not guaranteed
Direct construction does not validate even the syntax of the supplied
records. The code alone does not determine the missing complement decoration
or certify a closed Heegaard diagram. It is not a resolved Lickorish-based HC.

# Sources
cf: lib/topology.jl/docs/specs/module-structure.md §2 — legacy combinatorial code as distinct from resolved Lickorish-based HC.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]);

julia> d = H.CombinatorialHeegaardCode(c.alpha, c.beta, c.intersections, c.metadata);

julia> (H.validate(d).valid, any(i -> i.code == :missing_decoration, H.validate(d).issues))
(false, true)
```

See also [`heegaard_code`](@ref), [`decorate`](@ref), [`validate`](@ref).

Role: convert
"""
struct CombinatorialHeegaardCode
    alpha::Vector{Curve}
    beta::Vector{Curve}
    intersections::Vector{Intersection}
    metadata::Dict{String,Any}
end

"""
    RegionPartition(blocks::AbstractVector) -> RegionPartition

Store the grouping of ribbon-boundary cycles into complement regions.

# Arguments
`blocks` is a vector whose entries can each be collected and converted
elementwise to `Int`. Iteration or conversion errors propagate. No check
is made that the entries form a partition of any diagram's boundary cycles.

# Returns
`RegionPartition` with `boundary_cycle_blocks::Vector{Vector{Int}}`.
The meaning of complement regions belongs to
[legacy region records](../specs/heegaard-terminology.md#legacy-cellular-records).

# Guarantees
1. The order within and between supplied blocks is retained after conversion.
   The vectors are copied and the inputs are not modified.
   [proved/agent: by construction]
2. The value records a proposed complement decoration for the combinatorial
   diagram path; it does not compute a complement from curve words.
   [proved/agent: by construction]

# Not guaranteed
Empty blocks, repeated entries and indices outside the eventual diagram's
range are not rejected here. Region topology and compatibility are not
certified; use [`validate`](@ref) on the decorated code.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §5 — complement regions as distinct from individual boundary circuits and reference-subdivided faces.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]);

julia> p = H.RegionPartition([[1]]); d = H.certify_heegaard(H.DecoratedHeegaardCode(c,p,Dict{String,Any}()));

julia> H.complement_regions(d)[1].boundary_cycles == [1]
true
```

See also [`decorate`](@ref), [`validate`](@ref), [`complement_regions`](@ref).

Role: construct
"""
struct RegionPartition
    boundary_cycle_blocks::Vector{Vector{Int}}
    function RegionPartition(blocks::AbstractVector)
        new([Int.(collect(block)) for block in blocks])
    end
end

"""
    DecoratedHeegaardCode(code, regions, metadata) -> DecoratedHeegaardCode

Store a combinatorial Heegaard code with its complement-region decoration.

# Arguments
`code::CombinatorialHeegaardCode`, `regions::RegionPartition`,
`metadata::Dict{String,Any}`, or values convertible to those field types.
Conversion errors propagate. Region and incidence compatibility is not checked.

# Returns
`DecoratedHeegaardCode` with the supplied `code`, region partition and
metadata. Complement regions are distinguished from subdivided faces in
[the terminology specification](../specs/heegaard-terminology.md#5-面細分した領域補領域).

# Guarantees
1. Fields preserve the supplied values after conversion. The supplied code,
   partition and typed dictionary are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. The value carries both incidence and region data for the cellular-diagram
   validation path. [`certify_heegaard`](@ref) checks the data before returning
   a `CellularHeegaardDiagram`. [proved/agent: by construction]

# Not guaranteed
The constructor does not ensure that the blocks partition the boundary
cycles, that the resulting surface is valid or that either curve family is
a complete cut system. A decoration alone is not a validation certificate.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §5 and lib/topology.jl/docs/specs/heegaard-diagram.md — complement-region data and the cellular representation boundary.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]);

julia> d = H.DecoratedHeegaardCode(c, H.RegionPartition([[1]]), Dict{String,Any}());

julia> H.validate(d).valid
true
```

See also [`decorate`](@ref), [`certify_heegaard`](@ref).

Role: convert
"""
struct DecoratedHeegaardCode
    code::CombinatorialHeegaardCode
    regions::RegionPartition
    metadata::Dict{String,Any}
end

"""
    ValidationIssue(level, code, message, path) -> ValidationIssue

Store one diagnostic about a Heegaard-code validation check.

# Arguments
`level, code::Symbol` and `message, path::String`, or convertible values.
Conversion errors propagate. The constructor does not check the tags or
whether the path identifies an element of any code.

# Returns
`ValidationIssue` with the check's level, reason code, explanatory message
and location path. The validation scope is described in
[the Heegaard interface](../specs/heegaard-diagram.md).

# Guarantees
1. Fields retain their supplied values after conversion; arguments are not
   modified. [proved/agent: by construction]
2. The record is a diagnostic and performs no mathematical validation of a
   curve arrangement. [proved/agent: by construction]

# Not guaranteed
The record does not establish that the described issue occurred. Its presence
alone neither certifies nor refutes a geometric realization.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-diagram.md — validation at the public Heegaard boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> i = ValidationIssue(:syntax, :invalid_visit, "positive id required", "curves.alpha");

julia> (i.level, i.code, i.path)
(:syntax, :invalid_visit, "curves.alpha")
```

See also [`ValidationReport`](@ref), [`validate`](@ref).

Role: construct
"""
struct ValidationIssue
    level::Symbol
    code::Symbol
    message::String
    path::String
end

"""
    ValidationReport(level, valid, issues, checks) -> ValidationReport

Store the declared outcome and diagnostics of a Heegaard-code validation.

# Arguments
`level::Symbol`, `valid::Bool`, `issues::Vector{ValidationIssue}` and
`checks::Dict{Symbol,Bool}`, or convertible values. Conversion errors
propagate. No consistency check relates the flag, issues and checks.

# Returns
`ValidationReport` with the requested validation level, overall validity
flag, issue records and named check results.

# Guarantees
1. Fields retain the supplied values after conversion; typed containers are
   shared without mutation. [proved/agent: by construction]
2. The record can carry an outcome from the
   [Heegaard validation interface](../specs/heegaard-diagram.md);
   its constructor performs no validation. [proved/agent: by construction]

# Not guaranteed
An arbitrary true `valid` field is not a certificate of a surface or
Heegaard diagram. Use `validate` to obtain a computed report and inspect its
requested scope; a successful partial check is not full certification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-diagram.md — validation and certification are distinct from record construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> r = ValidationReport(:syntax, true, ValidationIssue[], Dict(:syntax => true));

julia> (r.level, r.valid, isempty(r.issues), r.checks[:syntax])
(:syntax, true, true, true)
```

See also [`ValidationIssue`](@ref), [`ValidationError`](@ref), [`validate`](@ref).

Role: convert
"""
struct ValidationReport
    level::Symbol
    valid::Bool
    issues::Vector{ValidationIssue}
    checks::Dict{Symbol,Bool}
end

"""
    ValidationError(report) -> ValidationError

Store a Heegaard validation report as an exception value.

# Arguments
`report::ValidationReport`, or a value convertible to that type. Conversion
errors propagate. The constructor does not require `report.valid == false`.

# Returns
`ValidationError` with the report in `report`. `showerror` prints
"invalid Heegaard diagram" followed by each issue's level, code and message.

# Guarantees
1. The supplied report is retained without mutation or copying of its
   containers. [proved/agent: by construction]
2. Construction neither validates a diagram nor throws the resulting value;
   it packages the diagnostic scope of the
   [Heegaard interface](../specs/heegaard-diagram.md).
   [proved/agent: by construction]

# Not guaranteed
A directly constructed exception need not describe a failed validation.
The display does not include issue paths or the report's check dictionary.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-diagram.md — validation failure at the Heegaard boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> i = ValidationIssue(:syntax, :invalid_visit, "positive id required", "curves.alpha");

julia> e = ValidationError(ValidationReport(:syntax, false, [i], Dict(:syntax => false)));

julia> (e.report.valid, only(e.report.issues).code)
(false, :invalid_visit)
```

See also [`ValidationReport`](@ref), [`certify_heegaard`](@ref).

Role: convert
"""
struct ValidationError <: Exception
    report::ValidationReport
end

function Base.showerror(io::IO, err::ValidationError)
    print(io, "invalid Heegaard diagram")
    for issue in err.report.issues
        print(io, "; ", issue.level, "/", issue.code, ": ", issue.message)
    end
end

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/docs/DESIGN.md §判断 — セル分割で認証した図式は公開 HC ではない。公開名は Topology で HeegaardCode に束縛する
"""
    CellularHeegaardDiagram(code::DecoratedHeegaardCode, report::ValidationReport, ::Val{:certified}) -> CellularHeegaardDiagram

Store a decorated combinatorial diagram with a supplied successful validation report.

# Arguments
`code` and `report` should describe the same decorated diagram.
`report.valid` must be true; otherwise `ArgumentError` says that only a
valid report can certify a diagram. The last argument is `Val(:certified)`.
For validation of the actual code use [`certify_heegaard`](@ref).

# Returns
`CellularHeegaardDiagram` with fields `code` and `validation` retaining the
supplied objects. Its distinction from the resolved public diagram type is
specified in [the public boundary](../specs/heegaard-diagram.md).

# Guarantees
1. `validation.valid` is true when construction returns. Already-typed arguments are not modified. [proved/agent: by construction]
2. When `report` is a successful full validation of this exact `code`, the
   value carries the validated decorated diagram. Construction itself
   does not run that validation or modify either argument. [proved/agent: by construction]

# Not guaranteed
The constructor does not check the report's level or provenance. A stale or
manually fabricated report can be accepted. The type does not identify a
3-manifold by name and is not the resolved `Topology.HeegaardDiagram` type.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-diagram.md — cellular certification as a representation distinct from resolved HC.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> c = H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]);

julia> d = H.CellularHeegaardDiagram(c, H.validate(c; level=:full), Val(:certified));

julia> (H.is_valid(d), H.surface_genus(d), H.euler_characteristic(d))
(true, 1, 0)
```

See also [`certify_heegaard`](@ref), [`validate`](@ref).

Role: convert
"""
struct CellularHeegaardDiagram
    code::DecoratedHeegaardCode
    validation::ValidationReport
    function CellularHeegaardDiagram(code::DecoratedHeegaardCode, report::ValidationReport, ::Val{:certified})
        report.valid || throw(ArgumentError("only a valid report can certify a Heegaard diagram"))
        new(code, report)
    end
end

# Existing methods, including lescop_theta.jl, still annotate this name.
const HeegaardDiagram = CellularHeegaardDiagram

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

"""
    RibbonGraph(vertices, edges, darts, rotation, next_dart) -> RibbonGraph

Store the ribbon graph of a signed curve arrangement and its boundary traversal data.

# Arguments
`vertices::Vector{Int}`, `edges::Vector{RibbonEdge}`,
`darts::Vector{RibbonDart}`, `rotation::Dict{Int,Vector{Symbol}}`, and
`next_dart::Vector{Int}`, or convertible values. Conversion errors propagate.
Indices, incidences and the boundary permutation are not validated.

# Returns
`RibbonGraph` with vertex identifiers, edges, directed darts, local rotations
and boundary successors. `RibbonEdge` and `RibbonDart` are accessible in
`Topology.HeegaardCodes`. Its intended surface context is
[the cellular ribbon graph](../specs/heegaard-snappy.md#曲面の切開と残る接続).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed containers
   are retained and not modified. [proved/agent: by construction]
2. For output of [`ribbon_graph`](@ref), the data encode the signed curve
   arrangement's ribbon neighbourhood. Complement-region attachment is
   separate data and is not selected by this constructor.
   [proved/agent: by construction]

# Not guaranteed
Arbitrary fields need not define a graph or a permutation, and boundary
traversal can fail on invalid fields. The constructor does not certify a
closed surface, complete cut systems or a 3-manifold.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md §曲面の切開と残る接続 — signed cellular ribbon graph as input to surface cutting.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> g = H.ribbon_graph(d); h = H.RibbonGraph(g.vertices,g.edges,g.darts,g.rotation,g.next_dart);

julia> (length(h.darts) == 2length(h.edges), sort(vcat(H.boundary_cycles(h)...)) == collect(1:length(h.darts)))
(true, true)
```

See also [`ribbon_graph`](@ref), [`boundary_cycles`](@ref), [`RegionPartition`](@ref).

Role: construct
"""
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

"""
    GroupPresentation(generators, relators, beta_ids, metadata) -> GroupPresentation

Store a finite group presentation with relator identifiers and metadata.

# Arguments
`generators::Vector{String}`, `relators::Vector{Vector{Tuple{Int,Int}}}`,
`beta_ids::Vector{Int}`, `metadata::Dict{String,Any}`, or values convertible
to those types. Conversion errors propagate. Generator references and
relator identifiers are not validated.

# Returns
`GroupPresentation` with generator names, relator words as generator-index
and exponent pairs, beta identifiers associated with relators, and metadata.
Use the [presentation-word encoding](../specs/heegaard-terminology.md#presentation-word-encoding).

# Guarantees
1. Fields preserve the supplied values after conversion; typed containers
   are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. Well-formed fields encode the group of the supplied presentation.
   A presentation produced from a valid closed Heegaard diagram `H` with the
   required complete cut system is intended to present the fundamental
   group of `Mhat(H)` from
   [the Heegaard manifold specification](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
   [unverified: general manifold identification of the computed presentation]

# Not guaranteed
Direct construction does not simplify words, check index ranges or identify
the group with any manifold. Names and metadata are not a certificate of
such an identification. The general identification of the computed
presentation with the manifold group is not certified here. Abelianization
does not classify the group.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 and lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — Heegaard presentation, abelianization and the represented manifold Mhat(H).

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> p = H.fundamental_group_presentation(d); q = H.GroupPresentation(p.generators, p.relators, p.beta_ids, p.metadata);

julia> h = H.first_homology(q); (h.free_rank, Tuple(h.torsion))
(0, ())
```

See also [`fundamental_group_presentation`](@ref), [`Topology.HeegaardCodes.first_homology`](@ref).

Role: construct
"""
struct GroupPresentation
    generators::Vector{String}
    relators::Vector{Vector{Tuple{Int,Int}}}
    beta_ids::Vector{Int}
    metadata::Dict{String,Any}
end

"""
    Puncture(id, pair, side, center, radius, orientation) -> Puncture

Store a circular chart placement for one boundary of an alpha-cut surface.

# Arguments
`id::Int`, `pair::Int`, `side::Symbol`, `center::NTuple{2,Float64}`,
`radius::Float64`, `orientation::Int`, or convertible values. Conversion
errors propagate. No positivity, orientation, side or disjointness check
is performed.

# Returns
`Puncture` with boundary identifier, boundary-pair identifier, side tag,
planar center and radius, and orientation tag. Its chart role belongs to
[legacy cut and puncture conventions](../specs/heegaard-terminology.md#legacy-cut-indices).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. For output of [`default_puncture_placement`](@ref), the record supplies a
   drawing placement for a cut boundary, not a new hole in the original
   closed surface. [proved/agent: by construction]

# Not guaranteed
The record does not establish a valid cut system, nonintersecting placement,
or correspondence with any `CutSurface`. Arbitrary radii and coordinates
are accepted, including nonfinite values.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §7 — cut boundaries and their displayed holes are distinct objects.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> c = H.cut_along_alpha(d); ps = [H.Puncture(p.id,p.pair,p.side,p.center,p.radius,p.orientation) for p in H.default_puncture_placement(c)];

julia> (length(ps), length(H.place_ports(c,ps)))
(2, 2)
```

See also [`default_puncture_placement`](@ref), [`place_ports`](@ref).

Role: construct
"""
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

"""
    CutSurface(diagram_hash, genus, boundaries, alpha_orders, beta_orders, intersection_signs, intersection_betas) -> CutSurface

Store the paired-boundary and incidence data of an alpha-cut Heegaard surface.

# Arguments
`diagram_hash::String`, `genus::Int`, `boundaries::Vector{CutBoundary}`,
`alpha_orders::Dict{Int,Vector{Int}}`, `beta_orders::Dict{Int,Vector{Int}}`,
`intersection_signs::Dict{Int,Int}`, `intersection_betas::Dict{Int,Int}`,
or values convertible to those types. Conversion errors propagate. Obtain
computed records from [`cut_along_alpha`](@ref).

# Returns
`CutSurface` with the supplied diagram identity, genus, boundary records,
ordered alpha and beta visits, crossing signs and crossing-to-beta lookup.
`CutBoundary` is available in `Topology.HeegaardCodes`. Boundary IDs and
port correspondence follow [legacy cut conventions](../specs/heegaard-terminology.md#legacy-cut-indices).

# Guarantees
1. Fields equal the supplied values after conversion; typed containers are
   retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. For output of `cut_along_alpha`, the data describe the alpha-cut surface
   of that input diagram in the scope of
   [cutting and drawing](../specs/heegaard-terminology.md#7-切断と描画). [proved/agent: by construction]

# Not guaranteed
Direct construction does not verify a complete cut system, the diagram
identity, the number of boundaries or any incidence consistency. It does
not supply geometric coordinates or a new manifold certificate.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §§6–7 — complete cut systems, paired cut boundaries and their drawing data.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> c = H.cut_along_alpha(d);

julia> copied = H.CutSurface(c.diagram_hash, c.genus, c.boundaries, c.alpha_orders, c.beta_orders, c.intersection_signs, c.intersection_betas);

julia> (length(H.default_puncture_placement(copied)), length(H.place_ports(copied)))
(2, 2)
```

See also [`cut_along_alpha`](@ref), [`place_ports`](@ref).

Role: construct
"""
struct CutSurface
    diagram_hash::String
    genus::Int
    boundaries::Vector{CutBoundary}
    alpha_orders::Dict{Int,Vector{Int}}
    beta_orders::Dict{Int,Vector{Int}}
    intersection_signs::Dict{Int,Int}
    intersection_betas::Dict{Int,Int}
end

"""
    Port(id, intersection_id, boundary_id, beta_id, parameter, point) -> Port

Store one beta-curve endpoint on a displayed boundary of an alpha-cut surface.

# Arguments
`id`, `intersection_id`, `boundary_id`, `beta_id` are `Int`;
`parameter::Float64` and `point::NTuple{2,Float64}`, or convertible values.
Conversion errors propagate. Identifiers, parameter range and coordinate
consistency are not checked.

# Returns
`Port` with port, intersection, cut-boundary and beta-curve identifiers,
the boundary parameter and planar coordinate. Its role is specified in
[legacy cut and port conventions](../specs/heegaard-terminology.md#legacy-cut-indices).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are not modified. [proved/agent: by construction]
2. For records from [`place_ports`](@ref), the endpoint refers to a beta
   half-edge at an alpha cut of the supplied diagram, with a position in
   the supplied puncture placement. [proved/agent: by construction]

# Not guaranteed
Construction does not verify any incidence or placement, and an arbitrary
point need not lie on the indicated boundary. Port identifiers alone do
not specify the beta arcs joining them.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §7 — ports and wiring on the alpha-cut surface.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> c = H.cut_along_alpha(d); ps = [H.Port(p.id,p.intersection_id,p.boundary_id,p.beta_id,p.parameter,p.point) for p in H.place_ports(c)];

julia> H.build_wiring_problem(c,ps).wires == [(beta_id=1, from_port=1, to_port=2)]
true
```

See also [`place_ports`](@ref), [`build_wiring_problem`](@ref).

Role: construct
"""
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

"""
    MoveCertificate(move_type, attachment_data, label_correspondence, before_hash, after_hash, inverse_data, preserved_invariants) -> MoveCertificate

Store a combinatorial Heegaard move's attachment, correspondence and inverse data.

# Arguments
`move_type::Symbol`; `attachment_data`, `label_correspondence`, `inverse_data`
and `preserved_invariants` are `Dict{String,Any}`; `before_hash` and
`after_hash` are strings. Convertible values are accepted; conversion errors
propagate. No hash, correspondence or preservation claim is checked.

# Returns
`MoveCertificate` with the move tag, attachment choices, label correspondence,
input and output identity strings, inverse-operation data and recorded
invariant data. Use the distinction between a move and its witness in
[the move-witness boundary](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed dictionaries
   are retained and not modified. [proved/agent: by construction]
2. The record can accompany a move on a decorated combinatorial diagram;
   construction performs no move and computes no invariant.
   [proved/agent: by construction]

# Not guaranteed
The type name does not certify geometric preservation or invertibility.
The stored invariant data need not have been independently compared on
both endpoints. Manually supplied hashes and inverse data may be stale.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §§9, 9.2 — moves, correspondence and the separate scope of their evidence.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> _, c = H.reverse_curve_orientation(d, (:alpha,1)); e = H.MoveCertificate(c.move_type,c.attachment_data,c.label_correspondence,c.before_hash,c.after_hash,c.inverse_data,c.preserved_invariants);

julia> H.inverse_move(e)["move"]
"reverse_curve_orientation"
```

See also [`apply_move`](@ref), [`inverse_move`](@ref).

Role: construct
"""
struct MoveCertificate
    move_type::Symbol
    attachment_data::Dict{String,Any}
    label_correspondence::Dict{String,Any}
    before_hash::String
    after_hash::String
    inverse_data::Dict{String,Any}
    preserved_invariants::Dict{String,Any}
end

"""
    ConversionCertificate(strategy, source_hash, joining_disk, isotopy_moves, branch_convention, label_map, generated_hash, comparison_invariants, evidence_status) -> ConversionCertificate

Store correspondence and evidence data for a Heegaard-to-spine conversion.

# Arguments
`strategy::Symbol`, `source_hash::String`, `joining_disk::Dict{String,Any}`,
`isotopy_moves::Vector{Any}`, `branch_convention::String`,
`label_map::Dict{String,Any}`, `generated_hash::Union{Nothing,String}`,
`comparison_invariants::Dict{String,Any}`, `evidence_status::Symbol`, or
values convertible to those field types. Conversion errors propagate.

# Returns
`ConversionCertificate` with the supplied strategy, source identity, joining
disk data, move records, convention name, label correspondence, optional
output identity, comparison results and evidence-status tag.

# Guarantees
1. Fields equal the supplied values after conversion. Typed containers are
   retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. The record can carry the source-to-spine correspondence of
   [the conversion specification](../specs/heegaard-koda-conversion.md#64-返す対応).
   Its contents are supplied evidence, not an additional computation. [proved/agent: by construction]

# Not guaranteed
No hash, correspondence, joining disk or evidence tag is verified. Merely
constructing this record does not prove that a conversion preserves the
manifold, branching, combing or framing. Mutable dictionaries can invalidate
previously recorded evidence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.4, 10.2–10.3 — correspondence and the separate preservation claims recorded by a conversion.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "HeegaardAbalone.jl"));

julia> e = HeegaardAbaloneExample.fixture(); c = H.to_odata(e.diagram; joining_disk=e.disk).certificate;

julia> d = H.ConversionCertificate(c.strategy, c.source_hash, c.joining_disk, c.isotopy_moves, c.branch_convention, c.label_map, c.generated_hash, c.comparison_invariants, c.evidence_status);

julia> (Tuple(d.label_map["smoothed_intersections"]), Tuple(d.label_map["remaining_intersections"]))
((1, 2), (3,))
```

See also [`make_punctured`](@ref), [`to_odata`](@ref).

Role: construct
"""
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

"""
    JoiningDisk(labelled_diagram_hash, region, boundary_darts) -> JoiningDisk

Store a selected joining-region boundary bound to a labelled legacy diagram identity.

# Arguments
`labelled_diagram_hash::String`, `region::Int`, and
`boundary_darts::Vector{Int}`, or convertible values. Conversion errors
propagate. This constructor does not check the hash, region or boundary
against the [joining-disk domain](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).

# Returns
A `JoiningDisk` with the supplied identity string, region index and boundary
dart sequence in [legacy coordinates](../specs/heegaard-terminology.md#legacy-cellular-records).

# Guarantees
1. Fields preserve supplied values after conversion; a typed boundary
   vector is retained without mutation. [proved/agent: by construction]
2. The record can identify the selected disk of the punctured construction;
   it does not establish that the supplied region is such a disk.
   [proved/agent: by construction]

# Not guaranteed
A raw value need not be embedded, joining, or attached to the intended
diagram. Use `joining_disk` or `validate_joining_disk` for the checked path.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.2 — joining disk for a punctured Heegaard diagram.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = JoiningDisk("example", 1, [1,8]); (d.region, Tuple(d.boundary_darts))
(1, (1, 8))
```

See also [`joining_disk`](@ref), [`validate_joining_disk`](@ref).

Role: construct
"""
struct JoiningDisk
    labelled_diagram_hash::String
    region::Int
    boundary_darts::Vector{Int}
end

"""
    PuncturedHeegaardDiagram(diagram, joining_disk, cut, certificate) -> PuncturedHeegaardDiagram

Store a cellular Heegaard diagram with a selected joining disk and conversion data.

# Arguments
`diagram::CellularHeegaardDiagram`, `joining_disk::JoiningDisk`,
`cut::CutSurface`, `certificate::ConversionCertificate`, or convertible
values. Conversion errors propagate. Direct construction does not check
[the joining conditions](../specs/heegaard-koda-conversion.md#5-heegaard-条件と-joining-disk-の入力検査).
Use [`make_punctured`](@ref) for the checked path.

# Returns
`PuncturedHeegaardDiagram` with the retained diagram, joining disk, alpha-cut
data and conversion certificate. This record uses the cellular diagram
type, as distinguished by [the public boundary](../specs/heegaard-diagram.md).

# Guarantees
1. Already-typed arguments are retained without copying or modification.
   [proved/agent: by construction]
2. When produced by `make_punctured`, the selected disk is attached to that
   input diagram's data. The intended punctured diagram belongs to
   [the Koda conversion domain](../specs/heegaard-koda-conversion.md).
   [proved/agent: by construction]

# Not guaranteed
Direct construction does not validate the joining disk, source hashes,
cut data or certificate, and does not construct a spine. A stale disk or
certificate can be stored. The constructor makes no manifold-preservation claim.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§5–6 — joining data and the punctured-diagram domain of the spine construction.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "HeegaardAbalone.jl"));

julia> e = HeegaardAbaloneExample.fixture(); p = H.make_punctured(e.diagram; joining_disk=e.disk); q = H.PuncturedHeegaardDiagram(p.diagram,p.joining_disk,p.cut,p.certificate);

julia> b = H.to_odata(q); (Tuple(only(b.ograph.gauss)), Tuple(b.ograph.orientation))
((1, -1), (1,))
```

See also [`make_punctured`](@ref), [`validate_joining_disk`](@ref), [`to_odata`](@ref).

Role: convert
"""
struct PuncturedHeegaardDiagram
    diagram::HeegaardDiagram
    joining_disk::JoiningDisk
    cut::CutSurface
    certificate::ConversionCertificate
end

"""
    SpineEdge(source_vertex, source_family, target_vertex, target_family, source_darts) -> SpineEdge

Store a directed singular strand of the Koda spine with its source darts.

# Arguments
Integer source and target vertex IDs, symbolic source and target curve
families, and `source_darts::Vector{Int}`, or convertible values. Conversion
errors propagate. The caller supplies the valid strand data of
[Koda conversion §6.3](../specs/heegaard-koda-conversion.md#63-特異辺と-sector);
the constructor does not check IDs, families or incidence.

# Returns
A `SpineEdge` retaining its two endpoint/family pairs and source-dart list.

# Guarantees
1. Fields preserve supplied values after conversion; a typed dart vector
   is shared and not modified. [proved/agent: by construction]
2. The record can carry a directed singular strand of the specified spine
   construction, without constructing its ambient thickening.
   [proved/agent: by construction]

# Not guaranteed
Arbitrary edge records need not form closed strand cycles or a branched spine.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.3, 7.1 — singular-edge and encoding correspondence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> e = SpineEdge(1, :alpha, 1, :beta, Int[]); (e.source_vertex, e.target_family)
(1, :beta)
```

Role: construct
"""
struct SpineEdge
    source_vertex::Int
    source_family::Symbol
    target_vertex::Int
    target_family::Symbol
    source_darts::Vector{Int}
end

"""
    SpineRegion(source_cells, source_darts, boundary) -> SpineRegion

Store a spine sector's source cells and oriented attaching data.

# Arguments
`source_cells::Vector{String}`, `source_darts::Vector{Int}` and
`boundary::Vector{Int}`, or convertible values. Conversion errors propagate.
The caller supplies consistent [sector and cellular coordinates](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates);
no disk, boundary or incidence check is performed by this constructor.

# Returns
A `SpineRegion` retaining source-cell labels, source darts and the signed
edge-attaching word in the linked convention.

# Guarantees
1. Fields preserve supplied values after conversion; typed vectors are
   retained without mutation. [proved/agent: by construction]
2. Repeated incidences are retained. The record can carry a sector of the
   Koda complex without constructing it geometrically.
   [proved/agent: by construction]

# Not guaranteed
A raw attaching word does not certify a disk sector or a valid branched complex.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.3–6.5 — sectors, repeated incidences and orientations.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> r = SpineRegion(["cell"], Int[], [1,-1,1]); (Tuple(r.source_cells), Tuple(r.boundary))
(("cell",), (1, -1, 1))
```

Role: construct
"""
struct SpineRegion
    source_cells::Vector{String}
    source_darts::Vector{Int}
    boundary::Vector{Int}
end

"""
    BranchedSpineComplex(regions, gluings, edge_incidence, checks) -> BranchedSpineComplex

Store the regions, gluing correspondence and incidence data of a branched-spine complex.

# Arguments
`regions::Vector{SpineRegion}`, `gluings::Vector{Dict{String,Any}}`,
`edge_incidence::Vector{Vector{Tuple{Int,Int}}}` and `checks::Dict{String,Any}`,
or convertible values. Conversion errors propagate. The caller supplies the
consistent complex of [Koda conversion §6](../specs/heegaard-koda-conversion.md#6-heegaard-diagramから-branched-spine-を作る);
the constructor does not validate it.

# Returns
A `BranchedSpineComplex` retaining sector records, source-cell gluing records,
signed region incidences at each edge and diagnostic checks.

# Guarantees
1. The fields preserve supplied values after conversion, including repeated
   incidences. Already-typed containers are retained without mutation.
   [proved/agent: by construction]
2. The record can carry the finite complex and correspondence of the
   specified punctured-diagram construction. It does not itself construct
   a thickening. [proved/agent: by construction]

# Not guaranteed
Arbitrary field values need not form a branched spine or represent a
3-manifold. Diagnostic entries are supplied assertions, not fresh checks.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.3–6.5 — sectors, gluing correspondence and cellular coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> c = BranchedSpineComplex(SpineRegion[], Dict{String,Any}[], Vector{Tuple{Int,Int}}[], Dict{String,Any}()); (length(c.regions), length(c.edge_incidence))
(0, 0)
```

See also [`heegaard_to_spine`](@ref).

Role: convert
"""
struct BranchedSpineComplex
    regions::Vector{SpineRegion}
    gluings::Vector{Dict{String,Any}}
    edge_incidence::Vector{Vector{Tuple{Int,Int}}}
    checks::Dict{String,Any}
end

"""
    FlowSpine(genus, edges, vertex_signs, source, certificate, complex) -> FlowSpine

Store a singular-graph presentation and its branched-spine complex and correspondence data.

# Arguments
`genus::Int`, `edges::Vector{SpineEdge}`, `vertex_signs::Dict{Int,Int}`,
arbitrary `source`, `certificate::ConversionCertificate`, and
`complex::BranchedSpineComplex`, or values convertible to those types.
Conversion errors propagate. No geometric or incidence validation is made.

# Returns
`FlowSpine` with the supplied source genus, directed singular edges, vertex
signs, retained source, conversion certificate and sector/gluing complex.
The represented object and its conversion domain are specified in
[the spine construction](../specs/heegaard-koda-conversion.md#6-heegaard-diagramから-branched-spine-を作る).

# Guarantees
1. Fields preserve the supplied values after conversion; typed containers
   and supplied objects are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. When supplied by `to_flow_spine` on its validated punctured-diagram
   domain, the record describes that diagram's branched spine, with the
   correspondence and conventions of the linked construction.
   [unverified: general correctness of the implemented source-to-spine correspondence]

# Conventions
Use [Heegaard and Koda orientation conventions](../conventions/HeegaardCodes.md).

# Not guaranteed
The constructor does not certify a spine, a flow, a closed 3-manifold or a
framing. `source` can be arbitrary. The presence of a certificate does not
check its content or prove that it corresponds to `complex` and `edges`.
The general correctness of the computed source-to-spine correspondence is
not certified by this storage constructor.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.1–6.4 — punctured-diagram spine construction and the scope of Koda Lemma 5.3.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "HeegaardAbalone.jl"));

julia> e = HeegaardAbaloneExample.fixture(); s = H.to_odata(e.diagram; joining_disk=e.disk).spine;

julia> t = H.FlowSpine(s.genus, s.edges, s.vertex_signs, s.source, s.certificate, s.complex);

julia> Tuple(Tuple(edge.source_darts) for edge in t.edges)
((5, 12), (10, 3))
```

See also [`to_flow_spine`](@ref), [`spine_to_odata`](@ref).

Role: convert
"""
struct FlowSpine
    genus::Int
    edges::Vector{SpineEdge}
    vertex_signs::Dict{Int,Int}
    source::Any
    certificate::ConversionCertificate
    complex::BranchedSpineComplex
end

"""
    OData(gauss, orientation, metadata) -> OData

Store o-graph component words, vertex signs and metadata for the Koda conversion boundary.

# Arguments
`gauss::Vector{Vector{Int}}`, `orientation::Vector{Int}`,
`metadata::Dict{String,Any}`, or convertible values. Conversion errors
propagate. No word, sign or closed-normal condition is checked.

# Returns
`OData` with component words, vertex-sign data and metadata. The meaning of
o-data is owned by [the o-graph terminology](../specs/ograph-terminology.md)
and [the conversion boundary](../specs/heegaard-diagram.md).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed containers
   are retained and not modified. [proved/agent: by construction]
2. The result stores an o-graph presentation, not Heegaard curve words
   interpreted as such. The constructor performs no conversion between the
   two representations. [proved/agent: by construction]

# Not guaranteed
Arbitrary stored data need not be a valid o-graph or represent a closed
3-manifold. Construction does not establish a spine, branching, combing,
framing or equivalence to an input Heegaard diagram. Metadata is unchecked.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-diagram.md — o-data as the VirtualLink representation at the Koda conversion boundary.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; g = H.OData([[1,-1]], [1], Dict{String,Any}());

julia> s = H.odata_to_spine(g); (Tuple(only(s.gauss)), s.vertex_signs[1])
((1, -1), 1)
```

See also [`odata_to_spine`](@ref), [`compare_odata_words`](@ref), [`spine_to_odata`](@ref).

Role: construct
"""
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
StructTypes.StructType(::Type{CombinatorialHeegaardCode}) = StructTypes.Struct()
