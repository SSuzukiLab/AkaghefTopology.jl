# cf:[[BP97#^bp97-sec-0549]] Definition 2.1.4 and Fig. 2.5 — O0–O3 and the edge chart (1 2 3)^α; adjustment: upright legs (upper-right, up, lower-left), upside-down legs (lower-left, down, upper-right); the colour is the undirected Z_3 label and is not negated when the drawing is reversed; unverified

"""Coloured o-graphs in the sense of Definition 2.1.4.

The general object is a half-edge multigraph. Loops and parallel edges stay.
It has no vertex sign, no edge direction, and no Gauss word.

Each vertex stores its incident half-edges in cyclic order. On a quadrivalent
vertex the opposite positions are `(1, 3)` and `(2, 4)`. `over_positions[v]`
lists the local positions marked over. Exactly one of those two pairs must be
marked.

Each edge carries one colour, an element of Z_3, stored on either of its
half-edges. Integers are read modulo 3, with the representative in
`{-1, 0, 1}`: `2`, `-2`, and `4` are colours, not malformed data. The printed
Fig. 2.5 rule, once an end has been drawn on the right, glues right leg `i` to
left leg `τ^α(i)`. Reversing the drawing reverses the transverse frame by
`ρ = (1 3)` and keeps `α`. Freezing the leg numbers and only reversing the
bijection produces `-α`; that is a different calculation, and it is not how a
colour is stored. Status of this transport: draft.

`bp97-fig2.5-rl-tau/draft1` names the chart. Figure 2.5 (left) specifies the
local polyhedron. The over pair is the upright T and the under pair is the
upside-down T, and leg 2 meets the opposite half-edge. Figure 2.7 fixes the
local embedding once disc D has the usual orientation of the horizontal plane
and the screw is left-handed. Which side sheet of that embedding is leg 1 is
not registered onto the counter-clockwise slot. Candidate B is the calibration
derived from Proposition 3.2.1 together with the §2.1 reconstruction and the
§3.2 construction, under the draft dependencies on the convention page.
Candidate A on the closed presentation `s2xs1_2vertex` gives five discs and
Euler characteristic 3. `face_cycles` builds the discs. The default
calibration is `:b`. `OGRAPH_FACE_CYCLE_STATUS` records that derivation.
The status stays draft.

The reading check is `threads/PJ261120/sources/pro/core-20261121-235637.md`,
harvested at `threads/PJ261120/sources/pro/core-20261121-235637-answer.md`.
The colour sits on the undirected edge, as Definition 2.1.4 states. The frame
reversal `ρ = (1 3)` and the frozen-label inverse `-α` are inferences from
Figure 2.5. The twelve pairings are derived from the four panels of Figure 2.5
(right) and the gluing rule `τ^α`. The PDF does not print those twelve rows.
The condition-id check is `threads/PJ261120/sources/pro/core-20261122-012651.md`.
The slot-transcription check is `threads/PJ261120/sources/pro/core-20261122-033946.md`.
It leaves both charts of leg 1 open. The later calibration check is
`threads/PJ261120/sources/pro/core-20261122-045101.md`. Candidate B is the
calibration derived from Proposition 3.2.1, and the status stays draft.
"""

using SHA

const OGRAPH_COLOUR_CONVENTION_ID = "bp97-fig2.5-rl-tau/draft1"
const OGRAPH_FACE_CYCLE_STATUS = "prop321-calibration-b/draft"
const _COLOUR_CANON_MAX_VERTICES = 6

"""
    DiscGerm(halfedge::Integer, leg::Integer) -> DiscGerm

Construct a numbered disc-germ label at a half-edge.

# Arguments
- `halfedge`: a positive integer representable as `Int`; a nonpositive value
  throws `ArgumentError` identifying the half-edge, and an out-of-range
  conversion can throw `InexactError`.
- `leg`: one of `1:3`; otherwise `ArgumentError` identifies the leg range.
  The numbered leg is interpreted by the
  [edge chart](../specs/ograph-colours.md#chart).

# Returns
`DiscGerm` with `halfedge::Int` and `leg::Int` equal to the supplied indices.

# Guarantees
1. `halfedge > 0` and `leg in 1:3` for every returned value. [proved/agent: by construction]
2. The value labels a local disc germ in the Figure 2.5 chart when its
   half-edge belongs to the intended coloured graph. It does not select
   a global disc without that graph and its gluing data. [proved/agent: by construction]

# Conventions
Use the [numbered-leg convention](../conventions/OGraphColours.md).

# Not guaranteed
The constructor has no graph argument and does not check that the half-edge
exists. It does not settle the leg-to-slot calibration.

# Examples
```jldoctest
julia> germs = edge_disc_germs(1, 4, (2, 1, 4, 3));

julia> germs[1] == (DiscGerm(4, 1), DiscGerm(3, 2))
true
```

See also [`edge_disc_germs`](@ref), [`face_cycles`](@ref).

Role: construct
"""
struct DiscGerm
    halfedge::Int
    leg::Int

    function DiscGerm(halfedge::Integer, leg::Integer)
        halfedge > 0 || throw(ArgumentError("disc-germ half-edge must be positive"))
        leg in 1:3 || throw(ArgumentError("disc-germ leg must be 1, 2, or 3"))
        new(Int(halfedge), Int(leg))
    end
end

"""
    EdgeColour(halfedge, colour) -> EdgeColour

Store an integer colour representative at a half-edge.

# Arguments
- `halfedge::Int` and `colour::Int`, or values convertible to `Int`.
  Failed conversions propagate their Julia exceptions, including
  `InexactError` for an unrepresentable integer.
- Edge membership and colour consistency are checked by O3 in
  [`validate_coloured_ograph`](@ref), not by this constructor.

# Returns
`EdgeColour` with `halfedge` the supplied edge-end index and `colour` the
supplied integer representative. The mathematical colour belongs to the
undirected edge under [O3](../specs/ograph-colours.md#condition-identifiers).

# Guarantees
1. Fields equal the supplied integers after conversion to `Int`. [proved/agent: by construction]
2. In a graph satisfying O3, this record contributes its colour to the
   coloured o-graph of Definition 2.1.4 under the
   [stored-colour convention](../specs/ograph-colours.md#chart).
   [proved/agent: by construction]

# Not guaranteed
Construction does not check positivity or existence of the half-edge,
normalize the colour, or validate the containing graph.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> records = [EdgeColour(r.halfedge, r.colour) for r in g.edge_colours];

julia> Tuple(sort([r.colour for r in records]))
(-1, 1)
```

See also [`ColouredOGraph`](@ref), [`validate_coloured_ograph`](@ref).

Role: construct
"""
struct EdgeColour
    halfedge::Int
    colour::Int
end

"""
    OGraphColourCheck(id, passed, message) -> OGraphColourCheck

Store the reported outcome of one coloured-o-graph condition.

# Arguments
`id::String`, `passed::Bool`, and `message::String`, or values convertible
to those types. Failed conversions propagate their Julia exceptions.
No consistency check is made between these arguments.

# Returns
`OGraphColourCheck` with `id` the condition label, `passed` its supplied
truth value and `message` the supplied explanation. Condition meanings
are owned by the [condition table](../specs/ograph-colours.md#condition-identifiers).

# Guarantees
1. Fields equal the supplied values after conversion. [proved/agent: by construction]
2. A record returned by `validate_coloured_ograph` reports one of O0–O3
   for the input graph; one returned by `validate_normal_colouring` reports
   one of N1, N2, C1–C3, P1 or P2 within that validator's stated scope.
   [proved/agent: by construction]

# Not guaranteed
An arbitrary `id` is accepted. Constructing a record with `passed=true`
does not prove a graph condition or certify a spine or manifold.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> c = first(validate_coloured_ograph(g)); d = OGraphColourCheck(c.id, c.passed, c.message);

julia> (d.id, d.passed, occursin(d.id, d.message))
("O0", true, true)
```

See also [`validate_coloured_ograph`](@ref), [`validate_normal_colouring`](@ref).

Role: construct
"""
struct OGraphColourCheck
    id::String
    passed::Bool
    message::String
end

"""
    ColouredOGraph(vertex_halfedges, mate, over_positions, edge_colours; convention_id=OGRAPH_COLOUR_CONVENTION_ID) -> ColouredOGraph

Store a half-edge presentation of a coloured o-graph for subsequent validation.

# Arguments
- `vertex_halfedges`, `mate`, `over_positions`: iterable integer data for the
  fields below; entries are converted to `Int`. There must be one
  `over_positions` entry per vertex, otherwise `ArgumentError` identifies
  that mismatch. Unrepresentable integer conversions can throw `InexactError`.
- `edge_colours`: iterable of `EdgeColour` values or two-integer tuples.
  Other record shapes throw `ArgumentError` describing the accepted forms.
- `convention_id::AbstractString`: must be nonempty after stripping surrounding
  whitespace; otherwise `ArgumentError` identifies the empty identifier.
  An unfamiliar nonempty identifier is stored without checking its meaning.

# Returns
`ColouredOGraph` with `vertex_halfedges` the tuple of cyclic half-edge tuples,
`mate` the partner tuple, `over_positions` the tuple of marked-position tuples,
`edge_colours` a tuple of `EdgeColour` records, and `convention_id` the stripped
string. Their mathematical interpretation is owned by the
[coloured o-graph specification](../specs/ograph-colours.md#condition-identifiers).

# Guarantees
1. The input sequences are stored as tuples with the declared integer or
   record element types; the input containers are not mutated. [proved/agent: by construction]
2. Integer colours are retained as supplied. Validation interprets them using
   O3; construction does not normalize them. [proved/agent: by construction]
3. If O0–O3 pass, the data constitute a coloured o-graph of Definition 2.1.4.
   The associated oriented standard spine uses the
   [specified reconstruction](../specs/ograph-colours.md#face-cycles).
   [proved/agent: by construction]

# Conventions
The edge chart and colour interpretation are in
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
Construction does not enforce O0–O3, closedness, or the existence of a
compatible normal orientation. It does not remove the reconstruction
dependencies recorded on the convention page.

# Sources
cf:[[BP97#^bp97-sec-0551]] Proposition 2.1.6 and Remark 2.1.7 — correspondence with oriented standard spines and validity of the reconstruction for one vertex.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> h = ColouredOGraph(g.vertex_halfedges, g.mate, g.over_positions, g.edge_colours; convention_id=g.convention_id);

julia> all(check.passed for check in validate_coloured_ograph(h))
true
```

See also [`validate_coloured_ograph`](@ref), [`canonical_coloured_ograph`](@ref).

Role: convert
"""
struct ColouredOGraph
    vertex_halfedges::Tuple{Vararg{Tuple{Vararg{Int}}}}
    mate::Tuple{Vararg{Int}}
    over_positions::Tuple{Vararg{Tuple{Vararg{Int}}}}
    edge_colours::Tuple{Vararg{EdgeColour}}
    convention_id::String

    function ColouredOGraph(vertex_halfedges, mate, over_positions, edge_colours;
                            convention_id::AbstractString=OGRAPH_COLOUR_CONVENTION_ID)
        convention = strip(String(convention_id))
        isempty(convention) && throw(ArgumentError("convention_id must not be empty"))
        vertices = Tuple(Tuple(Int.(collect(cyclic))) for cyclic in vertex_halfedges)
        partners = Tuple(Int.(collect(mate)))
        overs = Tuple(Tuple(Int.(collect(positions))) for positions in over_positions)
        stored = Tuple(_edge_colour(record) for record in edge_colours)
        length(overs) == length(vertices) || throw(ArgumentError(
            "over_positions must contain one entry per vertex"))
        new(vertices, partners, overs, stored, convention)
    end
end

_edge_colour(record::EdgeColour) = record

function _edge_colour(record::Tuple{Integer,Integer})
    EdgeColour(Int(record[1]), Int(record[2]))
end

function _edge_colour(record)
    throw(ArgumentError("an edge colour must be an EdgeColour or a (halfedge, colour) pair"))
end

function _colour_mod3(value::Integer)
    reduced = mod(Int(value), 3)
    reduced == 2 ? -1 : reduced
end

"""
    fig25_tau(alpha::Integer, leg::Integer) -> Int

Return the image of a numbered leg under the Figure 2.5 cyclic permutation.

# Arguments
- `alpha`: an integer exponent convertible to `Int`; conversion can throw
  `InexactError`. It is interpreted modulo three.
- `leg`: a leg in the [edge chart](../specs/ograph-colours.md#chart).
  Checked: an index outside `1:3` throws `ArgumentError` naming the leg range.

# Returns
The resulting leg number as an `Int`.

# Guarantees
1. The result is the image under the chart's permutation with exponent `alpha`. [proved/agent: by construction]
2. Congruent exponents give identical images. [proved/agent: by construction]
3. Arguments are not modified. [proved/agent: by construction]

# Conventions
Use the [Figure 2.5 leg numbering](../conventions/OGraphColours.md).

# Examples
```jldoctest
julia> Tuple(fig25_tau(1, leg) for leg in 1:3)
(2, 3, 1)
```

Role: convert
"""
function fig25_tau(alpha::Integer, leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    shift = mod(Int(alpha), 3)
    shift == 0 && return Int(leg)
    shift == 1 && return leg == 3 ? 1 : Int(leg) + 1
    return leg == 1 ? 3 : Int(leg) - 1
end

"""
    fig25_leg_direction(upright::Bool, leg::Integer) -> Tuple{Int,Int}

Return the planar direction of one numbered leg in the Figure 2.5 drawing.

# Arguments
- `upright`: selects the upright (`true`) or upside-down (`false`) chart.
- `leg`: a leg in the [edge chart](../specs/ograph-colours.md#chart).
  Checked: an index outside `1:3` throws `ArgumentError` naming the leg range.

# Returns
The pair of integer direction coordinates assigned to that leg by the chart.

# Guarantees
1. The coordinates agree with the selected chart's leg direction. [proved/agent: by construction]
2. Arguments are not modified. [proved/agent: by construction]

# Conventions
The [Figure 2.5 coordinate chart](../conventions/OGraphColours.md) fixes the
drawing directions and has draft geometric status.

# Not guaranteed
The pair is a direction label, not a unit vector or a three-dimensional embedding.

# Examples
```jldoctest
julia> fig25_leg_direction(false, 2)
(0, -1)
```

Role: access
"""
function fig25_leg_direction(upright::Bool, leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    upright && return ((1, 1), (0, 1), (-1, -1))[leg]
    return ((-1, -1), (0, -1), (1, 1))[leg]
end

"""
    fig25_leg_pairing(left_up::Bool, right_up::Bool, alpha::Integer) -> NTuple{3,Tuple{Tuple{Int,Int},Tuple{Int,Int}}}

Return the three pairs of planar leg directions in a coloured edge drawing.

# Arguments
- `left_up`, `right_up`: choose the two upright or upside-down charts.
- `alpha`: an integer colour in the [edge chart](../specs/ograph-colours.md#chart),
  convertible to `Int`; conversion can throw `InexactError`.

# Returns
Three `(right_direction, left_direction)` pairs in right leg order. Each
direction is the coordinate pair returned by [`fig25_leg_direction`](@ref).

# Guarantees
1. The pairs realize the linked chart's gluing for the supplied colour and drawing frames. [proved/agent: by construction]
2. Congruent integer colours give identical pairings. [proved/agent: by construction]
3. Arguments are not modified. [proved/agent: by construction]

# Conventions
Use the [Figure 2.5 pairing table](../conventions/OGraphColours.md), whose
geometric reading remains draft.

# Not guaranteed
These direction pairs alone do not reconstruct an entire polyhedron.

# Examples
```jldoctest
julia> fig25_leg_pairing(true, true, 0)
(((1, 1), (1, 1)), ((0, 1), (0, 1)), ((-1, -1), (-1, -1)))
```

Role: construct
"""
function fig25_leg_pairing(left_up::Bool, right_up::Bool, alpha::Integer)
    colour = _colour_mod3(alpha)
    ntuple(leg -> (fig25_leg_direction(right_up, leg),
                   fig25_leg_direction(left_up, fig25_tau(colour, leg))), 3)
end

"""
    edge_gluing_from_colour(alpha::Integer) -> NTuple{3,Int}

Return the leg permutation of a colour in the Figure 2.5 edge chart.

# Arguments
- `alpha`: a colour in the [edge chart](../specs/ograph-colours.md#chart).
  Conversion to `Int` is required and can throw `InexactError`; the colour
  is then read modulo three.

# Returns
`NTuple{3,Int}` giving the left leg paired with each numbered right leg.

# Guarantees
1. The tuple is the gluing assigned to `alpha` by the linked chart. [proved/agent: by construction]
2. Congruent integer colours give the same tuple. [proved/agent: by construction]
3. The argument is not modified. [proved/agent: by construction]

# Conventions
The [Figure 2.5 chart](../conventions/OGraphColours.md) fixes the leg numbering
and drawing direction; its geometric reading remains draft.

# Not guaranteed
This permutation alone does not specify an entire coloured o-graph or certify
its reconstruction.

# Examples
```jldoctest
julia> edge_gluing_from_colour(1)
(2, 3, 1)
```

See also [`colour_from_gluing`](@ref).

Role: convert
"""
function edge_gluing_from_colour(alpha::Integer)
    colour = _colour_mod3(alpha)
    ntuple(leg -> fig25_tau(colour, leg), 3)
end

"""
    colour_from_gluing(left_legs::NTuple{3,Integer}) -> Int

Return the colour represented by a leg permutation in the Figure 2.5 chart.

# Arguments
- `left_legs`: the images of the numbered right legs in the
  [edge chart](../specs/ograph-colours.md#chart). Each entry must convert to
  `Int`; conversion can throw `InexactError`. Checked: a tuple outside the
  chart's cyclic permutation subgroup throws `ArgumentError` naming that chart.

# Returns
`Int`, the normalized colour representative of the supplied gluing.

# Guarantees
1. `edge_gluing_from_colour(result) == left_legs`. [proved/agent: by construction]
2. For a normalized colour `a`, `colour_from_gluing(edge_gluing_from_colour(a)) == a`. [proved/agent: by construction]
3. The argument is not modified. [proved/agent: by construction]

# Conventions
Use the [Figure 2.5 chart](../conventions/OGraphColours.md), including its
draft geometric interpretation.

# Examples
```jldoctest
julia> colour_from_gluing(edge_gluing_from_colour(-1))
-1
```

Role: convert
"""
function colour_from_gluing(left_legs::NTuple{3,Integer})
    images = ntuple(leg -> Int(left_legs[leg]), 3)
    for alpha in (-1, 0, 1)
        edge_gluing_from_colour(alpha) == images && return alpha
    end
    throw(ArgumentError("gluing is not a power of (1 2 3) under the Fig. 2.5 chart"))
end

"""
    fig25_gauged_colour(alpha::Integer, right_shift::Integer, left_shift::Integer) -> Int

Return the colour read after cyclic changes of the two edge-leg numberings.

# Arguments
- `alpha`, `right_shift`, `left_shift`: the colour and numbering shifts in the
  [chart's gauge convention](../specs/ograph-colours.md#chart). Conversion to
  `Int` is required and can throw `InexactError`. The caller must ensure the
  intermediate sum and difference fit `Int`; overflow is not checked and can
  invalidate the gauge interpretation.

# Returns
`Int`, the normalized representative of the colour in the shifted chart.

# Guarantees
1. Under the arithmetic precondition, the result represents the same physical gluing in the shifted leg numberings. [proved/agent: by construction]
2. Zero shifts preserve a normalized input colour. [proved/agent: by construction]
3. Arguments are not modified. [proved/agent: by construction]

# Conventions
Use the [Figure 2.5 gauge convention](../conventions/OGraphColours.md).

# Not guaranteed
This does not apply a move or a gauge quotient to a stored coloured graph.

# Examples
```jldoctest
julia> Tuple(fig25_gauged_colour(a, 0, 0) for a in -1:1)
(-1, 0, 1)
```

Role: convert
"""
function fig25_gauged_colour(alpha::Integer, right_shift::Integer, left_shift::Integer)
    _colour_mod3(Int(alpha) + Int(left_shift) - Int(right_shift))
end

"""
    fig25_frame_reversal(leg::Integer) -> Int

Return the leg number in the reversed transverse frame of an edge drawing.

# Arguments
- `leg`: a leg in the [edge chart](../specs/ograph-colours.md#chart).
  Checked: an index outside `1:3` throws `ArgumentError` naming the leg range.

# Returns
The leg number after the chart's transverse frame reversal, as an `Int`.

# Guarantees
1. The returned label is the image under the frame reversal fixed by the linked chart. [proved/agent: by construction]
2. Applying the function twice returns the original leg. [proved/agent: by construction]
3. The argument is not modified. [proved/agent: by construction]

# Conventions
Use the [Figure 2.5 frame transport](../conventions/OGraphColours.md), whose
geometric interpretation remains draft.

# Not guaranteed
Reversing this label alone does not reverse a graph or alter its stored colour.

# Examples
```jldoctest
julia> all(fig25_frame_reversal(fig25_tau(-a, fig25_frame_reversal(i))) == fig25_tau(a, i) for a in -1:1 for i in 1:3)
true
```

Role: convert
"""
function fig25_frame_reversal(leg::Integer)
    leg in 1:3 || throw(ArgumentError("leg must be 1, 2, or 3"))
    leg == 2 ? 2 : 4 - Int(leg)
end

"""
    edge_disc_germs(colour::Integer, right_halfedge::Integer, mate::Tuple{Vararg{Int}}) -> NTuple{3,Tuple{DiscGerm,DiscGerm}}

Return the three paired disc germs of one coloured edge drawing.

# Arguments
- `colour`: an integer colour in the [edge chart](../specs/ograph-colours.md#chart),
  convertible to `Int`; conversion can throw `InexactError`.
- `right_halfedge`: checked to index `mate`; `ArgumentError` identifies an
  endpoint outside the half-edge set.
- `mate`: the caller supplies the edge involution of the same graph.
  The involution and the range of the left endpoint are not validated here;
  a nonpositive left endpoint fails the `DiscGerm` constructor with
  `ArgumentError`. Without valid edge data the result is only a list of labels.

# Returns
Three `(right_germ, left_germ)` pairs, ordered by right leg number, with each
germ represented by a `DiscGerm`.

# Guarantees
1. The pairs describe the gluing of the chosen edge in the linked chart. [proved/agent: by construction]
2. Every right germ has the supplied half-edge and every left germ has its mate. [proved/agent: by construction]
3. Arguments are not modified. [proved/agent: by construction]

# Conventions
The [Figure 2.5 drawing and frame transport](../conventions/OGraphColours.md)
apply; the geometric interpretation remains draft.

# Not guaranteed
The rest of the graph and O0–O3 are not certified.

# Examples
```jldoctest
julia> pairs = edge_disc_germs(1, 4, (2, 1, 4, 3));

julia> Tuple((a.leg, b.leg) for (a, b) in pairs)
((1, 2), (2, 3), (3, 1))
```

Role: construct
"""
function edge_disc_germs(colour::Integer, right_halfedge::Integer, mate::Tuple{Vararg{Int}})
    1 <= right_halfedge <= length(mate) || throw(ArgumentError(
        "drawing right end is outside the half-edge set"))
    left = mate[right_halfedge]
    alpha = _colour_mod3(colour)
    ntuple(leg -> (DiscGerm(right_halfedge, leg),
                   DiscGerm(left, fig25_tau(alpha, leg))), 3)
end

function _colour_check(id::AbstractString, passed::Bool, detail::AbstractString)
    message = passed ? "$id passed: $detail" : "$id failed: $detail"
    occursin(id, message) || error("colour check message lost its condition id")
    OGraphColourCheck(String(id), passed, message)
end

function _listed_halfedges(graph::ColouredOGraph)
    listed = Int[]
    for cyclic in graph.vertex_halfedges
        append!(listed, cyclic)
    end
    listed
end

function _colour_edges(graph::ColouredOGraph)
    edges = Tuple{Int,Int}[]
    reasons = String[]
    n = length(graph.mate)
    for (halfedge, other) in enumerate(graph.mate)
        if other < 1 || other > n || graph.mate[other] != halfedge || other == halfedge
            push!(reasons, "half-edge $halfedge is not paired")
            continue
        end
        halfedge < other && push!(edges, (halfedge, other))
    end
    edges, reasons
end

function _colour_components(graph::ColouredOGraph, edges::Vector{Tuple{Int,Int}})
    n = length(graph.vertex_halfedges)
    parent = collect(1:n)
    function find(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        vertex
    end
    owner = Dict{Int,Int}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for halfedge in cyclic
            owner[halfedge] = vertex
        end
    end
    for (left, right) in edges
        haskey(owner, left) || continue
        haskey(owner, right) || continue
        a = find(owner[left])
        b = find(owner[right])
        a == b || (parent[b] = a)
    end
    n == 0 ? 0 : length(unique(find(vertex) for vertex in 1:n))
end

function _check_o0(graph::ColouredOGraph)
    reasons = String[]
    listed = _listed_halfedges(graph)
    n = length(graph.mate)
    if sort(listed) != collect(1:n) || !allunique(listed)
        push!(reasons, "half-edges are not exactly the set 1:$n")
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        length(cyclic) == 4 || push!(reasons, "vertex $vertex is not quadrivalent")
    end
    edges, pair_reasons = _colour_edges(graph)
    append!(reasons, unique(pair_reasons))
    if isempty(reasons) && _colour_components(graph, edges) != 1
        push!(reasons, "graph is not connected")
    end
    isempty(reasons) && return _colour_check("O0", true, "connected finite quadrivalent graph")
    _colour_check("O0", false, join(reasons, "; "))
end

function _check_o1(graph::ColouredOGraph)
    if length(graph.vertex_halfedges) >= 1
        return _colour_check("O1", true, "at least one vertex")
    end
    _colour_check("O1", false, "an o-graph has at least one vertex")
end

function _o2_detail(cyclic::Tuple, positions::Tuple)
    length(cyclic) == 4 || return "missing an over pair"
    isempty(positions) && return "missing an over pair"
    any(position -> position < 1 || position > length(cyclic), positions) &&
        return "missing an over pair"
    marked = Set(positions)
    (marked == Set((1, 3)) || marked == Set((2, 4))) && return ""
    marked == Set((1, 2, 3, 4)) && return "more than one over pair"
    length(marked) == 2 && return "over branches are not opposite"
    return "missing an over pair"
end

function _check_o2(graph::ColouredOGraph)
    reasons = String[]
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        detail = _o2_detail(cyclic, graph.over_positions[vertex])
        isempty(detail) || push!(reasons, "vertex $vertex has $detail")
    end
    isempty(reasons) && return _colour_check("O2", true,
        "every vertex has one marked opposite over pair")
    _colour_check("O2", false, join(reasons, "; "))
end

function _recorded_colours(graph::ColouredOGraph)
    recorded = Dict{Tuple{Int,Int},Int}()
    reasons = String[]
    for record in graph.edge_colours
        right = record.halfedge
        if right < 1 || right > length(graph.mate)
            push!(reasons, "colour record does not select an edge")
            continue
        end
        left = graph.mate[right]
        if left < 1 || left > length(graph.mate) || graph.mate[left] != right || left == right
            push!(reasons, "colour record does not select an edge")
            continue
        end
        key = minmax(left, right)
        reduced = _colour_mod3(record.colour)
        if haskey(recorded, key) && recorded[key] != reduced
            push!(reasons, "edge $key has more than one colour")
        else
            recorded[key] = reduced
        end
    end
    recorded, reasons
end

function _check_o3(graph::ColouredOGraph)
    edges, _ = _colour_edges(graph)
    recorded, reasons = _recorded_colours(graph)
    for edge in edges
        haskey(recorded, edge) || push!(reasons, "edge $edge has no colour")
    end
    isempty(reasons) && return _colour_check("O3", true, "every edge has one colour in Z_3")
    _colour_check("O3", false, join(reasons, "; "))
end

"""
    validate_coloured_ograph(graph::ColouredOGraph) -> NTuple{4,OGraphColourCheck}

Report the O0–O3 conditions for a stored coloured-o-graph presentation.

# Arguments
`graph` is a value made by the `ColouredOGraph` constructor. It need not
satisfy O0–O3: those failures are returned as checks, not raised as exceptions.
The condition meanings are owned by the
[condition table](../specs/ograph-colours.md#condition-identifiers).

# Returns
Four `OGraphColourCheck` records in order `("O0", "O1", "O2", "O3")`.
Each has `id::String` the condition label, `passed::Bool` its outcome,
and `message::String` an explanation containing that label.

# Guarantees
1. Every record reports its named condition on the input data. Failure of
   one condition does not suppress the other three checks. [proved/agent: by construction]
2. All four passing means the data meet the repository's Definition 2.1.4
   coloured-o-graph conditions. This is a check of that combinatorial
   structure, not of a normal branching or closed manifold. [proved/agent: by construction]
3. The input, including its stored integer colours, is unchanged. [proved/agent: by construction]

# Conventions
The stored-colour interpretation is owned by
[the chart specification](../specs/ograph-colours.md#chart).

# Not guaranteed
This does not check N1–N2, C1–C3 or P1–P2, validate `convention_id` against
a known registry, or certify the face-cycle reconstruction or boundary flags.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> checks = validate_coloured_ograph(g);

julia> (map(c -> c.id, checks), all(c.passed for c in checks))
(("O0", "O1", "O2", "O3"), true)
```

See also [`ColouredOGraph`](@ref), [`validate_normal_colouring`](@ref).

Role: condition
"""
function validate_coloured_ograph(graph::ColouredOGraph)
    (_check_o0(graph), _check_o1(graph), _check_o2(graph), _check_o3(graph))
end

function _require_coloured_ograph(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "canonical form requires O0-O3: " * join((check.message for check in failed), "; ")))
    length(graph.vertex_halfedges) <= _COLOUR_CANON_MAX_VERTICES || throw(ArgumentError(
        "coloured o-graph canonical search supports at most $(_COLOUR_CANON_MAX_VERTICES) vertices"))
    graph
end

function _rotated_over_code(positions::Tuple, rotation::Integer)
    moved = Set(mod1(position - rotation, 4) for position in positions)
    moved == Set((1, 3)) && return 1
    moved == Set((2, 4)) && return 2
    throw(ArgumentError("over pair is not one opposite pair"))
end

function _permutations(n::Integer)
    results = Vector{Int}[]
    used = falses(n)
    current = Int[]
    function walk()
        if length(current) == n
            push!(results, copy(current))
            return
        end
        for vertex in 1:n
            used[vertex] && continue
            used[vertex] = true
            push!(current, vertex)
            walk()
            pop!(current)
            used[vertex] = false
        end
    end
    walk()
    results
end

function _rotations(n::Integer)
    results = Vector{Int}[]
    current = zeros(Int, n)
    function walk(index)
        if index > n
            push!(results, copy(current))
            return
        end
        for rotation in 0:3
            current[index] = rotation
            walk(index + 1)
        end
    end
    walk(1)
    results
end

function _canonical_candidate(graph::ColouredOGraph, permutation::Vector{Int},
                              rotations::Vector{Int}, colour_of_edge::Dict{Tuple{Int,Int},Int})
    n = length(permutation)
    new_of_old = Dict{Int,Int}()
    new_over = Int[]
    next_id = 1
    for new_vertex in 1:n
        old_vertex = permutation[new_vertex]
        cyclic = graph.vertex_halfedges[old_vertex]
        rotation = rotations[new_vertex]
        rotated = (cyclic[rotation + 1:end]..., cyclic[1:rotation]...)
        for halfedge in rotated
            new_of_old[halfedge] = next_id
            next_id += 1
        end
        push!(new_over, _rotated_over_code(graph.over_positions[old_vertex], rotation))
    end
    mate = zeros(Int, 4n)
    colour_by_end = Dict{Int,Int}()
    for (left, right) in keys(colour_of_edge)
        new_left = new_of_old[left]
        new_right = new_of_old[right]
        mate[new_left] = new_right
        mate[new_right] = new_left
        smaller = min(new_left, new_right)
        colour_by_end[smaller] = colour_of_edge[(left, right)]
    end
    ordered_colours = Int[]
    for halfedge in 1:4n
        halfedge < mate[halfedge] || continue
        push!(ordered_colours, colour_by_end[halfedge])
    end
    (graph.convention_id, Tuple(new_over), Tuple(mate), Tuple(ordered_colours))
end

"""
    canonical_coloured_ograph(graph::ColouredOGraph) -> ColouredOGraph

Return a canonical representative of the coloured graph under relabelling and cyclic-origin changes.

# Arguments
`graph` must satisfy O0–O3 in the
[condition table](../specs/ograph-colours.md#condition-identifiers) and have
at most six vertices. A failed condition throws `ArgumentError` listing
the condition failures; exceeding the bound throws `ArgumentError` naming
the six-vertex search limit.

# Returns
`ColouredOGraph` with the same field meanings as its constructor, representing
the same coloured graph with canonical labels and normalized colour records.

# Guarantees
1. The result satisfies O0–O3 and has the same vertex and edge counts and
   `convention_id` as the input. [proved/agent: by construction]
2. Applying this function again gives an equal result. Relabelling half-edges
   or vertices, rotating cyclic origins with their over marks, changing
   which endpoint stores a colour, and replacing a colour by a congruent
   integer give equal results when the other data are transported consistently.
   [proved/agent: by construction]
3. The result represents the same coloured graph, hence the same selected
   polyhedron under the [reconstruction conventions](../conventions/OGraphColours.md).
   The input is unchanged. [proved/agent: by construction]

# Conventions
Representative selection and the permitted changes are specified in
[Canonical form and serialization](../specs/ograph-colours.md#canonical-form-and-serialization).

# Not guaranteed
This does not classify spines or manifolds up to C- or MP-moves, quotient
colour gauge, or identify reflected cyclic orders. Exact canonical labels
are not a version-independent serialization contract.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> c = canonical_coloured_ograph(g);

julia> canonical_coloured_ograph(c) == c
true
```

See also [`coloured_ograph_hash`](@ref), [`coloured_ograph_serialization`](@ref).

Role: convert
"""
function canonical_coloured_ograph(graph::ColouredOGraph)
    _require_coloured_ograph(graph)
    n = length(graph.vertex_halfedges)
    colour_of_edge, _ = _recorded_colours(graph)
    best = nothing
    for permutation in _permutations(n)
        for rotations in _rotations(n)
            candidate = _canonical_candidate(graph, permutation, rotations, colour_of_edge)
            if best === nothing || candidate < best
                best = candidate
            end
        end
    end
    convention, over_codes, mate, colours = best
    vertices = Tuple(Tuple((4 * (vertex - 1) + 1):(4 * vertex)) for vertex in 1:n)
    overs = Tuple(code == 1 ? (1, 3) : (2, 4) for code in over_codes)
    records = EdgeColour[]
    colour_index = 1
    for halfedge in 1:length(mate)
        halfedge < mate[halfedge] || continue
        push!(records, EdgeColour(halfedge, colours[colour_index]))
        colour_index += 1
    end
    ColouredOGraph(vertices, mate, overs, records; convention_id=convention)
end

"""
    coloured_ograph_serialization(graph::ColouredOGraph) -> String

Return the canonical text serialization of a coloured-graph presentation.

# Arguments
`graph` must meet the O0–O3 and six-vertex preconditions of
[`canonical_coloured_ograph`](@ref). Their failures propagate as `ArgumentError`
identifying the condition failures or search bound.

# Returns
A `String` encoding the convention identifier, vertex count, over positions,
mate data and colour records of the canonical coloured graph.
The grammar, escaping scope and compatibility policy are in
[Canonical form and serialization](../specs/ograph-colours.md#canonical-form-and-serialization).

# Guarantees
1. Equal canonical representatives have equal serialization strings. [proved/agent: by construction]
2. Serializing a graph and its canonical representative gives the same text;
   repeated calls on the same input give the same text and leave the input unchanged.
   [proved/agent: by construction]
3. The encoded mathematical object is the coloured graph under the
   [chart conventions](../specs/ograph-colours.md#chart), with the quotient
   scope of `canonical_coloured_ograph`; it is not a manifold classification.
   [proved/agent: by construction]

# Not guaranteed
There is no public inverse parser or version-independent wire-format promise.
Arbitrary convention identifiers are not escaped, so this contract does not
promise unambiguous decoding for every user-supplied identifier.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> a = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> b = canonical_coloured_ograph(a);

julia> coloured_ograph_serialization(a) == coloured_ograph_serialization(b)
true
```

See also [`canonical_coloured_ograph`](@ref), [`coloured_ograph_hash`](@ref).

Role: convert
"""
function coloured_ograph_serialization(graph::ColouredOGraph)
    canonical = canonical_coloured_ograph(graph)
    colours = join((string(record.halfedge, ":", record.colour) for record in canonical.edge_colours), ",")
    mates = join(canonical.mate, ",")
    overs = join((join(positions, "+") for positions in canonical.over_positions), ",")
    "ColouredOGraph;convention=$(canonical.convention_id);n=$(length(canonical.vertex_halfedges));" *
        "over_positions=$overs;mate=$mates;colours=$colours"
end

"""
    coloured_ograph_hash(graph::ColouredOGraph) -> String

Return a SHA-256 fingerprint of the canonical coloured-graph presentation.

# Arguments
`graph` must meet the O0–O3 and six-vertex preconditions of
[`canonical_coloured_ograph`](@ref). Their failures propagate as `ArgumentError`
identifying the condition failures or search bound.

# Returns
A 64-character lowercase hexadecimal `String`, the SHA-256 digest of the
UTF-8 code units of [`coloured_ograph_serialization`](@ref).
Its compatibility scope is specified in
[Canonical form and serialization](../specs/ograph-colours.md#canonical-form-and-serialization).

# Guarantees
1. Equal serializations have equal fingerprints. Repeated calls are
   deterministic and leave the input unchanged. [proved/agent: by construction]
2. The fingerprint is invariant under the relabellings and cyclic-origin
   changes identified by `canonical_coloured_ograph`. It refers to the
   coloured graph, with the [chart scope](../specs/ograph-colours.md#chart).
   [proved/agent: by construction]

# Not guaranteed
Hash equality is not a collision-free proof of graph equality or a test of
homeomorphism of spines or manifolds. This does not quotient C- or MP-moves
or colour gauge. Changes to the serialization can change the fingerprint.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> a = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> b = canonical_coloured_ograph(a);

julia> coloured_ograph_hash(a) == coloured_ograph_hash(b)
true
```

See also [`coloured_ograph_serialization`](@ref), [`canonical_coloured_ograph`](@ref).

Role: invariant
"""
function coloured_ograph_hash(graph::ColouredOGraph)
    bytes2hex(SHA.sha256(codeunits(coloured_ograph_serialization(graph))))
end
