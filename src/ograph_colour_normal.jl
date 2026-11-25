# cf:[[BP97#^bp97-sec-0566]] Proposition 3.2.1 — colours and the orientation witness of a closed normal o-graph; adjustment: above = positive passage joins the repository Gauss code to the over branch; unverified

"""
    NormalEdgeOrientation(tail, head) -> NormalEdgeOrientation

Store a direction for one edge of a normal o-graph.

# Arguments
`tail::Int` and `head::Int`, or values convertible to `Int`. Conversion
failures propagate their Julia exceptions. This constructor has no graph
argument and performs no edge-membership check.

# Returns
`NormalEdgeOrientation` with `tail` the supplied first half-edge and `head`
the supplied second half-edge. Their use is the
[normal adapter](../specs/ograph-colours.md#normal-adapter).

# Guarantees
1. Both fields equal the supplied values after conversion. [proved/agent: by construction]
2. When the endpoints are mates in the intended graph, this records an
   orientation of that edge for the N2/P1 and P2 conditions in the
   [condition table](../specs/ograph-colours.md#condition-identifiers).
   [proved/agent: by construction]

# Not guaranteed
The constructor does not check that the endpoints are distinct, in range,
or paired. A collection of these records need not satisfy N2 or P2 and
does not itself construct the B7 disc orientations in the
[branching specification](../specs/ograph-colour-branching.md#condition-identifiers).

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")));

julia> edges = Tuple(NormalEdgeOrientation(e.tail, e.head) for e in c.edges);

julia> any(found -> Set((e.tail, e.head) for e in found) == Set((e.tail, e.head) for e in edges), prop321_orientations(c.graph))
true
```

See also [`normal_orientations`](@ref), [`prop321_orientations`](@ref).

Role: construct
"""
struct NormalEdgeOrientation
    tail::Int
    head::Int
end

"""
    ClosedNormalColouring(graph::ColouredOGraph, edges, colours, witness::ClosednessWitness, state::OGraphState) -> ClosedNormalColouring

Store a coloured graph together with a normal orientation, colour sequence and closedness witness.

# Arguments
- `graph::ColouredOGraph`, `edges::Tuple{Vararg{NormalEdgeOrientation}}`,
  `colours::Tuple{Vararg{Int}}`, `witness::ClosednessWitness`, `state::OGraphState`:
  the fields to store; the converting constructor accepts values convertible
  to these field types. Failed conversions propagate their Julia exceptions.
- For a certified conversion use [`normal_coloured_ograph`](@ref). Direct
  construction does not check N1, N2, C1–C3 or P1–P2 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers).

# Returns
`ClosedNormalColouring` with `graph` the coloured graph, `edges` its supplied
orientation, `colours` the supplied colour sequence, `witness` the supplied
closedness evidence, and `state` the supplied normal presentation.
The sequence correspondence is the [normal adapter](../specs/ograph-colours.md#normal-adapter).

# Guarantees
1. Each field contains the supplied value after conversion to its declared type.
   [proved/agent: by construction]
2. When the fields come from `normal_coloured_ograph`, they record the
   Proposition 3.2.1 colouring of the represented normal spine, under the
   [adopted adapter conventions](../conventions/OGraphColours.md).
   The constructor alone does not establish that correspondence.
   [proved/agent: by construction]

# Conventions
Use the normal adapter and colour ordering linked above.

# Not guaranteed
The constructor does not verify that `graph`, `edges`, `colours`, `state`
and `witness` describe the same input. It does not construct disc orientations
or strengthen the geometric evidence status on the convention page.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")));

julia> d = ClosedNormalColouring(c.graph, c.edges, c.colours, c.witness, c.state);

julia> (d.colours, all(check.passed for check in validate_normal_colouring(d)))
((1, -1), true)
```

See also [`normal_coloured_ograph`](@ref), [`validate_normal_colouring`](@ref).

Role: convert
"""
struct ClosedNormalColouring
    graph::ColouredOGraph
    edges::Tuple{Vararg{NormalEdgeOrientation}}
    colours::Tuple{Vararg{Int}}
    witness::ClosednessWitness
    state::OGraphState
end

"""
    prop321_colour(first_above::Bool, second_above::Bool) -> Int

Return the Proposition 3.2.1 colour determined by two endpoint branch marks.

# Arguments
`first_above` and `second_above` are the endpoint marks in the
[normal-adapter convention](../specs/ograph-colours.md#normal-adapter).
Their correspondence to the intended directed edge is the caller's
responsibility; the function only receives Boolean values.

# Returns
The normalized integer colour representative prescribed by P2 in the
[condition table](../specs/ograph-colours.md#condition-identifiers).

# Guarantees
1. The result is the P2 colour for the supplied ordered marks. [proved/agent: by construction]
2. Arguments are not modified. [proved/agent: by construction]

# Conventions
Use [OGraphColours](../conventions/OGraphColours.md) for the adapter's reading
of above/below marks and its draft source correspondence.

# Not guaranteed
No graph, normal orientation or closedness condition is checked.

# Examples
```jldoctest
julia> (prop321_colour(true,true), prop321_colour(false,true), prop321_colour(true,false))
(0, 1, -1)
```

Role: convert
"""
function prop321_colour(first_above::Bool, second_above::Bool)
    first = first_above ? 0 : 1
    second = second_above ? 0 : 1
    _colour_mod3(first - second)
end

function _normal_port_slot(vertex_sign::Int, port::Symbol)
    vertex_sign == 1 && port === :mi && return 1
    vertex_sign == 1 && port === :po && return 2
    vertex_sign == 1 && port === :mo && return 3
    vertex_sign == 1 && port === :pi && return 4
    vertex_sign == -1 && port === :mi && return 1
    vertex_sign == -1 && port === :pi && return 2
    vertex_sign == -1 && port === :mo && return 3
    vertex_sign == -1 && port === :po && return 4
    throw(ArgumentError("vertex sign and port do not name a Fig. 3.4 or Fig. 3.5 slot"))
end

_normal_halfedge(vertex::Int, slot::Int) = (vertex - 1) * 4 + slot

function _normal_pieces(state::OGraphState)
    n = true_vertex_count(state)
    n >= 1 || throw(ArgumentError("N1 failed: an o-graph has at least one vertex"))
    link = _ograph_virtual_link(state)
    table = real_edge_table(link)
    length(table) == length(state.gauss) || throw(ArgumentError(
        "real edge table does not have one edge per Gauss entry"))
    vertices = Tuple(ntuple(slot -> _normal_halfedge(vertex, slot), 4) for vertex in 1:n)
    overs = Tuple((2, 4) for _ in 1:n)
    mate = zeros(Int, 4n)
    edges = NormalEdgeOrientation[]
    colours = Int[]
    for edge in table
        tail_label, head_label = edge.crossing
        (tail_label == 0 || head_label == 0) && throw(ArgumentError(
            "a closed normal colouring has no cut endpoint"))
        tail_vertex = abs(tail_label)
        head_vertex = abs(head_label)
        tail_port = tail_label > 0 ? :po : :mo
        head_port = head_label > 0 ? :pi : :mi
        tail = _normal_halfedge(tail_vertex, _normal_port_slot(state.vertex_signs[tail_vertex], tail_port))
        head = _normal_halfedge(head_vertex, _normal_port_slot(state.vertex_signs[head_vertex], head_port))
        (tail == head || mate[tail] != 0 || mate[head] != 0) && throw(ArgumentError(
            "Gauss passages do not give four distinct half-edges at each vertex"))
        mate[tail] = head
        mate[head] = tail
        push!(edges, NormalEdgeOrientation(tail, head))
        push!(colours, prop321_colour(tail_label > 0, head_label > 0))
    end
    all(!iszero, mate) || throw(ArgumentError("a half-edge was left unpaired"))
    records = Tuple(EdgeColour(edge.tail, colour) for (edge, colour) in zip(edges, colours))
    graph = ColouredOGraph(vertices, Tuple(mate), overs, records)
    graph, Tuple(edges), Tuple(colours)
end

"""
    normal_coloured_ograph(state::OGraphState) -> ClosedNormalColouring

Return the Proposition 3.2.1 colouring and normal orientation of a closed normal presentation.

# Arguments
`state` must satisfy N1–N2 and C1–C3 in the
[condition table](../specs/ograph-colours.md#condition-identifiers).
N1 and operational closedness are checked; `ArgumentError` identifies the
failed condition. The converted data must also pass O0–O3 and the normal
colour checks; otherwise `ArgumentError` reports those failed checks.

# Returns
`ClosedNormalColouring` containing the coloured `graph`, directed `edges`,
edge `colours`, replayed closedness `witness`, and the supplied `state`.

# Guarantees
1. The returned graph and orientation pass `validate_coloured_ograph` and `validate_normal_colouring`. [proved/agent: by construction]
2. The colour data express P2 in the normal adapter's chart; their reconstruction meaning is conditional on the draft dependencies recorded in Slot transcription. [proved/agent: by construction]
3. The input state is not modified; its labelled presentation is retained. [proved/agent: by construction]

# Conventions
Use the [normal adapter](../specs/ograph-colours.md#normal-adapter) and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
This does not construct B7 disc orientations or discharge the reconstruction
and C3 dependencies listed in
[Slot transcription](../specs/ograph-colours.md#slot-transcription).
The crossing-free extension is rejected by N1.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")));

julia> all(check.passed for check in validate_normal_colouring(c))
true
```

See also [`face_cycles`](@ref), [`thickening_boundary`](@ref).

Role: convert
"""
function normal_coloured_ograph(state::OGraphState)
    true_vertex_count(state) >= 1 || throw(ArgumentError(
        "N1 failed: an o-graph has at least one vertex"))
    witness = closedness_witness(state)
    all(witness.criteria) || throw(ArgumentError(
        "C1-C3 failed under $(witness.verifier_version): $(witness.criteria)"))
    graph, edges, colours = _normal_pieces(state)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "converted o-graph failed Definition 2.1.4: " *
        join((check.message for check in failed), "; ")))
    colouring = ClosedNormalColouring(graph, edges, colours, witness, state)
    checks = validate_normal_colouring(colouring)
    failed_normal = filter(check -> !check.passed, checks)
    isempty(failed_normal) || throw(ArgumentError(
        "converted o-graph failed the normal colour checks: " *
        join((check.message for check in failed_normal), "; ")))
    colouring
end

"""
    normal_coloured_ograph(closed::ClosedNormalOGraph) -> ClosedNormalColouring

Return the Proposition 3.2.1 colouring of a certified closed normal presentation.

# Arguments
`closed` supplies the stored state. The N1, operational C1–C3, O0–O3 and
normal colour checks of [`normal_coloured_ograph(::OGraphState)`](@ref normal_coloured_ograph) are
performed anew; failures throw `ArgumentError` naming the failed conditions.

# Returns
A `ClosedNormalColouring` with the same field meanings as the state overload.

# Guarantees
1. The result is the normal adapter's colouring of `closed.state`, with fresh validation results. [proved/agent: by construction]
2. `closed` is not modified. [proved/agent: by construction]

# Conventions
Use the [normal adapter](../specs/ograph-colours.md#normal-adapter) and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
The [draft reconstruction dependencies](../specs/ograph-colours.md#slot-transcription)
remain; no B7 disc orientations are returned.

# Examples
```jldoctest
julia> c = normal_coloured_ograph(certify_closed_ograph(lens_space_link(2, 1)));

julia> all(check.passed for check in validate_normal_colouring(c))
true
```

Role: convert
"""
normal_coloured_ograph(closed::ClosedNormalOGraph) = normal_coloured_ograph(closed.state)

"""
    normal_coloured_ograph(link::VirtualLink) -> ClosedNormalColouring

Return the Proposition 3.2.1 colouring of a closed normal link presentation.

# Arguments
`link` must be accepted by [`ograph_state`](@ref), including its unweighted,
closed, single-circuit and coefficient-mode restrictions. The resulting state
must satisfy the checked conditions of the state overload. Their
`ArgumentError` domain failures and representation errors propagate.

# Returns
A `ClosedNormalColouring` with the graph, orientation, colours, witness and
state of the link's real-Gauss presentation.

# Guarantees
1. The result is the normal adapter's colouring of `ograph_state(link)`. [proved/agent: by construction]
2. `link` is not modified. [proved/agent: by construction]

# Conventions
Use the [normal adapter](../specs/ograph-colours.md#normal-adapter) and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
The [draft reconstruction dependencies](../specs/ograph-colours.md#slot-transcription)
remain; weights and B7 disc orientations are not transported.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(vl_example("s3_abalone"));

julia> all(check.passed for check in validate_normal_colouring(c))
true
```

Role: convert
"""
function normal_coloured_ograph(link::VirtualLink)
    normal_coloured_ograph(ograph_state(link))
end

# cf:[[BP97#^bp97-sec-0536]] N1 and N2 — opposite-edge directions, checked apart from closedness_witness; unverified

function _check_n1(graph::ColouredOGraph)
    reasons = String[]
    if isempty(graph.vertex_halfedges)
        push!(reasons, "an o-graph has at least one vertex")
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        length(cyclic) == 4 || push!(reasons, "vertex $vertex is not a normal crossing")
        detail = _o2_detail(cyclic, graph.over_positions[vertex])
        isempty(detail) || push!(reasons, "vertex $vertex has $detail")
    end
    isempty(reasons) && return _colour_check("N1", true,
        "at least one vertex, each a marked normal crossing")
    _colour_check("N1", false, join(reasons, "; "))
end

function _above_halfedges(graph::ColouredOGraph)
    above = falses(length(graph.mate))
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        overs = Set(graph.over_positions[vertex])
        for (position, halfedge) in enumerate(cyclic)
            1 <= halfedge <= length(above) || continue
            above[halfedge] = position in overs
        end
    end
    above
end

function _orientation_is_pairing(graph::ColouredOGraph, edges)
    seen = Set{Tuple{Int,Int}}()
    for edge in edges
        tail, head = edge.tail, edge.head
        if tail < 1 || head < 1 || tail > length(graph.mate) || head > length(graph.mate)
            return false
        end
        if graph.mate[tail] != head || graph.mate[head] != tail || tail == head
            return false
        end
        push!(seen, minmax(tail, head))
    end
    length(seen) == length(edges) || return false
    for halfedge in 1:length(graph.mate)
        other = graph.mate[halfedge]
        halfedge < other || continue
        minmax(halfedge, other) in seen || return false
    end
    true
end

function _directions_match(graph::ColouredOGraph, edges)
    _orientation_is_pairing(graph, edges) || return false
    tail_at = falses(length(graph.mate))
    for edge in edges
        tail_at[edge.tail] = true
    end
    for cyclic in graph.vertex_halfedges
        length(cyclic) == 4 || return false
        for pair in ((1, 3), (2, 4))
            count(position -> tail_at[cyclic[position]], pair) == 1 || return false
        end
    end
    true
end

function _check_n2(graph::ColouredOGraph, edges)
    _directions_match(graph, edges) && return _colour_check("N2", true,
        "directions of opposite edges match through the vertices")
    _colour_check("N2", false, "directions of opposite edges do not match through the vertices")
end

function _check_p1(graph::ColouredOGraph, edges)
    _directions_match(graph, edges) && return _colour_check("P1", true,
        "edge orientations match across marked vertices")
    _colour_check("P1", false, "edge orientations do not match across marked vertices")
end

function _check_p2(graph::ColouredOGraph, edges)
    recorded, reasons = _recorded_colours(graph)
    isempty(reasons) || return _colour_check("P2", false, join(reasons, "; "))
    above = _above_halfedges(graph)
    mismatches = String[]
    for edge in edges
        key = minmax(edge.tail, edge.head)
        if !haskey(recorded, key)
            push!(mismatches, "edge $key has no colour")
            continue
        end
        expected = prop321_colour(above[edge.tail], above[edge.head])
        expected == recorded[key] || push!(mismatches, "edge $key colour does not match its ends")
    end
    isempty(mismatches) && return _colour_check("P2", true,
        "every edge colour agrees with its first and second ends")
    _colour_check("P2", false, join(mismatches, "; "))
end

function _check_closed(witness::ClosednessWitness, index::Int, id::String)
    passed = witness.criteria[index]
    passed && return _colour_check(id, true, "closedness witness criterion $(id)")
    _colour_check(id, false, "closedness witness rejected $(id)")
end

"""
    validate_normal_colouring(colouring::ClosedNormalColouring) -> NTuple{7,OGraphColourCheck}

Report normal-orientation and colour compatibility together with the stored closedness criteria.

# Arguments
`colouring` should come from [`normal_coloured_ograph`](@ref), or contain
mutually corresponding fields with in-range half-edge references. It may
fail the named conditions. Such failures are returned as checks. Arbitrarily
malformed fields from the unchecked constructor can instead throw indexing
or lookup errors; this is not a general structural validator.

# Returns
Seven `OGraphColourCheck` records in the order
`("N1", "N2", "C1", "C2", "C3", "P1", "P2")`. Each record's `id` names
the condition, `passed` is the reported outcome, and `message` contains the id.
The conditions are owned by the
[condition table](../specs/ograph-colours.md#condition-identifiers).

# Guarantees
1. N1, N2, P1 and P2 concern `colouring.graph` and `colouring.edges`.
   The N2 and P1 results agree. The C1, C2 and C3 outcomes equal the three
   entries of `colouring.witness.criteria` in order. [proved/agent: by construction]
2. On a coherent value from `normal_coloured_ograph`, the checks concern
   the Proposition 3.2.1 colouring of that normal presentation, with the
   [adapter interpretation](../specs/ograph-colours.md#normal-adapter).
   [proved/agent: by construction]
3. The supplied colouring is unchanged. [proved/agent: by construction]

# Not guaranteed
C1–C3 are not recomputed. No correspondence of `witness` or `state` to
`graph` is checked. `colouring.colours` is not compared with graph colours
or the Gauss word. All seven passing on a manually assembled value therefore
does not certify that its fields come from one closed normal presentation.
This function does not construct B7 disc orientations or verify the
geometric reconstruction.

# Sources
cf:[[BP97#^bp97-sec-0536]] N1–N2 and C1–C3, PDF p.15 — the named normal and closed-normal conditions; this function reports stored C1–C3 outcomes rather than proving them again.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")));

julia> checks = validate_normal_colouring(c);

julia> (map(c -> c.id, checks), all(c.passed for c in checks))
(("N1", "N2", "C1", "C2", "C3", "P1", "P2"), true)
```

See also [`normal_coloured_ograph`](@ref), [`validate_coloured_ograph`](@ref).

Role: condition
"""
function validate_normal_colouring(colouring::ClosedNormalColouring)
    graph = colouring.graph
    edges = colouring.edges
    (
        _check_n1(graph),
        _check_n2(graph, edges),
        _check_closed(colouring.witness, 1, "C1"),
        _check_closed(colouring.witness, 2, "C2"),
        _check_closed(colouring.witness, 3, "C3"),
        _check_p1(graph, edges),
        _check_p2(graph, edges),
    )
end

const _NORMAL_ORIENTATION_MAX_EDGES = 12

function _paired_edges(graph::ColouredOGraph)
    edges = Tuple{Int,Int}[]
    for halfedge in 1:length(graph.mate)
        other = graph.mate[halfedge]
        halfedge < other && push!(edges, (halfedge, other))
    end
    edges
end

"""
    normal_orientations(graph::ColouredOGraph) -> Vector

Return every N2-compatible orientation of the coloured graph's underlying edges.

# Arguments
`graph` must satisfy O0, O1 and O2 in the
[condition table](../specs/ograph-colours.md#condition-identifiers) and have
at most twelve edges. Failed conditions throw `ArgumentError` naming O0–O2;
an excessive edge count throws `ArgumentError` naming the search bound.
O3 is not a precondition of this orientation search.

# Returns
A vector of tuples of `NormalEdgeOrientation`; each tuple has one entry per
edge. Each entry has the `tail` and `head` fields described by that type.

# Guarantees
1. Each tuple orients every underlying edge exactly once and satisfies N2.
   [proved/agent: by construction]
2. Every orientation satisfying N2 is returned exactly once, for inputs
   within the supported bound. The result concerns the graph-edge directions
   in [the normal adapter](../specs/ograph-colours.md#normal-adapter).
   [proved/agent: by construction]
3. The result does not depend on stored colours. The input graph is unchanged.
   [proved/agent: by construction]

# Not guaranteed
These orientations need not satisfy P2 or closed-normal conditions. They are
not constructed B7 disc orientations of the spine; see
[the branching scope](../specs/ograph-colour-branching.md#what-is-recovered).
The enumeration order is not a canonical orientation identifier.

# Sources
cf:[[BP97#^bp97-sec-0536]] N2, PDF p.15 — normal edge directions; this function enumerates directions, not disc branchings.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> length(normal_orientations(g))
2
```

See also [`prop321_orientations`](@ref), [`validate_normal_colouring`](@ref).

Role: convert
"""
function normal_orientations(graph::ColouredOGraph)
    checks = validate_coloured_ograph(graph)
    passed = Dict(check.id => check.passed for check in checks)
    (passed["O0"] && passed["O1"] && passed["O2"]) || throw(ArgumentError(
        "normal orientations require O0, O1, and O2"))
    edges = _paired_edges(graph)
    m = length(edges)
    m <= _NORMAL_ORIENTATION_MAX_EDGES || throw(ArgumentError(
        "normal orientation search supports at most $(_NORMAL_ORIENTATION_MAX_EDGES) edges"))
    found = Vector{NTuple{m,NormalEdgeOrientation}}()
    for mask in 0:(1 << m) - 1
        oriented = Vector{NormalEdgeOrientation}(undef, m)
        for index in 1:m
            left, right = edges[index]
            if (mask >> (index - 1)) & 1 == 1
                oriented[index] = NormalEdgeOrientation(left, right)
            else
                oriented[index] = NormalEdgeOrientation(right, left)
            end
        end
        _directions_match(graph, oriented) || continue
        push!(found, Tuple(oriented))
    end
    found
end

"""
    prop321_orientations(graph::ColouredOGraph) -> Vector

Return the N2-compatible edge orientations whose P2 colours agree with the graph.

# Arguments
`graph` has the O0–O2 and twelve-edge preconditions of
[`normal_orientations`](@ref); violations propagate its `ArgumentError`.
For the Proposition 3.2.1 interpretation, require O3 as well; see the
[condition table](../specs/ograph-colours.md#condition-identifiers).
O3 is not separately rejected here: missing or conflicting colour records
cause the candidate orientations to fail the colour comparison.

# Returns
A vector of tuples of `NormalEdgeOrientation`, one record per edge in each
tuple. An empty vector means no enumerated orientation passed P2.

# Guarantees
1. The result is precisely the subset of `normal_orientations(graph)` that
   satisfies the stored-colour comparison P2; each is returned once.
   [proved/agent: by construction]
2. On O0–O3 inputs, these are the graph-edge witnesses for Proposition 3.2.1
   under the [branching interpretation](../specs/ograph-colour-branching.md#what-is-recovered).
   The input graph is unchanged. [proved/agent: by construction]

# Conventions
P2 and its adapter to over passages are owned by the
[condition table](../specs/ograph-colours.md#condition-identifiers) and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
This does not construct B7 disc orientations, impose C1–C3 or return a unique
normal presentation. On a graph failing O3, an empty result is not evidence
that a valid coloured spine has no compatible branching.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> length(prop321_orientations(g))
1
```

See also [`normal_orientations`](@ref), [`validate_normal_colouring`](@ref).

Role: condition
"""
function prop321_orientations(graph::ColouredOGraph)
    filter(edges -> _check_p2(graph, edges).passed, normal_orientations(graph))
end

"""
    coloured_gc_display(colouring::ClosedNormalColouring) -> String

Return the coloured Gauss-code display of a stored normal presentation.

# Arguments
`colouring` supplies the normal state and colour sequence. Their lengths
must agree; otherwise `ArgumentError` states that one colour per Gauss entry
is required. No graph, orientation or closedness checks are run here.

# Returns
A `String` containing the stored Gauss word, vertex signs and colour sequence
in the [normal-adapter display](../specs/ograph-colours.md#normal-adapter).

# Guarantees
1. Repeated calls on the same fields give the same string and leave the
   colouring unchanged. [proved/agent: by construction]
2. For a value from `normal_coloured_ograph`, the text displays that normal
   presentation together with its Proposition 3.2.1 colours. It describes
   the stored presentation, not an equivalence class of spines or manifolds.
   [proved/agent: by construction]

# Conventions
The display's ordering is owned by the normal adapter linked above and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
This is not a canonical identifier or a validation result. Inconsistent
fields passed through the direct `ClosedNormalColouring` constructor can
be displayed when the single length check passes.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> c = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")));

julia> coloured_gc_display(c)
"coloured GC [[[+1,-1]];[+];[1,-1]]"
```

See also [`normal_coloured_ograph`](@ref), [`validate_normal_colouring`](@ref).

Role: convert
"""
function coloured_gc_display(colouring::ClosedNormalColouring)
    state = colouring.state
    length(colouring.colours) == length(state.gauss) || throw(ArgumentError(
        "coloured GC display needs one colour per Gauss entry"))
    passages = join(map(label -> label > 0 ? "+$(label)" : string(label), state.gauss), ",")
    signs = join(map(sign -> sign > 0 ? "+" : "−", state.vertex_signs), ",")
    colours = join(map(string, colouring.colours), ",")
    "coloured GC [[[$(passages)]];[$(signs)];[$(colours)]]"
end
