# cf:[[BP97#^bp97-sec-0549]] Fig. 2.5 — discs pasted after the edge gluings; adjustment: leg 1 is candidate B, derived from Proposition 3.2.1, draft
# cf:[[BP97#^bp97-sec-0566]] Proposition 3.2.1 — the P2 colouring represents the same spine; candidate A is refuted by the closed-spine count

import AbstractAlgebra

"""
    FaceCycle(germs, corners, edge_letters) -> FaceCycle

Store the germ, corner and signed-edge data of a disc boundary.

# Arguments
- `germs::Vector{DiscGerm}`, `corners::Vector{NTuple{3,Int}}`, and
  `edge_letters::Vector{Int}`, or values convertible to those types.
  Failed conversions propagate their Julia exceptions.
- The constructor does not check that the three sequences describe one
  disc boundary; obtain computed cycles from [`face_cycles`](@ref).

# Returns
`FaceCycle` with `germs` the departing sheet-end labels, `corners` the
arriving vertex-corner labels and `edge_letters` the signed singular-edge
indices, in boundary traversal order for a cycle returned by `face_cycles`.
The corner representation is in
[Face cycles](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. Fields preserve the supplied values after conversion. [proved/agent: by construction]
2. A record supplied by `face_cycles` describes one attaching boundary of
   the selected face-cycle polyhedron, with the geometric interpretation
   subject to that operation's [calibration](../conventions/OGraphColours.md).
   [proved/agent: by construction]
3. Input vectors of the declared types are retained without copying. [proved/agent: by construction]

# Conventions
Signed-edge numbering is owned by
[the invariant specification](../specs/ograph-colour-invariants.md#closed-examples).

# Not guaranteed
Direct construction does not check equal lengths, closure, graph membership
or a choice of branching on the disc. It does not establish a preferred
cyclic origin or orientation for comparing records from different graphs.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> p = face_cycles(normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph);

julia> copies = [FaceCycle(c.germs, c.corners, c.edge_letters) for c in p.cycles];

julia> Tuple(sort([length(c.germs) for c in copies]))
(1, 5)
```

See also [`face_cycles`](@ref).

Role: convert
"""
struct FaceCycle
    germs::Vector{DiscGerm}
    corners::Vector{NTuple{3,Int}}
    edge_letters::Vector{Int}
end

"""
    ColouredPolyhedron(vertex_count, edge_count, disc_count, euler, halfedge_count, germ_count, boundary_length, cycles, calibration, status) -> ColouredPolyhedron

Store face-cycle data and cell counts for a coloured-o-graph polyhedron.

# Arguments
- The first seven arguments are `Int` counts; `cycles::Vector{FaceCycle}`,
  `calibration::Symbol`, and `status::String` supply the remaining fields.
  Values convertible to those field types are also accepted; conversion
  failures propagate their Julia exceptions.
- For computed counts use [`face_cycles`](@ref). Direct construction does not
  check the counts, cycles, calibration or status against each other.

# Returns
`ColouredPolyhedron` with `vertex_count`, `edge_count`, `disc_count` the
supplied cell counts, `euler` the supplied Euler characteristic,
`halfedge_count` the half-edge count, `germ_count` the counted germs,
`boundary_length` the total disc-boundary length, `cycles` the face-cycle
records, `calibration` the chart selection and `status` its evidence label.
These quantities refer to [face cycles](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. Fields preserve the supplied values after field-type conversion. [proved/agent: by construction]
2. When supplied by `face_cycles`, the record describes the selected
   face-cycle polyhedron. Its interpretation as the represented standard
   spine has the [calibration dependencies](../conventions/OGraphColours.md)
   of that operation. [proved/agent: by construction]
3. A supplied `Vector{FaceCycle}` is retained without copying. [proved/agent: by construction]

# Not guaranteed
Direct construction does not establish an Euler relation, coverage of germs,
realizability as a spine, or the validity of the supplied status label.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> p = face_cycles(normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph);

julia> q = ColouredPolyhedron(p.vertex_count, p.edge_count, p.disc_count, p.euler, p.halfedge_count, p.germ_count, p.boundary_length, p.cycles, p.calibration, p.status);

julia> (q.vertex_count, q.edge_count, q.disc_count, q.euler, q.germ_count)
(1, 2, 2, 1, 6)
```

See also [`face_cycles`](@ref).

Role: convert
"""
struct ColouredPolyhedron
    vertex_count::Int
    edge_count::Int
    disc_count::Int
    euler::Int
    halfedge_count::Int
    germ_count::Int
    boundary_length::Int
    cycles::Vector{FaceCycle}
    calibration::Symbol
    status::String
end

"""
    ColourHomology(free_rank, torsion_factors) -> ColourHomology

Store the free rank and torsion invariant factors of an abelian group.

# Arguments
- `free_rank::Int` and `torsion_factors::Vector{BigInt}`, or values convertible
  to these field types. Conversion failures propagate their Julia exceptions.
  The constructor does not check nonnegativity or the invariant-factor conditions.

# Returns
`ColourHomology` with `free_rank` the supplied free rank and `torsion_factors`
the supplied cyclic-summand orders. For a homology computation use
[`cellular_homology`](@ref).

# Guarantees
1. The fields equal the supplied values after field-type conversion. [proved/agent: by construction]
2. For data returned by `cellular_homology`, this record describes the first
   integral homology of the selected face-cycle polyhedron, with the
   [invariant-factor convention](../specs/ograph-colour-invariants.md#closed-examples).
   Direct construction supplies data, not a homology computation. [proved/agent: by construction]
3. An input `Vector{BigInt}` is retained, so later mutations are visible in the record.
   [proved/agent: by construction]

# Not guaranteed
Arbitrary fields need not be invariant factors of a homology group. The
record neither identifies a manifold nor certifies a geometric realization.

# Sources
cf:[[BP97#^bp97-sec-0703]] §10.1 — homology of the spine and its associated manifold; storage of invariant factors is a repository representation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> h = cellular_homology(g); stored = ColourHomology(h.free_rank, h.torsion_factors);

julia> (stored.free_rank, Tuple(stored.torsion_factors))
(0, ())
```

See also [`cellular_homology`](@ref).

Role: construct
"""
struct ColourHomology
    free_rank::Int
    torsion_factors::Vector{BigInt}
end

"""
    DualFacePairing(source_vertex, source_slot, target_vertex, target_slot, orientation_reversing) -> DualFacePairing

Store the two faces and orientation flag of a dual face pair.

# Arguments
- The four face indices are `Int` and `orientation_reversing` is `Bool`,
  or values convertible to those types. Failed conversions propagate their
  Julia exceptions. No range or geometric check is performed.

# Returns
`DualFacePairing` with `source_vertex`, `source_slot`, `target_vertex`,
`target_slot` the supplied face labels, and `orientation_reversing` the
supplied reversal flag. Face labels and the flag use
[Face cycles](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. Fields equal the supplied values after conversion. [proved/agent: by construction]
2. For records produced by [`dual_face_pairings`](@ref), the labels specify
   a pair of faces in the slot model of the dual ideal triangulation of the
   selected polyhedron. Direct construction does not compute that gluing.
   [proved/agent: by construction]

# Conventions
The face orientation and endpoint order are in
[OGraphColours](../conventions/OGraphColours.md#四面体と-link-triangle-の向き).

# Not guaranteed
The constructor does not order endpoints, compute the reversal flag, or
check that a face pairing exists. The record does not contain the vertex
bijection of the face gluing and cannot by itself reconstruct it.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — ideal triangulations and orientation-reversing face identifications; slot labels are the repository representation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> p = first(dual_face_pairings(g));

julia> q = DualFacePairing(p.source_vertex, p.source_slot, p.target_vertex, p.target_slot, p.orientation_reversing);

julia> q.orientation_reversing
true
```

See also [`dual_face_pairings`](@ref).

Role: construct
"""
struct DualFacePairing
    source_vertex::Int
    source_slot::Int
    target_vertex::Int
    target_slot::Int
    orientation_reversing::Bool
end

"""
    SpineBoundaryComponent(euler, face_count, edge_count, vertex_count, orientable, closed_surface) -> SpineBoundaryComponent

Store the counts and flags reported for one boundary component.

# Arguments
`euler`, `face_count`, `edge_count` and `vertex_count` are `Int` values;
`orientable` and `closed_surface` are `Bool` values. Convertible arguments
are also accepted; conversion failures propagate their Julia exceptions.
The constructor performs no surface or consistency checks.

# Returns
`SpineBoundaryComponent` with `euler` the supplied Euler characteristic,
`face_count`, `edge_count`, `vertex_count` the supplied cell counts, and
`orientable`, `closed_surface` the supplied boundary flags. Computed records
come from [`thickening_boundary`](@ref), whose space is described in
[Face cycles](../specs/ograph-colours.md#face-cycles). The meanings and scope
of the flags are owned by
[Boundary diagnostic flags](../specs/ograph-colours.md#boundary-diagnostic-flags).

# Guarantees
1. Each field equals the supplied value after conversion. [proved/agent: by construction]
2. For output of `thickening_boundary`, a record corresponds to one connected
   component of its corner-link model. Interpreting it as a component of
   the polyhedron's regular-neighbourhood boundary requires that operation's
   geometric hypotheses. [proved/agent: by construction]

# Not guaranteed
Neither the Euler relation nor either Boolean flag is checked on direct
construction. A supplied `orientable=true` is not a proof of orientability.
The record contains counts and flags, not the surface gluing itself.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — dual ideal triangulations; the corner-link boundary interpretation is reconstructed.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> b = only(thickening_boundary(g).components);

julia> c = SpineBoundaryComponent(b.euler, b.face_count, b.edge_count, b.vertex_count, b.orientable, b.closed_surface);

julia> (c.euler, c.face_count, c.edge_count, c.vertex_count, c.orientable, c.closed_surface)
(2, 4, 6, 4, true, true)
```

See also [`thickening_boundary`](@ref), [`SpineBoundary`](@ref).

Role: construct
"""
struct SpineBoundaryComponent
    euler::Int
    face_count::Int
    edge_count::Int
    vertex_count::Int
    orientable::Bool
    closed_surface::Bool
end

"""
    SpineBoundary(components, euler_sum) -> SpineBoundary

Store boundary-component records and their reported total Euler characteristic.

# Arguments
`components::Vector{SpineBoundaryComponent}` and `euler_sum::Int`, or values
convertible to these types. Conversion failures propagate their Julia
exceptions. The constructor does not recompute the sum or validate components.

# Returns
`SpineBoundary` with `components` the supplied component records and
`euler_sum` the supplied total Euler characteristic.
Component flag meanings are in
[Boundary diagnostic flags](../specs/ograph-colours.md#boundary-diagnostic-flags).

# Guarantees
1. Fields equal the supplied values after conversion; an input vector of
   the declared type is retained without copying. [proved/agent: by construction]
2. For data computed by [`thickening_boundary`](@ref), this describes the
   corner-link boundary model of the selected polyhedron. Its interpretation
   as the boundary of a regular neighbourhood is subject to the
   [geometric interpretation](../specs/ograph-colours.md#face-cycles).
   [proved/agent: by construction]

# Not guaranteed
Direct construction does not establish a surface, orientability, an Euler
sum, or a closed normal spine. It does not certify a closed 3-manifold.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — dual ideal triangulations; corner links as a boundary model are reconstructed.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> b = thickening_boundary(g); stored = SpineBoundary(b.components, b.euler_sum);

julia> (length(stored.components), stored.euler_sum)
(1, 2)
```

See also [`thickening_boundary`](@ref).

Role: convert
"""
struct SpineBoundary
    components::Vector{SpineBoundaryComponent}
    euler_sum::Int
end

"""
    DiscPartitionComparison(equal, ours, pipeline, germ_count, unique_germs) -> DiscPartitionComparison

Store a comparison of two disc partitions in the shared corner coordinates.

# Arguments
- `equal::Bool`, `ours::Vector{Vector{NTuple{3,Int}}}`,
  `pipeline::Vector{Vector{NTuple{3,Int}}}`, `germ_count::Int`,
  `unique_germs::Int`, or values convertible to those types.
  Conversion failures propagate their Julia exceptions.
- Obtain a computed comparison from [`compare_disc_partition`](@ref).
  This constructor does not recompute any field.

# Returns
`DiscPartitionComparison` with `equal` the supplied comparison flag,
`ours` the face-cycle partition, `pipeline` the disk-table partition,
`germ_count` the supplied total corner count and `unique_germs` the supplied
distinct-corner count for `ours`. The comparison's meaning is specified in
[Face cycles](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. Fields preserve the supplied values after conversion. [proved/agent: by construction]
2. For data returned by `compare_disc_partition`, the record compares disc
   partitions of the selected face-cycle model and the normal disk-table
   model on the same presentation; it is a finite comparison of those models.
   [proved/agent: by construction]
3. Input vectors of the declared types are retained without copying. [proved/agent: by construction]

# Not guaranteed
Direct construction does not establish equality, sorting, counts, or that
the two partitions refer to the same presentation. Agreement is not a proof
of the general geometric reconstruction.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> link = vl_example("s3_abalone"); c = normal_coloured_ograph(ograph_state(link));

julia> p = compare_disc_partition(link, c);

julia> q = DiscPartitionComparison(p.equal, p.ours, p.pipeline, p.germ_count, p.unique_germs);

julia> (q.equal, q.germ_count, q.unique_germs)
(true, 6, 6)
```

See also [`compare_disc_partition`](@ref).

Role: construct
"""
struct DiscPartitionComparison
    equal::Bool
    ours::Vector{Vector{NTuple{3,Int}}}
    pipeline::Vector{Vector{NTuple{3,Int}}}
    germ_count::Int
    unique_germs::Int
end

function _face_owner(graph::ColouredOGraph)
    owner = Dict{Int,Tuple{Int,Int}}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for (slot, halfedge) in enumerate(cyclic)
            owner[Int(halfedge)] = (vertex, slot)
        end
    end
    owner
end

function _face_colours(graph::ColouredOGraph)
    colour = Dict{Tuple{Int,Int},Int}()
    for record in graph.edge_colours
        other = Int(graph.mate[record.halfedge])
        colour[minmax(Int(record.halfedge), other)] = _colour_mod3(record.colour)
    end
    colour
end

function _face_edges(graph::ColouredOGraph)
    ids = zeros(Int, length(graph.mate))
    edges = Tuple{Int,Int}[]
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        push!(edges, (halfedge, Int(other)))
        ids[halfedge] = length(edges)
        ids[Int(other)] = length(edges)
    end
    edges, ids
end

"""Leg chart. `:a` puts leg 1 on the counter-clockwise side of an over slot.
`:b` swaps legs 1 and 3. Leg 2 stays the opposite half-edge.
"""
function _leg_chart(graph::ColouredOGraph, calibration::Symbol)
    calibration in (:a, :b) || throw(ArgumentError(
        "leg calibration must be :a or :b"))
    charts = Vector{NTuple{3,Int}}(undef, length(graph.mate))
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        over = graph.over_positions[vertex]
        for (slot, halfedge) in enumerate(cyclic)
            opposite = Int(cyclic[mod1(slot + 2, 4)])
            ccw = Int(cyclic[mod1(slot + 1, 4)])
            cw = Int(cyclic[mod1(slot - 1, 4)])
            upright = slot in over
            adjacent = upright ? (ccw, cw) : (cw, ccw)
            calibration == :b && (adjacent = (adjacent[2], adjacent[1]))
            charts[Int(halfedge)] = (adjacent[1], opposite, adjacent[2])
        end
    end
    charts
end

"Transport of one sheet across an edge: ρ∘τ^α, with the stored colour."
function _sheet_transport(colour::Integer, leg::Integer)
    fig25_frame_reversal(fig25_tau(Int(colour), leg))
end

function _corner(owner, halfedge::Int, partner::Int)
    vertex, slot_a = owner[halfedge]
    other_vertex, slot_b = owner[partner]
    other_vertex == vertex || throw(ArgumentError(
        "a vertex corner joins half-edges at two vertices"))
    lo, hi = minmax(slot_a, slot_b)
    (vertex, lo, hi)
end

"""
    face_cycles(graph::ColouredOGraph; calibration::Symbol=:b) -> ColouredPolyhedron

Return the disc boundaries and cell counts of the selected coloured-o-graph polyhedron.

# Arguments
- `graph` must satisfy O0–O3 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers).
  Failed checks throw `ArgumentError` identifying them.
- `calibration::Symbol` must be `:b` or `:a`; otherwise `ArgumentError` names
  those choices. Use the [adopted calibration](../conventions/OGraphColours.md)
  for the normal-spine interpretation; `:a` is the documented comparison.
  A reused or missing sheet end or an inconsistent corner or coverage
  count throws `ArgumentError` identifying that failure.

# Returns
`ColouredPolyhedron`: `vertex_count`, `edge_count`, `disc_count`, `euler`
are the cell counts and Euler characteristic; `halfedge_count`, `germ_count`
and `boundary_length` are the half-edge, counted-germ and total boundary-length
counts. `cycles::Vector{FaceCycle}` gives the disc boundaries; `calibration`
records the selected symbol and `status` is `OGRAPH_FACE_CYCLE_STATUS`.
See [Face cycles](../specs/ograph-colours.md#face-cycles) for the geometric objects.

# Guarantees
1. For `n = length(graph.vertex_halfedges)`, the result has `vertex_count == n`,
   `edge_count == 2n`, `halfedge_count == 4n`, and
   `germ_count == boundary_length == 6n`. [proved/agent: by construction]
2. `disc_count == length(cycles)` and `euler == vertex_count - edge_count + disc_count`.
   Each cycle's three sequences have equal lengths; their total length is `6n`.
   [proved/agent: by construction]
3. The returned corners partition the vertex germs of the selected gluing.
   The cycles represent its attached disc boundaries. Identification with
   the BP97 standard spine is under the reconstruction assumptions on the
   [convention page](../conventions/OGraphColours.md).
   [proved/agent: by construction]
4. The input graph is unchanged. `calibration` records the requested choice;
   `status` is a library evidence label, including when `calibration=:a`.
   [proved/agent: by construction]

# Conventions
Use the leg calibration, corner coordinates and signed-edge convention linked
above and in [the invariant specification](../specs/ograph-colour-invariants.md#closed-examples).

# Not guaranteed
O0–O3 do not establish closed-normal conditions, a single spherical boundary,
or `disc_count == n + 1`. Germ coverage alone does not prove the geometric
reconstruction. Calibration `:a` is not promised to recover the intended
normal spine, and the status string is not an input-specific certificate.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> p = face_cycles(g);

julia> (p.vertex_count, p.edge_count, p.disc_count, p.euler, p.germ_count, Tuple(sort([length(c.germs) for c in p.cycles])))
(1, 2, 2, 1, 6, (1, 5))
```

See also [`cellular_homology`](@ref), [`compare_disc_partition`](@ref), [`thickening_boundary`](@ref).

Role: convert
"""
function face_cycles(graph::ColouredOGraph; calibration::Symbol=:b)
    _require_coloured_ograph_shape(graph)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    _, edge_id = _face_edges(graph)
    n = length(graph.mate)
    edge_count = count(i -> i < graph.mate[i], 1:n)
    used = falses(n, 3)
    cycles = FaceCycle[]
    for start_halfedge in 1:n, start_leg in 1:3
        used[start_halfedge, start_leg] && continue
        germs = DiscGerm[]
        corners = NTuple{3,Int}[]
        letters = Int[]
        halfedge, leg = start_halfedge, start_leg
        while true
            used[halfedge, leg] && throw(ArgumentError(
                "sheet end ($halfedge, $leg) is reused"))
            used[halfedge, leg] = true
            other = Int(graph.mate[halfedge])
            alpha = get(colours, minmax(halfedge, other), nothing)
            alpha === nothing && throw(ArgumentError(
                "edge ($halfedge, $other) has no colour"))
            other_leg = _sheet_transport(alpha, leg)
            used[other, other_leg] && throw(ArgumentError(
                "arrival sheet end ($other, $other_leg) is reused"))
            used[other, other_leg] = true
            partner = charts[other][other_leg]
            back = findfirst(==(other), charts[partner])
            back === nothing && throw(ArgumentError(
                "half-edge $partner has no leg returning to $other"))
            push!(germs, DiscGerm(halfedge, leg))
            push!(corners, _corner(owner, other, partner))
            letter = edge_id[halfedge]
            push!(letters, halfedge < other ? letter : -letter)
            if (partner, back) == (start_halfedge, start_leg)
                break
            end
            halfedge, leg = partner, back
        end
        push!(cycles, FaceCycle(germs, corners, letters))
    end
    all(used) || throw(ArgumentError("a sheet end was not used"))
    vertex_count = length(graph.vertex_halfedges)
    disc_count = length(cycles)
    germ_count = sum(length, (cycle.germs for cycle in cycles))
    germ_count == 6 * vertex_count || throw(ArgumentError(
        "face cycles cover $germ_count germs, not $(6 * vertex_count)"))
    ColouredPolyhedron(
        vertex_count, edge_count, disc_count,
        vertex_count - edge_count + disc_count,
        n, germ_count, germ_count, cycles, calibration, OGRAPH_FACE_CYCLE_STATUS)
end

function _require_coloured_ograph_shape(graph::ColouredOGraph)
    failed = filter(check -> !check.passed, validate_coloured_ograph(graph))
    isempty(failed) || throw(ArgumentError(
        "face cycles require O0-O3: " * join((check.message for check in failed), "; ")))
    graph
end

function _spanning_tree(edges::Vector{Tuple{Int,Int}}, owner, vertex_count::Int)
    parent = collect(1:vertex_count)
    function find(vertex)
        while parent[vertex] != vertex
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        end
        vertex
    end
    tree = Int[]
    for (index, (left, right)) in enumerate(edges)
        a = find(owner[left][1])
        b = find(owner[right][1])
        a == b && continue
        parent[b] = a
        push!(tree, index)
    end
    length(tree) == vertex_count - 1 || throw(ArgumentError(
        "the 1-skeleton is not connected"))
    tree
end

"""
    cellular_homology(graph::ColouredOGraph; calibration::Symbol=:b) -> ColourHomology

Return the first integral homology of the selected face-cycle polyhedron.

# Arguments
- `graph` must satisfy O0–O3 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers).
  A failed condition throws `ArgumentError` identifying the failing checks.
- `calibration` selects the [leg calibration](../conventions/OGraphColours.md).
  It must be `:b` or `:a`; another symbol throws `ArgumentError` naming those
  choices. `:b` is the default; `:a` is retained for the documented comparison.
  An inconsistent face-cycle or disconnected-skeleton computation throws
  `ArgumentError` naming the failed consistency check.

# Returns
`ColourHomology`: `free_rank::Int` is the free abelian rank and
`torsion_factors::Vector{BigInt}` contains the nontrivial cyclic-summand orders,
using the [invariant-factor convention](../specs/ograph-colour-invariants.md#closed-examples).

# Guarantees
1. `free_rank >= 0`. Every torsion factor exceeds one and divides its
   successor; the factors and free rank specify the computed abelian group.
   [proved/agent: by construction]
2. This is cellular first homology with integer coefficients for the disc
   boundaries returned by [`face_cycles`](@ref), with the signed-edge
   convention linked above. The input graph is unchanged. [proved/agent: by construction]
3. Under the [geometric reconstruction](../specs/ograph-colour-invariants.md#closed-examples),
   this is also the first homology of the represented spine's thickening.
   When that thickening has one spherical boundary component, capping it
   with a ball leaves this first homology unchanged.
   [proved/agent: ../specs/ograph-colour-invariants.md#closed-examples, conditional on the stated reconstruction]

# Conventions
Use the selected calibration and the signed-edge and invariant-factor
conventions linked above.

# Not guaranteed
Identification with the BP97 spine retains the calibration dependencies on
the convention page. The `:a` comparison is not promised to recover the
intended normal spine. This function does not certify closedness, identify
a manifold, compute its fundamental group, or establish the general
raw-edge-presentation identification discussed in [`ograph_h1`](@ref).

# Sources
cf:[[BP97#^bp97-sec-0703]] §10.1 — cellular homology of a spine and its thickening, and the additional 3-cell in the closed case; invariant factors use the repository convention.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> h = cellular_homology(g);

julia> (h.free_rank, Tuple(h.torsion_factors))
(0, ())
```

See also [`face_cycles`](@ref), [`thickening_boundary`](@ref).

Role: invariant
"""
function cellular_homology(graph::ColouredOGraph; calibration::Symbol=:b)
    poly = face_cycles(graph; calibration)
    owner = _face_owner(graph)
    edges, _ = _face_edges(graph)
    tree = _spanning_tree(edges, owner, poly.vertex_count)
    generator_of = zeros(Int, poly.edge_count)
    generator_count = 0
    tree_set = Set(tree)
    for edge in 1:poly.edge_count
        edge in tree_set && continue
        generator_count += 1
        generator_of[edge] = generator_count
    end
    rows = length(poly.cycles)
    matrix = zeros(BigInt, rows, generator_count)
    for (row, cycle) in enumerate(poly.cycles)
        for letter in cycle.edge_letters
            generator = generator_of[abs(letter)]
            generator == 0 && continue
            matrix[row, generator] += sign(letter)
        end
    end
    diagonal = BigInt[]
    if rows > 0 && generator_count > 0
        smith = AbstractAlgebra.snf(AbstractAlgebra.matrix(AbstractAlgebra.ZZ, matrix))
        diagonal = [abs(BigInt(smith[index, index])) for index in 1:min(size(smith)...)]
    end
    nonzero = filter(!iszero, diagonal)
    all(value > 0 for value in nonzero) || throw(ArgumentError(
        "Smith diagonal must be nonnegative"))
    all(nonzero[index + 1] % nonzero[index] == 0 for index in 1:max(length(nonzero) - 1, 0)) ||
        throw(ArgumentError("Smith invariant factors are not divisibility ordered"))
    free_rank = generator_count - length(nonzero)
    free_rank >= 0 || throw(ArgumentError("Smith rank exceeds the generator count"))
    factors = BigInt[value for value in nonzero if value > 1]
    ColourHomology(free_rank, factors)
end

function _perm_sign(values)
    sign = 1
    sequence = collect(Int.(values))
    for left in 1:length(sequence) - 1, right in left + 1:length(sequence)
        sequence[left] == sequence[right] && return 0
        sequence[left] > sequence[right] && (sign = -sign)
    end
    sign
end

"""Induced orientation of the face opposite slot `omitted`.

The positive orientation of the tetrahedron is the cyclic order of the four
slots, begun at either slot of the over pair. The two over slots are
opposite, so the two beginnings differ by a double transposition and are the
same orientation. In that order, the face opposite the i-th slot carries
(−1)^i times the other three slots. The result is three slot indices. It
does not read half-edge ids.
"""
function _positive_face_order(over, omitted::Int)
    slots = Tuple(sort(Int.(over)))
    length(slots) == 2 && slots[2] - slots[1] == 2 ||
        throw(ArgumentError("the over pair must be opposite slots"))
    positive = ntuple(step -> mod1(slots[1] + step - 1, 4), 4)
    index = findfirst(==(omitted), positive)
    index === nothing && throw(ArgumentError("slot $omitted is outside 1:4"))
    others = [positive[step] for step in 1:4 if step != index]
    isodd(index) ? (others[1], others[2], others[3]) : (others[1], others[3], others[2])
end

function _order_sign(reference::NTuple{3,Int}, image::NTuple{3,Int})
    position = Dict(reference[index] => index for index in 1:3)
    _perm_sign(ntuple(index -> position[image[index]], 3))
end

"""Face identifications of the dual ideal triangulation.

Each record sends the three slots of the face opposite `source_slot` to the
three slots of the face opposite `target_slot`. The reversal flag compares
the induced orders of those faces. The half-edge comparison only selects one
end of the undirected edge. `dual_face_pairings` lists the smaller
`(vertex, slot)` first, and exchanging the ends keeps the flag.
"""
function _slot_gluings(graph::ColouredOGraph, calibration::Symbol)
    _require_coloured_ograph_shape(graph)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    gluings = NamedTuple[]
    seen = Set{Tuple{Int,Int}}()
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        source_vertex, source_slot = owner[halfedge]
        target_vertex, target_slot = owner[Int(other)]
        alpha = colours[minmax(halfedge, Int(other))]
        images = Dict{Int,Int}()
        for leg in 1:3
            neighbour = charts[halfedge][leg]
            _, neighbour_slot = owner[neighbour]
            image_leg = _sheet_transport(alpha, leg)
            image_halfedge = charts[Int(other)][image_leg]
            _, image_slot = owner[image_halfedge]
            images[neighbour_slot] = image_slot
        end
        source_positive = _positive_face_order(graph.over_positions[source_vertex], source_slot)
        image_positive = ntuple(index -> images[source_positive[index]], 3)
        target_positive = _positive_face_order(graph.over_positions[target_vertex], target_slot)
        target_vertices = Set(slot for slot in 1:4 if slot != target_slot)
        Set(values(images)) == target_vertices && length(images) == 3 ||
            throw(ArgumentError(
                "face ($source_vertex, $source_slot) does not match three slots"))
        reversing = _order_sign(target_positive, image_positive) == -1
        push!(gluings, (; source_vertex, source_slot, target_vertex, target_slot,
            images, reversing))
        push!(seen, (source_vertex, source_slot))
        push!(seen, (target_vertex, target_slot))
    end
    vertex_count = length(graph.vertex_halfedges)
    expected = Set((vertex, slot) for vertex in 1:vertex_count for slot in 1:4)
    seen == expected || throw(ArgumentError("a tetrahedron face is unpaired"))
    length(gluings) == 2 * vertex_count || throw(ArgumentError(
        "a quadrivalent graph has 2n face pairs"))
    gluings
end

"""
    dual_face_pairings(graph::ColouredOGraph; calibration::Symbol=:b) -> Vector{DualFacePairing}

Return the paired faces and orientation flags of the selected dual slot model.

# Arguments
- `graph` must satisfy O0–O3 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers);
  failed conditions throw `ArgumentError` identifying the failing checks.
- `calibration::Symbol` must be `:b` or `:a`; otherwise `ArgumentError` names
  the choices. The [calibration](../conventions/OGraphColours.md) selects the
  geometric interpretation. Inconsistent face coverage or slot matching
  throws `ArgumentError` naming that consistency failure.

# Returns
`Vector{DualFacePairing}`. Each record gives `source_vertex`, `source_slot`,
`target_vertex`, `target_slot` and `orientation_reversing` with the meanings
in [Face cycles](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. For `n` vertices there are `2n` records; every face label occurs once,
   all vertex labels lie in `1:n`, and all slots lie in `1:4`.
   [proved/agent: by construction]
2. Endpoints and reversal flags obey the
   [face-orientation convention](../conventions/OGraphColours.md#四面体と-link-triangle-の向き).
   Changing which end of an edge is considered first does not change its flag.
   [proved/agent: by construction]
3. Under the [geometric interpretation](../specs/ograph-colour-invariants.md#dual-ideal-triangulation),
   these are face pairs of the dual ideal triangulation of the selected
   polyhedron. The input graph is unchanged. [proved/agent: by construction]

# Not guaranteed
Successful return does not imply that every reversal flag is true. O0–O3
alone are not a certificate of the desired oriented triangulation. The
record omits the face's vertex bijection, and the specification establishes
no isomorphism to the separate `disk_tetrahedra` certificate.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — ideal triangulations use orientation-reversing face identifications; the repository slot reconstruction supplies the flags.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> pairs = dual_face_pairings(g);

julia> (length(pairs), all(p.orientation_reversing for p in pairs))
(2, true)
```

See also [`dual_edge_orbits`](@ref), [`thickening_boundary`](@ref).

Role: access
"""
function dual_face_pairings(graph::ColouredOGraph; calibration::Symbol=:b)
    rows = DualFacePairing[]
    for gluing in _slot_gluings(graph, calibration)
        source = (gluing.source_vertex, gluing.source_slot)
        target = (gluing.target_vertex, gluing.target_slot)
        # The reversal flag is a property of the unordered pair. The record
        # lists the end with the smaller (vertex, slot) first.
        if target < source
            source, target = target, source
        end
        push!(rows, DualFacePairing(source[1], source[2], target[1], target[2],
            gluing.reversing))
    end
    rows
end

"""
    dual_edge_orbits(graph::ColouredOGraph; calibration::Symbol=:b) -> Vector{Vector{NTuple{3,Int}}}

Return the corner-germ classes representing dual edges of the selected polyhedron.

# Arguments
- `graph` must satisfy O0–O3 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers).
  Call [`validate_coloured_ograph`](@ref) first on untrusted data: this
  function does not run that validator. Malformed data can throw `KeyError`,
  `BoundsError` or `ArgumentError` instead of a condition-labelled rejection.
- `calibration::Symbol` must be `:b` or `:a`; an unsupported symbol throws
  `ArgumentError` naming these choices. Use the
  [adopted calibration](../conventions/OGraphColours.md) for the normal-spine interpretation.

# Returns
A vector of classes, each a vector of corner triples in the
[face-cycle coordinates](../specs/ograph-colours.md#face-cycles).
Neither class order nor order within a class is part of the contract.

# Guarantees
1. For `n` vertices the classes partition all `6n` corner germs: each occurs
   in exactly one class, without repetition. [proved/agent: by construction]
2. The equivalence classes are those of the selected sheet identifications.
   Under the [dual ideal triangulation interpretation](../specs/ograph-colour-invariants.md#dual-ideal-triangulation),
   they represent its edges, one edge per class. [proved/agent: by construction]
3. `sorted_germ_partition` of this result equals that of the corners of
   `face_cycles(graph; calibration)` for the same valid graph and calibration.
   The input graph is unchanged. [proved/agent: by construction]

# Conventions
Use the calibration and corner coordinates linked above.

# Not guaranteed
The result does not certify an ideal triangulation, its orientation, or its
identification with the separate `disk_tetrahedra` certificate. The geometric
interpretation retains the reconstruction dependencies on the convention page.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — dual ideal triangulations; edge classes from the repository slot gluings are reconstructed.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> orbits = dual_edge_orbits(g);

julia> (length(orbits), sum(length, orbits))
(2, 6)
```

See also [`dual_face_pairings`](@ref), [`face_cycles`](@ref), [`sorted_germ_partition`](@ref).

Role: access
"""
function dual_edge_orbits(graph::ColouredOGraph; calibration::Symbol=:b)
    charts = _leg_chart(graph, calibration)
    owner = _face_owner(graph)
    colours = _face_colours(graph)
    parent = Dict{NTuple{3,Int},NTuple{3,Int}}()
    function find(germ)
        while parent[germ] != germ
            parent[germ] = parent[parent[germ]]
            germ = parent[germ]
        end
        germ
    end
    function unite(left, right)
        left_root = find(left)
        right_root = find(right)
        left_root == right_root || (parent[right_root] = left_root)
    end
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for left in 1:3, right in left + 1:4
            germ = (vertex, left, right)
            parent[germ] = germ
        end
    end
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        alpha = colours[minmax(halfedge, Int(other))]
        for leg in 1:3
            tail_partner = charts[halfedge][leg]
            other_leg = _sheet_transport(alpha, leg)
            head_partner = charts[Int(other)][other_leg]
            unite(_corner(owner, halfedge, tail_partner),
                  _corner(owner, Int(other), head_partner))
        end
    end
    orbits = Dict{NTuple{3,Int},Vector{NTuple{3,Int}}}()
    for germ in keys(parent)
        push!(get!(orbits, find(germ), NTuple{3,Int}[]), germ)
    end
    collect(values(orbits))
end

"""
    sorted_germ_partition(discs) -> Vector{Vector{NTuple{3,Int}}}

Return a sorted comparison form of a collection of disc-corner lists.

# Arguments
`discs` is an iterable of iterable, sortable corner triples convertible to
`NTuple{3,Int}` in the [face-cycle coordinates](../specs/ograph-colours.md#face-cycles).
Incompatible elements propagate sorting or conversion errors. No validity
check of discs, corners or a partition is performed.

# Returns
A fresh vector of sorted corner vectors, itself sorted lexicographically.

# Guarantees
1. The result represents the same multiset of corner multisets, retaining
   repeated corners and repeated lists. [proved/agent: by construction]
2. It forgets the order of lists and their internal cyclic order, so equality
   compares the corner partitions in the linked face-cycle comparison.
   The supplied lists are unchanged. [proved/agent: by construction]

# Not guaranteed
Equality does not identify cyclic attaching maps, orientations, coloured
o-graphs or their thickenings. Arbitrary input lists are not certified as
the regions of a special polyhedron.

# Sources
cf:docs/specs/ograph-colours.md#face-cycles — comparison of disc-corner partitions.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); c = normal_coloured_ograph(v); p = face_cycles(c.graph);

julia> sorted_germ_partition(cycle.corners for cycle in p.cycles) == sorted_germ_partition(disk_table_corners(v, c.state))
true
```

See also [`compare_disc_partition`](@ref), [`disk_table_corners`](@ref).

Role: convert
"""
function sorted_germ_partition(discs)
    parts = Vector{NTuple{3,Int}}[sort!(collect(disc)) for disc in discs]
    sort!(parts)
    parts
end

function _port_order(sign::Integer)
    sign == 1 && return (:mi, :po, :mo, :pi)
    sign == -1 && return (:mi, :pi, :mo, :po)
    throw(ArgumentError("vertex sign must be +1 or -1"))
end

const _DISK_PORT_LANES = Dict(
    :mi => (1, 2, 3),
    :pi => (4, 5, 6),
    :mo => (7, 8, 9),
    :po => (10, 11, 12),
)

"""
    disk_table_corners(link::VirtualLink, state::OGraphState) -> Vector{Vector{NTuple{3,Int}}}

Return the disk-table circuits of a normal presentation in vertex-corner coordinates.

# Arguments
- `link`, `state`: the caller supplies the same labelled closed normal
  presentation, with `state = ograph_state(link)`. The correspondence is not
  checked; mismatches can give incorrect corners or propagate `KeyError`,
  `BoundsError` or `ArgumentError` from the disk and corner data.
  Use the [normal adapter](../specs/ograph-colours.md#normal-adapter) and
  [face-cycle comparison](../specs/ograph-colours.md#face-cycles).

# Returns
One vector of `(vertex, slot_lo, slot_hi)` corners per disk-table path, with
path order and multiplicity retained, in the linked corner coordinates.

# Guarantees
1. Each returned vector has the length of the corresponding `disk_table(link).paths` entry. [proved/agent: by construction]
2. Under the matching-presentation precondition, the result expresses the same disk-table circuit partition in the normal adapter's corner coordinates. [proved/agent: by construction]
3. Neither argument is modified. [proved/agent: by construction]

# Conventions
Use the [normal adapter and slot transcription](../conventions/OGraphColours.md).

# Not guaranteed
This conversion alone neither verifies C1–C3 nor identifies the circuit
partition with a spine; the dependencies of that identification are recorded
in [Slot transcription](../specs/ograph-colours.md#slot-transcription).

# Sources
cf:docs/specs/ograph-colours.md#face-cycles — corner coordinates and comparison of the two circuit partitions.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone");

julia> sum(length, disk_table_corners(v, ograph_state(v)))
6
```

See also [`compare_disc_partition`](@ref).

Role: convert
"""
function disk_table_corners(link::VirtualLink, state::OGraphState)
    table = disk_table(link)
    dictionary = _disk_node_dictionary()
    slot_of_port = Dict{Tuple{Int,Symbol},Int}()
    node_ports = Dict{Tuple{Int,Int},NTuple{2,Symbol}}()
    for (vertex, sign) in enumerate(state.vertex_signs)
        for (slot, port) in enumerate(_port_order(sign))
            slot_of_port[(vertex, port)] = slot
        end
        slots = dictionary[sign]
        owners = Dict{Int,Vector{Symbol}}()
        for (port, indices) in _DISK_PORT_LANES
            for index in indices
                push!(get!(owners, slots[index], Symbol[]), port)
            end
        end
        for (node, ports) in owners
            unique_ports = unique(ports)
            length(unique_ports) == 2 || throw(ArgumentError(
                "disk node $node at vertex $vertex meets ports $ports"))
            a, b = sort(collect(unique_ports))
            node_ports[(vertex, node)] = (a, b)
        end
    end
    discs = Vector{NTuple{3,Int}}[]
    for path in table.paths
        germs = NTuple{3,Int}[]
        for node in path
            vertex = (node - 1) ÷ 6 + 1
            local_node = mod1(node, 6)
            port_a, port_b = node_ports[(vertex, local_node)]
            slot_a = slot_of_port[(vertex, port_a)]
            slot_b = slot_of_port[(vertex, port_b)]
            push!(germs, (vertex, min(slot_a, slot_b), max(slot_a, slot_b)))
        end
        push!(discs, germs)
    end
    discs
end

"""
    compare_disc_partition(link::VirtualLink, colouring::ClosedNormalColouring; calibration::Symbol=:b) -> DiscPartitionComparison

Return a comparison record of the face-cycle and disk-table disc partitions of a normal presentation.

# Arguments
- `link` and `colouring` must describe the same labelled normal presentation;
  obtain `colouring` with `normal_coloured_ograph(ograph_state(link))`.
  This correspondence is a caller precondition and is not checked here.
- `colouring.graph` must satisfy O0–O3 and `calibration` must be `:b` or `:a`,
  as required by [`face_cycles`](@ref). Its `ArgumentError` failures propagate.
- The disk-table and slot conversion must be defined for `link` and
  `colouring.state`. Inconsistent manually assembled inputs can propagate
  `KeyError`, `BoundsError` or `ArgumentError`; there is no dedicated
  mismatched-presentation exception.

# Returns
`DiscPartitionComparison` with `ours` and `pipeline` the two sorted corner
partitions, `equal::Bool` their equality, `germ_count::Int` the total number
of corners in `ours`, and `unique_germs::Int` its number of distinct corners.
See the [partition comparison](../specs/ograph-colours.md#face-cycles).

# Guarantees
1. `equal == (ours == pipeline)`. The comparison retains corner multiplicities
   and disregards the cyclic order of corners within a disc. [proved/agent: by construction]
2. The counts refer to `ours`, even when `equal` is false. [proved/agent: by construction]
3. On corresponding inputs the result compares two combinatorial models
   of the normal spine's discs under the selected calibration. It is
   evidence about this input, with the scope specified on the linked page.
   [proved/agent: by construction]

# Conventions
The corner coordinates and calibration are those of
[Face cycles](../specs/ograph-colours.md#face-cycles) and
[OGraphColours](../conventions/OGraphColours.md).

# Not guaranteed
Partition equality does not prove the general BP97 reconstruction, identify
manifolds, or choose a calibration. No equality between the two presentations
is inferred from the supplied `ClosedNormalColouring` type alone.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> link = vl_example("s3_abalone"); c = normal_coloured_ograph(ograph_state(link));

julia> r = compare_disc_partition(link, c);

julia> (r.equal, r.germ_count, r.unique_germs)
(true, 6, 6)
```

See also [`face_cycles`](@ref), [`disk_table_corners`](@ref).

Role: condition
"""
function compare_disc_partition(link::VirtualLink, colouring::ClosedNormalColouring;
                                calibration::Symbol=:b)
    poly = face_cycles(colouring.graph; calibration)
    ours = [collect(cycle.corners) for cycle in poly.cycles]
    pipeline = disk_table_corners(link, colouring.state)
    left = sorted_germ_partition(ours)
    right = sorted_germ_partition(pipeline)
    flat = collect(Iterators.flatten(left))
    DiscPartitionComparison(left == right, left, right, length(flat), length(unique(flat)))
end

"""
    thickening_boundary(graph::ColouredOGraph; calibration::Symbol=:b) -> SpineBoundary

Return the component counts and Euler characteristics of the selected corner-link boundary model.

# Arguments
- `graph` must satisfy O0–O3 in the
  [condition table](../specs/ograph-colours.md#condition-identifiers).
  Failed conditions throw `ArgumentError` identifying the failing checks.
- `calibration::Symbol` must be `:b` or `:a`; otherwise `ArgumentError` names
  the choices. Inconsistent face or side pairings, corner membership or
  boundary walks throw `ArgumentError` naming the failed check.
- For the oriented regular-neighbourhood interpretation, require reversing
  dual face pairs under the [face-orientation convention](../conventions/OGraphColours.md#四面体と-link-triangle-の向き)
  and the geometric reconstruction assumptions linked below. Nonreversing
  pairs are reported by the component flags, not rejected by this function.

# Returns
`SpineBoundary` with `components::Vector{SpineBoundaryComponent}` and
`euler_sum::Int`. Each component records `euler`, `face_count`, `edge_count`,
`vertex_count`, `orientable` and `closed_surface`; the counts refer to its
corner-link cellulation. The flags are diagnostics for that model.
Their exact scope is specified in
[Boundary diagnostic flags](../specs/ograph-colours.md#boundary-diagnostic-flags).

# Guarantees
1. `euler_sum` equals the sum of the component Euler characteristics, and
   each `euler == vertex_count - edge_count + face_count`.
   For `n` graph vertices, the component face counts sum to `4n` and the
   component edge counts sum to `6n`. [proved/agent: by construction]
2. Every returned component has `closed_surface == true`. If every dual
   face pair reverses orientation, every component has `orientable == true`.
   [proved/agent: by construction]
3. Under the [reconstruction interpretation](../specs/ograph-colours.md#face-cycles),
   these components are the boundary of a regular neighbourhood of the
   selected polyhedron. They need not form a single sphere.
   [proved/agent: by construction]
4. Components are sorted by `(euler, face_count, edge_count, vertex_count)`.
   The input graph is unchanged. [proved/agent: by construction]

# Conventions
Use the selected leg calibration and the corner-link orientations on the
convention page linked above. The separate ribbon construction is owned by
[the boundary specification](../specs/ograph-colour-boundary.md).

# Not guaranteed
O0–O3 alone do not certify the geometric interpretation or the Euler relation
between this boundary and the face-cycle polyhedron; the
[documented move failures](../specs/ograph-colour-invariants.md#c-move-and-mp-move)
exhibit that distinction. A false `orientable` flag is not a complete
nonorientability certificate. A spherical boundary alone does not certify
the closed-normal conditions or identify the capped 3-manifold.

# Sources
cf:[[BP97#^bp97-sec-0557]] §2.3 — the dual ideal triangulation; the corner-link boundary model is reconstructed.
cf:[[BP97#^bp97-sec-0553]] Proposition 2.2.1 — the boundary surface associated with an o-graph; this function uses the separate slot model, not the ribbon construction of that proposition.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = normal_coloured_ograph(ograph_state(vl_example("s3_abalone"))).graph;

julia> b = thickening_boundary(g); c = only(b.components);

julia> (b.euler_sum, c.euler, c.face_count, c.edge_count, c.vertex_count, c.orientable, c.closed_surface)
(2, 2, 4, 6, 4, true, true)
```

See also [`dual_face_pairings`](@ref), [`face_cycles`](@ref).

Role: invariant
"""
function thickening_boundary(graph::ColouredOGraph; calibration::Symbol=:b)
    gluings = _slot_gluings(graph, calibration)
    partner = Dict{NTuple{3,Int},NTuple{3,Int}}()
    reversing_of = Dict{NTuple{3,Int},Bool}()
    for gluing in gluings
        for (source_slot, target_slot) in gluing.images
            left = (gluing.source_vertex, source_slot, gluing.source_slot)
            right = (gluing.target_vertex, target_slot, gluing.target_slot)
            partner[left] = right
            partner[right] = left
            reversing_of[left] = gluing.reversing
            reversing_of[right] = gluing.reversing
        end
    end
    vertex_count = length(graph.vertex_halfedges)
    triangles = Vector{NTuple{3,Int}}[]
    for vertex in 1:vertex_count, slot in 1:4
        order = _positive_face_order(graph.over_positions[vertex], slot)
        push!(triangles, [(vertex, slot, order[index]) for index in 1:3])
    end
    next_of = Dict{NTuple{3,Int},NTuple{3,Int}}()
    owner = Dict{NTuple{3,Int},Int}()
    for (index, triangle) in enumerate(triangles)
        for position in 1:3
            side = triangle[position]
            haskey(partner, side) || throw(ArgumentError("link side $side is unpaired"))
            partner[partner[side]] == side ||
                throw(ArgumentError("link side pairing is not an involution"))
            next_of[side] = triangle[mod1(position + 1, 3)]
            owner[side] = index
        end
    end
    parent = collect(1:length(triangles))
    function find(index)
        while parent[index] != index
            parent[index] = parent[parent[index]]
            index = parent[index]
        end
        index
    end
    for side in keys(partner)
        haskey(owner, side) || throw(ArgumentError("link side $side is not on a corner"))
        left = find(owner[side])
        right = find(owner[partner[side]])
        left == right || (parent[right] = left)
    end
    seen = Set{NTuple{3,Int}}()
    cycle_of = Dict{NTuple{3,Int},Int}()
    cycle_count = 0
    for side in keys(next_of)
        side in seen && continue
        cycle_count += 1
        current = side
        while true
            current in seen && throw(ArgumentError("link vertex walk is not a cycle"))
            push!(seen, current)
            cycle_of[current] = cycle_count
            current = next_of[partner[current]]
            current == side && break
        end
    end
    groups = Dict{Int,Vector{Int}}()
    for index in 1:length(triangles)
        push!(get!(groups, find(index), Int[]), index)
    end
    components = SpineBoundaryComponent[]
    for indices in values(groups)
        sides = NTuple{3,Int}[]
        orientable = true
        for index in indices, side in triangles[index]
            push!(sides, side)
            reversing_of[side] || (orientable = false)
        end
        faces = length(indices)
        edge_count = length(sides) ÷ 2
        vertices = length(unique(cycle_of[side] for side in sides))
        push!(components, SpineBoundaryComponent(
            vertices - edge_count + faces, faces, edge_count, vertices, orientable, true))
    end
    sort!(components, by=component ->
        (component.euler, component.face_count, component.edge_count, component.vertex_count))
    SpineBoundary(components, sum(component.euler for component in components; init=0))
end
