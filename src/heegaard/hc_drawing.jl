# cf: lib/topology.jl/docs/specs/drawing-objects.md §3–§6, §8 item 1
# Coordinates are copied from plot(hc). Ordinary edges run between drawn events
# (alpha–beta crossings and ports on l_i). Reference and marked edges are the
# LC intervals from _hp_mark_regions, and only that path has faces.

const _hd_layout_calls = Ref(0)

"""
    DrawingPrimitive(id, shape, points, closed) -> DrawingPrimitive

Store one geometric part of a Heegaard drawing.

# Arguments
`id::Int`, `shape::Symbol`, `points::Tuple{Vararg{NTuple{2,Float64}}}`,
`closed::Bool`, or values convertible to those types. Conversion errors
propagate. No shape, coordinate-finiteness or identifier check is performed.

# Returns
`DrawingPrimitive` with the supplied part identifier, shape tag, point
coordinates and closure flag, for the
[drawing representation](../specs/drawing-objects.md#31-要素が持つもの).

# Guarantees
1. Fields preserve the supplied values after conversion; coordinates are
   stored as a tuple of two-coordinate tuples. Already-typed arguments are not modified. [proved/agent: by construction]
2. For parts returned by [`drawing`](@ref), the coordinates belong to that
   drawing's planar chart. They describe a displayed part rather than a
   complete topological curve on the glued surface. [proved/agent: by construction]

# Not guaranteed
The constructor does not place or validate a diagram, establish the meaning
of `shape`, or associate the part with a `DrawingElement`.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261111-plot-acceptance", "package-run-5", "inputs", "lens21-m", "diagram.hc");

julia> d = H.drawing(H.read_hc(path)); c = first(e for e in d.elements if e.kind == :crossing);

julia> p = d.primitives[first(c.primitives)]; q = H.DrawingPrimitive(p.id, p.shape, p.points, p.closed);

julia> length(q.points)
1
```

See also [`DrawingElement`](@ref), [`drawing`](@ref).

Role: construct
"""
struct DrawingPrimitive
    id::Int
    shape::Symbol
    points::Tuple{Vararg{NTuple{2,Float64}}}
    closed::Bool
end

"""
    DrawingElement(kind, key, input, drawn, primitives) -> DrawingElement

Store one semantic element and its part references in a Heegaard drawing.

# Arguments
`kind::Symbol`, an arbitrary `key`, `input::NamedTuple`, `drawn::NamedTuple`,
and `primitives::Tuple{Vararg{Int}}`, or values convertible to those field
types. Conversion errors propagate. References and named-tuple schemas are
not checked against a drawing.

# Returns
`DrawingElement` with the supplied element kind, semantic key, input-derived
attributes, drawing-derived attributes and primitive references. Their
meaning is owned by [drawing elements](../specs/drawing-objects.md#3-意味の要素).

# Guarantees
1. Fields preserve the supplied values after conversion. The constructor
   performs no layout or geometric computation. Already-typed arguments are not modified. [proved/agent: by construction]
2. For elements returned by [`drawing`](@ref), the element refers to the
   indicated feature of that input diagram, in the scope of the linked
   specification. One element can refer to several displayed parts. [proved/agent: by construction]

# Not guaranteed
The constructor does not ensure that the key exists in an input diagram,
that primitive references are valid, or that nested attribute values are
immutable. It does not certify a geometric embedding or a 3-manifold.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261111-plot-acceptance", "package-run-5", "inputs", "lens21-m", "diagram.hc");

julia> d = H.drawing(H.read_hc(path)); c = first(e for e in d.elements if e.kind == :crossing);

julia> q = H.DrawingElement(c.kind, c.key, c.input, c.drawn, c.primitives);

julia> H.crossing_curves(d, q.key) == (alpha=q.input.alpha, beta=q.input.beta)
true
```

See also [`DrawingPrimitive`](@ref), [`drawing`](@ref), [`element_at`](@ref).

Role: construct
"""
struct DrawingElement
    kind::Symbol
    key::Any
    input::NamedTuple
    drawn::NamedTuple
    primitives::Tuple{Vararg{Int}}
end

mutable struct DrawingStyleState
    rules::Vector{NamedTuple}
    listeners::Vector{Any}
end

"""
    HeegaardDrawing(primitives, elements, owner, bigons, path, layout_calls, style) -> HeegaardDrawing

Store a Heegaard drawing's geometric parts, semantic elements and shared style state.

# Arguments
`primitives::Tuple{Vararg{DrawingPrimitive}}`,
`elements::Tuple{Vararg{DrawingElement}}`, `owner::Tuple{Vararg{Int}}`,
`bigons::Tuple{Vararg{Tuple{Int,Int}}}`, `path::Symbol`, `layout_calls::Int`,
and `style::DrawingStyleState`, or convertible values. Conversion errors
propagate. The constructor does not check cross-references or drawing validity.

# Returns
`HeegaardDrawing` with geometric primitives, semantic elements, a primitive-to-
element `owner` map, bigon corner pairs, drawing-path tag, layout count and
mutable style state. Use [the drawing object](../specs/drawing-objects.md)
for their geometric scope.

# Guarantees
1. Fields preserve supplied values after conversion; already-typed arguments
   are retained and not modified. In particular, `style` is shared.
   [proved/agent: by construction]
2. When returned by [`drawing`](@ref), this is a planar rendering of the
   input diagram with semantic references to its curves and intersections.
   Construction from fields performs no layout.
   [proved/agent: by construction]

# Not guaranteed
Arbitrary fields need not describe a consistent drawing. Shared mutable style
or nested attributes are not frozen by the outer immutable struct. Planar
coordinates alone do not encode the full glued surface or a 3-manifold.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261111-plot-acceptance", "package-run-5", "inputs", "lens21-m", "diagram.hc");

julia> d = H.drawing(H.read_hc(path)); e = H.HeegaardDrawing(d.primitives, d.elements, d.owner, d.bigons, d.path, d.layout_calls, d.style);

julia> e.path
:reference
```

See also [`drawing`](@ref), [`element_at`](@ref), [`style!`](@ref).

Role: convert
"""
struct HeegaardDrawing
    primitives::Tuple{Vararg{DrawingPrimitive}}
    elements::Tuple{Vararg{DrawingElement}}
    owner::Tuple{Vararg{Int}}
    bigons::Tuple{Vararg{Tuple{Int,Int}}}
    path::Symbol
    layout_calls::Int
    style::DrawingStyleState
end

function Base.show(io::IO, d::HeegaardDrawing)
    print(io, "HeegaardDrawing(", d.path, ", ", length(d.elements), " elements, ", length(d.primitives), " parts)")
end

const _HD_ATTRS = (:visible, :color, :linewidth, :linestyle, :fill, :marker, :markersize, :label, :labelsize, :zorder)
const _HD_SELS = (:name, :role, :face, :sign, :bigon, :upper_alpha, :chart, :curve, :alpha, :beta)

"""
    drawingplot(d::HeegaardDrawing) -> Makie.Figure

Render the placed semantic drawing into a new Makie figure.

# Arguments
Load a Makie backend, such as `CairoMakie`, to activate `TopologyMakieExt`; `using Topology` alone supplies the generic binding without rendering methods. `d` must have valid primitive ownership, style and element records in [the drawing interface](../specs/drawing-objects.md#6-makie-との接続). Without the extension, calling the binding raises `MethodError`; backend and malformed-record errors propagate.

# Returns
A `Makie.Figure` containing a drawing axis and data inspector.

# Guarantees
1. The extension renders the drawing's existing primitives with resolved semantic styles and the specified screen orientation. It does not recompute mathematical layout. [proved/agent: by construction]
2. Primitive coordinates are retained; a listener is appended to the drawing's mutable style state for subsequent style updates. [proved/agent: by construction]

# Not guaranteed
This generic binding does not install or activate a backend. Rendering does not validate the source diagram or save a file; use Makie's display/save interface.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md §6 — optional Makie extension and rendering boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> parentmodule(H.drawingplot) === H
true
```

Role: draw
"""
function drawingplot end
"""
    drawingplot!(ax, d::HeegaardDrawing) -> Makie.Axis

Add a placed semantic drawing to an existing Makie axis.

# Arguments
The optional Makie extension must be loaded. `ax` is a `Makie.Axis`; `d` has `drawingplot`'s record and style requirements. Missing methods, incompatible axes and malformed drawing data raise the corresponding backend or access errors.

# Returns
The supplied axis after adding the drawing's plots.

# Guarantees
1. Existing drawing coordinates are rendered with resolved styles; the axis uses the screen orientation of [the Makie interface](../specs/drawing-objects.md#6-makie-との接続). [proved/agent: by construction]
2. The axis is mutated and a listener is appended to `d.style.listeners`; mathematical input and primitive coordinates are not recomputed. [proved/agent: by construction]

# Not guaranteed
Existing axis plots are not cleared. Repeated calls add plots and listeners. This operation neither verifies the mathematical input nor saves a file.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md §6 — rendering existing parts and style updates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> parentmodule(H.drawingplot!) === H
true
```

Role: draw
"""
function drawingplot! end

_hd_in(;role=:none, curve=nothing, sign=nothing, bigon=false, faces=(),
        upper_alpha=false, chart=:none, alpha=nothing, beta=nothing, boundary=()) =
    (;role, curve, sign, bigon, faces=Tuple(faces), upper_alpha, chart, alpha, beta, boundary=Tuple(boundary))

_hd_drawn(;angle=nothing, labelpoint=nothing) = (;angle, labelpoint)

_hd_tup(pts) = Tuple((Float64(p[1]), Float64(p[2])) for p in pts)

function _hd_add!(prims, owners, elements, kind, key, input, drawn, parts; shape=:polyline, closed=false)
    el = length(elements) + 1
    ids = Int[]
    for pts in parts
        isempty(pts) && continue
        t = _hd_tup(pts)
        length(t) < (shape == :point ? 1 : 2) && continue
        id = length(prims) + 1
        flag = closed && !(length(t) >= 2 && _hp_dist(t[1], t[end]) < 1e-8)
        push!(prims, DrawingPrimitive(id, shape, t, flag))
        push!(owners, el)
        push!(ids, id)
    end
    isempty(ids) && throw(ArgumentError("drawing element $kind $key has no part"))
    push!(elements, DrawingElement(kind, key, input, drawn, Tuple(ids)))
    elements[end]
end

function _hd_piece(r, a, b)
    m = r.closed ? length(r.points) : length(r.points) - 1
    !r.closed && b < a - 1e-9 && throw(ArgumentError("open piece runs backward on $(r.name)"))
    r.closed && b <= a + 1e-10 && (b += m)
    at(s) = begin
        !r.closed && s >= m - 1e-10 && return last(r.points)
        u = mod(s, m)
        i = floor(Int, u) + 1
        _hp_lerp(r.points[i], r.points[mod1(i + 1, length(r.points))], u - floor(u))
    end
    q = HCPlotPoint[at(a)]
    for k in floor(Int, a) + 1:ceil(Int, b) - 1
        _hp_dist(last(q), at(k)) > 1e-10 && push!(q, at(k))
    end
    _hp_dist(last(q), at(b)) > 1e-10 && push!(q, at(b))
    q
end

_hd_run_index(r) = begin
    i = findlast('#', r.name)
    i === nothing ? 0 : parse(Int, r.name[nextind(r.name, i):end])
end

function _hd_ordered_runs(p, contacts, component)
    runs = [r for r in p.runs if r.component == component && !startswith(r.name, "@K:")]
    sort!(runs; by=_hd_run_index)
    isempty(runs) && throw(ArgumentError("curve $component has no placed run"))
    n = length(runs)
    if n > 1 || !runs[1].closed
        for t in 1:n
            a = runs[t]
            b = runs[mod1(t + 1, n)]
            ca = get(contacts, (a.name, :end), nothing)
            cb = get(contacts, (b.name, :start), nothing)
            ca == cb || throw(ArgumentError("contacts do not chain $component: $(a.name) end=$ca, $(b.name) start=$cb"))
        end
    end
    runs
end

function _hd_span_pieces(runs, depart, arrive)
    i, j = depart.pos, arrive.pos
    r = runs[i]
    ahead = arrive.station > depart.station + 1e-9
    if i == j && (ahead || r.closed)
        return [_hd_piece(r, depart.station, arrive.station)]
    end
    n = length(runs)
    i == j && n == 1 && throw(ArgumentError("open run $(r.name) has a drawn edge against its direction"))
    seq = Int[i]
    k = i
    while true
        k = mod1(k + 1, n)
        push!(seq, k)
        k == j && break
        length(seq) > n + 1 && throw(ArgumentError("run cycle does not reach the next drawn event on $(r.component)"))
    end
    pieces = Vector{HCPlotPoint}[]
    for (t, pos) in enumerate(seq)
        rr = runs[pos]
        a = t == 1 ? depart.station : 0.0
        mend = Float64(rr.closed ? length(rr.points) : length(rr.points) - 1)
        b = t == length(seq) ? arrive.station : mend
        q = _hd_piece(rr, a, b)
        length(q) >= 2 && push!(pieces, q)
    end
    isempty(pieces) && throw(ArgumentError("drawn edge on $(r.component) has no polyline"))
    pieces
end

_hd_dedup(pts; tol=1e-5) = begin
    out = HCPlotPoint[]
    for p in pts
        any(_hp_dist(p, q) <= tol for q in out) || push!(out, p)
    end
    out
end

function _hd_same_poly(a, b; tol=1e-5)
    length(a) == length(b) || return false
    all(_hp_dist(a[i], b[i]) <= tol for i in eachindex(a)) && return true
    all(_hp_dist(a[i], b[end - i + 1]) <= tol for i in eachindex(a))
end

function _hd_event_sets(hc)
    g = hc.lc.genus
    lset = Set{Int}(abs(x) for i in 1:g for x in get(hc.lc.words, "l$i", Int[]))
    alphas = hc_curve_names(hc, :alpha)
    betas = hc_curve_names(hc, :beta)
    aset = Set{Int}(abs(x) for n in alphas for x in hc.lc.words[n])
    bset = Set{Int}(abs(x) for n in betas for x in hc.lc.words[n])
    labels = Dict{String,Set{Int}}()
    for n in alphas
        labels[n] = union(bset, lset)
    end
    for n in betas
        labels[n] = union(aset, lset)
    end
    (labels=labels, crossings=sort!(collect(intersect(aset, bset))), alphas=alphas, betas=betas, lset=lset)
end

function _hd_sites(p, contacts, audit, ordered, curve, label, crossings)
    if abs(label) in crossings
        hits = NamedTuple[]
        for (i, c) in enumerate(audit.crossings)
            get(audit.labels, i, 0) == abs(label) || continue
            for (ri, st) in ((c.run_a, c.station_a), (c.run_b, c.station_b))
                r = p.runs[ri]
                r.component == curve || continue
                pos = findfirst(x -> x.name == r.name, ordered)
                pos === nothing && continue
                push!(hits, (pos=pos, station=Float64(st), point=c.point))
            end
        end
        length(hits) == 1 || throw(ArgumentError("drawn crossing ($curve, $label) has $(length(hits)) stations"))
        return (arrive=hits[1], depart=hits[1])
    end
    arrive = nothing
    depart = nothing
    for (pos, r) in enumerate(ordered)
        if get(contacts, (r.name, :end), 0) == abs(label)
            arrive !== nothing && throw(ArgumentError("duplicate arrival port ($curve, $label)"))
            arrive = (pos=pos, station=Float64(length(r.points) - 1), point=last(r.points))
        end
        if get(contacts, (r.name, :start), 0) == abs(label)
            depart !== nothing && throw(ArgumentError("duplicate departure port ($curve, $label)"))
            depart = (pos=pos, station=0.0, point=first(r.points))
        end
    end
    (arrive === nothing || depart === nothing) &&
        throw(ArgumentError("drawn port ($curve, $label) has no contact on the placed runs"))
    (arrive=arrive, depart=depart)
end

function _hd_ordinary_edges(hc, result)
    info = _hd_event_sets(hc)
    p = result.placement
    edges = Tuple{Tuple{String,Int},Vector{Vector{HCPlotPoint}},String,Symbol}[]
    for name in vcat(info.alphas, info.betas)
        role = hc.roles[name]
        word = hc.lc.words[name]
        events = Int[]
        for (i, x) in enumerate(word)
            abs(x) in info.labels[name] && push!(events, i)
        end
        if isempty(events)
            closed = any(r.component == name && r.closed for r in p.runs)
            closed && throw(ArgumentError("closed curve $name has no drawn event; the specification gives the edge and not its key"))
            throw(ArgumentError("curve $name has no drawn event"))
        end
        ordered = _hd_ordered_runs(p, result.contacts, name)
        sites = map(events) do i
            _hd_sites(p, result.contacts, result.audit, ordered, name, word[i], info.crossings)
        end
        for (t, i) in enumerate(events)
            nxt = sites[mod1(t + 1, length(sites))]
            pieces = _hd_span_pieces(ordered, sites[t].depart, nxt.arrive)
            push!(edges, ((name, word[i]), pieces, name, role))
        end
    end
    (edges=edges, info=info, faces=nothing)
end

function _hd_reference_edges(hc, result)
    regions = _hp_mark_regions(hc, result.placement, result.contacts; geometry=result.audit)
    info = _hd_event_sets(hc)
    face_of = Dict{Any,Vector{String}}()
    faces = NamedTuple[]
    for id in sort!(collect(keys(regions.faces)))
        arcs = regions.faces[id]
        boundary = [(edge=a.edge, orientation=Int(a.orientation)) for a in arcs]
        for a in arcs
            push!(get!(face_of, a.edge, String[]), id)
        end
        push!(faces, (id=id, polygon=_hp_mark_polygon(arcs), boundary=boundary))
    end
    edges = Tuple{Tuple{String,Int},Vector{Vector{HCPlotPoint}},String,Symbol,Tuple}[]
    for name in vcat(info.alphas, info.betas)
        role = hc.roles[name]
        for x in hc.lc.words[name]
            key = (name, x)
            parts = if haskey(regions.intervals, key)
                [regions.intervals[key]]
            elseif haskey(regions.cutcopies, key)
                raw = [copy.points for copy in regions.cutcopies[key]]
                kept = Vector{HCPlotPoint}[]
                for q in raw
                    any(_hd_same_poly(q, u) for u in kept) || push!(kept, q)
                end
                kept
            else
                throw(ArgumentError("layout has no interval for edge $key"))
            end
            push!(edges, (key, parts, name, role, Tuple(get(face_of, key, String[]))))
        end
    end
    (edges=edges, info=info, faces=faces, regions=regions, areas=_hd_areas(hc, faces))
end

# cf: lib/topology.jl/docs/specs/drawing-objects.md §3.3.1 — a region is a component of Σ ∖ (α ∪ β): the faces glued
# across the K curves without a role, taken from forget_reference (heegaard-koda-conversion.md §4.1).
function _hd_areas(hc, faces)
    polygon = Dict(f.id => f.polygon for f in faces)
    areas = NamedTuple[]
    for r in forget_reference(realize_hc(hc)).regions
        missing_faces = [f for f in r.face_roots if !haskey(polygon, f)]
        isempty(missing_faces) || throw(ArgumentError("region $(r.id) uses faces without a polygon: $(join(missing_faces, ", "))"))
        push!(areas, (id=r.id, faces=Tuple(sort(r.face_roots; by=f -> parse(Int, f[2:end]))), parts=[polygon[f] for f in r.face_roots]))
    end
    covered = sort!([f for a in areas for f in a.faces])
    covered == sort!(collect(keys(polygon))) || throw(ArgumentError("regions do not partition the faces"))
    areas
end

function _hd_circle(center, radius)
    n = max(16, ceil(Int, 2π * radius / 10))
    [(center[1] + radius * cos(2π * k / n), center[2] + radius * sin(2π * k / n)) for k in 0:n-1]
end

function _hd_curve_parts(p, name)
    plain = [r for r in p.runs if r.component == name && !startswith(r.name, "@K:")]
    guides = [r for r in p.runs if r.component == name && startswith(r.name, "@K:")]
    !isempty(plain) && return [r.points for r in sort(plain; by=_hd_run_index)]
    !isempty(guides) && return [r.points for r in sort(guides; by=_hd_run_index)]
    if !p.square && startswith(name, "l")
        index = name[2:end]
        holes = [h for h in p.holes if h.first == "top" * index]
        isempty(holes) && throw(ArgumentError("curve $name has no positive opening"))
        center, radius = holes[1].second
        return [_hd_circle(center, radius)]
    end
    p.square || throw(ArgumentError("curve $name has no placed run"))
    # In screen coordinates l⁺ is the lower edge, at y = H.
    H, W = p.height, p.width
    if startswith(name, "m")
        return [[(0.0, 0.0), (0.0, H)], [(W, 0.0), (W, H)]]
    elseif startswith(name, "l")
        return [[(0.0, H), (W, H)]]
    end
    throw(ArgumentError("curve $name has no placed run on the square chart"))
end

function _hd_opening_key(name)
    m = match(r"^(top|bot)(\d+)$", name)
    m === nothing && throw(ArgumentError("opening name $name is outside the chart convention"))
    ("l$(m.captures[2])", m.captures[1] == "top" ? (:+) : (:-))
end

function _hd_angles(result, label)
    angles = Float64[]
    audit = result.audit
    for (i, c) in enumerate(audit.crossings)
        get(audit.labels, i, 0) == abs(label) && push!(angles, c.angle)
    end
    if haskey(audit, :reference)
        for row in audit.reference.angles
            row.label == abs(label) && push!(angles, row.degrees)
        end
    end
    angles
end

function _hd_crossing_points(p, contacts, audit, label, curves)
    pts = HCPlotPoint[]
    for (i, c) in enumerate(audit.crossings)
        get(audit.labels, i, 0) == abs(label) && push!(pts, c.point)
    end
    for curve in curves
        for r in p.runs
            r.component == curve && !r.closed || continue
            for side in (:start, :end)
                get(contacts, (r.name, side), 0) == abs(label) || continue
                push!(pts, side == :start ? first(r.points) : last(r.points))
            end
        end
    end
    got = _hd_dedup(pts)
    isempty(got) && throw(ArgumentError("crossing $label has no placed point"))
    got
end

function _hd_sign_on(word, label)
    x = only(v for v in word if abs(v) == abs(label))
    x > 0 ? (:+) : (:-)
end

_hd_role_curves(hc) = vcat(hc_curve_names(hc, :alpha), hc_curve_names(hc, :beta))

function _hd_theme()
    NamedTuple[
        (kind=:curve, selector=(;role=:alpha), attrs=(;color="#ba433f", linewidth=2.5, visible=true)),
        (kind=:curve, selector=(;role=:beta), attrs=(;color="#255dba", linewidth=2.5, visible=true)),
        (kind=:opening, selector=(;upper_alpha=false), attrs=(;color="#000000", linewidth=2.5, visible=true)),
        (kind=:opening, selector=(;upper_alpha=true), attrs=(;linewidth=2.5, visible=true)),
        (kind=:frame, selector=(;chart=:square), attrs=(;visible=true, color="#777777", linewidth=1.2)),
        (kind=:point_mark, selector=NamedTuple(), attrs=(;label=true, labelsize=15)),
        (kind=:disk_mark, selector=NamedTuple(), attrs=(;label=true, labelsize=14)),
    ]
end

function _hd_finish(hc, result, calls, path, built)
    p = result.placement
    prims = DrawingPrimitive[]
    owners = Int[]
    elements = DrawingElement[]
    info = built.info
    for (key, parts, curve, role, faces...) in built.edges
        facekeys = isempty(faces) ? () : faces[1]
        _hd_add!(prims, owners, elements, :edge, key,
            _hd_in(;role, curve, faces=facekeys), _hd_drawn(), parts)
    end
    if built.faces !== nothing
        for face in built.faces
            length(face.polygon) >= 3 || throw(ArgumentError("face $(face.id) has no polygon"))
            _hd_add!(prims, owners, elements, :face, face.id,
                _hd_in(;boundary=face.boundary, faces=(face.id,)), _hd_drawn(), [face.polygon]; shape=:polygon, closed=true)
        end
        for area in built.areas
            _hd_add!(prims, owners, elements, :region, area.id,
                _hd_in(;faces=area.faces), _hd_drawn(), area.parts; shape=:polygon, closed=true)
        end
    end
    for name in _hd_role_curves(hc)
        role = hc.roles[name]
        parts = _hd_curve_parts(p, name)
        _hd_add!(prims, owners, elements, :curve, name, _hd_in(;role, curve=name), _hd_drawn(), parts; closed=any(r.component == name && r.closed for r in p.runs))
    end
    if p.square
        H, W = p.height, p.width
        specs = [(("l1", :+), [(0.0, H), (W, H)]), (("l1", :-), [(0.0, 0.0), (W, 0.0)])]
        for (key, pts) in specs
            upper = key[2] == :+ && get(hc.roles, key[1], :none) == :alpha
            _hd_add!(prims, owners, elements, :opening, key, _hd_in(;curve=key[1], upper_alpha=upper), _hd_drawn(), [pts])
        end
        if !any(startswith(r.name, "@K:m") for r in p.runs)
            sides = [("m1", 1, [(0.0, 0.0), (0.0, H)]), ("m1", 2, [(W, 0.0), (W, H)])]
            for (curve, index, pts) in sides
                _hd_add!(prims, owners, elements, :reference, (curve, index), _hd_in(;curve), _hd_drawn(), [pts])
            end
        end
        _hd_add!(prims, owners, elements, :frame, nothing, _hd_in(;chart=:square), _hd_drawn(),
            [[(0.0, 0.0), (W, 0.0), (W, H), (0.0, H)]]; closed=true)
    else
        for (name, (center, radius)) in p.holes
            key = _hd_opening_key(name)
            upper = key[2] == :+ && get(hc.roles, key[1], :none) == :alpha
            _hd_add!(prims, owners, elements, :opening, key, _hd_in(;curve=key[1], upper_alpha=upper), _hd_drawn(), [_hd_circle(center, radius)]; closed=true)
        end
        center, radius = _hp_frame(p)
        _hd_add!(prims, owners, elements, :frame, nothing, _hd_in(;chart=:disk), _hd_drawn(), [_hd_circle(center, radius)]; closed=true)
    end
    for r in p.runs
        startswith(r.name, "@K:") || continue
        parts = split(r.name, ':')
        key = (String(parts[2]), parse(Int, parts[3]))
        _hd_add!(prims, owners, elements, :reference, key, _hd_in(;curve=key[1]), _hd_drawn(), [r.points])
    end
    groups = Dict{Tuple{String,Int},Vector{HCPlotPoint}}()
    for r in p.runs
        get(hc.roles, r.component, :none) in (:alpha, :beta) || continue
        r.closed && continue
        for side in (:start, :end)
            haskey(result.contacts, (r.name, side)) || continue
            label = result.contacts[(r.name, side)]
            push!(get!(groups, (r.component, label), HCPlotPoint[]), side == :start ? first(r.points) : last(r.points))
        end
    end
    for key in sort!(collect(keys(groups)); by=string)
        pts = _hd_dedup(groups[key])
        role = hc.roles[key[1]]
        _hd_add!(prims, owners, elements, :port, key, _hd_in(;role, curve=key[1]), _hd_drawn(), [[q] for q in pts]; shape=:point)
    end
    pairs = _hp_bigon_pairs(p.runs, result.audit.crossings)
    bigons = Tuple{Int,Int}[]
    bigon_labels = Set{Int}()
    for (a, b) in pairs
        la = get(result.audit.labels, a, nothing)
        lb = get(result.audit.labels, b, nothing)
        (la === nothing || lb === nothing) && continue
        push!(bigons, la < lb ? (la, lb) : (lb, la))
        push!(bigon_labels, la, lb)
    end
    for label in info.crossings
        alpha = nothing
        beta = nothing
        for name in _hd_role_curves(hc)
            any(abs(x) == label for x in hc.lc.words[name]) || continue
            hc.roles[name] == :alpha ? (alpha = name) : (beta = name)
        end
        (alpha === nothing || beta === nothing) && throw(ArgumentError("crossing $label is not an alpha–beta pair"))
        pts = _hd_crossing_points(p, result.contacts, result.audit, label, (alpha, beta))
        angles = _hd_angles(result, label)
        isempty(angles) && throw(ArgumentError("crossing $label has no placed angle"))
        _hd_add!(prims, owners, elements, :crossing, label,
            _hd_in(;sign=_hd_sign_on(hc.lc.words[alpha], label), bigon=label in bigon_labels, alpha, beta),
            _hd_drawn(;angle=minimum(angles)), [[q] for q in pts]; shape=:point)
    end
    if haskey(result, :marks)
        for point in result.marks.points
            _hd_add!(prims, owners, elements, :point_mark, point.name, _hd_in(;curve=point.support),
                _hd_drawn(;labelpoint=point.point), [[point.point]]; shape=:point)
        end
        for disk in result.marks.disks
            parts = Any[q for q in disk.paths]
            shapes = fill(:polygon, length(disk.paths))
            append!(parts, disk.boundary)
            append!(shapes, fill(:polyline, length(disk.boundary)))
            el = length(elements) + 1
            ids = Int[]
            for (pts, shape) in zip(parts, shapes)
                length(pts) < 2 && continue
                id = length(prims) + 1
                push!(prims, DrawingPrimitive(id, shape, _hd_tup(pts), shape == :polygon))
                push!(owners, el)
                push!(ids, id)
            end
            isempty(ids) && throw(ArgumentError("disk mark $(disk.name) has no part"))
            push!(elements, DrawingElement(:disk_mark, disk.name, _hd_in(), _hd_drawn(;labelpoint=disk.labelpoint), Tuple(ids)))
        end
    end
    HeegaardDrawing(Tuple(prims), Tuple(elements), Tuple(owners), Tuple(sort!(unique(bigons))), path, calls, DrawingStyleState(_hd_theme(), Any[]))
end

"""
    drawing(hc::HeegaardCode; infinity_face=nothing) -> HeegaardDrawing

Return the semantic elements and geometric parts of a planar Heegaard drawing.

# Arguments
`hc` must satisfy the [drawing domain](../specs/heegaard-plot.md).
`infinity_face` is passed to chart selection. Validation, placement and
element-construction errors propagate, including the `ArgumentError` for the
[single-crossing fixture](../specs/heegaard-plot.md#single-crossing-drawing-failure).

# Returns
`HeegaardDrawing` containing `primitives`, semantic `elements`, their `owner`
indices, `bigons`, the `path` tag, a layout-call serial and mutable `style`.
Their meaning is owned by [drawing objects](../specs/drawing-objects.md).

# Guarantees
1. Returned coordinates and semantic parts refer to the input diagram's
   planar chart. [proved/agent: by construction]
2. The input is not modified. A layout-call counter advances and the result
   has its own initially empty style rules. [proved/agent: by construction]

# Not guaranteed
A successful `plot` does not ensure successful `drawing`; the linked bug
affects the single-crossing fixture. Angles and distances are chart geometry,
not topological invariants.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md §§3–6 — semantic drawing representation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-2], "m1"=>[1,3], "c1"=>[-4,-3], "c2"=>[4,2]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("c1"=>:alpha, "c2"=>:beta));

julia> try H.drawing(hc) catch e; occursin("open run c1#0 has a drawn edge against its direction", sprint(showerror,e)) end
true
```

Role: draw
"""
function drawing(hc::HeegaardCode; infinity_face=nothing)
    _hd_layout_calls[] += 1
    calls = _hd_layout_calls[]
    result = plot(hc; infinity_face)
    placing = haskey(result, :reference_alpha) || haskey(result, :mark_contacts)
    built = placing ? _hd_reference_edges(hc, result) : _hd_ordinary_edges(hc, result)
    _hd_finish(hc, result, calls, placing ? :reference : :ordinary, built)
end

function _hd_match(el, sel)
    for (k, v) in pairs(sel)
        if k == :name
            keys = v isa AbstractVector ? v : (v,)
            el.key in keys || return false
        elseif k == :face
            wanted = v isa AbstractVector ? Tuple(v) : (v,)
            any(f -> f in el.input.faces, wanted) || return false
        elseif k in _HD_SELS
            haskey(el.input, k) || return false
            getproperty(el.input, k) == v || return false
        else
            throw(ArgumentError("unknown drawing selector $k"))
        end
    end
    true
end

function _hd_base(el)
    kind = el.kind
    z = kind == :frame ? 0 : kind == :face || kind == :region || kind == :disk_mark ? 1 : kind == :reference ? 2 :
        kind == :curve ? 3 : kind == :edge ? 4 : kind == :opening ? 5 : kind == :port ? 6 :
        kind == :crossing ? 7 : 8
    attrs = Dict{Symbol,Any}(:visible => true, :color => "#555555", :linewidth => 2.5, :linestyle => nothing,
        :fill => ("#000000", 0.0), :marker => :circle, :markersize => 3.5, :label => false, :labelsize => 12, :zorder => z)
    if kind == :crossing
        attrs[:visible] = true
        attrs[:color] = (:black, 0.0)
        attrs[:markersize] = 8
    elseif kind == :port
        attrs[:visible] = false
        attrs[:markersize] = 2.8
        attrs[:color] = "#255dba"
    elseif kind == :reference
        attrs[:visible] = false
        attrs[:linewidth] = 1.5
    elseif kind == :frame
        attrs[:visible] = false
        attrs[:color] = "#777777"
        attrs[:linewidth] = 1.2
    elseif kind == :face || kind == :region
        attrs[:linewidth] = 0
    elseif kind == :point_mark
        attrs[:color] = "#252525"
        attrs[:markersize] = 3.5
    elseif kind == :disk_mark
        attrs[:fill] = "#eee2b8"
        attrs[:color] = "#9a813c"
        attrs[:linewidth] = 1.7
    elseif kind == :opening
        attrs[:color] = "#000000"
    end
    attrs
end

function _hd_find(d, kind, key)
    for el in d.elements
        el.kind == kind && el.key == key && return el
    end
    throw(ArgumentError("no $kind element $key"))
end

"""
    resolved_style(d::HeegaardDrawing, el::DrawingElement) -> Dict{Symbol,Any}

Resolve the current visual attributes of one semantic drawing element.

# Arguments
`el` must carry the fields required by its kind and any selectors in `d.style`. Parent-curve references of edges and openings must exist; missing parents raise `ArgumentError`. Malformed fields or style data can raise access errors. Attribute semantics follow [drawing styles](../specs/drawing-objects.md).

# Returns
A new dictionary of effective attributes, including visibility, color, line and marker properties, labels and drawing order.

# Guarantees
1. Matching later style rules override earlier ones; inherited curve attributes apply to child edges. The result describes appearance, without changing the mathematical input or placement. [proved/agent: by construction]
2. The drawing and element are not modified. [proved/agent: by construction]

# Not guaranteed
This does not validate backend-specific attribute values or redraw a figure. Nested mutable attribute values need not be copied.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md §5 — style precedence and semantic selectors.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> e = H.DrawingElement(:curve, "c", (;role=:alpha), (;), ()); d = H.HeegaardDrawing((), (e,), (), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> H.style!(d, :curve; name="c", color="#ff00aa"); H.resolved_style(d,e)[:color]
"#ff00aa"
```

Role: access
"""
function resolved_style(d::HeegaardDrawing, el::DrawingElement)
    attrs = _hd_base(el)
    if el.kind == :edge && el.input.curve !== nothing
        parent = _hd_find(d, :curve, el.input.curve)
        pa = resolved_style(d, parent)
        for k in (:visible, :color, :linewidth, :linestyle)
            attrs[k] = pa[k]
        end
    end
    specified = Set{Symbol}()
    for rule in d.style.rules
        rule.kind == el.kind || continue
        _hd_match(el, rule.selector) || continue
        for (k, v) in pairs(rule.attrs)
            attrs[k] = v
            push!(specified, k)
        end
    end
    if el.kind == :opening && el.input.upper_alpha && !(:color in specified)
        parent = _hd_find(d, :curve, el.input.curve)
        pa = resolved_style(d, parent)
        pa[:visible] && (attrs[:color] = pa[:color])
    end
    attrs
end

"""
    style!(d::HeegaardDrawing, kind::Symbol; kwargs...) -> HeegaardDrawing

Append a visual rule to a semantic Heegaard drawing.

# Arguments
Keyword names must be attributes or selectors from [drawing styles](../specs/drawing-objects.md#5-見た目の規則); unknown names raise `ArgumentError`. Attribute values and the kind are not fully validated here. Existing style listeners may throw, in which case the appended rule remains stored.

# Returns
The same drawing object `d`.

# Guarantees
1. The new rule participates after earlier rules in `resolved_style`; already attached backend listeners are notified. [proved/agent: by construction]
2. The style rule list is mutated, while primitive geometry and mathematical input data are not replaced or recomputed. [proved/agent: by construction]

# Not guaranteed
No matching element is required, and backend compatibility of arbitrary values is not guaranteed. A listener failure does not roll back the rule.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md §5 — semantic style rules.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> e = H.DrawingElement(:curve, "c", (;role=:alpha), (;), ()); d = H.HeegaardDrawing((), (e,), (), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> H.style!(d,:curve; name="c",color="#ff00aa") === d
true
```

Role: other
"""
function style!(d::HeegaardDrawing, kind::Symbol; kwargs...)
    sel = Pair{Symbol,Any}[]
    attrs = Pair{Symbol,Any}[]
    for (k, v) in kwargs
        if k in _HD_ATTRS
            push!(attrs, k => v)
        elseif k in _HD_SELS
            push!(sel, k => v)
        else
            throw(ArgumentError("unknown drawing attribute or selector $k"))
        end
    end
    push!(d.style.rules, (kind=kind, selector=NamedTuple(sel), attrs=NamedTuple(attrs)))
    for listener in d.style.listeners
        listener()
    end
    d
end

"""
    curve_incidents(d::HeegaardDrawing, name::AbstractString) -> NamedTuple

Return the edge and crossing keys incident to a named curve in a Heegaard drawing.

# Arguments
`d` is normally produced by `drawing`; `name` is a stored curve name.
No existence check is made: an absent name returns empty lists. The caller
must ensure that manually constructed edge and crossing elements have the
input attributes required by [drawing relations](../specs/drawing-objects.md#35-関係);
missing attributes can raise property-access errors.

# Returns
A named tuple with `edges` and `crossings`, each a fresh vector of element
keys in their stored order.

# Guarantees
1. The lists contain exactly the stored edge keys assigned to `name` and
   crossing keys whose alpha or beta curve is `name`.
   [proved/agent: by construction]
2. For a drawing of an HC, these are incidences in its displayed curve
   arrangement. The drawing is not modified. [proved/agent: by construction]

# Not guaranteed
This performs no geometric intersection test. Empty output need not mean an
embedded curve is disjoint; the name may be absent from the stored elements.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.DrawingElement(:crossing, 4, (alpha="c1", beta="c2"), (;), ());

julia> d = H.HeegaardDrawing((), (c,), (), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> Tuple(H.curve_incidents(d, "c1").crossings)
(4,)
```

See also [`crossing_curves`](@ref), [`drawing`](@ref).

Role: access
"""
function curve_incidents(d::HeegaardDrawing, name::AbstractString)
    edges = [el.key for el in d.elements if el.kind == :edge && el.input.curve == name]
    crossings = [el.key for el in d.elements if el.kind == :crossing && (el.input.alpha == name || el.input.beta == name)]
    (;edges, crossings)
end

"""
    face_boundary(d::HeegaardDrawing, name::AbstractString) -> Tuple

Return the stored oriented edge boundary of a named face in a Heegaard drawing.

# Arguments
`d` normally comes from `drawing`, and `name` is a face key. A missing face
throws `ArgumentError` naming it. A manually constructed face must carry
the tuple-valued boundary attribute of [drawing relations](../specs/drawing-objects.md#35-関係);
this is not validated, and missing attributes cause property-access errors.

# Returns
The face's stored boundary tuple, with edge keys and orientation values.

# Guarantees
1. For a drawing that records faces, the result is the indicated face's
   stored boundary in the input curve arrangement.
   [proved/agent: by construction]
2. No data are modified or recomputed. [proved/agent: by construction]

# Not guaranteed
Not every drawing path records faces. This is not a boundary of a whole
complement region after forgetting K, nor a fresh geometric validation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> f = H.DrawingElement(:face, "F1", (boundary=((edge=("a",1), orientation=1),),), (;), ());

julia> d = H.HeegaardDrawing((), (f,), (), (), :reference, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> !isempty(H.face_boundary(d, "F1"))
true
```

See also [`drawing`](@ref), [`curve_incidents`](@ref).

Role: access
"""
function face_boundary(d::HeegaardDrawing, name::AbstractString)
    el = _hd_find(d, :face, String(name))
    el.input.boundary
end

"""
    crossing_curves(d::HeegaardDrawing, label::Integer) -> NamedTuple

Return the alpha and beta curve names recorded at a crossing of a Heegaard drawing.

# Arguments
`d` is normally produced by `drawing`; `label` is an intersection key
convertible to `Int`. Conversion errors propagate. A missing crossing throws
`ArgumentError` identifying the absent element. For manually constructed
drawings, the caller must provide alpha and beta input attributes on that
element; otherwise property-access errors propagate.

# Returns
A named tuple `(alpha=..., beta=...)` with the selected element's curve names,
in the sense of [drawing relations](../specs/drawing-objects.md#35-関係).

# Guarantees
1. The names equal the selected crossing's stored alpha and beta attributes;
   for output of `drawing`, they identify the incident input curves.
   [proved/agent: by construction]
2. The drawing and its style are not modified. [proved/agent: by construction]

# Not guaranteed
The accessor does not compute a crossing sign, angle or geometric validation.
A manually assembled record need not describe an actual crossing.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.DrawingElement(:crossing, 4, (alpha="c1", beta="c2"), (;), ());

julia> d = H.HeegaardDrawing((), (c,), (), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> H.crossing_curves(d, 4)
(alpha = "c1", beta = "c2")
```

See also [`curve_incidents`](@ref), [`drawing`](@ref).

Role: access
"""
function crossing_curves(d::HeegaardDrawing, label::Integer)
    el = _hd_find(d, :crossing, Int(label))
    (alpha=el.input.alpha, beta=el.input.beta)
end

"""
    bigon_corners(d::HeegaardDrawing) -> Tuple{Vararg{Tuple{Int,Int}}}

Return the recorded crossing-label pairs at the corners of a Heegaard drawing's bigons.

# Arguments
`d` is normally produced by `drawing`. No consistency check is performed on
a manually constructed drawing; its interpretation is the caller's responsibility.

# Returns
The tuple stored in `d.bigons`, interpreted as the
[bigon relation](../specs/drawing-objects.md#35-関係) of the drawing.

# Guarantees
1. The result equals the stored bigon pairs and does not change the drawing,
   its style or its geometry. [proved/agent: by construction]
2. For output of `drawing`, these are the recorded bigon corners of that
   displayed curve arrangement, rather than a manifold invariant.
   [proved/agent: by construction]

# Not guaranteed
This accessor performs no geometric search or validation. The pairs alone
do not certify applicability of a bigon-removal move.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.HeegaardDrawing((), (), (), ((4,5), (5,6)), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> length(H.bigon_corners(d))
2
```

See also [`drawing`](@ref), [`crossing_curves`](@ref).

Role: access
"""
bigon_corners(d::HeegaardDrawing) = d.bigons

"""
    element_at(d::HeegaardDrawing, primitive_id::Integer, point_index::Integer) -> DrawingElement

Return the semantic drawing element owning a selected primitive point.

# Arguments
`primitive_id` and `point_index` are 1-based positions in `d.primitives` and
that primitive's points. Each range is checked; `ArgumentError` identifies
an unknown primitive or out-of-range point. The caller must supply consistent
owner references in a manually constructed drawing; malformed references can
raise `BoundsError`.

# Returns
The `DrawingElement` selected by the primitive's owner entry, as specified
by [drawing reverse lookup](../specs/drawing-objects.md#36-逆引き).

# Guarantees
1. All valid point indices of one primitive return its stored owner element.
   [proved/agent: by construction]
2. The drawing, style and coordinates are not modified.
   [proved/agent: by construction]

# Not guaranteed
This does not perform nearest-point selection or validate a geometric hit.
Positions are array indices, not an arbitrary primitive record's `id` field.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> p = H.DrawingPrimitive(1, :point, ((0.0,0.0),), false); c = H.DrawingElement(:crossing, 4, (alpha="c1", beta="c2"), (;), (1,));

julia> d = H.HeegaardDrawing((p,), (c,), (1,), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[]));

julia> H.element_at(d, 1, 1).kind
:crossing
```

See also [`drawing_label`](@ref), [`drawing_inspector_label`](@ref).

Role: access
"""
function element_at(d::HeegaardDrawing, primitive_id::Integer, point_index::Integer)
    1 <= primitive_id <= length(d.primitives) || throw(ArgumentError("unknown primitive $primitive_id"))
    prim = d.primitives[primitive_id]
    1 <= point_index <= length(prim.points) || throw(ArgumentError("point index $point_index is outside primitive $primitive_id"))
    d.elements[d.owner[primitive_id]]
end

function _hd_curve_text(name)
    m = match(r"^([A-Za-z]+)(\d*)$", name)
    m === nothing && return name
    head = m.captures[1] == "alpha" ? "α" : m.captures[1] == "beta" ? "β" : m.captures[1]
    head * m.captures[2]
end

"""
    drawing_label(el::DrawingElement, spec) -> Union{Nothing,String}

Return the display label requested for a semantic Heegaard drawing element.

# Arguments
`spec=false` suppresses the label, `true` selects its key-based label,
a `Function` supplies `spec(el)`, and any other value is converted with
`string`. Callback and string-conversion errors propagate. No drawing-validity
check is performed. A callback is responsible for its own side effects.

# Returns
`nothing` for suppression or a `String` for display; the true case uses
the crossing-key label for crossings and the stringified key otherwise.

# Guarantees
1. The result follows the supplied label selector for the
   [drawing element](../specs/drawing-objects.md#3-意味の要素).
   [proved/agent: by construction]
2. The operation itself does not change `el`; any mutation by a supplied
   callback or custom string conversion is outside this guarantee.
   [proved/agent: by construction]

# Not guaranteed
Labels need not be unique or geometrically meaningful. This supplies no
layout, escaping for a target renderer or mathematical identification.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.DrawingElement(:crossing, 4, (alpha="c1", beta="c2"), (;), ());

julia> H.drawing_label(c, true)
"V4"
```

See also [`element_at`](@ref), [`drawing_inspector_label`](@ref).

Role: draw
"""
function drawing_label(el::DrawingElement, spec)
    spec === false && return nothing
    spec === true && return el.kind == :crossing ? "V$(el.key)" : string(el.key)
    spec isa Function && return string(spec(el))
    string(spec)
end

"""
    drawing_inspector_label(d::HeegaardDrawing, el::DrawingElement) -> String

Return inspector text for one semantic element of a Heegaard drawing.

# Arguments
`d` is a drawing and `el` is an element in the [drawing representation](../specs/drawing-objects.md).
Membership is not checked. Crossing elements must supply `input.sign`,
`input.alpha`, `input.beta` and a finite `drawn.angle` convertible to a rounded
`Int`; missing fields and conversion errors propagate.

# Returns
A string naming the element; crossing text also carries its curve names,
sign and rounded chart angle in degrees.

# Guarantees
1. Noncrossing elements are labeled by their kind and key; crossing labels
   refer to the supplied element's chart data. [proved/agent: by construction]
2. Neither argument is modified. [proved/agent: by construction]

# Not guaranteed
The label does not validate membership, crossing data or a topological
intersection number. It is presentation text, not a stable serialization.

# Sources
cf: lib/topology.jl/docs/specs/drawing-objects.md — semantic element labels and chart attributes.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.HeegaardDrawing((), (), (), (), :ordinary, 0, H.DrawingStyleState(NamedTuple[], Any[])); el = H.DrawingElement(:edge, 7, (;), (;), ());

julia> H.drawing_inspector_label(d, el)
"edge 7"
```

Role: draw
"""
function drawing_inspector_label(d::HeegaardDrawing, el::DrawingElement)
    if el.kind == :crossing
        sig = el.input.sign == :+ ? "+" : "-"
        ang = round(Int, el.drawn.angle)
        return string("V", el.key, "、", _hd_curve_text(el.input.beta), " × ", _hd_curve_text(el.input.alpha), "、符号 ", sig, "、交差角 ", ang, "°")
    end
    string(el.kind, " ", el.key)
end

export drawing, style!, drawingplot, drawingplot!, element_at, resolved_style,
    drawing_inspector_label, drawing_label, curve_incidents, face_boundary, crossing_curves, bigon_corners,
    HeegaardDrawing, DrawingElement, DrawingPrimitive
