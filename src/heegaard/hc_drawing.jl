# cf: lib/topology.jl/docs/specs/drawing-objects.md §3–§6, §8 item 1
# Coordinates are copied from plot(hc). Ordinary edges run between drawn events
# (alpha–beta crossings and ports on l_i). Reference and marked edges are the
# LC intervals from _hp_mark_regions, and only that path has faces.

const _hd_layout_calls = Ref(0)

struct DrawingPrimitive
    id::Int
    shape::Symbol
    points::Tuple{Vararg{NTuple{2,Float64}}}
    closed::Bool
end

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

function drawingplot end
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

"""Read-only drawing of a resolved Heegaard diagram. Uses the same placement as `plot`."""
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

"""Resolved attributes of one element. Later `style!` rules win."""
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

"""Append one style rule. Does not place or smooth the diagram again."""
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

function curve_incidents(d::HeegaardDrawing, name::AbstractString)
    edges = [el.key for el in d.elements if el.kind == :edge && el.input.curve == name]
    crossings = [el.key for el in d.elements if el.kind == :crossing && (el.input.alpha == name || el.input.beta == name)]
    (;edges, crossings)
end

function face_boundary(d::HeegaardDrawing, name::AbstractString)
    el = _hd_find(d, :face, String(name))
    el.input.boundary
end

function crossing_curves(d::HeegaardDrawing, label::Integer)
    el = _hd_find(d, :crossing, Int(label))
    (alpha=el.input.alpha, beta=el.input.beta)
end

bigon_corners(d::HeegaardDrawing) = d.bigons

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

function drawing_label(el::DrawingElement, spec)
    spec === false && return nothing
    spec === true && return el.kind == :crossing ? "V$(el.key)" : string(el.key)
    spec isa Function && return string(spec(el))
    string(spec)
end

"""Text for Makie's `inspector_label`. Crossings read 「V17、β2 × α1、符号 +、交差角 51°」."""
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
