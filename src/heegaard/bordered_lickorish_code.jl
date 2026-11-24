# Bordered Lickorish code, docs/specs/bordered-lickorish-code.md §2–§6 and §7.7.
# Face orbits are ρ ∘ τ, the same edge convention as lickorish-code.md §3.1 (D5).

struct BLCLabel
    sign::Int
    n::Int
end

struct BLCOrigin
    parent::String
    component::String
    start::BLCLabel
    cut::Union{Nothing, BLCLabel}
end

struct BLCBoundary
    name::String
    base::String
    points::Vector{String}
end

struct BLCArc
    kind::Symbol
    name::String
    start::String
    stop::String
    word::Vector{BLCLabel}
    origin::Union{Nothing, BLCOrigin}
end

struct BLCLoop
    kind::Symbol
    name::String
    word::Vector{BLCLabel}
    origin::Union{Nothing, BLCOrigin}
end

struct BLCPlace
    components::Vector{String}
    host::String
    outside::String
end

struct BorderedLickorishCode
    name::String
    genus::Int
    nboundaries::Int
    boundaries::Vector{BLCBoundary}
    k_arcs::Vector{BLCArc}
    k_loops::Vector{BLCLoop}
    q_arcs::Vector{BLCArc}
    q_loops::Vector{BLCLoop}
    arcs::Vector{BLCArc}
    loops::Vector{BLCLoop}
    places::Vector{BLCPlace}
end

Base.:(==)(a::BLCLabel, b::BLCLabel) = a.sign == b.sign && a.n == b.n
Base.hash(a::BLCLabel, h::UInt) = hash((a.sign, a.n), h)

Base.:(==)(a::BLCOrigin, b::BLCOrigin) =
    a.parent == b.parent && a.component == b.component && a.start == b.start && a.cut == b.cut
Base.hash(a::BLCOrigin, h::UInt) = hash((a.parent, a.component, a.start, a.cut), h)

Base.:(==)(a::BLCBoundary, b::BLCBoundary) =
    a.name == b.name && a.base == b.base && a.points == b.points
Base.hash(a::BLCBoundary, h::UInt) = hash((a.name, a.base, a.points), h)

Base.:(==)(a::BLCArc, b::BLCArc) =
    a.kind == b.kind && a.name == b.name && a.start == b.start && a.stop == b.stop &&
    a.word == b.word && a.origin == b.origin
Base.hash(a::BLCArc, h::UInt) = hash((a.kind, a.name, a.start, a.stop, a.word, a.origin), h)

Base.:(==)(a::BLCLoop, b::BLCLoop) =
    a.kind == b.kind && a.name == b.name && a.word == b.word && a.origin == b.origin
Base.hash(a::BLCLoop, h::UInt) = hash((a.kind, a.name, a.word, a.origin), h)

Base.:(==)(a::BLCPlace, b::BLCPlace) =
    a.components == b.components && a.host == b.host && a.outside == b.outside
Base.hash(a::BLCPlace, h::UInt) = hash((a.components, a.host, a.outside), h)

Base.:(==)(a::BorderedLickorishCode, b::BorderedLickorishCode) =
    a.name == b.name && a.genus == b.genus && a.nboundaries == b.nboundaries &&
    a.boundaries == b.boundaries && a.k_arcs == b.k_arcs && a.k_loops == b.k_loops &&
    a.q_arcs == b.q_arcs && a.q_loops == b.q_loops && a.arcs == b.arcs &&
    a.loops == b.loops && a.places == b.places
Base.hash(a::BorderedLickorishCode, h::UInt) =
    hash((a.name, a.genus, a.nboundaries, a.boundaries, a.k_arcs, a.k_loops,
          a.q_arcs, a.q_loops, a.arcs, a.loops, a.places), h)

struct BorderedLCError <: Exception
    section::String
    reason::Symbol
    detail::String
end

Base.showerror(io::IO, e::BorderedLCError) =
    print(io, "BorderedLCError(", e.section, ", ", e.reason, "): ", e.detail)

struct BorderedLCDiagnosis
    V::Int
    E::Int
    F::Int
    boundary_cycles::Int
    chi::Int
    connected::Bool
    boundary_ports::Vector{Vector{String}}
end

struct BorderedCut
    piece::BorderedLickorishCode
    cutter::String
    genus::Int
    places::Dict{String, String}
    disks::Dict{String, Vector{String}}
    points::Dict{String, String}
    weights::Dict{String, Int}
end

struct _BLCDart
    kind::Int
    name::String
    label::Int
    sign::Int
    port::String
    dir::Int
end

struct _BLCCurve
    name::String
    ref::Bool
    isloop::Bool
    startp::String
    endp::String
    word::Vector{BLCLabel}
end

function _blc_fail(section, reason, detail)
    throw(BorderedLCError(section, reason, detail))
end

_blc_fmt(lab::BLCLabel) = (lab.sign < 0 ? "-" : "+") * string(lab.n)

function _blc_parse_label(tok::AbstractString)
    m = match(r"^([+-]?)(\d+)$", tok)
    m === nothing && _blc_fail("3", :word, "bad label $tok")
    m.captures[2] == "0" && _blc_fail("3", :word, "label 0")
    BLCLabel(m.captures[1] == "-" ? -1 : 1, parse(Int, m.captures[2]))
end

function _blc_parse_word(raw::AbstractString)
    s = strip(raw)
    isempty(s) && return BLCLabel[]
    [_blc_parse_label(tok) for tok in split(s)]
end

function _blc_parse_origin(text::AbstractString, component::AbstractString)
    m = match(r"^([^.]+)\.([^@]+)@([+-]?\d+)(?::([+-]?\d+))?$", text)
    m === nothing && _blc_fail("3", :origin, "bad origin $text")
    m.captures[2] == component || _blc_fail("3", :origin, "origin component $(m.captures[2])")
    start = _blc_parse_label(m.captures[3])
    cut = m.captures[4] === nothing ? nothing : _blc_parse_label(m.captures[4])
    BLCOrigin(m.captures[1], component, start, cut)
end

function _blc_ints(word::Vector{BLCLabel})
    [lab.sign * lab.n for lab in word]
end

function _blc_from_ints(word::Vector{Int})
    out = BLCLabel[]
    for x in word
        x == 0 && _blc_fail("6.1", :word, "label 0")
        push!(out, BLCLabel(x < 0 ? -1 : 1, abs(x)))
    end
    out
end

function _blc_curves(piece::BorderedLickorishCode)
    out = _BLCCurve[]
    for a in piece.k_arcs
        push!(out, _BLCCurve(a.name, true, false, a.start, a.stop, a.word))
    end
    for a in piece.k_loops
        push!(out, _BLCCurve(a.name, true, true, "", "", a.word))
    end
    for a in piece.q_arcs
        push!(out, _BLCCurve(a.name, true, false, a.start, a.stop, a.word))
    end
    for a in piece.q_loops
        push!(out, _BLCCurve(a.name, true, true, "", "", a.word))
    end
    for a in piece.arcs
        push!(out, _BLCCurve(a.name, false, false, a.start, a.stop, a.word))
    end
    for a in piece.loops
        push!(out, _BLCCurve(a.name, false, true, "", "", a.word))
    end
    out
end

function _blc_names(piece::BorderedLickorishCode)
    names = String[]
    for c in _blc_curves(piece)
        push!(names, c.name)
    end
    names
end

function _blc_is_path(input::AbstractString)
    (isempty(input) || '\n' in input || '\r' in input) && return false
    try
        return isfile(input)
    catch
        return false
    end
end

function read_bordered_lc(input::AbstractString)
    text = _blc_is_path(input) ? read(input, String) : String(input)
    _blc_parse(text)
end

function _blc_parse(text::AbstractString)
    name = ""
    genus = -1
    nboundaries = -1
    boundaries = BLCBoundary[]
    k_arcs = BLCArc[]
    k_loops = BLCLoop[]
    q_arcs = BLCArc[]
    q_loops = BLCLoop[]
    arcs = BLCArc[]
    loops = BLCLoop[]
    places = BLCPlace[]
    seen_arc = Set{String}()
    seen_loop = Set{String}()
    rank = 0
    started = false
    for raw in split(text, '\n')
        line = strip(replace(raw, '\r' => ""))
        line = strip(split(line, '#'; limit=2)[1])
        isempty(line) && continue
        started = true
        if (m = match(r"^piece\s+(\S+)\s*:\s*genus=(-?\d+)\s+boundaries=(-?\d+)\s*$", line)) !== nothing
            rank > 0 && _blc_fail("3", :order, "piece is not first")
            rank = 1
            name = String(m.captures[1])
            genus = parse(Int, m.captures[2])
            nboundaries = parse(Int, m.captures[3])
            continue
        end
        if match(r"^marking\s*:\s*reference\s*$", line) !== nothing
            rank <= 1 || _blc_fail("3", :order, "marking")
            rank = 2
            continue
        end
        if (m = match(r"^boundary\s+(\S+)\s+base=(\S+)\s*:\s*(.*)$", line)) !== nothing
            rank <= 3 || _blc_fail("3", :order, "boundary")
            rank = 3
            pts = String.(split(strip(m.captures[3])))
            isempty(pts) && _blc_fail("3", :boundary, "empty boundary")
            String(m.captures[2]) in pts || _blc_fail("3", :boundary, "base not on $(m.captures[1])")
            push!(boundaries, BLCBoundary(String(m.captures[1]), String(m.captures[2]), pts))
            continue
        end
        if (m = match(r"^place\s+\{([^}]*)\}\s*:\s*(\S+)\s*;\s*outside=(\S+)\s*$", line)) !== nothing
            rank <= 10 || _blc_fail("3", :order, "place")
            rank = 10
            comps = String.(split(strip(m.captures[1])))
            isempty(comps) && _blc_fail("3", :place, "empty place")
            push!(places, BLCPlace(comps, String(m.captures[2]), String(m.captures[3])))
            continue
        end
        m = match(r"^(?:(K|Q)\s+)?(arc|loop)\s+(\S+)(?:\s+(\S+)\s+->\s+(\S+))?(?:\s+origin=(\S+))?\s*:\s*(.*)$", line)
        m === nothing && _blc_fail("3", :line, line)
        prefix = m.captures[1]
        shape = m.captures[2]
        cname = String(m.captures[3])
        kind = prefix === nothing ? :input : (prefix == "K" ? :K : :Q)
        origin = m.captures[6] === nothing ? nothing : _blc_parse_origin(m.captures[6], cname)
        if shape == "arc"
            (m.captures[4] === nothing || m.captures[5] === nothing) &&
                _blc_fail("3", :arc, cname)
            body = strip(m.captures[7])
            (startswith(body, "[") && endswith(body, "]")) || _blc_fail("3", :word, cname)
            word = _blc_parse_word(body[2:end-1])
            cname in seen_loop && _blc_fail("3", :duplicate, cname)
            push!(seen_arc, cname)
            arc = BLCArc(kind, cname, String(m.captures[4]), String(m.captures[5]), word, origin)
            cat = kind == :K ? 4 : kind == :Q ? 6 : 8
            rank <= cat || _blc_fail("3", :order, cname)
            rank = cat
            kind == :K ? push!(k_arcs, arc) : kind == :Q ? push!(q_arcs, arc) : push!(arcs, arc)
        else
            m.captures[4] === nothing || _blc_fail("3", :loop, cname)
            body = strip(m.captures[7])
            (startswith(body, "(") && endswith(body, ")")) || _blc_fail("3", :word, cname)
            word = _blc_parse_word(body[2:end-1])
            (cname in seen_loop || cname in seen_arc) && _blc_fail("3", :duplicate, cname)
            push!(seen_loop, cname)
            loop = BLCLoop(kind, cname, word, origin)
            cat = kind == :K ? 5 : kind == :Q ? 7 : 9
            rank <= cat || _blc_fail("3", :order, cname)
            rank = cat
            kind == :K ? push!(k_loops, loop) : kind == :Q ? push!(q_loops, loop) : push!(loops, loop)
        end
    end
    started || _blc_fail("3", :empty, "no piece")
    isempty(name) && _blc_fail("3", :piece, "missing piece")
    genus < 0 && _blc_fail("3", :genus, string(genus))
    nboundaries < 1 && _blc_fail("3", :boundaries, string(nboundaries))
    length(boundaries) == nboundaries ||
        _blc_fail("3", :boundaries, "declared $nboundaries, wrote $(length(boundaries))")
    BorderedLickorishCode(name, genus, nboundaries, boundaries, k_arcs, k_loops,
        q_arcs, q_loops, arcs, loops, places)
end

function _blc_word_text(word::Vector{BLCLabel}, bracket::Tuple{Char, Char})
    string(bracket[1], join((_blc_fmt(lab) for lab in word), " "), bracket[2])
end

function _blc_origin_text(origin::BLCOrigin)
    text = " origin=$(origin.parent).$(origin.component)@$(_blc_fmt(origin.start))"
    origin.cut === nothing ? text : text * ":" * _blc_fmt(origin.cut)
end

function bordered_lc_text(piece::BorderedLickorishCode)
    lines = String[
        "piece $(piece.name) : genus=$(piece.genus) boundaries=$(piece.nboundaries)",
        "marking : reference",
    ]
    for b in piece.boundaries
        push!(lines, "boundary $(b.name) base=$(b.base) : $(join(b.points, " "))")
    end
    function arc_line(prefix, arc)
        head = isempty(prefix) ? "arc" : prefix * " arc"
        origin = arc.origin === nothing ? "" : _blc_origin_text(arc.origin)
        "$head $(arc.name) $(arc.start) -> $(arc.stop)$origin : $(_blc_word_text(arc.word, ('[', ']')))"
    end
    function loop_line(prefix, loop)
        head = isempty(prefix) ? "loop" : prefix * " loop"
        origin = loop.origin === nothing ? "" : _blc_origin_text(loop.origin)
        "$head $(loop.name)$origin : $(_blc_word_text(loop.word, ('(', ')')))"
    end
    for a in piece.k_arcs
        push!(lines, arc_line("K", a))
    end
    for a in piece.k_loops
        push!(lines, loop_line("K", a))
    end
    for a in piece.q_arcs
        push!(lines, arc_line("Q", a))
    end
    for a in piece.q_loops
        push!(lines, loop_line("Q", a))
    end
    for a in piece.arcs
        push!(lines, arc_line("", a))
    end
    for a in piece.loops
        push!(lines, loop_line("", a))
    end
    for place in piece.places
        push!(lines, "place {$(join(place.components, " "))} : $(place.host) ; outside=$(place.outside)")
    end
    join(lines, "\n") * "\n"
end

function _blc_find(parent, x)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    x
end

function _blc_cellulation(curves::Vector{_BLCCurve}, boundaries::Vector{BLCBoundary})
    own = Dict{Tuple{Int, Int}, String}()
    for c in curves
        for lab in c.word
            key = (lab.n, lab.sign)
            haskey(own, key) && _blc_fail("4.1", :label, "repeated $( _blc_fmt(lab) )")
            own[key] = c.name
        end
    end
    labels = sort!(unique(lab.n for c in curves for lab in c.word))
    for n in labels
        haskey(own, (n, 1)) && haskey(own, (n, -1)) ||
            _blc_fail("4.1", :label, "label $n needs both signs")
    end
    vindex = Dict{Tuple, Int}()
    function vid(key)
        get!(vindex, key) do
            length(vindex) + 1
        end
    end
    for b in boundaries
        for p in b.points
            vid((:b, b.name, p))
        end
    end
    for n in labels
        vid((:x, n))
    end
    for c in curves
        c.isloop && isempty(c.word) && vid((:e, c.name))
    end
    nv = length(vindex)
    parent = collect(1:max(nv, 1))
    nedge = 0
    function link(a, b)
        nedge += 1
        nv == 0 && return
        ra, rb = _blc_find(parent, a), _blc_find(parent, b)
        parent[ra] = rb
    end
    for b in boundaries
        pts = b.points
        for i in eachindex(pts)
            j = mod1(i + 1, length(pts))
            link(vindex[(:b, b.name, pts[i])], vindex[(:b, b.name, pts[j])])
        end
    end
    for c in curves
        if c.isloop && isempty(c.word)
            continue
        end
        verts = Tuple[]
        if !c.isloop
            push!(verts, (:b, _boundary_of(boundaries, c.startp), c.startp))
        end
        for lab in c.word
            push!(verts, (:x, lab.n))
        end
        if !c.isloop
            push!(verts, (:b, _boundary_of(boundaries, c.endp), c.endp))
        end
        if c.isloop
            for i in eachindex(verts)
                link(vindex[verts[i]], vindex[verts[mod1(i + 1, length(verts))]])
            end
        else
            for i in 1:(length(verts) - 1)
                link(vindex[verts[i]], vindex[verts[i + 1]])
            end
        end
    end
    V = length(labels) + sum(length(b.points) for b in boundaries; init=0)
    E = nedge
    twin = Dict{_BLCDart, _BLCDart}()
    function join_darts(x::_BLCDart, y::_BLCDart)
        twin[x] = y
        twin[y] = x
    end
    for b in boundaries
        pts = b.points
        for i in eachindex(pts)
            j = mod1(i + 1, length(pts))
            join_darts(_BLCDart(3, b.name, 0, 0, pts[i], 1), _BLCDart(3, b.name, 0, 0, pts[j], -1))
        end
    end
    for c in curves
        (c.isloop && isempty(c.word)) && continue
        seq = Tuple{Int, Int, String, Int}[]
        if !c.isloop
            push!(seq, (2, 0, c.startp, 1))
        end
        for lab in c.word
            push!(seq, (1, lab.n, "", lab.sign))
        end
        if !c.isloop
            push!(seq, (2, 0, c.endp, -1))
        end
        m = length(seq)
        span = c.isloop ? m : m - 1
        for i in 1:span
            j = c.isloop ? mod1(i + 1, m) : i + 1
            a, b = seq[i], seq[j]
            left = a[1] == 1 ? _BLCDart(1, c.name, a[2], a[4], "", 1) :
                               _BLCDart(2, c.name, 0, 0, a[3], a[4])
            right = b[1] == 1 ? _BLCDart(1, c.name, b[2], b[4], "", -1) :
                                _BLCDart(2, c.name, 0, 0, b[3], b[4])
            join_darts(left, right)
        end
    end
    groups = Dict{Tuple, Vector{_BLCDart}}()
    for d in keys(twin)
        key = d.kind == 1 ? (:x, d.label) : (:b, d.kind == 3 ? d.name : _boundary_of(boundaries, d.port), d.kind == 3 ? d.port : d.port)
        # boundary name for an arc-end dart:
        if d.kind == 2
            key = (:b, _boundary_of(boundaries, d.port), d.port)
        elseif d.kind == 3
            key = (:b, d.name, d.port)
        end
        push!(get!(groups, key, _BLCDart[]), d)
    end
    rho = Dict{_BLCDart, _BLCDart}()
    for (key, darts) in groups
        ordered = _blc_order_vertex(key, darts, own)
        for i in eachindex(ordered)
            rho[ordered[i]] = ordered[mod1(i + 1, length(ordered))]
        end
    end
    length(rho) == length(twin) || _blc_fail("4.2", :rotation, "darts $(length(twin)) rho $(length(rho))")
    length(twin) == 2 * E || _blc_fail("4.2", :rotation, "darts $(length(twin)) edges $(2E)")
    unseen = Set(keys(twin))
    faces = Vector{_BLCDart}[]
    while !isempty(unseen)
        start = _blc_min_dart(unseen)
        z = start
        face = _BLCDart[]
        guard = 0
        while true
            z in unseen || _blc_fail("4.2", :face, string(z))
            delete!(unseen, z)
            push!(face, z)
            z = rho[twin[z]]
            guard += 1
            guard > length(twin) + 1 && _blc_fail("4.2", :face, "did not close")
            z == start && break
        end
        push!(faces, face)
    end
    pure = [face for face in faces if all(d -> d.kind == 3, face)]
    interior = [face for face in faces if any(d -> d.kind != 3, face)]
    curve_edges = E - sum(length(b.points) for b in boundaries; init=0)
    roots = nv == 0 ? 0 : length(unique(_blc_find(parent, i) for i in 1:nv))
    connected = roots <= 1
    use_disk = curve_edges == 0 && length(boundaries) == 1 && connected
    boundary_faces = pure
    F = length(interior)
    if use_disk
        forward = [face for face in pure if all(d -> d.dir == 1, face)]
        reverse = [face for face in pure if all(d -> d.dir == -1, face)]
        F = length(interior) + length(forward)
        boundary_faces = reverse
    end
    ports = Vector{String}[]
    assigned = Dict{String, Vector{String}}()
    for face in boundary_faces
        bname = ""
        seq = String[]
        for d in face
            bname = d.name
            if isempty(seq) || seq[end] != d.port
                push!(seq, d.port)
            end
        end
        if length(seq) >= 2 && seq[1] == seq[end]
            pop!(seq)
        end
        assigned[bname] = seq
    end
    for b in boundaries
        push!(ports, get(assigned, b.name, String[]))
    end
    chi = V - E + F
    counted = interior
    if use_disk
        append!(counted, [face for face in pure if all(d -> d.dir == 1, face)])
    end
    (
        diagnosis = BorderedLCDiagnosis(V, E, F, length(boundary_faces), chi, connected, ports),
        counted = counted,
        twin = twin,
    )
end

function _boundary_of(boundaries, port)
    hits = [b.name for b in boundaries if port in b.points]
    length(hits) == 1 || _blc_fail("4.1", :port, "port $port is on $(length(hits)) boundaries")
    hits[1]
end

function _blc_order_vertex(key, darts, own)
    if key[1] == :x
        label = key[2]
        minus = own[(label, -1)]
        plus = own[(label, 1)]
        wanted = [
            _BLCDart(1, minus, label, -1, "", 1),
            _BLCDart(1, plus, label, 1, "", 1),
            _BLCDart(1, minus, label, -1, "", -1),
            _BLCDart(1, plus, label, 1, "", -1),
        ]
        for d in wanted
            d in darts || _blc_fail("4.2", :rotation, "missing dart at $label")
        end
        return wanted
    end
    port = key[3]
    outb = _BLCDart(3, key[2], 0, 0, port, 1)
    inb = _BLCDart(3, key[2], 0, 0, port, -1)
    interior = [d for d in darts if d.kind == 2]
    length(interior) <= 1 || _blc_fail("4.1", :port, "port $port has $(length(interior)) arcs")
    ordered = [outb]
    if length(interior) == 1
        push!(ordered, interior[1])
    end
    push!(ordered, inb)
    for d in ordered
        d in darts || _blc_fail("4.2", :rotation, "missing boundary dart $port")
    end
    length(ordered) == length(darts) || _blc_fail("4.2", :rotation, "extra darts at $port")
    ordered
end

function _blc_min_dart(unseen)
    best = nothing
    for d in unseen
        if best === nothing || _blc_dart_key(d) < _blc_dart_key(best)
            best = d
        end
    end
    best
end

_blc_dart_key(d::_BLCDart) = (d.kind, d.name, d.label, d.sign, d.port, d.dir)

function _blc_selected(piece::BorderedLickorishCode; reference::Bool, drop::Set{String})
    curves = [c for c in _blc_curves(piece) if !(c.name in drop)]
    if !reference
        return curves
    end
    own = Dict{Tuple{Int, Int}, String}()
    refnames = Set(c.name for c in curves if c.ref)
    for c in _blc_curves(piece)
        for lab in c.word
            own[(lab.n, lab.sign)] = c.name
        end
    end
    out = _BLCCurve[]
    for c in curves
        c.ref || continue
        word = BLCLabel[]
        for lab in c.word
            other = lab.sign == 1 ? -1 : 1
            haskey(own, (lab.n, other)) || continue
            own[(lab.n, other)] in refnames && own[(lab.n, lab.sign)] in refnames &&
                push!(word, lab)
        end
        push!(out, _BLCCurve(c.name, true, c.isloop, c.startp, c.endp, word))
    end
    out
end

function diagnose_bordered_lc(piece::BorderedLickorishCode; reference::Bool=false)
    drop = Set{String}()
    if !reference
        for place in piece.places
            for name in place.components
                push!(drop, name)
            end
        end
    end
    curves = _blc_selected(piece; reference=reference, drop=drop)
    _blc_cellulation(curves, piece.boundaries).diagnosis
end

function _blc_check_endpoints(piece::BorderedLickorishCode)
    used = Dict{String, String}()
    for arc in Iterators.flatten((piece.k_arcs, piece.q_arcs, piece.arcs))
        for port in (arc.start, arc.stop)
            haskey(used, port) && _blc_fail("4.1", :port, "port $port is shared")
            used[port] = arc.name
            _boundary_of(piece.boundaries, port)
        end
        arc.start == arc.stop && _blc_fail("4.1", :port, "$(arc.name) has equal ends")
    end
    listed = Set{String}()
    for b in piece.boundaries
        for p in b.points
            p in listed && _blc_fail("4.1", :port, "point $p repeats")
            push!(listed, p)
        end
    end
    for port in keys(used)
        port in listed || _blc_fail("4.1", :port, "port $port is not on a boundary")
    end
end

function _blc_floating(piece::BorderedLickorishCode)
    inputs = _BLCCurve[]
    refs = Set{String}()
    for c in _blc_curves(piece)
        c.ref ? push!(refs, c.name) : push!(inputs, c)
    end
    own = Dict{Tuple{Int, Int}, String}()
    for c in _blc_curves(piece)
        for lab in c.word
            own[(lab.n, lab.sign)] = c.name
        end
    end
    parent = Dict(c.name => c.name for c in inputs)
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    touches = Set{String}()
    for c in inputs
        if !c.isloop
            push!(touches, c.name)
        end
        for lab in c.word
            other = own[(lab.n, lab.sign == 1 ? -1 : 1)]
            if other in refs
                push!(touches, c.name)
            elseif any(d.name == other for d in inputs)
                ra, rb = find(c.name), find(other)
                parent[ra] = rb
            end
        end
    end
    groups = Dict{String, Vector{String}}()
    for c in inputs
        find(c.name) in touches && continue
        # a component touches if any member of its class touches
        push!(get!(groups, find(c.name), String[]), c.name)
    end
    floating = Vector{String}[]
    touched_roots = Set(find(name) for name in touches)
    for (root, members) in groups
        root in touched_roots && continue
        push!(floating, sort(members))
    end
    sort!(floating; by=first)
    floating
end

function _blc_place_token_ok(token::String, face_count::Int, components::Set{String})
    if (m = match(r"^F([1-9]\d*)$", token)) !== nothing
        n = parse(Int, m.captures[1])
        1 <= n <= face_count || _blc_fail("4.5", :host, token)
        return
    end
    if (m = match(r"^face\(([^/]+)/([+-]\d+)/out\)$", token)) !== nothing
        m.captures[1] in components || _blc_fail("4.5", :host, token)
        return
    end
    if (m = match(r"^(.+)([+-])$", token)) !== nothing && m.captures[1] in components
        return
    end
    _blc_fail("4.5", :host, token)
end

function _blc_check_places(piece::BorderedLickorishCode, ::BorderedLCDiagnosis)
    floating = _blc_floating(piece)
    declared = [sort(collect(place.components)) for place in piece.places]
    sort!(declared; by=first)
    floating == declared || _blc_fail("4.5", :place,
        "floating $(floating) place $(declared)")
    isempty(piece.places) && return
    names = Set(_blc_names(piece))
    parent = Dict{Int, Int}()
    comp_group = Dict{String, Int}()
    for (i, place) in enumerate(piece.places)
        parent[i] = i
        for name in place.components
            name in names || _blc_fail("4.5", :place, name)
            haskey(comp_group, name) && _blc_fail("4.5", :place, "repeat $name")
            comp_group[name] = i
        end
    end
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    function link(a, b)
        ra, rb = find(a), find(b)
        ra == rb && _blc_fail("4.5", :place, "cycle")
        parent[ra] = rb
    end
    for (i, place) in enumerate(piece.places)
        drop = Set{String}()
        for later in piece.places[i:end]
            for name in later.components
                push!(drop, name)
            end
        end
        faces = bordered_lc_faces(piece; drop=drop)
        ids = Set(face.id for face in faces)
        members = Set(place.components)
        _blc_accept_token(piece, place.host, drop, ids, members, true)
        _blc_accept_token(piece, place.outside, drop, ids, members, false)
        for token in (place.host, place.outside)
            m = match(r"^(.+)([+-])$", token)
            m === nothing && continue
            startswith(token, "face(") && continue
            cname = String(m.captures[1])
            haskey(comp_group, cname) || continue
            other = comp_group[cname]
            other == i || link(i, other)
        end
    end
end

function validate_bordered_lc(piece::BorderedLickorishCode)
    _blc_check_endpoints(piece)
    for c in _blc_curves(piece)
        for lab in c.word
            lab.n > 0 || _blc_fail("4.1", :label, _blc_fmt(lab))
        end
    end
    diagnose_bordered_lc(piece; reference=false)
    reference = diagnose_bordered_lc(piece; reference=true)
    expected = 2 - 2 * piece.genus - piece.nboundaries
    problems = String[]
    reference.connected || push!(problems, "disconnected")
    reference.boundary_cycles == piece.nboundaries ||
        push!(problems, "boundary_cycles=$(reference.boundary_cycles)")
    reference.chi == expected || push!(problems, "chi=$(reference.chi)")
    isempty(problems) || _blc_fail("4.3", :reference_graph, join(problems, " "))
    full = diagnose_bordered_lc(piece; reference=false)
    full.chi == reference.chi || _blc_fail("4.4", :input_chi,
        "full chi $(full.chi), reference chi $(reference.chi)")
    _blc_check_per_face(piece)
    _blc_check_places(piece, full)
    true
end

function _blc_is_k(name::AbstractString, genus::Int)
    m = match(r"^([lmn])(\d+)$", name)
    m === nothing && return false
    i = parse(Int, m.captures[2])
    letter = m.captures[1]
    letter == "l" && return 1 <= i <= genus
    letter == "m" && return 1 <= i <= genus
    letter == "n" && return 1 <= i <= genus - 1
    false
end

function _blc_comp_key(name::AbstractString)
    m = match(r"^([lmn])(\d+)$", name)
    m === nothing && return (9, name, 0)
    letter = m.captures[1]
    idx = letter == "l" ? 1 : letter == "m" ? 2 : 3
    (idx, "", parse(Int, m.captures[2]))
end

function cut_bordered_lc(lc::LickorishCode, curve::AbstractString; parent::AbstractString="P")
    haskey(lc.words, curve) || _blc_fail("6.1", :curve, curve)
    cutter_word = lc.words[curve]
    isempty(cutter_word) && _blc_fail("6.2", :separating, "$curve has no crossing")
    seen = Dict{Int, Vector{Tuple{String, Int}}}()
    for (name, word) in lc.words
        for (i, x) in enumerate(word)
            push!(get!(seen, abs(x), Tuple{String, Int}[]), (name, x))
        end
    end
    cut_labels = Set{Int}()
    for x in cutter_word
        lab = abs(x)
        push!(cut_labels, lab)
        sites = get(seen, lab, Tuple{String, Int}[])
        length(sites) == 2 || _blc_fail("6.1", :label, "label $lab appears $(length(sites)) times")
        owners = Set(site[1] for site in sites)
        curve in owners || _blc_fail("6.1", :label, "label $lab is not on $curve")
        length(owners) == 2 || _blc_fail("6.2", :not_simple, "$curve crosses itself at $lab")
    end
    k_arcs = BLCArc[]
    k_loops = BLCLoop[]
    arcs = BLCArc[]
    loops = BLCLoop[]
    for name in sort!(collect(keys(lc.words)); by=_blc_comp_key)
        name == curve && continue
        word = lc.words[name]
        kind = _blc_is_k(name, lc.genus) ? :K : :input
        hits = [i for i in eachindex(word) if abs(word[i]) in cut_labels]
        origin_start = isempty(word) ? nothing : _blc_from_ints([word[1]])[1]
        if isempty(hits)
            isempty(word) && continue
            origin = BLCOrigin(parent, name, origin_start, nothing)
            loop = BLCLoop(kind, name, _blc_from_ints(word), origin)
            kind == :K ? push!(k_loops, loop) : push!(loops, loop)
            continue
        end
        n = length(word)
        for t in eachindex(hits)
            i = hits[t]
            j = hits[mod1(t + 1, length(hits))]
            frag = Int[]
            p = mod1(i + 1, n)
            if p != j
                while p != j
                    push!(frag, word[p])
                    p = mod1(p + 1, n)
                end
            end
            start_shore = word[i] > 0 ? "L" : "R"
            end_shore = word[j] > 0 ? "R" : "L"
            start_port = string(abs(word[i]), start_shore)
            end_port = string(abs(word[j]), end_shore)
            origin = BLCOrigin(parent, name, origin_start, _blc_from_ints([word[i]])[1])
            arc = BLCArc(kind, name, start_port, end_port, _blc_from_ints(frag), origin)
            kind == :K ? push!(k_arcs, arc) : push!(arcs, arc)
        end
    end
    sort!(k_arcs; by=a -> (a.origin.cut.n, a.name))
    sort!(arcs; by=a -> (a.origin.cut.n, a.name))
    sort!(k_loops; by=a -> _blc_comp_key(a.name))
    sort!(loops; by=a -> _blc_comp_key(a.name))
    left = String["oL"]
    right_cross = String[]
    for x in cutter_word
        push!(left, string(abs(x), "L"))
        push!(right_cross, string(abs(x), "R"))
    end
    right = String["oR"; reverse(right_cross)]
    boundaries = [
        BLCBoundary("BL", "oL", left),
        BLCBoundary("BR", "oR", right),
    ]
    piece = BorderedLickorishCode("cut-$curve", lc.genus - 1, 2, boundaries,
        k_arcs, k_loops, BLCArc[], BLCLoop[], arcs, loops, BLCPlace[])
    full = diagnose_bordered_lc(piece; reference=false)
    full.connected || _blc_fail("6.2", :separating, curve)
    BorderedCut(piece, String(curve), lc.genus,
        copy(lc.places), copy(lc.disks), copy(lc.points), copy(lc.weights))
end

function _blc_port_label(port::String)
    m = match(r"^(\d+)([LR])$", port)
    m === nothing && _blc_fail("6.3", :port, port)
    parse(Int, m.captures[1]), m.captures[2]
end

function _blc_rotate(word::Vector{Int}, start::Int)
    i = findfirst(==(start), word)
    i === nothing && _blc_fail("6.3", :origin, "start $start not in $(word)")
    [word[i:end]; word[1:i - 1]]
end

function glue_bordered_lc(cut::BorderedCut)
    piece = cut.piece
    length(piece.boundaries) == 2 || _blc_fail("6.3", :boundary, "need two shores")
    left = piece.boundaries[1]
    right = piece.boundaries[2]
    left.base == "oL" && right.base == "oR" || _blc_fail("6.3", :boundary, "shores")
    left_labs = Int[_blc_port_label(p)[1] for p in left.points if p != "oL"]
    right_labs = Int[_blc_port_label(p)[1] for p in right.points if p != "oR"]
    right_labs == reverse(left_labs) || _blc_fail("6.3", :pairing, "shores disagree")
    words = Dict{String, Vector{Int}}()
    groups = Dict{String, Vector{BLCArc}}()
    for arc in Iterators.flatten((piece.k_arcs, piece.q_arcs, piece.arcs))
        push!(get!(groups, arc.name, BLCArc[]), arc)
    end
    for (name, arcs) in groups
        by = Dict{Int, BLCArc}()
        for arc in arcs
            arc.origin === nothing && _blc_fail("6.3", :origin, name)
            by[arc.origin.cut.n] = arc
        end
        seq = Int[]
        arc = arcs[1]
        seen = Set{Int}()
        for _ in eachindex(arcs)
            cut_at = arc.origin.cut.n
            cut_at in seen && _blc_fail("6.3", :arc, name)
            push!(seen, cut_at)
            append!(seq, _blc_ints(arc.word))
            lab, shore = _blc_port_label(arc.stop)
            push!(seq, (shore == "R" ? 1 : -1) * lab)
            arc = by[lab]
        end
        arc.origin.cut.n == arcs[1].origin.cut.n || _blc_fail("6.3", :arc, name)
        start = arcs[1].origin.start
        words[name] = _blc_rotate(seq, start.sign * start.n)
    end
    for loop in Iterators.flatten((piece.k_loops, piece.q_loops, piece.loops))
        words[loop.name] = _blc_ints(loop.word)
    end
    all_arcs = BLCArc[]
    for arcs in values(groups)
        append!(all_arcs, arcs)
    end
    cutter = Int[]
    for lab in left_labs
        sign = 0
        for arc in all_arcs
            if arc.start == string(lab, "L")
                sign = 1
            elseif arc.start == string(lab, "R")
                sign = -1
            end
        end
        sign == 0 && _blc_fail("6.3", :pairing, "no arc leaves $lab")
        push!(cutter, -sign * lab)
    end
    words[cut.cutter] = cutter
    LickorishCode(cut.genus, words, cut.places, cut.disks, cut.points, cut.weights, nothing)
end

struct BLCFace
    id::Int
    ports::Vector{String}
    by_label::Bool
end

Base.:(==)(a::BLCFace, b::BLCFace) = a.id == b.id && a.ports == b.ports && a.by_label == b.by_label
Base.hash(a::BLCFace, h::UInt) = hash((a.id, a.ports, a.by_label), h)

struct BLCPair
    name::String
    left::BorderedLickorishCode
    right::BorderedLickorishCode
    corner_labels::Vector{Int}
    pin_below::Vector{Int}
    above_k::Bool
    arc_order::Vector{String}
end

Base.:(==)(a::BLCPair, b::BLCPair) =
    a.name == b.name && a.left == b.left && a.right == b.right &&
    a.corner_labels == b.corner_labels && a.pin_below == b.pin_below &&
    a.above_k == b.above_k && a.arc_order == b.arc_order
Base.hash(a::BLCPair, h::UInt) =
    hash((a.name, a.left, a.right, a.corner_labels, a.pin_below, a.above_k, a.arc_order), h)

struct BLCPairHit
    bindings::Vector{Pair{String,String}}
    starts::Vector{Int}
    labels::Vector{Int}
    output::LickorishCode
end

function _blc_is_empty_loop(piece::BorderedLickorishCode, name::AbstractString)
    for loop in Iterators.flatten((piece.k_loops, piece.q_loops, piece.loops))
        loop.name == name && return isempty(loop.word)
    end
    false
end

function _blc_passage(piece::BorderedLickorishCode, drop::Set{String}, cname::AbstractString, sign::Int, n::Int)
    cname in drop && return false
    for c in _blc_curves(piece)
        c.name == cname || continue
        for lab in c.word
            lab.sign == sign && lab.n == n && return true
        end
    end
    false
end

function _blc_accept_token(piece, token, drop, ids, members, is_host)
    if (m = match(r"^F([1-9]\d*)$", token)) !== nothing
        is_host || _blc_fail("4.5", :host, token)
        parse(Int, m.captures[1]) in ids || _blc_fail("4.5", :host, token)
        return
    end
    if (m = match(r"^face\(([^/]+)/([+-])(\d+)/out\)$", token)) !== nothing
        cname = String(m.captures[1])
        sign = m.captures[2] == "-" ? -1 : 1
        n = parse(Int, m.captures[3])
        if is_host
            _blc_passage(piece, drop, cname, sign, n) || _blc_fail("4.5", :host, token)
        else
            (cname in members && !_blc_is_empty_loop(piece, cname) &&
                _blc_passage(piece, Set{String}(), cname, sign, n)) ||
                _blc_fail("4.5", :host, token)
        end
        return
    end
    if (m = match(r"^(.+)([+-])$", token)) !== nothing
        cname = String(m.captures[1])
        if is_host
            (cname in drop || !_blc_is_empty_loop(piece, cname)) && _blc_fail("4.5", :host, token)
        else
            (cname in members && _blc_is_empty_loop(piece, cname)) || _blc_fail("4.5", :host, token)
        end
        return
    end
    _blc_fail("4.5", :host, token)
end

function _blc_isolated_chi(piece::BorderedLickorishCode, group::Vector{String})
    wanted = Set(group)
    curves = [c for c in _blc_curves(piece) if c.name in wanted]
    if all(c.isloop && isempty(c.word) for c in curves)
        return 2
    end
    _blc_cellulation(curves, BLCBoundary[]).diagnosis.chi
end

_blc_dart_vertex(d::_BLCDart) = d.kind == 1 ? (:x, d.label) : (:b, d.port)

function _blc_edge_key(d::_BLCDart, tw::_BLCDart)
    a, b = _blc_dart_key(d), _blc_dart_key(tw)
    a <= b ? (a, b) : (b, a)
end

function _blc_region_chi(cell, face_ids::Vector{Int})
    verts = Set{Tuple}()
    edges = Set{Tuple}()
    for i in face_ids
        for d in cell.counted[i]
            tw = cell.twin[d]
            push!(verts, _blc_dart_vertex(d))
            push!(verts, _blc_dart_vertex(tw))
            push!(edges, _blc_edge_key(d, tw))
        end
    end
    length(verts) - length(edges) + length(face_ids)
end

function _blc_input_edge(d::_BLCDart, tw::_BLCDart, refnames::Set{String})
    for x in (d, tw)
        x.kind == 3 && return false
        x.name in refnames && return false
    end
    true
end

function _blc_place_drop(piece::BorderedLickorishCode)
    drop = Set{String}()
    for place in piece.places
        for name in place.components
            push!(drop, name)
        end
    end
    drop
end

function _blc_check_per_face(piece::BorderedLickorishCode)
    for group in _blc_floating(piece)
        chi = _blc_isolated_chi(piece, group)
        chi == 2 || _blc_fail("4.4", :isolated, "isolated $(join(group, " ")) chi $chi")
    end
    drop = _blc_place_drop(piece)
    full_curves = _blc_selected(piece; reference=false, drop=drop)
    ref_curves = _blc_selected(piece; reference=true, drop=drop)
    # No reference curve: §4.3 counts the forward boundary as the one
    # reference face. Full interior faces meet the opposite boundary darts
    # (dir differs, and the dart key includes dir), so the seed below stays
    # empty and line "reference face is empty" fires. Spec §7.7.1 says this
    # input passes §4.4 because the two cellulations have the same chi,
    # which validate_bordered_lc already compared. Isolated components were
    # checked above.
    if isempty(ref_curves)
        return
    end
    full = _blc_cellulation(full_curves, piece.boundaries)
    ref = _blc_cellulation(ref_curves, piece.boundaries)
    refnames = Set(c.name for c in ref_curves)
    ref_of = Dict{Any, Int}()
    for (i, face) in enumerate(ref.counted)
        for d in face
            ref_of[_blc_dart_key(d)] = i
        end
    end
    face_of = Dict{Any, Int}()
    for (i, face) in enumerate(full.counted)
        for d in face
            face_of[_blc_dart_key(d)] = i
        end
    end
    seed = Dict{Int, Int}()
    for (i, face) in enumerate(full.counted)
        curve_hits = Int[]
        boundary_hits = Int[]
        for d in face
            k = _blc_dart_key(d)
            haskey(ref_of, k) || continue
            push!(d.kind == 3 ? boundary_hits : curve_hits, ref_of[k])
        end
        hits = isempty(curve_hits) ? boundary_hits : curve_hits
        isempty(hits) && continue
        chosen = unique(hits)
        length(chosen) == 1 || _blc_fail("4.4", :disk,
            "full face $i meets reference faces $(chosen)")
        seed[i] = chosen[1]
    end
    n = length(full.counted)
    parent = collect(1:max(n, 1))
    function find(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        x
    end
    changed = true
    while changed
        changed = false
        for i in 1:n
            for d in full.counted[i]
                tw = full.twin[d]
                _blc_input_edge(d, tw, refnames) || continue
                j = get(face_of, _blc_dart_key(tw), 0)
                (j == 0 || j == i) && continue
                if haskey(seed, i) && haskey(seed, j)
                    seed[i] == seed[j] || _blc_fail("4.4", :disk,
                        "reference face $(seed[i]) meets $(seed[j])")
                elseif haskey(seed, i)
                    seed[j] = seed[i]
                    changed = true
                elseif haskey(seed, j)
                    seed[i] = seed[j]
                    changed = true
                else
                    find(i) == find(j) || (parent[find(i)] = find(j); changed = true)
                end
            end
        end
    end
    groups = Dict{Int, Vector{Int}}()
    for i in 1:n
        haskey(seed, i) && push!(get!(groups, seed[i], Int[]), i)
        haskey(seed, i) && continue
        r = find(i)
        donor = 0
        for j in 1:n
            find(j) == r && haskey(seed, j) && (donor = seed[j])
        end
        donor == 0 && continue
        push!(get!(groups, donor, Int[]), i)
    end
    for ri in eachindex(ref.counted)
        faces = unique(get(groups, ri, Int[]))
        isempty(faces) && _blc_fail("4.4", :disk, "reference face $ri is empty")
        # A reference face may repeat a vertex on its boundary walk, so the
        # carrier's V-E+F is the surface piece, not 1. Filling that face by
        # disks keeps this number; the disk itself is the difference plus one.
        carrier = _blc_region_chi(ref, [ri])
        chi = _blc_region_chi(full, faces) - carrier + 1
        chi == 1 || _blc_fail("4.4", :disk, "reference face $ri chi $chi")
    end
end

function _blc_port_key(face::Vector{_BLCDart})
    ports = String[]
    for d in face
        isempty(d.port) || push!(ports, d.port)
    end
    isempty(ports) && return String[]
    best = ports
    for i in 2:length(ports)
        rot = vcat(ports[i:end], ports[1:i - 1])
        rot < best && (best = rot)
    end
    best
end

function _blc_face_tie(face::Vector{_BLCDart})
    best = _blc_dart_key(face[1])
    for d in face
        k = _blc_dart_key(d)
        k < best && (best = k)
    end
    best
end

function bordered_lc_faces(piece::BorderedLickorishCode; drop::Union{Nothing,Set{String}}=nothing)
    if drop === nothing
        drop = _blc_place_drop(piece)
    end
    curves = _blc_selected(piece; reference=false, drop=drop)
    cell = _blc_cellulation(curves, piece.boundaries)
    n = length(cell.counted)
    face_of = Dict{Any, Int}()
    for (i, face) in enumerate(cell.counted)
        for d in face
            face_of[_blc_dart_key(d)] = i
        end
    end
    own = Dict{Tuple{Int,Int}, String}()
    labels = Int[]
    for c in curves
        for lab in c.word
            own[(lab.n, lab.sign)] = c.name
            push!(labels, lab.n)
        end
    end
    sort!(unique!(labels))
    ids = zeros(Int, n)
    by_label = falses(n)
    next_id = 1
    for label in labels
        for (sign, dir) in ((-1, 1), (1, 1), (-1, -1), (1, -1))
            haskey(own, (label, sign)) || continue
            dart = _BLCDart(1, own[(label, sign)], label, sign, "", dir)
            i = get(face_of, _blc_dart_key(dart), 0)
            (i == 0 || ids[i] != 0) && continue
            ids[i] = next_id
            by_label[i] = true
            next_id += 1
        end
    end
    pending = [(i, _blc_port_key(cell.counted[i]), _blc_face_tie(cell.counted[i]))
               for i in 1:n if ids[i] == 0]
    sort!(pending; by=row -> (row[2], row[3]))
    for row in pending
        ids[row[1]] = next_id
        next_id += 1
    end
    out = BLCFace[]
    for i in 1:n
        push!(out, BLCFace(ids[i], _blc_port_key(cell.counted[i]), by_label[i]))
    end
    sort!(out; by=face -> face.id)
    out
end

function _blc_arc_table(piece::BorderedLickorishCode)
    table = Dict{String, BLCArc}()
    for arc in Iterators.flatten((piece.k_arcs, piece.q_arcs, piece.arcs))
        table[arc.name] = arc
    end
    table
end

function _blc_make_pair(name, left, right, corners, pin, above, order)
    lt = _blc_arc_table(left)
    rt = _blc_arc_table(right)
    for arc in order
        (haskey(lt, arc) && haskey(rt, arc)) || _blc_fail("7", :pair, arc)
        lt[arc].kind == rt[arc].kind || _blc_fail("7", :pair, arc)
    end
    BLCPair(name, left, right, corners, pin, above, String.(order))
end

function vertex_bordered_pair()
    left = read_bordered_lc("""
piece vertex-L : genus=0 boundaries=1
marking : reference
boundary d base=Bm : Bm Am Cp Bp Ap Cm
K arc A Am -> Ap : [+2 +1]
K arc B Bm -> Bp : [+3 -1]
arc C Cm -> Cp : [-3 -2]
""")
    right = read_bordered_lc("""
piece vertex-R : genus=0 boundaries=1
marking : reference
boundary d base=Bm : Bm Am Cp Bp Ap Cm
K arc A Am -> Ap : [+1 +2]
K arc B Bm -> Bp : [-1 +3]
arc C Cm -> Cp : [-2 -3]
""")
    validate_bordered_lc(left)
    validate_bordered_lc(right)
    _blc_make_pair("vertex", left, right, [1, 2, 3], [1], false, ["A", "B", "C"])
end

function digon_bordered_pair()
    left = read_bordered_lc("""
piece digon-L : genus=0 boundaries=1
marking : reference
boundary d base=Cm : Cm Km Kp Cp
K arc K Km -> Kp : [+1 -2]
arc C Cm -> Cp : [-1 +2]
""")
    right = read_bordered_lc("""
piece digon-R : genus=0 boundaries=1
marking : reference
boundary d base=Cm : Cm Km Kp Cp
K arc K Km -> Kp : []
arc C Cm -> Cp : []
""")
    validate_bordered_lc(left)
    validate_bordered_lc(right)
    _blc_make_pair("digon", left, right, [1, 2], Int[], true, ["C", "K"])
end

function swap_bordered_pair(pair::BLCPair)
    BLCPair(pair.name, pair.right, pair.left, pair.corner_labels, pair.pin_below,
        pair.above_k, pair.arc_order)
end

function _blc_windows(word::Vector{Int}, pattern::Vector{Int}, subst::Dict{Int,Int})
    hits = Tuple{Int,Dict{Int,Int}}[]
    (isempty(word) || isempty(pattern)) && return hits
    n = length(word)
    taken = Dict(v => k for (k, v) in subst)
    for start in 1:n
        next = copy(subst)
        used = copy(taken)
        ok = true
        for offset in 0:(length(pattern) - 1)
            got = word[mod1(start + offset, n)]
            want = pattern[offset + 1]
            if sign(got) != sign(want)
                ok = false
                break
            end
            local_label = abs(want)
            global_label = abs(got)
            if get(next, local_label, global_label) != global_label
                ok = false
                break
            end
            if haskey(used, global_label) && used[global_label] != local_label
                ok = false
                break
            end
            next[local_label] = global_label
            used[global_label] = local_label
        end
        ok && push!(hits, (start, next))
    end
    hits
end

function _blc_subst_word(pattern::Vector{Int}, subst::Dict{Int,Int})
    out = Int[]
    for entry in pattern
        haskey(subst, abs(entry)) || return nothing
        push!(out, sign(entry) * subst[abs(entry)])
    end
    out
end

function _blc_occurrences(word::Vector{Int}, pattern::Vector{Int})
    (isempty(pattern) || length(pattern) > length(word)) && return Int[]
    n = length(word)
    found = Int[]
    for start in 1:n
        all(offset -> word[mod1(start + offset, n)] == pattern[offset + 1], 0:(length(pattern) - 1)) ||
            continue
        push!(found, start)
    end
    found
end

function _blc_is_rotation(seq::Vector{Int}, target::Vector{Int})
    length(seq) == length(target) || return false
    n = length(seq)
    for shift in 0:(n - 1)
        all(j -> seq[mod1(j + shift, n)] == target[j], 1:n) && return true
    end
    false
end

function _blc_face_ok(lc::LickorishCode, corner_globals::Vector{Int})
    active = Dict(k => v for (k, v) in lc.words if !isempty(v))
    faces, _ = _lc_faces(active)
    n = length(corner_globals)
    placed = Set(values(lc.places))
    for face in faces
        length(face.corners) == n || continue
        face.id in placed && continue
        seq = Int[corner[1] for corner in face.corners]
        _blc_is_rotation(seq, corner_globals) || continue
        return true
    end
    false
end

function _blc_locations_fit(locations, lengths, word_lengths)
    occupied = Set{Tuple{Int,Int}}()
    for (i, (comp, start)) in enumerate(locations)
        n = word_lengths[comp]
        n == 0 && return false
        for offset in 0:(lengths[i] - 1)
            key = (comp, mod1(start + offset, n))
            key in occupied && return false
            push!(occupied, key)
        end
    end
    true
end

function _blc_location_products(options, lengths, word_lengths)
    found = Vector{Vector{Tuple{Int,Int}}}()
    current = Tuple{Int,Int}[]
    function visit(i)
        if i > length(options)
            _blc_locations_fit(current, lengths, word_lengths) && push!(found, copy(current))
            return
        end
        for loc in options[i]
            push!(current, loc)
            visit(i + 1)
            pop!(current)
        end
    end
    visit(1)
    found
end

function _blc_apply_at(lc::LickorishCode, names, before_frags, after_frags, locations, vertices)
    words = Dict(name => copy(get(lc.words, name, Int[])) for name in names)
    capacity = collect(1:length(before_frags))
    result, old_to_new, _, _, _, _, _ = _kernel_apply(
        names, words, before_frags, after_frags, locations, vertices, capacity)
    removed = setdiff(vertices, keys(old_to_new))
    _pin_ok(lc, old_to_new, removed)
    preview = Dict(names[i] => result[i] for i in eachindex(names))
    assigned = Dict{String,String}()
    emptied = String[]
    for (i, name) in enumerate(names)
        isempty(result[i]) || continue
        haskey(lc.places, name) && continue
        isempty(get(lc.words, name, Int[])) && continue
        push!(emptied, name)
    end
    if !isempty(emptied)
        face = _merged_place(lc, names, locations, Set(removed), old_to_new, preview)
        for name in emptied
            assigned[name] = face
        end
    end
    new_words, new_places, _ = _install_words(lc, names, result, lc.places, assigned)
    new_places = _move_transport_places(lc, new_words, new_places, old_to_new;
        digon_create_untouched_corners=Dict{String,Vector{Tuple{Int,String,String}}}())
    _build_output(lc, lc.genus, new_words, new_places)
end

function apply_bordered_pair(lc::LickorishCode, pair::BLCPair)
    validate_lc_input(lc)
    left = _blc_arc_table(pair.left)
    right = _blc_arc_table(pair.right)
    arcs = [left[name] for name in pair.arc_order]
    any(isempty(arc.word) for arc in arcs) && return BLCPairHit[]
    patterns = [_blc_ints(arc.word) for arc in arcs]
    after_frags = [_blc_ints(right[name].word) for name in pair.arc_order]
    max_local = maximum(abs(x) for word in patterns for x in word)
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    words = Dict(name => get(lc.words, name, Int[]) for name in names)
    word_lengths = [length(words[name]) for name in names]
    pools = Vector{Vector{String}}()
    for arc in arcs
        if arc.kind == :K
            push!(pools, [name for name in names if _blc_is_k(name, lc.genus)])
        else
            push!(pools, [name for name in names if !_blc_is_k(name, lc.genus)])
        end
    end
    limit = _k_limit(lc.genus)
    hits = BLCPairHit[]
    seen = Set{Tuple}()
    used = Set{String}()
    assignment = fill("", length(arcs))
    function leaf(subst)
        all(haskey(subst, k) for k in 1:max_local) || return
        if pair.above_k
            all(subst[k] > limit for k in 1:max_local) || return
        end
        for local_label in pair.pin_below
            subst[local_label] <= limit || return
        end
        corner = Int[subst[k] for k in pair.corner_labels]
        _blc_face_ok(lc, corner) || return
        concrete = Vector{Vector{Int}}()
        options = Vector{Vector{Tuple{Int,Int}}}()
        for i in eachindex(arcs)
            word = _blc_subst_word(patterns[i], subst)
            word === nothing && return
            push!(concrete, word)
            found = _blc_occurrences(words[assignment[i]], word)
            isempty(found) && return
            push!(options, [(index[assignment[i]], hit) for hit in found])
        end
        labels = Int[subst[k] for k in 1:max_local]
        curves = copy(assignment)
        products = _blc_location_products(options, [length(word) for word in concrete], word_lengths)
        for locations in products
            starts = [loc[2] for loc in locations]
            key = (tuple(curves...), tuple(starts...), tuple(labels...))
            key in seen && continue
            push!(seen, key)
            output = _blc_apply_at(lc, names, patterns, after_frags, locations, labels)
            bindings = Pair{String,String}[pair.arc_order[i] => curves[i] for i in eachindex(arcs)]
            push!(hits, BLCPairHit(bindings, starts, labels, output))
        end
    end
    function visit(i, subst)
        if i > length(arcs)
            leaf(subst)
            return
        end
        for curve in pools[i]
            curve in used && continue
            push!(used, curve)
            assignment[i] = curve
            for (_, next) in _blc_windows(words[curve], patterns[i], subst)
                visit(i + 1, next)
            end
            delete!(used, curve)
        end
    end
    visit(1, Dict{Int,Int}())
    hits
end

# §7.7.1. Rows are the shared digon_left_face table, with C pushing and D crossed.
const _BLC_DIGON_ROWS = Dict{Tuple{Int,Int},NamedTuple}(
    (1, 1) => (points=("Cm", "Cp", "Dm", "Dp"), c=[-1, 2], d=[-2, 1]),
    (1, -1) => (points=("Cm", "Cp", "Dp", "Dm"), c=[1, -2], d=[-1, 2]),
    (-1, 1) => (points=("Cm", "Dm", "Dp", "Cp"), c=[-1, 2], d=[1, -2]),
    (-1, -1) => (points=("Cm", "Dp", "Dm", "Cp"), c=[1, -2], d=[2, -1]),
)

function _blc_kind_ok(kind::Symbol)
    kind in (:K, :Q, :input)
end

function _blc_arc_head(kind::Symbol)
    kind == :K ? "K arc" : kind == :Q ? "Q arc" : "arc"
end

function _blc_signed_text(word::Vector{Int})
    join([entry > 0 ? "+" * string(entry) : string(entry) for entry in word], " ")
end

function _blc_host_matches(kind::Symbol, name::AbstractString, genus::Int)
    _blc_kind_ok(kind) || return false
    _blc_is_k(name, genus) == (kind == :K)
end

function _blc_digon_text(title::AbstractString, points, c_kind, d_kind, c_word, d_word)
    lines = String[
        "piece $title : genus=0 boundaries=1",
        "marking : reference",
        "boundary d base=Cm : " * join(points, " "),
    ]
    rows = [
        (c_kind == :K ? 4 : c_kind == :Q ? 6 : 8, "C", c_kind, c_word),
        (d_kind == :K ? 4 : d_kind == :Q ? 6 : 8, "D", d_kind, d_word),
    ]
    sort!(rows; by=row -> (row[1], row[2]))
    for (_, name, kind, word) in rows
        startp = name == "C" ? "Cm" : "Dm"
        stopp = name == "C" ? "Cp" : "Dp"
        push!(lines, _blc_arc_head(kind) * " $name $startp -> $stopp : [" * _blc_signed_text(word) * "]")
    end
    join(lines, "\n")
end

"""Removal pair for one sign row. Creation is this pair swapped."""
function digon_bordered_family(eps_c::Integer, eps_d::Integer;
        create::Bool=true, c_kind::Symbol=:input, d_kind::Symbol=:K)
    (eps_c in (1, -1) && eps_d in (1, -1)) ||
        _blc_fail("7.7.1", :sign, "signs must be +1 or -1")
    (_blc_kind_ok(c_kind) && _blc_kind_ok(d_kind)) ||
        _blc_fail("7.7.1", :kind, "component kind is not K, Q, or input")
    row = _BLC_DIGON_ROWS[(Int(eps_c), Int(eps_d))]
    left = read_bordered_lc(_blc_digon_text("digon-L", row.points, c_kind, d_kind, row.c, row.d))
    right = read_bordered_lc(_blc_digon_text("digon-R", row.points, c_kind, d_kind, Int[], Int[]))
    validate_bordered_lc(left)
    validate_bordered_lc(right)
    pair = _blc_make_pair("digon", left, right, [1, 2], Int[], true, ["C", "D"])
    create ? swap_bordered_pair(pair) : pair
end

function _blc_apply_insert(lc::LickorishCode, names, after_frags, locations, parking)
    before = [Int[] for _ in after_frags]
    words = Dict(name => copy(get(lc.words, name, Int[])) for name in names)
    capacity = collect(1:length(before))
    result, old_to_new, _, _, _, _, _ = _kernel_apply(
        names, words, before, after_frags, locations, Int[], capacity)
    _pin_ok(lc, old_to_new, Int[])
    new_words, new_places, _ = _install_words(lc, names, result, lc.places, Dict{String,String}())
    new_places = _move_transport_places(lc, new_words, new_places, old_to_new;
        digon_create_untouched_corners=parking)
    _build_output(lc, lc.genus, new_words, new_places)
end

function _blc_create_hits(lc::LickorishCode, pair::BLCPair, arcs, after_frags)
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    active = Dict(k => copy(v) for (k, v) in lc.words if !isempty(v))
    faces, _ = _lc_faces(active)
    occ = _lc_rotation_occ(active)
    hits = BLCPairHit[]
    seen = Set{Tuple}()
    for face in faces
        edges = _move_left_boundary(active, occ, face)
        for push_edge in edges, cross_edge in edges
            push_edge.component == cross_edge.component && continue
            _blc_host_matches(arcs[1].kind, push_edge.component, lc.genus) || continue
            _blc_host_matches(arcs[2].kind, cross_edge.component, lc.genus) || continue
            c_pat, d_pat = _derived_digon(push_edge.forward, cross_edge.forward)
            (c_pat == after_frags[1] && d_pat == after_frags[2]) || continue
            c_len = length(lc.words[push_edge.component])
            d_len = length(lc.words[cross_edge.component])
            c_at = push_edge.forward ? mod1(push_edge.position + 1, c_len) : push_edge.position
            d_at = cross_edge.forward ? mod1(cross_edge.position + 1, d_len) : cross_edge.position
            key = (push_edge.component, c_at, cross_edge.component, d_at)
            key in seen && continue
            push!(seen, key)
            corners = Tuple{Int,String,String}[
                (Int(c[1]), String(c[2]), String(c[3])) for c in face.corners]
            parking = Dict{String,Vector{Tuple{Int,String,String}}}(String(face.id) => corners)
            locations = [(index[push_edge.component], c_at), (index[cross_edge.component], d_at)]
            output = _blc_apply_insert(lc, names, after_frags, locations, parking)
            bindings = Pair{String,String}[
                pair.arc_order[1] => push_edge.component,
                pair.arc_order[2] => cross_edge.component,
            ]
            push!(hits, BLCPairHit(bindings, [c_at, d_at], Int[], output))
        end
    end
    hits
end

"""Apply a family value. Empty left words insert the right-hand digon at a face-edge site."""
function apply_bordered_family(lc::LickorishCode, pair::BLCPair)
    validate_lc_input(lc)
    left = _blc_arc_table(pair.left)
    right = _blc_arc_table(pair.right)
    arcs = BLCArc[]
    patterns = Vector{Int}[]
    after_frags = Vector{Int}[]
    for name in pair.arc_order
        (haskey(left, name) && haskey(right, name)) || _blc_fail("7.7", :pair, name)
        push!(arcs, left[name])
        push!(patterns, _blc_ints(left[name].word))
        push!(after_frags, _blc_ints(right[name].word))
    end
    if all(isempty, patterns)
        all(isempty, after_frags) && return BLCPairHit[]
        return _blc_create_hits(lc, pair, arcs, after_frags)
    end
    apply_bordered_pair(lc, pair)
end

function _blc_loop_text(name::AbstractString, word::Vector{Int})
    "loop $name : (" * _blc_signed_text(word) * ")"
end

function _blc_slide_text(title::AbstractString, q_word, r_word, first_name, first_word, second_word)
    join([
        "piece $title : genus=0 boundaries=3",
        "marking : reference",
        "boundary A base=oA : oA qA rA",
        "boundary B base=oB : oB qB",
        "boundary C base=oC : oC rC",
        "Q arc q qA -> qB : [" * _blc_signed_text(q_word) * "]",
        "Q arc r rA -> rC : [" * _blc_signed_text(r_word) * "]",
        _blc_loop_text(first_name, first_word),
        _blc_loop_text("b", second_word),
    ], "\n")
end

"""Pants pair for orientations of a and b. `eps_a` and `eps_b` are the signs on q and r."""
function slide_bordered_family(eps_a::Integer, eps_b::Integer)
    (eps_a in (1, -1) && eps_b in (1, -1)) ||
        _blc_fail("7.7.2", :sign, "signs must be +1 or -1")
    sa = Int(eps_a)
    sb = Int(eps_b)
    # Label 3 carries the same sign as label 1. The numbers stay 1, 2, 3.
    left = read_bordered_lc(_blc_slide_text(
        "slide-L", [sa], [sb * 2], "a", [-sa], [-sb * 2]))
    right = read_bordered_lc(_blc_slide_text(
        "slide-R", [sa], [sa * 3, sb * 2], "a_new", [-sa, -sa * 3], [-sb * 2]))
    validate_bordered_lc(left)
    validate_bordered_lc(right)
    _blc_make_pair("slide", left, right, Int[], Int[], false, ["q", "r"])
end

function _blc_slide_stop(err)
    kind = err.kind
    reason = kind == :path_self_intersection ? :band :
        kind == :path_meets_family ? :winding :
        kind == :disconnected_path ? :rectangle :
        kind == :path_order_required ? :rectangle :
        :band
    _blc_fail("7.7.2", reason, err.message)
end

function _blc_slide_route(lc::LickorishCode, a, b, band, avoid)
    band isa AbstractDict || _blc_fail("7.7.2", :band, "band needs start and finish")
    haskey(band, "start") && haskey(band, "finish") ||
        _blc_fail("7.7.2", :band, "band needs start and finish")
    route = try
        _move_certify_route(lc, band; avoid=avoid)
    catch err
        err isa LCMoveError || rethrow()
        _blc_slide_stop(err)
    end
    (route.start.component, route.finish.component) == (String(a), String(b)) ||
        _blc_fail("7.7.2", :band, "band ends on another curve")
    route
end

"""Build the pants pair from two curves and a band. Stops are separate reasons."""
function slide_bordered_family(lc::LickorishCode, a::AbstractString, b::AbstractString, band;
        family_a::AbstractString, family_b::AbstractString, avoid::Vector{String}=String[])
    validate_lc_input(lc)
    (haskey(lc.words, a) && haskey(lc.words, b)) ||
        _blc_fail("7.7.2", :curve, "missing curve")
    a == b && _blc_fail("7.7.2", :distinct, "a slide needs two curves")
    String(family_a) == String(family_b) ||
        _blc_fail("7.7.2", :family, "the two curves are not one family")
    for name in (a, b)
        word = lc.words[name]
        isempty(word) && _blc_fail("7.7.2", :empty, name)
        length(unique(abs.(word))) == length(word) || _blc_fail("7.7.2", :simple, name)
    end
    isempty(intersect(abs.(lc.words[a]), abs.(lc.words[b]))) ||
        _blc_fail("7.7.2", :intersect, "$a meets $b")
    avoided = isempty(avoid) ? String[a, b] : avoid
    route = _blc_slide_route(lc, a, b, band, avoided)
    si = -_move_dir(route.start.dart)
    sj = -_move_dir(route.finish.dart)
    slide_bordered_family(-si, -sj)
end

const _BLC_SLIDE_ROW = Dict{String,Any}(
    "id" => "bordered.slide",
    "inverse_id" => "bordered.slide",
    "digest" => "bordered",
    "relation" => "handle_slide",
    "layer" => "LC",
    "label_policy" => "pin_k",
    "genus_delta" => 0,
)

function _blc_slide_surgery(lc::LickorishCode, a, b, route)
    si = -_move_dir(route.start.dart)
    sj = -_move_dir(route.finish.dart)
    cycle = _move_barbell_boundary(si, sj)
    runs = Any[]
    for dart in cycle
        if abs(dart) == 1
            push!(runs, (kind=:curve, component=a, edge=route.start.edge, side=si,
                direction=sign(dart), survives=true))
        elseif abs(dart) == 2
            push!(runs, (kind=:curve, component=b, edge=route.finish.edge, side=sj,
                direction=sign(dart), survives=false))
        else
            push!(runs, (kind=:bank, route=route, side=si * sign(dart), direction=sign(dart)))
        end
    end
    witness = _move_compile_boundary(lc, _BLC_SLIDE_ROW, a, runs, Dict(
        "kind" => "pants", "route" => route, "boundary_cycle" => cycle,
        "coefficient" => si * sj, "stationary" => b))
    _apply_witness(lc, witness)
end

"""Replace the pants neighborhood of `band` by the other side of the pair it determines."""
function apply_slide_bordered(lc::LickorishCode, a::AbstractString, b::AbstractString, band;
        family_a::AbstractString, family_b::AbstractString, avoid::Vector{String}=String[])
    pair = slide_bordered_family(lc, a, b, band; family_a=family_a, family_b=family_b, avoid=avoid)
    avoided = isempty(avoid) ? String[a, b] : avoid
    route = _blc_slide_route(lc, a, b, band, avoided)
    result = _blc_slide_surgery(lc, a, b, route)
    (pair=pair, output=result.output, inverse=result.inverse)
end

"""Apply the swapped pants by the reverse fragment rewrite of a forward slide."""
function apply_slide_inverse(lc::LickorishCode, inverse)
    _apply_witness(lc, inverse).output
end

function _blc_flip_dart(dart)
    way = String(dart[3])
    Any[dart[1], dart[2], way == "out" ? "in" : "out"]
end

function _blc_band_end(raw; flip::Bool)
    dart = flip ? _blc_flip_dart(raw["dart"]) : Any[raw["dart"]...]
    Dict{String,Any}("dart" => dart, "u" => raw["u"])
end

"""
Inverse band. The parallel arc along the stationary curve is the opposite
half-edge at the finish, and the return crossings are that side in reverse
order. The replaced curve keeps the forward start half-edge, which is the
attachment of a_new.
"""
function reverse_slide_band(band::AbstractDict)
    raw = get(band, "crossings", Any[])
    crossings = Any[_blc_band_end(p; flip=true) for p in reverse(collect(raw))]
    Dict{String,Any}(
        "start" => _blc_band_end(band["start"]; flip=false),
        "finish" => _blc_band_end(band["finish"]; flip=true),
        "crossings" => crossings,
    )
end

# The forward pants word is the rotated source curve, then one fresh label per
# crossing in each bank direction, then one fresh label per label of b.
# Those fresh labels are the parallel copy. Deleting them and unrotating at
# the attachment edge is the geometric inverse. The flipped finish dart is not
# on one face with the start dart after the surgery, so certify rejects it.
function _blc_slide_mint(lc::LickorishCode, a, b, route)
    si = -_move_dir(route.start.dart)
    sj = -_move_dir(route.finish.dart)
    cycle = _move_barbell_boundary(si, sj)
    cycle[1] == 1 || _blc_fail("7.7.3", :band, "pants cycle does not start at +1")
    labels = Set(abs.(reduce(vcat, values(lc.words); init=Int[])))
    fresh = maximum(vcat([0], collect(labels)))
    copied = Int[]
    minted = Int[]
    for dart in cycle
        direction = sign(dart)
        if abs(dart) == 1
            word = lc.words[a]
            n = length(word)
            inds = direction == 1 ? [mod1(route.start.edge + k, n) for k in 1:n] :
                [mod1(route.start.edge - k, n) for k in 0:n - 1]
            for pos in inds
                push!(copied, direction * word[pos])
            end
        elseif abs(dart) == 2
            word = lc.words[b]
            n = length(word)
            inds = direction == 1 ? [mod1(route.finish.edge + k, n) for k in 1:n] :
                [mod1(route.finish.edge - k, n) for k in 0:n - 1]
            for _ in inds
                fresh += 1
                push!(minted, fresh)
            end
        else
            events = direction == 1 ? route.crossings : reverse(route.crossings)
            for _ in events
                fresh += 1
                push!(minted, fresh)
            end
        end
    end
    (copied=copied, minted=minted, edge=route.start.edge)
end

"""
Restore the diagram of a forward slide by deleting the parallel copy.

`source` is the diagram the band was drawn on. `output` is the diagram after
`apply_slide_bordered` of that same band. The attachment edge is the forward
start edge, so the copied labels unrotate back to the source word.
"""
function apply_reversed_slide(output::LickorishCode, source::LickorishCode,
        a::AbstractString, b::AbstractString, band; avoid::Vector{String}=String[])
    avoided = isempty(avoid) ? String[a, b] : avoid
    route = _blc_slide_route(source, a, b, band, avoided)
    info = _blc_slide_mint(source, a, b, route)
    n = length(info.copied)
    word = output.words[a]
    length(word) >= n || _blc_fail("7.7.3", :band, "replaced curve lost the attachment")
    prefix = word[1:n]
    prefix == info.copied || _blc_fail("7.7.3", :band, "attachment prefix is not the source curve")
    drop = Set(info.minted)
    restored = Int[prefix[mod1(j - info.edge, n)] for j in 1:n]
    words = Dict{String,Vector{Int}}()
    for (name, entry) in output.words
        words[name] = name == a ? restored : Int[x for x in entry if !(abs(x) in drop)]
    end
    _build_output(output, output.genus, words, copy(output.places))
end

# Host darts of the reversed band on the diagram after surgery.
# The start is the in half-edge of the first letter of the surviving complement.
# The finish is the forward finish with in and out exchanged.
function _blc_geometric_inverse_band(first_letter::Integer, band::AbstractDict)
    start = Dict{String,Any}(
        "dart" => Any[abs(first_letter), first_letter > 0 ? "+" : "-", "in"],
        "u" => "1/2",
    )
    raw = get(band, "crossings", Any[])
    crossings = Any[_blc_band_end(p; flip=true) for p in reverse(collect(raw))]
    Dict{String,Any}(
        "start" => start,
        "finish" => _blc_band_end(band["finish"]; flip=true),
        "crossings" => crossings,
    )
end

"""
Inverse of one slide, certified on the diagram after the surgery.

`output` is the result of `apply_slide_bordered` for `band` on `source`.
The rectangle is the route of the reversed band on `output`. Restored words
keep their signs and return to the source up to a cyclic shift. The record
holds the host darts of that band and the pushed-curve correspondence.
"""
function apply_geometric_slide_inverse(output::LickorishCode, source::LickorishCode,
        a::AbstractString, b::AbstractString, band; avoid::Vector{String}=String[])
    avoided = isempty(avoid) ? String[a, b] : avoid
    haskey(output.words, String(a)) || _blc_fail("7.7.3", :curve, "replaced curve is missing after surgery")
    # The forward route names the attachment. The rectangle of the inverse is
    # the following certification on `output`.
    forward = _blc_slide_route(source, a, b, band, avoided)
    info = _blc_slide_mint(source, a, b, forward)
    word = output.words[String(a)]
    n = length(info.copied)
    length(word) >= n + 1 || _blc_fail("7.7.3", :band, "replaced curve has no fresh label")
    prefix = word[1:n]
    prefix == info.copied || _blc_fail("7.7.3", :band, "attachment prefix is not the source curve")
    isempty(info.minted) && _blc_fail("7.7.3", :band, "the parallel copy has no fresh label")
    abs(word[n + 1]) == info.minted[1] ||
        _blc_fail("7.7.3", :band, "corridor does not start at the attachment index")
    inverse_band = _blc_geometric_inverse_band(word[1], band)
    certified = _blc_slide_route(output, a, b, inverse_band, avoided)
    drop = Set(info.minted)
    restored = Int[prefix[mod1(j - info.edge, n)] for j in 1:n]
    words = Dict{String,Vector{Int}}()
    for (name, entry) in output.words
        words[name] = name == String(a) ? restored : Int[x for x in entry if !(abs(x) in drop)]
    end
    diagram = _build_output(output, output.genus, words, copy(output.places))
    record = Dict{String,Any}(
        "inverse_start" => inverse_band["start"]["dart"],
        "inverse_finish" => inverse_band["finish"]["dart"],
        "inverse_start_edge" => certified.start.edge,
        "inverse_finish_edge" => certified.finish.edge,
        "attachment_index" => n + 1,
        "forward_edge" => info.edge,
        "replaced" => String(a),
        "stationary" => String(b),
        "minted" => info.minted,
        "copied" => info.copied,
    )
    (diagram=diagram, record=record)
end
