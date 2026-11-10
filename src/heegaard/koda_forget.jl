# Source: docs/specs/heegaard-koda-conversion.md §4.1, §5.2.
"""Combinatorial K removal, complement regions, §5.2 joining disk tests."""

struct RegionPort
    curve::String
    family::Symbol
    curve_id::Int
    from_label::Int
    to_label::Int
    # +1 if the ∂J walk leaves this corner along the LC word, -1 if against it.
    # LCCorner way "out" steps forward in the word; "in" steps backward
    # (lc_combinatorial.jl `_lc_along`). Labels do not determine this when the
    # word has at most two visits.
    word_step::Int
end

struct CombinatorialRegion
    id::Int
    face_roots::Vector{String}
    chi::Int
    genus::Int
    boundary_components::Int
    boundary_ports::Vector{RegionPort}
    joining_disk::Bool
end

struct ForgottenHeegaard
    genus::Int
    realization::Union{Nothing, LCRealization}
    alpha_words::Dict{String, Vector{Int}}
    beta_words::Dict{String, Vector{Int}}
    crossings::Vector{Dict{String, Any}}
    regions::Vector{CombinatorialRegion}
    diagram::Union{Nothing, HeegaardDiagram}
    correspondence::Dict{String, Any}
end

function _uf_find!(parent::Dict{String, String}, a::String)
    while parent[a] != a
        parent[a] = parent[parent[a]]
        a = parent[a]
    end
    a
end

function _uf_union!(parent::Dict{String, String}, a::String, b::String)
    ra, rb = _uf_find!(parent, a), _uf_find!(parent, b)
    parent[rb] = ra
end

function _curve_family(name::AbstractString, roles::Dict{String, Symbol})
    get(roles, name, :none)
end

function _ab_labels(words::AbstractDict, roles::Dict{String, Symbol})
    labels = Set{Int}()
    for (name, word) in words
        fam = _curve_family(name, roles)
        fam in (:alpha, :beta) || continue
        for x in word
            push!(labels, abs(x))
        end
    end
    alpha = Set{Int}()
    beta = Set{Int}()
    for (name, word) in words
        fam = _curve_family(name, roles)
        fam == :alpha && (for x in word; push!(alpha, abs(x)); end)
        fam == :beta && (for x in word; push!(beta, abs(x)); end)
    end
    sort(collect(intersect(alpha, beta)))
end

function _crossing_records(words::AbstractDict, roles::Dict{String, Symbol}, occ)
    out = Dict{String, Any}[]
    for label in _ab_labels(words, roles)
        plus = occ[label]["+"]
        minus = occ[label]["-"]
        ap, an = plus[1], minus[1]
        ra, rb = _curve_family(ap, roles), _curve_family(an, roles)
        sign = 1
        ai, bi = 0, 0
        if ra == :alpha && rb == :beta
            sign = 1
            ai = _curve_index(words, roles, :alpha, ap)
            bi = _curve_index(words, roles, :beta, an)
        elseif ra == :beta && rb == :alpha
            sign = -1
            ai = _curve_index(words, roles, :alpha, an)
            bi = _curve_index(words, roles, :beta, ap)
        else
            throw(KodaImplementationError(:forget_reference, "non-transverse crossing", Dict("label" => label)))
        end
        push!(out, Dict("id" => label, "sign" => sign, "alpha" => ai, "beta" => bi))
    end
    sort!(out, by = x -> x["id"])
    out
end

function _curve_index(words, roles, fam::Symbol, name::AbstractString)
    ids = sort([n for (n, _) in words if _curve_family(n, roles) == fam])
    idx = findfirst(==(name), ids)
    idx === nothing && throw(KodaImplementationError(:forget_reference, "unknown curve", Dict("name" => name)))
    idx
end

function _corner_component(words, c::LCCorner)
    label, sg, = c
    for (name, word) in words
        for (pos, x) in enumerate(word)
            abs(x) == label && _lc_sign_char(x) == sg && return (name, pos)
        end
    end
    return (nothing, 0)
end

function _lc_face_across(words, occ, corner_face::Dict{LCCorner, String}, c::LCCorner)
    other, s2, w2 = _lc_tau(words, occ, c)
    s3, w3 = _lc_rot_prev(s2, w2)
    get(corner_face, (other, s3, w3), nothing)
end

"""Ordered boundary ports on ∂J: one entry per maximal α/β arc (§5.2)."""
function _region_boundary_ports(
    words::AbstractDict,
    roles::Dict{String, Symbol},
    face_ids::Set{String},
    faces::Vector{LCFace},
    corner_face::Dict{LCCorner, String},
    occ,
    k_names::Set{String},
)
    face_set = face_ids
    edges = Tuple{LCCorner, LCCorner, String}[]
    for f in faces
        f.id in face_set || continue
        n = length(f.corners)
        for i in 1:n
            c1 = f.corners[i]
            c2 = f.corners[mod1(i + 1, n)]
            across = _lc_face_across(words, occ, corner_face, c1)
            across === nothing && continue
            across in face_set && continue
            comp, = _corner_component(words, c1)
            comp === nothing && continue
            comp in k_names && continue
            fam = _curve_family(comp, roles)
            fam in (:alpha, :beta) || continue
            push!(edges, (c1, c2, comp))
        end
    end
    isempty(edges) && return RegionPort[]
    forward = Dict{LCCorner, Tuple{LCCorner, String}}()
    for (c1, c2, comp) in edges
        !haskey(forward, c1) && (forward[c1] = (c2, comp))
    end
    ordered = Tuple{String, Int, Int, Int}[]
    if !isempty(forward)
        start = sort(collect(keys(forward)))[1]
        c = start
        for _ in 1:length(forward)
            haskey(forward, c) || break
            c2, comp = forward[c]
            step = c[3] == "out" ? 1 : -1
            push!(ordered, (comp, c[1], c2[1], step))
            c = c2
            c == start && break
        end
    end
    if length(ordered) != length(forward)
        ordered = Tuple{String, Int, Int, Int}[
            (comp, c1[1], c2[1], c1[3] == "out" ? 1 : -1) for (c1, c2, comp) in edges
        ]
    end
    darts = ordered
    # merge consecutive darts on the same curve into arcs
    ports = RegionPort[]
    for (comp, a, b, step) in darts
        fam = _curve_family(comp, roles)
        if !isempty(ports) && ports[end].curve == comp
            # Later corners on the same boundary arc are arrival corners. Their
            # way field is the encoding of the arrival, not a second measurement
            # of the arc. The arc direction is the first corner's word step.
            ports[end] = RegionPort(ports[end].curve, ports[end].family, ports[end].curve_id, ports[end].from_label, b, ports[end].word_step)
        else
            cid = _curve_index(words, roles, fam, comp)
            push!(ports, RegionPort(string(comp), fam, cid, a, b, step))
        end
    end
    # rotate so boundary cycle starts deterministically
    if length(ports) > 1
        start = argmin(i -> (ports[i].family == :alpha ? 0 : 1, ports[i].curve), eachindex(ports))
        ports = circshift(ports, 1 - start)
    end
    ports
end

function _region_topology(face_count::Int, ports::Vector{RegionPort})
    isempty(ports) && return (0, 0, 0)
    verts = Set{Int}()
    for p in ports
        push!(verts, p.from_label, p.to_label)
    end
    v = length(verts)
    e = length(ports)
    chi = v - e + face_count
    b = 1
    twog = 2 - chi - b
    g = twog >= 0 && iseven(twog) ? twog ÷ 2 : max(0, twog ÷ 2)
    chi, g, b
end

function _contiguous_on_cycle(indices::Vector{Int}, n::Int)
    isempty(indices) && return false
    length(unique(indices)) == length(indices) || return false
    length(indices) == 1 && return true
    s = sort(indices)
    gaps = [s[i + 1] - s[i] for i in 1:(length(s) - 1)]
    if all(==(1), gaps)
        return true
    end
    # wrap-around contiguous
    return s[1] + (n - s[end]) == 1
end

function _ab_crossing_labels(words::AbstractDict, roles::Dict{String, Symbol})
    α = Set{Int}()
    β = Set{Int}()
    for (name, word) in words
        fam = _curve_family(name, roles)
        fam == :alpha && union!(α, abs.(word))
        fam == :beta && union!(β, abs.(word))
    end
    intersect(α, β)
end

"""Steps along `word` from `from_l` to `to_l` in `word_step` direction. `n` means the whole curve."""
function _arc_step_count(word::Vector{Int}, from_l::Int, to_l::Int, step::Int)
    n = length(word)
    n == 0 && return 0
    step in (1, -1) || return -1
    start = findfirst(x -> abs(x) == from_l, word)
    start === nothing && return -1
    pos = start
    for k in 1:n
        pos = mod1(pos + step, n)
        abs(word[pos]) == to_l && return k
    end
    -1
end

"""§5.2 failure, or `nothing` when the region is a joining disk.

A contact that is the whole closed curve is not an arc. Endpoints must be α∩β corners.
cf: docs/specs/heegaard-koda-conversion.md §5.2
"""
function joining_disk_failure(region::CombinatorialRegion, words::AbstractDict, roles::Dict{String, Symbol})
    region.chi == 1 && region.genus == 0 && region.boundary_components == 1 || return :not_disk
    bd = region.boundary_ports
    length(bd) >= 2 || return :short_boundary
    αnames = sort([n for (n, _) in words if _curve_family(n, roles) == :alpha])
    βnames = sort([n for (n, _) in words if _curve_family(n, roles) == :beta])
    length(αnames) == length(βnames) || return :unequal_family_count
    ab = _ab_crossing_labels(words, roles)
    for name in vcat(αnames, βnames)
        hits = [i for (i, p) in enumerate(bd) if p.curve == name]
        length(hits) == 1 || return :not_one_arc
        _contiguous_on_cycle(hits, length(bd)) || return :not_one_arc
        word = Int[abs(x) for x in words[name]]
        ports = [bd[i] for i in hits]
        steps = 0
        for p in ports
            k = _arc_step_count(word, p.from_label, p.to_label, p.word_step)
            k < 0 && return :arc_not_on_word
            steps += k
        end
        steps >= length(word) && return :whole_curve
        (ports[1].from_label in ab && ports[end].to_label in ab) || return :no_ab_corner
    end
    nothing
end

function _spec_joining_disk(
    region::CombinatorialRegion,
    words::AbstractDict,
    roles::Dict{String, Symbol},
)
    joining_disk_failure(region, words, roles) === nothing
end

function _merge_faces_forget(realization::LCRealization, k_names::Vector{String})
    words = _lc_nonempty_words(realization.lc.words)
    isempty(words) && return realization.faces, Dict{String, String}()
    parent = Dict(f.id => f.id for f in realization.faces)
    corner_face = realization.corner_face
    kset = Set(k_names)
    for c in k_names
        haskey(words, c) || continue
        word = words[c]
        length(word) < 2 && continue
        n = length(word)
        for i in 1:n
            c_out = (abs(word[i]), _lc_sign_char(word[i]), "out")
            c_in = (abs(word[mod1(i + 1, n)]), _lc_sign_char(word[mod1(i + 1, n)]), "in")
            f1 = get(corner_face, c_out, nothing)
            f2 = get(corner_face, c_in, nothing)
            (f1 === nothing || f2 === nothing) && continue
            _uf_union!(parent, f1, f2)
        end
    end
    groups = Dict{String, Vector{String}}()
    for f in realization.faces
        r = _uf_find!(parent, f.id)
        push!(get!(groups, r, String[]), f.id)
    end
    groups
end

# Computed region blocks for a non-cellular complement. Not used: certify_heegaard
# rejects them with closed_heegaard/not_closed (validation.jl). See the round-6 report.
function _cycle_region_blocks(forgot::ForgottenHeegaard, code::HeegaardCode)
    graph = ribbon_graph(code)
    cycles = boundary_cycles(graph)
    nreg = length(forgot.regions)
    length(cycles) == nreg && return [[i] for i in eachindex(cycles)]
    forgot.realization === nothing && return nothing
    xs = forgot.crossings
    pids = Dict(x["id"] => i for (i, x) in enumerate(xs))
    anames = sort(collect(keys(forgot.alpha_words)))
    bnames = sort(collect(keys(forgot.beta_words)))
    fam_id = Dict{String, Tuple{Symbol, Int}}()
    word_of = Dict{String, Vector{Int}}()
    for (i, name) in enumerate(anames)
        fam_id[name] = (:alpha, i)
        word_of[name] = forgot.alpha_words[name]
    end
    for (i, name) in enumerate(bnames)
        fam_id[name] = (:beta, i)
        word_of[name] = forgot.beta_words[name]
    end
    by_curve = Dict{Tuple{Symbol, Int}, Vector{Int}}()
    for e in graph.edges
        push!(get!(by_curve, (e.family, e.curve_id), Int[]), e.id)
    end
    cover = Dict{Tuple{Symbol, Int, Int}, Int}()
    for (name, word) in word_of
        fam, cid = fam_id[name]
        n = length(word)
        visits = Int[i for i in 1:n if haskey(pids, abs(word[i]))]
        eids = get(by_curve, (fam, cid), Int[])
        length(eids) == length(visits) || return nothing
        for k in eachindex(visits)
            s = visits[k]
            stop = visits[mod1(k + 1, length(visits))]
            guard = 0
            while s != stop && guard <= n
                cover[(fam, cid, s)] = eids[k]
                s = mod1(s + 1, n)
                guard += 1
            end
        end
    end
    real = forgot.realization
    face_region = Dict{String, Int}()
    for r in forgot.regions, fid in r.face_roots
        face_region[fid] = r.id
    end
    region_edges = Dict(r.id => Set{Int}() for r in forgot.regions)
    for (corner, fid) in real.corner_face
        rid = get(face_region, fid, 0)
        rid == 0 && continue
        comp, pos = _corner_component(real.lc.words, corner)
        comp === nothing && continue
        haskey(fam_id, string(comp)) || continue
        fam, cid = fam_id[string(comp)]
        step = corner[3] == "out" ? pos : mod1(pos - 1, length(word_of[string(comp)]))
        eid = get(cover, (fam, cid, step), 0)
        eid == 0 && continue
        push!(region_edges[rid], eid)
    end
    owned = Dict(r.id => Int[] for r in forgot.regions)
    for (ci, cyc) in enumerate(cycles)
        eids = Set(graph.darts[d].edge for d in cyc)
        hits = [rid for (rid, es) in region_edges if eids ⊆ es]
        length(hits) == 1 || return nothing
        push!(owned[hits[1]], ci)
    end
    blocks = Vector{Int}[]
    for r in forgot.regions
        isempty(owned[r.id]) && return nothing
        push!(blocks, sort(owned[r.id]))
    end
    blocks
end

function _build_heegaard_from_combinatorial(forgot::ForgottenHeegaard)
    xs = forgot.crossings
    labels = [x["id"] for x in xs]
    pids = Dict(x => i for (i, x) in enumerate(labels))
    anames = sort(collect(keys(forgot.alpha_words)))
    bnames = sort(collect(keys(forgot.beta_words)))
    α_curves = [
        Curve(:alpha, i, Int[pids[abs(x)] for x in forgot.alpha_words[name] if haskey(pids, abs(x))]; label = name)
        for (i, name) in enumerate(anames)
    ]
    β_curves = [
        Curve(:beta, i, Int[pids[abs(x)] for x in forgot.beta_words[name] if haskey(pids, abs(x))]; label = name)
        for (i, name) in enumerate(bnames)
    ]
    inter = [
        Intersection(pids[x["id"]], x["alpha"], x["beta"], x["sign"]; source_label = string(x["id"]))
        for x in xs
    ]
    certify_heegaard(decorate(heegaard_code(α_curves, β_curves, inter); cellular = true))
end

function _ribbon_joining_indices(h::HeegaardDiagram)
    idxs = Int[]
    for i in eachindex(h.code.regions.boundary_cycle_blocks)
        try
            joining_disk(h, i)
            push!(idxs, i)
        catch
        end
    end
    idxs
end

function _ribbon_boundary_arcs(h::HeegaardDiagram, cycle::Vector{Int})
    graph = ribbon_graph(h)
    [(graph.edges[graph.darts[d].edge].family, graph.edges[graph.darts[d].edge].curve_id) for d in cycle]
end

"""Match §5.2 combinatorial boundary ports to a ribbon complement region."""
function _ribbon_cycle_matches_ports(h::HeegaardDiagram, cycle::Vector{Int}, ports::Vector{RegionPort})
    arcs = _ribbon_boundary_arcs(h, cycle)
    length(arcs) == length(ports) || return false
    expected = [(p.family, p.curve_id) for p in ports]
    for k in 0:(length(arcs) - 1)
        rotated = circshift(arcs, -k)
        rotated == expected && return true
    end
    port_keys = [(p.family, p.curve_id) for p in ports]
    length(arcs) == length(port_keys) || return false
    for key in unique(port_keys)
        hits = [i for (i, a) in enumerate(arcs) if a == key]
        length(hits) == 1 || return false
        !_contiguous_on_cycle(hits, length(arcs)) && return false
    end
    sort(arcs) == sort(port_keys)
end

function _find_ribbon_joining_disk(h::HeegaardDiagram, region::CombinatorialRegion)
    rb = boundary_cycles(h)
    blocks = h.code.regions.boundary_cycle_blocks
    for (ri, block) in enumerate(blocks)
        length(block) != 1 && continue
        cycle = rb[only(block)]
        _ribbon_cycle_matches_ports(h, cycle, region.boundary_ports) || continue
        try
            return joining_disk(h, ri)
        catch
        end
    end
    nothing
end

function _map_ribbon_joining_regions!(regions::Vector{CombinatorialRegion}, h::HeegaardDiagram)
    out = Dict{Int, Int}()
    for r in regions
        r.joining_disk || continue
        jd = _find_ribbon_joining_disk(h, r)
        jd === nothing && continue
        out[r.id] = jd.region
    end
    out
end

function forget_reference(realized::KodaConversionResult)
    doc = realized.value
    lc = doc["lc"]::LickorishConfiguration
    realization = doc["realization"]::LCRealization
    roles = get(doc, "roles_applied", Dict{String, Symbol}())
    isempty(roles) && (roles = get(realized.correspondence, "roles_applied", Dict{String, Symbol}()))
    g = lc.genus
    k_names = _lc_k_names(g)
    words = lc.words
    occ = realization.occurrence
    groups = _merge_faces_forget(realization, k_names)
    regions = CombinatorialRegion[]
    rid = 0
    words_active = _lc_nonempty_words(words)
    for (root, face_ids) in sort(collect(groups))
        rid += 1
        ports = _region_boundary_ports(
            words_active,
            roles,
            Set(face_ids),
            realization.faces,
            realization.corner_face,
            realization.occurrence,
            Set(k_names),
        )
        chi, h, b = _region_topology(length(face_ids), ports)
        region = CombinatorialRegion(rid, face_ids, chi, h, b, ports, false)
        region = CombinatorialRegion(
            rid,
            face_ids,
            chi,
            h,
            b,
            ports,
            _spec_joining_disk(region, words_active, roles),
        )
        push!(regions, region)
    end
    αw = Dict(n => copy(words[n]) for n in keys(words) if _curve_family(n, roles) == :alpha)
    βw = Dict(n => copy(words[n]) for n in keys(words) if _curve_family(n, roles) == :beta)
    crossings = _crossing_records(words_active, roles, occ)
    forgot = ForgottenHeegaard(g, realization, αw, βw, crossings, regions, nothing, Dict{String, Any}())
    ribbon_map = Dict{Int, Int}()
    diagram_build_error = nothing
    try
        hdiag = _build_heegaard_from_combinatorial(forgot)
        ribbon_map = _map_ribbon_joining_regions!(regions, hdiag)
        forgot = ForgottenHeegaard(g, realization, αw, βw, crossings, regions, hdiag, forgot.correspondence)
    catch err
        diagram_build_error = sprint(showerror, err)
    end
    corr = Dict(
        "regions" => [
            Dict(
                "id" => r.id,
                "chi" => r.chi,
                "genus" => r.genus,
                "boundary_components" => r.boundary_components,
                "joining_disk" => r.joining_disk,
                "faces" => r.face_roots,
            ) for r in regions
        ],
        "crossings" => length(crossings),
        "ribbon_joining_map" => ribbon_map,
        "diagram_build_error" => diagram_build_error,
    )
    ForgottenHeegaard(forgot.genus, forgot.realization, forgot.alpha_words, forgot.beta_words, forgot.crossings, forgot.regions, forgot.diagram, corr)
end


function enumerate_joining_disk_candidates(forgot::ForgottenHeegaard)
    [Dict("region_id" => r.id, "joining_disk" => r.joining_disk, "chi" => r.chi, "genus" => r.genus, "boundary_components" => r.boundary_components) for r in forgot.regions if r.joining_disk]
end

function _resolve_joining_region(forgot::ForgottenHeegaard, selection)
    selection === nothing && return nothing
    if selection isa Integer
        return only(r for r in forgot.regions if r.id == selection)
    end
    key = string(selection)
    for r in forgot.regions
        string(r.id) == key && return r
        r.joining_disk && string(r.id) == key && return r
    end
    throw(KodaInputError(:as_punctured_heegaard, "unknown joining disk selection", Dict("selection" => selection)))
end

"""Punctured diagram for Koda spine when no full alpha-cut chart is recorded (§5–§6 words-only path)."""
function make_punctured_combinatorial(H::HeegaardDiagram, jd::JoiningDisk)
    validate_joining_disk(H, jd)
    # Labelled content hash only. canonical_hash searches curve relabelings and
    # is not part of the spine; on a genus-2 diagram with long words it dominates.
    identity = _labelled_hash(H)
    cut = CutSurface(
        string(identity),
        0,
        CutBoundary[],
        Dict{Int, Vector{Int}}(),
        Dict{Int, Vector{Int}}(),
        Dict{Int, Int}(),
        Dict{Int, Int}(),
    )
    cert = ConversionCertificate(
        :koda,
        identity,
        _disk_json(jd),
        Any[],
        KODA_CONVENTION,
        Dict{String, Any}("alpha_cut_deferred" => true),
        nothing,
        Dict{String, Any}(),
        :unverified,
    )
    PuncturedHeegaardDiagram(H, deepcopy(jd), cut, cert)
end

function as_punctured_heegaard(forgot::ForgottenHeegaard, selection; choices = Dict{String, Any}())
    region = _resolve_joining_region(forgot, selection)
    region === nothing && throw(KodaInputError(:as_punctured_heegaard, "joining disk selection required", Dict()))
    region.joining_disk ||
        throw(KodaInputError(:as_punctured_heegaard, "§5.2 joining disk conditions fail", Dict("region" => region.id)))
    forgot.diagram === nothing &&
        throw(KodaImplementationError(:as_punctured_heegaard, "Heegaard diagram not built", Dict("region" => region.id)))
    core = _undecorated(forgot.diagram)
    length(core.alpha) == length(core.beta) ||
        throw(KodaInputError(:as_punctured_heegaard, "#alpha must equal #beta", Dict()))
    jd = _find_ribbon_joining_disk(forgot.diagram, region)
    jd === nothing && throw(KodaImplementationError(
        :as_punctured_heegaard,
        "no ribbon region for combinatorial joining disk",
        Dict("region" => region.id, "ports" => length(region.boundary_ports)),
    ))
    ph = make_punctured_combinatorial(forgot.diagram, jd)
    merge!(choices, Dict("joining_region" => region.id))
    KodaConversionResult(ph, choices, Dict{String, Any}("joining_region" => region.id), :constructed)
end
