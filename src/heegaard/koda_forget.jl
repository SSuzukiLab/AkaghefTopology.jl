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

"""
    ForgottenHeegaard(genus, realization, alpha_words, beta_words, crossings, regions, diagram, correspondence) -> ForgottenHeegaard

Store the Heegaard arrangement and complement regions after forgetting auxiliary reference curves.

# Arguments
`genus::Int`, `realization::Union{Nothing,LCRealization}`,
`alpha_words::Dict{String,Vector{Int}}`, `beta_words::Dict{String,Vector{Int}}`,
`crossings::Vector{Dict{String,Any}}`, `regions::Vector{CombinatorialRegion}`,
`diagram::Union{Nothing,CellularHeegaardDiagram}`,
`correspondence::Dict{String,Any}`, or values convertible to those types.
Conversion errors propagate. This constructor does not perform forgetting.

# Returns
`ForgottenHeegaard` with the supplied surface genus, optional LC realization,
two curve-word dictionaries, crossing records, complement regions, optional
cellular diagram and correspondence data. Helper types are available in
`Topology.HeegaardCodes`. The object is specified in
[reference forgetting](../specs/heegaard-koda-conversion.md#41-忘却の計算).

# Guarantees
1. Fields preserve the supplied values after conversion; typed containers
   and supplied objects are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. For output of [`forget_reference`](@ref), the regions refer to the
   complement of the retained alpha and beta curves on the realized surface.
   The correspondence records which reference curves were forgotten or retained. [proved/agent: by construction]

# Not guaranteed
Direct construction does not verify the surface, incidence, region data or
joining-disk flags. `diagram` may be `nothing` even when forgetting returned
region data. Its presence is not checked against the other fields. No
3-manifold, combing or framing is identified by construction.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§4.1–4.2 — complement regions and transported marks after forgetting auxiliary K curves.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-13", "case-01.hc");

julia> f = H.forget_reference(H.realize_hc(H.read_hc(path)));

julia> g = H.ForgottenHeegaard(f.genus, f.realization, f.alpha_words, f.beta_words, f.crossings, f.regions, f.diagram, f.correspondence);

julia> only(r for r in g.regions if r.id == 3).joining_disk
false
```

See also [`forget_reference`](@ref), [`enumerate_joining_disk_candidates`](@ref).

Role: convert
"""
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

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/overlay/lib/topology.jl/docs/specs/hc-standard.md §計算と忘却 — K_forget is the reference curves with no α/β role
function _koda_reference_names_to_forget(genus::Integer, roles)
    [name for name in _lc_k_names(genus) if !haskey(roles, name)]
end

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/overlay/lib/topology.jl/docs/specs/hc-standard.md §計算と忘却 — reference names, then the other component names in dictionary order
function _koda_component_order(words, genus::Integer)
    reference = _lc_k_names(genus)
    vcat(reference, sort([name for name in keys(words) if !(name in reference)]))
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

function _crossing_records(words::AbstractDict, roles::Dict{String, Symbol}, occ; order=nothing)
    out = Dict{String, Any}[]
    for label in _ab_labels(words, roles)
        plus = occ[label]["+"]
        minus = occ[label]["-"]
        ap, an = plus[1], minus[1]
        ra, rb = _curve_family(ap, roles), _curve_family(an, roles)
        sign = 1
        ai, bi = 0, 0
        if ra == :alpha && rb == :beta
            # cf: docs/specs/heegaard-koda-conversion.md §7.3.1: ε is det(tα,tβ).
            sign = -1
            ai = _curve_index(words, roles, :alpha, ap; order=order)
            bi = _curve_index(words, roles, :beta, an; order=order)
        elseif ra == :beta && rb == :alpha
            sign = 1
            ai = _curve_index(words, roles, :alpha, an; order=order)
            bi = _curve_index(words, roles, :beta, ap; order=order)
        else
            throw(KodaImplementationError(:forget_reference, "non-transverse crossing", Dict("label" => label)))
        end
        push!(out, Dict("id" => label, "sign" => sign, "alpha" => ai, "beta" => bi))
    end
    sort!(out, by = x -> x["id"])
    out
end

function _curve_index(words, roles, fam::Symbol, name::AbstractString; order=nothing)
    ids = if order === nothing
        sort([n for (n, _) in words if _curve_family(n, roles) == fam])
    else
        [n for n in order if haskey(words, n) && _curve_family(n, roles) == fam]
    end
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

"""Normalized face-polygon quotient and maximal alpha/beta boundary arcs.
cf: lib/topology.jl/src/heegaard/lc_combinatorial.jl::_lc_faces; docs/specs/heegaard-koda-conversion.md §5.2 — glue only forgotten-reference sides [reconstructed]
"""
function _region_cell_data(words, roles, face_ids, realization, k_names; order=nothing)
    fail(msg, data=Dict{String,Any}()) =
        throw(KodaImplementationError(:forget_reference, msg, data))
    fset = Set(string(x) for x in face_ids)
    kset = Set(string(x) for x in k_names)
    faces = [f for f in realization.faces if f.id in fset]
    length(faces) == length(fset) ||
        fail("unknown LC face in region", Dict("faces" => collect(fset)))
    isempty(faces) && fail("empty LC face union")

    occ = realization.occurrence
    darts = LCCorner[]
    foreach(f -> append!(darts, f.corners), faces)
    length(unique(darts)) == length(darts) || fail("duplicate face dart")
    node = Dict{LCCorner,Int}(c => i for (i,c) in enumerate(darts))

    tau(c) = _lc_tau(words, occ, c)
    function nxt(c)
        t = tau(c)
        s,w = _lc_rot_next(t[2], t[3])
        (t[1], s, w)
    end
    dkey(c) = (c[1], c[2] == "-" ? 0 : 1, c[3] == "out" ? 0 : 1)

    comp = Dict{LCCorner,String}()
    for c in darts
        name, = _corner_component(words, c)
        name === nothing && fail("face dart has no component", Dict("corner" => c))
        comp[c] = string(name)
    end

    parent = collect(eachindex(darts))
    function root(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        i
    end
    function unite(i,j)
        i,j = root(i),root(j)
        i == j || (parent[j] = i)
    end

    # Quotient the disjoint selected face polygons only across forgotten-K sides.
    internal = Set{LCCorner}()
    paired = Set{LCCorner}()
    boundary = LCCorner[]
    E = 0
    for c in darts
        t = tau(c)
        if comp[c] in kset
            haskey(node,t) || fail("forgotten-K face union is not closed under tau",
                Dict("curve"=>comp[c], "corner"=>c, "mate"=>t))
            push!(internal,c)
            if !(c in paired)
                push!(paired,c); push!(paired,t); E += 1
                unite(node[c], node[nxt(t)])
                unite(node[nxt(c)], node[t])
            end
        else
            E += 1
            push!(boundary,c)
        end
    end
    V = length(Set(root(i) for i in eachindex(parent)))
    F = length(faces)
    chi = V - E + F

    # After a boundary side, rotate past every glued K side at that vertex.
    function bnext(c)
        d = nxt(c)
        for _ in 0:length(darts)
            !(d in internal) && return d
            d = nxt(tau(d))
        end
        fail("boundary walk is trapped at a K vertex", Dict("corner" => c))
    end

    cycles = Vector{Vector{LCCorner}}()
    if !isempty(boundary)
        bset = Set(boundary)
        succ = Dict(c => bnext(c) for c in boundary)
        all(d -> d in bset, values(succ)) || fail("bad boundary successor")
        length(Set(values(succ))) == length(boundary) ||
            fail("normalized boundary is not a 1-manifold")
        seen = Set{LCCorner}()
        for start in sort(boundary; by=dkey)
            start in seen && continue
            cyc = LCCorner[]
            local_seen = Set{LCCorner}()
            c = start
            while !(c in local_seen)
                c in seen && fail("boundary cycles merge", Dict("corner" => c))
                push!(local_seen,c); push!(cyc,c)
                c = succ[c]
            end
            c == start || fail("boundary closes at the wrong dart")
            union!(seen,local_seen)
            push!(cycles,cyc)
        end
    end

    b = length(cycles)
    twog = 2 - b - chi
    twog >= 0 && iseven(twog) ||
        fail("region is not a connected orientable surface",
             Dict("chi"=>chi, "boundary_components"=>b, "V"=>V, "E"=>E, "F"=>F))
    genus = twog ÷ 2

    frank(f) = f == :alpha ? 0 : f == :beta ? 1 : 2
    function mk_run(c)
        name = comp[c]
        fam = _curve_family(name,roles)
        fam in (:alpha,:beta) ||
            fail("non-alpha/beta boundary edge", Dict("curve"=>name, "corner"=>c))
        t = tau(c)
        (curve=name, family=fam, from=c[1], to=t[1],
         step=(c[3] == "out" ? 1 : -1), first=c, last=c)
    end
    function follows(a,z)
        t = tau(a.last)
        a.curve == z.curve && a.step == z.step &&
            z.first == (t[1],t[2],a.last[3])
    end
    merge_run(a,z) = (curve=a.curve, family=a.family, from=a.from, to=z.to,
                 step=a.step, first=a.first, last=z.last)

    # Compress only after topology is fixed; compression never enters V-E+F.
    grouped = Tuple{NTuple{9,Int},Vector{RegionPort}}[]
    for cyc in cycles
        runs = NamedTuple[]
        for c in cyc
            r = mk_run(c)
            if !isempty(runs) && follows(runs[end],r)
                runs[end] = merge_run(runs[end],r)
            else
                push!(runs,r)
            end
        end
        if length(runs) > 1 && follows(runs[end],runs[1])
            runs[1] = merge_run(runs[end],runs[1])
            pop!(runs)
        end

        entries = Tuple{RegionPort,LCCorner}[]
        # A whole retained circle stays as one closed port. It is rejected by
        # joining_disk_failure rather than lost from the boundary certificate.
        for r in runs
            cid = _curve_index(words,roles,r.family,r.curve; order=order)
            p = RegionPort(r.curve,r.family,cid,r.from,r.to,r.step)
            push!(entries,(p,r.first))
        end

        if isempty(entries)
            q = dkey(cyc[1])
            push!(grouped,((1,2,0,0,0,0,q[1],q[2],q[3]),RegionPort[]))
            continue
        end
        i0 = argmin(i -> begin
            p,c = entries[i]; q = dkey(c)
            (frank(p.family),p.curve_id,p.from_label,p.to_label,
             p.word_step,q[1],q[2],q[3])
        end, eachindex(entries))
        entries = circshift(entries,1-i0)
        ps = RegionPort[x[1] for x in entries]
        p,c = entries[1]; q = dkey(c)
        key = (0,frank(p.family),p.curve_id,p.from_label,p.to_label,
               p.word_step,q[1],q[2],q[3])
        # cf: heegaard-koda-conversion.md §§3.1, 7.3.1. The LC face walk
        # keeps the face on the right; positive boundary ports keep it on the left.
        boundary = [RegionPort(p.curve,p.family,p.curve_id,p.to_label,p.from_label,-p.word_step)
                    for p in Iterators.reverse(ps)]
        push!(grouped,(key,boundary))
    end

    sort!(grouped; by=x -> x[1])
    ports = RegionPort[]
    foreach(x -> append!(ports,x[2]), grouped)
    (chi,genus,b,ports)
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

"""
    joining_disk_failure(region::Topology.HeegaardCodes.CombinatorialRegion, words::AbstractDict, roles::Dict{String,Symbol}) -> Union{Nothing,Symbol}

Diagnose a forgotten region's implemented joining-disk checks.

# Arguments
`region`, `words` and `roles` must describe the same forgotten curve arrangement. The caller supplies consistent region topology and ordered boundary ports, as produced by `forget_reference`; this function does not reconstruct them. Word entries must convert to integer labels. Conversion and malformed-record errors propagate.

# Returns
`nothing` on success, otherwise a diagnostic symbol such as `:arc_not_on_word`, `:whole_curve`, `:no_ab_corner`, `:not_disk`, `:short_boundary`, `:unequal_family_count` or `:not_one_arc`.

# Guarantees
1. The result diagnoses the supplied data against the implemented [joining-disk domain](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk). In particular, whole-curve contact is rejected. [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Not guaranteed
Success alone does not independently certify an embedded disk closure, the consistency of hand-built region records, or complete Heegaard cut systems. Diagnostics do not enumerate all simultaneous failures.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.2 — joining disk, nonempty contact arcs and whole-curve exclusion.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> r = H.CombinatorialRegion(1, String[], 1, 0, 1, H.RegionPort[], false); H.joining_disk_failure(r, Dict{String,Vector{Int}}(), Dict{String,Symbol}())
:no_ab_corner
```

Role: condition
"""
function joining_disk_failure(region::CombinatorialRegion, words::AbstractDict, roles::Dict{String, Symbol})
    bd = region.boundary_ports
    # Diagnose whole-curve contact before topology, including annular regions.
    for p in bd
        haskey(words,p.curve) || return :arc_not_on_word
        word = Int[abs(x) for x in words[p.curve]]
        isempty(word) && continue
        steps = _arc_step_count(word,p.from_label,p.to_label,p.word_step)
        steps < 0 && return :arc_not_on_word
        steps >= length(word) && return :whole_curve
    end
    isempty(_ab_crossing_labels(words,roles)) && return :no_ab_corner
    region.chi == 1 && region.genus == 0 && region.boundary_components == 1 || return :not_disk
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

function _merge_faces_forget(realization::LCRealization, k_names::Vector{String}; include_single_visit::Bool=true)
    words = _lc_nonempty_words(realization.lc.words)
    isempty(words) && return realization.faces, Dict{String, String}()
    parent = Dict(f.id => f.id for f in realization.faces)
    corner_face = realization.corner_face
    kset = Set(k_names)
    for c in k_names
        haskey(words, c) || continue
        word = words[c]
        isempty(word) && continue
        !include_single_visit && length(word) == 1 && continue
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
function _cycle_region_blocks(forgot::ForgottenHeegaard, code::CombinatorialHeegaardCode)
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
    anames = get(forgot.correspondence, "alpha_order", sort(collect(keys(forgot.alpha_words))))
    bnames = get(forgot.correspondence, "beta_order", sort(collect(keys(forgot.beta_words))))
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

function _source_label_ids(h::HeegaardDiagram)
    ids = Dict{Int, Int}()
    for point in h.code.code.intersections
        point.source_label === nothing && continue
        label = tryparse(Int, point.source_label)
        label === nothing && continue
        ids[label] = point.id
    end
    ids
end

"""Consecutive darts of one curve are one boundary arc. A cycle that starts
mid-arc splits that arc across the two ends; those two runs are the same arc."""
function _ribbon_boundary_runs(h::HeegaardDiagram, cycle::Vector{Int})
    graph = ribbon_graph(h)
    runs = Tuple{Symbol, Int, Int, Int}[]
    for dart_id in cycle
        dart = graph.darts[dart_id]
        edge = graph.edges[dart.edge]
        if !isempty(runs) && runs[end][1] == edge.family && runs[end][2] == edge.curve_id
            prev = runs[end]
            runs[end] = (prev[1], prev[2], prev[3], dart.target)
        else
            push!(runs, (edge.family, edge.curve_id, dart.origin, dart.target))
        end
    end
    if length(runs) > 1 &&
       runs[1][1] == runs[end][1] && runs[1][2] == runs[end][2] &&
       runs[end][4] == runs[1][3]
        last = runs[end]
        first = runs[1]
        runs[1] = (first[1], first[2], last[3], first[4])
        pop!(runs)
    end
    runs
end

"""Same curve is not the same arc. The corners have to be the same α∩β points."""
function _ribbon_cycle_matches_ports(h::HeegaardDiagram, cycle::Vector{Int}, ports::Vector{RegionPort})
    runs = _ribbon_boundary_runs(h, cycle)
    length(runs) == length(ports) || return false
    labels = _source_label_ids(h)
    expected = Tuple{Symbol, Int, Int, Int}[]
    for port in ports
        from_id = get(labels, port.from_label, nothing)
        to_id = get(labels, port.to_label, nothing)
        (from_id === nothing || to_id === nothing) && return false
        push!(expected, (port.family, port.curve_id, from_id, to_id))
    end
    any(k -> circshift(runs, k) == expected, 0:(length(runs) - 1))
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

"""
    forget_reference(realized::KodaConversionResult) -> ForgottenHeegaard

Return the retained Heegaard arrangement and complement regions after forgetting auxiliary K curves.

# Arguments
`realized` normally comes from `realize_hc`; it must carry a compatible LC,
its realization and applied roles under
[the forgetting domain](../specs/heegaard-koda-conversion.md#41-忘却の計算).
The caller is responsible for manually assembled results. Missing keys,
wrong field types or inconsistent combinatorial data can raise lookup,
type or Koda conversion errors. Full input validity is not rechecked.

# Returns
`ForgottenHeegaard` with genus, the retained realization, copies of the
alpha/beta input words, crossing records, complement-region records, an
optional cellular `diagram`, and `correspondence`. The correspondence lists
forgotten/retained reference curves, region data, ribbon-region matches and
`diagram_build_error`. The family words retain their input labels.

# Guarantees
1. On a valid realized HC, the region data describe the complement of the
   retained alpha and beta curves under
   [reference forgetting](../specs/heegaard-koda-conversion.md#41-忘却の計算).
   [proved/agent: by construction]
2. Region joining flags record the current checks of
   [the joining domain](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).
   [proved/agent: by construction]
3. The input is not modified. The realization is shared, while the family
   word vectors are copied. [proved/agent: by construction]

# Not guaranteed
Successful forgetting need not produce a cellular diagram: a failure in that
construction is retained as `diagram_build_error` with `diagram === nothing`.
Region IDs can have gaps and are not array positions. This does not choose a
joining disk or construct a spine, combing or framing.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§4.1–4.2, 5.2 — reference forgetting, complement regions and joining checks.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc");

julia> f = forget_reference(realize_hc(read_hc(path)));

julia> length(enumerate_joining_disk_candidates(f))
2
```

See also [`enumerate_joining_disk_candidates`](@ref), [`as_punctured_heegaard`](@ref), [`heegaard_to_spine`](@ref).

Role: convert
"""
function forget_reference(realized::KodaConversionResult)
    doc = realized.value
    lc = doc["lc"]::LickorishCode
    realization = doc["realization"]::LCRealization
    roles = get(doc, "roles_applied", Dict{String, Symbol}())
    isempty(roles) && (roles = get(realized.correspondence, "roles_applied", Dict{String, Symbol}()))
    g = lc.genus
    k_names = _koda_reference_names_to_forget(g, roles)
    words = lc.words
    order = _koda_component_order(words, g)
    occ = realization.occurrence
    groups = _merge_faces_forget(realization, k_names)
    # Preserve established region labels under newly included one-visit gluings.
    # IDs are provenance labels, not array positions; a merged region inherits
    # the minimum pre-merge ID and labels can consequently have gaps.
    legacy = _merge_faces_forget(realization, k_names; include_single_visit=false)
    face_label = Dict(fid => i for (i, (_, fs)) in enumerate(sort(collect(legacy))) for fid in fs)
    regions = CombinatorialRegion[]
    words_active = _lc_nonempty_words(words)
    for (root, face_ids) in sort(collect(groups))
        rid = minimum(face_label[fid] for fid in face_ids)
        chi, h, b, ports = _region_cell_data(words_active, roles, Set(face_ids), realization, Set(k_names); order=order)
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
    crossings = _crossing_records(words_active, roles, occ; order=order)
    role_map = Dict{String, Any}(
        "alpha_order" => [n for n in order if get(roles, n, :none) == :alpha],
        "beta_order" => [n for n in order if get(roles, n, :none) == :beta],
        "forgotten_reference" => k_names,
        "retained_reference" => [n for n in _lc_k_names(g) if haskey(roles, n)],
    )
    forgot = ForgottenHeegaard(g, realization, αw, βw, crossings, regions, nothing, role_map)
    ribbon_map = Dict{Int, Int}()
    diagram_build_error = nothing
    try
        hdiag = _build_heegaard_from_combinatorial(forgot)
        ribbon_map = _map_ribbon_joining_regions!(regions, hdiag)
        forgot = ForgottenHeegaard(g, realization, αw, βw, crossings, regions, hdiag, forgot.correspondence)
    catch err
        diagram_build_error = sprint(showerror, err)
    end
    corr = merge(copy(role_map), Dict(
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
    ))
    ForgottenHeegaard(forgot.genus, forgot.realization, forgot.alpha_words, forgot.beta_words, forgot.crossings, forgot.regions, forgot.diagram, corr)
end


"""
    enumerate_joining_disk_candidates(forgot::ForgottenHeegaard) -> Vector

Return the complement regions recorded as joining-disk candidates after forgetting K.

# Arguments
`forgot` normally comes from `forget_reference`. Its region records must be
consistent with [joining conditions §5.2](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).
The caller is responsible for manually constructed records: stored flags
are used without recomputing or validating the conditions.

# Returns
A fresh vector of dictionaries with `region_id`, `joining_disk`, `chi`,
`genus` and `boundary_components`, in the order of the selected region records.

# Guarantees
1. Exactly the records with a true `joining_disk` flag are included.
   For output of `forget_reference`, these flags are its recorded checks of
   the joining domain. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The result does not select a disk, construct a punctured diagram or certify
that a cellular counterpart is available. IDs need not be contiguous indices
into the region vector.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.2 — joining union and joining-disk candidate conditions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc");

julia> f = forget_reference(realize_hc(read_hc(path)));

julia> length(enumerate_joining_disk_candidates(f))
2
```

See also [`as_punctured_heegaard`](@ref), [`heegaard_to_spine`](@ref).

Role: access
"""
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

"""
    as_punctured_heegaard(forgot::ForgottenHeegaard, selection; choices=Dict{String,Any}()) -> KodaConversionResult{PuncturedHeegaardDiagram}

Return the punctured cellular diagram associated with a selected joining region after forgetting K.

# Arguments
`forgot` is obtained from `forget_reference`. `selection` is a region ID or
its string representation and must select a joining disk satisfying
[Koda conditions §5.2](../specs/heegaard-koda-conversion.md#52-joining-union-と-joining-disk).
A missing selection, a non-joining region or unequal alpha/beta counts
throws `KodaInputError` describing the failed requirement. An unknown string
selection throws `KodaInputError`; an integer selection without a unique
matching region throws `ArgumentError`. A missing cellular diagram or ribbon
counterpart throws `KodaImplementationError`. Joining-disk validation errors
propagate. `choices` must be a mutable string-keyed dictionary accepting the
selected region ID; the default has type `Dict{String,Any}`.

# Returns
A `KodaConversionResult` with the punctured diagram in `value`, the updated
`choices`, a `correspondence` containing `joining_region`, and
`source_status == :constructed`. The diagram retains the cellular source
and selected joining disk; its full alpha-cut chart is deferred.

# Guarantees
1. On success, the selected region has a joining-disk counterpart in the
   returned cellular diagram under the
   [forgetful correspondence](../specs/heegaard-koda-conversion.md#42-印の輸送).
   [proved/agent: by construction]
2. `forgot` is not modified. The supplied `choices` dictionary is updated at
   `"joining_region"` and is retained in the result when already correctly
   typed. [proved/agent: by construction]

# Not guaranteed
This operation does not compute a full alpha-cut chart, a spine or o-data.
The embedded conversion certificate retains its unverified evidence status;
`:constructed` alone is not a literature proof of the entire Koda correspondence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§4.2, 5.2, 6 — joining-region transport and the punctured input to the spine construction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc");

julia> f = forget_reference(realize_hc(read_hc(path)));

julia> as_punctured_heegaard(f, 2).value.joining_disk.region
2
```

See also [`enumerate_joining_disk_candidates`](@ref), [`heegaard_to_spine`](@ref).

Role: convert
"""
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
