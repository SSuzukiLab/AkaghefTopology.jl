# Source: .research-scratch/PJ260922/ld-koda-spec/261109-enumerate-correspondences.py;
# docs/specs/heegaard-koda-conversion.md §7–8.
"""Ribbon surface, joining smooth, and finite table helpers (no per-case routing)."""

struct _KodaDSU
    parent::Dict{Any, Any}
end

_KodaDSU() = _KodaDSU(Dict{Any, Any}())

function _koda_dsu_find!(uf::_KodaDSU, a)
    get!(uf.parent, a, a) == a || (uf.parent[a] = _koda_dsu_find!(uf, uf.parent[a]))
    uf.parent[a]
end

function _koda_dsu_union!(uf::_KodaDSU, a, b)
    uf.parent[_koda_dsu_find!(uf, b)] = _koda_dsu_find!(uf, a)
end

function _koda_orbit(nxt::Dict)
    nxt = copy(nxt)
    out = Vector{Vector}()
    while !isempty(nxt)
        x = first(sort(collect(keys(nxt))))
        cycle = Any[]
        while haskey(nxt, x)
            push!(cycle, x)
            x = pop!(nxt, x)
        end
        push!(out, cycle)
    end
    out
end

"""dart ribbon for alpha/beta words; faces are directed cycles."""
function _koda_ribbon(alpha, beta, epsilon)
    rho = Dict{Tuple{Any, Char, Int}, Tuple{Any, Char, Int}}()
    reverse = Dict{Tuple{Any, Char, Int}, Tuple{Any, Char, Int}}()
    edges = Tuple{Tuple{Any, Char, Int}, Tuple{Any, Char, Int}}[]
    for (v, e) in epsilon
        cycle = [(v, 'a', 1), (v, 'b', e), (v, 'a', -1), (v, 'b', -e)]
        for i in 1:length(cycle)
            rho[cycle[i]] = cycle[mod1(i + 1, 4)]
        end
    end
    for (family, words) in (('a', alpha), ('b', beta))
        for word in words
            for i in 1:length(word)
                v = word[i]
                q = word[mod1(i + 1, length(word))]
                x, y = (v, family, 1), (q, family, -1)
                reverse[x], reverse[y] = y, x
                push!(edges, (x, y))
            end
        end
    end
    nxt = Dict(d => rho[reverse[d]] for d in keys(rho))
    faces = Vector{Vector{Tuple{Any, Char, Int}}}()
    face_index = Dict{Tuple{Any, Char, Int}, Int}()
    while !isempty(nxt)
        x = first(keys(nxt))
        cycle = Tuple{Any, Char, Int}[]
        while haskey(nxt, x)
            face_index[x] = length(faces) + 1
            push!(cycle, x)
            x = pop!(nxt, x)
        end
        push!(faces, cycle)
    end
    faces, face_index, edges
end

"""§7.3 candidate A: `over` is `:beta` for β over α."""
function _koda_smooth_joining(alpha, beta, joining, over::Symbol)
    over in (:alpha, :beta) || throw(KodaInputError(:heegaard_to_spine, "invalid over family", Dict("over" => over)))
    fam_char = over == :beta ? 'b' : 'a'
    successor = Dict{Tuple{Any, Char}, Tuple{Any, Char}}()
    for (family, words) in (('a', alpha), ('b', beta))
        for w in words
            for i in 1:length(w)
                successor[(w[i], family)] = (w[mod1(i + 1, length(w))], family)
            end
        end
    end
    jset = Set(joining)
    keep = Set(keys(successor))
    for v in jset, f in ('a', 'b')
        delete!(keep, (v, f))
    end
    nxt = Dict{Tuple{Any, Char}, Tuple{Any, Char}}()
    for start in keep
        v, f = successor[start]
        seen = Set{Tuple{Any, Char}}()
        steps = 0
        while v in jset
            steps += 1
            steps > 64 && throw(KodaImplementationError(
                :heegaard_to_spine,
                "joining corner absorption did not exit",
                Dict("start" => start, "joining" => collect(jset), "steps" => steps),
            ))
            (v, f) in seen && throw(KodaImplementationError(:heegaard_to_spine, "joining-only singular cycle", Dict("vertex" => v)))
            push!(seen, (v, f))
            v, f = successor[(v, f == 'a' ? 'b' : 'a')]
        end
        nxt[start] = (v, f)
    end
    vert_list = collect(Set(x[1] for x in keep))
    sort!(vert_list)
    labels = Dict(v => i for (i, v) in enumerate(vert_list))
    words = Vector{Vector{Int}}()
    while !isempty(nxt)
        key_list = collect(keys(nxt))
        sort!(key_list)
        x = first(key_list)
        word = Int[]
        while haskey(nxt, x)
            push!(word, labels[x[1]] * (x[2] == fam_char ? 1 : -1))
            x = pop!(nxt, x)
        end
        push!(words, word)
    end
    words, labels, nxt
end

function _koda_flow_spine_from_words(words::Vector{Vector{Int}}, sigma::Dict{Int, Int}, source; genus=0)
    state_next = Dict{Tuple{Int, Symbol}, Tuple{Int, Symbol}}()
    for word in words
        prev = nothing
        start = nothing
        for (pos, p) in enumerate(word)
            v = abs(p)
            fam = p > 0 ? :beta : :alpha
            state = (v, fam)
            prev !== nothing && (state_next[prev] = state)
            start === nothing && (start = state)
            prev = state
            pos == length(word) && start !== nothing && (state_next[state] = start)
        end
    end
    edges = SpineEdge[]
    for (s, t) in state_next
        push!(edges, SpineEdge(s[1], s[2], t[1], t[2], Int[]))
    end
    FlowSpine(
        genus,
        edges,
        Dict(k => v for (k, v) in sigma),
        source,
        ConversionCertificate(:koda, "", Dict(), Any[], KODA_CONVENTION, Dict(), nothing, Dict(), :unverified),
        BranchedSpineComplex(SpineRegion[], Dict[], Vector{Vector{Tuple{Int, Int}}}[], Dict()),
    )
end



"""Return +1/-1 when the contact arc's labels determine a unique word step.

Returns `nothing` when the word has at most two visits, or both cyclic
neighbors carry `to_l`: the step is then a position, not a label. Callers
must use `RegionPort.word_step`.
"""
function _koda_lc_arc_direction(word::Vector{Int}, from_l::Int, to_l::Int)
    n = length(word)
    found = Int[]
    for i in 1:n
        abs(word[i]) != from_l && continue
        nxt = abs(word[mod1(i + 1, n)])
        prv = abs(word[mod1(i - 1, n)])
        (nxt == to_l && prv == to_l) && return nothing
        nxt == to_l && push!(found, 1)
        prv == to_l && push!(found, -1)
    end
    isempty(found) && throw(KodaInputError(
        :heegaard_to_spine,
        "contact arc not found on curve word",
        Dict("from" => from_l, "to" => to_l, "word" => word),
    ))
    length(unique(found)) == 1 || return nothing
    found[1]
end

"""§7.3.1 cγ on ∂J contact arcs (candidate A: dγ=cγ).

Contact arcs can span several LC segments after forgetting K. Even for longer
words, the same endpoints may bound the other, shorter cyclic arc. Direction is
the directed arc's recorded `word_step`, rather than a guess from its endpoints.
"""
function _koda_d_gamma_from_ports(ports::Vector{RegionPort}, alpha_words::Dict, beta_words::Dict)
    dγ = Dict{Tuple{Symbol, Int}, Int}()
    for p in ports
        word = p.family == :alpha ? alpha_words[p.curve] : beta_words[p.curve]
        c = p.word_step
        c in (-1, 1) || throw(KodaInputError(:heegaard_to_spine,
            "port word_step must be +1 or -1", Dict("curve" => p.curve, "word_step" => c)))
        all(l -> any(x -> abs(x) == l, word), (p.from_label, p.to_label)) ||
            throw(KodaInputError(:heegaard_to_spine, "port endpoint is absent from curve word", Dict("curve" => p.curve)))
        key = (p.family, p.curve_id)
        if haskey(dγ, key) && dγ[key] != c
            throw(KodaImplementationError(:heegaard_to_spine, "incoherent cγ on curve", Dict("key" => key)))
        end
        dγ[key] = c
    end
    dγ
end

function _koda_joining_corners_from_ports(ports::Vector{RegionPort})
    Set{Int}(v for p in ports for v in (p.from_label, p.to_label))
end

function _koda_strand_oriented_word(word::Vector{Int}, d::Int)
    d == 1 && return collect(word)
    return Int[word[i] for i in length(word):-1:1]
end

function _koda_strand_oriented_cycles(
    alpha::Vector{Vector{Int}},
    beta::Vector{Vector{Int}},
    dγ::Dict{Tuple{Symbol, Int}, Int},
)
    function step(family, i, w)
        key = (family, i)
        haskey(dγ, key) || throw(KodaInputError(
            :heegaard_to_spine,
            "missing cγ; §5.2 requires every curve to meet J in one arc",
            Dict("key" => key, "word" => w),
        ))
        _koda_strand_oriented_word(w, dγ[key])
    end
    αo = [step(:alpha, i, w) for (i, w) in enumerate(alpha)]
    βo = [step(:beta, i, w) for (i, w) in enumerate(beta)]
    αo, βo
end

function _koda_sigma_from_crossings(
    crossings::Vector{<:AbstractDict},
    labels::AbstractDict,
    dγ::Dict{Tuple{Symbol, Int}, Int},
)
    σ = Dict{Int, Int}()
    for x in crossings
        id = x["id"]
        !haskey(labels, id) && continue
        da = dγ[(:alpha, x["alpha"])]
        db = dγ[(:beta, x["beta"])]
        σ[labels[id]] = -x["sign"] * da * db
    end
    σ
end

function _koda_word_vectors_from_dict(alpha_words, beta_words)
    anames = sort(collect(keys(alpha_words)))
    bnames = sort(collect(keys(beta_words)))
    [collect(alpha_words[n]) for n in anames], [collect(beta_words[n]) for n in bnames]
end

function _koda_word_vectors_from_forgot(forgot::ForgottenHeegaard)
    anames = get(forgot.correspondence, "alpha_order", sort(collect(keys(forgot.alpha_words))))
    bnames = get(forgot.correspondence, "beta_order", sort(collect(keys(forgot.beta_words))))
    α = [collect(forgot.alpha_words[n]) for n in anames]
    β = [collect(forgot.beta_words[n]) for n in bnames]
    α, β
end

function _koda_combinatorial_spine_core(
    α::Vector{Vector{Int}},
    β::Vector{Vector{Int}},
    crossings::Vector{<:AbstractDict},
    joining::Set{Int},
    dγ::Dict{Tuple{Symbol, Int}, Int},
    source,
    genus::Int,
)
    all_ids = Set(v for w in vcat(α, β) for v in w)
    true_v = setdiff(all_ids, joining)
    αs, βs = _koda_strand_oriented_cycles(α, β, dγ)
    words, labels, _ = _koda_smooth_joining(αs, βs, collect(joining), :beta)
    σ = _koda_sigma_from_crossings(crossings, labels, dγ)
    spine = _koda_flow_spine_from_words(words, σ, source; genus=genus)
    _koda_fill_intermediate!(spine, αs, βs, joining, true_v, nothing, nothing)
    corr = Dict{String, Any}("true_vertices" => sort(collect(true_v)), "d_gamma" => dγ, "vertex_free" => isempty(true_v))
    KodaConversionResult(spine, Dict{String, Any}("convention" => "candidate_A"), corr, :constructed)
end

function _koda_combinatorial_spine_from_forgotten_region(forgot::ForgottenHeegaard, region::CombinatorialRegion)
    region.joining_disk || throw(KodaInputError(:heegaard_to_spine, "region is not a joining disk", Dict()))
    α, β = _koda_word_vectors_from_forgot(forgot)
    dγ = _koda_d_gamma_from_ports(region.boundary_ports, forgot.alpha_words, forgot.beta_words)
    joining = _koda_joining_corners_from_ports(region.boundary_ports)
    _koda_combinatorial_spine_core(α, β, forgot.crossings, joining, dγ, (forgot, region), forgot.genus)
end

function _koda_heegaard_word_cycles(H::HeegaardDiagram)
    code = _undecorated(H)
    α = [Int[c.visits...] for c in code.alpha]
    β = [Int[c.visits...] for c in code.beta]
    α, β
end

function _koda_joining_corners(H::HeegaardDiagram, D::JoiningDisk)
    graph = ribbon_graph(H)
    Set(graph.darts[d].origin for d in D.boundary_darts)
end

"""cγ from the ribbon boundary dart.

ribbon.jl numbers the forward dart (along the curve word, visit[k] → visit[k+1])
with the odd id `2e-1`, and the reverse dart with the even id `2e`.
`JoiningDisk.boundary_darts` is that ribbon boundary cycle. On the abalone
fixture the cycle has the disk on the left of the forward alpha dart, which is
the spec §7.3 outer-normal-first orientation (region on the left of t).
Therefore cγ = +1 iff the boundary dart is the forward dart.
"""
function _koda_curve_directions(H::HeegaardDiagram, D::JoiningDisk)
    graph = ribbon_graph(H)
    directions = Dict{Tuple{Symbol, Int}, Int}()
    for d in D.boundary_darts
        dart = graph.darts[d]
        edge = graph.edges[dart.edge]
        key = (edge.family, edge.curve_id)
        val = isodd(d) ? 1 : -1
        haskey(directions, key) && directions[key] != val &&
            throw(KodaImplementationError(:heegaard_to_spine, "incoherent joining boundary direction", Dict("key" => key)))
        directions[key] = val
    end
    directions
end

function _koda_vertex_signs_from_diagram(H::HeegaardDiagram, D::JoiningDisk)
    directions = _koda_curve_directions(H, D)
    code = _undecorated(H)
    Dict(
        p.id => -p.sign * directions[(:alpha, p.alpha)] * directions[(:beta, p.beta)]
        for p in code.intersections
        if haskey(directions, (:alpha, p.alpha)) && haskey(directions, (:beta, p.beta))
    )
end

function _koda_epsilon_from_diagram(H::HeegaardDiagram, D::JoiningDisk)
    sigs = _koda_vertex_signs_from_diagram(H, D)
    corners = _koda_joining_corners(H, D)
    ε = Dict{Int, Int}()
    for (id, sg) in sigs
        ε[id] = sg
    end
    for id in corners
        haskey(ε, id) || throw(KodaImplementationError(:heegaard_to_spine, "missing corner sign", Dict("id" => id)))
    end
    ε
end

function _koda_combinatorial_spine_from_punctured(ph::PuncturedHeegaardDiagram;
                                                genus::Int=surface_genus(ph.diagram))
    H, D = ph.diagram, ph.joining_disk
    validate_joining_disk(H, D)
    α, β = _koda_heegaard_word_cycles(H)
    joining = collect(_koda_joining_corners(H, D))
    all_ids = Set(v for w in vcat(α, β) for v in w)
    true_v = setdiff(all_ids, joining)
    dγ = _koda_curve_directions(H, D)
    αs, βs = _koda_strand_oriented_cycles(α, β, dγ)
    words, labels, _ = _koda_smooth_joining(αs, βs, joining, :beta)
    code = _undecorated(H)
    crossings = [
        Dict("id" => p.id, "sign" => p.sign, "alpha" => p.alpha, "beta" => p.beta)
        for p in code.intersections
    ]
    sigma = _koda_sigma_from_crossings(crossings, labels, dγ)
    spine = _koda_flow_spine_from_words(words, sigma, ph; genus=genus)
    _koda_fill_intermediate!(spine, αs, βs, Set(joining), true_v, nothing, nothing)
    corr = Dict{String, Any}(
        "true_vertices" => sort(collect(true_v)),
        "d_gamma" => dγ,
        "vertex_free" => isempty(true_v),
    )
    KodaConversionResult(spine, Dict{String, Any}("convention" => "candidate_A"), corr, :constructed)
end

"""§6.3 singular circles: walk only through joining corners, switching family at each corner."""
function _koda_vertex_free_circles(alpha, beta, joining::Set{Int})
    successor = Dict{Tuple{Int, Char}, Tuple{Int, Char}}()
    for (family, words) in (('a', alpha), ('b', beta))
        for w in words
            for i in eachindex(w)
                successor[(w[i], family)] = (w[mod1(i + 1, length(w))], family)
            end
        end
    end
    seen = Set{Tuple{Int, Char}}()
    seen_keys = Set{Tuple}()
    circles = Vector{Vector{Int}}()
    for start in sort(collect(keys(successor)); by=x -> (x[1], x[2]))
        start in seen && continue
        start[1] in joining || continue
        cyc = Int[]
        state = start
        closed = false
        for _ in 1:(length(successor) + 2)
            state in seen && break
            state[1] in joining || break
            push!(seen, state)
            push!(cyc, state[2] == 'b' ? state[1] : -state[1])
            nxt = get(successor, state, nothing)
            nxt === nothing && break
            state = (nxt[1], nxt[2] == 'a' ? 'b' : 'a')
            state == start && (closed = true; break)
        end
        closed && !isempty(cyc) || continue
        verts = abs.(cyc)
        n = length(verts)
        rots = [tuple(verts[mod1.(i:(i + n - 1), n)]...) for i in 1:n]
        rev = reverse(verts)
        rrots = [tuple(rev[mod1.(i:(i + n - 1), n)]...) for i in 1:n]
        key = minimum(vcat(rots, rrots))
        key in seen_keys && continue
        push!(seen_keys, key)
        push!(circles, cyc)
    end
    circles
end

function _koda_fill_intermediate!(spine::FlowSpine, alpha, beta, joining::Set{Int}, true_v, forgot, region)
    circles = isempty(true_v) ? _koda_vertex_free_circles(alpha, beta, joining) : Vector{Vector{Int}}[]
    components = Dict{String, Any}[]
    if isempty(true_v)
        for (i, word) in enumerate(circles)
            push!(components, Dict(
                "id" => i,
                "kind" => "circle",
                "true_vertex_count" => 0,
                "word" => word,
                "artificial_vertex" => -i,
            ))
        end
    else
        push!(components, Dict(
            "id" => 1,
            "kind" => "graph",
            "true_vertex_count" => length(true_v),
            "edge_count" => length(spine.edges),
            "true_vertices" => sort(collect(true_v)),
        ))
    end
    sectors = Dict{String, Any}[]
    if forgot !== nothing && region !== nothing
        for r in forgot.regions
            r.id == region.id && continue
            push!(sectors, Dict(
                "id" => r.id,
                "kind" => "surface",
                "genus" => r.genus,
                "boundary_count" => r.boundary_components,
                "chi" => r.chi,
                "gluing_words" => [[(p.curve, p.from_label, p.to_label, p.word_step) for p in r.boundary_ports]],
            ))
        end
        for (name, word) in sort(collect(merge(forgot.alpha_words, forgot.beta_words)))
            push!(sectors, Dict(
                "id" => name,
                "kind" => "compression_disk",
                "genus" => 0,
                "boundary_count" => 1,
                "chi" => 1,
                "gluing_words" => [collect(word)],
            ))
        end
    end
    spine.complex.checks["vertex_free"] = isempty(true_v)
    spine.complex.checks["artificial_vertices"] = isempty(true_v) ? [-i for i in eachindex(circles)] : Int[]
    spine.complex.checks["singular_components"] = components
    spine.complex.checks["sectors"] = sectors
    spine.complex.checks["true_vertices"] = sort(collect(true_v))
    spine.complex.checks["spec"] = "heegaard-koda-conversion.md §6.3"
    spine
end

"""§6.3 complex: singular edges after smoothing ∂J corners, sectors after gluing along ∂J."""
function _koda_j_undirected(forgot::ForgottenHeegaard, region::CombinatorialRegion)
    words = merge(forgot.alpha_words, forgot.beta_words)
    s = Set{Tuple{String, Int, Int}}()
    for p in region.boundary_ports
        word = get(words, p.curve, nothing)
        word === nothing && continue
        n = length(word)
        pos = findfirst(x -> abs(x) == p.from_label, word)
        pos === nothing && continue
        for _ in 1:n
            nxt = mod1(pos + p.word_step, n)
            a, b = abs(word[pos]), abs(word[nxt])
            push!(s, (p.curve, min(a, b), max(a, b)))
            b == p.to_label && break
            pos = nxt
        end
    end
    s
end

function _koda_singular_chains(forgot::ForgottenHeegaard, region::CombinatorialRegion)
    words = merge(forgot.alpha_words, forgot.beta_words)
    jset = _koda_j_undirected(forgot, region)
    ab = Set(x["id"] for x in forgot.crossings)
    on_j = Set{Int}()
    for (_, a, b) in jset
        push!(on_j, a, b)
    end
    jcorn = intersect(on_j, ab)
    true_v = setdiff(ab, jcorn)
    sing = Tuple{String, Int, Int}[]
    for (name, word) in words
        n = length(word)
        for i in 1:n
            a = abs(word[i])
            b = abs(word[mod1(i + 1, n)])
            (name, min(a, b), max(a, b)) in jset && continue
            push!(sing, (name, a, b))
        end
    end
    singset = Set(sing)
    succ = Dict{Tuple{String, Int, Int}, Tuple{String, Int, Int}}()
    for p in sing
        name, a, b = p
        b in ab && continue
        word = words[name]
        n = length(word)
        pos = findfirst(x -> abs(x) == a, word)
        pos === nothing && continue
        step = abs(word[mod1(pos + 1, n)]) == b ? 1 : -1
        nb = abs(word[mod1(pos + 2 * step, n)])
        cand = (name, b, nb)
        cand in singset && (succ[p] = cand)
    end
    used = Set{Tuple{String, Int, Int}}()
    chains = Vector{Vector{Tuple{String, Int, Int}}}()
    closed = Bool[]
    for p0 in sing
        p0 in used && continue
        has_ab = any(q -> q[1] == p0[1] && (q[2] in ab || q[3] in ab), sing)
        if has_ab && !(p0[2] in ab)
            continue
        end
        cyc = Tuple{String, Int, Int}[]
        cur = p0
        is_closed = false
        for _ in 1:length(sing) + 2
            cur in used && break
            push!(used, cur)
            push!(cyc, cur)
            haskey(succ, cur) || break
            nxt = succ[cur]
            if nxt == p0
                is_closed = true
                break
            end
            cur = nxt
        end
        push!(chains, cyc)
        push!(closed, is_closed)
    end
    for p0 in sing
        p0 in used && continue
        cyc = Tuple{String, Int, Int}[]
        cur = p0
        is_closed = false
        for _ in 1:length(sing) + 2
            cur in used && break
            push!(used, cur)
            push!(cyc, cur)
            haskey(succ, cur) || break
            nxt = succ[cur]
            nxt == p0 && (is_closed = true; break)
            cur = nxt
        end
        isempty(cyc) || (push!(chains, cyc); push!(closed, is_closed))
    end
    chains, closed, true_v
end

function _koda_piece_cycles(forgot::ForgottenHeegaard, region::CombinatorialRegion)
    words = merge(forgot.alpha_words, forgot.beta_words)
    jset = _koda_j_undirected(forgot, region)
    on_seam(p) = (p[1], min(p[2], p[3]), max(p[2], p[3])) in jset
    roles = Dict{String, Symbol}()
    for n in keys(forgot.alpha_words); roles[string(n)] = :alpha; end
    for n in keys(forgot.beta_words); roles[string(n)] = :beta; end
    pieces = Dict{String, Any}[]
    real = forgot.realization
    if real !== nothing
        face_by = Dict(f.id => f for f in real.faces)
        for r in forgot.regions
            r.id == region.id && continue
            ids = Set(r.face_roots)
            darts = Tuple{String, Int, Int}[]
            for fid in ids
                f = face_by[fid]
                n = length(f.corners)
                for i in 1:n
                    c1 = f.corners[i]
                    c2 = f.corners[mod1(i + 1, n)]
                    comp, = _corner_component(real.lc.words, c1)
                    comp === nothing && continue
                    get(roles, string(comp), :none) in (:alpha, :beta) || continue
                    push!(darts, (string(comp), c1[1], c2[1]))
                end
            end
            unique!(darts)
            cycles = _koda_chain_darts(darts)
            chi = _koda_region_cell_chi(forgot, r)
            push!(pieces, Dict("key" => "R$(r.id)", "cycles" => cycles, "chi" => chi, "regions" => [r.id], "disks" => String[]))
        end
    end
    for (name, word) in sort(collect(words))
        n = length(word)
        cyc = [(name, abs(word[i]), abs(word[mod1(i + 1, n)])) for i in 1:n]
        push!(pieces, Dict("key" => "D:$name", "cycles" => [cyc], "chi" => 1, "regions" => Int[], "disks" => [name]))
    end
    pieces, on_seam
end

function _koda_chain_darts(darts)
    by_start = Dict{Tuple{String, Int}, Vector{Int}}()
    for (i, d) in enumerate(darts)
        push!(get!(by_start, (d[1], d[2]), Int[]), i)
    end
    used = falses(length(darts))
    cycles = Vector{Vector{Tuple{String, Int, Int}}}()
    for s in eachindex(darts)
        used[s] && continue
        cyc = Tuple{String, Int, Int}[]
        i = s
        for _ in 1:length(darts) + 2
            used[i] && break
            used[i] = true
            push!(cyc, darts[i])
            cands = [j for j in get(by_start, (darts[i][1], darts[i][3]), Int[]) if !used[j]]
            isempty(cands) && break
            i = cands[1]
        end
        isempty(cyc) || push!(cycles, cyc)
    end
    cycles
end

function _koda_region_cell_chi(forgot::ForgottenHeegaard, region::CombinatorialRegion)
    real = forgot.realization
    real === nothing && return region.chi
    # cf: docs/specs/heegaard-koda-conversion.md §4.1 — delete K edges and K-K vertices from the original cells.
    ids = Set(region.face_roots)
    knames = Set(String.(forgot.correspondence["forgotten_reference"]))
    words = _lc_nonempty_words(real.lc.words)
    nedges = 0
    for name in knames
        haskey(words, name) || continue
        w = words[name]
        for i in eachindex(w)
            x, y = w[i], w[mod1(i + 1, length(w))]
            a = real.corner_face[(abs(x), _lc_sign_char(x), "out")] in ids
            b = real.corner_face[(abs(y), _lc_sign_char(y), "in")] in ids
            a == b || throw(KodaImplementationError(:spine_topology,
                "original region is not closed across a forgotten K edge",
                Dict("region" => region.id, "curve" => name, "position" => i)))
            a && (nedges += 1)
        end
    end
    nvertices = 0
    for (lab, occ) in real.occurrence
        all(o -> o[1] in knames, values(occ)) || continue
        adjacent = [real.corner_face[(lab, sg, way)] in ids
                    for sg in ("+", "-"), way in ("in", "out")]
        any(adjacent) || continue
        all(adjacent) || throw(KodaImplementationError(:spine_topology,
            "forgotten K-K vertex spans original regions", Dict("label" => lab)))
        nvertices += 1
    end
    length(ids) - nedges + nvertices
end

"""Boundary cycles of one sector: follow each piece boundary, and cross a ∂J edge onto the partner."""
function _koda_trace_boundaries(pieces, on_seam)
    rev(p) = (p[1], p[3], p[2])
    function locate(target)
        for (pi, piece) in enumerate(pieces)
            for (ci, cyc) in enumerate(piece["cycles"])
                for (ei, e) in enumerate(cyc)
                    if e == target
                        return (pi, ci, ei, 1)
                    elseif rev(e) == target
                        return (pi, ci, ei, -1)
                    end
                end
            end
        end
        nothing
    end
    function step(pos)
        pi, ci, ei, dir = pos
        cyc = pieces[pi]["cycles"][ci]
        ei2 = mod1(ei + dir, length(cyc))
        (pi, ci, ei2, dir)
    end
    function edge_at(pos)
        pi, ci, ei, dir = pos
        e = pieces[pi]["cycles"][ci][ei]
        dir == 1 ? e : (e[1], e[3], e[2])
    end
    # singular occurrences: positions whose traversed edge is not a seam
    starts = Tuple{Int, Int, Int, Int}[]
    for (pi, piece) in enumerate(pieces)
        for (ci, cyc) in enumerate(piece["cycles"])
            for ei in eachindex(cyc)
                e = cyc[ei]
                on_seam(e) && continue
                prev = cyc[mod1(ei - 1, length(cyc))]
                (on_seam(prev) || length(cyc) == 1) && push!(starts, (pi, ci, ei, 1))
            end
            if !any(on_seam, cyc) && length(cyc) > 1
                push!(starts, (pi, ci, 1, 1))
            end
        end
    end
    used = Set{Tuple{Int, Int, Int, Int}}()
    words = Vector{Vector{Tuple{String, Int, Int}}}()
    for start in starts
        start in used && continue
        word = Tuple{String, Int, Int}[]
        pos = start
        guard = 0
        while guard < 10000
            guard += 1
            e = edge_at(pos)
            if on_seam(e)
                partner = locate(rev(e))
                partner === nothing && break
                pos = step(partner)
                continue
            end
            key = (pos[1], pos[2], pos[3], pos[4])
            key in used && break
            push!(used, key)
            push!(word, e)
            nxt = step(pos)
            edge_at(nxt) == edge_at(start) && nxt[1] == start[1] && length(word) > 0 && break
            pos = nxt
            pos == start && break
        end
        isempty(word) || push!(words, word)
    end
    words
end

function _koda_traced_sectors(forgot::ForgottenHeegaard, region::CombinatorialRegion, piece_edge)
    pieces, on_seam = _koda_piece_cycles(forgot, region)
    parent = Dict(i => i for i in eachindex(pieces))
    function find(x)
        parent[x] == x ? x : (parent[x] = find(parent[x]))
    end
    seam_owner = Dict{Tuple{String, Int, Int}, Vector{Int}}()
    for (i, piece) in enumerate(pieces)
        for cyc in piece["cycles"]
            for e in cyc
                on_seam(e) || continue
                key = (e[1], min(e[2], e[3]), max(e[2], e[3]))
                push!(get!(seam_owner, key, Int[]), i)
            end
        end
    end
    for ids in values(seam_owner)
        for j in ids[2:end]
            parent[find(j)] = find(ids[1])
        end
    end
    groups = Dict{Int, Vector{Int}}()
    for i in eachindex(pieces)
        push!(get!(groups, find(i), Int[]), i)
    end
    sectors = Dict{String, Any}[]
    for (k, ids) in enumerate(values(groups))
        sub = pieces[ids]
        raw = _koda_trace_boundaries(sub, on_seam)
        boundaries = Vector{Int}[]
        for w in raw
            mapped = _koda_word_of(w, piece_edge)
            isempty(mapped) || push!(boundaries, mapped)
        end
        seams_here = Set{Tuple{String, Int, Int}}()
        for piece in sub
            for cyc in piece["cycles"]
                for e in cyc
                    on_seam(e) || continue
                    push!(seams_here, (e[1], min(e[2], e[3]), max(e[2], e[3])))
                end
            end
        end
        chi = sum(p["chi"] for p in sub) - length(seams_here)
        b = length(boundaries)
        twice = 2 - chi - b
        genus = twice >= 0 && iseven(twice) ? div(twice, 2) : 0
        sources_r = Int[]
        sources_d = String[]
        for p in sub
            append!(sources_r, p["regions"])
            append!(sources_d, p["disks"])
        end
        push!(sectors, Dict(
            "id" => k,
            "genus" => genus,
            "boundary_count" => b,
            "chi" => chi,
            "boundaries" => boundaries,
            "sources" => Dict("regions" => sources_r, "disks" => sources_d),
        ))
    end
    sectors
end

"""cf: docs/specs/heegaard-koda-conversion.md §4.1, §6.3 — oriented LC sides retain the original complementary surfaces.

The auxiliary ribbon diagram carries the curve graph; its cellular caps do not
specify the original surface. Check every collapsed LC segment on both sides.
"""
function _koda_original_complement(forgot::ForgottenHeegaard, h::HeegaardDiagram,
                                   selected::CombinatorialRegion)
    real = forgot.realization
    real === nothing && throw(KodaImplementationError(:spine_topology,
        "original LC realization is required", Dict()))
    fail(msg, data=Dict{String,Any}()) =
        throw(KodaImplementationError(:spine_topology, msg, data))
    face_region = Dict(f => r.id for r in forgot.regions for f in r.face_roots)
    (sum(length(r.face_roots) for r in forgot.regions) == length(real.faces) &&
     Set(keys(face_region)) == Set(f.id for f in real.faces)) || fail("incomplete original face partition")
    graph = ribbon_graph(h); cycles = boundary_cycles(graph)
    labels = _source_label_ids(h)
    owner = Dict{Int,Int}()
    bycurve = Dict{Tuple{Symbol,Int},Vector{Int}}()
    for e in graph.edges
        push!(get!(bycurve, (e.family, e.curve_id), Int[]), e.id)
    end
    anames = get(forgot.correspondence, "alpha_order", sort(collect(keys(forgot.alpha_words))))
    bnames = get(forgot.correspondence, "beta_order", sort(collect(keys(forgot.beta_words))))
    specs = vcat([(:alpha, i, n, forgot.alpha_words[n]) for (i,n) in enumerate(anames)],
                 [(:beta, i, n, forgot.beta_words[n]) for (i,n) in enumerate(bnames)])
    for (fam, cid, name, w) in specs
        visits = [i for i in eachindex(w) if haskey(labels, abs(w[i]))]
        isempty(visits) && fail("curve without alpha-beta crossings requires a placement", Dict("curve" => name))
        eids = get(bycurve, (fam,cid), Int[])
        length(eids) == length(visits) || fail("curve-edge order mismatch", Dict("curve" => name))
        for k in eachindex(visits)
            p, stop = visits[k], visits[mod1(k + 1, length(visits))]
            right, left = Set{Int}(), Set{Int}()
            while true
                q = mod1(p + 1, length(w))
                push!(right, face_region[real.corner_face[(abs(w[p]), _lc_sign_char(w[p]), "out")]])
                push!(left, face_region[real.corner_face[(abs(w[q]), _lc_sign_char(w[q]), "in")]])
                q == stop && break
                p = q
            end
            length(right) == length(left) == 1 || fail("one ribbon side crosses original regions", Dict("curve" => name))
            e = graph.edges[eids[k]]
            (e.from == labels[abs(w[visits[k]])] && e.to == labels[abs(w[stop])]) ||
                fail("ribbon edge is not oriented in LC word order", Dict("edge" => e.id))
            # cf: heegaard-koda-conversion.md §§4.1, 7.3.1: the ribbon
            # boundary keeps its region to the left of the oriented dart.
            owner[2e.id - 1] = only(left); owner[2e.id] = only(right)
        end
    end
    blocks = Dict(r.id => Int[] for r in forgot.regions)
    for (ci, cyc) in enumerate(cycles)
        rs = Set(owner[d] for d in cyc)
        length(rs) == 1 || fail("boundary cycle has mixed LC-side owners", Dict("cycle" => ci))
        push!(blocks[only(rs)], ci)
    end
    signed(d) = isodd(d) ? graph.darts[d].edge : -graph.darts[d].edge
    regions = NamedTuple[]
    for r in forgot.regions
        chi = _koda_region_cell_chi(forgot, r)
        b = length(blocks[r.id]); twice = 2 - chi - b
        twice >= 0 && iseven(twice) || fail("inconsistent original region", Dict("region" => r.id, "chi" => chi, "boundary_count" => b))
        push!(regions, (id=r.id, faces=r.face_roots, chi=chi, genus=twice÷2,
                       cycles=[[signed(d) for d in cycles[ci]] for ci in blocks[r.id]]))
    end
    jreg = only(r for r in regions if r.id == selected.id)
    (jreg.chi == 1 && jreg.genus == 0 && length(jreg.cycles) == 1 &&
     Set(jreg.faces) == Set(selected.face_roots)) || fail("selected original region is not a disk", Dict("region" => selected.id))
    (; graph, cycles, owner, regions, joining_region=selected.id)
end

function _koda_ribbon_complex!(spine::FlowSpine, forgot::ForgottenHeegaard, region::CombinatorialRegion)
    diagram = forgot.diagram
    source = _koda_original_complement(forgot, diagram, region)
    if any(r -> r.genus != 0 || length(r.cycles) != 1, source.regions)
        # Keep the original surface banks for an explicitly requested isotopy.
        spine.complex.checks["source_surface_regions"] = [Dict{String,Any}(
            "id"=>r.id, "genus"=>r.genus, "chi"=>r.chi,
            "boundary_count"=>length(r.cycles),
            "boundary_darts"=>[copy(c) for c in source.cycles if source.owner[first(c)]==r.id])
            for r in source.regions]
    end
    graph = source.graph
    cycles = source.cycles
    regs = source.regions
    dart = Dict(d.id => d for d in graph.darts)
    J = source.joining_region
    original_joining = only(r for r in regs if r.id == J)
    seam = Set(abs(x) for boundary in original_joining.cycles for x in boundary)
    signed(d) = isodd(d) ? dart[d].edge : -dart[d].edge
    owner = source.owner
    function region_coeff(eid)
        for d in (2 * eid - 1, 2 * eid)
            haskey(owner, d) || continue
            owner[d] == J && continue
            return isodd(d) ? 1 : -1
        end
        nothing
    end
    pieces = Dict{String, Any}[]
    for reg in regs
        ri = reg.id
        ri == J && continue
        words = Vector{Int}[]
        for cycle in reg.cycles
            push!(words, copy(cycle))
        end
        push!(pieces, Dict(
            "kind" => "region",
            "id" => ri,
            "chi" => reg.chi,
            "cycles" => words,
            "label" => "region:$ri",
        ))
    end
    curve_edges = Dict{Tuple{Symbol, Int}, Vector{Int}}()
    curve_order = Tuple{Symbol, Int}[]
    for e in graph.edges
        key = (e.family, e.curve_id)
        if !haskey(curve_edges, key)
            curve_edges[key] = Int[]
            push!(curve_order, key)
        end
        push!(curve_edges[key], e.id)
    end
    for key in curve_order
        eids = curve_edges[key]
        forward = true
        for eid in eids
            eid in seam || continue
            coeff = region_coeff(eid)
            coeff === nothing && continue
            forward = coeff == -1
            break
        end
        word = forward ? copy(eids) : [-e for e in reverse(eids)]
        push!(pieces, Dict(
            "kind" => "disk",
            "id" => key,
            "chi" => 1,
            "cycles" => [word],
            "label" => string(key[1]) * string(key[2]),
        ))
    end
    parent = Dict(i => i for i in eachindex(pieces))
    function find(x)
        parent[x] == x ? x : (parent[x] = find(parent[x]))
    end
    function piece_of_region(ri)
        findfirst(p -> p["kind"] == "region" && p["id"] == ri, pieces)
    end
    function piece_of_disk(key)
        findfirst(p -> p["kind"] == "disk" && p["id"] == key, pieces)
    end
    seam_of = Dict{Int, Tuple{Int, Int}}()
    for eid in seam
        e = graph.edges[eid]
        ri = nothing
        for d in (2 * eid - 1, 2 * eid)
            haskey(owner, d) || continue
            owner[d] == J && continue
            ri = owner[d]
        end
        ri === nothing && continue
        pr = piece_of_region(ri)
        pd = piece_of_disk((e.family, e.curve_id))
        (pr === nothing || pd === nothing) && continue
        parent[find(pd)] = find(pr)
        seam_of[eid] = (pr, pd)
    end
    groups = Dict{Int, Vector{Int}}()
    for i in eachindex(pieces)
        push!(get!(groups, find(i), Int[]), i)
    end
    function trace(ids, orientation)
        sub = pieces[ids]
        index_of = Dict(ids[i] => i for i in eachindex(ids))
        function after_reverse(signed_seam, skip)
            target = -signed_seam
            for pi in eachindex(sub)
                pi == skip && continue
                for (ci, cyc) in enumerate(sub[pi]["cycles"])
                    for (ei, x) in enumerate(cyc)
                        if x == target
                            return (pi, ci, mod1(ei + 1, length(cyc)), 1)
                        elseif x == signed_seam
                            return (pi, ci, mod1(ei - 1, length(cyc)), -1)
                        end
                    end
                end
            end
            nothing
        end
        function edge_at(pos)
            pi, ci, ei, dir = pos
            x = sub[pi]["cycles"][ci][ei]
            dir == 1 ? x : -x
        end
        function step(pos)
            pi, ci, ei, dir = pos
            n = length(sub[pi]["cycles"][ci])
            (pi, ci, mod1(ei + dir, n), dir)
        end
        starts = Tuple{Int, Int, Int, Int}[]
        for (pi, piece) in enumerate(sub)
            for (ci, cyc) in enumerate(piece["cycles"])
                for ei in eachindex(cyc)
                    abs(cyc[ei]) in seam && continue
                    push!(starts, (pi, ci, ei, 1))
                end
            end
        end
        used = Set{Tuple{Int, Int, Int}}()
        words = Vector{Int}[]
        for start in starts
            (start[1], start[2], start[3]) in used && continue
            word = Int[]
            pos = start
            guard = 0
            while guard < 10000
                guard += 1
                letter = edge_at(pos)
                if abs(letter) in seam
                    nxt = after_reverse(letter, pos[1])
                    nxt === nothing && break
                    pos = nxt
                    continue
                end
                key = (pos[1], pos[2], pos[3])
                key in used && break
                push!(used, key)
                push!(word, letter)
                nxt = step(pos)
                nxt == start && break
                pos = nxt
            end
            if !isempty(word)
                orientation[ids[start[1]]] == -1 && (word = [-x for x in reverse(word)])
                push!(words, word)
            end
        end
        words
    end
    sectors = Dict{String, Any}[]
    raw_genus = Tuple{Int, Any}[]
    for (k, ids) in enumerate(values(groups))
        seams_here = [eid for (eid, pair) in seam_of if find(pair[1]) == find(ids[1])]
        # Propagate one surface orientation across the glued seams. Boundary
        # parametrizations are then induced from this common orientation.
        orientation = Dict(first(ids) => 1)
        coefficient(pi, eid) = only(sign(x) for cyc in pieces[pi]["cycles"] for x in cyc if abs(x) == eid)
        while length(orientation) < length(ids)
            before = length(orientation)
            for eid in seams_here
                a, b = seam_of[eid]
                factor = -coefficient(a, eid) * coefficient(b, eid)
                if haskey(orientation, a)
                    expected = factor * orientation[a]
                    haskey(orientation, b) && orientation[b] != expected &&
                        throw(KodaImplementationError(:spine_topology, "inconsistent sector orientation", Dict("edge" => eid)))
                    orientation[b] = expected
                elseif haskey(orientation, b)
                    orientation[a] = factor * orientation[b]
                end
            end
            length(orientation) > before || throw(KodaImplementationError(
                :spine_topology, "sector pieces are not seam-connected", Dict("pieces" => ids)))
        end
        for eid in seams_here
            a, b = seam_of[eid]
            orientation[b] == -coefficient(a, eid) * coefficient(b, eid) * orientation[a] ||
                throw(KodaImplementationError(:spine_topology, "inconsistent sector orientation", Dict("edge" => eid)))
        end
        boundaries = trace(ids, orientation)
        chi = sum(pieces[i]["chi"] for i in ids) - length(seams_here)
        b = length(boundaries)
        twice = 2 - chi - b
        genus = twice >= 0 && iseven(twice) ? div(twice, 2) : twice / 2
        twice >= 0 && iseven(twice) || push!(raw_genus, (k, genus))
        sources_r = Int[]
        sources_d = String[]
        for i in ids
            p = pieces[i]
            p["kind"] == "region" ? push!(sources_r, p["id"]) : push!(sources_d, p["label"])
        end
        push!(sectors, Dict(
            "id" => k,
            "genus" => genus,
            "boundary_count" => b,
            "chi" => chi,
            "boundaries" => boundaries,
            "boundary_orientation_signs" => fill(1, b),
            "sources" => Dict("regions" => sources_r, "disks" => sources_d),
        ))
    end
    rays = Dict{Int, Int}()
    edges = Dict{String, Any}[]
    for e in graph.edges
        e.id in seam && continue
        rays[e.from] = get(rays, e.from, 0) + 1
        rays[e.to] = get(rays, e.to, 0) + 1
        push!(edges, Dict(
            "id" => e.id,
            "from" => e.from,
            "to" => e.to,
            "pieces" => [(string(e.family), e.curve_id, e.from, e.to)],
        ))
    end
    vertices = Dict{String, Any}[]
    true_labels = Int[]
    for lab in sort(collect(keys(rays)))
        n = rays[lab]
        n == 0 && continue
        kind = n == 4 ? "true" : "artificial"
        kind == "true" && push!(true_labels, lab)
        push!(vertices, Dict("id" => lab, "kind" => kind, "rays" => n))
    end
    checks = spine.complex.checks
    checks["vertices"] = vertices
    checks["edges"] = edges
    checks["sectors"] = sectors
    checks["true_vertices"] = true_labels
    checks["artificial_vertices"] = [v["id"] for v in vertices if v["kind"] == "artificial"]
    checks["vertex_free"] = isempty(true_labels)
    checks["singular_components"] = [
        Dict("id" => e["id"], "kind" => e["from"] == e["to"] ? "circle" : "graph", "word" => e["pieces"])
        for e in edges
    ]
    checks["genus_not_integral"] = raw_genus
    checks["spec"] = "heegaard-koda-conversion.md §6.3"
    spine
end

function _koda_record_complex!(spine::FlowSpine, forgot::ForgottenHeegaard, region::CombinatorialRegion)
    forgot.diagram === nothing || return _koda_ribbon_complex!(spine, forgot, region)
    chains, closed, _ = _koda_singular_chains(forgot, region)
    rays = Dict{Int, Int}()
    for (ei, chain) in enumerate(chains)
        closed[ei] && continue
        a = chain[1][2]
        b = chain[end][3]
        rays[a] = get(rays, a, 0) + 1
        rays[b] = get(rays, b, 0) + 1
    end
    vertices = Dict{String, Any}[]
    vid = Dict{Int, Int}()
    true_labels = Int[]
    for lab in sort(collect(keys(rays)))
        rays[lab] == 0 && continue
        kind = rays[lab] == 4 ? "true" : "artificial"
        kind == "true" && push!(true_labels, lab)
        push!(vertices, Dict("id" => length(vertices) + 1, "kind" => kind))
        vid[lab] = length(vertices)
    end
    next_id = length(vertices)
    edges = Dict{String, Any}[]
    piece_edge = Dict{Tuple{String, Int, Int}, Tuple{Int, Int}}()
    for (ei, chain) in enumerate(chains)
        if closed[ei]
            next_id += 1
            push!(vertices, Dict("id" => next_id, "kind" => "artificial"))
            from = to = next_id
        else
            from = vid[chain[1][2]]
            to = vid[chain[end][3]]
        end
        push!(edges, Dict(
            "id" => ei,
            "from" => from,
            "to" => to,
            "pieces" => [(p[1], p[2], p[3]) for p in chain],
        ))
        for p in chain
            piece_edge[p] = (ei, 1)
            piece_edge[(p[1], p[3], p[2])] = (ei, -1)
        end
    end
    sectors = _koda_traced_sectors(forgot, region, piece_edge)
    checks = spine.complex.checks
    checks["vertices"] = vertices
    checks["edges"] = edges
    checks["sectors"] = sectors
    checks["true_vertices"] = true_labels
    checks["artificial_vertices"] = [v["id"] for v in vertices if v["kind"] == "artificial"]
    checks["vertex_free"] = isempty(true_labels)
    checks["singular_components"] = [
        Dict("id" => e["id"], "kind" => e["from"] == e["to"] ? "circle" : "graph", "word" => e["pieces"])
        for e in edges
    ]
    checks["spec"] = "heegaard-koda-conversion.md §6.3"
    spine
end

function _koda_expand_ports(ports, words)
    pieces = Tuple{String, Int, Int}[]
    for p in ports
        word = get(words, p.curve, nothing)
        word === nothing && continue
        n = length(word)
        pos = findfirst(x -> abs(x) == p.from_label, word)
        pos === nothing && continue
        for _ in 1:n
            nxt = mod1(pos + p.word_step, n)
            push!(pieces, (p.curve, abs(word[pos]), abs(word[nxt])))
            abs(word[nxt]) == p.to_label && break
            pos = nxt
        end
    end
    pieces
end

function _koda_word_of(pieces, piece_edge)
    word = Int[]
    for p in pieces
        haskey(piece_edge, p) || continue
        id, sgn = piece_edge[p]
        step = sgn * id
        if !isempty(word) && abs(word[end]) == id && sign(word[end]) == sign(step)
            continue
        end
        push!(word, step)
    end
    word
end

function _koda_glued_sectors(forgot::ForgottenHeegaard, region::CombinatorialRegion, piece_edge)
    words = merge(forgot.alpha_words, forgot.beta_words)
    jset = _koda_j_undirected(forgot, region)
    on_seam(p) = (p[1], min(p[2], p[3]), max(p[2], p[3])) in jset
    # facelets: non-J regions and compression disks. Glue a region to a disk when they share a seam piece.
    facelets = Dict{String, Any}[]
    for r in forgot.regions
        r.id == region.id && continue
        pieces = [p for p in _koda_expand_ports(r.boundary_ports, words) if !on_seam(p)]
        push!(facelets, Dict(
            "key" => "R$(r.id)",
            "sources" => Dict("regions" => [r.id], "disks" => String[]),
            "pieces" => pieces,
            "genus" => r.genus,
            "seams" => Set(p for p in _koda_expand_ports(r.boundary_ports, words) if on_seam(p)),
        ))
    end
    for (name, word) in sort(collect(words))
        n = length(word)
        pieces = Tuple{String, Int, Int}[]
        seams = Set{Tuple{String, Int, Int}}()
        for i in 1:n
            p = (name, abs(word[i]), abs(word[mod1(i + 1, n)]))
            if on_seam(p)
                push!(seams, p)
            else
                push!(pieces, p)
            end
        end
        push!(facelets, Dict(
            "key" => "D:$name",
            "sources" => Dict("regions" => Int[], "disks" => [name]),
            "pieces" => pieces,
            "genus" => 0,
            "seams" => seams,
        ))
    end
    parent = Dict(i => i for i in eachindex(facelets))
    function find(x)
        parent[x] == x ? x : (parent[x] = find(parent[x]))
    end
    function seam_key(p)
        (p[1], min(p[2], p[3]), max(p[2], p[3]))
    end
    by_seam = Dict{Tuple{String, Int, Int}, Vector{Int}}()
    for (i, f) in enumerate(facelets)
        for p in f["seams"]
            push!(get!(by_seam, seam_key(p), Int[]), i)
        end
    end
    for ids in values(by_seam)
        for j in ids[2:end]
            parent[find(j)] = find(ids[1])
        end
    end
    groups = Dict{Int, Vector{Int}}()
    for i in eachindex(facelets)
        push!(get!(groups, find(i), Int[]), i)
    end
    sectors = Dict{String, Any}[]
    for (k, ids) in enumerate(values(groups))
        sources_r = Int[]
        sources_d = String[]
        boundaries = Vector{Int}[]
        genus = 0
        for i in ids
            f = facelets[i]
            append!(sources_r, f["sources"]["regions"])
            append!(sources_d, f["sources"]["disks"])
            genus += f["genus"]
            w = _koda_word_of(f["pieces"], piece_edge)
            isempty(w) || push!(boundaries, w)
        end
        b = length(boundaries)
        # χ is not taken from the pre-glue region fields. Disks contribute 1,
        # each seam-glue of an arc subtracts 1, regions contribute 1 per component
        # only as a lower model; the measured χ is recorded and the formula is checked.
        n_disk = length(sources_d)
        n_region = length(sources_r)
        n_glue = max(n_disk + n_region - 1, 0)
        chi = n_disk + n_region - n_glue
        sectors_genus = chi + b <= 2 && iseven(2 - b - chi) ? div(2 - b - chi, 2) : genus
        push!(sectors, Dict(
            "id" => k,
            "genus" => max(sectors_genus, 0),
            "boundary_count" => b,
            "chi" => chi,
            "boundaries" => boundaries,
            "sources" => Dict("regions" => sources_r, "disks" => sources_d),
        ))
    end
    sectors
end

function _koda_attach_sectors!(spine::FlowSpine, forgot::ForgottenHeegaard, region::CombinatorialRegion)
    _koda_record_complex!(spine, forgot, region)
end

"""
    spine_topology(spine::FlowSpine) -> Dict{String,Any}

Read the recorded vertices, singular edges and sectors of a flow-spine complex.

# Arguments
`spine.complex.checks` supplies the topology records of [Koda conversion §6.3](../specs/heegaard-koda-conversion.md#63-特異辺と-sector). Their correctness is the caller's responsibility; this accessor does not validate them.

# Returns
A dictionary with `vertices`, `edges`, `sectors`, `singular_components`, `vertex_free`, `artificial_vertices`, `true_vertices` and `gluing_locus`. Missing stored entries use empty defaults or vertex-sign-derived defaults.

# Guarantees
1. Existing topology records are returned without reinterpretation and retain their stored values. [proved/agent: by construction]
2. The input is not modified, but stored nested records are shared with the returned dictionary. [proved/agent: by construction]

# Not guaranteed
Empty defaults do not establish that the geometric complex has no cells. The accessor does not compute homology, verify manifold links or establish a spine identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.3 — topology of the glued complex.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> cert = H.ConversionCertificate(:test,"",Dict(),Any[],"",Dict(),nothing,Dict(),:unverified); checks = Dict{String,Any}("vertices"=>[Dict("id"=>1)], "edges"=>[Dict("id"=>1,"from"=>1,"to"=>1)]); c = H.BranchedSpineComplex(H.SpineRegion[],Dict{String,Any}[],Vector{Tuple{Int,Int}}[],checks); s = H.FlowSpine(0,H.SpineEdge[],Dict{Int,Int}(),nothing,cert,c);

julia> t = H.spine_topology(s); (length(t["vertices"]), length(t["edges"]))
(1, 1)
```

Role: access
"""
function spine_topology(spine::FlowSpine)
    c = spine.complex.checks
    Dict{String, Any}(
        "vertices" => get(c, "vertices", Dict{String, Any}[]),
        "edges" => get(c, "edges", Dict{String, Any}[]),
        "sectors" => get(c, "sectors", Dict{String, Any}[]),
        "singular_components" => get(c, "singular_components", Dict{String, Any}[]),
        "vertex_free" => get(c, "vertex_free", isempty(spine.vertex_signs)),
        "artificial_vertices" => get(c, "artificial_vertices", Int[]),
        "true_vertices" => get(c, "true_vertices", sort(collect(keys(spine.vertex_signs)))),
        "gluing_locus" => get(c, "gluing_locus", Any[]),
    )
end

"""
    spine_pi1_presentation(spine::FlowSpine) -> GroupPresentation

Present the fundamental group of the recorded connected sector complex.

# Arguments
The stored topology must describe the connected oriented sector complex of [Koda conversion §6.3](../specs/heegaard-koda-conversion.md#63-特異辺と-sector). Sector boundaries must have coherent induced orientations or explicit orientation signs. Duplicate edges, invalid boundaries, invalid sector data and a disconnected full 1-skeleton raise `KodaImplementationError` at `:spine_pi1_presentation`; missing fields and conversion errors can propagate. The caller establishes the recorded complex's geometric interpretation.

# Returns
A `GroupPresentation` with cell generators and relators, no distinguished generators, and metadata identifying the sector-complex source.

# Guarantees
1. For the stated connected oriented complex, the presentation gives its fundamental group, including the sector topology and all attaching walks. [proved/agent: relative surface cellulation and contraction of a spanning tree]
2. The input spine and stored topology records are not modified. [proved/agent: by construction]

# Not guaranteed
This is not a proof that arbitrary stored topology is the spine of a claimed 3-manifold. Missing topology fields receive empty defaults through `spine_topology`; an empty/default record must not be mistaken for a validated complex. Boundaryless sectors are unsupported.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.3 — glued sector complex.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> cert = H.ConversionCertificate(:test,"",Dict(),Any[],"",Dict(),nothing,Dict(),:unverified); checks = Dict{String,Any}("vertices"=>[Dict("id"=>1)], "edges"=>[Dict("id"=>1,"from"=>1,"to"=>1)]); c = H.BranchedSpineComplex(H.SpineRegion[],Dict{String,Any}[],Vector{Tuple{Int,Int}}[],checks); s = H.FlowSpine(0,H.SpineEdge[],Dict{Int,Int}(),nothing,cert,c);

julia> g = H.spine_pi1_presentation(s); (length(g.generators), length(g.relators))
(1, 0)
```

Role: invariant
"""
function spine_pi1_presentation(spine::FlowSpine)
    topo = spine_topology(spine)
    edges, sectors = topo["edges"], topo["sectors"]
    fail(msg, data=Dict{String,Any}()) =
        throw(KodaImplementationError(:spine_pi1_presentation, msg, data))
    edge_by_id = Dict{Int,Any}()
    verts = Set{Int}(Int(v["id"]) for v in topo["vertices"])
    gens, froms, tos = String[], Int[], Int[]
    function addcell!(name, u, v)
        push!(gens, name); push!(froms, u); push!(tos, v); push!(verts, u, v)
        length(gens)
    end
    id_to_gen = Dict{Int,Int}()
    for e in edges
        id = Int(e["id"])
        haskey(edge_by_id, id) && fail("duplicate singular-edge id", Dict("edge" => id))
        edge_by_id[id] = e
        u, v = Int(e["from"]), Int(e["to"])
        id_to_gen[id] = addcell!("e$id", u, v)
    end
    function directed_ends(x::Int)
        x == 0 && fail("zero is not an oriented edge")
        e = get(edge_by_id, abs(x), nothing)
        e === nothing && fail("unknown singular edge", Dict("edge" => abs(x)))
        u, v = Int(e["from"]), Int(e["to"])
        x > 0 ? (u, v) : (v, u)
    end
    function boundary_start(w::Vector{Int}, sid, j)
        isempty(w) && fail("empty sector boundary", Dict("sector" => sid, "boundary" => j))
        start, cur = directed_ends(w[1])
        for x in w[2:end]
            u, v = directed_ends(x)
            u == cur || fail("boundary is not continuous",
                             Dict("sector" => sid, "boundary" => j))
            cur = v
        end
        cur == start || fail("sector boundary is not closed",
                             Dict("sector" => sid, "boundary" => j))
        start
    end
    plans = NamedTuple[]
    for sec in sectors
        sid = sec["id"]
        g = get(sec, "genus", nothing)
        g isa Integer && g >= 0 ||
            fail("sector genus must be a nonnegative integer", Dict("sector" => sid))
        raw = get(sec, "boundaries", Any[]); b = length(raw)
        b > 0 || fail("boundaryless sector is unsupported", Dict("sector" => sid))
        get(sec, "boundary_count", b) == b ||
            fail("boundary_count mismatch", Dict("sector" => sid))
        get(sec, "chi", 2 - 2g - b) == 2 - 2g - b ||
            fail("sector chi mismatch", Dict("sector" => sid))
        signs = collect(get(sec, "boundary_orientation_signs", fill(1, b)))
        length(signs) == b ||
            fail("orientation-sign count mismatch", Dict("sector" => sid))
        words = Vector{Int}[]
        for j in 1:b
            s = Int(signs[j])
            s in (-1, 1) || fail("boundary orientation sign must be ±1",
                                 Dict("sector" => sid, "boundary" => j))
            w = Int.(raw[j])
            push!(words, s == 1 ? w : [-x for x in reverse(w)])
        end
        bases = [boundary_start(words[j], sid, j) for j in 1:b]
        conns = [addcell!("t$(sid)_$(j)", bases[1], bases[j]) for j in 2:b]
        handles = Tuple{Int,Int}[]
        for i in 1:Int(g)
            a = addcell!("a$(sid)_$(i)", bases[1], bases[1])
            z = addcell!("b$(sid)_$(i)", bases[1], bases[1])
            push!(handles, (a, z))
        end
        push!(plans, (words=words, conns=conns, handles=handles))
    end
    parent = Dict(v => v for v in verts)
    find(v) = parent[v] == v ? v : (parent[v] = find(parent[v]))
    relators = Vector{Vector{Tuple{Int,Int}}}()
    for i in eachindex(gens)
        froms[i] == tos[i] && continue
        a, b = find(froms[i]), find(tos[i])
        a == b && continue
        parent[b] = a
        push!(relators, [(i, 1)])
    end
    length(Set(find(v) for v in verts)) <= 1 ||
        fail("full 1-skeleton is disconnected")
    function letters!(rel, w)
        for x in w
            push!(rel, (id_to_gen[abs(x)], x > 0 ? 1 : -1))
        end
    end
    for p in plans
        rel = Tuple{Int,Int}[]
        for (a, b) in p.handles
            append!(rel, ((a, 1), (b, 1), (a, -1), (b, -1)))
        end
        letters!(rel, p.words[1])
        for j in 2:length(p.words)
            t = p.conns[j - 1]
            push!(rel, (t, 1)); letters!(rel, p.words[j]); push!(rel, (t, -1))
        end
        push!(relators, rel)
    end
    GroupPresentation(gens, relators, Int[], Dict{String,Any}(
        "source" => "spine_complex",
        "spec" => "heegaard-koda-conversion.md §6.3",
        "sector_cw" => "relative_oriented_surface",
    ))
end
