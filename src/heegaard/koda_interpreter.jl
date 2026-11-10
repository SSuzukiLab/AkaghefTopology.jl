# Source: .research-scratch/PJ260922/ld-koda-spec/261109-replay-local-tables.py;
# docs/specs/heegaard-koda-conversion.md §8, §13.1–13.3 (finite table interpreter).
"""Decode o-data words through local JSON tables; build Heegaard diagrams from BW boundary."""

function _koda_tables_pack(tables=nothing)
    tables === nothing && (tables = load_koda_local_tables())
    V = tables.vertex["models"]
    EG = tables.edge["models"]
    BW = tables.bw_boundary
    C = tables.c_selection
    CHART = tables.vertex["chart_port_vectors"]
    (; V, EG, BW, C, CHART)
end

function _koda_port(x, out::Bool)
    name = (x > 0 ? "+" : "-") * (out ? "o" : "i")
    return name
end

function _koda_pushgroup!(groups::Dict, key, val)
    push!(get!(groups, key, Any[]), val)
end

function koda_table_decode(words::Vector{Vector{Int}}, sigma::Dict{Int, Int}; tables=nothing)
    pack = _koda_tables_pack(tables)
    V, EG, BW, C = pack.V, pack.EG, pack.BW, pack.C
    models = Dict(v => V[string(s)] for (v, s) in sigma)
    G = Dict{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}()
    edge_at = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    edges = Vector{Dict{String, Any}}()
    lift_edges = Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int}, Tuple{Int, Int}}}()
    port = (x, out) -> begin
        name = _koda_port(x, out)
        findfirst(==(name), models[abs(x)]["pd_ports_ccw"])
    end
    for word in words
        ext = vcat(word, word[1:1])
        n = length(word)
        for i in 1:n
            x, y = ext[i], ext[i + 1]
            v, p = abs(x), port(x, true)
            w, t = abs(y), port(y, false)
            key = string(sigma[v], ":", models[v]["ports"][string(p)]["name"], "->", sigma[w], ":", models[w]["ports"][string(t)]["name"])
            model = EG[key]
            iedge = length(edges)
            push!(edges, Dict("id" => iedge, "source" => [v, p], "target" => [w, t], "model" => key))
            edge_at[(v, p)] = (iedge, 0)
            edge_at[(w, t)] = (iedge, 1)
            for (ak, bv) in model["sheet_map"]
                qa = [parse(Int, string(c)) for c in string(ak)]
                qb = [parse(Int, string(c)) for c in string(bv)]
                q = only(z for z in qa if z != p)
                r = only(z for z in qb if z != t)
                G[(v, p, q)] = (w, t, r)
                G[(w, t, r)] = (v, p, q)
            end
            for (rank, (r, s)) in enumerate(model["boundary_lift_map"])
                push!(lift_edges, ((iedge, rank), (v, r), (w, s)))
            end
        end
    end
    sectors = _KodaDSU()
    for a in keys(G)
        _koda_dsu_union!(sectors, a, G[a])
        _koda_dsu_union!(sectors, a, (a[1], a[3], a[2]))
    end
    groups = Dict{Any, Vector{Any}}()
    for a in keys(G)
        _koda_pushgroup!(groups, _koda_dsu_find!(sectors, a), a)
    end
    cycles = Vector{Vector{Tuple{Int, Int, Int}}}()
    sector_flags = Vector{Vector{Tuple{Int, Int, Int}}}()
    group_sizes = [length(v) for v in values(groups)]
    for flags in values(groups)
        start = minimum(flags)
        x = start
        cycle = Tuple{Int, Int, Int}[]
        while isempty(cycle) || x != start
            push!(cycle, x)
            x = G[(x[1], x[3], x[2])]
        end
        length(cycle) * 2 == length(flags) ||
            throw(KodaImplementationError(:koda_table_decode, "sector cycle length mismatch", Dict("cycle" => length(cycle), "flags" => length(flags), "start" => start, "G" => length(G), "group_sizes" => group_sizes)))
        push!(cycles, cycle)
        push!(sector_flags, sort(flags))
    end
    conn = _KodaDSU()
    for (_, a, b) in lift_edges
        _koda_dsu_union!(conn, a, b)
    end
    faces = Vector{Dict{String, Any}}()
    edge_faces = Dict{Any, Vector{Any}}()
    for (sector, walk) in enumerate(cycles)
        for color in BW["sector_interior"]["lifts"]
            vs = Tuple{Int, Int}[]
            es = Tuple{Int, Int}[]
            for (v, p, q) in walk
                ps = sort([p, q])
                sheet = string(ps[1], ps[2])
                r = models[v]["sheets"][sheet]["boundary_lifts"][color]
                push!(vs, (v, r))
                eid = edge_at[(v, q)][1]
                rank = findfirst(==(r), models[v]["ports"][string(q)]["face_vertices"]) - 1
                push!(es, (eid, rank))
            end
            fid = length(faces) + 1
            for e in es
                _koda_pushgroup!(edge_faces, e, fid)
            end
            push!(faces, Dict(
                "id" => fid,
                "sector" => sector,
                "color" => color,
                "vertices" => vs,
                "edges" => es,
            ))
        end
    end
    all(length(fs) == 2 for fs in values(edge_faces)) ||
        throw(KodaImplementationError(:koda_table_decode, "edge face incidence not 2", Dict()))
    regions = _KodaDSU()
    for f in faces
        _koda_dsu_find!(regions, f["id"])
    end
    for ((_, rank), fs) in edge_faces
        colors = Set(faces[f]["color"] for f in fs)  # f is 1-based face index
        if rank == 1
            colors == Set(["W", "B"]) ||
                throw(KodaImplementationError(:koda_table_decode, "rank-1 edge not W/B", Dict()))
        else
            length(colors) == 1 || throw(KodaImplementationError(:koda_table_decode, "rank≠1 edge mixed colors", Dict()))
            for f in fs[2:end]
                _koda_dsu_union!(regions, fs[1], f)
            end
        end
    end
    function collect_regions(grouping, closed::Bool)
        out = Dict{String, Any}[]
        for (key, inds) in sort(collect(grouping); by=first)
            vertices = Set(v for i in inds for (v, _) in faces[i]["vertices"])
            es = Set(e for i in inds for e in faces[i]["edges"])
            chi = length(vertices) - length(es) + length(inds)
            item = Dict(
                "id" => length(out),
                "root" => key isa Tuple ? collect(key) : key,
                "V" => length(vertices),
                "E" => length(es),
                "F" => length(inds),
                "chi" => chi,
                "colors" => sort(unique(faces[i]["color"] for i in inds)),
                "face_ids" => inds,
                "vertices" => sort(collect(vertices)),
                "edges" => sort(collect(es)),
            )
            closed && (item["genus"] = (2 - chi) ÷ 2)
            push!(out, item)
        end
        out
    end
    bgroups = Dict{Any, Vector{Int}}()
    rgroups = Dict{Any, Vector{Int}}()
    for f in faces
        _koda_pushgroup!(bgroups, _koda_dsu_find!(conn, f["vertices"][1]), f["id"])
        _koda_pushgroup!(rgroups, _koda_dsu_find!(regions, f["id"]), f["id"])
    end
    boundary = collect_regions(bgroups, true)
    bw = collect_regions(rgroups, false)
    rho = Dict{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}()
    tau = Dict{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}()
    for (v, model) in models
        for (r, row) in model["boundary_vertices"]
            ports = row["ports_ccw"]
            ports_ext = vcat(ports, ports[1:1])
            ri = parse(Int, string(r))
            for i in 1:length(ports)
                a, b = ports_ext[i], ports_ext[i + 1]
                rho[(v, a, ri)] = (v, b, ri)
            end
        end
    end
    for edge in edges
        v, p = edge["source"][1], edge["source"][2]
        w, t = edge["target"][1], edge["target"][2]
        model = EG[edge["model"]]
        for (r, s) in model["boundary_lift_map"]
            tau[(v, p, r)] = (w, t, s)
            tau[(w, t, s)] = (v, p, r)
        end
    end
    ribbon = _koda_orbit(Dict(h => rho[tau[h]] for h in keys(rho)))
    length(ribbon) == length(faces) ||
        throw(KodaImplementationError(:koda_table_decode, "ribbon face count mismatch", Dict()))
    enext = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    eout = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()
    for (v, model) in models
        for (r, row) in model["boundary_vertices"]
            row["stratum"] == "e" || continue
            p = only([p for p in row["e_ports"] if endswith(models[v]["ports"][string(p)]["name"], "o")])
            eid, tend = edge_at[(v, p)]
            edge = edges[eid + 1]
            w, t = edge["target"][1], edge["target"][2]
            tend == 0 || throw(KodaImplementationError(:koda_table_decode, "unexpected edge end", Dict()))
            rr = models[v]["ports"][string(p)]["lifts"]["e"]
            emodel = EG[edge["model"]]
            blmap = Dict{Any, Any}(r => s for (r, s) in emodel["boundary_lift_map"])
            ss = blmap[rr]
            enext[(v, parse(Int, string(r)))] = (w, ss)
            eout[(v, parse(Int, string(r)))] = (eid, 1)
        end
    end
    ecycles = _koda_orbit(enext)
    ekey(point) = begin
        v, r = point
        row = models[v]["boundary_vertices"][string(r)]
        p = only([p for p in row["e_ports"] if endswith(models[v]["ports"][string(p)]["name"], "o")])
        (v, models[v]["ports"][string(p)]["name"], edge_at[(v, p)][1])
    end
    ecycles = [
        begin
            idx = argmin(i -> ekey(cycle[i]), 1:length(cycle))
            vcat(cycle[idx:end], cycle[1:(idx - 1)])
        end for cycle in ecycles
    ]
    sort!(ecycles, by=c -> ekey(c[1]))
    for (i, cycle) in enumerate(ecycles)
        edge = eout[cycle[1]]
        wf = [f for f in get(edge_faces, edge, Int[]) if faces[f]["color"] == "W"]
        isempty(wf) && continue
        for region in bw
            if wf[1] in region["face_ids"]
                haskey(region, "e_components") || (region["e_components"] = Int[])
                push!(region["e_components"], i - 1)
                break
            end
        end
    end
    for region in bw
        b = length(get(region, "e_components", Int[]))
        region["boundary_count"] = b
        region["genus"] = (2 - b - region["chi"]) ÷ 2
    end
    Dict(
        "input_vertices" => length(sigma),
        "input_edges" => length(edges),
        "sector_count" => length(cycles),
        "spine_chi" => length(sigma) - length(edges) + length(cycles),
        "sector_cycles" => cycles,
        "sector_flags" => sector_flags,
        "edges" => edges,
        "boundary" => boundary,
        "BW" => bw,
        "e_cycles" => ecycles,
        "boundary_ribbon_face_lengths" => sort([length(r) for r in ribbon]),
    )
end

function _koda_rational(t)
    Rational(t[1], t[2])
end

function _koda_polyline_crossings(polyline, source::Rational, target::Rational)
    values = Dict("source" => source, "target" => target, "opposite-source" => -source)
    points = [(_koda_rational(t), values[string(value)]) for (t, value) in polyline]
    all(points[i][1] < points[i + 1][1] for i in 1:(length(points) - 1)) ||
        throw(KodaImplementationError(:koda_hc_from_boundary, "polyline not monotone", Dict()))
    out = Tuple{String, String, Vector{Int}}[]
    for i in 1:(length(points) - 1)
        (t, y) = points[i]
        (u, z) = points[i + 1]
        y != 0 && z != 0 || throw(KodaImplementationError(:koda_hc_from_boundary, "polyline hits zero", Dict()))
        if y * z < 0
            crossing = t + (u - t) * (-y) // (z - y)
            push!(out, (y > 0 ? "B" : "W", z > 0 ? "B" : "W", [numerator(crossing), denominator(crossing)]))
        end
    end
    out
end

function koda_hc_from_boundary(spine::Dict{String, Any}, sigma::Dict{Int, Int}; tables=nothing)
    pack = _koda_tables_pack(tables)
    V, C, CHART = pack.V, pack.C, pack.CHART
    rule = Dict(r["spur_side"] => r for r in C["vertex_rules"])
    events = Vector{Vector{Dict{String, Any}}}()
    corners = Dict{String, Dict{String, Any}}()
    choices = Dict{String, Any}[]
    serial = 0
    positions = Dict{Int, Tuple{Int, Int}}(parse(Int, string(k)) => (Int(v[1]), Int(v[2])) for (k, v) in CHART)
    for (ci, cycle) in enumerate(spine["e_cycles"])
        states = [rule[V[string(sigma[v])]["boundary_vertices"][string(r)]["spur_side"]]["C_side"] for (v, r) in cycle]
        homogeneous = length(unique(states)) == 1
        seq = Dict{String, Any}[]
        coordinates = Dict{String, Any}[]
        for (i, ((v, r), side)) in enumerate(zip(cycle, states))
            push!(seq, Dict("type" => "true", "id" => string(v), "side" => side, "lift" => [v, r]))
            target = states[mod1(i + 1, length(states))]
            spur = V[string(sigma[v])]["boundary_vertices"][string(r)]["spur_side"]
            off = rule[spur]["offset"]
            source_height = Rational(Int(off[1]), Int(off[2]))
            (tv, tr) = cycle[mod1(i + 1, length(cycle))]
            spur2 = V[string(sigma[tv])]["boundary_vertices"][string(tr)]["spur_side"]
            off2 = rule[spur2]["offset"]
            target_height = Rational(Int(off2[1]), Int(off2[2]))
            geometry = if homogeneous && i == 1
                C["edge_rules"]["no_crossing_component"]["polyline"]
            else
                C["edge_rules"][side != target ? "opposite_sides" : "same_side"]
            end
            changes = _koda_polyline_crossings(geometry, source_height, target_height)
            push!(coordinates, Dict(
                "vertex" => [v, r],
                "offset" => [numerator(source_height), denominator(source_height)],
                "edge_after" => i - 1,
                "crossings" => [c[3] for c in changes],
            ))
            for (before, after, param) in changes
                serial += 1
                name = "j$serial"
                corners[name] = Dict("before" => before, "after" => after, "e_component" => ci - 1, "edge_after" => i - 1, "param" => param)
                push!(seq, Dict("type" => "corner", "id" => name, "before" => before, "after" => after))
            end
        end
        any(x -> x["type"] == "corner", seq) || koda_not_computed(
            :spine_to_heegaard_koda,
            :vertex_free_odata_not_encoded,
            Dict{String, Any}("e_component" => ci - 1, "spec" => "§8.2 / §13.3"),
        )
        push!(events, seq)
        push!(choices, Dict("e_component" => ci - 1, "collar_vertices" => coordinates, "added_pair_for_no_switch" => homogeneous))
    end
    alpha = Vector{Vector{String}}()
    beta = Vector{Vector{String}}()
    epsilon = Dict{String, Int}()
    local_rotations = Dict{String, Dict{String, Any}}()
    for (v, s) in sigma
        model = V[string(s)]
        out = Dict{String, Tuple{Int, Int}}()
        for (r, row) in model["boundary_vertices"]
            row["stratum"] == "e" || continue
            side = rule[row["spur_side"]]["C_side"]
            family = C["disk_types"]["C_in_" * side]["family"]
            p = only([p for p in row["e_ports"] if endswith(model["ports"][string(p)]["name"], "o")])
            out[family] = positions[p]
        end
        a, b = out["alpha"], out["beta"]
        epsilon[string(v)] = a[1] * b[2] - a[2] * b[1]
        abs(epsilon[string(v)]) == 1 || throw(KodaImplementationError(:koda_hc_from_boundary, "bad epsilon", Dict("vertex" => v)))
    end
    jcycle = [(part, direction) for (part, direction) in C["joining_corner"]["ports_ccw"]]
    for (name, row) in corners
        rays = Dict{String, Tuple{Int, Int}}()
        for (i, (part, direction)) in enumerate(jcycle)
            side = row[part]
            family = C["disk_types"]["C_in_" * side]["family"]
            if endswith(direction, "out")
                rays[family] = positions[i]
            end
        end
        a, b = rays["alpha"], rays["beta"]
        epsilon[name] = a[1] * b[2] - a[2] * b[1]
        local_rotations[name] = Dict("ports_ccw" => jcycle, "epsilon" => epsilon[name])
    end
    for seq in events
        n = length(seq)
        for (i, item) in enumerate(seq)
            item["type"] == "corner" || continue
            side = item["after"]
            ids = [item["id"]]
            j = mod1(i + 1, n)
            while true
                push!(ids, seq[j]["id"])
                seq[j]["type"] == "corner" && break
                seq[j]["side"] == side || throw(KodaImplementationError(:koda_hc_from_boundary, "incident side mismatch", Dict()))
                j = mod1(j + 1, n)
            end
            (side == "B" ? alpha : beta) |> (curves -> push!(curves, ids))
        end
    end
    length(alpha) == length(beta) == length(corners) ÷ 2 ||
        throw(KodaImplementationError(:koda_hc_from_boundary, "alpha/beta curve count", Dict()))
    _, _, _ = _koda_ribbon(alpha, beta, epsilon)
    js = Int[]
    for (ci, seq) in enumerate(events)
        boundary_darts = Tuple{String, Char, Int}[]
        for (family, curves) in (('a', alpha), ('b', beta))
            for curve in curves
                corners[curve[1]]["e_component"] == ci - 1 || continue
                push!(boundary_darts, (curve[1], family, -1))
            end
        end
        faces, _, _ = _koda_ribbon(alpha, beta, epsilon)
        bd = Set(boundary_darts)
        candidates = [i for (i, f) in enumerate(faces) if Set(f) == bd]
        length(candidates) == 1 || throw(KodaImplementationError(:koda_hc_from_boundary, "joining face ambiguous", Dict("component" => ci - 1)))
        push!(js, candidates[1] - 1)
    end
    joining = Dict{String, Any}[]
    for region in spine["BW"]
        colors = region["colors"]
        ecomp = get(region, "e_components", nothing)
        if colors == ["W"] && ecomp !== nothing && !isempty(ecomp)
            push!(joining, Dict(
                "genus" => 0,
                "boundary_count" => length(ecomp),
                "e_components" => ecomp,
                "boundary_face_orbits" => [js[i + 1] for i in ecomp],
            ))
        end
    end
    chi = -length(epsilon) + length(_koda_ribbon(alpha, beta, epsilon)[1]) + 2 * (length(joining) - length(js))
    predicted = spine["spine_chi"] - length(alpha) - length(beta) + 2 * length(joining) - length(js)
    chi == predicted || throw(KodaImplementationError(:koda_hc_from_boundary, "chi mismatch", Dict("chi" => chi, "predicted" => predicted)))
    (2 - chi) % 2 == 0 || throw(KodaImplementationError(:koda_hc_from_boundary, "odd chi", Dict()))
    Dict(
        "genus" => (2 - chi) ÷ 2,
        "chi" => chi,
        "alpha" => alpha,
        "beta" => beta,
        "epsilon" => epsilon,
        "intersections" => length(epsilon),
        "true_intersections" => length(sigma),
        "joining_intersections" => length(corners),
        "joining" => joining,
        "J_vertices" => collect(keys(corners)),
        "C_choices" => choices,
        "corner_local_models" => local_rotations,
        "face_lengths" => sort([length(f) for f in _koda_ribbon(alpha, beta, epsilon)[1]]),
    )
end

function _koda_key_order(k::AbstractString)
    startswith(k, "j") ? (1, parse(Int, k[2:end])) : (0, parse(Int, k))
end

function koda_heegaard_from_hc(hc::Dict{String, Any})
    ε = hc["epsilon"]
    keys_sorted = sort(collect(keys(ε)); by=_koda_key_order)
    idmap = Dict{String, Int}(k => i for (i, k) in enumerate(keys_sorted))
    αraw = hc["alpha"]
    βraw = hc["beta"]
    α_curves = [Curve(:alpha, i, Int[idmap[s] for s in w]; label = "α$i") for (i, w) in enumerate(αraw)]
    β_curves = [Curve(:beta, i, Int[idmap[s] for s in w]; label = "β$i") for (i, w) in enumerate(βraw)]
    inter = Intersection[]
    for k in keys_sorted
        ai = findfirst(w -> k in w, αraw)
        bi = findfirst(w -> k in w, βraw)
        (ai === nothing || bi === nothing) &&
            throw(KodaImplementationError(:spine_to_heegaard_koda, "intersection missing from curve", Dict("id" => k)))
        push!(inter, Intersection(idmap[k], ai, bi, Int(ε[k])))
    end
    certify_heegaard(decorate(heegaard_code(α_curves, β_curves, inter); cellular = true))
end

function _koda_select_joining_region(h::HeegaardDiagram, hc::Dict{String, Any})
    hints = Int[]
    for block in hc["joining"]
        for x in get(block, "boundary_face_orbits", [])
            push!(hints, Int(x))
        end
    end
    cands = _ribbon_joining_indices(h)
    scored = Tuple{Int, Int}[]
    for r in cands
        joining_disk(h, r)
        score = count(h0 -> (h0 + 1) == r, hints)
        push!(scored, (score, r))
    end
    isempty(scored) && throw(KodaImplementationError(:spine_to_heegaard_koda, "no joining disk candidate", Dict("hints" => hints)))
    sort!(scored; rev=true)
    best = scored[1][1]
    tied = [r for (s, r) in scored if s == best]
    return best > 0 ? only(tied) : minimum(tied)
end

function koda_forward_gauss_from_hc(hc::Dict{String, Any}, sigma::Dict{Int, Int})
    # Curve words in the reverse HC are parameter words, not strand words.
    # §7.3: strand = dγ × parameter word, and σ = −ε dα dβ. Do not assume dγ = +1.
    h = koda_heegaard_from_hc(hc)
    region = _koda_select_joining_region(h, hc)
    dγ = _koda_curve_directions(h, joining_disk(h, region))
    function oriented(words, family)
        out = Vector{Vector{String}}()
        for (i, w) in enumerate(words)
            key = (family, i)
            haskey(dγ, key) || throw(KodaInputError(
                :heegaard_to_spine,
                "missing cγ; §5.2 requires every curve to meet J in one arc",
                Dict("key" => key),
            ))
            push!(out, dγ[key] == 1 ? collect(w) : reverse(collect(w)))
        end
        out
    end
    α = oriented(hc["alpha"], :alpha)
    β = oriented(hc["beta"], :beta)
    back, labels = _koda_smooth_joining(α, β, hc["J_vertices"], :beta)
    lbl_rev = Dict(v => k for (k, v) in labels)
    signs = Dict{Int, Int}()
    for v in keys(sigma)
        lbl = get(lbl_rev, string(v), nothing)
        lbl === nothing && continue
        ai = findfirst(w -> string(v) in w, hc["alpha"])
        bi = findfirst(w -> string(v) in w, hc["beta"])
        (ai === nothing || bi === nothing) && throw(KodaImplementationError(
            :heegaard_to_spine, "true vertex missing from a curve", Dict("id" => v),
        ))
        signs[lbl] = -Int(hc["epsilon"][string(v)]) * dγ[(:alpha, ai)] * dγ[(:beta, bi)]
    end
    (back, signs)
end
