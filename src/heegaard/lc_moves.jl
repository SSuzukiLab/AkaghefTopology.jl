# cf: lib/topology.jl/docs/specs/heegaard-moves.md §3–§8 — catalog rows and one apply path
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3–§4 — words, passage sign, place
# cf: lib/topology.jl/src/moves.jl — _cyclic_occurrences, _location_products, _vertex_maps, _map_fragment, _simple_pattern_rewrite, _reconnect_pattern

# Parent of Heegaard. Unexported Topology move kernels are called through this binding.
const _LC_MV = parentmodule(@__MODULE__)

struct LCMoveError <: Exception
    kind::Symbol
    message::String
    detail::Dict{String,Any}
end

function Base.showerror(io::IO, e::LCMoveError)
    print(io, "LCMoveError(", e.kind, "): ", e.message)
    isempty(e.detail) || print(io, " ", e.detail)
end

struct LCMoveWitness
    row_id::String
    inverse_id::String
    digest::String
    relation::String
    layer::String
    names::Vector{String}
    locations::Vector{Tuple{Int,Int}}
    vertices::Vector{Int}
    before_frags::Vector{Vector{Int}}
    after_frags::Vector{Vector{Int}}
    capacity::Vector{Int}
    fingerprint::String
    label_policy::String
    genus_delta::Int
    extra::Dict{String,Any}
end

struct LCMoveResult
    output
    inverse::LCMoveWitness
    correspondence::Dict{String,Any}
end

const _LC_CATALOG = Ref{Union{Nothing,Dict{String,Any}}}(nothing)

_plain(x::AbstractDict) = Dict{String,Any}(String(k) => _plain(v) for (k, v) in pairs(x))
_plain(x::AbstractVector) = Any[_plain(v) for v in x]
_plain(x::Integer) = x isa Bool ? x : Int(x)
_plain(x) = x

function _moves_dir()
    normpath(joinpath(@__DIR__, "..", "..", "data", "heegaard", "moves"))
end

function _load_json(path)
    _plain(JSON3.read(read(path, String)))
end

function _row_digest(row)
    bytes2hex(sha256(JSON3.write(row)))
end

function load_heegaard_moves(; force::Bool=false)
    _LC_CATALOG[] !== nothing && !force && return _LC_CATALOG[]
    root = _moves_dir()
    catalog = _load_json(joinpath(root, "catalog.json"))
    rows = Dict{String,Any}()
    for rel in catalog["files"]
        doc = _load_json(joinpath(root, rel))
        for row in doc["rows"]
            row["digest"] = _row_digest(row)
            row["file"] = rel
            haskey(rows, row["id"]) && throw(LCMoveError(:catalog, "duplicate row id", Dict("id" => row["id"])))
            rows[row["id"]] = row
        end
    end
    loaded = Dict{String,Any}("schema_version" => catalog["schema_version"], "rows" => rows, "dir" => root,
        "shared" => _load_json(joinpath(root, "shared.json")))
    _LC_CATALOG[] = loaded
    loaded
end

function heegaard_move_row(id::AbstractString)
    rows = load_heegaard_moves()["rows"]
    haskey(rows, id) || throw(LCMoveError(:catalog, "unknown row", Dict("id" => id)))
    rows[id]
end

heegaard_move_ids() = sort(collect(keys(load_heegaard_moves()["rows"])))

include("lc_moves_support.jl")

_k_limit(g::Integer) = g <= 0 ? 0 : 3 * g - 2
_lc_of(lc::LickorishConfiguration) = lc
_lc_of(hc::HeegaardLayer) = hc.lc
_roles_of(::LickorishConfiguration) = Dict{String,Symbol}()
_roles_of(hc::HeegaardLayer) = hc.roles
_marks(lc::LickorishConfiguration) = (disks=lc.disks, points=lc.points, weights=lc.weights)
_marks(hc::HeegaardLayer) = (disks=hc.disks, points=hc.points, weights=hc.weights)

function _component_names(lc::LickorishConfiguration)
    k = lc.genus > 0 ? _lc_k_names(lc.genus) : String[]
    rest = sort([n for n in keys(lc.words) if !(n in k)])
    String[vcat(k, rest)...]
end

function _fingerprint(lc::LickorishConfiguration, roles=Dict{String,Symbol}())
    parts = String["g=$(lc.genus)"]
    for name in _component_names(lc)
        push!(parts, name * ":" * join(get(lc.words, name, Int[]), ","))
        haskey(lc.places, name) && push!(parts, "place:" * name * "=" * lc.places[name])
    end
    for name in sort(collect(keys(roles)))
        push!(parts, "role:" * name * "=" * string(roles[name]))
    end
    bytes2hex(sha256(join(parts, "|")))
end

function _missing_marks(input)
    marks = _marks(input)
    missing = String[]
    isempty(marks.disks) || push!(missing, "disk")
    isempty(marks.points) || push!(missing, "point")
    isempty(marks.weights) || push!(missing, "weight")
    missing
end

function _refuse_marks(input, row)
    get(row, "transport", "") == "refuse_if_marks" || return
    missing = _missing_marks(input)
    isempty(missing) || throw(LCMoveError(
        :missing_transport,
        "印の輸送指定が足りない",
        Dict("missing" => missing, "row" => row["id"]),
    ))
end

_frag_word(fragment) = Int[Int(x) for x in fragment["word"]]
_frags(side) = side["fragments"]
_capacity_of(side, n) = haskey(side, "capacity") && !isempty(side["capacity"]) ? Int[Int(x) for x in side["capacity"]] : collect(1:n)

function _binding_roles(row, side)
    bindings = get(row, "bindings", Dict{String,Any}())
    roles = String[]
    kinds = Dict{String,String}()
    for fragment in _frags(side)
        role = String(fragment["role"])
        push!(roles, role)
        kind = get(bindings, role, nothing)
        kind isa AbstractString || throw(LCMoveError(:bindings, "role has no binding", Dict("role" => role, "row" => row["id"])))
        kinds[role] = String(kind)
    end
    roles, kinds
end

function _is_k(lc, name)
    name in _lc_k_names(lc.genus)
end

function _candidates_for(lc, kind)
    names = _component_names(lc)
    if kind == "k"
        return [n for n in names if _is_k(lc, n)]
    elseif kind == "input"
        return [n for n in names if !_is_k(lc, n)]
    elseif kind == "any"
        return [n for n in names if !isempty(get(lc.words, n, Int[]))]
    end
    throw(LCMoveError(:bindings, "unknown binding kind", Dict("kind" => kind)))
end

function _assignments(lc, roles, kinds)
    pools = [_candidates_for(lc, kinds[role]) for role in roles]
    found = Vector{Dict{String,String}}()
    function visit(i, used, acc)
        if i > length(roles)
            push!(found, copy(acc))
            return
        end
        for name in pools[i]
            name in used && continue
            push!(used, name)
            acc[roles[i]] = name
            visit(i + 1, used, acc)
            delete!(acc, roles[i])
            pop!(used)
        end
    end
    visit(1, Set{String}(), Dict{String,String}())
    found
end

function _window(word, start, len)
    n = length(word)
    [word[mod1(start + offset, n)] for offset in 0:(len - 1)]
end

function _consistent_windows(word, pattern, subst)
    hits = Tuple{Int,Dict{Int,Int}}[]
    isempty(word) && return hits
    isempty(pattern) && return hits
    for start in 1:length(word)
        window = _window(word, start, length(pattern))
        next = copy(subst)
        ok = true
        for (got, want) in zip(window, pattern)
            sign(got) == sign(want) || (ok = false; break)
            local_label = abs(want)
            global_label = abs(got)
            if haskey(next, local_label) && next[local_label] != global_label
                ok = false
                break
            end
            next[local_label] = global_label
        end
        ok && push!(hits, (start, next))
    end
    hits
end

function _apply_subst(pattern, subst)
    out = Int[]
    for entry in pattern
        haskey(subst, abs(entry)) || return nothing
        push!(out, sign(entry) * subst[abs(entry)])
    end
    out
end

function _face_sets(words)
    active = Dict(k => v for (k, v) in words if !isempty(v))
    faces, _ = _lc_faces(active)
    [(id=face.id, n=length(face.corners), labels=Set(c[1] for c in face.corners)) for face in faces]
end

function _applicability_ok(lc, row, words, labels)
    checks = String[string(x) for x in get(row, "applicability", Any[])]
    faces = _face_sets(words)
    limit = _k_limit(lc.genus)
    label_set = Set(values(labels))
    for check in checks
        if check == "digon_face"
            any(face.n == 2 && face.labels == label_set for face in faces) || return false
        elseif check == "empty_triangle"
            any(face.n == 3 && face.labels == label_set for face in faces) || return false
        elseif check == "labels_above_k"
            all(label > limit for label in values(labels)) || return false
        elseif check == "pin_is_k"
            pin = get(get(row, "bindings", Dict()), "pin_local", Any[])
            for local_label in pin
                global_label = get(labels, Int(local_label), nothing)
                global_label === nothing && return false
                global_label <= limit || return false
            end
        elseif check in ("k_edge_fixed", "gaps_on_common_face", "same_family", "band_avoids_family",
                         "support_supplied", "q_integer", "disk_supplied", "genus_plus_one")
            continue
        else
            throw(LCMoveError(:applicability, "unknown applicability term", Dict("term" => check, "row" => row["id"])))
        end
    end
    true
end

function _witness_from(input, row, names, locations, vertices, before_frags, after_frags, capacity; extra=Dict{String,Any}())
    lc = _lc_of(input)
    LCMoveWitness(
        String(row["id"]),
        row["inverse_id"] === nothing ? "" : String(row["inverse_id"]),
        String(row["digest"]),
        String(row["relation"]),
        String(row["layer"]),
        names,
        locations,
        vertices,
        before_frags,
        after_frags,
        capacity,
        _fingerprint(lc, _roles_of(input)),
        String(get(row, "label_policy", "pin_k")),
        Int(get(row, "genus_delta", 0)),
        extra,
    )
end

function _search_literal(input, row)
    lc = _lc_of(input)
    op = String(get(get(row, "expansion", Dict()), "op", "literal"))
    op == "literal" || return LCMoveWitness[]
    String(get(row, "status", "")) == "unresolved" && throw(LCMoveError(
        :unresolved_row, "関係が未指定の行は適用しない", Dict("row" => row["id"])))
    side = row["source"]
    fragments = _frags(side)
    isempty(fragments) && return LCMoveWitness[]
    any(isempty(_frag_word(fragment)) for fragment in fragments) && return LCMoveWitness[]
    roles, kinds = _binding_roles(row, side)
    patterns = Dict(String(fragment["role"]) => _frag_word(fragment) for fragment in fragments)
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    words = Dict(name => get(lc.words, name, Int[]) for name in names)
    capacity = _capacity_of(side, length(fragments))
    after_frags = [_frag_word(fragment) for fragment in _frags(row["target"])]
    length(after_frags) == length(fragments) || throw(LCMoveError(
        :catalog, "source and target fragment counts differ", Dict("row" => row["id"])))
    witnesses = LCMoveWitness[]
    seen = Set{Any}()
    for assignment in _assignments(lc, roles, kinds)
        function visit(role_index, subst)
            if role_index > length(roles)
                substituted = Vector{Vector{Int}}()
                options = Vector{Vector{Tuple{Int,Int}}}()
                for role in roles
                    concrete = _apply_subst(patterns[role], subst)
                    concrete === nothing && return
                    push!(substituted, concrete)
                    hits = _LC_MV._cyclic_occurrences(words[assignment[role]], concrete)
                    isempty(hits) && return
                    push!(options, [(index[assignment[role]], hit) for hit in hits])
                end
                products = _LC_MV._location_products(options, substituted, [words[name] for name in names])
                for locations in products
                    key = (copy(assignment), locations, copy(subst))
                    key in seen && continue
                    push!(seen, key)
                    _applicability_ok(lc, row, words, subst) || continue
                    max_local = maximum(abs.(reduce(vcat, values(patterns); init=Int[])))
                    vertices = Int[]
                    for local_label in 1:max_local
                        haskey(subst, local_label) || throw(LCMoveError(
                            :catalog, "local label missing from substitution", Dict("label" => local_label, "row" => row["id"])))
                        push!(vertices, subst[local_label])
                    end
                    before_local = [patterns[role] for role in roles]
                    push!(witnesses, _witness_from(input, row, names, locations, vertices, before_local, after_frags, capacity))
                end
                return
            end
            role = roles[role_index]
            for (_, next) in _consistent_windows(words[assignment[role]], patterns[role], subst)
                visit(role_index + 1, next)
            end
        end
        visit(1, Dict{Int,Int}())
    end
    witnesses
end

function _selection_dict(; kwargs...)
    Dict{String,Any}(String(k) => v for (k, v) in kwargs)
end

function _require_selection(row, selections)
    required = String[string(x) for x in get(row, "requires", Any[])]
    missing = [key for key in required if !haskey(selections, key) || selections[key] === nothing]
    isempty(missing) || throw(LCMoveError(
        :missing_selection, "支持または選択が足りない", Dict("missing" => missing, "row" => row["id"])))
end

function _rot_left(word, k)
    isempty(word) && return Int[]
    r = mod(k, length(word))
    r == 0 ? copy(word) : vcat(word[r + 1:end], word[1:r])
end

function _rot_right(word, k)
    isempty(word) && return Int[]
    _rot_left(word, -k)
end

function _vmax(frags)
    labels = abs.(reduce(vcat, frags; init=Int[]))
    isempty(labels) ? 0 : maximum(labels)
end

function _kernel_apply(names, words, before_frags, after_frags, locations, vertices, capacity;
        label_map=nothing)
    work = [copy(get(words, name, Int[])) for name in names]
    locs = Tuple{Int,Int}[]
    rotations = zeros(Int, length(names))
    for (comp, start) in locations
        push!(locs, (comp, start))
    end
    for (i, (comp, start)) in enumerate(locs)
        frag = before_frags[i]
        isempty(frag) && continue
        n = length(work[comp])
        if start + length(frag) - 1 > n
            r = start - 1
            work[comp] = _rot_left(work[comp], r)
            locs[i] = (comp, 1)
            rotations[comp] = r
        end
    end
    before = _LC_MV.MovePattern("lc", before_frags, ones(Int, _vmax(before_frags)), capacity)
    after = _LC_MV.MovePattern("lc", after_frags, ones(Int, _vmax(after_frags)), capacity)
    max_label = maximum(vcat([0], abs.(reduce(vcat, work; init=Int[]))))
    orientation = ones(Int, max(max_label, 0))
    if label_map === nothing
        length(vertices) == before.vertex_count || throw(LCMoveError(
            :witness, "vertex vector does not match the source pattern",
            Dict("vertices" => length(vertices), "expected" => before.vertex_count)))
        old_to_new, after_vertices, new_orientation = _LC_MV._vertex_maps(orientation, vertices, before, after)
    else
        # Fragments are already global labels. The row supplies the surviving map.
        old_to_new = label_map
        after_vertices = Int[]
        new_orientation = Int[]
    end
    mapped_before = label_map === nothing ?
        [_LC_MV._map_fragment(frag, vertices) for frag in before_frags] : before_frags
    mapped_after = label_map === nothing ?
        [_LC_MV._map_fragment(frag, after_vertices) for frag in after_frags] : after_frags
    for (i, (comp, start)) in enumerate(locs)
        frag = mapped_before[i]
        isempty(frag) && continue
        start in _LC_MV._cyclic_occurrences(work[comp], frag) || throw(LCMoveError(
            :witness_mismatch, "witness fragment is not at the recorded location",
            Dict("component" => names[comp], "start" => start, "fragment" => frag)))
    end
    candidate = _LC_MV.MoveCandidate("LC", "lc", vertices, NTuple{2,Int}[], :forward, before, after, locs)
    weights = [zeros(Float64, length(word)) for word in work]
    same = length(before_frags) == length(after_frags) && before.capacity == after.capacity
    converted = same ? _LC_MV._simple_pattern_rewrite(
        work, weights, candidate, mapped_before, mapped_after, old_to_new, new_orientation; weighted=false) : nothing
    if converted === nothing
        label_map === nothing || throw(LCMoveError(
            :witness, "global fragment rewrite did not apply", Dict("components" => names)))
        converted = _LC_MV._reconnect_pattern(
            work, weights, candidate, mapped_before, mapped_after, old_to_new, new_orientation; weighted=false)
    end
    result = converted[1]
    # A wrap rotation only linearizes the fragment. Survivors are already in the
    # original cyclic order. Unrotating by the old shift modulo the new length
    # scrambles that order (p0254 digon.remove.11).
    result, old_to_new, after_vertices, mapped_before, mapped_after, locs, rotations
end

function _pin_ok(lc, old_to_new, removed)
    limit = _k_limit(lc.genus)
    for label in 1:limit
        label in removed && throw(LCMoveError(:label_policy, "K∩K label would be removed", Dict("label" => label)))
        haskey(old_to_new, label) && old_to_new[label] != label && throw(LCMoveError(
            :label_policy,
            "src/moves.jl::_vertex_maps densely renumbered a K∩K label; pin policy is not in that function",
            Dict("label" => label, "image" => old_to_new[label]),
        ))
    end
end

function _map_entry(entry, old_to_new)
    label = abs(entry)
    haskey(old_to_new, label) ? sign(entry) * old_to_new[label] : entry
end

function _fragment_neighbors(word, start, len)
    n = length(word)
    (n == 0 || len >= n) && return nothing
    prev = word[mod1(start - 1, n)]
    next = word[mod1(start + len, n)]
    (prev, next)
end

function _insert_before(word, prev_signed, next_signed)
    isempty(word) && return 1
    prev_signed == next_signed && return 1
    prev_at = findfirst(==(prev_signed), word)
    next_at = findfirst(==(next_signed), word)
    (prev_at === nothing || next_at === nothing) && throw(LCMoveError(
        :inverse, "digon neighbor is not in the output",
        Dict("prev" => prev_signed, "next" => next_signed, "word" => word)))
    mod1(prev_at + 1, length(word)) == next_at || throw(LCMoveError(
        :inverse, "digon neighbors are not adjacent after removal",
        Dict("prev" => prev_signed, "next" => next_signed, "word" => word)))
    next_at
end

function _adjacent_gap(word, prev_signed, next_signed)
    isempty(word) && return 1
    prev_signed == next_signed && return 1
    prev_at = findfirst(==(prev_signed), word)
    next_at = findfirst(==(next_signed), word)
    (prev_at === nothing || next_at === nothing) && return nothing
    mod1(prev_at + 1, length(word)) == next_at || return nothing
    next_at
end

# cf: lib/topology.jl/docs/specs/heegaard-moves.md §7 — 隣接して消した交点の再挿入空位は、その位置より前に残る文字の個数; reconstructed
function _output_gap_after_deletions(names, locations, before_frags, after_frags, input_words, comp, start)
    word = get(input_words, names[comp], Int[])
    n = length(word)
    removed = falses(n)
    for (i, (c, s)) in enumerate(locations)
        c == comp || continue
        isempty(after_frags[i]) || continue
        for k in 0:(length(before_frags[i]) - 1)
            p = s + k
            (1 <= p <= n) || return nothing
            removed[p] && return nothing
            removed[p] = true
        end
    end
    (1 <= start <= n + 1) || return nothing
    count(p -> !removed[p], 1:(start - 1)) + 1
end

function _reinsert_locations(names, locations, before_frags, after_frags, input_words, output_words, old_to_new)
    found = Tuple{Int,Int}[]
    for (i, (comp, start)) in enumerate(locations)
        frag = after_frags[i]
        if !isempty(frag)
            hits = _LC_MV._cyclic_occurrences(output_words[names[comp]], frag)
            isempty(hits) && throw(LCMoveError(
                :inverse, "after fragment is not in the output", Dict("component" => names[comp], "fragment" => frag)))
            length(hits) == 1 || throw(LCMoveError(
                :inverse, "after fragment is not unique; refusing to pick another site",
                Dict("component" => names[comp], "hits" => hits)))
            push!(found, (comp, hits[1]))
            continue
        end
        word = get(input_words, names[comp], Int[])
        neighbors = _fragment_neighbors(word, start, length(before_frags[i]))
        if neighbors === nothing
            push!(found, (comp, 1))
            continue
        end
        prev, next = neighbors
        out = get(output_words, names[comp], Int[])
        mapped_prev = _map_entry(prev, old_to_new)
        mapped_next = _map_entry(next, old_to_new)
        gap = _adjacent_gap(out, mapped_prev, mapped_next)
        if gap === nothing
            gap = _output_gap_after_deletions(
                names, locations, before_frags, after_frags, input_words, comp, start)
        end
        if gap === nothing
            gap = _insert_before(out, mapped_prev, mapped_next)
        end
        push!(found, (comp, gap))
    end
    found
end

function _edge_face_ids(words, corner_face, name, index)
    word = words[name]
    n = length(word)
    n < 2 && return String[]
    leaving = (abs(word[index]), word[index] > 0 ? "+" : "-", "out")
    arriving_at = mod1(index + 1, n)
    arriving = (abs(word[arriving_at]), word[arriving_at] > 0 ? "+" : "-", "in")
    ids = String[]
    haskey(corner_face, leaving) && push!(ids, corner_face[leaving])
    haskey(corner_face, arriving) && push!(ids, corner_face[arriving])
    unique(ids)
end

function _crossed_index(lc, names, locations)
    for i in eachindex(locations)
        _is_k(lc, names[locations[i][1]]) && return i
    end
    length(locations) >= 2 ? 2 : 1
end

function _merged_place(lc, names, locations, removed, old_to_new, after_words)
    isempty(removed) && return nothing
    k_index = _crossed_index(lc, names, locations)
    comp, start = locations[k_index]
    k_name = names[comp]
    haskey(lc.words, k_name) || return nothing
    active = Dict(k => v for (k, v) in lc.words if !isempty(v))
    faces, corner_face = _lc_faces(active)
    incident = _edge_face_ids(active, corner_face, k_name, start)
    digons = Set(face.id for face in faces if length(face.corners) == 2 && all(corner[1] in removed for corner in face.corners))
    dest_ids = [id for id in incident if !(id in digons)]
    if length(dest_ids) != 1
        # Both sides of the crossed edge are digons. The interior digon is the
        # one on the pushed curve's first visit; it merges into the non-digon
        # face that holds the opposite corner.
        push_index = only(i for i in eachindex(locations) if i != k_index)
        pcomp, pstart = locations[push_index]
        pname = names[pcomp]
        pword = lc.words[pname]
        entry = pword[pstart]
        sign = entry > 0 ? "+" : "-"
        interior = [get(corner_face, (abs(entry), sign, way), "") for way in ("in", "out")]
        interior = [id for id in interior if id in digons]
        length(interior) == 1 || throw(LCMoveError(
            :missing_place,
            "越えた辺の向こう側の面が一意でない",
            Dict("component" => k_name, "incident" => incident, "digons" => sort(collect(digons)), "interior" => interior)))
        dest_face = only(face for face in faces if face.id == interior[1])
        targets = String[]
        for (label, sg, way) in dest_face.corners
            other = way == "out" ? "in" : "out"
            id = get(corner_face, (label, sg, other), "")
            (id == "" || id in digons) && continue
            push!(targets, id)
        end
        dest_ids = unique(targets)
        length(dest_ids) == 1 || throw(LCMoveError(
            :missing_place,
            "二角形が合併する面が一意でない",
            Dict("component" => pname, "targets" => dest_ids, "interior" => interior)))
    end
    dest = only(face for face in faces if face.id == dest_ids[1])
    anchors = LCCorner[]
    for (label, sign, way) in dest.corners
        label in removed && continue
        haskey(old_to_new, label) || continue
        push!(anchors, (old_to_new[label], sign, way))
    end
    isempty(anchors) && throw(LCMoveError(
        :missing_place,
        "合併後の面を特定する残る角がない",
        Dict("component" => k_name, "dest" => dest.id)))
    after_active = Dict(k => v for (k, v) in after_words if !isempty(v))
    after_faces, _ = _lc_faces(after_active)
    hits = unique([face.id for face in after_faces if any(corner in face.corners for corner in anchors)])
    length(hits) == 1 || throw(LCMoveError(
        :missing_place,
        "除去の後に合併してできる面が一意でない",
        Dict("faces" => hits, "anchors" => anchors)))
    hits[1]
end

function _install_words(lc, names, result, places, assigned=Dict{String,String}())
    words = Dict{String,Vector{Int}}()
    new_names = String[]
    for (i, name) in enumerate(names)
        i > length(result) && break
        words[name] = result[i]
    end
    for i in (length(names) + 1):length(result)
        throw(LCMoveError(:component_name, "reconnect added a component with no name in the row", Dict("index" => i)))
    end
    for (name, word) in lc.words
        haskey(words, name) || (words[name] = copy(word))
    end
    new_places = copy(places)
    for (name, face) in assigned
        new_places[name] = face
    end
    for (name, word) in words
        isempty(word) || delete!(new_places, name)
    end
    active = Dict(k => v for (k, v) in words if !isempty(v))
    face_ids = String[]
    if !isempty(active)
        faces, _ = _lc_faces(active)
        face_ids = [face.id for face in faces]
    end
    for (name, word) in words
        isempty(word) || continue
        haskey(new_places, name) && continue
        length(face_ids) == 1 || throw(LCMoveError(
            :missing_place,
            "空語成分の place を一意に決められない",
            Dict("component" => name, "faces" => face_ids),
        ))
        new_places[name] = face_ids[1]
    end
    words, new_places, new_names
end

function _build_output(input, genus, words, places)
    lc0 = _lc_of(input)
    lc = LickorishConfiguration(genus, words, places, copy(lc0.disks), copy(lc0.points), copy(lc0.weights), lc0.source_path)
    validate_lc_input(lc)
    if any(!isempty(word) for word in values(words))
        _lc_validate_surface!(lc)
    end
    if input isa HeegaardLayer
        roles = copy(input.roles)
        return HeegaardLayer(lc, roles, copy(input.disks), copy(input.points), copy(input.weights), input.hc_path)
    end
    lc
end

function _inverse_locations(names, output_words, mapped_after, locations, rotations)
    found = Tuple{Int,Int}[]
    for (i, (comp, _)) in enumerate(locations)
        frag = mapped_after[i]
        if isempty(frag)
            word = output_words[names[comp]]
            # Reinsertion point recorded by the forward splice: the linearized start, unrotated.
            start = locations[i][2]
            r = rotations[comp]
            if r != 0 && !isempty(word)
                start = mod1(1 + r, length(word) + 1)
                start == length(word) + 1 && (start = length(word) + 1)
            elseif isempty(word)
                start = 1
            end
            push!(found, (comp, start))
        else
            hits = _LC_MV._cyclic_occurrences(output_words[names[comp]], frag)
            isempty(hits) && throw(LCMoveError(
                :inverse, "after fragment is not in the output", Dict("component" => names[comp], "fragment" => frag)))
            length(hits) == 1 || throw(LCMoveError(
                :inverse, "after fragment is not unique; refusing to pick another site",
                Dict("component" => names[comp], "hits" => hits)))
            push!(found, (comp, hits[1]))
        end
    end
    found
end

function _label_map_of(extra)
    raw = get(extra, "label_map", nothing)
    raw isa AbstractDict ? raw : nothing
end

function _copy_extra(extra)
    Dict{String,Any}(String(k) => v for (k, v) in extra)
end

function _attach_roles(output, witness::LCMoveWitness, drop_names)
    output isa HeegaardLayer || return output
    roles = copy(output.roles)
    for name in drop_names
        delete!(roles, name)
    end
    if witness.genus_delta > 0
        for (name, role) in get(witness.extra, "new_roles", Dict{String,Any}())
            roles[String(name)] = role isa Symbol ? role : Symbol(String(role))
        end
    end
    HeegaardLayer(output.lc, roles, output.disks, output.points, output.weights, output.hc_path)
end

function _apply_witness(input, witness::LCMoveWitness)
    if get(witness.extra, "identity", false)
        return LCMoveResult(input, witness, Dict{String,Any}("relation" => witness.relation, "identity" => true))
    end
    lc = _lc_of(input)
    _fingerprint(lc, _roles_of(input)) == witness.fingerprint || throw(LCMoveError(
        :witness_mismatch, "witness was not built for this input", Dict("row" => witness.row_id)))
    words = Dict(name => copy(get(lc.words, name, Int[])) for name in witness.names)
    result, old_to_new, after_vertices, mapped_before, mapped_after, locs, rotations = _kernel_apply(
        witness.names, words, witness.before_frags, witness.after_frags,
        witness.locations, witness.vertices, witness.capacity;
        label_map=_label_map_of(witness.extra))
    removed = setdiff(witness.vertices, keys(old_to_new))
    witness.label_policy == "pin_k" && _pin_ok(lc, old_to_new, removed)
    preview = Dict(witness.names[i] => result[i] for i in eachindex(witness.names) if i <= length(result))
    for (name, word) in lc.words
        haskey(preview, name) || (preview[name] = copy(word))
    end
    assigned = Dict{String,String}()
    created = String[name for name in get(witness.extra, "created_components", String[])]
    emptied = [name for (name, word) in preview if isempty(word) && !haskey(lc.places, name) && !isempty(get(lc.words, name, Int[])) && !(name in created)]
    if !isempty(emptied)
        face = _merged_place(lc, witness.names, witness.locations, Set(removed), old_to_new, preview)
        for name in emptied
            assigned[name] = face
        end
    end
    drop_names = [name for name in created if isempty(get(preview, name, Int[0]))]
    if !isempty(drop_names)
        new_words = Dict{String,Vector{Int}}()
        for (i, name) in enumerate(witness.names)
            i > length(result) && break
            name in drop_names && continue
            new_words[name] = result[i]
        end
        for (name, word) in lc.words
            name in drop_names && continue
            haskey(new_words, name) || (new_words[name] = copy(word))
        end
        new_places = copy(lc.places)
        for name in drop_names
            delete!(new_places, name)
        end
        for (name, face) in assigned
            haskey(new_words, name) && (new_places[name] = face)
        end
    else
        new_words, new_places, _ = _install_words(lc, witness.names, result, lc.places, assigned)
    end
    output = _attach_roles(
        _build_output(input, lc.genus + witness.genus_delta, new_words, new_places),
        witness, drop_names)
    out_lc = _lc_of(output)
    out_words = Dict(name => get(out_lc.words, name, Int[]) for name in witness.names)
    inverse_locations = any(isempty, mapped_after) ?
        _reinsert_locations(witness.names, witness.locations, witness.before_frags, mapped_after, words, out_words, old_to_new) :
        _inverse_locations(witness.names, out_words, mapped_after, witness.locations, rotations)
    inverse_vertices = isempty(witness.after_frags) || _vmax(witness.after_frags) == 0 ? Int[] : after_vertices
    # Inverse source is the forward image. Local labels of that image are 1..n in after_vertices order.
    inverse_before = witness.after_frags
    inverse_after = witness.before_frags
    if _vmax(inverse_before) == 0
        inverse_vertices = Int[]
    end
    inverse_extra = _copy_extra(witness.extra)
    if haskey(inverse_extra, "label_map") && haskey(inverse_extra, "inverse_label_map")
        forward_map = inverse_extra["label_map"]
        inverse_extra["label_map"] = inverse_extra["inverse_label_map"]
        inverse_extra["inverse_label_map"] = forward_map
    end
    # Added crossings are fresh locals above the surviving ones. _vertex_maps sorts
    # the selected globals and would hand the inverse the wrong surviving labels.
    # Carry the output words themselves and skip that renumbering.
    if get(witness.extra, "global_inverse", false)
        inverse_before = mapped_after
        inverse_after = [Int[x == 0 ? 0 : sign(x) * get(old_to_new, abs(x), abs(x)) for x in frag] for frag in mapped_before]
        inverse_vertices = Int[]
        inverse_extra["label_map"] = Dict{Int,Int}()
    end
    inverse = LCMoveWitness(
        witness.inverse_id == "" ? witness.row_id : witness.inverse_id,
        witness.row_id,
        witness.digest,
        witness.relation,
        witness.layer,
        witness.names,
        inverse_locations,
        inverse_vertices,
        inverse_before,
        inverse_after,
        witness.capacity,
        _fingerprint(out_lc, _roles_of(output)),
        witness.label_policy,
        -witness.genus_delta,
        inverse_extra,
    )
    correspondence = Dict{String,Any}(
        "label_map" => old_to_new,
        "removed_labels" => sort(collect(removed)),
        "created_labels" => after_vertices,
        "rotations" => rotations,
        "place" => assigned,
        "relation" => witness.relation,
        "row" => witness.row_id,
    )
    LCMoveResult(output, inverse, correspondence)
end

function _expand_gaps(input, row, selections)
    _require_selection(row, selections)
    gaps = selections["gaps"]
    lc = _lc_of(input)
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    target = _frags(row["target"])
    locations = Tuple{Int,Int}[]
    for fragment in target
        role = String(fragment["role"])
        haskey(gaps, role) || throw(LCMoveError(:missing_selection, "gap missing for role", Dict("role" => role)))
        pair = gaps[role]
        comp = String(pair[1])
        start = Int(pair[2])
        haskey(index, comp) || throw(LCMoveError(:missing_selection, "gap component is not in the LC", Dict("component" => comp)))
        push!(locations, (index[comp], start))
    end
    # Gaps alone are not a support. Select only embeddings certified by a face.
    [w for w in _expand_face_edges(input, row, selections) if w.locations == locations]
end

function _owner_of(words, label, sign_char)
    for (name, word) in words
        for entry in word
            abs(entry) == label && (entry > 0 ? "+" : "-") == sign_char && return name
        end
    end
    nothing
end

# Equal-gap inserts are applied from the left, so the partner word comes out reversed.
# Emit each gap last-to-first. The new curve's order is the one that keeps χ = 2-2g.
function _emit_gap_inserts!(locations, before_frags, after_frags, index, ops)
    grouped = Dict{Tuple{String,Int},Vector{Int}}()
    keys = Tuple{String,Int}[]
    for (partner, gap, signed_local) in ops
        key = (partner, gap)
        if !haskey(grouped, key)
            push!(keys, key)
            grouped[key] = Int[]
        end
        push!(grouped[key], signed_local)
    end
    for key in keys
        partner, gap = key
        for signed_local in reverse(grouped[key])
            push!(locations, (index[partner], gap))
            push!(before_frags, Int[])
            push!(after_frags, [-signed_local])
        end
    end
end

# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.2 — F_q(s,t)=(s+q*t mod 1, t); reconstructed; q>0 is one right parallel of the support at every visit, q<0 one left parallel
function _expand_affine(input, row, selections)
    _require_selection(row, selections)
    q = Int(selections["q"])
    acted = String[string(x) for x in selections["acted_components"]]
    support = selections["support"]
    lc = _lc_of(input)
    core = support isa AbstractString ? String(support) : throw(LCMoveError(
        :missing_selection,
        "support must be a component name or an oriented edge list; edge-list extraction is the next concrete case",
        Dict("support_type" => string(typeof(support))),
    ))
    haskey(lc.words, core) || throw(LCMoveError(:missing_selection, "support component is missing", Dict("support" => core)))
    q == 0 && return [_witness_from(input, row, _component_names(lc), Tuple{Int,Int}[], Int[], Vector{Vector{Int}}(), Vector{Vector{Int}}(), Int[];
        extra=Dict{String,Any}("identity" => true))]
    # F_q shears the acted curve along the support. Each other support crossing is copied
    # onto that crossing's partner, |q| times. The partner is the opposite passage, not the
    # support itself: a visit on the support does not split a face. Every visit uses the same
    # bank. q>0 uses the right bank and q<0 the left. A forward reading of that bank is placed
    # before the visit; a backward reading is placed after it. Splitting the visits across
    # both banks is not this twist, and the rotation system then has the wrong Euler
    # characteristic for this genus.
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    words = lc.words
    next_local = 1
    local_of = Dict{Int,Int}()
    function local_label(global_label)
        haskey(local_of, global_label) && return local_of[global_label]
        local_of[global_label] = next_local
        next_local += 1
        local_of[global_label]
    end
    before_frags = Vector{Vector{Int}}()
    after_frags = Vector{Vector{Int}}()
    locations = Tuple{Int,Int}[]
    # Partner inserts are emitted after every visit is known. An after-gap stack
    # is reversed so the first visit stays next to the support crossing; the kernel
    # writes equal gaps from the left.
    visits = NamedTuple{(:ops, :before, :after, :pos, :curve), Tuple{Vector{Tuple{String,Int,Int,Bool}}, Vector{Int}, Vector{Int}, Int, String}}[]
    core_word = words[core]
    # Existing visits take the low local ids. Fresh ids then sit strictly above
    # before.vertex_count, which is the range _vertex_maps treats as new labels.
    for curve in acted
        haskey(words, curve) || throw(LCMoveError(:missing_selection, "acted component is missing", Dict("component" => curve)))
        for entry in words[curve]
            _owner_of(words, abs(entry), entry > 0 ? "-" : "+") == core || continue
            local_label(abs(entry))
        end
    end
    for curve in acted
        haskey(words, curve) || throw(LCMoveError(:missing_selection, "acted component is missing", Dict("component" => curve)))
        word = words[curve]
        for (pos, entry) in enumerate(word)
            partner = _owner_of(words, abs(entry), entry > 0 ? "-" : "+")
            partner == core || continue
            core_pos = findfirst(x -> abs(x) == abs(entry), core_word)
            core_pos === nothing && throw(LCMoveError(:expansion, "core visit missing", Dict("label" => abs(entry))))
            rest = Int[]
            for offset in 1:(length(core_word) - 1)
                rest_entry = core_word[mod1(core_pos + offset, length(core_word))]
                other = _owner_of(words, abs(rest_entry), rest_entry > 0 ? "-" : "+")
                other === nothing && throw(LCMoveError(
                    :expansion, "support visit has no partner", Dict("visit" => rest_entry)))
                other == curve && continue
                push!(rest, rest_entry)
            end
            direction = sign(q) * (entry > 0 ? 1 : -1)
            side = q > 0 ? "right" : "left"
            insert_word = Int[]
            for _ in 1:abs(q)
                append!(insert_word, direction > 0 ? rest : reverse(-rest))
            end
            inserted_signed = Int[]
            ops = Tuple{String,Int,Int,Bool}[]
            for inserted in insert_word
                core_letter = direction > 0 ? inserted : -inserted
                partner_name = _owner_of(words, abs(core_letter), core_letter > 0 ? "-" : "+")
                partner_name === nothing && throw(LCMoveError(
                    :expansion, "support visit has no partner", Dict("visit" => core_letter)))
                ppos = findfirst(x -> abs(x) == abs(core_letter), words[partner_name])
                ppos === nothing && throw(LCMoveError(
                    :expansion, "partner visit missing", Dict("visit" => core_letter, "partner" => partner_name)))
                fresh = next_local
                next_local += 1
                acted_signed = sign(inserted) * fresh
                push!(inserted_signed, acted_signed)
                gap = _parallel_after(core_letter > 0, side) ? ppos + 1 : ppos
                push!(ops, (partner_name, gap, acted_signed, gap != ppos))
            end
            visit = [sign(entry) * local_label(abs(entry))]
            # Forward along the bank, the parallel is read before the visit.
            # Backward, it is read after the visit. Both readings stay on `side`.
            after = direction > 0 ? vcat(inserted_signed, visit) : vcat(visit, inserted_signed)
            push!(visits, (ops=ops, before=visit, after=after, pos=pos, curve=curve))
        end
    end
    for after_bank in (false, true)
        seq = after_bank ? reverse(visits) : visits
        for rec in seq
            chosen = Tuple{String,Int,Int}[(op[1], op[2], op[3]) for op in rec.ops if op[4] == after_bank]
            isempty(chosen) && continue
            _emit_gap_inserts!(locations, before_frags, after_frags, index, chosen)
        end
    end
    for rec in visits
        push!(locations, (index[rec.curve], rec.pos))
        push!(before_frags, rec.before)
        push!(after_frags, rec.after)
    end
    isempty(locations) && throw(LCMoveError(
        :expansion, "acted components do not meet the support", Dict("support" => core, "acted" => acted)))
    vertices = zeros(Int, next_local - 1)
    for (global_label, local_label) in local_of
        vertices[local_label] = global_label
    end
    # Fresh locals are not in vertices; _vertex_maps assigns them only if they exceed before.vertex_count.
    # Here before fragments use only existing locals, so before.vertex_count is the max existing local,
    # and fresh locals in after fragments are assigned by _vertex_maps if they are > that count.
    capacity = collect(1:length(before_frags))
    [_witness_from(input, row, names, locations, vertices[1:_vmax(before_frags)], before_frags, after_frags, capacity;
        extra=Dict{String,Any}("q" => q, "support" => core, "acted" => acted, "global_inverse" => true))]
end

function _shift_label(label, genus)
    limit = _k_limit(genus)
    label <= genus && return label
    label <= limit && return label + 1
    label + 3
end

# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.3 — empty band is a left parallel of a_j; reconstructed
function _expand_glue(input, row, selections)
    _require_selection(row, selections)
    input isa HeegaardLayer || throw(LCMoveError(:layer, "handle slide requires an HC", Dict("row" => row["id"])))
    curves = selections["curves"]
    length(curves) == 2 || throw(LCMoveError(:missing_selection, "curves must be the pair (a_i, a_j)", Dict()))
    a_i, a_j = String(curves[1]), String(curves[2])
    a_i == a_j && throw(LCMoveError(:applicability, "a slide needs two distinct curves", Dict()))
    role_i = get(input.roles, a_i, nothing)
    role_j = get(input.roles, a_j, nothing)
    (role_i === nothing || role_j === nothing || role_i != role_j) && throw(LCMoveError(
        :applicability, "handle slide requires two curves of the same family",
        Dict("a_i" => a_i, "a_j" => a_j, "roles" => Dict(a_i => role_i, a_j => role_j))))
    band = selections["band"]
    for step in band
        crossed = String(step["component"])
        get(input.roles, crossed, nothing) == role_i && throw(LCMoveError(
            :applicability, "band meets another curve of the same family", Dict("component" => crossed)))
    end
    lc = input.lc
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    words = lc.words
    next_local = 1
    local_of = Dict{Int,Int}()
    function use_local(global_label)
        haskey(local_of, global_label) && return local_of[global_label]
        local_of[global_label] = next_local
        next_local += 1
        local_of[global_label]
    end
    before_frags = Vector{Vector{Int}}()
    after_frags = Vector{Vector{Int}}()
    locations = Tuple{Int,Int}[]
    # a_i is replaced by its own visits followed by a parallel copy of a_j and the band crossings.
    before_i = [sign(entry) * use_local(abs(entry)) for entry in words[a_i]]
    after_i = copy(before_i)
    function add_fresh(partner, gap, sign_on_new)
        fresh = next_local
        next_local += 1
        push!(locations, (index[partner], gap))
        push!(before_frags, Int[])
        push!(after_frags, [-sign_on_new * fresh])
        fresh
    end
    # Same passage signs as a_j, on its left. The opposite signs sit off the surface:
    # the rotation system then has the wrong Euler characteristic for this genus.
    ops = Tuple{String,Int,Int}[]
    for entry in words[a_j]
        partner = _owner_of(words, abs(entry), entry > 0 ? "-" : "+")
        partner === nothing && throw(LCMoveError(:expansion, "a_j visit has no partner", Dict("visit" => entry)))
        ppos = findfirst(x -> abs(x) == abs(entry), words[partner])
        ppos === nothing && throw(LCMoveError(
            :expansion, "partner visit missing", Dict("visit" => entry, "partner" => partner)))
        sign_on_new = sign(entry)
        gap = _parallel_after(entry > 0, "left") ? ppos + 1 : ppos
        fresh = next_local
        next_local += 1
        signed = sign_on_new * fresh
        push!(after_i, signed)
        push!(ops, (partner, gap, signed))
    end
    _emit_gap_inserts!(locations, before_frags, after_frags, index, ops)
    for step in band
        crossed = String(step["component"])
        gap = Int(step["at"]) + 1
        approach = String(get(step, "approach", "right"))
        sign_on_new = approach == "right" ? 1 : -1
        fresh = add_fresh(crossed, gap, sign_on_new)
        push!(after_i, sign_on_new * fresh)
    end
    push!(locations, (index[a_i], 1))
    push!(before_frags, before_i)
    push!(after_frags, after_i)
    vertices = zeros(Int, max(1, next_local - 1))
    for (global_label, local_label) in local_of
        vertices[local_label] = global_label
    end
    used = _vmax(before_frags)
    capacity = collect(1:length(before_frags))
    [_witness_from(input, row, names, locations, vertices[1:used], before_frags, after_frags, capacity;
        extra=Dict{String,Any}("a_i" => a_i, "a_j" => a_j, "global_inverse" => true))]
end

# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.4 — ladder_end crossing tokens; reconstructed
function _ladder_letter(token, env)
    text = String(token)
    sig = 1
    if startswith(text, "-")
        sig = -1
        text = text[2:end]
    elseif startswith(text, "+")
        text = text[2:end]
    end
    haskey(env, text) || throw(LCMoveError(:expansion, "unknown ladder token", Dict("token" => text)))
    sig * env[text]
end

_ladder_word(tokens, env) = Int[_ladder_letter(token, env) for token in tokens]

function _expand_k_inclusion(input, row, selections)
    _require_selection(row, selections)
    delta = Int(get(get(row, "expansion", Dict()), "genus_delta", row["genus_delta"]))
    delta == 1 || throw(LCMoveError(
        :expansion,
        "逆安定化は、順方向の witness が残す一穴トーラスの支持でだけ適用する",
        Dict("genus_delta" => delta, "row" => row["id"]),
    ))
    input isa HeegaardLayer || throw(LCMoveError(:layer, "stabilization requires an HC", Dict()))
    String(selections["k_identification"]) == "ladder_end" || throw(LCMoveError(
        :missing_selection, "only the ladder_end identification is in the row data",
        Dict("k_identification" => selections["k_identification"])))
    pattern = get(get(row, "expansion", Dict()), "pattern", nothing)
    pattern isa AbstractDict || throw(LCMoveError(
        :expansion, "ladder_end pattern is not in the row", Dict("row" => row["id"])))
    lc = input.lc
    g = lc.genus
    g >= 1 || throw(LCMoveError(:applicability, "ladder_end stabilization starts at genus >= 1", Dict("genus" => g)))
    disk = String(selections["disk"])
    active = Dict(k => copy(v) for (k, v) in lc.words if !isempty(v))
    faces, _ = _lc_faces(active)
    any(face -> face.id == disk, faces) || throw(LCMoveError(
        :missing_selection, "disk is not a face of this diagram", Dict("disk" => disk)))
    lname = "l$g"
    mname = "m$g"
    haskey(lc.words, lname) && haskey(lc.words, mname) || throw(LCMoveError(
        :expansion, "ladder end is missing l_g or m_g", Dict("genus" => g)))
    pos = findfirst(==(-g), lc.words[lname])
    pos === nothing && throw(LCMoveError(:expansion, "lm visit not found on l_g", Dict("component" => lname)))
    limit = _k_limit(g)
    non_k = [entry for entry in lc.words[mname] if abs(entry) > limit]
    length(non_k) == 1 || throw(LCMoveError(
        :expansion, "ladder_end needs one non-K visit on m_g",
        Dict("component" => mname, "visits" => non_k)))
    entry = only(non_k)
    partner = _owner_of(lc.words, abs(entry), entry > 0 ? "-" : "+")
    (partner === nothing || partner == lname) && throw(LCMoveError(
        :expansion, "m_g visit has no ladder partner", Dict("visit" => entry)))
    ref = _lc_k_constants(g + 1)
    old_k = Int[sign(e) * _shift_label(abs(e), g) for e in lc.words[lname] if abs(e) <= limit]
    fresh_k = [letter for letter in ref[lname] if !(abs(letter) in Set(abs.(old_k)))]
    length(fresh_k) == 1 || throw(LCMoveError(
        :expansion, "new l_g K letter is not unique", Dict("letters" => fresh_k)))
    new_signed = only(fresh_k)
    labels = unique(abs(e) for word in values(lc.words) for e in word)
    forward = Dict{Int,Int}()
    backward = Dict{Int,Int}()
    shifted_max = 0
    for label in labels
        image = _shift_label(label, g)
        forward[label] = image
        backward[image] = label
        shifted_max = max(shifted_max, image)
    end
    nxt = max(shifted_max, _k_limit(g + 1)) + 1
    env = Dict("A" => nxt, "B" => nxt + 1, "C" => nxt + 2, "D" => nxt + 3, "E" => nxt + 4)
    l_new = "l$(g + 1)"
    m_new = "m$(g + 1)"
    n_name = "n$g"
    aname = "alpha$(g + 1)"
    bname = "beta$(g + 1)"
    for name in (l_new, m_new, n_name, aname, bname)
        haskey(lc.words, name) && throw(LCMoveError(
            :component_name, "stabilization name is already used", Dict("name" => name)))
    end
    nref = ref[n_name]
    length(nref) == 2 || throw(LCMoveError(
        :expansion, "n_g constant is not a pair", Dict("word" => nref)))
    created_words = Dict(
        l_new => vcat(ref[l_new], _ladder_word(pattern["l_suffix"], env)),
        m_new => vcat(ref[m_new], _ladder_word(pattern["m_suffix"], env)),
        n_name => vcat(nref[1:1], _ladder_word(pattern["n_between"], env), nref[2:2]),
        aname => _ladder_word(pattern["alpha0"], env),
        bname => _ladder_word(pattern["beta0"], env),
    )
    touch = only(_ladder_word([pattern["old_touch"]], env))
    ppos = findfirst(x -> abs(x) == abs(entry), lc.words[partner])
    ppos === nothing && throw(LCMoveError(
        :expansion, "partner visit missing", Dict("visit" => entry, "partner" => partner)))
    # cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.4 — n meets the partner on the right bank of the m_g visit; reconstructed
    gap = _parallel_after(entry > 0, "right") ? ppos + 1 : ppos
    names = vcat(_component_names(lc), [l_new, m_new, n_name, aname, bname])
    index = Dict(name => i for (i, name) in enumerate(names))
    before_frags = Vector{Vector{Int}}()
    after_frags = Vector{Vector{Int}}()
    locations = Tuple{Int,Int}[]
    function push_frag(name, start, pre, post)
        push!(locations, (index[name], start))
        push!(before_frags, pre)
        push!(after_frags, post)
    end
    push_frag(lname, pos, Int[-g], Int[-g, new_signed])
    push_frag(partner, gap, Int[], Int[touch])
    for name in (l_new, m_new, n_name, aname, bname)
        push_frag(name, 1, Int[], created_words[name])
    end
    extra = Dict{String,Any}(
        "label_map" => forward,
        "inverse_label_map" => backward,
        "created_components" => [l_new, m_new, n_name, aname, bname],
        "new_roles" => Dict{String,Symbol}(aname => :alpha, bname => :beta),
        "disk" => disk,
        "acted" => [lname, partner],
    )
    [_witness_from(input, row, names, locations, Int[], before_frags, after_frags,
        collect(1:length(before_frags)); extra=extra)]
end

function _kind_ok(lc, kind, name)
    kind == "any" && return !isempty(get(lc.words, name, Int[]))
    kind == "k" && return _is_k(lc, name)
    kind == "input" && return !_is_k(lc, name) && !isempty(get(lc.words, name, Int[]))
    throw(LCMoveError(:bindings, "unknown binding kind", Dict("kind" => kind)))
end

function _derived_digon(push_forward, cross_forward)
    # cf: .research-scratch/PJ260922/261111-gpt-pro-review.md §4 — left-face table;
    # checked against all 384 insertions on g1-base (48 opposite-sign passes).
    key = string(push_forward ? "+" : "-", cross_forward ? "+" : "-")
    entry = load_heegaard_moves()["shared"]["digon_left_face"][key]
    Int.(entry["push"]), Int.(entry["cross"])
end

function _expand_face_edges(input, row, selections)
    lc = _lc_of(input)
    names = _component_names(lc)
    index = Dict(name => i for (i, name) in enumerate(names))
    active = Dict(k => copy(v) for (k, v) in lc.words if !isempty(v))
    faces, _ = _lc_faces(active)
    occ = _lc_rotation_occ(active)
    target = _frags(get(row, "target", Dict("fragments" => [])))
    fixed = any(!isempty(_frag_word(fragment)) for fragment in target)
    fixed_c = fixed ? _frag_word(only(fragment for fragment in target if String(fragment["role"]) == "C")) : Int[]
    fixed_k = fixed ? _frag_word(only(fragment for fragment in target if String(fragment["role"]) == "K")) : Int[]
    bindings = get(row, "bindings", Dict{String,Any}())
    c_kind = String(get(bindings, "C", "any"))
    k_kind = String(get(bindings, "K", "any"))
    witnesses = LCMoveWitness[]
    for face in faces
        # A place-only circle makes this region punctured; its side and the
        # path around it are not encoded by two edge occurrences alone.
        face.id in values(lc.places) && continue
        edges = _move_left_boundary(active, occ, face)
        for push_edge in edges, cross_edge in edges
            push_edge.component == cross_edge.component && continue
            _kind_ok(lc, c_kind, push_edge.component) || continue
            _kind_ok(lc, k_kind, cross_edge.component) || continue
            if input isa HeegaardLayer
                a = get(input.roles, push_edge.component, nothing)
                b = get(input.roles, cross_edge.component, nothing)
                a !== nothing && a == b && continue
            end
            c_pat, k_pat = _derived_digon(push_edge.forward, cross_edge.forward)
            fixed && (c_pat != fixed_c || k_pat != fixed_k) && continue
            c_len = length(lc.words[push_edge.component])
            k_len = length(lc.words[cross_edge.component])
            c_at = push_edge.forward ? mod1(push_edge.position + 1, c_len) : push_edge.position
            k_at = cross_edge.forward ? mod1(cross_edge.position + 1, k_len) : cross_edge.position
            extra = Dict{String,Any}(
                "range" => "one face, two distinct boundary edges",
                "face" => face.id,
                "push" => push_edge.component,
                "cross" => cross_edge.component,
                "support" => Dict("kind" => "disk_band", "face" => face.id,
                    "start" => collect(push_edge.corner), "finish" => collect(cross_edge.corner)),
            )
            push!(witnesses, _witness_from(
                input, row, names,
                [(index[push_edge.component], c_at), (index[cross_edge.component], k_at)],
                Int[], [Int[], Int[]], [c_pat, k_pat], [1, 2]; extra=extra))
        end
    end
    witnesses
end

function _parallel_after(c_sign_positive, side)
    seq = [("-", "out"), ("+", "out"), ("-", "in"), ("+", "in")]
    c_sign = c_sign_positive ? "+" : "-"
    i = findfirst(==( (c_sign, "out") ), seq)
    step = side == "left" ? 1 : -1
    seq[mod1(i + step, 4)][2] == "out"
end

function _expand_parallel(input, row, selections)
    lc = _lc_of(input)
    curve = get(selections, "curve", nothing)
    side = get(selections, "side", nothing)
    names0 = _component_names(lc)
    curves = curve === nothing ? [n for n in names0 if !isempty(get(lc.words, n, Int[]))] : [String(curve)]
    sides = side === nothing ? ["left", "right"] : [String(side)]
    max_label = maximum(vcat([0], [abs(x) for w in values(lc.words) for x in w]))
    witnesses = LCMoveWitness[]
    for c in curves, s in sides
        word = get(lc.words, c, Int[])
        isempty(word) && continue
        s in ("left", "right") || throw(LCMoveError(:missing_selection, "side must be left or right", Dict("side" => s)))
        copy_name = c * (s == "left" ? "_l" : "_r")
        haskey(lc.words, copy_name) && throw(LCMoveError(:component_name, "parallel copy name is already used", Dict("name" => copy_name)))
        partners = Dict{String,Vector{Int}}()
        new_local = max_label
        copy_word = Int[]
        for entry in word
            new_local += 1
            sign = entry > 0 ? 1 : -1
            push!(copy_word, sign * new_local)
            placed = false
            for (name, partner) in lc.words
                name == c && continue
                for (i, y) in enumerate(partner)
                    abs(y) == abs(entry) || continue
                    w = get!(partners, name, copy(partner))
                    insert_at = _parallel_after(entry > 0, s) ? i + 1 : i
                    # later inserts on the same word shift; record and apply from the right below
                    placed = true
                end
            end
            placed || throw(LCMoveError(:expansion, "crossing has no partner", Dict("label" => abs(entry), "curve" => c)))
        end
        # rebuild partner inserts from the right so indices stay valid
        jobs = Dict{String,Vector{Tuple{Int,Int,Bool}}}()
        new_local = max_label
        for entry in word
            new_local += 1
            sign = entry > 0 ? 1 : -1
            partner_entry = -sign * new_local
            for (name, partner) in lc.words
                name == c && continue
                for (i, y) in enumerate(partner)
                    abs(y) == abs(entry) || continue
                    push!(get!(jobs, name, Tuple{Int,Int,Bool}[]), (i, partner_entry, _parallel_after(entry > 0, s)))
                end
            end
        end
        after_of = Dict{String,Vector{Int}}()
        for (name, items) in jobs
            w = copy(lc.words[name])
            for (i, signed_new, after) in sort(items; by=first, rev=true)
                insert!(w, after ? i + 1 : i, signed_new)
            end
            after_of[name] = w
        end
        names = vcat(names0, [copy_name])
        before = Vector{Vector{Int}}()
        after = Vector{Vector{Int}}()
        locations = Tuple{Int,Int}[]
        for (i, name) in enumerate(names0)
            old = get(lc.words, name, Int[])
            push!(before, copy(old))
            push!(after, get(after_of, name, copy(old)))
            push!(locations, (i, 1))
        end
        push!(before, Int[])
        push!(after, copy_word)
        push!(locations, (length(names), 1))
        extra = Dict{String,Any}(
            "range" => "every nonempty curve, both sides",
            "curve" => c,
            "side" => s,
            "created_components" => [copy_name],
        )
        push!(witnesses, _witness_from(
            input, row, names, locations, collect(1:max_label), before, after,
            collect(1:length(names)); extra=extra))
    end
    witnesses
end

function _expand(input, row, selections)
    op = String(get(get(row, "expansion", Dict()), "op", "literal"))
    if op == "literal"
        return _search_literal(input, row)
    elseif op == "gaps"
        if haskey(selections, "gaps")
            return _expand_gaps(input, row, selections)
        end
        return _expand_face_edges(input, row, selections)
    elseif op == "face_edges"
        return _expand_face_edges(input, row, selections)
    elseif op == "parallel"
        return _expand_parallel(input, row, selections)
    elseif op == "affine"
        haskey(selections, "q") || return LCMoveWitness[]
        return _expand_affine(input, row, selections)
    elseif op == "glue"
        haskey(selections, "curves") || return LCMoveWitness[]
        return _expand_glue(input, row, selections)
    elseif op == "k_inclusion"
        haskey(selections, "disk") || haskey(selections, "support") || return LCMoveWitness[]
        return _expand_k_inclusion(input, row, selections)
    end
    throw(LCMoveError(:expansion, "unknown expansion op", Dict("op" => op, "row" => row["id"])))
end

function movable(input::Union{LickorishConfiguration,HeegaardLayer}, row::AbstractDict; kwargs...)
    _refuse_marks(input, row)
    if String(row["layer"]) == "HC" && !(input isa HeegaardLayer)
        throw(LCMoveError(:layer, "HC の行には役割を持つ HC が要る", Dict("row" => row["id"])))
    end
    selections = _selection_dict(; kwargs...)
    _expand(input, row, selections)
end

function movable(input::Union{LickorishConfiguration,HeegaardLayer}, row_id::AbstractString; kwargs...)
    movable(input, heegaard_move_row(row_id); kwargs...)
end

function move(input::Union{LickorishConfiguration,HeegaardLayer}, witness::LCMoveWitness)
    row = heegaard_move_row(witness.row_id)
    _refuse_marks(input, row)
    _apply_witness(input, witness)
end

function lc_text(lc::LickorishConfiguration)
    lines = String["genus: $(lc.genus)"]
    for name in _component_names(lc)
        haskey(lc.words, name) || continue
        body = join([entry > 0 ? "+$entry" : string(entry) for entry in lc.words[name]], " ")
        push!(lines, "$name : $body")
    end
    for name in sort(collect(keys(lc.places)))
        push!(lines, "place $name : $(lc.places[name])")
    end
    join(lines, "\n") * "\n"
end

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/overlay/lib/topology.jl/docs/specs/heegaard-code.md §5.1 — encoding はテキストの選択であり、役割も語も変えない
function hc_text(hc::HeegaardLayer; lc_filename::AbstractString, encoding::Symbol=:auto)
    validate_hc_input(hc)
    encoding in (:auto, :standard, :explicit) || throw(ArgumentError("unknown HC encoding: $encoding"))
    (isempty(lc_filename) || occursin(r"[\r\n#]", lc_filename) || isabspath(lc_filename)) &&
        throw(ArgumentError("lc_filename must be a nonempty relative HC-safe path"))
    aligned = is_standard(hc)
    encoding == :standard && !aligned &&
        throw(ArgumentError("cannot serialize a nonstandard binding as standard"))
    compact = encoding != :explicit && aligned
    lines = String["format : hc/2", "lc : $lc_filename"]
    compact && push!(lines, "standard")
    for name in hc_component_names(hc.lc)
        haskey(hc.roles, name) || continue
        compact && hc.roles[name] == :alpha && continue
        role = hc.roles[name] == :alpha ? "alpha" : hc.roles[name] == :beta ? "beta" :
            throw(ArgumentError("role is not alpha or beta"))
        push!(lines, "role $name : $role")
    end
    for name in sort(collect(keys(hc.disks)))
        push!(lines, "disk $name : " * join(hc.disks[name], " "))
    end
    for name in sort(collect(keys(hc.points)))
        push!(lines, "point $name : $(hc.points[name])")
    end
    for name in sort(collect(keys(hc.weights)))
        push!(lines, "weight $name : $(hc.weights[name])")
    end
    join(lines, "\n") * "\n"
end

function write_lc(path::AbstractString, lc::LickorishConfiguration)
    open(path, "w") do io
        write(io, lc_text(lc))
    end
    path
end

function write_hc(path::AbstractString, hc::HeegaardLayer; lc_filename::AbstractString, encoding::Symbol=:auto)
    text = hc_text(hc; lc_filename=lc_filename, encoding=encoding)
    open(path, "w") do io
        write(io, text)
    end
    path
end

export LCMoveError, LCMoveWitness, LCMoveResult,
       load_heegaard_moves, heegaard_move_row, heegaard_move_ids,
       movable, move, lc_text, hc_text, write_lc, write_hc
