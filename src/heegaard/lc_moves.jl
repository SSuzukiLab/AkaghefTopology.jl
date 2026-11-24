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
include("lc_moves_ribbon.jl")
include("lc_moves_annulus.jl")
include("lc_moves_pl.jl")
include("lc_moves_surface.jl")
include("lc_moves_roles.jl")
include("lc_moves_marks.jl")
include("lc_moves_mark_edit.jl")
include("lc_moves_sequence.jl")
include("lc_moves_prepare.jl")
include("lc_moves_certificates.jl")
# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.2–§6.4 — candidates when a selection is omitted
include("lc_moves_enumerate.jl")

_k_limit(g::Integer) = g <= 0 ? 0 : 3 * g - 2
_lc_of(lc::LickorishCode) = lc
_lc_of(hc::HeegaardCode) = hc.lc
_roles_of(::LickorishCode) = Dict{String,Symbol}()
_roles_of(hc::HeegaardCode) = hc.roles
_marks(lc::LickorishCode) = (disks=lc.disks, points=lc.points, weights=lc.weights)
_marks(hc::HeegaardCode) = (disks=hc.disks, points=hc.points, weights=hc.weights)

function _component_names(lc::LickorishCode)
    k = lc.genus > 0 ? _lc_k_names(lc.genus) : String[]
    rest = sort([n for n in keys(lc.words) if !(n in k)])
    String[vcat(k, rest)...]
end

function _fingerprint(lc::LickorishCode, roles=Dict{String,Symbol}())
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

_input_fingerprint(input) = _fingerprint(_lc_of(input),_roles_of(input))

function _missing_marks(input)
    marks = _marks(input)
    missing = String[]
    isempty(marks.disks) || push!(missing, "disk")
    isempty(marks.points) || push!(missing, "point")
    isempty(marks.weights) || push!(missing, "weight")
    missing
end

# Free marks never gate a geometric move (owner decision 2026-11-12).
_refuse_marks(input,row;policy=nothing) = nothing

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
            delete!(used, name)
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
            any(face.n == 2 && face.labels == label_set && !(face.id in values(lc.places)) for face in faces) || return false
        elseif check == "empty_triangle"
            any(face.n == 3 && face.labels == label_set && !(face.id in values(lc.places)) for face in faces) || return false
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
    extra["origin_row_id"]=String(row["id"])
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
        _input_fingerprint(input),
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
    # Join bindings through already matched labels instead of materializing
    # the Cartesian product of all curves (which explodes at genus 28).
    assignment = Dict{String,String}()
    used_names = Set{String}()
    pools = Dict(role => _candidates_for(lc,kinds[role]) for role in roles)
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
            for name in pools[role]
                name in used_names && continue
                assignment[role] = name
                push!(used_names,name)
                for (_, next) in _consistent_windows(words[name], patterns[role], subst)
                    visit(role_index + 1, next)
                end
                delete!(used_names,name)
                delete!(assignment,role)
            end
    end
    visit(1, Dict{Int,Int}())
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
        elseif start > 1 && start + length(frag) - 1 == n && isempty(after_frags[i])
            # An empty image is deleted from the end, so its next neighbor is
            # the first letter. Reinsertion uses that cyclic gap and writes the
            # fragment at the front. A non-empty image stays where the splice
            # wrote it. Rotating that word (vertex.100 on a suffix) shifts a
            # restoration that was already linear.
            rotations[comp] = start - 1
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
    lc = LickorishCode(genus, words, places, copy(lc0.disks), copy(lc0.points), copy(lc0.weights), lc0.source_path)
    validate_lc_input(lc)
    if any(!isempty(word) for word in values(words))
        _lc_validate_surface!(lc)
    end
    if input isa HeegaardCode
        roles = copy(input.roles)
        return HeegaardCode(lc, roles, copy(input.disks), copy(input.points), copy(input.weights), input.hc_path)
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
    output isa HeegaardCode || return output
    roles = copy(output.roles)
    for name in drop_names
        delete!(roles, name)
    end
    if witness.genus_delta > 0
        for (name, role) in get(witness.extra, "new_roles", Dict{String,Any}())
            roles[String(name)] = role isa Symbol ? role : Symbol(String(role))
        end
    end
    for t in get(witness.extra,"role_transfers",Any[])
        get(output.roles,t.from,nothing)==t.role || throw(LCMoveError(:role_transport_conflict,"source role differs",Dict("component"=>t.from)))
        haskey(roles,t.to) && throw(LCMoveError(:role_transport_conflict,"target already has a role",Dict("component"=>t.to)))
        delete!(roles,t.from);roles[t.to]=t.role
    end
    HeegaardCode(output.lc, roles, output.disks, output.points, output.weights, output.hc_path)
end

function _apply_witness(input, witness::LCMoveWitness)
    _input_fingerprint(input) == witness.fingerprint || throw(LCMoveError(
        :witness_mismatch, "witness was not built for these words, places and roles", Dict("row"=>witness.row_id)))
    haskey(witness.extra,"mark_patch") && return _move_apply_mark_patch(input,witness)
    if get(witness.extra, "identity", false)
        return LCMoveResult(input, witness, Dict{String,Any}("relation" => witness.relation, "identity" => true))
    end
    lc = _lc_of(input)
    if get(witness.extra,"annulus_recompute",false)
        inverse_plan=_move_affine_witness(input,heegaard_move_row(witness.row_id),
            witness.extra["annulus_atlas"],witness.extra["acted"],witness.extra["q"])
        return _restore_recorded_marks(_apply_witness(input,inverse_plan), witness)
    end
    if get(witness.extra,"embedded_annulus_recompute",false)
        inverse_plan=_move_embedded_affine(input,heegaard_move_row(witness.row_id),
            witness.extra["embedded_annulus_atlas"],witness.extra["acted"],witness.extra["q"])
        return _restore_recorded_marks(_apply_witness(input,inverse_plan), witness)
    end
    words = Dict(name => copy(get(lc.words, name, Int[])) for name in witness.names)
    result, old_to_new, after_vertices, mapped_before, mapped_after, locs, rotations = _kernel_apply(
        witness.names, words, witness.before_frags, witness.after_frags,
        witness.locations, witness.vertices, witness.capacity;
        label_map=_label_map_of(witness.extra))
    # The forward kernel may rotate a word so a wrapped fragment starts at 1,
    # and it leaves that rotation in the output. The inverse has to turn the
    # same word back, or the previous witness sees a different linear order.
    undo = get(witness.extra, "undo_rotation", nothing)
    if undo isa AbstractDict
        for (comp, name) in enumerate(witness.names)
            shift = Int(get(undo, String(name), 0))
            shift == 0 && continue
            comp > length(result) && continue
            result[comp] = _rot_right(result[comp], shift)
        end
    end
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
    parking=get(witness.extra,"support_direction",1)==1 ?
        get(witness.extra,"creation_parking",Dict{String,Vector{Tuple{Int,String,String}}}()) :
        Dict{String,Vector{Tuple{Int,String,String}}}()
    new_places = _move_transport_places(lc,new_words,new_places,old_to_new;
        digon_create_untouched_corners=parking)
    output = _attach_roles(
        _build_output(input, lc.genus + witness.genus_delta, new_words, new_places),
        witness, drop_names)
    output,mark_cells = _move_transport_marks(input,output,witness,old_to_new)
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
    if haskey(inverse_extra,"role_transfers")
        inverse_extra["role_transfers"]=[(from=t.to,to=t.from,role=t.role) for t in reverse(witness.extra["role_transfers"])]
    end
    inverse_extra["support_direction"]=-get(witness.extra,"support_direction",1)
    if haskey(inverse_extra, "label_map") && haskey(inverse_extra, "inverse_label_map")
        forward_map = inverse_extra["label_map"]
        inverse_extra["label_map"] = inverse_extra["inverse_label_map"]
        inverse_extra["inverse_label_map"] = forward_map
    end
    # A removal densely renumbers every higher survivor and forgets the deleted
    # labels. The local inverse pattern would mint new high labels, so the words
    # no longer match the previous witness. Carry the global fragments and the
    # reverse numbering whenever labels were removed or survivors moved.
    survivor_shifted = any(pair -> pair[1] != pair[2], pairs(old_to_new))
    removed_labels = setdiff(
        Set(abs.(reduce(vcat, values(lc.words); init=Int[]))),
        Set(keys(old_to_new)),
    )
    if get(witness.extra, "global_inverse", false) || survivor_shifted || !isempty(removed_labels)
        inverse_before = mapped_after
        inverse_after = mapped_before
        inverse_vertices = Int[]
        inverse_extra["label_map"] = Dict(v => k for (k, v) in old_to_new)
        inverse_extra["global_inverse"] = true
    end
    delete!(inverse_extra, "undo_rotation")
    undo_rotation = Dict{String,Int}()
    for (comp, shift) in enumerate(rotations)
        shift == 0 && continue
        undo_rotation[String(witness.names[comp])] = shift
    end
    isempty(undo_rotation) || (inverse_extra["undo_rotation"] = undo_rotation)
    if haskey(witness.extra,"annulus_after")
        inverse_extra["annulus_recompute"]=true
        inverse_extra["annulus_atlas"]=witness.extra["annulus_after"]
        inverse_extra["q"]=-witness.extra["q"]
    end
    if haskey(witness.extra,"embedded_annulus_after")
        inverse_extra["embedded_annulus_recompute"]=true
        inverse_extra["embedded_annulus_atlas"]=witness.extra["embedded_annulus_after"]
        inverse_extra["q"]=-witness.extra["q"]
    end
    inverse_extra["mark_before"] = _recorded_marks(input)
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
        _input_fingerprint(output),
        witness.label_policy,
        -witness.genus_delta,
        inverse_extra,
    )
    correspondence = Dict{String,Any}(
        "label_map" => old_to_new,
        "removed_labels" => sort(collect(setdiff(Set(abs.(reduce(vcat,values(lc.words);init=Int[]))),Set(keys(old_to_new))))),
        "created_labels" => sort(collect(setdiff(Set(abs.(reduce(vcat,values(out_lc.words);init=Int[]))),Set(values(old_to_new))))),
        "rotations" => rotations,
        "place" => assigned,
        "place_images" => Dict(n=>(old=f,new=new_places[n]) for(n,f) in lc.places if haskey(new_places,n)),
        "relation" => witness.relation,
        "row" => witness.row_id,
        "support" => get(witness.extra,"support",nothing),
        "support_direction" => get(witness.extra,"support_direction",1),
        "k_identification" => get(witness.extra,"k_identification",nothing),
        "role_transfers" => get(witness.extra,"role_transfers",Any[]),
        "mark_cells" => mark_cells,
        "attachment" => get(witness.extra,"attachment",nothing),
    )
    _restore_recorded_marks(LCMoveResult(output, inverse, correspondence), witness)
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
    acted = unique(String[string(x) for x in selections["acted_components"]])
    if q!=0 && input isa HeegaardCode && any(_is_k(input.lc,n) && haskey(input.roles,n) for n in acted)
        return _move_prepare_affine(input,row,selections,acted)
    end
    q != 0 && _move_require_separated(input, acted)
    support = selections["support"]
    lc = _lc_of(input)
    if support isa AbstractDict
        return [_move_embedded_affine(input,row,support,acted,q)]
    end
    core = support isa AbstractString ? String(support) : throw(LCMoveError(
        :missing_selection,
        "support must be a component name or an oriented edge list; edge-list extraction is the next concrete case",
        Dict("support_type" => string(typeof(support))),
    ))
    haskey(lc.words, core) || throw(LCMoveError(:missing_selection, "support component is missing", Dict("support" => core)))
    q == 0 && return [_witness_from(input, row, _component_names(lc), Tuple{Int,Int}[], Int[], Vector{Vector{Int}}(), Vector{Vector{Int}}(), Int[];
        extra=Dict{String,Any}("identity" => true))]
    [_move_affine_witness(input,row,_move_thin_annulus(lc,core),acted,q)]
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
    input isa HeegaardCode || throw(LCMoveError(:layer, "handle slide requires an HC", Dict("row" => row["id"])))
    curves = selections["curves"]
    length(curves) == 2 || throw(LCMoveError(:missing_selection, "curves must be the pair (a_i, a_j)", Dict()))
    a_i, a_j = String(curves[1]), String(curves[2])
    if _is_k(input.lc,a_i) && haskey(input.roles,a_i)
        return _move_prepare_slide(input,row,selections)
    end
    a_i == a_j && throw(LCMoveError(:applicability, "a slide needs two distinct curves", Dict()))
    role_i = get(input.roles, a_i, nothing)
    role_j = get(input.roles, a_j, nothing)
    (role_i === nothing || role_j === nothing || role_i != role_j) && throw(LCMoveError(
        :applicability, "handle slide requires two curves of the same family",
        Dict("a_i" => a_i, "a_j" => a_j, "roles" => Dict(a_i => role_i, a_j => role_j))))
    band = selections["band"]
    band isa AbstractDict && haskey(band,"route") && (band=band["route"])
    bands = band isa AbstractDict ? [band] :
        band isa AbstractVector && isempty(band) ? _move_face_bands(input.lc,a_i,a_j) :
        throw(LCMoveError(:path_order_required,"band needs start/finish darts and ordered crossing points",Dict()))
    return [_move_pants_expand(input,row,[a_i,a_j],b) for b in bands]
end

function _expand_k_inclusion(input, row, selections)
    _require_selection(row, selections)
    delta = Int(get(get(row, "expansion", Dict()), "genus_delta", row["genus_delta"]))
    delta == 1 || throw(LCMoveError(
        :expansion,
        "逆安定化は、順方向の witness が残す一穴トーラスの支持でだけ適用する",
        Dict("genus_delta" => delta, "row" => row["id"]),
    ))
    input isa HeegaardCode || throw(LCMoveError(:layer, "stabilization requires an HC", Dict()))
    String(selections["k_identification"])=="ladder_end" || throw(LCMoveError(:missing_selection,"K identification is not in the catalogue",Dict()))
    if get(row["expansion"],"support_policy",nothing)=="free_attachment_disk"
        return _move_prepare_stabilization(input,row,selections)
    end
    [_move_surface_inclusion(input,row,String(selections["disk"]))]
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
        # Creation uses interior edge ports and a narrow embedded band disjoint
        # from old vertices. Actual placed circles can first be parked near an
        # untouched old corner; their nesting and oriented sides are preserved.
        edges = _move_left_boundary(active, occ, face)
        for push_edge in edges, cross_edge in edges
            push_edge.component == cross_edge.component && continue
            _kind_ok(lc, c_kind, push_edge.component) || continue
            _kind_ok(lc, k_kind, cross_edge.component) || continue
            if input isa HeegaardCode
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
                "creation_parking" => Dict(face.id=>copy(face.corners)),
                "support" => Dict("kind" => "disk_band", "face" => face.id,
                    "port_u"=>1//2,"closed_band_avoids_old_vertices"=>true,
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
    # s*sigma=+1 is before, s*sigma=-1 is after.
    side in ("left","right") || throw(LCMoveError(:support,"side must be left or right",Dict("side"=>side)))
    (side == "left") != c_sign_positive
end

function _expand_parallel(input, row, selections)
    input isa HeegaardCode && throw(LCMoveError(:parallel_role_required,
        "parallel copy is an LC construction; a new HC component needs an explicit role and support transport",
        Dict("row"=>row["id"])))
    lc = _lc_of(input)
    isempty(lc.places) || throw(LCMoveError(:missing_place_transport,
        "parallel collar needs the image of every place-only curve",Dict("places"=>copy(lc.places))))
    curve = get(selections, "curve", nothing)
    side = get(selections, "side", nothing)
    names0 = _component_names(lc)
    curves = curve === nothing ? [n for n in names0 if !isempty(get(lc.words, n, Int[]))] : [String(curve)]
    sides = side === nothing ? ["left", "right"] : [String(side)]
    witnesses = LCMoveWitness[]
    for c in curves, s in sides
        word = get(lc.words,c,Int[])
        isempty(word) && continue
        length(unique(abs.(word)))==length(word) || throw(LCMoveError(:support,"parallel core has a self crossing",Dict("curve"=>c)))
        s in ("left","right") || throw(LCMoveError(:missing_selection,"side must be left or right",Dict()))
        target=c*(s=="left" ? "_l" : "_r")
        haskey(lc.words,target) && throw(LCMoveError(:component_name,"parallel copy name is already used",Dict("name"=>target)))
        run=(kind=:curve,component=c,edge=length(word),side=s=="left" ? 1 : -1,direction=1,survives=false)
        w=_move_compile_boundary(input,row,target,[run],Dict("kind"=>"collar","curve"=>c,"side"=>s))
        w.extra["curve"]=c;w.extra["side"]=s;w.extra["range"]="every nonempty simple curve, both sides"
        push!(witnesses,w)
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
    elseif op == "role_parallel"
        haskey(selections, "curve") || return _enumerate_role_parallel(input, row, selections)
        return _move_separate_role(input,row,selections)
    elseif op == "mark_edit"
        (haskey(selections, "name") && haskey(selections, "target")) || return _enumerate_mark_edit(input, row, selections)
        return [_move_mark_edit(input,row,selections)]
    elseif op == "affine"
        _affine_selection_complete(selections) || return _enumerate_affine(input, row, selections)
        return _expand_affine(input, row, selections)
    elseif op == "glue"
        haskey(selections, "curves") || return _enumerate_glue(input, row, selections)
        return _expand_glue(input, row, selections)
    elseif op == "k_inclusion"
        (haskey(selections, "disk") || haskey(selections, "support")) || return _enumerate_k_inclusion(input, row, selections)
        return _expand_k_inclusion(input, row, selections)
    end
    throw(LCMoveError(:expansion, "unknown expansion op", Dict("op" => op, "row" => row["id"])))
end

function movable(input::Union{LickorishCode,HeegaardCode}, row::AbstractDict; kwargs...)
    input isa HeegaardCode ? validate_hc_input(input) : validate_lc_input(input)
    _move_validate_input(_lc_of(input))
    _refuse_marks(input, row;policy=get(kwargs,:mark_transport,nothing))
    if String(row["layer"]) == "HC" && !(input isa HeegaardCode)
        throw(LCMoveError(:layer, "HC の行には役割を持つ HC が要る", Dict("row" => row["id"])))
    end
    selections = _selection_dict(; kwargs...)
    ws=_expand(input, row, selections)
    for w in ws
        haskey(selections,"mark_transport") && (w.extra["mark_transport"]=selections["mark_transport"])
    end
    ws
end

function movable(input::Union{LickorishCode,HeegaardCode}, row_id::AbstractString; kwargs...)
    movable(input, heegaard_move_row(row_id); kwargs...)
end

function move(input::Union{LickorishCode,HeegaardCode}, witness::LCMoveWitness)
    row = heegaard_move_row(witness.row_id)
    _refuse_marks(input, row;policy=get(witness.extra,"mark_transport",nothing))
    _move_checked_apply(input, witness)
end

function lc_text(lc::LickorishCode)
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
function hc_text(hc::HeegaardCode; lc_filename::AbstractString, encoding::Symbol=:auto)
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

function write_lc(path::AbstractString, lc::LickorishCode)
    open(path, "w") do io
        write(io, lc_text(lc))
    end
    path
end

function write_hc(path::AbstractString, hc::HeegaardCode; lc_filename::AbstractString, encoding::Symbol=:auto)
    text = hc_text(hc; lc_filename=lc_filename, encoding=encoding)
    open(path, "w") do io
        write(io, text)
    end
    path
end

export LCMoveError, LCMoveWitness, LCMoveResult,
       load_heegaard_moves, heegaard_move_row, heegaard_move_ids,
       movable, move, lc_text, hc_text, write_lc, write_hc
