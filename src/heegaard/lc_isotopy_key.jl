# Isotopy representative of one curve against the fixed K curves.
# Moves are applied through movable/move. This file does not name a diagram type
# above LickorishCode.

const _LC_ISO_DIGON_ROWS = (
    "lc.isotopy.digon.remove.00",
    "lc.isotopy.digon.remove.01",
    "lc.isotopy.digon.remove.10",
    "lc.isotopy.digon.remove.11",
)
const _LC_ISO_VERTEX_ROWS = (
    "lc.isotopy.vertex.000",
    "lc.isotopy.vertex.001",
    "lc.isotopy.vertex.010",
    "lc.isotopy.vertex.011",
    "lc.isotopy.vertex.100",
    "lc.isotopy.vertex.101",
    "lc.isotopy.vertex.110",
    "lc.isotopy.vertex.111",
)
const _LC_ISO_DIGON_CAP = 4000
const _LC_ISO_WORD_CAP = 800
const _LC_ISO_ORBIT_CAP = 4000
const _LC_ISO_CANON_CAP = 20000

"""A search stopped before a key existed. `kind` is `:algebra` or `:cap`."""
struct LCIsotopyKeyFailure
    kind::Symbol
    row::String
    cap::Symbol
    before::Vector{Tuple{String,String,Int}}
    after::Vector{Tuple{String,String,Int}}
end

"""Representative of one curve in minimal position with K.

`key` is `nothing` when `failure` is set. `words` is empty in that case.
"""
struct LCIsotopyKey
    curve::String
    names::Vector{String}
    words::Vector{Vector{Int}}
    key::Union{Nothing,String}
    orbit::Int
    digons_after_isolation::Int
    digons_after_vertex::Int
    algebra_ok::Bool
    failure::Union{Nothing,LCIsotopyKeyFailure}
end

function _lc_iso_owners(lc::LickorishCode)
    plus = Dict{Int,String}()
    minus = Dict{Int,String}()
    for (name, word) in lc.words
        for x in word
            if x > 0
                plus[abs(x)] = String(name)
            else
                minus[abs(x)] = String(name)
            end
        end
    end
    plus, minus
end

function _lc_iso_algebraic(lc::LickorishCode, a::AbstractString, b::AbstractString)
    plus, minus = _lc_iso_owners(lc)
    n = 0
    for lab in keys(plus)
        pair = (plus[lab], minus[lab])
        if pair == (String(a), String(b))
            n += 1
        elseif pair == (String(b), String(a))
            n -= 1
        end
    end
    n
end

function _lc_iso_k_algebra(lc::LickorishCode)
    kn = String.(_lc_k_names(lc.genus))
    rows = Tuple{String,String,Int}[]
    for a in sort!(collect(keys(lc.words))), b in kn
        String(a) == b && continue
        push!(rows, (String(a), b, _lc_iso_algebraic(lc, String(a), b)))
    end
    rows
end

function _lc_iso_blank(genus::Integer, words::Dict{String,Vector{Int}})
    LickorishCode(
        Int(genus),
        words,
        Dict{String,String}(),
        Dict{String,Vector{String}}(),
        Dict{String,String}(),
        Dict{String,Int}(),
        nothing,
    )
end

function _lc_iso_keep_words(lc::LickorishCode, keep::Set{String})
    plus, minus = _lc_iso_owners(lc)
    words = Dict{String,Vector{Int}}()
    for name in keep
        words[name] = Int[
            x for x in lc.words[name] if
            plus[abs(x)] in keep && minus[abs(x)] in keep
        ]
    end
    words
end

function _lc_iso_isolate(lc::LickorishCode, curve::AbstractString)
    kn = Set(String.(_lc_k_names(lc.genus)))
    keep = push!(copy(kn), String(curve))
    _lc_iso_blank(lc.genus, _lc_iso_keep_words(lc, keep))
end

function _lc_iso_maxlen(lc::LickorishCode)
    isempty(lc.words) ? 0 : maximum(length, values(lc.words))
end

function _lc_iso_one_digon(lc::LickorishCode)
    for row in _LC_ISO_DIGON_ROWS
        hits = movable(lc, row)
        isempty(hits) && continue
        return row, first(hits)
    end
    return nothing
end

function _lc_iso_algebra_failure(row, before, after)
    LCIsotopyKeyFailure(:algebra, String(row), :none, before, after)
end

function _lc_iso_cap_failure(cap::Symbol, row, before, after)
    LCIsotopyKeyFailure(:cap, String(row), cap, before, after)
end

function _lc_iso_reduce_digons(lc::LickorishCode)
    state = lc
    removed = String[]
    while true
        found = _lc_iso_one_digon(state)
        found === nothing && return state, removed, nothing
        row, hit = found
        before = _lc_iso_k_algebra(state)
        state = move(state, hit).output
        after = _lc_iso_k_algebra(state)
        if before != after
            return state, removed, _lc_iso_algebra_failure(row, before, after)
        end
        push!(removed, row)
        if length(removed) > _LC_ISO_DIGON_CAP
            return state, removed, _lc_iso_cap_failure(:digons, row, before, after)
        end
        if _lc_iso_maxlen(state) > _LC_ISO_WORD_CAP
            return state, removed, _lc_iso_cap_failure(:word_length, row, before, after)
        end
    end
end

function _lc_iso_rotations(word::Vector{Int})
    n = length(word)
    n == 0 && return [Int[]]
    Vector{Int}[
        Int[word[mod1(i + d, n)] for d in 0:(n - 1)]
        for i in 1:n
    ]
end

function _lc_iso_product(dims, cap::Int)
    p = 1
    for d in dims
        d < 1 && return cap + 1
        p > cap ÷ d && return cap + 1
        p *= d
    end
    p
end

function _lc_iso_canon_words(lc::LickorishCode)
    kn = String.(_lc_k_names(lc.genus))
    limit = _k_limit(lc.genus)
    names = String[]
    for name in kn
        haskey(lc.words, name) && push!(names, name)
    end
    for name in sort!(collect(setdiff(keys(lc.words), Set(kn))))
        push!(names, String(name))
    end
    rots = [_lc_iso_rotations(Int.(lc.words[name])) for name in names]
    dims = [length(r) for r in rots]
    prodlen = isempty(dims) ? 1 : _lc_iso_product(dims, _LC_ISO_CANON_CAP)
    if prodlen > _LC_ISO_CANON_CAP
        return names, Vector{Vector{Int}}(), nothing
    end
    idx = ones(Int, length(dims))
    best = nothing
    for _ in 1:prodlen
        choice = [rots[i][idx[i]] for i in eachindex(idx)]
        next_id = Ref(limit)
        maplab = Dict{Int,Int}()
        function mapx(x)
            a = abs(x)
            a <= limit && return sign(x) * a
            if !haskey(maplab, a)
                next_id[] += 1
                maplab[a] = next_id[]
            end
            sign(x) * maplab[a]
        end
        key = [Int[mapx(x) for x in word] for word in choice]
        (best === nothing || key < best) && (best = key)
        j = length(idx)
        while j >= 1
            idx[j] += 1
            if idx[j] <= dims[j]
                break
            end
            idx[j] = 1
            j -= 1
        end
    end
    names, best, true
end

function _lc_iso_form_id(names, words)
    join((names[i] * ":" * join(words[i], ",") for i in eachindex(names)), "|")
end

function _lc_iso_shorter(words, best)
    sw, sb = sum(length, words), sum(length, best)
    sw < sb || (sw == sb && words < best)
end

function _lc_iso_vertex_orbit(lc::LickorishCode)
    names, words, canon_ok = _lc_iso_canon_words(lc)
    if canon_ok === nothing
        return (
            names = String[],
            words = Vector{Vector{Int}}(),
            key = nothing,
            orbit = 0,
            digons_after_vertex = 0,
            failure = _lc_iso_cap_failure(:canonical, "", _lc_iso_k_algebra(lc), _lc_iso_k_algebra(lc)),
        )
    end
    start_id = _lc_iso_form_id(names, words)
    best_words = words
    best_names = names
    seen = Set{String}([start_id])
    queue = [lc]
    head = 1
    exposed = 0
    while head <= length(queue)
        cur = queue[head]
        head += 1
        for row in _LC_ISO_VERTEX_ROWS
            hits = movable(cur, row)
            for hit in hits
                before = _lc_iso_k_algebra(cur)
                nxt = move(cur, hit).output
                mid = _lc_iso_k_algebra(nxt)
                if before != mid
                    return (
                        names = String[],
                        words = Vector{Vector{Int}}(),
                        key = nothing,
                        orbit = length(seen),
                        digons_after_vertex = exposed,
                        failure = _lc_iso_algebra_failure(row, before, mid),
                    )
                end
                nxt, removed, info = _lc_iso_reduce_digons(nxt)
                exposed += length(removed)
                if info !== nothing
                    return (
                        names = String[],
                        words = Vector{Vector{Int}}(),
                        key = nothing,
                        orbit = length(seen),
                        digons_after_vertex = exposed,
                        failure = info,
                    )
                end
                nn, ww, ok2 = _lc_iso_canon_words(nxt)
                if ok2 === nothing
                    return (
                        names = String[],
                        words = Vector{Vector{Int}}(),
                        key = nothing,
                        orbit = length(seen),
                        digons_after_vertex = exposed,
                        failure = _lc_iso_cap_failure(:canonical, row, before, mid),
                    )
                end
                id = _lc_iso_form_id(nn, ww)
                if !(id in seen)
                    push!(seen, id)
                    if _lc_iso_shorter(ww, best_words)
                        best_words = ww
                        best_names = nn
                    end
                    push!(queue, nxt)
                end
                if length(seen) > _LC_ISO_ORBIT_CAP
                    return (
                        names = String[],
                        words = Vector{Vector{Int}}(),
                        key = nothing,
                        orbit = length(seen),
                        digons_after_vertex = exposed,
                        failure = _lc_iso_cap_failure(:orbit, row, before, mid),
                    )
                end
            end
        end
    end
    (
        names = best_names,
        words = best_words,
        key = _lc_iso_form_id(best_names, best_words),
        orbit = length(seen),
        digons_after_vertex = exposed,
        failure = nothing,
    )
end

function _lc_iso_stopped(curve, isolated_digons, orbit, vertex_digons, failure)
    LCIsotopyKey(
        String(curve),
        String[],
        Vector{Vector{Int}}(),
        nothing,
        orbit,
        isolated_digons,
        vertex_digons,
        failure.kind != :algebra,
        failure,
    )
end

"""
    lc_isotopy_key(lc::LickorishCode, curve::AbstractString) -> LCIsotopyKey
    lc_isotopy_key(lc::LickorishCode, curves) -> Vector{LCIsotopyKey}

Representative of `curve` in minimal position with the K curves of `lc`.

The curve is separated together with K. Digons are removed by the four
`lc.isotopy.digon.remove` rows, taking the first witness of the first row
that applies. The vertex orbit uses `lc.isotopy.vertex.000` through `.111`.
Each exchange is followed by the same digon removal. The representative is a
shortest word, then the lexicographically least cyclic relabeling. Labels of
K-K crossings stay fixed. A word with reversed orientation is a different key.

Algebraic intersections with K are compared around every removal and every
vertex exchange. When one changes, `key` is `nothing` and `failure` holds the
row and the intersections before and after. The search also stops without a
key after 4000 digon removals, when a word exceeds length 800, when the
vertex orbit exceeds 4000 states, or when the product of the cyclic shifts
exceeds 20000. The input `lc` is left unchanged.
"""
function lc_isotopy_key(lc::LickorishCode, curve::AbstractString)
    name = String(curve)
    haskey(lc.words, name) || throw(KodaInputError(
        :lc_isotopy_key, "missing curve", Dict("curve" => name),
    ))
    iso = _lc_iso_isolate(lc, name)
    validate_lc_input(iso)
    reduced, steps, info = _lc_iso_reduce_digons(iso)
    if info !== nothing
        return _lc_iso_stopped(name, length(steps), 0, 0, info)
    end
    rep = _lc_iso_vertex_orbit(reduced)
    if rep.failure !== nothing
        return _lc_iso_stopped(name, length(steps), rep.orbit, rep.digons_after_vertex, rep.failure)
    end
    LCIsotopyKey(
        name,
        rep.names,
        rep.words,
        rep.key,
        rep.orbit,
        length(steps),
        rep.digons_after_vertex,
        true,
        nothing,
    )
end

function lc_isotopy_key(lc::LickorishCode, curves::AbstractVector{<:AbstractString})
    LCIsotopyKey[lc_isotopy_key(lc, curve) for curve in curves]
end
