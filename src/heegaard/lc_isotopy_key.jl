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

"""
    LCIsotopyKey(curve, names, words, key, orbit, digons_after_isolation, digons_after_vertex, algebra_ok, failure) -> LCIsotopyKey

Store a curve representative and diagnostics from a bounded isotopy-key search.

# Arguments
`curve::String`, `names::Vector{String}`, `words::Vector{Vector{Int}}`,
`key::Union{Nothing,String}`, three integer search counts, `algebra_ok::Bool`
and an optional `LCIsotopyKeyFailure`, or convertible values. Conversion
errors propagate. The constructor does not check consistency of a key,
word data, counts or failure record.

# Returns
An `LCIsotopyKey` retaining the curve name, representative words, optional
key, search counts, algebra-check flag and optional failure.

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers are
   shared and not modified. [proved/agent: by construction]
2. The record can carry a search result for the fixed-surface curve-isotopy
   problem in [the move specification](../specs/heegaard-moves.md).
   [proved/agent: by construction]

# Conventions
LC key normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/lickorish-code.md#lc-isotopy-key-normalization).

# Not guaranteed
Construction does not certify minimal position or isotopy. A stopped bounded
search does not prove that no equivalence exists. The [key specification](../specs/lickorish-code.md#lc-isotopy-key-normalization)
states the normalization and bounded-search scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> k = LCIsotopyKey("c", String[], Vector{Int}[], nothing, 0, 0, 0, false, nothing); (k.curve, k.key, k.algebra_ok)
("c", nothing, false)
```

See also [`lc_isotopy_key`](@ref).

Role: construct
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

Compute a bounded representative key for one oriented curve against the fixed LC reference curves.

# Arguments
`curve` names an existing component of `lc`; a missing name raises `KodaInputError` at `:lc_isotopy_key`. The curve together with the fixed reference curves must pass `validate_lc_input`. Move-domain and validation errors propagate. Search caps and detected changes of algebraic intersections are reported in the result rather than as a successful key.

# Returns
An `LCIsotopyKey` with representative names and words, search counts and diagnostics. `key === nothing` with a failure record means no key was certified by this run; a string key accompanies `algebra_ok=true` and `failure=nothing`.

# Guarantees
1. A successful record describes a representative reached by the supported local isotopy moves in the [LC move domain](../specs/heegaard-moves.md#6-各-move-のデータ). Its checked algebraic intersections with the reference curves agree with those of the isolated input. [proved/agent: by construction]
2. The input LC is not modified. [proved/agent: by construction]

# Conventions
LC key normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/lickorish-code.md#lc-isotopy-key-normalization).

# Not guaranteed
A successful key is not a proved complete invariant of surface isotopy, and a failed or unequal key does not prove nonisotopy. Global minimality over all geometric representatives and orientation-reversal invariance are not guaranteed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §6 — local isotopy moves on reference-marked curves.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-2,-1], "m1"=>[1,3], "a"=>[-4,-3], "b"=>[2,4]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> r = H.lc_isotopy_key(lc, "a"); (r.key isa String, r.algebra_ok, r.failure === nothing)
(true, true, true)
```

Role: invariant
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

"""
    lc_isotopy_key(lc::LickorishCode, curves::AbstractVector{<:AbstractString}) -> Vector{LCIsotopyKey}

Compute bounded representative keys for a list of oriented curves.

# Arguments
Each requested name has the single-curve `lc_isotopy_key` overload's preconditions and exceptions. A thrown failure interrupts the call; diagnostic result values are retained.

# Returns
One `LCIsotopyKey` per requested name, in the input order, including repeated names.

# Guarantees
1. Each element has the single-curve overload's representative and algebraic-intersection guarantees. [proved/agent: by construction]
2. The LC and the name vector are not modified. [proved/agent: by construction]

# Conventions
LC key normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/lickorish-code.md#lc-isotopy-key-normalization).

# Not guaranteed
This does not compute a simultaneous isotopy of the whole curve system or certify that independently chosen representatives are mutually compatible. The single-curve normalization and completeness limitations apply.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §6 — supported local moves for each selected curve.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-2,-1], "m1"=>[1,3], "a"=>[-4,-3], "b"=>[2,4]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> r = H.lc_isotopy_key(lc, ["a", "b"]); (length(r), all(x -> x.key isa String, r))
(2, true)
```

Role: invariant
"""
function lc_isotopy_key(lc::LickorishCode, curves::AbstractVector{<:AbstractString})
    LCIsotopyKey[lc_isotopy_key(lc, curve) for curve in curves]
end
