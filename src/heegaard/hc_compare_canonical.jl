# The stored start of a word is not a tie-break. Rotations that emit the same
# sequence can still assign different labels; later curves select the lex-min map.
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3 — a cyclic shift is not data

"""
    HCCanonicalWitness(label_map, name_map, rotation_starts) -> HCCanonicalWitness

Store a correspondence from the labels and words of an HC to its canonical representation.

# Arguments
`label_map::Dict{Int,Int}`, `name_map::Dict{String,String}`,
`rotation_starts::Dict{String,Int}`, or values convertible to those types.
Conversion errors propagate. No source diagram is supplied or checked.

# Returns
`HCCanonicalWitness` with the supplied label correspondence, component-name
correspondence and word-start data. Obtain computed witnesses from
[`hc_canonical_witness`](@ref). Dictionary domains and the correspondence
direction are specified in [comparison records](../specs/heegaard-code.md#comparison-records).

# Guarantees
1. Fields equal the supplied values after conversion; typed dictionaries
   are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. For a computed witness, its data describe the correspondence to
   `hc_canonical` for the same input. The cyclic-word objects are those of
   [the LC specification](../specs/lickorish-code.md). [proved/agent: by construction]

# Not guaranteed
The constructor does not check completeness, bijectivity or compatibility
with words and annotations. A record alone does not prove identity, isotopy
or manifold equivalence and does not carry its source diagram.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-13", "case-01.hc");

julia> h = H.read_hc(path); w = H.hc_canonical_witness(h);

julia> copied = H.HCCanonicalWitness(w.label_map, w.name_map, w.rotation_starts);

julia> corr = H.HCIdentityCorrespondence(copied.label_map, copied.name_map, copied.rotation_starts);

julia> H.hc_identical(h, H.hc_apply_correspondence(h, corr))
true
```

See also [`hc_canonical_witness`](@ref), [`HCIdentityCorrespondence`](@ref).

Role: construct
"""
struct HCCanonicalWitness
    label_map::Dict{Int,Int}
    name_map::Dict{String,String}
    rotation_starts::Dict{String,Int}
end

"""
    HCIdentityCorrespondence(label_map, name_map, rotation_starts) -> HCIdentityCorrespondence

Store label, component-name and word-start data for a correspondence between two HCs.

# Arguments
`label_map::Dict{Int,Int}`, `name_map::Dict{String,String}`,
`rotation_starts::Dict{String,Int}`, or values convertible to those types.
Conversion errors propagate. Neither source nor target is an argument.

# Returns
`HCIdentityCorrespondence` with the supplied three dictionaries, for use
with [`hc_apply_correspondence`](@ref). Computed records are obtained from
[`hc_identity_correspondence`](@ref). Dictionary domains and the correspondence
direction are specified in [comparison records](../specs/heegaard-code.md#comparison-records).

# Guarantees
1. Fields preserve supplied values after conversion; typed dictionaries
   are retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. A computed correspondence describes the identity relation tested by
   `hc_identity_correspondence` on the two inputs' stored words and
   annotations. It acts on the cyclic-word data of
   [the LC specification](../specs/lickorish-code.md). [proved/agent: by construction]

# Not guaranteed
The constructor does not check map domains, bijectivity, word positions,
annotation transport or preservation of curve families. A manually
assembled record is not a checked identity, isotopy or manifold certificate.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-13", "case-01.hc");

julia> h = H.read_hc(path); c = H.hc_identity_correspondence(h, h);

julia> d = H.HCIdentityCorrespondence(c.label_map, c.name_map, c.rotation_starts);

julia> (H.hc_apply_correspondence(h, d).lc.words == h.lc.words, H.hc_identical(H.hc_apply_correspondence(h, d), h))
(true, true)
```

See also [`hc_identity_correspondence`](@ref), [`hc_apply_correspondence`](@ref).

Role: construct
"""
struct HCIdentityCorrespondence
    label_map::Dict{Int,Int}
    name_map::Dict{String,String}
    rotation_starts::Dict{String,Int}
end

const _HC_WITNESS_CAP = 256

function _hc_form_less(a::HCCanonical, b::HCCanonical)
    for i in 1:fieldcount(HCCanonical)
        left = getfield(a, i)
        right = getfield(b, i)
        left < right && return true
        right < left && return false
    end
    false
end

function _hc_rotate_word(word::Vector{Int}, start::Int)
    n = length(word)
    n == 0 && return Int[]
    Int[word[mod1(start + i - 1, n)] for i in 1:n]
end

function _hc_map_key(label_map::Dict{Int,Int})
    Tuple(sort!(collect(pairs(label_map))))
end

function _hc_future_labels(names, words)
    seen = Set{Int}()
    labels = Int[]
    for name in names
        haskey(words, name) || continue
        for x in words[name]
            lab = abs(x)
            lab in seen && continue
            push!(seen, lab)
            push!(labels, lab)
        end
    end
    sort!(labels)
end

function _hc_min_rotations(word::Vector{Int}, label_map::Dict{Int,Int}, next_id::Int, future::Vector{Int})
    n = length(word)
    n == 0 && return [(Int[], copy(label_map), next_id, 1)]
    best = nothing
    rows = Tuple{Vector{Int},Dict{Int,Int},Int,Int}[]
    for start in 1:n
        emitted, local_map, nid = _hc_emit_rotation(word, start, label_map, next_id)
        if best === nothing || emitted < best
            best = emitted
        end
        push!(rows, (emitted, local_map, nid, start))
    end
    seen = Dict{Tuple,Tuple{Vector{Int},Dict{Int,Int},Int,Int}}()
    for (emitted, local_map, nid, start) in rows
        emitted == best || continue
        key = (nid, Tuple(lab => local_map[lab] for lab in future if haskey(local_map, lab)))
        if !haskey(seen, key) || start < seen[key][4]
            seen[key] = (emitted, local_map, nid, start)
        end
    end
    chosen = collect(values(seen))
    sort!(chosen, by=row -> row[4])
    chosen
end

function _hc_guard_tie_groups(groups)
    complexity = 1
    for group in groups
        n = length(group)
        n <= 1 && continue
        n > 8 && throw(HCNotComputed(
            :canonical_tie_group,
            Dict{String,Any}("n" => n, "limit" => 8),
        ))
        complexity *= factorial(n)
        if complexity > _HC_TIE_LIMIT
            throw(HCNotComputed(
                :canonical_tie_group,
                Dict{String,Any}("complexity" => complexity, "limit" => _HC_TIE_LIMIT),
            ))
        end
    end
    nothing
end

function _hc_witness_rank(witness::HCCanonicalWitness)
    (
        sort!(collect(pairs(witness.rotation_starts))),
        sort!(collect(pairs(witness.label_map))),
        sort!(collect(pairs(witness.name_map))),
    )
end

function _hc_build_form(
    hc::HeegaardCode,
    k_order, alpha_order, beta_order, other_order,
    k_words, alpha_words, beta_words, other_emitted,
    alpha_empty, beta_empty, other_empty,
    label_map, starts,
)
    owners = isempty(hc.lc.words) ? Dict{Int,Tuple{String,String}}() : _hc_owners(hc.lc.words)
    for lab in keys(owners)
        haskey(label_map, lab) || throw(KodaImplementationError(
            :hc_canonical, "label was not assigned", Dict("label" => lab)))
    end
    face_map = _hc_has_annotations(hc) && !isempty(label_map) ?
        _hc_face_map(hc.lc.words, label_map) : Dict{String,String}()
    name_map = _hc_name_map(k_order, alpha_order, beta_order, other_order)
    if _hc_has_annotations(hc)
        alpha_empty_order = _hc_best_empty_order(alpha_empty, hc, name_map, face_map, label_map, "a")
        beta_empty_order = _hc_best_empty_order(beta_empty, hc, name_map, face_map, label_map, "b")
        other_empty_order = _hc_best_empty_order(other_empty, hc, name_map, face_map, label_map, "o")
    else
        alpha_empty_order = alpha_empty
        beta_empty_order = beta_empty
        other_empty_order = other_empty
    end
    alpha_order = vcat(alpha_order, alpha_empty_order)
    beta_order = vcat(beta_order, beta_empty_order)
    other_order = vcat(other_order, other_empty_order)
    alpha_words = vcat(alpha_words, [Int[] for _ in alpha_empty_order])
    beta_words = vcat(beta_words, [Int[] for _ in beta_empty_order])
    other_emitted = vcat(other_emitted, [Int[] for _ in other_empty_order])
    full_starts = copy(starts)
    for name in vcat(alpha_empty_order, beta_empty_order, other_empty_order)
        full_starts[name] = 1
    end
    name_map = _hc_name_map(k_order, alpha_order, beta_order, other_order)
    for name in keys(hc.lc.words)
        sname = String(name)
        haskey(name_map, sname) || throw(KodaImplementationError(
            :hc_canonical, "curve was not named", Dict("curve" => sname)))
        haskey(full_starts, sname) || throw(KodaImplementationError(
            :hc_canonical, "curve has no rotation start", Dict("curve" => sname)))
    end
    places, lc_disks, hc_disks, lc_points, hc_points, lc_weights, hc_weights =
        _hc_annotation(hc, name_map, face_map, label_map)
    form = HCCanonical(
        hc.lc.genus,
        [name => word for (name, word) in zip(k_order, k_words)],
        alpha_words,
        beta_words,
        [name => word for (name, word) in zip(other_order, other_emitted)],
        places, lc_disks, hc_disks, lc_points, hc_points, lc_weights, hc_weights,
    )
    form, HCCanonicalWitness(copy(label_map), name_map, full_starts)
end

function _hc_remember!(best_form, best_witnesses, form, witness)
    if best_form[] === nothing || _hc_form_less(form, best_form[])
        best_form[] = form
        empty!(best_witnesses)
        push!(best_witnesses, witness)
        return
    end
    form == best_form[] || return
    rank = _hc_witness_rank(witness)
    any(existing -> _hc_witness_rank(existing) == rank, best_witnesses) && return
    if length(best_witnesses) >= _HC_WITNESS_CAP
        throw(HCNotComputed(
            :canonical_search_budget,
            Dict{String,Any}("budget" => _HC_WITNESS_CAP, "reason" => "witness_cap"),
        ))
    end
    push!(best_witnesses, witness)
end

function _hc_canonicalize(hc::HeegaardCode)
    words = hc.lc.words
    isempty(words) || _hc_owners(words)
    alphas, betas, others = _hc_curve_lists(hc)
    comp_colors, _ = _hc_colors(words, hc.roles, hc.lc.genus)
    alpha_present, alpha_empty = _hc_split_empty(alphas, words)
    beta_present, beta_empty = _hc_split_empty(betas, words)
    other_present, other_empty = _hc_split_empty(others, words)
    k_present = String[]
    for name in _lc_k_names(hc.lc.genus)
        haskey(words, name) || continue
        push!(k_present, name)
    end
    groups = Vector{String}[]
    for name in k_present
        push!(groups, String[name])
    end
    for group in _hc_tie_groups(alpha_present, comp_colors)
        push!(groups, group.second)
    end
    for group in _hc_tie_groups(beta_present, comp_colors)
        push!(groups, group.second)
    end
    for group in _hc_tie_groups(other_present, comp_colors)
        push!(groups, group.second)
    end
    _hc_guard_tie_groups(groups)
    n_k = length(k_present)
    n_alpha = length(alpha_present)
    n_beta = length(beta_present)
    n_fixed = n_k + n_alpha + n_beta
    annotated = _hc_has_annotations(hc)
    all_labels = annotated ? _hc_future_labels(collect(String.(keys(words))), words) : Int[]
    best_form = Ref{Union{Nothing,HCCanonical}}(nothing)
    best_code = Ref{Union{Nothing,Vector{Vector{Int}}}}(nothing)
    best_witnesses = HCCanonicalWitness[]
    nodes = Ref(0)
    function note!()
        nodes[] += 1
        nodes[] > _HC_TIE_LIMIT && throw(HCNotComputed(
            :canonical_search_budget,
            Dict{String,Any}("budget" => _HC_TIE_LIMIT),
        ))
    end
    function worse_prefix(partial, complete)
        m = min(length(partial), length(complete))
        for i in 1:m
            partial[i] < complete[i] && return false
            partial[i] > complete[i] && return true
        end
        false
    end
    function finish(order, emitted, label_map, starts)
        note!()
        form, witness = _hc_build_form(
            hc,
            order[1:n_k],
            order[n_k+1:n_k+n_alpha],
            order[n_k+n_alpha+1:n_k+n_alpha+n_beta],
            order[n_k+n_alpha+n_beta+1:end],
            emitted[1:n_k],
            emitted[n_k+1:n_k+n_alpha],
            emitted[n_k+n_alpha+1:n_k+n_alpha+n_beta],
            emitted[n_k+n_alpha+n_beta+1:end],
            alpha_empty, beta_empty, other_empty, label_map, starts,
        )
        before = best_form[]
        _hc_remember!(best_form, best_witnesses, form, witness)
        if best_form[] !== before
            best_code[] = emitted
        end
    end
    # α と β は、出した語の列が正規形の比較順そのものである。
    # 今の位置で語がすでに大きい曲線を先に置く並べ方は、辞書式最小にはならない。
    word_ordered = falses(length(groups))
    for (i, members) in enumerate(groups)
        length(members) > 1 || continue
        role = get(hc.roles, members[1], nothing)
        (role == :alpha || role == :beta) || continue
        word_ordered[i] = all(name -> get(hc.roles, name, nothing) == role, members)
    end
    function future_of(index, remaining, name)
        annotated && return all_labels
        later = String[]
        for other in remaining
            other == name && continue
            push!(later, other)
        end
        for group in groups[index+1:end]
            append!(later, group)
        end
        _hc_future_labels(later, words)
    end
    function take(index, remaining, name, rots, map, nid, code, order, starts)
        next_remaining = String[]
        for other in remaining
            other == name && continue
            push!(next_remaining, other)
        end
        for (emitted, local_map, nn, start) in rots
            note!()
            new_code = vcat(code, [emitted])
            if length(new_code) <= n_fixed && best_code[] !== nothing && worse_prefix(new_code, best_code[])
                continue
            end
            new_starts = copy(starts)
            new_starts[name] = start
            descend(index, next_remaining, local_map, nn, new_code, vcat(order, [name]), new_starts)
        end
    end
    function descend(index, remaining, map, nid, code, order, starts)
        if isempty(remaining)
            rec(index + 1, map, nid, code, order, starts)
            return
        end
        if word_ordered[index]
            best_word = nothing
            selected = Any[]
            for name in remaining
                rots = _hc_min_rotations(words[name], map, nid, future_of(index, remaining, name))
                word = rots[1][1]
                if best_word === nothing || word < best_word
                    best_word = word
                    empty!(selected)
                    push!(selected, (name, rots))
                elseif word == best_word
                    push!(selected, (name, rots))
                end
            end
            if best_code[] !== nothing && length(code) + 1 <= n_fixed &&
                    worse_prefix(vcat(code, [best_word]), best_code[])
                return
            end
            for (name, rots) in selected
                take(index, remaining, name, rots, map, nid, code, order, starts)
            end
        else
            for name in remaining
                rots = _hc_min_rotations(words[name], map, nid, future_of(index, remaining, name))
                take(index, remaining, name, rots, map, nid, code, order, starts)
            end
        end
    end
    function rec(index, map, nid, code, order, starts)
        if index > length(groups)
            finish(order, code, map, starts)
            return
        end
        descend(index, groups[index], map, nid, code, order, starts)
    end
    rec(1, Dict{Int,Int}(), 1, Vector{Vector{Int}}(), String[], Dict{String,Int}())
    best_form[] === nothing && throw(KodaImplementationError(
        :hc_canonical, "search produced no form", Dict{String,Any}()))
    sort!(best_witnesses, by=_hc_witness_rank)
    best_form[], best_witnesses
end

"""
    hc_canonical(hc::HeegaardCode) -> HCCanonical

Return the canonical comparison record of an HC's oriented stored data.

# Arguments
`hc` must satisfy the [HC input and comparison domain](../specs/heegaard-code.md#comparison-records).
Input and annotation errors propagate. A bounded search that cannot finish
raises `HCNotComputed` with its reason rather than returning a difference.

# Returns
`HCCanonical` with genus, reference words, alpha and beta words, other
components and transported annotations in the comparison format.

# Guarantees
1. The returned record uses the identity relation and ordering of the linked
   comparison convention. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is not a canonical form for 3-manifolds or an isotopy decision, and does
not identify general orientation reversals or mirrors. Search completion
is not guaranteed on every input.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — canonical comparison records and ordering.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta, "m1"=>:alpha));

julia> c = H.hc_canonical(hc); (c.genus, Tuple(p.first for p in c.k_words))
(1, ("l1", "m1"))
```

Role: invariant
"""
function hc_canonical(hc::HeegaardCode)
    form, _ = _hc_canonicalize(hc)
    form
end

"""
    hc_canonical_witness(hc::HeegaardCode) -> HCCanonicalWitness

Return a correspondence from an HC's stored data to its canonical comparison record.

# Arguments
`hc` has the domain and failure conditions of `hc_canonical` in
[comparison records](../specs/heegaard-code.md#comparison-records).
Invalid input and annotation errors propagate; bounded search exhaustion
raises `HCNotComputed`.

# Returns
`HCCanonicalWitness` with `label_map`, `name_map` and `rotation_starts` in
the linked correspondence convention.

# Guarantees
1. The witness carries the input's data to its computed canonical record
   under the stated dictionary interpretation. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The witness need not preserve fixed K names when roles are assigned to K
curves; it is not by itself a valid resolved-HC construction or a manifold
homeomorphism certificate. There may be several witnessing symmetries.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — witness direction, ordering and K-name limitation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta, "m1"=>:alpha));

julia> w = H.hc_canonical_witness(hc); (w.name_map["l1"], w.name_map["m1"])
("b1", "a1")
```

Role: access
"""
function hc_canonical_witness(hc::HeegaardCode)
    _, witnesses = _hc_canonicalize(hc)
    witnesses[1]
end

function _hc_transport_disk(dict, transport)
    out = Dict{String,Vector{String}}()
    for (name, refs) in dict
        out[String(name)] = String[transport(ref) for ref in refs]
    end
    out
end

function _hc_transport_point(dict, transport)
    Dict{String,String}(String(name) => transport(ref) for (name, ref) in dict)
end

function _hc_transport_weight(dict, transport)
    Dict{String,Int}(transport(String(name)) => Int(value) for (name, value) in dict)
end

"""
    hc_apply_correspondence(hc::HeegaardCode, corr::HCIdentityCorrespondence) -> HeegaardCode

Transport an HC's words, roles and annotations along a supplied identity correspondence.

# Arguments
`corr` must obey the [correspondence dictionary convention](../specs/heegaard-code.md#comparison-records).
Missing curve or crossing mappings raise `KodaInputError` identifying the
missing item; absent rotation keys and annotation-transport failures propagate.
The caller supplies complete bijective, role-compatible maps; this entry does
not independently validate a manually assembled correspondence.

# Returns
A new `HeegaardCode` with transported words, roles, places, disks, points and
weights; both source paths are cleared.

# Guarantees
1. For a valid correspondence, the output represents the same oriented
   cyclic-word and annotation data in the target notation.
   [proved/agent: by construction]
2. The input and correspondence are not modified. [proved/agent: by construction]

# Not guaranteed
Arbitrary dictionaries can collide or destroy information without a dedicated
bijection error. A canonical witness is not automatically a valid resolved-HC
identity map; use `hc_identity_correspondence` between stored HCs.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — dictionary direction and canonical-versus-identity distinction.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta, "m1"=>:alpha));

julia> c = H.hc_identity_correspondence(hc,hc); out = H.hc_apply_correspondence(hc,c);

julia> (out.lc.words == hc.lc.words, out.roles == hc.roles)
(true, true)
```

Role: convert
"""
function hc_apply_correspondence(hc::HeegaardCode, corr::HCIdentityCorrespondence)
    words = Dict{String,Vector{Int}}()
    for (name, word) in hc.lc.words
        sname = String(name)
        haskey(corr.name_map, sname) || throw(KodaInputError(
            :hc_apply_correspondence, "missing curve", Dict("curve" => sname)))
        rotated = _hc_rotate_word(word, corr.rotation_starts[sname])
        mapped = Int[]
        for x in rotated
            lab = abs(x)
            haskey(corr.label_map, lab) || throw(KodaInputError(
                :hc_apply_correspondence, "missing label", Dict("label" => lab)))
            push!(mapped, sign(x) * corr.label_map[lab])
        end
        words[corr.name_map[sname]] = mapped
    end
    roles = Dict{String,Symbol}()
    for (name, role) in hc.roles
        sname = String(name)
        haskey(corr.name_map, sname) || throw(KodaInputError(
            :hc_apply_correspondence, "missing role curve", Dict("curve" => sname)))
        roles[corr.name_map[sname]] = role
    end
    face_map = _hc_has_annotations(hc) && !isempty(corr.label_map) ?
        _hc_face_map(hc.lc.words, corr.label_map) : Dict{String,String}()
    transport = ref -> _hc_transport_ref(ref, corr.name_map, face_map, corr.label_map)
    places = Dict{String,String}()
    for (name, ref) in hc.lc.places
        sname = String(name)
        haskey(corr.name_map, sname) || throw(KodaInputError(
            :hc_apply_correspondence, "missing place curve", Dict("curve" => sname)))
        places[corr.name_map[sname]] = transport(ref)
    end
    lc = LickorishCode(
        hc.lc.genus, words, places,
        _hc_transport_disk(hc.lc.disks, transport),
        _hc_transport_point(hc.lc.points, transport),
        _hc_transport_weight(hc.lc.weights, transport),
        nothing,
    )
    HeegaardCode(
        lc, roles,
        _hc_transport_disk(hc.disks, transport),
        _hc_transport_point(hc.points, transport),
        _hc_transport_weight(hc.weights, transport),
        nothing,
    )
end

function _hc_compose_witness(a::HeegaardCode, left::HCCanonicalWitness, right::HCCanonicalWitness)
    inv_label = Dict{Int,Int}()
    for (old, canonical) in right.label_map
        haskey(inv_label, canonical) && throw(KodaImplementationError(
            :hc_identity_correspondence, "canonical label is not unique", Dict("label" => canonical)))
        inv_label[canonical] = old
    end
    inv_name = Dict{String,String}()
    for (old, canonical) in right.name_map
        haskey(inv_name, canonical) && throw(KodaImplementationError(
            :hc_identity_correspondence, "canonical curve name is not unique", Dict("curve" => canonical)))
        inv_name[canonical] = old
    end
    label_map = Dict{Int,Int}()
    for (old, canonical) in left.label_map
        haskey(inv_label, canonical) || throw(KodaImplementationError(
            :hc_identity_correspondence, "canonical label is missing on the target", Dict("label" => canonical)))
        label_map[old] = inv_label[canonical]
    end
    name_map = Dict{String,String}()
    starts = Dict{String,Int}()
    for name in keys(a.lc.words)
        sname = String(name)
        haskey(left.name_map, sname) || throw(KodaImplementationError(
            :hc_identity_correspondence, "source curve was not named", Dict("curve" => sname)))
        canonical = left.name_map[sname]
        haskey(inv_name, canonical) || throw(KodaImplementationError(
            :hc_identity_correspondence, "canonical curve is missing on the target", Dict("curve" => canonical)))
        target = inv_name[canonical]
        name_map[sname] = target
        n = length(a.lc.words[name])
        if n == 0
            starts[sname] = 1
        else
            starts[sname] = mod1(left.rotation_starts[sname] - right.rotation_starts[target] + 1, n)
        end
    end
    HCIdentityCorrespondence(label_map, name_map, starts)
end

function _hc_stored_words(hc::HeegaardCode)
    Dict{String,Vector{Int}}(String(name) => copy(word) for (name, word) in hc.lc.words)
end

"""
    hc_identity_correspondence(a::HeegaardCode, b::HeegaardCode) -> HCIdentityCorrespondence

Return a checked correspondence between two HCs with the same canonical comparison record.

# Arguments
The inputs have the domain of `hc_canonical` under
[comparison records](../specs/heegaard-code.md#comparison-records).
Different forms raise `HCNotComputed(:not_identical)`; equal forms without a
checked transporting correspondence raise `HCNotComputed(:correspondence_inconsistent)`.
Input and bounded-search failures propagate.

# Returns
`HCIdentityCorrespondence` from `a` to `b`, with label, name and rotation-start
dictionaries in the linked convention.

# Guarantees
1. Applying the result to `a` gives exactly `b`'s stored curve words and
   the same canonical comparison record as `b`.
   [proved/agent: by construction]
2. Neither input is modified. [proved/agent: by construction]

# Not guaranteed
This is a correspondence of the encoded diagram data, not a general manifold
homeomorphism or isotopy search. Exact identity of source paths or arbitrary
metadata is not asserted.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#comparison-records — identity correspondences and checked transport.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> c = H.hc_identity_correspondence(hc,hc); H.hc_apply_correspondence(hc,c).lc.words == hc.lc.words
true
```

Role: access
"""
function hc_identity_correspondence(a::HeegaardCode, b::HeegaardCode)
    form_a, witnesses_a = _hc_canonicalize(a)
    form_b, witnesses_b = _hc_canonicalize(b)
    form_a == form_b || throw(HCNotComputed(
        :not_identical,
        Dict{String,Any}("reason" => "canonical forms differ"),
    ))
    target_words = _hc_stored_words(b)
    for left in witnesses_a, right in witnesses_b
        corr = _hc_compose_witness(a, left, right)
        applied = hc_apply_correspondence(a, corr)
        _hc_stored_words(applied) == target_words || continue
        hc_canonical(applied) == form_b || continue
        return corr
    end
    throw(HCNotComputed(
        :correspondence_inconsistent,
        Dict{String,Any}("witnesses" => length(witnesses_a) * length(witnesses_b)),
    ))
end
