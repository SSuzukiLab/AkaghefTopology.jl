# The stored start of a word is not a tie-break. Rotations that emit the same
# sequence can still assign different labels; later curves select the lex-min map.
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3 — a cyclic shift is not data

struct HCCanonicalWitness
    label_map::Dict{Int,Int}
    name_map::Dict{String,String}
    rotation_starts::Dict{String,Int}
end

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
    hc::HeegaardLayer,
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

function _hc_canonicalize(hc::HeegaardLayer)
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
    function emit_group(index, perm_i, perm, map, nid, code, order, starts)
        if perm_i > length(perm)
            rec(index + 1, map, nid, code, order, starts)
            return
        end
        name = perm[perm_i]
        later = String[]
        perm_i < length(perm) && append!(later, perm[perm_i+1:end])
        for group in groups[index+1:end]
            append!(later, group)
        end
        future = annotated ? all_labels : _hc_future_labels(later, words)
        for (emitted, local_map, nn, start) in _hc_min_rotations(words[name], map, nid, future)
            note!()
            new_code = vcat(code, [emitted])
            if length(new_code) <= n_fixed && best_code[] !== nothing && worse_prefix(new_code, best_code[])
                continue
            end
            new_starts = copy(starts)
            new_starts[name] = start
            emit_group(
                index, perm_i + 1, perm, local_map, nn,
                new_code, vcat(order, [name]), new_starts,
            )
        end
    end
    function rec(index, map, nid, code, order, starts)
        if index > length(groups)
            finish(order, code, map, starts)
            return
        end
        members = groups[index]
        perms = length(members) == 1 ? [members] : _hc_permutations(members)
        for perm in perms
            emit_group(index, 1, perm, map, nid, code, order, starts)
        end
    end
    rec(1, Dict{Int,Int}(), 1, Vector{Vector{Int}}(), String[], Dict{String,Int}())
    best_form[] === nothing && throw(KodaImplementationError(
        :hc_canonical, "search produced no form", Dict{String,Any}()))
    sort!(best_witnesses, by=_hc_witness_rank)
    best_form[], best_witnesses
end

"""
同一性の正規形を返す。許すのは交点ラベルの付け替え、各語の巡回、role の中での成分名の付け替えだけである。

曲線の向きの取り替えと鏡像は正規化しない。K の成分名は基準曲線の名前として固定する。
成分の順序は、K を名前順にしたあと、色細分で分かれた α、β の順である。同じ色の組は全部の並べ方を見る。
語の巡回で、出した列は同じでもラベルの写像が違うものは、後続の曲線まで含めて辞書式最小になる方を残す。
保存された語の開始位置では選ばない。探索が予算を超えたときは HCNotComputed を投げ、
同一かどうか分からないのに different とは返さない。
"""
function hc_canonical(hc::HeegaardLayer)
    form, _ = _hc_canonicalize(hc)
    form
end

"""
入力から正規形への対応を返す。

label_map は元の交点ラベルから正規形のラベル、name_map は元の成分名から正規形の成分名である。
rotation_starts は元の語の 1 始まりの位置で、そこを先頭にすると正規形の語になる。
"""
function hc_canonical_witness(hc::HeegaardLayer)
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
許した変換を施す。rotation_starts は元の語の 1 始まりの位置で、
そこへ回し、label_map と name_map を施すと、相手の保存された語になる。
"""
function hc_apply_correspondence(hc::HeegaardLayer, corr::HCIdentityCorrespondence)
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
    lc = LickorishConfiguration(
        hc.lc.genus, words, places,
        _hc_transport_disk(hc.lc.disks, transport),
        _hc_transport_point(hc.lc.points, transport),
        _hc_transport_weight(hc.lc.weights, transport),
        nothing,
    )
    HeegaardLayer(
        lc, roles,
        _hc_transport_disk(hc.disks, transport),
        _hc_transport_point(hc.points, transport),
        _hc_transport_weight(hc.weights, transport),
        nothing,
    )
end

function _hc_compose_witness(a::HeegaardLayer, left::HCCanonicalWitness, right::HCCanonicalWitness)
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

function _hc_stored_words(hc::HeegaardLayer)
    Dict{String,Vector{Int}}(String(name) => copy(word) for (name, word) in hc.lc.words)
end

"""
同一な 2 つの HC の間の、許した変換を返す。同一でないときは HCNotComputed(:not_identical) を投げる。

label_map は a の交点ラベルから b の交点ラベル、name_map は a の成分名から b の成分名である。
rotation_starts は a の語の 1 始まりの位置で、そこへ回して 2 つの写像を施すと b の保存された語になる。
正規形が一致しても、取り出した対応が a を b の保存された語へ移さないときは
HCNotComputed(:correspondence_inconsistent) を投げる。
"""
function hc_identity_correspondence(a::HeegaardLayer, b::HeegaardLayer)
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
