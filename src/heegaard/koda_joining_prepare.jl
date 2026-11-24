# cf: docs/specs/heegaard-koda-conversion.md §5.2 — a joining disk meets every α and every β in one proper arc.
# The rows are lc.isotopy.digon.create.ab and hc.slide.pants. Both preserve the closed manifold:
# the digon is a curve isotopy, and a same-family slide is a handle slide.
# An annulus left in the complement is cut by a finger across both of its boundaries.
# The finger is written back into the HeegaardCode. It is not a catalogued move row.

"""Prepare `hc` so `enumerate_joining_disk_candidates` is nonempty and no
complementary annulus remains.

Return `(hc=prepared, witnesses=applied)`. `witnesses` is the `movable` /
`move` sequence that produces a joining disk. An input that already has a
joining disk contributes an empty sequence. The sequence is deterministic:

1. While a β and an α share no label and bound one face, push that β across
   that α with the first `lc.isotopy.digon.create.ab` witness.
2. Otherwise try same-family `hc.slide.pants` witnesses, β pairs before α
   pairs, and repeat step 1 on that result. Keep the first result that has a
   joining disk.
3. If the curves already intersect and neither step applies, push the first
   β–α pair once. This is the genus-1 bigon.

After that search, while a complementary region is an annulus, take the first
opposite-bank pair whose β edge and α edge are each one boundary step of a
single annulus face. If no such pair exists, take the first pair whose edges
each meet one annulus face along their whole length. Push that β across the
curve that separates those two faces, then across α on the other bank. The
push starts at the boundary step when that step is a single corner. Each push
is one digon written into the LC words. One pass cuts one annulus, so two
annuli take two passes.
`hc.stabilization.once` is not used. A diagram where these rows still produce
no joining disk raises `KodaInputError`.

Pass `trace` as a `Vector` to record `(stage, seconds)` for each search step.
"""
function hc_prepare_joining_disk(hc::HeegaardCode; trace::Union{Nothing,Vector}=nothing)
    validate_hc_input(hc)
    found = _prepare_search_joining(hc, trace)
    t0 = time_ns()
    cut = _prepare_cut_annuli(found.hc)
    _prepare_mark(trace, :annulus_finger, t0)
    (hc=cut, witnesses=found.witnesses)
end

function _prepare_mark(trace, stage, t0)
    trace === nothing && return
    seconds = (time_ns() - t0) / 1.0e9
    push!(trace, (stage=stage, seconds=seconds))
end

function _prepare_search_joining(hc::HeegaardCode, trace)
    _prepare_has_joining(hc; trace=trace) && return (hc=hc, witnesses=LCMoveWitness[])
    trial = LCMoveWitness[]
    tried = _prepare_zero_digons(hc, trial, trace)
    _prepare_has_joining(tried; trace=trace) && return (hc=tried, witnesses=trial)
    for curves in _prepare_slide_pairs(hc)
        t0 = time_ns()
        slides = try
            movable(hc, "hc.slide.pants"; curves=curves, band=Any[])
        catch
            LCMoveWitness[]
        end
        _prepare_mark(trace, :slide_movable, t0)
        isempty(slides) && continue
        applied = LCMoveWitness[first(slides)]
        slid = _prepare_zero_digons(move(hc, first(slides)).output, applied, trace)
        _prepare_has_joining(slid; trace=trace) && return (hc=slid, witnesses=applied)
    end
    applied = LCMoveWitness[]
    once = _prepare_one_digon(hc, applied, trace)
    _prepare_has_joining(once; trace=trace) && return (hc=once, witnesses=applied)
    throw(KodaInputError(
        :hc_prepare_joining_disk,
        "digon creation and same-family slides did not produce a joining disk",
        Dict(
            "missing_move" => "finger of one beta along a face path to an alpha that shares neither a label nor a face",
            "rows" => ["lc.isotopy.digon.create.ab", "hc.slide.pants"],
        ),
    ))
end

function _prepare_has_joining(hc::HeegaardCode; trace=nothing)
    t0 = time_ns()
    ok = !isempty(enumerate_joining_disk_candidates(forget_reference(realize_hc(hc))))
    _prepare_mark(trace, :has_joining, t0)
    ok
end

function _prepare_shared(lc::LickorishCode, left::AbstractString, right::AbstractString)
    labels = Set(abs.(get(lc.words, left, Int[])))
    count(entry -> abs(entry) in labels, get(lc.words, right, Int[]))
end

function _prepare_digon_index(hc::HeegaardCode; trace=nothing)
    t0 = time_ns()
    alphas = Set(string.(hc_curve_names(hc, :alpha)))
    betas = Set(string.(hc_curve_names(hc, :beta)))
    found = Dict{Tuple{String,String}, LCMoveWitness}()
    order = Tuple{String,String}[]
    for witness in movable(hc, "lc.isotopy.digon.create.ab")
        key = (string(witness.extra["push"]), string(witness.extra["cross"]))
        (key[1] in betas && key[2] in alphas) || continue
        haskey(found, key) && continue
        found[key] = witness
        push!(order, key)
    end
    _prepare_mark(trace, :digon_index, t0)
    (order, found)
end

function _prepare_zero_digons(hc::HeegaardCode, witnesses::Vector{LCMoveWitness}, trace=nothing)
    steps = max(1, hc.lc.genus^2)
    for _ in 1:steps
        _prepare_has_joining(hc; trace=trace) && return hc
        order, found = _prepare_digon_index(hc; trace=trace)
        zeros = [key for key in order if _prepare_shared(hc.lc, key[1], key[2]) == 0]
        isempty(zeros) && return hc
        key = zeros[1]
        before = _prepare_shared(hc.lc, key[1], key[2])
        witness = found[key]
        hc = move(hc, witness).output
        push!(witnesses, witness)
        _prepare_shared(hc.lc, key[1], key[2]) > before || return hc
    end
    hc
end

function _prepare_one_digon(hc::HeegaardCode, witnesses::Vector{LCMoveWitness}, trace=nothing)
    _prepare_has_joining(hc; trace=trace) && return hc
    order, found = _prepare_digon_index(hc; trace=trace)
    isempty(order) && return hc
    witness = found[order[1]]
    push!(witnesses, witness)
    move(hc, witness).output
end

function _prepare_slide_pairs(hc::HeegaardCode)
    jobs = Vector{Vector{String}}()
    for family in (hc_curve_names(hc, :beta), hc_curve_names(hc, :alpha))
        length(family) < 2 && continue
        for i in eachindex(family), j in eachindex(family)
            i == j && continue
            push!(jobs, String[family[i], family[j]])
        end
    end
    jobs
end

function _prepare_fail(message, detail=Dict{String,Any}())
    throw(KodaImplementationError(:hc_prepare_joining_disk, message, detail))
end

function _complement_annulus_count(hc::HeegaardCode)
    realized = realize_hc(hc)
    forgot = forget_reference(realized)
    forgot.diagram === nothing && return 0
    cands = enumerate_joining_disk_candidates(forgot)
    isempty(cands) && return 0
    rid = minimum(c["region_id"] for c in cands)
    region = only(r for r in forgot.regions if r.id == rid)
    source = _koda_original_complement(forgot, forgot.diagram, region)
    count(r -> r.genus == 0 && r.chi == 0 && length(r.cycles) == 2, source.regions)
end

function _prepare_cut_annuli(hc::HeegaardCode)
    for step in 1:8
        before = _complement_annulus_count(hc)
        before == 0 && return hc
        nxt = try
            _cut_one_annulus(hc)
        catch err
            detail = hasproperty(err, :detail) ? copy(err.detail) : Dict{String,Any}()
            detail["finger_step"] = step
            detail["annuli_before"] = before
            _prepare_fail(sprint(showerror, err), detail)
        end
        nxt === hc && _prepare_fail("annulus finger made no change",
            Dict("annuli" => before, "step" => step))
        after = _complement_annulus_count(nxt)
        after < before || _prepare_fail("annulus finger kept the complementary annulus",
            Dict("before" => before, "after" => after, "step" => step))
        hc = nxt
    end
    _prepare_fail("annulus finger did not finish", Dict("steps" => 8))
end

function _curve_named(curves, id::Int)
    hit = [curve for curve in curves if curve.id == id]
    length(hit) == 1 || return nothing
    hit[1].label === nothing ? nothing : String(hit[1].label)
end

"""Cut one complementary annulus, or return `hc` unchanged.

The pass starts at the lowest-id annulus. The new β arc then selects the
annulus it bounds on exactly one side, rather than the lowest remaining id.
"""
function _cut_one_annulus(hc::HeegaardCode)
    realized = realize_hc(hc)
    forgot = forget_reference(realized)
    forgot.diagram === nothing && return hc
    h = forgot.diagram
    cands = enumerate_joining_disk_candidates(forgot)
    isempty(cands) && return hc
    rid = minimum(c["region_id"] for c in cands)
    region = only(r for r in forgot.regions if r.id == rid)
    source = _koda_original_complement(forgot, h, region)
    annuli = [r for r in source.regions if r.genus == 0 && r.chi == 0 && length(r.cycles) == 2]
    isempty(annuli) && return hc
    ann = first(sort!(annuli, by=r -> r.id))
    graph = source.graph
    banks = [copy(cycle) for cycle in source.cycles if source.owner[first(cycle)] == ann.id]
    length(banks) == 2 || _prepare_fail("annulus does not have two boundary components",
        Dict("region" => ann.id, "banks" => length(banks)))
    jd = _find_ribbon_joining_disk(h, region)
    jd === nothing && _prepare_fail("joining disk is missing on the ribbon", Dict("region" => rid))
    old = _undecorated(h)
    sumchi = sum(r.chi for r in source.regions)
    length(old.intersections) - length(graph.edges) + sumchi == 2 - 2 * hc.lc.genus ||
        _prepare_fail("surface Euler characteristic does not match the genus",
            Dict("genus" => hc.lc.genus, "chi_sum" => sumchi))
    removed = Set(graph.darts[dart_id].edge for dart_id in jd.boundary_darts)
    pairs = Tuple{Int,Int}[]
    for side in 1:2, push_dart in banks[side], cross_dart in banks[3 - side]
        push_edge = graph.edges[graph.darts[push_dart].edge]
        cross_edge = graph.edges[graph.darts[cross_dart].edge]
        push_edge.family == :beta && cross_edge.family == :alpha || continue
        push_edge.id in removed && continue
        cross_edge.id in removed && continue
        push!(pairs, (push_dart, cross_dart))
    end
    isempty(pairs) && _prepare_fail("annulus has no beta-alpha finger across both boundaries",
        Dict("region" => ann.id))
    real = forgot.realization
    ann_faces = Set(ann.faces)
    pids = _source_label_ids(h)
    sort!(pairs)
    single = nothing
    broad = nothing
    for (pd, cd) in pairs
        push_edge = graph.edges[graph.darts[pd].edge]
        cross_edge = graph.edges[graph.darts[cd].edge]
        beta_name = _curve_named(old.beta, push_edge.curve_id)
        alpha_name = _curve_named(old.alpha, cross_edge.curve_id)
        (beta_name === nothing || alpha_name === nothing) && continue
        beta_step = _boundary_step(hc, real, pids, ann_faces, beta_name, push_edge)
        alpha_step = _boundary_step(hc, real, pids, ann_faces, alpha_name, cross_edge)
        if beta_step !== nothing && alpha_step !== nothing
            single = (beta_name, alpha_name, beta_step.face, alpha_step.face, beta_step.corner)
            break
        end
        broad !== nothing && continue
        beta_face = _unique_annulus_face(hc, real, pids, ann_faces, beta_name, push_edge)
        alpha_face = _unique_annulus_face(hc, real, pids, ann_faces, alpha_name, cross_edge)
        (beta_face === nothing || alpha_face === nothing) && continue
        broad = (beta_name, alpha_name, beta_face, alpha_face, nothing)
    end
    selected = single === nothing ? broad : single
    selected === nothing && _prepare_fail("annulus has no single-face beta-alpha pair",
        Dict("region" => ann.id, "pairs" => length(pairs)))
    beta_name, alpha_name, beta_face, alpha_face, beta_corner = selected
    before = _label_set(hc)
    if beta_face == alpha_face
        return _apply_digon(hc, beta_name, alpha_name, beta_face; corner=beta_corner)
    end
    wall = _separating_curve(hc, real, ann_faces, beta_face, alpha_face, Set((beta_name, alpha_name)))
    crossed = _apply_digon(hc, beta_name, wall, beta_face; corner=beta_corner)
    _apply_new_arc_digon(crossed, beta_name, alpha_name, before)
end

function _label_set(hc::HeegaardCode)
    Set(abs(letter) for word in values(hc.lc.words) for letter in word)
end

function _step_faces(real, x::Int, y::Int)
    right = real.corner_face[(abs(x), _lc_sign_char(x), "out")]
    left = real.corner_face[(abs(y), _lc_sign_char(y), "in")]
    (right, left)
end

function _ribbon_steps(word::Vector{Int}, pids, from_pid::Int, to_pid::Int)
    visits = [i for i in eachindex(word) if haskey(pids, abs(word[i]))]
    for k in eachindex(visits)
        start = visits[k]
        stop = visits[mod1(k + 1, length(visits))]
        pids[abs(word[start])] == from_pid || continue
        pids[abs(word[stop])] == to_pid || continue
        steps = Tuple{Int,Int}[]
        i = start
        while true
            j = mod1(i + 1, length(word))
            push!(steps, (i, j))
            j == stop && break
            i = j
        end
        return steps
    end
    nothing
end

function _boundary_step(hc::HeegaardCode, real, pids, ann_faces, name::String, edge)
    steps = _ribbon_steps(hc.lc.words[name], pids, edge.from, edge.to)
    (steps === nothing || length(steps) != 1) && return nothing
    i, j = only(steps)
    word = hc.lc.words[name]
    right, left = _step_faces(real, word[i], word[j])
    if right in ann_faces && !(left in ann_faces)
        return (face=right, corner=(abs(word[i]), _lc_sign_char(word[i]), "out"))
    elseif left in ann_faces && !(right in ann_faces)
        return (face=left, corner=(abs(word[j]), _lc_sign_char(word[j]), "in"))
    end
    nothing
end

function _unique_annulus_face(hc::HeegaardCode, real, pids, ann_faces, name::String, edge)
    steps = _ribbon_steps(hc.lc.words[name], pids, edge.from, edge.to)
    steps === nothing && return nothing
    face = nothing
    for (i, j) in steps
        right, left = _step_faces(real, hc.lc.words[name][i], hc.lc.words[name][j])
        on_right = right in ann_faces
        on_left = left in ann_faces
        on_right == on_left && return nothing
        side = on_right ? right : left
        face === nothing && (face = side)
        face == side || return nothing
    end
    face
end

function _separating_curve(hc::HeegaardCode, real, ann_faces, beta_face, alpha_face, skip)
    wanted = Set((beta_face, alpha_face))
    hits = Tuple{String,Int}[]
    for name in sort(collect(keys(hc.lc.words)))
        name in skip && continue
        word = hc.lc.words[name]
        for i in eachindex(word)
            right, left = _step_faces(real, word[i], word[mod1(i + 1, length(word))])
            right in ann_faces && left in ann_faces && Set((right, left)) == wanted || continue
            push!(hits, (name, i))
        end
    end
    isempty(hits) && _prepare_fail("annulus banks do not share a separating curve",
        Dict("beta_face" => beta_face, "alpha_face" => alpha_face))
    first(hits)[1]
end

function _corner_matches(spec, corner)
    length(spec) >= 3 && spec[1] == corner[1] && string(spec[2]) == corner[2] && string(spec[3]) == corner[3]
end

function _apply_digon(hc::HeegaardCode, push_name::String, cross_name::String, face; corner=nothing)
    # A digon witness stores tau of the face corner, the far end of that boundary edge.
    active = _lc_nonempty_words(hc.lc.words)
    occ = _lc_rotation_occ(active)
    expected = corner === nothing ? nothing :
        _lc_tau(active, occ, (Int(corner[1]), String(corner[2]), String(corner[3])))
    hits = LCMoveWitness[]
    for witness in movable(hc, "lc.isotopy.digon.create.ab")
        string(witness.extra["push"]) == push_name || continue
        string(witness.extra["cross"]) == cross_name || continue
        string(witness.extra["face"]) == string(face) || continue
        if expected !== nothing
            start = get(get(witness.extra, "support", Dict()), "start", Any[])
            _corner_matches(start, expected) || continue
        end
        push!(hits, witness)
    end
    isempty(hits) && _prepare_fail("annulus finger has no digon on the selected face",
        Dict("push" => push_name, "cross" => cross_name, "face" => string(face),
             "corner" => expected === nothing ? nothing : collect(expected)))
    move(hc, first(hits)).output
end

function _apply_new_arc_digon(hc::HeegaardCode, beta_name::String, alpha_name::String, previous)
    fresh = setdiff(_label_set(hc), previous)
    length(fresh) == 2 || _prepare_fail("wall crossing did not add two intersections",
        Dict("fresh" => sort(collect(fresh))))
    realized = realize_hc(hc)
    forgot = forget_reference(realized)
    forgot.diagram === nothing && _prepare_fail("wall crossing has no diagram", Dict())
    real = forgot.realization
    cands = enumerate_joining_disk_candidates(forgot)
    isempty(cands) && _prepare_fail("wall crossing removed the joining disk", Dict())
    rid = minimum(c["region_id"] for c in cands)
    region = only(r for r in forgot.regions if r.id == rid)
    source = _koda_original_complement(forgot, forgot.diagram, region)
    annuli = [r for r in source.regions if r.genus == 0 && r.chi == 0 && length(r.cycles) == 2]
    isempty(annuli) && return hc
    word = hc.lc.words[beta_name]
    slot = nothing
    for i in eachindex(word)
        j = mod1(i + 1, length(word))
        abs(word[i]) in fresh && abs(word[j]) in fresh || continue
        slot = (i, j)
        break
    end
    slot === nothing && _prepare_fail("new beta arc is not a single edge", Dict("curve" => beta_name))
    i, j = slot
    right, left = _step_faces(real, word[i], word[j])
    # One annulus, the one this edge bounds on one side. Do not pick by id.
    hits = Tuple{Any,Tuple{Int,String,String}}[]
    for ann in annuli
        faces = Set(ann.faces)
        on_right = right in faces
        on_left = left in faces
        if on_right && !on_left
            push!(hits, (right, (abs(word[i]), _lc_sign_char(word[i]), "out")))
        elseif on_left && !on_right
            push!(hits, (left, (abs(word[j]), _lc_sign_char(word[j]), "in")))
        end
    end
    length(hits) == 1 || _prepare_fail("new beta arc does not bound the annulus on one side",
        Dict("right" => right, "left" => left, "one_sided" => length(hits)))
    face, corner = only(hits)
    _apply_digon(hc, beta_name, alpha_name, face; corner=corner)
end

export hc_prepare_joining_disk
