# cf: lib/topology.jl/docs/specs/heegaard-diagram.md — `standard` binds each α to an existing meridian m_i
# cf: lib/topology.jl/docs/specs/lickorish-code.md §3 — a component word is a cyclic sequence of signed passages
# cf: lib/topology.jl/docs/specs/heegaard-moves.md §6.1 — K-edge bigon removal, K-vertex passage, input-crossing slide

"""One α that this check does not bind to a meridian.

`homology_class` is a proof that the curve is not isotopic to that meridian.
`undecided` means the isotopy search stopped without that proof.
"""
struct HCStandardFailure
    curve::String
    reason::Symbol
    detail::Dict{String,Any}
end

"""`standard` judgment.

A true result binds α to `m1,…,mg`. `witness` lists the isotopy row ids that
took a different word to that meridian. An empty witness means no isotopy row
was applied. Reversing an α is not a row id.
"""
struct HCStandardCheck
    standard::Bool
    diagram::Union{Nothing,HeegaardCode}
    correspondence::Vector{Pair{String,String}}
    failures::Vector{HCStandardFailure}
    witness::Vector{String}
end

# Distinct diagrams expanded while searching the three isotopy moves.
# Fixed before any example is run; a later failure does not raise it.
const _HC_STANDARD_STATE_LIMIT = 10_000

function _hc_signed_partners(curve::AbstractString, word::Vector{Int}, owners::Dict{Int,Tuple{String,String}})
    partners = Tuple{String,Int}[]
    for x in word
        plus_c, minus_c = owners[abs(x)]
        if x > 0
            plus_c == curve || throw(KodaInputError(:hc_standard_check, "plus passage is not on its curve",
                Dict{String,Any}("curve" => String(curve), "label" => x)))
            push!(partners, (minus_c, 1))
        else
            minus_c == curve || throw(KodaInputError(:hc_standard_check, "minus passage is not on its curve",
                Dict{String,Any}("curve" => String(curve), "label" => x)))
            push!(partners, (plus_c, -1))
        end
    end
    partners
end

function _hc_cycles_equal(left::Vector{Tuple{String,Int}}, right::Vector{Tuple{String,Int}})
    length(left) == length(right) || return false
    n = length(left)
    n == 0 && return true
    for shift in 0:(n - 1)
        if all(i -> left[i] == right[mod1(i + shift, n)], 1:n)
            return true
        end
    end
    false
end

function _hc_matchings(graph::Vector{Vector{Int}}, limit::Int)
    n = length(graph)
    used = falses(n)
    found = Vector{Vector{Int}}()
    nodes = Ref(0)
    assignment = zeros(Int, n)
    function search(i::Int)
        length(found) >= limit && return
        if i > n
            push!(found, copy(assignment))
            return
        end
        nodes[] += 1
        nodes[] > 100_000 && return
        for m in graph[i]
            used[m] && continue
            used[m] = true
            assignment[i] = m
            search(i + 1)
            used[m] = false
            length(found) >= limit && return
        end
    end
    search(1)
    found, nodes[] > 100_000 && length(found) < limit
end

function _hc_decoration_mentions(text::AbstractString, curve::AbstractString)
    occursin(Regex("(^|[^A-Za-z0-9_])" * escape_string(curve) * "([^A-Za-z0-9_]|\$)"), text)
end

function _hc_decoration_block(hc::HeegaardCode, dropped::Vector{String})
    isempty(dropped) && return nothing
    texts = String[]
    for (name, ref) in hc.lc.places
        push!(texts, name * " " * ref)
    end
    for (name, refs) in hc.lc.disks
        push!(texts, name * " " * join(refs, " "))
    end
    for (name, ref) in hc.lc.points
        push!(texts, name * " " * ref)
    end
    for (name, weight) in hc.lc.weights
        push!(texts, name * " " * string(weight))
    end
    for (name, refs) in hc.disks
        push!(texts, name * " " * join(refs, " "))
    end
    for (name, ref) in hc.points
        push!(texts, name * " " * ref)
    end
    for (name, weight) in hc.weights
        push!(texts, name * " " * string(weight))
    end
    for curve in dropped, text in texts
        if curve == text || _hc_decoration_mentions(text, curve)
            return HCStandardFailure(curve, :decoration_blocks_removal,
                Dict{String,Any}("support" => text))
        end
    end
    nothing
end

function _hc_without_curves(lc::LickorishCode, dropped::Vector{String})
    drop = Set(dropped)
    labels = Set{Int}()
    for name in dropped
        for x in lc.words[name]
            push!(labels, abs(x))
        end
    end
    words = Dict{String,Vector{Int}}()
    for (name, word) in lc.words
        sname = String(name)
        sname in drop && continue
        words[sname] = Int[x for x in word if abs(x) ∉ labels]
    end
    places = Dict{String,String}(String(k) => String(v) for (k, v) in lc.places if !(String(k) in drop))
    disks = Dict{String,Vector{String}}(
        String(k) => String.(v) for (k, v) in lc.disks if !(String(k) in drop))
    points = Dict{String,String}(String(k) => String(v) for (k, v) in lc.points if !(String(k) in drop))
    weights = Dict{String,Int}(String(k) => Int(v) for (k, v) in lc.weights if !(String(k) in drop))
    LickorishCode(lc.genus, words, places, disks, points, weights, nothing)
end

function _hc_check(standard::Bool, diagram, correspondence, failures; witness::Vector{String}=String[])
    HCStandardCheck(standard, diagram, correspondence, failures, witness)
end

function _hc_failure(curve::AbstractString, reason::Symbol, detail::Dict{String,Any}=Dict{String,Any}())
    _hc_check(false, nothing, Pair{String,String}[],
        [HCStandardFailure(String(curve), reason, detail)])
end

# + on a and − on b means the frame (b, a) is positive, so i(a, b) = −1.
function _hc_algebraic_intersection(curve::AbstractString, other::AbstractString, owners::Dict{Int,Tuple{String,String}})
    curve == other && return 0
    count = 0
    for (_, (plus_c, minus_c)) in owners
        if plus_c == curve && minus_c == other
            count -= 1
        elseif minus_c == curve && plus_c == other
            count += 1
        end
    end
    count
end

# Coefficients of [m_1], [l_1], …, [m_g], [l_g]. x_j = −i(α, l_j) and y_j = i(α, m_j).
function _hc_homology_class(curve::AbstractString, genus::Integer, owners::Dict{Int,Tuple{String,String}})
    class = Int[]
    for j in 1:genus
        push!(class, -_hc_algebraic_intersection(curve, "l$j", owners))
        push!(class, _hc_algebraic_intersection(curve, "m$j", owners))
    end
    class
end

_hc_opposite_class(class::Vector{Int}) = .-class

function _hc_same_word(hc::HeegaardCode)
    validate_hc_input(hc)
    genus = hc.lc.genus
    owners = _hc_owners(hc.lc.words)
    alphas = hc_curve_names(hc, :alpha)
    meridians = String["m$i" for i in 1:genus]
    if length(alphas) != genus
        return _hc_failure("", :count_mismatch,
            Dict{String,Any}("alpha_count" => length(alphas), "genus" => genus))
    end
    cycles = Dict(name => _hc_signed_partners(name, hc.lc.words[name], owners) for name in alphas)
    meridian_cycles = Dict{String,Vector{Tuple{String,Int}}}()
    for name in meridians
        haskey(hc.lc.words, name) || return _hc_failure(name, :missing_meridian, Dict{String,Any}())
        meridian_cycles[name] = _hc_signed_partners(name, hc.lc.words[name], owners)
    end
    graph = Vector{Vector{Int}}()
    for alpha in alphas
        hits = Int[j for (j, meridian) in enumerate(meridians) if _hc_cycles_equal(cycles[alpha], meridian_cycles[meridian])]
        push!(graph, hits)
    end
    if any(isempty, graph)
        return nothing
    end
    found, exhausted = _hc_matchings(graph, 2)
    if exhausted
        return _hc_failure("", :correspondence_unresolved, Dict{String,Any}("search_nodes" => 100_000))
    end
    if isempty(found)
        return _hc_failure("", :no_bijection, Dict{String,Any}("alphas" => alphas, "meridians" => meridians))
    end
    if length(found) > 1
        return _hc_failure("", :ambiguous_correspondence, Dict{String,Any}("alphas" => alphas))
    end
    match = found[1]
    correspondence = Pair{String,String}[alphas[i] => meridians[match[i]] for i in eachindex(alphas)]
    for (alpha, meridian) in correspondence
        role = get(hc.roles, meridian, :none)
        if role == :beta || (role == :alpha && alpha != meridian)
            return _hc_failure(alpha, :meridian_occupied,
                Dict{String,Any}("meridian" => meridian, "role" => String(role)))
        end
    end
    dropped = String[alpha for (alpha, meridian) in correspondence if alpha != meridian]
    meridian_labels = Set{Int}()
    for meridian in meridians, x in hc.lc.words[meridian]
        push!(meridian_labels, abs(x))
    end
    for name in dropped, x in hc.lc.words[name]
        if abs(x) in meridian_labels
            return _hc_failure(name, :intersects_meridian,
                Dict{String,Any}("label" => abs(x)))
        end
    end
    blocked = _hc_decoration_block(hc, dropped)
    blocked !== nothing && return _hc_check(false, nothing, Pair{String,String}[], [blocked])
    if isempty(dropped)
        return _hc_check(true, hc, correspondence, HCStandardFailure[])
    end
    roles = Dict{String,Symbol}()
    for (name, role) in hc.roles
        sname = String(name)
        sname in dropped && continue
        roles[sname] = role
    end
    for meridian in meridians
        roles[meridian] = :alpha
    end
    reduced = resolve_hc(_hc_without_curves(hc.lc, dropped), roles;
        disks=Dict{String,Vector{String}}(String(k) => String.(v) for (k, v) in hc.disks),
        points=Dict{String,String}(String(k) => String(v) for (k, v) in hc.points),
        weights=Dict{String,Int}(String(k) => Int(v) for (k, v) in hc.weights))
    is_standard(reduced) || return _hc_failure("", :reduction_not_standard, Dict{String,Any}())
    _hc_check(true, reduced, correspondence, HCStandardFailure[])
end

function _hc_unmatched_alphas(hc::HeegaardCode)
    genus = hc.lc.genus
    owners = _hc_owners(hc.lc.words)
    alphas = hc_curve_names(hc, :alpha)
    meridians = String["m$i" for i in 1:genus]
    cycles = Dict(name => _hc_signed_partners(name, hc.lc.words[name], owners) for name in alphas)
    unmatched = String[]
    for alpha in alphas
        hit = any(meridian -> _hc_cycles_equal(cycles[alpha],
            _hc_signed_partners(meridian, hc.lc.words[meridian], owners)), meridians)
        hit || push!(unmatched, alpha)
    end
    unmatched, owners, meridians
end

function _hc_negative_meridian_alphas(hc::HeegaardCode, unmatched::Vector{String}, owners, meridians)
    classes = Dict(name => _hc_homology_class(name, hc.lc.genus, owners) for name in vcat(unmatched, meridians))
    chosen = String[]
    for alpha in unmatched
        aclass = classes[alpha]
        positive = any(meridian -> aclass == classes[meridian], meridians)
        negative = any(meridian -> aclass == _hc_opposite_class(classes[meridian]), meridians)
        negative && !positive && push!(chosen, alpha)
    end
    chosen
end

function _hc_homology_failures(hc::HeegaardCode, unmatched::Vector{String}, owners, meridians)
    classes = Dict(name => _hc_homology_class(name, hc.lc.genus, owners) for name in vcat(unmatched, meridians))
    failures = HCStandardFailure[]
    for alpha in unmatched
        aclass = classes[alpha]
        hits = [meridian for meridian in meridians if aclass == classes[meridian] || aclass == _hc_opposite_class(classes[meridian])]
        isempty(hits) || continue
        detail = Dict{String,Any}("alpha_class" => aclass)
        if length(meridians) == 1
            detail["meridian"] = meridians[1]
            detail["meridian_class"] = classes[meridians[1]]
        else
            detail["meridian_classes"] = Dict(meridian => classes[meridian] for meridian in meridians)
        end
        push!(failures, HCStandardFailure(alpha, :homology_class, detail))
    end
    failures
end

function _hc_isotopy_row_ids()
    ids = String[]
    for id in heegaard_move_ids()
        row = String(id)
        if startswith(row, "lc.isotopy.digon.remove.") && !occursin(".inverse.", row)
            push!(ids, row)
        elseif startswith(row, "lc.isotopy.vertex.") || startswith(row, "lc.isotopy.crossing.")
            push!(ids, row)
        end
    end
    ids
end

function _hc_move_touches(witness::LCMoveWitness, curves)
    involved = Set(witness.names[comp] for (comp, _) in witness.locations)
    any(in(involved), curves)
end

function _hc_beta_sequence(hc::HeegaardCode, curve::AbstractString, owners)
    betas = Set(hc_curve_names(hc, :beta))
    Tuple{String,Int}[item for item in _hc_signed_partners(curve, hc.lc.words[curve], owners) if first(item) in betas]
end

function _hc_beta_preserved(before::HeegaardCode, after::HeegaardCode, alphas)
    before_owners = _hc_owners(before.lc.words)
    after_owners = _hc_owners(after.lc.words)
    for alpha in alphas
        haskey(after.lc.words, alpha) || return false
        _hc_cycles_equal(_hc_beta_sequence(before, alpha, before_owners),
            _hc_beta_sequence(after, alpha, after_owners)) || return false
    end
    true
end

function _hc_state_key(hc::HeegaardCode)
    join((string(name, "=", _lc_canon_cyclic(hc.lc.words[name]))
        for name in sort!(collect(keys(hc.lc.words)))), ";")
end

function _hc_undecided(alphas, seen, stopped)
    failures = [HCStandardFailure(alpha, :undecided, Dict{String,Any}(
        "scope" => "isotopy_search",
        "states" => length(seen),
        "limit" => _HC_STANDARD_STATE_LIMIT,
        "stopped" => stopped,
    )) for alpha in alphas]
    _hc_check(false, nothing, Pair{String,String}[], failures)
end

function _hc_search_isotopy(hc::HeegaardCode, alphas)
    rows = _hc_isotopy_row_ids()
    seen = Set([_hc_state_key(hc)])
    queue = Tuple{HeegaardCode,Vector{String}}[(hc, String[])]
    head = 1
    while head <= length(queue)
        current, path = queue[head]
        head += 1
        for row in rows
            hits = try
                movable(current, row)
            catch
                LCMoveWitness[]
            end
            for witness in hits
                _hc_move_touches(witness, alphas) || continue
                moved = try
                    move(current, witness)
                catch
                    nothing
                end
                moved === nothing && continue
                nxt = moved.output
                nxt isa HeegaardCode || continue
                _hc_beta_preserved(hc, nxt, alphas) || continue
                key = _hc_state_key(nxt)
                key in seen && continue
                if length(seen) >= _HC_STANDARD_STATE_LIMIT
                    return _hc_undecided(alphas, seen, "limit")
                end
                push!(seen, key)
                next_path = vcat(path, [row])
                checked = _hc_same_word(nxt)
                # A bound meridian word ends the search. A child that matches
                # the word but cannot be bound is not a proof about the start.
                if checked !== nothing
                    checked.standard || continue
                    return _hc_check(true, checked.diagram, checked.correspondence,
                        HCStandardFailure[]; witness=next_path)
                end
                push!(queue, (nxt, next_path))
            end
        end
    end
    _hc_undecided(alphas, seen, "exhausted")
end

"""Judge whether each α can be written as a meridian of K.

The same cyclic signed word is `standard`. A homology class other than
±[m_i] is `homology_class`. A class equal to −[m_i] and to no +[m_j] is
reversed with `hc_reverse_orientation` before the same judgment: the signed
word, or the three-move search. A class equal to +[m_i] is searched without
reversing. Passing the state limit, or finishing that search without the
meridian word, is `undecided`.
"""
function hc_standard_check(hc::HeegaardCode)
    checked = _hc_same_word(hc)
    checked !== nothing && return checked
    unmatched, owners, meridians = _hc_unmatched_alphas(hc)
    proved = _hc_homology_failures(hc, unmatched, owners, meridians)
    isempty(proved) || return _hc_check(false, nothing, Pair{String,String}[], proved)
    to_reverse = _hc_negative_meridian_alphas(hc, unmatched, owners, meridians)
    oriented = isempty(to_reverse) ? hc : hc_reverse_orientation(hc, to_reverse)
    if !isempty(to_reverse)
        checked = _hc_same_word(oriented)
        checked !== nothing && return checked
        unmatched, _, _ = _hc_unmatched_alphas(oriented)
    end
    _hc_search_isotopy(oriented, unmatched)
end

export HCStandardFailure, HCStandardCheck, hc_standard_check
