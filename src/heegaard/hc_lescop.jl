# Source: docs/specs/lickorish-code.md §3 (passage sign, lines 83–86);
# docs/specs/heegaard-koda-conversion.md §3.1 (line 122) and §6.1 (line 212);
# docs/specs/heegaard-terminology.md §4 (lines 61–62).
# Lescop σ(c) = sign det(t_α, t_β) on ∂H_α. The LC frame (b, a) is positive
# when b carries −k and a carries +k, and S is oriented as ∂H_α, so
# σ(c) = −(passage sign on α) = (passage sign on β).
"""Forget roleless K and evaluate Lescop's Θ on the resulting cellular diagram."""

struct HCLescopView
    diagram::HeegaardDiagram
    face_region::Dict{String,Int}
    label_of_id::Dict{Int,Int}
end

function _lescop_check_reference!(lc::LickorishCode)
    g = lc.genus
    g < 1 && return
    limit = 3g - 2
    ref = _lc_k_constants(g)
    for name in _lc_k_names(g)
        haskey(lc.words, name) || throw(ArgumentError("missing reference curve $name"))
        kept = Int[x for x in lc.words[name] if abs(x) <= limit]
        _lc_canon_cyclic(kept) == _lc_canon_cyclic(ref[name]) ||
            throw(ArgumentError("reference curve $name does not restrict to the Lickorish constants"))
    end
end

function _lescop_passage(word::Vector{Int}, label::Int)
    hits = [x for x in word if abs(x) == label]
    length(hits) == 1 || throw(ArgumentError("label $label appears $(length(hits)) times on one curve"))
    sign(hits[1])
end

"""α–β crossings with Lescop signs. Roleless endpoints are not crossings of the diagram."""
function _lescop_crossing_table(words, roles, alpha_names, beta_names)
    alpha_of = Dict(name => i for (i, name) in enumerate(alpha_names))
    beta_of = Dict(name => i for (i, name) in enumerate(beta_names))
    role = Dict(name => get(roles, name, :none) for name in keys(words))
    labels = sort(collect(keys(_lc_rotation_occ(words))))
    rows = NamedTuple{(:label, :sigma, :alpha, :beta, :alpha_id, :beta_id),
                       Tuple{Int,Int,String,String,Int,Int}}[]
    for label in labels
        owners = Pair{String,Int}[]
        for (name, word) in words
            passage = 0
            for x in word
                abs(x) == label || continue
                passage = sign(x)
                break
            end
            passage == 0 && continue
            push!(owners, name => passage)
        end
        length(owners) == 2 || throw(ArgumentError("label $label does not have two owners"))
        fam = [get(role, name, :none) for (name, _) in owners]
        count(in((:alpha, :beta)), fam) == 2 || continue
        fam[1] == fam[2] && throw(ArgumentError("label $label joins two $(fam[1]) curves"))
        alpha_at = findfirst(==(:alpha), fam)
        beta_at = findfirst(==(:beta), fam)
        alpha_name, alpha_passage = owners[alpha_at]
        beta_name, beta_passage = owners[beta_at]
        sigma = -alpha_passage
        sigma == beta_passage ||
            throw(ArgumentError("passage signs of label $label are not opposite"))
        push!(rows, (label=label, sigma=sigma, alpha=alpha_name, beta=beta_name,
                     alpha_id=alpha_of[alpha_name], beta_id=beta_of[beta_name]))
    end
    rows
end

function _lescop_merge_roleless!(parent, words, roles, corner_face)
    function find(a)
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        a
    end
    for (name, word) in words
        get(roles, name, :none) in (:alpha, :beta) && continue
        n = length(word)
        n == 0 && continue
        for i in 1:n
            j = mod1(i + 1, n)
            left = (abs(word[i]), _lc_sign_char(word[i]), "out")
            right = (abs(word[j]), _lc_sign_char(word[j]), "in")
            f1 = get(corner_face, left, nothing)
            f2 = get(corner_face, right, nothing)
            (f1 === nothing || f2 === nothing) &&
                throw(ArgumentError("roleless arc of $name is missing a face"))
            parent[find(f2)] = find(f1)
        end
    end
    find
end

function _lescop_face_regions(words, roles)
    active = _lc_nonempty_words(words)
    faces, corner_face = _lc_faces(active)
    isempty(faces) && throw(ArgumentError("LC has no face to carry an exterior point"))
    parent = Dict(face.id => face.id for face in faces)
    find = _lescop_merge_roleless!(parent, active, roles, corner_face)
    groups = Dict{String,Vector{String}}()
    for face in faces
        push!(get!(groups, find(face.id), String[]), face.id)
    end
    number(id) = parse(Int, id[2:end])
    roots = sort(collect(keys(groups)); by=root -> minimum(number, groups[root]))
    face_region = Dict{String,Int}()
    for (index, root) in enumerate(roots)
        for id in groups[root]
            face_region[id] = index
        end
    end
    face_region, corner_face, length(roots)
end

"""Faces along one side of the arc from an α–β crossing to the next one."""
function _lescop_arc_regions(word, start_at, stop_at, ab, corner_face, face_region)
    n = length(word)
    left = Set{Int}()
    right = Set{Int}()
    index = start_at
    for _ in 1:n
        nxt = mod1(index + 1, n)
        cout = (abs(word[index]), _lc_sign_char(word[index]), "out")
        cin = (abs(word[nxt]), _lc_sign_char(word[nxt]), "in")
        push!(left, face_region[corner_face[cout]])
        push!(right, face_region[corner_face[cin]])
        index = nxt
        index == stop_at && break
    end
    index == stop_at || throw(ArgumentError("alpha-beta arc did not close"))
    (length(left) == 1 && length(right) == 1) ||
        throw(ArgumentError("forgetting K separated the two sides of one alpha-beta arc"))
    only(left), only(right)
end

function _lescop_visits(word, ab, id_of)
    Int[id_of[abs(x)] for x in word if abs(x) in ab]
end

"""
    hc_to_cellular(hc::HeegaardCode) -> HCLescopView

Carry crossings, cyclic order, Lescop signs, a matching's labels, and
complement regions from a resolved HC to the cellular diagram read by
`lescop_theta`. Roleless K arcs are deleted and the LC faces they separate
are one complement region. A K curve with an alpha or beta role stays.
`face_region` sends an LC face id to that complement-region index.
"""
function hc_to_cellular(hc::HeegaardCode)
    validate_hc_input(hc)
    _lc_validate_surface!(hc.lc)
    _lescop_check_reference!(hc.lc)
    words = hc.lc.words
    alpha_names = hc_curve_names(hc, :alpha)
    beta_names = hc_curve_names(hc, :beta)
    (isempty(alpha_names) || isempty(beta_names)) &&
        throw(ArgumentError("Lescop's formula needs both curve families"))
    rows = _lescop_crossing_table(words, hc.roles, alpha_names, beta_names)
    isempty(rows) && throw(ArgumentError("no alpha-beta crossing"))
    ab = Set(row.label for row in rows)
    id_of = Dict(row.label => index for (index, row) in enumerate(rows))
    alpha = [Curve(:alpha, index, _lescop_visits(words[name], ab, id_of); label=name)
             for (index, name) in enumerate(alpha_names)]
    beta = [Curve(:beta, index, _lescop_visits(words[name], ab, id_of); label=name)
            for (index, name) in enumerate(beta_names)]
    any(curve -> isempty(curve.visits), vcat(alpha, beta)) &&
        throw(ArgumentError("a role curve meets the other family nowhere"))
    points = [Intersection(id_of[row.label], row.alpha_id, row.beta_id, row.sigma;
                           source_label=string(row.label)) for row in rows]
    code = heegaard_code(alpha, beta, points)
    face_region, corner_face, nregion = _lescop_face_regions(words, hc.roles)
    side = Dict{Tuple{Symbol,Int,Int,Int},Tuple{Int,Int}}()
    for (family, curves) in ((:alpha, alpha_names), (:beta, beta_names))
        for (index, name) in enumerate(curves)
            word = words[name]
            positions = [i for (i, x) in enumerate(word) if abs(x) in ab]
            for k in eachindex(positions)
                start_at = positions[k]
                stop_at = positions[mod1(k + 1, length(positions))]
                left, right = _lescop_arc_regions(word, start_at, stop_at, ab, corner_face, face_region)
                from = id_of[abs(word[start_at])]
                to = id_of[abs(word[stop_at])]
                side[(family, index, from, to)] = (left, right)
            end
        end
    end
    graph = ribbon_graph(code)
    cycles = boundary_cycles(graph)
    blocks = [Int[] for _ in 1:nregion]
    for (cycle_index, cycle) in enumerate(cycles)
        regions = Set{Int}()
        for dart_id in cycle
            dart = graph.darts[dart_id]
            edge = graph.edges[dart.edge]
            left, right = side[(edge.family, edge.curve_id, edge.from, edge.to)]
            # The ribbon walk turns clockwise, so the interior lies to the right
            # of a forward dart. That is the face of the backward LC corner.
            push!(regions, isodd(dart_id) ? right : left)
        end
        length(regions) == 1 ||
            throw(ArgumentError("one boundary cycle meets more than one complement region"))
        push!(blocks[only(regions)], cycle_index)
    end
    for (index, block) in enumerate(blocks)
        isempty(block) &&
            throw(ArgumentError("forgetting K dropped every boundary of complement region $index"))
    end
    diagram = certify_heegaard(decorate(code, blocks))
    HCLescopView(diagram, face_region, Dict(id => row.label for (id, row) in enumerate(rows)))
end

function _lescop_exterior_face(hc::HeegaardCode, exterior)
    exterior isa Integer && throw(ArgumentError(
        "the exterior point is an LC face id; that number is not the complement-region index"))
    key = String(exterior)
    if haskey(hc.points, key)
        key = hc.points[key]
        (startswith(key, "F") && all(isdigit, SubString(key, 2))) ||
            throw(ArgumentError("the exterior point must lie in an LC face"))
    end
    (startswith(key, "F") && all(isdigit, SubString(key, 2))) ||
        throw(ArgumentError("exterior face must be an LC face id such as F1"))
    key
end

"""
    hc_lescop_theta(hc, matching_labels, exterior_face)

Evaluate `(ell2, sl, e, theta)` from a resolved HC. `matching_labels` are the
LC labels of one crossing on each alpha and each beta. `exterior_face` is the
id of an LC face in the complement region of the exterior point, or the name
of an HC point supported on such a face. The LC face number is mapped to the
complement region and is not passed through.
"""
function hc_lescop_theta(hc::HeegaardCode, matching_labels, exterior)
    view = hc_to_cellular(hc)
    face = _lescop_exterior_face(hc, exterior)
    haskey(view.face_region, face) ||
        throw(ArgumentError("exterior face $face is not a face of this diagram"))
    wanted = Set(Int(label) for label in matching_labels)
    ids = Int[id for (id, label) in view.label_of_id if label in wanted]
    length(ids) == length(wanted) ||
        throw(ArgumentError("a matching label is not an alpha-beta crossing"))
    lescop_theta(view.diagram, ids, view.face_region[face])
end

"""
    hc_lescop_p1_difference(hc, matching_labels, from_label, to_label)

Return `(lk, delta_p1)` for the oriented beta edge from `from_label` to
`to_label` on the cellular diagram of a resolved HC. `delta_p1` is the
difference of p₁ between the two sides of that edge. The labels are LC
labels of alpha-beta crossings. Theta is not consulted.
"""
function hc_lescop_p1_difference(hc::HeegaardCode, matching_labels, from_label, to_label)
    view = hc_to_cellular(hc)
    id_of = Dict(label => id for (id, label) in view.label_of_id)
    function crossing(label)
        id = get(id_of, Int(label), nothing)
        id === nothing && throw(ArgumentError("label $label is not an alpha-beta crossing"))
        id
    end
    wanted = [crossing(label) for label in matching_labels]
    lescop_p1_difference(view.diagram, wanted, crossing(from_label), crossing(to_label))
end

function _lescop_matching_ids(hc::HeegaardCode, matching_labels)
    view = hc_to_cellular(hc)
    id_of = Dict(label => id for (id, label) in view.label_of_id)
    ids = Int[]
    for label in matching_labels
        id = get(id_of, Int(label), nothing)
        id === nothing && throw(ArgumentError("label $label is not an alpha-beta crossing"))
        push!(ids, id)
    end
    view, ids
end

"""
    hc_lescop_epsilon(hc, matching_labels, other_labels)

Return ε of two matchings on the cellular diagram of a resolved HC. The
labels are LC labels of one alpha-beta crossing on each curve. Θ and p₁
are not consulted.
"""
function hc_lescop_epsilon(hc::HeegaardCode, matching_labels, other_labels)
    view, first = _lescop_matching_ids(hc, matching_labels)
    _, second = _lescop_matching_ids(hc, other_labels)
    lescop_epsilon(view.diagram, first, second)
end

"""
    hc_lescop_spinc_label(hc, matching_labels, exterior)

Return the Spin^c label of one matching and one exterior point on the
cellular diagram of a resolved HC. `exterior` is an LC face id or an HC
point supported on such a face. Θ and p₁ are not consulted.
"""
function hc_lescop_spinc_label(hc::HeegaardCode, matching_labels, exterior)
    view, ids = _lescop_matching_ids(hc, matching_labels)
    face = _lescop_exterior_face(hc, exterior)
    haskey(view.face_region, face) ||
        throw(ArgumentError("exterior face $face is not a face of this diagram"))
    lescop_spinc_label(view.diagram, ids, view.face_region[face])
end
