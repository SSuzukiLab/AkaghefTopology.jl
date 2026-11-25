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
    hc_to_cellular(hc::HeegaardCode) -> Topology.HeegaardCodes.HCLescopView

Return the cellular Heegaard view used by the HC Lescop calculations.

# Arguments
`hc` satisfies the [HC and reference-surface conditions](../specs/heegaard-code.md).
Input and surface validation run. Both curve families must exist and every
role curve must meet the other family; reference-word, passage, region and
cellular certification failures raise `ArgumentError` or their existing
validation exceptions. These checks do not certify complete cut systems.

# Returns
A `Topology.HeegaardCodes.HCLescopView` with the legacy `diagram`, an LC-face-to-
complement-region `face_region` map, and cellular-crossing-to-LC-label `label_of_id`.

# Guarantees
1. The view retains alpha/beta crossing order and converts signs under
   [the ordered-pairing translation](../specs/heegaard-code.md#intersection-sign-translation).
   [proved/agent: by construction]
2. Its region map groups LC faces across forgotten roleless reference arcs;
   K curves carrying an alpha/beta role remain. [proved/agent: by construction]
3. The input is not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
The view does not by itself certify a rational homology sphere, a complete
cut system, or the Lescop infinity-combing normalization. An LC face ID is
not its complement-region number.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §§5, 9.1 — LC faces, complement regions and exterior data.
cf: lib/topology.jl/docs/specs/heegaard-code.md#intersection-sign-translation — Lescop ordered signs.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1,Dict("l1"=>[-1,-4],"m1"=>[1,-2,-3],"b"=>[4,2,3]),Dict(),Dict(),Dict(),Dict(),nothing); v = H.hc_to_cellular(H.resolve_hc(lc,Dict("m1"=>:alpha,"b"=>:beta)));

julia> Tuple(sort!(collect(values(v.label_of_id))))
(2, 3)
```

Role: convert
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
    hc_lescop_theta(hc::HeegaardCode, matching_labels, exterior) -> NamedTuple

Return Lescop's matching- and exterior-dependent theta calculation on a resolved HC.

# Arguments
`hc` must be accepted by `hc_to_cellular`. LC crossing labels specify a
[matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語);
repeated labels are collapsed by this entry. `exterior` names an LC face or
an HC point supported on one; integer complement-region indices are rejected.
Unknown labels or faces, incomplete cut systems, a singular intersection
matrix and invalid rotation data raise `ArgumentError`. Cellular-conversion
and exact-arithmetic errors propagate.

# Returns
A named tuple of exact rational `ell2`, `sl`, `e` and `theta` values for
the converted cellular diagram, matching and exterior region.

# Guarantees
1. The calculation uses the complement region containing the specified
   LC face; `theta == ell2 + sl - e`. [proved/agent: by construction]
2. On the registered Lescop D1 HC fixture the result is
   `(ell2=0, sl=0, e=-1//4, theta=1//4)`.
   [verified-computationally: test/heegaard/lescop_theta.jl, D1 HC conversion]
3. Neither input nor matching labels are modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
This value depends on the matching and exterior choice and is not an
unqualified manifold invariant. The [Lescop coordinate specification](../specs/heegaard-terminology.md#lescop-coordinates)
states the combing interpretation; a general comparison with the geometric
Lescop normalization is not established by finite examples.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — Lescop theta and the exterior region.
cf: [[Les14]] §§1.3–1.5 — intended diagram formula.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-4], "m1"=>[1,-2,-3], "b"=>[4,2,3]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc,Dict("m1"=>:alpha,"b"=>:beta));

julia> H.hc_lescop_theta(hc,[2],"F3")
(ell2 = 0//1, sl = 0//1, e = -1//4, theta = 1//4)
```

Role: invariant
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
    hc_lescop_p1_difference(hc::HeegaardCode, matching_labels, from_label, to_label) -> NamedTuple

Return the linking-number and p₁-difference calculation across an oriented beta edge.

# Arguments
`hc` must be accepted by `hc_to_cellular`; labels refer to LC crossings.
`matching_labels` specifies a [matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語).
Checked: labels are alpha–beta crossings, endpoints are successive along one
oriented beta curve, both cut systems are complete and the intersection
matrix is invertible. Violations raise `ArgumentError` describing the failed
condition; cellular-conversion and arithmetic errors propagate.

# Returns
A named tuple `lk` and `delta_p1`, both exact rational values, for the
transported cellular diagram, matching and edge.

# Guarantees
1. The output is the cellular beta-edge calculation with the specified
   LC labels, and `delta_p1 == 8lk`. [proved/agent: by construction]
2. Neither the HC nor the matching list is modified.
   [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
An absolute p₁ is not returned. The [Lescop coordinate specification](../specs/heegaard-terminology.md#lescop-coordinates)
states the side correspondence and combing interpretation.
This wrapper does not establish the general geometric p₁ identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — matching-dependent diagram invariants.
cf: [[Les14]] Theorem 3.4 and Proposition 5.1 — intended p₁ difference and linking relation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-4], "m1"=>[1,-2,-3], "b"=>[4,2,3]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc,Dict("m1"=>:alpha,"b"=>:beta));

julia> r = H.hc_lescop_p1_difference(hc,[2],2,3); r.delta_p1 == 8r.lk
true
```

Role: invariant
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
    hc_lescop_epsilon(hc::HeegaardCode, matching_labels, other_labels) -> Vector{Int}

Return the homology-coordinate difference of two matchings on an HC's cellular diagram.

# Arguments
`hc` must be accepted by `hc_to_cellular`. The two label lists specify
[matchings](../specs/heegaard-terminology.md#91-群と不変量で使う用語) using LC crossing labels.
Unknown or non-alpha–beta labels, invalid matchings, singular intersection
matrices and oversized homology-coordinate searches raise `ArgumentError`.
Other cellular-conversion errors propagate. The caller supplies a valid
rational-homology-sphere Heegaard diagram for a manifold interpretation.

# Returns
An integer coordinate vector for the difference class on the converted
cellular diagram.

# Guarantees
1. The result equals the cellular matching-difference calculation after
   transporting LC crossing labels. [proved/agent: by construction]
2. Equal matchings give the zero class on the tested genus-one diagrams.
   [verified-computationally: test/heegaard/lescop_theta.jl, genus-one matching comparisons]
3. The input and label lists are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
The [Lescop coordinate specification](../specs/heegaard-terminology.md#lescop-coordinates)
states the representative and affine reference choices. Equality of this vector across unrelated
coordinate bases is not asserted. This wrapper does not prove a general
geometric identification with the corresponding manifold class.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — matchings and diagram invariants.
cf: [[OS04a]] subsection subsec:Disks — intended matching-difference class.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-4], "m1"=>[1,-2,-3], "b"=>[4,2,3]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc,Dict("m1"=>:alpha,"b"=>:beta));

julia> Tuple(H.hc_lescop_epsilon(hc,[2],[2]))
(0,)
```

Role: invariant
"""
function hc_lescop_epsilon(hc::HeegaardCode, matching_labels, other_labels)
    view, first = _lescop_matching_ids(hc, matching_labels)
    _, second = _lescop_matching_ids(hc, other_labels)
    lescop_epsilon(view.diagram, first, second)
end

"""
    hc_lescop_spinc_label(hc::HeegaardCode, matching_labels, exterior) -> Vector{Int}

Return the cellular Spinᶜ-coordinate label associated with a matching and exterior point.

# Arguments
`hc` must be accepted by `hc_to_cellular`. LC crossing labels specify a
[matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語).
`exterior` names an LC face or an HC point supported on an LC face; integer
region indices are rejected with `ArgumentError`. Unknown labels, invalid
matchings or faces, singular intersection data, oversized homology searches
and nonunique basepoint shifts raise `ArgumentError`; other conversion
errors propagate. A valid manifold diagram is the caller's responsibility.

# Returns
A vector of integer coordinates for the converted cellular matching and
exterior complement region.

# Guarantees
1. The exterior face is transported to its actual cellular complement
   region, and matching labels to cellular crossings.
   [proved/agent: by construction]
2. The result is that cellular record's Spinᶜ-label calculation.
   [proved/agent: by construction]
3. Neither input data nor the matching list is modified.
   [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
The [Lescop coordinate specification](../specs/heegaard-terminology.md#lescop-coordinates)
states the affine origin and basepoint transport. A vector alone is not an absolute,
basis-independent Spinᶜ structure or a general geometric certification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — exterior region versus LC face.
cf: [[OS04a]] subsection subsec:SpinCStructures, lemma VarySpinC — intended basepoint-dependent Spinᶜ interpretation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-4], "m1"=>[1,-2,-3], "b"=>[4,2,3]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc,Dict("m1"=>:alpha,"b"=>:beta));

julia> try H.hc_lescop_spinc_label(hc,[2],3) catch e; e isa ArgumentError end
true
```

Role: invariant
"""
function hc_lescop_spinc_label(hc::HeegaardCode, matching_labels, exterior)
    view, ids = _lescop_matching_ids(hc, matching_labels)
    face = _lescop_exterior_face(hc, exterior)
    haskey(view.face_region, face) ||
        throw(ArgumentError("exterior face $face is not a face of this diagram"))
    lescop_spinc_label(view.diagram, ids, view.face_region[face])
end
