# Source: Lescop, https://arxiv.org/abs/1402.2261, §§1.3–1.5.
# groups.jl::intersection_matrix has alpha rows and beta columns.
function _lescop_data(H::HeegaardDiagram, matching)
    report = compression_boundary_report(H)
    report.alpha_complete && report.beta_complete ||
        throw(ArgumentError("Lescop's formula requires two complete cut systems"))
    code = _undecorated(H)
    points = _intersection_map(code)
    alpha_match = Dict{Int,Int}()
    beta_match = Dict{Int,Int}()
    for id in matching
        point = get(points, id, nothing)
        point === nothing && throw(ArgumentError("matching crossing $id is absent"))
        (haskey(alpha_match, point.alpha) || haskey(beta_match, point.beta)) &&
            throw(ArgumentError("matching must meet each alpha and beta exactly once"))
        alpha_match[point.alpha] = id
        beta_match[point.beta] = id
    end
    length(alpha_match) == length(code.alpha) ||
        throw(ArgumentError("matching must meet each alpha and beta exactly once"))
    J = try
        inv(Rational{BigInt}.(intersection_matrix(H)))
    catch err
        err isa SingularException || rethrow()
        throw(ArgumentError("Lescop's formula requires a rational homology sphere"))
    end

    intersections = code.intersections
    index = Dict(point.id => c for (c, point) in enumerate(intersections))
    alpha_index = Dict(curve.id => i for (i, curve) in enumerate(code.alpha))
    beta_index = Dict(curve.id => j for (j, curve) in enumerate(code.beta))
    a = [alpha_index[p.alpha] for p in intersections]
    b = [beta_index[p.beta] for p in intersections]
    signs = BigInt[p.sign for p in intersections]
    u = [J[b[c],a[c]] * signs[c] for c in eachindex(intersections)]
    (; code, alpha_match, beta_match, J, index, a, b, signs, u)
end

# Entries are twice the arc weights. Both endpoints have half weight,
# and an arc from a matching crossing to itself is empty (Lescop §1.3).
function _lescop_arcs(curves, matching, point_index)
    n = length(point_index)
    arcs = zeros(Int, n, n)
    for curve in curves
        visits = curve.visits
        start = findfirst(==(matching[curve.id]), visits)
        first_point = point_index[visits[start]]
        previous = Int[]
        for offset in 1:(length(visits) - 1)
            current = point_index[visits[mod1(start + offset, length(visits))]]
            arcs[first_point, current] = 1
            arcs[current, current] = 1
            arcs[previous, current] .= 2
            push!(previous, current)
        end
    end
    arcs
end

function _lescop_intersections(data)
    (; code, alpha_match, beta_match, J, index, a, b, signs, u) = data
    n = length(signs)
    A = _lescop_arcs(code.alpha, alpha_match, index)
    B = _lescop_arcs(code.beta, beta_match, index)
    alpha_beta = zeros(BigInt, n, length(code.beta))
    beta_alpha = zeros(BigInt, length(code.alpha), n)
    for c in 1:n, d in 1:n
        alpha_beta[c,b[d]] += A[d,c] * signs[d]
        beta_alpha[a[d],c] += B[d,c] * signs[d]
    end
    correction = (alpha_beta * J * beta_alpha) ./ 4
    raw = transpose(A) * (signs .* B)
    ell2 = sl = zero(eltype(J))
    for c in 1:n, d in 1:n
        kernel = raw[c,d] // 4 - correction[c,d]
        ell2 += (J[b[c],a[d]] * J[b[d],a[c]] * signs[c] * signs[d]) * kernel
        c == d && (ell2 -= u[c] * kernel)
        sl += u[c] * u[d] * kernel
    end
    (; ell2, sl, beta_alpha)
end

"""
    lescop_intersection_terms(H::CellularHeegaardDiagram, matching) -> NamedTuple

Compute the exact intersection terms of the implemented Lescop diagram formula.

# Arguments
`H` is a legacy cellular diagram for a rational homology sphere and `matching` is a [matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語) of crossing IDs. Both complete cut systems, the matching and nonsingularity are checked; failed checks raise `ArgumentError`. Incidence and exact-linear-algebra errors propagate.

# Returns
A named tuple `(ell2, sl)` of exact rational numbers.

# Guarantees
1. The result is the pair of diagram-and-matching intersection terms used by `lescop_theta`. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
A general identification with the geometric invariant requires the matching, homology-coordinate and combing conventions to agree with the source. The current specification does not establish that full identification; the result is the implemented diagram quantity, not an independently certified manifold invariant.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — Lescop's diagram-dependent quantity and matching.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> D = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2]); cellular=true));

julia> r = H.lescop_intersection_terms(D, [1]); (r.ell2, r.sl)
(0//1, 0//1)
```

Role: invariant
"""
function lescop_intersection_terms(H::HeegaardDiagram, matching)
    terms = _lescop_intersections(_lescop_data(H, matching))
    (; terms.ell2, terms.sl)
end

# Reconstruction from Lescop §1.5, Fig. 2, using planar turning number.
# Place the matching port of each alpha circle at the upward tangent and
# cluster the others, in cyclic order, near the downward tangent. In the
# cluster limit each alpha interval turns a half-turn per matched endpoint,
# and all other intervals turn zero. Each region traverses the holes
# clockwise, with interior on the left. Each crossing corner turns +1/4.
# The total turn of the oriented boundary equals the Euler characteristic
# after removing the disk at w. Its new outer boundary contributes +1.
# The cluster limit preserves the sum of endpoint and arc turns.
function _lescop_rotation(H, data, terms, exterior_region)
    (; code, alpha_match, beta_match, J, index, u) = data
    blocks = _surface_partition(H).boundary_cycle_blocks
    exterior_region in eachindex(blocks) ||
        throw(ArgumentError("exterior point must specify a complement region index"))
    graph = ribbon_graph(H)
    cycles = boundary_cycles(graph)
    owner = _ribbon_region_owners(graph, blocks)
    matched = Set(values(alpha_match))
    b = [(2 - length(block) - 2(r == exterior_region)) // big(1)
         for (r, block) in enumerate(blocks)]
    for (r, block) in enumerate(blocks), c in block
        b[r] -= length(cycles[c]) // big(4)
        for id in cycles[c]
            dart = graph.darts[id]
            graph.edges[dart.edge].family == :alpha || continue
            b[r] += ((dart.origin in matched) + (dart.target in matched)) // big(2)
        end
    end

    # Coefficient of each oriented beta edge turn in Lescop's e. Its
    # prefix contribution uses subsequent crossings before returning to
    # the matching. The full beta degree enters every correction term.
    correction = J * (terms.beta_alpha * u) ./ 2
    weight = Dict{Int,Rational{BigInt}}()
    for (j, curve) in enumerate(code.beta)
        visits = curve.visits
        start = findfirst(==(beta_match[curve.id]), visits)
        suffix = -correction[j]
        for offset in (length(visits) - 1):-1:0
            id = visits[mod1(start + offset, length(visits))]
            weight[id] = suffix
            suffix += u[index[id]]
        end
    end

    # The dual incidence column is +1 on the left and -1 on the right.
    # The clockwise predecessor in ribbon.jl traces the left-hand region.
    # A dual potential certifies that the e coefficient vector annihilates
    # the kernel of the face equations, so e needs no chosen individual turns.
    T = Rational{BigInt}
    adjacency = [Tuple{Int,T}[] for _ in blocks]
    for edge in graph.edges
        edge.family == :beta || continue
        left = owner[2edge.id - 1]
        right = owner[2edge.id]
        w = weight[edge.from]
        push!(adjacency[left], (right, -w))
        push!(adjacency[right], (left, w))
    end
    potential = zeros(T, length(blocks))
    seen = falses(length(blocks))
    for root in eachindex(blocks)
        seen[root] && continue
        queue = [root]
        seen[root] = true
        for r in queue
            for (s, delta) in adjacency[r]
                value = potential[r] + delta
                if seen[s]
                    potential[s] == value ||
                        throw(ArgumentError("face turns do not determine Lescop's rotation term"))
                else
                    potential[s] = value
                    seen[s] = true
                    push!(queue, s)
                end
            end
        end
        iszero(sum(b[queue])) || throw(ArgumentError("inconsistent face-turn equations"))
    end
    sum(potential .* b)
end

"""
    lescop_theta(H::CellularHeegaardDiagram, matching, exterior_region::Integer) -> NamedTuple

Evaluate the exact implemented Lescop diagram expression at an exterior region.

# Arguments
`H` is a legacy cellular diagram for a rational homology sphere and `matching` is a [matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語) of crossing IDs. Both complete cut systems, the matching and nonsingularity are checked; failed checks raise `ArgumentError`. Incidence and exact-linear-algebra errors propagate. `exterior_region` indexes a complement region; an invalid index or inconsistent face-turn equations raises `ArgumentError`.

# Returns
A named tuple `(ell2, sl, e, theta)` of exact rational numbers.

# Guarantees
1. The intersection terms agree with `lescop_intersection_terms` and `theta == ell2 + sl - e`; `e` is the implemented exterior-region rotation correction. [proved/agent: by construction]
2. The value retains its diagram, matching and exterior-region dependence. Inputs are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
A general identification with the geometric invariant requires the matching, homology-coordinate and combing conventions to agree with the source. The current specification does not establish that full identification; the result is the implemented diagram quantity, not an independently certified manifold invariant. Independence of the exterior region or matching is not asserted.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — Lescop's Θ and its selected diagram data.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> D = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2]); cellular=true));

julia> r = H.lescop_theta(D, [1], 2); (r.ell2, r.sl, r.e, r.theta)
(0//1, 0//1, -1//4, 1//4)
```

Role: invariant
"""
function lescop_theta(H::HeegaardDiagram, matching, exterior_region::Integer)
    data = _lescop_data(H, matching)
    terms = _lescop_intersections(data)
    (; ell2, sl) = terms
    e = _lescop_rotation(H, data, terms, exterior_region)
    (; ell2, sl, e, theta=ell2 + sl - e)
end

# Les14 Theorem 3.4 (thmbousillete): the p1 difference of two infinity
# combings that differ by the exterior point is 8 times the linking
# number of Proposition 5.1 (propevallinkbis). The arc weights are the
# half-bracket of section 1.3. This function does not read Theta.
function _lescop_half_arc(visits, base, target)
    weights = Dict{Int,Rational{BigInt}}()
    edges = Tuple{Int,Int}[]
    target == base && return weights, edges
    start = findfirst(==(base), visits)
    start === nothing && throw(ArgumentError("crossing $base is not on its curve"))
    n = length(visits)
    prev = base
    for _ in 1:n
        nxt = visits[mod1(start + (length(edges) + 1), n)]
        push!(edges, (prev, nxt))
        if nxt == target
            weights[base] = 1//2
            weights[target] = get(weights, target, 0//1) + 1//2
            return weights, edges
        end
        weights[nxt] = get(weights, nxt, 0//1) + 1
        prev = nxt
    end
    throw(ArgumentError("crossing $target is not on the curve of crossing $base"))
end

"""
    lescop_p1_difference(H::CellularHeegaardDiagram, matching, from_id::Integer, to_id::Integer) -> NamedTuple

Compute the implemented linking term and p₁ change across an oriented beta edge.

# Arguments
`H` is a legacy cellular diagram for a rational homology sphere and `matching` is a [matching](../specs/heegaard-terminology.md#91-群と不変量で使う用語) of crossing IDs. Both complete cut systems, the matching and nonsingularity are checked; failed checks raise `ArgumentError`. Incidence and exact-linear-algebra errors propagate. `from_id` and `to_id` must be consecutive crossing IDs in the directed beta word; absence, different beta curves or a wrong successor raises `ArgumentError`.

# Returns
A named tuple `(lk, delta_p1)` of exact rational numbers.

# Guarantees
1. The returned `delta_p1` is eight times the returned `lk`, for the specified directed edge and matching. [proved/agent: by construction]
2. The quantity is the implemented edge-crossing contribution intended for the infinity-combing p₁ difference; inputs are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
A general identification with the geometric invariant requires the matching, homology-coordinate and combing conventions to agree with the source. The current specification does not establish that full identification; the result is the implemented diagram quantity, not an independently certified manifold invariant. The edge-side translation to a geometric left-to-right basepoint motion is not separately certified by this call.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — diagram and matching dependence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> D = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2]); cellular=true));

julia> r = H.lescop_p1_difference(D, [1], 1, 2); (r.lk, r.delta_p1)
(-1//4, -2//1)
```

Role: invariant
"""
function lescop_p1_difference(H::HeegaardDiagram, matching, from_id::Integer, to_id::Integer)
    data = _lescop_data(H, matching)
    (; code, beta_match, J, index, u) = data
    points = Dict(p.id => p for p in code.intersections)
    (haskey(points, from_id) && haskey(points, to_id)) ||
        throw(ArgumentError("the beta edge must join two crossings of the diagram"))
    points[from_id].beta == points[to_id].beta ||
        throw(ArgumentError("the two crossings do not lie on one beta curve"))
    owner = only(c for c in code.beta if c.id == points[from_id].beta)
    start = findfirst(==(from_id), owner.visits)
    start === nothing && throw(ArgumentError("crossing $from_id is not on its beta curve"))
    owner.visits[mod1(start + 1, length(owner.visits))] == to_id ||
        throw(ArgumentError("the oriented beta edge does not run from $from_id to $to_id"))

    alpha_index = Dict(c.id => i for (i, c) in enumerate(code.alpha))
    beta_index = Dict(c.id => j for (j, c) in enumerate(code.beta))
    match = Set(Int.(matching))
    k = Dict(id => (id in match ? 1//1 : 0//1) - u[index[id]] for id in keys(points))
    g = length(code.alpha)
    P = Dict(id => 0//1 for id in keys(points))
    A = [Dict(id => 0//1 for id in keys(points)) for _ in 1:g]
    for id in keys(points)
        curve = only(c for c in code.beta if c.id == points[id].beta)
        weights, edges = _lescop_half_arc(curve.visits, beta_match[curve.id], id)
        curve.id == owner.id && (from_id, to_id) in edges && (P[id] = 1//1)
        for (x, w) in weights
            A[alpha_index[points[x].alpha]][id] += w * points[x].sign
        end
    end
    lk = sum(k[id] * P[id] for id in keys(P); init=0//1)
    jowner = beta_index[owner.id]
    for i in 1:g
        acc = sum(k[id] * A[i][id] for id in keys(P); init=0//1)
        lk -= J[jowner, i] * acc
    end
    (; lk, delta_p1=8lk)
end

# OS04a, subsection subsec:Disks (anchor anonnd72-blk-0753) and
# subsection subsec:SpinCStructures, lemma:VarySpinC. The beta-edge
# class is the coboundary forced by that lemma's alpha formula.
# Theta and p1 are not read.

function _lescop_intersection_order(I::AbstractMatrix{<:Integer})
    size(I, 1) == size(I, 2) ||
        throw(ArgumentError("Spin^c labels require a rational homology sphere"))
    order = abs(det(BigInt.(I)))
    order == 0 && throw(ArgumentError("Spin^c labels require a rational homology sphere"))
    order
end

function _lescop_in_lattice(vector, I::AbstractMatrix{<:Integer})
    solution = Rational{BigInt}.(I) \ Rational{BigInt}.(collect(vector))
    all(value -> isone(denominator(value)), solution)
end

function _lescop_canonical_class(vector, I::AbstractMatrix{<:Integer})
    order = _lescop_intersection_order(I)
    genus = size(I, 1)
    genus == 1 && return Int[mod(vector[1], order)]
    box = order^genus
    box > 20000 && throw(ArgumentError("the homology box is too large for this label"))
    best = nothing
    coordinates = zeros(Int, genus)
    for n in 0:box-1
        remainder = n
        for index in genus:-1:1
            coordinates[index] = remainder % order
            remainder ÷= order
        end
        _lescop_in_lattice(coordinates .- vector, I) || continue
        if best === nothing || coordinates < best
            best = copy(coordinates)
        end
    end
    best === nothing && throw(ArgumentError("the homology class has no representative"))
    best
end

function _lescop_classes(I::AbstractMatrix{<:Integer})
    order = _lescop_intersection_order(I)
    genus = size(I, 1)
    box = order^genus
    box > 20000 && throw(ArgumentError("the homology box is too large for this label"))
    classes = Vector{Int}[]
    coordinates = zeros(Int, genus)
    for n in 0:box-1
        remainder = n
        for index in genus:-1:1
            coordinates[index] = remainder % order
            remainder ÷= order
        end
        canonical = _lescop_canonical_class(coordinates, I)
        canonical == coordinates && push!(classes, copy(coordinates))
    end
    length(classes) == order ||
        throw(ArgumentError("the homology representatives do not cover H1"))
    classes
end

function _lescop_match_maps(code, matching)
    points = _intersection_map(code)
    alpha_match = Dict{Int,Int}()
    beta_match = Dict{Int,Int}()
    for id in matching
        point = get(points, Int(id), nothing)
        point === nothing && throw(ArgumentError("matching crossing $id is absent"))
        (haskey(alpha_match, point.alpha) || haskey(beta_match, point.beta)) &&
            throw(ArgumentError("matching must meet each alpha and beta exactly once"))
        alpha_match[point.alpha] = Int(id)
        beta_match[point.beta] = Int(id)
    end
    (length(alpha_match) == length(code.alpha) && length(beta_match) == length(code.beta)) ||
        throw(ArgumentError("matching must meet each alpha and beta exactly once"))
    (; alpha_match, beta_match)
end

function _lescop_matchings(code)
    by_alpha = Dict(curve.id => Intersection[] for curve in code.alpha)
    for point in code.intersections
        haskey(by_alpha, point.alpha) && push!(by_alpha[point.alpha], point)
    end
    results = Vector{Int}[]
    chosen = Int[]
    used = Set{Int}()
    function walk(index)
        if index > length(code.alpha)
            length(used) == length(code.beta) && push!(results, copy(chosen))
            return
        end
        for point in by_alpha[code.alpha[index].id]
            point.beta in used && continue
            push!(used, point.beta)
            push!(chosen, point.id)
            walk(index + 1)
            pop!(chosen)
            pop!(used, point.beta)
        end
    end
    walk(1)
    isempty(results) && throw(ArgumentError("the diagram has no matching"))
    results
end

function _lescop_arrival_vector(code, first_matching, second_matching)
    points = _intersection_map(code)
    first = _lescop_match_maps(code, first_matching)
    second = _lescop_match_maps(code, second_matching)
    alpha_index = Dict(curve.id => index for (index, curve) in enumerate(code.alpha))
    arrival = zeros(Int, length(code.alpha))
    for curve in code.beta
        start_id = first.beta_match[curve.id]
        stop_id = second.beta_match[curve.id]
        start_id == stop_id && continue
        start = findfirst(==(start_id), curve.visits)
        start === nothing && throw(ArgumentError("crossing $start_id is not on its beta curve"))
        arrived = false
        for step in 1:(length(curve.visits) - 1)
            id = curve.visits[mod1(start + step, length(curve.visits))]
            point = points[id]
            arrival[alpha_index[point.alpha]] += point.sign
            if id == stop_id
                arrived = true
                break
            end
        end
        arrived || throw(ArgumentError("crossing $stop_id is not on its beta curve"))
    end
    arrival
end

"""
    lescop_epsilon(H::CellularHeegaardDiagram, matching, other) -> Vector{Int}

Compute the implemented homology-class difference of two matchings.

# Arguments
`matching` and `other` are matchings of crossing IDs in the [legacy diagram](../specs/heegaard-terminology.md#legacy-cellular-records). The caller supplies a rational homology sphere diagram. Matching validity, a nonzero determinant and the finite homology-label search bound are checked, with `ArgumentError` on failure. Arithmetic and incidence errors propagate.

# Returns
An integer coordinate vector for the class difference in the implemented intersection-matrix quotient.

# Guarantees
1. Equal matchings return the zero class whenever the call succeeds; the vector represents the matching-cycle class in the implemented lattice coordinates. [proved/agent: by construction]
2. The diagram and matching inputs are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
A general identification with the geometric invariant requires the matching, homology-coordinate and combing conventions to agree with the source. The current specification does not establish that full identification; the result is the implemented diagram quantity, not an independently certified manifold invariant.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — matching and diagram-dependent invariant inputs.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> D = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2]); cellular=true));

julia> Tuple(H.lescop_epsilon(D, [1], [1]))
(0,)
```

Role: invariant
"""
function lescop_epsilon(H::HeegaardDiagram, matching, other)
    code = _undecorated(H)
    I = intersection_matrix(H)
    arrival = _lescop_arrival_vector(code, matching, other)
    _lescop_canonical_class(-arrival, I)
end

function _lescop_edge_sides(H::HeegaardDiagram)
    graph = ribbon_graph(H)
    blocks = _surface_partition(H).boundary_cycle_blocks
    owner = _ribbon_region_owners(graph, blocks)
    graph, owner
end

function _lescop_region_shifts(H::HeegaardDiagram)
    code = _undecorated(H)
    I = intersection_matrix(H)
    genus = size(I, 1)
    graph, owner = _lescop_edge_sides(H)
    regions = length(_surface_partition(H).boundary_cycle_blocks)
    alpha_index = Dict(curve.id => index for (index, curve) in enumerate(code.alpha))
    beta_index = Dict(curve.id => index for (index, curve) in enumerate(code.beta))
    classes = _lescop_classes(I)
    beta_count = length(code.beta)
    solutions = Tuple{Vector{Vector{Int}},Vector{Vector{Int}}}[]
    candidates = collect(Iterators.product((classes for _ in 1:beta_count)...))
    for choice in candidates
        beta_class = Vector{Int}[collect(item) for item in choice]
        function edge_class(edge)
            if edge.family == :alpha
                class = zeros(Int, genus)
                class[alpha_index[edge.curve_id]] = -1
                return class
            end
            edge.family == :beta ||
                throw(ArgumentError("a ribbon edge is neither alpha nor beta"))
            beta_class[beta_index[edge.curve_id]]
        end
        potential = Union{Nothing,Vector{Int}}[nothing for _ in 1:regions]
        potential[1] = zeros(Int, genus)
        queue = [1]
        consistent = true
        for region in queue
            for edge in graph.edges
                left = owner[2edge.id - 1]
                right = owner[2edge.id]
                (left == region || right == region) || continue
                class = edge_class(edge)
                if left == right
                    _lescop_in_lattice(class, I) || (consistent = false)
                    continue
                end
                parent = region
                child = parent == left ? right : left
                step = parent == left ? class : -class
                expected = potential[parent] + step
                if potential[child] === nothing
                    potential[child] = expected
                    push!(queue, child)
                elseif !_lescop_in_lattice(potential[child] .- expected, I)
                    consistent = false
                end
            end
            consistent || break
        end
        consistent || continue
        any(value -> value === nothing, potential) && continue
        shifts = [_lescop_canonical_class(value, I) for value in potential]
        push!(solutions, (shifts, [_lescop_canonical_class(class, I) for class in beta_class]))
    end
    length(solutions) == 1 ||
        throw(ArgumentError("the basepoint class is not a unique coboundary"))
    shifts, beta_class = only(solutions)
    (; shifts, beta_class)
end

"""
    lescop_spinc_label(H::CellularHeegaardDiagram, matching, region::Integer) -> Vector{Int}

Compute the implemented affine homology label of a matching and exterior region.

# Arguments
`H` is a legacy cellular rational homology sphere diagram; `matching` uses crossing IDs and `region` indexes its complement partition. Invalid region or matching, singular homology data, exceeded finite-label bounds, or a nonunique basepoint coboundary raises `ArgumentError`. The caller establishes the complete Heegaard interpretation; this routine does not independently check it.

# Returns
An integer vector in the implemented finite homology quotient, with the routine's affine reference choices.

# Guarantees
1. The result combines the matching-class difference with the selected region's computed basepoint shift in the same homology quotient. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Conventions
Lescop and affine homology coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-terminology.md#lescop-coordinates).

# Not guaranteed
A general identification with the geometric invariant requires the matching, homology-coordinate and combing conventions to agree with the source. The current specification does not establish that full identification; the result is the implemented diagram quantity, not an independently certified manifold invariant. A coordinate label alone does not identify Spinᶜ structures across different diagrams or reference choices.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — matching and exterior-region inputs.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> D = H.certify_heegaard(H.decorate(H.heegaard_code([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2]); cellular=true));

julia> Tuple(H.lescop_spinc_label(D, [1], 2))
(1,)
```

Role: invariant
"""
function lescop_spinc_label(H::HeegaardDiagram, matching, region::Integer)
    code = _undecorated(H)
    blocks = _surface_partition(H).boundary_cycle_blocks
    1 <= region <= length(blocks) ||
        throw(ArgumentError("exterior point must specify a complement region index"))
    I = intersection_matrix(H)
    reference = first(sort(_lescop_matchings(code)))
    arrival = _lescop_arrival_vector(code, reference, matching)
    shifts, _ = _lescop_region_shifts(H)
    _lescop_canonical_class(-arrival .+ shifts[region], I)
end
