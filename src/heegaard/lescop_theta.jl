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
    lescop_intersection_terms(H::HeegaardDiagram, matching)

Return the exact intersection terms (ell2, sl) for a certified rational
homology sphere diagram. The matching gives the crossing IDs of a perfect
matching between the alpha and beta systems.
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
    lescop_theta(H::HeegaardDiagram, matching, exterior_region::Integer)

Evaluate (ell2, sl, e, theta) exactly from a rational homology sphere
Heegaard diagram, a matching of crossing IDs, and the index of the
complement region containing the exterior point. The rotation term is
recovered from planar face turning numbers and a dual potential. The
result uses Lescop's infinity combing convention, with the Heegaard
surface oriented as the boundary of the alpha handlebody.
"""
function lescop_theta(H::HeegaardDiagram, matching, exterior_region::Integer)
    data = _lescop_data(H, matching)
    terms = _lescop_intersections(data)
    (; ell2, sl) = terms
    e = _lescop_rotation(H, data, terms, exterior_region)
    (; ell2, sl, e, theta=ell2 + sl - e)
end
