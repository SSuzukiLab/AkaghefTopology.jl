# cf:[[Kup91#^kup91-sec-0105]] §The invariant — Z is the contraction of the two circle families; a negative tangent basis interposes the antipode; adjustment: α carries Δ^(k-1)(right_cointegral), β carries right_integral of the word product, the antipode output is the β leg and the input is the α leg, and λ(Λ)≠1 is rescaled by pairing^{-g} rather than by (dim H)^{g-n_u-n_l}. eq26–eq28 unverified.
# cf:[[Kup91#^kup91-blk-0138]] — a recorded group-algebra consequence, not a code path: those weights are diagonal coproducts and identity-word products. This file uses only the structure constants.
const _KUPERBERG_INTEGRAL_ATOL = 1e-10

mutable struct _KuperbergRNG
    state::UInt64
end

function _kuperberg_rand!(rng::_KuperbergRNG)
    x = rng.state
    x = xor(x, x << 13)
    x = xor(x, x >> 7)
    x = xor(x, x << 17)
    rng.state = x == 0 ? 0x9e3779b97f4a7c15 : x
    rng.state
end

function _hc_role_curves(hc::Heegaard.HeegaardCode, role::Symbol)
    [name for (name, family) in hc.roles if family == role]
end

function _lc_label_owners(words::Dict{String, Vector{Int}})
    owners = Dict{Int, Dict{Symbol, String}}()
    for (name, word) in words
        for x in word
            k = abs(x)
            side = x > 0 ? :+ : :-
            d = get!(owners, k, Dict{Symbol, String}())
            haskey(d, side) && d[side] != name &&
                throw(ArgumentError("label $k has conflicting $(side) owners"))
            d[side] = name
        end
    end
    owners
end

_curve_role(roles::Dict{String, Symbol}, name::AbstractString) = get(roles, name, :none)

function _is_alpha_beta_label(owners::Dict, label::Int, roles::Dict{String, Symbol})
    haskey(owners, label) || return false
    sides = owners[label]
    haskey(sides, :+) && haskey(sides, :-) || return false
    ra, rb = _curve_role(roles, sides[:+]), _curve_role(roles, sides[:-])
    (ra == :alpha && rb == :beta) || (ra == :beta && rb == :alpha)
end

"""Return whether this code treats the crossing as positive (no antipode).

`docs/specs/lickorish-code.md` §3 defines passage sign only: the component
holding `−k` and the component holding `+k`, in that order, form a positive
frame. `docs/specs/heegaard-code.md` does not define a Kuperberg crossing
sign; it keeps the LC word signs. Those two specs therefore do not decide
which passage orientation inserts the antipode. eq26–eq28 unverified.
The current 10 algebras do not distinguish a global flip of this choice.
This predicate is the code's choice, not a value check: true when the `+`
owner is α and the `−` owner is β.
"""
function _ab_crossing_passage_positive(hc::Heegaard.HeegaardCode, label::Int)
    owners = _lc_label_owners(hc.lc.words)
    sides = owners[label]
    get(hc.roles, sides[:+], :none) == :alpha && get(hc.roles, sides[:-], :none) == :beta
end

function _filter_ab_word(word::Vector{Int}, owners, roles)
    [x for x in word if _is_alpha_beta_label(owners, abs(x), roles)]
end

"""α then β components, α∩β crossings only, orientation = §3 passage sign."""
function _heegaard_kuperberg_link(hc::Heegaard.HeegaardCode)
    owners = _lc_label_owners(hc.lc.words)
    alphas = sort(_hc_role_curves(hc, :alpha))
    betas = sort(_hc_role_curves(hc, :beta))
    curves = vcat(alphas, betas)
    gauss = [_filter_ab_word(hc.lc.words[name], owners, hc.roles) for name in curves]
    raw_labels = sort(unique(abs(x) for word in gauss for x in word))
    remap = Dict(label => index for (index, label) in enumerate(raw_labels))
    gauss = [[sign(x) * remap[abs(x)] for x in word] for word in gauss]
    orientation = [_ab_crossing_passage_positive(hc, label) ? 1 : -1 for label in raw_labels]
    link = VirtualLink()
    set_data!(link; rgauss=gauss, orientation=orientation)
    link
end

function _kuperberg_right_integral_pairing(algebra::FiniteHopfAlgebra)
    total = 0.0 + 0im
    for ((index,), integral_value) in algebra.right_integral.terms
        haskey(algebra.right_cointegral.terms, (index,)) || continue
        total += integral_value * algebra.right_cointegral.terms[(index,)]
    end
    total
end

function _require_involutory_kuperberg_algebra(algebra::FiniteHopfAlgebra)
    is_involutory(algebra) ||
        throw(DomainError(algebra, "Kuperberg on HC requires an involutory Hopf algebra"))
    pairing = _kuperberg_right_integral_pairing(algebra)
    abs(pairing) <= _KUPERBERG_INTEGRAL_ATOL &&
        throw(DomainError(algebra,
            "Kuperberg convention requires λ(Λ) ≠ 0 so the integral can be normalized"))
    algebra
end

"""
    heegaard_pi1_presentation(hc::HeegaardCodes.HeegaardCode) -> HeegaardCodes.GroupPresentation

Return the alpha-generator and beta-relator group associated with the Heegaard code's words.

# Arguments
Requires H8 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
The caller supplies coherent curve roles, words and passages. Conflicting
same-sign passage owners raise `ArgumentError` naming the label; missing
curve words or indices can raise `KeyError`. No HC validity, complete-family
count or cutting-system check occurs here. Non-alpha/beta passages are omitted.

# Returns
A `HeegaardCodes.GroupPresentation` with `generators`, `relators`, empty
`beta_ids`, and `metadata` identifying `lc_beta_words` and its passage
convention. The group and exponent-pair meaning are the
[Heegaard word group](../specs/hopf-invariants.md#heegaard-word-group-and-antipode-labels).

# Guarantees
1. The result presents the coordinate group defined in the linked specification, retaining the stored beta-word order without relator reduction. [proved/agent: by construction]
2. The intended identification is the fundamental group of the closed 3-manifold of the complete Heegaard splitting. [unverified: identify these alpha generators and signed beta relators with the complete-cut-system manifold of heegaard-koda-conversion.md §5.1 and the group-state constraints of Kup91 ^kup91-blk-0140, in both local passage orientations and arbitrary genus]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This routine does not certify the geometric identification in guarantee 2.
The missing comparison must identify the alpha generators with handlebody
generators and each beta relator with its attaching circle, including the
orientation of each crossing and the order along each beta circle. It must
cover complete cutting systems in arbitrary genus; the registered S3 and T3
homomorphism-count comparisons do not prove that identification. Malformed
inputs need not be rejected, and the presentation need not be minimal or simplified.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Heegaard word group and antipode labels — coordinate presentation and geometric evidence boundary.
cf:[[Kup91#^kup91-blk-0140]] group-state constraints on lower circles and upper attaching words — comparison target, not a certification of this HC coordinate map.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology),"..","..","threads","PJ260922","data","261107-lickorish-data","group-15","case-01.hc"));

julia> p = heegaard_pi1_presentation(Topology.HeegaardCodes.read_hc(path)); Tuple(only(p.relators))
((1, -1), (1, 1), (1, -1))
```

Role: convert
"""
function heegaard_pi1_presentation(hc::Heegaard.HeegaardCode)
    owners = _lc_label_owners(hc.lc.words)
    alphas = sort(_hc_role_curves(hc, :alpha))
    betas = sort(_hc_role_curves(hc, :beta))
    alpha_index = Dict(name => i for (i, name) in enumerate(alphas))
    relators = Vector{Vector{Tuple{Int, Int}}}()
    for beta_name in betas
        relator = Tuple{Int, Int}[]
        for x in hc.lc.words[beta_name]
            label = abs(x)
            _is_alpha_beta_label(owners, label, hc.roles) || continue
            sides = owners[label]
            if get(hc.roles, sides[:+], :none) == :beta
                push!(relator, (alpha_index[sides[:-]], 1))
            elseif get(hc.roles, sides[:-], :none) == :beta
                push!(relator, (alpha_index[sides[:+]], -1))
            end
        end
        push!(relators, relator)
    end
    Heegaard.GroupPresentation(
        ["x$name" for name in alphas],
        relators,
        Int[],
        Dict{String, Any}(
            "source" => "lc_beta_words",
            "convention" => "relator walks each beta word; passage sign §3",
        ),
    )
end

function _require_heegaard_cutting_domain(hc::Heegaard.HeegaardCode)
    Heegaard.validate_hc_input(hc)
    g = hc.lc.genus
    alphas = _hc_role_curves(hc, :alpha)
    betas = _hc_role_curves(hc, :beta)
    length(alphas) == g || throw(DomainError(hc,
        "complete Heegaard diagram requires genus-many alpha curves"))
    length(betas) == g || throw(DomainError(hc,
        "complete Heegaard diagram requires genus-many beta curves"))
    foreach(name -> haskey(hc.lc.words, name) || throw(ArgumentError("missing curve word $name")),
        vcat(alphas, betas))
    owners = _lc_label_owners(hc.lc.words)
    foreach(label -> begin
        sides = owners[label]
        (haskey(sides, :+) && haskey(sides, :-)) ||
            throw(ArgumentError("label $label is not a valid passage"))
    end, keys(owners))
    hc
end


function _kuperberg_genus_normalization(pairing, genus::Int)
    isapprox(pairing, 1; atol=_KUPERBERG_INTEGRAL_ATOL, rtol=_KUPERBERG_INTEGRAL_ATOL) &&
        return one(ComplexF64)
    ComplexF64(pairing)^(-genus)
end

_kuperberg_genus_normalization(algebra::FiniteHopfAlgebra, genus::Int) =
    _kuperberg_genus_normalization(_kuperberg_right_integral_pairing(algebra), genus)

_kuperberg_genus_normalization(algebra::FiniteHopfAlgebra, hc::Heegaard.HeegaardCode) =
    _kuperberg_genus_normalization(algebra, hc.lc.genus)


function _port_find!(parent::Dict{Int, Int}, v::Int)
    p = get(parent, v, v)
    if p != v
        parent[v] = _port_find!(parent, p)
    elseif !haskey(parent, v)
        parent[v] = v
    end
    parent[v]
end

function _port_union!(parent::Dict{Int, Int}, a::Int, b::Int)
    ra, rb = _port_find!(parent, a), _port_find!(parent, b)
    ra != rb && (parent[rb] = ra)
end

function _push_coproduct_chain!(factors, algebra::FiniteHopfAlgebra, legs::Vector{Int}, fresh)
    k = length(legs)
    if k == 0
        wire = fresh()
        push!(factors, label_tensor(algebra.counit, [wire]))
        push!(factors, label_tensor(algebra.right_cointegral, [wire]))
        return
    end
    if k == 1
        push!(factors, label_tensor(algebra.right_cointegral, legs))
        return
    end
    wire = fresh()
    push!(factors, label_tensor(algebra.right_cointegral, [wire]))
    for j in 1:(k - 1)
        if j < k - 1
            nxt = fresh()
            push!(factors, label_tensor(algebra.coproduct, [wire, legs[j], nxt]))
            wire = nxt
        else
            push!(factors, label_tensor(algebra.coproduct, [wire, legs[k - 1], legs[k]]))
        end
    end
end

function _push_product_chain!(factors, algebra::FiniteHopfAlgebra, legs::Vector{Int}, fresh)
    k = length(legs)
    if k == 0
        wire = fresh()
        push!(factors, label_tensor(algebra.right_integral, [wire]))
        push!(factors, label_tensor(algebra.unit, [wire]))
        return
    end
    if k == 1
        push!(factors, label_tensor(algebra.right_integral, legs))
        return
    end
    wire = fresh()
    push!(factors, label_tensor(algebra.multiplication, [legs[1], legs[2], wire]))
    for j in 3:k
        nxt = fresh()
        push!(factors, label_tensor(algebra.multiplication, [wire, legs[j], nxt]))
        wire = nxt
    end
    push!(factors, label_tensor(algebra.right_integral, [wire]))
end

function _ab_crossing_visits(gauss::Vector{Vector{Int}}, alpha_count::Int)
    alpha = Dict{Int, Tuple{Int, Int}}()
    beta = Dict{Int, Tuple{Int, Int}}()
    for strand in 1:alpha_count
        for (idx, c) in enumerate(gauss[strand])
            alpha[abs(c)] = (strand, idx)
        end
    end
    for strand in (alpha_count + 1):length(gauss)
        for (idx, c) in enumerate(gauss[strand])
            beta[abs(c)] = (strand, idx)
        end
    end
    alpha, beta
end

"""Coproduct chain of Λ on α, product chain of λ on β, antipode at negative crossings.

Antipode terms are stored as `(output, input)`. The contracted edge runs from
the α leg `a` to the β leg `b`, so the β index is the output and the α index
is the input.
"""
function _kuperberg_factors_from_words(alpha_words::Vector{Vector{Int}},
                                       beta_words::Vector{Vector{Int}},
                                       sign_of::Dict{Int, Int},
                                       algebra::FiniteHopfAlgebra)
    alpha_count = length(alpha_words)
    alpha_legs = Dict{Tuple{Int, Int}, Int}()
    beta_legs = Dict{Tuple{Int, Int}, Int}()
    var_id = 0
    for (strand, word) in enumerate(alpha_words)
        for idx in eachindex(word)
            var_id += 1
            alpha_legs[(strand, idx)] = var_id
        end
    end
    for (strand, word) in enumerate(beta_words)
        for idx in eachindex(word)
            var_id += 1
            beta_legs[(strand, idx)] = var_id
        end
    end
    alpha_at = Dict{Int, Int}()
    beta_at = Dict{Int, Int}()
    for (strand, word) in enumerate(alpha_words), (idx, xid) in enumerate(word)
        alpha_at[xid] = alpha_legs[(strand, idx)]
    end
    for (strand, word) in enumerate(beta_words), (idx, xid) in enumerate(word)
        beta_at[xid] = beta_legs[(strand, idx)]
    end
    parent = Dict{Int, Int}()
    negative = Tuple{Int, Int}[]
    for xid in keys(sign_of)
        haskey(alpha_at, xid) && haskey(beta_at, xid) || continue
        a, b = alpha_at[xid], beta_at[xid]
        if sign_of[xid] > 0
            _port_union!(parent, a, b)
        else
            push!(negative, (a, b))
        end
    end
    resolve(v) = _port_find!(parent, v)
    fresh() = (var_id += 1)
    factors = SparseFactor[]
    for strand in 1:alpha_count
        legs = [resolve(alpha_legs[(strand, j)]) for j in eachindex(alpha_words[strand])]
        _push_coproduct_chain!(factors, algebra, legs, fresh)
    end
    for strand in eachindex(beta_words)
        legs = [resolve(beta_legs[(strand, j)]) for j in eachindex(beta_words[strand])]
        _push_product_chain!(factors, algebra, legs, fresh)
    end
    for (a, b) in negative
        push!(factors, label_tensor(algebra.antipode, [resolve(b), resolve(a)]))
    end
    factors
end

function _heegaard_kuperberg_factors(hc::Heegaard.HeegaardCode, algebra::FiniteHopfAlgebra)
    link = _heegaard_kuperberg_link(hc)
    gauss, orientation = real_gauss_code(link)
    alpha_count = length(_hc_role_curves(hc, :alpha))
    alpha_words = [abs.(word) for word in gauss[1:alpha_count]]
    beta_words = [abs.(word) for word in gauss[(alpha_count + 1):end]]
    sign_of = Dict(index => orientation[index] for index in eachindex(orientation))
    _kuperberg_factors_from_words(alpha_words, beta_words, sign_of, algebra)
end

function _diagram_crossing_words(diagram::Heegaard.HeegaardDiagram)
    code = diagram.code.code
    alphas = sort(code.alpha; by=curve -> curve.id)
    betas = sort(code.beta; by=curve -> curve.id)
    alpha_words = [Int.(curve.visits) for curve in alphas]
    beta_words = [Int.(curve.visits) for curve in betas]
    sign_of = Dict(point.id => point.sign for point in code.intersections)
    alpha_words, beta_words, sign_of
end

function _dense_rank_cap(dimension::Int)
    rank = 0
    acc = 1.0
    while acc * dimension <= 2.0e7
        acc *= dimension
        rank += 1
    end
    rank
end

function _size_within(dimension::Int, rank::Int, limit::Float64)
    acc = 1.0
    for _ in 1:rank
        acc *= dimension
        acc > limit && return false
    end
    true
end

function _candidate_pairs(leg_sets::Vector{Set{Int}})
    pairs = Tuple{Int, Int, Int}[]
    for i in 1:(length(leg_sets) - 1)
        for j in (i + 1):length(leg_sets)
            shared = 0
            for leg in leg_sets[i]
                leg in leg_sets[j] && (shared += 1)
            end
            shared == 0 && continue
            push!(pairs, (length(leg_sets[i]) + length(leg_sets[j]) - 2 * shared, i, j))
        end
    end
    pairs
end

function _record_contraction(leg_lists::Vector{Vector{Int}}, pick)
    sets = [Set(legs) for legs in leg_lists]
    order = Tuple{Int, Int}[]
    peak = maximum(length, sets; init=0)
    while true
        pairs = _candidate_pairs(sets)
        isempty(pairs) && break
        _, i, j = pick(pairs)
        shared = intersect(sets[i], sets[j])
        merged = setdiff(union(sets[i], sets[j]), shared)
        peak = max(peak, length(merged))
        push!(order, (i, j))
        deleteat!(sets, max(i, j))
        deleteat!(sets, min(i, j))
        push!(sets, merged)
    end
    order, peak
end

function _pick_contraction_pair(pairs, rng::_KuperbergRNG)
    best = minimum(pair -> pair[1], pairs)
    tight = [pair for pair in pairs if pair[1] == best]
    if _kuperberg_rand!(rng) % 4 == 0
        wider = [pair for pair in pairs if pair[1] <= best + 1]
        return wider[Int(_kuperberg_rand!(rng) % length(wider)) + 1]
    end
    tight[Int(_kuperberg_rand!(rng) % length(tight)) + 1]
end

function _best_contraction_order(leg_lists::Vector{Vector{Int}}, dimension::Int)
    cap = _dense_rank_cap(dimension)
    best_order, best_peak = _record_contraction(leg_lists, pairs -> first(sort!(copy(pairs))))
    _size_within(dimension, best_peak, 1.0e6) && return best_order, best_peak
    # Trial count is fixed by the network, not by the clock, so a slow machine
    # reaches the same trial as a fast one. Seeds are the trial numbers.
    trials = _size_within(dimension, best_peak, 2.0e7) ? 40 : 800
    for trial in 1:trials
        rng = _KuperbergRNG(UInt64(trial) + 0x9e3779b97f4a7c15)
        order, peak = _record_contraction(leg_lists, pairs -> _pick_contraction_pair(pairs, rng))
        if peak < best_peak
            best_peak = peak
            best_order = order
        end
        best_peak <= cap && _size_within(dimension, best_peak, 1.0e6) && break
    end
    best_peak > cap && throw(ArgumentError(
        "no dense contraction order within rank $cap (best peak $best_peak)"))
    best_order, best_peak
end

function _dense_labeled(factor::SparseFactor, dimension::Int)
    legs = Int.(factor.variables)
    rank = length(legs)
    if rank == 0
        return legs, fill(ComplexF64(get(factor.terms, (), 0.0 + 0im)))
    end
    data = zeros(ComplexF64, ntuple(_ -> dimension, rank))
    for (key, value) in factor.terms
        data[key...] = ComplexF64(value)
    end
    legs, data
end

function _move_axes(data::Array{ComplexF64}, legs::Vector{Int}, order::Vector{Int})
    isempty(order) && return data
    perm = [findfirst(==(leg), legs)::Int for leg in order]
    ndims(data) == 0 ? data : permutedims(data, perm)
end

function _contract_dense_pair(legs_a, data_a, legs_b, data_b, dimension::Int)
    shared = Int[leg for leg in legs_a if leg in legs_b]
    free_a = Int[leg for leg in legs_a if leg ∉ shared]
    free_b = Int[leg for leg in legs_b if leg ∉ shared]
    left = _move_axes(data_a, legs_a, vcat(free_a, shared))
    right = _move_axes(data_b, legs_b, vcat(shared, free_b))
    na = dimension^length(free_a)
    ns = dimension^length(shared)
    nb = dimension^length(free_b)
    product = reshape(left, na, ns) * reshape(right, ns, nb)
    out_legs = vcat(free_a, free_b)
    out = isempty(out_legs) ? fill(product[1]) : reshape(product, ntuple(_ -> dimension, length(out_legs)))
    out_legs, out
end

function _contract_kuperberg(factors::Vector{SparseFactor}, dimension::Int)
    isempty(factors) && return 1.0 + 0im
    tensors = Tuple{Vector{Int}, Array{ComplexF64}}[_dense_labeled(factor, dimension) for factor in factors]
    leg_lists = [copy(legs) for (legs, _) in tensors]
    order, _ = _best_contraction_order(leg_lists, dimension)
    for (i, j) in order
        legs, data = _contract_dense_pair(tensors[i][1], tensors[i][2], tensors[j][1], tensors[j][2], dimension)
        deleteat!(tensors, max(i, j))
        deleteat!(tensors, min(i, j))
        push!(tensors, (legs, data))
    end
    total = 1.0 + 0im
    for (legs, data) in tensors
        isempty(legs) || throw(ArgumentError("Kuperberg network has an open leg"))
        total *= data[]
    end
    total
end

"""    KuperbergTensorExpression(hc::HeegaardCode, link::VirtualLink, factors::Vector{SparseFactor}, normalization::ComplexF64, dimension::Int) -> KuperbergTensorExpression

Store the finite tensor contraction associated with a Heegaard-code calculation.

# Arguments
The fields use the [Kuperberg coordinate target](../specs/hopf-invariants.md#kuperberg-coordinate-target).
H6 ([conditions](../specs/hopf-invariants.md#numbered-conditions)) is the
caller's responsibility: the constructor does not validate the algebra,
diagram, closed legs or normalization. The general field-converting method
converts arguments to the declared types; unsupported conversions propagate
Julia conversion errors.

# Returns
A `KuperbergTensorExpression` with input `hc`, incidence `link`, tensor
`factors`, scalar `normalization` and common `dimension`, with meanings
specified in the linked coordinate-target entry.

# Guarantees
1. Fields equal the supplied values after conversion. Under H6 they specify the normalized finite coordinate contraction of the associated Heegaard calculation, subject to the [open HC/cellular sign convention](../conventions/Algebra.md#kuperberg-sign-contradiction) when interpreted as a Kuperberg invariant. [proved/agent: by construction]
2. Construction does not modify its arguments; already correctly typed mutable fields remain shared with the result. [proved/agent: by construction]

# Conventions
[Kuperberg normalization](../specs/hopf-invariants.md#kuperberg-coordinate-target)
and [tensor ports](../conventions/Algebra.md#tensor-ports).

# Not guaranteed
A raw record does not certify H1–H6 or invariance of a represented
3-manifold. The source crossing/dual/normalization identification is
unverified as recorded in the specification and inherited source references.
The [intended HC/cellular sign convention](../conventions/Algebra.md#kuperberg-sign-contradiction)
is undecided; equality of the two entrances' values is not established.
Numerical evaluation is not an exact algebraic certificate.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md, Kuperberg coordinate target and H6 — expression record meaning and validation boundary.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc"));

julia> expr = invariant_expression(Heegaard.read_hc(path), cyclic_group_algebra(2));

julia> isapprox(evaluate_invariant(expr), 1; atol=1e-8)
true
```

See also [`invariant_expression`](@ref), `evaluate_invariant`.

Role: convert
"""
struct KuperbergTensorExpression
    hc::Heegaard.HeegaardCode
    link::VirtualLink
    factors::Vector{SparseFactor}
    normalization::ComplexF64
    dimension::Int
end

"""
    invariant_expression(hc::HeegaardCodes.HeegaardCode, algebra::FiniteHopfAlgebra) -> KuperbergTensorExpression

Return the normalized finite Kuperberg coordinate network of a Heegaard code.

# Arguments
Requires H1–H5 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
Numerical involutivity and nonzero distinguished pairing are checked with
`DomainError` identifying the failed algebra requirement. HC validation
errors propagate. Wrong complete-family counts raise `DomainError` naming
alpha or beta; missing words, conflicting owners or invalid passages raise
`ArgumentError`. Full geometric cutting-system and source-port hypotheses
remain caller obligations. Tensor shape and arithmetic errors propagate.

# Returns
A `KuperbergTensorExpression` containing the input `hc`, a new alpha/beta
incidence `link`, coordinate `factors`, scalar `normalization` and algebra
`dimension`, as defined by the
[Kuperberg target](../specs/hopf-invariants.md#kuperberg-coordinate-target).

# Guarantees
1. On the supported coordinate domain, the expression denotes that normalized finite contraction with the HC entrance's current sign assignment; its intended relation to cellular signs is the [open sign point](../conventions/Algebra.md#kuperberg-sign-contradiction). [proved/agent: by construction]
2. Inputs are not modified. The returned `hc` field shares the input HC object. [proved/agent: by construction]

# Not guaranteed
The input guards do not certify H1–H3, complete cutting systems, or the
general source invariant identification. The source antipode and
dimension-versus-pairing normalization comparisons remain unverified.
The [intended HC/cellular sign convention](../conventions/Algebra.md#kuperberg-sign-contradiction)
is undecided; equality with the cellular entrance's contraction is not established.
Numeric tensor formation provides no exact-arithmetic or error certificate.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Kuperberg coordinate target and H1–H6 — finite network and geometric boundary.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology),"..","..","threads","PJ260922","data","261107-lickorish-data","group-15","case-01.hc"));

julia> x = invariant_expression(Topology.HeegaardCodes.read_hc(path),cyclic_group_algebra(2));

julia> isapprox(evaluate_invariant(x),1; atol=1e-8)
true
```

Role: convert
"""
function invariant_expression(hc::Heegaard.HeegaardCode, algebra::FiniteHopfAlgebra)
    _require_involutory_kuperberg_algebra(algebra)
    _require_heegaard_cutting_domain(hc)
    link = _heegaard_kuperberg_link(hc)
    KuperbergTensorExpression(
        hc,
        link,
        _heegaard_kuperberg_factors(hc, algebra),
        _kuperberg_genus_normalization(algebra, hc),
        algebra.dimension,
    )
end

"""
    invariant_expression(path::AbstractString, algebra::FiniteHopfAlgebra) -> KuperbergTensorExpression

Read a Heegaard code and return its normalized finite Kuperberg coordinate network.

# Arguments
`path` is a readable file accepted by `HeegaardCodes.read_hc`; I/O, parsing
and HC validation errors propagate. The parsed code and `algebra` require
H1–H5 ([conditions](../specs/hopf-invariants.md#numbered-conditions)), with
the checks and failures of `invariant_expression(::HeegaardCode,algebra)`.
Complete cutting systems and the geometric/source-port hypotheses remain
caller obligations.

# Returns
A `KuperbergTensorExpression` containing the parsed HC, incidence link,
coordinate factors, normalization and dimension, as defined by the
[Kuperberg target](../specs/hopf-invariants.md#kuperberg-coordinate-target).

# Guarantees
1. The result is the coordinate expression of the parsed HC and supplied algebra, with the HC method's supported finite meaning and the same [open HC/cellular sign point](../conventions/Algebra.md#kuperberg-sign-contradiction). [proved/agent: by construction]
2. The input file and algebra are not modified. [proved/agent: by construction]

# Not guaranteed
Parsing and algebra guards do not certify H1–H3 or the full manifold-invariant
identification. The [intended HC/cellular sign convention](../conventions/Algebra.md#kuperberg-sign-contradiction)
is undecided; equality with the cellular entrance's contraction is not established.
Numeric tensor formation has no exact-arithmetic certificate.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Kuperberg coordinate target — finite network and geometric boundary.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology),"..","..","threads","PJ260922","data","261107-lickorish-data","group-15","case-01.hc"));

julia> isapprox(evaluate_invariant(invariant_expression(path,cyclic_group_algebra(2))),1; atol=1e-8)
true
```

Role: convert
"""
invariant_expression(path::AbstractString, algebra::FiniteHopfAlgebra) =
    invariant_expression(Heegaard.read_hc(path), algebra)

"""Copy the factors. Every accepted `strategy` performs this same contraction.

`Z_MST` for a `FiniteHopfAlgebra` accepts the same names on
`simplify_invariant(::MSTTensorExpression)` and does not rewrite the network
either. Taft's symbolic reducer is a different method and is not used here.
`:auto`, `:qa`, `:network`, `:none`, and `:numeric` therefore leave the value
unchanged. An unknown name still throws.
"""
function simplify_invariant(expr::KuperbergTensorExpression; strategy=:auto)
    strategy in (:auto, :qa, :network, :none, :numeric) ||
        throw(ArgumentError("unknown Kuperberg strategy $(strategy)"))
    KuperbergTensorExpression(
        expr.hc,
        expr.link,
        deepcopy(expr.factors),
        expr.normalization,
        expr.dimension,
    )
end

function evaluate_invariant(expr::KuperbergTensorExpression)
    isempty(expr.factors) ? expr.normalization :
        _contract_kuperberg(copy(expr.factors), expr.dimension) * expr.normalization
end

function _kuperberg_value(alpha_words, beta_words, sign_of, algebra, genus::Int)
    _require_involutory_kuperberg_algebra(algebra)
    factors = _kuperberg_factors_from_words(alpha_words, beta_words, sign_of, algebra)
    _contract_kuperberg(factors, algebra.dimension) *
        _kuperberg_genus_normalization(algebra, genus)
end

function Z_Kuperberg(hc::Heegaard.HeegaardCode, algebra::FiniteHopfAlgebra;
                     strategy=:auto, kwargs...)
    isempty(kwargs) || throw(ArgumentError("Z_Kuperberg does not accept extra keywords"))
    evaluate_invariant(simplify_invariant(invariant_expression(hc, algebra); strategy=strategy))
end

Z_Kuperberg(path::AbstractString, algebra::FiniteHopfAlgebra; kwargs...) =
    Z_Kuperberg(Heegaard.read_hc(path), algebra; kwargs...)

function Z_Kuperberg(diagram::Heegaard.HeegaardDiagram, algebra::FiniteHopfAlgebra; kwargs...)
    isempty(kwargs) || throw(ArgumentError("Z_Kuperberg does not accept extra keywords"))
    alpha_words, beta_words, sign_of = _diagram_crossing_words(diagram)
    _kuperberg_value(alpha_words, beta_words, sign_of, algebra, Heegaard.surface_genus(diagram))
end

Z_Kuperberg(punctured::Heegaard.PuncturedHeegaardDiagram, algebra::FiniteHopfAlgebra; kwargs...) =
    Z_Kuperberg(punctured.diagram, algebra; kwargs...)

"""
    kuperberg_hc_negative_labels(hc::HeegaardCodes.HeegaardCode) -> Set{Int}

Return the alpha/beta passage labels assigned an antipode by the Kuperberg coordinate convention.

# Arguments
Requires H8 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
Conflicting same-sign owners raise `ArgumentError` naming the label.
Other malformed data can raise indexing or conversion errors; HC validity
and cutting-system completeness are not checked. One-sided and non-alpha/beta
passages are omitted rather than certified as valid crossings.

# Returns
A set of integer labels with the meaning in
[Heegaard antipode labels](../specs/hopf-invariants.md#heegaard-word-group-and-antipode-labels).

# Guarantees
1. The set contains exactly the alpha/beta labels selected by that coordinate convention. [proved/agent: by construction]
2. On the three convertible fixtures listed in the [HC/cellular sign contradiction](../conventions/Algebra.md#kuperberg-sign-contradiction), the set converted to strings equals the `source_label` set of positive-sign intersections of `forget_reference(realize_hc(hc)).diagram`, and differs from its negative-sign set. [verified-computationally: three fixtures of test/kuperberg_invariant.jl, against the assertion of that test]
3. The input is unchanged. [proved/agent: by construction]

# Conventions
[Kuperberg and cellular crossing signs](../conventions/Algebra.md#kuperberg-and-cellular-crossing-signs).

# Not guaranteed
The intended HC/cellular sign convention is undecided: the
[recorded bug](../conventions/Algebra.md#kuperberg-sign-contradiction) is a
contradiction between the current code and its integration test.
The returned set is neither an invariant nor an HC validity certificate. Source agreement remains
unverified: identify the HC alpha/beta families with Kup91's lower/upper
circles and match their ordered tangent bases and antipode input/output at
both local crossing orientations in equations 26–28. The coordinate sign
translation alone does not establish that source identification or equality
of the two Kuperberg evaluator entrances.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Heegaard word group and antipode labels — coordinate selection and unresolved sign comparison.

# Examples
```jldoctest
julia> path = normpath(joinpath(pkgdir(Topology),"..","..","threads","PJ260922","data","261107-lickorish-data","group-15","case-01.hc"));

julia> hc = Topology.HeegaardCodes.read_hc(path); diagram = Topology.HeegaardCodes.forget_reference(Topology.HeegaardCodes.realize_hc(hc)).diagram;

julia> labels = Set(string(k) for k in kuperberg_hc_negative_labels(hc));

julia> labels == Set(p.source_label for p in diagram.code.code.intersections if p.sign > 0 && p.source_label !== nothing)
true

julia> (Tuple(sort!(collect(labels))), Tuple(sort!([p.source_label for p in diagram.code.code.intersections if p.sign < 0 && p.source_label !== nothing])))
(("5",), ("4", "6"))
```

Role: access
"""
function kuperberg_hc_negative_labels(hc::Heegaard.HeegaardCode)
    owners = _lc_label_owners(hc.lc.words)
    Set(label for label in keys(owners)
        if _is_alpha_beta_label(owners, label, hc.roles) &&
           !_ab_crossing_passage_positive(hc, label))
end

export Z_Kuperberg, KuperbergTensorExpression, invariant_expression,
       simplify_invariant, evaluate_invariant, heegaard_pi1_presentation,
       kuperberg_hc_negative_labels
