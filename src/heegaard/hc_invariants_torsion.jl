# cf: Massuyeau, An introduction to the abelian Reidemeister torsion of three-dimensional manifolds, arXiv:1003.2517, Lemma 5.3 — abelian Reidemeister torsion from the Heegaard Fox matrix. τ_0 is taken to be +1. β^♯ is the surface class meeting that β curve once and the other β curves zero times, pushed to H₁. reconstructed. No registered Work key.
import AkaghefAlgebra

"""Abelian Reidemeister torsion of a closed HC with b₁ = 0, evaluated at characters of H₁.

`order` is the exponent of H₁ and the cyclotomic order. `invariant_factors` is the
Smith diagonal, including 1s, so a character key has one entry per factor.
`keys` lists the nontrivial characters. `values[i]` is τ at `keys[i]`, or `nothing`
when the twisted complex is not acyclic. `acyclic` is true when at least one
nontrivial character is acyclic.
"""
struct HCReidemeisterTorsion{N}
    order::Int
    invariant_factors::Vector{Int}
    keys::Vector{Vector{Int}}
    values::Vector{Union{Nothing,AkaghefAlgebra.TaftExact.Cyclotomic{N}}}
    acyclic::Bool
end

struct HCTorsionError <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::HCTorsionError) = print(io, "HCTorsionError(", e.reason, "): ", e.detail)
_hct_require(ok, reason, detail) = ok || throw(HCTorsionError(reason, detail))

function _hct_smith(matrix::AbstractMatrix{<:Integer})
    A = Int.(matrix)
    rows, cols = size(A)
    U = Matrix{Int}(I, rows, rows)
    V = Matrix{Int}(I, cols, cols)
    k = 1
    while k <= min(rows, cols)
        location = nothing
        for i in k:rows, j in k:cols
            if A[i, j] != 0
                location = (i, j)
                break
            end
        end
        location === nothing && break
        if location[1] != k
            A[[k, location[1]], :] = A[[location[1], k], :]
            U[[k, location[1]], :] = U[[location[1], k], :]
        end
        if location[2] != k
            A[:, [k, location[2]]] = A[:, [location[2], k]]
            V[:, [k, location[2]]] = V[:, [location[2], k]]
        end
        while true
            reduced = false
            for i in k+1:rows
                A[i, k] == 0 && continue
                q = div(A[i, k], A[k, k])
                A[i, :] .-= q .* A[k, :]
                U[i, :] .-= q .* U[k, :]
                if A[i, k] != 0 && abs(A[i, k]) < abs(A[k, k])
                    A[[k, i], :] = A[[i, k], :]
                    U[[k, i], :] = U[[i, k], :]
                end
                reduced = true
                break
            end
            reduced && continue
            for j in k+1:cols
                A[k, j] == 0 && continue
                q = div(A[k, j], A[k, k])
                A[:, j] .-= q .* A[:, k]
                V[:, j] .-= q .* V[:, k]
                if A[k, j] != 0 && abs(A[k, j]) < abs(A[k, k])
                    A[:, [k, j]] = A[:, [j, k]]
                    V[:, [k, j]] = V[:, [j, k]]
                end
                reduced = true
                break
            end
            reduced && continue
            bad = nothing
            for i in k+1:rows, j in k+1:cols
                if mod(A[i, j], A[k, k]) != 0
                    bad = (i, j)
                    break
                end
            end
            bad === nothing && break
            A[k, :] .+= A[bad[1], :]
            U[k, :] .+= U[bad[1], :]
        end
        if A[k, k] < 0
            A[k, :] .*= -1
            U[k, :] .*= -1
        end
        k += 1
    end
    A, U, V
end

function _hct_solve(W::AbstractMatrix{<:Integer}, rhs::AbstractVector{<:Integer})
    D, U, V = _hct_smith(W)
    image = U * Int.(rhs)
    z = zeros(Int, size(W, 2))
    for i in 1:size(W, 1)
        pivot = i <= size(D, 2) ? D[i, i] : 0
        if pivot == 0
            image[i] == 0 || throw(HCTorsionError(:dual_class, "the dual system has no integer solution"))
        else
            image[i] % pivot == 0 || throw(HCTorsionError(
                :dual_class, "the dual class is not integral ($(image[i]) / $pivot)"))
            z[i] = div(image[i], pivot)
        end
    end
    v = V * z
    W * v == Int.(rhs) || throw(HCTorsionError(:dual_class, "the dual class does not meet the beta curves"))
    v
end

_hct_dot(left, right) = sum(left[i] * right[i] for i in eachindex(left); init=0)

function _hct_omega(owners, c, d)
    s = 0
    for (_, (plus_c, minus_c)) in owners
        if minus_c == c && plus_c == d
            s += 1
        elseif minus_c == d && plus_c == c
            s -= 1
        end
    end
    s
end

function _hct_coords(hc, owners, name)
    genus = hc.lc.genus
    for i in 1:genus
        haskey(hc.lc.words, "l$i") && haskey(hc.lc.words, "m$i") || throw(HCTorsionError(
            :missing_reference, "l$i or m$i is missing"))
    end
    A = [_hct_omega(owners, name, "m$i") for i in 1:genus]
    B = [-_hct_omega(owners, name, "l$i") for i in 1:genus]
    vcat(A, B)
end

function _hct_dual_exponents(hc, presentation)
    genus = length(presentation.generators)
    genus == 0 && return Vector{Int}[]
    owners = _hc_owners(hc.lc.words)
    alphas = presentation.metadata["generator_names"]
    betas = presentation.metadata["relator_names"]
    alpha = hcat([_hct_coords(hc, owners, name) for name in alphas]...)
    beta = hcat([_hct_coords(hc, owners, name) for name in betas]...)
    # ω((A,B),(A',B')) = A·B' − B·A', so J*(A,B) = (B, −A).
    JB = vcat(beta[genus+1:end, :], .-beta[1:genus, :])
    exponents = Vector{Vector{Int}}()
    for j in 1:genus
        rhs = zeros(Int, genus)
        rhs[j] = 1
        class = _hct_solve(transpose(JB), rhs)
        push!(exponents, [_hct_dot(class, vcat(alpha[genus+1:end, i], .-alpha[1:genus, i])) for i in 1:genus])
    end
    exponents
end

function _hct_fox_entry(relator, generator, genus)
    poly = Dict{Vector{Int},Int}()
    prefix = zeros(Int, genus)
    for (index, exponent) in relator
        if index == generator
            if exponent == 1
                key = copy(prefix)
                poly[key] = get(poly, key, 0) + 1
                prefix[index] += 1
            elseif exponent == -1
                prefix[index] -= 1
                key = copy(prefix)
                poly[key] = get(poly, key, 0) - 1
            else
                throw(HCTorsionError(:bad_letter, "a Fox letter is not a unit exponent"))
            end
        else
            prefix[index] += exponent
        end
    end
    poly
end

function _hct_eval_poly(poly, phi)
    sample = first(phi)
    acc = zero(sample)
    for (exponents, coefficient) in poly
        monomial = one(sample)
        for index in eachindex(phi)
            exponents[index] == 0 && continue
            monomial *= phi[index]^exponents[index]
        end
        acc += coefficient * monomial
    end
    acc
end

function _hct_det(matrix)
    n = size(matrix, 1)
    n == 0 && throw(ArgumentError("the empty minor is handled by the caller"))
    work = copy(matrix)
    determinant = one(work[1, 1])
    for col in 1:n
        pivot = findfirst(i -> !iszero(work[i, col]), col:n)
        pivot === nothing && return zero(determinant)
        pivot += col - 1
        if pivot != col
            work[col, :], work[pivot, :] = work[pivot, :], work[col, :]
            determinant = -determinant
        end
        pivot_value = work[col, col]
        determinant *= pivot_value
        for row in col+1:n
            iszero(work[row, col]) && continue
            work[row, :] .-= (work[row, col] / pivot_value) .* work[col, :]
        end
    end
    determinant
end

function _hct_partial3(exponents, phi)
    one_value = one(first(phi))
    map(row -> begin
        value = one_value
        for index in eachindex(phi)
            row[index] == 0 && continue
            value *= phi[index]^row[index]
        end
        value - one_value
    end, exponents)
end

function _hct_sign_patterns(genus::Int)
    patterns = Vector{Vector{Int}}()
    for mask in 0:(1 << genus) - 1
        push!(patterns, [((mask >> (j - 1)) & 1) == 1 ? -1 : 1 for j in 1:genus])
    end
    # Fewer minus signs first, and a minus on a later curve before a minus on an earlier one.
    sort!(patterns, by=signs -> (count(==(-1), signs), [sign == -1 ? 1 : 0 for sign in signs]))
    patterns
end

function _hct_partial3_signed(exponents, phi, signs)
    raw = _hct_partial3(exponents, phi)
    [signs[index] * raw[index] for index in eachindex(raw)]
end

function _hct_accept_boundary(fox, right, factors, zeta, exponents, signs, fail_key)
    status, key = _hct_boundary_status(fox, right, factors, zeta, exponents, signs)
    status == :image1 && throw(HCTorsionError(:bad_differential, "the Fox boundary does not square to zero at $key"))
    status == :ok && return (true, exponents, signs, fail_key)
    (false, exponents, signs, isempty(fail_key) ? key : fail_key)
end

function _hct_boundary_status(fox, right, factors, zeta, exponents, signs)
    genus = length(signs)
    for ks in _hct_characters(factors)
        all(iszero, ks) && continue
        phi = _hct_phi(right, factors, ks, zeta)
        matrix = [_hct_eval_poly(fox[row, column], phi) for row in 1:genus, column in 1:genus]
        partial1 = [phi[index] - 1 for index in 1:genus]
        zero_value = zero(zeta)
        image1 = [sum(matrix[row, column] * partial1[column] for column in 1:genus; init=zero_value) for row in 1:genus]
        all(iszero, image1) || return (:image1, ks)
        image3 = _hct_image3(matrix, _hct_partial3_signed(exponents, phi, signs))
        all(iszero, image3) || return (:image3, ks)
    end
    (:ok, Int[])
end

function _hct_image3(matrix, partial3)
    rows, cols = size(matrix)
    zero_value = zero(first(partial3))
    [sum(matrix[row, column] * partial3[row] for row in 1:rows; init=zero_value) for column in 1:cols]
end

function _hct_rank(matrix)
    rows, cols = size(matrix)
    work = copy(matrix)
    rank = 0
    row = 1
    for col in 1:cols
        row > rows && break
        pivot = findfirst(i -> !iszero(work[i, col]), row:rows)
        pivot === nothing && continue
        pivot += row - 1
        work[row, :], work[pivot, :] = work[pivot, :], work[row, :]
        scale = work[row, col]
        for i in row+1:rows
            iszero(work[i, col]) && continue
            work[i, :] .-= (work[i, col] / scale) .* work[row, :]
        end
        row += 1
        rank += 1
    end
    rank
end

function _hct_characters(factors)
    keys = Vector{Int}[[]]
    for factor in factors
        grown = Vector{Int}[]
        for prefix in keys, k in 0:(factor - 1)
            push!(grown, vcat(prefix, k))
        end
        keys = grown
    end
    keys
end

function _hct_phi(V, factors, ks, zeta)
    λ = 1
    for factor in factors
        λ = lcm(λ, factor)
    end
    group = [zeta^((λ ÷ factor) * ks[index]) for (index, factor) in enumerate(factors)]
    generators = size(V, 1)
    phi = typeof(zeta)[]
    for j in 1:generators
        value = one(zeta)
        for i in eachindex(group)
            V[j, i] == 0 && continue
            value *= group[i]^V[j, i]
        end
        push!(phi, value)
    end
    phi
end

function _hct_cyclic_slot(factors)
    isempty(factors) && return nothing
    order = prod(factors)
    order == 1 && return nothing
    slots = findall(==(order), factors)
    length(slots) == 1 || return nothing
    all(factors[i] == 1 || i == only(slots) for i in eachindex(factors)) || return nothing
    only(slots)
end

function _hct_character_values(fox, right, factors, zeta, exponents, signs)
    genus = length(signs)
    values = Any[]
    for ks in _hct_characters(factors)
        all(iszero, ks) && continue
        phi = _hct_phi(right, factors, ks, zeta)
        matrix = [_hct_eval_poly(fox[row, column], phi) for row in 1:genus, column in 1:genus]
        partial1 = [phi[index] - 1 for index in 1:genus]
        partial3 = _hct_partial3_signed(exponents, phi, signs)
        if _hct_rank(matrix) != genus - 1 || all(iszero, partial1) || all(iszero, partial3)
            push!(values, nothing)
            continue
        end
        got = nothing
        for column in 1:genus, row in 1:genus
            (iszero(partial1[column]) || iszero(partial3[row])) && continue
            rows = [index for index in 1:genus if index != row]
            cols = [index for index in 1:genus if index != column]
            minor = isempty(rows) ? one(zeta) : _hct_det(matrix[rows, cols])
            sign = isodd(column + row + genus + 1) ? -one(zeta) : one(zeta)
            torsion = sign * minor / (partial1[column] * partial3[row])
            if got === nothing
                got = torsion
            elseif got != torsion
                return :mismatch
            end
        end
        got === nothing && return :mismatch
        push!(values, got)
    end
    values
end

function _hct_character_unit(order, factors, ks, exponents, zeta)
    value = one(zeta)
    for index in eachindex(factors)
        exponents[index] == 0 && continue
        value *= zeta^((order ÷ factors[index]) * ks[index] * exponents[index])
    end
    value
end

function _hct_same_torsion_class(order, factors, keys, left, right)
    length(left) == length(right) || return false
    for index in eachindex(left)
        (left[index] === nothing) == (right[index] === nothing) || return false
    end
    order <= 1 && return true
    all(value -> value === nothing, left) && return true
    zeta = AkaghefAlgebra.TaftExact.zeta(order)
    for sgn in (1, -1), exponents in _hct_characters(factors)
        good = true
        for index in eachindex(left)
            left[index] === nothing && continue
            unit = _hct_character_unit(order, factors, keys[index], exponents, zeta)
            if left[index] != sgn * unit * right[index]
                good = false
                break
            end
        end
        good && return true
    end
    false
end

function _hct_require_unique_class(order, factors, keys, samples)
    length(samples) <= 1 && return nothing
    for sample in samples[2:end]
        _hct_same_torsion_class(order, factors, keys, samples[1], sample) && continue
        throw(HCTorsionError(
            :ambiguous_orientation,
            "more than one 2-cell orientation squares the Fox boundary to zero and the torsions are not the same class in F^x / ±φ(H1)",
        ))
    end
    nothing
end

function _hct_orientation_tag(exponents, dual, signs)
    (exponents == dual ? "+1" : "inverse") * ":" * string(signs)
end

function _hct_passing_orientations(fox, right, factors, zeta, dual, inverse)
    passing = Tuple{Vector{Vector{Int}}, Vector{Int}}[]
    fail_key = Int[]
    for signs in _hct_sign_patterns(length(dual))
        for exponents in (dual, inverse)
            status, key = _hct_boundary_status(fox, right, factors, zeta, exponents, signs)
            if status == :image1
                throw(HCTorsionError(:bad_differential, "the Fox boundary does not square to zero at $key"))
            elseif status == :ok
                push!(passing, (exponents, copy(signs)))
            elseif isempty(fail_key)
                fail_key = key
            end
        end
    end
    (passing, fail_key)
end

function _hct_select_orientation(fox, right, factors, order, zeta, dual)
    inverse = [.-exponents for exponents in dual]
    passing, fail_key = _hct_passing_orientations(fox, right, factors, zeta, dual, inverse)
    isempty(passing) && throw(HCTorsionError(
        :bad_differential, "the Fox boundary does not square to zero at $fail_key"))
    keys = [ks for ks in _hct_characters(factors) if !all(iszero, ks)]
    chosen, chosen_signs = passing[1]
    chosen_values = _hct_character_values(fox, right, factors, zeta, chosen, chosen_signs)
    if chosen_values !== :mismatch
        tags = String[_hct_orientation_tag(chosen, dual, chosen_signs)]
        samples = Any[chosen_values]
        for index in 2:length(passing)
            exponents, signs = passing[index]
            values = _hct_character_values(fox, right, factors, zeta, exponents, signs)
            values === :mismatch && continue
            _hct_same_torsion_class(order, factors, keys, chosen_values, values) && continue
            push!(samples, values)
            push!(tags, _hct_orientation_tag(exponents, dual, signs))
        end
        if length(samples) > 1
            throw(HCTorsionError(
                :ambiguous_orientation,
                "more than one 2-cell orientation squares the Fox boundary to zero and the torsions are not the same class in F^x / ±φ(H1): " * join(tags, ", "),
            ))
        end
    end
    (chosen, chosen_signs)
end

"""
Reidemeister torsion of a closed HC at the abelian characters of H₁, for b₁ = 0.

Generators are the α curves and relators are the β curves, as in `hc_presentation`.
The Fox matrix is (∂β_k/∂α_j). For a character φ with acyclic twisted homology,

    τ^φ · (φ(α_i) − 1) · (ε_j (φ(β_j^♯) − 1)) = (−1)^{i+j+g+1} det φ(minor_{j,i})

in F^× / ±φ(H₁). Each ε_j is ±1. The representative is the first pattern, in the
order all signs +1 on the +1 class, the same signs on the inverse class, then the
remaining patterns, for which ∂₂∂₃ = 0. Another pattern is kept in that choice
only when its torsion lies in the same class. Two classes raise
`ambiguous_orientation`. τ_0 of Lemma 5.3 stays +1. A character whose twisted
complex is not acyclic is stored as `nothing`. b₁ > 0 raises
`HCTorsionError(:positive_betti)`.
"""
function hc_reidemeister_torsion(hc::HeegaardCode)
    alphas, betas, _ = _hc_curve_lists(hc)
    genus = hc.lc.genus
    _hct_require(length(alphas) == genus && length(betas) == genus, :not_closed,
        "alpha count $(length(alphas)), beta count $(length(betas)), genus $genus")
    presented = hc_presentation(hc).presentation
    relators = presented.relators
    relation = zeros(Int, genus, genus)
    for (row, relator) in enumerate(relators), (column, exponent) in relator
        relation[row, column] += exponent
    end
    diagonal, _, right = _hct_smith(relation)
    factors = genus == 0 ? Int[] : [diagonal[i, i] for i in 1:genus]
    _hct_require(all(!iszero, factors), :positive_betti, "H1 has positive free rank")
    order = 1
    for factor in factors
        order = lcm(order, factor)
    end
    fox = [_hct_fox_entry(relators[row], column, genus) for row in 1:genus, column in 1:genus]
    dual = _hct_dual_exponents(hc, presented)
    # Passage sign: minus on b and plus on a means the frame (b, a) is positive on S.
    # S = ∂Hα with the outward normal of Hα first, and that normal points into Hβ, so
    # M is oriented by (n_out, positive frame of S). ω(c, d) = +1 for that frame.
    # Walking the stored β word, a crossing with the plus passage on α is the Fox
    # letter +α, and that crossing has ω(β, α) = +1. The +1 dual class solves
    # ω(β♯, β_j) = +1. The 2-cell is the co-core disk in Hβ, with Stokes boundary
    # equal to the stored β: its outward-in-disk normal is −n_out. Massuyeau
    # orients e² so that e² · β♯ = +1 and writes ∂₃ = Σ ([β♯] − 1) ê². Genus-one
    # diagrams reproduce the minor formula with τ_0 = +1 and every ε_j = +1, so
    # this reading inserts no further sign. The standard-form dual insertion writes
    # the l word with +ε and the β word with −ε, which flips every l–β intersection
    # and replaces the +1 class by its inverse; it does not single out one 2-cell.
    # Where the Fox augmentation equals the relation matrix, a per-curve minus is
    # not a crossing sign. It is the coefficient in ∂₃ = Σ ε_j ([β_j♯] − 1) ê²_j.
    # When β is Lagrangian, ω(x, β_k) = δ_jk has one solution in H₁(M), and reversing
    # the surface sends it to the class with every exponent negated. Those are the
    # only two solutions of ω(x, β_k) = ±δ_jk. A character of order 2 does not
    # separate them. Alexander still reads the +1 class from `_hct_dual_exponents`.
    # The representative is the first pattern that squares the boundary to zero.
    # Another passing pattern is accepted only when its torsion lies in the same
    # class of F^×/±φ(H₁); two classes raise ambiguous_orientation.
    # The reverse conversion writes vertex sign −S, which cancels the global
    # flip coming from the l = +ε, β = −ε insertion. The five prism diagrams
    # and four lens diagrams rebuilt under 261114-invariants/rebuilt/ then
    # square ∂₂∂₃ on the +1 class with every ε_j = +1, and the inverse class
    # does not. The scan is the fallback for a diagram on which that
    # normalization does not square.
    Z = AkaghefAlgebra.TaftExact
    zeta = order == 1 ? nothing : Z.zeta(order)
    chosen, chosen_signs = _hct_select_orientation(fox, right, factors, order, zeta, dual)
    keys = Vector{Int}[]
    values = Union{Nothing,Z.Cyclotomic{order}}[]
    for ks in _hct_characters(factors)
        all(iszero, ks) && continue
        push!(keys, ks)
        phi = _hct_phi(right, factors, ks, zeta)
        matrix = [_hct_eval_poly(fox[row, column], phi) for row in 1:genus, column in 1:genus]
        partial1 = [phi[index] - 1 for index in 1:genus]
        partial3 = _hct_partial3_signed(chosen, phi, chosen_signs)
        zero_value = zero(zeta)
        image1 = [sum(matrix[row, column] * partial1[column] for column in 1:genus; init=zero_value) for row in 1:genus]
        image3 = _hct_image3(matrix, partial3)
        all(iszero, image1) && all(iszero, image3) || throw(HCTorsionError(
            :bad_differential, "the Fox boundary does not square to zero at $ks"))
        if _hct_rank(matrix) != genus - 1 || all(iszero, partial1) || all(iszero, partial3)
            push!(values, nothing)
            continue
        end
        got = nothing
        for column in 1:genus, row in 1:genus
            (iszero(partial1[column]) || iszero(partial3[row])) && continue
            rows = [index for index in 1:genus if index != row]
            cols = [index for index in 1:genus if index != column]
            minor = isempty(rows) ? one(zeta) : _hct_det(matrix[rows, cols])
            sign = isodd(column + row + genus + 1) ? -one(zeta) : one(zeta)
            torsion = sign * minor / (partial1[column] * partial3[row])
            if got === nothing
                got = torsion
            elseif got != torsion
                throw(HCTorsionError(:minor_mismatch, "Fox minors disagree at $ks"))
            end
        end
        got === nothing && throw(HCTorsionError(:minor_mismatch, "no Fox minor is defined at $ks"))
        push!(values, got)
    end
    HCReidemeisterTorsion{order}(order, factors, keys, values, any(!isnothing, values))
end

function hc_torsion_character(torsion::HCReidemeisterTorsion, exponents::AbstractVector{<:Integer})
    index = findfirst(==(Int.(collect(exponents))), torsion.keys)
    index === nothing && throw(HCTorsionError(
        :character, "no character with exponents $(collect(exponents))"))
    value = torsion.values[index]
    value === nothing && throw(HCTorsionError(:not_acyclic, "the twisted complex is not acyclic"))
    value
end

function hc_torsion_character(torsion::HCReidemeisterTorsion, k::Integer)
    slot = _hct_cyclic_slot(torsion.invariant_factors)
    slot === nothing && throw(HCTorsionError(:character, "H1 is not cyclic of order $(torsion.order)"))
    0 <= k < torsion.order || throw(HCTorsionError(:character, "k=$k is outside 0:$(torsion.order - 1)"))
    key = [index == slot ? Int(k) : 0 for index in eachindex(torsion.invariant_factors)]
    hc_torsion_character(torsion, key)
end

export HCReidemeisterTorsion, HCTorsionError, hc_reidemeister_torsion, hc_torsion_character
