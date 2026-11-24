# Massuyeau, arXiv:1003.2517, Theorem 4.4, PDF p. 19: for a closed 3-manifold,
# β₁ ≥ 2 gives ±τ^μ = Δ in Z[G]/±G, and β₁ = 1 gives ±τ^μ = Δ/(t−1)^2.
# Lemma 4.3, PDF pp. 17–18, with τ₀ = +1, is the Heegaard Fox quotient below.
# b₁ = 0 stays on hc_reidemeister_torsion and is not handled here.
import AbstractAlgebra

"""Reidemeister–Turaev torsion for a closed HC with b₁ ≥ 1.

`free_basis` has one column per generator of the free part of H₁. Column k is the
exponent vector of that generator on the α curves, in the order of the presentation.
`numerator / denominator` lies in the fraction field of Z[t₁^{±1}, …, t_b^{±1}].
The class is taken in that field modulo ± monomials. b₁ = 0 raises
`HCTorsionError(:finite_homology)`.
"""
struct HCGroupRingFraction
    betti::Int
    free_basis::Matrix{Int}
    numerator
    denominator
end

"""Multivariable Alexander polynomial of a closed HC with b₁ ≥ 2.

The polynomial is the torsion of Theorem 4.4, as a Laurent polynomial on the same
`free_basis` as `hc_positive_torsion`. Comparison with another normalization is
modulo GL(b₁, Z) and multiplication by ± a monomial. b₁ < 2 raises `HCTorsionError`.
"""
struct HCMultivariableAlexander
    betti::Int
    free_basis::Matrix{Int}
    polynomial
end

function Base.show(io::IO, value::HCGroupRingFraction)
    print(io, "HCGroupRingFraction(b1=", value.betti, ", ", value.numerator, " / ", value.denominator, ")")
end

function Base.show(io::IO, value::HCMultivariableAlexander)
    print(io, "HCMultivariableAlexander(b1=", value.betti, ", ", value.polynomial, ")")
end

function _hcm_ring(betti::Int)
    names = [Symbol("t", i) for i in 1:betti]
    AbstractAlgebra.laurent_polynomial_ring(AbstractAlgebra.ZZ, names)
end

function _hcm_monomial(ring, variables, exponents)
    monomial = one(ring)
    for index in eachindex(variables)
        exponents[index] == 0 && continue
        monomial *= variables[index]^exponents[index]
    end
    monomial
end

function _hcm_terms(value)
    # Copy each exponent vector as it is produced. The Laurent iterator reuses its buffer.
    pairs = Tuple{Vector{Int},Int}[]
    exponents = AbstractAlgebra.exponent_vectors(value)
    coefficients = AbstractAlgebra.coefficients(value)
    exp_state = iterate(exponents)
    coeff_state = iterate(coefficients)
    while exp_state !== nothing || coeff_state !== nothing
        (exp_state === nothing || coeff_state === nothing) && throw(HCTorsionError(
            :bad_differential, "Laurent term iterators disagree"))
        push!(pairs, (Int.(exp_state[1]), Int(coeff_state[1])))
        exp_state = iterate(exponents, exp_state[2])
        coeff_state = iterate(coefficients, coeff_state[2])
    end
    pairs
end

function _hcm_is_unit(value)
    iszero(value) && return false
    terms = _hcm_terms(value)
    length(terms) == 1 && (terms[1][2] == 1 || terms[1][2] == -1)
end

function _hcm_unit_multiple(left, right)
    iszero(left) && return iszero(right)
    iszero(right) && return false
    ok, quotient = AbstractAlgebra.divides(left, right)
    ok && _hcm_is_unit(quotient)
end

function _hcm_eval(ring, variables, poly, basis)
    acc = zero(ring)
    betti = length(variables)
    for (exponents, coefficient) in poly
        power = zeros(Int, betti)
        for index in eachindex(exponents)
            exponents[index] == 0 && continue
            for slot in 1:betti
                power[slot] += exponents[index] * basis[index, slot]
            end
        end
        acc += coefficient * _hcm_monomial(ring, variables, power)
    end
    acc
end

function _hcm_det_laplace(matrix)
    ring = parent(matrix[1, 1])
    n = size(matrix, 1)
    n == 1 && return matrix[1, 1]
    total = zero(ring)
    for column in 1:n
        iszero(matrix[1, column]) && continue
        minor = AbstractAlgebra.zero_matrix(ring, n - 1, n - 1)
        for row in 2:n
            dest = 1
            for source in 1:n
                source == column && continue
                minor[row - 1, dest] = matrix[row, source]
                dest += 1
            end
        end
        sign = isodd(column) ? 1 : -1
        total += sign * matrix[1, column] * _hcm_det_laplace(minor)
    end
    total
end

# Same value as the Laplace expansion. Fraction-free elimination keeps the
# intermediate polynomials smaller than the expanded product of a large minor.
function _hcm_det(matrix)
    n = size(matrix, 1)
    n == 1 && return matrix[1, 1]
    ring = parent(matrix[1, 1])
    entries = [matrix[row, column] for row in 1:n, column in 1:n]
    sign = one(ring)
    previous = one(ring)
    for pivot_index in 1:(n - 1)
        pivot_row = 0
        for row in pivot_index:n
            iszero(entries[row, pivot_index]) && continue
            pivot_row = row
            break
        end
        pivot_row == 0 && return zero(ring)
        if pivot_row != pivot_index
            for column in 1:n
                entries[pivot_index, column], entries[pivot_row, column] =
                    entries[pivot_row, column], entries[pivot_index, column]
            end
            sign = -sign
        end
        pivot = entries[pivot_index, pivot_index]
        for row in (pivot_index + 1):n
            for column in (pivot_index + 1):n
                numerator = entries[row, column] * pivot -
                    entries[row, pivot_index] * entries[pivot_index, column]
                entries[row, column] = AbstractAlgebra.divexact(numerator, previous)
            end
        end
        previous = pivot
    end
    sign * entries[n, n]
end

function _hcm_reduce(numerator, denominator)
    iszero(denominator) && throw(HCTorsionError(:not_laurent, "the Fox denominator is zero"))
    iszero(numerator) && return zero(numerator), one(parent(numerator))
    common = AbstractAlgebra.gcd(numerator, denominator)
    reduced_num = AbstractAlgebra.divexact(numerator, common)
    reduced_den = AbstractAlgebra.divexact(denominator, common)
    leading = _hcm_terms(reduced_den)[1][2]
    if leading < 0
        reduced_num = -reduced_num
        reduced_den = -reduced_den
    end
    reduced_num, reduced_den
end

function _hcm_free_basis(right, free, genus)
    basis = zeros(Int, genus, length(free))
    for (slot, column) in enumerate(free)
        weights = [right[row, column] for row in 1:genus]
        sign_at = findfirst(!iszero, weights)
        _hct_require(sign_at !== nothing, :bad_differential, "a free generator is zero on every alpha")
        if weights[sign_at] < 0
            weights .*= -1
        end
        basis[:, slot] = weights
    end
    basis
end

"""
    hc_positive_torsion(hc::HeegaardCode) -> HCGroupRingFraction

The representative is the first defined Fox minor

    τ = (−1)^{i+j+g+1} det(minor_{j,i}) / ((φ(α_i) − 1)(φ(β_j^♯) − 1)),

with τ₀ = +1. Later minors are kept only when they differ from it by ± a monomial.
"""
function hc_positive_torsion(hc::HeegaardCode)
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
    free = findall(iszero, factors)
    _hct_require(!isempty(free), :finite_homology, "H1 is finite")
    basis = _hcm_free_basis(right, free, genus)
    betti = length(free)
    ring, variables = _hcm_ring(betti)
    fox = [_hct_fox_entry(relators[row], column, genus) for row in 1:genus, column in 1:genus]
    entries = [_hcm_eval(ring, variables, fox[row, column], basis) for row in 1:genus, column in 1:genus]
    dual = _hct_dual_exponents(hc, presented)
    partial1 = [_hcm_monomial(ring, variables, [basis[index, slot] for slot in 1:betti]) - 1 for index in 1:genus]
    partial3 = map(exponents -> begin
        power = zeros(Int, betti)
        for index in eachindex(exponents)
            exponents[index] == 0 && continue
            for slot in 1:betti
                power[slot] += exponents[index] * basis[index, slot]
            end
        end
        _hcm_monomial(ring, variables, power) - 1
    end, dual)
    got_num = nothing
    got_den = nothing
    for column in 1:genus, row in 1:genus
        (iszero(partial1[column]) || iszero(partial3[row])) && continue
        keep_rows = [index for index in 1:genus if index != row]
        keep_cols = [index for index in 1:genus if index != column]
        minor = if isempty(keep_rows)
            one(ring)
        else
            block = AbstractAlgebra.zero_matrix(ring, length(keep_rows), length(keep_cols))
            for (i, source_row) in enumerate(keep_rows), (j, source_col) in enumerate(keep_cols)
                block[i, j] = entries[source_row, source_col]
            end
            _hcm_det(block)
        end
        sign = isodd(column + row + genus + 1) ? -one(ring) : one(ring)
        numerator = sign * minor
        denominator = partial1[column] * partial3[row]
        if got_num === nothing
            got_num, got_den = numerator, denominator
        else
            _hct_require(_hcm_same_fraction(got_num, got_den, numerator, denominator), :minor_mismatch,
                "Fox minors disagree at row $row, column $column")
        end
    end
    got_num === nothing && throw(HCTorsionError(:minor_mismatch, "no Fox minor is defined"))
    reduced_num, reduced_den = _hcm_reduce(got_num, got_den)
    HCGroupRingFraction(betti, basis, reduced_num, reduced_den)
end

_hcm_same_fraction(left_num, left_den, right_num, right_den) =
    _hcm_unit_multiple(left_num * right_den, right_num * left_den)

"""
    hc_multivariable_alexander(hc::HeegaardCode) -> HCMultivariableAlexander
"""
function hc_multivariable_alexander(hc::HeegaardCode)
    torsion = hc_positive_torsion(hc)
    _hct_require(torsion.betti >= 2, :univariate, "b1 = $(torsion.betti)")
    _hct_require(iszero(torsion.numerator) || _hcm_is_unit(torsion.denominator), :not_laurent,
        "the b1 >= 2 torsion is not a Laurent polynomial")
    polynomial = iszero(torsion.numerator) ? torsion.numerator :
        AbstractAlgebra.divexact(torsion.numerator, torsion.denominator)
    HCMultivariableAlexander(torsion.betti, torsion.free_basis, polynomial)
end

function _hcm_polynomial(ring, terms)
    variables = [AbstractAlgebra.gen(ring, index) for index in 1:AbstractAlgebra.nvars(ring)]
    acc = zero(ring)
    for (exponents, coefficient) in terms
        acc += coefficient * _hcm_monomial(ring, variables, exponents)
    end
    acc
end

function _hcm_apply(value, matrix::Matrix{Int})
    ring = parent(value)
    variables = [AbstractAlgebra.gen(ring, index) for index in 1:size(matrix, 1)]
    acc = zero(ring)
    for (exponents, coefficient) in _hcm_terms(value)
        acc += coefficient * _hcm_monomial(ring, variables, matrix * exponents)
    end
    acc
end

function _hcm_det_int(matrix::Matrix{Int})
    n = size(matrix, 1)
    n == 1 && return matrix[1, 1]
    if n == 2
        return matrix[1, 1] * matrix[2, 2] - matrix[1, 2] * matrix[2, 1]
    end
    matrix[1, 1] * (matrix[2, 2] * matrix[3, 3] - matrix[2, 3] * matrix[3, 2]) -
        matrix[1, 2] * (matrix[2, 1] * matrix[3, 3] - matrix[2, 3] * matrix[3, 1]) +
        matrix[1, 3] * (matrix[2, 1] * matrix[3, 2] - matrix[2, 2] * matrix[3, 1])
end

"""
    hc_laurent_gl_matrix(polynomial, terms; bound::Int=2) -> Union{Matrix{Int}, Nothing}

Search GL(b, Z) with entries in `-bound:bound`. The result is a matrix A such that
the substituted polynomial and the term list differ by ± a monomial, where column
j of A is the image of the j-th variable. `nothing` means no such matrix was found
inside the bound.
"""
function hc_laurent_gl_matrix(polynomial, terms; bound::Int=2)
    ring = parent(polynomial)
    betti = AbstractAlgebra.nvars(ring)
    expected = _hcm_polynomial(ring, terms)
    found = Ref{Union{Nothing,Matrix{Int}}}(nothing)
    current = zeros(Int, betti, betti)
    function walk(column)
        found[] !== nothing && return
        if column > betti
            abs(_hcm_det_int(current)) == 1 || return
            if _hcm_unit_multiple(_hcm_apply(polynomial, current), expected)
                found[] = copy(current)
            end
            return
        end
        for entries in Iterators.product(ntuple(_ -> -bound:bound, betti)...)
            current[:, column] = collect(entries)
            walk(column + 1)
        end
    end
    walk(1)
    found[]
end

"""
    hc_turaev_matches(torsion::HCGroupRingFraction, alexander) -> Bool

`alexander` is a Laurent polynomial on the same basis. b₁ ≥ 2 requires
τ ∼ Δ, and b₁ = 1 requires τ ∼ Δ/(t−1)^2, both modulo ± monomials.
"""
function hc_turaev_matches(torsion::HCGroupRingFraction, alexander)
    parent(alexander) == parent(torsion.numerator) || return false
    if torsion.betti == 1
        variable = AbstractAlgebra.gen(parent(alexander), 1)
        return _hcm_unit_multiple(torsion.numerator * (variable - 1)^2, alexander * torsion.denominator)
    end
    torsion.betti >= 2 && _hcm_unit_multiple(torsion.numerator, alexander * torsion.denominator)
end

export HCGroupRingFraction, HCMultivariableAlexander, hc_positive_torsion,
    hc_multivariable_alexander, hc_laurent_gl_matrix, hc_turaev_matches
