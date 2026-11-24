# cf: Nicolaescu, Notes on the Reidemeister Torsion, January 2002, Theorem 2.37, printed p. 69 — for b₁ = 1, T^aug ∼ Δ/(t−1)^2. The minor identity is Massuyeau, arXiv:1003.2517, Lemma 5.3, with τ_0 = +1. No registered Work key.
import AbstractAlgebra

"""Alexander polynomial of a closed HC with b₁ = 1, from the Heegaard Fox matrix.

`coefficients` lists the terms of degrees 0 through deg, lowest degree shifted to 0,
with positive leading coefficient. For the value Δ₀ produced by the formula below,
this representative ρ satisfies Δ₀(t) = ± t^k ρ(t). Against another normalization Δ_*,

    ρ(t) = ± t^k Δ_*(t^ε),   ε ∈ {±1}.

The formula, for a column i and a row j whose factors are nonzero, is

    Δ₀(t) = (−1)^{i+j+g+1} det(minor_{j,i}) (t − 1)^2 / ((φ(α_i) − 1)(φ(β_j^♯) − 1)).

φ is the projection of H₁ onto its free summand, sent to ⟨t⟩ so that the first nonzero
exponent among the α curves is positive. b₁ = 0 raises `HCTorsionError(:finite_homology)`.
b₁ > 1 raises `HCTorsionError(:multivariable)`.
"""
struct HCAlexanderPolynomial
    betti::Int
    coefficients::Vector{Int}
end

function _hca_eval(ring, variable, poly, weights)
    acc = zero(ring)
    for (exponents, coefficient) in poly
        power = 0
        for index in eachindex(exponents)
            power += exponents[index] * weights[index]
        end
        acc += coefficient * variable^power
    end
    acc
end

function _hca_bounds(value)
    iszero(value) && return nothing
    relative_high = AbstractAlgebra.degree(value.poly)
    relative_low = 0
    while relative_low <= relative_high && iszero(AbstractAlgebra.coeff(value.poly, relative_low))
        relative_low += 1
    end
    (value.mindeg + relative_low, value.mindeg + relative_high)
end

function _hca_same_class(left, right)
    iszero(left) && return iszero(right)
    iszero(right) && return false
    ok, quotient = AbstractAlgebra.divides(left, right)
    ok || return false
    bounds = _hca_bounds(quotient)
    bounds === nothing && return false
    low, high = bounds
    low == high || return false
    coefficient = AbstractAlgebra.coeff(quotient, low)
    coefficient == 1 || coefficient == -1
end

function _hca_coefficients(value)
    bounds = _hca_bounds(value)
    bounds === nothing && return Int[]
    low, high = bounds
    coefficients = Int[AbstractAlgebra.coeff(value, degree) for degree in low:high]
    if coefficients[end] < 0
        coefficients .*= -1
    end
    coefficients
end

"""
    hc_alexander_polynomial(hc::HeegaardCode) -> HCAlexanderPolynomial
"""
function hc_alexander_polynomial(hc::HeegaardCode)
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
    _hct_require(length(free) == 1, :multivariable, "b1 = $(length(free))")
    slot = only(free)
    weights = [right[row, slot] for row in 1:genus]
    sign_at = findfirst(!iszero, weights)
    _hct_require(sign_at !== nothing, :bad_differential, "the free generator is zero on every alpha")
    if weights[sign_at] < 0
        weights .*= -1
    end
    ring, variable = AbstractAlgebra.laurent_polynomial_ring(AbstractAlgebra.ZZ, :t)
    fox = [_hct_fox_entry(relators[row], column, genus) for row in 1:genus, column in 1:genus]
    entries = [_hca_eval(ring, variable, fox[row, column], weights) for row in 1:genus, column in 1:genus]
    dual = _hct_dual_exponents(hc, presented)
    partial1 = [variable^weights[index] - 1 for index in 1:genus]
    partial3 = map(exponents -> variable^_hct_dot(exponents, weights) - 1, dual)
    twist = (variable - 1)^2
    got = nothing
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
            AbstractAlgebra.det(block)
        end
        sign = isodd(column + row + genus + 1) ? -one(ring) : one(ring)
        numerator = sign * minor * twist
        denominator = partial1[column] * partial3[row]
        ok, value = AbstractAlgebra.divides(numerator, denominator)
        _hct_require(ok, :not_laurent, "the Fox quotient at row $row, column $column is not a Laurent polynomial")
        if got === nothing
            got = value
        else
            _hct_require(_hca_same_class(got, value), :minor_mismatch,
                "Fox minors disagree at row $row, column $column")
        end
    end
    got === nothing && throw(HCTorsionError(:minor_mismatch, "no Fox minor is defined"))
    HCAlexanderPolynomial(1, _hca_coefficients(got))
end

export HCAlexanderPolynomial, hc_alexander_polynomial
