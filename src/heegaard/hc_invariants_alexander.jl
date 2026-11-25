# cf: Nicolaescu, Notes on the Reidemeister Torsion, January 2002, Theorem 2.37, printed p. 69 — for b₁ = 1, T^aug ∼ Δ/(t−1)^2. The minor identity is Massuyeau, arXiv:1003.2517, Lemma 5.3, with τ_0 = +1. No registered Work key.
import AbstractAlgebra

"""
    HCAlexanderPolynomial(betti, coefficients) -> HCAlexanderPolynomial

Store a first Betti number and a univariate Alexander-polynomial coefficient record.

# Arguments
`betti::Int` and `coefficients::Vector{Int}`, or convertible values.
Conversion errors propagate. The constructor does not check the Betti
number, normalization or whether these coefficients arise from an HC.

# Returns
An `HCAlexanderPolynomial` retaining the supplied Betti number and
coefficient vector, in the format consumed by `hc_alexander_polynomial`.

# Guarantees
1. Fields preserve supplied values after conversion. A typed coefficient
   vector is shared and not modified. [proved/agent: by construction]
2. Construction stores the claimed polynomial data without performing an
   invariant computation. [proved/agent: by construction]

# Conventions
Alexander coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
Arbitrary coefficients need not be a manifold's Alexander polynomial.
The producer's general manifold identification remains subject to the
source-alignment limits of the linked specification; this constructor does
not establish it.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> a = HCAlexanderPolynomial(1, [1,-1,1]); (a.betti, Tuple(a.coefficients))
(1, (1, -1, 1))
```

See also [`hc_alexander_polynomial`](@ref).

Role: construct
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

Return the univariate Alexander-polynomial record computed from a closed HC presentation.

# Arguments
`hc` must have one alpha and one beta curve per handle and first Betti
number one. These are checked with `HCTorsionError` conditions `:not_closed`,
`:finite_homology` or `:multivariable`. Complete cut systems for the intended
3-manifold are the caller's responsibility under
[the Heegaard conditions](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
Invalid dual data, a non-Laurent quotient or disagreeing minors raise
`HCTorsionError` naming the failed condition; arithmetic errors propagate.

# Returns
`HCAlexanderPolynomial` with `betti=1` and an integer `coefficients` vector.

# Guarantees
1. The record has first Betti number one and a coefficient vector agreed
   by all defined Fox-minor comparisons. [proved/agent: by construction]
2. For the saved S¹×S² fixture, the polynomial is equivalent to 1.
   [verified-computationally: test/heegaard/hc_invariants_extra.jl, Alexander polynomial of S1xS2]
3. The input is not modified. [proved/agent: by construction]

# Conventions
Alexander coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
The general identification with the represented manifold's Alexander
polynomial remains subject to the source-alignment limits in
[the Alexander specification](../specs/heegaard-code.md#alexander-laurent-coordinates). This entry does not handle
first Betti number zero or greater than one; finite checks do not establish
the general identification.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#alexander-laurent-coordinates — Alexander coordinates, normalization and geometric-identification boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261114-usecases", "g1-S1xS2", "diagram.hc"); a = H.hc_alexander_polynomial(H.read_hc(path));

julia> (a.betti, Tuple(a.coefficients))
(1, (1,))
```

Role: invariant
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
