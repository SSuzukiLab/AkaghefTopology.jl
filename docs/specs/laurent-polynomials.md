# Laurent polynomials and rational functions

## Definitions

### Polynomial data

For an ordered list of variables `x_1,...,x_m`, a Laurent polynomial over the
rationals is a finite sum of monomials with integer exponent vectors and
rational coefficients. `LaurentPolynomial` stores that ordered list as
strings and the coefficient map as `Dict{Tuple{Vararg{Int}},Rational{BigInt}}`.
The zero polynomial has an empty map. Normalized data omit zero coefficients.
These are local coordinate definitions in [laurent.jl](../../src/laurent.jl),
not an identification with a knot polynomial or a literature invariant.

### Fraction data and equality

`LaurentRational` is an ordered numerator/denominator pair in the same Laurent
polynomial ring. A mathematical fraction requires a nonzero denominator.
The raw two-field constructor does not impose these requirements.
Arithmetic normalization checks the variable lists and denominator,
cancels univariate common factors, and makes the denominator monic after
exponent adjustment. Multivariate normalization performs exact divisibility
and supported binomial cancellation, then scales the leading denominator
coefficient to one. It is not a complete multivariate factorization.
Stored equality is equality of these coordinates, not a decision procedure
for arbitrary presentations of the same rational function.

### Evaluation and zero limit

Evaluation substitutes a value for every variable in the stored order, with
exact `Rational{BigInt}` coercion. A zero limit is the finite univariate limit
at zero. Its algebraic criterion is the lowest nonzero exponent: a negative
valuation is a pole, a positive one gives zero, and zero valuation gives the
lowest coefficient. For a fraction, subtract denominator valuation and use
the ratio of lowest coefficients. This definition applies to valid nonzero
denominators, including a zero numerator with limit zero.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| L1 | Ordered variable names are distinct; each exponent tuple has the same length as the variable list. | Raw fields are unchecked; callers must preserve this representation. |
| L2 | Numerator and denominator use identical ordered variables, with denominator nonzero. | Arithmetic normalization checks; raw fraction construction does not. |
| L3 | A substitution supplies one value per variable and all evaluated negative powers/divisions are defined. | Count is checked; coefficient arithmetic determines singular behavior. |
| L4 | A zero-limit input is univariate and has finite valuation at zero (nonnegative valuation difference for fractions). | Variable count and poles are checked; raw malformed data are outside the guarantee. |

## Coordinate choices and evidence

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.
Variable order is the caller's order; two permutations of variable names are
different coordinate rings unless translated. The map between such orders
permutes the exponent-vector coordinates in the same way. No automatic
variable permutation is performed. Named variables use the supplied prefix
and increasing index. Univariate normalization and the current multivariate
identities are embodied by [topology_core.jl](../../test/topology_core.jl),
`C251020 rank-one Laurent workflow`.

question: open — Full multivariate canonical equality is not established.
The current cancellations establish the tested identities, but completeness
would require a canonical multivariate gcd/normal-form argument. No such
argument or general equality algorithm is supplied.

## Operations

`LaurentPolynomial`, `LaurentRational`, `laurent_variable`,
`laurent_variables`, `evaluate`, `limit_zero`, `+`, `-`, `*`, `/`, `^`,
`zero`, `one`, `iszero`, `==`.
