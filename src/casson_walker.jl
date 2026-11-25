# Casson–Walker λ from surgery data. These functions do not read a Heegaard code.
# Normalization: [[Les07#^anonnd70-blk-0138]]. Draft: docs/conventions/CassonWalker.md.

"""
    sawtooth(x::Rational) -> Rational

Return the rational sawtooth value used in the Dedekind sum.

# Arguments
`x` satisfies CW1 ([conditions](../specs/casson-walker.md#numbered-conditions)) in the source's Dedekind-sum definition,
referenced by [the Casson–Walker convention](../conventions/CassonWalker.md).
Finiteness and representable arithmetic are caller obligations; exceptional
division, flooring or overflow behavior is not validated separately.

# Returns
A rational value with coefficient type determined by the input arithmetic.

# Guarantees
1. For finite input and supported arithmetic, the result is the source
   sawtooth function referenced by the convention.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
No floating-input method, overflow check or analytic continuation to an
infinite rational value is provided. The input domain is CW1.

# Sources
cf:[[Sav12#^sav12-sec-0540]] Dedekind-sum definition and the preceding sawtooth definition.

# Examples
```jldoctest
julia> sawtooth(1//3)^2 + sawtooth(2//3)^2
1//18
```

See also [`dedekind_sum`](@ref).

Role: invariant
"""
function sawtooth(x::Rational)
    isinteger(x) && return zero(x)
    return x - floor(x) - 1//2
end

"""
    dedekind_sum(b::Integer, a::Integer) -> Rational

Return the rational Dedekind-sum expression for the supplied integer pair.

# Arguments
`a,b` satisfy CW1–CW2 ([conditions](../specs/casson-walker.md#numbered-conditions)) and have the meaning in the
[Dedekind convention](../conventions/CassonWalker.md).
Positive `a` is checked; failure raises `ArgumentError` stating that the
second argument must be positive. Coprimality is not checked: noncoprime
inputs still evaluate the same finite sum. Arithmetic must be representable;
fixed-width intermediate products can overflow, so callers requiring large
exact values use `BigInt`.

# Returns
A rational scalar. The empty-sum case has value `0//1`.

# Guarantees
1. With supported arithmetic, the value is the finite sum in the linked
   convention, extended by that same expression to noncoprime pairs.
   [proved/agent: by construction]
2. The first argument may be replaced by a congruent value modulo the
   second without changing the result when arithmetic is representable.
   [proved/agent: by construction]
3. The input arguments are not modified. [proved/agent: by construction]

# Not guaranteed
This does not certify coprimality, a lens-space presentation, or a reciprocity
theorem outside its required source hypotheses. There is no fixed-width
overflow check or sublinear run-time promise.

# Sources
cf:[[Sav12#^sav12-sec-0540]] finite Dedekind sum, on the source's coprime domain; the implementation also evaluates the expression without that restriction.

# Examples
```jldoctest
julia> dedekind_sum(2, 7)
1//14
```

See also [`sawtooth`](@ref), [`casson_walker_lens`](@ref).

Role: invariant
"""
function dedekind_sum(b::Integer, a::Integer)
    a > 0 || throw(ArgumentError("the second argument of a Dedekind sum is positive"))
    bb = mod(b, a)
    s = 0//1
    for k in 1:(a - 1)
        s += sawtooth(k//a) * sawtooth((k * bb)//a)
    end
    return s
end

"""
    casson_walker_lens(p::Integer, q::Integer) -> Rational

Return the Casson–Walker lens-space expression in the repository's surgery convention.

# Arguments
`p,q` satisfy CW1 and CW3 in the
[Casson–Walker specification §1](../specs/casson-walker.md#1-関数).
The positivity and coprimality checks raise `ArgumentError` naming the
failed condition. Arithmetic must remain representable; fixed-width
overflow is not checked, so callers needing large exact values use `BigInt`.

# Returns
A rational scalar in the
[lens naming and normalization convention](../conventions/CassonWalker.md).

# Guarantees
1. On supported arithmetic, the result is the lens expression specified
   in that convention for the supplied surgery data.
   [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker normalization and lens naming](../conventions/CassonWalker.md).

# Not guaranteed
The function does not identify a lens space from another representation or
translate orientation conventions used by other modules. The owning
convention page remains draft. Overflowed arithmetic is not an exact value.

# Sources
cf:[[Les07#^anonnd70-blk-0314]] surgery formula specialized to the unknot; the repository uses the linked lens naming.
cf:[[Les07#^anonnd70-blk-0138]] normalization relative to the Walker invariant.

# Examples
```jldoctest
julia> casson_walker_lens(3, 1)
-1//36
```

See also [`casson_walker_knot_surgery`](@ref), [`dedekind_sum`](@ref).

Role: invariant
"""
function casson_walker_lens(p::Integer, q::Integer)
    p > 0 || throw(ArgumentError("lens order p is positive in this naming"))
    gcd(p, abs(q)) == 1 || throw(ArgumentError("p and q are coprime"))
    return -dedekind_sum(q, p) // 2
end

"""
    casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Union{Integer,Rational}) -> Rational

Return the Casson–Walker surgery expression for a knot in the three-sphere.

# Arguments
- `p,q` satisfy CW1 and CW4 ([conditions](../specs/casson-walker.md#numbered-conditions)) in the
  [surgery convention](../conventions/CassonWalker.md). A zero `q` or
  noncoprime pair raises `ArgumentError` identifying that condition.
  A zero `p` is outside the finite rational-homology-sphere interpretation;
  the caller must exclude it. Its eventual failure is not a dedicated
  zero-slope diagnostic. Negative `q` is accepted.
- `delta2` satisfies CW5 in the
  [surgery derivative normalization](../specs/knot-classical-invariants.md#surgery-derivative).
  Its knot origin and normalization are caller obligations, not checked.
- All integer/rational arithmetic must be representable. Overflow is not
  checked; use `BigInt` data for large exact calculations.

# Returns
A rational scalar in the linked Casson–Walker convention.

# Guarantees
1. On supported arithmetic, the result equals the surgery expression in
   [Casson–Walker §1–§2](../specs/casson-walker.md) for the supplied data.
   With the source's knot and normalization hypotheses it has the stated
   surgery interpretation. [proved/agent: by construction]
2. Simultaneously negating `p,q` leaves the result unchanged when the
   negations and subsequent arithmetic are representable.
   [proved/agent: by construction]
3. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker normalization and surgery slope](../conventions/CassonWalker.md).

# Not guaranteed
This call neither obtains `delta2` from a knot nor checks that a supplied
polynomial has the source normalization. Differentiating the stored output
of `alexander_polynomial` without the linked translation is not justified.
No convention translation for an arbitrary triangulation or Heegaard
presentation is supplied. The owning convention page remains draft.

# Sources
cf:[[Les07#^anonnd70-blk-0313]] Alexander normalization required by the surgery formula.
cf:[[Les07#^anonnd70-blk-0314]] knot-surgery formula and its S3 specialization.
cf:[[Les07#^anonnd70-blk-0138]] Casson–Walker normalization.

# Examples
```jldoctest
julia> casson_walker_knot_surgery(7, 1, -2)
-9//28
```

See also [`casson_walker_lens`](@ref), [`casson_walker_surgery_pup`](@ref).

Role: invariant
"""
casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Integer) =
    casson_walker_knot_surgery(p, q, delta2//1)

"""
    casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Rational) -> Rational

Return the Casson–Walker surgery expression for a knot in the three-sphere.

# Arguments
- `p,q` satisfy CW1 and CW4 ([conditions](../specs/casson-walker.md#numbered-conditions)) in the
  [surgery convention](../conventions/CassonWalker.md). A zero `q` or
  noncoprime pair raises `ArgumentError` identifying that condition.
  A zero `p` is outside the finite rational-homology-sphere interpretation;
  the caller must exclude it. Its eventual failure is not a dedicated
  zero-slope diagnostic. Negative `q` is accepted.
- `delta2` satisfies CW5 in the
  [surgery derivative normalization](../specs/knot-classical-invariants.md#surgery-derivative).
  Its knot origin and normalization are caller obligations, not checked.
- All integer/rational arithmetic must be representable. Overflow is not
  checked; use `BigInt` data for large exact calculations.

# Returns
A rational scalar in the linked Casson–Walker convention.

# Guarantees
1. On supported arithmetic, the result equals the surgery expression in
   [Casson–Walker §1–§2](../specs/casson-walker.md) for the supplied data.
   With the source's knot and normalization hypotheses it has the stated
   surgery interpretation. [proved/agent: by construction]
2. Simultaneously negating `p,q` leaves the result unchanged when the
   negations and subsequent arithmetic are representable.
   [proved/agent: by construction]
3. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker normalization and surgery slope](../conventions/CassonWalker.md).

# Not guaranteed
This call neither obtains `delta2` from a knot nor checks that a supplied
polynomial has the source normalization. Differentiating the stored output
of `alexander_polynomial` without the linked translation is not justified.
No convention translation for an arbitrary triangulation or Heegaard
presentation is supplied. The owning convention page remains draft.

# Sources
cf:[[Les07#^anonnd70-blk-0313]] Alexander normalization required by the surgery formula.
cf:[[Les07#^anonnd70-blk-0314]] knot-surgery formula and its S3 specialization.
cf:[[Les07#^anonnd70-blk-0138]] Casson–Walker normalization.

# Examples
```jldoctest
julia> casson_walker_knot_surgery(7, 1, -2//1)
-9//28
```

See also [`casson_walker_lens`](@ref), [`casson_walker_surgery_pup`](@ref).

Role: invariant
"""
function casson_walker_knot_surgery(p::Integer, q::Integer, delta2::Rational)
    q == 0 && throw(ArgumentError("q = 0 is not a finite surgery coefficient"))
    if q < 0
        p = -p
        q = -q
    end
    gcd(abs(p), q) == 1 || throw(ArgumentError("p and q are coprime"))
    knot = (q//p) * (delta2 // 2)
    lens = p > 0 ? casson_walker_lens(p, q) : -casson_walker_lens(-p, q)
    return knot + lens
end

"""
    casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Union{Integer,Rational}) -> Rational

Return the positive-denominator expansion of the Casson–Walker knot-surgery expression.

# Arguments
`p,q,delta2` satisfy CW1 and CW4–CW6 ([conditions](../specs/casson-walker.md#numbered-conditions)), with the meaning and caller obligations of
[`casson_walker_knot_surgery`](@ref), in the
[surgery convention](../conventions/CassonWalker.md).
This expansion additionally checks positive `q`; failure raises
`ArgumentError` stating that restriction. Noncoprime inputs raise
`ArgumentError`. A nonzero `p`, correctly normalized derivative and
representable arithmetic are caller obligations; division or overflow on
invalid inputs is not a valid surgery result.

# Returns
A rational scalar in the linked Casson–Walker normalization.

# Guarantees
1. For supported arithmetic, the value is the positive-denominator
   expansion identified in [Casson–Walker §2](../specs/casson-walker.md#2-一致).
   [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker surgery and normalization](../conventions/CassonWalker.md).

# Not guaranteed
The function does not validate the knot or derivative normalization, accept
a zero surgery numerator as a rational-homology-sphere calculation, or
translate another module's orientation convention. Fixed-width overflow is
unchecked, and the owning convention remains draft.

# Sources
cf:[[Les07#^anonnd70-blk-0314]] positive-denominator expansion in the proof of the surgery lemma.

# Examples
```jldoctest
julia> casson_walker_surgery_pup(-7, 1, -2)
9//28
```

See also [`casson_walker_knot_surgery`](@ref).

Role: invariant
"""
casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Integer) =
    casson_walker_surgery_pup(p, q, delta2//1)

"""
    casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Rational) -> Rational

Return the positive-denominator expansion of the Casson–Walker knot-surgery expression.

# Arguments
`p,q,delta2` satisfy CW1 and CW4–CW6 ([conditions](../specs/casson-walker.md#numbered-conditions)), with the meaning and caller obligations of
[`casson_walker_knot_surgery`](@ref), in the
[surgery convention](../conventions/CassonWalker.md).
This expansion additionally checks positive `q`; failure raises
`ArgumentError` stating that restriction. Noncoprime inputs raise
`ArgumentError`. A nonzero `p`, correctly normalized derivative and
representable arithmetic are caller obligations; division or overflow on
invalid inputs is not a valid surgery result.

# Returns
A rational scalar in the linked Casson–Walker normalization.

# Guarantees
1. For supported arithmetic, the value is the positive-denominator
   expansion identified in [Casson–Walker §2](../specs/casson-walker.md#2-一致).
   [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker surgery and normalization](../conventions/CassonWalker.md).

# Not guaranteed
The function does not validate the knot or derivative normalization, accept
a zero surgery numerator as a rational-homology-sphere calculation, or
translate another module's orientation convention. Fixed-width overflow is
unchecked, and the owning convention remains draft.

# Sources
cf:[[Les07#^anonnd70-blk-0314]] positive-denominator expansion in the proof of the surgery lemma.

# Examples
```jldoctest
julia> casson_walker_surgery_pup(-7, 1, -2//1)
9//28
```

See also [`casson_walker_knot_surgery`](@ref).

Role: invariant
"""
function casson_walker_surgery_pup(p::Integer, q::Integer, delta2::Rational)
    q > 0 || throw(ArgumentError("this expansion is stated for q > 0"))
    gcd(abs(p), q) == 1 || throw(ArgumentError("p and q are coprime"))
    sign = p > 0 ? 1 : -1
    return (q//p) * (delta2 // 2 - 1//24 - (p * p + 1)//(24 * q * q)) +
           sign//8 + dedekind_sum(p, q) // 2
end

"""
    casson_walker_brieskorn(p::Integer, q::Integer, r::Integer) -> Rational

Return the Casson expression for the Brieskorn sphere with the supplied multiplicities.

# Arguments
`p,q,r` satisfy CW1 and CW7 ([conditions](../specs/casson-walker.md#numbered-conditions)). Positivity is checked; failure raises
`ArgumentError` with `Seifert invariants are positive`. For the Brieskorn
homology-sphere interpretation, the caller must supply the hypotheses of
the theorem cited by [Casson–Walker §1](../specs/casson-walker.md#1-関数);
pairwise coprimality and that geometric interpretation are not checked.
Arithmetic must be representable: fixed-width powers/products can overflow
without a diagnostic. Use `BigInt` for large exact calculations.

# Returns
A rational scalar in the source orientation and the repository's
[Casson–Walker normalization](../conventions/CassonWalker.md).

# Guarantees
1. On supported arithmetic, the result is the expression of Sav12
   Theorem 19.4 referenced in the specification. Its interpretation as
   the Casson invariant is conditional on that theorem's input domain.
   [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Conventions
[Casson–Walker normalization](../conventions/CassonWalker.md).

# Not guaranteed
Positive inputs alone do not certify a Brieskorn homology sphere. The value
outside the source domain is only the evaluated formula, not a proved
invariant of a space. This call does not translate orientations from another
manifold representation, and the owning convention remains draft.

# Sources
cf:[[Sav12#^sav12-sec-0544]] Theorem 19.4 — the Brieskorn homology-sphere expression.
cf:[[Sav12#^sav12-sec-0540]] Dedekind-sum normalization entering that expression.

# Examples
```jldoctest
julia> casson_walker_brieskorn(2, 3, 5)
-1//1
```

See also [`casson_walker_seifert`](@ref).

Role: invariant
"""
function casson_walker_brieskorn(p::Integer, q::Integer, r::Integer)
    (p > 0 && q > 0 && r > 0) || throw(ArgumentError("Seifert invariants are positive"))
    term = (1 - (p * q * r)^2 + (p * q)^2 + (q * r)^2 + (p * r)^2)//(3 * p * q * r)
    d(a, b) = -4 * dedekind_sum(b, a)
    bracket = 1 - term - d(p, q * r) - d(q, p * r) - d(r, p * q)
    return -bracket // 8
end

"""
    casson_walker_seifert(b::Integer, fibers::Vector{Tuple{Int,Int}}) -> Union{Nothing,NamedTuple}

Return the repository's Casson–Walker and homology-order expressions for Seifert data.

# Arguments
`b,fibers` satisfy CW1 and CW8 ([conditions](../specs/casson-walker.md#numbered-conditions)) in
[Casson–Walker §1](../specs/casson-walker.md#1-関数).
Each fiber multiplicity must be nonzero and each pair coprime; these are
checked, with `ArgumentError` naming the failed condition. A geometric
Seifert presentation in the owning convention is the caller's responsibility.
The vector element type is specifically `Tuple{Int,Int}`. Products of its
entries use fixed-width arithmetic; preventing overflow remains the caller's
responsibility even when `b` is a `BigInt`.

# Returns
`nothing` for the zero-Euler case. Otherwise a named tuple `(lambda,h1,E)`
of rational values, with the field meanings and formulas defined in
[Casson–Walker §1](../specs/casson-walker.md#1-関数).

# Guarantees
1. On representable input arithmetic, the returned fields equal those
   specified algebraic expressions. [proved/agent: by construction]
2. The input vector and its entries are not modified.
   [proved/agent: by construction]

# Conventions
[Seifert expression and normalization](../conventions/CassonWalker.md).

# Not guaranteed
The general geometric identification of this Seifert formula is unverified
in the owning convention page; its listed calibration examples do not
establish it for every Seifert presentation. `nothing` is a domain outcome,
not a claim that the invariant is zero. The `h1` field is an order expression,
not the full first-homology group or a linking form.

# Sources
cf:lib/topology.jl/docs/specs/casson-walker.md §1 — the current Seifert expression and field definitions.
cf:[[Sav12#^sav12-sec-0544]] Brieskorn calibration values; not a source for the unrestricted Seifert formula.

# Examples
```jldoctest
julia> r = casson_walker_seifert(-1, [(2,1),(3,1),(5,1)]);

julia> (r.lambda, r.h1, r.E)
(-1//1, 1//1, -1//30)
```

See also [`casson_walker_brieskorn`](@ref), [`casson_walker_lens`](@ref).

Role: invariant
"""
function casson_walker_seifert(b::Integer, fibers::Vector{Tuple{Int,Int}})
    E = -b//1
    ssum = 0//1
    quad = (2 - length(fibers))//1
    order = 1
    for (alpha, beta) in fibers
        alpha == 0 && throw(ArgumentError("a fiber multiplicity is nonzero"))
        gcd(abs(alpha), abs(beta)) == 1 || throw(ArgumentError("a fiber pair is coprime"))
        E -= beta//alpha
        quad += 1//(alpha * alpha)
        order *= abs(alpha)
        s = dedekind_sum(beta, abs(alpha))
        ssum += 12 * (alpha > 0 ? 1 : -1) * s
    end
    E == 0 && return nothing
    inside = 3 * (E > 0 ? 1 : -1) - E - ssum - quad / E
    return (lambda = inside / 24, h1 = abs(E) * order, E = E)
end
