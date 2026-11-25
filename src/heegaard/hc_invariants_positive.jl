# Massuyeau, arXiv:1003.2517, Theorem 4.4, PDF p. 19: for a closed 3-manifold,
# β₁ ≥ 2 gives ±τ^μ = Δ in Z[G]/±G, and β₁ = 1 gives ±τ^μ = Δ/(t−1)^2.
# Lemma 4.3, PDF pp. 17–18, with τ₀ = +1, is the Heegaard Fox quotient below.
# b₁ = 0 stays on hc_reidemeister_torsion and is not handled here.
import AbstractAlgebra

"""
    HCGroupRingFraction(betti, free_basis, numerator, denominator) -> HCGroupRingFraction

Store a rational Laurent-expression record for positive-Betti Heegaard torsion.

# Arguments
`betti::Int`, `free_basis::Matrix{Int}`, and arbitrary numerator and denominator
values, with field conversions where applicable. Conversion errors propagate.
The caller supplies consistent coordinates and compatible ring elements;
the constructor checks neither dimensions nor a nonzero denominator.

# Returns
An `HCGroupRingFraction` with the supplied Betti number, free-coordinate
matrix and fraction components.

# Guarantees
1. Fields preserve supplied values after conversion. Typed arrays and ring
   values are retained without modification. [proved/agent: by construction]
2. Construction stores claimed torsion data without computing a Fox complex
   or proving a manifold identification. [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
A hand-built fraction need not define a valid torsion. The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the producer's basis and normalization; comparing raw fractions across HC presentations is not certified here.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> f = HCGroupRingFraction(1, ones(Int,1,1), 1, 2); (f.betti, f.numerator, f.denominator)
(1, 1, 2)
```

See also [`hc_positive_torsion`](@ref).

Role: construct
"""
struct HCGroupRingFraction
    betti::Int
    free_basis::Matrix{Int}
    numerator
    denominator
end

"""
    HCMultivariableAlexander(betti, free_basis, polynomial) -> HCMultivariableAlexander

Store a multivariable Alexander-polynomial record with its free-coordinate matrix.

# Arguments
`betti::Int`, `free_basis::Matrix{Int}` and an arbitrary polynomial value.
Field conversions are applied and conversion errors propagate. The caller
supplies compatible coordinates and a suitable polynomial; the constructor
checks neither the Betti number nor the coefficient ring.

# Returns
An `HCMultivariableAlexander` retaining the supplied Betti number,
free-coordinate matrix and polynomial.

# Guarantees
1. Fields preserve supplied values after conversion. Typed containers and
   the polynomial are retained without mutation. [proved/agent: by construction]
2. No invariant computation or manifold certification is performed.
   [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
Arbitrary fields need not describe a manifold's Alexander polynomial.
The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the producer's coordinates and normalization; raw polynomial equality across presentations
is not an invariant comparison supplied by this constructor.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> a = HCMultivariableAlexander(2, [1 0; 0 1], 1); (a.betti, size(a.free_basis), a.polynomial)
(2, (2, 2), 1)
```

See also [`hc_multivariable_alexander`](@ref).

Role: construct
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

Return a Laurent-fraction torsion record for the free abelian quotient of an HC presentation.

# Arguments
`hc` has the [Heegaard input domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
Closed curve counts and positive first Betti number are checked with
`HCTorsionError(:not_closed)` or `HCTorsionError(:finite_homology)`.
Invalid dual data, absent or disagreeing minors and a zero denominator raise
stage-specific `HCTorsionError`; algebra errors propagate. The caller supplies
complete cut systems for a geometric manifold interpretation.

# Returns
`HCGroupRingFraction` with `betti`, integer `free_basis`, and reduced Laurent
`numerator` and `denominator`.

# Guarantees
1. All defined Fox-minor fractions agree with the returned quotient up to
   a signed Laurent monomial. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the positive-Betti basis, augmentation and normalization. The record's general identification with the represented
3-manifold's Reidemeister torsion is not independently certified here. It is
not an absolute value invariant under arbitrary basis changes.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#reidemeister-character-coordinates — Fox-complex identification boundary; positive-Betti normalization requires its own entry.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> try H.hc_positive_torsion(hc) catch e; e isa H.HCTorsionError && e.reason == :finite_homology end
true
```

Role: invariant
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

Return the multivariable Laurent-polynomial record associated with a positive-Betti HC.

# Arguments
`hc` has the domain of `hc_positive_torsion`, including closed curve counts
and positive first Betti number. Its errors propagate. Checked: Betti number
is at least two, otherwise `HCTorsionError(:univariate)`; nonzero torsion must
have a Laurent-unit denominator, otherwise `HCTorsionError(:not_laurent)`.
Complete cut systems are the caller's responsibility in the
[Heegaard domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).

# Returns
`HCMultivariableAlexander` with `betti`, the producer's `free_basis` and a
Laurent `polynomial` in that basis.

# Guarantees
1. The polynomial equals the positive-torsion quotient on the reported
   free coordinates, including a zero numerator. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the free-basis and unit convention. General identification
with the represented manifold's Alexander polynomial is not certified by this
wrapper, and literal coefficient comparison across bases is not justified.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#reidemeister-character-coordinates — related Fox-complex identification boundary; positive-Betti conventions remain separate.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> try H.hc_multivariable_alexander(hc) catch e; e isa H.HCTorsionError && e.reason == :finite_homology end
true
```

Role: invariant
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
    hc_laurent_gl_matrix(polynomial, terms; bound::Int=2) -> Union{Matrix{Int},Nothing}

Search for a bounded unimodular change of Laurent-polynomial coordinates.

# Arguments
`polynomial` is an AbstractAlgebra integer Laurent polynomial and `terms`
is an iterable of integer exponent-vector/coefficient pairs in the same
number of variables. The caller must restrict the number of variables to
one, two or three; that determinant-domain restriction is not checked.
A nonnegative `bound` is the caller's responsibility. Integer arithmetic must
fit in `Int`. Ring, exponent and conversion failures propagate.

# Returns
A unimodular integer matrix witnessing unit-equivalence to the term-list
polynomial, or `nothing` if none was found inside the supplied entry bound.

# Guarantees
1. Within the stated dimension and arithmetic domain, a returned matrix
   witnesses equality up to a signed Laurent monomial after substitution.
   [proved/agent: by construction]
2. Neither input is modified. [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
A failed bounded search does not rule out an unbounded change of basis.
The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the exponent transport and unit comparison;
this numerical witness alone does not identify manifold invariants.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#alexander-laurent-coordinates — exponent transport and Laurent-unit comparison.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes; A = H.AbstractAlgebra;

julia> r, v = A.laurent_polynomial_ring(A.ZZ, [:u,:v]); p = one(r) + v[2]^2;

julia> H.hc_laurent_gl_matrix(p, [(Int[0,0],1),(Int[2,0],1)]; bound=1) !== nothing
true
```

Role: other
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

Test the supported torsion–Alexander identity on common Laurent coordinates.

# Arguments
`alexander` must be a Laurent polynomial; `torsion` must have consistent
positive-Betti fraction fields. A different parent ring returns false.
The caller ensures that both objects use the same free basis and that the
denominator is nonzero; malformed ring data and exact-division errors propagate.

# Returns
A Boolean reporting equality in the supported torsion–Alexander unit class.

# Guarantees
1. In the stated domain the result tests the Betti-one or higher-Betti
   torsion–Alexander relation implemented for these records.
   [proved/agent: by construction]
2. Neither input is modified. [proved/agent: by construction]

# Conventions
Laurent coordinates and normalization: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-code.md#alexander-laurent-coordinates).

# Not guaranteed
Matching is not a manifold identification or a proof that either input is
its claimed invariant. The [Laurent-coordinate specification](../specs/heegaard-code.md#alexander-laurent-coordinates)
states the normalization, basis transport and comparison relation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md#reidemeister-character-coordinates — related torsion-identification boundary; positive-Betti normalization remains separate.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> A = H.AbstractAlgebra; r, v = A.laurent_polynomial_ring(A.ZZ, [:t]); t = only(v); f = H.HCGroupRingFraction(1, ones(Int,1,1), one(r), (t-1)^2);

julia> H.hc_turaev_matches(f,one(r))
true
```

Role: condition
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
