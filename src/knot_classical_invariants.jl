"""Exact invariants and PD bridges for one-component classical knot diagrams.

This file is intentionally self-contained at the public boundary: a caller may
include it after `knot_gauss.jl` (or after `Topology` has been loaded), and it
does not change the package module.  The two invariant implementations share
only the input data, not a KnotInfo value or a table lookup.

Conventions
-----------

`KnotGaussDiagram.gauss` is a cyclic traversal word.  A positive visit is an
over-crossing and a negative visit is an under-crossing.  `orientation[c]` is
the oriented crossing sign, +1 or -1.  The Gauss-to-PD bridge uses one edge
between consecutive visits, with the last edge closed back to the first.  A
PD row is cyclically ordered counter-clockwise and has the two strands in
slots (1,3) and (2,4); (1,3) is under and (2,4) is over.  For a crossing sign
`-1` the row is `[u_in,o_in,u_out,o_out]`; for `+1` it is
`[u_in,o_out,u_out,o_in]`.

For the Alexander polynomial, the Wirtinger/Fox row for crossing `c` is
`(1-t^s)x_over + t^s*x_under_in - x_under_out`, where `s` is its crossing
sign.  A last row and column are deleted and the result is normalized up to
`+/- t^k` by shifting its least exponent to zero and making the leading
coefficient positive.

For Jones, the PD state sum uses A-smoothing `(1,2)(3,4)`, B-smoothing
`(2,3)(4,1)`, `delta=-A^2-A^-2`, and
`V_D(t)=(-A^3)^(-w(D))*<D>` with `t=A^-4`, where
`w(D)=sum(orientation)`.  The implementation rejects a residual A exponent;
for a valid oriented one-component diagram all normalized exponents are in a
single class modulo four.
"""

if !isdefined(@__MODULE__, :KnotGaussDiagram)
    include(joinpath(@__DIR__, "knot_gauss.jl"))
end

"""    ExactKnotLaurentPolynomial(variable::Symbol, terms::AbstractDict) -> ExactKnotLaurentPolynomial

Construct an exact integer Laurent-polynomial record.

# Arguments
`terms` supplies exponent/coefficient pairs convertible to `Int`/`BigInt`;
conversion failures propagate (`InexactError`, `MethodError` or other numeric
conversion errors). `variable` is a `Symbol`. The one-dictionary overload accepts keyword `variable=:t`. The representation is defined
in [Exact polynomial records](../specs/knot-classical-invariants.md#exact-polynomial-records).
For a normalized record the caller must omit zero coefficients when passing
an already typed `Dict{Int,BigInt}`: that field constructor does not clean it.

# Returns
An `ExactKnotLaurentPolynomial` with the supplied `variable` and a `terms`
dictionary denoting the linked integer Laurent polynomial.

# Guarantees
1. For the general dictionary conversion, coefficients are converted exactly to `BigInt`, exponents to `Int`, coincident converted exponents are combined, and zero coefficients are omitted. [proved/agent: by construction]
2. The `Dict{Int,BigInt}` field-constructor method retains the supplied dictionary unchanged; other dictionary inputs are converted to a new dictionary. Inputs are not modified during construction. [proved/agent: by construction]

# Not guaranteed
No knot, Alexander or Jones interpretation is assigned by construction. A
variable symbol is not checked for a printable mathematical name. The typed
field constructor does not promise a canonical zero-free representation.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, Exact polynomial records — local coefficient-ring representation.

# Examples
```jldoctest
julia> p = ExactKnotLaurentPolynomial(Dict(0 => 1));

julia> p == alexander_polynomial(pd_to_gauss(Int[]))
true
```

See also [`polynomial_terms`](@ref), [`knot_polynomial_string`](@ref).

Role: construct
"""
struct ExactKnotLaurentPolynomial
    variable::Symbol
    terms::Dict{Int,BigInt}
end

function ExactKnotLaurentPolynomial(variable::Symbol,
                                    terms::AbstractDict)
    clean = Dict{Int,BigInt}()
    for (exponent, coefficient) in terms
        c = BigInt(coefficient)
        c == 0 && continue
        e = Int(exponent)
        clean[e] = get(clean, e, big(0)) + c
        clean[e] == 0 && delete!(clean, e)
    end
    ExactKnotLaurentPolynomial(variable, clean)
end

ExactKnotLaurentPolynomial(terms::AbstractDict; variable::Symbol=:t) =
    ExactKnotLaurentPolynomial(variable, terms)

Base.copy(p::ExactKnotLaurentPolynomial) =
    ExactKnotLaurentPolynomial(p.variable, copy(p.terms))
Base.:(==)(a::ExactKnotLaurentPolynomial, b::ExactKnotLaurentPolynomial) =
    a.variable == b.variable && a.terms == b.terms
Base.hash(p::ExactKnotLaurentPolynomial, h::UInt) = hash((p.variable, p.terms), h)
Base.iszero(p::ExactKnotLaurentPolynomial) = isempty(p.terms)

"""
    knot_polynomial_string(p::ExactKnotLaurentPolynomial) -> String

Return a readable expression for an exact integer Laurent polynomial.

# Arguments
`p` is an `ExactKnotLaurentPolynomial`. For the usual polynomial display,
its terms should omit zero coefficients; this function does not normalize
manually assembled data or validate the variable name.

# Returns
A `String` expressing the stored coefficients and exponents in `p.variable`.
The polynomial interface belongs to
[KnotGauss](../specs/knot-classical-invariants.md#definitions).

# Guarantees
1. Nonzero stored terms are rendered with their exact integer coefficients
   and exponents, ordered by increasing exponent. The empty dictionary
   renders as `"0"`. [proved/agent: by construction]
2. The expression denotes the stored Laurent polynomial when its variable
   name has an unambiguous mathematical interpretation. The input is unchanged. [proved/agent: by construction]

# Not guaranteed
No Alexander or Jones normalization is performed here. Text output is not
a versioned interchange format or an identity/classification certificate;
no escaping or parser round trip is promised for arbitrary variable names.

# Sources
cf:lib/topology.jl/src/knot_classical_invariants.jl, knot_polynomial_string — stored polynomial and PD data access; no independent knot classification.

# Examples
```jldoctest
julia> p = alexander_polynomial(pd_to_gauss(Int[]));

julia> knot_polynomial_string(p)
"1"
```

See also [`polynomial_terms`](@ref).

Role: convert
"""
function knot_polynomial_string(p::ExactKnotLaurentPolynomial)
    isempty(p.terms) && return "0"
    variable = String(p.variable)
    pieces = String[]
    for exponent in sort(collect(keys(p.terms)))
        coefficient = p.terms[exponent]
        negative = coefficient < 0
        magnitude = abs(coefficient)
        body = if exponent == 0
            string(magnitude)
        elseif exponent == 1
            magnitude == 1 ? variable : "$(magnitude)*$(variable)"
        else
            magnitude == 1 ? "$(variable)^$(exponent)" :
                "$(magnitude)*$(variable)^$(exponent)"
        end
        if isempty(pieces)
            push!(pieces, negative ? "-$(body)" : body)
        else
            push!(pieces, negative ? "- $(body)" : "+ $(body)")
        end
    end
    join(pieces, " ")
end

Base.show(io::IO, p::ExactKnotLaurentPolynomial) =
    print(io, knot_polynomial_string(p))

"""
    polynomial_terms(p::ExactKnotLaurentPolynomial) -> Dict{Int,BigInt}

Return a copy of the exponent-to-coefficient data of an integer Laurent polynomial.

# Arguments
`p` is an `ExactKnotLaurentPolynomial`; no additional validation is performed.

# Returns
A new `Dict{Int,BigInt}` with the stored exponents as keys and exact integer
coefficients as values, in the [KnotGauss API](../specs/knot-classical-invariants.md#definitions).

# Guarantees
1. The result equals `p.terms` and describes the same stored Laurent polynomial. [proved/agent: by construction]
2. Adding, replacing or deleting dictionary entries in the result does not
   alter the input dictionary. The call itself leaves `p` unchanged. [proved/agent: by construction]

# Not guaranteed
Dictionary iteration order is not a term ordering. Stored zero coefficients
are not removed, and no normalization or invariant computation is performed.
The type of `p` alone does not establish an Alexander or Jones interpretation.

# Sources
cf:lib/topology.jl/src/knot_classical_invariants.jl, polynomial_terms — stored polynomial and PD data access; no independent knot classification.

# Examples
```jldoctest
julia> p = alexander_polynomial(pd_to_gauss(Int[]));

julia> polynomial_terms(p) == Dict(0=>big(1))
true
```

See also [`knot_polynomial_string`](@ref).

Role: access
"""
polynomial_terms(p::ExactKnotLaurentPolynomial) = copy(p.terms)

_lp_zero() = Dict{Int,BigInt}()
_lp_one() = Dict(0 => big(1))
_lp_monomial(exponent::Integer, coefficient::Integer=1) =
    coefficient == 0 ? _lp_zero() : Dict(Int(exponent) => BigInt(coefficient))

function _lp_add(a::AbstractDict{Int,BigInt}, b::AbstractDict{Int,BigInt})
    result = copy(a)
    for (exponent, coefficient) in b
        result[exponent] = get(result, exponent, big(0)) + coefficient
        result[exponent] == 0 && delete!(result, exponent)
    end
    result
end

function _lp_scale(a::AbstractDict{Int,BigInt}, coefficient::Integer)
    c = BigInt(coefficient)
    c == 0 && return _lp_zero()
    Dict(exponent => value * c for (exponent, value) in a)
end

function _lp_mul(a::AbstractDict{Int,BigInt}, b::AbstractDict{Int,BigInt})
    result = _lp_zero()
    for (ea, ca) in a, (eb, cb) in b
        exponent = ea + eb
        result[exponent] = get(result, exponent, big(0)) + ca * cb
        result[exponent] == 0 && delete!(result, exponent)
    end
    result
end

function _lp_pow(a::AbstractDict{Int,BigInt}, exponent::Integer)
    exponent >= 0 || throw(ArgumentError("Laurent powers must be nonnegative"))
    result = _lp_one()
    for _ in 1:Int(exponent)
        result = _lp_mul(result, a)
    end
    result
end

function _lp_det(matrix::Vector{Vector{Dict{Int,BigInt}}})
    n = length(matrix)
    n == 0 && return _lp_one()
    n == 1 && return copy(matrix[1][1])
    result = _lp_zero()
    for column in 1:n
        minor = Vector{Vector{Dict{Int,BigInt}}}()
        for row in 2:n
            push!(minor, [copy(matrix[row][j]) for j in 1:n if j != column])
        end
        term = _lp_mul(matrix[1][column], _lp_det(minor))
        isodd(column + 1) && (term = _lp_scale(term, -1))
        result = _lp_add(result, term)
    end
    result
end

function _normalized_alexander(terms::AbstractDict{Int,BigInt})
    isempty(terms) && throw(ArgumentError("Alexander determinant is zero"))
    least = minimum(keys(terms))
    shifted = Dict(exponent - least => coefficient for
                   (exponent, coefficient) in terms)
    leading = shifted[maximum(keys(shifted))]
    leading < 0 && (shifted = Dict(exponent => -coefficient
                                   for (exponent, coefficient) in shifted))
    ExactKnotLaurentPolynomial(:t, shifted)
end

function _knot_gauss_arcs(d::KnotGaussDiagram)
    n = length(d.orientation)
    n == 0 && return Int[], Int[], Int[]
    over = zeros(Int, n)
    under_in = zeros(Int, n)
    under_out = zeros(Int, n)
    current = 1
    for value in d.gauss
        crossing = abs(value)
        if value > 0
            over[crossing] = current
        else
            under_in[crossing] = current
            current = current == n ? 1 : current + 1
            under_out[crossing] = current
        end
    end
    all(!=(0), over) && all(!=(0), under_in) && all(!=(0), under_out) ||
        throw(ArgumentError("Gauss code must have one over and one under visit per crossing"))
    over, under_in, under_out
end

"""    alexander_polynomial(d::KnotGaussDiagram) -> ExactKnotLaurentPolynomial

Return the repository-normalized Alexander minor polynomial of a signed knot diagram.

# Arguments
`d` satisfies K1 ([Gauss conditions](../specs/knot-gauss.md#numbered-conditions))
by its checked constructor. P3 ([polynomial conditions](../specs/knot-classical-invariants.md#numbered-conditions))
is the caller's responsibility for a classical knot interpretation; this
function does not check classicality. A zero determinant raises `ArgumentError`
with message `Alexander determinant is zero`; invalid visit/sign data detected
by the calculation also raise `ArgumentError` describing that data.

# Returns
An `ExactKnotLaurentPolynomial` with variable `:t` and exact integer terms,
in [Alexander coordinates](../specs/knot-classical-invariants.md#alexander-coordinates).

# Guarantees
1. The polynomial is the normalized minor defined by the linked Alexander coordinates; the crossing-free circle returns one. [proved/agent: by construction]
2. The input diagram is unchanged. [proved/agent: by construction]

# Conventions
[Alexander normalization](../specs/knot-classical-invariants.md#alexander-coordinates).

# Not guaranteed
The general identification with the Alexander polynomial of a classical knot
is unverified at the repository source-comparison level; the specification
records the tested fixtures and missing source/convention bridge. No knot-type
classification follows from equality. This normalization is not the symmetric
surgery normalization: use the [surgery derivative translation](../specs/knot-classical-invariants.md#surgery-derivative)
before supplying a second derivative to the Casson–Walker surgery API.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, Alexander coordinates and Surgery derivative — local polynomial and normalization meanings.

# Examples
```jldoctest
julia> knot_polynomial_string(alexander_polynomial(pd_to_gauss(Int[])))
"1"
```

See also [`polynomial_terms`](@ref), [`jones_polynomial`](@ref).

Role: invariant
"""
function alexander_polynomial(d::KnotGaussDiagram)
    n = length(d.orientation)
    n == 0 && return ExactKnotLaurentPolynomial(:t, _lp_one())
    over, under_in, under_out = _knot_gauss_arcs(d)
    matrix = Vector{Vector{Dict{Int,BigInt}}}()
    for crossing in 1:n
        row = [_lp_zero() for _ in 1:n]
        s = d.orientation[crossing]
        s == 1 || s == -1 || throw(ArgumentError("crossing signs must be +/-1"))
        row[over[crossing]] = _lp_add(row[over[crossing]],
                                      _lp_add(_lp_monomial(0),
                                              _lp_monomial(s, -1)))
        row[under_in[crossing]] = _lp_add(row[under_in[crossing]],
                                          _lp_monomial(s))
        row[under_out[crossing]] = _lp_add(row[under_out[crossing]],
                                           _lp_monomial(0, -1))
        push!(matrix, row)
    end
    # The standard Alexander matrix has rank n-1 for a knot.  This minor is
    # exact and is defined for n=1 as the empty determinant 1.
    minor = Vector{Vector{Dict{Int,BigInt}}}()
    for row in 1:max(n - 1, 0)
        push!(minor, [copy(matrix[row][column]) for column in 1:max(n - 1, 0)])
    end
    _normalized_alexander(_lp_det(minor))
end

"""    KnotPD(rows, crossing_signs::AbstractVector{<:Integer}, component_count::Integer=1) -> KnotPD

Construct a record of oriented planar-diagram crossing coordinates.

# Arguments
`rows`, `crossing_signs` and `component_count` use
[PD coordinates](../specs/knot-classical-invariants.md#pd-coordinates).
P1–P2 ([conditions](../specs/knot-classical-invariants.md#numbered-conditions))
are the caller's responsibility for a knot diagram; construction is not full
validation. The general converting method checks row width, sign-count
agreement and positive component count, raising `ArgumentError` naming the
failure; integer conversion errors propagate. The more specific three-field
method with `Vector{NTuple{4,Int}}`, `Vector{Int}`, `Int` bypasses these checks.

# Returns
A `KnotPD` with `rows`, `crossing_signs` and `component_count` as defined in
the PD coordinate entry. Use `verify_pd` before relying on P1–P2 for raw data.

# Guarantees
1. Successful conversion retains the supplied coordinate values in the declared field types. With P1–P2 these are the crossing coordinates of the specified one-component traversal. [proved/agent: by construction]
2. The general converting method copies row/sign data; the typed field constructor shares its arrays. Neither changes an input array during construction. [proved/agent: by construction]

# Conventions
[PD slot and sign coordinates](../specs/knot-classical-invariants.md#pd-coordinates).

# Not guaranteed
A successful constructor does not certify planarity, component connectivity,
local-frame signs or a knot type. The general literature identification of
the coordinate convention remains unverified as stated in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, PD coordinates — repository crossing representation.

# Examples
```jldoctest
julia> pd = KnotPD(NTuple{4,Int}[], Int[], 1);

julia> verify_pd(pd)
true
```

See also [`verify_pd`](@ref), [`pd_to_gauss`](@ref).

Role: construct
"""
struct KnotPD
    rows::Vector{NTuple{4,Int}}
    crossing_signs::Vector{Int}
    component_count::Int
end

function _normalize_pd_rows(rows)
    normalized = NTuple{4,Int}[]
    for row in rows
        values = Int.(collect(row))
        length(values) == 4 || throw(ArgumentError("each PD row needs four edge labels"))
        push!(normalized, (values[1], values[2], values[3], values[4]))
    end
    normalized
end

function KnotPD(rows, crossing_signs::AbstractVector{<:Integer},
                component_count::Integer=1)
    normalized = _normalize_pd_rows(rows)
    signs = Int.(collect(crossing_signs))
    length(normalized) == length(signs) ||
        throw(ArgumentError("PD sign count must equal row count"))
    component_count > 0 || throw(ArgumentError("component_count must be positive"))
    KnotPD(normalized, signs, Int(component_count))
end

Base.length(pd::KnotPD) = length(pd.rows)

"""
    pd_matrix(pd::KnotPD) -> Matrix{Int}

Return the crossing rows of a planar-diagram record as an integer matrix.

# Arguments
`pd` is a `KnotPD`. No edge-incidence, crossing-sign, component or planarity
validation is performed by this accessor.

# Returns
A newly allocated `Matrix{Int}` with one row per entry of `pd.rows` and
four columns, preserving the record's row and slot orders. This belongs to
the [KnotGauss API](../specs/knot-classical-invariants.md#definitions).

# Guarantees
1. Matrix row `i` contains the four entries of `pd.rows[i]` in order. [proved/agent: by construction]
2. The result carries the same PD incidence rows as the input. Changes to
   the matrix do not change `pd.rows`, and the input is unchanged by the call. [proved/agent: by construction]
3. An empty row list returns a matrix with zero rows and four columns. [proved/agent: by construction]

# Not guaranteed
The matrix omits `crossing_signs` and `component_count` as separate fields.
This operation does not certify a planar embedding or a knot type, choose
a PD convention, or relabel the edges.

# Sources
cf:lib/topology.jl/src/knot_classical_invariants.jl, pd_matrix — stored polynomial and PD data access; no independent knot classification.

# Examples
```jldoctest
julia> pd = pd_from_gauss(pd_to_gauss(Int[]));

julia> size(pd_matrix(pd), 1)
0
```

See also [`pd_from_matrix`](@ref), [`verify_pd`](@ref).

Role: access
"""
function pd_matrix(pd::KnotPD)
    matrix = zeros(Int, length(pd.rows), 4)
    for (row_index, row) in enumerate(pd.rows)
        matrix[row_index, :] = collect(row)
    end
    matrix
end

function _pd_occurrences(pd::KnotPD)
    occurrences = Dict{Int,Vector{Tuple{Int,Int}}}()
    for (row_index, row) in enumerate(pd.rows), slot in 1:4
        edge = row[slot]
        edge > 0 || throw(ArgumentError("PD edge labels must be positive"))
        push!(get!(occurrences, edge, Tuple{Int,Int}[]), (row_index, slot))
    end
    all(length(value) == 2 for value in values(occurrences)) ||
        throw(ArgumentError("every PD edge must occur exactly twice"))
    occurrences
end

function _pd_other_halfedge(pd::KnotPD, occurrences, halfedge::Tuple{Int,Int})
    row, slot = halfedge
    edge = pd.rows[row][slot]
    pair = occurrences[edge]
    pair[1] == halfedge ? pair[2] : pair[1]
end

function _pd_transition(pd::KnotPD, occurrences, halfedge::Tuple{Int,Int})
    row, slot = halfedge
    opposite = (row, mod1(slot + 2, 4))
    _pd_other_halfedge(pd, occurrences, opposite)
end

function _pd_cycles(pd::KnotPD, occurrences)
    seen = Set{Tuple{Int,Int}}()
    cycles = Vector{Vector{Tuple{Int,Int}}}()
    for row in 1:length(pd.rows), slot in 1:4
        start = (row, slot)
        start in seen && continue
        cycle = Tuple{Int,Int}[]
        current = start
        while !(current in seen)
            push!(seen, current)
            push!(cycle, current)
            current = _pd_transition(pd, occurrences, current)
        end
        current == start || throw(ArgumentError("PD transition is not cyclic"))
        push!(cycles, cycle)
    end
    cycles
end

const _PD_SLOT_VECTOR = ((1, 0), (0, 1), (-1, 0), (0, -1))

function _pd_code_and_signs(pd::KnotPD, occurrences, cycle)
    n = length(pd.rows)
    length(cycle) == 2n || throw(ArgumentError(
        "selected PD orientation does not visit one component"))
    code = Int[]
    directions = [Dict{Symbol,Tuple{Int,Int}}() for _ in 1:n]
    for halfedge in cycle
        row, slot = halfedge
        opposite = (row, mod1(slot + 2, 4))
        endpoint = slot in (1, 3) ? :under : :over
        haskey(directions[row], endpoint) && throw(ArgumentError(
            "selected PD orientation visits a crossing twice on one strand"))
        directions[row][endpoint] = (slot, opposite[2])
        push!(code, endpoint === :over ? row : -row)
    end
    signs = zeros(Int, n)
    for row in 1:n
        haskey(directions[row], :under) && haskey(directions[row], :over) ||
            throw(ArgumentError("selected PD orientation misses a crossing"))
        under = directions[row][:under]
        over = directions[row][:over]
        uv = (_PD_SLOT_VECTOR[under[2]][1] - _PD_SLOT_VECTOR[under[1]][1],
              _PD_SLOT_VECTOR[under[2]][2] - _PD_SLOT_VECTOR[under[1]][2])
        ov = (_PD_SLOT_VECTOR[over[2]][1] - _PD_SLOT_VECTOR[over[1]][1],
              _PD_SLOT_VECTOR[over[2]][2] - _PD_SLOT_VECTOR[over[1]][2])
        determinant = ov[1] * uv[2] - ov[2] * uv[1]
        determinant == 0 && throw(ArgumentError("degenerate PD crossing frame"))
        signs[row] = determinant > 0 ? 1 : -1
    end
    code, signs
end

"""
    verify_pd(pd::KnotPD; require_one_component::Bool=true, check_signs::Bool=true) -> Bool

Check the combinatorial incidence, traversal and optional local signs of a PD record.

# Arguments
`pd` is a raw [PD record](../specs/knot-classical-invariants.md#pd-coordinates).
P1–P2 ([conditions](../specs/knot-classical-invariants.md#numbered-conditions))
are checked in the supported one-component domain. Invalid sign count or
values, labels, occurrences, traversal cycles or component count raise
`ArgumentError` identifying the failure. `check_signs=false` omits comparison
with reconstructed local signs, but still checks that stored signs are `±1`.
`require_one_component=false` omits the explicit one-component guard; other
traversal-length checks remain active.

# Returns
`true` when the enabled checks succeed. Invalid supported data raise an
exception rather than returning `false`.

# Guarantees
1. Successful return certifies the checked combinatorial PD conditions; with default options, the reconstructed local signs agree with the stored signs and the component count is one. [proved/agent: by construction]
2. The input record and its arrays are unchanged. [proved/agent: by construction]

# Conventions
[PD frame](../specs/knot-classical-invariants.md#pd-coordinates).

# Not guaranteed
This does not certify a plane embedding, P3 or a knot type. Disabling the
one-component guard is not a general link-validation API: traversal checks
still impose the current knot-oriented restrictions. Disabling sign comparison
does not make the stored signs the reconstructed signs.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, PD coordinates and P1–P3 — checked representation versus classical interpretation.

# Examples
```jldoctest
julia> verify_pd(pd_from_gauss(pd_to_gauss(Int[])))
true
```

See also [`pd_to_gauss`](@ref).

Role: condition
"""
function verify_pd(pd::KnotPD; require_one_component::Bool=true,
                   check_signs::Bool=true)
    n = length(pd.rows)
    length(pd.crossing_signs) == n || throw(ArgumentError("PD sign count mismatch"))
    all(sign in (-1, 1) for sign in pd.crossing_signs) ||
        throw(ArgumentError("PD crossing signs must be +/-1"))
    occurrences = _pd_occurrences(pd)
    cycles = _pd_cycles(pd, occurrences)
    if n == 0
        pd.component_count == 1 || throw(ArgumentError(
            "the empty PD represents exactly one unknot component"))
        isempty(cycles) || throw(ArgumentError("empty PD has half-edge cycles"))
        return true
    end
    length(cycles) == 2 * pd.component_count || throw(ArgumentError(
        "PD rotation has $(length(cycles)) oriented cycles, expected $(2 * pd.component_count)"))
    all(length(cycle) == 2n for cycle in cycles) || throw(ArgumentError(
        "PD component cycles have inconsistent lengths"))
    if require_one_component
        pd.component_count == 1 || throw(ArgumentError("links are outside this knot API"))
    end
    code, signs = _pd_code_and_signs(pd, occurrences, cycles[1])
    check_signs && signs != pd.crossing_signs && throw(ArgumentError(
        "PD local frame signs disagree with stored crossing signs"))
    length(code) == 2n || throw(ArgumentError("PD Gauss traversal has wrong length"))
    true
end

"""    pd_from_gauss(d::KnotGaussDiagram) -> KnotPD

Convert signed Gauss data to the repository's oriented PD crossing coordinates.

# Arguments
`d` satisfies K1 ([Gauss conditions](../specs/knot-gauss.md#numbered-conditions)).
The result is checked by `verify_pd` under P1–P2
([PD conditions](../specs/knot-classical-invariants.md#numbered-conditions));
failed incidence, traversal or local-sign checks raise `ArgumentError` naming
the failure. A planar-knot interpretation additionally requires caller-supplied
P3; no independent classicality check is performed here.

# Returns
A `KnotPD` with converted rows and crossing signs, with component count one,
in [PD coordinates](../specs/knot-classical-invariants.md#pd-coordinates).

# Guarantees
1. Successful return satisfies `verify_pd(result)` and represents the local crossing/arc data of the input in the specified PD frame. [proved/agent: by construction]
2. The output owns new row and sign arrays; the input is unchanged. [proved/agent: by construction]

# Conventions
[PD coordinates](../specs/knot-classical-invariants.md#pd-coordinates).

# Not guaranteed
A checked combinatorial conversion does not supply a plane embedding or
identify a knot type. The general source comparison for the PD frame remains
unverified as recorded in the specification. Conversion back from raw rows
need not preserve the original traversal root or direction.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, PD coordinates — local Gauss/PD dictionary.

# Examples
```jldoctest
julia> d = pd_to_gauss([[1,5,2,4], [3,1,4,6], [5,3,6,2]]);

julia> verify_pd(pd_from_gauss(d))
true
```

See also [`pd_matrix`](@ref), [`pd_to_gauss`](@ref).

Role: convert
"""
function pd_from_gauss(d::KnotGaussDiagram)
    n = length(d.orientation)
    n == 0 && return KnotPD(NTuple{4,Int}[], Int[], 1)
    length(d.gauss) == 2n || throw(ArgumentError("Gauss code length mismatch"))
    visits = Dict{Int,Tuple{Int,Int}}()
    for (position, value) in enumerate(d.gauss)
        visits[value] = (position, mod1(position + 1, 2n))
    end
    rows = NTuple{4,Int}[]
    for crossing in 1:n
        under = visits[-crossing]
        over = visits[crossing]
        row = if d.orientation[crossing] == -1
            (under[1], over[1], under[2], over[2])
        elseif d.orientation[crossing] == 1
            (under[1], over[2], under[2], over[1])
        else
            throw(ArgumentError("crossing signs must be +/-1"))
        end
        push!(rows, row)
    end
    pd = KnotPD(rows, collect(d.orientation), 1)
    verify_pd(pd)
    pd
end

"""
    pd_to_gauss(rows) -> KnotGaussDiagram

Convert raw PD crossing rows to a signed one-component Gauss diagram.

# Arguments
`rows` uses the raw-row interface in
[PD coordinates](../specs/knot-classical-invariants.md#pd-coordinates): pass
`eachrow(M)` for a Julia matrix, or `.rows` for a `KnotPD`. Row entries must
convert to `Int`; conversion errors propagate. P1–P2
([conditions](../specs/knot-classical-invariants.md#numbered-conditions)) are
checked, with `ArgumentError` describing invalid row width, labels, incidence,
traversal or component structure. P3 for a classical-knot interpretation is
the caller's responsibility, not independently checked here.

# Returns
A `KnotGaussDiagram` with visit and orientation tuples reconstructed from the PD traversal.

# Guarantees
1. On successful return, the converted data use the specified PD/Gauss local crossing correspondence; empty input represents the crossing-free circle. [proved/agent: by construction]
2. The result satisfies K1. [proved/agent: by construction]
3. The input rows are unchanged and output data are newly allocated. [proved/agent: by construction]

# Conventions
[PD frame and traversal coordinates](../specs/knot-classical-invariants.md#pd-coordinates).

# Not guaranteed
No plane embedding, classicality proof or knot-type identification is supplied.
The general literature convention comparison remains unverified as stated in
the specification. Input row order and local traversal choice can affect the
returned representative; no canonical-root or traversal-reversal quotient
is promised.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, PD coordinates — one-component coordinate conversion and source boundary.

# Examples
```jldoctest
julia> crossing_count(pd_to_gauss(Int[]))
0
```

Role: construct
"""
function pd_to_gauss(rows)
    pd = pd_from_matrix(rows)
    n = length(pd.rows)
    n == 0 && return KnotGaussDiagram(Int[], Int[])
    occurrences = _pd_occurrences(pd)
    cycles = _pd_cycles(pd, occurrences)
    code, signs = _pd_code_and_signs(pd, occurrences, cycles[1])
    KnotGaussDiagram(code, signs)
end

"""
    pd_from_matrix(rows) -> KnotPD

Convert raw crossing rows to a sign-checked one-component PD record.

# Arguments
`rows` uses the raw-row interface in
[PD coordinates](../specs/knot-classical-invariants.md#pd-coordinates): pass
`eachrow(M)` for a Julia matrix, or `.rows` for a `KnotPD`. Row entries must
convert to `Int`; conversion errors propagate. P1–P2
([conditions](../specs/knot-classical-invariants.md#numbered-conditions)) are
checked, with `ArgumentError` describing invalid row width, labels, incidence,
traversal or component structure. P3 for a classical-knot interpretation is
the caller's responsibility, not independently checked here.

# Returns
A `KnotPD` with normalized rows, crossing signs reconstructed from the local frame and component count one.

# Guarantees
1. On successful return, the converted data use the specified PD/Gauss local crossing correspondence; empty input represents the crossing-free circle. [proved/agent: by construction]
2. `verify_pd(result)` is true. [proved/agent: by construction]
3. The input rows are unchanged and output data are newly allocated. [proved/agent: by construction]

# Conventions
[PD frame and traversal coordinates](../specs/knot-classical-invariants.md#pd-coordinates).

# Not guaranteed
No plane embedding, classicality proof or knot-type identification is supplied.
The general literature convention comparison remains unverified as stated in
the specification. Input row order and local traversal choice can affect the
returned representative; no canonical-root or traversal-reversal quotient
is promised.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, PD coordinates — one-component coordinate conversion and source boundary.

# Examples
```jldoctest
julia> verify_pd(pd_from_matrix(Int[]))
true
```

Role: construct
"""
function pd_from_matrix(rows)
    normalized = _normalize_pd_rows(rows)
    n = length(normalized)
    n == 0 && return KnotPD(NTuple{4,Int}[], Int[], 1)
    provisional = KnotPD(normalized, fill(1, n), 1)
    verify_pd(provisional; check_signs=false)
    occurrences = _pd_occurrences(provisional)
    cycles = _pd_cycles(provisional, occurrences)
    _, signs = _pd_code_and_signs(provisional, occurrences, cycles[1])
    pd = KnotPD(normalized, signs, 1)
    verify_pd(pd)
    pd
end

function _union_find(n::Integer)
    parent = collect(1:Int(n))
    rank = zeros(Int, Int(n))
    parent, rank
end

function _uf_find!(parent::Vector{Int}, node::Int)
    root = node
    while parent[root] != root
        root = parent[root]
    end
    while parent[node] != node
        next = parent[node]
        parent[node] = root
        node = next
    end
    root
end

function _uf_union!(parent::Vector{Int}, rank::Vector{Int}, a::Int, b::Int)
    ra = _uf_find!(parent, a)
    rb = _uf_find!(parent, b)
    ra == rb && return
    if rank[ra] < rank[rb]
        parent[ra] = rb
    elseif rank[ra] > rank[rb]
        parent[rb] = ra
    else
        parent[rb] = ra
        rank[ra] += 1
    end
end

function _state_loops(pd::KnotPD, state::Integer)
    n = length(pd.rows)
    parent, rank = _union_find(4n)
    occurrences = _pd_occurrences(pd)
    for edge_occurrences in values(occurrences)
        a = (edge_occurrences[1][1] - 1) * 4 + edge_occurrences[1][2]
        b = (edge_occurrences[2][1] - 1) * 4 + edge_occurrences[2][2]
        _uf_union!(parent, rank, a, b)
    end
    for row in 1:n
        is_b = ((Int(state) >> (row - 1)) & 1) == 1
        pairs = is_b ? ((2, 3), (4, 1)) : ((1, 2), (3, 4))
        for (a_slot, b_slot) in pairs
            a = (row - 1) * 4 + a_slot
            b = (row - 1) * 4 + b_slot
            _uf_union!(parent, rank, a, b)
        end
    end
    length(Set(_uf_find!.(Ref(parent), 1:4n)))
end

function _bracket_polynomial(pd::KnotPD)
    n = length(pd.rows)
    n < 62 || throw(ArgumentError("state sum is bounded to fewer than 62 crossings"))
    n == 0 && return _lp_one()
    delta = Dict(2 => big(-1), -2 => big(-1))
    total = _lp_zero()
    for state in 0:(2^n - 1)
        b_count = count_ones(UInt(state))
        exponent = n - 2 * b_count
        loops = _state_loops(pd, state)
        contribution = _lp_shift(_lp_pow(delta, loops - 1), exponent)
        total = _lp_add(total, contribution)
    end
    total
end

function _lp_shift(terms::AbstractDict{Int,BigInt}, shift::Integer)
    Dict(exponent + Int(shift) => coefficient for
         (exponent, coefficient) in terms)
end

function _jones_t_terms(normalized_a::AbstractDict{Int,BigInt})
    isempty(normalized_a) && throw(ArgumentError("normalized Jones polynomial is zero"))
    residues = Set(mod(exponent, 4) for exponent in keys(normalized_a))
    residues == Set([0]) || throw(ArgumentError(
        "Jones state sum has a nonzero A exponent residue: $residues"))
    Dict(-div(exponent, 4) => coefficient
         for (exponent, coefficient) in normalized_a)
end

"""    jones_polynomial(d::KnotGaussDiagram) -> ExactKnotLaurentPolynomial

Return the repository-normalized Jones bracket polynomial of a signed knot diagram.

# Arguments
`d` satisfies K1 ([Gauss conditions](../specs/knot-gauss.md#numbered-conditions)).
P3 ([polynomial conditions](../specs/knot-classical-invariants.md#numbered-conditions))
is the caller's responsibility for a classical-knot interpretation. The
PD conversion checks its combinatorial conditions; invalid data raise
`ArgumentError` describing the failed check. P4 is checked: too many crossings
or residual nonintegral powers of the Jones variable raise `ArgumentError`.

# Returns
An `ExactKnotLaurentPolynomial` with variable `:t` and integer terms in the
[Jones coordinates](../specs/knot-classical-invariants.md#jones-coordinates).

# Guarantees
1. The result is the exact state-sum polynomial with the writhe and variable normalization defined by those coordinates; the crossing-free circle returns one. [proved/agent: by construction]
2. The input diagram is unchanged. [proved/agent: by construction]

# Conventions
[Jones smoothing and variable normalization](../specs/knot-classical-invariants.md#jones-coordinates).

# Not guaranteed
The general identification with the Jones polynomial of a classical knot is
unverified at the repository source-comparison level; the specification states
the tested fixtures and missing convention bridge. The crossing guard is not
a practical runtime bound: the state sum grows exponentially. Polynomial
equality does not classify knots, and no classicality certificate is returned.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, Jones coordinates — exact local bracket definition and evidence boundary.

# Examples
```jldoctest
julia> knot_polynomial_string(jones_polynomial(pd_to_gauss(Int[])))
"1"
```

See also [`invert_knot_variable`](@ref), [`alexander_polynomial`](@ref).

Role: invariant
"""
function jones_polynomial(d::KnotGaussDiagram)
    pd = pd_from_gauss(d)
    verify_pd(pd)
    bracket = _bracket_polynomial(pd)
    writhe = sum(pd.crossing_signs)
    normalized = _lp_shift(bracket, -3 * writhe)
    if isodd(abs(writhe))
        normalized = _lp_scale(normalized, -1)
    end
    ExactKnotLaurentPolynomial(:t, _jones_t_terms(normalized))
end

"""    invert_knot_variable(p::ExactKnotLaurentPolynomial) -> ExactKnotLaurentPolynomial

Return the Laurent polynomial obtained by inverting its formal variable.

# Arguments
`p` is an [exact polynomial record](../specs/knot-classical-invariants.md#exact-polynomial-records).
Its exponent negations must fit `Int` (P5 in the
[conditions](../specs/knot-classical-invariants.md#numbered-conditions)); this
is the caller's responsibility and overflow is not detected. Invalid raw
record data are outside the normalized-polynomial interpretation.

# Returns
An `ExactKnotLaurentPolynomial` with the same variable symbol and a new term
dictionary representing `p(t^-1)`.

# Guarantees
1. Under P5, the coefficient at exponent `k` in the result is the input coefficient at `-k`; no coefficient conjugation is performed. [proved/agent: by construction]
2. The input record and its dictionary are unchanged; the result owns a new dictionary. [proved/agent: by construction]

# Not guaranteed
No knot is constructed or checked. The mirror interpretation for Jones data
has the fixture scope and missing general source comparison stated in the
polynomial specification. The operation does not impose Alexander normalization.

# Sources
cf:lib/topology.jl/docs/specs/knot-classical-invariants.md, Exact polynomial records and Coordinate choices — variable inversion and its finite mirror checks.

# Examples
```jldoctest
julia> d = pd_to_gauss([[1,5,2,4], [3,1,4,6], [5,3,6,2]]);

julia> jones_polynomial(mirror(d)) == invert_knot_variable(jones_polynomial(d))
true
```

Role: move
"""
function invert_knot_variable(p::ExactKnotLaurentPolynomial)
    ExactKnotLaurentPolynomial(p.variable,
                               Dict(-exponent => coefficient
                                    for (exponent, coefficient) in p.terms))
end

export ExactKnotLaurentPolynomial, KnotPD, pd_matrix, pd_from_gauss, pd_to_gauss,
       pd_from_matrix, verify_pd, polynomial_terms, knot_polynomial_string,
       alexander_polynomial, jones_polynomial, invert_knot_variable
