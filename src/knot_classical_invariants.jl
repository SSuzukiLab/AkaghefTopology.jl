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

"""A sparse exact Laurent polynomial with integer coefficients."""
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

"""Compute the normalized Alexander polynomial by a Wirtinger/Fox minor."""
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

"""PD rows and their independently reconstructed oriented crossing signs."""
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

"""Validate a PD/rotation representation; links are rejected by default."""
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

"""Convert a signed Gauss diagram to the explicit local PD convention above."""
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

"""Reconstruct one signed one-component Gauss diagram from PD notation.

The selected orientation starts at row 1, slot 1.  Its reversal is equivalent
and is identified later by the Gauss canonicalizer; no invariant value is used
to choose the signed decoration.
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

"""Normalize and independently sign-check a raw PD matrix."""
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

"""Compute the normalized Jones polynomial by an exact bracket state sum."""
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

"""Substitute `t -> t^-1`, used for the mirror cross-check."""
function invert_knot_variable(p::ExactKnotLaurentPolynomial)
    ExactKnotLaurentPolynomial(p.variable,
                               Dict(-exponent => coefficient
                                    for (exponent, coefficient) in p.terms))
end

export ExactKnotLaurentPolynomial, KnotPD, pd_matrix, pd_from_gauss, pd_to_gauss,
       pd_from_matrix, verify_pd, polynomial_terms, knot_polynomial_string,
       alexander_polynomial, jones_polynomial, invert_knot_variable
