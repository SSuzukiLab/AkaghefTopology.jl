"""Validated signed Gauss representations for one-component diagrams.

The code is a cyclic vector of signed crossing labels.  Each absolute label
occurs exactly twice, with opposite visit signs.  `orientation` stores one
independent crossing sign (`+1` or `-1`) for every label.  This combinatorial
constructor does *not* certify realizability by a classical planar diagram;
callers needing that property must run a separately justified predicate.
Neither reversing the traversal orientation nor mirroring is identified by
canonicalization.
"""

using SHA

const KNOT_GAUSS_CANON_VERSION = 1

"""This module enumerates signed one-component Gauss representations."""
const KNOT_GAUSS_CENSUS_KIND = :signed_gauss_representation

"""No classical-planarity/realizability oracle is bundled in this lane."""
classical_knot_census_available() = false

struct KnotGaussDiagram
    gauss::Tuple{Vararg{Int}}
    orientation::Tuple{Vararg{Int}}

    function KnotGaussDiagram(gauss::AbstractVector{<:Integer},
                              orientation::AbstractVector{<:Integer})
        code = Int.(collect(gauss))
        signs = Int.(collect(orientation))
        n = length(signs)
        iseven(length(code)) || throw(ArgumentError("Gauss code must have even length"))
        all(x -> x != 0, code) || throw(ArgumentError("Gauss code labels cannot be zero"))
        all(x -> x == 1 || x == -1, signs) ||
            throw(ArgumentError("crossing orientation signs must be +1 or -1"))
        n == length(code) ÷ 2 ||
            throw(ArgumentError("orientation must contain one sign per crossing"))
        for label in 1:n
            occurrences = findall(x -> abs(x) == label, code)
            length(occurrences) == 2 ||
                throw(ArgumentError("crossing $label must occur exactly twice"))
            code[occurrences[1]] * code[occurrences[2]] < 0 ||
                throw(ArgumentError("crossing $label must occur once positive and once negative"))
        end
        new(Tuple(code), Tuple(signs))
    end
end

KnotGaussDiagram(; gauss, orientation) = KnotGaussDiagram(gauss, orientation)
KnotGaussDiagram(code::Tuple, orientation) =
    KnotGaussDiagram(collect(code), collect(orientation))

"""Fail explicitly rather than treating a representation as a classical knot.

A separately justified realizability algorithm must be supplied before a
classical-knot census can be requested.
"""
function require_classical_realizability(::KnotGaussDiagram)
    throw(ArgumentError("classical realizability is not implemented; this is a signed Gauss representation only"))
end

Base.copy(d::KnotGaussDiagram) = KnotGaussDiagram(d.gauss, d.orientation)
Base.:(==)(a::KnotGaussDiagram, b::KnotGaussDiagram) =
    a.gauss == b.gauss && a.orientation == b.orientation
Base.hash(d::KnotGaussDiagram, h::UInt) = hash((d.gauss, d.orientation), h)

function Base.getproperty(d::KnotGaussDiagram, name::Symbol)
    name === :code && return getfield(d, :gauss)
    name === :crossing_sign && return getfield(d, :orientation)
    name === :signs && return getfield(d, :orientation)
    getfield(d, name)
end

crossing_count(d::KnotGaussDiagram) = length(d.orientation)

"""Return a mutable copy of the visit code for APIs that require vectors."""
gauss_vector(d::KnotGaussDiagram) = collect(d.gauss)

"""Return a mutable copy of the crossing signs for APIs that require vectors."""
orientation_vector(d::KnotGaussDiagram) = collect(d.orientation)

"""Rotate the rooted representative by `amount` positions."""
function cyclic_shift(d::KnotGaussDiagram, amount::Integer)
    n = length(d.gauss)
    n == 0 && return copy(d)
    KnotGaussDiagram(circshift(d.gauss, mod(Int(amount), n)), d.orientation)
end

cyclic_shift(d::KnotGaussDiagram) = cyclic_shift(d, 1)

"""Relabel crossings by an old-label => new-label permutation."""
function reindex(d::KnotGaussDiagram, permutation::AbstractVector{<:Integer})
    n = crossing_count(d)
    p = Int.(collect(permutation))
    length(p) == n || throw(DimensionMismatch("reindex permutation has the wrong length"))
    sort(p) == collect(1:n) || throw(ArgumentError("reindex must be a permutation of 1:$n"))
    code = [sign(x) * p[abs(x)] for x in d.gauss]
    orientation = fill(0, n)
    for old in 1:n
        orientation[p[old]] = d.orientation[old]
    end
    KnotGaussDiagram(code, orientation)
end

function reindex(d::KnotGaussDiagram, permutation::AbstractDict)
    values = Int[]
    for i in 1:crossing_count(d)
        haskey(permutation, i) || throw(ArgumentError("reindex map is missing label $i"))
        push!(values, permutation[i])
    end
    reindex(d, values)
end

reindex(d::KnotGaussDiagram, f::Function) =
    reindex(d, [f(i) for i in 1:crossing_count(d)])

function _rotation_signature(d::KnotGaussDiagram, start::Int)
    length(d.gauss) == 0 && return Tuple{Int,Int,Int}[]
    l = length(d.gauss)
    positions = Dict{Int,Vector{Int}}()
    for (position, label) in enumerate(d.gauss)
        push!(get!(positions, abs(label), Int[]), position)
    end
    result = Tuple{Int,Int,Int}[]
    for j in 0:l-1
        position = mod1(start + j, l)
        pair = positions[abs(d.gauss[position])]
        other = pair[1] == position ? pair[2] : pair[1]
        offset = mod(other - position, l)
        push!(result, (sign(d.gauss[position]), offset, d.orientation[abs(d.gauss[position])]))
    end
    result
end

function _canonical_from_signature(signature::Vector{<:NTuple{3,Int}})
    labels = Dict{Int,Int}()
    code = Int[]
    orientation = Int[]
    next_label = 1
    for position in eachindex(signature)
        visit_sign, offset, crossing_sign = signature[position]
        partner = mod1(position + offset, length(signature))
        pair_key = min(position, partner)
        if !haskey(labels, pair_key)
            labels[pair_key] = next_label
            push!(orientation, crossing_sign)
            next_label += 1
        end
        label = labels[pair_key]
        push!(code, visit_sign * label)
    end
    KnotGaussDiagram(code, orientation)
end

function _diagram_key(d::KnotGaussDiagram)
    Tuple(vcat(d.gauss, [0], d.orientation))
end

"""Oracle canonicalization: minimum of all rotations after first-occurrence relabeling."""
function canonicalize_slow(d::KnotGaussDiagram)
    l = length(d.gauss)
    l == 0 && return copy(d)
    best = nothing
    for start in 1:l
        code = Int[]
        orientation = Int[]
        labels = Dict{Int,Int}()
        next_label = 1
        for j in 0:l-1
            value = d.gauss[mod1(start + j, l)]
            old = abs(value)
            if !haskey(labels, old)
                labels[old] = next_label
                push!(orientation, d.orientation[old])
                next_label += 1
            end
            push!(code, sign(value) * labels[old])
        end
        candidate = KnotGaussDiagram(code, orientation)
        (best === nothing || _diagram_key(candidate) < _diagram_key(best)) && (best = candidate)
    end
    best
end

"""Production canonicalization from label-free partner-offset signatures."""
function canonicalize(d::KnotGaussDiagram)
    l = length(d.gauss)
    l == 0 && return copy(d)
    best = nothing
    for start in 1:l
        candidate = _canonical_from_signature(_rotation_signature(d, start))
        (best === nothing || _diagram_key(candidate) < _diagram_key(best)) && (best = candidate)
    end
    best
end

canonicalize_fast(d::KnotGaussDiagram) = canonicalize(d)

function canonical_serialization(d::KnotGaussDiagram;
                                canon_version::Integer=KNOT_GAUSS_CANON_VERSION)
    canon_version == KNOT_GAUSS_CANON_VERSION ||
        throw(ArgumentError("unsupported KnotGauss canonicalization version $canon_version"))
    c = canonicalize(d)
    "KnotGaussDiagram;canon=v$(canon_version);n=$(crossing_count(c));" *
        "gauss=$(join(c.gauss, ','));orientation=$(join(c.orientation, ','))"
end

canonical_serialize(d::KnotGaussDiagram; kwargs...) = canonical_serialization(d; kwargs...)
serialize_canonical(d::KnotGaussDiagram; kwargs...) = canonical_serialization(d; kwargs...)

function canonical_hash(d::KnotGaussDiagram;
                       canon_version::Integer=KNOT_GAUSS_CANON_VERSION)
    bytes2hex(SHA.sha256(codeunits(canonical_serialization(d; canon_version=canon_version))))
end

function rotation_automorphism_count(d::KnotGaussDiagram)
    c = canonicalize(d)
    count(start -> _canonical_from_signature(_rotation_signature(d, start)) == c,
          1:length(d.gauss))
end

symmetry_rotation_count(d::KnotGaussDiagram) = rotation_automorphism_count(d)
automorphism_rotation_count(d::KnotGaussDiagram) = rotation_automorphism_count(d)
symmetry_count(d::KnotGaussDiagram) = rotation_automorphism_count(d)

"""Reverse the orientation of the parametrized knot (not a quotient operation)."""
orientation_reverse(d::KnotGaussDiagram) =
    KnotGaussDiagram(reverse(d.gauss), d.orientation)

reverse_orientation(d::KnotGaussDiagram) = orientation_reverse(d)

"""Mirror the diagram under the stored crossing-sign convention."""
mirror(d::KnotGaussDiagram) =
    KnotGaussDiagram([-x for x in d.gauss], [-x for x in d.orientation])
mirror_diagram(d::KnotGaussDiagram) = mirror(d)

function _rooted_signed_count(n::Int)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    n <= 6 || throw(ArgumentError("deterministic enumeration is bounded by n=6"))
    # Pairings of the 2n rooted positions contribute (2n)!/(2^n n!),
    # while choosing which visit is positive contributes 2^n.
    factorial(big(2n)) ÷ factorial(big(n))
end

rooted_normalized_signed_gc_count(n::Integer) = Int(_rooted_signed_count(Int(n)))
rooted_normalized_signed_gauss_count(n::Integer) = rooted_normalized_signed_gc_count(n)

"""A deterministic lazy stream of rooted signed representations.

This is a representation census, not a census of classical knots: entries are
not filtered for planar realizability, and `length` is available without
materializing the stream (including at the n=6 boundary).
"""
struct RootedSignedGaussEnumeration
    n::Int
end

Base.IteratorSize(::Type{RootedSignedGaussEnumeration}) = Base.HasLength()
Base.eltype(::Type{RootedSignedGaussEnumeration}) = Vector{Int}
Base.length(it::RootedSignedGaussEnumeration) = rooted_normalized_signed_gc_count(it.n)

function _unrank_rooted_signed(n::Int, rank::Int)
    remaining = 2n
    counts = zeros(Int, n)
    first_sign = zeros(Int, n)
    next_label = 1
    result = Int[]
    completion(r, k, open) = _rooted_completion_count(n, r, k, open)
    for position in 1:2n
        open_labels = [label for label in 1:n if counts[label] == 1]
        chosen = false
        for label in open_labels
            branch = Int(completion(remaining - 1, next_label - 1, length(open_labels) - 1))
            if rank >= branch
                rank -= branch
            else
                counts[label] = 2
                push!(result, -first_sign[label] * label)
                chosen = true
                break
            end
        end
        chosen && (remaining -= 1; continue)
        next_label <= n || error("internal rooted Gauss-code unranking error")
        branch_each_sign = Int(completion(remaining - 1, next_label, length(open_labels) + 1))
        rank < 2 * branch_each_sign || error("internal rooted Gauss-code rank error")
        first_sign[next_label] = rank < branch_each_sign ? 1 : -1
        rank < branch_each_sign || (rank -= branch_each_sign)
        counts[next_label] = 1
        push!(result, first_sign[next_label] * next_label)
        next_label += 1
        remaining -= 1
    end
    result
end

const _ROOTED_COMPLETION_CACHE = Dict{NTuple{4,Int},BigInt}()

function _rooted_completion_count(n::Int, r::Int, k::Int, open::Int)
    key = (n, r, k, open)
    haskey(_ROOTED_COMPLETION_CACHE, key) && return _ROOTED_COMPLETION_CACHE[key]
    value = if r == 0
        open == 0 ? big(1) : big(0)
    else
        result = big(0)
        open > 0 && (result += open * _rooted_completion_count(n, r - 1, k, open - 1))
        k < n && (result += 2 * _rooted_completion_count(n, r - 1, k + 1, open + 1))
        result
    end
    _ROOTED_COMPLETION_CACHE[key] = value
end

function Base.iterate(it::RootedSignedGaussEnumeration, state::Int=0)
    state >= length(it) && return nothing
    (_unrank_rooted_signed(it.n, state), state + 1)
end

enumerate_rooted_normalized_signed_gc(n::Integer) = RootedSignedGaussEnumeration(Int(n))
enumerate_rooted_normalized_signed_gauss(n::Integer) = enumerate_rooted_normalized_signed_gc(n)

struct KnotGaussEnumeration
    n::Int
    orientations::Bool
end

Base.IteratorSize(::Type{KnotGaussEnumeration}) = Base.HasLength()
Base.eltype(::Type{KnotGaussEnumeration}) = KnotGaussDiagram
function Base.length(it::KnotGaussEnumeration)
    base = rooted_normalized_signed_gc_count(it.n)
    it.orientations ? base * 2^it.n : base
end

function Base.iterate(it::KnotGaussEnumeration, state::Int=0)
    state >= length(it) && return nothing
    assignments = it.orientations ? 2^it.n : 1
    code_index = state ÷ assignments
    orientation_index = state % assignments
    orientation = [isodd(orientation_index >>> (label - 1)) ? -1 : 1
                   for label in 1:it.n]
    ((KnotGaussDiagram(_unrank_rooted_signed(it.n, code_index), orientation)), state + 1)
end

enumerate_knot_gauss(n::Integer; orientations::Bool=true) =
    KnotGaussEnumeration(Int(n), orientations)
enumerate_knot_gauss_diagrams(n::Integer; orientations::Bool=true) =
    enumerate_knot_gauss(n; orientations=orientations)
enumerate_rooted_normalized_signed_gc_with_orientations(n::Integer) =
    enumerate_knot_gauss(n; orientations=true)
