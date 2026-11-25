"""Exact classical-realizability checks for one-component Gauss data.

There are two deliberately separate routes in this file.

* `classical_gauss_certificate(code, orientation)` uses the Carter surface of
  the signed Gauss paragraph.  A positive repository visit is the overcrossing
  (the arrow tail), a negative visit is the undercrossing (the arrow head), and
  `orientation[i]` is the `+`/`-` exponent carried by the undercrossing.  The
  Carter surface is the least-genus oriented surface carrying that diagram, so
  genus zero is an exact classical-realizability decision.
* `unsigned_gauss_certificate(word)` forgets over/under and crossing signs.
  It exhausts the two alternating cyclic orders at every 4-valent crossing in
  the arc graph of the double-occurrence word.  A genus-zero orientable ribbon
  map is exactly a spherical immersed-circle realization.  The enumeration is
  finite and constructive, and does not use `planarize_virtual_gauss` as a
  realizability oracle.

The signed route follows V. Kurlin, *Gauss paragraphs of classical links and a
characterization of virtual link groups*, arXiv:math/0610929v2,
Definition 3.1, Lemma 3.2, and Theorem 3.6:
https://arxiv.org/abs/math/0610929

The unsigned route is the rotation-system form of the graph criterion in
M. L. Marx, *The Gauss realizability problem*, Proc. Amer. Math. Soc. 22
(1969), Theorem 2 (DOI 10.1090/S0002-9939-1969-0244984-6):
https://doi.org/10.1090/S0002-9939-1969-0244984-6
The graph has one vertex per crossing and one edge per arc between consecutive
visits.  The prescribed smoothing pairs at a vertex must be opposite in a
plane rotation; enumerating the two possible orientations of that alternating
order gives all admissible plane maps.  The face permutation is
`rotation * edge-reversal`, so `V - E + F = 2 - 2g` is the certificate.

This file is intentionally standalone: tests and scripts include it after
`knot_gauss.jl` and the package's `Topology.jl` is not modified here.
"""

const CLASSICAL_GAUSS_REALIZABILITY_VERSION = 1
const KURLIN_GAUSS_PARAGRAPH_URL = "https://arxiv.org/abs/math/0610929"
const MARX_GAUSS_REALIZABILITY_URL =
    "https://doi.org/10.1090/S0002-9939-1969-0244984-6"

"""    ClassicalGaussCertificate(code, orientation, route, realizable, genus, vertices, edges, faces, rotation_bits, rotation, face_cycles) -> ClassicalGaussCertificate

Store the coordinate and Euler data of an oriented ribbon-map certificate.

# Arguments
The fields use the types and meanings in [Ribbon records](../specs/knot-gauss.md#ribbon-records).
K6 ([conditions](../specs/knot-gauss.md#numbered-conditions)) is the caller's
responsibility: this field constructor does not validate a certificate.
Julia converts arguments to the declared field types; unsupported conversions
raise `MethodError` or an appropriate conversion exception.

# Returns
A `ClassicalGaussCertificate` containing `code`, `orientation`, `route`,
`realizable`, `genus`, `vertices`, `edges`, `faces`, `rotation_bits`,
`rotation` and `face_cycles`, as defined in the linked record entry.

# Guarantees
1. The returned fields equal the supplied values after field-type conversion. For K6 inputs they describe the capped oriented ribbon surface of the recorded rotation system. [proved/agent: by construction]
2. Construction does not modify the supplied arrays; arrays already of the declared type are retained and can subsequently be mutated through aliases. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
An arbitrary constructed record need not satisfy K6. Use
`classical_gauss_certificate` to obtain consistent data. The general literature
identification of the signed genus-zero test with classical realizability is
unverified here; the missing source and rotation-convention comparison are
recorded in the specification. A stored flag is not itself that proof.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — local surface and record definitions; unresolved literature comparison is stated there.

# Examples
```jldoctest
julia> c = classical_gauss_certificate(KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, 1, 1]));

julia> (c isa ClassicalGaussCertificate, c.genus, c.faces)
(true, 0, 5)
```

See also [`classical_gauss_certificate`](@ref).

Role: construct
"""
struct ClassicalGaussCertificate
    code::Tuple{Vararg{Int}}
    orientation::Tuple{Vararg{Int}}
    route::Symbol
    realizable::Bool
    genus::Int
    vertices::Int
    edges::Int
    faces::Int
    rotation_bits::Tuple{Vararg{Int}}
    rotation::Vector{NTuple{4,Int}}
    face_cycles::Vector{Vector{Int}}
end

"""    ClassicalGaussImmersion(word, rotation_bits, rotation, face_cycles, genus) -> ClassicalGaussImmersion

Store a rotation-map representation of an unsigned immersed-circle shadow.

# Arguments
The five fields have the meanings in [Ribbon records](../specs/knot-gauss.md#ribbon-records).
K5 ([conditions](../specs/knot-gauss.md#numbered-conditions)) is the caller's
responsibility for a spherical immersion: no consistency or genus check is
performed. Julia field conversions may raise `MethodError` or conversion errors.

# Returns
A `ClassicalGaussImmersion` with `word`, `rotation_bits`, `rotation`,
`face_cycles` and `genus` as specified by the linked record definition.

# Guarantees
1. The fields retain the supplied data after conversion. For K5 inputs they encode a spherical rotation map with the prescribed circle traversal. [proved/agent: by construction]
2. Input arrays are not modified; correctly typed arrays remain shared with the result. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
This constructor does not prove K5, embed a curve in Euclidean coordinates,
or choose over/under decorations. The specification records the unresolved
literature comparison for the general planar-realizability criterion.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — local rotation-map meaning.

# Examples
```jldoctest
julia> immersion = Topology.KnotGauss.construct_classical_immersion([1, 1]);

julia> (immersion isa ClassicalGaussImmersion, immersion.genus)
(true, 0)
```

See also `construct_classical_immersion`.

Role: construct
"""
struct ClassicalGaussImmersion
    word::Tuple{Vararg{Int}}
    rotation_bits::Tuple{Vararg{Int}}
    rotation::Vector{NTuple{4,Int}}
    face_cycles::Vector{Vector{Int}}
    genus::Int
end

function _cg_validate_unsigned(code)
    raw = Int[x for x in code]
    iseven(length(raw)) ||
        throw(ArgumentError("an unsigned Gauss word must have even length"))
    all(!=(0), raw) ||
        throw(ArgumentError("Gauss-word labels cannot be zero"))
    word = abs.(raw)
    n = length(word) ÷ 2
    n == 0 && return word
    sort!(unique(word)) == collect(1:n) ||
        throw(ArgumentError("Gauss-word labels must be exactly 1:$n"))
    all(label -> count(==(label), word) == 2, 1:n) ||
        throw(ArgumentError("every Gauss-word label must occur exactly twice"))
    word
end

function _cg_validate_signed(code, orientation)
    raw = Int[x for x in code]
    signs = Int[x for x in orientation]
    iseven(length(raw)) ||
        throw(ArgumentError("a signed Gauss code must have even length"))
    n = length(raw) ÷ 2
    length(signs) == n ||
        throw(DimensionMismatch("orientation must contain one sign per crossing"))
    all(x -> x == 1 || x == -1, signs) ||
        throw(ArgumentError("crossing orientation signs must be +1 or -1"))
    all(!=(0), raw) ||
        throw(ArgumentError("Gauss-code labels cannot be zero"))
    n == 0 && return raw, signs
    sort!(unique(abs.(raw))) == collect(1:n) ||
        throw(ArgumentError("Gauss-code labels must be exactly 1:$n"))
    for label in 1:n
        count(==(label), raw) == 1 && count(==(-label), raw) == 1 ||
            throw(ArgumentError(
                "crossing $label must occur once as over (+$label) and once as under (-$label)"))
    end
    raw, signs
end

"""Return the dart leaving word position `p` and the dart entering it."""
_cg_leave_dart(p::Int) = 2 * (p - 1) + 1
_cg_predecessor_enter_dart(p::Int, length_word::Int) =
    2 * (mod1(p - 1, length_word) - 1) + 2

function _cg_edge_reversal(n::Int)
    reversal = Vector{Int}(undef, 4n)
    for edge in 1:2n
        source = 2 * (edge - 1) + 1
        target = source + 1
        reversal[source] = target
        reversal[target] = source
    end
    reversal
end

"""Compute the boundary cycles of a cyclic rotation system."""
function _cg_face_cycles(rotation::Vector{NTuple{4,Int}}, n::Int)
    n == 0 && return Vector{Vector{Int}}()
    dart_count = 4n
    rotation_map = fill(0, dart_count)
    for order in rotation
        length(unique(order)) == 4 ||
            throw(ArgumentError("a crossing rotation must contain four distinct darts"))
        all(1 <= dart <= dart_count for dart in order) ||
            throw(ArgumentError("rotation contains an out-of-range dart"))
        for index in 1:4
            rotation_map[order[index]] = order[mod1(index + 1, 4)]
        end
    end
    all(!=(0), rotation_map) ||
        throw(ArgumentError("rotation does not cover every dart exactly once"))
    reversal = _cg_edge_reversal(n)
    face_next = rotation_map[reversal]
    seen = falses(dart_count)
    cycles = Vector{Vector{Int}}()
    for start in 1:dart_count
        seen[start] && continue
        cycle = Int[]
        dart = start
        while !seen[dart]
            push!(cycle, dart)
            seen[dart] = true
            dart = face_next[dart]
        end
        dart == start ||
            throw(ArgumentError("face permutation did not close a cycle"))
        push!(cycles, cycle)
    end
    cycles
end

"""Count Carter boundary components without materializing their dart lists."""
function _cg_face_count(rotation::Vector{NTuple{4,Int}}, n::Int)
    n == 0 && return 2
    dart_count = 4n
    rotation_map = fill(0, dart_count)
    for order in rotation
        for index in 1:4
            rotation_map[order[index]] = order[mod1(index + 1, 4)]
        end
    end
    reversal = _cg_edge_reversal(n)
    face_next = rotation_map[reversal]
    seen = falses(dart_count)
    faces = 0
    for start in 1:dart_count
        seen[start] && continue
        faces += 1
        dart = start
        while !seen[dart]
            seen[dart] = true
            dart = face_next[dart]
        end
    end
    faces
end

function _cg_genus(vertices::Int, edges::Int, faces::Int)
    chi = vertices - edges + faces
    iseven(2 - chi) ||
        throw(ArgumentError("ribbon-map Euler characteristic has odd parity"))
    genus = (2 - chi) ÷ 2
    genus >= 0 || throw(ArgumentError("ribbon map has negative genus"))
    genus
end

function _cg_unsigned_positions(word::Vector{Int})
    n = length(word) ÷ 2
    positions = [Int[] for _ in 1:n]
    for (position, label) in enumerate(word)
        push!(positions[label], position)
    end
    positions
end

"""Build one of the two admissible alternating rotations at every crossing."""
function _cg_unsigned_rotation(word::Vector{Int}, bits)
    n = length(word) ÷ 2
    length(bits) == n || throw(DimensionMismatch("one rotation bit per crossing"))
    length_word = length(word)
    positions = _cg_unsigned_positions(word)
    rotation = Vector{NTuple{4,Int}}(undef, n)
    for label in 1:n
        p, q = positions[label]
        in_p = _cg_predecessor_enter_dart(p, length_word)
        out_p = _cg_leave_dart(p)
        in_q = _cg_predecessor_enter_dart(q, length_word)
        out_q = _cg_leave_dart(q)
        bit = Int(bits[label])
        bit == 0 || bit == 1 ||
            throw(ArgumentError("rotation bits must be 0 or 1"))
        rotation[label] = bit == 0 ? (in_p, in_q, out_p, out_q) :
            (in_p, out_q, out_p, in_q)
    end
    rotation
end

"""Build the Carter-surface rotation for a repository signed code.

For `orientation[label] == +1`, the counter-clockwise order is
`in_over, in_under, out_over, out_under`; for `-1` it is its reverse.  This is
the local `always turn left` convention of Kurlin Definition 3.1 after mapping
the repository's positive visit to the arrow tail and negative visit to the
arrow head.
"""
function _cg_carter_rotation(code::Vector{Int}, orientation::Vector{Int})
    n = length(orientation)
    length_word = length(code)
    rotation = Vector{NTuple{4,Int}}(undef, n)
    for label in 1:n
        over_position = findfirst(==(label), code)
        under_position = findfirst(==(-label), code)
        over_position === nothing && error("internal code validation failure")
        under_position === nothing && error("internal code validation failure")
        in_over = _cg_predecessor_enter_dart(over_position, length_word)
        out_over = _cg_leave_dart(over_position)
        in_under = _cg_predecessor_enter_dart(under_position, length_word)
        out_under = _cg_leave_dart(under_position)
        rotation[label] = orientation[label] == 1 ?
            (in_over, in_under, out_over, out_under) :
            (in_over, out_under, out_over, in_under)
    end
    rotation
end

function _cg_make_certificate(code::Vector{Int}, orientation::Vector{Int},
                              route::Symbol, rotation_bits::Tuple,
                              rotation::Vector{NTuple{4,Int}})
    n = route == :unsigned_ribbon ? length(code) ÷ 2 : length(orientation)
    cycles = _cg_face_cycles(rotation, n)
    faces = n == 0 ? 2 : length(cycles)
    genus = _cg_genus(n, 2n, faces)
    ClassicalGaussCertificate(Tuple(code), Tuple(orientation), route,
                              genus == 0, genus, n, 2n, faces,
                              rotation_bits, rotation, cycles)
end

"""    classical_gauss_certificate(code, orientation) -> ClassicalGaussCertificate

Return the oriented ribbon-surface certificate of signed Gauss data or an unsigned word.

# Arguments
With two arguments, `code` and `orientation` are integer-convertible sequences
subject to K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)).
The one-argument form accepts an object with `gauss` and `orientation`
properties and checks those data, or an `AbstractVector`/`Tuple` subject to K8
for the unsigned route. These conditions are checked: a sign-count mismatch
raises `DimensionMismatch`; invalid word length, labels, visits or signs raise
`ArgumentError` naming the failure. An unsupported one-argument object raises
`ArgumentError` requesting signed data. Integer conversion failures propagate;
values must fit `Int`, including their absolute labels.

# Returns
A `ClassicalGaussCertificate` with every field interpreted by
[Ribbon records](../specs/knot-gauss.md#ribbon-records). The signed route uses
the supplied crossing data; the unsigned route returns a genus-zero witness
if one exists, otherwise a least-genus witness among its permitted rotations.

# Guarantees
1. The certificate satisfies K6; its Euler data describe the capped oriented ribbon surface of the stored rotation system. [proved/agent: by construction]
2. In the signed route, `realizable` is true exactly when that ribbon genus is zero. In the unsigned route it is true exactly when an admissible rotation has genus zero. [proved/agent: by construction]
3. The input containers are unchanged; the certificate owns newly created rotation and face-cycle arrays. [proved/agent: by construction]

# Conventions
[Ribbon coordinates and empty circle](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
The general literature identification with classical realizability is
unverified at the source-comparison level stated in the specification. A
ribbon certificate does not provide a coordinate drawing or a knot-type
classification. Unsigned search can be exponential in the crossing count.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — combinatorial surface and certificate semantics; source comparison remains recorded there.

# Examples
```jldoctest
julia> c = classical_gauss_certificate(KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, 1, 1]));

julia> (c.realizable, c.genus, c.faces)
(true, 0, 5)
```

See also [`classical_gauss_realizable`](@ref).

Role: convert
"""
function classical_gauss_certificate(code, orientation)
    signed_code, signs = _cg_validate_signed(code, orientation)
    n = length(signs)
    n == 0 && return _cg_make_certificate(signed_code, signs, :carter_surface,
                                           (), NTuple{4,Int}[])
    rotation = _cg_carter_rotation(signed_code, signs)
    # The orientation signs are the fixed local choices in the Carter map.
    _cg_make_certificate(signed_code, signs, :carter_surface, Tuple(signs), rotation)
end

"""    classical_gauss_certificate(diagram) -> ClassicalGaussCertificate

Return the oriented ribbon-surface certificate of signed Gauss data or an unsigned word.

# Arguments
With two arguments, `code` and `orientation` are integer-convertible sequences
subject to K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)).
The one-argument form accepts an object with `gauss` and `orientation`
properties and checks those data, or an `AbstractVector`/`Tuple` subject to K8
for the unsigned route. These conditions are checked: a sign-count mismatch
raises `DimensionMismatch`; invalid word length, labels, visits or signs raise
`ArgumentError` naming the failure. An unsupported one-argument object raises
`ArgumentError` requesting signed data. Integer conversion failures propagate;
values must fit `Int`, including their absolute labels.

# Returns
A `ClassicalGaussCertificate` with every field interpreted by
[Ribbon records](../specs/knot-gauss.md#ribbon-records). The signed route uses
the supplied crossing data; the unsigned route returns a genus-zero witness
if one exists, otherwise a least-genus witness among its permitted rotations.

# Guarantees
1. The certificate satisfies K6; its Euler data describe the capped oriented ribbon surface of the stored rotation system. [proved/agent: by construction]
2. In the signed route, `realizable` is true exactly when that ribbon genus is zero. In the unsigned route it is true exactly when an admissible rotation has genus zero. [proved/agent: by construction]
3. The input containers are unchanged; the certificate owns newly created rotation and face-cycle arrays. [proved/agent: by construction]

# Conventions
[Ribbon coordinates and empty circle](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
The general literature identification with classical realizability is
unverified at the source-comparison level stated in the specification. A
ribbon certificate does not provide a coordinate drawing or a knot-type
classification. Unsigned search can be exponential in the crossing count.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — combinatorial surface and certificate semantics; source comparison remains recorded there.

# Examples
```jldoctest
julia> c = classical_gauss_certificate(KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, 1, 1]));

julia> (c.realizable, c.genus, c.faces)
(true, 0, 5)
```

See also [`classical_gauss_realizable`](@ref).

Role: convert
"""
function classical_gauss_certificate(diagram)
    if hasproperty(diagram, :gauss) && hasproperty(diagram, :orientation)
        return classical_gauss_certificate(getproperty(diagram, :gauss),
                                           getproperty(diagram, :orientation))
    end
    # A one-argument call is the explicitly unsigned route.
    (diagram isa AbstractVector || diagram isa Tuple) &&
        return unsigned_gauss_certificate(diagram)
    throw(ArgumentError("expected a signed Gauss code and orientation"))
end

"""Exhaustive rotation-system certificate for an unsigned double-occurrence word."""
function unsigned_gauss_certificate(code)
    word = _cg_validate_unsigned(code)
    n = length(word) ÷ 2
    n == 0 && return _cg_make_certificate(word, Int[], :unsigned_ribbon,
                                           (), NTuple{4,Int}[])
    best = nothing
    # BigInt keeps the decision exact for arbitrary finite words; the N<=5
    # campaign only visits at most 2^5 local rotation choices per word.
    for mask in BigInt(0):((BigInt(1) << n) - 1)
        bits = ntuple(label -> Int((mask >> (label - 1)) & 1), n)
        rotation = _cg_unsigned_rotation(word, bits)
        candidate = _cg_make_certificate(word, Int[], :unsigned_ribbon,
                                         bits, rotation)
        best === nothing && (best = candidate; continue)
        candidate.genus < best.genus ||
            (candidate.genus == best.genus && candidate.rotation_bits < best.rotation_bits) ||
            continue
        best = candidate
    end
    best
end

"""Return whether the unsigned word has a spherical immersed-circle realization."""
unsigned_gauss_realizable(code) =
    unsigned_gauss_certificate(code).realizable

is_unsigned_gauss_realizable(code) =
    unsigned_gauss_realizable(code)

"""
    classical_gauss_realizable(code, orientation) -> Bool

Return whether the signed Gauss rotation surface has genus zero.

# Arguments
The signed forms use K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)):
integer-convertible sequences, or an object with `gauss` and `orientation`
properties. Validation is performed; wrong sign count raises `DimensionMismatch`,
and invalid lengths, labels, visits or signs raise `ArgumentError` naming the
failure. Conversion errors propagate; absolute labels must fit `Int`.
A single `AbstractVector` or `Tuple` instead uses the unsigned K8 route, checked with `ArgumentError` on invalid data. Other unsupported objects raise `ArgumentError`.

# Returns
`Bool` for the genus-zero decision of the corresponding
[ribbon surface](../specs/knot-gauss.md#ribbon-records).

# Guarantees
1. For signed inputs the result equals `classical_gauss_certificate(code, orientation).realizable`. [proved/agent: by construction]
2. The input data are unchanged. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
No witness, coordinate embedding or knot-type classification is returned.
The general literature identification with classical realizability retains
the unresolved source/convention comparison recorded in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — genus-zero predicate and its mathematical scope.

# Examples
```jldoctest
julia> d = KnotGaussDiagram([1, -2, -1, 2], [-1, 1]);

julia> classical_gauss_realizable(d)
false
```

See also [`classical_gauss_certificate`](@ref).

Role: condition
"""
classical_gauss_realizable(code, orientation) =
    classical_gauss_certificate(code, orientation).realizable

"""
    classical_gauss_realizable(diagram) -> Bool

Return whether the signed Gauss rotation surface has genus zero.

# Arguments
The signed forms use K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)):
integer-convertible sequences, or an object with `gauss` and `orientation`
properties. Validation is performed; wrong sign count raises `DimensionMismatch`,
and invalid lengths, labels, visits or signs raise `ArgumentError` naming the
failure. Conversion errors propagate; absolute labels must fit `Int`.
A single `AbstractVector` or `Tuple` instead uses the unsigned K8 route, checked with `ArgumentError` on invalid data. Other unsupported objects raise `ArgumentError`.

# Returns
`Bool` for the genus-zero decision of the corresponding
[ribbon surface](../specs/knot-gauss.md#ribbon-records).

# Guarantees
1. The result equals `classical_gauss_certificate(diagram).realizable`. [proved/agent: by construction]
2. The input data are unchanged. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
No witness, coordinate embedding or knot-type classification is returned.
The general literature identification with classical realizability retains
the unresolved source/convention comparison recorded in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — genus-zero predicate and its mathematical scope.

# Examples
```jldoctest
julia> d = KnotGaussDiagram([1, -2, -1, 2], [-1, 1]);

julia> classical_gauss_realizable(d)
false
```

See also [`classical_gauss_certificate`](@ref).

Role: condition
"""
classical_gauss_realizable(diagram) = classical_gauss_certificate(diagram).realizable
is_classical_gauss_realizable(args...) = classical_gauss_realizable(args...)

"""
    classical_gauss_realizable_fast(code, orientation) -> Bool

Return whether the signed Gauss rotation surface has genus zero.

# Arguments
The signed forms use K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)):
integer-convertible sequences, or an object with `gauss` and `orientation`
properties. Validation is performed; wrong sign count raises `DimensionMismatch`,
and invalid lengths, labels, visits or signs raise `ArgumentError` naming the
failure. Conversion errors propagate; absolute labels must fit `Int`.
The one-argument form requires both properties and rejects other objects with `ArgumentError`.

# Returns
`Bool` for the genus-zero decision of the corresponding
[ribbon surface](../specs/knot-gauss.md#ribbon-records).

# Guarantees
1. For signed inputs the result equals `classical_gauss_certificate(code, orientation).realizable`. [proved/agent: by construction]
2. The input data are unchanged. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
No witness, coordinate embedding or knot-type classification is returned.
The general literature identification with classical realizability retains
the unresolved source/convention comparison recorded in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — genus-zero predicate and its mathematical scope.

# Examples
```jldoctest
julia> d = KnotGaussDiagram([1, -2, -1, 2], [-1, 1]);

julia> classical_gauss_realizable_fast(d)
false
```

See also [`classical_gauss_certificate`](@ref).

Role: condition
"""
function classical_gauss_realizable_fast(code, orientation)
    signed_code, signs = _cg_validate_signed(code, orientation)
    n = length(signs)
    n == 0 && return true
    rotation = _cg_carter_rotation(signed_code, signs)
    _cg_face_count(rotation, n) == n + 2
end

"""
    classical_gauss_realizable_fast(diagram) -> Bool

Return whether the signed Gauss rotation surface has genus zero.

# Arguments
The signed forms use K1 ([conditions](../specs/knot-gauss.md#numbered-conditions)):
integer-convertible sequences, or an object with `gauss` and `orientation`
properties. Validation is performed; wrong sign count raises `DimensionMismatch`,
and invalid lengths, labels, visits or signs raise `ArgumentError` naming the
failure. Conversion errors propagate; absolute labels must fit `Int`.
The one-argument form requires both properties and rejects other objects with `ArgumentError`.

# Returns
`Bool` for the genus-zero decision of the corresponding
[ribbon surface](../specs/knot-gauss.md#ribbon-records).

# Guarantees
1. The result equals `classical_gauss_certificate(diagram).realizable`. [proved/agent: by construction]
2. The input data are unchanged. [proved/agent: by construction]

# Conventions
[Ribbon coordinates](../specs/knot-gauss.md#ribbon-records).

# Not guaranteed
No witness, coordinate embedding or knot-type classification is returned.
The general literature identification with classical realizability retains
the unresolved source/convention comparison recorded in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Ribbon records — genus-zero predicate and its mathematical scope.

# Examples
```jldoctest
julia> d = KnotGaussDiagram([1, -2, -1, 2], [-1, 1]);

julia> classical_gauss_realizable_fast(d)
false
```

See also [`classical_gauss_certificate`](@ref).

Role: condition
"""
function classical_gauss_realizable_fast(diagram)
    hasproperty(diagram, :gauss) && hasproperty(diagram, :orientation) ||
        throw(ArgumentError("expected a signed Gauss code and orientation"))
    classical_gauss_realizable_fast(getproperty(diagram, :gauss),
                                    getproperty(diagram, :orientation))
end

"""Construct a combinatorial plane immersion from an unsigned word.

The returned rotation system and face cycles are the construction; no
coordinate embedding or knot-type identification is asserted.
"""
function construct_classical_immersion(code)
    certificate = unsigned_gauss_certificate(code)
    certificate.realizable ||
        throw(ArgumentError("unsigned Gauss word is not classically realizable"))
    ClassicalGaussImmersion(certificate.code, certificate.rotation_bits,
                            certificate.rotation, certificate.face_cycles,
                            certificate.genus)
end

classical_gauss_immersion(code) =
    construct_classical_immersion(code)

"""Attach one of the `2^n` over/under choices to a chosen shadow.

`over_first[label]` says that the first occurrence of `label` in `word` is the
overcrossing.  The local rotation supplied by the chosen shadow determines the
crossing sign: swapping which branch is over swaps that sign.  Thus every
over/under decoration of a fixed planar immersion is classical, while an
arbitrary independent repository `orientation` vector need not be compatible
with that immersion (the Carter certificate decides compatibility).

If `KnotGaussDiagram` is loaded in the including module, a diagram value is
returned; otherwise a named tuple with the same two fields is returned.
"""
function decorate_classical_immersion(code, over_first;
                                       immersion=nothing)
    word = _cg_validate_unsigned(code)
    n = length(word) ÷ 2
    length(over_first) == n ||
        throw(DimensionMismatch("one over/under choice per crossing"))
    chosen = immersion === nothing ? construct_classical_immersion(word) : immersion
    chosen.word == Tuple(word) ||
        throw(ArgumentError("the supplied immersion belongs to another word"))
    bits = chosen.rotation_bits
    local_sign = [bits[label] == 0 ? 1 : -1 for label in 1:n]
    seen = falses(n)
    signed_code = Int[]
    for label in word
        first = !seen[label]
        first && (seen[label] = true)
        over = first ? over_first[label] : !over_first[label]
        push!(signed_code, over ? label : -label)
    end
    orientation = [over_first[label] ? local_sign[label] : -local_sign[label]
                   for label in 1:n]
    if isdefined(@__MODULE__, :KnotGaussDiagram)
        return KnotGaussDiagram(signed_code, orientation)
    end
    (gauss=signed_code, orientation=orientation)
end

"""Enumerate all classical over/under decorations of a chosen immersion."""
function enumerate_classical_decorations(code;
                                          immersion=nothing)
    word = _cg_validate_unsigned(code)
    n = length(word) ÷ 2
    result = Any[]
    for mask in BigInt(0):((BigInt(1) << n) - 1)
        over_first = [isodd(mask >> (label - 1)) for label in 1:n]
        push!(result, decorate_classical_immersion(word, over_first;
                                                   immersion=immersion))
    end
    result
end
