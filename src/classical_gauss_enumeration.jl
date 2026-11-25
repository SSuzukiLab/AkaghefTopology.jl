"""Direct enumeration of classical rooted signed one-component Gauss data.

This lane is a generator, not a knot-type classifier.  It avoids scanning the
42,577,920 arbitrary signed inputs at `n == 6` by generating exactly the
classical inputs from their planar shadows.

Completeness is the following finite rotation-system argument.  A rooted
unsigned double-occurrence word has one graph vertex for each crossing and
one edge for each arc between consecutive visits.  At each crossing the two
visits are the prescribed opposite smoothing pairs, so an oriented planar
rotation has exactly two possible alternating cyclic orders.  Enumerating all
`2^n` local choices therefore includes every planar immersion of the word.
For a connected map with `V = n`, `E = 2n`, the face permutation has a sphere
exactly when `F = n + 2`; smoothing the opposite pairs then recovers the
original one-component traversal.  Conversely, every planar immersed-circle
realization induces one of those local choices.  Each planar immersion has
all `2^n` choices of which visit is over.  With the repository convention
`+i = over/tail`, `-i = under/head`, the crossing orientation is the Kurlin
undercrossing exponent determined by the selected local rotation.  Thus every
classical signed code occurs in the generated set, and every generated code
has a genus-zero Carter rotation.

The signed convention follows V. Kurlin, *Gauss paragraphs of classical links
and a characterization of virtual link groups*, arXiv:math/0610929v2,
Definition 3.1 and Lemma 3.2:
https://arxiv.org/abs/math/0610929
The unsigned graph criterion is the rotation-system form of M. L. Marx,
*The Gauss realizability problem*, Proc. Amer. Math. Soc. 22 (1969),
Theorem 2:
https://doi.org/10.1090/S0002-9939-1969-0244984-6

The implementation is self-contained apart from the optional
`KnotGaussDiagram` constructor used by `classical_gauss_key_to_diagram` when
`knot_gauss.jl` has already been included.
"""

const CLASSICAL_GAUSS_ENUMERATION_VERSION = 1

"""    ClassicalGaussEnumerationResult(n, unsigned_words, rotations_checked, spherical_immersions, generated_decorations, raw_duplicates, raw_keys) -> ClassicalGaussEnumerationResult

Store the finite counts and signed-word keys of a direct Gauss enumeration.

# Arguments
The counters are `Int` fields and `raw_keys` is a `Set{UInt64}` with the meanings
in [Enumeration records](../specs/knot-gauss.md#enumeration-records).
K7 ([conditions](../specs/knot-gauss.md#numbered-conditions)) is the caller's
responsibility; no run-consistency checks are performed. Incompatible field
conversions raise Julia conversion errors.

# Returns
A `ClassicalGaussEnumerationResult` storing `n`, `unsigned_words`,
`rotations_checked`, `spherical_immersions`, `generated_decorations`,
`raw_duplicates` and `raw_keys` as defined in the linked record entry.

# Guarantees
1. Each field equals its supplied value after conversion; under K7 these are counts and keys of the stated finite rooted signed-diagram population. [proved/agent: by construction]
2. The input set is not modified; a `Set{UInt64}` is retained by reference. [proved/agent: by construction]

# Conventions
[Enumeration keys](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
Construction does not enumerate anything, certify the counters, or classify
knot isotopy types. The general literature comparison for classicality is
unverified as recorded in the specification.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration records — local finite-population and key definitions.

# Examples
```jldoctest
julia> r = direct_classical_gauss_enumeration(1);

julia> (r.unsigned_words, r.rotations_checked, r.spherical_immersions, length(r.raw_keys))
(1, 2, 2, 4)
```

See also [`direct_classical_gauss_enumeration`](@ref).

Role: construct
"""
struct ClassicalGaussEnumerationResult
    n::Int
    unsigned_words::Int
    rotations_checked::Int
    spherical_immersions::Int
    generated_decorations::Int
    raw_duplicates::Int
    raw_keys::Set{UInt64}
end

function _cge_validate_n(n::Integer)
    value = Int(n)
    0 <= value <= 6 ||
        throw(ArgumentError("direct classical Gauss enumeration is bounded by n=6"))
    value
end

"""    rooted_unsigned_gauss_count(n::Integer) -> Int

Return the number of rooted unsigned Gauss words with the requested crossing count.

# Arguments
`n` satisfies K4 ([conditions](../specs/knot-gauss.md#numbered-conditions)).
The range is checked with `ArgumentError`; integer conversion failures
propagate.

# Returns
An `Int` giving the size of the
[rooted unsigned population](../specs/knot-gauss.md#enumeration-and-move-witnesses).

# Guarantees
1. The result equals `length(enumerate_rooted_unsigned_gauss_words(n))` in the supported range, including the single empty word at zero crossings. [proved/agent: by construction]
2. The input is unchanged and no word collection is returned. [proved/agent: by construction]

# Not guaranteed
This count includes nonplanar words. It is neither a count of canonical
cyclic-root classes nor a count of knot types.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration and move witnesses — finite rooted unsigned population.

# Examples
```jldoctest
julia> Tuple(rooted_unsigned_gauss_count(n) for n in 0:4)
(1, 1, 3, 15, 105)
```

Role: access
"""
function rooted_unsigned_gauss_count(n::Integer)
    value = _cge_validate_n(n)
    result = big(1)
    for crossing in 1:value
        result *= 2crossing - 1
    end
    Int(result)
end

"""    foreach_rooted_unsigned_gauss_word(f::Function, n::Integer) -> Nothing

Visit every rooted unsigned Gauss word with the requested crossing count.

# Arguments
`n` satisfies K4 and `f` satisfies K10
([conditions](../specs/knot-gauss.md#numbered-conditions)) for one word vector.
The reversed positional order `(n, f)` is also supported. An invalid range
raises `ArgumentError`; integer conversion and callback failures propagate.

# Returns
`nothing` after all callbacks return normally.

# Guarantees
1. Each word of the specified rooted unsigned population is passed to `f` exactly once, in the same order as `enumerate_rooted_unsigned_gauss_words(n)`. [proved/agent: by construction]
2. Callback words are independently allocated vectors; a callback can retain or modify its vector without changing later words. No input container is modified by the enumerator. [proved/agent: by construction]

# Not guaranteed
Callback results are discarded and do not stop enumeration. Callback side
effects are not undone after an exception. Words are not certified planar
and do not carry crossing-sign or over/under data.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration and move witnesses — rooted unsigned population and ordering.

# Examples
```jldoctest
julia> words = Vector{Int}[];

julia> foreach_rooted_unsigned_gauss_word(w -> push!(words, w), 2)

julia> length(words)
3
```

See also [`enumerate_rooted_unsigned_gauss_words`](@ref).

Role: other
"""
function foreach_rooted_unsigned_gauss_word(f::Function, n::Integer)
    crossings = _cge_validate_n(n)
    length_word = 2crossings
    if crossings == 0
        f(Int[])
        return nothing
    end

    word = fill(0, length_word)
    occurrence_count = zeros(Int, crossings)
    word[1] = 1
    occurrence_count[1] = 1

    function extend!(position::Int, next_label::Int)
        if position > length_word
            next_label == crossings + 1 && f(copy(word))
            return nothing
        end

        # Closing an open label preserves the restricted-growth condition and
        # is ordered to make the stream deterministic.
        for label in 1:(next_label - 1)
            occurrence_count[label] == 1 || continue
            word[position] = label
            occurrence_count[label] = 2
            extend!(position + 1, next_label)
            occurrence_count[label] = 1
        end

        # A fresh label can be opened only after all earlier labels have been
        # introduced, and at most once for each label.
        if next_label <= crossings
            word[position] = next_label
            occurrence_count[next_label] = 1
            extend!(position + 1, next_label + 1)
            occurrence_count[next_label] = 0
        end
        nothing
    end

    extend!(2, 2)
    nothing
end

"""
    foreach_rooted_unsigned_gauss_word(n::Integer, f::Function) -> Nothing

Visit every rooted unsigned Gauss word with the requested crossing count.

# Arguments
`n` satisfies K4 and `f` satisfies K10
([conditions](../specs/knot-gauss.md#numbered-conditions)) for one word vector.
The reversed positional order `(n, f)` is also supported. An invalid range
raises `ArgumentError`; integer conversion and callback failures propagate.

# Returns
`nothing` after all callbacks return normally.

# Guarantees
1. Each word of the specified rooted unsigned population is passed to `f` exactly once, in the same order as `enumerate_rooted_unsigned_gauss_words(n)`. [proved/agent: by construction]
2. Callback words are independently allocated vectors; a callback can retain or modify its vector without changing later words. No input container is modified by the enumerator. [proved/agent: by construction]

# Not guaranteed
Callback results are discarded and do not stop enumeration. Callback side
effects are not undone after an exception. Words are not certified planar
and do not carry crossing-sign or over/under data.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration and move witnesses — rooted unsigned population and ordering.

# Examples
```jldoctest
julia> words = Vector{Int}[];

julia> foreach_rooted_unsigned_gauss_word(2, w -> push!(words, w))

julia> length(words)
3
```

See also [`enumerate_rooted_unsigned_gauss_words`](@ref).

Role: other
"""
foreach_rooted_unsigned_gauss_word(n::Integer, f::Function) =
    foreach_rooted_unsigned_gauss_word(f, n)

"""    enumerate_rooted_unsigned_gauss_words(n::Integer) -> Vector{Vector{Int}}

Return all rooted unsigned Gauss words with the requested crossing count.

# Arguments
`n` satisfies K4 ([conditions](../specs/knot-gauss.md#numbered-conditions)).
The range is checked with `ArgumentError`; integer conversion failures
propagate. No classical-realizability restriction is imposed.

# Returns
A `Vector{Vector{Int}}` of the rooted unsigned words defined in
[Enumeration and move witnesses](../specs/knot-gauss.md#enumeration-and-move-witnesses).

# Guarantees
1. Every word in that finite population occurs once, in the specified order; its size is `rooted_unsigned_gauss_count(n)`. [proved/agent: by construction]
2. Each word has its own mutable vector; modifying one result word does not modify another. No input is changed. [proved/agent: by construction]

# Not guaranteed
Unsigned words need not be realizable by plane curves. Root changes remain
distinct here; no over/under choices or crossing signs are supplied.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration and move witnesses — rooted unsigned population.

# Examples
```jldoctest
julia> length(enumerate_rooted_unsigned_gauss_words(2))
3
```

See also [`foreach_rooted_unsigned_gauss_word`](@ref).

Role: construct
"""
function enumerate_rooted_unsigned_gauss_words(n::Integer)
    words = Vector{Vector{Int}}()
    foreach_rooted_unsigned_gauss_word(n) do word
        push!(words, word)
    end
    words
end

rooted_unsigned_gauss_words(n::Integer) =
    enumerate_rooted_unsigned_gauss_words(n)

enumerate_rooted_unsigned_gauss(n::Integer) =
    enumerate_rooted_unsigned_gauss_words(n)

"""Precompute the four darts attached to every crossing visit pair."""
function _cge_word_arcs(word::Vector{Int}, n::Int)
    first_position = zeros(Int, n)
    second_position = zeros(Int, n)
    for (position, label) in enumerate(word)
        if first_position[label] == 0
            first_position[label] = position
        else
            second_position[label] = position
        end
    end

    in_first = Vector{Int}(undef, n)
    out_first = Vector{Int}(undef, n)
    in_second = Vector{Int}(undef, n)
    out_second = Vector{Int}(undef, n)
    for label in 1:n
        p = first_position[label]
        q = second_position[label]
        in_first[label] = p == 1 ? 2length(word) : 2p - 2
        out_first[label] = 2p - 1
        in_second[label] = q == 1 ? 2length(word) : 2q - 2
        out_second[label] = 2q - 1
    end
    (in_first=in_first, out_first=out_first,
     in_second=in_second, out_second=out_second)
end

"""Count faces for one local rotation bit mask using reusable buffers."""
function _cge_face_count(arcs, n::Int, rotation_mask::Int,
                          rotation_map::Vector{Int}, seen::AbstractVector{Bool})
    n == 0 && return 2
    for label in 1:n
        in_first = arcs.in_first[label]
        out_first = arcs.out_first[label]
        in_second = arcs.in_second[label]
        out_second = arcs.out_second[label]
        if ((rotation_mask >> (label - 1)) & 1) == 0
            # Opposite smoothing pairs: (in_first,out_first) and
            # (in_second,out_second).
            a, b, c, d = in_first, in_second, out_first, out_second
        else
            a, b, c, d = in_first, out_second, out_first, in_second
        end
        rotation_map[a] = b
        rotation_map[b] = c
        rotation_map[c] = d
        rotation_map[d] = a
    end

    fill!(seen, false)
    faces = 0
    for start in eachindex(seen)
        seen[start] && continue
        faces += 1
        dart = start
        while !seen[dart]
            seen[dart] = true
            opposite = isodd(dart) ? dart + 1 : dart - 1
            dart = rotation_map[opposite]
        end
    end
    faces
end

"""Pack one signed rooted code into a compact, deterministic `UInt64` key.

The key stores `n` in a four-bit prefix, then each visit as four bits
(`2label` for positive and `2label+1` for negative), followed by one bit per
crossing orientation (`0` for `+1`, `1` for `-1`).  The bounded N<=6 lane
therefore uses at most 62 bits and avoids retaining hundreds of thousands of
small tuples in the direct census.
"""
function _cge_pack_signed(code, orientation)
    n = length(orientation)
    _cge_validate_n(n)
    length(code) == 2n || throw(DimensionMismatch("signed code length is not 2n"))
    key = UInt64(n)
    for value_raw in code
        value = Int(value_raw)
        label = abs(value)
        1 <= label <= n || throw(ArgumentError("signed key label out of range"))
        token = UInt64(2label + (value < 0 ? 1 : 0))
        key = (key << 4) | token
    end
    for sign_raw in orientation
        sign = Int(sign_raw)
        sign == 1 || sign == -1 ||
            throw(ArgumentError("orientation must contain only +/-1"))
        key = (key << 1) | UInt64(sign == -1)
    end
    key
end

"""
    classical_gauss_key_to_code(key::UInt64, n::Integer) -> Tuple{Vector{Int},Vector{Int}}

Decode a packed signed-Gauss key into its visit and crossing-sign vectors.

# Arguments
`n` must satisfy K4 and `key` must satisfy K9
([conditions](../specs/knot-gauss.md#numbered-conditions)). The crossing-count
range, packed prefix and label ranges are checked, with `ArgumentError`
identifying the failed check; conversion of `n` to `Int` may fail.
Full decoded K1 validity is the caller's responsibility: this vector decoder does not check the signed occurrence counts.

# Returns
A pair of new `Vector{Int}` values `(code, orientation)` in the signed Gauss representation.

# Guarantees
1. The output has exactly the signed coordinates encoded by a K9 input key in the specified packed format. [proved/agent: by construction]
2. Input arguments are unchanged. [proved/agent: by construction]

# Conventions
[Packed-key format](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No classicality or knot-type decision is made. A structurally valid packed
key need not have arisen from a classical census. Cross-version format
compatibility is not certified by these routines.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration records and K9 — local packed-coordinate representation.

# Examples
```jldoctest
julia> r = direct_classical_gauss_enumeration(0); key = only(r.raw_keys);

julia> code, signs = classical_gauss_key_to_code(key, 0);

julia> (Tuple(code), Tuple(signs))
((), ())
```

Role: construct
"""
function classical_gauss_key_to_code(key::UInt64, n::Integer)
    crossings = _cge_validate_n(n)
    remainder = key
    orientation = Vector{Int}(undef, crossings)
    for label in crossings:-1:1
        orientation[label] = (remainder & UInt64(1)) == UInt64(1) ? -1 : 1
        remainder >>= 1
    end
    code = Vector{Int}(undef, 2crossings)
    for position in 2crossings:-1:1
        token = remainder & UInt64(0x0f)
        label = Int(token >> 1)
        bit = token & UInt64(1)
        1 <= label <= crossings || throw(ArgumentError("invalid packed label"))
        code[position] = bit == UInt64(1) ? -label : label
        remainder >>= 4
    end
    remainder == UInt64(crossings) ||
        throw(ArgumentError("packed key prefix does not match n"))
    code, orientation
end

"""
    classical_gauss_key_to_diagram(key::UInt64, n::Integer) -> KnotGaussDiagram

Decode a packed signed-Gauss key into a knot-diagram record.

# Arguments
`n` must satisfy K4 and `key` must satisfy K9
([conditions](../specs/knot-gauss.md#numbered-conditions)). The crossing-count
range, packed prefix and label ranges are checked, with `ArgumentError`
identifying the failed check; conversion of `n` to `Int` may fail.
Decoded K1 validity is also checked by `KnotGaussDiagram`; invalid signed pairing raises `ArgumentError`.

# Returns
A `KnotGaussDiagram` with the decoded `gauss` and `orientation` fields in the loaded library.

# Guarantees
1. The output has exactly the signed coordinates encoded by a K9 input key in the specified packed format. [proved/agent: by construction]
2. Input arguments are unchanged. [proved/agent: by construction]

# Conventions
[Packed-key format](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
No classicality or knot-type decision is made. A structurally valid packed
key need not have arisen from a classical census. Cross-version format
compatibility is not certified by these routines.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration records and K9 — local packed-coordinate representation.

# Examples
```jldoctest
julia> r = direct_classical_gauss_enumeration(0); key = only(r.raw_keys);

julia> d = classical_gauss_key_to_diagram(key, 0);

julia> (d.gauss, d.orientation)
((), ())
```

Role: construct
"""
function classical_gauss_key_to_diagram(key::UInt64, n::Integer)
    code, orientation = classical_gauss_key_to_code(key, n)
    if isdefined(@__MODULE__, :KnotGaussDiagram)
        return KnotGaussDiagram(code, orientation)
    end
    (gauss=code, orientation=orientation)
end

"""    direct_classical_gauss_enumeration(n::Integer; progress=nothing) -> ClassicalGaussEnumerationResult

Enumerate the rooted signed Gauss data obtained from spherical admissible rotations.

# Arguments
`n` must satisfy K4 ([conditions](../specs/knot-gauss.md#numbered-conditions));
the checked range failure raises `ArgumentError`, and `Int` conversion errors
propagate. An optional `progress` callback satisfies K10 for three cumulative
integer counts `(unsigned_words, spherical_immersions, generated_decorations)`;
it is called once per unsigned word. Unsupported calls and callback errors
propagate; callback side effects are the caller's responsibility.

# Returns
A `ClassicalGaussEnumerationResult` with counters and `raw_keys` as defined in
[Enumeration records](../specs/knot-gauss.md#enumeration-records).

# Guarantees
1. `raw_keys` is exactly the set of encoded signed decorations from all spherical admissible rotations of the rooted unsigned words with `n` crossings. Each key decodes to K1 data. [proved/agent: by construction]
2. The counters refer to that finite enumeration, and `generated_decorations == length(raw_keys) + raw_duplicates`. [proved/agent: by construction]
3. The routine creates a new result and key set and changes no input container; a supplied callback may perform its own side effects. [proved/agent: by construction]

# Conventions
[Rooted population](../specs/knot-gauss.md#enumeration-and-move-witnesses),
[ribbon choices](../specs/knot-gauss.md#ribbon-records) and
[packed keys](../specs/knot-gauss.md#enumeration-records).

# Not guaranteed
This is neither a canonical-root quotient nor a knot-isotopy census.
The general source identification with all classically realizable diagrams
retains the specification's unverified literature comparison. Memory and
run time can grow rapidly toward the supported upper crossing bound.

# Sources
cf:lib/topology.jl/docs/specs/knot-gauss.md, Enumeration and ribbon records — finite rotation/decorated-word population.

# Examples
```jldoctest
julia> r = direct_classical_gauss_enumeration(1);

julia> (r.spherical_immersions, r.generated_decorations, r.raw_duplicates, length(r.raw_keys))
(2, 4, 0, 4)
```

See also [`classical_gauss_key_to_diagram`](@ref).

Role: construct
"""
function direct_classical_gauss_enumeration(n::Integer; progress=nothing)
    crossings = _cge_validate_n(n)
    masks = 1 << crossings
    raw_keys = Set{UInt64}()
    unsigned_words = 0
    rotations_checked = 0
    spherical_immersions = 0
    generated_decorations = 0
    raw_duplicates = 0

    foreach_rooted_unsigned_gauss_word(crossings) do word
        unsigned_words += 1
        arcs = _cge_word_arcs(word, crossings)
        rotation_map = zeros(Int, 4crossings)
        seen_darts = falses(4crossings)
        code_buffer = Vector{Int}(undef, 2crossings)
        orientation_buffer = Vector{Int}(undef, crossings)
        seen_labels = falses(crossings)

        for rotation_mask in 0:(masks - 1)
            rotations_checked += 1
            faces = _cge_face_count(arcs, crossings, rotation_mask,
                                    rotation_map, seen_darts)
            faces == crossings + 2 || continue
            spherical_immersions += 1

            # Every over/under assignment on this fixed spherical immersion
            # is classical after deriving its crossing signs from the local
            # rotation.  The raw Set remains the final duplicate guard.
            for over_mask in 0:(masks - 1)
                generated_decorations += 1
                fill!(seen_labels, false)
                for (position, label) in enumerate(word)
                    first_visit = !seen_labels[label]
                    seen_labels[label] = true
                    first_over = ((over_mask >> (label - 1)) & 1) == 1
                    over = first_visit ? first_over : !first_over
                    code_buffer[position] = over ? label : -label
                end
                for label in 1:crossings
                    first_over = ((over_mask >> (label - 1)) & 1) == 1
                    rotation_reversal = ((rotation_mask >> (label - 1)) & 1) == 1
                    # local_sign is +1 for the first rotation and -1 for its
                    # reverse; switching the over branch flips that sign.
                    orientation_buffer[label] = first_over != rotation_reversal ?
                        1 : -1
                end
                key = _cge_pack_signed(code_buffer, orientation_buffer)
                if key in raw_keys
                    raw_duplicates += 1
                else
                    push!(raw_keys, key)
                end
            end
        end
        progress === nothing || progress(unsigned_words, spherical_immersions,
                                         generated_decorations)
    end

    ClassicalGaussEnumerationResult(crossings, unsigned_words,
                                    rotations_checked, spherical_immersions,
                                    generated_decorations, raw_duplicates,
                                    raw_keys)
end
