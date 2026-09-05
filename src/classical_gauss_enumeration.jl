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

"""Summary of one direct rooted unsigned-word enumeration."""
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

"""Number `(2n-1)!!` of rooted first-occurrence restricted-growth words."""
function rooted_unsigned_gauss_count(n::Integer)
    value = _cge_validate_n(n)
    result = big(1)
    for crossing in 1:value
        result *= 2crossing - 1
    end
    Int(result)
end

"""Call `f` once for every rooted unsigned double-occurrence word.

The recursion keeps labels in first-occurrence restricted-growth order.  At a
position it closes each currently open label in increasing label order, then
opens the next fresh label.  Hence each word is reached by one and only one
branch and no all-pairing superset is constructed.
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

foreach_rooted_unsigned_gauss_word(n::Integer, f::Function) =
    foreach_rooted_unsigned_gauss_word(f, n)

"""Materialized convenience view of the rooted unsigned word stream."""
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

"""Decode a compact key into mutable signed code and orientation vectors."""
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

"""Return a repository diagram when available, otherwise a two-field tuple."""
function classical_gauss_key_to_diagram(key::UInt64, n::Integer)
    code, orientation = classical_gauss_key_to_code(key, n)
    if isdefined(@__MODULE__, :KnotGaussDiagram)
        return KnotGaussDiagram(code, orientation)
    end
    (gauss=code, orientation=orientation)
end

"""Enumerate the complete direct classical raw set for one crossing number."""
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
