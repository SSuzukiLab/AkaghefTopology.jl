# Source: MST25 Appendix A, definitions and Proposition A.1,
# [[MST25#^mst25-blk-0186]], [[MST25#^mst25-blk-0188]].
# Encoding and scope: docs/specs/ograph-symmetries.md.
# No integer/Z2 framing transport is asserted by these operations.
function _check_diagram_symmetry(v::VirtualLink; reflection=false)
    _validate_cut_markers(v.gauss)
    v.is_weighted && throw(ArgumentError("symmetry transport of framing weights is not specified"))
    reflection && any(iszero,v.orientation) &&
        throw(ArgumentError("plane reflection of stored virtual-crossing cyclic orders is not supported"))
end

"""
    switch_crossings!(v::VirtualLink) -> VirtualLink

Exchange over and under at every real crossing of a link presentation in place.

# Arguments
`v` must be a coherent unweighted diagram under the
[symmetry specification](../specs/ograph-symmetries.md#source-and-encoding).
Assigned weights and malformed cut markers throw `ArgumentError`.
Incoherent crossing labels can raise indexing errors. Virtual crossings and
cut strands are permitted.

# Returns
The same `VirtualLink`, with the crossing switch applied.

# Guarantees
1. The result has the crossing-switch action of the linked symmetry table,
   retaining virtual visits and component traversal directions.
   [proved/agent: by construction]
2. The link is modified in place; applying this operation twice restores
   its stored presentation. [proved/agent: by construction]

# Not guaranteed
This primitive does not check closed-normality or certify a combed manifold.
Use `mst_sign_change!` for the source-backed closed-normal interpretation.
It supplies no framing lift and does not perform plane reflection.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — the diagram crossing-switch operation.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); expected = switch_crossings(v); result = switch_crossings!(v);

julia> (result === v, v.gauss == expected.gauss, v.orientation == expected.orientation)
(true, true, true)
```

See also [`switch_crossings`](@ref), [`mst_sign_change!`](@ref).

Role: move
"""
function switch_crossings!(v::VirtualLink)
    _check_diagram_symmetry(v)
    v.gauss=[[x==0 || v.orientation[abs(x)]==0 ? x : -x for x in w] for w in v.gauss]
    v.orientation .*= -1
    v
end
"""
    switch_crossings(v::VirtualLink) -> VirtualLink

Return the diagram with over and under exchanged at every real crossing.

# Arguments
`v` must be a coherent unweighted diagram under the
[symmetry specification](../specs/ograph-symmetries.md#source-and-encoding).
Assigned weights throw `ArgumentError` identifying unspecified framing
transport. Malformed cut markers throw `ArgumentError`; incoherent labels can
raise indexing errors. Virtual crossings and cut strands are permitted.

# Returns
A fresh `VirtualLink` with independent mutable data and the named crossing
switch applied.

# Guarantees
1. The result has the crossing-switch action of the linked symmetry table;
   virtual visits and component traversal directions are retained.
   [proved/agent: by construction]
2. Applying the operation twice restores the stored presentation. The input
   is unchanged. [proved/agent: by construction]

# Not guaranteed
This primitive does not check closed-normality or certify a combed manifold.
Use `mst_sign_change` for the source-backed closed-normal interpretation.
It supplies no integer or Z2 framing lift and is distinct from plane reflection.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — the diagram crossing-switch operation.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); w = switch_crossings(v);

julia> (w.gauss == [-word for word in v.gauss], w.orientation == -v.orientation)
(true, true)
```

See also [`switch_crossings!`](@ref), [`mst_sign_change`](@ref).

Role: move
"""
switch_crossings(v::VirtualLink)=switch_crossings!(deepcopy(v))

"""
    reflect_diagram!(v::VirtualLink) -> VirtualLink

Reflect the plane of a real-crossing diagram in place.

# Arguments
`v` must be coherent and unweighted with no explicit virtual crossings.
Weights, virtual crossings and malformed cut markers throw `ArgumentError`.
Closed-normality is not checked; cut strands are permitted.

# Returns
The same `VirtualLink` with the plane-reflection action of the
[symmetry table](../specs/ograph-symmetries.md#source-and-encoding).

# Guarantees
1. The diagram is reflected with passage signs and component traversal retained. [proved/agent: by construction]
2. The input is modified in place; applying reflection twice restores its presentation data. [proved/agent: by construction]

# Not guaranteed
This primitive does not certify a combed-manifold interpretation or provide
a framing lift; use [`mst_reflect_reverse!`](@ref) for the meaning-bearing entrypoint.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — diagram reflection, with repository encoding in the linked symmetry specification.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); expected = reflect_diagram(v);

julia> (reflect_diagram!(v) === v, v.gauss == expected.gauss, v.orientation == expected.orientation)
(true, true, true)
```

Role: move
"""
function reflect_diagram!(v::VirtualLink)
    _check_diagram_symmetry(v;reflection=true)
    v.orientation .*= -1
    v
end
"""
    reflect_diagram(v::VirtualLink) -> VirtualLink

Return the plane-reflected real-crossing diagram.

# Arguments
`v` must be coherent and unweighted, with no explicit virtual crossings,
under the [diagram-symmetry scope](../specs/ograph-symmetries.md#source-and-encoding).
Weights or virtual crossings throw `ArgumentError` identifying the unsupported
transport. Malformed cut markers throw `ArgumentError`; other representation
errors can propagate. Closed-normality is not checked.

# Returns
A fresh `VirtualLink` carrying the named reflection, with independent mutable
storage and the same field meanings as the input.

# Guarantees
1. The result has the reflection action in the
   [source and encoding table](../specs/ograph-symmetries.md#source-and-encoding).
   [proved/agent: by construction]
2. Reflection is an involution on the stored presentation. The input is unchanged.
   [proved/agent: by construction]

# Not guaranteed
This diagram operation alone asserts no combed-manifold interpretation and
supplies no framing lift. For the source-backed reflected-and-reversed
closed-normal operation, use `mst_reflect_reverse`.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — reflection of the diagram, distinct from exchanging over and under.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); w = reflect_diagram(v);

julia> (w.gauss == v.gauss, w.orientation == -v.orientation)
(true, true)
```

See also [`reflect_diagram!`](@ref), [`mst_reflect_reverse`](@ref).

Role: move
"""
reflect_diagram(v::VirtualLink)=reflect_diagram!(deepcopy(v))

"""
    reverse_strands!(v::VirtualLink) -> VirtualLink

Reverse every strand's direction in a link presentation in place.

# Arguments
`v` must be coherent and unweighted. Assigned weights or malformed cut
markers throw `ArgumentError`; virtual crossings and cut strands are accepted.
See the [symmetry scope](../specs/ograph-symmetries.md#source-and-encoding).

# Returns
The same `VirtualLink` after reversal of each strand traversal.

# Guarantees
1. The result follows the linked reversal convention, including endpoint exchange on cut strands, with crossing signs retained. [proved/agent: by construction]
2. The argument is modified in place; applying reversal twice restores its presentation data. [proved/agent: by construction]

# Not guaranteed
This does not certify a combed manifold or define integer/Z2 framing
transport or a new direction-reversed 0–2 move.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — strand reversal and its repository encoding.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); expected = reverse_strands(v);

julia> (reverse_strands!(v) === v, v.gauss == expected.gauss, v.orientation == expected.orientation)
(true, true, true)
```

Role: move
"""
function reverse_strands!(v::VirtualLink)
    _check_diagram_symmetry(v)
    v.gauss=[_is_cut_component(w) ? vcat(reverse(w[1:end-1]),0) : reverse(w) for w in v.gauss]
    v
end
"""
    reverse_strands(v::VirtualLink) -> VirtualLink

Return the diagram with every strand's direction reversed.

# Arguments
`v` must be a coherent unweighted `VirtualLink` under the
[symmetry specification](../specs/ograph-symmetries.md#source-and-encoding).
Weighted inputs throw `ArgumentError` because framing transport is unspecified.
Malformed cut markers throw `ArgumentError`; other invalid data can propagate
representation errors. Virtual crossings and cut strands are permitted.

# Returns
A fresh `VirtualLink` with independent mutable data and the named reversal
of all its strands.

# Guarantees
1. The result obeys the reversal convention in the linked symmetry table,
   including endpoint exchange for cut strands. [proved/agent: by construction]
2. Applying reversal twice restores the stored presentation. The input and
   its crossing signs are unchanged. [proved/agent: by construction]

# Not guaranteed
This does not assert a closed-normal or framed interpretation, introduce a
direction-reversed 0–2 move, or imply preservation of a combing. For the
source-backed combined operation use `mst_reflect_reverse`.

# Sources
cf:[[MST25#^mst25-blk-0186]] Appendix A — reversal of all strand orientations.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); w = reverse_strands(v);

julia> (w.gauss == reverse.(v.gauss), w.orientation == v.orientation)
(true, true)
```

See also [`reverse_strands!`](@ref), [`mst_reflect_reverse`](@ref).

Role: move
"""
reverse_strands(v::VirtualLink)=reverse_strands!(deepcopy(v))

function _check_mst_symmetry(v::VirtualLink)
    _check_diagram_symmetry(v;reflection=true)
    first(is_closed(v)) || throw(ArgumentError("MST semantic methods require a closed normal o-graph"))
end

"""
    mst_sign_change!(v::VirtualLink) -> VirtualLink

Replace a closed normal o-graph by the MST sign-change presentation of the opposite manifold and negated transported combing.

# Arguments
`v` must be an unweighted real-crossing presentation accepted by `is_closed`,
under the [meaning-bearing symmetry API](../specs/ograph-symmetries.md#meaning-bearing-api).
Weighted input, virtual crossings or failed closedness throw `ArgumentError`
identifying the restriction. Malformed representation errors propagate.

# Returns
The same `VirtualLink` with the MST sign change applied.

# Guarantees
1. For an input representing `(M, v)`, the output represents `(M^op, -v^op)` under the linked interpretation. [proved/literature: [[MST25#^mst25-blk-0188]] Proposition A.1]
2. The result has the sign-change action of the linked symmetry table; applying the operation twice restores the stored presentation. [proved/agent: by construction]
3. The argument is modified in place. [proved/agent: by construction]

# Conventions
Use the [source and encoding table](../specs/ograph-symmetries.md#source-and-encoding).

# Not guaranteed
No integer or Z2 framing transport is provided. Passing the operational
closedness guard does not prove its general source correspondence.

# Sources
cf:[[MST25#^mst25-blk-0188]] Proposition A.1 — opposite manifold and negated transported combing.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); expected = mst_sign_change(v);

julia> (mst_sign_change!(v) === v, v.gauss == expected.gauss, v.orientation == expected.orientation)
(true, true, true)
```

Role: move
"""
function mst_sign_change!(v::VirtualLink)
    _check_mst_symmetry(v)
    switch_crossings!(v)
end
"""
    mst_sign_change(v::VirtualLink) -> VirtualLink

Return the MST sign-change presentation of the oppositely oriented manifold with negated transported combing.

# Arguments
`v` must be an unweighted real-crossing presentation accepted by `is_closed`,
as required by the [meaning-bearing symmetry API](../specs/ograph-symmetries.md#meaning-bearing-api).
Weighted input, explicit virtual crossings or a failed closed-normal gate
throw `ArgumentError` identifying that restriction. Malformed cut markers
or other incoherent data can propagate validation errors.

# Returns
A fresh `VirtualLink` with the same field meanings as the input and the
MST sign-change operation applied. Its nested mutable data are independent
of the input.

# Guarantees
1. The array action is the named MST sign change in the
   [source and encoding table](../specs/ograph-symmetries.md#source-and-encoding).
   [proved/agent: by construction]
2. For a closed normal input representing `(M, v)`, the output represents
   `(M^op, -v^op)` in the sense of Proposition A.1 and the linked specification.
   [proved/literature: [[MST25#^mst25-blk-0188]] Proposition A.1]
3. The input is unchanged; applying the operation twice restores its stored
   presentation data. [proved/agent: by construction]

# Not guaranteed
The computational closedness guard is not a formal proof of its own
correctness. This operation supplies no integer or Z2 framing lift and
does not assert preservation of the original combed structure.

# Sources
cf:[[MST25#^mst25-blk-0188]] Appendix A, Proposition A.1 — the combed 3-manifold represented by the sign-change operation.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); w = mst_sign_change(v);

julia> (first(is_closed(w)), mst_sign_change(w).gauss == v.gauss)
(true, true)
```

See also [`mst_reflect_reverse`](@ref), [`switch_crossings`](@ref).

Role: convert
"""
mst_sign_change(v::VirtualLink)=mst_sign_change!(deepcopy(v))

"""
    mst_reflect_reverse!(v::VirtualLink) -> VirtualLink

Replace a closed normal o-graph by the reflected and reversed presentation of the opposite manifold and transported combing.

# Arguments
`v` must be an unweighted real-crossing presentation accepted by `is_closed`,
under the [meaning-bearing symmetry API](../specs/ograph-symmetries.md#meaning-bearing-api).
Weighted input, virtual crossings or failed closedness throw `ArgumentError`
identifying the restriction. Malformed representation errors propagate.

# Returns
The same `VirtualLink` with reflection and strand reversal applied.

# Guarantees
1. For an input representing `(M, v)`, the output represents `(M^op, v^op)` under the linked interpretation. [proved/literature: [[MST25#^mst25-blk-0188]] Proposition A.1]
2. The result has the reflected-and-reversed action of the linked symmetry table; applying the operation twice restores the stored presentation. [proved/agent: by construction]
3. The argument is modified in place. [proved/agent: by construction]

# Conventions
Use the [source and encoding table](../specs/ograph-symmetries.md#source-and-encoding).

# Not guaranteed
No integer or Z2 framing transport is provided. Passing the operational
closedness guard does not prove its general source correspondence.

# Sources
cf:[[MST25#^mst25-blk-0188]] Proposition A.1 — opposite manifold and transported combing.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); expected = mst_reflect_reverse(v);

julia> (mst_reflect_reverse!(v) === v, v.gauss == expected.gauss, v.orientation == expected.orientation)
(true, true, true)
```

Role: move
"""
function mst_reflect_reverse!(v::VirtualLink)
    _check_mst_symmetry(v)
    reflect_diagram!(v)
    reverse_strands!(v)
end
"""
    mst_reflect_reverse(v::VirtualLink) -> VirtualLink

Return the MST reflected and reversed presentation of the oppositely oriented manifold with transported combing.

# Arguments
`v` must be an unweighted real-crossing presentation accepted by `is_closed`,
under the [meaning-bearing symmetry API](../specs/ograph-symmetries.md#meaning-bearing-api).
Weighted input, explicit virtual crossings or a failed closed-normal gate
throw `ArgumentError` identifying the restriction. Malformed data may propagate
validation, indexing or diagram errors.

# Returns
A fresh `VirtualLink` carrying the reflected and strand-reversed presentation,
with nested mutable data independent of `v`.

# Guarantees
1. The output has the array action assigned to this operation by the
   [symmetry specification](../specs/ograph-symmetries.md#source-and-encoding).
   [proved/agent: by construction]
2. For a closed normal input representing `(M, v)`, the output represents
   `(M^op, v^op)` as stated by Proposition A.1 and the linked interpretation.
   [proved/literature: [[MST25#^mst25-blk-0188]] Proposition A.1]
3. The input is unchanged; applying the operation twice restores its stored
   presentation data. [proved/agent: by construction]

# Not guaranteed
No integer or Z2 framing transport is defined. Passing the existing
computational closedness guard is not a formal proof of its correctness.
The result is distinct from the combing negation in `mst_sign_change`.

# Sources
cf:[[MST25#^mst25-blk-0188]] Appendix A, Proposition A.1 — the combed 3-manifold represented by reflection followed by orientation reversal.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); w = mst_reflect_reverse(v);

julia> (first(is_closed(w)), w.gauss == reverse.(v.gauss))
(true, true)
```

See also [`mst_sign_change`](@ref), [`reflect_diagram`](@ref), [`reverse_strands`](@ref).

Role: convert
"""
mst_reflect_reverse(v::VirtualLink)=mst_reflect_reverse!(deepcopy(v))

"""
    mirror_manifold!(v::VirtualLink) -> VirtualLink

Apply the legacy crossing-switch operation to a link presentation.

# Arguments
`v` must satisfy the unweighted diagram preconditions of
[`switch_crossings!`](@ref). Assigned weights throw `ArgumentError` naming
unspecified framing transport; malformed cut markers and label errors
propagate. Closed-normality is not checked.

# Returns
The same `VirtualLink` with the crossing-switch action.

# Guarantees
1. The result is exactly the diagram operation `switch_crossings!(v)`. [proved/agent: by construction]
2. The argument is modified in place. [proved/agent: by construction]

# Conventions
The [symmetry table](../specs/ograph-symmetries.md#source-and-encoding)
distinguishes crossing switch from plane reflection.

# Not guaranteed
The legacy name does not certify a manifold mirror or combing transport.
Use [`mst_sign_change!`](@ref) for the checked closed-normal entrypoint.

# Sources
cf:docs/specs/ograph-symmetries.md — legacy alias and the distinction between diagram and meaning-bearing operations.

# Examples
```jldoctest
julia> v = lens_space_link(2, 1); set_weight!(v, zeros(4));

julia> try mirror_manifold!(v); false catch e; e isa ArgumentError end
true
```

Role: move
"""
mirror_manifold!(v::VirtualLink)=switch_crossings!(v)
