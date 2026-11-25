"""Immutable one-circuit presentation.

A nonempty value is the Gauss code of an unweighted normal o-graph.  The empty
word is the MST crossing-free circuit, outside that class.  [BP97, §1.4] N1
excludes it from `ClosedNormalOGraph`.

This value records representation-level data only.  Raw equality and `hash`
compare the stored root and labels exactly.  `canonicalize_ograph` and
`ograph_canonical_hash` identify cyclic re-rooting and crossing relabeling; they
do not identify MP-, PS-, CP-, mirror-, or orientation-reversal equivalence.

N2 ([BP97, §1.4]) holds by construction of one oriented cyclic word: each edge
is directed from one entry to the next, so opposite directions match through
each vertex.  A head-to-head mismatch cannot be stored.  The passage-sign rule
is the one in `virtual-link-rgauss/1`: a positive entry is the over passage.
The one-over-one-under check below is N1's passage mark, not N2.
"""

using SHA

const OGRAPH_PRESENTATION_SCHEMA = "akaghef.ograph/1"
const OGRAPH_PRESENTATION_CANON_VERSION = 1
const OGRAPH_PRESENTATION_CANONICALIZATION = "ograph-presentation/1"
const OGRAPH_DEFAULT_CONVENTION = "virtual-link-rgauss/1"
# Persisted evidence still uses this version.  N1 is a gate on
# `ClosedNormalOGraph`; it does not change the stored C2/C3 numbers.
# See docs/specs/ograph-terminology.md, "Persisted closedness evidence".
const OGRAPH_CLOSEDNESS_VERIFIER_VERSION = "virtual-link-c1-c3/2"
const N2_ENFORCEMENT_SITE = "OGraphState"   # N2 [BP97, §1.4], by construction
const C1_ENFORCEMENT_SITE = "ograph_state"  # C1 [BP97, §1.4], one cyclic word
const C2_GRAPH_SPEC_VERSION = "bp97-fig1.2-pd-slots/1"

"""
    OGraphState(gauss::AbstractVector{<:Integer}, vertex_signs::AbstractVector{<:Integer}; convention_id=OGRAPH_DEFAULT_CONVENTION) -> OGraphState

Construct an immutable one-circuit normal o-graph presentation or the empty MST circuit state.

# Arguments
`gauss` and `vertex_signs` must satisfy the representation checks assigned
to `OGraphState` in [Enforcement](../specs/ograph-terminology.md#enforcement),
using the [passage and vertex signs](../specs/ograph-terminology.md#4-passage-sign-and-vertex-sign).
Failed shape, passage or sign checks throw `ArgumentError` naming the failed
requirement. Entries are converted to `Int`; failed conversions propagate.
The tuple overload collects its inputs and uses the same vector constructor.
`convention_id::AbstractString` must be nonempty after stripping whitespace,
otherwise `ArgumentError` names that field. An unknown nonempty identifier
is accepted without interpreting it.

# Returns
`OGraphState` with `gauss::Tuple{Vararg{Int}}`,
`vertex_signs::Tuple{Vararg{Int}}` and the stripped `convention_id::String`.
These are the supplied presentation's word, vertex signs and convention label.

# Guarantees
1. The supplied order, labels and signs are retained after integer conversion;
   later mutation of input vectors cannot change the stored tuples.
   [proved/agent: by construction]
2. For a nonempty state, the result is the one-circuit normal presentation
   described in [Enforcement](../specs/ograph-terminology.md#enforcement),
   with N2 enforced by that representation. It has not been certified for C2–C3.
   [proved/agent: by construction]
3. Empty input vectors are accepted as the
   [MST crossing-free circuit state](../specs/ograph-terminology.md#mst-crossing-free-circuit),
   which is outside the closed-normal N1 class. [proved/agent: by construction]

# Conventions
The word and vertex-sign encoding is owned by
[OGraph terminology](../specs/ograph-terminology.md#2-standard-gauss-code-displays).

# Not guaranteed
The constructor does not canonicalize roots or labels, check C2–C3, interpret
arbitrary convention identifiers, or encode weights or cut boundaries.
Use `certify_closed_ograph` for the closed-normal gate. Raw equality is not
equivalence of represented spines or manifolds.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4, N1–N2 and C1–C3 — normal and closed-normal conditions; the one-word enforcement and empty MST extension are repository encoding.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> code, signs = real_gauss_code(vl_example("s3_abalone"));

julia> state = OGraphState(only(code), signs);

julia> (!Base.ismutabletype(typeof(state)), certify_closed_ograph(state).witness.criteria)
(true, (true, true, true))
```

See also [`ograph_state`](@ref), [`canonicalize_ograph`](@ref), [`certify_closed_ograph`](@ref).

Role: construct
"""
struct OGraphState
    gauss::Tuple{Vararg{Int}}
    vertex_signs::Tuple{Vararg{Int}}
    convention_id::String

    function OGraphState(gauss::AbstractVector{<:Integer},
                         vertex_signs::AbstractVector{<:Integer};
                         convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION)
        code = Int.(collect(gauss))
        signs = Int.(collect(vertex_signs))
        convention = strip(String(convention_id))
        isempty(convention) && throw(ArgumentError("convention_id must not be empty"))
        iseven(length(code)) || throw(ArgumentError(
            "N1 preamble: an o-graph Gauss code has even length (quadrivalence) [BP97, §1.4]"))
        all(!=(0), code) || throw(ArgumentError("o-graph real Gauss labels cannot be zero"))
        all(sign -> sign in (-1, 1), signs) || throw(ArgumentError(
            "o-graph vertex signs must be +1 or -1"))
        length(signs) == length(code) ÷ 2 || throw(ArgumentError(
            "vertex_signs must contain one sign per true vertex"))
        for label in eachindex(signs)
            occurrences = filter(value -> abs(value) == label, code)
            length(occurrences) == 2 || throw(ArgumentError(
                "N1 preamble: true vertex $label must occur exactly twice (quadrivalence) [BP97, §1.4]"))
            occurrences[1] * occurrences[2] < 0 || throw(ArgumentError(
                "N1: true vertex $label must occur once as an over passage and once as an under passage [BP97, §1.4]"))
        end
        new(Tuple(code), Tuple(signs), convention)
    end
end

"""
    OGraphState(gauss::Tuple, vertex_signs::Tuple; convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) -> OGraphState

Construct an immutable one-circuit presentation from tuple data.

# Arguments
The tuple entries and convention identifier must pass the vector constructor's
checks in [Enforcement](../specs/ograph-terminology.md#enforcement).
Invalid incidence, signs or an empty identifier throw `ArgumentError`;
conversion and dispatch errors propagate for unsuitable tuple element types.

# Returns
An `OGraphState` with the supplied word, signs and stripped convention label.

# Guarantees
1. Successful construction has exactly the vector constructor's representation meaning and field values. [proved/agent: by construction]
2. The input tuples are not modified; C2–C3 are not certified. [proved/agent: by construction]

# Not guaranteed
The empty word is the separate MST extension, not an N1 closed-normal graph.
Construction does not canonicalize or identify a manifold.

# Sources
cf:docs/specs/ograph-terminology.md#enforcement — one-circuit representation and normal-condition gates.

# Examples
```jldoctest
julia> state = OGraphState((-2,-1,1,2), (1,-1)); state.gauss
(-2, -1, 1, 2)
```

Role: construct
"""
OGraphState(gauss::Tuple, vertex_signs::Tuple;
            convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) =
    OGraphState(collect(gauss), collect(vertex_signs); convention_id=convention_id)

Base.copy(state::OGraphState) = OGraphState(
    state.gauss, state.vertex_signs; convention_id=state.convention_id)
Base.:(==)(left::OGraphState, right::OGraphState) =
    left.gauss == right.gauss &&
    left.vertex_signs == right.vertex_signs &&
    left.convention_id == right.convention_id
Base.hash(state::OGraphState, seed::UInt) =
    hash((state.gauss, state.vertex_signs, state.convention_id), seed)

"""
    true_vertex_count(state::OGraphState) -> Int

Return the number of true vertices in an immutable o-graph presentation.

# Arguments
`state` is an `OGraphState` satisfying its checked representation contract.
No further [closed-normal conditions](../specs/ograph-terminology.md#enforcement)
are required or tested.

# Returns
An `Int` equal to the number of stored vertex signs.

# Guarantees
1. The result counts true vertices of this presentation, with the meaning in
   [the terminology specification](../specs/ograph-terminology.md). It is
   nonnegative and equals half the Gauss-word length.
   [proved/agent: by construction]
2. The state is unchanged. [proved/agent: by construction]

# Not guaranteed
This is not a minimal complexity over moves or a manifold invariant.
A zero value does not certify the state as a BP97 normal o-graph.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §Enforcement — paired real visits and the true-vertex count.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> true_vertex_count(ograph_state(vl_example("s3_abalone")))
1
```

See also [`closedness_witness`](@ref), [`certify_closed_ograph`](@ref).

Role: access
"""
true_vertex_count(state::OGraphState) = length(state.vertex_signs)

function _ograph_relabel_from_root(state::OGraphState, root::Int)
    isempty(state.gauss) && return copy(state)
    labels = Dict{Int,Int}()
    signs = Int[]
    code = Int[]
    next_label = 1
    for offset in 0:length(state.gauss)-1
        value = state.gauss[mod1(root + offset, length(state.gauss))]
        old = abs(value)
        if !haskey(labels, old)
            labels[old] = next_label
            push!(signs, state.vertex_signs[old])
            next_label += 1
        end
        push!(code, sign(value) * labels[old])
    end
    OGraphState(code, signs; convention_id=state.convention_id)
end

_ograph_state_key(state::OGraphState) =
    Tuple(vcat(collect(state.gauss), [0], collect(state.vertex_signs)))

"""
    canonicalize_ograph(state::OGraphState) -> OGraphState

Return a representative of an o-graph presentation modulo cyclic re-rooting and crossing relabelling.

# Arguments
`state` must be a valid `OGraphState` under its constructor's
[representation contract](../specs/ograph-terminology.md#enforcement).
The empty state is accepted. No closedness or geometric-realizability test is made.

# Returns
An immutable `OGraphState` with the same convention identifier and the selected
word and vertex signs. Root and crossing labels may differ from the input.

# Guarantees
1. Applying canonicalization twice gives an equal state. Consistent crossing
   relabelling and cyclic re-rooting of the input give equal results.
   [proved/agent: by construction]
2. The result retains the same signed oriented cyclic presentation and vertex
   count; input data are unchanged. This is
   [presentation identity](../specs/ograph-terminology.md#persisted-closedness-evidence),
   not a quotient by geometric moves. [proved/agent: by construction]

# Not guaranteed
Reflection, reversal, MP, PS and CP are not identified. No manifold or combing
classification follows. The
[representative and serialization specification](../specs/ograph-terminology.md#presentation-canonicalization-and-serialization)
defines the coordinates and their compatibility boundary.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §Persisted closedness evidence — canonical presentation identity used by stored witnesses.

# Examples
```jldoctest
julia> a = OGraphState([-2,-1,1,2], [1,-1]); b = OGraphState([1,2,-2,-1], [1,-1]);

julia> (canonicalize_ograph(a) == canonicalize_ograph(b), canonicalize_ograph(canonicalize_ograph(a)) == canonicalize_ograph(a))
(true, true)
```

See also [`ograph_canonical_hash`](@ref), [`ograph_canonical_serialization`](@ref).

Role: move
"""
function canonicalize_ograph(state::OGraphState)
    isempty(state.gauss) && return copy(state)
    best = _ograph_relabel_from_root(state, 1)
    for root in 2:length(state.gauss)
        candidate = _ograph_relabel_from_root(state, root)
        _ograph_state_key(candidate) < _ograph_state_key(best) && (best = candidate)
    end
    best
end

"""
    ograph_canonical_serialization(state::OGraphState; canon_version::Integer=OGRAPH_PRESENTATION_CANON_VERSION) -> String

Return the domain-qualified text identity of a canonically rooted o-graph presentation.

# Arguments
`state` is a valid presentation. `canon_version` must equal the supported
canonicalization version; otherwise `ArgumentError` identifies the unsupported
version. C2–C3 are not checked.

# Returns
A `String` in the specified
[canonical presentation format](../specs/ograph-terminology.md#presentation-canonicalization-and-serialization).

# Guarantees
1. Re-rootings and consistent crossing relabellings identified by `canonicalize_ograph` have identical text. [proved/agent: by construction]
2. The text identifies the canonical presentation data, retaining its convention identifier, not a manifold-equivalence class. [proved/agent: by construction]
3. The state is not modified and repeated calls with the same supported version return identical text. [proved/agent: by construction]

# Not guaranteed
There is no public inverse parser or cross-version compatibility promise.
This does not quotient moves, reflection or direction reversal. The empty
state has the separate [MST circuit meaning](../specs/ograph-terminology.md#mst-crossing-free-circuit).

# Sources
cf:docs/specs/ograph-terminology.md — canonical presentation identity and persisted closedness evidence.

# Examples
```jldoctest
julia> state = OGraphState([-2,-1,1,2], [1,-1]);

julia> startswith(ograph_canonical_serialization(state), "OGraphState;schema=akaghef.ograph/1;canonicalization=ograph-presentation/1;")
true
```

See also [`ograph_canonical_hash`](@ref).

Role: convert
"""
function ograph_canonical_serialization(state::OGraphState;
        canon_version::Integer=OGRAPH_PRESENTATION_CANON_VERSION)
    canon_version == OGRAPH_PRESENTATION_CANON_VERSION || throw(ArgumentError(
        "unsupported o-graph canonicalization version $canon_version"))
    canonical = canonicalize_ograph(state)
    convention = canonical.convention_id
    "OGraphState;schema=$(OGRAPH_PRESENTATION_SCHEMA);" *
        "canonicalization=$(OGRAPH_PRESENTATION_CANONICALIZATION);" *
        "canon=v$(canon_version);" *
        "convention-bytes=$(ncodeunits(convention)):$(convention);" *
        "n=$(true_vertex_count(canonical));" *
        "gauss=$(join(canonical.gauss, ','));" *
        "vertex_signs=$(join(canonical.vertex_signs, ','))"
end

"""
    ograph_canonical_hash(state::OGraphState; canon_version=OGRAPH_PRESENTATION_CANON_VERSION) -> String

Return the fingerprint of the presentation modulo cyclic re-rooting and crossing relabelling.

# Arguments
`state::OGraphState` is a valid state. `canon_version` must be the supported
serialization version; an unsupported version throws `ArgumentError` naming
it. Other keywords are forwarded to `ograph_canonical_serialization` and
are unsupported unless accepted there.

# Returns
A 64-character lowercase hexadecimal `String`, the SHA-256 digest of the
UTF-8 canonical serialization of `state`.

# Guarantees
1. Cyclic re-rooting and consistent crossing relabelling leave the result
   unchanged. Repeated calls in this implementation are deterministic.
   [proved/agent: by construction]
2. The fingerprint concerns [presentation identity](../specs/ograph-terminology.md#enforcement),
   retaining vertex signs and the convention identifier; it does not
   identify manifolds or combed structures. The input is unchanged.
   [proved/agent: by construction]

# Not guaranteed
Hash equality is not a collision-free proof. The call does not quotient
MP, PS, CP, reflection or strand reversal, and does not check C2–C3.
It is distinct from the raw identity used for move replay.

# Sources
cf: lib/topology.jl/docs/specs/ograph-terminology.md §Persisted closedness evidence — the presentation hash used by replayable closedness records.

# Examples
```jldoctest
julia> code = [-2, -1, 1, 2]; signs = [1, -1];

julia> a = OGraphState(code, signs); b = OGraphState([1, 2, -2, -1], signs);

julia> ograph_canonical_hash(a) == ograph_canonical_hash(b)
true
```

See also [`ograph_canonical_serialization`](@ref), [`ograph_raw_hash`](@ref).

Role: invariant
"""
ograph_canonical_hash(state::OGraphState; kwargs...) =
    bytes2hex(SHA.sha256(codeunits(ograph_canonical_serialization(state; kwargs...))))

"""Reject raw cut-boundary inputs at closed o-graph API boundaries.

The raw `VirtualLink.gauss` representation keeps a cut as a trailing zero.
This check deliberately runs before `real_gauss_code`, whose normalized view
may otherwise hide that boundary from a closed-only caller.
"""
function _require_closed_ograph_input(link::VirtualLink, api::AbstractString)
    any(word -> !isempty(word) && last(word) == 0, link.gauss) &&
        throw(ArgumentError(
            "$api requires a closed normal o-graph; trailing cut boundary markers are unsupported"))
    nothing
end

"""
    ograph_state(link::VirtualLink; convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION) -> OGraphState

Return the immutable one-circuit normal presentation underlying an unweighted link.

# Arguments
`link` must be unweighted, in integer mode, with one closed real-Gauss component.
Checked: cuts, Z2 mode, assigned weights or failure of C1 throw `ArgumentError`
naming the restriction. The real-Gauss projection must pass the `OGraphState`
representation checks; their errors propagate. `convention_id` is a nonempty
string after stripping, checked by that constructor. C2–C3 are not required.
See [Enforcement](../specs/ograph-terminology.md#enforcement).

# Returns
An `OGraphState` with the projected Gauss word and vertex signs stored as
immutable integer tuples, and the supplied convention identifier.

# Guarantees
1. The stored data equal the link's real-Gauss projection with its existing root and crossing labels. [proved/agent: by construction]
2. The result satisfies the representation-enforced conditions of the linked specification; the empty word has the separate MST circuit meaning. [proved/agent: by construction]
3. The input and its mutable data are unchanged and are not retained by the state. [proved/agent: by construction]

# Not guaranteed
The call does not certify C2–C3 or preserve explicit virtual crossings.
Arbitrary convention identifiers are stored rather than interpreted.

# Sources
cf:docs/specs/ograph-terminology.md#enforcement — representation-level normal conditions and the C1 adapter gate.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [-2,1]);

julia> try ograph_state(v); false catch e; e isa ArgumentError end
true
```

See also [`certify_closed_ograph`](@ref), [`normal_coloured_ograph`](@ref).

Role: convert
"""
function ograph_state(link::VirtualLink;
                      convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION)
    _require_closed_ograph_input(link, "ograph_state")
    weight_mode(link)==:integer || throw(ArgumentError(
        "Z2 mode requires a typed weighted state adapter; it cannot be silently projected to an unweighted o-graph"))
    link.is_weighted && throw(ArgumentError(
        "weighted VirtualLink input must not be silently projected to an unweighted o-graph"))
    code, signs = real_gauss_code(link)
    # C1 ([BP97, §1.4]): deleting vertices and joining opposite edges leaves
    # one oriented circuit.  The Gauss components are that circuit decomposition.
    # An OGraphState holds one word, so the condition cannot fail inside the type.
    length(code) == 1 || throw(ArgumentError(
        "C1: deleting the vertices and joining opposite edges must leave one oriented circuit [BP97, §1.4]"))
    OGraphState(only(code), signs; convention_id=convention_id)
end

"""
    ClosednessWitness(canonical_hash, verifier_version, criteria, c1_oriented_circuits, c2_trivalent_components, c3_method_circuits, true_vertices, method_circuit_lengths) -> ClosednessWitness

Store closed-normal criterion data and representation-enforcement metadata.

# Arguments
Supply string-convertible hash and verifier version, a three-Boolean criterion
tuple, integer-convertible circuit/component and vertex counts, and an iterable
of integer-convertible circuit lengths. Conversion errors propagate. Agreement
of these values with one presentation is the caller's responsibility; this
constructor does not recompute the criteria.

# Returns
A `ClosednessWitness` retaining the supplied fields after conversion and
adding `n1`, `n2_enforcement` and `c1_enforcement` with the meanings in
[Enforcement](../specs/ograph-terminology.md#enforcement). The criterion order
and distinction between a word count and a smoothing count are owned there.

# Guarantees
1. Counts and lengths are stored as `Int` values and a tuple; mutable input
   containers are not retained. `n1` is derived from the supplied vertex count.
   [proved/agent: by construction]
2. When computed by [`closedness_witness`](@ref), the record concerns the
   normal presentation and the [MST zero-vertex boundary](../specs/ograph-terminology.md#mst-crossing-free-circuit)
   of that operation. Construction alone supplies no geometric certificate.
   [proved/agent: by construction]

# Not guaranteed
Negative counts, inconsistent flags, stale fingerprints and unknown verifier
versions are not rejected here. Even three true criterion bits do not suffice
for N1. `ClosedNormalOGraph` checks the record against its state.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4 — named closed-normal conditions; the MST zero-vertex branch is the repository extension.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> w = closedness_witness(ograph_state(vl_example("s3_abalone")));

julia> copy = ClosednessWitness(w.canonical_hash, w.verifier_version, w.criteria, w.c1_oriented_circuits, w.c2_trivalent_components, w.c3_method_circuits, w.true_vertices, w.method_circuit_lengths);

julia> (copy.criteria, copy.true_vertices, copy.method_circuit_lengths)
((true, true, true), 1, (5, 1))
```

See also [`closedness_witness`](@ref), [`ClosedNormalOGraph`](@ref).

Role: construct
"""
struct ClosednessWitness
    canonical_hash::String
    verifier_version::String
    criteria::NTuple{3,Bool}
    c1_oriented_circuits::Int
    c2_trivalent_components::Int
    c3_method_circuits::Int
    true_vertices::Int
    method_circuit_lengths::Tuple{Vararg{Int}}
    n1::Bool
    n2_enforcement::String
    c1_enforcement::String

    function ClosednessWitness(canonical_hash, verifier_version, criteria,
            c1_oriented_circuits, c2_trivalent_components, c3_method_circuits,
            true_vertices, method_circuit_lengths)
        vertices = Int(true_vertices)
        new(String(canonical_hash), String(verifier_version),
            NTuple{3,Bool}(criteria),
            Int(c1_oriented_circuits), Int(c2_trivalent_components),
            Int(c3_method_circuits), vertices,
            Tuple(Int(n) for n in method_circuit_lengths),
            vertices >= 1, N2_ENFORCEMENT_SITE, C1_ENFORCEMENT_SITE)
    end
end

"""
    C2Graph(node_count::Integer, edges::Vector{NTuple{2,Int}}, origins::Vector{NTuple{3,Int}}) -> C2Graph

Store an occurrence-preserving C2 multigraph with edge provenance.

# Arguments
`node_count` is nonnegative and every edge endpoint lies in its declared node
set; violations throw `ArgumentError` identifying the count or endpoint.
`edges` and `origins` must have equal length, otherwise `DimensionMismatch`
requires one origin per edge occurrence. Origin membership and the BP C2
construction are not checked; the caller supplies those meanings.

# Returns
A `C2Graph` with integer `node_count`, the supplied `edges` and row-aligned
`origins`. Their coordinates and lane provenance are owned by the
[C2 specification](../specs/ograph-c2.md#source-to-index-trace).

# Guarantees
1. Every stored endpoint is in the declared node set and the vectors have
   equal length. Loops, duplicate edges and input order are retained.
   The vectors are shared, not copied. [proved/agent: by construction]
2. A value produced by `c2_graph` is the BP Figure 1.2 graph of its normal
   presentation, under the [repository frame](../specs/ograph-c2.md#repository-frame).
   This constructor alone only validates graph storage.
   [proved/agent: by construction]

# Not guaranteed
Direct construction does not check trivalence, connectedness, origin ranges,
or correspondence to an o-graph, spine or manifold. Later input-vector mutation
can invalidate the checked storage conditions.

# Sources
cf:[[BP97#^bp97-sec-0016]] Figure 1.2 — the graph used by C2; the specification owns its index correspondence.

# Examples
```jldoctest
julia> g = c2_graph([[1,-1]], [1]); copy = C2Graph(g.node_count, g.edges, g.origins);

julia> copy.node_count
4
```

See also [`c2_graph`](@ref), [`closedness_witness`](@ref).

Role: construct
"""
struct C2Graph
    node_count::Int
    edges::Vector{NTuple{2,Int}}
    origins::Vector{NTuple{3,Int}}

    function C2Graph(node_count::Integer,
                     edges::Vector{NTuple{2,Int}},
                     origins::Vector{NTuple{3,Int}})
        node_count >= 0 || throw(ArgumentError("C2 node count must be non-negative"))
        length(edges) == length(origins) || throw(DimensionMismatch(
            "each C2 edge occurrence must have one origin"))
        all(edge -> all(node -> 1 <= node <= node_count, edge), edges) ||
            throw(ArgumentError("C2 edge endpoint is outside the declared node set"))
        new(Int(node_count), edges, origins)
    end
end

# Source and index trace: docs/specs/ograph-c2.md.  Port columns are
# (+i, -i, +o, -o); the two rows are vertex signs +1 and -1.
const _C2_LOCAL_TRIPLES = (
    ((3, 4, 1), (4, 1, 2), (1, 2, 3), (2, 3, 4)),
    ((1, 2, 3), (4, 1, 2), (3, 4, 1), (2, 3, 4)),
)
const _C2_PAIR_PERMUTATIONS = ((3, 2, 1), (2, 3, 1)) # same, mixed

function _validate_c2_gc(component_words, vertex_signs)
    words = [Int.(collect(word)) for word in component_words]
    signs = Int.(collect(vertex_signs))
    n = length(signs)
    n >= 1 || throw(ArgumentError(
        "N1: the Figure 1.2 graph is stated for a normal o-graph with at least one vertex [BP97, §1.4]"))
    !isempty(words) && all(!isempty, words) || throw(ArgumentError(
        "C2 component words must be nonempty cyclic words"))
    all(sign -> sign in (-1, 1), signs) || throw(ArgumentError(
        "C2 vertex signs must be +1 or -1"))
    entries = reduce(vcat, words; init=Int[])
    length(entries) == 2n || throw(ArgumentError(
        "C2 GC must contain exactly two entries per vertex"))
    all(entry -> entry != 0 && abs(entry) <= n, entries) || throw(ArgumentError(
        "C2 GC labels must be nonzero and lie in 1:n"))
    for label in 1:n
        count(==(label), entries) == 1 && count(==(-label), entries) == 1 ||
            throw(ArgumentError(
                "N1: C2 vertex $label must occur once as an over passage and once as an under passage [BP97, §1.4]"))
    end
    words, signs
end

"""
    c2_graph(component_words, vertex_signs) -> C2Graph

Return the BP Figure 1.2 multigraph of a complete normal-o-graph Gauss presentation.

# Arguments
Inputs must be integer-convertible component words and vertex signs satisfying
the [C2 input domain](../specs/ograph-c2.md) and N1 representation requirements
in [Enforcement](../specs/ograph-terminology.md#enforcement).
The function checks nonempty cyclic components, nonzero real signs, complete
signed visits and their label ranges. Failure throws `ArgumentError` identifying
the count, sign, label or visit condition; conversion errors propagate.
C1, C2 connectedness and C3 are not preconditions.

# Returns
A `C2Graph` with node count, edge-occurrence pairs and aligned provenance
triples in the [repository frame](../specs/ograph-c2.md#repository-frame).

# Guarantees
1. For `n` true vertices, there are `4n` nodes and `6n` edge occurrences.
   Loops and parallel edges are retained, with one origin per occurrence.
   Inputs are unchanged. [proved/agent: by construction]
2. The result is the source-indexed BP C2 multigraph of this presentation;
   its source correspondence is recorded in the [C2 specification](../specs/ograph-c2.md#source-to-index-trace).
   It need not be connected. [proved/agent: by construction]

# Not guaranteed
This does not certify C1–C3, infer a manifold, transport weights or produce a
simple graph. The zero-vertex MST extension is outside this construction.

# Sources
cf:[[BP97#^bp97-sec-0016]] Figure 1.2 — the multigraph used in the C2 condition.

# Examples
```jldoctest
julia> g = c2_graph([[1,-1]], [1]);

julia> (g.node_count, length(g.edges), Tuple(g.edges[1:3]))
(4, 6, ((1, 1), (2, 2), (3, 4)))
```

See also [`C2Graph`](@ref), [`ograph_graph_a`](@ref), [`closedness_witness`](@ref).

Role: construct
"""
function c2_graph(component_words, vertex_signs)
    words, signs = _validate_c2_gc(component_words, vertex_signs)
    n = length(signs)
    edge_count = 3sum(length, words)
    edges = Vector{NTuple{2,Int}}(undef, edge_count)
    origins = Vector{NTuple{3,Int}}(undef, edge_count)
    row = 0

    for (component, word) in enumerate(words), tail_occurrence in eachindex(word)
        a = word[tail_occurrence]
        b = word[mod1(tail_occurrence + 1, length(word))]
        left_port = a > 0 ? 3 : 4       # +o / -o
        right_port = b > 0 ? 1 : 2      # +i / -i
        left = _C2_LOCAL_TRIPLES[signs[abs(a)] > 0 ? 1 : 2][left_port]
        right = _C2_LOCAL_TRIPLES[signs[abs(b)] > 0 ? 1 : 2][right_port]
        permutation = _C2_PAIR_PERMUTATIONS[sign(a) == sign(b) ? 1 : 2]

        for lane in 1:3
            row += 1
            edges[row] = (
                4 * (abs(a) - 1) + left[lane],
                4 * (abs(b) - 1) + right[permutation[lane]],
            )
            origins[row] = (component, tail_occurrence, lane)
        end
    end
    C2Graph(4n, edges, origins)
end

"""
    c2_graph(state::OGraphState) -> C2Graph

Return the BP C2 multigraph of an immutable normal presentation.

# Arguments
`state` must satisfy N1 in [Enforcement](../specs/ograph-terminology.md#enforcement).
An empty state throws `ArgumentError` identifying N1. The component/sign
method's representation and conversion errors propagate; C2–C3 need not pass.

# Returns
A `C2Graph` carrying the node set, edge occurrences and provenance of this
state's single component, in the [C2 frame](../specs/ograph-c2.md#repository-frame).

# Guarantees
1. Its fields equal those of `c2_graph([collect(state.gauss)], collect(state.vertex_signs))`.
   The state is unchanged. [proved/agent: by construction]
2. It is the BP Figure 1.2 graph associated to that presentation, with loops
   and multiplicities retained; connectivity is a separate condition.
   [proved/agent: by construction]

# Not guaranteed
Construction does not certify C2 or C3, choose a canonical graph labeling, or
identify a 3-manifold. It does not accept the zero-vertex MST circuit.

# Sources
cf:[[BP97#^bp97-sec-0016]] Figure 1.2 — the C2 graph for a normal presentation.

# Examples
```jldoctest
julia> s = OGraphState([1,2,-1,-2], [1,1]);

julia> c2_graph(s).edges == c2_graph([[1,2,-1,-2]], [1,1]).edges
true
```

See also [`c2_graph`](@ref), [`closedness_witness`](@ref).

Role: convert
"""
c2_graph(state::OGraphState) =
    c2_graph([collect(state.gauss)], collect(state.vertex_signs))

"""
    c2_graph(link::VirtualLink) -> C2Graph

Return the BP C2 multigraph of a link's underlying real normal presentation.

# Arguments
`link` must have closed cyclic components, no explicit virtual crossings and
the complete nonempty real-Gauss data required by the component/sign method.
Cuts, virtual crossings or invalid representation data throw `ArgumentError`
identifying that restriction; conversion errors propagate. Weighted inputs
are accepted; C1 and C3 are not required.

# Returns
A `C2Graph` with node count, occurrence-preserving edges and provenance for
the real-Gauss components in the [C2 specification](../specs/ograph-c2.md).

# Guarantees
1. The fields equal the component/sign result on `real_gauss_code(link)`.
   The link and all its weights are unchanged. [proved/agent: by construction]
2. This is the BP Figure 1.2 graph of the underlying presentation; weights
   do not enter this graph or establish a weighted equivalence.
   [proved/agent: by construction]

# Not guaranteed
This is not a connectedness or closedness certificate, a framing transport,
or a manifold identifier. A crossing-free input is outside the C2 domain.

# Sources
cf:[[BP97#^bp97-sec-0016]] Figure 1.2 — the unweighted C2 construction.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> v = vl_example("s3_abalone"); set_weight!(v, [-2,1]);

julia> c2_graph(v).edges == c2_graph([[1,-1]], [1]).edges
true
```

See also [`ograph_graph_a`](@ref), [`closedness_witness`](@ref).

Role: convert
"""
function c2_graph(link::VirtualLink)
    _require_closed_ograph_input(link, "c2_graph")
    all(!=(0), link.orientation) || throw(ArgumentError(
        "C2 requires a real o-graph; virtual crossings are unsupported"))
    all(word -> all(!=(0), word), link.gauss) || throw(ArgumentError(
        "C2 requires closed cyclic words; cut markers are unsupported"))
    code, signs = real_gauss_code(link)
    c2_graph(code, signs)
end

"""
    ograph_graph_a(link::VirtualLink) -> SimpleGraph

Return the simple-graph projection of the presentation's C2 multigraph.

# Arguments
`link` must be a coherent closed-strand presentation. Cuts throw `ArgumentError`;
on a nonempty real projection, the restrictions and errors of [`c2_graph`](@ref)
also apply, including its rejection of explicit virtual crossings.
Weights and coefficient mode are ignored. C3 is not required.
See [C2 construction](../specs/ograph-c2.md).

# Returns
A `SimpleGraph` on the C2 nodes, with loops omitted and parallel edges merged;
the crossing-free presentation returns a graph with no vertices.

# Guarantees
1. On a nonempty presentation, the result has the same connected components as `c2_graph(link)`. [proved/agent: by construction]
2. For nonempty presentations its node set is the C2 node set, but it retains only simple adjacency. [proved/agent: by construction]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is a compatibility projection, not the authoritative C2 multigraph,
a spine, or a C3 certificate. The empty graph does not represent the special
MST circuit's operational trivalent-component count.

# Sources
cf:docs/specs/ograph-c2.md — C2 multigraph, source correspondence and compatibility projection.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> g = ograph_graph_a(v); (Topology.Graphs.nv(g), Topology.Graphs.ne(g))
(4, 3)
```

Role: convert
"""
function ograph_graph_a(link::VirtualLink)::SimpleGraph
    _require_closed_ograph_input(link, "ograph_graph_a")
    _, orientation = real_gauss_code(link)
    isempty(orientation) && return SimpleGraph(0)
    c2 = c2_graph(link)
    graph = SimpleGraph(c2.node_count)
    for (left, right) in c2.edges
        left == right || add_edge!(graph, left, right)
    end
    graph
end

"""
    ograph_graph_a(state::OGraphState) -> SimpleGraph

Return the simple-graph projection of a state's C2 multigraph.

# Arguments
`state` must satisfy its representation contract. The
[C2 construction](../specs/ograph-c2.md) is used without imposing C3;
its representation errors propagate.

# Returns
A fresh `SimpleGraph` on the C2 nodes, or `SimpleGraph(0)` for the empty state.

# Guarantees
1. For nonempty states, simple adjacency has the same connected components as the C2 multigraph, with loops and multiplicities discarded. [proved/agent: by construction]
2. The state is not modified. [proved/agent: by construction]

# Not guaranteed
The result is not the authoritative C2 multigraph or a C3 certificate.
The empty-state convention is not a reconstruction of the MST circuit.

# Sources
cf:docs/specs/ograph-c2.md — multigraph connectivity and the legacy simple projection.

# Examples
```jldoctest
julia> g = ograph_graph_a(OGraphState([1,-1], [1])); (Topology.Graphs.nv(g), Topology.Graphs.ne(g))
(4, 3)
```

Role: convert
"""
function ograph_graph_a(state::OGraphState)::SimpleGraph
    isempty(state.vertex_signs) && return SimpleGraph(0)
    c2 = c2_graph(state)
    graph = SimpleGraph(c2.node_count)
    for (left, right) in c2.edges
        left == right || add_edge!(graph, left, right)
    end
    graph
end

function _ograph_trivalent_component_count(link::VirtualLink)
    _, orientation = real_gauss_code(link)
    isempty(orientation) && return 1
    length(connected_components(ograph_graph_a(link)))
end

function _ograph_virtual_link(state::OGraphState)
    link = VirtualLink()
    # `collect(())` has element type `Union{}`.  Preserve the concrete integer
    # representation so the MST crossing-free state remains a valid
    # `OGraphState` rather than failing in `_normalize_gauss`.
    set_data!(link; rgauss=[Int[state.gauss...]],
              orientation=Int[state.vertex_signs...])
end

"""
    closedness_witness(state::OGraphState) -> ClosednessWitness

Return the closed-normal criterion record for a canonicalized o-graph presentation.

# Arguments
`state` must satisfy the `OGraphState` representation contract. C2 and C3 may
fail; they produce false criterion bits rather than rejection. Errors from
the [C2 construction](../specs/ograph-c2.md) or Method-B data propagate.
The empty state follows the separately specified
[MST boundary](../specs/ograph-terminology.md#mst-crossing-free-circuit).

# Returns
A `ClosednessWitness` with canonical fingerprint, current verifier version,
criteria, counts, circuit lengths and enforcement metadata, interpreted by
[Enforcement](../specs/ograph-terminology.md#enforcement).

# Guarantees
1. The fingerprint equals `ograph_canonical_hash(state)`. Consistent re-rooting
   and relabelling give equal witnesses. The state is unchanged.
   [proved/agent: by construction]
2. On nonempty states the computed C2/C3 data concern the canonical normal
   presentation, with the source-correspondence scope recorded in the
   [closed-normal specification](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3).
   [proved/agent: by construction]
3. On the empty state the result has `n1 == false`, `criteria == (true,true,true)`
   and `method_circuit_lengths == (0,)`, as specified for the MST extension.
   [proved/agent: by construction]

# Not guaranteed
The C1 slot is not a newly measured smoothing count. Passing criteria alone
does not imply N1 or create a `ClosedNormalOGraph`. This is not a formal proof
of the general C3 correspondence, a framing certificate or a manifold identifier.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4 — C2 and C3; the representation-enforced conditions and empty-state extension are owned by the repository specification.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> w = closedness_witness(ograph_state(vl_example("s3_abalone"))); (w.criteria, w.method_circuit_lengths)
((true, true, true), (5, 1))
```

See also [`certify_closed_ograph`](@ref), [`ClosednessWitness`](@ref).

Role: condition
"""
function closedness_witness(state::OGraphState)
    canonical = canonicalize_ograph(state)
    vertices = true_vertex_count(canonical)
    if vertices == 0
        # MST extension.  `disk_table` has no vertex-local nodes on this input.
        c2_count = 1
        c3_count = 1
        method_circuit_lengths = (0,)
    else
        link = _ograph_virtual_link(canonical)
        table = disk_table(link)
        c2_count = _ograph_trivalent_component_count(link)
        c3_count = size(table.delta, 1)
        method_circuit_lengths = Tuple(Int.(table.path_lengths))
    end
    # Derived from the single-word representation.  Persisted under
    # `virtual-link-c1-c3/2`.  Not a smoothing count.  C1 is enforced
    # by `ograph_state`.  On the empty word this integer is the MST contract.
    c1_word_count = 1
    criteria = (true, c2_count == 1, c3_count == vertices + 1)
    ClosednessWitness(
        ograph_canonical_hash(canonical),
        OGRAPH_CLOSEDNESS_VERIFIER_VERSION,
        criteria,
        c1_word_count,
        c2_count,
        c3_count,
        vertices,
        method_circuit_lengths,
    )
end

function _closedness_failure(witness::ClosednessWitness)
    failed = String[]
    witness.criteria[2] || push!(failed, "C2")
    witness.criteria[3] || push!(failed, "C3")
    "presentation fails $(join(failed, ", ")) [BP97, §1.4] under $(witness.verifier_version): $(witness.criteria)"
end

"""
    ClosedNormalOGraph(state::OGraphState, witness::ClosednessWitness) -> ClosedNormalOGraph

Return a canonical normal presentation accepted by the closed-normal checks.

# Arguments
`state` and `witness` must meet [N1–N2 and C1–C3](../specs/ograph-terminology.md#enforcement).
The constructor checks N1, recomputes the witness on the canonical state and
requires exact witness equality and passing criteria. Failure throws
`ArgumentError` identifying N1, a stale or mismatching witness, or failed C2/C3.
Errors from canonicalization or criterion computation propagate.

# Returns
A `ClosedNormalOGraph` with canonical `state::OGraphState` and the matching
`witness::ClosednessWitness`.

# Guarantees
1. The returned state is canonical, has at least one true vertex, and its
   witness equals `closedness_witness(state)` with all criteria true.
   Inputs are unchanged. [proved/agent: by construction]
2. The intended geometric object is the closed normal o-graph and associated
   branched spine of the [closed-normal specification](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3),
   within that specification's correspondence dependencies.
   [proved/agent: by construction]

# Not guaranteed
This is a computational certificate, not a formal proof of the C3 geometric
correspondence. It provides no framing, manifold recognition or group
identification. The MST crossing-free circuit is outside this type.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4 — the closed-normal conditions and their represented spine.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> s = ograph_state(vl_example("s3_abalone")); g = ClosedNormalOGraph(s, closedness_witness(s));

julia> (g.witness.criteria, g.witness.method_circuit_lengths)
((true, true, true), (5, 1))
```

See also [`certify_closed_ograph`](@ref), [`ograph_pi1_presentation`](@ref).

Role: convert
"""
struct ClosedNormalOGraph
    state::OGraphState
    witness::ClosednessWitness

    function ClosedNormalOGraph(state::OGraphState, witness::ClosednessWitness)
        canonical = canonicalize_ograph(state)
        true_vertex_count(canonical) >= 1 || throw(ArgumentError(
            "N1: a closed normal o-graph has at least one vertex [BP97, §1.4]"))
        expected = closedness_witness(canonical)
        witness == expected || throw(ArgumentError(
            "closedness witness is stale or does not match the presentation"))
        all(witness.criteria) || throw(ArgumentError(_closedness_failure(witness)))
        new(canonical, witness)
    end
end

"""
    certify_closed_ograph(state::OGraphState) -> ClosedNormalOGraph

Certify an immutable presentation with the current closed-normal o-graph checks.

# Arguments
`state` must satisfy [N1–N2 and C1–C3](../specs/ograph-terminology.md#enforcement).
Its representation enforces N2/C1; this call checks the remaining domain and
computed criteria. Empty input or failed criteria throw `ArgumentError`
identifying N1 or failed C2/C3, with the verifier diagnostic. Criterion
computation errors propagate.

# Returns
A `ClosedNormalOGraph` with canonical state and a freshly checked matching
closedness witness.

# Guarantees
1. The result satisfies the `ClosedNormalOGraph` constructor contract and
   retains the input's canonical presentation identity. The input is unchanged.
   [proved/agent: by construction]
2. The intended object is the closed normal presentation of the
   [branched spine](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3),
   subject to the specification's geometric correspondence dependencies.
   [proved/agent: by construction]

# Not guaranteed
Certification is computational, not a formal proof of the C3 correspondence.
It gives no framing, manifold identification or raw-edge-group identification.
The zero-vertex MST circuit uses a different entry point.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4 — the named conditions of the closed-normal class.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> certify_closed_ograph(ograph_state(vl_example("s3_abalone"))).witness.criteria
(true, true, true)
```

See also [`ClosedNormalOGraph`](@ref), [`ograph_pi1_presentation`](@ref).

Role: convert
"""
function certify_closed_ograph(state::OGraphState)
    canonical = canonicalize_ograph(state)
    ClosedNormalOGraph(canonical, closedness_witness(canonical))
end

"""MST/MATLAB crossing-free one-circuit state, outside [BP97, §1.4].

One oriented circuit, no vertex.  The fundamental group is infinite cyclic
and the first homology is Z.  This is not a closed normal o-graph.

cf: lib/topology.jl/src/ograph_invariants.jl — MATLAB crossing-free circuit; source MST, owner decision 2026-11-22
"""
struct MSTCrossingFreeCircuit
    state::OGraphState
    witness::ClosednessWitness

    function MSTCrossingFreeCircuit(state::OGraphState, witness::ClosednessWitness)
        canonical = canonicalize_ograph(state)
        true_vertex_count(canonical) == 0 || throw(ArgumentError(
            "MST crossing-free circuit has no vertices"))
        witness == closedness_witness(canonical) || throw(ArgumentError(
            "MST crossing-free witness does not match the empty presentation"))
        witness.n1 && throw(ArgumentError(
            "MST crossing-free circuit is outside N1 [BP97, §1.4]"))
        new(canonical, witness)
    end
end

function mst_crossing_free_circuit()
    state = OGraphState(Int[], Int[])
    MSTCrossingFreeCircuit(state, closedness_witness(state))
end

"""
    certify_closed_ograph(link::VirtualLink; convention_id=OGRAPH_DEFAULT_CONVENTION) -> ClosedNormalOGraph

Certify the real-Gauss normal presentation of an unweighted closed link input.

# Arguments
`link` and optional `convention_id` must meet [`ograph_state`](@ref)'s checked
domain, and the resulting state must pass [`certify_closed_ograph`](@ref).
Cuts, weighted or Z2 input, multiple real-Gauss components, invalid convention
or representation data, and failed N1/C2/C3 checks raise the delegated
`ArgumentError` diagnostics. Other errors propagate; unknown keywords are
forwarded to `ograph_state` and rejected there.

# Returns
A `ClosedNormalOGraph` for the canonical real-Gauss state with its freshly
checked witness.

# Guarantees
1. The result equals `certify_closed_ograph(ograph_state(link; convention_id))`
   in its state and witness fields. The link is unchanged.
   [proved/agent: by construction]
2. The geometric meaning is that of the accepted closed-normal presentation
   under [Enforcement](../specs/ograph-terminology.md#enforcement), retaining
   its source-correspondence dependencies. [proved/agent: by construction]

# Not guaranteed
The raw virtual-crossing layout is not retained in the real projection.
Certification does not identify a manifold or construct a framing, and is
not a formal proof of the C3 correspondence.

# Sources
cf:[[BP97#^bp97-sec-0536]] §1.4 — closed-normal conditions for the projected presentation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> g = certify_closed_ograph(vl_example("s3_abalone")); (g.witness.true_vertices, g.witness.c3_method_circuits)
(1, 2)
```

See also [`ograph_state`](@ref), [`ograph_pi1_presentation`](@ref), [`ograph_h1`](@ref).

Role: convert
"""
function certify_closed_ograph(link::VirtualLink; kwargs...)
    certify_closed_ograph(ograph_state(link; kwargs...))
end
