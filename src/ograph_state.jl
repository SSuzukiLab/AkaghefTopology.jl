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

"""Canonicalize only cyclic re-rooting and crossing labels of a presentation."""
function canonicalize_ograph(state::OGraphState)
    isempty(state.gauss) && return copy(state)
    best = _ograph_relabel_from_root(state, 1)
    for root in 2:length(state.gauss)
        candidate = _ograph_relabel_from_root(state, root)
        _ograph_state_key(candidate) < _ograph_state_key(best) && (best = candidate)
    end
    best
end

"""Deterministic, domain-qualified serialization for o-graph DB identity."""
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

"""Create an immutable unweighted state from a fresh real-Gauss view."""
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

"""C1–C3 record for one `OGraphState`, plus the N1 and N2 marks.

`criteria` is `(C1 slot, C2, C3)`.  The C1 slot is true because the state
holds one cyclic word.  It is not a count of smoothed circuits.  C1's
enforcement site is `c1_enforcement` (`ograph_state`).  N2's site is
`n2_enforcement` (`OGraphState`).  `n1` is the computed [BP97, §1.4] vertex
count check, true when `true_vertices >= 1`.

`c1_oriented_circuits` is the integer derived from the single-word
representation, persisted under `virtual-link-c1-c3/2`.  Stored evidence
compares that integer exactly.  It is not a smoothing count and not the C1
test.  On the empty word the same integer is the MST contract.  C2 and C3
are the computed fields.

cf:[[BP97]] §1.4 — N1, N2, C1, C2, C3; the zero-vertex branch is the MST extension
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

"""Duplicate-preserving BP Figure 1.2 graph and GC-edge provenance.

`origins[i] == (component, tail_occurrence, lane)` describes `edges[i]`.
`lane` is the source-triple position before the pair permutation is applied.
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

"""Construct the BP Figure 1.2 C2 multigraph directly from a complete GC.

The returned edge vector retains loops, parallel edges, and row-aligned origins.
The construction does not pass through PD.  See `docs/specs/ograph-c2.md`.
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

c2_graph(state::OGraphState) =
    c2_graph([collect(state.gauss)], collect(state.vertex_signs))

"""Construct C2 from the underlying o-graph of a `VirtualLink`.

Edge weights do not occur in BP condition C2.  A weighted input is therefore
read without mutation and with its GC/sign data unchanged; this function does
not claim or perform any weighted-move equivalence.
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

"""Lossy `SimpleGraph` projection retained for legacy connectivity callers.

The authoritative C2 representation is `c2_graph`; this projection drops loops
and parallel-edge multiplicity and must not be stored as the C2 graph.
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

"""Recompute the closedness record for one presentation.

C2 is the component count of the [BP97, §1.4] Figure 1.2 graph.  C3, when
there is at least one vertex, is the row count of `disk_table` (Method B;
the spec identifies this count with Figure 1.3).  C1 is not measured here.

The zero-vertex branch is the MST crossing-free circuit: one declared
circuit, length 0, and `n1 == false`.  It is not a Figure 1.3 count.
`ClosedNormalOGraph` refuses that branch under N1.  Callers that only read
`criteria`, including the move verifier on a zero-vertex raw target, still
see `(true, true, true)`.

cf:[[BP97]] §1.4 — C2 and C3; N1 and C1 are recorded, not counted here
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

"""A closed normal o-graph in the sense of [BP97, §1.4].

N1 requires at least one vertex.  N2 and C1 hold by construction of the
state.  C2 and C3 are the computed witness bits.

cf:[[BP97]] §1.4 — class of closed normal o-graphs
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

"""Certify a presentation as a closed normal o-graph.

The empty Gauss code is refused.  N1 ([BP97, §1.4]) requires at least one
vertex.  The MST crossing-free circuit is `mst_crossing_free_circuit`.

cf:[[BP97]] §1.4 — N1, N2, C1, C2, C3
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

function certify_closed_ograph(link::VirtualLink; kwargs...)
    certify_closed_ograph(ograph_state(link; kwargs...))
end
