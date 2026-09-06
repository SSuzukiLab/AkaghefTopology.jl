"""Immutable, one-circuit presentation state for an unweighted normal o-graph.

This value records representation-level data only.  Raw equality and `hash`
compare the stored root and labels exactly.  `canonicalize_ograph` and
`ograph_canonical_hash` identify cyclic re-rooting and crossing relabeling; they
do not identify MP-, PS-, CP-, mirror-, or orientation-reversal equivalence.
"""

using SHA

const OGRAPH_PRESENTATION_SCHEMA = "akaghef.ograph/1"
const OGRAPH_PRESENTATION_CANON_VERSION = 1
const OGRAPH_PRESENTATION_CANONICALIZATION = "ograph-presentation/1"
const OGRAPH_DEFAULT_CONVENTION = "virtual-link-rgauss/1"
const OGRAPH_CLOSEDNESS_VERIFIER_VERSION = "virtual-link-c1-c3/1"

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
        iseven(length(code)) || throw(ArgumentError("o-graph Gauss code must have even length"))
        all(!=(0), code) || throw(ArgumentError("o-graph real Gauss labels cannot be zero"))
        all(sign -> sign in (-1, 1), signs) || throw(ArgumentError(
            "o-graph vertex signs must be +1 or -1"))
        length(signs) == length(code) ÷ 2 || throw(ArgumentError(
            "vertex_signs must contain one sign per true vertex"))
        for label in eachindex(signs)
            occurrences = filter(value -> abs(value) == label, code)
            length(occurrences) == 2 || throw(ArgumentError(
                "true vertex $label must occur exactly twice"))
            occurrences[1] * occurrences[2] < 0 || throw(ArgumentError(
                "true vertex $label must occur once positive and once negative"))
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

"""Create an immutable unweighted state from a fresh real-Gauss view."""
function ograph_state(link::VirtualLink;
                      convention_id::AbstractString=OGRAPH_DEFAULT_CONVENTION)
    link.is_weighted && throw(ArgumentError(
        "weighted VirtualLink input must not be silently projected to an unweighted o-graph"))
    code, signs = real_gauss_code(link)
    length(code) == 1 || throw(ArgumentError(
        "closed normal o-graph state requires one oriented Gauss circuit"))
    OGraphState(only(code), signs; convention_id=convention_id)
end

struct ClosednessWitness
    canonical_hash::String
    verifier_version::String
    criteria::NTuple{3,Bool}
    c1_oriented_circuits::Int
    c2_trivalent_components::Int
    c3_method_circuits::Int
    true_vertices::Int
    method_circuit_lengths::Tuple{Vararg{Int}}
end

"""Construct the standard `Graphs.jl` object for BP Figure 1.2 (C2).

The vertex set is the local trivalent-node set of the replacement graph.  Each
PD edge occurs twice; the three local nodes at its two occurrences are joined
in order.  The returned object is the semantic graph used for C2.  Conversion
from the o-graph/PD representation is local glue; graph semantics are owned by
`Graphs.jl`, not by a private graph implementation.
"""
function ograph_graph_a(link::VirtualLink)::SimpleGraph
    _, orientation = real_gauss_code(link)
    vertex_count = length(orientation)
    vertex_count == 0 && return SimpleGraph(0)

    pd = pd_code(link)
    order = [4, 1, 2, 3, 4, 1]
    graph_edges = Tuple{Int, Int}[]
    used = Set{Int}()

    for edge in unique(vec(pd))
        occurrences = findall(==(edge), pd)
        length(occurrences) == 2 || throw(ArgumentError(
            "PD edge $edge must occur exactly twice"))

        endpoints = Vector{Vector{Int}}()
        for occurrence in occurrences
            crossing, slot = occurrence[1], occurrence[2]
            local_nodes = 3 * (crossing - 1) .+ order[slot:slot + 2]
            push!(endpoints, local_nodes)
            union!(used, local_nodes)
        end

        for index in 1:3
            push!(graph_edges, (endpoints[1][index], endpoints[2][index]))
        end
    end

    graph = SimpleGraph(maximum(used))
    for (left, right) in graph_edges
        add_edge!(graph, left, right)
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
    # representation so the crossing-free one-circuit state remains a valid
    # closed o-graph input rather than failing in `_normalize_gauss`.
    set_data!(link; rgauss=[Int[state.gauss...]],
              orientation=Int[state.vertex_signs...])
end

"""Recompute a finite, auditable C1-C3 implementation witness from immutable data.

This is an operational replay under `OGRAPH_CLOSEDNESS_VERIFIER_VERSION`, not a
standalone proof that the implementation exactly realizes the source theorem.
"""
function closedness_witness(state::OGraphState)
    canonical = canonicalize_ograph(state)
    link = _ograph_virtual_link(canonical)
    table = disk_table(link)
    c1_count = 1
    c2_count = _ograph_trivalent_component_count(link)
    vertices = true_vertex_count(canonical)
    # The crossing-free one-circuit o-graph has one Method circuit although
    # `disk_table` has no vertex-local nodes from which to reconstruct it.
    # Record that boundary case explicitly (with zero local-node length).
    c3_count = vertices == 0 ? 1 : size(table.delta, 1)
    method_circuit_lengths = vertices == 0 ? (0,) : Tuple(Int.(table.path_lengths))
    criteria = (c1_count == 1, c2_count == 1, c3_count == vertices + 1)
    ClosednessWitness(
        ograph_canonical_hash(canonical),
        OGRAPH_CLOSEDNESS_VERIFIER_VERSION,
        criteria,
        c1_count,
        c2_count,
        c3_count,
        vertices,
        method_circuit_lengths,
    )
end

struct ClosedNormalOGraph
    state::OGraphState
    witness::ClosednessWitness

    function ClosedNormalOGraph(state::OGraphState, witness::ClosednessWitness)
        canonical = canonicalize_ograph(state)
        expected = closedness_witness(canonical)
        witness == expected || throw(ArgumentError(
            "closedness witness is stale or does not match the presentation"))
        all(witness.criteria) || throw(ArgumentError(
            "presentation is not a closed normal o-graph under $(witness.verifier_version): $(witness.criteria)"))
        new(canonical, witness)
    end
end

"""Return an operationally C1-C3-validated wrapper, rejecting any failed check."""
function certify_closed_ograph(state::OGraphState)
    canonical = canonicalize_ograph(state)
    ClosedNormalOGraph(canonical, closedness_witness(canonical))
end

function certify_closed_ograph(link::VirtualLink; kwargs...)
    certify_closed_ograph(ograph_state(link; kwargs...))
end
