"""Derived π₁-presentation and exact H₁ cache for closed normal o-graphs.

The raw edge presentation follows the closed, ordinary-strand branch of
MATLAB `VirtualLink.calcSpec("pi1")` and Koda's relations `g=l`, `hg=k`.
Presentation identity is not group identity.  H₁ is a deterministic derived
invariant of that presentation; optional group backends remain a separate,
currently unwritable enrichment layer.
"""

using AbstractAlgebra
using JSON
using SHA
# The store is a child of the caller when a script includes it first, and a
# sibling once this file is included from VirtualLinks.
if isdefined(@__MODULE__, :MoveGraphStore) &&
        getfield(@__MODULE__, :MoveGraphStore) isa Module &&
        parentmodule(getfield(@__MODULE__, :MoveGraphStore)) === @__MODULE__
    @eval using .MoveGraphStore
else
    @eval using ..MoveGraphStore
end

const OGRAPH_PI1_PRESENTATION_SCHEMA = "akaghef.ograph-pi1-presentation/1"
const OGRAPH_PI1_ALGORITHM = "koda-edge-presentation"
const OGRAPH_PI1_ALGORITHM_VERSION = 1
const OGRAPH_PI1_CONVENTION = "pd-code-sage-matlab-closed-ordinary/1"
const OGRAPH_H1_SCHEMA = "akaghef.ograph-h1/1"
const OGRAPH_H1_ALGORITHM = "exponent-sum-smith-normal-form"
const OGRAPH_H1_ALGORITHM_VERSION = 1
const OGRAPH_INVARIANT_STORE_SCHEMA_VERSION = 1

"""
    OGraphPi1Presentation(presentation_id, state_id, presentation_hash, generator_count, relators, algorithm_id, algorithm_version, convention_id, serialization, evidence_json) -> OGraphPi1Presentation

Store a finite group presentation and its supplied identity and evidence fields.

# Arguments
`presentation_id`, `state_id`, `presentation_hash`, `algorithm_id`,
`convention_id`, `serialization` and `evidence_json` are `String` fields;
`generator_count` and `algorithm_version` are `Int`; `relators` is a `Tuple`.
Values convertible to these types are accepted and conversion errors propagate.
This field constructor does not validate the contents. Use the three-argument
constructor for checked signed words and derived identity fields.

# Returns
`OGraphPi1Presentation` with the supplied generator count and ordered relators,
state and presentation identifiers, presentation digest, algorithm name and
version, convention identifier, serialized presentation and JSON evidence.

# Guarantees
1. Every field contains the corresponding supplied value after conversion.
   [proved/agent: by construction]
2. When the fields come from [`ograph_pi1_presentation`](@ref), they describe
   its raw edge presentation. The intended space and the conditional
   spine interpretation are those of the
   [invariant specification](../specs/ograph-colour-invariants.md#closed-examples).
   Direct construction does not establish that interpretation.
   [proved/agent: by construction]

# Not guaranteed
This constructor does not check relator letters, counts, JSON, hashes or
agreement of the metadata with the relators. An arbitrary value need not be
usable by `ograph_h1`. Presentation identity is not abstract-group identity,
and the general raw-edge-presentation identification with the spine's
fundamental group remains unestablished.

# Sources
cf:[[MST25#^mst25-blk-0059]] Closed normal o-graphs and combed 3-manifolds — the represented branched polyhedron; this does not identify arbitrary stored relators with its fundamental group.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> p = ograph_pi1_presentation(certify_closed_ograph(vl_example("s3_abalone")));

julia> q = OGraphPi1Presentation(p.presentation_id, p.state_id, p.presentation_hash, p.generator_count, p.relators, p.algorithm_id, p.algorithm_version, p.convention_id, p.serialization, p.evidence_json);

julia> (q.generator_count, q.relators)
(2, ((1, -2), (2, 1, -1)))
```

See also [`ograph_pi1_presentation`](@ref), [`ograph_h1`](@ref).

Role: construct
"""
struct OGraphPi1Presentation
    presentation_id::String
    state_id::String
    presentation_hash::String
    generator_count::Int
    relators::Tuple
    algorithm_id::String
    algorithm_version::Int
    convention_id::String
    serialization::String
    evidence_json::String
end

_ograph_pi1_value_tuple(value::OGraphPi1Presentation) = (
    value.presentation_id, value.state_id, value.presentation_hash,
    value.generator_count, value.relators, value.algorithm_id,
    value.algorithm_version, value.convention_id, value.serialization,
    value.evidence_json,
)

Base.:(==)(left::OGraphPi1Presentation, right::OGraphPi1Presentation) =
    _ograph_pi1_value_tuple(left) == _ograph_pi1_value_tuple(right)
Base.isequal(left::OGraphPi1Presentation, right::OGraphPi1Presentation) =
    isequal(_ograph_pi1_value_tuple(left), _ograph_pi1_value_tuple(right))
Base.hash(value::OGraphPi1Presentation, seed::UInt) =
    hash(_ograph_pi1_value_tuple(value), seed)

function _normalize_signed_relators(generator_count::Integer, relators)
    count = Int(generator_count)
    count >= 0 || throw(ArgumentError("generator_count must be nonnegative"))
    normalized = Tuple(Tuple(begin
        letter isa Integer || throw(ArgumentError(
            "presentation letters must be signed integers"))
        value = big(letter)
        value != 0 || throw(ArgumentError("presentation letters cannot be zero"))
        abs(value) <= count || throw(ArgumentError(
            "presentation letter refers to a missing generator"))
        Int(value)
    end for letter in relator) for relator in relators)
    count, normalized
end

"""
    OGraphPi1Presentation(state_id::AbstractString, generator_count::Integer, relators; algorithm_id=OGRAPH_PI1_ALGORITHM, algorithm_version=OGRAPH_PI1_ALGORITHM_VERSION, convention_id=OGRAPH_PI1_CONVENTION, evidence=Dict{String,Any}()) -> OGraphPi1Presentation

Construct a checked finite signed-word group presentation with deterministic identity metadata.

# Arguments
`state_id`, `algorithm_id` and `convention_id` must be nonempty after stripping
whitespace. `generator_count` must convert to a nonnegative `Int`, and
`algorithm_version` to a positive `Int`. `relators` is an iterable of iterable
signed integer words accepted by this presentation type. Invalid counts,
zero or noninteger letters, out-of-range generator references, empty identifiers
and nonpositive versions throw `ArgumentError` naming the violation; failed
integer conversions propagate their conversion errors.
`evidence` must be supported by [`ograph_canonical_json`](@ref); its errors
propagate. A state identifier need not name a stored or certified graph.

# Returns
`OGraphPi1Presentation` with stripped identifiers, integer count and version,
tuple relators preserving the supplied word order, and derived `serialization`,
`presentation_hash` and `presentation_id`. `evidence_json` stores the supplied
evidence. The remaining field meanings are those of the type's field constructor.

# Guarantees
1. Returned words contain nonzero `Int` letters whose absolute values do not
   exceed `generator_count`; the constructor preserves letters and word order.
   [proved/agent: by construction]
2. Repeating a call with equal normalized data yields equal identity fields.
   Changing only `evidence` leaves the presentation identity unchanged.
   [proved/agent: by construction]
3. The relators specify an abstract finite presentation to which
   [`ograph_h1`](@ref) applies. The constructor makes no assertion that this
   group is the fundamental group of a given spine or 3-manifold.
   [proved/agent: by construction]

# Not guaranteed
No free reduction, relator simplification, abstract-group canonicalization,
state lookup or geometric certification is performed. The algorithm name
and convention are supplied metadata, not evidence of how the words arose.
For the closed-normal pipeline the
[space identification](../specs/ograph-colour-invariants.md#closed-examples)
retains the qualification in `ograph_pi1_presentation`.

# Examples
```jldoctest
julia> p = OGraphPi1Presentation("cyclic", 1, [fill(1, 12)]);

julia> h = ograph_h1(p);

julia> (h.free_rank, h.torsion_invariant_factors, h.h1_finite_order)
(0, (12,), 12)
```

See also [`ograph_pi1_presentation`](@ref), [`ograph_h1`](@ref).

Role: construct
"""
function OGraphPi1Presentation(state_id::AbstractString,
        generator_count::Integer, relators;
        algorithm_id::AbstractString=OGRAPH_PI1_ALGORITHM,
        algorithm_version::Integer=OGRAPH_PI1_ALGORITHM_VERSION,
        convention_id::AbstractString=OGRAPH_PI1_CONVENTION,
        evidence=Dict{String,Any}())
    state = strip(String(state_id))
    isempty(state) && throw(ArgumentError("state_id must not be empty"))
    count, words = _normalize_signed_relators(generator_count, relators)
    algorithm = strip(String(algorithm_id))
    convention = strip(String(convention_id))
    version = Int(algorithm_version)
    isempty(algorithm) && throw(ArgumentError("algorithm_id must not be empty"))
    isempty(convention) && throw(ArgumentError("convention_id must not be empty"))
    version > 0 || throw(ArgumentError("algorithm_version must be positive"))
    serialization = _ograph_canonical_json(Dict(
        "schema" => OGRAPH_PI1_PRESENTATION_SCHEMA,
        "state_id" => state,
        "algorithm_id" => algorithm,
        "algorithm_version" => version,
        "convention_id" => convention,
        "generator_count" => count,
        "ordered_relators" => [collect(word) for word in words],
    ))
    digest = _ograph_sha256(serialization)
    OGraphPi1Presentation("ograph-pi1-" * digest, state, digest, count, words,
        algorithm, version, convention, serialization,
        _ograph_canonical_json(evidence))
end

"""
    ograph_pi1_presentation(graph::ClosedNormalOGraph; algorithm_id=OGRAPH_PI1_ALGORITHM, algorithm_version=OGRAPH_PI1_ALGORITHM_VERSION, convention_id=OGRAPH_PI1_CONVENTION) -> OGraphPi1Presentation

Return the deterministic Koda/MATLAB raw edge presentation associated with the closed normal o-graph.

# Arguments
- `graph::ClosedNormalOGraph`: obtain this value with
  [`certify_closed_ograph`](@ref) applied to a `VirtualLink` or `OGraphState`.
  Requires N1, N2, C1, C2 and C3 in the
  [closed normal conditions](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3).
  Certification rejects a failed condition with `ArgumentError` naming it;
  this function rechecks the stored state. A raw `VirtualLink` or `OGraphState`
  is not an argument of this method.
  For the `VirtualLink` certification route, the link must be unweighted,
  have `weight_mode(link) == :integer`, and have no trailing cut-boundary
  markers; these exclusions throw `ArgumentError` identifying the unsupported
  input. See [cut markers](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
- `algorithm_id::AbstractString`, `algorithm_version::Integer`, and
  `convention_id::AbstractString`: metadata must equal the displayed defaults.
  Unsupported values throw `ArgumentError` with a message identifying
  unsupported presentation metadata. A version not representable as `Int`
  can throw `InexactError` during conversion.

# Returns
`OGraphPi1Presentation` with these fields:
- `generator_count::Int`: number of generators, indexed from 1.
- `relators::Tuple`: ordered tuple of signed-integer words; a positive index
  denotes a generator, a negative index its inverse, and each word equals 1.
- `state_id::String`: identity of the canonical input presentation.
- `presentation_hash::String`: SHA-256 digest of `serialization`.
- `presentation_id::String`: `"ograph-pi1-"` followed by `presentation_hash`.
- `algorithm_id::String`, `algorithm_version::Int`, `convention_id::String`:
  the supported algorithm and convention metadata.
- `serialization::String`: canonical JSON of the schema, state identity,
  algorithm metadata, generator count and ordered relators.
- `evidence_json::String`: JSON containing the PD rows, slot names, generator
  rule, relation order and source-rule description.

# Guarantees
1. If the graph has `n` true vertices, `generator_count == 2n` and
   `length(relators) == 2n`. The first `n` words have length 2 and the last
   `n` words have length 3; every letter has absolute value between 1 and
   `generator_count`. [proved/agent: by construction]
2. Repeated calls on the same certified state with the same supported
   metadata return equal `OGraphPi1Presentation` values. Cyclic re-rooting
   and crossing relabelling identified by [`canonicalize_ograph`](@ref)
   give the same presentation. [proved/agent: by construction]
3. [`ograph_h1`](@ref) accepts the result and returns the integral
   abelianization of the group defined by its generators and relators in the
   [closed-normal invariant pipeline](../specs/ograph-colour-invariants.md#closed-examples).
   [proved/agent: by construction]
4. The input `graph` is unchanged. [proved/agent: by construction]
5. The intended space is the closed oriented 3-manifold obtained by capping
   the spherical boundary of the thickening of the branched spine represented
   by `graph`, as in the [spine/manifold interpretation](../specs/ograph-colour-invariants.md#closed-examples).
   The spine, its thickening and that closed manifold have isomorphic
   fundamental groups. Thus, **conditional on identifying the returned raw
   edge presentation with the fundamental group of that spine**, the result
   presents the fundamental group of that closed manifold, and
   [`ograph_h1`](@ref) gives its first integral homology.
   [open: identification of the raw edge presentation with the spine's fundamental group]

# Conventions
- Passage sign, vertex sign and local ports:
  [o-graph terminology](../specs/ograph-terminology.md#4-passage-sign-and-vertex-sign)
  and [vertex conventions, rules](../conventions/OGraphVertex.md#規則).
- The algorithm metadata select `pd-code-sage-matlab-closed-ordinary/1`.
  The input is the canonical presentation, rather than the received labels.

# Not guaranteed
The result is not a reduced presentation or a canonical identifier of the
abstract group or manifold. Equality of presentation identities is stronger
than group isomorphism. Weights and framing are not part of this invariant.
The MST crossing-free circuit is outside N1 and uses a separate overload.
The general identification of these particular raw edge relators with the
spine's fundamental group is not established by the cited sources and current
specification. The specification records agreement of first homology on its
listed closed examples; that finite comparison does not prove the general
identification of fundamental groups. N1--C3 certification checks the input
conditions, not this identification.

# Sources
cf:[[MST25#^mst25-def-0042]] §2.1 — a spine is a collapse of its thickening.
cf:[[MST25#^mst25-blk-0057]] Closed normal o-graphs and combed 3-manifolds — the closed manifold is obtained by capping the spherical boundary of the thickening.
cf:[[MST25#^mst25-blk-0059]] Closed normal o-graphs and combed 3-manifolds — a closed normal o-graph represents a closed branched polyhedron; these citations identify the space, not the raw edge presentation theorem.

# Examples
```jldoctest
julia> graph = certify_closed_ograph(lens_space_link(3, 1));

julia> pres = ograph_pi1_presentation(graph);

julia> (pres.generator_count, length(pres.relators))
(6, 6)

julia> sort(unique(length.(pres.relators)))
2-element Vector{Int64}:
 2
 3
```

See also [`certify_closed_ograph`](@ref), [`ograph_h1`](@ref).

Role: invariant
"""
function ograph_pi1_presentation(graph::ClosedNormalOGraph;
        algorithm_id::AbstractString=OGRAPH_PI1_ALGORITHM,
        algorithm_version::Integer=OGRAPH_PI1_ALGORITHM_VERSION,
        convention_id::AbstractString=OGRAPH_PI1_CONVENTION)
    String(algorithm_id) == OGRAPH_PI1_ALGORITHM &&
        Int(algorithm_version) == OGRAPH_PI1_ALGORITHM_VERSION &&
        String(convention_id) == OGRAPH_PI1_CONVENTION || throw(ArgumentError(
            "unsupported o-graph pi1 presentation algorithm metadata"))
    canonical = certify_closed_ograph(graph.state)
    state = ograph_canonical_state(canonical)
    pd = pd_code(_ograph_virtual_link(canonical.state))
    # MATLAB uses `max(RPD) + Ncircle`: the crossing-free oriented circuit has
    # one edge generator and no relations, hence H1 = Z rather than 0.
    generator_count = isempty(pd) ? 1 : maximum(pd)
    isempty(pd) || Set(vec(pd)) == Set(1:generator_count) || throw(ArgumentError(
        "PD edge labels are not contiguous from one"))
    first_relations = Tuple((Int(row[1]), -Int(row[3])) for row in eachrow(pd))
    second_relations = Tuple((Int(row[4]), Int(row[1]), -Int(row[2]))
                             for row in eachrow(pd))
    relators = (first_relations..., second_relations...)
    evidence = Dict(
        "pd_rows" => [Int.(collect(row)) for row in eachrow(pd)],
        "pd_row_slots" => ["g", "k", "l", "h"],
        "generator_rule" =>
            "max PD edge label plus one generator per crossing-free circuit",
        "crossing_free_circuit_count" => isempty(pd) ? 1 : 0,
        "relation_order" => "all g*l^-1, then all h*g*k^-1",
        "matlab_scope" => "calcSpec(pi1), closed ordinary strands only",
        "source_rule" => "Koda/MST: generators are edges; g=l and hg=k at each vertex",
    )
    OGraphPi1Presentation(state.state_id, generator_count, relators;
        algorithm_id=algorithm_id, algorithm_version=algorithm_version,
        convention_id=convention_id, evidence=evidence)
end

"""
    ograph_pi1_presentation(circuit::MSTCrossingFreeCircuit; algorithm_id::AbstractString=OGRAPH_PI1_ALGORITHM, algorithm_version::Integer=OGRAPH_PI1_ALGORITHM_VERSION, convention_id::AbstractString=OGRAPH_PI1_CONVENTION) -> OGraphPi1Presentation

Return the infinite-cyclic group presentation assigned to the MST crossing-free circuit.

# Arguments
- `circuit`: the separately typed
  [MST crossing-free circuit](../specs/ograph-terminology.md#mst-crossing-free-circuit).
- `algorithm_id`, `algorithm_version`, `convention_id`: must equal the
  supported defaults. `ArgumentError` identifies unsupported metadata;
  conversion of the version to `Int` can throw `InexactError`.

# Returns
An `OGraphPi1Presentation` with one generator and no relators, state and
presentation identities, supported algorithm metadata, canonical JSON
serialization and evidence recording the crossing-free rule.

# Guarantees
1. `generator_count == 1` and `isempty(relators)`, so the presented group is infinite cyclic. [proved/agent: by construction]
2. Under the repository's specified identification of this circuit with S²×S¹, the result presents its fundamental group. [proved/agent: docs/specs/ograph-terminology.md#mst-crossing-free-circuit, using the stated identification]
3. The circuit is not modified; repeated calls with the same metadata return equal presentation data. [proved/agent: by construction]

# Not guaranteed
This separate extension is outside the BP97 N1 class. The function does not
derive its S²×S¹ identification from BP97 closed-normal certification.

# Sources
cf:docs/specs/ograph-terminology.md#mst-crossing-free-circuit — the repository's MST extension and its S²×S¹ identification.

# Examples
```jldoctest
julia> p = ograph_pi1_presentation(Topology.mst_crossing_free_circuit());

julia> (p.generator_count, isempty(p.relators))
(1, true)
```

See also [`ograph_h1`](@ref).

Role: invariant
"""
function ograph_pi1_presentation(circuit::MSTCrossingFreeCircuit;
        algorithm_id::AbstractString=OGRAPH_PI1_ALGORITHM,
        algorithm_version::Integer=OGRAPH_PI1_ALGORITHM_VERSION,
        convention_id::AbstractString=OGRAPH_PI1_CONVENTION)
    String(algorithm_id) == OGRAPH_PI1_ALGORITHM &&
        Int(algorithm_version) == OGRAPH_PI1_ALGORITHM_VERSION &&
        String(convention_id) == OGRAPH_PI1_CONVENTION || throw(ArgumentError(
            "unsupported o-graph pi1 presentation algorithm metadata"))
    serialization = ograph_canonical_serialization(circuit.state)
    state_id = bytes2hex(SHA.sha256(codeunits(
        "v$(OGRAPH_PRESENTATION_CANON_VERSION)\0$(serialization)")))
    evidence = Dict(
        "pd_rows" => Vector{Vector{Int}}(),
        "pd_row_slots" => ["g", "k", "l", "h"],
        "generator_rule" =>
            "max PD edge label plus one generator per crossing-free circuit",
        "crossing_free_circuit_count" => 1,
        "relation_order" => "all g*l^-1, then all h*g*k^-1",
        "matlab_scope" => "calcSpec(pi1), closed ordinary strands only",
        "source_rule" => "MST: one generator for the crossing-free circuit; not a closed normal o-graph",
    )
    OGraphPi1Presentation(state_id, 1, ();
        algorithm_id=algorithm_id, algorithm_version=algorithm_version,
        convention_id=convention_id, evidence=evidence)
end

"""
    OGraphH1Result(result_id, presentation_id, free_rank, torsion_invariant_factors, torsion_subgroup_order, h1_finite_order, signature, result_fingerprint, algorithm_id, algorithm_version, backend, backend_version, evidence_json) -> OGraphH1Result

Store integral abelianization data and its supplied computation metadata.

# Arguments
`free_rank` and `algorithm_version` are `Int`, `torsion_invariant_factors`
is a `Tuple`, `torsion_subgroup_order` is `BigInt`, and `h1_finite_order`
is `nothing` or `BigInt`. The remaining arguments are `String` fields.
Convertible values are accepted; conversion errors propagate. No consistency
checks are performed. Obtain a computed value with [`ograph_h1`](@ref).

# Returns
`OGraphH1Result` with the supplied free rank, torsion factors, torsion-subgroup
order and optional finite group order. `result_id`, `presentation_id`,
`result_fingerprint` and `signature` store identity data; `algorithm_id`,
`algorithm_version`, `backend`, `backend_version` and `evidence_json` store
the computation metadata and evidence.

# Guarantees
1. Each field equals the supplied value after conversion to its declared type.
   [proved/agent: by construction]
2. A value computed by `ograph_h1` describes the integral abelianization of
   its input presentation's group. Its interpretation as the first homology
   of a spine and capped thickening has the same conditional scope as that
   operation and the [invariant specification](../specs/ograph-colour-invariants.md#closed-examples).
   [proved/agent: by construction]

# Conventions
Use the invariant-factor convention on the specification page linked above.

# Not guaranteed
The constructor does not recompute a homology group, validate invariant factors
or verify the supplied fingerprints and JSON. Arbitrary fields need not
describe any valid computation or manifold. Equal homology groups need not
have equal result identities.

# Sources
cf:[[BP97#^bp97-sec-0703]] §10.1 — homology of the spine and closed case; the geometric interpretation additionally needs the raw-presentation identification stated by `ograph_h1`.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> h = ograph_h1(ograph_pi1_presentation(certify_closed_ograph(vl_example("s3_abalone"))));

julia> r = OGraphH1Result(h.result_id, h.presentation_id, h.free_rank, h.torsion_invariant_factors, h.torsion_subgroup_order, h.h1_finite_order, h.signature, h.result_fingerprint, h.algorithm_id, h.algorithm_version, h.backend, h.backend_version, h.evidence_json);

julia> (r.free_rank, r.torsion_invariant_factors, r.h1_finite_order)
(0, (), 1)
```

See also [`ograph_h1`](@ref).

Role: construct
"""
struct OGraphH1Result
    result_id::String
    presentation_id::String
    free_rank::Int
    torsion_invariant_factors::Tuple
    torsion_subgroup_order::BigInt
    h1_finite_order::Union{Nothing,BigInt}
    signature::String
    result_fingerprint::String
    algorithm_id::String
    algorithm_version::Int
    backend::String
    backend_version::String
    evidence_json::String
end

_ograph_h1_value_tuple(value::OGraphH1Result) = (
    value.result_id, value.presentation_id, value.free_rank,
    value.torsion_invariant_factors, value.torsion_subgroup_order,
    value.h1_finite_order, value.signature, value.result_fingerprint,
    value.algorithm_id, value.algorithm_version, value.backend,
    value.backend_version, value.evidence_json,
)

Base.:(==)(left::OGraphH1Result, right::OGraphH1Result) =
    _ograph_h1_value_tuple(left) == _ograph_h1_value_tuple(right)
Base.isequal(left::OGraphH1Result, right::OGraphH1Result) =
    isequal(_ograph_h1_value_tuple(left), _ograph_h1_value_tuple(right))
Base.hash(value::OGraphH1Result, seed::UInt) =
    hash(_ograph_h1_value_tuple(value), seed)

function _ograph_exponent_sum_matrix(presentation::OGraphPi1Presentation)
    matrix = zeros(BigInt, length(presentation.relators),
                   presentation.generator_count)
    for (row, relator) in enumerate(presentation.relators), letter in relator
        matrix[row, abs(letter)] += sign(letter)
    end
    matrix
end

function _ograph_h1_from_matrix(presentation_id::AbstractString,
                                exponent_matrix::AbstractMatrix{<:Integer};
        algorithm_id::AbstractString=OGRAPH_H1_ALGORITHM,
        algorithm_version::Integer=OGRAPH_H1_ALGORITHM_VERSION,
        backend::AbstractString="AbstractAlgebra",
        backend_version::AbstractString=string(Base.pkgversion(AbstractAlgebra)))
    matrix = BigInt.(exponent_matrix)
    relation_count, generator_count = size(matrix)
    diagonal = BigInt[]
    if relation_count > 0 && generator_count > 0
        smith = AbstractAlgebra.snf(AbstractAlgebra.matrix(AbstractAlgebra.ZZ, matrix))
        diagonal = [abs(BigInt(smith[index, index]))
                    for index in 1:min(size(smith)...)]
    end
    nonzero = filter(!iszero, diagonal)
    all(value > 0 for value in nonzero) || throw(ArgumentError(
        "Smith diagonal must be nonnegative"))
    all(nonzero[index + 1] % nonzero[index] == 0
        for index in 1:max(length(nonzero) - 1, 0)) || throw(ArgumentError(
            "Smith invariant factors are not divisibility ordered"))
    free_rank = generator_count - length(nonzero)
    free_rank >= 0 || throw(ArgumentError("Smith rank exceeds generator count"))
    factors = Tuple(value for value in nonzero if value > 1)
    torsion_order = prod(factors; init=BigInt(1))
    finite_order = free_rank == 0 ? torsion_order : nothing
    signature = "H1;schema=$(OGRAPH_H1_SCHEMA);free_rank=$(free_rank);" *
        "torsion=[$(join(string.(factors), ','))];finite_order=" *
        (finite_order === nothing ? "null" : string(finite_order))
    algorithm = strip(String(algorithm_id))
    version = Int(algorithm_version)
    backend_name = strip(String(backend))
    backend_release = strip(String(backend_version))
    isempty(algorithm) && throw(ArgumentError("H1 algorithm_id must not be empty"))
    version > 0 || throw(ArgumentError("H1 algorithm_version must be positive"))
    isempty(backend_name) && throw(ArgumentError("H1 backend must not be empty"))
    isempty(backend_release) && throw(ArgumentError("H1 backend_version must not be empty"))
    evidence = _ograph_canonical_json(Dict(
        "exponent_sum_matrix_decimal" =>
            [[string(matrix[row, column]) for column in axes(matrix, 2)]
             for row in axes(matrix, 1)],
        "smith_diagonal_decimal" => string.(diagonal),
        "relation_count" => relation_count,
        "generator_count" => generator_count,
    ))
    result_data = _ograph_canonical_json(Dict(
        "schema" => OGRAPH_H1_SCHEMA,
        "presentation_id" => String(presentation_id),
        "algorithm_id" => algorithm,
        "algorithm_version" => version,
        "backend" => backend_name,
        "backend_version" => backend_release,
        "signature" => signature,
        "evidence" => JSON.parse(evidence),
    ))
    fingerprint = _ograph_sha256(result_data)
    OGraphH1Result("ograph-h1-" * fingerprint, String(presentation_id),
        free_rank, factors, torsion_order, finite_order, signature, fingerprint,
        algorithm, version, backend_name, backend_release, evidence)
end

"""
    ograph_h1(presentation::OGraphPi1Presentation) -> OGraphH1Result

Return the exact integral abelianization of the group defined by the presentation.

# Arguments
- `presentation::OGraphPi1Presentation`: for the first homology of a closed
  normal o-graph, call [`certify_closed_ograph`](@ref), then
  [`ograph_pi1_presentation`](@ref), then this function. The graph must satisfy
  N1, N2, C1, C2 and C3 in the
  [closed normal conditions](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3).
  Graph certification throws `ArgumentError` naming failed conditions.
- A synthetic presentation may instead be made with
  `OGraphPi1Presentation(state_id, generator_count, relators)`. That constructor
  validates the count and signed letters and throws `ArgumentError` for
  negative counts representable as `Int`, zero letters, noninteger letters or
  missing generators. Integer conversions outside `Int` can throw `InexactError`.
  This function trusts the supplied presentation fields; malformed values
  made through the all-field constructor are not validated here and invalid
  generator indices can throw `BoundsError`.

# Returns
`OGraphH1Result` describing `Z^free_rank` plus the cyclic summands whose orders
are `torsion_invariant_factors`, with fields:
- `free_rank::Int`: the rank of the free abelian summand.
- `torsion_invariant_factors::Tuple`: `BigInt` orders greater than 1, in
  divisibility order.
- `torsion_subgroup_order::BigInt`: product of those factors, or 1 when empty.
- `h1_finite_order::Union{Nothing,BigInt}`: the group order when `free_rank == 0`,
  and `nothing` otherwise.
- `presentation_id::String`: the supplied presentation's identity.
- `signature::String`: schema, free rank, torsion factors and finite order.
- `algorithm_id::String`, `algorithm_version::Int`: the H1 algorithm metadata.
- `backend::String`, `backend_version::String`: the exact-arithmetic backend
  name and installed version.
- `evidence_json::String`: JSON with the exponent-sum matrix and Smith
  diagonal as decimal strings, and the relation and generator counts.
- `result_fingerprint::String`: SHA-256 digest of canonical result data,
  including the presentation identity, backend metadata, signature and evidence.
- `result_id::String`: `"ograph-h1-"` followed by `result_fingerprint`.

# Guarantees
1. `free_rank >= 0`; each torsion factor is greater than 1 and divides its
   successor. The factors and free rank describe the quotient of
   `Z^presentation.generator_count` by the subgroup generated by the
   exponent-sum rows of the relators. [proved/agent: by construction]
2. `torsion_subgroup_order == prod(torsion_invariant_factors; init=big(1))`.
   `h1_finite_order == torsion_subgroup_order` exactly when `free_rank == 0`;
   otherwise `h1_finite_order === nothing`. [proved/agent: by construction]
3. `presentation_id == presentation.presentation_id`, and the input is unchanged. [proved/agent: by construction]
4. Repeated calls on the same presentation in the same backend version
   return equal `OGraphH1Result` values. [proved/agent: by construction]
5. For a presentation returned by [`ograph_pi1_presentation`](@ref), the
   intended space is the closed oriented 3-manifold obtained by capping the
   spherical boundary of the thickening of the graph's branched spine;
   see the [spine/manifold interpretation](../specs/ograph-colour-invariants.md#closed-examples).
   **Conditional on that presentation presenting the spine's fundamental
   group**, the returned abelian group is the first integral homology of
   both the spine and this closed manifold. For a valid synthetic
   presentation, the unconditional guarantee is the integral abelianization
   of the group it defines, as described in guarantee 1.
   [open: identification of the raw edge presentation with the spine's fundamental group]

# Conventions
The invariant-factor convention, including omission of unit factors, is
recorded in [o-graph invariant specifications](../specs/ograph-colour-invariants.md#closed-examples).

# Not guaranteed
This function does not certify that a synthetic presentation comes from a
3-manifold or distinguish groups with the same abelianization. The result
identity includes backend version and presentation identity, so equal homology
groups need not have equal result identities. For an infinite group the finite
torsion-subgroup order is not the order of the whole group.
For the closed normal o-graph pipeline, the general identification of the raw
edge presentation with the spine's fundamental group remains unestablished in
the cited sources and specification. The finite homology comparisons listed
there do not establish that general identification. This function computes
abelianization and does not check that missing geometric identification.

# Sources
cf:[[MST25#^mst25-blk-0057]] Closed normal o-graphs and combed 3-manifolds — the space in the closed-manifold interpretation is the capped thickening of the branched spine.
cf:[[BP97#^bp97-sec-0703]] §10.1 — cellular homology of the spine and the extra 3-cell in the closed case; the raw-presentation identification is an additional hypothesis, not a theorem supplied here.

# Examples
```jldoctest
julia> pres = ograph_pi1_presentation(certify_closed_ograph(lens_space_link(3, 1)));

julia> h = ograph_h1(pres);

julia> (h.free_rank, h.torsion_invariant_factors)
(0, (3,))

julia> h.h1_finite_order
3
```

See also [`ograph_pi1_presentation`](@ref), [`certify_closed_ograph`](@ref).

Role: invariant
"""
function ograph_h1(presentation::OGraphPi1Presentation)
    _ograph_h1_from_matrix(presentation.presentation_id,
                           _ograph_exponent_sum_matrix(presentation))
end

function _ograph_invariant_schema_path()
    normpath(joinpath(@__DIR__, "..", "db", "ograph_invariant_schema.sql"))
end

function _ograph_invariant_existing_schema_version(db)
    meta = MoveGraphStore._one(db, """
        SELECT 1 AS present FROM sqlite_master
        WHERE type='table' AND name='ograph_invariant_schema_meta'
    """)
    if meta === nothing
        names = ("ograph_pi1_presentations", "ograph_h1_results",
                 "ograph_group_backend_runs", "ograph_finite_quotient_fingerprints",
                 "ograph_manifold_identifications")
        placeholders = join(fill("?", length(names)), ",")
        partial = MoveGraphStore._rows(db,
            "SELECT name FROM sqlite_master WHERE type='table' AND name IN ($placeholders)",
            names)
        isempty(partial) || throw(ArgumentError(
            "o-graph invariant tables exist without invariant schema metadata"))
        return nothing
    end
    row = MoveGraphStore._one(db, """
        SELECT value FROM ograph_invariant_schema_meta WHERE key='schema_version'
    """)
    row === nothing && throw(ArgumentError(
        "ograph_invariant_schema_meta has no schema_version"))
    try
        parse(Int, String(row.value))
    catch
        throw(ArgumentError("o-graph invariant schema_version is not an integer"))
    end
end

function _assert_ograph_invariant_schema(store::MoveGraphStoreDB)
    _ograph_invariant_existing_schema_version(store.db) ==
        OGRAPH_INVARIANT_STORE_SCHEMA_VERSION || throw(ArgumentError(
            "o-graph invariant schema migration did not reach version " *
            string(OGRAPH_INVARIANT_STORE_SCHEMA_VERSION)))
    isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check")) ||
        throw(ArgumentError("o-graph invariant store foreign-key check failed"))
    true
end

"""
    migrate_ograph_invariant_store!(store::MoveGraphStoreDB) -> MoveGraphStoreDB

Prepare append-only storage for derived edge-group presentations and their abelianizations.

# Arguments
`store` must be an open writable database. Base migration and verification
must succeed; their errors propagate. An unsupported invariant schema version,
invalid foreign keys, or failed replay of existing invariant evidence throws
`ArgumentError` identifying the failure; file and database errors propagate.

# Returns
The same `MoveGraphStoreDB`, with the current
[invariant storage schema](../specs/ograph-relation-db.md#invariant-storage-and-replay) installed.

# Guarantees
1. Successful return means the invariant schema and its append-only guards are installed and existing invariant rows have passed the current store verifier. [proved/agent: by construction]
2. Repeating the migration on a valid current store preserves its stored invariant records. [proved/agent: by construction]
3. The database is modified; no new mathematical invariant of a presentation is inferred by migration. [proved/agent: by construction]

# Not guaranteed
The base migration can remain committed if a later invariant migration fails.
Replay does not resolve the raw edge presentation's identification with a
spine fundamental group, as recorded in the
[closed examples](../specs/ograph-colour-invariants.md#closed-examples).

# Sources
cf:docs/specs/ograph-colour-invariants.md#closed-examples — mathematical scope of the stored closed-normal presentation and H1 comparison.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_invariant_store!(store); migrate_ograph_invariant_store!(store);

julia> Topology._ograph_invariant_existing_schema_version(store.db)
1

julia> close_store!(store);
```

See also [`put_ograph_invariants!`](@ref), [`verify_ograph_invariant_store`](@ref).

Role: other
"""
function migrate_ograph_invariant_store!(store::MoveGraphStoreDB)
    migrate_ograph_store!(store)
    verify_ograph_store(store)
    version = _ograph_invariant_existing_schema_version(store.db)
    version === nothing || version <= OGRAPH_INVARIANT_STORE_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph invariant schema version $version"))
    sql = read(_ograph_invariant_schema_path(), String)
    MoveGraphStore._transaction(store.db) do
        MoveGraphStore._execute_schema!(store.db, sql)
        _assert_ograph_invariant_schema(store)
        # A database can already advertise schema v1 while missing the
        # append-only triggers (for example, an interrupted/manual bootstrap).
        # Replay every existing row before committing the installed guards so
        # invalid legacy evidence cannot become permanently frozen in place.
        verify_ograph_invariant_store(store)
    end
    store
end

function _pi1_relators_json(presentation::OGraphPi1Presentation)
    _ograph_canonical_json([collect(word) for word in presentation.relators])
end

function _assert_pi1_row(row, presentation::OGraphPi1Presentation,
                         created_at_utc::AbstractString)
    row === nothing && throw(ArgumentError("π1 presentation is not stored"))
    actual = (String(row.presentation_id), String(row.state_id),
        String(row.presentation_hash), String(row.algorithm_id),
        Int(row.algorithm_version), String(row.convention_id),
        Int(row.generator_count), Int(row.relator_count),
        String(row.relators_json), String(row.presentation_serialization),
        String(row.evidence_json), String(row.created_at_utc))
    expected = (presentation.presentation_id, presentation.state_id,
        presentation.presentation_hash, presentation.algorithm_id,
        presentation.algorithm_version, presentation.convention_id,
        presentation.generator_count, length(presentation.relators),
        _pi1_relators_json(presentation), presentation.serialization,
        presentation.evidence_json, String(created_at_utc))
    actual == expected || throw(ArgumentError(
        "π1 presentation id already exists with different data"))
    true
end

"""
    put_ograph_pi1_presentation!(store::MoveGraphStoreDB, graph::ClosedNormalOGraph; created_at_utc::AbstractString) -> OGraphPi1Presentation

Store the raw edge presentation computed from a persisted closed normal o-graph.

# Arguments
`store` must have the invariant schema installed by
[`migrate_ograph_invariant_store!`](@ref). `graph` must have its exact canonical
state already registered by [`put_closed_ograph!`](@ref). These checks throw
`ArgumentError` identifying schema or persisted-state disagreement; database
errors can occur first for absent tables. `created_at_utc` is caller-supplied
text, not a validated timestamp. Reusing an identity requires all stored data,
including this timestamp, to agree, otherwise `ArgumentError` reports a conflict.

# Returns
The computed `OGraphPi1Presentation`, with fields and mathematical scope as
in [`ograph_pi1_presentation`](@ref).

# Guarantees
1. Successful return leaves an exact matching presentation row in the store;
   an identical repeated call returns equal data without replacing the row.
   The graph is unchanged. [proved/agent: by construction]
2. The stored group has the raw-edge-presentation interpretation of the
   [closed-normal pipeline](../specs/ograph-colour-invariants.md#closed-examples).
   Persistence does not strengthen its conditional identification with the
   spine's fundamental group. [proved/agent: by construction]

# Not guaranteed
The state is not inserted automatically. Storage does not prove group
triviality, identify a manifold or resolve the raw-presentation identification
gap. Arbitrary timestamp strings are not authenticated provenance.

# Sources
cf: lib/topology.jl/docs/specs/ograph-colour-invariants.md §Closed examples — the algebraic and geometric scope of the stored presentation.
cf: lib/topology.jl/db/ograph_invariant_schema.sql — append-only presentation records.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_invariant_store!(store); g = certify_closed_ograph(vl_example("s3_abalone")); put_closed_ograph!(store, g);

julia> p = put_ograph_pi1_presentation!(store, g; created_at_utc="2026-11-25T00:00:00Z");

julia> put_ograph_pi1_presentation!(store, g; created_at_utc="2026-11-25T00:00:00Z") == p
true

julia> close_store!(store);
```

See also [`put_ograph_h1!`](@ref), [`put_ograph_invariants!`](@ref).

Role: other
"""
function put_ograph_pi1_presentation!(store::MoveGraphStoreDB,
        graph::ClosedNormalOGraph; created_at_utc::AbstractString)
    _assert_ograph_invariant_schema(store)
    state = ograph_canonical_state(graph)
    presentation = ograph_pi1_presentation(graph)
    MoveGraphStore._transaction(store.db) do
        persisted = get_state(store, state.state_id)
        persisted !== nothing && _same_canonical_state(persisted, state) ||
            throw(ArgumentError("π1 input is not the exact persisted o-graph state"))
        MoveGraphStore._one(store.db,
            "SELECT state_id FROM ograph_states WHERE state_id=?", (state.state_id,)) !== nothing ||
            throw(ArgumentError("π1 input is not a registered o-graph state"))
        row = MoveGraphStore._one(store.db, """
            SELECT * FROM ograph_pi1_presentations WHERE presentation_id=?
        """, (presentation.presentation_id,))
        if row === nothing
            MoveGraphStore._dbexec(store.db, """
                INSERT INTO ograph_pi1_presentations(
                  presentation_id,state_id,presentation_hash,algorithm_id,
                  algorithm_version,convention_id,generator_count,relator_count,
                  relators_json,presentation_serialization,evidence_json,created_at_utc)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
            """, (presentation.presentation_id, presentation.state_id,
                    presentation.presentation_hash, presentation.algorithm_id,
                    presentation.algorithm_version, presentation.convention_id,
                    presentation.generator_count, length(presentation.relators),
                    _pi1_relators_json(presentation), presentation.serialization,
                    presentation.evidence_json, String(created_at_utc)))
            row = MoveGraphStore._one(store.db, """
                SELECT * FROM ograph_pi1_presentations WHERE presentation_id=?
            """, (presentation.presentation_id,))
        end
        _assert_pi1_row(row, presentation, created_at_utc)
        presentation
    end
end

function _h1_factors_json(result::OGraphH1Result)
    _ograph_canonical_json(string.(result.torsion_invariant_factors))
end

function _assert_h1_row(row, result::OGraphH1Result,
                        created_at_utc::AbstractString)
    row === nothing && throw(ArgumentError("H1 result is not stored"))
    actual = (String(row.h1_result_id), String(row.presentation_id),
        String(row.algorithm_id), Int(row.algorithm_version),
        String(row.backend), String(row.backend_version), Int(row.free_rank),
        String(row.torsion_invariant_factors_json),
        String(row.torsion_subgroup_order_decimal),
        _ograph_nullable_string(row.h1_finite_order_decimal),
        String(row.h1_signature), String(row.result_fingerprint),
        String(row.evidence_json), String(row.created_at_utc))
    expected = (result.result_id, result.presentation_id, result.algorithm_id,
        result.algorithm_version, result.backend, result.backend_version,
        result.free_rank, _h1_factors_json(result),
        string(result.torsion_subgroup_order),
        result.h1_finite_order === nothing ? nothing : string(result.h1_finite_order),
        result.signature, result.result_fingerprint, result.evidence_json,
        String(created_at_utc))
    actual == expected || throw(ArgumentError(
        "H1 result id already exists with different data"))
    true
end

"""
    put_ograph_h1!(store::MoveGraphStoreDB, presentation::OGraphPi1Presentation; created_at_utc::AbstractString) -> OGraphH1Result

Store the integral abelianization of an exactly persisted group presentation.

# Arguments
`store` must have the invariant schema. `presentation` must exactly match an
existing presentation row. Missing or disagreeing rows throw `ArgumentError`
identifying the input conflict. `created_at_utc` is stored text; an existing
result must agree in every field, including this timestamp, or `ArgumentError`
reports differing data. Schema, database and [`ograph_h1`](@ref) errors propagate.

# Returns
The `OGraphH1Result` computed by `ograph_h1`, with its invariant factors,
rank, orders and algorithm and presentation identities.

# Guarantees
1. The stored result is the integral abelianization of the input presentation's
   group, with the scope of [`ograph_h1`](@ref). The presentation is unchanged.
   [proved/agent: by construction]
2. Successful return leaves an exactly matching result row. An identical
   repeated call returns equal data without replacing that row.
   [proved/agent: by construction]

# Not guaranteed
The presentation is not automatically inserted. Identification with a spine's
or a capped thickening's first homology has the same conditional scope as the
[invariant specification](../specs/ograph-colour-invariants.md#closed-examples).
Storage neither identifies a manifold nor proves its fundamental group trivial.

# Sources
cf: lib/topology.jl/db/ograph_invariant_schema.sql — append-only H1 records linked to exact presentations.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_invariant_store!(store); g = certify_closed_ograph(vl_example("s3_abalone")); put_closed_ograph!(store, g);

julia> p = put_ograph_pi1_presentation!(store, g; created_at_utc="2026-11-25T00:00:00Z"); h = put_ograph_h1!(store, p; created_at_utc="2026-11-25T00:00:00Z");

julia> (h.free_rank, h.torsion_invariant_factors, h.h1_finite_order)
(0, (), 1)

julia> close_store!(store);
```

See also [`put_ograph_pi1_presentation!`](@ref), [`verify_ograph_invariant_store`](@ref).

Role: other
"""
function put_ograph_h1!(store::MoveGraphStoreDB,
                        presentation::OGraphPi1Presentation;
                        created_at_utc::AbstractString)
    _assert_ograph_invariant_schema(store)
    MoveGraphStore._transaction(store.db) do
        presentation_row = MoveGraphStore._one(store.db, """
            SELECT * FROM ograph_pi1_presentations WHERE presentation_id=?
        """, (presentation.presentation_id,))
        presentation_row === nothing && throw(ArgumentError(
            "H1 input presentation is not stored"))
        _assert_pi1_row(presentation_row, presentation,
                        String(presentation_row.created_at_utc))
        result = ograph_h1(presentation)
        row = MoveGraphStore._one(store.db,
            "SELECT * FROM ograph_h1_results WHERE h1_result_id=?",
            (result.result_id,))
        if row === nothing
            MoveGraphStore._dbexec(store.db, """
                INSERT INTO ograph_h1_results(
                  h1_result_id,presentation_id,algorithm_id,algorithm_version,
                  backend,backend_version,free_rank,
                  torsion_invariant_factors_json,torsion_subgroup_order_decimal,
                  h1_finite_order_decimal,h1_signature,result_fingerprint,
                  evidence_json,created_at_utc)
                VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """, (result.result_id, result.presentation_id, result.algorithm_id,
                    result.algorithm_version, result.backend, result.backend_version,
                    result.free_rank, _h1_factors_json(result),
                    string(result.torsion_subgroup_order),
                    result.h1_finite_order === nothing ? nothing :
                        string(result.h1_finite_order),
                    result.signature, result.result_fingerprint,
                    result.evidence_json, String(created_at_utc)))
            row = MoveGraphStore._one(store.db,
                "SELECT * FROM ograph_h1_results WHERE h1_result_id=?",
                (result.result_id,))
        end
        _assert_h1_row(row, result, created_at_utc)
        result
    end
end

"""
    put_ograph_invariants!(store::MoveGraphStoreDB, graph::ClosedNormalOGraph; created_at_utc::AbstractString) -> NamedTuple

Store a closed normal o-graph's raw group presentation and its integral abelianization.

# Arguments
`store`, `graph` and `created_at_utc` must meet the checked preconditions of
[`put_ograph_pi1_presentation!`](@ref); the resulting presentation must meet
those of [`put_ograph_h1!`](@ref). Their schema, identity-conflict, database and
computation errors propagate. In particular, register `graph` first using
`put_closed_ograph!` after migrating the invariant store.

# Returns
A `NamedTuple` with `presentation::OGraphPi1Presentation` and
`h1::OGraphH1Result` computed from that presentation.

# Guarantees
1. Successful return leaves exact matching rows for both records, and
   `h1.presentation_id == presentation.presentation_id`. The graph is unchanged.
   [proved/agent: by construction]
2. The group and abelian group have exactly the meanings and geometric
   qualifications of the [closed-normal invariant pipeline](../specs/ograph-colour-invariants.md#closed-examples).
   [proved/agent: by construction]

# Not guaranteed
The two insertions are not one atomic operation: if H1 storage fails, the
presentation may remain stored. This is not manifold recognition or a proof
that the raw edge group is the spine's fundamental group.

# Sources
cf: lib/topology.jl/docs/specs/ograph-colour-invariants.md §Closed examples — the two-stage invariant pipeline and its geometric limitation.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_invariant_store!(store); g = certify_closed_ograph(vl_example("s3_abalone")); put_closed_ograph!(store, g);

julia> pair = put_ograph_invariants!(store, g; created_at_utc="2026-11-25T00:00:00Z");

julia> put_ograph_h1!(store, pair.presentation; created_at_utc="2026-11-25T00:00:00Z") == pair.h1
true

julia> close_store!(store);
```

See also [`verify_ograph_invariant_store`](@ref), [`ograph_h1`](@ref).

Role: other
"""
function put_ograph_invariants!(store::MoveGraphStoreDB,
                                graph::ClosedNormalOGraph;
                                created_at_utc::AbstractString)
    presentation = put_ograph_pi1_presentation!(store, graph;
                                               created_at_utc=created_at_utc)
    result = put_ograph_h1!(store, presentation; created_at_utc=created_at_utc)
    (presentation=presentation, h1=result)
end

"""
    verify_ograph_invariant_store(store::MoveGraphStoreDB) -> NamedTuple

Verify stored o-graph presentations and first-homology results by recomputation.

# Arguments
`store` must be open with the current o-graph and invariant schemas installed.
The requirements of [`verify_ograph_store`](@ref) also apply. Missing or
unsupported schemas, inconsistent replay data, unsupported H1 algorithms or
backends, and stored external group-enrichment rows throw `ArgumentError`.
Database, decoding and underlying invariant errors propagate.

# Returns
A summary with `ograph_schema_version`, `invariant_schema_version`, `states`,
`presentations`, `h1_results`, `backend_runs` and `identifications`.
Successful verification currently requires the last two counts to be zero,
as well as the stored finite-quotient count.

# Guarantees
1. On successful return, the joined stored presentations agree with
   `ograph_pi1_presentation` on their certified state payloads, and every H1
   row agrees with the abelianization of its referenced presentation.
   [proved/agent: by construction]
2. The mathematical meaning is exactly that of [`ograph_pi1_presentation`](@ref)
   and [`ograph_h1`](@ref), including their limits on identifying a space.
   This verification adds consistency of the stored results, not a stronger
   topological identification. [proved/agent: by construction]
3. Verification leaves the store unchanged. [proved/agent: by construction]

# Not guaranteed
This is not a manifold recognizer, a proof of group triviality or finiteness,
or a replay engine for external group backends. Stored backend version strings
are retained metadata, not proof that that historical backend was executed.
No complete invariant-store format or compatibility policy is owned by a
specification page yet.

# Sources
cf:lib/topology.jl/src/ograph_invariants.jl ograph_pi1_presentation and ograph_h1 — presentation and abelianization contracts; invariant schema v1 replay scope.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_invariant_store!(store);

julia> for name in ("s3_abalone", "lens_2_1", "s2xs1")
           g = certify_closed_ograph(vl_example(name)); put_closed_ograph!(store, g); put_ograph_invariants!(store, g; created_at_utc="2026-09-04T00:00:00Z");
       end

julia> s = verify_ograph_invariant_store(store); (s.states, s.presentations, s.h1_results)
(3, 3, 3)

julia> close_store!(store);
```

See also [`put_ograph_invariants!`](@ref), [`verify_ograph_store`](@ref).

Role: condition
"""
function verify_ograph_invariant_store(store::MoveGraphStoreDB)
    core = verify_ograph_store(store)
    _assert_ograph_invariant_schema(store)
    presentation_rows = MoveGraphStore._rows(store.db, """
        SELECT presentation.*, state.payload
        FROM ograph_pi1_presentations presentation
        JOIN states state ON state.state_id=presentation.state_id
        ORDER BY presentation.presentation_id
    """)
    presentations = Dict{String,OGraphPi1Presentation}()
    for row in presentation_rows
        graph = certify_closed_ograph(_decode_ograph_state_payload(row.payload))
        expected = ograph_pi1_presentation(graph;
            algorithm_id=String(row.algorithm_id),
            algorithm_version=Int(row.algorithm_version),
            convention_id=String(row.convention_id))
        _assert_pi1_row(row, expected, String(row.created_at_utc))
        presentations[expected.presentation_id] = expected
    end
    h1_rows = MoveGraphStore._rows(store.db,
        "SELECT * FROM ograph_h1_results ORDER BY h1_result_id")
    for row in h1_rows
        presentation_id = String(row.presentation_id)
        haskey(presentations, presentation_id) || throw(ArgumentError(
            "H1 result references an unreplayable π1 presentation"))
        String(row.algorithm_id) == OGRAPH_H1_ALGORITHM &&
            Int(row.algorithm_version) == OGRAPH_H1_ALGORITHM_VERSION &&
            String(row.backend) == "AbstractAlgebra" || throw(ArgumentError(
                "unsupported stored H1 replay algorithm or backend"))
        expected = _ograph_h1_from_matrix(
            presentation_id,
            _ograph_exponent_sum_matrix(presentations[presentation_id]);
            algorithm_id=String(row.algorithm_id),
            algorithm_version=Int(row.algorithm_version),
            backend=String(row.backend),
            backend_version=String(row.backend_version))
        _assert_h1_row(row, expected, String(row.created_at_utc))
    end
    backend_runs = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM ograph_group_backend_runs").n)
    quotient_rows = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM ograph_finite_quotient_fingerprints").n)
    identifications = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM ograph_manifold_identifications").n)
    backend_runs == 0 && quotient_rows == 0 && identifications == 0 ||
        throw(ArgumentError(
            "external group-enrichment replay is not implemented in invariant schema v1"))
    (
        ograph_schema_version=core.ograph_schema_version,
        invariant_schema_version=OGRAPH_INVARIANT_STORE_SCHEMA_VERSION,
        states=core.states,
        presentations=length(presentation_rows),
        h1_results=length(h1_rows),
        backend_runs=backend_runs,
        identifications=identifications,
    )
end
