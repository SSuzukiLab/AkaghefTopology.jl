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
using .MoveGraphStore

const OGRAPH_PI1_PRESENTATION_SCHEMA = "akaghef.ograph-pi1-presentation/1"
const OGRAPH_PI1_ALGORITHM = "koda-edge-presentation"
const OGRAPH_PI1_ALGORITHM_VERSION = 1
const OGRAPH_PI1_CONVENTION = "pd-code-sage-matlab-closed-ordinary/1"
const OGRAPH_H1_SCHEMA = "akaghef.ograph-h1/1"
const OGRAPH_H1_ALGORITHM = "exponent-sum-smith-normal-form"
const OGRAPH_H1_ALGORITHM_VERSION = 1
const OGRAPH_INVARIANT_STORE_SCHEMA_VERSION = 1

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

"""Construct the deterministic Koda/MATLAB raw edge presentation.

For every PD row `[g,k,l,h]`, relations are emitted in MATLAB order: first all
`g*l^-1`, then all `h*g*k^-1`.  The mutable/cut `typeStrand == 'p'` branch is
outside this closed Area-3 adapter.
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

"""Compute exact integral H₁ from a raw presentation's exponent-sum matrix."""
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

"""Install the separately versioned, append-only derived-invariant layer."""
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

function put_ograph_invariants!(store::MoveGraphStoreDB,
                                graph::ClosedNormalOGraph;
                                created_at_utc::AbstractString)
    presentation = put_ograph_pi1_presentation!(store, graph;
                                               created_at_utc=created_at_utc)
    result = put_ograph_h1!(store, presentation; created_at_utc=created_at_utc)
    (presentation=presentation, h1=result)
end

"""Replay every implemented derived invariant from immutable o-graph states."""
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
