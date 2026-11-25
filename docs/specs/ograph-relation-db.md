# Closed o-graph relation and invariant storage

This page specifies the current Area-3 adapters in
[`ograph_move_db.jl`](../../src/ograph_move_db.jl) and
[`ograph_invariants.jl`](../../src/ograph_invariants.jl). The generic move graph
uses schema 2, the o-graph relation extension schema 2, explicitly enabled
typed transitions schema 1, and invariant storage schema 1. The authoritative
SQL field sets and constraints are in
[`ograph_relation_schema.sql`](../../db/ograph_relation_schema.sql) and
[`ograph_invariant_schema.sql`](../../db/ograph_invariant_schema.sql).
The working unweighted search projection has its own
[specification](ograph-search-db.md); its storage verifier has a narrower scope.

## Canonical JSON

Status: fixed 2026-11-25 Codex (rule a). `ograph_canonical_json` serializes
acyclic supported values with no optional whitespace. `nothing` becomes
`null`, actual Booleans become `true`/`false`, and integers must fit signed
64-bit range and are printed in decimal. Finite `AbstractFloat` values and
strings use `JSON.json`; nonfinite floats, other numeric types and `missing`
are rejected. Tuples and vectors become ordered arrays. A `NamedTuple` becomes
an object using its field names. Dictionary keys must be strings or symbols,
are converted to strings, must remain distinct, and are sorted in Julia string
order. Their values are serialized recursively. Other value types are rejected.

String escaping and floating-point spelling are delegated to the installed
JSON.jl version. This is the repository's canonical encoding, not a claim of
conformance to an external canonical-JSON standard. Integer `1` and float `1.0`
need not have the same bytes. The [database tests](../../test/ograph_move_db.jl)
exercise key-order independence, collisions after key conversion, nonfinite
numbers and the signed-64-bit boundary. Cross-version JSON.jl byte stability
remains unverified: version changes need comparison of these serialized
identities, including strings and floating-point values, before reuse.

## State and provenance identity

A closed state uses the [canonical presentation serialization](ograph-terminology.md#presentation-canonicalization-and-serialization),
Area `3.0`, representation `"o-graph"`, canonical version `1`, and true-vertex
complexity. The identity domain is
`area3/closed-normal-ograph/presentation-v1`. Its replay payload contains
`identity_domain`, `schema`, `canonicalization`, `gauss`, `vertex_signs`,
`convention_id`, and the historical `closedness` evidence object specified in
[terminology](ograph-terminology.md#persisted-closedness-evidence).
The generic descriptor's `state_hash` is lowercase SHA-256 of the serialization;
`state_id` is lowercase SHA-256 of `v1`, a NUL byte, and that serialization.
The payload is replay evidence, not an additional manifold identifier.

`OGraphProvenance(source_kind, locator; snapshot_sha256, metadata)` strips
outer whitespace from the two required nonempty strings. Source kinds are
caller labels, not an enumerated authority hierarchy. A supplied snapshot is
stripped, lowercased and required to be 64 hexadecimal characters. Metadata
uses canonical JSON. The identity object has keys `source_kind`, `locator`,
`snapshot_sha256` (possibly null), and parsed `metadata`; the source ID is
`ograph-source-` followed by its canonical-JSON SHA-256. Construction does not
read the locator, compute its snapshot or validate the mathematical claim.

State insertion accepts an empty provenance vector (the default), re-derives
each supplied provenance identity, and stores its source association.
Reusing an identity with conflicting data is
an error. A matching source record certifies reproducible metadata, not the
truth or current availability of the source. Raw field constructors of record
types do not bypass the insertion checks.

## Relation profiles

Profile names are stripped and lowercased; the positive integer version and
move set participate in identity. Supported unweighted move sets are:

| Profile | Moves | Parent | Metadata requirement |
|---|---|---|---|
| `presentation` | none | none | `moves_complete=true`, `evidence_status="identity-only"` |
| `combed` | `MP`, `PS` | registered presentation profile | `verified` with complete=true, or `unimplemented` with complete=false |
| `manifold`, version 1 | `CP`, `MP`, `PS` | registered combed profile | same paired status/complete rule |
| `manifold`, version >=2 | `CP`, `MP`, `PS`, `S` | registered combed profile | same paired status/complete rule |

Weighted profiles reserved in SQL are not supported by this unweighted
constructor. Version >=2 inserts a `move_capabilities` snapshot into the
conventions object only if that key is absent; a supplied snapshot is retained.
The snapshot describes type/direction support, separately from calculus
completeness. The move source gates are owned by
[move verification](ograph-move-verifier.md).

The fingerprint hashes canonical JSON of `profile_name`, `profile_version`,
`parent_profile_id`, ordered `moves`, `moves_complete`, `evidence_status`, and
parsed `conventions`. The profile ID is `ograph-profile-` plus that fingerprint.
`created_at_utc` is excluded from identity but must agree when an existing row
is reused. Registration checks the parent kind and the Area-3 move registry.
The status string `verified` is stored profile metadata; it is not an API
guarantee's evidence status and does not itself prove the calculus complete.

## Finite campaign metadata

`register_ograph_campaign!` requires distinct persisted closed-state IDs and
a registered profile. IDs are sorted before hashing and stored with zero-based
ordinals. `N_input` is nonnegative, `B_path >= N_input`, and input complexity
must be finite integral and at most `N_input`. These are declared campaign
bounds; registration does not perform a search or enforce all intermediate
path bounds by itself.

`caps` and `counts` use canonical JSON. The input fingerprint hashes the sorted
ID array. The campaign fingerprint hashes `campaign_id`, `profile_id`,
`N_input`, `B_path`, `input_state_ids`, parsed `caps`, `input_fingerprint`,
`engine_version`, `verifier_version`, `status`, `stop_reason`, and parsed
`counts`. The caller supplies a nonempty stripped campaign ID. SQL accepts
`planned`, `completed`, `budget-exhausted`, or `failed`. These statuses record
the declared finite run; `completed` is not a nonexistence theorem.
Reuse requires matching metadata, including the creation timestamp, despite
that timestamp's exclusion from the fingerprint. Timestamps are stored text;
the API does not validate UTC syntax or chronological truth.

## Typed transition evidence

Ordinary relation migration does not enable typed transition insertion; call
`migrate_ograph_transition_store!` explicitly. Evidence version is
`ograph-transition-evidence/1`. Its required fields are `evidence_version`,
`source_state_id`, `target_state_id`, `source_state_hash`, `target_state_hash`,
`move_type`, `direction`, `profile_id`, `parameters`, `location`, `witness`,
`replay`, and `operation_id`. Parameters must contain the same `operation_id`,
`parameter`, `vertices`, `source_locations`, and `target_locations`. PS also
requires `method_b_witness`, equal to the embedded verified witness's Method-B
record. Additional envelope fields are permitted and affect its fingerprint.
The witness's exact nested schema is specified
[separately](ograph-move-verifier.md#witness-coordinates-and-serialization).

Insertion requires a successful semantic verification report, freshly replays
its witness, and binds canonical endpoints, hashes, move type, direction,
parameter, vertices and occurrence locations to the database row. The profile
must allow the move, and endpoints must already be registered closed states.
The canonical evidence JSON is fingerprinted and stored together with separate
witness/replay JSON fields. An envelope success flag alone is insufficient.
The caller's replay payload is recorded metadata; acceptance comes from the
fresh typed replay, not from trusting that payload.

`verify_ograph_store` checks base/extension versions, SQLite integrity and
foreign keys, exact state/closedness replay, provenance identities, profiles,
campaign fingerprints and registered input records. When enabled, typed
transition evidence is decoded and semantically replayed, and stored relation
certificate steps are checked in order against their endpoints and profile.
This verifies stored evidence within the move gates' scope. It does not
establish completeness of a search or classify all manifolds represented by
unconnected states.

## Invariant storage and replay

The presentation table stores the state association, algorithm/version and
convention, generator/relator counts, relators, serialization, fingerprint,
evidence JSON and creation timestamp. The implemented presentation algorithm
is `koda-edge-presentation` version 1 with convention
`pd-code-sage-matlab-closed-ordinary/1`. H1 rows store the presentation
association, algorithm/backend versions, free rank, torsion invariant factors,
torsion subgroup order, optional finite total order, signature, fingerprint,
evidence and timestamp. Arbitrary-size order values are decimal strings;
finite total order is null exactly when free rank is positive.

H1 uses `exponent-sum-smith-normal-form` version 1 with `AbstractAlgebra`.
Its mathematical object is the abelianization of the stored edge-group
presentation. Identification with the spine's fundamental group, and hence
with the spine's first homology by this route, remains
[open: missing identification of the raw edge presentation with the spine's
fundamental group]. The matching finite comparisons in
[closed examples](ograph-colour-invariants.md#closed-examples) support that
identification on those inputs; a source proof with matching conventions or
a counterexample would settle or overturn it. Storage does not upgrade it.

Invariant rows are append-only under the installed update/delete triggers.
Reinserting the same identity checks all expected fields, including timestamps;
conflicts are errors. `put_ograph_invariants!` writes the presentation and then
H1 in separate operations, so a failed H1 write can leave the presentation.
Migration installs guards and replays existing rows before committing the
invariant extension. Its preceding base migration can already have committed.

`verify_ograph_invariant_store` first verifies the relation store, then
recomputes each presentation and its Smith-normal-form result. It requires
supported algorithm IDs/versions and backend name; the stored backend-version
string is reused as metadata rather than requiring the installed backend to
have that same version. Reserved external group-backend, finite-quotient and
manifold-identification tables have no supported enrichment replay here:
nonempty rows cause rejection. These tables are not an implemented manifold
recognition interface.

## Evidence and compatibility policy

Status: fixed 2026-11-25 Codex (rule a) for the byte formats, identity inputs,
ordinals and metadata conventions above. The implementation correspondence is
[proved/agent: by construction]. The
[database tests](../../test/ograph_move_db.jl) exercise provenance forgery,
profile hierarchy, campaign identity and typed transition acceptance/rejection;
[invariant tests](../../test/ograph_invariants.jl) exercise deterministic
presentations/H1, duplicate insertion, conflict rejection and append-only
guards. The distinction from search-store verification is embodied by the
[search tests](../../test/ograph_search_db.jl). These are repository regression
cases, not a claim that they were rerun while this entry was written.

Readers and migrations accept only the versions explicitly supported by the
implementation; old witness replay follows the
[versioned replay rule](ograph-move-verifier.md#persistence-and-replay).
No general forward compatibility, migration of unknown formats, immutable
physical database protection against schema tampering, or cross-library-version
JSON byte stability is established. Expected compatibility is limited to these
explicit paths; additional migration/replay tests are required to expand it.
