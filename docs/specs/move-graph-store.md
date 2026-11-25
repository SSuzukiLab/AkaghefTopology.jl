# Move graph storage

This page owns the common store model. The o-graph extension is specified in
[ograph-search-db.md](ograph-search-db.md).

## Definitions

### Domain vocabulary

An area is a stored domain tag: `2.0` for knot/link diagram moves, `2.5` for
mixed-link/lo-graph moves, `3.0` for 3-manifold moves. A representation is a
lowercase, whitespace-trimmed string naming the chosen encoding. The move
vocabulary is the whitelist in [domain_taxonomy.jl](../../src/domain_taxonomy.jl):
R1/R2/R3 and their Reidemeister aliases in area 2; mixed/lo-graph aliases in
area 2.5; MP/PS/CP/S/02, their listed inverses, BMP/B02/H and listed
3-manifold aliases in area 3. Representation changes and serialization are
not topological moves merely because they change stored data.

### Textual state identity

A `CanonicalState` contains caller-canonicalized serialization, a positive
canonical version, domain/representation metadata, optional complexity and
payload. Its `state_hash` is lowercase hexadecimal SHA-256 of the serialization
code units. Its `state_id` hashes the code units of `vVERSION`, followed by a
NUL byte and the serialization. The metadata and payload are not additional
digest inputs. Canonicality is a property of the caller's encoding; computing
these digests neither canonicalizes nor parses it.

### Graph and evidence records

A `MoveType` names a domain, inverse type and inverse direction. A transition
is a directed stored edge with source and target state IDs, move type,
direction, location, parameters, engine/verifier versions and verification
evidence. It stores a positive verifier verdict. A certificate-path record is
a path string attached to a transition; it is not the contents or validation
of the external certificate. An invariant record is an attributed value
attached to a state. A certificate path returned by graph search is a path
in the stored directed graph, not an independent proof of its edges.

### Schema and classifications

`MoveGraphStoreDB` wraps a SQLite connection (`Store` is its alias). The common
schema is version 2 in [knot_move_schema.sql](../../db/knot_move_schema.sql).
It has states, move types, transitions, certificate paths, invariant values,
KnotInfo snapshots/raw records, knot types, classification runs and per-state
classifications. `ClassificationRun` identifies a dataset/run and its
metadata. Per-state classification records attach labels and evidence to that
run; they do not execute the classification algorithm. A KnotInfo import
stores the provided snapshot identifier and first-sheet records. The supplied
identifier is provenance metadata, not a digest recomputed by this importer.

### Knot adapter

The signed Gauss adapter canonicalizes its diagram before recording area 2,
representation `knotgauss`, the canonical serialization, a Gauss/orientation
payload and crossing-count complexity. A stored Reidemeister witness uses
its source and target canonical identities. Deletion has stored direction
`:inverse`; insertion and R3 use `:forward`. This is a direction label, not
an exchange of stored source and target. Adapter evidence records the local
candidate and both serializations. The adapter's finite local verification
has the [Gauss source boundaries](knot-gauss.md#coordinate-choices-and-evidence).

### Search-run provenance

A search run records a presentation-specific experiment in the same SQLite
store: its label, caller-supplied scope, script path and content digest,
creation/completion times, status and reported counts. This records provenance;
it neither runs the script nor makes the scope an enforced search bound.
Heegaard and o-graph records use their respective migrated search tables.
The type-directed common interface selects the presentation after the store
argument, as specified in [MoveAPI](module-structure.md#31-moveapi).

An HC search operation records a candidate move-row identity, source and
optional target state identities, depth, outcome, inverse-check status and
optional error data. Outcome distinguishes success, rejected local moves and
errors. Inverse status distinguishes unchecked, successful, mismatching and
failed inverse attempts. Success and inverse success are separate facts.
Identity mode is `canonical` or `literal`, using the HC search representation
owned by the Heegaard interface. A diagnostic snapshot has homology and
linking result dictionaries, including their `ok` flags; comparison with a
starting snapshot is a finite diagnostic, not a general invariant proof.
These record meanings follow
[hc_search_db.jl](../../src/heegaard/hc_search_db.jl) and
[the HC search tests](../../test/heegaard/hc_search.jl).

### Stored search paths

A search path is a directed sequence of successful stored operation records.
Its length counts those edges, independently of each operation's recorded
exploration depth. The HC path model has no repeated state. The o-graph model
allows a return to the initial state only as the last step, as owned by the
[o-graph path restriction](ograph-search-db.md#exact-h-insertion-feasibility).
A geometric-hash target selects presentations with that stored geometry;
it does not identify their represented manifolds. A zero-edge HC path is a
formal identity at the supplied state ID, even when no state row exists.
The o-graph path target must be a stored search state.

These are repository graph-query definitions in
[HC search](../../src/heegaard/hc_search_db.jl) and
[o-graph search](../../src/ograph_search_db.jl), not additional move theorems.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| M1 | Hash and ID match the serialization/version; version is positive and area belongs to the vocabulary. | Keyword state construction checks all; the raw eight-field constructor checks digest consistency only. |
| M2 | Serialized text, payload and declared representation describe the same mathematical state and use the stated canonicalization. | Caller obligation; store cannot establish it. |
| M3 | A transition has a successful verifier verdict, compatible domains and registered move metadata, with required structured evidence. | API and schema validation; the verifier callback remains trusted. |
| M4 | A classification claim identifies its run/dataset and evidence, not only its label. | Metadata/schema checks; mathematical adequacy is external. |
| M5 | A schema migration starts from a supported schema and valid existing evidence fields. | Migration checks and transaction; no mathematical replay is implied. |
| M6 | A neighbor expansion has nonnegative `store_cap` and `max_crossings`, and an open writable supported store. | Negative bounds raise `ArgumentError`; database errors may be counted as rejected attempts by the expansion. |
| M7 | A search store is open and has the supported migrated tables for the selected presentation. | Presentation-specific schema checks; SQLite failures propagate. |
| M8 | Run registration has an existing readable script path, a text label and scope representable by the search JSON format. | Missing paths/schema or unsupported values raise the delegated validation errors; read/database errors propagate. |
| M9 | HC expansion uses an existing run in the same store, a coherent source, one identity mode, meaningful source depth and a remaining-state allowance. Optional diagnostic containers/callbacks accept the generated problem records and a starting snapshot uses the current snapshot schema. | Identity/schema checks are explicit; record coherence and diagnostic compatibility are caller obligations. Per-candidate failures may be recorded rather than thrown. |
| M10 | A terminal status is `completed`, `budget-exhausted` or `failed` for HC; `completed-found`, `completed-none`, `budget-exhausted` or `failed` for o-graphs. Reported counts fit the presentation's search JSON format. | Status and serialization are checked; truth/completeness of reported counts is external. |
| M11 | A stored search writing has the selected presentation's supported schema and valid reconstruction data. | Loading checks reconstruction, not a cryptographic reauthentication of the state identity or a geometric manifold certificate. |
| M12 | Path bounds satisfy `0 <= min_depth <= max_depth` and `limit > 0`; identifiers are store-produced strings, HC row filters name permitted catalog rows, and o-graph kinds/patterns use MP, 02 or PS. Exactly one o-graph target selector is supplied. | Bounds/target/kind checks are explicit; integer and string conversions and database failures propagate. A path does not replay the underlying move evidence. |
| M13 | Classification-run method, fingerprints and timestamp are nonempty stripped strings; `canonical_version > 0`, `0 <= max_n <= barrier`; counts and conventions are JSON-serializable. | Registration checks fields and immutable metadata. A repeated key compares encoded counts/conventions exactly and retains its original timestamp. The raw `ClassificationRun` constructor only converts declared field types. |
| M14 | A classification state is stored in area 2 and its method is registered. `planar` is Boolean or integer 0/1. Status is `classified`, `unresolved` or `nonclassical`; a classified planar state has a nonempty knot name, an unresolved planar state has none, and a nonplanar state has no knot name and status `nonclassical`. | Julia and SQL check consistency; reusing the state/method key requires identical values. Mathematical planarity and knot identity are caller assertions. |
| M15 | A common-store path bound is a nonnegative integer and endpoints convert to state-ID strings. | The shortest-path query does not validate the bound or endpoint existence. Equal endpoint strings give the formal one-state identity path. |

## Choices and evidence boundaries

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)** for domain tags,
normalization and digest bytes. These are repository storage conventions,
established by [move_graph_store.jl](../../src/move_graph_store.jl),
[store tests](../../test/move_graph_store.jl) and
[knot move integration](../../test/knot_move_db_integration.jl).
Classification metadata is tested in
[N5 database tests](../../test/knot_n5_classification_db.jl) and
[N6 database tests](../../test/knot_n6_classical_classification_db.jl).
These local records have no literature definition requiring a cite-key alias.

question: unverified — A stored verdict certifies the mathematical move only
if its external verifier is sound for the stated representation and hypotheses.
The schema checks metadata, so it cannot establish that implication. Expected
validity rests on the named verifier's evidence; replay against its source and
target is the check that can support or refute it.

## Operations

`CanonicalState`, `state_from_serialization`, `state_hash`, `state_id`,
`MoveType`, `MoveGraphStoreDB`, `open_store`, `close_store!`, `migrate!`,
`put_state!`, `get_state`, `register_move_type!`, `get_move_type`,
`insert_verified_transition!`, `record_transition!`, `add_certificate!`,
`add_invariant!`, `neighbors`, `bfs_certificate_path`, `certificate_path`,
`import_knotinfo!`, `knotinfo_records`, `ClassificationRun`,
`register_classification_run!`, `put_state_classification!`,
`record_state_classification!`, `classification_run`, `state_classification`,
`classification_rows`, `classification_counts`, `normalize_area`,
`normalize_representation`, `move_area`, `validate_topological_move`.
