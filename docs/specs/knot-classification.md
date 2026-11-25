# Bounded knot move classifications

## Definitions

### Bounded graph

The vertices are canonical [signed Gauss representations](knot-gauss.md)
retained by the chosen classicality predicate, with crossing number at most
the run's bound. An edge is a replay-verified local Reidemeister transition
between retained vertices. Its undirected interpretation includes the inverse
transition. R1/R2 deletion witnesses account for insertions whose two endpoints
lie in the graph. Connected components are components of this finite graph,
not asserted components of the unrestricted knot isotopy relation.

### Classification record

`KnotN5Classification` and `KnotN6Classification` record a finite run's
states, edges, adjacency, components, representative polynomial values,
known representative matches, counts and verification diagnostics.
Known representative labels come from the local KnotInfo PD fixtures and
their mirrors. A component is resolved when its verified connectivity reaches
the applicable known representative; unmatched components remain unresolved.
A polynomial match alone is not the definition of resolution.

The `max_n` parameter bounds the enumerated input states. The separate
`barrier` field does not silently authorize a larger census: the current N5
and N6 classifiers restrict their actual supported runs as their contracts
state. State and edge arrays are transient workspace data. Persistence and
classification metadata are [MoveGraphStore](move-graph-store.md) objects.

### Targeted closure

`KnotN6TargetedBridgeResult` describes a barrier-7 search seeded by unresolved
barrier-6 components. `KnotN6TargetedBarrier8Result` describes the subsequent
targeted barrier-8 calculation. Their closure graph has seed-component nodes,
new canonical-state nodes and known-component boundary nodes. A reached known
boundary need not be expanded. Thus a closure is relative to its seed set,
allowed moves and stopping policy; it is not the complete crossing-7 or
crossing-8 census. A stored reach distance is distance in that closure graph,
where component nodes can contract several Reidemeister edges.

A resource-limited stop is a partial exploration. Queue exhaustion certifies
only the selected closure rule, provided every candidate was handled and its
verification diagnostics succeeded. Neither implies unrestricted knot
classification completeness.

The barrier-7 workspace is the transient graph data supplied to its optional
callback for the barrier-8 entry point. It contains seed and known component
nodes, canonical crossing-7 states and their node maps, closure adjacency,
component partition data and known-label connections. The public result does
not retain that workspace. Reusing it requires the same input classification
and bridge run; no independent fingerprint comparison establishes this link.

### Result field meanings

All exported result types in this page have unchecked field constructors.
A stored integer, flag or fingerprint is evidence only when linked to the run
that produced it. Mutable vectors and dictionaries remain mutable through
the record; these are not sealed certificates.

`KnotN5Edge` stores source and target state indices, the move-type symbol, and
the candidate directed from source to target. `KnotN5ComponentInvariant`
stores a component number, its representative state index, and that state's
Alexander and Jones records. `KnotN5KnownRepresentative` stores the snapshot
name, display label, mirror flag, diagram, state index, component number and
the two polynomial records. Their numbering refers to their enclosing result,
not to a persistent database. A raw constructor does not verify the witness,
name or polynomial values.

For `KnotN5Classification` and `KnotN6Classification`, `states` and
`state_index` give the finite population and inverse index. `edges` contains
the retained witnesses; each `adjacency` entry lists incident edge indices.
`components` gives each state's component number and `component_roots` its
chosen state representative. `state_counts` counts canonical states by
crossing number beginning at zero; `expected_state_counts` is the N5
reference vector. `raw_rooted_total` counts all rooted inputs and
`raw_rooted_classical_total` those retained by the classicality predicate.
The N6 `raw_rooted_classical_counts_by_n` gives the latter by crossing count.
The type-count dictionaries count canonical or rooted multiplicities by
resolved label; the corresponding unresolved count is separate.
`component_type_labels` is the N6 resolved component-to-label map.
`component_invariants`, `known_representatives`, and `unresolved_components`
hold the associated records and unmatched component numbers.

`direct_transition_count` counts considered candidate transitions before
state-pair deduplication; `duplicate_edge_transition_count` counts repeated
state pairs. `nonclassical_target_count` counts rejected targets. N5
`path_replay_count`/`path_replay_failures` count the path replay checks and
failures; N6 `spanning_tree_edge_replay_count`/`spanning_tree_replay_failures`
count its tree-edge checks and failures. N5 `edge_invariant_checks` and N6
`edge_invariant_sample_count` count polynomial comparisons;
`edge_invariant_failures` contains diagnostic strings. `fingerprint` is the
run's deterministic graph digest and `elapsed_seconds` is timing data.

For targeted result records, `seed_components` and `seed_state_count`
identify the selected input population. `input_fingerprint` (barrier 7),
`n6_input_fingerprint` and `barrier7_fingerprint` (barrier 8) identify its
inputs. `n7_state_count`/`n8_state_count` count the retained new-level states;
`imported_n7_state_count` counts the imported barrier-7 population.
`candidate_count`, `adopted_transition_count`, `duplicate_transition_count`,
`reused_n6_candidate_count`, `excluded_over_barrier_count`, and
`imported_graph_edge_count` describe considered, adopted, repeated, reused,
excluded and imported transitions respectively. `n7_processed_count`,
`n8_processed_count` and `max_frontier` measure exploration.
`transition_verification_failures`, `n6_membership_failures`,
`classicality_checks` and `classicality_failures` are validation diagnostics.
`connected_known_labels` and `known_reach_distance` are indexed by seed
component, with `nothing` for no reached known boundary. `resolved_components`
maps seeds to resolved labels; `resolved_component_type_counts` and
`resolved_state_type_counts` count seeds and their input states by label.
`closure_groups` partitions seed components by closure connectivity;
`unresolved_components` lists those not resolved. Barrier 8 also records
`component_status`, global `status`, and `resource_reason`. `fingerprint`
and `elapsed_seconds` record the closure digest and duration. Fields absent
from a particular type are not implied by this shared glossary.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| Q1 | Every state satisfies K1 and the run's classicality predicate, and state identities use the same canonical version. | Enumeration, canonicalization and run validation. |
| Q2 | Every retained transition replays from its stated source to target; its inverse interpretation uses the corresponding inverse witness. | Transition/replay checks; inspect failure counters. |
| Q3 | A classification claim names its crossing bound, classicality criterion, canonical equivalences and move set. | Required interpretation of a result; not inferred from a filename. |
| Q4 | A targeted result/workspace is attached to the exact input classification and seed selection. | Caller responsibility; structural import checks are partial and do not compare the input fingerprints. Changing the source requires a new run. |
| Q5 | Exhaustiveness requires successful completion with no omitted in-domain candidates or failed verification. | Resource termination and diagnostics must be inspected separately. |
| Q6 | All fields of a stored result or subordinate record agree with its enclosing run, state numbering, witnesses, polynomial data and counts. | Caller responsibility for raw constructors; not checked by field conversion. |
| Q7 | N5 bounds satisfy `0 <= max_n <= barrier <= 5`; the enumerated population ends at `max_n`. | The builder checks the range. |
| Q8 | N6 bounds satisfy `0 <= max_n <= barrier <= 6`; the enumerated population ends at `max_n`. | The builder checks the range. |
| Q9 | The KnotInfo JSON snapshot has the builder's required names and `pd_notation` column, with parseable one-component PD data; names identify the intended seeds. | Parsing, column/name presence and representative membership are checked. Literature identification of each name remains the caller's evidence obligation. |
| Q10 | Nonempty targeted searches start from the complete supported crossing-6 input classification and its unresolved components; smaller bounds are supported for the empty-seed fixtures. | The entry checks only `0 <= max_n <= 6`; missing lower-level targets can raise consistency errors. |

## Evidence and boundaries

Definitions are established by
[knot_n5_classification.jl](../../src/knot_n5_classification.jl),
[knot_n6_classification.jl](../../src/knot_n6_classification.jl),
[knot_n6_targeted_bridge.jl](../../src/knot_n6_targeted_bridge.jl) and
[knot_n6_targeted_barrier8.jl](../../src/knot_n6_targeted_barrier8.jl).
Their tests are [N5](../../test/knot_n5_classification.jl),
[N6](../../test/knot_n6_classification.jl),
[bridge](../../test/knot_n6_targeted_bridge.jl) and
[barrier 8](../../test/knot_n6_targeted_barrier8.jl).
The strict label choice retains mirrors separately where the fixtures
distinguish them; traversal and mirror conventions are those of the
[polynomial page](knot-classical-invariants.md#coordinate-choices-and-evidence).

question: open — Separation of two bounded components need not survive a
higher barrier. The expected obstruction is the missing larger-crossing path,
not a new knot type; a verified path would merge the components. Source
identification and general polynomial/classicality claims retain the evidence
boundaries of [K2](knot-gauss.md#numbered-conditions) and
[P3](knot-classical-invariants.md#numbered-conditions).

## Operations

`build_knot_n5_classification`, `build_knot_n6_classification`,
`classify_knot_n5`, `classify_knot_n6`, `run_knot_n6_targeted_bridge`, `build_knot_n6_targeted_bridge`,
`targeted_knot_n6_bridge`, `build_knot_n6_targeted_barrier8`,
`targeted_knot_n6_barrier8`, `run_knot_n6_targeted_barrier8`.
