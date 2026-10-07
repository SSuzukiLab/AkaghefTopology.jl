# Unweighted o-graph search database

The o-graph search database separates geometric path discovery from weighted
move admissibility.

## Contract

- Search states are exact raw `VirtualLink` presentations with weights removed.
- Search edges are successful unweighted `MP`, literal `02`, or `PS` moves.
- `02` and `PS` are distinct move kinds. Neither is an alias of the other.
- `H` is never a search edge. It is inserted only while replaying a selected
  geometric path on a weighted source.
- A completed `(state, move kind, selector, enumerator version)` expansion is
  reusable by later runs. Empty neighborhoods and rejected applications remain
  stored results.
- Run provenance is the executing script path plus its automatically captured
  SHA-256 snapshot. Mathematical source provenance does not belong in this
  working projection.

The state identity is intentionally presentation-level. `geometry_hash` is a
secondary lookup key for the exact `real_gauss_code` projection; it is not a
manifold identifier or a canonical equivalence-class proof.

## Workflow

```julia
using Topology

store = open_store("move-graph.sqlite3")
migrate_ograph_search_store!(store)

run_id = begin_ograph_search_run!(store, @__FILE__;
    label="bounded MP/02 exploration",
    scope=Dict("max_depth" => 6, "weights" => "excluded"))

source_state = put_ograph_search_state!(store, source)
expand_ograph_search_state!(store, run_id, source, "02";
    selector=Dict("direction" => "forward"))

paths = query_ograph_search_paths(store, source_state.state_id;
    target_geometry_hash=target_geometry_hash,
    kind_pattern=["02", "MP", "MP", "02", "MP", "MP"])

trial = record_weighted_ograph_path_trial!(store, replay_run_id,
    weighted_source, first(paths);
    h_moves=[(after_step=2, vertex=4, amount=1)])
```

`ograph_search_expansions` records finite neighborhood completeness.
`ograph_search_operations` is the directed geometric graph.
`ograph_search_weighted_trials` records admissible/rejected weighted replays
without introducing weighted states into the graph.

<!-- XQQQ1588700 -->

## Exact H-insertion feasibility

After selecting a geometric MP/literal-02 path, use
`solve_ograph_path_weights(store, weighted_source, path; target_weights=...)`.
It extracts the fixed internal-weight conditions and integer affine transport
of each edge and solves all H coefficients at every stage, including the two
endpoints, in one Smith-normal-form system. `:inconsistent` and
`:rational_only` are obstructions for ALL integer H coefficients on that
specific geometric path, not a bounded coefficient search.

`target_generator=z` instead computes the possible shifts `w -> w+s*z`.
The result's `shift_lattice=(offset,step)` describes the affine integer lattice
of possible `s`. `target_permutation` is an explicit target edge-occurrence
permutation for an independently checked cyclic-root / crossing-label
identification. No general manifold equivalence is inferred from it.

Call `record_ograph_path_constraint_trial!` to retain the equation system and
an obstruction or witness in `ograph_search_constraint_trials`. Always use
`record_weighted_ograph_path_trial!` to replay a positive solution through the
actual weighted move API before reporting it as a move word. Intermediate
framed-weight residuals and endpoint identification need checking as well.

Exact-state closed loops are admitted by `query_ograph_search_paths`: the
source may recur as the final endpoint, but is never expanded again inside
that walk.

`ograph_path_endpoint_permutations` gives all explicit cyclic-root / crossing
label identifications. Such an identification can act on H1, so a nonzero weight
difference alone must not be treated as the desired shift.
`lift_ograph_search_path!` executes the same geometric word on an isomorphic
source, retaining vertex roles through every intermediate move and persisting
each lifted edge. This is the entrypoint for composing selected cycles.

Edge affine relations are cached in `ograph_search_edge_weight_relations`,
separately from geometric graph identity. An interrupted `building` expansion
can resume once every previous owning run has a terminal status. A still-running
owner prevents concurrent regeneration. Retained operations are not discarded.
`open_store` defaults to a 10-second SQLite busy timeout for brief overlapping
reads/writes; applications can override `busy_timeout_ms`.

<!-- XQQQ1000000 -->

Use `finish_ograph_search_run!` with `completed-found`, `completed-none`,
`budget-exhausted`, or `failed`. These statuses describe only the explicit
finite scope stored in `scope_json`; they are not nonexistence proofs.
