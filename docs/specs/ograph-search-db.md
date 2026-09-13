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

Use `finish_ograph_search_run!` with `completed-found`, `completed-none`,
`budget-exhausted`, or `failed`. These statuses describe only the explicit
finite scope stored in `scope_json`; they are not nonexistence proofs.
