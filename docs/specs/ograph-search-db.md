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
`query_ograph_path_constraint_trials` retrieves stored path/H recipes and
obstructions by source, fixed target weights, and equation status. It returns
composed walks with repeated geometric states as well as simple paths, and
reports whether that exact recipe has a separate weighted replay result.

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

## Selectors and query boundaries

Status: fixed 2026-11-25 Codex (rule a), recording
[`ograph_search_db.jl`](../../src/ograph_search_db.jl). Move names are
case-sensitive strings. Selector keys are converted to strings; unrecognized
keys are retained in expansion identity but ignored by candidate generation.

| Kind | Recognized selector | Default and domain |
|---|---|---|
| `MP` | `vertex_arities` | `[2,3]`; sorted distinct integer-converted entries must belong to `{2,3}`. Ordered distinct vertex tuples are considered. |
| `02` | `direction` | `"both"`; one of `"forward"`, `"inverse"`, `"both"`. Forward candidates use consecutive real-Gauss entries; inverse candidates use ordered distinct vertices. |
| `PS` | `parameters`, `direction` | Parameters `"1"`–`"4"`, default all four; direction as for `02`. Forward candidates use distinct ordered edge-occurrence locations, inverse candidates ordered distinct vertices. |

Candidate generation works on the unweighted projection and deduplicates
candidate descriptors. This is the enumerator's finite neighbourhood; a
completeness theorem for the mathematical calculus does not follow from it.
Unknown move names, including `H`, are rejected. Literal `02` and `PS` retain
their separate selector vocabularies.

Path queries require exactly one target: an exact `state_id` or an exact
real-Gauss `geometry_hash`. Default allowed kinds are `("MP","02")`, minimum
depth `1`, maximum depth `12`, and result limit `1000`. Depth bounds must be
nonnegative and ordered, and the limit positive. A non-`nothing` `kind_pattern`
supplies the exact sequence of `MP`, `02`, `PS` kinds and overrides both depths
with its length after validating the originally supplied bounds. Without a
pattern, `kinds` is the allowed set.

Queries use stored successful edges, sort by depth and then the serialized
operation-ID sequence, and return at most `limit` paths. Repeated states are
excluded except for a final return to the source, which is not expanded again.
With `min_depth=0`, a matching stored source can yield a zero-step path. This
does not enumerate every walk with repeated intermediate states.

The empty-input cases are not uniform: an ordinary no-match plural query gives
an empty vector and the singular query gives `nothing`; an empty allowed-kind
set with no pattern makes the plural query return `nothing`, and the singular
wrapper then errors at `isempty(nothing)`. An empty explicit pattern produces
an invalid SQL `CASE` expression and errors. These are current implementation
boundaries, not additional supported zero-step query forms. To request a
zero-step match, use ordinary nonempty kinds and `min_depth=0`.

## Storage checks and semantic replay

`verify_ograph_search_store` checks the supported search schema version and
SQLite foreign-key consistency, then counts states, runs, expansions,
operations, successful edges and weighted trials. It does not reconstruct
states, authenticate hashes, rerun moves, check enumerator completeness, or
replay weighted trials. Those stronger checks must not be inferred from its
name or from positive counts.

`replay_weighted_ograph_path` requires assigned source weights and, for a
nonempty path, checks its first source identity. H insertions are indexed by
`after_step=0:length(path)` and preserve caller order at each stage. It applies
the candidate descriptors and actual weighted moves to a copy, retaining
intermediate states. It does not authenticate each supplied row against the
database or compare every recorded target ID with the resulting presentation.
Thus callers supply an authentic composable path, and check the final target,
framed residuals and endpoint marking separately. A successful replay is
evidence for the applied word, not for all possible paths.

## Affine caches and endpoint coordinates

The affine edge relation uses real-edge occurrence order, flattened by
component: `C*w=d` is the internal-weight condition and `T*w+b` the transported
target weights. The current relation extractor supports closed MP/literal-02
paths. Integer-path callers must supply finite integral coefficients in the
integer mode; the solver's rounding/conversion operations do not constitute
a framed-weight or coefficient-domain validator.

The database relation cache is keyed by `(operation_id, algorithm_version)`.
A cached JSON relation is decoded without recomputing or proving it. The
optional in-memory `relation_cache` is keyed only by operation ID and trusted
as supplied; callers must keep its store, version and coordinate meanings
consistent. A cache hit is reuse of evidence, not independent verification.

For `ograph_path_endpoint_permutations(target, reference)`, a returned vector
`p` lists target occurrences in reference order: `target_weights[p]` are the
coordinates to compare with reference weights. It represents a cyclic root
and a crossing-label bijection preserving signed passages and vertex signs.
All matching roots are returned in increasing order; two empty words produce
no permutations. Reversal, arbitrary homeomorphisms and an induced homology
marking are outside this identification.

`solve_ograph_path_weights` uses `target_permutation` as a row selection and
does not check that it is a bijection or one of these identifications. The
caller supplies that proof, and a target vector/generator in the selected
coordinates. In contrast, `lift_ograph_search_path!` checks its initial source
permutation against the returned set, then transports vertex roles through
the geometric word. Lifting persists edges incrementally; a later failure can
leave an earlier prefix stored.

These implementation facts are [proved/agent: by construction]. The
[search tests](../../test/ograph_search_db.jl) cover weight-free identities,
literal-02 expansion/reuse, no-match queries and weighted replay; the
[constraint tests](../../test/ograph_path_constraints.jl) compare the affine
result with a weighted move and exercise endpoint permutation, lifting and
cache reuse. The empty-kind/empty-pattern failures and hostile cache contents
lack dedicated regression assertions in those tests; their general handling
remains unverified. The expected failures and trust boundary above follow the
current branches and would be overturned by contrary execution evidence.
