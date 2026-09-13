# Topology.jl migration workbench

This directory is an in-progress MATLAB-to-Julia migration. It is not yet a
complete or fidelity-verified port.

The dependency boundary is deliberately scoped: the maintained Julia
plot/layout pipeline must not require Sage or Python, but this is **not** a
project-wide policy to replace Sage. Mature general-purpose CAS operations,
especially group computations already provided by Sage/GAP, should be reused
through a thin Julia-facing boundary rather than reimplemented here.

## Algebra units

公開入口の契約は [MST インターフェース](docs/260909_mst_interface.md)、
共通の設計方針は [AkaghefAlgebra](../algebra.jl/docs/260909_interface_design.md) を参照。
以下は現在実装済みの API。

For Taft MST invariants, the exact entry is:

```julia
using Topology
# v is a closed VirtualLink with the intended integer edge weights assigned.
H = TaftAlgebra(9; M=2)
z = Z_MST(v, H; strategy=:qa) # Cyclotomic{9,2}, canonical ζ₉ coefficients
expr = invariant_expression(v, H)
r = simplify_invariant(expr; strategy=:qa)
base = evaluate_invariant(r, 9) # internal Cyclotomic{9}; no M
z1 = normalize_root(base, RootSpec(9,1))
z2 = normalize_root(base, RootSpec(9,2)) # reuse the same exact base value
@assert normalize_root(z2, RootSpec(9,11)) === z2
```

`canonical_value(z)` returns the internal canonical value; `root_spec(z)` and
`specified_root(z)` expose the selected root. `==` compares canonical values
across M for the same N. Arithmetic requires the same N and M; `with_root(z,
RootSpec(N,M))` explicitly relabels a value without changing its coefficients.
To compute a different Galois conjugate, normalize the saved `base` instead.
`TaftAlgebra(:N)` returns a reduced expression until N is fixed. Strategies
`:qa`, `:network`, and `:none` select how much reduction precedes exact evaluation.
`hopf_invariant`, `taft_reduce`, and `taft_evaluate` remain compatibility entrances.

This path uses Julia 1.10 or later and the sibling `AkaghefAlgebra` package
(`../algebra.jl`, recorded in the manifest). It supports `style=:L`; unweighted inputs retain zero edge
weights and the function does not solve or change weights. Unrecognized sums
remain in `r.expression` and are evaluated exactly with a bounded summand budget.
`Z_MST(v, H; details=true)` also returns the base value, reduced expression and timings.

Existing calls with `uqsl2_borel_small(M//N; style=:L)` select this path too.
The exact return type replaces the previous floating result; use `ComplexF64(z)`
at the numerical boundary. `method=:numeric` on the algebra-valued call requests
the legacy tensor contraction. Modified tensors and `style=:R` retain that
numeric path; `method=:exact` rejects them. No example-specific values are stored.

`disk_tetrahedra(v)` supplies the ordered tetrahedral edge data dual to a
nonempty closed normal o-graph, reusing `disk_table` and the reviewed face
recipe. See [the input, orientation and validation contract](docs/specs/ograph-tetrahedra.md).

Finite group algebras can be built without the topology pipeline:

```julia
H = group_algebra(table)
H = group_algebra(elements, multiply)
```

The basis order is the row/column order of `table` or the order of `elements`.
A table entry `table[i,j] == k` means `eᵢeⱼ=eₖ`; the constructor validates
closure, a two-sided identity, associativity, and unique two-sided inverses.
`cyclic_group_algebra(n)` is the compatible cyclic specialization. The stored
`right_integral` / `left_integral` fields follow the repository's existing
coefficient-vector convention; `is_unimodular(H)` compares their one-dimensional
spans rather than their arbitrary normalizations. `is_involutory(H)` checks
`S²=id` directly in the stored basis. These operations do not enumerate
`Hom(π₁,G)` and do not evaluate a topological invariant.

The authoritative scope and current evidence are in
[`LIVE_SCRIPT_SCOPE.md`](LIVE_SCRIPT_SCOPE.md) and
[`LIVESCRIPT_PORTS.md`](LIVESCRIPT_PORTS.md). The original five generated SVGs
are rejected prototypes: they omit the MATLAB/Sage o-graph layout and crossing
semantics and must not be treated as valid output.

Current smoke tests can be run with:

```sh
julia --startup-file=no --pkgimages=no --project=. -e 'using Pkg; Pkg.test(; coverage=false)'
```

On this Windows machine, Application Control enforces user-mode Code
Integrity and can block Julia 1.12 package-image DLLs under
`%USERPROFILE%\\.julia\\compiled\\v1.12`. Keep `--pkgimages=no` explicit in
test and script invocations; it prevents Julia from using or creating those
native package images without deleting the cache or weakening the OS policy.
The repository VS Code settings apply the same flag to Julia extension
sessions. `--compiled-modules=no` is an emergency fallback only: it avoids
incremental compiled-module caches too, but makes this dependency tree
substantially slower to start.

Passing these tests only validates the currently implemented subset. Migration
acceptance additionally requires all 17 topology/dependency Live Scripts and all
44 o-graph plot calls to match the MATLAB structural and visual baselines.
The full-resolution manual review of the 46 dynamic Julia plot instances is
recorded in [`VISUAL_VERIFICATION.md`](VISUAL_VERIFICATION.md).

[`PLOT_PIPELINE.md`](PLOT_PIPELINE.md) is the current state of the plot path:
what runs between a `VirtualLink` and a drawn o-graph on both the MATLAB/Sage
and the Julia side, which stages have been measured to agree, and two
retractions worth reading before trusting an earlier claim. Read it before
touching `src/orthogonal_layout.jl` or `src/svg_plot.jl`. Three things it
settles:

The normative API, call graph, stage contracts, SVG semantics, errors, and
P1–P7 acceptance rules are specified separately in
[`docs/specs/virtual-link-plot.md`](docs/specs/virtual-link-plot.md). New
requirements belong in that specification; audit observations and historical
retractions belong in `PLOT_PIPELINE.md`.

The source-retaining, real-Gauss-centred multi-representation model and the
provisional Julia direction—immutable mathematical values with optional mutable
memoisation—are recorded in
[`docs/specs/virtual-link-data-model.md`](docs/specs/virtual-link-data-model.md).
It is design rationale, not yet a normative Julia API specification.

The data-driven move architecture, MATLAB compatibility boundary, public
`movable`/`move`/`move!` API, capacity reconnection, vertex lifecycle, weight
transport, errors, and support matrix are specified in
[`docs/specs/virtual-link-moves.md`](docs/specs/virtual-link-moves.md).

Reusable move-sequence exploration uses the unweighted SQLite projection in
[`docs/specs/ograph-search-db.md`](docs/specs/ograph-search-db.md). It stores
`MP`, literal `02`, and `PS` as distinct geometric edge kinds; weights and `H`
moves enter only when a queried path is replayed.

- Plots are **no longer accepted by comparison with MATLAB**. The criterion is
  that the layout is a valid bend-minimal orthogonal embedding, checked by
  property. Differences caused by the degenerate minimum-bend MILP are
  accepted.
- `SageWrapper` **has no Windows form** as designed, though Sage is reachable
  there through the WSL bridge. The canonical fixture is obtainable after a
  stateless rewrite; the gates waiting on it are retired by decision, not by
  impossibility.
- The open plot defect is IS7 in `_route_edges`: the edge closing a cycle can
  terminate on the wrong crossing, silently. It affects 3 of the 19
  `VLExample` plot calls.

The parity harness is `scripts/sage/sage_layout_parity.py` together with
`scripts/sage/julia_layout_parity.jl`.

Canonical ownership and the personal/iPad return path are defined in
[`CANONICAL_BOUNDARY.md`](CANONICAL_BOUNDARY.md). In particular, Mac replay
evidence is not Windows Math Harness verification, and personal-surface edits
cannot overwrite either this library or canonical Harness state without an
explicit import, rerun, and review on Windows.
