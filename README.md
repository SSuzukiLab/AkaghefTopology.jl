# Topology.jl migration workbench

This directory is an in-progress MATLAB-to-Julia migration. It is not yet a
complete or fidelity-verified port.

The dependency boundary is deliberately scoped: the maintained Julia
plot/layout pipeline must not require Sage or Python, but this is **not** a
project-wide policy to replace Sage. Mature general-purpose CAS operations,
especially group computations already provided by Sage/GAP, should be reused
through a thin Julia-facing boundary rather than reimplemented here.

## Algebra units

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
