# Topology.jl migration workbench

This directory is an in-progress MATLAB-to-Julia migration. It is not yet a
complete or fidelity-verified port.

The authoritative scope and current evidence are in
[`LIVE_SCRIPT_SCOPE.md`](LIVE_SCRIPT_SCOPE.md) and
[`LIVESCRIPT_PORTS.md`](LIVESCRIPT_PORTS.md). The original five generated SVGs
are rejected prototypes: they omit the MATLAB/Sage o-graph layout and crossing
semantics and must not be treated as valid output.

Current smoke tests can be run with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

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
