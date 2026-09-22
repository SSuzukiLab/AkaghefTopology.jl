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

## Heegaard subsystem

`Topology.Heegaard` stores signed cyclic words and explicit region
decorations separately from layout data. The core path provides
immutable-by-convention code values, layered validation/certification, ribbon
boundary cycles, exact intersection and Smith-normal-form computations,
representation canonicalization, certificate-carrying local moves,
deterministic alpha-cut layout data, semantic SVG/TikZ output, and typed
punctured-diagram/flow-spine/o-data boundaries.

The pure-Julia core does not silently delegate topology to Catlab or a
renderer. `layout_heegaard(; backend=:catlab)` records its pure-Julia fallback
in the backend trace, `group_probe` reports when optional GAP is unavailable,
and non-canonical spine output is marked in its conversion certificate.

### Heegaard plots in QMD/IJulia

`HeegaardLayout` implements the `image/svg+xml` display MIME, so a layout at
the end of a QMD cell renders inline. `Heegaard.plot(diagram, feasible)` applies
the fixed production profile (four review settings fixed at 100). The joint
Julia solver replaces the retired browser-only preview; it does not claim
pixel identity with that preview. There are no public tuning sliders.

```julia
using Topology
const H = Topology.Heegaard

diagram = H.certify_heegaard(H.load_heegaard("diagram.hgd.json"))
feasible = H.load_layout("diagram.hgd.layout.json", diagram)
figure = H.plot(diagram, feasible)
figure
```

Persist the same verified result with the existing APIs:

```julia
H.render(figure; backend=:svg, path="diagram.svg")
H.save_layout("diagram.hgd.layout.json", figure)
```

[`examples/HeegaardPlot.qmd`](examples/HeegaardPlot.qmd) is the executable
Poincaré-sphere example. The MVP intentionally does not expose slider mapping,
presets, or a second figure-management abstraction.

### Local Heegaard application

Build the selected representative/comparison plots and serve the review app:

```powershell
julia --threads=2 --startup-file=no --compiled-modules=existing --pkgimages=no --project=. scripts/heegaard/build_app.jl
python scripts/heegaard/serve.py
```

Open <http://localhost:8731/> for a curated representative gallery. Selection and
the structural differences are explicit in `web/heegaard/gallery-selection.json`.
Stabilizations and connected sums are collapsed construction comparisons, not
evidence of structural diversity. Previously generated parameter files are
retained and their direct URLs still work, but the default manifest lists only
the selection. `--extended` explicitly reindexes the old sweep; its count is
not a success metric. The T³ example
is genuinely genus three: each beta meets two distinct alpha curves, with
geometric incidence `[2 0 2; 2 2 0; 0 2 2]` and computed H₁ = Z³.
Its source is Ozsváth–Stipsicz–Szabó, *Heegaard Floer homology*, Example 2.1.13,
Figure 2.5, pp. 34–35 (2024-08-20 author PDF). This is not a stabilized lens
diagram. The source octahedron is embedded explicitly, then an annular
half-twist and two inverse disk spins put its six punctures in the 2×3 grid.
`scripts/heegaard/build_t3_seed.jl` reproduces the saved feasible seed;
`examples/HeegaardT3.qmd` uses the ordinary `H.plot(diagram, feasible)` API.
The seed is source-specific; it does not make the bounded general router complete.

The builder defaults to the selected examples and reuses current saved plots.
It publishes completed entries incrementally and records failures in
`build-status.json`; `--fresh` recomputes the selected plots. Only the explicit
`--extended` option rebuilds the old parameter sweep.

Individual pages provide paired-puncture plots, synchronized
before/after comparison, pan/zoom, beta-curve and paired-intersection highlighting,
SVG/core/layout downloads, and local JSON import. No data leaves localhost.
The before/after comparison uses a common coordinate scale. The website and
QMD use `Heegaard.render` without an independent browser layout algorithm.

The punctures remain in two rows. The production curve solver uses length,
line-charge spacing and disk clearance, a whole-route preconditioned update,
and local contact sliding. Continuous motion is checked during relaxation;
the final cubic curves are adaptively flattened and validated. Fixed endpoint
contacts are exempt from disk clearance. Stopping is reported honestly; no
global optimality or universal automatic-routing guarantee is claimed.

Every port has a short, fixed outward-normal guard (8% of the disk radius when
the surrounding free space permits). A monotone radial collar moves nearby
strands together to initialize the boundary condition, before joint elastic
relaxation. Length, bending and repulsion then shape the outgoing curve
together; no straight stub is attached after relaxation. The final endpoint angle is
checked against a 15-degree normal cone; the constructed exits are exactly
normal. This constraint belongs to Julia geometry, not SVG styling.

Poincaré starts from the explicitly recorded, fixture-scoped feasible layout.
The original examples retain their saved finite-router layouts. Additional lens
and connected-sum examples use a constructive annulus-family seed from
`scripts/heegaard/gallery_examples.jl`, checked before the same production
relaxation. This family-specific construction is not a general routing algorithm.
The connected-sum decoration joins one complementary region from each summand;
stabilization adds S³ summands. Source transcription and constructed examples are
distinguished in each entry. For a new diagram,
`H.plot(diagram)` tries the bounded router; supply `H.plot(diagram, feasible)`
if a geometric representative is already available. An unsuccessful search
raises an error rather than displaying an invalid diagram. An imported layout
must match the core hash and the standard two-row disk placement.
