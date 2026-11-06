# Full alpha/beta surface atlas

## Mandatory display convention

Use the owner-specified [standard spatial display convention](C:/Users/Akaghef/data/21_研究/threads/PJ260922/docs/notes/2026-09-25-marked-paths-integration.md#標準立体表示規約)
for every app use, explanation, implementation change and acceptance check.
This is the normative source for the spatial interpretation; the dated
implementation descriptions below do not override it. An implementation gap
must be reported as such, not used to choose a different convention.

## Combinatorial generation (2026-09-22)

`H.layout_heegaard(h; backend=:pure_julia, router=:planar)` reads only the
diagram. It flag-refines the terminal ribbon map, solves its planar embedding,
truncates vertices to circles, realizes the exact polar-port parameters in
disjoint collars, and transports the entire plane to the usual 2×g placement.
An orientation-preserving global similarity first aligns the computed disks;
a final coordinatewise positive PL map compresses only the exterior, leaving
all disks/ports fixed so that a remote outer-face detour cannot dominate the
viewport. These are changes of the plane, not independently moved strands.
Every accepted placement/transport is geometrically checked. No example name,
saved coordinates, Catlab anchor, or route fixture enters this path. The
existing `router=:auto/:feasible` and saved-layout APIs retain their behavior.

The app's **組合せデータだけで再生成** sends only the current mathematical
input (including D) to `/api/render`. Surface imports with no supplied layout
use this same planar route. The resulting seed goes through the unchanged
fixed-100 physical parameters of `H.plot`; a saved layout is still optional.
Equivalence: the alpha-cut picture is identified up to ambient isotopy of S²
(labelled disks and routes, orientation preserved), which is the same as the
labelled diagram up to orientation-preserving homeomorphism. For a filling
diagram the plane graph of routes and circles is connected, so its rotation
system, fixed by the signed intersection orders, determines the picture; only
the face placed at infinity is an extra choice. Disk positions, disk rotations
and braiding of disks are therefore layout, not data. Marked-surface layers
(seam-word marking, whole-surface cellulation) were moved to the isolated branch
`codex/heegaard-cellulated`. `scripts/heegaard/check_planar_isotopy.py` checks a
saved layout against this statement without using `Topology.Heegaard`: routes
simple and pairwise disjoint outside the open disks, ports in data order, and
drawn faces equal to the regions traced from the signed intersections.
The normal-exit construction now splits route segments along an exact
piecewise-affine map and validates before relaxation. Its previous separate
smooth-map approximations could cross tightly packed routes. No failed plot
is replaced by the old example's coordinates.

`scripts/heegaard/build_generalized.jl` builds AbaloneAuto and PoincareDAuto
from the punctured inputs, with no layout-file reads. It preserves the old
example geometries, adds conversion metadata, and registers the successful
new plots separately. `--conversion-only` keeps Koda conversion independent
of plot generation. The same arithmetic/normal-exit checks apply to imports.

The implementation has no fixed genus. Current geometric scope is a connected
planar terminal ribbon map, or disconnected components sharing one specified
outer region. Other nested disconnected region decorations are rejected, not
silently rearranged. Floating-point degeneracy and transport/plot failure are
reported explicitly; this is not a total algorithm for arbitrary generalized
Heegaard systems. The atlas still requires alpha and beta counts equal to
genus and intersected beta curves. The focused genus-three test checks three
disconnected terminal components, exact port identities, and unchanged input.

The older fixture-scoped demonstrations below remain historical baselines.

Focused checks: `test/heegaard/automatic_layout.jl` (including the public
genus-three surface entrance), `test/heegaard/full_curves.jl`,
`test/heegaard/joining_disk_surface.jl`, and
`test/heegaard/koda_complex.jl`. Browser entry tests:
`scripts/heegaard/check_generalized.py` and `check_poincare_auto.py`, run with
`uv run --with playwright python ...` against port 8741. They submit no layout
to the real API, preserve the complete input/D, test Koda downloads and source
selection, and reject a malformed D. All new figures remain owner-reviewable;
the previous abalone acceptance does not transfer to a different layout.

The dedicated page is `/surface.html` (for example
`http://localhost:8741/surface_Poincare.html`). The original alpha-cut app
remains at `/`. Both display SVG computed by `Topology.Heegaard`.

## Julia / QMD

```julia
using Topology
H = Topology.Heegaard
diagram = H.certify_heegaard(H.load_heegaard("input.hgd.json"))
feasible = H.load_layout("input.hgd.layout.json", diagram)
figure = H.plot(diagram, feasible; view=:surface)
figure                        # image/svg+xml in IJulia / Quarto
H.render(figure; backend=:svg, path="surface.svg")
```

For an already-computed fixed-production plot, use
`H.surface_chart(diagram, fixed_plot)`; this preserves its beta geometry exactly.
The default `H.plot(diagram, feasible)` still returns the fixed-100 alpha-cut plot.
There are no new geometric tuning controls.

The same surface app now accepts `PuncturedHeegaardDiagram` and preserves its
joining disk. The `S³ · abalone + D` example is at `/surface_Abalone.html`;
`/abalone.html` redirects there. Use `H.plot(ph; view=:surface)` or
`H.surface_chart(ph, fixed_plot)`. The app's file-open dialog accepts both
ordinary and punctured data plus an optional saved layout. See
[joining-disk.md](joining-disk.md) for the data contract, chart interpretation,
and current finite verification scope. Abalone was accepted by the owner;
the subsequent `Poincaré + D₀` example at `/surface_PoincareD.html` remains
owner-review pending. It reuses the existing Poincare layout and displays the
referenced conversation's x/y labels and explicit disk boundary.

`scripts/heegaard/build_full_curves.jl` builds the seven baseline examples from
this checkout's versioned `artifacts/heegaard/app/`. It also writes a
checkout-local Jupyter kernel under
`artifacts/heegaard/full_curves/jupyter/`. Select that directory with
`JUPYTER_PATH` when rendering `examples/HeegaardFullCurves.qmd`; it launches
Julia with this checkout as `--project`, without changing any installed kernel.
The local notebook kernel uses `-O0` to fit Quarto's 60-second
kernel-start limit on this machine. The example reuses computed geometry;
fresh mechanical solves were checked with the normal Julia command.

## Shared marked-display fairing (2026-10-26)

The public `H.plot(view::AbstractDict)` overload accepts an explicit
`marked-cut-view/1` embedded layout, with source curve intervals, crossing marks,
optional joining disk, equal circular boundaries on a 2-by-genus grid, and the
native boundary-vertex pairing. `H.render(result)` returns its checked SVG.
This does not infer an absent marking or initialise arbitrary unlaid-out input.
The legacy alpha-cut overloads remain separate and unchanged.

The existing `/api/render` and `surface.html` use this same overload for marked
input. The viewer is no longer tied to six holes, eleven crossings, or a required
disk. The server and page structure were retained. Lickorish reference curves
for genera 1--4 are registered in the existing example list; l curves are shown
as paired cut boundaries, m/n as arcs whose source words close under the pairing.
They are not silently substituted for the Koda input's different marking.

The numerical backend reuses `fair_marked_paths.py`. Before fairing,
`marked_display_pipeline.pair_ports` transports the entire embedding in disjoint
annular collars to realise one reflected circle parameter on each upper/lower
pair. Already consistent ports are retained. Crossings and hole centres/radii
are fixed. D is transported with its source edges. This fixes the earlier
independent boundary-angle placement; only normal tangents were checked before.

Defaults are 12 SVG units for strand clearance and non-port hole clearance.
Only shared intended crossings have a radius-20 clearance exception. Only a
branch incident to a specified port has that port's radius-20 hole-clearance
exception; unrelated curves do not. All curve ends leave their circles normally.
The final Bezier control hulls, source interval order, signed crossings,
paired port positions, and D boundary are checked independently of optimisation.
Continuous-time quadratic segment events check the interpolation from the
transported seed. These are floating-point finite checks, not an optimality or
general completeness theorem. Failure never falls back to a named example.

Verification: Koda retains all 11 crossing coordinates and D boundary incidence;
minimum strand spacing is 14.4278 and non-port hole clearance is 13.4966 SVG units.
Lickorish genera 1, 2, 3, 4 ran the same method. The public Julia overload and the
actual browser import/regenerate/download path were exercised. A separate
diagnostic adds a closed, non-intersecting component with an interior disk; it is
not presented as a manifold example. Missing pairing and off-boundary ports are
negative controls. Repeated rendering must not emit roundoff-sized zero edges.

Entry points are `scripts/heegaard/check_common_marked_app.py` and
`scripts/heegaard/check_common_marked_geometry.py`. Existing fixed-grid examples
can use `marked_display_pipeline.py INPUT.json OUTPUT.json`; parameters and an
optional prepared-complex pairing are explicit CLI arguments. This numerical
backend uses cached NumPy/SciPy/Shapely/NetworkX through uv; it does not fetch
dependencies during a render. A joining disk split by the display cuts is
explicitly unsupported by this layout format. The native prepared complex still
retains such pieces, but this renderer does not discard or reinterpret them.

The owner regarded the current Koda shape as adoptable. The broad lower arc is
retained: the objective favours low bending with weak tension, and no claim is
made that this global bulge is topologically necessary or minimised.

## Meaning of the drawing

The base chart is a sphere minus 2g open disks; the unbounded region includes
infinity. For each alpha curve, a genuine annulus is attached along two
boundary circles. Its red core is one closed alpha curve. Beta arcs continue
from the base through the annuli. The annular beta pieces are radial geodesics,
and the base curves retain the elastic tension/bending/repulsion solution.
Normal endpoint tangents give smooth joins in the indicated collar charts;
no straight stubs are added to an already-relaxed base curve.

With t in cyclic alpha order, annulus angle is -2πt. The plus base boundary
has angle +2πt and attaches to the annulus outer boundary; the minus base
boundary has angle -2πt and attaches to the inner boundary. The base orientation
is -dx∧dy, while the annulus orientation is dx∧dy. These maps reverse the
induced boundary orientations. The Euler characteristic is 2-2g.

An alpha curve runs clockwise in its annulus. Positive beta runs inner to
outer, and negative beta outer to inner. Therefore the local tangent
determinant has the recorded sign. Black signed points are the mathematical
intersections; white points labelled x_i are chart attachments. A sphere arc
from x_i to x_j, with the two corresponding annular segments, implements the
beta successor i -> j. No projection crossings are introduced.

The displayed signed words are the input cyclic words. Curve and intersection
ids are preserved, including ids that differ from positional puncture indices.
The displayed surface is assembled from charts, not embedded as a genus-g
surface in the plane.

## Scope and evidence boundary

The current `HeegaardDiagram` software certification checks incidence,
ribbon/region consistency, and that the reconstructed genus equals the
number of alpha and beta curves. In this baseline,
`compression_boundary_report` uses curve counts; this renderer does not
upgrade it to an independent proof of meridian-system completeness.

The atlas requires a certified object, valid alpha-cut geometry with disjoint
round boundaries, exact polar ports, and smooth normal exits. All beta curves
must have an intersection. Raw or decorated code, failed routing, and
zero-intersection beta components are rejected. Curve systems with counts
different from genus are outside this API. No complete general routing
algorithm, arbitrary immersed diagram renderer, or 3D surface is supplied.

S³, L(5,2), and the inherited Poincaré example are demonstrated, together with
the remaining four baseline lens examples. Poincaré retains the source-code
and candidate-region-decoration provenance and its fixture-scoped feasible
layout; this is not general-router evidence.

The browser check exercises curve selection, signed-intersection tracking,
comparison, and exact SVG download. Julia checks trace beta successors through
both kinds of chart, alpha order, tangent signs, snapshot isolation, rejection
boundaries, and the public fresh-plot entrance.

Earlier verification in this checkout on 2026-09-22: the two focused Julia testsets passed
354 and 11 assertions. Quarto executed all six cells and generated
`examples/HeegaardFullCurves.html`. Its three embedded surface SVGs and the
three QMD-saved SVGs match the web artifacts byte for byte. The browser check
also verifies that a missing fixture shows an error rather than another diagram.
The Quarto checks above predate joining-disk integration; they were not rerun
for that addition. See joining-disk.md for the integration checks and the known
Catlab-dependent failure of the broader Heegaard suite.
