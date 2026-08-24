# O-graph visual verification log

Inspection date: 2026-08-13.  Status revised 2026-08-20, corrected 2026-08-21.

> **This log is a review record, not an acceptance gate.**  Two things changed
> on 2026-08-20.  Plot acceptance became property-based rather than a
> comparison against MATLAB
> ([`LIVESCRIPT_PORTS.md`](LIVESCRIPT_PORTS.md), criterion 1–2), and the
> Windows fixture this log waits on stopped being a gate
> ([`PLOT_PIPELINE.md`](PLOT_PIPELINE.md) section 7).  A 2026-08-20 revision
> called that fixture unobtainable; that was corrected on 2026-08-21 and the
> gate is retired by decision instead.
>
> The manual pass recorded below also missed real defects.  Of the 19
> `VLExample` instances marked inspected and passing, three contain an edge
> that terminates on the wrong crossing and consequently overlaps another edge
> (IS7).  A five-line endpoint assertion detects all three.  Read every "checked"
> entry here as a human screen for legibility, which it is good at, and not as
> evidence of geometric correctness, which it demonstrably is not.

The reviewer opened the full-resolution contact sheet for every generated plot,
not merely checked that an image file existed.  MATLAB reference sheets and the
corresponding Julia sheets were viewed in the same pass.  This was later found
insufficient: a same-sheet topology check did not prove geometric/directional
fidelity.  The current generated set contains 46 Julia plot instances: 19
`VLExample`, 26 other `VirtualLink` instances, and one `URDiagram` quiver.

| Source Live Script | Julia instances inspected | MATLAB reference inspected | Result |
|---|---:|---|---|
| `topology/Manifold/VLExample.mlx` | 19 | replay v5, 19 captured figures | topology, gaps, arrow direction and real-vertex numbering checked |
| `Docs/experiment/C260225M_Oinv.mlx` | 1 | replay v2 | topology and eight real-arc weights checked |
| `projects/invariant/C250811MSTinvUR.mlx` | 2 | replay v2 plus embedded output | Hopf checked; later Poincare call retained as source-call coverage |
| `projects/invariant/C250818MSTinv_uqsl2.mlx` | 6 | replay v5 | all six topologies and exact replay weights checked |
| `projects/invariant/C250823MSTinv_extalg.mlx` | 2 | embedded output | both source calls checked; MATLAB replay stops earlier in algebra verification |
| `projects/invariant/C251020KnotDsl2bfromUR.mlx` | 2 | replay v2 | **rejected pending redraw**: 2026-08-14 direct side-by-side inspection found both Borromean and Whitehead differ in outer routing, relative vertex placement, and arrow placement. **Rejection withdrawn 2026-08-20**: it was recorded against a MATLAB-comparison criterion that no longer applies, and it waited on a fixture that is no longer a gate (see the header). Both figures pass P1/P3/P5/P6. The routing difference from MATLAB is the degenerate minimum-bend MILP choosing another optimum (`PLOT_PIPELINE.md` IS1), which is now accepted. Whitehead was never explained by IS1 in any case: its bending vector equals the one the fork's default backend returns |
| `projects/invariant/C260125MSTinv_uqsl2.mlx` | 12 | replay v5 | ten call sites including the three-instance `arrayfun` call checked |
| `projects/spine/example/C260220LensSpaces.mlx` | 1 | replay v2 | topology, edge IDs and weights checked |
| `topology/URDiagram/dev_URDiagram.mlx` | 1 | MATLAB source implementation | one leftward horizontal quiver checked |

Checks performed on every sheet:

- no blank/collapsed output or missing component;
- all real crossings have a gap on the under-strand and one real-vertex marker;
- directed arrows follow the replayed Gauss traversal;
- virtual crossings are not rendered as real-vertex dots;
- weighted calls display one label per MATLAB `REdgeTable` real arc, not per
  planar PD segment;
- self-loop arcs use the two distinct slots at the same crossing.

None of these checks covers whether an edge terminates on its own crossing,
which is exactly the defect the pass missed.  Add P3 to any future manual list,
or better, let the automated check carry it.

The contact sheets live under `artifacts/julia_contact_sheets/`; individual SVG
and PNG files live under `artifacts/julia_plots/`.

`plot_svg(...; replay_positions=...)` bypasses the Julia layout solver and draws
a MATLAB `REdgeTable.Position` polyline verbatim.  It is retained as an
optional Mac-side comparison tool.  The sentence that previously stood here —
that a canonical Windows replay must export that field before any plot can be
promoted — is withdrawn, because plot promotion no longer depends on a MATLAB
comparison.  (A 2026-08-20 revision justified the withdrawal by claiming such a
replay was impossible; that was corrected on 2026-08-21.  It is producible once
`SageWrapper` is made stateless — see `PLOT_PIPELINE.md` section 7.)

The portable JSON representation is one N-by-2 `[real, imag]` row array per
`REdgeTable.Position` cell. The C251020 fixture importer verifies that its arc
order exactly matches the Julia `REdgeTable` traversal before it draws either
figure, so an incorrectly ordered replay cannot create a false visual pass.
No plot geometry is inferred from the MATLAB PNG: a temporary hand-traced
fallback was rejected and removed because it would replace canonical drawing
data with an unverified approximation.

## URDiagram diagnostic views (2026-08-14)

MATLAB `URDiagram.plot` renders only a one-arrow quiver.  It is retained as a
source-compatible output, but it is not accepted as evidence of a reduction.
`plot_urd_structure_svg` now emits a separate diagnostic SVG with signed
endpoints, edge pairings, and a bounded weight label.  The C240623 fourth-cell
reduction was regenerated and visually inspected: its remaining edge is `2`
to `-2`; its unsimplified symbolic coefficient is intentionally reported as a
269-character expression rather than expanded across the canvas.  This is a
Mac-local diagnostic check, not Windows canonical verification.
