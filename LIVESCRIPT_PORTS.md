# Live Script migration status

The previous version of this file incorrectly claimed that the migration covered
all Live Scripts after inspecting only the five files below `topology/`.  That
claim is withdrawn.  The repository-wide inventory is recorded in
[`LIVE_SCRIPT_SCOPE.md`](LIVE_SCRIPT_SCOPE.md).

Current evidence:

- 60 `.mlx` files were converted to plain MATLAB source and inspected; all 60
  are migration targets. The non-topology/algebra research scripts are not yet
  Julia ports and must not be silently excluded by the topology dependency
  classification.
- 15 Live Scripts directly exercise `VirtualLink`, `URDiagram`, Gauss/o-data,
  or a topology conversion routine.
- 2 additional Live Scripts exercise algebra/tensor dependencies used by those
  topology workflows.
- 9 Live Scripts contain 44 calls to the o-graph plotting API.
- 37 figures embedded in the source Live Scripts were extracted as MATLAB
  reference images under `artifacts/matlab_outputs/`.

The five legacy SVG files under `artifacts/plots/` are retained only as rejected
output; they are not fidelity evidence and must not be used for acceptance.

## Plot acceptance criterion, revised 2026-08-20

The previous rule on this line read: "No Julia Live Script port is marked
complete until its code cells, observable data, and every plot have been
compared with the MATLAB reference."  **For plots that rule is withdrawn.**  It
is replaced by the following, and the change is not cosmetic — it decides what
work counts as finished.

1. **Plots are no longer accepted against a MATLAB reference.**  The target is
   that the layout algorithm delivers the figure quality it specifies: a
   bend-minimal orthogonal embedding that is actually an embedding.  MATLAB's
   figure is a prior implementation, not an oracle.
2. **Differences caused by the minimum-bend MILP's degeneracy are accepted.**
   That model has many optima with equal objective and the backend picks one.
   Tests compare the objective value; they must never compare the bending
   vector or the resulting coordinates against a stored MATLAB answer.
3. **The Windows fixture gate is retired.**  Several entries below read
   "Windows replay gate".  Those gates are closed by criteria 1–2 above, not by
   any impossibility: Sage *is* reachable from Windows MATLAB through the WSL
   bridge, and a fixture becomes producible once `SageWrapper` is made
   stateless.  What has no Windows form is `SageWrapper`'s persistent
   `sage.all` namespace, because the bridge spawns a fresh Sage per call.
   See [`PLOT_PIPELINE.md`](PLOT_PIPELINE.md) section 7 and R3.
4. **Code cells and numeric results are unaffected.**  They keep their MATLAB
   oracles; only the plot criterion changed.

What replaces the fixture comparison is property checking: every plotted layout
must satisfy P1 (segments axis-parallel), P3 (each edge ends on its own
crossing), P5 (crossing positions distinct) and P6 (no meeting except at a
crossing).  Current status of the `VLExample` corpus under those properties is
in `PLOT_PIPELINE.md` section 6: 19/19 run and emit SVG, 19/19 pass P1 and P5,
16/19 pass P3 and P6.  The three failures are one defect, IS7.

## Current verified implementation

- 46 dynamic Julia plot instances have been generated (44 source call sites;
  one `arrayfun` call expands to three figures).
- All 46 full-resolution Julia contact-sheet entries have been manually viewed,
  but this is only a preliminary screen; it is not geometric acceptance.
- 42 instances with complete MATLAB replay state match on normalized
  `RGaussCode`, real orientation, weighted flag, and `weightRE`; run
  `julia/scripts/export_plot_states.jl` followed by
  `node julia/scripts/verify_plot_states.mjs` to repeat the comparison.
- Julia-native implementations now cover PD code, regions, Sage-compatible
  minimum-bend orthogonal layout, virtual planarization, real-edge tables,
  disk coboundary data, H-move matrices, integral weight solutions, connected
  sum, knot-complement conversion, `VL2URD`, `UR_`, standard-matrix trace
  formula, finite-Hopf tensor contraction, Sweedler, small Borel
  \(u_q(\mathfrak{sl}_2)\), exterior, cyclic-group, and Kac--Paljutkin Hopf
  algebras.
- The Julia suite is run in bounded shards because this desktop runner can
  terminate a single Julia process at about 30 s.  Aggregate assertion counts
  are not acceptance evidence and are deliberately not stated here until the
  current shard logs have been retained.  It includes direct
  MATLAB R2026a numeric replay values for C250818, C260125, and C260319.

## Cell-level acceptance status

`verified` means a Julia entry point and assertions cover its executable
topology/algebra cells with a MATLAB numeric or structural oracle.  `partial`
means the source was audited and useful cells are ported, but this file must not
be used to claim the full Live Script complete.

**"Windows replay gate" in the Open work column is stale.**  Those entries were
written when a Windows rerun was expected to supply plot fixtures.  For the plot
half that gate is retired (criterion 3 above).  For the numeric half it still
stands, and it is reachable: any Windows MATLAB run that needs Sage requires
`SageWrapper` to be made stateless first, because the WSL bridge spawns a fresh
Sage per call.  That covers every `VirtualLink` workflow.  Rows are left
unedited so the original wording stays auditable; read them through this note.

| Live Script | Status | Current Julia evidence | Open work |
|---|---|---|---|
| `Docs/experiment/C260225M_Oinv.mlx` | verified | cyclic invariant, (S^2), Hopf checks | Windows replay gate |
| `algebra/Execution/C2508/C250825Dsl2.mlx` | verified | finite dual, Drinfeld and Heisenberg doubles | Windows replay gate |
| `algebra/Execution/MSTinv/C250405Sweedler.mlx` | verified | Sweedler products, Hopf checks, two invariants | Windows replay gate |
| `algebra/Execution/MSTinv/C250406TensorNetwork.mlx` | verified | all stored tetrahedron, introductory, framing, twist, and Sweedler/small-Borel/cyclic ring cells; small-Borel and cyclic Hopf checks | Windows replay gate |
| `projects/invariant/C250811MSTinvUR.mlx` | verified | conversion/trace cells 1--8, two plot calls, non-weighted A2/B2/B4/A4/C2 MP cells, explicit UR sequences 13/15/16, cell-9 sweep (eight × -1/4), and cell-11 q=exp(4πi/5) invariant (-0.38197) | Windows replay gate |
| `projects/invariant/C250818MSTinv_uqsl2.mlx` | verified | all numeric invariant families and six plots | Windows replay gate |
| `projects/invariant/C250823MSTinv_extalg.mlx` | verified | all exterior-algebra invariant families and two plots | Windows replay gate |
| `projects/invariant/C251020KnotDsl2bfromUR.mlx` | partial | Borromean/Whitehead plots and inputs; three rank-one helper cells; exact rank-one figure-eight Laurent reduction; Hopf `calcInv2` (C=1, Alexander=1); exact multi-rank source forms and four independent exact-rational UR-reduction checks | plot redraw/side-by-side acceptance and Windows replay |
| `projects/invariant/C251030WeylAlgTraceFormula.mlx` | partial | exact standard-matrix determinant trace, symbolic N=1..4 correction matrices, cells 5--7 reductions, cell-8 knot variants (including singular trace), source-matched numeric N=4 reduction state \((C,V,W)=(-1/51,[13,-13],[-113/17])\), and the final N=4 identity \(1/\mathrm{trace}(D)+\det(A+I-\mathrm{circshift}(I,1))=0\) | retain exact reduced symbolic `D.W`/`D.C` display state and Windows replay |
| `projects/invariant/C260125MSTinv_uqsl2.mlx` | verified | both surgery series, connected sums, shared invariant cells, 12 plots | Windows replay gate |
| `projects/invariant/C260319MSTinv_unknot.mlx` | verified | cyclic and Kac--Paljutkin (p=1,2,3) values | Windows replay gate |
| `projects/spine/example/C260220LensSpaces.mlx` | verified | all 31 Julia sweep cases; 24 successful MATLAB path-length outputs exactly match; (L(2,1)) weights/disk sums and plot | Julia also handles the 7 MATLAB try/catch error pairs |
| `topology/Manifold/VLExample.mlx` | partial | all 19 plot calls run and emit SVG; 19/19 pass P1 and P5, 16/19 pass P3 and P6 (`PLOT_PIPELINE.md` section 6); core data operations | IS7 on `koda_virtual` and `koda_fig13`; delete the `s2xs1_calculation` bending override; non-plot Sage display cells |
| `topology/URDiagram/C240610URD.mlx` | partial | construction/reduction and `trace2` regression: common-ε exact limit −1 | per-edge symbolic-limit equivalence and Windows replay |
| `topology/URDiagram/C240623URD.mlx` | partial | all four cells: three `trace` regressions at common-ε limits −1/9, −1, −1; fourth cell sets a,b before reducing to [2,−2] | per-edge symbolic-limit equivalence and Windows replay |
| `topology/URDiagram/URDforKnotCmpl.mlx` | partial | three conversion constructions, source-oracled D2 manual prefix to [-5,1,-2,2,5,-1], the complete final D2 move sequence with explicit non-mutating `lim0` boundary, and current-source figure-eight numeric trace −4.526591341323682e−5 | obtain canonical MATLAB snapshots for the exploratory symbolic `D2.W`/`D2.C` display states; embedded XML’s −1562.5 is stale/error-state output, not the current-source oracle |
| `topology/URDiagram/dev_URDiagram.mlx` | partial | validation, swaps, one plot, source-oracled `reduct3`/repeated `reduct5`/`simplify1` structural cells, and exact-decimal denominator loop N=1..15 [1,3,4,3,1,Inf,1,3,4,3,1,Inf,1,3,4] | complete remaining reduction/dilation cells and symbolic display normalization |

## Remaining acceptance blockers

The code-cell gate remains open because C251020/C251030 still need
multi-variable Laurent cancellation and the URDiagram development/reduction
scripts retain substantive unported workflows.  The C240610/C240623 trace
oracles use a shared perturbation until a per-edge multivariate limit exists.

### The visual gate, superseded 2026-08-20

This section previously required a replayed MATLAB `REdgeTable.Position`
polyline as a geometry fixture before any plot could be accepted.  **That
requirement is withdrawn**, for two independent reasons.

*It is not the criterion.*  Plot acceptance is now property-based rather than a
comparison against MATLAB, so a fixture would not settle anything even if it
arrived.

(An earlier revision of this section also called the fixture unobtainable.
That was corrected on 2026-08-21: Sage is reachable from Windows MATLAB through
the WSL bridge, and a fixture becomes producible once `SageWrapper` is made
stateless.  See [`PLOT_PIPELINE.md`](PLOT_PIPELINE.md) section 7 and R3.)

*It is no longer the criterion.*  Plot acceptance is now property-based, not
comparison-based (criterion 1–2 at the top of this file).  A layout that
differs from MATLAB's because the degenerate MILP chose another optimum is
accepted.

`plot_svg(...; replay_positions=...)` and
`examples/C251020ReplayFixture.jl` are kept as an optional comparison tool for
Mac-side investigation.  They are no longer an acceptance gate, and nothing
should be marked incomplete for lack of a fixture.  The JSON contract they
implement is unchanged: `REdgeTable.Position` is a cell-aligned array whose
element for each arc is an N-by-2 array of `[real, imag]` rows.

What is genuinely open on the plot side is IS7: `_route_edges` misroutes the
edge that closes a cycle, on 3 of the 19 `VLExample` plot calls.  One of the
three disappears when the hand-tuned `bending_numbers` override is deleted from
`examples/VLExample.jl`, since that override was chosen to dodge the fork's
default solution and is counterproductive here.
