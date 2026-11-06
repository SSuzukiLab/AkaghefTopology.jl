# Joining disk and abalone smoke example

## Koda executable complex and app integration (2026-09-22)

`to_flow_spine(ph).complex` now contains oriented sector attaching words,
all source surface/compression-disk identities, and the D-edge gluing sequence.
Repeated incidences of one region at an edge are retained. This implements
B = S union the compression disks, followed by deletion of D, from
`K08#^k08-sec-0079` (p.210, §5). The normal-o-graph convention remains the
agent-derived Fig.12 adapter described below.

Each D-edge merges its two remaining oppositely oriented disk sheets.
Self-gluing and incomplete smoothed-edge paths are rejected. Checks require
closed attaching words, three sheets at every edge, induced signs (-,+,+),
and Euler characteristic 1. Supported input: cellular certified diagrams
with one validated joining disk; no implicit standardization is performed.

`to_odata` retains source-intersection/output-vertex maps, output-edge/source-
dart paths, and component runs with the switching-edge o-block count.
JSON bundles contain the intermediate complex, o-data, and certificate.
The common app's **Koda変換を確認** calls `/api/convert` on the displayed input.
No layout or target o-data is supplied. Vertex-source buttons select the
corresponding Heegaard intersection.

`test/heegaard/koda_complex.jl` checks the hand transcript in message
`7c91383a-4724-460a-b098-3107bc03d47d`: label order
`(x3,x2,x1,x7,x5,y5,y3,y2)`, word
`[-1,-2,-3,-4,-5,6,4,7,-7,-8,-6,5,8,3,2,1]`, types
`[r,r,r,l,l,l,r,l]`. Transcript epsilon_K is det(alpha,beta) after
reorientation; repository sigma is det(beta,alpha), so r maps to -1 and l
to +1. The comparison permits this explicit relabelling and cyclic rotation
only, and also checks the four source-curve runs.

Finite results: abalone (V,E,F)=(1,2,2), one block; Poincare (8,16,9), two
blocks; both pass repository C1/C2/C3. These are computational checks, not an
owner-checked proof of the general adapter. The hand word is an independent
oracle; the graph extraction and complex builder are not two independent
implementations of the entire Koda construction.

## Data contract

`JoiningDisk` stores an exact labelled-diagram hash, region number, and the
oriented boundary darts. It is part of `PuncturedHeegaardDiagram`, not a layout
annotation. `joining_disk(H, region)` checks that the region is a disk, its
closure has a simple boundary, and that boundary meets each alpha and beta
curve in exactly one connected arc. The region model has genus-zero pieces;
one boundary component therefore means a disk here.

`make_punctured(H; joining_disk=D)` requires an explicit disk. The previous
implicit region 1 and the curve-word-as-o-data placeholder have been removed.
Relabelling or changing the source invalidates D: its dart IDs are not bound
by the representation-invariant canonical hash. Reselect explicitly after a
transformation. `save_punctured` / `load_punctured` round-trip the core and disk
in `akaghef.heegaard.punctured/1`; existing `.hgd.json` version 1 is unchanged.

## Construction and conventions

The source construction is [K08, §3, p.206] and [K08, §5, pp.210,216–218]:
form B = S union the alpha and beta compressing disks, then remove D. Along
the boundary of D, a singular edge becomes smooth; each corner becomes a
bivalent point of the remaining singular locus and is smoothed.

The implementation records each resulting directed edge as an ordered list
of source darts. Away from D, continue through the opposite incidence of the
same curve family. At a smoothed corner, continue into the other family.
Tracing this permutation gives actual signed Gauss component words, not the
original alpha and beta words.

The local adapter is an agent-derived reading of [K08, Fig.12, p.216] under
K08's **left-hand** ambient convention (§3). Draw S with its given orientation
counter-clockwise; its alpha-side normal then points away from the viewer.
Thus beta is the over-passage and alpha the under-passage in the normal
o-graph. Singular-edge directions agree with the induced boundary of D
(Theorem 5.14(ii)); the oriented slopes of §5 run oppositely. Input curve
orientations are parameter choices and are not silently equated with either.

For input intersection sign epsilon = det(alpha,beta) and disk-induced
direction multipliers da, db, the repository vertex sign is
sigma = det(over,under) = -epsilon * da * db. This follows from the repository
CCW PD frame, starting at incoming-under. Passage signs (+over/-under) and
vertex signs are distinct. Do not reuse this adapter unchanged with the
right-hand convention of another source.

**Evidence.**
- literature: [K08, §3]
- literature: [K08, §5]
- literature: [K08, Fig. 12]

Managed source: `w-74eb606e-3061-4192-b7e6-22abd44ef542`, cite key K08.
Exact Work anchors: `K08#^k08-sec-0073` (definition),
`K08#^k08-sec-0079` (removal construction), `K08#^k08-sec-0085`
(Fig.12), `K08#^k08-sec-0086` and `K08#^k08-sec-0087` (Theorem 5.14).
The managed text remains OCR-unchecked/reading-aid. Load-bearing original
pages 206, 210, 216 and 217–218 are inspected visually for this implementation;
this does not upgrade the Work's corpus status or prove the general adapter.
The consumer resolves K08, but its human figure/section locators are not
mechanically mapped. Receipt creation is unavailable for this external library
worktree (`receipt requires a document inside the repository`); no successful
strict receipt is claimed.

## Finite example

`examples/HeegaardAbalone.jl` draws both curves on R²/Z². Alpha has class (1,0),
beta (0,1), with one cancelling pair; the input is a genus-one S³ splitting.
The simple bigon at intersections 1 and 2 can be cancelled to obtain the
meridian-longitude diagram with one intersection. This is the S³ argument,
not an inference from H1=0 alone.

Both cyclic words are (1 2 3); intersection signs are (+,-,+). The selected
disk has boundary darts [1,8]: 1 ->alpha 2 ->beta-reversed 1. Removing its
boundary edges and smoothing vertices 1 and 2 leaves vertex 3 and paths
[5,12] and [10,3]. The output is `GC [[[1,-1]];[+]]`.

The build script compares with the existing `abalone` object returned by
`examples/VLExample.jl`, without supplying target o-data to the converter.
It also passes the result to `VirtualLink` and checks C1/C2/C3. Source H1 is
not copied into the output as supposed independent evidence.

Run from this worktree (Julia 1.12.6):

```powershell
julia --startup-file=no --compiled-modules=existing --pkgimages=no -O0 --project=. scripts/heegaard/build_abalone.jl
julia --startup-file=no --compiled-modules=existing --pkgimages=no -O0 --project=. test/heegaard/joining_disk.jl
```

Outputs: `artifacts/heegaard/full_curves/abalone/`. Browser entry:
`http://localhost:8741/surface_Abalone.html`, in the existing Heegaard app's
example list. The old `/abalone.html` entry redirects here; it no longer serves
the standalone torus/conversion page. Earlier conversion artifacts are retained.

## Joining disk in the common app

`scripts/heegaard/build_full_curves.jl` loads the persisted punctured input,
routes it with the common bounded router, and uses the existing fixed-100 plot.
It registers `Abalone` alongside the baseline examples. Its cache is
`artifacts/heegaard/full_curves/Abalone.hgd.layout.json`; no torus coordinates
from `HeegaardAbaloneExample.render_svg` enter this rendering path.

```julia
ph = H.load_punctured("S3.punctured.hgd.json")
figure = H.plot(ph; view=:alpha_cut) # 2g puncture disks + beta arcs + D
# Or reuse a saved valid plot without recomputing its geometry:
svg = H.render(ph, fixed_plot; view=:alpha_cut)
# The sphere + annuli atlas remains an explicit alternative:
atlas = H.surface_chart(ph, fixed_plot)
```

`joining_disk_layout.jl` traces the stored boundary darts. Beta intervals reuse
the exact layout Bezier curves (with reversed controls for reverse darts).
Alpha intervals become boundary-circle arcs in the sphere chart and sectors
between the alpha core and the appropriate annular rim. These pieces glue to
the selected D. The source diagram, D, and input layout remain unchanged.
The drawing is derived data; D remains in the mathematical input, not merely
in the SVG. There is no abalone-specific geometry branch.

The app displays D in yellow, its oriented boundary, and a visibility toggle.
Boundary vertices participate in existing intersection tracking. Download
"数学データ" and "再読込用の配置", then use "データを開く" to reopen both.
The same `/api/render` endpoint accepts punctured data in alpha-cut (default)
or surface view and rejects invalid D. Neither view drops D. The browser
defaults to the 2g-puncture view; its atlas toggle reuses the same geometry.

`/surface_KodaGenus3Auto.html` is the genus-3, six-puncture entry. The old
`/surface.html?cellulated=KodaGenus3` URL also opens this entry, not the original
image projection. `scripts/heegaard/check_alpha_cut.py` exercises that URL,
SVG download, D visibility, fixture-free API regeneration and invalid-D refusal.
`refresh_cut_views.jl` only refreshes view payloads from validated saved layouts;
it must not be presented as a fixture-free routing test.

This standard alpha-cut view uses the ordinary input's declared equivalence:
labelled oriented diagrams up to homeomorphism, or ambient isotopy of the
whole picture in S². It does not recover a marking discarded by the earlier
migration or certify isotopy relative to pointwise fixed puncture disks.

The current atlas restrictions still apply: alpha and beta counts equal genus,
valid polar cut layout, and nonempty beta curves. This integration does not
claim arbitrary generalized Heegaard-system support; the Poincare case below
adds one genus-two D example, not an arbitrary-genus routing result.

The focused checks are `test/heegaard/joining_disk_surface.jl` and
`scripts/heegaard/check_abalone.py`. The latter requires a running local app.

## Poincare example from the referenced conversation

The app entry `/surface_PoincareD.html` adds `Poincaré + D₀` while retaining the
old no-D Poincare entry. `examples/HeegaardPoincareJoiningDisk.jl` consumes the
existing baseline and checks all four cyclic words and all intersection signs
against the final Poincare answer in the owner-specified
[HeegaardDiagramProgram conversation](chatgpt-conversation://6a8ece48-b94c-83ee-9018-e1181ed9134b),
message `7c91383a-4724-460a-b098-3107bc03d47d`, sections PG and PF.

Internal IDs 1 through 12 correspond respectively to the conversation labels
`x1,x2,x3,x4,x6,x5,x7,y1,y5,y2,y3,y4`. The numeric topology IDs and saved
layout incidences are not renumbered. The new input stores these names as
`Intersection.source_label`, with `display_source_labels=true` in metadata;
the common SVG and UI use them while selection remains bound to numeric IDs.
Other fixtures retain their existing numeric display. The baseline input and
its geometry are unchanged.

This uses PG's original alpha/beta parameter orientations, not the alpha
reversal performed later for the conversation's o-graph calculation. D₀ has
the same four edges as `x4 ->alpha1 x6 ->beta2^-1 y4 ->alpha2 y1 ->beta1^-1 x4`.
Under this implementation's surface-induced, disk-on-left boundary convention,
the stored boundary runs in reverse:
`x4 ->beta1 y1 ->alpha2^-1 y4 ->beta2 x6 ->alpha1^-1 x4`.
It is region 5 with darts `[31,24,43,8]`. The UI states this orientation
difference explicitly rather than changing signs to force the source sequence.
All four edges, not merely the four corners, are checked when selecting D₀.

The existing validated Poincare fixed-100 layout is reused with `load_layout`
validation. D₀ is rendered by the same `surface_chart`/joining-disk renderer as
Abalone: a bounded sphere-chart region and an inner sector in each of the two
annuli. No Poincare-specific D coordinates or drawing function are introduced.
The baseline layout is fixture-scoped; no new generic-router success is claimed.
The input retains the baseline's candidate-region/source-transcription status.
No new source-PDF certification or Heegaard-to-o-graph conversion is done here.

Build with `scripts/heegaard/build_full_curves.jl`. Reopen the generated
`PoincareD.punctured.hgd.json` and `PoincareD.hgd.layout.json` through the app's
file dialog. Focused checks are `test/heegaard/poincare_joining_disk.jl` and
`scripts/heegaard/check_poincare_disk.py`; the latter also compares the full
SVG after downloading and reopening both files. Only the object-key order
inside JSON metadata is normalized; all geometry and all metadata values must
match. Julia's dictionary insertion order can change when the file is parsed.
On 2026-09-22, the focused Julia check passed 31 assertions. The browser check
passed source-label display, both curve families, all three D chart pieces,
boundary tracing, D visibility, download/re-import SVG equivalence, and return
to the accepted Abalone example, with no page errors. Poincare file re-import
currently takes minutes for Julia's validation; the generated example tab
does not rerun that validation on each selection.

## Limits

Verification on 2026-09-22: the focused joining-disk test has 34 passes. The
build's persisted input produces an exact match to the registered abalone and
passes all three `VirtualLink.is_closed` conditions. The common-atlas disk test
passes 21 assertions, including both embedded bigons and the fresh plot entrance.
The browser check verifies alpha/beta/D rendering, D visibility, boundary and
curve tracing, downloads, exact displayed-SVG equality after file re-import,
invalid-D rejection, baseline navigation, and no page errors. The in-app tab
was visually inspected with D displayed. The owner accepted Abalone with
"OK 次はPoincare". The new Poincare D₀ display remains owner-review pending.
The broader `test/heegaard/runtests.jl` run passed its first four groups (41
checks), then stopped in the Catlab-dependent alpha-cut layout group (24 pass,
2 fail, 1 error: `pure-julia-fallback` instead of `catlab`, missing `catlab_layout`).
That separate renderer failure was not repaired here; the full suite is not green.

This is a one-example smoke comparison, not a proof of the converter on all
Heegaard diagrams. `FlowSpine` holds the singular graph and retained source,
not an independently certified spine CW-complex. Non-cellular inputs and
vertex-free singular components are rejected. The inherited Heegaard
validator's compression check is count-based; it is not a general disk-system
recognizer. General higher-genus geometric realizability, joining unions, weights,
framing, and Koda 2012 Fig.4 are outside this test. `to_heegaard(spine)` only
retrieves the retained source; it does not invert an arbitrary o-graph.
