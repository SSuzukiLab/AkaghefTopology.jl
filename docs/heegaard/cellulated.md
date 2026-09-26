# Marked whole-diagram input

The existing Heegaard app accepts the delivered `version: 1` marked cellulation
document through **データを開く** and `/api/render`. Its new disk-atlas view uses the
same Julia kernel exported as `Topology.Heegaard.MarkedSurfaceDiagrams`; it does
not invoke the old alpha-cut router or infer a route from endpoints.

## Contract

The mathematical document contains the reference surface's named disks and
orientation-reversing side gluings, a full triangulation of each disk, directed
edge words, crossing labels/signs, and optional region marks. A word is interpreted
in that environment, e.g. `P{v1_0 v2_1 v3_2 v4_3} | P{v0_3 v1_4}`. `|` must be
an actual paired boundary gate. Region seed triangle indices are **zero-based**.

Save **数学データ** to preserve that entire document. **図式＋配置** additionally
stores exact rational positions under `akaghef.heegaard.cellulated-project/1`.
The layout schema is `akaghef.heegaard.cellulated-layout/1`; its representation
SHA-256 must match the mathematical document. Loading checks every triangle and
mesh edge again, and requires the fixed canonical chart boundary. Interior vertex
adjustments can be made with `move_vertex` and imported as exact rational positions.
With no layout, `realize` determines coordinates for the saved triangulation,
without changing any path. No geometry is silently synthesized for old inputs.

The hash is a representation checksum, **not an isotopy-class canonical form**.
The intended conditional uniqueness argument is in the supplied
`vendor/marked_surface_diagrams/docs/PROOFS_JA.md`; integration tests do not upgrade
that LLM-authored argument to a checked theorem. Certified output means ideal
polyline centerlines after decimal conversion, not stroke widths or pixels.

## Joining disk and current boundaries

`embedded_closed_disk` and `koda_joining` are different predicates. The former
requires an injective boundary image, not merely a normalized disk complement.
The latter additionally checks the source kernel's connected-arc and outward
orientation conventions. The app never reverses a curve to satisfy that test.
The derived `torus_three_crossings_D` selects triangle `P:68` and reports geometric
disk true, Koda joining false. It is **not** the existing abalone presentation.
The Koda conversion endpoint rejects this geometric-only D. For closed alpha/beta
inputs with both complete cut systems, disk complements, and one oriented joining
disk, `cellulated_punctured` collapses every fine region boundary to a ribbon
boundary and checks a bijection before invoking the existing converter. It retains
the marked document's hashes and intersection labels in the conversion receipt.
Aux/open curves and non-cellular inputs are explicitly unsupported for conversion
(they remain supported by the renderer). It never chooses a disk implicitly from
multiple marks or reorients the curves to make a disk pass.

Zero-intersection curves, disconnected/nested systems and non-disk complement
regions are preserved. The supplied examples include empty genus 0–4 surfaces;
the input format is not genus-1-specific. The exact reference kernel has explicit
limits of 64 patches and 256 vertices per patch. No general-purpose refinement,
smooth atlas transport, or recovery of old Abalone/Poincaré routes is implemented.
Old marked seam-word projects retain their original semantics and renderer.

## Source and reproduction

Imported unchanged from owner-downloaded `heegaard_marked_paths_codex_bundle.zip`,
SHA-256 `59866f50af8107529dfda8f75fda9e58a49b846dfe1583bab721ffde7e6c13dd`.
Source conversation: `6ab50601-7edc-83ee-a004-10d8f1fbb761` (一意性を比較).
Kernel SHA-256: `db69514535674467194d81597db919ce6b4df5b9a11a8e07d02bedc422a66280`.
The vendor directory carries source, documentation, fixtures and original test
oracles; original generated artifacts/overview are retained in the preserved ZIP,
not copied into the app. The supplied README's “Julia not run” describes its
delivery environment, not the later native verification recorded below.

From this checkout with Julia 1.12.6:

```text
julia --startup-file=no -O0 vendor/marked_surface_diagrams/test/runtests.jl
julia --startup-file=no -O0 scripts/heegaard/build_cellulated.jl
julia --startup-file=no --project=. examples/HeegaardCellulated.jl
uv run --no-project --offline --python 3.12 --with playwright python scripts/heegaard/check_cellulated.py
```

The last command uses the running loopback app on port 8741 and headless Edge.
It verifies all supplied cases through project read-back, malformed inputs,
browser save/import and legacy/new-view switching. Generated evidence lives in
`artifacts/heegaard/cellulated/app-verification.json` when that run completes.
The original Python suite was independently rerun (36 tests). The unchanged
Julia suite passed 192 assertions, including exact comparison to the original
12-case Python oracle. No full Topology test-suite claim is made.

## Supplied genus-3 example

`vendor/koda_genus3/` is the unchanged owner-supplied
`koda_genus3_verified_bundle.zip` (SHA-256
`bf6b3c149923bb53e06819a3152d102aa14cfe5e93a56432168121ffda126417`).
All 33 member checksums were verified. Original Python 25 checks and Julia 11
assertions were rerun locally. The source's paper identification is not certified
by these checks. Its manual image transcription remains subject to owner review.

The app entry is `/surface.html?cellulated=KodaGenus3`. The default front/back
projection is a PL transcription witness bound to the whole marked representation,
not an automatically inferred embedding. `build_koda_projection.py` independently
reconstructs its cellulation and checks the full representation hash. The certified
rational disk atlas is available with the adjacent button. SVG downloads follow
the selected view; saved projects retain the mathematical document and exact atlas.
The optional source underlay is provenance, not part of the layout certificate.

The source traversal has D on the right. Repository singular edges follow D on the
left. Consequently the COMPLETE source word must be reversed (cyclic start arbitrary),
with vertex types unchanged. Source word and repository word, respectively:

```text
[-1,-2,3,1,4,-5,-4,5,-3,2]
[1,3,-2,-1,2,-3,5,-4,-5,4]
types = [-1,1,1,1,-1]
```

The repository output has 5 vertices, 10 edges, 6 disk sectors and all three
closedness criteria pass. The source word entered WITHOUT convention transport
has only 4 sectors: that failure was not a counterexample to the transcription.
The older user target instead changes the two visits at vertices 4 and 5 and is
retained as a negative control. The MATLAB plotting reproduction in `VLExample.jl`
and its manually weighted virtual version are not silently rewritten from this
different presentation. The working genus-3 fixture is `HeegaardKodaGenus3.jl`.

Marked conversion uses an explicitly prefixed labelled identity for its lossy
ribbon projection, and separately retains the full marked hashes. This avoids a
costly, unnecessary canonical relabelling search; it does not claim canonical
isotopy equality. Existing unmarked conversion keeps its canonical identity.

```text
uv run --no-project --offline --python 3.12 python scripts/heegaard/build_koda_projection.py
julia --startup-file=no --project=. scripts/heegaard/build_koda_genus3.jl
julia --startup-file=no --project=. test/heegaard/cellulated_conversion.jl
uv run --no-project --offline --python 3.12 --with playwright python scripts/heegaard/check_koda_genus3.py
```

Regenerate the genus-3 entry after rebuilding the original generic fixture gallery.
The actual app checker writes `KodaGenus3.koda.json`, the corrected
`KodaGenus3.repository.odata.json`, and `KodaGenus3.app-verification.json` under
`artifacts/heegaard/cellulated/`.
