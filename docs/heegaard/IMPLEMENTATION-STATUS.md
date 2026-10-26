# Heegaard subsystem implementation status

## Constructive marked-operation checkpoint (2026-10-26)

The native whole-cellulation API now constructs empty-triangle and exact
single-vertex isotopies, oriented simplicial actions, barycentric refinement,
quotient-vertex repatching, auxiliary parallel copies, band-sum handle slides
and connected-sum stabilization. Generated inverses reconstruct current routes,
neighborhoods and cell incidence and regenerate the forward result; they do not
return a saved source diagram. This is generated inverse support, not arbitrary
inverse recognition or a complete isotopy decision procedure.

The slide consumes two same-family closed curves and an embedded band word,
constructs a normalized pair-of-pants frontier and computes its crossings.
Default marked-region preservation rejects a changed D. Explicit `follow_seed`
may continue its protected point to a different region, records that change and
revalidates the original claim. Stabilization preserves protected embedded D
outside its support, but rejects an unchanged Koda joining disk whose boundary
cannot meet the newly added curves. Optional display-cut extension constructs a
distinct auxiliary meridian and verifies the extended cut system; it never
silently treats the new alpha as an auxiliary cut. Existing 64-patch and
256-vertex-per-patch bounds remain in force.

`marked_computation_report` binds matrices, generators, relators and curve words
to source elements. `marked_curve_words` additionally exposes oriented fine-edge
words and actual primal-tree basing paths. These words live in the recorded
surface presentation: the retained root-face boundary is a surface relator,
not a freely trivial word.

The independent raw-o-data reverse constructs a new marked Heegaard splitting
from actual oriented tetrahedron gluing and vertex links. It does not consume
retained source metadata or recover the original marking, D, genus or branching.
The retained-source route remains explicitly separate. Comparison distinguishes
full data, fixed-model isotopy, marked cell-model identification, underlying
Heegaard moves and represented-manifold evidence. Independently extracted raw
o-data relations check H1 consistency; matching H1 alone is not a homeomorphism
or group-isomorphism decision.

The lead independently verified finite public cases: moves 49, refinement and
repatching 57+36, parallel copies 14, slides 19+13 (including a nonseparating
genus-two slide on a complete splitting), stabilization 70, cellular computation
32+14+48, reverse construction 46, and comparison 154. Exact Polyhedra 0.8.3
polygon-overlay tests passed 34 and the rational annulus-product tests 21.
Those two helpers are prerequisites only: a first Dehn-twist construction still
exceeds the unchanged per-chart bound (991 versus 256), so no completed twist or
ordered mapping-class-word capability is claimed at this checkpoint.

Initial drawing robustness, complete application integration, mapping-class
operations and whole-book reproduction remain ongoing. The focused results do
not claim that the entire package suite or every valid marked input succeeds.

## Marked input integration (2026-10-26)

The public `load_heegaard` / `save_heegaard` entry now accepts the native marked
document and preserves its whole cellulation, directed words, auxiliary curves,
region marks and metadata. Invalid edits are rejected before opening the output
file. `examples/HeegaardInput.jl` demonstrates an actual curve/incidence correction;
`docs/heegaard/input.md` describes the construction and literature-input contract.
`validate(document; level=:surface)` checks the embedded object and stored D;
full validation separately requires actual closed compression. Structured reports
distinguish proper arcs, positive-genus compression and Koda adapter prerequisites.
The lead independently ran the focused public input and validity tests (48 and
29 checks respectively). These checks establish the stated finite cases, not
general drawing, marked operations or inverse-conversion completeness.

The active viewer checkout now reuses the existing whole-cellulation kernel,
reference constructor and cut preparation from the cellulated checkout.
`marked_topology_report(document)` and the `AbstractDict` overload of
`compression_boundary_report` compute the actual normalized cut components and
cap only the newly created boundaries. Auxiliary reference curves do not split
the mathematical joining disk. The focused public-API test distinguishes a
meridian from a contractible curve with the same genus and family count, handles
redundant parallel curves and original surface boundaries, and retains split D.

The decorated ribbon overload now also computes compression by gluing the
complementary regions across the other family's strips. A genus-one example
with one curve in each family but a sphere-plus-torus compression is rejected,
matching an independently built marked torus. Missing zero-intersection
embedding data is rejected instead of fabricating a genus from curve counts.
Initial geometry generation, full move/convert coverage and executable book
figures remain incomplete.

`surface_group_presentation`, `marked_curve_words` and marked
`fundamental_group_presentation` use a primal tree and dual triangle elimination
on the same cellulation. They retain the generator edges and each elimination.
The finite tests cover reference genera 1--4, disk/annulus, sphere, the meridian/
longitude torus and redundant parallel meridians. Surface group words are
derived computations, not a replacement for the saved whole embedding.
The reference constructor records the managed Work's actual OCR image asset and
hash for Fig. 12.3. The following section describes the earlier implementation.

The genus 1--4 saved reference views now carry that same library provenance.
The existing app's import, regeneration and JSON save preserve it. The accepted
coordinates and SVGs were not changed by this provenance refresh.

For cellular marked diagrams with two minimal complete cut systems, the reused
`cellulated_punctured` adapter checks a bijection from fine regions to ribbon
boundaries before calling the existing Koda converter. Full auxiliary data is
retained. `to_odata(::AbstractDict)` and the existing `/api/convert` also accept a
saved marked view after its source edge words, cut pairing, crossings and entire
oriented D boundary match the retained marked document. The actual Koda view
produced the independently supplied o-data (with the documented traversal
reversal), passed all three closedness checks, and selected the original crossing
from the resulting graph. Its SVG remained unchanged. This is not an inverse
construction, nor support for split D or arbitrary non-cellular conversion.

The old suite exposed a Julia 1.12 world-age failure in the late-loaded Catlab
adapter. Loading Catlab before entering the latest-world callback fixed it;
the existing Heegaard suite and the source-binding negative checks then passed.

Focused verification: `julia --startup-file=no --compiled-modules=existing
--pkgimages=no -O0 --project=. test/heegaard/marked_topology.jl`.
Direct analysis: `scripts/heegaard/analyze_marked.jl INPUT.json REPORT.json`.

## Earlier intersection-only implementation

This worktree implements the core contract from the supplied Heegaard
specification as `Topology.Heegaard`, not as a standalone package.

Implemented and tested:

- signed cyclic words, intersections, explicit region partitions, layered
  validation, certification, ribbon boundary cycles, surface Euler/genus,
  presentations, exact integer homology, and representation hashes;
- separate `.hgd.json`, `.hgd.layout.json`, and `.hgd2ograph.json`
  serializers, with core-hash verification on layout replay;
- epsilon/R/S moves and explicitly supplied combinatorial H moves returning
  fresh values with `MoveCertificate`s;
- Gauss and alpha-cut paired-puncture layouts, semantic SVG/TikZ output, and
  the typed punctured-diagram/flow-spine/OData boundary;
- Catlab as a direct dependency for ordered-port wiring seeds only. Puncture
  placement, terminal ribbon topology, geometric validation, and acceptance
  remain in `Topology.Heegaard`;
- a bounded `feasible_fixed_endpoint_routes` prototype. It consumes the fixed
  punctures, ordered ports, wiring, and optional Catlab anchor hints, searches
  a finite visibility graph, and accepts a result only after
  `validate_layout` reconstructs the beta successor relation and finds zero
  crossings, self-crossings, and puncture intrusions;
- generated review fixtures for `S3`, `L(2,1)`, `L(5,1)`, `L(5,2)`,
  `L(7,1)`, and `L(7,2)`. All six are accepted by the strict layout validator.
  In particular, touching at a shared interior routing-grid vertex is treated
  as an intersection, not as a legal shared endpoint;
- deterministic mechanical relaxation of verified PL routes: topology-preserving
  multiscale shortcuts, densification, length/bending/repulsion dynamics,
  arc-length redistribution, clearance-bounded accepted steps, knot
  reduction, G1 cubic Bezier output, serialization, and final layout validation;
- a Poincare homology-sphere source fixture and a separately certified
  paired-puncture layout fixture. The layout uses Catlab provenance plus a
  fixture-scoped terminal-DCEL corridor certificate and passes strict
  crossing, self-crossing, puncture, endpoint, beta-successor, and layout-hash
  checks.

The following boundaries remain explicit:

- The bounded feasible router is not complete: it uses a finite grid, a fixed
  family of wire orders, and a bounded number of shortest-path attempts.
  Failure is only
  failure to find a solution within those bounds. It is not a general
  fixed-boundary or arbitrary genus-zero terminal-DCEL router, and it uses no
  Tutte embedding, moving punctures, or radial postprocess.
- The Poincare route certificate remains a separate fixture artifact and is
  not evidence of general router completeness. The current bounded search
  does not find the Poincare layout even at 24 attempts, so that fixture still
  requires the checked certificate. This failure is explicit and is not
  replaced by Tutte embedding or puncture motion.
- A geometric band-slide search, VirtualLink standardization of OData, and a
  fully canonical flow-spine/o-graph conversion are not claimed. The H move
  requires replacement visits, a path record, and the post-slide region
  partition.
- The Poincare decorated region partition is a cellular reconstruction
  candidate, not a region decoration transcribed from Klebanov's figure.

The package's full test command is:

```text
julia --startup-file=no --pkgimages=no --project=. -e "using Pkg; Pkg.test(; coverage=false)"
```

Julia 1.12.6 is available on the current Windows machine. Julia 1.10 LTS was
not installed during this run, so that compatibility job remains unexecuted.

Current verification evidence is recorded in the Poincare fixture note and in
the test suite. The focused feasible-router suite passes 32/32 after adding
the shared-interior-point regression. Before/after smoothing artifacts are
generated for Poincare and all six S3/lens review fixtures with zero reported
route crossings. A full
`Pkg.test()` run reached the existing top-level fixture test and then failed
at `test/runtests.jl:158` on a Windows temporary filename with `Invalid
argument`; therefore this work does not claim that the unrelated full package
suite is green. GAP verification is preserved under
`experiments/heegaard/poincare_sphere/gap/`.
