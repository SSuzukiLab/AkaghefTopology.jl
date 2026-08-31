# Heegaard subsystem implementation status

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
the shared-interior-point regression. Candidate smoothing generation accepts
all three Poincare variants with zero flattened-geometry violations. A full
`Pkg.test()` run reached the existing top-level fixture test and then failed
at `test/runtests.jl:158` on a Windows temporary filename with `Invalid
argument`; therefore this work does not claim that the unrelated full package
suite is green. GAP verification is preserved under
`experiments/heegaard/poincare_sphere/gap/`.
