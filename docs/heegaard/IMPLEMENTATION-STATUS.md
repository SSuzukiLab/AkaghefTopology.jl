# Heegaard subsystem implementation status

This worktree implements the core contract from the supplied Heegaard
specification:

- `Topology.Heegaard` is a separate namespace under the existing `Topology`
  package.
- Signed cyclic words, intersections, explicit region partitions, layered
  validation, certification, ribbon boundary cycles, surface Euler/genus,
  presentations, exact integer homology, and representation hashes are live.
- `.hgd.json`, `.hgd.layout.json`, and `.hgd2ograph.json` have separate
  serializers and the layout loader checks the core hash.
- Epsilon/R/S moves and an explicitly supplied combinatorial H move return
  fresh values with `MoveCertificate`s.
- Gauss and alpha-cut layout data, semantic SVG/TikZ output, and the typed
  punctured-diagram/flow-spine/OData boundary are live.

The following acceptance layers remain intentionally explicit rather than
silently approximated:

- Catlab/Graphviz/CairoMakie and GAP are not direct core dependencies in this
  worktree. The layout path records a pure-Julia fallback, and `group_probe`
  reports optional GAP unavailability.
- A geometric band-slide search, VirtualLink standardization of OData, and
  fully canonical flow-spine/o-graph conversion are not claimed. The H move
  requires replacement visits, a path record, and the post-slide region
  partition.
- The supplied bundle contains specification and bibliography files, not the
  reviewed Poincare/Hempel fixture transcriptions. No fixture was invented;
  tests use checked synthetic S3, lens, and S2 x S1 records.

The package's full test command is:

```text
julia --startup-file=no --pkgimages=no --project=. -e "using Pkg; Pkg.test(; coverage=false)"
```

Julia 1.12.6 is available on the current Windows machine. Julia 1.10 LTS was
not installed during this run, so that compatibility job remains unexecuted.

Verification at commit `cf4dede`:

- the dedicated Heegaard suite passes all 60 tests;
- the complete package suite passed before rebasing onto upstream `eaec7a6`;
- the post-rebase complete suite is blocked by the pre-existing
  `svg_plot.jl` → `orthogonal_layout.jl` `planarization` keyword mismatch.
  The parent main worktree already contains an uncommitted compatibility
  change for that unrelated issue, and it was not copied into this worktree.
