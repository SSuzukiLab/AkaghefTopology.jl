# Hempel candidate: rendering intake

The input JSON is a byte-for-byte copy from the literature-gallery task:
`C:\Users\Akaghef\dev\topology-heegaard-literature-gallery\artifacts\heegaard\literature\klebanov-hempel-candidate.hgd.json`.
SHA-256: `3db46b211d410dec8f5a4902557eb31605b87061757ecdf32987756244946549`.
The originating task is `01a0c813-cd0e-7483-be11-7e929d0357f0`.
No curve word, sign, region partition, or source-worktree file was changed.

The supplied attribution is Anna Klebanov, *Gauss diagrams of 3-manifolds*,
[arXiv math/0308055v1](https://arxiv.org/pdf/math/0308055v1),
PDF page 20, Figure 2.16. This task has not re-audited the source figure.
The transcribing task reports source-matched signed words and a candidate
singleton region reconstruction. Those two evidence levels remain distinct.
The local corpus resolver reports this Work as unresolved. No corpus admission
or Inbox mutation is performed within this isolated rendering task.

Only the supplied combinatorial candidate is tested here. A successful layout
would not strengthen the inherited source correspondence or independently
certify the manifold identification. Source crops are not copied or distributed.

`hempel-routing-report.json` records the current bounded computation, with
the exact input hash, terminal ribbon report, finite visibility-grid settings,
route order, first failed wire, and feasible-prefix length for each attempted
configuration. An unsuccessful bounded search is not a nonexistence result.
No failed prefix is exported as a valid diagram.

Reproduction uses `scripts/heegaard/probe_full_curves_hempel.jl` in this
worktree's Julia project. The script calls the existing Heegaard router helpers,
keeps the default 2×g disks and polar endpoint positions fixed, and imposes a
600-second budget after construction. It tries all 52 distinct default
order/weight combinations before coarser-grid alternatives if time permits.
It introduces no new router and does not use Tutte or change the input.

## Observed result

All 80 recorded attempts failed (52 with step 0.14 / margin 1.5, 16 with
step 0.18 / margin 1.5, and 12 with step 0.18 / margin 2.2). The input passes
the existing software checks. Its terminal ribbon has V=4, E=13, F=11,
Euler=2, genus=0, with a boundary-cycle bijection.

The best prefix is 9 of 13 arcs. Default attempts 20 and 46 stop at beta 1
segment 3, from intersection 5 to 7 (port 8 to 26). This is failure of the
existing greedy finite-grid construction, before elastic plotting or atlas
construction. Required next owner: an alpha-cut feasible-seed/router task.
No feasible layout, fixed-100 SVG, full-curve SVG, or Hempel web page is asserted.
The original 7-example atlas is unaffected.

The corpus inventory command succeeded, but its receipt writer rejected this
external library-worktree note with `receipt requires a document inside the
repository`. No receipt was fabricated and no research-main file was created
to bypass that placement restriction. This does not obstruct computation on
the supplied JSON; source correspondence remains owned by the literature task.
