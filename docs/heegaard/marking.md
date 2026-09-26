# Marked Heegaard input

This page describes the earlier seam-word format, which still needs a geometric
witness. For whole-diagram triangulations, zero-intersection/nested curves and
word-driven realization without a saved layout, see [cellulated input](cellulated.md).
The two formats have distinct semantics and are not silently interconverted.

The original `HeegaardDiagram` remains **unmarked**. Its signed cyclic words
and region boundary-cycle partition do not identify a diagram on a fixed
reference surface modulo ambient isotopy. `canonical_hash` is deliberately
unchanged. A layout, its geometric validity, or its fixed-100 solver profile
is not a surface marking.

`MarkedHeegaardDiagram` adds a mathematical field separate from the layout:

```julia
m = H.mark_heegaard(punctured_or_plain_diagram, displayed_layout)
H.save_marked("example.marked.hgd.json", m)
H.marked_hash(m)                 # labelled marked-data identity
H.validate_marking(m, layout)    # throws on a different relative path class
H.plot(m, layout; view=:surface)
```

Selecting a marking is explicit. It selects the supplied representative; it
does **not** recover an unstored historic/source marking. There is no implicit
conversion from unmarked input, and no implicit re-marking when validation
fails. Joining disk D remains bound to the exact labelled underlying diagram.

## Reference surface and mathematical meaning

Version 1 uses the application's standard alpha-cut reference plane: disk
centres `(i,+1)` and `(i,-1)`, radius `.28`, the prescribed labelled polar
ports, and the existing radial annulus attachments (`radial-zero/1`). Alpha
cores, their parameterization/attachment and infinity are fixed. In particular,
moving disks around one another is not an innocent change of coordinates
unless the marking is transported too.

Seam `2i-1` is the upward vertical ray from `(i+1/25,1)` and seam `2i` is the
downward ray from `(i+1/25,-1)`. Ignore the initial pieces inside the disks.
These rays are disjoint and avoid all other disks. A positive letter records
a left-to-right crossing in the reference xy plane, independent of whether
the ray points upward or downward. The xy convention is distinct from the
atlas surface-orientation convention. For every oriented beta segment, save
its exact beta/segment/endpoint IDs and freely reduced signed crossing word.
Words are noncommutative: `[1,2,-1,-2]` is not erased merely because its signed
intersection totals vanish.

The argument for sufficiency here is an **agent derivation, not an
owner-checked proof or a literature-certified implementation theorem**:
cutting this plane-with-disks along all the rays gives a simply connected
fundamental domain. Successive crossings specify a lift through copies of
that domain; adjacent inverse crossings are exactly backtracking. Thus a
reduced word and the fixed endpoints specify the relative homotopy class.
For embedded disjoint proper arcs, the disk/bigon elimination argument gives
relative isotopy; for a system of disjoint labelled arcs it can be applied
successively while keeping already matched arcs fixed. The fixed annulus
attachments then give a definite labelled diagram on the closed reference
surface. D is the selected complementary disk of this same diagram.

This is a **sufficient encoding, not a canonical classification under every
weaker equivalence**. Fixing infinity and all boundary points is stronger
than unmarked isotopy on the closed surface. Distinct marked hashes need not
mean different closed-surface isotopy classes, let alone different manifolds.
The one-twist torus regression does distinguish genuinely different slopes.
There are `2g` rays here because infinity is an additional fixed puncture;
a sphere with `2g` boundary circles can instead be cut into a disk by a tree
of `2g-1` proper arcs. The referenced chat's `2g-3` disk-cutting count was
incorrect. Region/face incidence without a **fixed reference** still does
not remove mapping-class ambiguity.

## Data, realization and rejection

`akaghef.heegaard.marked/1` contains the underlying plain or punctured diagram
and the frame ID, radial-attachment convention, source binding and segment
words. It contains **no route coordinates**. Loading checks syntax, binding,
port incidence, every segment exactly once, and freely reduced valid letters.
It does not claim an arbitrary collection of words has a disjoint realization.

A separately validated `HeegaardLayout` is its realization witness. The app's
`akaghef.heegaard.project/1` download bundles the mathematical object and this
layout in separate fields for one-file re-opening. A supplied witness must
pass ordinary geometric checks **and** realize every specified word. Stored
`report.valid`, plot trace and unmarked core hash alone are not trusted.

Crossing predicates use exact rationals representing the stored coordinates.
For cubic Beziers, rational de Casteljau subdivision continues until each
control bounding box excludes every puncture reference point. Each cubic can
then be replaced by its chord within its convex hull without crossing a
puncture point. This certifies the word independently of sampled flatness.
Tangencies use a half-open crossing rule and cancel inverse letters. The
ordinary disjointness/clearance checks remain the existing numerical checks;
this is not a new exact-arithmetic geometric embedding prover.

Old automatic artifacts contain a few rounding ulps in transported disk and
port coordinates. Frame comparison permits at most `1e-10` coordinate error;
words include straight endpoint connectors to the prescribed standard ports.
This bounded near-identity identification is not permission to move holes or
to relabel a different placement as the same marking.

`plot(marked, witness)` checks the words before and after the common solver.
A different class is rejected, never silently accepted or written back.
`load_layout(path, marked)` also checks the marking. `surface_chart(marked,
layout)` preserves it in chart/SVG provenance. Koda conversion uses the
underlying `(H,D)` and records the marked source hash explicitly.

## Scope and remaining work

There is no fixed genus limit, but this version requires the application's
standard complete-alpha-cut presentation and intersections on every curve.
Arbitrary incomplete alpha systems, positive-genus complementary regions and
zero-intersection components are **not newly supported**. Existing region
decoration only records planar complementary regions.

Word-only routing is **not implemented**. A marked input without its witness
is rejected rather than passed to the unmarked automatic router. Arbitrary
puncture motion, nonzero annulus attachment conventions and changing the
reference frame likewise require explicit marking transport, currently not
implemented. These restrictions concern realization/transport, not loss of
the stored relative classes.

Tests exercise a one-intersection torus with/without a full twist (same
unmarked core, different words), subdivision invariance, noncommutative
reduction, schema round-trip, changed-word rejection, common plotting, and
the existing abalone/Poincare witnesses with unchanged D. The browser check
uses actual mark, download, import and redraw operations and tests rejection
of a tampered word or a missing witness. Owner visual review is separate.
