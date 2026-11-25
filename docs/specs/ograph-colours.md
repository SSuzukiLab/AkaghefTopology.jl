---
thread: PJ261120
date: 2026-11-22
agent: Grok
---

# Coloured o-graphs

<!-- XQQQ1500000 -->

Scope: the shared coloured o-graph, the Fig. 2.5 edge chart, the conversion of a closed normal o-graph by Proposition 3.2.1, and the face cycles of the standard polyhedron. Figure 2.5 (left) specifies the local polyhedron. The over pair is the upright T, the under pair is the upside-down T, and leg 2 meets the opposite half-edge. Figures 2.5 and 2.7 do not register which adjacent half-edge is leg 1. Candidate B is the calibration derived from Proposition 3.2.1 together with the §2.1 reconstruction and the §3.2 construction. The status is draft, `prop321-calibration-b/draft`.

The reading of the edge chart was checked against the PDF of [BP97](../../../../references/items/book/w-5dde91f6-e933-462c-8257-62ff4bb63ad8/sources/original/BP1997.pdf). A review of that reading is saved at `threads/PJ261120/sources/pro/core-20261121-235637-answer.md`. It is advice. The page `docs/conventions/OGraphColours.md` stays draft.

## Condition identifiers

| ID | BP97 locator | Spec | Function | Test | Review |
|---|---|---|---|---|---|
| O0 | Definition 2.1.4, preamble, PDF p.24 | connected finite quadrivalent graph | `validate_coloured_ograph` | `O0-O3 validator` | applied, ZA1, answer file above |
| O1 | Definition 2.1.4, item 1, PDF p.24 | at least one vertex | `validate_coloured_ograph` | empty graph in `O0-O3 validator` | applied, ZA1, ZA14 |
| O2 | Definition 2.1.4, item 2, PDF p.24 | one opposite over pair at each vertex | `validate_coloured_ograph` | missing, adjacent, and both pairs | applied, ZA1, ZA14 |
| O3 | Definition 2.1.4, item 3, PDF p.24 | one element of Z₃ on each edge; the integer 2 is the colour −1; two congruent records on one edge are that one colour | `validate_coloured_ograph`, `_colour_mod3` | colour 2 passes and canonicalises to −1; congruent duplicates pass; incongruent duplicates fail | applied, ZA1, ZA14 |
| N1 | PDF p.15, normal axioms | at least one vertex, each a marked normal crossing | `validate_normal_colouring` | `closed normal colouring of the example registry` | reading in the same answer; the function was written after it |
| N2 | PDF p.15, normal axioms | directions of opposite edges match through vertices | `validate_normal_colouring`, `normal_orientations` | registry; three small graphs have 4, 2, and 2 N2 orientations, of which 4, 0, and 1 satisfy P2 | same |
| C1, C2, C3 | PDF p.15 | the existing closedness witness, not renamed | `closedness_witness`, copied by `validate_normal_colouring` | registry passes; `koda_fig13` is rejected | not a new reading |
| P1 | Proposition 3.2.1, item 1, PDF p.35 | the same opposite-edge matching as N2, named for the proposition | `validate_normal_colouring` | registry | reading ZA7; implementation after the answer |
| P2 | Proposition 3.2.1, item 2, PDF p.35 | both ends equal gives 0; first below and second above gives 1; first above and second below gives −1 | `prop321_colour` | `Prop. 3.2.1 colour of one step`; registry | reading ZA8 for Figure 3.6 |

`closedness_witness` is the operational C1–C3 replay. It is not an N1 or N2 proof, and its empty-graph branch is not used as N1. `normal_coloured_ograph` rejects the empty presentation before calling it. `validate_normal_colouring` is also called on an empty half-edge graph, on a vertex whose cyclic order has length 3, on one reversed edge of `lens_2_1`, and on a witness whose criteria fail C1, C2, or C3 one at a time. A real `OGraphState` does not produce a false C1 slot, so that witness is built by hand.

## Chart

The twelve physical pairings are the table in `docs/conventions/OGraphColours.md`. They are derived from the four panels of Figure 2.5 (right) and the gluing rule τ^α, and they are tested by `Fig. 2.5 golden pairings`. The PDF does not print those twelve rows. The first reading, ZA3 in `threads/PJ261120/sources/pro/core-20261121-235637-answer.md`, calls the table a direct reading. The condition-id check in `threads/PJ261120/sources/pro/core-20261122-012651.md` calls it a derivation. The rows are unchanged. The convention page is still draft.

The stored colour is the undirected Z₃ label of Definition 2.1.4. Reversing the drawing applies ρ = (1 3) on both frames and keeps α. Freezing the printed numbers and reversing only the bijection gives −α. That calculation is tested, and it is not the stored colour. ρ and the frozen inverse are inferences from Figure 2.5. ZA5 and ZA24 of the first reading are the review of the stored colour. The code follows them.

A numbering shift τ^u on the right and τ^v on the left changes the read colour to α + v − u. The hash does not quotient that gauge, and it does not quotient C or MP.

## Normal adapter

`normal_coloured_ograph` reads an `OGraphState`, a `ClosedNormalOGraph`, or a closed `VirtualLink`. Slots run counter-clockwise from the incoming under-passage. Index +1 uses (−i, +o, −o, +i). Index −1 uses (−i, +i, −o, +o). In both orders the over slots are 2 and 4. The positive Gauss label is the over passage. That join is the repository adapter, still draft.

Colours are listed in `real_edge_table` order. The edge whose head is word position k has index k. The display is `coloured GC [component words; vertex signs; colours]`.

An undirected coloured graph may have 0, 1, or several N2 orientations whose P2 value equals the stored colour. `prop321_orientations` returns that list. On the one-vertex graph with colours 0 and 0, two N2 orientations exist and none satisfy P2. A single orientation failing P2 is not a report that no orientation works. Disc orientations of the spine are not constructed. `normal_orientations` enumerates at most 12 edges, and `coloured_ograph_hash` searches at most 6 vertices. The validators themselves have no such bound. `OGRAPH_FACE_CYCLE_STATUS` is `prop321-calibration-b/draft`.

Shifting one colour by +1 keeps O0–O3, changes the edge gluing, and makes the original orientation fail P2. Homology of the polyhedron is not part of that check.

The nine closed names in the registry test are `s3_abalone`, `lens_2_1`, `lens_3_1_mst2_example23`, `lens_3_2`, `poincare`, `quaternion`, `s2xs1`, `s2xs1_2vertex`, and `odata_variant`. `koda_fig13` is not closed. No sqlite census is included.

## Slot transcription

The crops are `threads/PJ261120/data/bp97-figures/fig2_5.png`, `fig3_4.png`, and `fig3_5.png`, taken because the Work's page images are too small to read arrowheads. Arrowhead direction is the narrow end of the ink span.

Figure 3.4's cross has both horizontal arrowheads pointing left and both vertical arrowheads pointing down. Figure 3.5's cross has both horizontal arrowheads pointing left and both vertical arrowheads pointing up. That is the BP97 row of `docs/specs/ograph-source-conventions.md` §8. The adapter does not store those compass directions. It stores the cyclic ports of `docs/specs/ograph-terminology.md` §6.

Figure 2.5 (left) draws the vertical wall on the left-right singular line and the downward sheet on the transverse line. The three panels of Figure 3.4 are one vertex, so the wall line is the over pair and the transverse line is the under pair. The over T is upright. The under T is upside-down. Leg 2 is the stem, and that sheet contains both opposite rays, so leg 2 meets the opposite half-edge. Legs 1 and 3 meet the two adjacent half-edges and exchange when the T is upside-down, as Figure 2.5 (right) prints.

Figures 2.5 and 2.7 do not register leg 1. Figure 2.5 (right) draws leg 1 upper-right when the edge lies horizontal in front of the reader. No sentence and no number on the left polyhedron identifies that stroke with the counter-clockwise neighbour, the clockwise neighbour, the far sheet, or the near sheet. The caption of Figure 2.7 in [[BP97#^bp97-sec-0550]] says that the orientation of the union of the four disc germs determines the embedding near the vertex. The paragraph that sends a 3-manifold to an o-graph states that, once that disc D lies on the horizontal plane with the usual orientation and the screws are left-handed, the embedding of a neighbourhood of the vertex is essentially unique. That sentence fixes the local embedding. It does not name the repository's counter-clockwise slot. Figure 2.7 prints no leg numbers.

Candidate A, with the repository order counter-clockwise from the incoming under-passage: an over slot has `(leg 1, leg 2, leg 3) = (ccw, opposite, cw)`, and an under slot has `(cw, opposite, ccw)`. Candidate B swaps legs 1 and 3 on every slot. Transport stays `ρ∘τ^α`. Because `ρτρ = τ^{−1}`, candidate B with that transport is the same gluing as candidate A with transport `ρ∘τ^{−α}`.

The adopted calibration is B, derived from [[BP97#^bp97-sec-0566]] Proposition 3.2.1 together with the reconstruction in §2.1 and the construction in §3.2. The proposition sentence by itself is not read as the whole spine claim. A spine of a closed manifold has `n + 1` discs, and the same count is `χ(P) = 1`. Built by `face_cycles` from the P2 colours under candidate A, the registered closed presentation `s2xs1_2vertex` has `n = 2`, five discs, and `χ = 3`. Its Gauss code is `[{-1},{-2},{+2},{+1}]` and its vertex signs are `[{-1},{+1}]`. The manifold name is not part of the refutation. All `6n = 12` germs lie on face cycles. That coverage checks the walk. It does not prove that the cycles are the discs of §2.1. The disc count differs from the closed spine when the walk counts the discs of candidate A's gluing. Among the registered closed examples with at most two vertices, this is the smallest failure of `d = n + 1`. `lens_3_1_mst2_example23` has the same vertex count, three discs, and `χ = 1`, so it is not the refuting example. At larger vertex counts, `lens_3_2` has two discs and `χ = −1`, and `poincare` has two discs and `χ = −3`. These counts are a recomputation by `face_cycles`.

Fixing leg 2 as the opposite half-edge does not by itself leave a global pair of calibrations. The dichotomy uses a further draft premise: one local choice of leg 1 propagates through the connections inside a vertex, both vertex indices, the reversal of the T, and the endpoint reversal ρ. Under that premise, the draft entries below, and the requirement that the input is a closed normal presentation independently of the calibration, candidate B is the remaining calibration for this stored colour and this transport `ρ∘τ^α`. A convention that also changes the colour or the transport lies outside this elimination. No sentence of Figure 2.5 or Figure 2.7 names leg 1 on a slot.

The derivation depends on these draft entries:

- P2 reads above as the over-passage.
- The Figure 3.4 and Figure 3.5 slot orders of the adapter.
- The twelve Figure 2.5 pairings and the transport `ρ∘τ^α` with the stored colour kept.
- The over pair is the upright T, the under pair is the upside-down T, and leg 2 meets the opposite half-edge.
- One local choice of leg 1 propagates as above.
- The face walk counts the discs of that calibration's gluing. Germ coverage is not that geometric identification.
- The input is a closed normal presentation independently of the calibration. When the closedness of `s2xs1_2vertex` is read from C1–C3 of `closedness_witness`, C3 uses the row count of `disk_table`, and `ograph-terminology.md` §11 does not re-derive the identification of those rows with the circuits of Figure 1.3. That identification is a dependency of the closedness claim. Equality of candidate B's corner partition with `disk_table` is a later check, and it is a different claim.
- The disc counts are the recomputation by `face_cycles`.

Agreement of H₁ with a known manifold, and equality with `disk_table`, are checks run after the calibration. They are not the reason B is adopted. `OGRAPH_FACE_CYCLE_STATUS` is `prop321-calibration-b/draft`. `face_cycles` builds the discs. The default calibration is `:b`.

The slot-transcription check is `threads/PJ261120/sources/pro/core-20261122-033946.md`. It is advice. ZA1 leaves both candidates open on the figures alone. ZA2 names the missing registration as one Z₂ choice: whether the side sheet of leg 1 follows J_D h or −J_D h, where J_D h is the counter-clockwise neighbour of half-edge h inside the oriented disc D. ZA16 confirms the algebra relating the two charts. ZA18 keeps the local embedding fixed by D and the left-handed screw. That review does not apply Proposition 3.2.1, and it is not the adoption above.

The calibration check is `threads/PJ261120/sources/pro/core-20261122-045101.md`. It is advice. The file records model GPT-6 Pro and thinking effort Pro, both verified. ZA1 reads the spine claim as strong enough to refute a leg calibration by a wrong disc count, once §2.1 and §3.2 are included, and reads `d = n + 1` and `χ = 1` as one count. ZA2 accepts the elimination of A when the correct calibration is already inside {A, B}. ZA3 finds no sentence that places leg 1 on the counter-clockwise slot. ZA5–ZA8 are the dependencies added above: closedness of the input, correctness of the disc count, coherent propagation of one leg-1 choice, and a recomputed count. ZA9–ZA12 are the limits on the wording: the spine sentence is the combined reading, the dichotomy needs the propagation premise, germ coverage is not geometric correctness, and uniqueness is for this colour and this transport.

## Face cycles

`face_cycles` partitions the `6n` sheet germs. The default calibration is B. `calibration=:a` reproduces the refuted chart. On every registered closed example the corner partition is compared with `disk_table` germ by germ, as a sorted partition of corners `(vertex, slot_lo, slot_hi)` with multiplicity kept. Under B the partitions are equal on all nine names. Under A they are equal on none. That equality is a finding after the adoption. It is not removed by editing B, and it is not the reason B is adopted. `thickening_boundary` builds the boundary of a regular neighbourhood from the corner link of the dual face gluings. `cellular_homology` is the cellular H₁ of the same face cycles. The twelve pairings stay as they are. The registry colour tuples are compared with `prop321_colour` on the same positive-label reading the adapter uses. C1–C3 in this layer copy `closedness_witness`. The separate negative test feeds one false criterion at a time. `koda_fig13` remains an aggregate closedness rejection before `validate_normal_colouring`.

The positive orientation of the tetrahedron at a vertex is the cyclic order of its four rays, begun at either ray of the over pair. The two over rays are opposite. The order begun at the other over ray exchanges the two over rays and the two under rays, and that double transposition is the same orientation. Write the positive order as `(v0, v1, v2, v3)`. The face opposite `vi` carries the induced orientation `(−1)^i` of the other three rays, taken in the order they have in this tuple. The link triangle of the corner `vi` is that induced order. A dual face pair transports one face's induced order by `ρ∘τ^α`. The pair reverses orientation when the image is the opposite cycle of the other face's induced order. [[BP97#^bp97-sec-0557]] states that the faces of an ideal triangulation are glued by orientation-reversing simplicial maps. The sign above is the orientation this repository defines from the cyclic order and the over pair. `DualFacePairing` writes the end with the smaller `(vertex, slot)` first. Exchanging the two ends keeps the reversal flag. A signed edge index in a face cycle numbers the singular edge from the smaller half-edge id. `cellular_homology` and `coloured_pi1` share that numbering, and the homology group is the Smith form of the matrix.

The test `face cycles of the adopted calibration` prints the mandatory counts. On each of the nine converted closed examples under B: `d = n + 1`, `χ(P) = 1`, one boundary component of Euler characteristic 2, `2n` orientation-reversing dual face pairs, dual edge orbits equal to the face-cycle partition, and H₁ equal to `ograph_h1` of the normal pipeline. The named groups are the trivial group for `s3_abalone` and `poincare`, Z/2 for `lens_2_1`, Z/3 for `lens_3_1_mst2_example23` and `lens_3_2`, Z/2×Z/2 for `quaternion`, and Z for `s2xs1` and `s2xs1_2vertex`. `odata_variant` matches the pipeline's trivial group and is not given a manifold name. Cyclic rotation of the vertex origins, together with the over slots, and the exchange of the two vertices of `lens_2_1`, leave these invariants unchanged. The test `half-edge renaming leaves the public surface unchanged` applies a bijection of the half-edge ids, a renumbering of the vertices, and a rotation of the slots, together, to the nine closed examples and to colourings outside the closed normal branch. The corner partition, transported by the vertex renumbering and the slot rotation, the boundary components, the reversal flags, transported the same way, and H₁ stay fixed. On each of those graphs the component Euler characteristics of `thickening_boundary` equal those of the Section 2.2 ribbon `OGraphColourBoundary.boundary_from_coloured` ([[BP97#^bp97-sec-0551]]). A uniform shift of the leg numbers that keeps the stored colour fixed is a different operation from the colour gauge `α' = α + v − u`.

Off the closed-normal branch, `d = n + 1` is not required. Under B, the one-vertex graph with mate `(2, 1, 4, 3)`, over pair `(2, 4)`, and colours `(0, 0)` has two discs, lengths 2 and 4, `χ = 1`, one boundary sphere, and H₁ = Z/4. Colours `(1, 1)` have the same counts and a different corner partition, with the same H₁. Setting the stored colour −1 of `s3_abalone` to +1 produces that second coloured graph. Disc orientations of the converse of Proposition 3.2.1 are not constructed.

## Canonical form and serialization

This section records the current representation of `canonical_coloured_ograph`,
`coloured_ograph_serialization` and `coloured_ograph_hash` in
[`ograph_colours.jl`](../../src/ograph_colours.jl). The domain is O0–O3 with
at most six vertices. The representation choices are fixed 2026-11-25 Codex
(rule a: the existing implementation and its canonical-hash tests).

For each permutation of the vertices and each independent cyclic rotation of
its four half-edges, assign half-edge labels `1:4n` consecutively in that order.
An opposite over pair `(1,3)` has over-code `1`, and `(2,4)` has over-code `2`.
Record the mate tuple in increasing half-edge order. Record one colour per
edge in increasing order of its smaller endpoint, using the representative
`-1`, `0` or `1` of its Z₃ value. Congruent duplicate colour records contribute
one value. The canonical representative is the least tuple
`(convention_id, over_codes, mate, colours)` in Julia's lexicographic tuple
order. Its vertex tuples are consecutive blocks of four labels, its over
positions are the pair encoded above, and its colour record is attached to
the smaller endpoint. Ties describe the same returned fields.

This fixes half-edge and vertex relabelling, cyclic origins, colour-record
endpoint choice and integer representatives of Z₃. It keeps cyclic orientation,
the convention identifier and the stored colour class. It does not quotient
the colour gauge described in [Chart](#chart), C-moves, MP-moves or reflection.
The reconstruction interpretation is the one in
[OGraphColours](../conventions/OGraphColours.md); choosing a representative adds
no geometric certification.

The current serialization is UTF-8 text, with exactly these literal separators
and no terminal newline:

```text
ColouredOGraph;convention=<identifier>;n=<n>;over_positions=<pairs>;mate=<mates>;colours=<records>
```

Here `<identifier>` is the canonical graph's stored `convention_id`, `<n>` is
the decimal vertex count, `<pairs>` is the comma-separated list of over pairs
written as `1+3` or `2+4`, `<mates>` is the comma-separated decimal mate tuple,
and `<records>` lists `halfedge:colour` pairs separated by commas in increasing
smaller-endpoint order. Decimal integers use Julia's ordinary base-10 string
form, with a minus sign on negative colours. The constructor strips surrounding
whitespace from the identifier; serialization inserts the stored identifier
verbatim, without escaping or length-prefixing. Arbitrary identifiers therefore
have no promised unambiguous inverse grammar. There is no public inverse parser.
The hash is the lowercase hexadecimal SHA-256 digest of these UTF-8 bytes.

Compatibility policy: these bytes describe the current implementation, not a
versioned interchange format. No cross-version byte or hash compatibility is
promised. A consumer needing persistent interchange must retain the serializer
revision and convention identifier; the identifier alone does not version this
grammar. The [canonical-hash tests](../../test/ograph_colours.jl) check the
identified changes and colour-sensitive distinctions. They do not supply a
cross-version golden-byte fixture or an arbitrary-identifier decoding guarantee.

## Boundary diagnostic flags

This section owns the flags of `SpineBoundaryComponent` returned by
[`thickening_boundary`](../../src/ograph_face_cycles.jl). They refer to the
corner-link model in [Face cycles](#face-cycles). The distinct ribbon model and
its component counts are specified in
[Boundary of a coloured o-graph](ograph-colour-boundary.md#what-is-counted).
The face orientation used here is fixed on
[OGraphColours](../conventions/OGraphColours.md#四面体と-link-triangle-の向き).

`orientable` is true exactly when every dual face gluing incident to that
corner-link component reverses the prescribed local orientation. This tests
compatibility of those prescribed orientations; it does not search for a new
choice of orientation on each triangle. A false flag is therefore not a
certificate that the underlying surface admits no orientation. The geometric
surface interpretation still requires the reconstruction assumptions above.

`closed_surface` is set to true for every returned component after the side
pairing, involution, corner-membership and vertex-walk consistency checks
succeed. It is not a separately computed recognition result for an arbitrary
surface. It is neither the closed-normal conditions C1–C3 nor a statement that
the corresponding 3-manifold has empty boundary. The public record constructor
checks neither flag and permits arbitrary stored Boolean values.

For the slot model, the component counts give `euler = vertex_count -
edge_count + face_count`; `euler_sum` is their sum. Equivalence with the
[ribbon model](ograph-colour-boundary.md#checks) is checked on the examples
listed in [Face cycles](#face-cycles), not established here for every input.
The [invariant specification](ograph-colour-invariants.md#c-move-and-mp-move)
records move images with nonreversing face pairs and failed boundary sanity
checks. These flags do not remove that limitation.

<!-- XQQQ1000000 -->
