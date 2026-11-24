---
thread: PJ261120
date: 2026-11-22
agent: Grok
---

# Oriented branchings and the converse of Proposition 3.2.1

<!-- XQQQ1500000 -->

Scope: from a coloured o-graph, decide whether some straight-through orientation satisfies the colour rule of Proposition 3.2.1, and recover the normal diagram when one does. Disc orientations of the spine are not constructed. The matching from the four ports of a vertex to the four T's in Figure 2.5 (left) is not printed.

The reading below was checked against the PDF of [BP97](../../../../references/items/book/w-5dde91f6-e933-462c-8257-62ff4bb63ad8/sources/original/BP1997.pdf). Printed page n sits near PDF page n+9. Each locator is an active block anchor in the Work. The Work face is an unchecked reading aid, so the anchor only names the block; the sentence was read on the PDF. A review of that reading is saved at `threads/PJ261120/sources/pro/branching-20261122-000145-answer.md`. It is advice. Each row below says whether that review was applied. A second review, of this table against the functions and the tests, is saved at `threads/PJ261120/sources/pro/branching-20261122-013056.md`. It is advice. Its model-selection and thinking-effort lines both read `verified=yes`. The locator column was confirmed there. The code and this note then took up the mismatches it named.

## Condition identifiers

| ID | Statement | Locator | How read | In code | In tests | Review |
|---|---|---|---|---|---|---|
| O0 | Connected, finite, quadrivalent | [[BP97#^bp97-sec-0549]], Definition 2.1.4 preamble, PDF p.24 | text | `ograph_conditions` | disconnected pair of loops; a missing port | applied, ZA6 |
| O1 | At least one vertex | [[BP97#^bp97-sec-0549]], Definition 2.1.4 item 1, PDF p.24 | text | `ograph_conditions` | empty graph | applied, ZA6 |
| O2 | Normal crossing with one opposite over pair | [[BP97#^bp97-sec-0549]], Definition 2.1.4 item 2, PDF p.24 | text | `ograph_conditions` | round trip; an adjacent over pair fails O2 | applied, ZA6 |
| O3 | One element of Z₃ on each edge | [[BP97#^bp97-sec-0549]], Definition 2.1.4 item 3, PDF p.24 | text | `colour_mod3`, `ograph_conditions` | the integer 2 is accepted; a missing colour is not | applied, ZA6 |
| P1 | Edge directions match across a marked vertex | [[BP97#^bp97-sec-0566]], Proposition 3.2.1 item 1, PDF p.35 | text | `p1_orientations` | round trip; two-circuit diagram | applied, ZA11, ZA14, ZA47 |
| P2a | Both ends above, or both below, gives colour 0 | [[BP97#^bp97-sec-0566]], Proposition 3.2.1 item 2(a), PDF p.35 | text | `p2_colour` | direct values | applied, ZA11, ZA23 |
| P2b | First below and second above gives colour 1 | [[BP97#^bp97-sec-0566]], Proposition 3.2.1 item 2(b), PDF p.35 | text | `p2_colour` | direct values; Figure 3.6 panel 3 | applied, ZA23 |
| P2c | First above and second below gives colour −1 | [[BP97#^bp97-sec-0566]], Proposition 3.2.1 item 2(c), PDF p.35 | text | `p2_colour` | direct values; Figure 3.6 panel 4 | applied, ZA23 |
| B1 | Oriented branching: disc orientations so that no edge is induced the same way by all three germs | [[BP97#^bp97-sec-0565]], Corollary 3.1.7, PDF p.35 | text | not constructed; `compatible_branching` records `discs = :unread_fig25` | the flag is asserted | applied, ZA9, ZA77; construction withheld |
| IDX | Index +1 is Figure 3.4; index −1 is Figure 3.5 | [[BP97#^bp97-sec-0567]], sentence after the proof of Proposition 3.2.1, PDF p.36 | text and figure | `vertex_index`, `port_slot` | round trip, including a vertex of index −1; two-circuit index change | applied, ZA16–ZA18, ZA49–ZA51 |

P1 is the same geometric matching as N2 in [[BP97#^bp97-sec-0536]], PDF p.15. The check is named P1. N1 and C1–C3 are not validity conditions for a general coloured o-graph. A broken opposite pairing is reported as O2, not as a failure of P1. The vertex signs use the BP97 entry of `docs/specs/ograph-source-conventions.md`: Figure 3.4 is type +, Figure 3.5 is type −, the strand direction matches MST3, and the letters l and r are not used.

## Colours

Above means the over passage of the marked crossing, and below means the under passage. It is not a vertical coordinate in the figure. The positive Gauss label is the repository's name for an over passage; that name is fixed in the o-graph terminology note, and the figures show the over strand as the horizontal one.

Definition 2.1.4 attaches one Z₃ element to the undirected edge. Proposition 3.2.1 then reads that element from an orientation: with the first and second ends named by the arrow, the three cases are P2a–P2c. Reversing the arrow swaps the ends, so the value required by P2 changes sign. The stored colour does not change with the candidate arrow. A nonzero colour therefore selects one direction of a mixed circuit. The negated colour tuple selects the other direction.

Evaluating a stored colour as an odd function of the candidate direction would make both directions pass. That reading was not kept. The same conclusion is ZA36–ZA39 of the saved review.

On one straight-through circuit the signed colours sum to 0 in Z₃. A one-circuit diagram meets each vertex once over and once under, so the all-zero tuple does not satisfy P2. An all-zero compatible colouring needs at least two circuits, and each homogeneous circuit may be reversed on its own.

For a fixed marked graph, let r be the number of straight-through circuits that use both an over passage and an under passage. The count of branching colourings is an inference from Proposition 3.2.1, not a sentence of the book. Along one circuit the P2 value is the difference of the above-indicators, so the colours sum to 0 in Z₃. Reversing the circuit negates every colour on it. A homogeneous circuit has every colour 0, and reversal leaves that tuple unchanged. A mixed circuit has a nonzero colour, so reversal produces a different tuple. The circuits use disjoint edges, so r mixed circuits give 2^r tuples. The other tuples in {−1, 0, 1}^E still satisfy O0–O3. Shifting one edge by +1 in Z₃, starting from a tuple that already satisfies P2, makes that circuit sum 1, so the new tuple comes from no branching. A shift of a tuple that does not satisfy P2 is not claimed here. These counts do not use Figure 2.5. `branching_colour_census` enumerates the tuples; it does not evaluate a closed formula, and it stops above 10 edges.

## Figures

Figure 3.4, middle cross, read from the figure: the horizontal strand points left and the vertical strand points down. Figure 3.5, middle cross: horizontal left, vertical up. The curves on the right of each figure are the boundary of a neighbourhood of the singular set, seen from above. Their vertical arrow is opposite the vertical strand. They are not passage signs.

Figure 3.6, read from the figure, singular-edge arrow upward, so the bottom end is first and the top end is second:

| Panel | Bottom | Top | Colour | Curves, from the viewer's left |
|---|---|---|---|---|
| 1 | + | + | 0 | three straight curves, down, down, up |
| 2 | − | − | 0 | same arrows; the middle curve is drawn dashed |
| 3 | − | + | 1 | left and middle curves cross; the right curve is straight and points up |
| 4 | + | − | −1 | same crossing pattern and the same three arrows |

The dashed segments are hidden parts of the projected boundary curves. They are not a mark of the preferred germ. The three arrows are not the three induced directions on the singular edge, so they are not fed into a majority vote. The edge direction of a normal o-graph is the sign of the sum of the three incidence signs of ∂β. That sum is not read off these arrows.

Figure 2.5 (right), read from the figure in [[BP97#^bp97-sec-0549]]. Right-way-up: legs (1, 2, 3) sit up-right, up, and down-left. Upside-down: down-left, down, and up-right. The four panels are the four combinations of those two T's. The glueing sentence in [[BP97#^bp97-sec-0550]] sends right leg i to left leg (1 2 3)^α(i). The horizontal segment in the figure has no arrow. Nothing in Figures 2.5, 2.7, or 3.4–3.6 labels which of those legs is which curve of Figure 3.6, so the choice of which end is the right end is not fixed here. Figures 3.4–3.6 are in [[BP97#^bp97-sec-0567]].

## What is recovered

`coloured_ograph_from_normal` builds the undirected coloured graph of a one-circuit normal presentation. Slots run counter-clockwise. For index +1 they are the incoming under passage, the outgoing over passage, the outgoing under passage, and the incoming over passage. For index −1 the two over slots are exchanged. The colour of each step is P2 in the direction of the Gauss word.

`compatible_orientations` tries both directions of every straight-through circuit and keeps those whose P2 values equal the stored colours. `recovered_normal_diagrams` returns one oriented word per circuit, plus the vertex indices of that orientation. `recovered_normal_states` keeps only the diagrams with a single circuit, because `OGraphState` is a one-circuit presentation. `compatible_branching` stores both: `diagrams` for every compatible orientation, and `states` for the one-circuit slice. A coloured graph whose compatible orientations all have two or more circuits has an empty `states` vector and a nonempty `diagrams` vector. `canonicalize_ograph` identifies cyclic re-rooting and relabeling. It does not identify reversal. For every registered closed example that `certify_closed_ograph` accepts, the stored colours recover the original presentation, and the negated colours recover the reversed presentation.

`coloured_ograph` accepts that presentation, a `BranchColourGraph`, or a `ColouredOGraph`. The shared type is read off its half-edges, mate, over positions, and undirected colours after its own O0–O3 checks pass. Its constructor from a normal o-graph is not called.

Loops, parallel edges, a single vertex, and more than one circuit are in the scope of Proposition 3.2.1. A non-orientable spine, as in Example 3.4.4 of [[BP97#^bp97-sec-0575]], is outside the oriented-spine convention of that chapter. Existence of a branching after a C-move is outside this note.

## Checks

`test/ograph_colour_branching.jl` includes the source directly. It was run with `julia --project=. test/ograph_colour_branching.jl` from `lib/topology.jl`, and the run passed.

The registered closed presentations that `certify_closed_ograph` accepts are `s3_abalone`, `lens_2_1`, `lens_3_1_mst2_example23`, `lens_3_2`, `poincare`, `quaternion`, `s2xs1`, `s2xs1_2vertex`, and `odata_variant`. `s2xs1_2vertex` and `odata_variant` each have a vertex of index −1. The same list refuses `koda_fig13`: its closedness triple is `(true, false, true)`, so the C2 bit fails. There is no separate closed census in the library. The full colour count is checked for the one-vertex example (9 tuples, 2 from a branching, so r = 1) and for two two-vertex one-circuit examples (81 tuples, 2 from a branching). A two-vertex graph whose over strand and under strand are distinct homogeneous circuits has four orientations, index patterns that are not a single `OGraphState`, and a census of 81 tuples of which 1 comes from a branching. A two-vertex graph with two mixed circuits has a census of 81 tuples of which 4 come from a branching.

<!-- XQQQ1000000 -->
