---
thread: PJ261120
date: 2026-11-22
agent: Grok
---

# Boundary of a coloured o-graph

<!-- XQQQ1500000 -->

Scope: from a coloured o-graph, build the trivalent graph Θ₀ of [[BP97#^bp97-sec-0554]] Figure 2.12, carry the Z₂ colours printed under that figure, and cap the ribbon of [[BP97#^bp97-sec-0552]] Figure 2.8. The result is one Euler characteristic and one component count per connected component of the capped surface. The polyhedron is not built. The C2 graph of Figure 1.2 is a different trivalent graph and is not used.

`reading` is `:provisional`. The strand permutations, the Z₂ digits, and the four trivalent vertices are read from the figures. Which arm of the Figure 2.12 cross is the incoming under-passage is a placement choice, recorded as DIR below. The figure reading and the condition check are advice. The review at the end records what that advice supports in the PDF.

Printed page n of [[BP97]] sits near PDF page n+9. Figure 2.8 is on PDF page 26, Figures 2.9 and 2.10 on PDF page 27, and Figures 2.11–2.13 on PDF page 28.

## Condition identifiers

| ID | Statement | Locator | How read | In code | In tests | Review |
|---|---|---|---|---|---|---|
| B1 | Each o-graph vertex contributes four trivalent vertices of Θ₀, with north (N, E, W), south (E, W, S), west (N, S, W), and east (N, S, E). There is no internal edge. | [[BP97#^bp97-sec-0554]] Fig. 2.12, and the same cross in [[BP97#^bp97-sec-0553]] Fig. 2.9 | figure | `THETA_VERTICES`, `OWNER`, `ROTATION` | counts `4n` vertices and `6n` edges | item check: Fig. 2.12 has four trivalent vertices and no internal edge. `OWNER` is unused; the partition is `ROTATION` |
| B2 | Six edge panels. Positions 0, 1, 2 are left to right. A tuple `(a,b,c)` sends the left, middle, and right strands to those ranks. The Z₂ triple is indexed by the arrival rank. | [[BP97#^bp97-sec-0554]] Fig. 2.12 | figure; the half-turn is reconstructed | `SAME_ROW`, `MIXED_ROW`, `_strand_table` | abalone colours; half-turn of a nonzero same-type edge | item check: the six drawn panels match; the index convention and the half-turn are repository encoding |
| B3 | Colour 0 is an untwisted band. Colour 1 is the half-twisted band of Fig. 2.8. The walk flips side on a colour-1 strand. | [[BP97#^bp97-sec-0552]] Fig. 2.8 | figure; the side flip is inferred | `twist` inside `_boundary_directed` | the side map is a permutation on every registered input | item check: the band is Fig. 2.8; the side flip implements that band and the test does not name an expected flip |
| B4 | At each trivalent vertex the three stubs are ordered counterclockwise on the page. | [[BP97#^bp97-sec-0554]] Fig. 2.12, local stars | figure | `ROTATION`, then `ccw` and `cw` | used by every surface in the test | item check: the cyclic order is the local star; the stub names are a transcription |
| B5 | The ribbon boundary is the side-walk. Each unoriented boundary circle is capped by one disc. The two directions of a circle are one circle. | [[BP97#^bp97-sec-0552]] paragraph before Fig. 2.8 | text, and the reversal orbit is inferred | `image`, then one orbit under reversal in `cycle_theta` | one-vertex graph has 4 circles; directed states are paired | applied. The item check and the orbit-count check both give the reversal identity |
| B6 | For one component, `euler = v - e + b`, with `b` the number of ribbon circles of that component. The number of boundary components of the 3-manifold is the number of these capped components. | [[BP97#^bp97-sec-0553]] Prop. 2.2.1 | inferred from the ribbon | `BoundaryComponent` inside `_boundary_directed` | the Euler relation; closed examples have one component of characteristic 2 | applied as an inference, not as a printed formula |
| DIR | PD slot 1, the incoming under-passage, is placed on the south arm. Counterclockwise slots are then south, east, north, west. The over pair is then the horizontal arms and the under pair is the vertical arms, which is the shared vertex model in `docs/conventions/OGraphColours.md`. | repository port table in `docs/specs/ograph-terminology.md` §6; the arm choice is not printed in Fig. 2.12 | inferred | `SLOT_DIRECTION`, `_ports` | gauss input and `ColouredOGraph` agree on the seven examples | stays inferred. Fig. 2.5 leg 1 is unregistered in that convention page, so a Fig. 2.5 leg is not copied onto an arm |

Proposition 2.2.1 names the graph of Figure 2.9 and the oriented surface obtained when every Z₂ colour is 0. The proof first obtains Θ₀ from Figures 2.10 and 2.11 and the rules of Figure 2.12, then applies Figure 2.13 and concludes that the result is that colour-0 graph. The code keeps Θ₀ and the printed Z₂ colours. The two descriptions are the same surface when the proof applies. The code does not insert the curls of Figure 2.13.

## Panels

Same-type ends use the top row of Figure 2.12. Let σ be +1 when both ends are above and −1 when both ends are below. The row is the one for σ times the Z₃ colour. Both 3-cycles are unchanged by a half-turn of the picture, so the same-type row is not reflected when the edge is reversed.

| Effective index | Permutation | Z₂ at arrivals 0, 1, 2 | Panel |
|---|---|---|---|
| 0 | (0, 1, 2) | 0, 0, 0 | top left |
| +1 | (2, 0, 1) | 1, 1, 0 | top middle |
| −1 | (1, 2, 0) | 0, 1, 1 | top right |

Mixed ends use the bottom row, which is drawn with the above end at the top of the panel and the below end at the bottom. The opposite mixed pattern is not drawn. Reversing the edge, keeping the Z₃ colour, replaces a permutation p by q(j) = 2 − p⁻¹(2 − j) and carries each Z₂ mark with its strand. That replacement is `_rotate180`.

| Colour | Drawn permutation, above end first | Z₂ | After the half-turn |
|---|---|---|---|
| 0 | (2, 1, 0) | 1, 1, 1 | the same permutation and the same marks |
| +1 | (0, 2, 1) | 1, 0, 0 | (1, 0, 2) with marks 0, 0, 1 |
| −1 | (1, 0, 2) | 0, 0, 1 | (0, 2, 1) with marks 1, 0, 0 |

The Z₂ digits were read from the glyphs under the braid ends on PDF page 28. The input colours of the six panels are 0 and 0, 1 and −1, −1 and 1 on the top row, and 0, 1, −1 on the bottom row. The minus glyphs on the three nonzero mixed and opposite-same panels were read from the page layout. They were not cut out as separate components.

Above means an over-passage. For a Gauss word, the tail is the entry g[k] and the head is g[k+1], cyclically. The tail port is outgoing and the head port is incoming. `prop321_colour` is the rule of Proposition 3.2.1: equal passage marks give 0, below then above gives +1, and above then below gives −1.

An edge-end uses the three stubs of its arm. After the arm has been turned so that it points down the page, left-to-right order is `DOWN_ORDER`: south (S0, S1, S2), north (N2, N1, N0), east (E2, E1, E0), west (W0, W1, W2). The head sits at the bottom of the panel, so its left-to-right order is the reverse of that list.

## Rotation

Page coordinates have y downward. Counterclockwise at the four vertices, named by stub:

| Vertex | Order |
|---|---|
| north | N1, W0, E0 |
| south | E2, W2, S1 |
| west | N0, W1, S0 |
| east | N2, S2, E1 |

Side 0 takes the counterclockwise successor of the arrival stub. Side 1 takes the clockwise successor. Crossing a colour-1 strand flips the side. The successor map on the 24n side-states is checked to be a permutation before the circles are counted.

## What is counted

For n o-graph vertices the intact graph has 4n vertices of Θ₀ and 6n edges. The ribbon retracts to that graph, so its characteristic is −2n. Capping one disc on each boundary circle gives characteristic −2n plus the number of circles. One sphere is one component with characteristic 2, hence 2n+2 circles. The circle count is not the number of boundary components of the 3-manifold.

The first walk listed every directed side-state cycle. Reversal R(h, s) = (partner(h), 1 ⊕ s ⊕ twist(h)) satisfies R² = 1 and (R T)² = 1, because the partner map is an involution, the Z₂ mark agrees on the two ends of a strand, and the two turns at a vertex are inverses. Each circle was stored twice. On the one-vertex o-graph the graph is V = 4, E = 6, and that walk found 8 directed cycles, characteristic 6. The geometric count is 4 circles, of degrees 1, 1, 5, 5, and the characteristic is 2. The same doubling turns the sphere count 2n+2 into the reported characteristics 2n+4. The panels and the arm placement do not enter that identity. `ribbon_cycles` is now the number of reversal orbits. A stored degree is the length of one directed traversal. The degrees of one component sum to twice its edge count. A degree of 1 is a monogon face of that cellulation.

`boundary_from_gauss` reads words and vertex signs. `alter_edge = 1` adds 1 in Z₃ to the first edge of the first word, with −1, 0, 1 in that cyclic order. `boundary_from_coloured` reads `vertex_halfedges`, `mate`, `over_positions`, and `edge_colours`. It accepts the `ColouredOGraph` constructor. The smaller half-edge is the tail. The stored colour is not negated. Over marks on slots (1, 3) and on slots (2, 4) are both accepted. On this adapter, `alter_edge = 1` shifts the edge whose smaller half-edge is first in increasing order. That edge need not be the first Gauss edge.

## What [[BP97]] predicts

The following is read from the text. After the circle count was corrected, the closed examples below are an acceptance check of that count. The arm placement was not chosen to satisfy it.

[[BP97#^bp97-sec-0624]] §5.2, PDF page 79: for a normal o-graph, N1 and N2 place it as the normal o-graph of a branched standard spine; C1 gives Euler characteristic 2 for the boundary; C2 gives a connected boundary, hence a 2-sphere; C3 gives the trivial bicoloration. C1 and C2 are the circuit axioms on [[BP97#^bp97-thm-0535]] PDF page 15, including the connectivity of the Figure 1.2 graph. They apply to the Proposition 3.2.1 colouring of a closed normal o-graph. They do not apply to a freely chosen Z₃ colouring.

Adding 1 in Z₃ to one edge leaves the Proposition 3.2.1 colouring. §5.2 does not say that the shifted colouring bounds a sphere. §2.2 still applies to it. Theorem 2.1.5 says that an o-graph with at least two vertices represents some compact connected oriented 3-manifold with nonempty boundary. [[BP97#^bp97-sec-0551]] Remark 2.1.7 says that the constructions of §2.1 remain valid for one vertex.

## Checks

`test/ograph_colour_boundary.jl` includes the source directly. It checks the Proposition 3.2.1 values, the Euler relation, the counts 4n and 6n, agreement of the Gauss path with a `ColouredOGraph` built from the same words, invariance of a both-above colour-+1 edge under reversal of its tail, and that an all-zero colouring still produces a surface. For the one-vertex graph it requires 4 circles and characteristic 2. For each closed name it requires one component of characteristic 2 under the Proposition 3.2.1 colouring. A colour shifted by +1 is printed and is not required to be a sphere. The panel tables and the DIR placement stay as read.

The seven names in the lane report are `s3_abalone`, `lens_2_1`, `lens_3_1_mst2_example23`, `poincare`, `quaternion`, `s2xs1`, and `s2xs1_2vertex`. The same sphere check covers `lens_3_2` and `odata_variant`. `koda_fig13` stays outside it.

## Review

The figure reading is `threads/PJ261120/sources/pro/boundary-20261122-003051-answer.md`. Advice. The six permutations and the six Z₂ triples agree with the tables above. The two strokes on a top-row panel are two single-edge inputs of that panel, which is the effective-index rule already in the table. The opposite mixed pattern is the half-turn implemented by `_rotate180`. The proof builds Θ₀ from Figure 2.12 and then applies Figure 2.13. That answer also reads the final Figure 2.9 graph as two monogon edges and two through-edges at each o-graph vertex, with the same counts V = 4n and E = 6n. Its ZA58 says the Θ₀ ribbon, with the colour-1 half-twist, has the same underlying surface. This code builds Θ₀, keeps the printed Z₂ colours, and leaves the curls of Figure 2.13 uninserted.

The condition check is `threads/PJ261120/sources/pro/boundary-20261122-013442.md`. Advice. It reads Figure 2.12 as four trivalent vertices and no internal edge, and it locates the excess characteristic in the directed double count of boundary circles. That account is the reversal identity above. The one-vertex graph is the hand check. DIR stays inferred. `docs/conventions/OGraphColours.md` still has `OGRAPH_FACE_CYCLE_STATUS` equal to `unresolved-fig25-vertex-chart/draft`, so this lane uses the shared placement of the over pair on the horizontal arms and leaves a Figure 2.5 leg off the arms of Figure 2.12.

The corrected count is checked in `threads/PJ261120/sources/pro/boundary-20261122-042259.md`. Advice. It reads §2.2 as V = 4n, E = 6n, and F equal to the number of unoriented ribbon circles, with F = 2n+2 when the boundary is one sphere. It recomputed the one-vertex side-walk as four circles of degrees 1, 1, 5, 5 and characteristic 2. A face of degree 1 is one monogon circle. The sphere check on a closed normal o-graph with the Proposition 3.2.1 colouring is the acceptance check of [[BP97#^bp97-sec-0624]] §5.2. That pass leaves DIR inferred. The Julia run gives one component of characteristic 2 on each closed name.

<!-- XQQQ1000000 -->
