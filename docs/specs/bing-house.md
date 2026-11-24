---
thread: PJ261120
date: 2026-11-22
agent: Grok
---

# Bing's house as a coloured o-graph

<!-- XQQQ1500000 -->

Scope: one coloured o-graph for Casler's Bing house, the standard polyhedron it determines, and the branching count of that polyhedron. The convention id stays `bp97-fig2.5-rl-tau/draft1`. The face-cycle status stays `prop321-calibration-b/draft`.

Section III opens in [[Cas65#^cas65-sec-0014]]. Figure 2, the collapse of the cube onto X, and the sentence that X is a standard spine of the cube are [[Cas65#^cas65-sec-0015]]. The sentence that X is a special spine of the cube is [[Mat73#^mat73-sec-0012]]. The local model that places the upper half-plane on one singular line and the lower half-plane on the other is [[BP97#^bp97-sec-0547]]. The o-graph, the gluing, and Figure 2.5 are [[BP97#^bp97-sec-0549]]. The passage from an embedded spine to an o-graph begins after Figure 2.7 in [[BP97#^bp97-sec-0550]] and concludes in the chunk [[BP97#^bp97-sec-0551]], before Proposition 2.1.6. Corollary 3.1.7 is [[BP97#^bp97-sec-0565]]. Proposition 3.2.1 is [[BP97#^bp97-sec-0566]]. Remark 3.2.2, on circuits of a normal o-graph, begins in [[BP97#^bp97-sec-0567]] and continues in [[BP97#^bp97-sec-0568]]. Section 2.2 and Figure 2.8 are [[BP97#^bp97-sec-0552]]. Proposition 2.2.1 is [[BP97#^bp97-sec-0553]]. The anchor [[BP97#^bp97-sec-0551]] also closes Proposition 2.1.6 and Remark 2.1.7.

Moment 1 is `threads/PJ261120/sources/pro/core-20261122-152820.md`. Moment 3 is `threads/PJ261120/sources/pro/core-20261122-165414.md`. Both are advice. The local model prints an upper half-plane and a lower half-plane. Calling the upper sheet the wall is inferred. Rows that name the lead rule are computed counts under the draft calibration. A passing count does not select the over pair. The book does not print the boundary words, the circuit lengths, or the statement that this graph is Bing's house.

## Polyhedron

The cube model is right-handed, with the partition horizontal. The upper membrane meets the left wall on the partition. The lower membrane meets the right wall on the partition. Those two points are the true vertices. The front and back arcs of the equator, the upper loop through the ceiling and the upper tube, and the lower loop through the floor and the lower tube, are the four triple lines. Casler states that the complement of the singular 1-skeleton is three open disks ([[Cas65#^cas65-sec-0016]]). The sheet trace closes as one 10-gon and two monogons, of lengths 10, 1, and 1.

The 10-gon, in the edge orientation front equator from the upper vertex to the lower vertex, back equator the same way, upper loop from the upward ray to the inward ray, and lower loop from the downward ray to the inward ray, has the boundary

(back −, front +, lower loop +, front −, upper loop −, front +, back −, upper loop +, back +, lower loop −).

The monogons are the upper membrane and the lower membrane.

## Coloured o-graph

`bing_house_coloured_ograph` in `examples/BingHouse.jl` stores slots counter-clockwise. At the upper vertex the order is up, back, in, front, and the over slots are back and front. At the lower vertex the order is front, in, back, down, and the over slots are in and down. Ordered by smaller half-edge, the colours are 0, +1, −1, 0, on the upper loop, the back equator, the front equator, and the lower loop. The constructor stores the records (4, −1), (2, +1), (1, 0), (6, 0), which is not the smaller-half-edge order.

On this graph, four orientations satisfy the colour rule of Proposition 3.2.1. The same trace on the other screw-consistent choices gives the same three boundary words, and the repository gauge splits those choices into more than one class. That enumeration is outside the unit test. The registered class is the one for which `prop321_orientations` is non-empty. In the library lettering the lower loop is positive from the inward ray to the downward ray, so a positive geometric step on that loop is the letter −4. With that sign, the geometric word is a rotation of the word locked by the test. The test identifies words up to rotation, reversal, and the order of the regions. The membrane names are the edge names above.

Each of those four orientations runs in three circuits, of lengths 1, 1, and 2. `OGraphState` stores one cyclic word, so this colouring has no one-component normal presentation in that type. Among the fixtures in `vl_virtual_gauss_inputs` for which `normal_coloured_ograph` returns, none has the same `canonical_coloured_ograph` form. A conversion that throws is skipped and is outside that claim.

## Checks

| ID | Locator | What is checked | Function | Test | Review |
|---|---|---|---|---|---|
| O0–O3 | [[BP97#^bp97-sec-0549]] | the registered graph passes the four checks of Definition 2.1.4 | `validate_coloured_ograph` | `Bing's house coloured o-graph` | moment 3, ZA21, in `threads/PJ261120/sources/pro/core-20261122-165414.md` |
| discs | [[Cas65#^cas65-sec-0016]] states three open disks; the words are the trace under the draft calibration | three cycles, lengths 1, 1, 10, up to rotation, reversal, and region order | `face_cycles` | same test, `_BING_WORDS` | moment 3, ZA22. The book does not print the words |
| Euler | the lead rule, not a printed sentence | n=2, E=4, 8 half-edges, 12 germs, d=3, χ=1, length sum 12 | `face_cycles` | same test | moment 3, ZA23. Computed under the draft calibration |
| boundary | [[BP97#^bp97-sec-0552]], [[BP97#^bp97-sec-0553]] for the constructions; [[Cas65#^cas65-sec-0015]] for the cube | `thickening_boundary` has one closed orientable component of Euler characteristic 2; the provisional Section 2.2 adapter has Euler list [2] | `thickening_boundary`; `OGraphColourBoundary.boundary_from_coloured` | same test. Closed and orientable are asserted for the first function. The second asserts the Euler list | moment 3, ZA24. The adapter's reading status stays provisional |
| dual | the lead rule, not a printed sentence | 4 orientation-reversing face pairs and 3 dual edges | `dual_face_pairings`, `dual_edge_orbits` | same test | moment 3, ZA25. Computed under the draft calibration |
| homology | the lead rule; contractibility of the cube is [[Cas65#^cas65-sec-0015]] | free rank 0 and no torsion | `cellular_homology` | same test | moment 3, ZA26. The Smith form is the draft face cycles |
| branch | [[BP97#^bp97-sec-0565]] | all 8 region orientations induce mixed directions on every triple line, counting each letter as a germ | `_branching_failures` | same test | moment 1, ZA4 and ZA10; moment 3, ZA19 and ZA27. The corollary is the orientation test. Definition 3.1.2 is not that test |
| P2 | [[BP97#^bp97-sec-0566]] | 4 orientations satisfy the colour rule | `prop321_orientations` | same test | moment 3, ZA28. The count is computed, not printed |
| circuits | [[BP97#^bp97-sec-0567]], [[BP97#^bp97-sec-0568]] | each of those orientations has circuit lengths 1, 1, 2 | `_strand_lengths` | same test | moment 3, ZA28. Split from the colour rule |
| gauge | the lead rule, not a printed sentence | slot rotation of the upper vertex preserves the canonical form | `canonical_coloured_ograph` | same test | moment 3, ZA29. Vertex exchange and the 21 classes are outside this test |
| registry | the fixture list, not a printed sentence | among fixtures for which conversion returns, none has this canonical form | `vl_virtual_gauss_inputs`, `normal_coloured_ograph`, `canonical_coloured_ograph` | same test. A thrown conversion is skipped | moment 3, ZA30 |

<!-- XQQQ1000000 -->
