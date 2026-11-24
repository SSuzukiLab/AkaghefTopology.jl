---
thread: PJ261120
date: 2026-11-22
agent: Grok
---

# Invariants of a coloured o-graph

<!-- XQQQ1500000 -->

The discs of a coloured o-graph are the cycles of `face_cycles`. The default calibration is B, and the status is `prop321-calibration-b/draft`, the value of `OGRAPH_FACE_CYCLE_STATUS`. The reading is the one in `docs/conventions/OGraphColours.md`, drawn from [[BP97#^bp97-sec-0566]] Proposition 3.2.1 together with §2.1 and §3.2. The dependencies are the draft bullets that follow the sentence 「この導出が依存する draft は次である。」 in that convention page. This spec uses that list. The ledger record is the durable claim of 2026-11-22 in `threads/PJ261120/LEDGER.md`.

`src/ograph_colour_invariants.jl` and `test/ograph_colour_invariants.jl` are the two files of this lane. `src/Topology.jl` does not include the source. A caller includes the source after `using Topology`.

## One function for each object

| Object | Function |
|---|---|
| Disc partition, cell counts, χ(P) | `face_cycles` |
| H₁ | `cellular_homology` |
| π₁ presentation | `coloured_pi1`, whose letters are `FaceCycle.edge_letters` |
| Boundary of the thickening | `thickening_boundary` |
| Dual face pairs | `dual_face_pairings` |
| Dual edges | `dual_edge_orbits` |
| Closed colouring | `normal_coloured_ograph` |
| Ordered tetrahedra | `disk_tetrahedra` on the `VirtualLink` |

`coloured_pi1` orders the singular edges by the smaller half-edge id and contracts the Kruskal tree in that order, the same order `cellular_homology` contracts. The test checks that the abelianization equals `cellular_homology`. The ordered-tetrahedron certificate is `disk_tetrahedra`. This lane has no wrapper around it.

The identifiers G1, G2, G3, T1, T2, T3, and CHI named checks of a chart and a dual walk that this lane used to build. The checks in the table above are the current ones.

The test file keeps one independent walk: chart A, with leg 1 on the counter-clockwise side of an over slot, and transport ρ∘τ^{−α}. The convention page records that this gluing is calibration B with transport ρ∘τ^α, because ρτρ = τ^{−1}. The test compares the sorted corner partitions. The walk is not a library function.

## Closed examples

The normal pipeline is `ograph_h1(ograph_pi1_presentation(certify_closed_ograph(...)))`, algorithm `exponent-sum-smith-normal-form`. Its presentation is the Koda/MATLAB edge presentation. The coloured column is `cellular_homology` of `normal_coloured_ograph`. The two columns agree on every row. `koda_fig13` is absent: `normal_coloured_ograph` throws, and the closedness witness reports criteria `(true, false, true)`.

| Registry name | n | d | χ(P) | Boundary | H₁ | Source string |
|---|---|---|---|---|---|---|
| s3_abalone | 1 | 2 | 1 | one component, χ = 2 | 0 | |
| lens_2_1 | 2 | 3 | 1 | one component, χ = 2 | Z/2 | |
| lens_3_1_mst2_example23 | 2 | 3 | 1 | one component, χ = 2 | Z/3 | MST2 Example 23 |
| lens_3_2 | 3 | 4 | 1 | one component, χ = 2 | Z/3 | |
| quaternion | 3 | 4 | 1 | one component, χ = 2 | Z/2 × Z/2 | K10 Fig. 14; positive trefoil o-graph |
| odata_variant | 4 | 5 | 1 | one component, χ = 2 | 0 | |
| s2xs1 | 4 | 5 | 1 | one component, χ = 2 | Z | |
| s2xs1_2vertex | 2 | 3 | 1 | one component, χ = 2 | Z | PJ260830: 261009_koda-bp-ograph-conventions |
| poincare | 5 | 6 | 1 | one component, χ = 2 | 0 | |

`odata_variant` has the pipeline's trivial group. The registry entry has no source string, and this row gives it no manifold name. The quaternion row keeps both the registry name and the stored source string.

On each of these nine graphs the counted objects are: `2n` edges, `4n` half-edges, `6n` germs, disc-boundary lengths summing to `6n`, `d = n + 1`, `χ(P) = 1`, one orientable closed boundary component of Euler characteristic 2, `χ(∂M) = 2`, `2n` dual face pairs, each pair orientation-reversing and using slots `1:4`, and `d` dual edges. The calibration field is `:b`.

[[BP97#^bp97-sec-0549]] Definition 2.1.4 is the coloured o-graph. The same section pastes discs after the Fig. 2.5 gluings. [[BP97#^bp97-sec-0551]] Proposition 2.1.6 is the bijection with oriented standard spines up to the C-move. For a closed manifold whose spine is P, the thickening retracts onto P and the boundary is one 2-sphere, so π₁(M) = π₁(P) and H₁(M) = H₁(P). [[BP97#^bp97-sec-0703]] §10.1 states cellular homology for a branched spine, with `n` vertices, `2n` edges, and `n + 1` discs, and for the closed case one 3-cell whose boundary is zero. Section 10.1 gives neither a maximal-tree presentation nor a Smith form. The positive letter is the sign `face_cycles` stores, the letter of the smaller half-edge. The Smith form drops invariant factors equal to 1. That is the cellular-chain convention in Hatcher's *Algebraic Topology*. Hatcher has no cite key in this repository. The same dropping convention is the one in `ograph_h1`.

## Colourings outside Proposition 3.2.1

The one-vertex graph with mate `(2, 1, 4, 3)` and over pair `(2, 4)` has two loops. Each loop joins an over slot to an under slot. The abalone colouring of this mate stores −1 and +1, so colours `(0, 0)` and `(1, 1)` are outside Proposition 3.2.1. The general counts hold on both. The closed-input requirements are a separate check and are not the acceptance condition for these two graphs.

| Colours | n, E, d | χ(P) | H₁ | Boundary |
|---|---|---|---|---|
| (0, 0) | 1, 2, 2 | 1 | Z/4 | one component, χ = 2, orientable, closed; 4 vertices, 6 edges, 4 faces |
| (1, 1) | 1, 2, 2 | 1 | Z/4 | one component, χ = 2, orientable, closed; 4 vertices, 6 edges, 4 faces |

The sorted disc partitions differ. The manifold is unidentified. A group Z/4 and one boundary sphere are the homology type of some lens space L(4, q), and q is undetermined.

The negative control starts from the Proposition 3.2.1 colouring of `s3_abalone` and sets the edge whose stored colour is −1 to +1. The resulting coloured graph is the `(1, 1)` row. Its H₁ is Z/4, with `d = 2` and `χ(P) = 1`. No manifold is expected.

## Gauge

On the nine closed graphs and on both loops, a cyclic shift of each vertex tuple by 1, 2, or 3 slots, carried with the over slots, leaves H₁, the disc count, χ(P), and χ(∂M) unchanged. Exchanging vertices 1 and 2 does the same on every graph with at least two vertices. `canonical_coloured_ograph` does the same on all of them. Every image passes `move_sanity`.

## Candidate A

`face_cycles(; calibration=:a)` on `s2xs1_2vertex` gives `n = 2`, five discs, `χ(P) = 3`, and 12 germs. The thickening has three components, each of Euler characteristic 2, and the boundary Euler sum is 6. The closed-input requirements `d = n + 1`, `χ(P) = 1`, and exactly one boundary sphere fail on this graph. Among the registered closed examples with at most two vertices, it is the smallest example of that failure. The adopted calibration remains B.

## C-move and MP-move

On the nine closed graphs, 54 C images and 2 MP images pass `move_sanity`, and H₁ is unchanged on each of them. The 54 C images are the four frames of `s3_abalone` and frames 3 and 4 at every vertex of the other eight graphs. The two MP images are the forward move of `lens_2_1` from vertex 2 to vertex 1 at frames `(2, 2)`, and the forward move of `s2xs1_2vertex` from vertex 2 to vertex 1 at frames `(4, 2)`. Each has 3 vertices, 4 discs, `χ(P) = 1`, and one boundary component of Euler characteristic 2. The first has H₁ = Z/2. The second has H₁ = Z. `s3_abalone`, `lens_3_1_mst2_example23`, `poincare`, and `quaternion` have no MP placement.

The other 64 images make `move_sanity` throw. Every such C image is frame 1 or frame 2. The smallest example is `lens_2_1`, vertex 1, frame 1 and frame 2: the image has `n = 2`, `d = 3`, and `χ(P) = 1`, and the throw names χ(∂M) = 2χ(P) (boundary Euler characteristic 0), a non-orientable boundary component, and a dual face pair that does not reverse orientation. The same two checks, orientability and orientation reversal, fail on every one of the 64. The move stays in `src/ograph_colour_moves.jl`. An image that fails the sanity rule is omitted from the homology table.

## Dual ideal triangulation

[[BP97#^bp97-sec-0557]] introduces ideal triangulations as tetrahedra with vertices removed, faces glued in pairs by orientation-reversing simplicial maps, and says the process is not formalized. `dual_face_pairings` and `dual_edge_orbits` are the pairs and the edges used here. `disk_tetrahedra` on `s3_abalone` returns one tetrahedron, with `face_orientations_verified` and `edge_orientations_verified` both true. That call reads the normal certificate. An isomorphism between those tetrahedra and the slot gluings is not claimed.

## Advice

Moment 1 is `threads/PJ261120/sources/pro/invariants-20261122-011544.md`. Moment 3 is `threads/PJ261120/sources/pro/invariants-20261122-022517.md`. Both are advice. Moment 3's IA9 says that Fig. 2.5 prints τ^α on one drawing and prints neither ρ nor a refusal of ρ∘τ^{−α}. Calibration B keeps CORE's transport ρ∘τ^α. The test-only walk is the other description of the same gluing.

<!-- XQQQ1000000 -->
