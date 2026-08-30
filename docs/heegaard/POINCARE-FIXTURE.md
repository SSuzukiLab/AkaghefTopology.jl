# Poincare homology-sphere fixture

`test/fixtures/poincare_homology_sphere_v1.hgd.json` records the exact
Heegaard/Gauss code recoverable from Anna Klebanov, *Gauss diagrams of
3-manifolds*, arXiv:math/0308055v1, Fig. 2.15 (PDF p. 20, source asset
`poincare.eps`). The source links are the [arXiv record](https://arxiv.org/abs/math/0308055),
[PDF](https://arxiv.org/pdf/math/0308055v1), and [source archive](https://arxiv.org/e-print/math/0308055v1).

With the source's `M+ = alpha` and `M- = beta` convention, the transcribed
cyclic words are

```text
alpha = [[1,2,3,4,5,6,7], [8,9,10,11,12]]
beta  = [[1,2,3,4,8,6,10], [7,11,12,5,9]]
```

The signs are negative on ids `1:4` and `11:12`, and positive on ids
`5:10`. Consequently the beta relators in the package convention are

```text
g1^-4 g2 g1 g2
g1 g2^-2 g1 g2
```

and their exponent/abelianization matrix is

```text
[-3  2
  2 -1]
```

The matrix has determinant `-1`, so the exact integer path gives
`SNF=diag(1,1)` and `H1=0`. The ribbon computation has 10 boundary cycles
and the fixture uses the singleton partition as a *cellular genus-2
reconstruction candidate*.

That partition is not directly transcribed from Fig. 2.15: the figure does
not expose boundary-cycle ids or a region partition. The fixture therefore
separates `code_transcribed` from `region_decoration_status =
candidate_not_source_transcribed`; it must not be described as a fully
source-transcribed decorated diagram.

The previous Kashaev--Virelizier convention is not used by this fixture and
must not be conflated with the Klebanov transcription.

## Paired-puncture layout certificate

`test/fixtures/poincare_homology_sphere_v1.hgd.layout.json` is a separate
derived artifact. Catlab constructs the ordered-port wiring seed; it does not
place punctures or certify topology. `Topology.Heegaard` independently obtains
the terminal ribbon map

```text
V=4, E=12, F=10, chi=2, genus=0,
face lengths = [2,2,2,2,2,2,3,3,3,3].
```

The ten terminal faces bijectively match the ten original ribbon boundary
cycles after alpha darts are removed. The layout generator then realizes the
outer triangular face `(e1,e12,e7)`, four triangular interior faces, and five
digon bands by fixed corridors. Inside each puncture collar it interpolates
the prescribed port angles to corridor slots without changing their cyclic
order. This is a fixture-scoped certificate, not a generic planar router.

The resulting 12 beta routes contain 346 vertices. Both the package validator
and an independent strict-sign audit check:

- zero pairwise intersections, including exact collinear contacts;
- zero self-intersections;
- zero puncture-interior intrusions;
- exact replay at every `Port.point` endpoint;
- the beta cyclic successor relation and one-time use of all 24 ports;
- core hash `691e24337be7458d8b8309229d22bdb1985c137c628312a26622a5f58625dcc0`;
- `catlab_used=true` and `post_router=certified-fixed-boundary-pl` provenance.

Artifact SHA-256 values:

```text
core fixture    6c72620416df3f001bb9d059240021a3e6ca98c4fd18de6f1cd9349292839719
route data      52eacb19222b2a2790680271197d83fb67b69ed996c0247b24301eaaa490c739
layout fixture  1fe0faf0d35592a1610b287fd86f70b2a6cc7f0cac7621d30277cdc074392286
```

## Independent group check

The GAP replay under `experiments/heegaard/poincare_sphere/gap/` checks the
two relators independently of the Smith-normal-form path. GAP 4.12.1 returns
group order 120, a perfect group, and empty abelian invariants. Together with
the exact integer calculation this verifies

```text
H1 = 0, |pi1| = 120, pi1_ab = 0.
```
