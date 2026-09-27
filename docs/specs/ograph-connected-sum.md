# Connected sum of closed normal o-graphs

Status: normative local-operation, real-Gauss, and integer-weight specification

Scope: closed normal o-graphs with optional integral disk-equation weights

Terminology: [`ograph-terminology.md`](ograph-terminology.md)

## Source statement

Koda defines the connected sum of two normal o-graphs by choosing one directed
edge in each disjoint diagram, requiring the two chosen edges to run around the
boundary of a rectangle in the same direction, and replacing the resulting
two-component normal o-tangle by the four-vertex tangle in Figure 6. Theorem
1.1 states that the result represents the oriented connected sum of the two
represented closed 3-manifolds.

- Work: `[[K10]]`, Yuya Koda, *O-graphic study of closed 3-manifolds*, Hakone
  Seminar 26 (2010), 51--61.
- Construction: `[[K10#^k10-sec-0032]]`, §1.2, PDF p.3 (printed p.53).
- Defining figure and theorem: `[[K10#^k10-sec-0033]]`, Fig. 6 and Theorem 1.1,
  PDF p.4 (printed p.54).
- Local expanded drawing: [MST1 Figure 13](../../../../threads/PJ260825/sources/MST1/connected_sum.pdf).

**Evidence.**

- literature: [K10, §1.2, Fig. 6, Theorem 1.1]
- literature: [MST24, Fig. 13]

The managed Koda Work is
`w-0be57e06-c1d1-47c3-bec1-7208afa6fe59`. Its text capability is currently
`ocr-unchecked`, so the construction, figure, and theorem above were also
checked visually against the
[source PDF](../../../../references/items/article/w-0be57e06-c1d1-47c3-bec1-7208afa6fe59/sources/original/Koda2010.pdf).

## Labelled local operation

The following is MST1 Figure 13 with the four inserted vertices labelled in
the order used by the implementation:

![Labelled connected-sum tangle](assets/ograph-connected-sum.svg)

The labels have the following fixed geometric roles.

| label | role | vertex sign |
|---|---|---:|
| `a` | lower self-loop vertex | `-1` |
| `b` | left splice vertex on the selected edge of `Γ_M` | `+1` |
| `c` | upper self-loop vertex | `-1` |
| `d` | right splice vertex on the selected edge of `Γ_N` | `+1` |

Thus the appended vertex-sign vector, in label order `(a,b,c,d)`, is

```text
(-1, +1, -1, +1).
```

## Real-Gauss splice

Let the input contain `n` real vertices and set

```text
a = n+1,  b = n+2,  c = n+3,  d = n+4.
```

Let `G_M` and `G_N` be the two selected cyclic component words after choosing
their cut seams. With `+v` denoting the over-passage and `-v` the
under-passage, the output component word is

```text
G_M, -b, -a, +a, +d, G_N, -d, -c, +c, +b.
```

The displayed word is cyclic. A cyclic rotation represents the same component
word, but changing any of the eight displayed passage signs does not implement
this local operation under the repository convention.

If an anchor occurrence is supplied for an input component, rotate that cyclic
word so the anchor is first. The cut edge is then the cyclic edge from the
predecessor of the anchor to the anchor. Without an anchor, the stored
last-to-first seam is used.

## Fixed-exterior integer-weight lift

Assume that both input diagrams carry integral weights satisfying the current
disk equation

```text
2 * delta * weight + cp2 = 0.
```

Let `u` and `v` be the weights of the two selected cut edges. Number the ten
new real edges along the two four-vertex bridge strands as follows:

```text
x1: M-out -> -b       x6: N-out -> -d
x2: -b -> -a          x7: -d -> -c
x3: -a -> +a          x8: -c -> +c
x4: +a -> +d          x9: +c -> +b
x5: +d -> N-in       x10: +b -> M-in
```

All input edge weights outside the two cut edges remain fixed. The six open
disk paths through the local tangle have exactly the same doubled boundary
contributions as the two removed edges if and only if

```text
x1 + x10 = u,
x5 + x6  = v.
```

The three internal disk equations are equivalent to

```text
x3 = 1,
x8 = 1,
x2 + x4 + x7 + x9 = -1.
```

Consequently an integral lift exists for every pair of integers `u,v`. The
implementation uses the fixed affine representative

```text
(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
    = (0,0,1,0,0,v,0,1,-1,u).
```

For the stored cyclic output word, the edge sequence begins with `x10`, then
the unchanged edges of `G_M`, followed by `x1,...,x5`, the unchanged edges of
`G_N`, and `x6,...,x9`.

This proves a universal integer lift for the current disk equation with fixed
external weights. It does not by itself identify this affine representative
as a canonical geometric connected sum of homotopy classes of framings; that
requires a separate comparison of the represented framings.

The derivation, implementation boundary, and bounded rerun are recorded in
[the PJ260829 note](../../../../threads/PJ260829/docs/notes/2026-09-27-102357-connected-sum-integer-lift.md).

## Preconditions and postconditions

Both selected components must be closed strands of closed normal o-graphs and
must be distinct. The operation merges them into one component, adds four real
vertices, and must return a closed normal o-graph. In particular, the result
must satisfy C1, C2, and C3; for `N` output vertices the Method-B circuit count
must be `N+1`.

When both operands of `disjoint_sum!` are weighted in the same mode, their
weights are retained. `connected_sum!` transports assigned integral weights by
the lift above and verifies the output disk equation. Z2 transport is not
defined here. `transport_weights=false` explicitly requests an unweighted
output.
