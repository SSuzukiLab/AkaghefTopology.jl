# Knot PD coordinates and polynomial data

## Definitions

### PD coordinates

For the crossing visits of [signed Gauss data](knot-gauss.md), a PD row lists
the four incident arc labels in the local cyclic slot order. Slots 1 and 3
form the under strand; 2 and 4 form the over strand. The local slot vectors
are `(1,0),(0,1),(-1,0),(0,-1)`. For oriented strand vectors from incoming to
outgoing slot, the crossing sign is the sign of `det(over,under)`.
Equal arc labels join two slots. `KnotPD` stores rows, signs and component
count. The empty PD has exactly one crossing-free component.

Raw PD input is an iterable of four-entry rows. Despite the name
`pd_from_matrix`, an ordinary Julia matrix must be passed as `eachrow(M)`;
iterating the matrix itself supplies scalar entries. A `KnotPD` record is
passed through its `.rows` when a raw-row conversion is wanted. These
conversions reconstruct signs from the local frame rather than trusting
separate supplied signs.

For a Gauss traversal, arc labels are the positions before and after each
visit, cyclically in `1:2n`. If under and over visits have pairs
`(u_in,u_out)` and `(o_in,o_out)`, the PD row is
`(u_in,o_out,u_out,o_in)` for positive crossing sign and
`(u_in,o_in,u_out,o_out)` for negative sign.

### Exact polynomial records

`ExactKnotLaurentPolynomial` records a variable symbol and a finite integer
exponent to `BigInt` coefficient map, with zero coefficients removed by its
normalizing constructor. Equality compares that variable and map. It is a
distinct type from the rational-coefficient
[LaurentPolynomial](laurent-polynomials.md).

### Alexander coordinates

The Wirtinger arc indices change at under visits. For a crossing of sign `s`,
the row of the repository Alexander matrix adds `1-t^s` at the over arc,
`t^s` at the incoming under arc and `-1` at the outgoing under arc; coincident
arcs add contributions. The polynomial is the determinant with the last row
and column removed, shifted to minimum exponent zero and with positive
highest-degree coefficient. The crossing-free case is one. A zero determinant
is rejected. This fixes a representative of the usual `±t^k` ambiguity.

### Jones coordinates

For a PD row, A smoothing pairs slots `(1,2),(3,4)` and B smoothing pairs
`(2,3),(4,1)`. For a state with `b` B choices and `l` circles, its bracket
contribution is `A^(n-2b)*(-A^2-A^-2)^(l-1)`. The writhe correction is
`(-1)^w*A^(-3w)`, where `w` is the sum of crossing signs. Substitution
`t=A^-4` defines the output: all surviving A exponents must be divisible by
four. The crossing-free case is one.

These formulas define the repository's polynomial calculations on its
coordinates. They do not turn a nonplanar signed word into a classical knot.

### Surgery derivative

The derivative input to the [Casson–Walker surgery API](casson-walker.md)
uses CW5. The [Alexander normalization dictionary](../conventions/CassonWalker.md#alexander-normalization-dictionary)
owns the translation from the stored polynomial, its five-fixture check,
and the remaining source-identification boundary. There is no automatic
library operation for that translation.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| P1 | PD rows have four positive arc labels, each label occurs twice, signs are `±1`, and the chosen oriented traversal visits each crossing once on each strand. | PD validators check the combinatorial data. |
| P2 | Input represents one component; the empty input means one unknot component. | Knot PD conversion/checking rejects links. |
| P3 | The signed data are classically realizable (K2). | Polynomial evaluation does not supply an independent classicality proof. |
| P4 | The Jones state sum has fewer than 62 crossings and its normalized exponents are divisible by four. | Explicit guards; the upper bound is not a practical run-time promise. |
| P5 | Exponents of an exact polynomial and the negations requested by variable inversion fit `Int`; normalized records omit zero coefficients. | Raw typed field construction and integer negation do not detect all violations. |

## Coordinate choices and evidence

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.
Definitions follow [knot_classical_invariants.jl](../../src/knot_classical_invariants.jl).
[knot_classical_invariants tests](../../test/knot_classical_invariants.jl)
compare the five PD fixtures `0_1,3_1,4_1,5_1,5_2` with independent listed
polynomials; they also check selected R1/R2/R3, relabelling, cyclic shifts,
traversal reversal and mirror. In this convention mirror changes Jones by
`t ↦ t^-1`; the tested trefoil changes while its Alexander polynomial does not.
That finite scope is `verified-computationally`, not a proof for all diagrams.

question: unverified — Registered cite keys and exact source locators for
the PD framing, Alexander minor and bracket normalization have not been
resolved in this lane. The formulas and fixture values are established local
definitions; their intended classical invariant identification is supported
by those fixtures. A general comparison with the source crossing and
variable conventions is still required and could reveal mirror/variable
inversion. No placeholder literature key is introduced.

Registered surgery source:
[Les07](../../../../references/items/article/w-bb972e2b-4a1b-4ec9-b12d-e5e5d6add53a/Les07.md).

## Operations

`KnotPD`, `pd_matrix`, `pd_from_gauss`, `pd_from_matrix`, `pd_to_gauss`,
`verify_pd`, `ExactKnotLaurentPolynomial`, `polynomial_terms`,
`knot_polynomial_string`, `alexander_polynomial`, `jones_polynomial`,
`invert_knot_variable`.
