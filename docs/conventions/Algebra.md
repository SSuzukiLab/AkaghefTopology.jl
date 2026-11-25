# Algebra conventions

Scope: coordinate tensors, finite algebra models and the legacy tensor
expression bridge. Laurent coordinates are defined in
[their specification](../specs/laurent-polynomials.md); the geometric use of
these tensors is in [Hopf invariants](../specs/hopf-invariants.md).
Heegaard and o-graph orientations are owned by their module pages. The
Kuperberg/cellular comparison below records an unresolved crossing-sign
contradiction at the tensor interface.

The choices in Tensor ports through Convention dictionary have status **fixed 2026-11-25 Codex lane lib-api-rest
(rule a)**: it is the existing code/test choice. A fixed coordinate choice
does not upgrade an unverified mathematical identification.

## Tensor ports

In an ordered basis `e_i`, multiplication is `Mu[i,j,k]` for the coefficient
of `e_k` in `e_i e_j`; coproduct is `De[i,j,k]` for the coefficient of
`e_j ⊗ e_k` in `Δ(e_i)`. Antipode is `An[output,input]`.
The one-leg names are `Et` (unit), `Ep` (counit), `Ir`, `Il`, `Cr`, `Cl`.
`Tp[output,right,input,delta_input]` contracts the first coproduct output
with the second multiplication input; `Tm` inserts the antipode on that
contracted leg. `Un` is the identity matrix. These names and port choices
are fixed by [tensor_network.jl](../../src/tensor_network.jl),
[hopf_invariant.jl](../../src/hopf_invariant.jl) and the tests
[topology_core.jl](../../test/topology_core.jl),
[tensor_validation.jl](../../test/tensor_validation.jl).

Coordinates pair bilinearly, without complex conjugation. Source comparison:
[[Kup91#^kup91-blk-0023]] evaluation on a vector and
[[Kup91#^kup91-blk-0039]] multiplication tensor.

## Integral names

The stored `*_integral` fields are functional coordinates and the stored
`*_cointegral` fields are element coordinates for the group-algebra and
Kuperberg routes. For a group basis the former is the identity delta and the
latter is the sum of basis elements. Neither field name alone certifies its
axioms for an arbitrary record. The tests are
[group_algebra.jl](../../test/group_algebra.jl),
[hopf_properties.jl](../../test/hopf_properties.jl) and
[kuperberg_invariant.jl](../../test/kuperberg_invariant.jl).

The source distinguishes integrals on the algebra and its dual
([[Kup91#^kup91-blk-0079]]); translating that distinction requires specifying
the domain, rather than assuming the repository field names use a universal
integral/cointegral vocabulary.

## Bases and small-Borel styles

Group bases retain the supplied label/table order. Exterior bases are indexed
by bitmasks plus one, with generators in increasing bit order and Koszul
signs; use the graded bialgebra check. Kac–Paljutkin coordinates are ordered
`1,x,y,xy,z,xz,yz,xyz`. Small-Borel coordinates use `i+N*j+1` for `(i,j)` in
`0:N-1`, with `q=exp(2πim*M/N)` for the reduced input `M/N`.
The product coefficient of `(i1,j1),(i2,j2)` is `q^(-j1*i2)` when
`j1+j2<N`, with the first index reduced modulo `N`, and is zero otherwise.
Style `:L` has the opposite coproduct and inverse antipode of style `:R`.
Evidence: [hopf_invariant.jl](../../src/hopf_invariant.jl),
[finite_algebra.jl](../../src/finite_algebra.jl),
[topology_core.jl](../../test/topology_core.jl) (including exterior graded
checks and small-Borel style fixtures),
[hopf_properties.jl](../../test/hopf_properties.jl) (orders 2 and 5).

## Double coordinates

The basis pair `(f,a)` in `H* ⊗ H` has index `f+d*(a-1)`.
The product tensor is the six-leg expression below, with output order
`[f1,a1,f2,a2,fo,ao]=[1,2,3,4,5,6]`, combined into three basis-pair indices.

| Model | Defining coordinate expression |
|---|---|
| Heisenberg record | `De{5,1,7}De{2,9,8}Mu{7,9,3}Mu{8,4,6}` |
| Drinfeld record | `Mu{7,8,3}Mu{9,10,8}De{2,10,11}De{11,12,13}Si{7,13}Mu{12,4,6}De{5,1,9}`, where `Si=An^-1` |

The unit is `ε ⊗ 1`. Evidence:
[finite_algebra.jl](../../src/finite_algebra.jl),
[topology_core.jl](../../test/topology_core.jl).
question: unverified — identification of these products with a particular
literature opposite/copposite convention needs an explicit action/basis map.
Their names express the intended models; the coordinate definitions are the
established content.

## Convention dictionary

| Other convention and locator | Map to this convention | Check and boundary |
|---|---|---|
| Dual evaluation in [[Kup91#^kup91-blk-0023]] | Corresponding primal/dual bases; no conjugation | `hopf_pairing` and the basis-pairing checks in [topology_core.jl](../../test/topology_core.jl) |
| Legacy MATLAB expression ports, [tensor_network.jl](../../src/tensor_network.jl), `legacy_hopf_tensors` | Transpose `An`; legacy `Tp=De{1,4,0}Mu{3,4,2}` with output `[1,0,3,2]`; legacy `Tm=De{1,4,0}An{4,5}Mu{3,5,2}` with output `[0,1,2,3]` | [tensor_validation.jl](../../test/tensor_validation.jl), C250406 expressions; a total literature diagram correspondence is unverified |
| MATLAB exactly-twice index rule | Select `index_rule=:matlab`; the network extension retains all requested output variables | [tensor_validation.jl](../../test/tensor_validation.jl), multiply-occurring output fixture |
| Small-Borel `:R` coordinates | To `:L`: exchange coproduct outputs and invert antipode, retain multiplication/basis | [hopf_invariant.jl](../../src/hopf_invariant.jl), [topology_core.jl](../../test/topology_core.jl); no general source basis comparison claimed |
| Ungraded bialgebra identity vs exterior algebra | Insert the Koszul middle-factor sign using supplied degrees | [topology_core.jl](../../test/topology_core.jl), exterior(2) ordinary/graded checks |

## Kuperberg and cellular crossing signs

Status: **open** — the implementation and its integration test disagree.
The intended relative sign convention has not been decided.

The LC passage frame belongs to [LC §3](../specs/lickorish-code.md).
The cellular sign belongs to
[Koda conversion §7.3.1](../specs/heegaard-koda-conversion.md#731-dαdβ-の計算).
The comparison concerns coherent alpha/beta intersections retained by
`forget_reference(realize_hc(hc))` when its cellular diagram exists.

| Positive LC owner | Negative LC owner | HC Kuperberg crossing sign | Cellular `Intersection.sign` | HC antipode selected |
|---|---|---:|---:|---|
| alpha | beta | +1 | −1 | no |
| beta | alpha | −1 | +1 | yes |

The two cases follow exhaustively from `_ab_crossing_passage_positive` in
[kuperberg_invariant.jl](../../src/kuperberg_invariant.jl) and
`_crossing_records` in [koda_forget.jl](../../src/heegaard/koda_forget.jl).
`_build_heegaard_from_combinatorial` copies that sign and stores the original
integer label as its decimal `source_label`. The table describes the current
local branches; it does not settle which sign the Kuperberg entrance is meant
to use or identify its contraction with the source invariant.

A direct driver on 2026-11-26 used the actual functions and all five inputs
of the test set `HC and diagram entrances mark the same crossings negative`
in [kuperberg_invariant.jl](../../test/kuperberg_invariant.jl):

| Registered fixture | HC antipode labels | Cellular negative labels | Cellular positive labels |
|---|---|---|---|
| `261107-lickorish-data/group-15/case-01.hc` | 5 | 4, 6 | 5 |
| `261107-ld-draw-b/g12-c01.hc` | 7 | 6 | 7 |
| `261107-ld-draw-b/g15-c01.hc` | 5 | 4, 6 | 5 |

Both fixture directories are under `threads/PJ260922/data/` from the research
repository root. The remaining inputs, `g13-c01.hc` and `g31-c01.hc`, return
no cellular diagram and are skipped by that test. All three available
diagrams satisfy the observed positive-set equality. Both compared sets
in the test assertion already contain strings.

<a id="kuperberg-sign-contradiction"></a>

bug: The integration testset `HC and diagram entrances mark the same crossings negative`
at [test/kuperberg_invariant.jl:275](../../test/kuperberg_invariant.jl#L275)
asserts `hc_neg == diag_neg` at [line 293](../../test/kuperberg_invariant.jl#L293).
For `threads/PJ260922/data/261107-lickorish-data/group-15/case-01.hc`, the library returns
HC antipode labels `{5}`, cellular negative labels `{4,6}` and positive labels `{5}`;
for `threads/PJ260922/data/261107-ld-draw-b/g12-c01.hc`, these sets are `{7}`, `{6}` and `{7}`;
for `threads/PJ260922/data/261107-ld-draw-b/g15-c01.hc`, they are `{5}`, `{4,6}` and `{5}`.
Thus the assertion fails on all three convertible fixtures. This observation
is `verified-computationally` on those three fixtures, reproduced independently
by the lead; it establishes a contradiction between the code and test.

question: open — Which crossing-sign convention is the Kuperberg HC entrance
meant to use relative to the cellular diagram, and how would choosing it change
the values of `Z_Kuperberg` on these three inputs? The observed opposite sets
make an entrance-dependent antipode assignment the working explanation;
equal scalar values on a particular algebra would not decide the convention.
To decide it, match the HC alpha/beta families to the lower/upper circles of
[[Kup91#^kup91-blk-0111]] and their tensor ports in equations 26–28, including
both ordered tangent-basis orientations, and compare both entrances' contractions
on these fixtures with involutory Hopf algebras whose antipode is nontrivial.
An explicit family/port correspondence and the resulting value comparisons are
missing; they may explain a sign reversal or identify which entrance must change.
Neither equality nor inequality of the two entrances' invariant values follows
from the label-set observation alone.

Registered source: [Kup91](../../../../references/items/article/w-6c3fdcc0-6151-4b25-9a8f-444e375d62fd/Kup91.md).
