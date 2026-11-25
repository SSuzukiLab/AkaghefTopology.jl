# Hopf tensor contractions and their geometric targets

## Definitions

### Structural-tensor contraction

A finite Hopf record supplies the structural maps and crossing tensors of
[the algebra specification](algebra.md). A diagram tensor expression is a
finite collection of these coordinate tensors with identified leg labels;
its scalar contraction sums every internal basis index. Its mathematical
meaning without additional hypotheses is that finite multilinear contraction.
The coordinates, dual pairing and legacy ports are owned by
[Algebra conventions](../conventions/Algebra.md).

### MST geometric target

[[MST25#^mst25-blk-0025]] associates positive crossings with the canonical
Heisenberg-double tensor. [[MST25#^mst25-blk-0027]] states that the literature
construction is an invariant of framed 3-manifolds with vanishing first Betti
number; involutivity gives the cited combed specialization, and unimodularity
and counimodularity give the stated topological specialization. Those are
hypotheses on the represented geometric object and algebra, not consequences
of obtaining a scalar from a coordinate network.

The library's `VirtualLink` route associates crossing/antipode tensors to real
Gauss data and integer edge weights. The source closed-normal hypotheses are
[N1, N2, C1–C3](ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3);
the framing equation belongs to [weight transport](virtual-link-moves.md#12-weight-transport).
The evaluator's closed-strand and integer-weight guards do not certify all
these geometric hypotheses.

### Exact Taft route

`TaftAlgebra` selects a symbolic/exact tensor calculation at a specified root
of unity. A root specification uses the reduced rational phase and the
associated cyclotomic coefficient domain. A finite record bearing the Taft
kind tag is eligible for the exact legacy bridge only when its structural
tensors match the supported unmodified `:L` prototype. Arbitrary mutation of
a record is not ignored by the dispatch. Reduction strategies and term
budgets control algebraic expression work, not a wall-clock bound.
The field/basis choices are [small-Borel styles](../conventions/Algebra.md#bases-and-small-borel-styles).

The concrete root record has integer order `N >= 2` fitting `Int` and
primitive exponent `M`, reduced modulo `N`, with `gcd(N,M)=1`.
`TaftAlgebra` stores this pair and style `:L` or `:R`. A symbolic order is an
integer-variable term (a Julia `Symbol` is converted to one); it permits
only `M=1`. The current network uses its common symbolic order variable `N`
and supports style `:L` only. These are the coordinate domains of
[RootSpec](../../../algebra.jl/src/CyclotomicValues.jl) and
[TaftAlgebra](../../src/taft_invariant.jl), not a proof of a Hopf model.

An `MSTTensorExpression` stores a vector of `SparseFactor` coordinate tensors.
The generic numeric route evaluates their closed contraction. The exact Taft
expression is the `TaftExpression` record owned by
[TaftNetworks](../../../algebra.jl/src/Hopf/TaftNetworks.jl).

### Kuperberg coordinate target

[[Kup91#^kup91-sec-0105]] defines the invariant by contracting the two circle
families of a Heegaard diagram, inserting the antipode at negative tangent
bases. The repository Heegaard tensor route assigns iterated coproducts of
the right-cointegral element to alpha circles and the right-integral
functional on word products to beta circles. The antipode output is the beta
leg and input the alpha leg. For pairing `p=λ(Λ)`, its normalization is
`p^(-g)` times the raw contraction, except that a pairing approximately one
at `atol=rtol=1e-10` uses normalization exactly one. Contraction and pairing
are evaluated in complex floating arithmetic. This is the coordinate definition in
[kuperberg_invariant.jl](../../src/kuperberg_invariant.jl); the comparison to
the source's dimension normalization remains a separate obligation.

`KuperbergTensorExpression` stores the input Heegaard code `hc`, its
alpha/beta incidence `link`, a list of coordinate `factors`, a scalar
`normalization` and the common algebra `dimension`. Its finite meaning is
the closed contraction of those factors times that normalization. The raw
field constructor does not validate the diagram, tensor legs, dimension or
normalization, and correctly typed mutable fields remain shared.

### Heegaard word group and antipode labels

For a coherent HC, sort its alpha and beta curve names lexicographically.
The coordinate group has one generator `x` followed by the alpha name for
each alpha curve. Each beta word gives one relator in its stored cyclic
order, retaining only alpha/beta intersections. A passage whose positive
owner is beta contributes exponent `+1` of its alpha generator; a passage
whose negative owner is beta contributes exponent `-1`. There is no free or
cyclic reduction of these relators. This defines the group of the stored
presentation in [kuperberg_invariant.jl](../../src/kuperberg_invariant.jl),
`heegaard_pi1_presentation`; the abstract presentation record is owned by
[Heegaard terminology](heegaard-terminology.md#91-群と不変量で使う用語).

The negative-antipode label set selects alpha/beta intersections whose
positive owner is beta and negative owner is alpha. It omits other curve
families. This is the repository's Kuperberg passage convention, not the
sign of a Gauss visit alone. The ambient LC passage frame is owned by
[LC §3](lickorish-code.md); the source target is negative tangent bases in
[[Kup91#^kup91-sec-0105]]. The current source-sign comparison remains
`unverified` in H3. The [HC/cellular sign contradiction](../conventions/Algebra.md#kuperberg-sign-contradiction)
records that this antipode-label set equals the cellular positive set on the
three convertible fixtures, while the integration test asserts equality with
the negative set. The intended relative sign convention is open. Its resolution
requires a source family/port comparison and a comparison of `Z_Kuperberg`
values on those inputs; the observed relation is not a fixed convention.

The intended group is the fundamental group of the closed 3-manifold of
a complete Heegaard splitting. The coordinate-to-geometric identification
requires complete cutting systems and the passage-sign correspondence;
`heegaard_pi1_presentation` checks neither. The generic Heegaard meaning is
the linked terminology's one, and the group-algebra target is
[[Kup91#^kup91-blk-0140]]. The outstanding proof must identify the sorted
alpha generators with handlebody generators and the signed beta words with
the attaching circles of the complete-cut-system manifold in
[Koda conversion §5.1](heegaard-koda-conversion.md#51-曲線系と多様体),
for both local passage orientations, the stored order around each beta
circle, and arbitrary genus. The circle-family and sign assignment must
also be compared with the source's lower-circle states and upper-circle
group constraints. This general dictionary is `unverified`. Registered S3,
finite-group and T3 homomorphism-count comparisons support the intended
identification but do not prove it. A contradictory attaching-word or
generator correspondence would overturn it. The opposite crossing-sign
conventions alone do not refute a group isomorphism.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| H1 | Input structural tensors satisfy the algebra/Hopf axioms in the chosen (ordinary or graded) category and coordinate convention. | Record construction does not certify this; finite numerical predicates have their stated scope. |
| H2 | The diagram and assigned weights satisfy all geometric hypotheses of the intended literature invariant. | Scalar evaluators check only their stated input guards; callers supply the stronger hypotheses. |
| H3 | Every crossing and antipode leg matches the source diagram dictionary, including orientation, dual and normalization choices. | Repository port formulas are fixed; the complete literature comparison is unverified. |
| H4 | Exact Taft evaluation uses a valid root specification and supported coefficient/record model. | Exact route validation; a finite kind tag alone is insufficient. |
| H5 | A Kuperberg Heegaard calculation uses coherent intersection/circle data and nonzero integral/cointegral pairing. | The evaluator checks its supported record and pairing; invariant identification still requires H1–H3. |
| H6 | A stored Kuperberg expression's HC, incidence link, tensor factors, dimension and normalization agree with one supported calculation and its closed contraction. | Caller responsibility for raw field construction; not checked by field conversion. |
| H7 | A Taft selector has the concrete or symbolic order, primitive exponent and style described above. Evaluation uses a concrete order and a nonnegative `Int` term budget. | Root/selector and evaluator guards check these domains; network evaluation further restricts style to `:L`. |
| H8 | HC curve roles and words give coherent two-sided passages; alpha and beta families are complete cutting systems of the represented Heegaard splitting when a manifold group is claimed. | The word-group and negative-label queries reject conflicting same-sign owners but do not validate HC or cutting-system completeness; missing words and malformed data can fail or be omitted. |

## Evidence and unresolved identifications

The finite route is in [hopf_link.jl](../../src/hopf_link.jl) and
[hopf_invariant.jl](../../src/hopf_invariant.jl); exact evaluation is in
[taft_invariant.jl](../../src/taft_invariant.jl).
[taft_exact_route.jl](../../test/taft_exact_route.jl) tests the one-crossing
weighted fixture, root-dependent values, route agreement and mutation
fallback. [kuperberg_invariant.jl](../../test/kuperberg_invariant.jl) tests
the Heegaard route. These checks establish their finite cases, not all
manifold invariance claims.

question: unverified — The general identification of each implemented
contraction with its literature invariant requires a complete port/basis map
and normalization comparison. The intended correspondence is supported by
the structural construction and tested examples; missing antipode/dual
convention comparisons or a failed move check could invalidate it. Kuperberg
source equations 26–28 and the dimension-versus-pairing normalization are
specifically not certified by the current source comment.

question: unverified — Numerical contractions may prune coefficients at the
repository tolerance and do not certify exact algebraic identities. Exact
Taft arithmetic removes that numerical issue in its domain, but does not
by itself settle H2 or H3.

## Operations

`hopf_invariant`, `invariant_expression`, `simplify_invariant`,
`evaluate_invariant`, `Z_MST`, `taft_reduce`, `taft_evaluate`,
`Z_Kuperberg`, `kuperberg_tensor_expression`, `kuperberg_hc_negative_labels`.

## Sources

[MST25](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/MST25.md),
[Kup91](../../../../references/items/article/w-6c3fdcc0-6151-4b25-9a8f-444e375d62fd/Kup91.md).
