# Finite algebras and tensor expressions

This page owns the algebraic data definitions used by the `algebra` object.
Coordinate choices live in [Algebra conventions](../conventions/Algebra.md).
Contracts live in the [algebra index](../api/algebra/INDEX.md).

## Definitions

### Coordinate tensors

For an ordered basis of a finite-dimensional vector space, a coordinate tensor
is an array indexed by one basis index per leg. A `SparseFactor{T}` records a
finite map from index tuples to coefficients of type `T`, together with an
ordered list of integer leg labels. Unstored coefficients are zero. A scalar
has the empty index tuple. Leg labels are contraction variables, not basis
indices. Repeated labels identify the corresponding indices; labels absent
from the output are summed. This is the repository's coordinate realization
of the tensor contractions in [[Kup91#^kup91-sec-0018]] and
[[Kup91#^kup91-blk-0024]].

### Algebra records

A `FiniteAlgebra` records a dimension, multiplication tensor, unit tensor and
name. An `AlgebraElement` is a coefficient vector relative to that record's
basis. A `FiniteHopfAlgebra` records multiplication, unit, coproduct, counit,
antipode, positive/negative crossing tensors and four distinguished
integral/cointegral fields. These are containers for structural tensors;
constructing a container does not assert the algebra or Hopf axioms.
The mathematical structural maps are those of [[Kup91#^kup91-sec-0038]];
the field/port dictionary is [owned here](../conventions/Algebra.md#tensor-ports).

A dual record uses the dual basis. Evaluation of a dual vector on a vector is
bilinear, as in [[Kup91#^kup91-blk-0023]], not Hermitian. The two double
constructors return records on the tensor space `H* ⊗ H`, with the repository
products specified by [double coordinates](../conventions/Algebra.md#double-coordinates).

### Expression language

A tensor expression is a concatenation of factors `Name{i,j,...}` or
`(Name^p){i,j}`, with whitespace between factors allowed. Names consist of a
letter or underscore followed by word characters; leg labels are nonnegative
decimal integers. A decimal integer power, optionally prefixed by `-`, is a matrix power and requires two
legs. The unpowered parenthesized form `(Name){i,j}` is also accepted.
An output list orders the retained leg labels. This is a repository language,
not notation specified by Kup91; its checks are embodied by
[tensor_validation.jl](../../test/tensor_validation.jl).

`SymExpr` is a stored text expression with the repository's elementary
arithmetic formatting and simplifications. Text equality is not equality in
a formally specified symbolic field. It supplies formal coefficients for the
[U(R) diagram](ur-diagrams.md); it is not a computer algebra system.

### Algebra tests and named families

An involutory Hopf algebra has antipode square equal to the identity
([[Kup91#^kup91-blk-0055]]). The numerical predicate compares the supplied
antipode matrix with this identity. The repository's unimodularity predicate
instead compares its two stored integral vectors for nonzero proportionality;
the [field names](../conventions/Algebra.md#integral-names) determine which
side of the duality this tests.

The group-algebra family has the ordered group elements as basis and the
usual product, diagonal coproduct, unit, counit and inverse antipode
([[Kup91#^kup91-blk-0052]], [[Kup91#^kup91-blk-0054]]).
The exterior family uses bit-subsets of generators and the graded tensor
product. The Sweedler, Kac–Paljutkin and small-Borel families are the explicit
coordinate models in [hopf_invariant.jl](../../src/hopf_invariant.jl) and
[finite_algebra.jl](../../src/finite_algebra.jl).

question: unverified — A complete literature basis dictionary for the latter
three named models and the two doubles remains absent. Their intended
identification is plausible from their structural tensors and finite Hopf
checks, but a basis/antipode comparison with a registered source is missing;
a failed comparison would refute that identification.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| A1 | Each tensor key has one index per leg, within `1:dimension`; coefficients support the requested arithmetic. | Raw tensor/record constructors do not certify this. |
| A2 | Multiplication and unit satisfy associativity and the two unit laws. | `verify_hopf` checks the supplied Hopf record numerically; `FiniteAlgebra` is unchecked. |
| A3 | Coproduct/counit satisfy coassociativity and the two counit laws; antipode satisfies both antipode laws; multiplication/coproduct satisfy bialgebra compatibility. | `verify_hopf` checks these nine named equations together with A2. It does not check every unit/counit compatibility or distinguished-vector axiom. |
| A4 | A graded check has one integer degree per basis element and uses the Koszul sign for the middle tensor factors in bialgebra compatibility. | Caller supplies degrees; this does not validate a grading. |
| A5 | Expression syntax is valid, outputs are distinct, tensor names/ranks agree when a dictionary is supplied, and powers have two legs. | Expression validator; dimensions and invertibility need separate checks. |
| A6 | `:matlab`: every label occurs exactly twice in the expression plus output list. `:network`: each output occurs in the expression; other occurrence counts are unrestricted. | Selected explicitly; validator defaults to `:matlab`, evaluator to `:network`. |
| A7 | A group table is square, closed on its indices, associative, with an identity and two-sided inverses. | Group constructor checks; a callback must also close on distinct supplied labels. |

`verify_hopf` uses absolute coefficient tolerance `1e-10` by default.
`is_involutory` uses `atol=rtol=1e-10`; finite floating checks establish only
their tested equations at the chosen tolerance. Float contraction may discard
products of magnitude at most `1e-12`. Exact coefficient preservation is tested
separately in [tensor_coefficients.jl](../../test/tensor_coefficients.jl).

## Operations

`SparseFactor`, `FiniteAlgebra`, `FiniteHopfAlgebra`, `AlgebraElement`,
`basis_element`, `dense_tensor`, `label_tensor`, `contract_network`,
`validate_tensor_expression`, `evaluate_tensor_expression`, `hopf_tensors`,
`legacy_hopf_tensors`, `hopf_pairing`, `dual_algebra`, `heisenberg_double`,
`drinfeld_double`, `group_algebra`, `cyclic_group_algebra`, `exterior_algebra`,
`sweedler_algebra`, `kac_paljutkin_algebra`, `uqsl2_borel_small`, `verify_hopf`,
`is_involutory`, `is_unimodular`, `SymExpr`, `sym`.

## Source resolution

Registered texts: [Kup91](../../../../references/items/article/w-6c3fdcc0-6151-4b25-9a8f-444e375d62fd/Kup91.md),
[MST25](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/MST25.md).
Local representation definitions above are established by their linked code
and tests, rather than attributed to those papers.
