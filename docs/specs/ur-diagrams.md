# U(R) diagrams and trace expressions

## Definitions

### Weighted endpoint data

A legacy `URDiagram` is a coefficient `C`, a linear word `V` of signed edge
endpoints and zero separators, a five-row edge table `E`, and weights `W`.
Each nonzero absolute label names an edge with one positive and one negative
endpoint. Zero separators split the word into tensor sections; the rank is
one plus their number, even for an empty word. Row 1 of `E` lists edge labels;
the remaining rows are auxiliary data, not a separately validated geometric
structure. `W` is aligned with the columns of `E`.

This is the data model of the local source
[SRC-262, properties, constructor and validateData](../../../../references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m),
realized in [ur_diagram.jl](../../src/ur_diagram.jl).
Endpoint order, regulator names and trace signs belong to
[URDiagrams conventions](../conventions/URDiagrams.md).

### Formal operator word and regulators

The displayed weight word uses two formal slots of a symbol `WO_e`, one for
each endpoint of edge `e`, and tensor separators. It retains only the sign of
each weight, not its magnitude or `C`; it is not a faithful encoding of the
weighted diagram. `SymExpr` coefficients are
[formal stored expressions](algebra.md#expression-language).
A regulator is an additive perturbation of a base weight. Its symbolic name
does not specify a topology, analytic family, direction of approach or
limiting prescription. Those additional objects would be necessary to speak
of an analytic trace or a regularized invariant.

question: unverified — The intended operator interpretation is suggested by
the source's `WO` notation, but an operator space, domains and trace/limit
hypotheses are not established by the library or its tests. The repository
formulas below are established algebraic expressions; an operator-model
derivation with those hypotheses would establish the stronger interpretation.

### Legacy trace expressions

The single-edge expression, section matrix, rank-one block matrix and cyclic
adjustment are defined in
[trace conventions](../conventions/URDiagrams.md#legacy-trace-formulas).
`trace1`, `trace2` and `trace3` refer to these distinct expressions. A
singular denominator is a failure of finite evaluation, not permission to
assign a regularized value. `trace_value` uses a reduced legacy diagram and
the section-matrix expression.

### Exact MST expression

`URAlgebra` selects a different exact evaluator; `URMSTExpression` stores its
integer matrix. The matrix and the rational value `-1/det(M)` are defined by
the [exact MST convention](../conventions/URDiagrams.md#exact-mst-word-and-matrix).
It uses one closed, nonempty real-Gauss component and crossing signs. Assigned
edge weights do not enter this calculation. A singular matrix has no value
in this finite rational evaluator.

[[MST25#^mst25-blk-0025]] gives the Heisenberg-double setting and
[[MST25#^mst25-blk-0027]] states the framed/combed/topological invariant
hypotheses. Those passages do not identify this particular determinant with
that invariant for every input. The tests in
[ur_mst.jl](../../test/ur_mst.jl) compare absolute values with the independently
known manifold orders on five named examples, check cyclic shifts there,
and check singular behavior on three named examples. This is
`verified-computationally` in that finite scope.

### Standard matrix data and local reduction records

A standard matrix diagram for an `n×n` coefficient matrix has edge
`e=i+n(j-1)` from positive section `i` to negative section `j`, with weight
`M[i,j]`. Each section lists its positive endpoints and then its negative
endpoints in increasing label order. Consecutive sections are separated by
zero. Its coefficient is the supplied `C`. The empty matrix gives the empty
word. This is the coordinate construction in
[ur_diagram.jl](../../src/ur_diagram.jl), `set_data_cm!`, and
[C251030WeylAlgTraceFormula.jl](../../examples/C251030WeylAlgTraceFormula.jl).

A sign group is a maximal consecutive run of entries with the same sign
in the endpoint word; zero has its own sign. A zero-weight dilation edge is
an additional paired edge of weight zero. The operation named `put_zero!`
inserts this edge, **not a zero word separator**. Its positional convention
and the source difference are in the
[convention dictionary](../conventions/URDiagrams.md#convention-dictionary).

The legacy reduction names refer to the local weighted-endpoint relations
of [SRC-262](../../../../references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m),
`reduct3`–`reduct8`. `reduct3` merges pairs with both endpoint sign groups in
common. `reduct4` relates chords in adjacent opposite-sign groups to supplied
section dilations. `reduct5` is the rank-one dilation relation; `reduct8` is
the sectionwise counterpart. `reduct6` applies the first sign-inversion
exchange. `reduct7` is a bounded pass of these exchanges, with a Boolean
record for each encountered separator interval and the final visited interval.
Its bound is the initial word length, even if exchanges lengthen the word;
thus its record need not certify all sections of the final word. A `true`
entry means no exchange was made in that visited interval, not a global
normal-form certificate.

A reduction-history entry is `(iteration, vertices, complete)`: the round
number, copied endpoint word, and copied `reduct7` flags at that round's end.
It does not contain weights or coefficient, and precedes the final merge of
`simplify2!`. These are local reduction records, not proofs of operator or
manifold equivalence. The general trace-preservation status remains `open`
as stated below.

### Knot-complement word

The endpoint-word construction named `trans_knot_complement` is defined in
[the conversion convention](../conventions/URDiagrams.md#knot-complement-word).
It creates word data; it does not itself supply crossing signs, weights,
an embedding, a boundary triangulation or a homeomorphism certificate.

question: unverified — A general identification of that word with a spine or
diagram of the input knot's complement is intended by the source name but
not established by the structural/trace fixtures. A local replacement proof
and a boundary/meridian convention are missing. Such a proof, or a failed
complement comparison, would settle that identification.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| U1 | Every positive endpoint is unique, every negative endpoint is unique, and the two absolute-label sets agree. | `validate_data` checks; setters do not. |
| U2 | The positive labels in increasing order equal `E[1,:]`, with one weight per edge. | `validate_data` checks; auxiliary rows have no mathematical validity check. |
| U3 | Coefficients support the required arithmetic and every denominator used is nonsingular when a finite scalar is claimed. | Arithmetic may fail or return infinities; not an analytic test. |
| U4 | A swap acts on two existing adjacent positions of the linear word. | Checked unless `skip=true`; not cyclic adjacency. |
| U5 | An edge sum has existing adjacent positive endpoints and adjacent negative endpoints. Distinct edges and coherent data are caller obligations. | Local adjacency guard unless skipped. |
| U6 | A dilation and its target occupy three consecutive existing positions; two composed dilations occupy four consecutive existing positions. | Local block guard unless skipped; no global operator-equivalence certificate. |
| U7 | A rank-one block trace has no zero separators; its edge/word data are coherent. | Rank checked; other validity is caller-supplied. |
| U8 | Exact MST input has one nonempty closed real component, every signed crossing label once and signs `±1`. | Exact conversion checks its supported domain; full closed-normal conditions are separate. |
| U9 | Exact finite rational evaluation requires a nonsingular trace matrix. | Determinant zero raises `DomainError`. |
| U10 | Endpoint insertion positions lie in `0:length(V)` in the pre-insertion word; labels and new label fit `Int`, and the existing data satisfy U1–U2. | Position bounds are checked. Label overflow and structural consistency are caller obligations. |
| U11 | A standard coefficient matrix is square and entries and coefficient support the arithmetic in U3. | Squareness is checked. Entry arithmetic and finite-denominator conditions are caller obligations. |
| U12 | Local reductions use U1–U3 and the relevant U4–U6 placements. Supplied dilation labels identify existing section dilations in section order. Iteration bounds are nonnegative integers. | No global preflight check: local placement and missing-dilation checks can fail after mutation. Bounds and consistency are caller obligations. |

question: open — General trace preservation by every legacy rewrite, complete
reduction within a finite iteration bound, equality of the three unreduced
trace expressions, and equality with the exact MST evaluator are not proved.
The registered reduction fixtures support their tested cases; a symbolic
identity proof with U1–U7 and the chosen regulator domain is the missing
general check. In particular, section/block matrix entries overwrite repeated
positions whereas the exact MST crossing contributions add.

## Operations

`URDiagram`, `set_data_cvw!`, `set_data_cvs!`, `set_data_v!`, `set_data_cm!`,
`validate_data`, `check_data`, `swap!`, `delete_edge!`, `add_edges!`,
`dilation!`, `compose_dilation!`, `simplify!`, `simplify1!`, `simplify2!`,
`trace1`, `trace2`, `trace3`, `trace_value`, `standard_matrix_adjustment`,
`standard_matrix_trace`, `odata`, `weight_expression`, `unicode_art`,
`trans_knot_complement`, `vl_to_urdiagram`, `ur_invariant`,
`URAlgebra`, `URMSTExpression`, `invariant_expression`, `simplify_invariant`,
`evaluate_invariant`, `Z_MST`, `ur_prime_invariant`.

## Sources

Registered source: [MST25](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/MST25.md).
question: unverified — SRC-262 is a local program source, not a resolved
literature cite key. Its exact file and method names are linked above rather
than assigning it an invented key. The absent operator/limit theorem is not
supplied by the generic MST25 invariant statement.
