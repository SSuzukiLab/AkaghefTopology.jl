# U(R) diagram conventions

This page owns the legacy UR endpoint ordering, formal regulators, algebraic
trace formulas, exact MST word convention and knot-complement word.
The objects and conditions are in [ur-diagrams.md](../specs/ur-diagrams.md).
Every coordinate choice below is **fixed 2026-11-25 Codex lane lib-api-rest
(rule a)**, following the existing code and tests. Analytic and geometric
identifications retain the evidence statuses stated on the specification.

## Endpoints and legacy VL2URD ordering

Positive endpoint `+e` is displayed as slot 2 of `WO_e`, negative endpoint
`-e` as slot 1. A negative or zero weight selects the `bar(...)` display;
a positive weight selects the unbarred display. Zero in the word is `⊗`.
The sign of a crossing is separate from either endpoint sign.

Legacy `VL2URD` takes real-Gauss components in their stored component order,
reverses each component's word, and inserts one zero between consecutive
components. It retains the crossing-label order for the base weights, which
are the real crossing orientations; `C=1`. Explicit virtual crossings and
assigned edge weights are excluded. In particular it does not cyclically
rotate the reversed word to a preferred positive endpoint.
Evidence: [virtual_link.jl](../../src/virtual_link.jl), `vl_to_urdiagram`, and
[topology_core.jl](../../test/topology_core.jl), C251030 conversion fixtures.

## Regulators

`set_data_cvs!` adds one perturbation to each base weight. With no supplied
epsilon, edge label `e` uses formal symbol `eps` followed by its decimal label.
A supplied scalar is shared by value across all additions; a supplied array
is aligned with the increasing edge labels. `set_data_v!` instead initializes
unknown weights with `NaN`.

Julia's legacy trace and odata operations do not take epsilon limits. A
numerical choice such as `1e-5` is evaluation at that perturbation, not a
limit. Evidence: [ur_diagram.jl](../../src/ur_diagram.jl) and
[topology_core.jl](../../test/topology_core.jl), C251030 and figure-eight
complement fixtures. The source comparison is
[SRC-262, setDataFromCVS, trace1, trace2, getOdata](../../../../references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m).

## Legacy trace formulas

For a single stored weight `w`, the expression is `-C/w`.
For the section matrix, the section of an endpoint is one plus the number of
zeros up to its position. The entry at `(section(+e),section(-e))` is `W_e`.
If several edges have the same matrix position, the last edge in table order
supplies that entry. Other entries are zero. The section trace is `C/det(-M)`.

For a rank-one word, its block number starts at one and increments whenever
a positive endpoint starts a new positive run after a negative endpoint.
Negative endpoints retain the current block. The block matrix uses the same
positive-row/negative-column assignment and overwrite rule.
For an `n×n` matrix, the standard adjustment `A` is zero for `n=1`; for
`n>1`, `A[i,i]=1` and `A[i-1,i]=-1`, with cyclic row indices. It is empty
for `n=0`. The block/standard-matrix expression is `-C/det(M+A)`.
These formulas are algebraic definitions, not assertions of trace-classness.
Evidence: [ur_diagram.jl](../../src/ur_diagram.jl),
[topology_core.jl](../../test/topology_core.jl), symbolic rank-one,
C251030 standard-matrix and C250811 reduction fixtures.

## Exact MST word and matrix

Let the one closed real-Gauss word have each `+i,-i` once, for `i=1,...,n`,
and crossing signs `s_i`. Reverse it and cyclically rotate it so its first
positive entry is first. For each signed entry `x`, let `b(x)` be the number
of positive entries up to and including it in that word. Thus `b` takes values
in `1:n`. Define the integer matrix

`M = I - P + sum(s_i * E[b(+i),b(-i)], i=1:n)`,

where `P` shifts the rows of the identity one step forward cyclically, and
`E[a,b]` has one at row `a`, column `b`. Repeated crossing contributions add.
The exact finite value is `-1/det(M)` in `Rational{BigInt}`; zero determinant
is rejected. The supported simplification strategies leave this matrix
unchanged in value.

Evidence: [ur_mst_link.jl](../../src/ur_mst_link.jl),
[ur_mst.jl](../../src/ur_mst.jl), and
[ur_mst tests](../../test/ur_mst.jl). The latter includes five named
manifold-order fixtures, all their cyclic shifts, and a five-member bridge
family whose signs alternate. Consequently a universal positive sign or
orientation-independent signed manifold value is not established by these
tests. [[MST25#^mst25-blk-0027]] alone supplies no formula equating this
matrix with the literature invariant.

## Knot-complement word

For each signed input visit `v`, define two ordered pairs `p(v),q(v)`:

| Visit | `p(v)` | `q(v)` |
|---|---|---|
| `v>0` | `(4v-2,4v-3)` | `(4v,4v-1)` |
| `v<0` | `(4v+3,4v+1)` | `(4v+2,4v)` |

Concatenate the `p` pairs in input order, then one zero, then the reverse of
the concatenated `q` pairs. The result has length `4*length(code)+1`.
The input is normally a paired nonzero signed knot word; the conversion
itself only performs integer conversion and does not certify that domain.
The complement workflows use the four-weight pattern `(1,-1,-1,1)` per
crossing, with separately chosen perturbations.
Evidence: [ur_diagram.jl](../../src/ur_diagram.jl), `trans_knot_complement`,
[topology_core.jl](../../test/topology_core.jl), figure-eight complement and
C251030 knot variants. A geometric complement theorem remains unverified.

## Convention dictionary

| Other convention and locator | Map to this convention | Check and boundary |
|---|---|---|
| SRC-262 `put0(argP,argN)` endpoint insertion | Julia `put_zero!` uses the two original-word gap positions, inserts `+e` at the positive gap and `-e` at the negative gap; equal gaps put `+e` first. For `argP > argN`, the source puts the signs in the opposite order. Exchange the two new endpoint signs to compare that branch. Neither version inserts a zero separator. | Direct branch comparison of [SRC-262 put0](../../../../references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m) and [ur_diagram.jl put_zero!](../../src/ur_diagram.jl). The registered reduction fixtures use equal gaps; unequal-gap operator equivalence is unverified. |
| Real-Gauss component traversal, [VirtualLink data](../specs/virtual-link-data-model.md) | Reverse each word; preserve component order; insert zero separators | C251030 conversions in [topology_core.jl](../../test/topology_core.jl) |
| [SRC-262](../../../../references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m), trace1/trace2/getOdata epsilon limits | Julia returns the expression before those limit operations | Direct code comparison; numerical epsilon fixture in [topology_core.jl](../../test/topology_core.jl). No universal interchange-of-limits claim |
| SRC-262 trace2 comment: sign provisional | Preserve the current `C/det(-M)` formula; do not infer a source-certified topological sign | Registered rank-one reductions give `-1` and `+1`; [topology_core.jl](../../test/topology_core.jl) |
| Legacy block adjustment | Exact MST uses the transpose cyclic adjustment, crossing sections counted by each positive entry, and additive entries | Formulas above and the two checked-in source files; no general equality map has been established |
| Legacy regulated trace vs exact MST | A proposed comparison must specify the word rotation, reduction, entry convention and regulator limit | `open`: separate fixture suites pass their own definitions; equality for arbitrary input is not tested or proved |

The `put0` dictionary check was run on 2026-11-25 with initial word
`[1,-1]` and gap pairs `(0,2)`, `(2,0)`, `(2,2)`. Julia returned respectively
`[2,1,-1,-2]`, `[-2,1,-1,2]`, `[1,-1,2,-2]`. Literal evaluation of the
source's indexing formulas gave `[2,1,-1,-2]`, `[2,1,-1,-2]`,
`[1,-1,2,-2]`. Exchanging only the two new signs in the reversed-gap case
matched all three words. This verifies the coordinate translation on each
positional branch; it does not prove an operator identity.

Registered source: [MST25](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/MST25.md).
