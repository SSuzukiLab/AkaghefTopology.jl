# Closed and relative spin weights

## Definitions

### Closed spin normal o-graph

Start with a closed normal o-graph satisfying
[N1, N2, C1–C3](ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3).
Give its edges weights in `Z₂`. For each source circuit with solid dots,
read the edge weights along the circuit with multiplicity. Condition S in
[[Mih23#^m23-blk-0020]] requires their sum modulo two to equal one plus half
the number of solid dots. The spin structure is the extension over the dual
2-skeleton of the weighted trivialization over the dual 1-skeleton, as stated
in [[Mih23#^m23-blk-0021]]. This geometric meaning requires the source circuit
interpretation and closed-normal hypotheses in addition to an equation check.

### Repository coordinates

The ordered edge vector is `z`, with canonical entries 0 or 1 and weight mode
`:z2`. The disk table supplies incidence matrix `delta`, dot counts `Ndots`
and closed/open row flags. On closed rows the coordinate equation is
`delta*z = 1 + Ndots/2 (mod 2)`; division of the even integer dot count
precedes reduction modulo two. These are the closed coordinates of condition S.

### Relative boundary data

An open disk path has an ordered pair of boundary endpoint identities
`(boundary_start,boundary_end)`. Each endpoint includes its strand, incoming
or outgoing side, and circuit lane. A relative target is an explicitly
supplied element `t_r` of `Z₄` for every open row. Closed rows retain target 2.
The measured value is `m_r = 2*(delta*z)_r - Ndots_r (mod 4)`.
The relative equation is `m_r=t_r`. The divided binary system has matrix
`delta mod 2` and right side `(t_r+Ndots_r)/2 mod 2`; it exists only when
that numerator is even.

Targets may be row-aligned or keyed by row number or the directed endpoint
pair. Endpoint identities do not follow from a zero cut marker alone. Their
transport is owned by [weight transport](virtual-link-moves.md#121-supported-weighted-moves).
Measuring a boundary value does not choose an independent boundary condition.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| S1 | Assigned weights have mode `:z2`, one canonical 0/1 value per real edge, and valid Gauss data. | Spin validation checks supported data. |
| S2 | Closed rows obey condition S in the above coordinates. | Residual checked; closed-normality and source circuit identification are separate. |
| S3 | Every open row has exactly one explicit exact-integer target modulo four; no closed row is reassigned away from 2. | Target validation; booleans are not integer targets. |
| S4 | Each target plus its dot count is even before division by two. | Checked; odd parity is an obstruction, not rounded. |
| S5 | Relative measured values equal the specified targets modulo four. | Relative residual check. |

## Coordinate choices and evidence

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)** for canonical
0/1 storage, directed endpoint order, target 2 on closed rows and the doubled
relative equation. Evidence: [spin_weights.jl](../../src/spin_weights.jl),
[spin tests](../../test/spin_weights.jl),
[cut spin tests](../../test/cut_spin_weights.jl), including odd dot-count
parity rejection and the explicit target vector `[0,2,2,0]`.
Integer framing equations and their storage range remain owned by
[virtual-link-moves §12](virtual-link-moves.md#12-weight-transport).

question: open — The relative equation has an established algebraic meaning,
but its equivalence to extension of a specified geometric boundary spin
structure is not established. The expected interpretation is relative spin
extension, since doubling recovers the closed source equation when `t=2`;
the missing data are a boundary trivialization and its identification with
each directed `Z₄` target. A relative obstruction derivation could confirm
that interpretation or require a translation of targets.

## Operations

### Equation records

`spin_boundary_data` returns `mode=:z2`, `modulus=4`, the assigned edge vector,
`delta`, `Ndots`, `isCycle`, `cp2`, `measured4`, endpoint arrays and
`boundary_identities`. A closed row's boundary identity is `nothing`; an open
row's identity is its ordered endpoint pair. This is measurement, without
choosing a prescribed target.

`spin_weight_equation` uses a modulus-two record when all rows are closed:
`delta` is the binary matrix, `rhs` its binary right side, and `residual` is
the coordinate equation residual. When any row is open, its record has
`modulus=4`, integral disk incidence as `delta`, prescribed doubled targets
as both `rhs` and `boundary_targets`, the measured doubled residues as
`measured4`, and directed endpoint arrays. Its residual is the measured value
minus the target modulo four. The assigned edge vector is returned in both
cases. Neither calculation modifies its input. The shape distinction and the
relative target vector `[0,2,2,0]` are exercised in
[cut spin tests](../../test/cut_spin_weights.jl).

### Binary solution records

`solve_spin_weight_system(A,b)` accepts one-based matrices/vectors with exact
integer entries and one right side per row. Bool entries are accepted; entries
are reduced modulo two. This input convention is already embodied by the
integer type check in [spin_weights.jl](../../src/spin_weights.jl); floating
coefficients are rejected in [spin tests](../../test/spin_weights.jl).

The schema is `topology-spin-weight-solution-v1`, with `mode=:z2`, `modulus=2`,
reduced `A` and `b`, and status `:spin_solution` or `:inconsistent`.
`particular` is a binary vector in the feasible case and `nothing` otherwise.
The columns of `kernel` are a full binary kernel basis in both cases; all
solutions are `particular + kernel*z` modulo two, uniquely for binary `z`.
The certificate contains augmented `rref`, a binary `row_transform` and
one-based pivot columns. On failure, `obstruction.row_combination` annihilates
`A` modulo two and pairs with `b` to give `obstruction.rhs=1`.

`solve_spin_weights` supplies this system from the diagram and attaches
`boundary_targets`. Its assignment status is `:not_requested`,
`:not_spin_feasible`, or `:assigned`; assignment changes only weights and the
assigned flag after validation. Mode remains Z2. The matrix record has
algebraic semantics; geometric spin interpretation additionally requires the
closed-normal and circuit hypotheses above.

Status: **proved/agent** by exact row operations and the free-coordinate
parameterization in [spin_weights.jl](../../src/spin_weights.jl). Independent
finite solution-set and witness comparisons are in
[spin tests](../../test/spin_weights.jl).

### Spin CP transport

The implemented local coordinates follow the CP and H rules labelled
[[Mih23#^m23-fig-0094]] and [[Mih23#^m23-fig-0093]]. The two binary cases use
the incoming-edge occurrence convention of the move catalogue. The rule data,
including the RHS vertex-role order, are in
[spin_cp.jl](../../src/spin_cp.jl); the exact edge arrays and local boundary
equations are compared in [spin CP tests](../../test/spin_cp.jl).

| Condition | Current enforcement |
|---|---|
| S6 | Assigned canonical Z2 weights and a CP candidate for one source interval; geometric site/templates match the current source. |
| S7 | An uncut source and target pass closed-normal certification and the spin equation; explicit boundary targets on an uncut source are rejected. |
| S8 | A cut source has open disk paths and valid relative targets. A missing source prescription is captured from its measurements; an explicit prescription must agree. |
| S9 | Internal edge occurrences match one of the two local cases, either directly or after the requested source-local H normalization. |

`normalize=false` requires a direct case match. `normalize=true` selects the
first matching local H mask in increasing integer mask order, trying case zero
before case one for each mask. The operation is H followed by CP, so inverse
CP restores the normalized representative, not necessarily the original raw
weights. This selection is existing behavior, not a claim of minimal H length.

The returned `target` is the resulting Z2 diagram; `case` and `direction`
identify the rule; `target_vertices` and `target_locations` identify its
output roles and matching intervals. `source_h_vertices` records the toggles,
`normalized_source` their result, `source_already_normalized` whether that
list is empty, and `source_figure` is the source locator. `modulus=2` records
edge coefficients; `boundary_modulus` is four on cut inputs and `nothing`
otherwise. `source_boundary`, `target_boundary` and `boundary` retain the
prescriptions; `boundary_transport` is the directed endpoint map for a cut
input and `nothing` otherwise. A supplied target map is an additional check,
not a replacement for the transported source prescription.

Status: **proved/agent** for the equation and directed-target guarantees by
the validation performed before returning a target. Finite closed and relative
comparisons, including H gauges and inverse moves, are in
[spin CP tests](../../test/spin_cp.jl) and
[cut spin CP tests](../../test/cut_spin_cp.jl). The open geometric
boundary-extension identification remains the question in
[Coordinate choices and evidence](#coordinate-choices-and-evidence).

`spin_boundary_data`, `spin_weight_equation`, `validate_spin_weight`,
`is_spin_weight`, `solve_spin_weights`.

## Source

[Mih23](../../../../references/items/article/w-f2194518-7624-4b07-9eb5-3e04e0a7d036/Mih23.md),
definition of spin normal o-graph, condition S and its following extension
paragraph. The open-path target model is a repository extension of that
closed statement, not a theorem attributed to Mih23.
