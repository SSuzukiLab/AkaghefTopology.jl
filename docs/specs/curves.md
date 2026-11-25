# Polynomial curve coordinates

## Definitions

A degree-`n` Bernstein basis on the unit interval has entries
`B_i,n(t)=binomial(n,i)*(1-t)^(n-i)*t^i`, for `i=0,...,n`.
For ordered control points `P_0,...,P_n` in a common coordinate space, the
Bezier curve is `sum(P_i*B_i,n(t), i=0:n)`.
These formulas define the local coordinate operations in
[helpers.jl](../../src/helpers.jl); they carry no diagram embedding assertion.

The sample matrix has one parameter per column and one basis value or spatial
coordinate per row. Control points are columns, in increasing basis-index
order. A scalar parameter is one sample, hence still one output column.
Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.
The cubic basis partition-of-unity samples and quadratic control-point example
in [topology_core.jl](../../test/topology_core.jl) embody this choice.
For a row-point convention the translation is matrix transposition on the
control-point and result matrices; that alternate convention is not an API mode.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| B1 | Degree is a positive integer. | Checked; degree zero is outside the implemented domain. |
| B2 | Every parameter is in `[0,1]`. | Checked; extrapolation is rejected. |
| B3 | Control points have one common row count and at least two columns. | Matrix shape and column count; scalar arithmetic must be supported. |
| B4 | Arithmetic and binomial coefficients do not overflow for any exact-arithmetic claim. | Caller obligation; floating evaluation supplies no error bound. |

question: unverified — A literature key for the Bernstein/Bezier convention
has not been resolved here. The local formulas and fixture values establish
the operation independently; no literature attribution is made. They do not
establish collision freedom or isotopy of a curve drawn by a caller.

## Operations

`bernstein`, `bezier`.
