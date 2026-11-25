# VirtualLink example routing

Status: normative experiment-input routing

Scope: selecting concrete `VirtualLink` and normal-o-graph inputs for experiments

Implementation references:

- named fixtures: [`examples/VLExample.jl`](../../examples/VLExample.jl)
- parameterized lens spaces: [`src/virtual_link.jl`](../../src/virtual_link.jl)
- integer-framed lens spaces: [`src/framed_lens_space.jl`](../../src/framed_lens_space.jl)

## First entrypoint

When an experiment asks for a concrete manifold, inspect `VLExample.jl` and
this specification before constructing a replacement presentation.  Select the
input by mathematical role, not merely by the manifold name.

For an initial observation of a phenomenon, the owner's default comparison is
several lens spaces, the Poincaré sphere, quaternion space, and `S²×S¹`.
An agent should actively select these inputs when planning a small exploratory
experiment, rather than merely knowing that the constructors exist. An explicit
target or a mathematical precondition takes precedence; explain an omitted
example when it cannot support the requested operation. This is a starting
comparison, not a requirement to run an expensive invariant or a broad sweep.

Load the registry once (the functions are example helpers, not Topology exports):

```julia
using Topology
include(joinpath(pkgdir(Topology), "examples", "VLExample.jl"))
names = (:lens_2_1, :lens_3_1_mst2_example23, :lens_3_2,
         :poincare, :quaternion, :s2xs1_2vertex)
for name in names
    v = vl_example(name) # fresh input for each experiment
    # Check the requested operation's preconditions, then observe it on v.
end
```

| Manifold / role | Retrieval | Real vertices |
|---|---|---:|
| `L(2,1)` | `vl_example(:lens_2_1)` | 2 |
| `L(3,1)` (unoriented name), MST2 Example 23 | `vl_example(:lens_3_1_mst2_example23)` | 2 |
| `L(3,2)` | `vl_example(:lens_3_2)` | 3 |
| Additional lens parameters | `lens_space_link(p,q)` | `p` |
| Poincaré sphere | `vl_example(:poincare)` | 5 |
| Quaternion space | `vl_example(:quaternion)` | 3 |
| `S²×S¹`, mixed vertex signs | `vl_example(:s2xs1_2vertex)` | 2 |
| `S²×S¹`, existing MATLAB presentation | `vl_example(:s2xs1)` | 4 |

These named inputs have no assigned framing. For a framed experiment, use the
framed lens entrance below and check/solve the required weights on other named
presentations; do not silently treat unassigned zero weights as a framing.

The quaternion input is `GC [[[+1,-2,+3,-1,+2,-3]];[+,+,+]]`, the positive
trefoil-shaped **o-graph**. It is used as a spine presentation, not as the
trefoil knot complement. Its reduced edge-group presentation is
`<x,y,z | xy=z, yz=x, zx=y>`, as in
[K10, Fig. 14](../../../../references/items/article/w-0be57e06-c1d1-47c3-bec1-7208afa6fe59/K10.md#^k10-sec-0045).
The diagram identification is agent-checked; this does not select a framing or
assert an orientation-preserving identification with a separately oriented quotient.

**Evidence.**

- literature: [[K10#^k10-sec-0045]] — Fig. 14 and its quaternion group presentation.

The two-vertex `S²×S¹` input is `GC [[[-1,-2,+2,+1]];[-,+]]` in the BP/runtime
convention. Its disk-boundary proof and distinction from the nonclosed
`[+,-]` input are in the
[owning note](../../../../threads/PJ260830/docs/notes/2026-10-09-232224-koda-bp-ograph-conventions.md).
The existing `s2xs1` key retains its four-vertex presentation.

## Routing

### Stored source-specific presentation

Use the finite `VLExample` registry when the request names an existing example,
paper figure, or calculation whose exact presentation matters.

```julia
vl_example_names()
input = vl_example_input(name)
v = vl_example(name)
```

The registry includes distinctions such as the two-vertex MST2 Example 23
presentation of `L(3,1)`, which is not interchangeable with the generic
three-vertex lens-space-family presentation for a source-faithful experiment.
Lens-space names in this registry are unoriented; orientation relations between
the presentations are in [the lens-space naming conventions](../conventions/LensSpaceNames.md) (draft).

### Parameterized unweighted or combed lens space

Use the original family only when the experiment asks for its unweighted or
combing presentation.

```julia
input = lens_space_rgauss(p,q)
v = lens_space_link(p,q)
```

`lens_space_link(p,q)` does not provide an integer framing.

### One integer-framed lens space

When an invariant, connected sum, or other operation requires assigned integer
weights, use the deterministic framed family directly.

```julia
input = framed_lens_space_input(p,q) # raw word, signs, weights, metadata
v = framed_lens_space_link(p,q)      # fresh weighted VirtualLink
```

Do not construct `lens_space_link(p,q)` and then invoke a weight solver.  The
unweighted family can have a non-extendable combing even though the framed
family supplies a valid presentation of the same lens space.

#### Framed-family conditions and stored data

The following conditions name the current public domain; they do not select a
new framing convention. The presentation names and orientation are those of
[LensSpaceNames](../conventions/LensSpaceNames.md), and weight coordinates are
those of [the data model §2.4](virtual-link-data-model.md#24-weight-assignment-and-endpoint-coordinates).

1. **F1 — Scalar parameter domain (checked).** `p` and `q` are non-Boolean
   `Integer` values convertible to machine `Int`; after conversion they satisfy
   `p > q > 0` and `gcd(p,q) == 1`. Boolean and domain violations raise
   `ArgumentError`; an inexact `Int` conversion raises `InexactError`.
2. **F2 — Construction range (caller).** Every intermediate integer operation,
   index and allocation of the construction must be representable. Checking F1
   does not check this range. Fixed-width overflow is not generally detected and
   removes the mathematical construction guarantee. Exhausted memory can prevent
   construction even for valid parameters.
3. **F3 — Assigned-weight storage (caller).** For `framed_lens_space_link`, the
   generated integer weights must be exactly representable by the `Float64`
   storage conversion of `set_weight!`. The raw input record retains `Int`
   weights; conversion to a link does not prove integrality, closed-normality or
   framing independently. Outside this storage range, exact equality with the
   generated integer weights is not guaranteed.
4. **F4 — Series bounds (checked).** `minimum_p` and `maximum_p` are non-Boolean
   integers convertible to `Int` with `2 <= minimum_p <= maximum_p`. Boolean or
   ordering violations raise `ArgumentError`; inexact conversion raises
   `InexactError`. The series contains all F1 pairs for that inclusive range of
   `p`, in increasing `p` then increasing `q`, subject to F2–F3 for every row.

The input record stores `schema = "framed-lens-o-graph-v2"`, the converted `p,q`,
one-component `rgauss::Vector{Vector{Int}}`, its `orientation::Vector{Int}`, and
one component's `weights::Vector{Vector{Int}}`. Both `search_used` and
`linear_solver_used` are false. `orientation_convention` is
`"local_all_positive_lens_space_link"`; this names the reference family and
does not assert that every returned crossing sign is positive.
`edge_weight_convention` is `"edge entering visit; word read cyclically"`.

For `q == p-1`, `construction` is `"previous_q_minus_one_base"` and
`parameters` contains only `vertices=p`. Otherwise `construction` is
`"fixed_insertion_in_L_p_plus_q_p_then_P_r_and_mirror"`; `parameters` contains
`m=p+q`, `d=gcd(2,p)`, `r_modulus=p÷d`, the modularly selected `r`, the patch
coefficient `A`, `B=0`, `patch_vertices`, total `vertices`,
`congruence="2*q*r == q-1 (mod p)"`, and
`raw_seed_orientation="opposite; mirror already applied"`. The construction
and these parameters refer to the
[general construction record](../../../../threads/PJ260829/docs/notes/2026-09-27-034949-deterministic-framed-lens-intake.md#一般構成).
These are construction metadata, not a canonical framing or an equivalence
certificate. Each call allocates fresh word/sign/weight containers. A series
stores `(p,q,link)` rows with independent mutable link containers.

Status: fixed 2026-11-25 Codex (rule a), recording the existing
[constructor](../../src/framed_lens_space.jl) and
[framed-family tests](../../test/lens_space_vg2g.jl). Those tests embody both
branches at `(5,2)` and `(4,3)`, the ordered series through `p=3`, exact weight
storage for `(5,2)`, and selected domain rejections. Field and domain statements
are proved/agent by construction; the general framed-manifold interpretation
retains the [evidence boundary](#evidence-boundary) below. The raw-record and
`Float64` link representations differ only by the declared numeric storage
conversion on this interface; F3 is the condition for an exact translation.
Compatibility beyond this schema's current implementation is not promised.

### Finite parameter sweep

For experiments over many coprime pairs, use the series entrance instead of
maintaining a list of named fixtures.

```julia
rows = framed_lens_space_series(maximum_p; minimum_p=2)
for row in rows
    p, q, v = row.p, row.q, row.link
    # Run the requested experiment on v.
end
```

Each row owns a fresh `VirtualLink`; destructive experimental operations on one
row do not modify another.

## Invariant boundary

The example selector supplies the intended input but does not choose an
invariant algebra.  For example, an exact Taft MST calculation remains an
explicit second step:

```julia
v = framed_lens_space_link(5,2)
z = Z_MST(v,TaftAlgebra(5;M=1,style=:L);strategy=:qa)
```

Do not rerun a broad family verification merely because an experiment consumes
one of these examples.  Verify only the experiment-specific preconditions or a
concrete suspected regression.

## Evidence boundary

The deterministic construction for every coprime `p>q>0` is currently
`proved/agent`, not owner-checked.  The explicit `L(5,2)` witness and recorded
finite checks are `verified-computationally`.  The owning record is
[PJ260829 deterministic framed lens intake](../../../../threads/PJ260829/docs/notes/2026-09-27-034949-deterministic-framed-lens-intake.md).
