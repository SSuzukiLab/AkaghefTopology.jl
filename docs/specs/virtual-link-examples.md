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
| `L(3,1)`, MST2 Example 23 | `vl_example(:lens_3_1_mst2_example23)` | 2 |
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
