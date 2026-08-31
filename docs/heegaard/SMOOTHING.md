# Mechanical relaxation

`Topology.Heegaard.smooth_layout(diagram, layout)` now performs iterative
mechanical relaxation of an already-valid, alpha-cut, fixed-boundary PL layout.
The former one-pass corner fillet remains internal and is not the public result.

Before local relaxation, a multiscale nonlocal shortcut pass removes
low-frequency detours. A shortcut is accepted only when it shortens the route,
its swept region contains no puncture or other route, and the ordinary layout
validator preserves the certified route/intersection structure. The local
objective is total route length with weak bending and obstacle repulsion;
points are redistributed by arc length after each accepted step.
Shortcutting continues until a complete pass accepts nothing. No arbitrary
pass limit is needed because every accepted shortcut strictly removes route
vertices, so this stage terminates finitely.

```text
verified feasible PL routes
-> densify (endpoint guard edges remain exact)
-> overdamped tension + bending + route/self/puncture repulsion
-> clearance-bounded step with backtracking and validate_layout
-> RDP knot reduction with fail-closed validation
-> G1 cubic Bezier chain
-> adaptive flattening and final validate_layout
```

The trace records the energy history, accepted steps, dense/reduced knot counts,
Bezier segment/control counts, minimum clearance, endpoint-tangent check, and
all parameters. The algorithm is deterministic. Fixed port endpoints, fixed
punctures, route identity/successor data, crossing freedom, and puncture
exclusion are hard gates.

Generate the owner-review candidates with:

```text
julia --startup-file=no --pkgimages=no --project=. scripts/generate_poincare_smoothing_candidates.jl
```

The generated directory contains `poincare-baseline.svg` and
`poincare-relaxed.svg`. Web child pages must show feasible input and relaxed
output together at a common scale. Mechanical checks do not decide visual
quality; owner review does.
