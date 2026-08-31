# Mechanical smoothing

`Topology.Heegaard.smooth_layout(diagram, layout)` is a bounded post-process
for an already-valid, alpha-cut, fixed-boundary PL layout. It does not create
or optimize a layout. The input Poincare route certificate remains the
feasibility witness.

Each interior PL vertex `V` is replaced by a cubic segment with controls
`(A,V,V,B)`, where `A` and `B` are stepped along the two incident PL edges.
The step is bounded by the local edge-length bound, one quarter of the
distance to every other route, and half of the nonnegative puncture margin.
Thus every generated corner lies in the corresponding local tube around the
original feasible route. Ports and puncture centers are copied unchanged.

`Route.points` is the deterministic flattened cache of `Route.curve`. A
`BezierChain` stores its flattening tolerance and recursion depth. Validation,
JSON serialization, SVG `C` commands, and TikZ output read the curve field;
the validator rejects a stale cache, malformed curve, non-exact endpoint, or
any flattened crossing, self-crossing, or puncture intrusion.

Generate the owner-review candidates with:

```text
julia --startup-file=no --pkgimages=no --project=. scripts/generate_poincare_smoothing_candidates.jl
```

The generated `artifacts/heegaard/poincare_smoothing_candidates/summary.json`
records the three step fractions and their flattened-geometry metrics. The
candidate SVGs are intended for visual review and do not replace the
validator’s finite flattened checks.
