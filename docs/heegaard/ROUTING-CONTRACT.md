# Heegaard alpha-cut routing contract

The accepted pipeline is:

```text
HeegaardDiagram
  -> cut_along_alpha
  -> fixed deterministic 2 x g puncture placement
  -> ordered port placement
  -> Catlab/Graphviz seed layout
  -> fixed-endpoint fitting
  -> Heegaard-owned feasible topological router
  -> constraint-preserving mechanical smoothing
  -> semantic and geometric verification
  -> rendering
```

Catlab supplies ordered-port seed geometry only. It does not move punctures,
choose the topology, or certify a route. Tutte, FPP, and other free-vertex
planar embeddings are not substitutes for the fixed-endpoint router.

The router is fail-closed. A seed may pass through unchanged only when the
Heegaard verifier proves that it already satisfies every semantic and geometric
condition. Otherwise the implementation must produce a verified reroute or
return an explicit failure; it must not publish an invalid or merely attractive
diagram.

Smoothing starts from a feasible route and remains inside a certified clearance
corridor. It preserves endpoints, puncture placement, boundary cyclic order,
route connectivity, and zero non-endpoint crossings. Visual acceptance remains
an owner review decision rather than a machine-generated claim.

This contract restores the original specification boundary after commit
`c059f61` introduced an unapproved Tutte-based geometric realization.
