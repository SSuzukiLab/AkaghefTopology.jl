# Disk IDs as ordered tetrahedral edges

`disk_tetrahedra(v)` exposes the existing dual-tetrahedron recipe for a nonempty
closed normal o-graph. It first uses `certify_closed_ograph`; this is the current
C1–C3 implementation gate, not an additional manifold-recognition theorem.
The crossing-free input is rejected explicitly.

```julia
data = disk_tetrahedra(v)
data.tetrahedra[1].edges  # disk IDs for (01,02,03,12,13,23)
data.tetrahedra[1].sign   # tetrahedron orientation; minus vertex_sign
data.node_to_disk        # production local-node ID -> disk ID
```

The implementation derives its six-germ maps from `_disk_node_dictionary()`
and the reviewed face-to-ordered-vertex map already used by the
[CP recipe audit](ograph-move-verifier.md). It does not copy a separate disk
table or select a reading by matching invariant values.

| Repository vertex sign | Germs for (01,02,03,12,13,23) | Tetrahedron sign |
|---|---|---|
| +1 | (1,5,6,3,2,4) | -1 |
| -1 | (1,3,5,2,6,4) | +1 |

For every actual o-graph edge the converter checks the complete directed face
map, not just unordered incidence: ascending ordered-edge endpoints must map
to the ascending target endpoints, and the paired germs must have the same
disk ID. Every generator of the disk equivalence relation therefore preserves
direction. Its transitive closure gives each disk ID a consistent global edge
orientation. All returned edge IDs are positive in this explicitly checked
orientation gauge. A future different face convention must pass the same check
or fail; it must not silently turn an inverse edge into a direct label.

Paired faces must also carry opposite induced orientations after multiplication
by the tetrahedron signs. This is checked on every gluing. The output flags
`edge_orientations_verified` and `face_orientations_verified` report those
checks, while `closedness` carries the existing operational C1–C3 witness.

For a finite group coloring, every ascending face uses `g_ij*g_jk=g_ik`, in
that order even for a negative repository vertex. The converter returns
combinatorial data only: it does not choose a cocycle, normalization, trace,
state-sum value or MST equivalence theorem.
