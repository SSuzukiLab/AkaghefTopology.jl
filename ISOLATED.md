# Isolated: marked-surface work (cellulation, seam-word marking)

This branch snapshots the `codex/heegaard-full-curves` working tree of 2026-09-27,
including the supplied marked-cellulation kernel (`vendor/marked_surface_diagrams`),
the genus-3 Koda bundle in cellulated form (`vendor/koda_genus3`), the cellulated
app view, and the 2026-09-24 seam-word marking.

Reason for isolation: these encode curves up to isotopy on a *marked* surface.
The intended equivalence for 2D alpha-cut Heegaard pictures is coarser: labelled
diagrams up to orientation-preserving homeomorphism, i.e. ambient isotopy of the
picture in S^2. For filling diagrams the signed intersection orders already
determine that class (plus the choice of the face at infinity), so the marking
layers are not needed on the main line. Marked-surface curve data on fixed
triangulations/square complexes is also covered by existing tools (flipper,
curver, SnapPy's twister).

`scripts/heegaard/migrate_koda_genus3.jl` is the one-shot conversion of the
genus-3 example to the ordinary punctured input used on the main line.
