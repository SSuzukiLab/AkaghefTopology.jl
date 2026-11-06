# Evaluating Lescop's Θ from a Heegaard diagram

<!-- XQQQ1500000 -->
Use Topology.Heegaard.lescop_theta(diagram, matching, exterior_region).
The diagram is a certified HeegaardDiagram with two complete cut systems.
Its surface has the boundary orientation of the alpha handlebody. Give one
crossing ID on each alpha and each beta curve. The exterior region is the
index of the complementary component containing w, in the original labelled
diagram. The output is a named tuple of exact rationals, (ell2, sl, e, theta).
Use lescop_intersection_terms(diagram, matching) for the first two terms alone.

The definitions follow [Lescop §§1.3–1.5](https://arxiv.org/abs/1402.2261).
Each common endpoint has a product of two half weights. For the rotation
term the evaluator writes the planar face turning equations after cutting
along alpha and placing the exterior point at infinity. It constructs a
dual face potential whose differences give the coefficient of each
oriented beta edge turn in e. This eliminates the individual arc turns.
The code checks the potential and the face equations exactly and rejects
inputs where this certificate fails. The output refers to Lescop's infinity
combing for the supplied exterior point and matching.

## Conversion and retained data

- to_odata(diagram; joining_disk=D) uses the supplied labelled JoiningDisk.
  The punctured result and the branching data retain the disk boundary,
  source crossing identities and the orientation of each source passage.
- cellulated_punctured(document) checks complete cut systems and the oriented
  joining region in the existing marked document. Its returned object is
  a PuncturedHeegaardDiagram. The original document, region correspondence
  and auxiliary-curve information are retained as provenance. The curve
  data for the calculation is the object's diagram. Select the exterior
  region and matching in that diagram using the correspondence.
- to_heegaard(odata) constructs a new marked splitting from the dual
  tetrahedra. Its conversion witness describes the new cellulation. An
  original joining disk, matching, genus, marking and combing require
  additional correspondence data. The retained-source route explicitly
  retrieves an original attached source when the input carries one.

The implementation migrated the Heegaard and cellulated Julia sources from
the full-curves worktree at deb1bae, with the plotting dependencies required
by the public method. Catlab uses the source worktree's vendored
RuntimeGeneratedFunctions 0.5.25 compatibility fix for Julia 1.12.
Package tests in test/heegaard/lescop_theta.jl use Lescop's two worked
decorated diagrams and a bounded family of signed ribbon fixtures.
<!-- XQQQ1000000 -->
