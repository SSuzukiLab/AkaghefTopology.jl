# Evaluating Lescop's Θ from a Heegaard diagram

<!-- XQQQ1500000 -->
Use Topology.Heegaard.lescop_theta(diagram, matching, exterior_region).
The diagram is a certified HeegaardDiagram with two complete cut systems.
Its surface has the boundary orientation of the alpha handlebody. Give one
crossing ID on each alpha and each beta curve. The exterior region is the
index of the complementary component containing w, in the original labelled
diagram. The output is a named tuple of exact rationals, (ell2, sl, e, theta).
Use lescop_intersection_terms(diagram, matching) for the first two terms alone.
Use Topology.Heegaard.hc_lescop_theta(hc, matching_labels, exterior_face)
for the same tuple on a resolved HeegaardCode. The labels are the LC
labels of one α–β crossing on each curve. The exterior is an LC face id
in the complement region of w. Roleless K is deleted and those LC faces
become one complement region. A K curve with an α or β role stays. The
face number is mapped to the region and is not passed through. The
intersection sign is Lescop's σ, the opposite of the passage sign on α.

The definitions follow [Lescop §§1.3–1.5](https://arxiv.org/abs/1402.2261).
Each common endpoint has a product of two half weights. For the rotation
term the evaluator writes the planar face turning equations after cutting
along alpha and placing the exterior point at infinity. It constructs a
dual face potential whose differences give the coefficient of each
oriented beta edge turn in e. This eliminates the individual arc turns.
The code checks the potential and the face equations exactly and rejects
inputs where this certificate fails. The output refers to Lescop's infinity
combing for the supplied exterior point and matching.

## The p₁ difference across one beta edge

Use Topology.Heegaard.lescop_p1_difference(diagram, matching, from_id, to_id).
The two ids are successive crossings of one oriented beta curve. The path
meets that edge once, with algebraic intersection +1, and meets no alpha
curve, so it runs from the geometric left of the edge to the geometric
right. The function returns the exact pair (lk, delta_p1). Here lk is the
linking number of Lescop's Proposition 5.1, summed over every crossing,
and delta_p1 is eight times that number. Theorem 3.4 identifies delta_p1
with p₁ on the right of the edge minus p₁ on the left. The half-bracket
gives both endpoints weight 1/2, an interior crossing weight 1, and the
arc from a matching crossing to itself weight 0.

The same pair on a resolved HeegaardCode is
Topology.Heegaard.hc_lescop_p1_difference(hc, matching_labels, from_label, to_label).
The labels are the LC labels of alpha-beta crossings. Roleless K is deleted
before the edge is read.

This difference does not consult Θ. An additive constant remains for each
diagram and matching: the variation formulae do not fix p₁(X(m, R_w)).
The +1-surgery expressions (4g(F)-1)g(F) and 4g(F)(g(F)-1) are changes for
the normal field of a Seifert surface, and this function does not use them.
Corollary 5.3 says that four times the change in Θ across the same edge
equals delta_p1. That comparison is a separate check.

## Spin^c labels

Use Topology.Heegaard.lescop_epsilon(diagram, matching, other) for the
class ε of [[OS04a]] subsection subsec:Disks. The cycle runs along alpha
from the first matching to the second and returns along beta. Coordinates
are the intersections with the alpha curves, reduced modulo the columns of
the intersection matrix, and stored as the lexicographically smallest
vector with entries in 0:|det|−1. An arc from a matching point to itself
contributes 0. The same pair on a resolved HeegaardCode is
Topology.Heegaard.hc_lescop_epsilon(hc, matching_labels, other_labels).

Use Topology.Heegaard.lescop_spinc_label(diagram, matching, region) for the
label of s_region(matching) in that coordinate group. The reference is the
lexicographically smallest matching and complement region 1. At one region,
the difference of two labels is ε, and [[OS04a]] lemma lemma:VarySpinC
identifies the Spin^c difference with the Poincaré dual of ε. An alpha edge
changes the label by −e_i from its odd-dart side to its even-dart side.
The beta class is the unique coboundary compatible with those alpha
classes; the paper states the basepoint formula for alpha arcs. Θ and p₁
are not read. The HeegaardCode form is
Topology.Heegaard.hc_lescop_spinc_label(hc, matching_labels, exterior_face).
On a rational homology sphere the number of labels is |H₁|. The first
Chern class in the same subsection is c₁(s) = s − conjugate(s). A Koda
Euler representative is not identified with that class here.

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
