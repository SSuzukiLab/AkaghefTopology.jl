# Signed knot Gauss data

## Definitions

### Signed one-component representation

The source objects are a cyclic traversal of one oriented immersed circle,
its crossing visits, and a crossing sign at each double point. The repository
stores the visits as a tuple `gauss` and the independent crossing signs as a
tuple `orientation`. A positive visit denotes the over strand and a negative
visit the under strand; its absolute value is the crossing label. The empty
word represents the crossing-free circle. `KnotGaussDiagram` is the record
subject to K1 below. It does not assert planar realizability.

The precise local crossing frame is owned by
[PD coordinates](knot-classical-invariants.md#pd-coordinates), not by the sign
of a visit. This distinguishes the two independent signs.

### Canonical coordinates

The canonical representative is the least pair `(gauss,orientation)` in tuple
lexicographic order among cyclic shifts, after renumbering crossings by first
occurrence and moving their orientation entries with their labels. This
quotients root choice and crossing labels. Traversal reversal, mirror and
Reidemeister equivalence are separate operations.

Version 1 serialization is the literal string
`KnotGaussDiagram;canon=v1;n=N;gauss=G;orientation=S`, where `N` is the crossing
count and `G,S` are comma-separated canonical tuple entries without brackets.
The canonical hash is lowercase hexadecimal SHA-256 of its UTF-8 code units.
It is a coordinate digest, not a collision-free knot invariant.

### Mirror and traversal reversal

The stored mirror of `(gauss,orientation)` is `(-gauss,-orientation)`, with
the same visit order. Traversal reversal is `(reverse(gauss),orientation)`.
These operations act on representatives and remain outside the canonical
quotient. Their knot-polynomial comparison is recorded in
[the polynomial evidence](knot-classical-invariants.md#coordinate-choices-and-evidence).

### Enumeration and move witnesses

Rooted enumeration fixes first-occurrence crossing labels, varies the pairing
of visits, their over/under choices and crossing signs, and retains a chosen
root. Canonical enumeration identifies the root choices as above. Neither
notion means an enumeration of knot isotopy classes. `KnotGaussEnumeration`
records the requested crossing count and whether all crossing-sign assignments
are selected (`orientations=true`) or only the all-positive assignment.

A rooted unsigned word has positive labels, each twice, introduced in order
of first occurrence `1,...,n`. Its population has size `(2n-1)!!`, with one
empty word for `n=0`. Words are ordered lexicographically by their integer
tuples in the materialized and callback enumerations. Spherical rotations
and their compatible over/under decorations use the ribbon choices below.

A `ReidemeisterCandidate` is a local proposed R1, R2 or R3 replacement together
with its location and parameters in a specific source representation.
It becomes evidence of a transition only when application and replay
verification succeed. A candidate retained after its source changes need
not remain applicable. Bounded classifications are defined
[separately](knot-classification.md).

### Move selectors

For knot witnesses, the family is `R1`, `R2` or `R3`, interpreted without
case. R1/R2 directions are `insert`, `delete` or `both`; omission means both.
R3 uses `self_inverse` or `both`, with omission selecting that same family.
Direction names are also case insensitive. An empty or absent vertex
selector places no restriction; otherwise its absolute labels must equal
the witness's ordered vertex tuple. Edge selectors are ordered pairs of
cyclic visit positions, including `(0,0)` for the empty circle's gap, and
must equal the witness's complete ordered edge tuple. An empty or absent
edge selector places no restriction. The compatibility flag `strict` does
not alter these knot selectors.

The named record uses `param`, `vertices`, `edges`, `strict`, with aliases
`p`, `v`, `e`. A long form takes precedence over its alias. Other fields are
outside this selector model. The status of these repository coordinates is
**fixed 2026-11-25 Codex lane lib-api-rest (rule a)**, following
[reidemeister.jl](../../src/reidemeister.jl); the
[selector tests](../../test/reidemeister.jl) include ordered edge locations,
named aliases and invalid directions.

### Ribbon records

For a word with `n` crossings, the arc from visit position `p` to its cyclic
successor has darts `2p-1` and `2p`. A crossing rotation is a cyclic ordering
of its four incident darts with the two traversal strands opposite.
`rotation` lists these orders; `face_cycles` lists the cycles of crossing
rotation composed with arc reversal. The associated closed oriented ribbon
surface is obtained by capping these boundary cycles. For nonempty words its
Euler data are `V=n`, `E=2n`, `F=length(face_cycles)` and
`g=(2-V+E-F)/2`. The empty circle has the special record `V=E=g=0`, `F=2`,
with empty rotation and face arrays.

`ClassicalGaussCertificate` stores the word, crossing orientations, route,
realizability flag, Euler data, local choices, rotations and face cycles.
For `:carter_surface`, `code` is signed and `rotation_bits` stores the crossing
signs. For `:unsigned_ribbon`, `code` is positive, `orientation` is empty and
`rotation_bits` consists of binary local choices. `realizable` records `g==0`.
`ClassicalGaussImmersion` stores an unsigned word, binary local choices,
rotations, face cycles and genus; its intended spherical use requires K5.
Raw record constructors do not check consistency or sphericality.

The signed local order for positive crossing sign is
`(in_over,in_under,out_over,out_under)`; negative sign reverses that order
while retaining its first dart. For unsigned visits `p,q`, bit zero uses
`(in_p,in_q,out_p,out_q)` and bit one uses
`(in_p,out_q,out_p,in_q)`. These repository choices and reversal invariance of
face count/genus are checked by
[classical_gauss_realizability.jl](../../test/classical_gauss_realizability.jl).
Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**. The source
comparison remains the unresolved literature entry below.

### Enumeration records

`ClassicalGaussEnumerationResult` records a crossing count, numbers of rooted
unsigned words, rotations tested, spherical rotations, generated signed
decorations, duplicate decorations, and the set of distinct packed signed
keys. From most to least significant fields, a key has a four-bit crossing
count, a four-bit token per visit (`2label` for positive and `2label+1` for
negative visits), then one bit per crossing sign (zero for positive, one for
negative). Leading zero bits are implicit in the `UInt64` value. The definition
is local to [classical_gauss_enumeration.jl](../../src/classical_gauss_enumeration.jl)
and its [round-trip tests](../../test/classical_gauss_enumeration.jl).

`KnotGaussCanonicalStreamStats` records the requested crossing count, total
rooted population, visited and emitted counts, supplied-membership true,
false and unknown counts, elapsed seconds, stop flag and reason, sum and
maximum of emitted serialization byte lengths, and a lowercase SHA-256
fingerprint. The fingerprint hashes emitted serializations in emission order,
each followed by a newline; it is not an individual canonical hash.
Membership counts describe the supplied predicate/set rather than an
independent classicality decision. These conventions are exercised by
[knot_n6_all_state_enumeration.jl](../../test/knot_n6_all_state_enumeration.jl).
Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.

## Numbered conditions

| Condition | Meaning | Enforcement boundary |
|---|---|---|
| K1 | For `n=length(orientation)`, the word has length `2n`, contains each of `1,...,n` exactly once with each visit sign, and each orientation is `±1`. | `KnotGaussDiagram` checks; zeros are forbidden. |
| K2 | A classical input is realizable by a planar immersed circle with the supplied crossing data. | Separate realizability certificate/check; K1 alone does not imply K2. |
| K3 | A move witness refers to the supplied source and its full local cyclic neighborhood. | The move/replay APIs check their supported witnesses. |
| K4 | Enumeration crossing count is in the supported finite range `0:6`. | Direct/canonical enumeration entry points check; the raw `KnotGaussEnumeration` constructor stores it unchecked and `length`/iteration check it. This is a resource/domain bound, not a mathematical bound. |
| K5 | An immersion record has a valid unsigned double-occurrence word, matching binary choices and rotations, their face cycles, and genus zero. | `construct_classical_immersion` establishes this; raw record constructors do not validate it. |
| K6 | Certificate Euler data, flag, route and local choices agree with the recorded word, rotations and face cycles. | Certificate-producing functions establish this; raw record constructors do not validate it. |
| K7 | Enumeration or stream counters, key sets and fingerprint refer to one run with the stated crossing count, ordering, stopping and membership parameters. | Raw record constructors do not validate this. |
| K8 | An unsigned word has even length, contains no zero, and after taking absolute values has labels `1:n`, each exactly twice. | Unsigned certificate entry points check this; signs of supplied visits are discarded. |
| K9 | A packed key has the crossing-count prefix and token format defined above for the supplied `n`, and its decoded signed word satisfies K1. | The decoder checks prefix and label range; the diagram constructor additionally checks K1. Classicality is a separate condition. |
| K10 | A supplied callback accepts the documented arguments and returns normally. A canonical-stream callback is a `Function`; its accepted shapes are a record, diagram, `(diagram,planar)` or `(diagram,serialization,planar)`, selected in that order of applicability. | Type/applicability checks where stated; callback failures propagate. |
| K11 | Optional canonical-stream limits are positive `Int`-convertible `max_raw` and positive finite `Float64`-convertible `max_seconds`. | Checked; conversion failures propagate. Time limits apply between raw inputs, not inside callbacks. |
| K12 | At most one canonical membership argument is supplied; it is a set of canonical strings or a function on such strings returning a `Bool`-convertible value. | Exclusivity/type checks are explicit; truth of the membership claim is external. |
| K13 | Move family, direction and selectors use the Move selectors coordinates above; entries must support their integer/string conversions. | Family/direction and pair-shape checks; a well-formed unmatched selector gives an empty population. |

## Coordinate choices and evidence

Status: **fixed 2026-11-25 Codex lane lib-api-rest (rule a)**.
The coordinate, ordering and serialization definitions above are embodied by
[knot_gauss.jl](../../src/knot_gauss.jl),
[reidemeister.jl](../../src/reidemeister.jl),
[knot_gauss tests](../../test/knot_gauss.jl) and
[Reidemeister tests](../../test/reidemeister.jl).
They are repository definitions rather than quotations of a literature theorem.

question: unverified — The source files refer to a literature realizability
criterion, but this page has no resolved registered cite key and locator for
that criterion. The present `classical_gauss_realizable` implementation uses
a Carter-surface certificate; its connection to planar realizability needs
the source/criterion comparison recorded with that certificate. The expected
bridge is genus zero of the oriented ribbon surface; a mismatch of rotation
conventions would invalidate it. The older
`classical_knot_census_available`/`require_classical_realizability` entry points
retain their explicit unsupported behavior and do not describe the newer API.

question: unverified — The Reidemeister source file names Ganzell et al.,
section 2, but this lane has no resolved registered key and locator for that
source. Its internal head/tail names reverse the source arrow terminology,
as the file records. The expected translation is global arrow reversal;
the R3 equal-three-sign predicate should survive it. The
[arrow-reversal and replay tests](../../test/reidemeister.jl) check selected
triples, stale candidates and forged witnesses; the full local-rule/source
comparison is still missing. A discrepancy in the R3 direction factor would
invalidate a general equivalence claim.

## Operations

`KnotGaussDiagram`, `KnotGaussEnumeration`, `crossing_count`, `canonicalize`,
`canonicalize_fast`, `canonicalize_slow`, `canonical_serialization`,
`canonical_hash`, `reindex`, `cyclic_shift`, `orientation_reverse`, `mirror`,
`classical_gauss_realizable`,
`classical_gauss_certificate`, `ReidemeisterCandidate`, `movable`,
`apply_reidemeister`, `verify_reidemeister_transition`.
