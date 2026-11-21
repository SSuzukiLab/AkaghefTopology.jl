# o-graph terminology and notation

Status: normative terminology and notation specification<br>
Scope: o-graphs, their Gauss-code representations, local ports, and conditions N1, N2, C1, C2, C3<br>
Audience: authors, maintainers, reviewers, and AI agents<br>
Owner context: 東京科学大学 情報理工学院 数理・計算科学系 鈴木咲衣研究室<br>
Implementation authority for the local frame: `Manifold/VirtualLink/VirtualLink.m:36`<br>

## 1. Purpose

This document fixes the terms and displays used for o-graphs in this repository.
It distinguishes the represented o-graph from its Gauss-code representation and
keeps the two independent uses of `+` and `−` unambiguous.

Use **o-graph** as the standard spelling. Use **vertex** and **edge** for its
combinatorial constituents. In a formula, `v` denotes a vertex or its positive
integer label; it never abbreviates “visit.” Use **occurrence** or **entry** for
a position in a Gauss word.

## 2. Standard Gauss-code displays

The standard display of an o-graph by Gauss code is

```text
GC [[[+1,+2,-2,-1]];[+,+]]
```

The standard weighted display is

```text
weighted GC [[[+1,+2,-2,-1]];[+,+];[0,1,-1,-1]]
```

In general, the fields are

```text
GC [component words; vertex signs]
weighted GC [component words; vertex signs; edge weights]
```

- **component words** are the cyclic Gauss words of the oriented components;
- **vertex signs** are listed in increasing vertex-label order;
- **edge weights** are listed along the component words, including the cyclic
  edge from the last entry back to the first.

Every positive passage label MUST have an explicit leading `+`, just as every
negative passage label has `-` (or typographic `−`). Thus write
`GC [[[-1,-2,+2,+1]];[-,+]]`. A bare positive passage label is not permitted
in a GC display. This applies to chat, Markdown, TeX, figure labels, HTML,
and human-readable script output, including weighted GC and displayed input
literals. Vertex signs remain `+`/`−`; this rule does not require explicit `+`
on edge weights, unsigned vertex identifiers, or zero.

Format each passage label explicitly; do not use a container's default string
or JSON serialization as its GC display. Language-native numeric storage,
JSON, database keys/hashes, and API arrays retain their existing representation
(JSON numbers cannot carry a leading `+`). Add signs at the display boundary.

The geometric **Gauss diagram** display is separate from this text format. It
places the cyclic word on an oriented circle, draws one arrow from the `+v`
over-passage to the `−v` under-passage for each vertex, and labels that arrow
with the vertex label and vertex sign `σ(v)`. The endpoint signs are passage
signs; `σ(v)` is the independent vertex sign. This display is a presentation
view and does not certify closed-normality or classical realizability.

A `GC` includes its vertex signs by definition. Do not call it a **signed GC**.
Use **weighted GC** only when the edge-weight field is present.

The coloured display is a third field, and that field is not a weight:

```text
coloured GC [[[+1,+2,-2,-1]];[+,+];[0,1,0,-1]]
```

```text
coloured GC [component words; vertex signs; colours]
```

Each colour is one element of Z₃, written `-1`, `0`, or `1`. The list follows
`real_edge_table`: the entry at index `k` belongs to the edge whose head is
word position `k`. Positive passage labels keep the explicit `+`. Use
**coloured GC** for this display. Do not call it a weighted GC. The general
coloured o-graph is not itself a Gauss code; coloured GC is the display of a
closed normal presentation together with those colours.

`GC` and `weighted GC` name representation formats, not the represented
objects. Write, for example, “the o-graph represented by `GC ...`.” Use
**weighted o-graph** for an o-graph whose edges carry weights.

## 3. Input literals

In prose and figures, write vertex signs as `+` and `−`. In an input literal,
encode them by the integers `+1` and `-1`:

```text
[[[+1,+2,-2,-1]];[+1,+1]]
```

A weighted input literal appends the edge-weight field:

```text
[[[+1,+2,-2,-1]];[+1,+1];[0,1,-1,-1]]
```

The display notation and the input literal are distinct surfaces for the same
data. The punctuation above is the canonical documentation form; individual
language APIs may use their language-native container syntax at the program
boundary.

## 4. Passage sign and vertex sign

Two independent signs occur in a GC and must not be conflated.

### 4.1 Passage sign

The sign of a Gauss entry `±v` is its **passage sign**:

- `+v` is the **over-passage** at vertex `v`;
- `−v` is the **under-passage** at vertex `v`;
- `|±v| = v` is the vertex label.

The two entries `+v` and `−v` are the over- and under-occurrences of the same
vertex. Do not call the passage sign a vertex sign, crossing sign, or
orientation.

### 4.2 Vertex sign

The sign attached to the vertex itself is its **vertex sign**. Write it as `+`
or `−` in mathematical documentation and as `+1` or `-1` in program data. When
a symbol is useful, write it as `σ(v)`.

Do not use the unqualified word **sign** when either passage sign or vertex sign
could be intended. Avoid **orientation** as a prose synonym for vertex sign,
because orientation also describes directed components and edges.

Owner terminology decision (2026-11-06): the repository has no global **r/l
type**. The labels `l`/`r` (also `ℓ`, `L`/`R`) belong to individual sources. The
only statement to make is a correspondence for a named Work: whether vertex sign
`+`/`−` corresponds to that Work's `r`/`l`, with a locator and its check status.
A correspondence checked for one Work does not transfer to another, and an
unchecked one stays unresolved. The per-Work table is §4 of
[the source-convention dictionary](../../../../docs/specs/ograph-source-conventions.md).

## 5. Edges and weights

An **edge** is the oriented segment between consecutive entries of a component
word, including its cyclic last-to-first segment. An **edge weight** is attached
to such an edge, not to a vertex.

Use **o-graph edge** when it must be distinguished from an edge of a derived
graph. Use **real edge** only at a VirtualLink implementation boundary where a
real/virtual distinction is active.

## 6. Local ports

The four incidences at a vertex are its **local ports**:

| Port | Meaning |
|---|---|
| `+i` | incoming port of the over-passage |
| `−i` | incoming port of the under-passage |
| `+o` | outgoing port of the over-passage |
| `−o` | outgoing port of the under-passage |

Here `+`/`−` records over/under passage and `i`/`o` records
incoming/outgoing direction. The signs in the port names are therefore passage
signs, not vertex signs.

The authoritative local frame is `Manifold/VirtualLink/VirtualLink.m:36`: a PD
row is ordered counter-clockwise starting from the incoming undercrossing.
Consequently, PD slot 1 is `−i`. With the repository vertex-sign convention,
the full correspondence is

| Vertex sign | PD slots 1, 2, 3, 4, counter-clockwise |
|---|---|
| `+` | `−i`, `+o`, `−o`, `+i` |
| `−` | `−i`, `+i`, `−o`, `+o` |

Program code may represent the four slots numerically, but documentation should
name the corresponding local ports whenever their geometric role matters.

Do not introduce `germ`, “芽,” `oi`, `oo`, or “outside-looking” as synonyms for
these ports.

## 7. Condition C2

The five closed-normal identifiers, their enforcement sites, and the MST
crossing-free circuit are defined in
[§11](#11-closed-normal-conditions-n1-n2-c1-c2-c3). This section names only
the C2 graph.

This specification uses **C2 graph** for the graph constructed to test condition
C2. At each o-graph vertex, apply the prescribed **vertex replacement** and
connect the resulting **local trivalent vertices** along the o-graph edges. The
resulting trivalent graph is the C2 graph.

**Condition C2** holds exactly when the C2 graph is connected. A connected
component of this graph is a **C2 component**.

The legacy implementation names `getGraphA` and `ograph_graph_a` denote the C2
graph. Use **Graph A** only when referring to those APIs or to a source that uses
that name; use **C2 graph** in new documentation.

This document does not introduce a DS diagram or a `4n+k` numbering scheme.
Those belong to neither the terminology nor the explanation of C2 here.

The source figure, repository indices, local triples, and pair permutations
are specified in [`ograph-c2.md`](ograph-c2.md). This terminology file does
not duplicate that numerical construction.

## 8. Preferred vocabulary

| Prefer | Avoid or restrict |
|---|---|
| o-graph | `ograph` in prose |
| vertex | crossing, except at a VirtualLink/PD boundary |
| occurrence, entry | visit; `v` meaning visit |
| passage sign | Gauss-entry sign, over/under `+/-` sign |
| vertex sign | crossing orientation, unqualified sign; r/l type without a named Work |
| local port | germ, dart, half-edge when referring specifically to `±i`, `±o` |
| C2 graph | Graph A, except for source or API correspondence |
| edge weight | vertex weight |
| GC | signed GC |
| weighted GC | weighted signed GC |
| coloured GC | weighted GC, when the third field is a Z₃ colour |

## 9. Implementation correspondence

| Documentation term | Current implementation surface |
|---|---|
| component words | `GaussCode`, `RGaussCode`, Gauss-code arrays |
| passage sign | sign of an integer Gauss entry |
| vertex signs | `orientation`, `ROrientation`, `OGraphState.vertex_signs` |
| local-port frame | PD row convention in `VirtualLink.m:36` |
| C2 graph | `getGraphA`, `ograph_graph_a` |
| condition N1 | `certify_closed_ograph`, `ClosedNormalOGraph`; over pair and quadrivalence in `OGraphState` |
| condition N2 | `OGraphState`; `ClosednessWitness.n2_enforcement` |
| condition C1 | `ograph_state`; `ClosednessWitness.c1_enforcement` |
| condition C2 | the second closedness criterion; Figure 1.2 in [`ograph-c2.md`](ograph-c2.md) |
| condition C3 | the third closedness criterion; `disk_table` row count for a nonempty state |
| MST crossing-free circuit | `mst_crossing_free_circuit`, `MSTCrossingFreeCircuit` |
| geometric Gauss diagram | `gauss_diagram_svg`, `plot_gauss_diagram_svg` |

Implementation identifiers do not override the mathematical terminology fixed
by this document.

## 10. 整数重み非可解

Owner terminology decision (2026-09-09): an o-graph is **整数重み非可解**
when its specified integer weight equations have no integer solution.
The term applies to the o-graph with the stated equation conventions, not to
the underlying manifold. State those conventions when the context is ambiguous.

All edge weights are variables unless the object's definition specifies otherwise.
Failure with fixed exterior weights, a prescribed parity class, or a particular
transport formula is only failure under that restriction; it does not establish
整数重み非可解. A search timeout or failure to find a solution does not establish it either.

整数重み非可解 does not by itself assert absence of a framing or spin structure,
or absence of another integer-weighted o-graph representing the same manifold.
Whether rational solutions exist is a separate property and should be stated
when relevant. Use this term consistently in research records and reports.

## 11. Closed normal conditions N1, N2, C1, C2, C3

Source: [BP97, §1.4], PDF p.15 (printed p.6), `[[BP97#^bp97-sec-0015]]`.
Figures 1.2 and 1.3 are on PDF p.16 (printed p.7), `[[BP97#^bp97-sec-0016]]`.
The sentences below were read from those PDF pages. The OCR extract of the
same pages drops the N1 and N2 sentences, so it is not the source of this
section. The numerical C2 construction stays in
[`ograph-c2.md`](ograph-c2.md).

### Book statements

Read from the text of [BP97, §1.4]. 𝒩 is the set of finite connected
quadrivalent graphs Γ with:

- **N1.** Γ has at least one vertex. Every vertex has a neighbourhood embedded
  in the plane as a normal crossing, up to planar isotopy, and two opposite
  germs are marked over the other two.
- **N2.** Each edge has a direction, and the directions of opposite edges match
  through the vertices.

𝒢 is the subset of 𝒩 whose members also satisfy:

- **C1.** Removing the vertices and joining opposite edges leaves a unique
  oriented circuit.
- **C2.** The trivalent graph of Figure 1.2 is connected.
- **C3.** The oriented circuits of Figure 1.3 number exactly one more than the
  vertices of Γ.

Read from the caption of Figure 1.3: the orientations and the solid dots on
the circuits in that figure matter only for spin and framed normal o-graphs.
The caption does not say that the input edge directions are unused. Those
directions choose which local panel applies. The count in C3 forgets the
arrowheads and the extra dots and keeps the number of circuits. Read from
the caption of Figure 1.2: non-marked vertices are fake. Read from the
figure: Figure 1.2 has one vertex panel and four edge panels; Figure 1.3 has
two oriented vertex panels and four edge panels. The four edge panels of
Figure 1.2 are transcribed in [`ograph-c2.md`](ograph-c2.md). This section
does not transcribe Figure 1.3.

The words "finite connected quadrivalent" are the preamble of 𝒩. They are
outside the labeled N1 sentence. The owner decision of 2026-11-22 places
that preamble under the repository identifier N1, together with the labeled
sentence. Error text for even length and for a label that does not occur
twice says `N1 preamble`, because those checks are the Gauss-code form of
quadrivalence. The error for no vertex, and the error for a missing over
pair, name the labeled N1 sentence.

Read from the paragraph above the axioms on the same page: a whole graph is
drawn as a generic immersion, and two such drawings are identified up to
planar isotopy and the Reidemeister-type moves of Figure 1.1. That global
convention is not the N1 sentence. N1's planar isotopy is the local
embedding of one vertex neighbourhood. A Gauss code stores neither drawing.
This library tests neither isotopy.

### Enforcement

| ID | Locator | Check | Site | Failing input exercised by `test/ograph_state.jl` |
|---|---|---|---|---|
| N1, at least one vertex | [BP97, §1.4] | `true_vertices >= 1`, stored as `n1` | `certify_closed_ograph`, `ClosedNormalOGraph` | `OGraphState(Int[], Int[])` |
| N1, over pair | [BP97, §1.4] | each label once positive and once negative | `OGraphState` | `OGraphState([1, 1], [1])` |
| N1, quadrivalence in one word | [BP97, §1.4] preamble, owner bundle | even length; each label twice | `OGraphState` | `OGraphState([1, -1, 2], [1])`, `OGraphState([1, 2, -2, 2], [1, 1])` |
| N2 | [BP97, §1.4] | opposite directions match | by construction of `OGraphState` (`n2_enforcement`) | none: a direction mismatch cannot be stored in one cyclic word |
| C1 | [BP97, §1.4] | one oriented circuit | `ograph_state` (`c1_enforcement`) | a `VirtualLink` with two Gauss components |
| C2 | [BP97, §1.4], Fig. 1.2 | Figure 1.2 graph connected | `closedness_witness` | `OGraphState([1, 2, -1, -2], [1, -1])` |
| C3 | [BP97, §1.4], Fig. 1.3 | circuit count equals vertices plus one | `closedness_witness` | `OGraphState([1, 2, -1, -2], [1, 1])` |

Connectedness is a preamble condition of 𝒩, separate from C1 in the source.
For a nonempty one-word state the word passes through every vertex, so the
quadrivalent graph is connected. That implication is inferred. The library
does not run a second connectivity search on the quadrivalent graph.
Several Gauss components can still give a connected quadrivalent graph, and
`ograph_state` rejects that input under C1.

N2 holds because one oriented cyclic word directs every edge from one entry
to the next, under the passage-sign rule of `virtual-link-rgauss/1`, and the
constructor rejects any other shape. The two occurrences of a vertex are the
two opposite passages. A head-to-head pair of directions has no encoding in
this type. The witness field `n2_enforcement` names the constructor
`OGraphState`. It is not a Boolean computed from a second orientation. The
constructor does not branch on `convention_id`, so the site string does not
include a convention version.

A repeated passage sign is the labeled N1 over-pair failure. Directions can
still match on that input: the word is cyclic. An earlier coverage note
treated the repeated sign as N2. The PDF assigns the over pair to N1 and the
direction match to N2. This section follows the PDF. A direction mismatch
that keeps the over pair has no Gauss-code encoding here.

C1 for a nonempty `OGraphState` holds because the value contains one cyclic
word. `ograph_state` is the site that rejects several Gauss components. The
Gauss components are the circuits obtained by deleting vertices and joining
opposite edges. Inside a nonempty `OGraphState`, C1 has no failing input.
The empty word is one circuit only as the MST contract. [BP97, §1.4] excludes
it by the labeled N1 sentence, so this section does not use it as a C1
counterexample. For a member of 𝒢 with at least one vertex, C1's count is 1
and C3's count is one more than the number of vertices, hence at least 2.
The two counts are different conditions.

`criteria` is the triple `(C1 slot, C2, C3)`. The C1 slot records the
single-word representation. `c1_oriented_circuits` is the integer derived
from that representation and persisted under verifier version
`virtual-link-c1-c3/2`. On every `OGraphState` that integer is 1. Neither
field is a count of smoothed circuits. The C1 refusal is the `ArgumentError`
raised by `ograph_state`.

C2's component count is `length(connected_components(ograph_graph_a(link)))`
for a nonempty state. The loops and parallel edges dropped by that simple
graph do not change the number of components. The edge list itself is the
construction in [`ograph-c2.md`](ograph-c2.md).

C3's count, for a nonempty state, is the number of rows of `disk_table`
(Method B). The library already identifies that row count with the circuits
of Figure 1.3. This section does not re-derive the identification.

`n1` is the comparison `true_vertices >= 1`.

### MST crossing-free circuit

`mst_crossing_free_circuit()` returns an `MSTCrossingFreeCircuit`: one
oriented circuit and no vertex. The source of this object is the MATLAB/MST
crossing-free generator, by the owner decision of 2026-11-22. It lies outside
𝒢. The fundamental group is infinite cyclic and the first homology is ℤ.
The identification of that state with S²×S¹ is the same owner decision.
Those group claims were not read from [BP97, §1.4].

`certify_closed_ograph` and `ClosedNormalOGraph` refuse the empty code. The
error names N1 and cites [BP97, §1.4].

`closedness_witness` on the empty code returns `n1 == false`,
`criteria == (true, true, true)`, C2 count 1, C3 count 1, and
`method_circuit_lengths == (0,)`. Those C2 and C3 numbers are the MST
extension. `c2_graph` refuses the empty code under N1. `disk_table` returns
a table with zero rows, and the witness does not store that zero: it writes
the MST count 1 and the length `(0,)`. The numbers are not a Figure 1.2 or
Figure 1.3 count. The move verifier accepts a raw target when `all(criteria)`
holds. The inverse PS move has a zero-vertex raw target, and this triple
keeps that target available. Membership of 𝒢 is `certify_closed_ograph`.
The criteria triple does not place the empty code in 𝒢.

`is_closed` on the one-component crossing-free `VirtualLink` is
`(false, (true, true, false))`. Its third flag is the `disk_table` row count,
which is 0, and C3 asks for `vertices + 1`. The MST object is
`mst_crossing_free_circuit`, not the `is_closed` flag.

### Persisted closedness evidence

The verifier version string remains `virtual-link-c1-c3/2`. A stored move-database
row records `closedness_verifier_version` and the object written by
`_closedness_dict`: canonical hash, version, the three criteria bits,
`c1_oriented_circuits`, the C2 count, the C3 count, the vertex count, and the
method-circuit lengths. `verify_ograph_store` rebuilds the presentation from
the Gauss payload (`gauss`, `vertex_signs`, `convention_id`) and requires the
recomputed evidence object to equal the stored object.

N1 is an in-memory gate on `ClosedNormalOGraph`. The fields `n1`,
`n2_enforcement`, and `c1_enforcement` are filled by the eight-argument
constructor of `ClosednessWitness`. `_closedness_dict` does not write them.
Existing nonempty rows still match: the serialized fields and the version
string are unchanged, the C1 slot stays true, and `c1_oriented_circuits`
stays 1. Changing either value, or the version string, would make
`verify_ograph_store` reject those rows. The stored integer is derived from
the single-word representation.

The version string was not given a new value. `certify_closed_ograph` now
refuses the empty code, so the same string sits on a narrower certified
class. Stored nonempty Gauss payloads still certify, and their evidence JSON
still matches. A separate version for 𝒢-membership was considered and not
adopted, because a new string would fail the exact replay of those rows.

The move-witness deserializer reads those same eight fields and calls the
eight-argument constructor, which fills `n1` from the stored vertex count.
A stored nonempty witness therefore reconstructs with `n1 == true`.

### Review status

Pro advice was checked against PDF pp.15–16 and folded above. It is not
evidence.

The question "N1 reading check" is saved as
`threads/PJ261120/sources/pro/normal-20261122-011858.md`. The saved Answer
section of that file is empty; the advice used here is the text in the
sidecar log of the same run. Folded from that advice, and checked on the
PDF: the preamble of 𝒩 is outside the labeled N1 sentence; Figure 1.1 is the
global drawing convention and is not N1; N1's planar isotopy is local;
Figure 1.3's caption limits arrowheads and extra dots, not the input
directions; the empty code is outside 𝒢 by N1 and is not a C1 counterexample;
C1's count and C3's count differ for every member of 𝒢; the over pair can
fail while the cyclic word still directs the edges. The C2 and C3 numerical
examples in the table are library inputs. The book does not discuss those
codes. The identification of `disk_table` with Figure 1.3 is not re-derived
in this section.

One item was not adopted. The same advice says that `all(criteria)` on a
zero-vertex witness cannot be the membership test for 𝒢. This lane agrees,
and `certify_closed_ograph` is that membership test. The move verifier still
accepts a raw target by `all(criteria)`, including the zero-vertex target of
the inverse PS move. That verifier is outside this lane's files, and the
owner decision keeps that raw target working. The criteria triple on the
empty code stays the MST contract.

The question "ID function test check" is saved as
`threads/PJ261120/sources/pro/normal-20261122-013747.md`. That run closed
before an answer was captured, and the Answer section is empty. It was not
sent again. The item-by-item check of these rows against the named functions
and tests is therefore still the lane's own check.
