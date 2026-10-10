# o-graph terminology and notation

Status: normative terminology and notation specification<br>
Scope: o-graphs, their Gauss-code representations, local ports, and condition C2<br>
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
| vertex sign | crossing orientation, unqualified sign |
| local port | germ, dart, half-edge when referring specifically to `±i`, `±o` |
| C2 graph | Graph A, except for source or API correspondence |
| edge weight | vertex weight |
| GC | signed GC |
| weighted GC | weighted signed GC |

## 9. Implementation correspondence

| Documentation term | Current implementation surface |
|---|---|
| component words | `GaussCode`, `RGaussCode`, Gauss-code arrays |
| passage sign | sign of an integer Gauss entry |
| vertex signs | `orientation`, `ROrientation`, `OGraphState.vertex_signs` |
| local-port frame | PD row convention in `VirtualLink.m:36` |
| C2 graph | `getGraphA`, `ograph_graph_a` |
| condition C2 | the second closedness criterion |
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
