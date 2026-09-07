# Source-grounded o-graph move verification

Status: operational candidate specification; MP-D2 and the parameter-3/PS-III
row are source-audited within the explicit gauges and finite fixtures below.

This specification is adjacent to
[`virtual-link-moves.md`](virtual-link-moves.md). That document defines move
generation and raw rewriting; this document defines what may be accepted as a
verified Area-3 transition. Generator success and invariant equality are not
verification.

## Evidence source

The primary source is the managed Work
`w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa`, arXiv `2405.18743v3` (MST25).
Relevant exact locators are:

- `MST25.md` anchor `^mst25-fig-0102`, “Correspondence between D2 and primary
  MP move” (`submit.tex` 608–615);
- §2.4–2.5, Figures 2.13–2.16 (`submit.tex` 477–544), for Method B and PS;
- `submit.tex` 513–516: PS-I `3→2`, PS-II `3→1`, PS-III `1→2`, PS-IV `3→3`;
- `submit.tex` 574 for the strongly-regular/closedness qualification.

The Work text and figure assets are read in place. They are not copied into the
source tree, and their presence is not owner approval of a mathematical claim.

## Verified transition contract

`verify_ograph_move(source, candidate)` does not call `movable`. It must
recompute, in this order:

1. the exact catalogue row, ordered concrete vertices, and source-fragment
   occurrences frozen in `MoveCandidate`;
2. the raw rewrite by that exact candidate;
3. C1–C3 separately on source and raw target;
4. every occurrence outside the replaced fragments, the old-to-new vertex map,
   and every consecutive outside adjacency;
5. raw endpoint hashes and the exact root/relabel maps used to canonicalize
   both endpoints.

The canonical target is an endpoint identity only. Canonicalization is not part
of the move and must not erase the raw target, application locations, or
outside-region evidence.

The accepted shared-library scope is currently `MP`, parameter `D2`, and `PS`,
parameter `3`. Other MP catalogue rows fail with
`mp_parameter_not_source_audited`; PS rows 1, 2, and 4 fail with
`ps_parameter_not_source_audited`. These are evidence boundaries, not
statements that the moves are illegal.

### D2 regression fixture

The nontrivial inverse-D2 fixture is:

```text
source GC: [-3,1,3,2,-2,-1,4,5,6,-6,-5,-4]; vertex signs: ++++++
vertices:  [1,2,3]
locations: [(1,1),(1,3),(1,5)]
raw target GC: [2,-1,1,-2,3,4,5,-5,-4,-3]; vertex signs: +++++
```

The removed vertex is old `3`; old vertices `4,5,6` map to new `3,4,5`.
The five outside adjacencies
`[4,5]`, `[5,6]`, `[6,-6]`, `[-6,-5]`, `[-5,-4]` map respectively to
`[3,4]`, `[4,5]`, `[5,-5]`, `[-5,-4]`, `[-4,-3]`.

## PS gate

MST25 Figure 2.13 replaces each true-vertex neighborhood by six local curve
arcs. Figure 2.14 numbers the three lines around an upward-viewed o-graph edge
from left to right, with circuit orientations `↑,↑,↓`. Therefore a replayable
PS witness must retain, for every Method-B segment,

```text
(real-edge occurrence, lane 1..3, directed local-node endpoints, circuit)
```

and identify the two selected source fragments as the left and right edges in
the Figure 2.15 frame. Merely being on the same aggregated disk row is
insufficient.

The managed source figures and the implementation overlay are shown together
so that this index choice can be checked without reading the arrays as magic
numbers:

![MST25 Figure 2.13, Method B](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/assets/MST25-method-b.png)

![MST25 Figure 2.14, lane convention](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/assets/MST25-correspondence-between-configuration-of-2-cells-around-triple-line-and-lines-in-circuit-diagram.png)

![Indexed positive Method-B frame](assets/ograph-method-b-positive-overlay.svg)

The Figure 2.13 asset used for the independent trace has SHA-256
`741A4FF1C9D0DFD2189560F81F669FC4C4314270BA3AB6F6D06EDC3DB822F9F3`.
In its physical frame, the top stubs `U1,U2,U3` and bottom stubs `D1,D2,D3`
run west-to-east, while `L1,L2,L3` and `R1,R2,R3` run top-to-bottom. The six
curves are
`U1-D1`, `U2-L1`, `U3-R1`, `D2-L3`, `D3-R3`, and `L2-R2`;
unmarked crossings of curves are not nodes. With under direction `D→U` and
over direction `L→R`, the port identification is
`D=-i`, `U=-o`, `L=+i`, `R=+o`.

Fix the twelve stubs as four consecutive blocks
`(-i,+i,-o,+o)`, each with Method-B lane index `1,2,3`. The absolute local-node
numbers below are an implementation gauge. Tracing the six curves of Figure
2.13 in the positive-vertex port frame gives:

| local node | two incident `(port,lane)` stubs | array positions |
|---:|---|---|
| 1 | `(-i,1)`, `(-o,1)` | `1,7` |
| 2 | `(-i,2)`, `(+i,3)` | `2,6` |
| 3 | `(+i,1)`, `(-o,2)` | `4,8` |
| 4 | `(+i,2)`, `(+o,2)` | `5,11` |
| 5 | `(-o,3)`, `(+o,1)` | `9,10` |
| 6 | `(-i,3)`, `(+o,3)` | `3,12` |

Thus the positive table is `dp` below. In the repository's normalized negative
vertex frame, the under-passage ports are unchanged while the over-passage
direction reverses. Therefore exchange `+i` and `+o` and reverse its viewed
lane order by `k ↦ 4-k`. Applying that transformation to the six positive
pairs gives the negative pairing positions
`(1,7),(2,10),(3,4),(5,11),(6,9),(8,12)`, hence `dm`:

```text
dp = [1,2,6,3,4,2,1,3,5,5,4,6]
dm = [1,2,3,3,4,5,1,6,5,2,4,6]
```

where the four three-stub blocks are currently consumed as
`(-i,+i,-o,+o)`. The lane-3 Method-B segment is stored in the reverse circuit
direction, matching Figure 2.14's `↑,↑,↓`; this changes directed trace order but
not circuit membership.

The current template row 3 has fragments `[1,2]`, `[-2,-1]` and vertex signs
`(+,-)`. Overlaying their directions and the two signed vertex frames on Figure
2.15 identifies the row with PS-III; Figure 2.16 requires left lane 1 to connect
to right lane 2. Rows 1, 2, and 4 remain unaccepted pending equally explicit
left/right source-edge overlays and positive circuit fixtures.

![MST25 Figure 2.15, PS-I through PS-IV](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/assets/MST25-pure-sliding-move-ps-i-ps-iv.png)

![MST25 Figure 2.16, corresponding lane conditions](../../../../references/items/article/w-49d1ab60-89f9-43f8-99f6-839e1ae5b0aa/assets/MST25-conditions-ps-i-ps-iv-for-ps-i-ps-iv-respectively.png)

For PS-III the left edge in Figure 2.15 points down, so its displayed strands
appear `3,2,1` from screen-left to screen-right in Figure 2.16; the intrinsic
upward-viewed lane remains lane 1. The right edge points up and is displayed
`1,2,3`. Thus the dashed connection is intrinsically left lane 1 to right lane
2, exactly as stated in `submit.tex` 510–516.

There is an important direction reversal. `make_move_data()` stores only the
inverse PS direction: its two-vertex RHS is the verifier input and the
zero-vertex LHS is the raw target. MST25 states the global condition on the
diagram *before the forward move*, i.e. on that LHS. The verifier therefore
does **not** test Method B on the inverse source. For each of the two removed
RHS fragments it records the exact source positions and walks past the complete
consumed patch to the nearest outside predecessor and successor. Those signed
occurrences are transported through the old-to-new vertex map and must form an
actual consecutive edge occurrence in the raw target. Method B and the lane
condition are evaluated on that raw target; the witness fixes
`evaluated_state_side=target_raw` and its raw hash.

For the stored PS-III fixture
`[-4,-3,-2,-1,1,2,3,4]` with vertex signs `(+,-,-,-)`, the concrete local
locations are `(1,5)` and `(1,3)`. The two distinct fragment incidences occupy
positions `[5,6]` and `[3,4]`. Both gap walks reach source boundary
`-3@2 → 3@7`, which transports to the same raw-target edge occurrence
`-1@2 → 1@3` (edge 3). Keeping the fragment index and source positions is
essential: the two region incidences must not be deduplicated merely because
their transported target edge ID agrees. The occurrence-preserving Method-B
trace places `(edge 3,lane 1)` and `(edge 3,lane 2)` on the same circuit. The
verifier records the complete directed circuit, not only this Boolean result.

This explicit witness requirement does not assert that closedness is
mathematically insufficient in the strongly regular setting.

Before rows 1, 2, or 4 are enabled, the specification and tests must add:

- the geometry-derived mapping from implementation parameters `1..4` to paper
  types I, II, and IV, not merely their common row order;
- positive circuit traces for those rows and lane/occurrence/type mutations.

## Persistence and replay

`ograph_move_witness_data` is the canonical JSON-like surface.
`ograph_move_witness_from_data` requires the exact schema and field set and
rejects type drift. `replay_ograph_move` reconstructs the source and catalogue
candidate, reruns the semantic verifier, and compares every regenerated witness
field. A DB must call this replay gate; a stored `valid` flag, matching
invariants, or the mere presence of a witness object is not acceptance.
