# VirtualLink multi-representation data-model rationale

Status: §2.1–§2.6 specify current Julia behaviour; the other sections retain the design rationale and provisional target<br>
Scope: mathematical state, representation conversion, caching, and mutation semantics of MATLAB and Julia `VirtualLink`<br>
Audience: maintainers, port authors, reviewers, and computation authors<br>
Owner context: 東京科学大学 情報理工学院 数理・計算科学系 鈴木咲衣研究室<br>
Source inspected: `Manifold/VirtualLink/VirtualLink.m` and the current Julia migration worktree, 2026-08-26

## 1. Thesis

A `VirtualLink` should not be forced into one universal stored representation.
Gauss, real Gauss, PD, table, o-data, and text inputs are useful at different
boundaries and may carry source-specific labels or provenance. The design should
therefore preserve the received input, prefer real Gauss data for mathematical
computation, and derive other representations only when required.

The core policy is:

> Preserve the supplied presentation as input authority; use real Gauss plus
> crossing orientation as the preferred computational view; treat full Gauss,
> PD, tables, layout, and positions as reproducible views or caches.

This separates two meanings of “source of truth”:

- **input authority** — the exact representation and payload supplied;
- **calculation authority** — the representation used as the common interface by
  algorithms.

Real Gauss is the preferred computational view. It is not automatically the
original input and must not erase another supplied representation.

## 2. Roles of stored data

| Role | Meaning | Examples |
|---|---|---|
| Ingress source record | Exact received representation and payload | Gauss, RGauss, PD, o-data, table identity, text input |
| Preferred computational view | Common representation requested by moves, weights, and invariants | RGauss + ROrientation |
| Derived combinatorial view | Recomputable conversion or table | Gauss, PD, RPD, crossing/edge/strand tables |
| Derived geometric view | Recomputable drawing state | planarisation, bending numbers, positions, layout, SVG |
| Evidence/replay data | Captured external output used for comparison or replay | MATLAB `REdgeTable.Position`, exported JSON, figures |

Weights are mathematical annotations on real arcs. They may be part of a
weighted o-graph computation without becoming the topological identity of the
underlying unweighted virtual link.

### 2.1 Open strands and cut markers

The Julia `VirtualLink` representation records an open strand with one trailing
`0` in its raw Gauss component. The marker is metadata, not a crossing: a zero
is valid only as the final entry of a component, and `strand_cuts(link)` derives
one actual `Bool` per strand from that marker. `gauss_code` and
`real_gauss_code` preserve the marker; virtual crossing visits remain in the
raw code and are omitted from the real code, while the trailing marker remains.

For a cut strand with (m) real crossing visits, `real_edge_table` exposes
(m+1) directed arcs, including the distinct boundary edges `(0, first)` and
`(last, 0)`. A strand with no real visits still has the crossing-free boundary
arc `(0, 0)`. These boundary IDs are representation-local and must not be
identified with closed PD edge labels.

`set_cuts!(link, flags)` requires one actual `Bool` per stored strand. On an
unweighted link it adds or removes only the trailing markers and resizes the
derived zero weight vectors by prefix. Changing the cut status of a weighted
link is rejected; callers must clear or rebuild the weights and then assign the
explicit (m+1) open-arc vectors. This prevents a closure weight from being
silently copied to both boundary edges.

`pd_code` is a full/raw traversal view: its edge IDs include virtual visits and
open boundary IDs, so they are not the IDs of `real_edge_table`. PD decoding
uses the supplied crossing orientation to recover directed strands; an open
boundary is preserved as a trailing marker when the conversion is defined.

A mathematical virtual link is not identified by a unique Gauss or PD normal
form in the current implementation. “Source” therefore means an adopted
presentation, not a canonical representative of an equivalence class.

### 2.2 Current Julia storage and coherence

The current [`VirtualLink`](../../src/virtual_link.jl) stores mutable Gauss
components, crossing orientations, `Float64` weight vectors, `real_only` and
`is_weighted` flags, strand-type characters, and a weight-mode symbol. The
positional constructor converts these fields but does not validate their mutual
consistency; correctly typed mutable arrays can remain shared with the caller.
The keyword empty constructor validates the mode (`:integer` or `:z2`).

A coherent crossing representation uses labels `1:length(orientation)`, with
one positive and one negative occurrence of each crossing. Orientations are
`-1` or `1` for real crossings and `0` for explicit virtual crossings; a
real-only representation has no virtual crossings. Cut markers obey §2.1.
There is one strand type and one weight vector per component; the weight-vector
length is the real-arc count of §2.1 (the real-visit count for a closed component).
An empty closed component has zero stored weights, an empty cut component one.
These are caller coherence requirements, not a claim that every constructor
or reader enforces them. In particular, the Gauss ingress validator counts
absolute labels twice but does not require opposite signs, and its empty-label
branch accepts even an unused orientation vector. Direct field mutation can
invalidate any of these relationships.

`head_map` has one row per orientation entry and two columns: column 1 is the
successor of the negative occurrence and column 2 the successor of the positive
occurrence. Successors are signed entries of the raw Gauss code, including
virtual crossings. Closed components wrap; cut endpoints have successor `0`.
Unfilled cells remain zero; repeated occurrences of the same sign overwrite
the same cell. This array is a presentation coordinate, not a validity test.

### 2.3 Ingress precedence

The current `set_data!` routing is fixed here rather than inferred from the
order of keyword arguments. `odata` first supplies real Gauss data, orientation
and, for a three-entry payload, weights; otherwise `table=(3,1)` supplies the
supported trefoil PD data and orientation. Other table values are rejected.
The resulting non-`nothing` `pd` then takes precedence over Gauss data. Without
PD, `gauss` supplies the word when present, otherwise `rgauss` does; the
`real_only` flag records whether `rgauss` is non-`nothing`. Thus supplying both
Gauss keywords does not mean that the real word wins. `odata` suppresses the
table branch, but does not suppress an explicit PD payload.

After geometry assignment the assigned-weight flag is cleared, supplied cuts
are applied, strand types are assigned (default `'o'`), and weights from a
three-entry o-data payload are assigned or zero arc-weight vectors are built.
The operation is not transactional: a later error can leave earlier field
changes in place. Use one ingress representation when a mixed-input precedence
is not intended. The implementation does not preserve a separate immutable
ingress record of the kind proposed in §2's design table.

### 2.4 Weight assignment and endpoint coordinates

`set_weight!` accepts a flat numeric vector, a numeric matrix flattened in
Julia column-major order, or one iterable per component. Flat values are split
in component order using the real-arc counts. Shape mismatches are rejected.
Integer mode stores converted `Float64` values; this assignment does not itself
check integrality or finiteness and need not preserve large integers exactly.
Z2 mode requires finite integral real inputs and reduces them modulo two before
storage as `0.0` or `1.0`.

Changing an unweighted link's mode changes its tag. Changing assigned integer
weights to Z2 reduces their values. Changing assigned Z2 weights to integer
mode is rejected because no lift is selected. Setting the already selected
mode is a no-op, not a revalidation. Explicit `set_weight!` can assign new
integer weights with `mode=:integer`.

Real-edge IDs run through components and their arcs in order. On a closed
component the weight at visit `i` belongs to the arc from its cyclic predecessor
to that visit. On a cut component the order starts with the incoming boundary
arc, continues between visits, and ends with the outgoing boundary arc.
`set_weight_by_crossing!` uses ordered signed endpoint pairs in these coordinates.
Repeated assignment keys retain the last supplied value; every edge with a
matching pair receives it, and unmatched keys change no edge. Consequently a
pair is not an edge-occurrence identifier when several edges have that pair.
Other arc weights are retained. Pair assignment is followed by ordinary
weight assignment and its checks.

### 2.5 Evidence and compatibility of the current conventions

Status: fixed 2026-11-25 Codex (rule a). These entries record the existing
implementation's choices rather than the future storage redesign. The
comparison cases are in [core tests](../../test/topology_core.jl),
[cut-weight tests](../../test/cut_weights.jl),
[closed-boundary tests](../../test/cut_closed_boundaries.jl), and
[spin-weight tests](../../test/spin_weights.jl). They exercise closed versus cut
arcs, endpoint assignment, Z2 reduction and rejection of a canonical integer
lift. The translation from endpoint-pair coordinates to occurrence coordinates
is the matching row set of `real_edge_table`; it is one-to-many when pairs repeat.
Raw Gauss head coordinates retain virtual visits, whereas real-edge coordinates
omit them; these coordinate systems must not be substituted for each other.

The branch behaviour above is supported by inspection of `set_data!`,
`head_map`, `set_weight!` and `set_weight_by_crossing!` [proved/agent: by
construction]. A comprehensive regression check of mixed ingress, duplicate
endpoint keys and deliberately incoherent raw constructors remains unverified:
the existing tests do not establish that coverage. Expected behaviour is the
branch behaviour recorded above; a test reaching a different branch/output
would require correction of this entry. No compatibility promise for a future
immutable redesign follows from this snapshot.

### 2.6 Circuit endpoint selection

These conditions specify `circuit_b_at`'s selection of a directed Method-B
path endpoint. They use the [cut markers](#21-open-strands-and-cut-markers),
[Method-B path model](ograph-move-verifier.md#ps-gate), and
[boundary lanes](spin-weights.md#relative-boundary-data).

1. **E1 — Coherent circuit data (caller).** The presentation has coherent
   real-edge, cut-marker and Method-B path data as required by `disk_table`.
   Marker and incidence checks can reject malformed data; geometric validity
   is not certified by selecting an endpoint.
2. **E2 — Strand and side (checked).** A non-Boolean `strand::Integer` indexes a component of
   `v.gauss` and that component is cut. Its `side::Symbol` is `:in` for the
   component's incoming cut endpoint or `:out` for its outgoing cut endpoint.
   An out-of-range strand raises `BoundsError`; an unsupported side or uncut
   strand raises `ArgumentError`. `false` fails the range check with
   `BoundsError`; `true` fails with `BoundsError` when there is no component,
   and otherwise reaches Julia's scalar-array-index rejection (`ArgumentError`).
3. **E3 — Lane selection (checked).** `lane=nothing` selects lanes `1,2,3`, in
   that order. Otherwise `lane` is an `Integer` in `1:3` and selects one lane.
   `true` therefore selects lane 1; `false` and noninteger values are rejected.
   An invalid selector raises `ArgumentError` naming the accepted lanes.
4. **E4 — Unique directed incidence (checked).** Each selected endpoint
   `(strand=Int(strand), side=side, lane=Int(lane))` occurs at the start or stop
   of exactly one open Method-B path. Zero or multiple matches raise
   `ErrorException` stating that a boundary lane must belong to exactly one
   directed B path. A successful query returns one record per selected lane.

The returned record extends the corresponding `circuit_b_boundary` record by
`endpoint`, `direction`, and `other`. `direction=:outgoing` means the selected
endpoint is the path's `start`; otherwise it is `:incoming` and the endpoint is
its `stop`. `other` is the opposite end of that same directed path. Strand side
and path direction are distinct: at the incoming endpoint of a crossing-free
interval, lanes 1 and 2 are outgoing paths and lane 3 is incoming. No strand is
cut by this query, and the input presentation is unchanged.

Status: fixed 2026-11-25 Codex (rule a), recording
[`circuit_b_at`](../../src/circuit_b.jl) and
[endpoint-selection tests](../../test/circuit_b.jl). The latter specify the
crossing-free interval's lane order, directions, signed edges `[1,1,-1]`,
opposite endpoint, single-lane query, and invalid side/lane/uncut selections.
The selection and error statements above are proved/agent by construction.
The Boolean lane acceptance and Boolean strand rejection follow the existing
integer predicate and Julia indexing, respectively; they are not a claim of
dedicated library-test coverage. There is no translation from these endpoint
coordinates to a boundary trivialization supplied by this query; the relative
framing/spin interpretation retains the boundary model's evidence limits.

## 3. Current MATLAB model

### 3.1 Broad ingress

MATLAB `VirtualLink.setData` accepts several forms. Native Gauss, RGauss, and PD
inputs are retained in their corresponding fields; adapter inputs such as
`odata`, `table`, and `ostring` are first routed to a supported working form.

| Input argument | Initial stored form | Initial route/cache state |
|---|---|---|
| `GaussCode` | `GaussCode`, `orientation` | `G=true`, `virtualFlag=false` |
| `RGaussCode` | `RGaussCode`, `ROrientation` | `VG=true`, `virtualFlag=true` |
| `PDCode` | `PDCode`, `orientation` | `PD=true`, `virtualFlag=false` |
| `odata` | routed to `RGaussCode`; optional weights are stored separately | real-Gauss route |
| `table` | converted through Sage to PD, then stored as `PDCode` | PD route |
| `ostring` | parsed and routed to real-Gauss input | real-Gauss route |

This is intentionally not destructive normalization. A caller can supply a
presentation in the form natural to its source boundary.

### 3.2 Two kinds of state change

`setData` defaults to `reset=true`. This establishes a new mathematical-state
checkpoint:

```text
setData(input)  ≡  setData(input, reset=true)
    → clear representation/table/position validity flags
    → clear the expected Sage link state
    → install the new supplied representation
```

In ordinary use this is mutation of the mathematical object represented by the
MATLAB handle.

`reset=false` has a different meaning:

```text
setData(derived_view, reset=false)
    → retain previously available representations
    → add or refresh one derived representation
```

`calcVG2G` uses this second path. It reads retained RGauss data, inserts virtual
crossings for a visualisable full Gauss representation, and stores the result
without clearing the RGauss availability flag.

### 3.3 Cache guards

`formatFlag` is a ledger of available or valid representations and tables:

```text
G, PD, VG, RPD, TC, TRC, TE, TRE, TS, HM, RHM, pos
```

A getter or calculator checks its corresponding flag:

```text
flag true  → return/reuse the stored value
flag false → calculate the view, store it, and mark it available
```

Thus the flag blocks redundant recomputation. It is an availability/cache flag,
not necessarily an input-provenance record: a flag may become true either
because that form was supplied or because a getter derived it later.

### 3.4 Active route

`virtualFlag` selects whether generic access prefers the real-Gauss or full-Gauss
route. For example, `getGC` delegates to `getRGaussCode` when the flag is true
and to `getGaussCode` otherwise.

`virtualFlag` is therefore an active-route flag, not a sufficient record of the
original input kind. In the current implementation, materialising `GaussCode`
through `calcVG2G` calls the `GaussCode` branch of `setData(reset=false, ...)`,
which sets `virtualFlag=false`. A cache fill can consequently change later
routing even though the retained RGauss data remains available. This is a
current coupling, not a required part of the data-model rationale.

### 3.5 MATLAB state machine

```mermaid
flowchart TD
    I[Supported input form] --> R{reset?}
    R -- true / default --> N[Invalidate prior views and install new source state]
    R -- false --> A[Retain prior views and add one materialised view]

    N --> S{Stored input form}
    S -->|RGauss or o-data| VG[RGauss + ROrientation; VG available]
    S -->|Gauss| G[Gauss + orientation; G available]
    S -->|PD or table| PD[PD + orientation; PD available]

    VG --> C[Preferred real-Gauss computation]
    G --> RG[getRGaussCode: derive/cache RGauss]
    PD --> PG[via Gauss/Sage as required]
    PG --> RG
    RG --> C

    VG --> V2G[calcVG2G when full representation is required]
    V2G --> G2[Store Gauss with reset=false]
    G2 --> T[PD, tables, layout, positions]
    G --> T
    PD --> T
```

## 4. What is good about this design

### PR1 — Preserve source information

An eager conversion to one representation can erase input labels, presentation
choices, or comparison identity. Retaining the supplied form supports
reproducibility, debugging, and MATLAB/Sage/Julia parity work.

### PR2 — Use RGauss as a computational interface

Real Gauss data describes the real-crossing traversal without making inserted
virtual crossings part of the mathematical calculation boundary. It is a
natural interface for moves, real-edge weights, and invariants.

### PR3 — Materialise visual and tabular views lazily

Planarisation, PD construction, tables, and positions need not run for every
calculation. Lazy materialisation avoids unnecessary Sage or layout work and
allows several representations to coexist.

### PR4 — Distinguish topology replacement from cache filling

The `reset=true` versus `reset=false` distinction captures an important semantic
boundary: replace the represented object versus add another view of the same
object.

## 5. Current implementation risks

### IS1 — Source, availability, and route are not fully separated

The implementation has three distinct questions:

```text
What was supplied?       source/provenance
What is available now?   cache validity
What should be used now? active route
```

`formatFlag` and `virtualFlag` answer parts of all three. After lazy conversion,
they cannot by themselves reconstruct the original input provenance.

### IS2 — Cache invalidation is distributed

Representation dependencies are implicit in setter and calculator methods.
Direct property mutation or an incomplete `reset=false` path can leave tables or
positions inconsistent with the mathematical state.

### IS3 — Read operations may have hidden effects

Getters can invoke conversion, Sage, cache writes, and route changes. A call that
looks observational may therefore change subsequent behavior.

### IS4 — Representation-local identity is unstable

PD edge labels, table rows, and coordinates are not global mathematical IDs.
Data attached to them must use a documented real-arc traversal or explicit key,
and must be remapped after a representation-changing conversion.

## 6. Julia target principle

Julia should preserve the MATLAB design intent without copying its handle-state
mechanism literally:

> The mathematical value is immutable; representation conversion is pure;
> optional memoisation is mutable but semantically invisible.

The target split is:

```text
immutable VirtualLink value
    = retained source presentation
    + mathematical annotations

pure representation functions
    = rgauss(link), gauss(link), pd_code(link), real_edges(link), ...

optional mutable workspace/cache
    = previously computed views for one exact immutable value
```

### 6.1 Represent ingress forms explicitly

A provisional shape is:

```julia
abstract type LinkSource end

struct RGaussSource <: LinkSource
    components::Vector{Vector{Int}}
    orientation::Vector{Int}
end

struct GaussSource <: LinkSource
    components::Vector{Vector{Int}}
    orientation::Vector{Int}
end

struct PDSource <: LinkSource
    code::Matrix{Int}
    orientation::Vector{Int}
end

struct VirtualLink{S<:LinkSource,W}
    source::S
    weights::W
    strand_types::Vector{Char}
end
```

The exact field containers remain provisional. Julia `struct` is only shallowly
immutable when it contains mutable arrays; constructors must take ownership by
copying or otherwise prevent external mutation, and public APIs must not expose
unguarded mutable internals.

The source type records the supplied representation directly. A cache flag is
not needed to infer provenance.

### 6.2 Use dispatch instead of an active-route flag

```julia
rgauss(link::VirtualLink{<:RGaussSource}) = link.source
rgauss(link::VirtualLink{<:GaussSource}) = remove_virtual_crossings(link.source)
rgauss(link::VirtualLink{<:PDSource}) = pd_to_rgauss(link.source)
```

Algorithms request the representation they require explicitly. Materialising a
Gauss view must not change what later generic operations mean.

### 6.3 Return new values for mathematical transformations

```julia
mirror_manifold(link) -> new_link
cyclic_shift(link, strand, amount) -> new_link
mp_move(link, parameter, vertices) -> new_link
connected_sum(left, right) -> new_link
```

A transformation creates a new mathematical value and leaves the original
available for comparison or reuse. Mutating `!` APIs, if retained, should belong
to a builder or interactive workspace rather than define the core value
semantics.

### 6.4 Keep memoisation separate

Begin with uncached pure conversions unless profiling establishes a material
cost. When repeated interactive work requires caching, use a separate typed
cache associated with one exact immutable value:

```julia
mutable struct LinkCache
    rgauss::Union{Nothing,RGaussData}
    gauss::Union{Nothing,GaussData}
    pd::Union{Nothing,PDData}
    real_edges::Union{Nothing,Vector{RealEdge}}
end

struct LinkWorkspace{L}
    value::L
    cache::LinkCache
end
```

Creating a transformed `VirtualLink` creates a new workspace or empty cache.
There is no invalidation graph across mathematical mutations because an old
cache remains attached only to the old value.

An internal mutable cache inside an otherwise immutable `VirtualLink` is also
possible, but then equality, hashing, serialization, and mathematical identity
must ignore the cache, and concurrent access requires explicit synchronization.
A separate workspace is the clearer default.

## 7. Required semantic invariants

Any Julia design derived from this rationale should preserve the following:

1. **Input retention.** The supplied representation and its provenance remain
   recoverable without inferring them from a derived-view flag.
2. **Preferred computation.** Algorithms may use RGauss as their common view
   without rewriting the retained source.
3. **Pure observation.** Requesting another representation does not change the
   mathematical value or later dispatch meaning.
4. **Explicit transformation.** A mathematical move returns a new value unless a
   clearly named builder/workspace mutation API is requested.
5. **Disposable cache.** Deleting every derived cache changes performance only,
   not results or identity.
6. **No unstable external identity.** PD labels, table rows, and positions are
   local views unless a separate stable-key contract is defined.
7. **Reproducible conversion.** Derived representations can be rebuilt from the
   retained source plus explicit algorithm/version parameters.
8. **Separated evidence.** Replay coordinates and rendered figures do not become
   mathematical source data merely because they are stored.

## 8. Migration direction

The current Julia `mutable struct VirtualLink` and `set_data!` reproduce part of
the MATLAB surface but collapse source provenance, active representation, and
working state into one mutable object. A future refactor should proceed in this
order:

1. add explicit source-representation types and preserve received payloads;
2. implement pure `rgauss`, `gauss`, and PD conversion functions with tests;
3. make mathematical transformations return new values;
4. move optional memoisation into a separate workspace/cache;
5. port tables, geometry, and replay paths as derived views;
6. retire `real_only` and mutation-based route selection only after equivalent
   behavior and source retention are verified.

This refactor must not be justified as “making RGauss the sole canonical data.”
Its purpose is to separate input authority, calculation authority, and derived
state while retaining broad input support.

## 9. Non-goals

This document does not yet specify:

- a canonical normal form for virtual-link equivalence classes;
- equality or hashing across different source presentations;
- a stable global edge/crossing identity scheme;
- the final Julia public API or exact container types;
- cache eviction, serialization, or thread-safety details;
- performance acceptance thresholds;
- changes to MATLAB behavior.

Those require separate decisions and executable tests.

## 10. Source traceability

| Concern | Current source |
|---|---|
| MATLAB state fields and flags | `Manifold/VirtualLink/VirtualLink.m`, properties and `formatFlag0` |
| Input handling and reset semantics | `VirtualLink.setData` |
| Active representation selection | `VirtualLink.getGC` and `virtualFlag` |
| Real-Gauss derivation/cache | `VirtualLink.getRGaussCode` |
| Full-Gauss materialisation | `VirtualLink.calcVG2G` and `vg2g` |
| PD derivation | `VirtualLink.getPDCode`, `getRPDCode` |
| Weights and real-edge ordering | `VirtualLink.setWeight`, `calcREdgeTable` |
| MATLAB plotting/cache behavior | `docs/specs/virtual-link-plot.md`, §§18–24 |
| Current Julia state model | `src/virtual_link.jl` |

## 11. Decision status

Accepted as design direction in the owner conversation:

- broad ingress support is intentional;
- received input should be retained;
- RGauss is the preferred computational view, not the sole input authority;
- flags currently guard cached recomputation;
- default MATLAB `setData` is mathematical-object mutation;
- Julia should separate immutable mathematical values from optional mutable
  memoisation.

The exact Julia types, equality semantics, and cache API remain open.
