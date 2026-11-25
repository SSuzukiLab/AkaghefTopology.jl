"""Generic names for moves and for the search database. Presentations add methods; this module has no implementation."""
module MoveAPI

"""
    movable(input::HeegaardCode, row; kwargs...) -> Vector{Topology.HeegaardCodes.LCMoveWitness}

Enumerate supported Heegaard move witnesses through the common presentation interface.

# Arguments
`input` satisfies the [Heegaard move domain](../specs/heegaard-moves.md).
`row` is a catalog dictionary or its string identifier; row-dependent
selection keywords have the domain of `Topology.HeegaardCodes.movable`.
Unknown identifiers raise `LCMoveError` naming the catalog failure. Input,
mark-transport, support and missing-selection errors propagate from the
selected query. Arbitrary dictionaries are not certified catalog rows.
The one-argument overload also accepts `LickorishCode` and selects its default
finite witness pool. Named-row LC queries use `Topology.HeegaardCodes.movable`.
Knot and VirtualLink methods have their separate contracts and selectors.

# Returns
A vector of the [LC move witness records](../specs/heegaard-moves.md#comparison-and-move-witnesses)
returned by the corresponding Heegaard interface. Each record refers to its
source presentation, selected support and catalog relation.

# Guarantees
1. A named-row HC query has the same result as `Topology.HeegaardCodes.movable(input,row;kwargs...)`. [proved/agent: by construction]
2. The one-argument LC/HC query returns that interface's default finite pool; it omits the inverse-stabilization row and suppresses row-level `LCMoveError` failures. [proved/agent: by construction]
3. The input presentation is not modified. [proved/agent: by construction]

# Not guaranteed
A returned relation tag does not itself prove preservation of a represented
manifold. The default pool does not enumerate arbitrary twist powers or all
geometric supports. Empty output establishes only absence from the selected
implemented search range. A one-argument query is not provided for every
registered presentation type.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — common move interface.
cf:lib/topology.jl/docs/specs/heegaard-moves.md §§4 and 6 — witness meaning and finite selection ranges.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> length(movable(hc, "lc.twist.annulus"))
4
```

Role: access
"""
function movable end

"""
    move(input::HeegaardCode, witness::LCMoveWitness) -> Topology.HeegaardCodes.LCMoveResult

Apply a supported Heegaard move witness through the common presentation interface.

# Arguments
The source-binding and supported transport conditions of
[Heegaard move witnesses](../specs/heegaard-moves.md#comparison-and-move-witnesses)
apply. Use a witness from `movable(input,row)` for the same source and an
unmodified catalog. The selected Heegaard operation checks its supported
binding, support and transport data; failures raise `LCMoveError` with the
reason. Validation, catalog and conversion failures propagate. Hand-built
lower-level witnesses retain the specification's trusted-constructor boundary.

# Returns
An `LCMoveResult` with transformed `output`, its `inverse` witness, and a
`correspondence` dictionary, with the meanings in the
[witness specification](../specs/heegaard-moves.md#comparison-and-move-witnesses).
Named-row LC operations remain available through `Topology.HeegaardCodes`.

# Guarantees
1. The result equals `Topology.HeegaardCodes.move(input,witness)` and has that selected operation's supported presentation meaning. [proved/agent: by construction]
2. The source presentation is not modified by this call. [proved/agent: by construction]

# Not guaranteed
General manifold preservation does not follow from a successful operation:
the relation is row-dependent and a partial-family twist can change the
manifold. The output need not be a canonical or minimal representative, nor
need an inverse recover the exact original text. Other presentation types
use their own method contracts.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — common move interface.
cf:lib/topology.jl/docs/specs/heegaard-moves.md §4.1 — checked witness application and relation scope.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> r = move(hc, first(movable(hc, "lc.twist.annulus"))); back = move(r.output, r.inverse);

julia> H.hc_search_state(back.output; identity=:canonical).state_id == H.hc_search_state(hc; identity=:canonical).state_id
true
```

Role: move
"""
function move end

"""
    move!(v::VirtualLink, type, arg::NamedTuple) -> VirtualLink

Replace a VirtualLink by the result of a move selected through a named record.

# Arguments
The [record selector coordinates](../specs/virtual-link-moves.md#41-namedtuple-argument-aliases)
and keyword `move!` source, local-placement and coefficient conditions apply.
`edge_occurrences` and `boundary` are also recognized; `strict` defaults to
`true`. Unknown fields and conflicting aliases are not diagnosed as schema
errors. Selection, unsupported-family and application failures propagate from
the selected overload. Other presentation types have no mutating move method.

# Returns
The same `VirtualLink` object, with the selected successful result installed.

# Guarantees
1. Successful return has the local presentation and supported transport meaning of the keyword overload for the recognized fields. [proved/agent: by construction]
2. Success replaces the source's presentation data in place; selection or application failure leaves those data unchanged. The selector record is unchanged. [proved/agent: by construction]

# Not guaranteed
The record is not a schema-validated certificate. The result is not
canonicalized and no general manifold, combing or framing equivalence is
certified. Replaced containers need not update references to old field arrays.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §§4, 12 and 13 — selectors, coefficient transport and mutation boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,2,-2,3,-3,1,0]], orientation=[1,1,1]);

julia> move!(v, "MP", (p="C3", v=[1,2], strict=true));

julia> count(identity, Topology.strand_cuts(v))
1
```

Role: move
"""
function move! end

"""
    search_state(x::Union{HeegaardCode,VirtualLink}; kwargs...) -> CanonicalState

Return the presentation-specific identity record used as a search vertex.

# Arguments
For HC, the input satisfies the [HC data model](../specs/heegaard-code.md),
and `identity=:canonical` or `:literal` selects its search writing. Invalid
identity names raise `ArgumentError`; serialization and canonicalization
errors propagate. For VirtualLink, use the [data model](../specs/virtual-link-data-model.md)
and the [search domain](../specs/ograph-search-db.md#contract); no keywords
are accepted. Input coherence is the caller's obligation in both domains.

# Returns
The `CanonicalState` of `hc_search_state` or `ograph_search_state`, with
serialized presentation, hash, versioned ID, raw writing payload, domain,
representation and complexity. Field meanings are in
[Textual state identity](../specs/move-graph-store.md#textual-state-identity).

# Guarantees
1. Each result equals its presentation-specific search-state record for the supplied arguments. [proved/agent: by construction]
2. Neither the input nor a database is modified. [proved/agent: by construction]

# Not guaranteed
Search-state identity does not classify manifolds. HC canonicality has the
scope of `hc_canonical`; o-graph identity retains the exact unweighted writing
and does not identify arbitrary relabellings or move-equivalent diagrams.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — domain dispatch of the search interface.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = H.resolve_hc(lc, Dict("l1"=>:beta); standard=true);

julia> Topology.MoveAPI.search_state(h; identity=:literal).representation
"hc-search"
```

Role: convert
"""
function search_state end

"""
    put_search_state!(store::MoveGraphStoreDB, hc::HeegaardCode; identity=:canonical) -> CanonicalState

Store a Heegaard search state and its retained writing through the common interface.

# Arguments
M7 in [the store specification](../specs/move-graph-store.md#numbered-conditions)
and the input requirements of `Topology.HeegaardCodes.hc_search_state` apply.
The HC schema and identity mode are checked. `identity` accepts `canonical`
or `literal` as a symbol or string; unsupported modes or schemas raise
`ArgumentError` naming the requirement. HC validation, canonicalization,
serialization and SQLite errors propagate. The store must be writable.

# Returns
A `CanonicalState` for the input's selected HC identity, with the
[identity and writing meanings](../specs/heegaard-code.md#comparison-records).
It contains state ID/hash, serialization, version, domain/representation,
genus complexity and the input's writing payload.

# Guarantees
1. Successful return stores the computed state and an HC writing record under its ID, with the same effects as `Topology.HeegaardCodes.put_hc_search_state!`. An existing ID retains its previous HC writing. [proved/agent: by construction]
2. The identity is a key for the selected HC comparison, not a classification of represented manifolds. [proved/agent: by construction]
3. The database can change; the HC input is unchanged. [proved/agent: by construction]

# Not guaranteed
For an existing identity, the returned input writing can differ from the
retained writing later loaded from the store. The writes are not one
whole-call transaction. This operation validates no move and certifies no
manifold equivalence. The VirtualLink overload has its own contract and
stores an unweighted search presentation.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Textual state identity — persistence and mathematical evidence boundary.
cf:lib/topology.jl/docs/specs/heegaard-code.md, Comparison records — HC identity meaning.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> store = open_store(":memory:"); H.migrate_hc_search_store!(store);

julia> same = try
           s = Topology.MoveAPI.put_search_state!(store, hc; identity=:canonical)
           loaded = H.load_hc_search_state(store, s.state_id)
           H.hc_search_state(loaded; identity=:canonical).state_id == s.state_id
       finally
           close_store!(store)
       end;

julia> same
true
```

Role: other
"""
function put_search_state! end

"""
    load_search_state(store::MoveGraphStoreDB, presentation::Union{Type{HeegaardCode},Type{VirtualLink}}, state_id::AbstractString) -> Union{Nothing,HeegaardCode,VirtualLink}

Reconstruct a stored search presentation of the selected type.

# Arguments
`presentation` is `HeegaardCode` or `VirtualLink`; omission selects the
o-graph interface. M7 and M11
([store conditions](../specs/move-graph-store.md#numbered-conditions)) apply.
Unknown state identity is accepted. HC schema and reconstruction validation
can raise `ArgumentError` or HC validation exceptions; malformed stored JSON,
fields and database operations can raise parsing, lookup, conversion or SQLite
errors. The o-graph method propagates the corresponding search reconstruction
errors. Unsupported presentation types have no method.

# Returns
`nothing` if absent from the selected search tables. Otherwise a fresh
`HeegaardCode` with saved words, roles and annotations and no source paths,
or an unweighted integer-mode `VirtualLink` with the saved search presentation.

# Guarantees
1. The result equals the selected `load_hc_search_state` or `load_ograph_search_state` reconstruction, with that presentation's supported record meaning. [proved/agent: by construction]
2. Loading does not change stored rows; modifying the returned presentation does not update the database. [proved/agent: by construction]

# Not guaranteed
This does not recompute/authenticate the stored state identity, replay move
history or certify a represented manifold. It returns the retained writing,
not a requested equivalent representative. An o-graph search record does not
recover edge weights removed from an earlier input.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, M11 — stored writing versus mathematical identity.
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — presentation-directed common search interface.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> store = open_store(":memory:"); H.migrate_hc_search_store!(store);

julia> same = try
           s = H.put_hc_search_state!(store, hc; identity=:canonical)
           loaded = Topology.MoveAPI.load_search_state(store, H.HeegaardCode, s.state_id)
           H.hc_search_state(loaded; identity=:canonical).state_id == s.state_id
       finally
           close_store!(store)
       end;

julia> same
true
```

Role: access
"""
function load_search_state end

"""
    begin_search_run!(store::MoveGraphStoreDB, presentation::Union{Type{HeegaardCode},Type{VirtualLink}}, script_path::AbstractString; label::AbstractString, scope=Dict{String,Any}()) -> String

Register the provenance of a presentation-specific move-search run and return its run identifier.

# Arguments
`presentation` is `HeegaardCode` or `VirtualLink`. M7–M8
([store conditions](../specs/move-graph-store.md#numbered-conditions)) apply.
The selected method delegates to `HeegaardCodes.begin_hc_search_run!` or
`begin_ograph_search_run!`; unsupported schema or absent files raise
`ArgumentError`, and read, serialization and database failures propagate.
`label` is required. Omitting the presentation type selects the existing
o-graph interface. Unsupported presentation types have no method.

# Returns
A `String` identifying the new running record in the selected search tables.
Use that presentation's run accessor or `finish_search_run!` with the same
type argument to inspect or finish it.

# Guarantees
1. Successful registration stores the supplied label/scope and the script path and content digest with running status, using the selected presentation's run interface. [proved/agent: by construction]
2. The result is computational provenance for the selected finite move-search workflow, as defined in the specification. It certifies no mathematical equivalence. [proved/agent: by construction]
3. The store is modified; the script and scope are unchanged and the script is not executed. [proved/agent: by construction]

# Not guaranteed
No search is started or automatically completed. Scope is caller metadata,
not a bound enforced by registration. Script bytes are not archived and
run identity is not a mathematical state identity. The delegated timestamp-based
identifier is not promised unique for identical same-timestamp registrations.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Search-run provenance — presentation-specific experiment record.

# Examples
```jldoctest
julia> store = open_store(":memory:"); migrate_ograph_search_store!(store); path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl");

julia> count = try
           Topology.MoveAPI.begin_search_run!(store, VirtualLink, path; label="literal-02-smoke")
           Topology.MoveAPI.begin_search_run!(store, VirtualLink, path; label="literal-02-reuse")
           verify_ograph_search_store(store).runs
       finally
           close_store!(store)
       end;

julia> count
2
```

Role: other
"""
function begin_search_run! end

"""
    finish_search_run!(store::MoveGraphStoreDB, presentation::Union{Type{HeegaardCode},Type{VirtualLink}}, run_id::AbstractString; status::AbstractString, counts=Dict{String,Any}()) -> Union{Nothing,NamedTuple}

Record the terminal status of a presentation-specific move-search run and return its summary, or nothing for an absent run.

# Arguments
`presentation` is `HeegaardCode` or `VirtualLink`; omitting it selects the
o-graph interface. M7 and M10
([store conditions](../specs/move-graph-store.md#numbered-conditions)) apply.
Invalid terminal status raises `ArgumentError`; serialization, schema and
database errors propagate from the selected `finish_hc_search_run!` or
`finish_ograph_search_run!`. An absent or already terminal run is accepted.

# Returns
`nothing` for an absent run, otherwise a named tuple with `run_id`,
`script_path`, `script_sha256`, `label`, decoded `scope`, `status` and reported
`counts`. The HC result adds stored `operations`, `successes`, `rejects` and
`errors`; the o-graph result adds linked `expansions`, `candidates`, `successes`,
`rejects` and `errors`. These are integer execution totals, distinct from the
caller-reported counts, with the [provenance meanings](../specs/move-graph-store.md#search-run-provenance).

# Guarantees
1. A running record is updated to the requested status/counts with a finish time; an already terminal record is unchanged. The selected presentation's run summary is returned. [proved/agent: by construction]
2. The operation records a caller-reported finite experiment outcome and does not verify mathematical search completeness. [proved/agent: by construction]
3. The store can change; the supplied counts are unchanged. [proved/agent: by construction]

# Not guaranteed
This does not cancel a running computation, compare the supplied counts with
the database totals, or prove manifold equivalence. Terminal status is not
an exhaustiveness certificate. Missing or malformed tables remain errors,
even though an absent run in valid tables returns `nothing`.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Search-run provenance and M10 — recorded outcome and presentation-specific terminal statuses.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; store = open_store(":memory:"); H.migrate_hc_search_store!(store);

julia> status = try
           id = H.begin_hc_search_run!(store, joinpath(pkgdir(Topology), "test", "heegaard", "hc_search.jl"); label="round-trip")
           Topology.MoveAPI.finish_search_run!(store, H.HeegaardCode, id; status="completed").status
       finally
           close_store!(store)
       end;

julia> status
"completed"
```

Role: other
"""
function finish_search_run! end

"""
    expand_search_state!(store::MoveGraphStoreDB, run_id::AbstractString, source::HeegaardCode; depth::Integer=0, identity=:canonical, room::Integer=typemax(Int), check_inverse::Bool=true, check_invariants::Bool=false, start_snap=nothing, problems=nothing, on_problem=nothing) -> NamedTuple

Expand one Heegaard search state and record its candidate outcomes in the move store.

# Arguments
M7 and M9 ([store conditions](../specs/move-graph-store.md#numbered-conditions))
apply. The HC schema and identity mode are checked; unsupported values raise
`ArgumentError`. `depth` converts to `Int`; `room` is the allowance for newly
stored states. Nonpositive room stops at the next candidate needing a new
state. Source-state validation/key failures and database errors outside the
candidate handler propagate. Within a candidate, local-move errors are recorded
as rejections and other errors as errors; interrupts are rethrown.
`check_inverse` enables inverse diagnostics. Invariant diagnostics require all
of `check_invariants=true`, a compatible `start_snap` and a mutable `problems`
collection. Optional `on_problem` accepts each diagnostic record; callback
failures can propagate. The source's mathematical domain is defined by
[Heegaard moves](../specs/heegaard-moves.md).

# Returns
A named tuple with `records`, `created`, `stopped`, `source_state_id` and
`candidate_count`. `created` counts new target states, `stopped` records the
room stop, and `candidate_count` is the entire enumerated candidate count.
Each record has `operation_id`, move `row_id`, `outcome`, `inverse_status`,
optional `target_state_id` and `target`, and target `depth`. These statuses
have the [search operation meanings](../specs/move-graph-store.md#search-run-provenance).

# Guarantees
1. This HC method has exactly the results and writes of `HeegaardCodes.expand_hc_search_state!` with these arguments; records distinguish successful local applications from recorded failures. [proved/agent: by construction]
2. A room stop leaves a prefix of candidate outcomes, not a complete expansion. A reported successful inverse check compares the returned inverse's state identity with the source identity. [proved/agent: by construction]
3. The store and supplied diagnostic collection can change; the source is unchanged by the standard move operations. Callback side effects remain the caller's responsibility. [proved/agent: by construction]

# Not guaranteed
A successful forward record need not have a successful inverse check.
Recorded failures do not imply an exception from this call, and prior writes
are not rolled back as one whole expansion. The room allowance is not a time
limit; it can require applying the next candidate before deciding to stop.
This finite search evidence does not establish unrestricted manifold
equivalence or completeness of the global move calculus.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Search-run provenance and M9 — finite candidate outcomes and diagnostics.
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — HC forwarding through the common move interface.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes; lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("l1" => :beta); standard=true);

julia> store = open_store(":memory:"); H.migrate_hc_search_store!(store);

julia> count = try
           id = H.begin_hc_search_run!(store, joinpath(pkgdir(Topology), "test", "heegaard", "hc_search.jl"); label="round-trip")
           Topology.MoveAPI.expand_search_state!(store, id, hc; room=0, check_inverse=false).candidate_count
       finally
           close_store!(store)
       end;

julia> count
9
```

Role: move
"""
function expand_search_state! end

"""
    search_neighbors(store::MoveGraphStoreDB, presentation::Union{Type{HeegaardCode},Type{VirtualLink}}, state_id::AbstractString; kwargs...) -> Vector

Read successful outgoing transitions from a presentation-specific stored search graph.

# Arguments
`presentation` is `HeegaardCode` or `VirtualLink`; omission selects o-graphs.
Requires M7 ([conditions](../specs/move-graph-store.md#numbered-conditions)).
HC accepts `row_ids=nothing` or an iterable of names convertible to strings;
its schema is checked. VirtualLink accepts `kinds=("MP","02")`; unsupported
kinds raise `ArgumentError` naming MP, 02 and PS. Empty filters give no rows.
The state ID need not exist. Conversion, decoding and database errors propagate.

# Returns
HC records have operation/run/source/target IDs, depth, row ID, inverse
status and decoded candidate. O-graph rows are the successful-edge records
of `ograph_search_neighbors`, including kind, direction, engine version and
encoded candidate. See [Search-run provenance](../specs/move-graph-store.md#search-run-provenance).

# Guarantees
1. Results equal the selected `hc_search_neighbors` or `ograph_search_neighbors` query for the supplied source and filters. [proved/agent: by construction]
2. HC results are ordered by depth, row ID and operation ID; o-graph results by kind and operation ID. [proved/agent: by construction]
3. Neither store nor arguments are modified. [proved/agent: by construction]

# Not guaranteed
The query neither generates moves nor replays stored edges. An empty result
does not imply that the mathematical state has no neighbors.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md Search-run provenance — recorded directed transitions.

# Examples
```jldoctest
julia> s = open_store(); migrate_ograph_search_store!(s);

julia> isempty(Topology.MoveAPI.search_neighbors(s, VirtualLink, "absent"; kinds=String[]))
true

julia> close_store!(s);
```

Role: access
"""
function search_neighbors end

"""
    query_search_paths(store::MoveGraphStoreDB, presentation::Union{Type{HeegaardCode},Type{VirtualLink}}, source_state_id::AbstractString; kwargs...) -> Union{Nothing,Vector}

Read bounded directed paths in a presentation-specific stored move graph.

# Arguments
`presentation` is `HeegaardCode` or `VirtualLink`; omission selects o-graphs.
M7 and M12 ([store conditions](../specs/move-graph-store.md#numbered-conditions))
apply. HC keywords are required `target_state_id`, optional `row_ids=nothing`,
`min_depth=1`, `max_depth=12`, `limit=100`. Its schema is checked.
VirtualLink keywords select exactly one of `target_state_id` and
`target_geometry_hash`, with `kinds=("MP","02")`, `kind_pattern=nothing`,
`min_depth=1`, `max_depth=12`, `limit=1000` by default.
Target/bound violations and unsupported selected move kinds raise
`ArgumentError`. A supplied o-graph pattern sets both depths to its length
after checking the supplied bounds. An empty pattern is not supported by the
current query. Conversion, stored-data and SQLite failures propagate.

# Returns
For HC, a vector of `(depth,steps)` records: path depth counts edges; each
step has operation/source/target IDs, stored exploration depth, move `row_id`,
`inverse_status` and decoded `candidate`. For VirtualLink, a vector of step
vectors with operation/source/target IDs, `move_kind`, `direction`,
`engine_version` and decoded `candidate`; an empty `kinds` selection without
a pattern returns `nothing`. No matching paths otherwise gives an empty vector.
These are the [stored path records](../specs/move-graph-store.md#stored-search-paths).

# Guarantees
1. The result equals the selected `query_hc_search_paths` or `query_ograph_search_paths` query: only successful stored transitions satisfying the chosen filters contribute. [proved/agent: by construction]
2. Results are ordered by path depth and operation-ID sequence and truncated to `limit`. HC paths have no repeated state; o-graph paths permit the source to recur only as their final endpoint. [proved/agent: by construction]
3. The query neither generates moves nor changes stored rows. [proved/agent: by construction]

# Not guaranteed
No edge is replayed and no manifold equivalence is independently certified.
A missing bounded stored path is not evidence of mathematical inequivalence.
A geometry hash is not a manifold identity; an unweighted o-graph path need
not admit the requested later weight transport. Path shape differs between
the two presentation types.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md, Stored search paths — finite stored graph and evidence boundary.
cf:lib/topology.jl/docs/specs/ograph-search-db.md, Exact H-insertion feasibility — final-source return and geometric path selection.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; gauss=[[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]], orientation=[1,1,1,1,0,0]);

julia> store = open_store(":memory:"); migrate_ograph_search_store!(store);

julia> same = try
           s = put_ograph_search_state!(store,v)
           r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); label="literal-02-smoke")
           edges = expand_ograph_search_state!(store,r,v,"02"; selector=Dict("direction"=>"forward"))
           target = String(first(edges).target_state_id)
           path = query_ograph_search_path(store,s.state_id; target_state_id=target,kinds=("02",),max_depth=1)
           Topology.MoveAPI.query_search_paths(store,VirtualLink,s.state_id; target_state_id=target,kinds=("02",),max_depth=1,limit=10) == [path]
       finally
           close_store!(store)
       end;

julia> same
true
```

Role: access
"""
function query_search_paths end

"""
    verify_search_store(store::MoveGraphStoreDB) -> NamedTuple

Check o-graph search-schema metadata and foreign keys and count its stored records.

# Arguments
Requires an open store with readable o-graph search tables under M7
([conditions](../specs/move-graph-store.md#numbered-conditions)). Missing or
unsupported version metadata and foreign-key violations raise `ArgumentError`
identifying the problem. Missing tables, malformed version text and database
failures can raise their parsing or database errors first.

# Returns
Integer fields `schema_version`, `states`, `runs`, `expansions`, `operations`,
`successful_edges` and `weighted_trials`, with the record-family meanings in
the [search specification](../specs/ograph-search-db.md#workflow).

# Guarantees
1. Results and validation equal `verify_ograph_search_store(store)`. Success means the schema version matches and the foreign-key check found no violations. [proved/agent: by construction]
2. The database is not modified. [proved/agent: by construction]

# Not guaranteed
No move or weighted trial is replayed, identity recomputed, or search
completeness established. There is no HC method at this generic entry point.

# Sources
cf:lib/topology.jl/docs/specs/module-structure.md §3.1 — the currently implemented o-graph verification adapter.

# Examples
```jldoctest
julia> s = open_store(); migrate_ograph_search_store!(s);

julia> path = joinpath(pkgdir(Topology), "test", "ograph_search_db.jl"); begin_ograph_search_run!(s,path; label="literal-02-smoke"); begin_ograph_search_run!(s,path; label="literal-02-reuse");

julia> Topology.MoveAPI.verify_search_store(s).runs
2

julia> close_store!(s);
```

Role: condition
"""
function verify_search_store end

export movable, move, move!,
       search_state, put_search_state!, load_search_state,
       begin_search_run!, finish_search_run!, expand_search_state!,
       search_neighbors, query_search_paths, verify_search_store

end
