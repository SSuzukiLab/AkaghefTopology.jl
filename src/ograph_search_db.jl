using Dates
using JSON
using SHA

const OGRAPH_SEARCH_SCHEMA_VERSION = 1
const OGRAPH_SEARCH_STATE_VERSION = "ograph-search-state/1"
const OGRAPH_SEARCH_ENUMERATOR_VERSION = "ograph-search-enumerator/1"
const OGRAPH_SEARCH_ENGINE_VERSION = "ograph-search-engine/1"

function _ograph_search_schema_path()
    normpath(joinpath(@__DIR__, "..", "db", "ograph_search_schema.sql"))
end

function migrate_ograph_search_store!(store::MoveGraphStoreDB)
    MoveGraphStore.migrate!(store)
    MoveGraphStore._transaction(store.db) do
        MoveGraphStore._execute_schema!(store.db, read(_ograph_search_schema_path(), String))
    end
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph search schema version"))
    store
end

function _require_ograph_search_store(store::MoveGraphStoreDB)
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("o-graph search schema is absent; call migrate_ograph_search_store!"))
    store
end

function _ograph_search_unweighted(link::VirtualLink)
    projected = VirtualLink(weight_mode=:integer)
    kwargs = (orientation=copy(link.orientation), strand_types=copy(link.strand_types))
    if link.real_only
        set_data!(projected; rgauss=deepcopy(link.gauss), kwargs...)
    else
        set_data!(projected; gauss=deepcopy(link.gauss), kwargs...)
    end
    projected
end

function _ograph_search_state_data(link::VirtualLink)
    projected = _ograph_search_unweighted(link)
    real_gauss, real_orientation = real_gauss_code(projected)
    raw = Dict(
        "schema" => OGRAPH_SEARCH_STATE_VERSION,
        "gauss" => projected.gauss,
        "orientation" => projected.orientation,
        "real_only" => projected.real_only,
        "strand_types" => string.(projected.strand_types),
    )
    geometry = Dict(
        "schema" => "ograph-search-geometry/1",
        "real_gauss" => real_gauss,
        "real_orientation" => real_orientation,
    )
    serialization = _ograph_canonical_json(raw)
    geometry_serialization = _ograph_canonical_json(geometry)
    (projected=projected, raw=raw, serialization=serialization,
     raw_hash=_ograph_sha256(serialization),
     geometry_hash=_ograph_sha256(geometry_serialization),
     real_gauss=real_gauss, real_orientation=real_orientation)
end

"""Return the unweighted exact-presentation state used by the search graph.

Weights are deliberately projected out. Weighted validity belongs to
`replay_weighted_ograph_path`, after a geometric path has been selected.
"""
function ograph_search_state(link::VirtualLink)
    data = _ograph_search_state_data(link)
    CanonicalState(data.serialization; area=AREA_3,
        representation="o-graph-search", payload=data.raw,
        canonical_version=OGRAPH_SEARCH_SCHEMA_VERSION,
        complexity=length(data.real_orientation))
end

function put_ograph_search_state!(store::MoveGraphStoreDB, link::VirtualLink)
    _require_ograph_search_store(store)
    data = _ograph_search_state_data(link)
    state = put_state!(store, ograph_search_state(link))
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_states(
          state_id,state_version,raw_hash,geometry_hash,raw_gauss_json,
          orientation_json,real_gauss_json,real_orientation_json,real_vertex_count)
        VALUES (?,?,?,?,?,?,?,?,?)
    """, (state.state_id, OGRAPH_SEARCH_STATE_VERSION, data.raw_hash,
          data.geometry_hash, _ograph_canonical_json(data.projected.gauss),
          _ograph_canonical_json(data.projected.orientation),
          _ograph_canonical_json(data.real_gauss),
          _ograph_canonical_json(data.real_orientation),
          length(data.real_orientation)))
    state
end

function load_ograph_search_state(store::MoveGraphStoreDB, state_id::AbstractString)
    row = MoveGraphStore._one(store.db, """
        SELECT state.serialization FROM states AS state
        JOIN ograph_search_states AS search_state USING(state_id)
        WHERE state_id=?
    """, (String(state_id),))
    row === nothing && return nothing
    data = JSON.parse(String(row.serialization))
    link = VirtualLink()
    gauss = [Int.(component) for component in data["gauss"]]
    orientation = Int.(data["orientation"])
    strand_types = only.(String.(data["strand_types"]))
    if Bool(data["real_only"])
        set_data!(link; rgauss=gauss, orientation=orientation,
                  strand_types=strand_types)
    else
        set_data!(link; gauss=gauss, orientation=orientation,
                  strand_types=strand_types)
    end
    link
end

function _ograph_search_script_snapshot(script_path::AbstractString)
    path = abspath(String(script_path))
    isfile(path) || throw(ArgumentError("provenance script does not exist: $path"))
    path, bytes2hex(SHA.sha256(read(path)))
end

function begin_ograph_search_run!(store::MoveGraphStoreDB,
                                  script_path::AbstractString;
                                  label::AbstractString,
                                  scope=Dict{String,Any}())
    _require_ograph_search_store(store)
    path, digest = _ograph_search_script_snapshot(script_path)
    started = string(Dates.now(Dates.UTC))
    run_id = "ograph-search-run-" * _ograph_sha256(
        join((replace(path, '\\' => '/'), digest, String(label), started), '\0'))
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO ograph_search_runs(
          run_id,script_path,script_sha256,label,scope_json,status,counts_json,
          started_at_utc)
        VALUES (?,?,?,?,?,'running','{}',?)
    """, (run_id, replace(path, '\\' => '/'), digest, String(label),
          _ograph_canonical_json(scope), started))
    run_id
end


function finish_ograph_search_run!(store::MoveGraphStoreDB,
                                   run_id::AbstractString;
                                   status::AbstractString,
                                   counts=Dict{String,Any}())
    status in ("completed-found", "completed-none", "budget-exhausted", "failed") ||
        throw(ArgumentError("invalid terminal o-graph search status"))
    MoveGraphStore._dbexec(store.db, """
        UPDATE ograph_search_runs
        SET status=?,counts_json=?,finished_at_utc=?
        WHERE run_id=? AND status='running'
    """, (String(status), _ograph_canonical_json(counts),
          string(Dates.now(Dates.UTC)), String(run_id)))
    ograph_search_run(store, run_id)
end

function ograph_search_run(store::MoveGraphStoreDB, run_id::AbstractString)
    row = MoveGraphStore._one(store.db, """
        SELECT run.*,
               COUNT(DISTINCT link.expansion_id) AS expansions,
               COALESCE(SUM(expansion.candidate_count),0) AS candidates,
               COALESCE(SUM(expansion.success_count),0) AS successes,
               COALESCE(SUM(expansion.reject_count),0) AS rejects,
               COALESCE(SUM(expansion.error_count),0) AS errors
        FROM ograph_search_runs AS run
        LEFT JOIN ograph_search_run_expansions AS link ON link.run_id=run.run_id
        LEFT JOIN ograph_search_expansions AS expansion
          ON expansion.expansion_id=link.expansion_id
        WHERE run.run_id=? GROUP BY run.run_id
    """, (String(run_id),))
    row === nothing && return nothing
    (run_id=String(row.run_id), script_path=String(row.script_path),
     script_sha256=String(row.script_sha256), label=String(row.label),
     scope=JSON.parse(String(row.scope_json)), status=String(row.status),
     counts=JSON.parse(String(row.counts_json)), expansions=Int(row.expansions),
     candidates=Int(row.candidates), successes=Int(row.successes),
     rejects=Int(row.rejects), errors=Int(row.errors))
end

function _ograph_search_candidate_data(candidate::MoveCandidate)
    Dict(
        "type" => candidate.type,
        "parameter" => candidate.parameter,
        "vertices" => candidate.vertices,
        "edges" => [collect(edge) for edge in candidate.edges],
        "direction" => String(candidate.direction),
        "locations" => [collect(location) for location in candidate.locations],
        "ps_split_k" => candidate.ps_split_k === nothing ? nothing :
                        string(candidate.ps_split_k),
    )
end

_ograph_search_candidate_key(candidate::MoveCandidate) =
    _ograph_canonical_json(_ograph_search_candidate_data(candidate))

function _ograph_search_unique_candidates(candidates)
    seen = Set{String}()
    result = MoveCandidate[]
    for candidate in candidates
        key = _ograph_search_candidate_key(candidate)
        key in seen && continue
        push!(seen, key)
        push!(result, candidate)
    end
    result
end

function _ograph_search_mp_candidates(link::VirtualLink, selector)
    arities = haskey(selector, "vertex_arities") ?
        sort!(unique(Int.(selector["vertex_arities"]))) : [2, 3]
    all(arity -> arity in (2, 3), arities) ||
        throw(ArgumentError("MP vertex_arities must be drawn from [2,3]"))
    result = MoveCandidate[]
    n = length(link.orientation)
    for arity in arities,
        vertices in Iterators.product(ntuple(_ -> 1:n, arity)...)
        length(unique(vertices)) == arity || continue
        append!(result, movable(link, "MP"; vertices=collect(vertices)))
    end
    _ograph_search_unique_candidates(result)
end

function _ograph_search_02_candidates(link::VirtualLink, selector)
    direction = get(selector, "direction", "both")
    direction in ("forward", "inverse", "both") ||
        throw(ArgumentError("02 direction must be forward, inverse, or both"))
    result = MoveCandidate[]
    if direction in ("forward", "both")
        code, _ = real_gauss_code(link)
        for word in code, index in eachindex(word)
            isempty(word) && continue
            vertices = [word[index], word[mod1(index + 1, length(word))]]
            append!(result, movable(link, "02"; param="dir", vertices=vertices))
        end
    end
    if direction in ("inverse", "both")
        n = length(link.orientation)
        for first_vertex in 1:n, second_vertex in 1:n
            first_vertex == second_vertex && continue
            append!(result, movable(link, "02"; param="inv",
                                    vertices=[first_vertex, second_vertex]))
        end
    end
    _ograph_search_unique_candidates(result)
end

function _ograph_search_ps_candidates(link::VirtualLink, selector)
    parameters = haskey(selector, "parameters") ?
        String.(selector["parameters"]) : ["1", "2", "3", "4"]
    all(parameter -> parameter in ("1", "2", "3", "4"), parameters) ||
        throw(ArgumentError("PS parameters must be drawn from 1,2,3,4"))
    direction = get(selector, "direction", "both")
    direction in ("forward", "inverse", "both") ||
        throw(ArgumentError("PS direction must be forward, inverse, or both"))
    result = MoveCandidate[]
    if direction in ("inverse", "both")
        n = length(link.orientation)
        for parameter in parameters, first_vertex in 1:n, second_vertex in 1:n
            first_vertex == second_vertex && continue
            append!(result, movable(link, "PS"; param=parameter,
                                    vertices=[first_vertex, second_vertex]))
        end
    end
    if direction in ("forward", "both")
        code, _ = real_gauss_code(link)
        locations = [(component, position) for (component, word) in enumerate(code)
                     for position in eachindex(word)]
        for parameter in parameters, first_location in locations,
            second_location in locations
            first_location == second_location && continue
            append!(result, movable(link, "PS"; param=parameter,
                                    edge_occurrences=[first_location, second_location]))
        end
    end
    _ograph_search_unique_candidates(result)
end

function ograph_search_candidates(link::VirtualLink, move_kind::AbstractString;
                                  selector=Dict{String,Any}())
    projected = _ograph_search_unweighted(link)
    kind = String(move_kind)
    normalized = Dict(String(key) => value for (key, value) in pairs(selector))
    kind == "MP" && return _ograph_search_mp_candidates(projected, normalized)
    kind == "02" && return _ograph_search_02_candidates(projected, normalized)
    kind == "PS" && return _ograph_search_ps_candidates(projected, normalized)
    kind == "H" && throw(ArgumentError(
        "H is not a geometric search edge; insert it only during weighted replay"))
    throw(ArgumentError("o-graph search supports MP, literal 02, and PS"))
end

function _ograph_search_link_expansion!(store, run_id, expansion_id, reused)
    ordinal = Int(MoveGraphStore._one(store.db, """
        SELECT COALESCE(MAX(ordinal),0)+1 AS ordinal
        FROM ograph_search_run_expansions WHERE run_id=?
    """, (String(run_id),)).ordinal)
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_run_expansions(
          run_id,expansion_id,ordinal,reused) VALUES (?,?,?,?)
    """, (String(run_id), String(expansion_id), ordinal, reused ? 1 : 0))
end

function _ograph_search_operations(store, expansion_id)
    MoveGraphStore._rows(store.db, """
        SELECT operation.* FROM ograph_search_expansion_operations AS link
        JOIN ograph_search_operations AS operation USING(operation_id)
        WHERE link.expansion_id=? ORDER BY link.ordinal
    """, (String(expansion_id),))
end

"""Exhaust one reusable, unweighted move-neighborhood and persist its outcome."""
function expand_ograph_search_state!(store::MoveGraphStoreDB,
                                     run_id::AbstractString,
                                     source::VirtualLink,
                                     move_kind::AbstractString;
                                     selector=Dict{String,Any}())
    kind = String(move_kind)
    kind == "H" && throw(ArgumentError(
        "H is replay-only and cannot enter the unweighted graph"))
    kind in ("MP", "02", "PS") || throw(ArgumentError(
        "o-graph search supports MP, literal 02, and PS"))
    source_link = _ograph_search_unweighted(source)
    source_state = put_ograph_search_state!(store, source_link)
    selector_json = _ograph_canonical_json(
        Dict(String(key) => value for (key, value) in pairs(selector)))
    expansion_id = "ograph-search-expansion-" * _ograph_sha256(join((
        source_state.state_id, kind, OGRAPH_SEARCH_ENUMERATOR_VERSION,
        selector_json), '\0'))
    existing = MoveGraphStore._one(store.db,
        "SELECT * FROM ograph_search_expansions WHERE expansion_id=?",
        (expansion_id,))
    if existing !== nothing && String(existing.status) == "completed"
        _ograph_search_link_expansion!(store, run_id, expansion_id, true)
        return _ograph_search_operations(store, expansion_id)
    elseif existing !== nothing
        String(existing.status)=="building" || throw(ArgumentError(
            "existing search expansion failed: $expansion_id"))
        active=MoveGraphStore._one(store.db,"""
            SELECT COUNT(*) AS n FROM ograph_search_run_expansions AS link
            JOIN ograph_search_runs AS run USING(run_id)
            WHERE link.expansion_id=? AND run.status='running' AND run.run_id<>?
        """,(expansion_id,String(run_id)))
        Int(active.n)==0 || throw(ArgumentError(
            "search expansion is owned by another running run: $expansion_id"))
        # A terminated run can leave the row created before enumeration. The
        # actual edge batch is transactional; resume without discarding rows.
    end

    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO ograph_search_expansions(
          expansion_id,source_state_id,move_kind,enumerator_version,selector_json,
          status,candidate_count,success_count,reject_count,error_count)
        VALUES (?,?,?,?,?,'building',0,0,0,0)
    """, (expansion_id, source_state.state_id, kind,
          OGRAPH_SEARCH_ENUMERATOR_VERSION, selector_json))
    _ograph_search_link_expansion!(store, run_id, expansion_id, false)

    candidates = try
        ograph_search_candidates(source_link, kind; selector=selector)
    catch error_value
        error_json = _ograph_canonical_json(Dict(
            "type" => string(typeof(error_value)),
            "message" => sprint(showerror, error_value)))
        MoveGraphStore._dbexec(store.db, """
            UPDATE ograph_search_expansions
            SET status='failed',error_count=1,error_json=? WHERE expansion_id=?
        """, (error_json, expansion_id))
        rethrow()
    end

    success_count = 0
    reject_count = 0
    error_count = 0
    MoveGraphStore._transaction(store.db) do
        for (ordinal, candidate) in enumerate(candidates)
            candidate_data = _ograph_search_candidate_data(candidate)
            candidate_json = _ograph_canonical_json(candidate_data)
            operation_id = "ograph-search-operation-" * _ograph_sha256(join((
                source_state.state_id, OGRAPH_SEARCH_ENGINE_VERSION,
                candidate_json), '\0'))
            target_state_id = nothing
            outcome = "success"
            error_json = nothing
            try
                target = move(source_link, candidate)
                target_state_id = put_ograph_search_state!(store, target).state_id
                success_count += 1
            catch error_value
                expected = error_value isa Union{ArgumentError,DimensionMismatch,BoundsError}
                outcome = expected ? "rejected" : "error"
                error_json = _ograph_canonical_json(Dict(
                    "type" => string(typeof(error_value)),
                    "message" => sprint(showerror, error_value)))
                expected ? (reject_count += 1) : (error_count += 1)
            end
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_search_operations(
                  operation_id,source_state_id,target_state_id,move_kind,direction,
                  engine_version,candidate_json,outcome,error_json)
                VALUES (?,?,?,?,?,?,?,?,?)
            """, (operation_id, source_state.state_id,
                  target_state_id === nothing ? missing : target_state_id,
                  kind, String(candidate.direction), OGRAPH_SEARCH_ENGINE_VERSION,
                  candidate_json, outcome,
                  error_json === nothing ? missing : error_json))
            MoveGraphStore._dbexec(store.db, """
                INSERT OR IGNORE INTO ograph_search_expansion_operations(
                  expansion_id,ordinal,operation_id) VALUES (?,?,?)
            """, (expansion_id, ordinal, operation_id))
        end
        MoveGraphStore._dbexec(store.db, """
            UPDATE ograph_search_expansions
            SET status='completed',candidate_count=?,success_count=?,reject_count=?,
                error_count=?,error_json=NULL WHERE expansion_id=?
        """, (length(candidates), success_count, reject_count, error_count,
              expansion_id))
    end
    _ograph_search_operations(store, expansion_id)
end

function ograph_search_neighbors(store::MoveGraphStoreDB,
                                 state_id::AbstractString;
                                 kinds=("MP", "02"))
    allowed = unique(String.(collect(kinds)))
    isempty(allowed) && return NamedTuple[]
    all(kind -> kind in ("MP", "02", "PS"), allowed) ||
        throw(ArgumentError("search-neighbor kinds must be MP, 02, or PS"))
    placeholders = join(fill("?", length(allowed)), ",")
    MoveGraphStore._rows(store.db, """
        SELECT * FROM successful_ograph_search_edges
        WHERE source_state_id=? AND move_kind IN ($placeholders)
        ORDER BY move_kind,operation_id
    """, Tuple(Any[String(state_id); allowed]))
end

function _ograph_search_path_rows(store, operation_ids)
    result = NamedTuple[]
    for operation_id in operation_ids
        row = MoveGraphStore._one(store.db, """
            SELECT * FROM successful_ograph_search_edges WHERE operation_id=?
        """, (String(operation_id),))
        row === nothing && throw(ArgumentError("search path operation is absent"))
        push!(result, (operation_id=String(row.operation_id),
                       source_state_id=String(row.source_state_id),
                       target_state_id=String(row.target_state_id),
                       move_kind=String(row.move_kind),
                       direction=String(row.direction),
                       engine_version=String(row.engine_version),
                       candidate=JSON.parse(String(row.candidate_json))))
    end
    result
end

"""Query simple paths, allowing the source to recur only as a final endpoint."""
function query_ograph_search_paths(store::MoveGraphStoreDB,
                                   source_state_id::AbstractString;
                                   target_state_id=nothing,
                                   target_geometry_hash=nothing,
                                   kinds=("MP", "02"),
                                   kind_pattern=nothing,
                                   min_depth::Integer=1,
                                   max_depth::Integer=12,
                                   limit::Integer=1000)
    (target_state_id === nothing) == (target_geometry_hash === nothing) &&
        throw(ArgumentError("specify exactly one path target"))
    min_depth >= 0 && max_depth >= min_depth ||
        throw(ArgumentError("invalid path depth bounds"))
    limit > 0 || throw(ArgumentError("path query limit must be positive"))
    allowed = unique(String.(collect(kinds)))
    pattern = kind_pattern === nothing ? nothing : String.(collect(kind_pattern))
    if pattern !== nothing
        all(kind -> kind in ("MP", "02", "PS"), pattern) ||
            throw(ArgumentError("path pattern contains an unsupported move"))
        max_depth = length(pattern)
        min_depth = length(pattern)
    else
        all(kind -> kind in ("MP", "02", "PS"), allowed) ||
            throw(ArgumentError("path kinds must be MP, 02, or PS"))
        isempty(allowed) && return nothing
    end
    edge_filter = if pattern === nothing
        "edge.move_kind IN (" * join(fill("?", length(allowed)), ",") * ")"
    else
        clauses = ["WHEN $(index - 1) THEN ?" for index in eachindex(pattern)]
        "edge.move_kind = CASE walk.depth " * join(clauses, " ") * " END"
    end
    goal_filter = target_state_id === nothing ?
        "target.geometry_hash=?" : "walk.state_id=?"
    sql = """
        WITH RECURSIVE walk(state_id,depth,operation_ids,visited) AS (
          SELECT ?,0,'','|' || ? || '|'
          UNION ALL
          SELECT edge.target_state_id,walk.depth+1,
                 CASE WHEN walk.operation_ids=''
                      THEN edge.operation_id
                      ELSE walk.operation_ids || ',' || edge.operation_id END,
                 walk.visited || edge.target_state_id || '|'
          FROM walk
          JOIN successful_ograph_search_edges AS edge
            ON edge.source_state_id=walk.state_id
          WHERE walk.depth < ? AND $edge_filter
            AND (walk.depth=0 OR walk.state_id<>?)
            AND (instr(walk.visited,'|' || edge.target_state_id || '|')=0
                 OR (edge.target_state_id=? AND walk.depth+1>=?))
        )
        SELECT walk.operation_ids,walk.depth FROM walk
        JOIN ograph_search_states AS target ON target.state_id=walk.state_id
        WHERE walk.depth>=? AND $goal_filter
        ORDER BY walk.depth,walk.operation_ids LIMIT ?
    """
    params = Any[String(source_state_id), String(source_state_id), max_depth]
    append!(params, pattern === nothing ? allowed : pattern)
    append!(params, [String(source_state_id),String(source_state_id),min_depth])
    push!(params, min_depth)
    push!(params, target_state_id === nothing ? String(target_geometry_hash) :
          String(target_state_id))
    push!(params, Int(limit))
    rows = MoveGraphStore._rows(store.db, sql, Tuple(params))
    [_ograph_search_path_rows(store,
        isempty(String(row.operation_ids)) ? String[] :
        split(String(row.operation_ids), ',')) for row in rows]
end

function query_ograph_search_path(store::MoveGraphStoreDB,
                                  source_state_id::AbstractString; kwargs...)
    paths = query_ograph_search_paths(store, source_state_id; limit=1, kwargs...)
    isempty(paths) ? nothing : first(paths)
end

function _ograph_search_match_candidate(candidate::MoveCandidate, data)
    candidate.type == String(data["type"]) || return false
    candidate.parameter == String(data["parameter"]) || return false
    candidate.vertices == Int.(data["vertices"]) || return false
    [collect(edge) for edge in candidate.edges] ==
        [Int.(edge) for edge in data["edges"]] || return false
    String(candidate.direction) == String(data["direction"]) || return false
    [collect(location) for location in candidate.locations] ==
        [Int.(location) for location in data["locations"]] || return false
    true
end

function _ograph_search_weighted_candidate(link, data)
    kind = String(data["type"])
    parameter = String(data["parameter"])
    vertices = Int.(data["vertices"])
    candidates = if kind == "PS" && String(data["direction"]) == "forward"
        locations = [Tuple(Int.(location)) for location in data["locations"]]
        movable(link, kind; param=parameter, edge_occurrences=locations, strict=true)
    else
        movable(link, kind; param=parameter, vertices=vertices, strict=true)
    end
    matching = filter(candidate -> _ograph_search_match_candidate(candidate, data),
                      candidates)
    length(matching) == 1 ? only(matching) : nothing
end

"""Replay a DB-selected geometric path on weights, inserting H moves only here."""
function replay_weighted_ograph_path(store::MoveGraphStoreDB,
                                     weighted_source::VirtualLink,
                                     path;
                                     h_moves=NamedTuple[])
    weighted_source.is_weighted || throw(ArgumentError(
        "weighted replay requires a weighted source"))
    steps = collect(path)
    source_state = ograph_search_state(weighted_source)
    if !isempty(steps) && String(first(steps).source_state_id) != source_state.state_id
        throw(ArgumentError("weighted source geometry does not match path source"))
    end
    grouped = Dict{Int,Vector{Any}}()
    for hmove in h_moves
        step = Int(hmove.after_step)
        0 <= step <= length(steps) || throw(ArgumentError("H after_step is outside path"))
        push!(get!(grouped, step, Any[]), hmove)
    end
    current = deepcopy(weighted_source)
    states = VirtualLink[deepcopy(current)]
    function apply_h!(step)
        for hmove in get(grouped, step, Any[])
            try
                move!(current, "H"; param=Int(hmove.amount),
                      vertices=[Int(hmove.vertex)], strict=true)
            catch error_value
                return (valid=false, failure_kind="H", failure_step=step,
                        message=sprint(showerror, error_value), target=nothing,
                        states=states)
            end
            push!(states, deepcopy(current))
        end
        nothing
    end
    failure = apply_h!(0)
    failure === nothing || return failure
    for (index, step) in enumerate(steps)
        candidate = _ograph_search_weighted_candidate(current, step.candidate)
        candidate === nothing && return (
            valid=false, failure_kind=String(step.move_kind), failure_step=index,
            message="unweighted path edge is not weight-admissible",
            target=nothing, states=states)
        try
            move!(current, candidate)
        catch error_value
            return (valid=false, failure_kind=String(step.move_kind),
                    failure_step=index, message=sprint(showerror, error_value),
                    target=nothing, states=states)
        end
        push!(states, deepcopy(current))
        failure = apply_h!(index)
        failure === nothing || return failure
    end
    (valid=true, failure_kind=nothing, failure_step=nothing, message=nothing,
     target=current, states=states)
end

function _ograph_search_weight_payload(link::VirtualLink)
    [[isinteger(value) ? Int(round(value)) : value for value in component]
     for component in link.weights]
end

"""Persist a weighted admissibility trial without adding weights to graph state."""
function record_weighted_ograph_path_trial!(store::MoveGraphStoreDB,
                                            run_id::AbstractString,
                                            weighted_source::VirtualLink,
                                            path;
                                            h_moves=NamedTuple[])
    weighted_source.is_weighted || throw(ArgumentError(
        "weighted path trial requires a weighted source"))
    steps = collect(path)
    source_state = put_ograph_search_state!(store, weighted_source)
    operation_ids = String[String(step.operation_id) for step in steps]
    operations_json = _ograph_canonical_json(operation_ids)
    source_weights_json = _ograph_canonical_json(
        _ograph_search_weight_payload(weighted_source))
    h_data = [Dict("after_step" => Int(move.after_step),
                   "vertex" => Int(move.vertex),
                   "amount" => Int(move.amount)) for move in h_moves]
    h_json = _ograph_canonical_json(h_data)
    path_fingerprint = _ograph_sha256(operations_json)
    trial_id = "ograph-search-weighted-trial-" * _ograph_sha256(join((
        String(run_id), path_fingerprint, source_weights_json, h_json), '\0'))
    replay = replay_weighted_ograph_path(store, weighted_source, steps;
                                         h_moves=h_moves)
    outcome = replay.valid ? "admissible" : "rejected"
    target_weights_json = replay.valid ? _ograph_canonical_json(
        _ograph_search_weight_payload(replay.target)) : nothing
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO ograph_search_weighted_trials(
          trial_id,run_id,source_state_id,path_fingerprint,operation_ids_json,
          source_weights_json,h_moves_json,outcome,target_weights_json,
          failure_kind,failure_step,message,observed_count)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,1)
        ON CONFLICT(run_id,path_fingerprint,source_weights_json,h_moves_json)
        DO UPDATE SET observed_count=observed_count+1
    """, (trial_id, String(run_id), source_state.state_id, path_fingerprint,
          operations_json, source_weights_json, h_json, outcome,
          target_weights_json === nothing ? missing : target_weights_json,
          replay.failure_kind === nothing ? missing : String(replay.failure_kind),
          replay.failure_step === nothing ? missing : Int(replay.failure_step),
          replay.message === nothing ? missing : String(replay.message)))
    (trial_id=trial_id, path_fingerprint=path_fingerprint, replay=replay)
end

function verify_ograph_search_store(store::MoveGraphStoreDB)
    version = MoveGraphStore._one(store.db,
        "SELECT value FROM ograph_search_schema_meta WHERE key='schema_version'")
    version !== nothing || throw(ArgumentError("o-graph search schema is absent"))
    parse(Int, String(version.value)) == OGRAPH_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported o-graph search schema version"))
    isempty(MoveGraphStore._rows(store.db, "PRAGMA foreign_key_check")) ||
        throw(ArgumentError("o-graph search foreign-key check failed"))
    count(table) = Int(MoveGraphStore._one(store.db,
        "SELECT COUNT(*) AS n FROM $table").n)
    (schema_version=OGRAPH_SEARCH_SCHEMA_VERSION,
     states=count("ograph_search_states"),
     runs=count("ograph_search_runs"),
     expansions=count("ograph_search_expansions"),
     operations=count("ograph_search_operations"),
     successful_edges=count("successful_ograph_search_edges"),
     weighted_trials=count("ograph_search_weighted_trials"))
end
