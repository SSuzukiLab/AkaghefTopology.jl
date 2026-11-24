# cf: lib/topology.jl/docs/specs/module-structure.md §3.1 — hc_search_* mirrors ograph_search_* on MoveGraphStoreDB
# cf: lib/topology.jl/src/ograph_search_db.jl — run, state, and transition tables over the same store
# cf: lib/topology.jl/src/heegaard/hc_compare_canonical.jl — canonical key, or the literal writing when that key is not used

import JSON
import Dates
import ..MoveGraphStore

const HC_SEARCH_SCHEMA_VERSION = 1
const HC_SEARCH_STATE_VERSION = "hc-search-state/1"
const HC_SEARCH_ENGINE_VERSION = "hc-search-engine/1"
const _HC_SEARCH_IDENTITIES = ("canonical", "literal")
const _HC_SEARCH_TERMINAL = ("completed", "budget-exhausted", "failed")

function _hc_search_schema_path()
    normpath(joinpath(@__DIR__, "..", "..", "db", "hc_search_schema.sql"))
end

function migrate_hc_search_store!(store::MoveGraphStore.MoveGraphStoreDB)
    MoveGraphStore.migrate!(store)
    MoveGraphStore._transaction(store.db) do
        MoveGraphStore._execute_schema!(store.db, read(_hc_search_schema_path(), String))
    end
    _require_hc_search_store(store)
end

function _require_hc_search_store(store::MoveGraphStore.MoveGraphStoreDB)
    present = MoveGraphStore._one(store.db, """
        SELECT 1 AS present FROM sqlite_master
        WHERE type='table' AND name='hc_search_schema_meta'
    """)
    present === nothing && throw(ArgumentError(
        "HC search schema is absent; call migrate_hc_search_store!"))
    row = MoveGraphStore._one(store.db,
        "SELECT value FROM hc_search_schema_meta WHERE key='schema_version'")
    row !== nothing && parse(Int, String(row.value)) == HC_SEARCH_SCHEMA_VERSION ||
        throw(ArgumentError("unsupported HC search schema version"))
    store
end

_hc_identity(identity::AbstractString) = identity in _HC_SEARCH_IDENTITIES ?
    String(identity) : throw(ArgumentError("identity must be canonical or literal"))
_hc_identity(identity::Symbol) = _hc_identity(String(identity))

function _hc_stable_json(value)
    value === nothing && return "null"
    value isa Bool && return value ? "true" : "false"
    if value isa Integer
        integer = BigInt(value)
        typemin(Int64) <= integer <= typemax(Int64) ||
            throw(ArgumentError("search JSON integer does not fit Int64"))
        return string(Int64(integer))
    end
    value isa AbstractString && return JSON.json(String(value))
    value isa Symbol && return JSON.json(String(value))
    value isa Pair && return _hc_stable_json([first(value), last(value)])
    if value isa NamedTuple
        return _hc_stable_json(Dict{String,Any}(String(key) => getfield(value, key) for key in keys(value)))
    end
    if value isa AbstractDict
        parts = String[]
        for key in sort!(collect(keys(value)); by=string)
            push!(parts, JSON.json(string(key)) * ":" * _hc_stable_json(value[key]))
        end
        return "{" * join(parts, ",") * "}"
    end
    if value isa Union{AbstractVector,Tuple}
        return "[" * join((_hc_stable_json(item) for item in value), ",") * "]"
    end
    throw(ArgumentError("not a search JSON value: $(typeof(value))"))
end

function _hc_pair_list(dict)
    [[string(key), dict[key]] for key in sort!(collect(keys(dict)); by=string)]
end

function _hc_writing(hc::HeegaardCode)
    lc = hc.lc
    order = hc_component_names(lc)
    Dict{String,Any}(
        "schema" => "hc-search-writing/1",
        "genus" => lc.genus,
        "words" => [Dict{String,Any}("name" => name, "word" => copy(lc.words[name]))
                    for name in order if haskey(lc.words, name)],
        "places" => _hc_pair_list(lc.places),
        "lc_disks" => _hc_pair_list(lc.disks),
        "lc_points" => _hc_pair_list(lc.points),
        "lc_weights" => _hc_pair_list(lc.weights),
        "roles" => [Dict{String,Any}("name" => name, "role" => String(hc.roles[name]))
                    for name in order if haskey(hc.roles, name)],
        "disks" => _hc_pair_list(hc.disks),
        "points" => _hc_pair_list(hc.points),
        "weights" => _hc_pair_list(hc.weights),
    )
end

function _hc_disk_dict(items)
    Dict{String,Vector{String}}(String(item[1]) => String.(item[2]) for item in items)
end

function _hc_string_dict(items)
    Dict{String,String}(String(item[1]) => String(item[2]) for item in items)
end

function _hc_int_dict(items)
    Dict{String,Int}(String(item[1]) => Int(item[2]) for item in items)
end

function _hc_from_writing(data)
    String(data["schema"]) == "hc-search-writing/1" ||
        throw(ArgumentError("unknown HC search writing"))
    words = Dict{String,Vector{Int}}(
        String(item["name"]) => Int.(item["word"]) for item in data["words"])
    lc = LickorishCode(
        Int(data["genus"]), words,
        _hc_string_dict(data["places"]),
        _hc_disk_dict(data["lc_disks"]),
        _hc_string_dict(data["lc_points"]),
        _hc_int_dict(data["lc_weights"]),
        nothing)
    roles = Dict{String,Symbol}(
        String(item["name"]) => Symbol(String(item["role"])) for item in data["roles"])
    hc = HeegaardCode(lc, roles,
        _hc_disk_dict(data["disks"]),
        _hc_string_dict(data["points"]),
        _hc_int_dict(data["weights"]),
        nothing)
    validate_hc_input(hc)
    hc
end

function _hc_form_data(form::HCCanonical)
    Dict{String,Any}(
        "genus" => form.genus,
        "k_words" => form.k_words,
        "alpha_words" => form.alpha_words,
        "beta_words" => form.beta_words,
        "other_words" => form.other_words,
        "places" => form.places,
        "lc_disks" => form.lc_disks,
        "hc_disks" => form.hc_disks,
        "lc_points" => form.lc_points,
        "hc_points" => form.hc_points,
        "lc_weights" => form.lc_weights,
        "hc_weights" => form.hc_weights,
    )
end

"""探索の鍵。`identity=:canonical` は `hc_canonical` の列、`:literal` は表記そのものである。"""
function hc_search_state(hc::HeegaardCode; identity=:canonical)
    mode = _hc_identity(identity)
    writing = _hc_writing(hc)
    body = if mode == "literal"
        Dict{String,Any}("mode" => "literal", "writing" => writing)
    else
        Dict{String,Any}("mode" => "canonical", "form" => _hc_form_data(hc_canonical(hc)))
    end
    MoveGraphStore.CanonicalState(_hc_stable_json(body); area=MoveGraphStore.AREA_3,
        representation="hc-search", payload=_hc_stable_json(writing),
        canonical_version=HC_SEARCH_SCHEMA_VERSION,
        complexity=Float64(hc.lc.genus))
end

function _store_hc_state!(store, state, hc, mode)
    MoveGraphStore.put_state!(store, state)
    MoveGraphStore._dbexec(store.db, """
        INSERT OR IGNORE INTO hc_search_states(
          state_id, state_version, identity_mode, raw_json)
        VALUES (?,?,?,?)
    """, (state.state_id, HC_SEARCH_STATE_VERSION, mode, _hc_stable_json(_hc_writing(hc))))
    state
end

function put_hc_search_state!(store::MoveGraphStore.MoveGraphStoreDB, hc::HeegaardCode; identity=:canonical)
    _require_hc_search_store(store)
    mode = _hc_identity(identity)
    _store_hc_state!(store, hc_search_state(hc; identity=mode), hc, mode)
end

function load_hc_search_state(store::MoveGraphStore.MoveGraphStoreDB, state_id::AbstractString)
    _require_hc_search_store(store)
    row = MoveGraphStore._one(store.db, """
        SELECT raw_json FROM hc_search_states WHERE state_id=?
    """, (String(state_id),))
    row === nothing && return nothing
    _hc_from_writing(JSON.parse(String(row.raw_json)))
end

function _hc_script_snapshot(script_path::AbstractString)
    path = abspath(String(script_path))
    isfile(path) || throw(ArgumentError("provenance script does not exist: $path"))
    path, bytes2hex(sha256(read(path)))
end

function begin_hc_search_run!(store::MoveGraphStore.MoveGraphStoreDB, script_path::AbstractString;
                              label::AbstractString, scope=Dict{String,Any}())
    _require_hc_search_store(store)
    path, digest = _hc_script_snapshot(script_path)
    started = string(Dates.now(Dates.UTC))
    run_id = "hc-search-run-" * bytes2hex(sha256(join(
        (replace(path, '\\' => '/'), digest, String(label), started), '\0')))
    MoveGraphStore._dbexec(store.db, """
        INSERT INTO hc_search_runs(
          run_id, script_path, script_sha256, label, scope_json, status, counts_json,
          started_at_utc)
        VALUES (?,?,?,?,?,'running','{}',?)
    """, (run_id, replace(path, '\\' => '/'), digest, String(label),
          _hc_stable_json(scope), started))
    run_id
end

function hc_search_run(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString)
    _require_hc_search_store(store)
    row = MoveGraphStore._one(store.db, """
        SELECT run.*,
          (SELECT COUNT(*) FROM hc_search_operations WHERE run_id=run.run_id) AS operations,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='success') AS successes,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='rejected') AS rejects,
          (SELECT COUNT(*) FROM hc_search_operations
            WHERE run_id=run.run_id AND outcome='error') AS errors
        FROM hc_search_runs AS run WHERE run.run_id=?
    """, (String(run_id),))
    row === nothing && return nothing
    (run_id=String(row.run_id), script_path=String(row.script_path),
     script_sha256=String(row.script_sha256), label=String(row.label),
     scope=JSON.parse(String(row.scope_json)), status=String(row.status),
     counts=JSON.parse(String(row.counts_json)),
     operations=Int(row.operations), successes=Int(row.successes),
     rejects=Int(row.rejects), errors=Int(row.errors))
end

function finish_hc_search_run!(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString;
                               status::AbstractString, counts=Dict{String,Any}())
    status in _HC_SEARCH_TERMINAL || throw(ArgumentError("invalid terminal HC search status"))
    MoveGraphStore._dbexec(store.db, """
        UPDATE hc_search_runs
        SET status=?, counts_json=?, finished_at_utc=?
        WHERE run_id=? AND status='running'
    """, (String(status), _hc_stable_json(counts), string(Dates.now(Dates.UTC)), String(run_id)))
    hc_search_run(store, run_id)
end

function _hc_jsonable(value)
    value === nothing && return nothing
    value isa Bool && return value
    if value isa Integer
        return typemin(Int) <= value <= typemax(Int) ? Int(value) : string(value)
    end
    value isa AbstractString && return String(value)
    value isa Symbol && return String(value)
    value isa Rational && return string(value)
    value isa Pair && return [_hc_jsonable(first(value)), _hc_jsonable(last(value))]
    value isa NamedTuple && return Dict{String,Any}(String(key) => _hc_jsonable(getfield(value, key)) for key in keys(value))
    value isa AbstractDict && return Dict{String,Any}(string(key) => _hc_jsonable(item) for (key, item) in value)
    value isa Union{AbstractVector,Tuple} && return [_hc_jsonable(item) for item in value]
    string(value)
end

function _hc_support_brief(value)
    value isa AbstractDict || return _hc_jsonable(value)
    brief = Dict{String,Any}()
    for key in ("kind", "core", "curve", "side")
        haskey(value, key) || continue
        item = value[key]
        brief[key] = item isa AbstractDict ? _hc_support_brief(item) : _hc_jsonable(item)
    end
    brief
end

function _hc_step_brief(step::LCMoveWitness)
    data = Dict{String,Any}(
        "row_id" => step.row_id,
        "locations" => [collect(Int, location) for location in step.locations],
        "before" => step.before_frags,
        "after" => step.after_frags,
    )
    for key in ("q", "curve", "side", "acted", "acted_components", "disk")
        haskey(step.extra, key) || continue
        data[key] = _hc_jsonable(step.extra[key])
    end
    haskey(step.extra, "support") && (data["support"] = _hc_support_brief(step.extra["support"]))
    data
end

function _hc_candidate_data(witness::LCMoveWitness)
    data = Dict{String,Any}(
        "row_id" => witness.row_id,
        "names" => witness.names,
        "locations" => [collect(Int, location) for location in witness.locations],
        "vertices" => witness.vertices,
        "before" => witness.before_frags,
        "after" => witness.after_frags,
        "genus_delta" => witness.genus_delta,
    )
    for key in ("q", "support", "acted", "acted_components", "disk", "curve", "side",
                "k_identification", "name", "target", "curves", "band")
        haskey(witness.extra, key) || continue
        value = witness.extra[key]
        data[key] = key == "support" ? _hc_support_brief(value) : _hc_jsonable(value)
    end
    # Prepared twists and stabilizations keep the choice inside the sequence, not on the outer witness.
    if haskey(witness.extra, "preparation")
        data["preparation"] = _hc_jsonable(witness.extra["preparation"])
    end
    if haskey(witness.extra, "sequence")
        steps = witness.extra["sequence"]
        data["sequence_length"] = length(steps)
        data["sequence"] = [_hc_step_brief(step) for step in steps if step isa LCMoveWitness]
    end
    if haskey(witness.extra, "created_components")
        data["created_components"] = _hc_jsonable(witness.extra["created_components"])
    end
    attachment = get(witness.extra, "attachment", nothing)
    if attachment isa AbstractDict
        brief = Dict{String,Any}()
        for key in ("requested_disk", "actual_face", "reassigned")
            haskey(attachment, key) && (brief[key] = _hc_jsonable(attachment[key]))
        end
        data["attachment"] = brief
    end
    data
end

_hc_candidate_json(witness::LCMoveWitness) = _hc_stable_json(_hc_candidate_data(witness))

"""`movable` が返す候補を、選択の安定な JSON で重複なく並べる。"""
function hc_search_candidates(hc::HeegaardCode)
    seen = Set{String}()
    result = LCMoveWitness[]
    for witness in movable(hc)
        key = _hc_candidate_json(witness)
        key in seen && continue
        push!(seen, key)
        push!(result, witness)
    end
    result
end

function _hc_err_json(err)
    message = sprint(showerror, err)
    length(message) > 800 && (message = first(message, 800))
    _hc_stable_json(Dict{String,Any}("type" => string(typeof(err)), "message" => message))
end

function _hc_known_state(store, state_id)
    MoveGraphStore._one(store.db,
        "SELECT 1 AS present FROM hc_search_states WHERE state_id=?", (state_id,)) !== nothing
end

function _hc_bundle_text(hc::HeegaardCode)
    try
        lc_text(hc.lc) * "\n--- hc ---\n" * hc_text(hc; lc_filename="diagram.lc", encoding=:explicit)
    catch err
        "text-failed: " * sprint(showerror, err)
    end
end

function _hc_snapshot(hc::HeegaardCode)
    homology = try
        value = hc_first_homology(hc)
        Dict{String,Any}(
            "ok" => true,
            "free_rank" => value.free_rank,
            "torsion" => copy(value.torsion),
            "invariant_factors" => copy(value.invariant_factors),
        )
    catch err
        Dict{String,Any}("ok" => false, "error" => sprint(showerror, err))
    end
    linking = try
        form = hc_linking_form(hc)
        rows, cols = size(form.matrix)
        entries = String[]
        for col in 1:cols, row in 1:rows
            push!(entries, string(form.matrix[row, col]))
        end
        Dict{String,Any}(
            "ok" => true,
            "free_rank" => form.free_rank,
            "torsion" => [Int(order) for order in form.torsion],
            "matrix" => entries,
        )
    catch err
        Dict{String,Any}("ok" => false, "error" => sprint(showerror, err))
    end
    (homology=homology, linking=linking)
end

function hc_search_same_homology(left::AbstractDict, right::AbstractDict)
    # `invariant_factors` in this code keeps the Smith 1s, so it changes with the presentation.
    # The group is the free rank together with the torsion coefficients above 1.
    get(left, "ok", false) && get(right, "ok", false) || return false
    left["free_rank"] == right["free_rank"] && left["torsion"] == right["torsion"]
end

function _hc_rational_matrix(entries)
    n = isqrt(length(entries))
    n * n == length(entries) || throw(ArgumentError("linking matrix is not square"))
    values = [parse(Rational{BigInt}, String(entry)) for entry in entries]
    reshape(values, n, n)
end

"""出発点の絡み形式と同じか。巡回成分は向きを保つ自己同型、それ以外は行列そのものを比べる。"""
function hc_search_same_linking(left::AbstractDict, right::AbstractDict)
    get(left, "ok", false) && get(right, "ok", false) || return false
    left["free_rank"] == right["free_rank"] || return false
    torsion = Int.(left["torsion"])
    torsion == Int.(right["torsion"]) || return false
    isempty(torsion) && return true
    left_matrix = _hc_rational_matrix(left["matrix"])
    right_matrix = _hc_rational_matrix(right["matrix"])
    if length(torsion) == 1
        order = torsion[1]
        order > 10_000 && return left_matrix == right_matrix
        return any(unit -> gcd(unit, order) == 1 &&
            isinteger(left_matrix[1, 1] - right_matrix[1, 1] * unit * unit), 1:(order - 1))
    end
    left_matrix == right_matrix
end

function _hc_note!(problems, kind, row_id, depth, detail, on_problem)
    key = (kind, row_id)
    if !haskey(problems, key) || depth < problems[key]["depth"]
        record = Dict{String,Any}(
            "kind" => kind,
            "row_id" => row_id,
            "depth" => depth,
        )
        for (name, value) in detail
            record[String(name)] = value
        end
        problems[key] = record
        on_problem === nothing || on_problem(record)
    end
end

function _hc_invariant_notes!(problems, start_snap, snap, row_id, depth, detail, on_problem)
    if !snap.homology["ok"]
        _hc_note!(problems, "homology_error", row_id, depth,
            merge(detail, Dict("value" => snap.homology)), on_problem)
    elseif start_snap.homology["ok"] && !hc_search_same_homology(start_snap.homology, snap.homology)
        _hc_note!(problems, "homology", row_id, depth,
            merge(detail, Dict("value" => snap.homology, "start" => start_snap.homology)), on_problem)
    end
    if !snap.linking["ok"]
        _hc_note!(problems, "linking_error", row_id, depth,
            merge(detail, Dict("value" => snap.linking)), on_problem)
    elseif start_snap.linking["ok"] && !hc_search_same_linking(start_snap.linking, snap.linking)
        _hc_note!(problems, "linking", row_id, depth,
            merge(detail, Dict("value" => snap.linking, "start" => start_snap.linking)), on_problem)
    end
end

"""
一状態の候補をすべて適用し、遷移を保存する。

`depth` は出発状態までの辺の数。新しい状態をあと `room` 個しか置けないとき、
新しい状態を作る候補に達した時点で、その候補以降は適用しない。
"""
function expand_hc_search_state!(store::MoveGraphStore.MoveGraphStoreDB, run_id::AbstractString,
                                 source::HeegaardCode;
                                 depth::Integer=0, identity=:canonical, room::Integer=typemax(Int),
                                 check_inverse::Bool=true, check_invariants::Bool=false,
                                 start_snap=nothing, problems=nothing, on_problem=nothing)
    _require_hc_search_store(store)
    mode = _hc_identity(identity)
    source_depth = Int(depth)
    source_state = put_hc_search_state!(store, source; identity=mode)
    records = NamedTuple[]
    created = 0
    stopped = false
    candidates = hc_search_candidates(source)
    for witness in candidates
        candidate = _hc_candidate_data(witness)
        candidate_json = _hc_stable_json(candidate)
        outcome = "error"
        inverse_status = "unchecked"
        target = nothing
        target_state_id = nothing
        error_json = nothing
        local result
        try
            result = move(source, witness)
            target = result.output
            target isa HeegaardCode || throw(LCMoveError(:expansion, "move did not return an HC",
                Dict("row" => witness.row_id, "type" => string(typeof(target)))))
            target_state = hc_search_state(target; identity=mode)
            known = _hc_known_state(store, target_state.state_id)
            if !known && room - created <= 0
                stopped = true
                target = nothing
                break
            end
            _store_hc_state!(store, target_state, target, mode)
            known || (created += 1)
            target_state_id = target_state.state_id
            outcome = "success"
            if check_inverse
                try
                    back = move(target, result.inverse)
                    back.output isa HeegaardCode || error("inverse did not return an HC")
                    back_state = hc_search_state(back.output; identity=mode)
                    inverse_status = back_state.state_id == source_state.state_id ? "ok" : "mismatch"
                    if inverse_status == "mismatch" && problems !== nothing
                        _hc_note!(problems, "inverse_mismatch", witness.row_id, source_depth + 1,
                            Dict("source_text" => _hc_bundle_text(source),
                                 "target_text" => _hc_bundle_text(target),
                                 "back_text" => _hc_bundle_text(back.output),
                                 "witness" => candidate,
                                 "source_state_id" => source_state.state_id,
                                 "back_state_id" => back_state.state_id), on_problem)
                    end
                catch err
                    inverse_status = "error"
                    error_json = _hc_err_json(err)
                    problems === nothing || _hc_note!(problems,
                        err isa HCNotComputed ? "key" : "inverse_error",
                        witness.row_id, source_depth + 1,
                        Dict("source_text" => _hc_bundle_text(source),
                             "target_text" => _hc_bundle_text(target),
                             "witness" => candidate,
                             "value" => JSON.parse(error_json)), on_problem)
                end
            end
            if check_invariants && start_snap !== nothing && problems !== nothing && !known
                snap = _hc_snapshot(target)
                _hc_invariant_notes!(problems, start_snap, snap, witness.row_id, source_depth + 1,
                    Dict("source_text" => _hc_bundle_text(source),
                         "target_text" => _hc_bundle_text(target),
                         "witness" => candidate), on_problem)
            end
        catch err
            err isa InterruptException && rethrow()
            target = nothing
            target_state_id = nothing
            outcome = err isa LCMoveError ? "rejected" : "error"
            inverse_status = "unchecked"
            error_json = _hc_err_json(err)
            kind = if err isa HCNotComputed
                "key"
            elseif err isa LCMoveError
                "move_rejected"
            else
                "move_error"
            end
            problems === nothing || _hc_note!(problems, kind, witness.row_id, source_depth + 1,
                Dict("source_text" => _hc_bundle_text(source),
                     "witness" => candidate,
                     "value" => JSON.parse(error_json)), on_problem)
        end
        operation_id = "hc-op-" * bytes2hex(sha256(join(
            (String(run_id), source_state.state_id, HC_SEARCH_ENGINE_VERSION, candidate_json), '\0')))
        MoveGraphStore._dbexec(store.db, """
            INSERT OR IGNORE INTO hc_search_operations(
              operation_id, run_id, source_state_id, target_state_id, depth, row_id,
              engine_version, candidate_json, outcome, inverse_status, error_json)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, (operation_id, String(run_id), source_state.state_id,
              target_state_id === nothing ? missing : target_state_id,
              source_depth + 1, witness.row_id, HC_SEARCH_ENGINE_VERSION, candidate_json,
              outcome, inverse_status, error_json === nothing ? missing : error_json))
        push!(records, (
            operation_id=operation_id, row_id=witness.row_id, outcome=outcome,
            inverse_status=inverse_status, target_state_id=target_state_id,
            target=target, depth=source_depth + 1,
        ))
        stopped && break
    end
    (records=records, created=created, stopped=stopped, source_state_id=source_state.state_id,
     candidate_count=length(candidates))
end

function hc_search_neighbors(store::MoveGraphStore.MoveGraphStoreDB, state_id::AbstractString; row_ids=nothing)
    _require_hc_search_store(store)
    if row_ids === nothing
        rows = MoveGraphStore._rows(store.db, """
            SELECT * FROM successful_hc_search_edges
            WHERE source_state_id=? ORDER BY depth, row_id, operation_id
        """, (String(state_id),))
    else
        allowed = unique(String.(collect(row_ids)))
        isempty(allowed) && return NamedTuple[]
        marks = join(fill("?", length(allowed)), ",")
        rows = MoveGraphStore._rows(store.db, """
            SELECT * FROM successful_hc_search_edges
            WHERE source_state_id=? AND row_id IN ($marks)
            ORDER BY depth, row_id, operation_id
        """, Tuple(Any[String(state_id); allowed]))
    end
    [(operation_id=String(row.operation_id), run_id=String(row.run_id),
      source_state_id=String(row.source_state_id), target_state_id=String(row.target_state_id),
      depth=Int(row.depth), row_id=String(row.row_id),
      inverse_status=String(row.inverse_status),
      candidate=JSON.parse(String(row.candidate_json))) for row in rows]
end

function query_hc_search_paths(store::MoveGraphStore.MoveGraphStoreDB, source_state_id::AbstractString;
                               target_state_id::AbstractString, min_depth::Integer=1,
                               max_depth::Integer=12, limit::Integer=100, row_ids=nothing)
    _require_hc_search_store(store)
    min_depth >= 0 && max_depth >= min_depth || throw(ArgumentError("invalid path depth bounds"))
    limit > 0 || throw(ArgumentError("path query limit must be positive"))
    allowed = row_ids === nothing ? nothing : unique(String.(collect(row_ids)))
    allowed !== nothing && isempty(allowed) && return NamedTuple[]
    edge_filter = allowed === nothing ? "1=1" :
        "edge.row_id IN (" * join(fill("?", length(allowed)), ",") * ")"
    sql = """
        WITH RECURSIVE walk(state_id, depth, operation_ids, visited) AS (
          SELECT ?, 0, '', '|' || ? || '|'
          UNION ALL
          SELECT edge.target_state_id, walk.depth + 1,
                 CASE WHEN walk.operation_ids='' THEN edge.operation_id
                      ELSE walk.operation_ids || ',' || edge.operation_id END,
                 walk.visited || edge.target_state_id || '|'
          FROM walk
          JOIN successful_hc_search_edges AS edge ON edge.source_state_id=walk.state_id
          WHERE walk.depth < ? AND $edge_filter
            AND instr(walk.visited, '|' || edge.target_state_id || '|')=0
        )
        SELECT operation_ids, depth FROM walk
        WHERE state_id=? AND depth>=? AND depth<=?
        ORDER BY depth, operation_ids LIMIT ?
    """
    params = Any[String(source_state_id), String(source_state_id), Int(max_depth)]
    allowed === nothing || append!(params, allowed)
    append!(params, [String(target_state_id), Int(min_depth), Int(max_depth), Int(limit)])
    rows = MoveGraphStore._rows(store.db, sql, Tuple(params))
    paths = NamedTuple[]
    for row in rows
        ids = isempty(String(row.operation_ids)) ? String[] : split(String(row.operation_ids), ',')
        steps = NamedTuple[]
        for operation_id in ids
            step = MoveGraphStore._one(store.db, """
                SELECT * FROM successful_hc_search_edges WHERE operation_id=?
            """, (String(operation_id),))
            step === nothing && throw(ArgumentError("search path operation is absent"))
            push!(steps, (
                operation_id=String(step.operation_id),
                source_state_id=String(step.source_state_id),
                target_state_id=String(step.target_state_id),
                depth=Int(step.depth), row_id=String(step.row_id),
                inverse_status=String(step.inverse_status),
                candidate=JSON.parse(String(step.candidate_json)),
            ))
        end
        push!(paths, (depth=Int(row.depth), steps=steps))
    end
    paths
end

"""
深さ `max_depth` 未満の状態だけを広げ、状態数が `max_states` に達したら止める。

深さは辺の数である。出発点は深さ 0 で、状態数に含まれる。
`identity` はこの探索の間は変えない。正規形が途中で予算を超えたときは、その遷移を鍵の失敗として残す。
"""
function hc_search_explore(root::HeegaardCode, store::MoveGraphStore.MoveGraphStoreDB,
                           script_path::AbstractString;
                           label::AbstractString, identity=:canonical,
                           max_depth::Integer=3, max_states::Integer=2000,
                           check_invariants::Bool=true, progress::Bool=false,
                           on_problem=nothing)
    mode = _hc_identity(identity)
    max_depth >= 0 && max_states >= 1 || throw(ArgumentError("invalid search bounds"))
    migrate_hc_search_store!(store)
    scope = Dict{String,Any}(
        "identity" => mode, "max_depth" => Int(max_depth), "max_states" => Int(max_states),
        "label" => String(label),
    )
    run_id = begin_hc_search_run!(store, script_path; label=String(label), scope=scope)
    problems = Dict{Tuple{String,String},Dict{String,Any}}()
    status = "completed"
    counts = Dict{String,Any}()
    try
        root_state = put_hc_search_state!(store, root; identity=mode)
        start_snap = check_invariants ? _hc_snapshot(root) : nothing
        if start_snap !== nothing
            _hc_invariant_notes!(problems, start_snap, start_snap, "(start)", 0,
                Dict("source_text" => _hc_bundle_text(root), "witness" => nothing), on_problem)
        end
        queue = [(root, 0)]
        depths = Dict(root_state.state_id => 0)
        created_total = 1
        transitions = 0
        while !isempty(queue) && status == "completed"
            hc, depth = popfirst!(queue)
            depth >= max_depth && continue
            if progress
                println("hc-search label=$label depth=$depth states=$created_total queue=$(length(queue)) transitions=$transitions")
                flush(stdout)
            end
            room = max_states - created_total
            expansion = expand_hc_search_state!(store, run_id, hc;
                depth=depth, identity=mode, room=room, check_inverse=true,
                check_invariants=check_invariants, start_snap=start_snap,
                problems=problems, on_problem=on_problem)
            transitions += length(expansion.records)
            for record in expansion.records
                record.outcome == "success" || continue
                record.target_state_id in keys(depths) && continue
                depths[record.target_state_id] = record.depth
                created_total += 1
                record.depth < max_depth && push!(queue, (record.target, record.depth))
            end
            expansion.stopped && (status = "budget-exhausted")
        end
        counts = Dict{String,Any}(
            "states" => created_total,
            "transitions" => transitions,
            "max_depth_reached" => isempty(depths) ? 0 : maximum(values(depths)),
            "identity" => mode,
            "problems" => length(problems),
        )
        finished = finish_hc_search_run!(store, run_id; status=status, counts=counts)
        return (
            run_id=run_id, status=status, identity=mode, counts=counts, run=finished,
            states=created_total, transitions=transitions,
            depth=counts["max_depth_reached"],
            start=start_snap, problems=collect(values(problems)),
        )
    catch err
        err isa InterruptException && rethrow()
        counts["error"] = sprint(showerror, err)
        try
            finish_hc_search_run!(store, run_id; status="failed", counts=counts)
        catch
        end
        rethrow()
    end
end
