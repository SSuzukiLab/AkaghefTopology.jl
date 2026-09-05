"""Build the full canonical signed-Gauss N<=6 classification database.

This is the full-census counterpart of the N6 classical-subset fallback.  The
all-state lane is deliberately streamed: the only retained per-state lookup is
the classical canonical serialization => final label map (or `nothing` for an
unresolved classical component).  Every emitted canonical state gets one
`states` row and one five-column `state_classifications` row.  The remaining
schema tables are left empty.

The production entry point requires the final barrier-8 result explicitly.  A
barrier-8 object is not fabricated from the barrier-7 result.  The CLI obtains
the three production inputs through `_n6full_production_inputs`, validates the
barrier-8 gate, and only then starts the full streamed database build.
"""

using JSON
using SHA
using SQLite
using DBInterface

const N6FULL_REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const N6FULL_PROJECT_ROOT = joinpath(N6FULL_REPO_ROOT, "threads", "PJ260830")
const N6FULL_SNAPSHOT_PATH = joinpath(N6FULL_PROJECT_ROOT, "data",
                                      "knotinfo-n6.json")
const N6FULL_DB_OUTPUT = joinpath(N6FULL_REPO_ROOT, "references", "items",
                                  "database", "knotinfo",
                                  "knot-n6-full-classification.sqlite3")
const N6FULL_JSON_OUTPUT = joinpath(N6FULL_PROJECT_ROOT, "data",
                                    "knot-n6-full-classification.json")
const N6FULL_MD_OUTPUT = joinpath(N6FULL_PROJECT_ROOT, "reports",
                                  "knot-n6-full-classification.md")

const N6FULL_METHOD_VERSION = "knot-n6-full-classification-v1"
# Captured once from the operating system while this method release was
# prepared; it is intentionally not regenerated on every build.
const N6FULL_METHOD_RELEASED_AT_UTC = "2026-08-31T16:43:02Z"
const N6FULL_MAX_N = 6
const N6FULL_BARRIER6 = 6
const N6FULL_BARRIER7 = 7
const N6FULL_BARRIER8 = 8
const N6FULL_BARRIER6_FINGERPRINT =
    "ae58e6b0c9e4c155d20b415cff0904f7babcd7e1b70a86f40d960c2e048170fb"
const N6FULL_BARRIER7_FINGERPRINT =
    "204a5b1f6d4f1df4aa9907c6f2b171c57f4903c070c179a427bfff50c908f839"
const N6FULL_BARRIER8_FINGERPRINT =
    "ceb102736747324058b0b538e7edf7a56a59cc2d85272f59a839b83b23c9986c"
const N6FULL_BARRIER8_STATUS = "barrier8_unresolved"
const N6FULL_BARRIER8_UNRESOLVED_COMPONENTS = collect(9:18)
const N6FULL_CANONICAL_STATE_COUNTS = [1, 2, 14, 164, 3388, 96776, 3548876]
const N6FULL_RAW_STATE_COUNTS = [1, 4, 48, 960, 26880, 967680, 42577920]
const N6FULL_FULL_CANONICAL_TOTAL = 3649221
const N6FULL_CLASSICAL_CANONICAL_TOTAL = 77989
const N6FULL_NONCLASSICAL_CANONICAL_TOTAL = 3571232
const N6FULL_TIME_LIMIT_SECONDS = 5 * 60 * 60
const N6FULL_RSS_LIMIT_BYTES = 22 * 1024^3
const N6FULL_FREE_DISK_LIMIT_BYTES = 10 * 1024^3
const N6FULL_PROGRESS_INTERVAL = 100_000
const N6FULL_HASH_CHUNK_BYTES = 1024^2

if !isdefined(@__MODULE__, :MoveGraphStoreDB)
    include(joinpath(@__DIR__, "..", "src", "move_graph_store.jl"))
end
if !isdefined(@__MODULE__, :KnotN6Classification)
    include(joinpath(@__DIR__, "..", "src", "knot_n6_classification.jl"))
end
if !isdefined(@__MODULE__, :KnotN6TargetedBridgeResult)
    include(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_bridge.jl"))
end
if !isdefined(@__MODULE__, :KnotN6TargetedBarrier8Result) &&
   isfile(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_barrier8.jl"))
    include(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_barrier8.jl"))
end
if !isdefined(@__MODULE__, :KnotGaussCanonicalStreamStats)
    include(joinpath(@__DIR__, "..", "src", "knot_n6_all_state_enumeration.jl"))
end
using .MoveGraphStore

function _n6full_sha256_file(path::AbstractString)
    context = SHA.SHA256_CTX()
    buffer = Vector{UInt8}(undef, N6FULL_HASH_CHUNK_BYTES)
    total = 0
    open(path, "r") do io
        while !eof(io)
            count = readbytes!(io, buffer)
            count == 0 && break
            SHA.update!(context, view(buffer, 1:count))
            total += count
        end
    end
    bytes2hex(SHA.digest!(context)), total
end

function _n6full_sha256_lines(lines)
    context = SHA.SHA256_CTX()
    for line in lines
        SHA.update!(context, codeunits(string(line, '\n')))
    end
    bytes2hex(SHA.digest!(context))
end

function _n6full_getproperty(value, name::Symbol, default=nothing)
    value === nothing && return default
    if value isa AbstractDict
        haskey(value, name) && return value[name]
        text = String(name)
        haskey(value, text) && return value[text]
        return default
    end
    hasproperty(value, name) ? getproperty(value, name) : default
end

function _n6full_component_key(value)
    value isa Integer && return Int(value)
    parse(Int, String(value))
end

function _n6full_label_name(label::AbstractString)
    text = strip(String(label))
    isempty(text) && throw(ArgumentError("empty component label"))
    first(split(text, '^'; limit=2))
end

function _n6full_label_entries(value, field::AbstractString)
    value === nothing && return Dict{Int,String}()
    entries = Dict{Int,String}()
    if value isa AbstractDict
        for (component, label) in value
            label === nothing && continue
            entries[_n6full_component_key(component)] = String(label)
        end
    elseif value isa AbstractVector
        for item in value
            if item isa Pair
                component, label = item
                label === nothing && continue
                entries[_n6full_component_key(component)] = String(label)
            elseif hasproperty(item, :component) && hasproperty(item, :label)
                label = getproperty(item, :label)
                label === nothing && continue
                entries[_n6full_component_key(getproperty(item, :component))] =
                    String(label)
            else
                throw(ArgumentError("unsupported $field entry representation"))
            end
        end
    else
        throw(ArgumentError("unsupported $field representation"))
    end
    for (component, label) in entries
        isempty(strip(label)) && throw(ArgumentError(
            "empty $field label for component $component"))
    end
    entries
end

"""Extract final component labels from either a bridge result or a map.

The barrier-8 implementation is intentionally allowed to introduce a compact
result type later.  This adapter accepts the existing `resolved_components`
field, a `component_type_labels` field, or a plain component=>label map, while
rejecting unknown shapes instead of silently classifying states.
"""
function _n6full_bridge_component_labels(bridge, field::AbstractString)
    bridge === nothing && return Dict{Int,String}()
    if bridge isa AbstractDict
        return _n6full_label_entries(bridge, field)
    end
    resolved = _n6full_getproperty(bridge, :resolved_components, nothing)
    resolved !== nothing && return _n6full_label_entries(resolved,
                                                         "$field.resolved_components")
    labels = _n6full_getproperty(bridge, :component_type_labels, nothing)
    labels !== nothing && return _n6full_label_entries(labels,
                                                       "$field.component_type_labels")
    throw(ArgumentError("$field result has no recognized component-label field"))
end

function _n6full_connection_labels(bridge)
    raw = _n6full_getproperty(bridge, :connected_known_labels, nothing)
    raw === nothing && return Dict{Int,Vector{String}}()
    result = Dict{Int,Vector{String}}()
    if !(raw isa AbstractDict)
        throw(ArgumentError("connected_known_labels must be a dictionary"))
    end
    for (component, labels) in raw
        labels isa AbstractVector || throw(ArgumentError(
            "connected_known_labels values must be vectors"))
        result[_n6full_component_key(component)] = String.(labels)
    end
    result
end

function _n6full_validate_barrier8_inputs(
        classification::KnotN6Classification, bridge7, barrier8)
    barrier8 === nothing && return nothing
    barrier = _n6full_getproperty(barrier8, :barrier, nothing)
    barrier === nothing || Int(barrier) == N6FULL_BARRIER8 || throw(ErrorException(
        "barrier8 result has the wrong barrier value"))
    n6_fingerprint = _n6full_getproperty(barrier8, :n6_input_fingerprint, nothing)
    n6_fingerprint === nothing || String(n6_fingerprint) == classification.fingerprint ||
        throw(ErrorException("barrier8 result does not match barrier-6 classification"))
    bridge7_fingerprint = _n6full_getproperty(barrier8, :barrier7_fingerprint, nothing)
    expected_bridge7 = _n6full_getproperty(bridge7, :fingerprint, nothing)
    bridge7_fingerprint === nothing || expected_bridge7 === nothing ||
        String(bridge7_fingerprint) == String(expected_bridge7) ||
        throw(ErrorException("barrier8 result does not match barrier-7 bridge"))
    seeds8 = _n6full_getproperty(barrier8, :seed_components, nothing)
    seeds7 = _n6full_getproperty(bridge7, :unresolved_components, nothing)
    if seeds8 !== nothing && seeds7 !== nothing
        Set(_n6full_component_key.(seeds8)) == Set(_n6full_component_key.(seeds7)) ||
            throw(ErrorException("barrier8 seed components disagree with barrier7"))
    end
    nothing
end

"""Merge one bridge's unambiguous labels at the persistence boundary."""
function _n6full_merge_bridge!(labels::Dict{Int,String},
                               classification::KnotN6Classification,
                               bridge,
                               field::AbstractString)
    bridge === nothing && return labels
    components = Set(Int.(classification.components))
    merged = _n6full_bridge_component_labels(bridge, field)
    connections = _n6full_connection_labels(bridge)
    seeds = _n6full_getproperty(bridge, :seed_components, nothing)
    if seeds !== nothing
        for component_value in seeds
            component = _n6full_component_key(component_value)
            component in components || throw(ErrorException(
                "$field references unknown component $component"))
            connected = get(connections, component, String[])
            length(connected) <= 1 || throw(ErrorException(
                "$field component $component has multiple known labels"))
            if !isempty(connected)
                connected_label = only(connected)
                if haskey(merged, component) && merged[component] != connected_label
                    throw(ErrorException(
                        "$field resolution mismatch for component $component"))
                end
                merged[component] = connected_label
            elseif haskey(merged, component)
                throw(ErrorException(
                    "$field marks component $component resolved without a known label"))
            end
        end
    end
    for (component, label) in merged
        component in components || throw(ErrorException(
            "$field labels unknown component $component"))
        name = _n6full_label_name(label)
        name in KNOT_N6_KNOWN_NAMES || throw(ErrorException(
            "$field contains unsupported knot name $name"))
        if haskey(labels, component)
            labels[component] == label || throw(ErrorException(
                "$field relabels component $component from $(labels[component]) to $label"))
        else
            labels[component] = label
        end
    end
    labels
end

"""Return final labels without assuming that barrier 7 resolves barrier 6."""
function _n6full_component_labels(classification::KnotN6Classification,
                                  bridge7;
                                  barrier8=nothing)
    _n6full_validate_barrier8_inputs(classification, bridge7, barrier8)
    labels = Dict{Int,String}(classification.component_type_labels)
    for (component, label) in labels
        _n6full_label_name(label) in KNOT_N6_KNOWN_NAMES || throw(ErrorException(
            "barrier-6 contains unsupported knot name $(_n6full_label_name(label))"))
    end
    _n6full_merge_bridge!(labels, classification, bridge7, "barrier7")
    barrier8 === nothing || _n6full_merge_bridge!(labels, classification, barrier8,
                                                  "barrier8")
    labels
end

"""Retain only canonical serialization => final label for classical states."""
function _n6full_classical_map(classification::KnotN6Classification,
                               labels::Dict{Int,String})
    result = Dict{String,Union{Nothing,String}}()
    for (index, diagram) in enumerate(classification.states)
        serialization = canonical_serialization(diagram)
        component = classification.components[index]
        label = get(labels, component, nothing)
        if haskey(result, serialization)
            result[serialization] == label || throw(ErrorException(
                "classical serialization receives conflicting labels"))
        else
            result[serialization] = label
        end
    end
    length(result) == length(classification.states) || throw(ErrorException(
        "classical classification states are not serialization-unique"))
    result
end

function _n6full_classical_counts(classification::KnotN6Classification,
                                  classical_map)
    counts = zeros(Int, classification.max_n + 1)
    for diagram in classification.states
        n = crossing_count(diagram)
        counts[n + 1] += 1
    end
    counts == classification.state_counts || throw(ErrorException(
        "classification state_counts disagree with retained states"))
    counts
end

function _n6full_final_classification_fingerprint(classification,
                                                  labels,
                                                  bridge7,
                                                  barrier8)
    lines = String[
        "classification_version=$(KNOT_N6_CLASSIFICATION_VERSION)",
        "barrier6=$(classification.fingerprint)",
        "barrier7=$(_n6full_getproperty(bridge7, :fingerprint, "none"))",
        "barrier8=$(_n6full_getproperty(barrier8, :fingerprint, "none"))"]
    append!(lines, "component|$(component)|$(label)" for
                   (component, label) in sort!(collect(labels); by=first))
    _n6full_sha256_lines(lines)
end

function _n6full_conventions()
    Dict{String,Any}(
        "scope" => "full canonical signed one-component Gauss census; N<=6",
        "area" => "Area 2 only",
        "canonicalization" =>
            "KnotGaussDiagram v$(KNOT_GAUSS_CANON_VERSION); cyclic root shift + first-occurrence relabel only",
        "mirror_and_traversal_reversal" => "distinct",
        "strict_label_side" =>
            "+ is local KnotInfo PD representative; - is its mirror; no right/left-handed assertion",
        "classification" =>
            "classical membership from the direct census; barrier-6 graph plus barrier-7 and final barrier-8 closures",
        "nonclassical_rows" =>
            "persisted with planar=false, knot_type=NULL, classification_status=nonclassical",
        "persistence" =>
            "states and planar/knot_type/classification_status/method_version only; no transitions, paths, certificates, Carter data, or per-state invariants")
end

function _n6full_free_disk_bytes(path::AbstractString)
    Sys.iswindows() || return -1
    absolute = abspath(path)
    isempty(absolute) && return -1
    drive = uppercase(string(first(absolute)))
    try
        text = read(`powershell.exe -NoProfile -NonInteractive -Command "(Get-PSDrive -Name $drive).Free"`,
                    String)
        parse(Int, strip(text))
    catch
        -1
    end
end

function _n6full_rss_bytes()
    try
        Int(Sys.maxrss())
    catch
        -1
    end
end

function _n6full_guard!(started_ns::UInt64, db_path::AbstractString)
    elapsed = (time_ns() - started_ns) / 1e9
    elapsed < N6FULL_TIME_LIMIT_SECONDS || throw(ErrorException(
        "N6 full census exceeded $(N6FULL_TIME_LIMIT_SECONDS)s"))
    rss = _n6full_rss_bytes()
    rss < 0 || rss <= N6FULL_RSS_LIMIT_BYTES || throw(ErrorException(
        "N6 full census RSS exceeded $(N6FULL_RSS_LIMIT_BYTES) bytes"))
    free = _n6full_free_disk_bytes(db_path)
    free < 0 || free >= N6FULL_FREE_DISK_LIMIT_BYTES || throw(ErrorException(
        "N6 full census free disk fell below $(N6FULL_FREE_DISK_LIMIT_BYTES) bytes"))
    nothing
end

function _n6full_progress!(;
        n::Int,
        raw_seen::Int,
        canonical::Int,
        stored::Int,
        started_ns::UInt64,
        db_path::AbstractString,
        force::Bool=false)
    force || (stored == 0 || stored % N6FULL_PROGRESS_INTERVAL == 0 || return)
    elapsed = (time_ns() - started_ns) / 1e9
    rate = elapsed > 0 ? stored / elapsed : 0.0
    bytes = isfile(db_path) ? filesize(db_path) : 0
    free = _n6full_free_disk_bytes(db_path)
    println(stderr, "[knot-n6-full] n=$n raw_seen=$raw_seen canonical=$canonical " *
                    "stored=$stored elapsed=$(round(elapsed; digits=2))s " *
                    "rate=$(round(rate; digits=1))/s db_bytes=$bytes " *
                    "free_disk_bytes=$free")
    _n6full_guard!(started_ns, db_path)
    nothing
end

function _n6full_checkpoint_counts(completed_n::Int,
                                   n_stats,
                                   raw_seen::Int,
                                   stored::Int)
    Dict{String,Any}(
        "checkpoint_version" => 1,
        "completed_n" => completed_n,
        "stream_n_stats" => n_stats,
        "raw_seen" => raw_seen,
        "canonical_total" => stored,
        "stored_total" => stored)
end

function _n6full_pending_run_sql(db, method::String, max_n::Int;
                                 completed_n::Int=-1,
                                 n_stats=Any[],
                                 raw_seen::Int=0,
                                 stored::Int=0)
    MoveGraphStore._dbexec(db, """
        INSERT INTO classification_runs(
          method_version,canonical_version,max_n,barrier,census_fingerprint,
          classification_fingerprint,counts_json,conventions_json,created_at_utc)
        VALUES (?,?,?,?,?,?,?,?,?)
    """, (method, KNOT_GAUSS_CANON_VERSION, max_n, max(N6FULL_BARRIER8, max_n),
           "pending-census-$method", "pending-classification-$method",
           JSON.json(_n6full_checkpoint_counts(completed_n, n_stats,
                                                raw_seen, stored)),
           JSON.json(_n6full_conventions()), N6FULL_METHOD_RELEASED_AT_UTC))
end

function _n6full_ensure_pending_run!(db, method::String, max_n::Int)
    row = MoveGraphStore._one(db,
                              "SELECT * FROM classification_runs WHERE method_version=?",
                              (method,))
    if row === nothing
        _n6full_pending_run_sql(db, method, max_n)
        return MoveGraphStore._one(db,
                                   "SELECT * FROM classification_runs WHERE method_version=?",
                                   (method,))
    end
    Int(row.canonical_version) == KNOT_GAUSS_CANON_VERSION || throw(ErrorException(
        "resume classification run has the wrong canonical version"))
    Int(row.max_n) == max_n || throw(ErrorException(
        "resume classification run has the wrong max_n"))
    Int(row.barrier) == max(N6FULL_BARRIER8, max_n) || throw(ErrorException(
        "resume classification run has the wrong barrier"))
    row
end

function _n6full_json_int(value, field::AbstractString)
    value === nothing && throw(ErrorException("missing checkpoint field $field"))
    value isa Integer && return Int(value)
    value isa Real && isfinite(value) && return Int(value)
    parse(Int, String(value))
end

function _n6full_checkpoint_stats(run)
    result = Dict{Int,Any}()
    run === nothing && return result
    # `_n6full_ensure_pending_run!` returns a materialized SQLite row, whose
    # JSON field is `counts_json`; callers may also pass the public
    # `ClassificationRun`, which exposes the decoded `counts` property.
    counts = _n6full_getproperty(run, :counts, nothing)
    if counts === nothing
        encoded = _n6full_getproperty(run, :counts_json, nothing)
        encoded === nothing && throw(ErrorException(
            "classification checkpoint has no counts_json"))
        counts = JSON.parse(String(encoded))
    end
    raw = _n6full_getproperty(counts, :stream_n_stats, nothing)
    raw === nothing && throw(ErrorException(
        "classification checkpoint has no stream_n_stats"))
    raw isa AbstractVector || throw(ErrorException(
        "classification checkpoint stream_n_stats is not an array"))
    for item in raw
        n = _n6full_json_int(_n6full_getproperty(item, :n, nothing),
                            "stream_n_stats.n")
        n >= 0 || throw(ErrorException(
            "classification checkpoint has a negative stream n"))
        haskey(result, n) && throw(ErrorException(
            "duplicate checkpoint stream statistics for n=$n"))
        result[n] = item
    end
    result
end

function _n6full_checkpoint_counts_value(run)
    run === nothing && throw(ErrorException("classification checkpoint is missing"))
    counts = _n6full_getproperty(run, :counts, nothing)
    if counts === nothing
        encoded = _n6full_getproperty(run, :counts_json, nothing)
        encoded === nothing && throw(ErrorException(
            "classification checkpoint has no counts_json"))
        counts = JSON.parse(String(encoded))
    end
    counts isa AbstractDict || throw(ErrorException(
        "classification checkpoint counts are not an object"))
    counts
end

"""Validate aggregate checkpoint metadata and return its completed prefix."""
function _n6full_validate_checkpoint_metadata!(run, stats, max_n::Int)
    counts = _n6full_checkpoint_counts_value(run)
    completed_n = _n6full_json_int(
        _n6full_getproperty(counts, :completed_n, nothing), "completed_n")
    -1 <= completed_n <= max_n || throw(ErrorException(
        "classification checkpoint completed_n is outside -1:max_n"))
    keys_sorted = sort!(collect(keys(stats)))
    expected_keys = completed_n < 0 ? Int[] : collect(0:completed_n)
    keys_sorted == expected_keys || throw(ErrorException(
        "classification checkpoint stream_n_stats does not match completed_n"))
    raw_seen = _n6full_json_int(
        _n6full_getproperty(counts, :raw_seen, nothing), "raw_seen")
    canonical_total = _n6full_json_int(
        _n6full_getproperty(counts, :canonical_total, nothing),
        "canonical_total")
    stored_total = _n6full_json_int(
        _n6full_getproperty(counts, :stored_total, nothing), "stored_total")
    expected_raw = completed_n < 0 ? 0 :
        sum(N6FULL_RAW_STATE_COUNTS[1:completed_n + 1])
    expected_canonical = completed_n < 0 ? 0 :
        sum(N6FULL_CANONICAL_STATE_COUNTS[1:completed_n + 1])
    raw_seen == expected_raw || throw(ErrorException(
        "classification checkpoint raw_seen does not match completed_n"))
    canonical_total == expected_canonical || throw(ErrorException(
        "classification checkpoint canonical_total does not match completed_n"))
    stored_total == canonical_total || throw(ErrorException(
        "classification checkpoint stored_total does not match canonical_total"))
    version = _n6full_json_int(
        _n6full_getproperty(counts, :checkpoint_version, nothing),
        "checkpoint_version")
    version == 1 || throw(ErrorException(
        "unsupported classification checkpoint version $version"))
    completed_n
end

"""Expected planar/status partitions of the all-state stream by crossing n."""
function _n6full_expected_status_counts(
        classification::KnotN6Classification,
        classical_map::AbstractDict)
    max_n = classification.max_n
    planar = zeros(Int, max_n + 1)
    classified = zeros(Int, max_n + 1)
    unresolved = zeros(Int, max_n + 1)
    for diagram in classification.states
        n = crossing_count(diagram)
        serialization = canonical_serialization(diagram)
        haskey(classical_map, serialization) || throw(ErrorException(
            "classical map is missing a classification state"))
        planar[n + 1] += 1
        classical_map[serialization] === nothing ?
            unresolved[n + 1] += 1 : classified[n + 1] += 1
    end
    sum(planar) == length(classical_map) || throw(ErrorException(
        "classical map contains states outside the classification set"))
    planar, classified, unresolved
end

function _n6full_db_n_summary(db, n::Int, method::String)
    rows = MoveGraphStore._rows(db, """
        SELECT COUNT(*) AS total,
               COUNT(c.state_id) AS classified_rows,
               COALESCE(SUM(CASE WHEN c.planar=1 THEN 1 ELSE 0 END),0) AS planar_true,
               COALESCE(SUM(CASE WHEN c.planar=0 THEN 1 ELSE 0 END),0) AS planar_false,
               COALESCE(SUM(CASE WHEN c.classification_status='classified' THEN 1 ELSE 0 END),0) AS status_classified,
               COALESCE(SUM(CASE WHEN c.classification_status='unresolved' THEN 1 ELSE 0 END),0) AS status_unresolved,
               COALESCE(SUM(CASE WHEN c.classification_status='nonclassical' THEN 1 ELSE 0 END),0) AS status_nonclassical,
               COALESCE(SUM(CASE WHEN s.area=2.0 THEN 1 ELSE 0 END),0) AS area_two,
               COALESCE(SUM(CASE WHEN s.representation='knotgauss' THEN 1 ELSE 0 END),0) AS knot_gauss,
               COALESCE(SUM(CASE WHEN s.payload IS NULL THEN 1 ELSE 0 END),0) AS null_payload,
               COALESCE(SUM(CASE WHEN c.state_id IS NULL OR NOT (
                   (c.planar=0 AND c.knot_type IS NULL AND
                    c.classification_status='nonclassical')
                   OR
                   (c.planar=1 AND c.knot_type IS NULL AND
                    c.classification_status='unresolved')
                   OR
                   (c.planar=1 AND c.knot_type IS NOT NULL AND
                    c.classification_status='classified')
               ) THEN 1 ELSE 0 END),0) AS invalid_classification,
               COALESCE(SUM(CASE WHEN s.canonical_version=? THEN 0 ELSE 1 END),0) AS invalid_canonical_version
        FROM states AS s
        LEFT JOIN state_classifications AS c
          ON c.state_id=s.state_id AND c.method_version=?
        WHERE s.complexity=?
    """, (KNOT_GAUSS_CANON_VERSION, method, Float64(n)))
    isempty(rows) && throw(ErrorException("resume summary query returned no row"))
    row = first(rows)
    Dict{String,Int}(String(name) => Int(getproperty(row, name)) for name in
                     (:total, :classified_rows, :planar_true, :planar_false,
                      :status_classified, :status_unresolved,
                      :status_nonclassical, :area_two, :knot_gauss,
                      :null_payload, :invalid_classification,
                      :invalid_canonical_version))
end

function _n6full_db_serialization_fingerprint(db; n::Union{Nothing,Int}=nothing)
    sql = n === nothing ?
        "SELECT serialization FROM states ORDER BY id" :
        "SELECT serialization FROM states WHERE complexity=? ORDER BY id"
    params = n === nothing ? () : (Float64(n),)
    context = SHA.SHA256_CTX()
    count = 0
    DBInterface.execute(db, sql, params) do query
        for row in query
            text = String(row[1])
            SHA.update!(context, codeunits(text))
            SHA.update!(context, UInt8[0x0a])
            count += 1
        end
    end
    bytes2hex(SHA.digest!(context)), count
end

function _n6full_validate_checkpoint_stat!(stat, n::Int, summary,
                                           fingerprint::AbstractString)
    stat === nothing && throw(ErrorException(
        "completed n=$n has no persisted stream fingerprint checkpoint"))
    _n6full_json_int(_n6full_getproperty(stat, :raw_total, nothing),
                     "stream_n_stats.raw_total") == N6FULL_RAW_STATE_COUNTS[n + 1] ||
        throw(ErrorException("checkpoint raw_total mismatch at n=$n"))
    _n6full_json_int(_n6full_getproperty(stat, :raw_seen, nothing),
                     "stream_n_stats.raw_seen") == N6FULL_RAW_STATE_COUNTS[n + 1] ||
        throw(ErrorException("checkpoint raw_seen mismatch at n=$n"))
    _n6full_json_int(_n6full_getproperty(stat, :canonical_emitted, nothing),
                     "stream_n_stats.canonical_emitted") == summary["total"] ||
        throw(ErrorException("checkpoint canonical count mismatch at n=$n"))
    _n6full_json_int(_n6full_getproperty(stat, :planar_emitted, nothing),
                     "stream_n_stats.planar_emitted") == summary["planar_true"] ||
        throw(ErrorException("checkpoint planar count mismatch at n=$n"))
    _n6full_json_int(_n6full_getproperty(stat, :nonclassical_emitted, nothing),
                     "stream_n_stats.nonclassical_emitted") == summary["planar_false"] ||
        throw(ErrorException("checkpoint nonclassical count mismatch at n=$n"))
    String(_n6full_getproperty(stat, :emitted_fingerprint, "")) == fingerprint ||
        throw(ErrorException("checkpoint fingerprint mismatch at n=$n"))
    summary["classified_rows"] == summary["total"] ||
        throw(ErrorException("classification row count mismatch at n=$n"))
    summary["area_two"] == summary["total"] ||
        throw(ErrorException("non-Area-2 state at n=$n"))
    summary["knot_gauss"] == summary["total"] ||
        throw(ErrorException("non-KnotGauss state at n=$n"))
    summary["null_payload"] == summary["total"] ||
        throw(ErrorException("state payload is not NULL at n=$n"))
    nothing
end

"""Validate committed n-prefix rows and remove only an incomplete next n.

The stream checkpoint is deliberately kept in the dataset-level run row.  A
completed crossing number is accepted only when its rows, partition, state
metadata, and serialization fingerprint agree with both the stream constants
and the in-memory classical census.  The first incomplete crossing number is
the only prefix boundary that may be restarted.
"""
function _n6full_resume_prefix!(store::MoveGraphStore.MoveGraphStoreDB,
                                classification::KnotN6Classification,
                                method::String,
                                max_n::Int;
                                expected_planar=nothing,
                                expected_classified=nothing,
                                expected_unresolved=nothing)
    db = store.db
    run = _n6full_ensure_pending_run!(db, method, max_n)
    _n6full_query_count(db, "SELECT COUNT(*) AS n FROM classification_runs") == 1 ||
        throw(ErrorException("resume DB must contain exactly one classification run"))
    other_methods = _n6full_query_count(db,
        "SELECT COUNT(*) AS n FROM state_classifications WHERE method_version<>?",
        (method,))
    other_methods == 0 || throw(ErrorException(
        "resume DB contains classifications from another method"))
    for table in ("transitions", "certificate_paths", "invariants",
                  "knotinfo_snapshots", "knotinfo_raw_records", "knot_types",
                  "move_types")
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM $table") == 0 ||
            throw(ErrorException("resume DB contains forbidden table rows: $table"))
    end
    invalid_complexity = _n6full_query_count(db, """
        SELECT COUNT(*) AS n FROM states
        WHERE complexity IS NULL OR complexity < 0 OR complexity > ?
              OR complexity != CAST(complexity AS INTEGER)
    """, (max_n,))
    invalid_complexity == 0 || throw(ErrorException(
        "resume DB contains a state with invalid crossing complexity"))
    stats = _n6full_checkpoint_stats(run)
    checkpoint_completed = _n6full_validate_checkpoint_metadata!(run, stats, max_n)
    expected_planar === nothing && (expected_planar = copy(classification.state_counts))
    for expected in (expected_planar, expected_classified, expected_unresolved)
        expected === nothing && continue
        length(expected) == max_n + 1 || throw(ErrorException(
            "resume expected status counts have the wrong length"))
    end
    completed_stats = Any[]
    first_incomplete = max_n + 1
    for n in 0:max_n
        summary = _n6full_db_n_summary(db, n, method)
        expected_total = N6FULL_CANONICAL_STATE_COUNTS[n + 1]
        if summary["total"] == 0 && summary["classified_rows"] == 0
            first_incomplete = n
            break
        elseif summary["total"] != expected_total
            first_incomplete = n
            break
        elseif summary["planar_true"] < 0 ||
               summary["planar_true"] + summary["planar_false"] != expected_total ||
               summary["status_classified"] + summary["status_unresolved"] +
               summary["status_nonclassical"] != expected_total
            throw(ErrorException("resume DB has inconsistent status counts at n=$n"))
        end
        expected_planar === nothing ||
            summary["planar_true"] == expected_planar[n + 1] ||
            throw(ErrorException("resume planar partition mismatch at n=$n"))
        expected_classified === nothing ||
            summary["status_classified"] == expected_classified[n + 1] ||
            throw(ErrorException("resume classified partition mismatch at n=$n"))
        expected_unresolved === nothing ||
            summary["status_unresolved"] == expected_unresolved[n + 1] ||
            throw(ErrorException("resume unresolved partition mismatch at n=$n"))
        summary["planar_false"] == expected_total - summary["planar_true"] ||
            throw(ErrorException("resume nonplanar partition mismatch at n=$n"))
        summary["invalid_classification"] == 0 || throw(ErrorException(
            "resume DB has an invalid classification row at n=$n"))
        summary["invalid_canonical_version"] == 0 || throw(ErrorException(
            "resume DB has an invalid canonical version at n=$n"))
        summary["area_two"] == expected_total || throw(ErrorException(
            "resume DB has a non-Area-2 state at n=$n"))
        summary["knot_gauss"] == expected_total || throw(ErrorException(
            "resume DB has a non-KnotGauss state at n=$n"))
        summary["null_payload"] == expected_total || throw(ErrorException(
            "resume DB has a non-NULL payload at n=$n"))
        fingerprint, fingerprint_count = _n6full_db_serialization_fingerprint(db; n=n)
        fingerprint_count == expected_total || throw(ErrorException(
            "resume DB fingerprint row count mismatch at n=$n"))
        _n6full_validate_checkpoint_stat!(get(stats, n, nothing), n, summary,
                                          fingerprint)
        push!(completed_stats, stats[n])
    end
    checkpoint_completed == first_incomplete - 1 || throw(ErrorException(
        "resume checkpoint completed_n does not match persisted prefix"))
    expected_raw = first_incomplete == 0 ? 0 :
        sum(N6FULL_RAW_STATE_COUNTS[1:first_incomplete])
    expected_canonical = first_incomplete == 0 ? 0 :
        sum(N6FULL_CANONICAL_STATE_COUNTS[1:first_incomplete])
    counts = _n6full_checkpoint_counts_value(run)
    _n6full_json_int(_n6full_getproperty(counts, :raw_seen, nothing), "raw_seen") ==
        expected_raw || throw(ErrorException(
            "resume checkpoint raw_seen does not match persisted prefix"))
    _n6full_json_int(_n6full_getproperty(counts, :canonical_total, nothing),
                     "canonical_total") ==
        expected_canonical || throw(ErrorException(
            "resume checkpoint canonical_total does not match persisted prefix"))
    _n6full_json_int(_n6full_getproperty(counts, :stored_total, nothing),
                     "stored_total") ==
        expected_canonical || throw(ErrorException(
            "resume checkpoint stored_total does not match persisted prefix"))
    if first_incomplete <= max_n
        if first_incomplete < max_n
            for n in first_incomplete + 1:max_n
                summary = _n6full_db_n_summary(db, n, method)
                summary["total"] == 0 && summary["classified_rows"] == 0 ||
                    throw(ErrorException(
                        "resume DB has rows after incomplete n=$first_incomplete"))
            end
        end
        _n6full_transaction = getfield(MoveGraphStore, Symbol("_transaction"))
        _n6full_transaction(db) do
            MoveGraphStore._dbexec(db,
                "DELETE FROM states WHERE complexity>=? AND complexity<=?",
                (Float64(first_incomplete), Float64(max_n)))
            MoveGraphStore._dbexec(db,
                "UPDATE classification_runs SET counts_json=? WHERE method_version=?",
                (JSON.json(_n6full_checkpoint_counts(first_incomplete - 1,
                                                      completed_stats,
                                                      expected_raw,
                                                      expected_canonical)),
                 method))
        end
    end
    first_incomplete, completed_stats
end

function _n6full_counts(classification, classical_map, labels,
                        stream_stats, raw_total::Int;
                        snapshot_sha256::AbstractString="unknown",
                        bridge7=nothing,
                        barrier8=nothing)
    max_n = classification.max_n
    all_counts = N6FULL_CANONICAL_STATE_COUNTS[1:max_n + 1]
    all_total = sum(all_counts)
    classical_total = length(classical_map)
    nonclassical_total = all_total - classical_total
    nonclassical_total >= 0 || throw(ErrorException("negative nonclassical total"))
    classified = count(value -> value !== nothing, values(classical_map))
    unresolved = classical_total - classified
    raw_expected = sum(N6FULL_RAW_STATE_COUNTS[1:max_n + 1])
    raw_total == raw_expected || throw(ErrorException(
        "raw all-GC total mismatch: got $raw_total expected $raw_expected"))
    type_counts = Dict{String,Int}()
    for value in values(classical_map)
        value === nothing || (type_counts[value] = get(type_counts, value, 0) + 1)
    end
    Dict{String,Any}(
        "canonical_state_counts" => copy(all_counts),
        "canonical_state_total" => all_total,
        "planar_true" => classical_total,
        "planar_false" => nonclassical_total,
        "classified" => classified,
        "unresolved" => unresolved,
        "nonclassical" => nonclassical_total,
        "transitions" => 0,
        "classical_canonical_state_counts" => copy(classification.state_counts),
        "nonclassical_canonical_state_counts" =>
            [all_counts[i] - classification.state_counts[i]
             for i in eachindex(all_counts)],
        "canonical_knot_type_counts" => type_counts,
        "raw_rooted_total" => raw_total,
        "raw_rooted_classical_total" => classification.raw_rooted_classical_total,
        "barrier6_component_count" => length(classification.component_roots),
        "barrier6_unresolved_component_count" =>
            length(classification.unresolved_components),
        "barrier7_resolved_component_count" =>
            length(_n6full_getproperty(bridge7, :resolved_components, Dict())),
        "barrier7_unresolved_component_count" =>
            length(_n6full_getproperty(bridge7, :unresolved_components, Int[])),
        "barrier8_resolved_component_count" =>
            length(_n6full_getproperty(barrier8, :resolved_components, Dict())),
        "barrier8_unresolved_component_count" =>
            length(_n6full_getproperty(barrier8, :unresolved_components, Int[])),
        "barrier8_final_label_count" => length(labels),
        "knotinfo_snapshot_sha256" => String(snapshot_sha256),
        "barrier6_graph_fingerprint_sha256" => classification.fingerprint,
        "barrier7_bridge_fingerprint_sha256" =>
            _n6full_getproperty(bridge7, :fingerprint, "none"),
        "barrier8_final_fingerprint_sha256" =>
            _n6full_getproperty(barrier8, :fingerprint, "none"),
        "barrier8_status" =>
            _n6full_getproperty(barrier8, :status, "not_supplied"),
        "stream_n_stats" => stream_stats)
end

"""Insert every streamed state without retaining the all-state census.

Each crossing number is an independent commit unit.  The `stop_after_n` hook
is intentionally only a fixture fault-injection hook: it throws after the
requested crossing-number transaction has committed, leaving a valid work DB
for the next invocation to resume.
"""
function _n6full_persist_stream(path::AbstractString,
                                classification::KnotN6Classification,
                                bridge7;
                                barrier8=nothing,
                                classical_map=nothing,
                                labels=nothing,
                                max_n::Integer=classification.max_n,
                                progress::Bool=true,
                                method_version::AbstractString=N6FULL_METHOD_VERSION,
                                snapshot_sha256::AbstractString="unknown",
                                stop_after_n=nothing)
    max_n = Int(max_n)
    max_n == classification.max_n || throw(ArgumentError(
        "stream max_n must equal classification.max_n"))
    0 <= max_n <= N6FULL_MAX_N || throw(ArgumentError("max_n must lie in 0:6"))
    stop_after_n = stop_after_n === nothing ? nothing : Int(stop_after_n)
    stop_after_n === nothing || 0 <= stop_after_n <= max_n || throw(ArgumentError(
        "stop_after_n must lie in 0:max_n"))
    labels = labels === nothing ?
        _n6full_component_labels(classification, bridge7; barrier8=barrier8) : labels
    classical_map = classical_map === nothing ?
        _n6full_classical_map(classification, labels) : classical_map
    length(classical_map) == length(classification.states) || throw(ArgumentError(
        "classical_map does not cover exactly the classification states"))
    expected_classical_counts = _n6full_classical_counts(classification, classical_map)
    expected_planar, expected_classified, expected_unresolved =
        _n6full_expected_status_counts(classification, classical_map)
    expected_planar == expected_classical_counts || throw(ErrorException(
        "expected planar partition disagrees with classical state counts"))
    started_ns = time_ns()
    store = MoveGraphStore.open_store(path)
    result = nothing
    try
        method = String(method_version)
        _n6full_ensure_pending_run!(store.db, method, max_n)
        first_incomplete, completed_stats = _n6full_resume_prefix!(store,
            classification, method, max_n;
            expected_planar=expected_planar,
            expected_classified=expected_classified,
            expected_unresolved=expected_unresolved)
        n_stats = Any[completed_stats...]
        raw_seen_total = first_incomplete == 0 ? 0 :
            sum(_n6full_json_int(_n6full_getproperty(stat, :raw_seen, nothing),
                                 "stream_n_stats.raw_seen")
                for stat in completed_stats)
        canonical_total = first_incomplete == 0 ? 0 :
            sum(N6FULL_CANONICAL_STATE_COUNTS[1:first_incomplete])
        stored_total = canonical_total
        state_stmt_sql = """
            INSERT INTO states(
              state_id,state_hash,canonical_version,serialization,payload,
              area,representation,complexity)
            VALUES (?,?,?,?,?,?,?,?)
        """
        class_stmt_sql = """
            INSERT INTO state_classifications(
              state_id,planar,knot_type,classification_status,method_version)
            VALUES (?,?,?,?,?)
        """
        if first_incomplete <= max_n
            for n in first_incomplete:max_n
                stream = nothing
                stat = nothing
                MoveGraphStore._transaction(store.db) do
                    state_stmt = SQLite.Stmt(store.db, state_stmt_sql)
                    class_stmt = SQLite.Stmt(store.db, class_stmt_sql)
                    try
                        stream = stream_canonical_knot_gauss(n;
                            planar_membership=serialization -> haskey(classical_map,
                                                                      serialization),
                            callback=record -> begin
                                record.planar isa Bool || throw(ErrorException(
                                    "all-state stream returned unknown planarity"))
                                in_classical = haskey(classical_map,
                                                       record.serialization)
                                in_classical == record.planar || throw(ErrorException(
                                    "stream/classical census membership mismatch"))
                                label = in_classical ?
                                    classical_map[record.serialization] : nothing
                                status = in_classical ?
                                    (label === nothing ? "unresolved" : "classified") :
                                    "nonclassical"
                                state = CanonicalState(record.serialization;
                                    area=AREA_2, representation="KnotGauss",
                                    complexity=Float64(crossing_count(record.diagram)))
                                DBInterface.execute(state_stmt, (
                                    state.state_id, state.state_hash,
                                    state.canonical_version, state.serialization,
                                    nothing, state.area, state.representation,
                                    state.complexity))
                                DBInterface.execute(class_stmt, (
                                    state.state_id, record.planar ? 1 : 0, label,
                                    status, method))
                                nothing
                            end)
                        stream.stopped && throw(ErrorException(
                            "canonical stream stopped at n=$n: $(stream.stop_reason)"))
                        stream.raw_seen == stream.raw_total || throw(ErrorException(
                            "canonical stream did not exhaust n=$n"))
                        stream.canonical_emitted ==
                            N6FULL_CANONICAL_STATE_COUNTS[n + 1] || throw(ErrorException(
                                "canonical count mismatch at n=$n"))
                        stream.planar_emitted == expected_planar[n + 1] ||
                            throw(ErrorException("classical count mismatch at n=$n"))
                        stream.nonplanar_emitted ==
                            N6FULL_CANONICAL_STATE_COUNTS[n + 1] - expected_planar[n + 1] ||
                            throw(ErrorException("nonclassical count mismatch at n=$n"))
                        stream.planar_unknown == 0 || throw(ErrorException(
                            "stream emitted unknown planarity at n=$n"))
                        persisted_fingerprint, persisted_count =
                            _n6full_db_serialization_fingerprint(store.db; n=n)
                        persisted_count == stream.canonical_emitted ||
                            throw(ErrorException("persisted row count mismatch at n=$n"))
                        persisted_fingerprint == stream.emitted_fingerprint ||
                            throw(ErrorException(
                                "persisted stream fingerprint mismatch at n=$n"))
                        stat = Dict{String,Any}(
                            "n" => n,
                            "raw_total" => stream.raw_total,
                            "raw_seen" => stream.raw_seen,
                            "canonical_emitted" => stream.canonical_emitted,
                            "planar_emitted" => stream.planar_emitted,
                            "nonclassical_emitted" => stream.nonplanar_emitted,
                            "emitted_fingerprint" => stream.emitted_fingerprint,
                            "elapsed_seconds" => stream.elapsed_seconds)
                        next_stats = Any[n_stats..., stat]
                        next_raw_seen = raw_seen_total + stream.raw_seen
                        next_canonical = canonical_total + stream.canonical_emitted
                        MoveGraphStore._dbexec(store.db, """
                            UPDATE classification_runs
                            SET counts_json=?
                            WHERE method_version=?
                        """, (JSON.json(_n6full_checkpoint_counts(
                            n, next_stats, next_raw_seen, next_canonical)), method))
                    finally
                        SQLite.finalize(state_stmt)
                        SQLite.finalize(class_stmt)
                    end
                end
                push!(n_stats, stat)
                raw_seen_total += stream.raw_seen
                canonical_total += stream.canonical_emitted
                stored_total += stream.canonical_emitted
                if progress
                    _n6full_progress!(n=n, raw_seen=raw_seen_total,
                        canonical=canonical_total, stored=stored_total,
                        started_ns=started_ns, db_path=path, force=true)
                end
                stop_after_n === nothing || n != stop_after_n || throw(ErrorException(
                    "simulated interruption after committed checkpoint n=$n"))
            end
        end
        all_counts = N6FULL_CANONICAL_STATE_COUNTS[1:max_n + 1]
        canonical_total == sum(all_counts) || throw(ErrorException(
            "all canonical row total mismatch"))
        stored_total == canonical_total || throw(ErrorException(
            "stored row counter mismatch"))
        census_fingerprint, census_count =
            _n6full_db_serialization_fingerprint(store.db)
        census_count == canonical_total || throw(ErrorException(
            "global persisted row count mismatch"))
        counts = _n6full_counts(classification, classical_map, labels,
                                n_stats, raw_seen_total;
                                snapshot_sha256=snapshot_sha256,
                                bridge7=bridge7, barrier8=barrier8)
        counts["stream_canonical_fingerprint_sha256"] = census_fingerprint
        counts["stored_state_total"] = stored_total
        counts["stored_classification_total"] = stored_total
        merge!(counts, _n6full_checkpoint_counts(max_n, n_stats,
                                                  raw_seen_total, stored_total))
        final_classification_fingerprint =
            _n6full_final_classification_fingerprint(classification, labels,
                                                     bridge7, barrier8)
        MoveGraphStore._transaction(store.db) do
            MoveGraphStore._dbexec(store.db, """
                UPDATE classification_runs
                SET census_fingerprint=?, classification_fingerprint=?, counts_json=?
                WHERE method_version=?
            """, (census_fingerprint, final_classification_fingerprint,
                   JSON.json(counts), method))
        end
        result = (counts=counts,
                  census_fingerprint=census_fingerprint,
                  classification_fingerprint=final_classification_fingerprint,
                  n_stats=n_stats,
                  raw_seen=raw_seen_total,
                  canonical_total=canonical_total,
                  stored_total=stored_total,
                  elapsed_seconds=(time_ns() - started_ns) / 1e9)
    finally
        MoveGraphStore.close_store!(store)
    end
    result
end

function _n6full_query_count(db, sql::AbstractString, params=())
    rows = MoveGraphStore._rows(db, sql, params)
    isempty(rows) && throw(ErrorException("count query returned no row"))
    value = first(rows).n
    value === nothing ? 0 : Int(value)
end

"""Check rows, public classification aggregation, hashes, and empty tables."""
function _n6full_verify_db(path::AbstractString, expected_counts::Dict{String,Any};
                           method_version::AbstractString=N6FULL_METHOD_VERSION)
    store = MoveGraphStore.open_store(path)
    try
        db = store.db
        total = expected_counts["canonical_state_total"]
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM states") == total ||
            throw(ErrorException("stored states count mismatch"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM state_classifications") == total ||
            throw(ErrorException("stored classifications count mismatch"))
        _n6full_query_count(db, "SELECT SUM(planar=1) AS n FROM state_classifications") ==
            expected_counts["planar_true"] || throw(ErrorException("planar=true count mismatch"))
        _n6full_query_count(db, "SELECT SUM(planar=0) AS n FROM state_classifications") ==
            expected_counts["planar_false"] || throw(ErrorException("planar=false count mismatch"))
        _n6full_query_count(db, "SELECT SUM(classification_status='classified') AS n FROM state_classifications") ==
            expected_counts["classified"] || throw(ErrorException("classified count mismatch"))
        _n6full_query_count(db, "SELECT SUM(classification_status='unresolved') AS n FROM state_classifications") ==
            expected_counts["unresolved"] || throw(ErrorException("unresolved count mismatch"))
        _n6full_query_count(db, "SELECT SUM(classification_status='nonclassical') AS n FROM state_classifications") ==
            expected_counts["nonclassical"] || throw(ErrorException("nonclassical count mismatch"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM transitions") == 0 ||
            throw(ErrorException("unexpected transitions"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM certificate_paths") == 0 ||
            throw(ErrorException("unexpected certificate paths"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM invariants") == 0 ||
            throw(ErrorException("unexpected invariants"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM knotinfo_snapshots") == 0 ||
            throw(ErrorException("unexpected KnotInfo snapshot rows"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM knot_types") == 0 ||
            throw(ErrorException("unexpected knot type rows"))
        _n6full_query_count(db, "SELECT COUNT(*) AS n FROM states WHERE payload IS NOT NULL") == 0 ||
            throw(ErrorException("unexpected state payload"))
        columns = String.(getproperty.(MoveGraphStore._rows(
            db, "PRAGMA table_info(state_classifications)"), :name))
        columns == ["state_id", "planar", "knot_type", "classification_status",
                    "method_version"] || throw(ErrorException(
            "unexpected state_classifications columns"))
        api_counts = MoveGraphStore.classification_counts(store, method_version)
        api_counts["total"] == total || throw(ErrorException("API total mismatch"))
        api_counts["planar_true"] == expected_counts["planar_true"] ||
            throw(ErrorException("API planar count mismatch"))
        api_counts["planar_false"] == expected_counts["planar_false"] ||
            throw(ErrorException("API nonclassical count mismatch"))
        api_counts["classified"] == expected_counts["classified"] ||
            throw(ErrorException("API classified count mismatch"))
        api_counts["unresolved"] == expected_counts["unresolved"] ||
            throw(ErrorException("API unresolved count mismatch"))
        api_counts["nonclassical"] == expected_counts["nonclassical"] ||
            throw(ErrorException("API nonclassical status mismatch"))
        run = MoveGraphStore.get_classification_run(store, method_version)
        run === nothing && throw(ErrorException("classification run is missing"))
        first_row = first(MoveGraphStore._rows(db,
            "SELECT state_id,serialization FROM states ORDER BY id LIMIT 1"))
        first_state = MoveGraphStore.get_state(store, String(first_row.state_id))
        first_state === nothing && throw(ErrorException("public state lookup failed"))
        first_state.serialization == String(first_row.serialization) ||
            throw(ErrorException("public state serialization mismatch"))
        MoveGraphStore.get_state_classification(store, first_state;
            method_version=method_version) === nothing &&
            throw(ErrorException("public classification lookup failed"))
        true
    finally
        MoveGraphStore.close_store!(store)
    end
end

function _n6full_barrier8_summary(barrier8)
    barrier8 === nothing && return nothing
    if isdefined(@__MODULE__, :knot_n6_targeted_barrier8_summary) &&
       applicable(knot_n6_targeted_barrier8_summary, barrier8)
        return knot_n6_targeted_barrier8_summary(barrier8)
    end
    resolved = _n6full_bridge_component_labels(barrier8, "barrier8")
    Dict{String,Any}(
        "barrier" => _n6full_getproperty(barrier8, :barrier, N6FULL_BARRIER8),
        "fingerprint" => _n6full_getproperty(barrier8, :fingerprint, "none"),
        "status" => _n6full_getproperty(barrier8, :status, "unknown"),
        "unresolved_components" =>
            copy(_n6full_getproperty(barrier8, :unresolved_components, Int[])),
        "resolved_components" => Dict(string(component) => label
                                       for (component, label) in
                                       sort!(collect(resolved); by=first)))
end

function _n6full_prepared(classification::KnotN6Classification,
                          bridge7;
                          barrier8=nothing)
    labels = _n6full_component_labels(classification, bridge7; barrier8=barrier8)
    classical_map = _n6full_classical_map(classification, labels)
    (classification=classification, bridge7=bridge7, barrier8=barrier8,
     labels=labels, classical_map=classical_map)
end

"""Prepare a bounded classification/bridge input for the focused fixture."""
function _n6full_run_in_memory(max_n::Integer;
                               snapshot_path::AbstractString=N6FULL_SNAPSHOT_PATH,
                               progress::Bool=false,
                               barrier8=nothing)
    max_n = Int(max_n)
    0 <= max_n <= N6FULL_MAX_N || throw(ArgumentError("max_n must lie in 0:6"))
    classification = build_knot_n6_classification(max_n; barrier=max_n,
                                                   snapshot_path=snapshot_path,
                                                   progress=progress)
    bridge7 = run_knot_n6_targeted_bridge(classification;
                                           barrier=N6FULL_BARRIER7,
                                           progress=progress)
    _n6full_prepared(classification, bridge7; barrier8=barrier8)
end

function _n6full_sha256_fingerprint(value, field::AbstractString)
    text = String(value)
    occursin(r"^[0-9a-f]{64}$", text) || throw(ArgumentError(
        "$field must be a lowercase SHA-256 fingerprint"))
    text
end

"""Validate the non-fixture production barrier-8 contract."""
function _n6full_validate_production_barrier8!(
        classification::KnotN6Classification,
        bridge7::KnotN6TargetedBridgeResult,
        barrier8::KnotN6TargetedBarrier8Result)
    classification.max_n == N6FULL_MAX_N || throw(ArgumentError(
        "production barrier-6 classification must have max_n=6"))
    classification.barrier == N6FULL_BARRIER6 || throw(ArgumentError(
        "production classification must be a barrier-6 result"))
    bridge7.barrier == N6FULL_BARRIER7 || throw(ArgumentError(
        "production bridge must be a barrier-7 result"))
    barrier8.barrier == N6FULL_BARRIER8 || throw(ArgumentError(
        "production final result must be a barrier-8 result"))
    _n6full_sha256_fingerprint(classification.fingerprint,
                               "barrier-6 fingerprint") ==
        N6FULL_BARRIER6_FINGERPRINT || throw(ArgumentError(
            "unexpected production barrier-6 fingerprint"))
    _n6full_sha256_fingerprint(bridge7.fingerprint,
                               "barrier-7 fingerprint") ==
        N6FULL_BARRIER7_FINGERPRINT || throw(ArgumentError(
            "unexpected production barrier-7 fingerprint"))
    _n6full_sha256_fingerprint(barrier8.fingerprint,
                               "barrier-8 fingerprint") ==
        N6FULL_BARRIER8_FINGERPRINT || throw(ArgumentError(
            "unexpected production barrier-8 fingerprint"))
    barrier8.n6_input_fingerprint == classification.fingerprint ||
        throw(ArgumentError("barrier-8 N6 input fingerprint mismatch"))
    barrier8.barrier7_fingerprint == bridge7.fingerprint ||
        throw(ArgumentError("barrier-8 barrier-7 fingerprint mismatch"))

    seed_components = sort!(copy(bridge7.unresolved_components))
    seed_components == N6FULL_BARRIER8_UNRESOLVED_COMPONENTS || throw(ArgumentError(
        "production barrier-7 unresolved components do not equal 9:18"))
    length(unique(seed_components)) == length(seed_components) ||
        throw(ArgumentError("production barrier-7 unresolved components are duplicated"))
    barrier8.seed_components == seed_components || throw(ArgumentError(
        "barrier-8 seed components do not match the ten barrier-7 unresolved components"))
    barrier8.unresolved_components == seed_components || throw(ArgumentError(
        "barrier-8 unresolved components do not match the ten barrier-7 seeds"))
    barrier8.status == N6FULL_BARRIER8_STATUS ||
        throw(ArgumentError("unexpected production barrier-8 status: $(barrier8.status)"))
    isempty(strip(String(barrier8.resource_reason))) || throw(ArgumentError(
        "barrier-8 reported a resource limit despite its resolved status"))
    nothing
end

"""Run the canonical production barrier-6 -> barrier-7 -> barrier-8 chain."""
function _n6full_production_inputs(;
        snapshot_path::AbstractString=N6FULL_SNAPSHOT_PATH,
        progress::Bool=true)
    classification = build_knot_n6_classification(
        N6FULL_MAX_N; barrier=N6FULL_BARRIER6,
        snapshot_path=snapshot_path, progress=progress)
    workspace_ref = Ref{Any}(nothing)
    bridge7 = run_knot_n6_targeted_bridge(
        classification; barrier=N6FULL_BARRIER7, progress=progress,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    workspace_ref[] === nothing && throw(ArgumentError(
        "barrier-7 did not provide the workspace required by barrier-8"))
    barrier8 = run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[];
        barrier=N6FULL_BARRIER8, progress=progress)
    _n6full_validate_production_barrier8!(classification, bridge7, barrier8)
    (classification=classification, bridge7=bridge7, barrier8=barrier8)
end

function _n6full_payload(prepared, result, snapshot_sha::AbstractString,
                         snapshot_bytes::Integer, db_sha::AbstractString,
                         db_bytes::Integer, elapsed_seconds::Real)
    classification = prepared.classification
    bridge7 = prepared.bridge7
    barrier8 = prepared.barrier8
    counts = result.counts
    runtime = Dict{String,Any}(
        "julia_version" => string(VERSION),
        "command" => "julia --startup-file=no --pkgimages=no --project=. scripts/build_knot_n6_full_classification_db.jl",
        "os" => string(Sys.KERNEL),
        "machine" => string(Sys.MACHINE),
        "elapsed_seconds" => elapsed_seconds)
    input = Dict{String,Any}(
        "snapshot" => "threads/PJ260830/data/knotinfo-n6.json",
        "snapshot_sha256" => snapshot_sha,
        "snapshot_bytes" => snapshot_bytes,
        "classical_canonical_fingerprint_sha256" =>
            _n6full_sha256_lines(sort!(collect(keys(prepared.classical_map)))),
        "barrier6_graph_fingerprint_sha256" => classification.fingerprint,
        "barrier7_bridge_fingerprint_sha256" =>
            _n6full_getproperty(bridge7, :fingerprint, "none"),
        "barrier8_final_fingerprint_sha256" =>
            _n6full_getproperty(barrier8, :fingerprint, "none"))
    method = Dict{String,Any}(
        "classification_version" => KNOT_N6_CLASSIFICATION_VERSION,
        "enumeration_version" => KNOT_N6_ALL_STATE_ENUMERATION_VERSION,
        "barrier6" => N6FULL_BARRIER6,
        "barrier7" => N6FULL_BARRIER7,
        "barrier8" => N6FULL_BARRIER8,
        "conventions" => _n6full_conventions())
    component_labels = Dict(string(component) => label
                            for (component, label) in
                            sort!(collect(prepared.labels); by=first))
    barrier8_summary = _n6full_barrier8_summary(barrier8)
    result_payload = Dict{String,Any}(
        "counts" => counts,
        "census_fingerprint_sha256" => result.census_fingerprint,
        "classification_fingerprint_sha256" => result.classification_fingerprint,
        "component_labels" => component_labels,
        "barrier7_summary" => bridge7 isa KnotN6TargetedBridgeResult ?
            knot_n6_targeted_bridge_summary(bridge7) : Dict{String,Any}(),
        "barrier8_summary" => barrier8_summary)
    database = Dict{String,Any}(
        "path" => "references/items/database/knotinfo/knot-n6-full-classification.sqlite3",
        "sha256" => db_sha,
        "bytes" => db_bytes)
    Dict{String,Any}(
        "schema_version" => 2,
        "method_version" => N6FULL_METHOD_VERSION,
        "method_released_at_utc" => N6FULL_METHOD_RELEASED_AT_UTC,
        "scope" => "full canonical census",
        "runtime" => runtime,
        "input" => input,
        "method" => method,
        "result" => result_payload,
        "database" => database)
end

function _n6full_write_json(path::AbstractString, payload)
    open(path, "w") do io
        JSON.print(io, payload)
        write(io, '\n')
    end
end

function _n6full_write_markdown(path::AbstractString, payload)
    counts = payload["result"]["counts"]
    input = payload["input"]
    db = payload["database"]
    runtime = payload["runtime"]
    open(path, "w") do io
        println(io, "# Knot N<=6 full canonical classification DB")
        println(io)
        println(io, "- scope: `full canonical signed one-component Gauss census`")
        println(io, "- method_version: `$(payload["method_version"])`")
        println(io, "- fixed method release UTC: `$(payload["method_released_at_utc"])`")
        println(io, "- KnotInfo snapshot SHA-256: `$(input["snapshot_sha256"])`")
        println(io, "- barrier-6 graph fingerprint: `$(input["barrier6_graph_fingerprint_sha256"])`")
        println(io, "- barrier-7 bridge fingerprint: `$(input["barrier7_bridge_fingerprint_sha256"])`")
        println(io, "- barrier-8 final fingerprint: `$(input["barrier8_final_fingerprint_sha256"])`")
        println(io, "- canonical census fingerprint: `$(payload["result"]["census_fingerprint_sha256"])`")
        println(io)
        println(io, "## Counts")
        println(io)
        println(io, "- canonical states: `$(counts["canonical_state_total"])`")
        println(io, "- planar=true/classical: `$(counts["planar_true"])`")
        println(io, "- planar=false/nonclassical: `$(counts["planar_false"])`")
        println(io, "- classified: `$(counts["classified"])`; unresolved: `$(counts["unresolved"])`; nonclassical: `$(counts["nonclassical"])`")
        println(io, "- transitions: `$(counts["transitions"])`")
        println(io)
        println(io, "| n | canonical states |")
        println(io, "|---:|---:|")
        for (n, count) in enumerate(counts["canonical_state_counts"])
            println(io, "| $(n - 1) | $count |")
        end
        println(io)
        println(io, "Only `states` and the five `state_classifications` columns are populated.")
        println(io, "Transitions, certificates, Carter data, per-state invariants, and paths are not persisted.")
        println(io, "The `+`/`-` labels indicate local-PD/mirror sides only; no handedness assertion is made.")
        println(io)
        println(io, "- database: `$(db["path"])`")
        println(io, "- database SHA-256: `$(db["sha256"])` ($(db["bytes"]) bytes)")
        println(io, "- elapsed seconds: `$(runtime["elapsed_seconds"])`")
        println(io, "- command: `$(runtime["command"])`")
    end
end

"""Run the production build; barrier-8 is a mandatory explicit dependency."""
function build_knot_n6_full_classification_db(; 
        snapshot_path::AbstractString=N6FULL_SNAPSHOT_PATH,
        classification=nothing,
        bridge7=nothing,
        barrier8=nothing,
        progress::Bool=true,
        resume_db_path::Union{Nothing,AbstractString}=nothing)
    barrier8 === nothing && throw(ArgumentError(
        "final barrier-8 result is required; production full run is disabled until supplied"))
    started = time_ns()
    classification = classification === nothing ?
        build_knot_n6_classification(N6FULL_MAX_N; barrier=N6FULL_BARRIER6,
                                     snapshot_path=snapshot_path, progress=progress) :
        classification
    classification.max_n == N6FULL_MAX_N || throw(ArgumentError(
        "production classification must have max_n=6"))
    bridge7 = bridge7 === nothing ?
        run_knot_n6_targeted_bridge(classification; barrier=N6FULL_BARRIER7,
                                    progress=progress) : bridge7
    prepared = _n6full_prepared(classification, bridge7; barrier8=barrier8)
    snapshot_sha, snapshot_bytes = _n6full_sha256_file(snapshot_path)
    length(prepared.classical_map) == N6FULL_CLASSICAL_CANONICAL_TOTAL ||
        throw(ErrorException("production classical canonical total mismatch"))
    mkpath(dirname(N6FULL_DB_OUTPUT))
    mkpath(dirname(N6FULL_JSON_OUTPUT))
    mkpath(dirname(N6FULL_MD_OUTPUT))
    tmp_db = resume_db_path === nothing ?
        tempname(dirname(N6FULL_DB_OUTPUT)) * ".sqlite3" :
        normpath(abspath(String(resume_db_path)))
    owns_tmp_db = resume_db_path === nothing
    tmp_json = tempname(dirname(N6FULL_JSON_OUTPUT)) * ".json"
    tmp_md = tempname(dirname(N6FULL_MD_OUTPUT)) * ".md"
    try
        result = _n6full_persist_stream(tmp_db, prepared.classification,
            prepared.bridge7; barrier8=prepared.barrier8,
            classical_map=prepared.classical_map, labels=prepared.labels,
            max_n=N6FULL_MAX_N, progress=progress,
            snapshot_sha256=snapshot_sha)
        result.counts["canonical_state_total"] == N6FULL_FULL_CANONICAL_TOTAL ||
            throw(ErrorException("production all canonical total mismatch"))
        result.counts["planar_false"] == N6FULL_NONCLASSICAL_CANONICAL_TOTAL ||
            throw(ErrorException("production nonclassical total mismatch"))
        _n6full_verify_db(tmp_db, result.counts)
        db_sha, db_bytes = _n6full_sha256_file(tmp_db)
        elapsed = (time_ns() - started) / 1e9
        payload = _n6full_payload(prepared, result, snapshot_sha, snapshot_bytes,
                                   db_sha, db_bytes, elapsed)
        _n6full_write_json(tmp_json, payload)
        _n6full_write_markdown(tmp_md, payload)
        JSON.parsefile(tmp_json)
        mv(tmp_db, N6FULL_DB_OUTPUT; force=true)
        mv(tmp_json, N6FULL_JSON_OUTPUT; force=true)
        mv(tmp_md, N6FULL_MD_OUTPUT; force=true)
        payload
    catch
        for path in (owns_tmp_db ? tmp_db : "", tmp_json, tmp_md)
            isfile(path) && rm(path; force=true)
        end
        rethrow()
    end
end

function _n6full_parse_args(args)
    progress = true
    resume_db_path = nothing
    i = 1
    while i <= length(args)
        arg = String(args[i])
        if arg == "--quiet"
            progress = false
            i += 1
        elseif arg == "--resume"
            i < length(args) || throw(ArgumentError(
                "--resume requires a work DB path"))
            resume_db_path = String(args[i + 1])
            isempty(strip(resume_db_path)) && throw(ArgumentError(
                "--resume requires a non-empty work DB path"))
            i += 2
        elseif startswith(arg, "--resume=")
            resume_db_path = split(arg, "="; limit=2)[2]
            isempty(strip(resume_db_path)) && throw(ArgumentError(
                "--resume requires a non-empty work DB path"))
            i += 1
        else
            throw(ArgumentError("unsupported argument $arg; use --quiet and/or --resume PATH"))
        end
    end
    (progress=progress, resume_db_path=resume_db_path)
end

"""Resolve the CLI route without executing the expensive census."""
function _n6full_cli_dispatch(args=ARGS)
    options = _n6full_parse_args(args)
    (mode=:production_full, progress=options.progress,
     resume_db_path=options.resume_db_path)
end

function main(args=ARGS)
    dispatch = _n6full_cli_dispatch(args)
    inputs = _n6full_production_inputs(; progress=dispatch.progress)
    payload = build_knot_n6_full_classification_db(
        ; classification=inputs.classification, bridge7=inputs.bridge7,
        barrier8=inputs.barrier8, progress=dispatch.progress,
        resume_db_path=dispatch.resume_db_path)
    counts = payload["result"]["counts"]
    println("method_version=$(payload["method_version"])")
    println("canonical_state_total=$(counts["canonical_state_total"]);planar_true=$(counts["planar_true"]);planar_false=$(counts["planar_false"])")
    println("classified=$(counts["classified"]);unresolved=$(counts["unresolved"]);nonclassical=$(counts["nonclassical"]);transitions=$(counts["transitions"])")
    println("census_fingerprint=$(payload["result"]["census_fingerprint_sha256"])")
    println("classification_fingerprint=$(payload["result"]["classification_fingerprint_sha256"])")
    println("database_bytes=$(payload["database"]["bytes"]);database_sha256=$(payload["database"]["sha256"])")
    println("elapsed_seconds=$(payload["runtime"]["elapsed_seconds"])")
    payload
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
