"""Build the minimal v2 persistence artifact for the barrier-5 N<=5 run.

The classification engine remains the source of truth.  This builder stores
one row for each canonical GC state and one five-column classification row;
it deliberately does not serialize Carter data, certificates, paths,
invariants, rooted records, or direct transitions.  The all-state census is
kept in memory only long enough to write the two required tables.
"""

using JSON
using SHA
using SQLite
using DBInterface

include(joinpath(@__DIR__, "..", "src", "move_graph_store.jl"))
include(joinpath(@__DIR__, "..", "src", "knot_n5_classification.jl"))
using .MoveGraphStore

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PJ260830_ROOT = joinpath(REPO_ROOT, "threads", "PJ260830")
const KNOTINFO_SNAPSHOT = joinpath(PJ260830_ROOT, "data", "knotinfo-n6.json")
const DB_OUTPUT = joinpath(REPO_ROOT, "references", "items", "database", "knotinfo",
                           "knot-n5-classification.sqlite3")
const JSON_OUTPUT = joinpath(PJ260830_ROOT, "data", "knot-n5-classification.json")
const MD_OUTPUT = joinpath(PJ260830_ROOT, "reports", "knot-n5-classification.md")

const METHOD_VERSION = "knot-n5-classification-v1"
# Captured from the operating system's UTC clock when this method release was
# stowed (2026-09-01 JST).  This identifies the method release, not a fresh
# generation time, so rebuilding the same method keeps the DB hash stable.
const METHOD_RELEASED_AT_UTC = "2026-08-31T15:26:06Z"
const CENSUS_FINGERPRINT =
    "85f8bf2cc5cf7b932350b62f67e0e539b0235cbbddc1a1c4ae54c59d5b73bdc9"
const CLASSIFICATION_FINGERPRINT =
    "35c4baa4d1800266b5334e64fc85476a49f9b600347e551176ff6563ee861047"
const ALL_CANONICAL_STATE_COUNTS = [1, 2, 14, 164, 3388, 96776]

function _sha256_file(path::AbstractString)
    bytes = read(path)
    bytes2hex(SHA.sha256(bytes)), length(bytes)
end

"Enumerate all rooted inputs once and retain canonical states only." 
function _all_canonical_states(max_n::Int)
    states = KnotGaussDiagram[]
    state_index = Dict{KnotGaussDiagram,Int}()
    state_counts = zeros(Int, max_n + 1)
    raw_total = 0
    for n in 0:max_n
        for diagram in enumerate_knot_gauss(n; orientations=true)
            raw_total += 1
            canonical = canonicalize(diagram)
            haskey(state_index, canonical) && continue
            push!(states, canonical)
            state_index[canonical] = length(states)
            state_counts[n + 1] += 1
        end
    end
    states, state_index, state_counts, raw_total
end

"Assign only the strict local-PD/mirror labels used by this N5 dataset." 
function _classify_all_states(states, state_index, n5)
    component_types = _n5_component_knot_types(n5.known_representatives)
    strict_types = Set(KNOT_N5_STRICT_KNOT_TYPES)
    labels = Vector{Union{Nothing,String}}(undef, length(states))
    planar = falses(length(states))
    canonical_counts = Dict{String,Int}(type => 0 for type in KNOT_N5_STRICT_KNOT_TYPES)
    for (index, diagram) in enumerate(states)
        is_planar = classical_gauss_realizable_fast(diagram)
        planar[index] = is_planar
        if is_planar
            n5_index = get(n5.state_index, diagram, 0)
            n5_index > 0 || throw(ErrorException(
                "planar canonical state is absent from N5 classification"))
            type = get(component_types, n5.components[n5_index], nothing)
            type === nothing && throw(ErrorException(
                "planar canonical state has no known N5 component label"))
            type in strict_types || throw(ErrorException("unexpected N5 strict label $type"))
            labels[index] = type
            canonical_counts[type] += 1
        else
            labels[index] = nothing
        end
    end
    planar, labels, canonical_counts
end

function _summary_counts(all_state_counts, raw_total, planar, labels, n5)
    planar_count = count(identity, planar)
    nonclassical_count = length(planar) - planar_count
    strict_counts = Dict{String,Int}(type => 0 for type in KNOT_N5_STRICT_KNOT_TYPES)
    for label in labels
        label === nothing || (strict_counts[label] += 1)
    end
    strict_counts == n5.canonical_knot_type_counts || throw(ErrorException(
        "strict canonical type counts disagree with N5 classification"))
    sum(values(strict_counts)) == planar_count == 6361 ||
        throw(ErrorException("strict canonical type totals do not equal planar count"))
    sum(values(n5.raw_rooted_classical_counts)) == n5.raw_rooted_classical_total == 62005 ||
        throw(ErrorException("strict raw rooted type totals do not equal classical total"))
    n5.raw_rooted_classical_unresolved_count == 0 ||
        throw(ErrorException("raw rooted classical input has unresolved labels"))
    all_state_counts == ALL_CANONICAL_STATE_COUNTS || throw(ErrorException(
        "all canonical state count mismatch: got $all_state_counts"))
    raw_total == 995573 || throw(ErrorException("raw rooted total mismatch: $raw_total"))
    Dict{String,Any}(
        "canonical_state_counts" => copy(all_state_counts),
        "canonical_state_total" => length(planar),
        "planar_true" => planar_count,
        "planar_false" => nonclassical_count,
        "classified" => planar_count,
        "unresolved" => 0,
        "nonclassical" => nonclassical_count,
        "component_count" => length(n5.component_roots),
        "known_representative_count" => length(n5.known_representatives),
        "transitions" => 0,
        "canonical_knot_type_counts" => strict_counts,
        "raw_rooted_total" => raw_total,
        "raw_rooted_classical_total" => n5.raw_rooted_classical_total,
        "raw_rooted_classical_counts" => copy(n5.raw_rooted_classical_counts),
        "raw_rooted_classical_unresolved" => n5.raw_rooted_classical_unresolved_count)
end

function _persist_db(path, states, planar, labels, counts, input_sha)
    store = open_store(path)
    try
        conventions = Dict{String,Any}(
            "canonicalization" =>
                "cyclic root shift + first-occurrence relabel only",
            "mirror_and_traversal_reversal" => "distinct",
            "strict_knot_type_labels" =>
                "0_1, 3_1^+, 3_1^-, 4_1, 5_1^+, 5_1^-, 5_2^+, 5_2^-",
            "strict_label_side" =>
                "+ is the local KnotInfo PD representative side; - is its mirror; " *
                "no right/left-handed identification",
            "classification" =>
                "Carter genus 0 plus barrier-5 R1/R2/R3 connectivity",
            "persistence" =>
                "canonical state and planar/knot_type/status only; no paths, " *
                "Carter data, invariants, or transitions")
        run_counts = copy(counts)
        run_counts["census_fingerprint"] = CENSUS_FINGERPRINT
        run_counts["classification_fingerprint"] = CLASSIFICATION_FINGERPRINT
        run_counts["input_knotinfo_sha256"] = input_sha
        register_classification_run!(store; method_version=METHOD_VERSION,
            canonical_version=KNOT_GAUSS_CANON_VERSION, max_n=5, barrier=5,
            census_fingerprint=CENSUS_FINGERPRINT,
            classification_fingerprint=CLASSIFICATION_FINGERPRINT,
            counts=run_counts, conventions=conventions,
            created_at_utc=METHOD_RELEASED_AT_UTC)

        # One transaction covers every state and classification row.  The
        # public one-row API performs the same Julia validation; this bulk
        # path avoids 200k SQLite BEGIN/COMMIT cycles while SQL CHECKs and the
        # Area-2 trigger remain active for every insert.
        MoveGraphStore._transaction(store.db) do
            for (index, diagram) in enumerate(states)
                state = CanonicalState(canonical_serialization(diagram);
                    area=AREA_2, representation="KnotGauss",
                    complexity=Float64(crossing_count(diagram)))
                put_state!(store, state)
                is_planar = planar[index]
                type = labels[index]
                status = is_planar ? "classified" : "nonclassical"
                is_planar && type === nothing && throw(ErrorException(
                    "classified state has no strict knot type"))
                !is_planar && type !== nothing && throw(ErrorException(
                    "nonclassical state has a knot type"))
                MoveGraphStore._dbexec(store.db, """
                    INSERT INTO state_classifications(
                      state_id,planar,knot_type,classification_status,method_version)
                    VALUES (?,?,?,?,?)
                """, (state.state_id, is_planar ? 1 : 0, type, status,
                       METHOD_VERSION))
            end
        end
        db_counts = classification_counts(store, METHOD_VERSION)
        db_counts["total"] == length(states) || throw(ErrorException(
            "stored classification row count mismatch"))
        db_counts["planar_true"] == counts["planar_true"] || throw(ErrorException(
            "stored planar count mismatch"))
        db_counts["planar_false"] == counts["planar_false"] || throw(ErrorException(
            "stored nonclassical count mismatch"))
        db_counts["classified"] == counts["classified"] || throw(ErrorException(
            "stored classified count mismatch"))
        db_counts["unresolved"] == 0 || throw(ErrorException(
            "stored unresolved count is nonzero"))
        db_counts["nonclassical"] == counts["nonclassical"] || throw(ErrorException(
            "stored nonclassical status count mismatch"))
        db_counts["knot_type_counts"] == counts["canonical_knot_type_counts"] ||
            throw(ErrorException("stored strict type counts mismatch"))
        first(MoveGraphStore._rows(store.db,
            "SELECT COUNT(*) AS n FROM transitions")).n == 0 ||
            throw(ErrorException("classification DB unexpectedly contains transitions"))
    finally
        close_store!(store)
    end
end

function _write_json(path, payload)
    open(path, "w") do io
        JSON.print(io, payload)
        write(io, '\n')
    end
end

function _write_markdown(path, payload)
    counts = payload["result"]["counts"]
    types = counts["canonical_knot_type_counts"]
    raw_types = counts["raw_rooted_classical_counts"]
    open(path, "w") do io
        println(io, "# Knot N<=5 classification DB")
        println(io)
        println(io, "- method_version: `$(payload["method_version"])`")
        println(io, "- schema_version: `$(payload["schema_version"])`")
        println(io, "- Julia: `$(payload["runtime"]["julia_version"])`")
        println(io, "- command: `$(payload["runtime"]["command"])`")
        println(io, "- input SHA-256: `$(payload["input"]["sha256"])`")
        println(io, "- census fingerprint: `$(payload["result"]["census_fingerprint"])`")
        println(io, "- classification fingerprint: `$(payload["result"]["classification_fingerprint"])`")
        println(io)
        println(io, "## Counts")
        println(io)
        println(io, "- canonical states: `$(counts["canonical_state_total"])` " *
                  "(planar `$(counts["planar_true"])`, nonclassical `$(counts["planar_false"])`)")
        println(io, "- statuses: classified `$(counts["classified"])`, " *
                  "unresolved `$(counts["unresolved"])`, nonclassical `$(counts["nonclassical"])`")
        println(io, "- components: `$(counts["component_count"])`; transitions: `$(counts["transitions"])`")
        println(io, "- raw rooted inputs: `$(counts["raw_rooted_total"])`; " *
                  "classical: `$(counts["raw_rooted_classical_total"])`")
        println(io)
        println(io, "| strict knot_type | canonical states | raw rooted classical |")
        println(io, "|---|---:|---:|")
        for type in KNOT_N5_STRICT_KNOT_TYPES
            println(io, "| `$(type)` | $(types[type]) | $(raw_types[type]) |")
        end
        println(io)
        println(io, "## Database")
        println(io)
        println(io, "- path: `$(payload["database"]["path"])`")
        println(io, "- bytes: `$(payload["database"]["bytes"])`")
        println(io, "- SHA-256: `$(payload["database"]["sha256"])`")
        println(io)
        println(io, "Reproduce with:")
        println(io)
        println(io, "```text")
        println(io, payload["runtime"]["command"])
        println(io, "```")
    end
end

function main(args=ARGS)
    isempty(args) || throw(ArgumentError("this builder accepts no arguments"))
    started = time_ns()
    input_sha, input_bytes = _sha256_file(KNOTINFO_SNAPSHOT)
    n5 = build_knot_n5_classification(5; barrier=5,
                                       snapshot_path=KNOTINFO_SNAPSHOT)
    states, state_index, all_state_counts, raw_total = _all_canonical_states(5)
    planar, labels, _ = _classify_all_states(states, state_index, n5)
    counts = _summary_counts(all_state_counts, raw_total, planar, labels, n5)

    for (type, count) in n5.canonical_knot_type_counts
        count == get(counts["canonical_knot_type_counts"], type, -1) ||
            throw(ErrorException("N5 strict count mismatch for $type"))
    end
    mkpath(dirname(DB_OUTPUT)); mkpath(dirname(JSON_OUTPUT)); mkpath(dirname(MD_OUTPUT))
    tmp_db = tempname(dirname(DB_OUTPUT)) * ".sqlite3"
    tmp_json = tempname(dirname(JSON_OUTPUT)) * ".json"
    tmp_md = tempname(dirname(MD_OUTPUT)) * ".md"
    try
        _persist_db(tmp_db, states, planar, labels, counts, input_sha)
        db_sha, db_bytes = _sha256_file(tmp_db)
        payload = Dict{String,Any}(
            "schema_version" => 2,
            "method_version" => METHOD_VERSION,
            "method_released_at_utc" => METHOD_RELEASED_AT_UTC,
            "runtime" => Dict{String,Any}(
                "julia_version" => string(VERSION),
                "command" =>
                    "julia --startup-file=no --pkgimages=no --project=. " *
                    "scripts/build_knot_n5_classification_db.jl"),
            "input" => Dict{String,Any}(
                "snapshot" => "threads/PJ260830/data/knotinfo-n6.json",
                "bytes" => input_bytes, "sha256" => input_sha,
                "census_fingerprint" => CENSUS_FINGERPRINT),
            "result" => Dict{String,Any}(
                "census_fingerprint" => CENSUS_FINGERPRINT,
                "classification_fingerprint" => CLASSIFICATION_FINGERPRINT,
                "counts" => counts),
            "database" => Dict{String,Any}(
                "path" => "references/items/database/knotinfo/knot-n5-classification.sqlite3",
                "bytes" => db_bytes, "sha256" => db_sha))
        _write_json(tmp_json, payload)
        _write_markdown(tmp_md, payload)
        # Each artifact is atomically replaced only after all generation and
        # validation above have succeeded; failures remove temporary files.
        mv(tmp_db, DB_OUTPUT; force=true)
        mv(tmp_json, JSON_OUTPUT; force=true)
        mv(tmp_md, MD_OUTPUT; force=true)
        tmp_db = tmp_json = tmp_md = ""
        elapsed = (time_ns() - started) / 1e9
        println("method_version=$(METHOD_VERSION)")
        println("canonical_state_total=$(counts["canonical_state_total"]);planar_true=$(counts["planar_true"]);" *
                "planar_false=$(counts["planar_false"])")
        println("classified=$(counts["classified"]);unresolved=$(counts["unresolved"]);" *
                "nonclassical=$(counts["nonclassical"]);component_count=$(counts["component_count"]);" *
                "transitions=$(counts["transitions"])")
        println("raw_rooted_total=$(counts["raw_rooted_total"]);" *
                "raw_rooted_classical_total=$(counts["raw_rooted_classical_total"])")
        for type in KNOT_N5_STRICT_KNOT_TYPES
            println("knot_type=$(type);canonical=$(counts["canonical_knot_type_counts"][type]);" *
                    "raw_rooted_classical=$(counts["raw_rooted_classical_counts"][type])")
        end
        println("census_fingerprint=$(CENSUS_FINGERPRINT)")
        println("classification_fingerprint=$(CLASSIFICATION_FINGERPRINT)")
        println("database_bytes=$(db_bytes);database_sha256=$(db_sha)")
        println("elapsed_seconds=$(elapsed)")
        println("unresolved_status=barrier5_resolved_for_enumerated_classical_states")
        println("command=$(payload["runtime"]["command"])")
        payload
    catch
        for path in (tmp_db, tmp_json, tmp_md)
            isempty(path) || (isfile(path) && rm(path))
        end
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
