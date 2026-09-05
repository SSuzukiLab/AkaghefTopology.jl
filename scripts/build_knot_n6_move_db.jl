# Build the bounded Area-2 one-step move cache for the 281 campaign seeds.
using Topology
using JSON
using SHA
using Dates

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PROJECT_ROOT = joinpath(ROOT, "threads", "PJ260830")
const CAMPAIGN_DEFAULT = joinpath(PROJECT_ROOT, "data", "knot-n6-campaign.json")
const KNOTINFO_DEFAULT = joinpath(PROJECT_ROOT, "data", "knotinfo-n6.json")
const DB_DEFAULT = joinpath(ROOT, "references", "items", "database", "knotinfo", "knot-n6-moves.sqlite3")
const JSON_DEFAULT = joinpath(PROJECT_ROOT, "data", "knot-n6-move-db.json")
const MD_DEFAULT = joinpath(PROJECT_ROOT, "reports", "knot-n6-move-db.md")
const COMMAND = "julia --startup-file=no --pkgimages=no --project=. scripts/build_knot_n6_move_db.jl"

_hash(path) = bytes2hex(SHA.sha256(read(path)))

function _parse_gauss(x)
    text = strip(String(x)); isempty(text) && return Int[]
    body = strip(text[2:end-1]); isempty(body) && return Int[]
    parse.(Int, strip.(split(body, ',')))
end

function _seed_diagrams(campaign)
    result = Tuple{String,KnotGaussDiagram}[]
    for seed in campaign["official_seed_decorations"]
        code = _parse_gauss(seed["gauss_notation"])
        for decoration in seed["decorations"]
            push!(result, (String(seed["name"]),
                          KnotGaussDiagram(code, Int[Int(x) for x in decoration["orientation"]])))
        end
    end
    result
end

function build(; campaign_path=CAMPAIGN_DEFAULT, knotinfo_path=KNOTINFO_DEFAULT,
               db_path=DB_DEFAULT, json_path=JSON_DEFAULT, markdown_path=MD_DEFAULT)
    started = time_ns()
    campaign = JSON.parsefile(campaign_path)
    seeds = _seed_diagrams(campaign)
    length(seeds) == 281 || error("expected exactly 281 campaign decorations")
    mkpath(dirname(db_path)); mkpath(dirname(json_path)); mkpath(dirname(markdown_path))
    temporary_db = tempname(dirname(db_path)) * ".sqlite3"
    temporary_json = tempname(dirname(json_path)) * ".json"
    temporary_markdown = tempname(dirname(markdown_path)) * ".md"
    store = nothing
    try
        store = Topology.MoveGraphStore.open_store(temporary_db)
        imported_knotinfo_rows = Topology.MoveGraphStore.import_knotinfo!(store, knotinfo_path)
        seed_ids = Set{String}()
        seed_dedupes = 0
        for (_, diagram) in seeds
            state = knot_canonical_state(diagram)
            state.state_id in seed_ids && (seed_dedupes += 1)
            push!(seed_ids, state.state_id)
            Topology.MoveGraphStore.put_state!(store, state)
        end
        length(seed_ids) == 215 || error("expected exactly 215 unique seed states")
        seed_dedupes == 66 || error("expected exactly 66 seed deduplications")

        candidates = persisted = rejected = over_limit = 0
        errors = Dict{String,Int}()
        for (seed_index, (_, diagram)) in enumerate(seeds)
            seed_index % 25 == 0 && println("[knot-n6] seed $(seed_index)/$(length(seeds)); candidates=$(candidates); db_open=$(isopen(store.db))")
            isopen(store.db) || error("SQLite store closed before seed $(seed_index)")
            for typ in ("R1", "R2", "R3")
                for candidate in movable(diagram, typ)
                    candidates += 1
                    after = apply_reidemeister(diagram, candidate)
                    if crossing_count(after) > 6
                        over_limit += 1
                        continue
                    end
                    persist_reidemeister_transition!(store, diagram, candidate)
                    persisted += 1
                end
            end
        end
        rejected == 0 || error("unexpected candidate rejections: $rejected")
        isempty(errors) || error("unexpected candidate errors: $errors")
        candidates > 12 || error("candidate census unexpectedly small: $candidates")
        state_count = Int(first(Topology.MoveGraphStore._rows(store.db, "SELECT count(*) AS n FROM states")).n)
        edge_count = Int(first(Topology.MoveGraphStore._rows(store.db, "SELECT count(*) AS n FROM transitions")).n)
        sqlite_ver = String(first(Topology.MoveGraphStore._rows(store.db, "SELECT sqlite_version() AS version")).version)
        elapsed = (time_ns() - started) / 1e9
        Topology.MoveGraphStore.close_store!(store)
        store = nothing
        edge_count > 0 || error("successful campaign produced no transitions")
        db_hash = _hash(temporary_db)
        payload = Dict{String,Any}(
        "schema_version" => 1,
        "campaign" => "Area-2 KnotInfo N<=6 one-step Reidemeister cache",
        "scope" => Dict("included" => ["281 official decorated representations", "R1/R2/R3 one-step candidates with target crossing count <=6", "verified transition cache"],
                        "excluded" => ["Area 2.5", "Area 3", "exhaustive classical-knot/isotopy graph", "KnotInfo-to-state association"]),
        "input" => Dict("campaign_path" => relpath(campaign_path, ROOT), "campaign_sha256" => _hash(campaign_path),
                        "knotinfo_path" => relpath(knotinfo_path, ROOT), "knotinfo_sha256" => _hash(knotinfo_path),
                        "campaign_input_snapshot_sha256" => String(campaign["input"]["snapshot_sha256"])),
        "versions" => Dict("canonicalization" => 1, "verifier" => 1, "schema" => 1),
        "counts" => Dict("seed_representations" => length(seeds), "unique_seed_states" => length(seed_ids),
                         "seed_dedupes" => seed_dedupes, "state_count" => state_count,
                         "candidate_count" => candidates, "over_limit_rejections" => over_limit,
                         "persisted_edge_count" => edge_count, "persist_successes" => persisted,
                         "rejection_count" => rejected, "error_counts" => errors,
                         "knotinfo_rows_imported" => imported_knotinfo_rows),
        "runtime" => Dict("julia_version" => string(VERSION), "sqlite_version" => sqlite_ver,
                          "command" => COMMAND, "elapsed_seconds" => elapsed,
                          "throughput_candidates_per_second" => candidates / max(elapsed, eps()),
                          "generated_at_utc" => string(Dates.now(Dates.UTC))),
        "database" => Dict("path" => relpath(db_path, ROOT), "sha256" => db_hash, "bytes" => filesize(temporary_db)),
        "oracle_boundary" => "KnotInfo names/invariants are oracle metadata; only independently verified Reidemeister paths are certificates.")
        open(temporary_json, "w") do io JSON.print(io, payload, 2); println(io) end
        open(temporary_markdown, "w") do io
        println(io, "# Knot n<=6 move database")
        println(io, "\n- Campaign input SHA-256: `$(payload["input"]["campaign_sha256"])`\n- KnotInfo input SHA-256: `$(payload["input"]["knotinfo_sha256"])`")
        println(io, "- Database: `$(payload["database"]["path"])`, SHA-256 `$(db_hash)`, $(payload["database"]["bytes"]) bytes")
        println(io, "\n## Scope\n\nExactly 281 official decorated representations were seeded and deduplicated. Every applicable one-step R1/R2/R3 candidate with target crossing count <=6 was independently verified and persisted. This is not an exhaustive knot or isotopy graph. Area 2.5 and Area 3 are excluded.")
        println(io, "\n## Counts\n\n```json\n", JSON.json(payload["counts"], 2), "\n```")
            println(io, "\n## Reproduction\n\n```text\n$COMMAND\n```")
        end
        mv(temporary_db, db_path; force=true)
        mv(temporary_json, json_path; force=true)
        mv(temporary_markdown, markdown_path; force=true)
        payload
    catch
        store !== nothing && isopen(store.db) && Topology.MoveGraphStore.close_store!(store)
        rm(temporary_db; force=true); rm(temporary_json; force=true); rm(temporary_markdown; force=true)
        rethrow()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    build()
end
