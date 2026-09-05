"""Area-2 combinatorial Gauss-code census and KnotInfo seed decoration run."""

using Dates
using JSON
using SHA

const KNOT_GAUSS_SOURCE = joinpath(@__DIR__, "..", "src", "knot_gauss.jl")
include(KNOT_GAUSS_SOURCE)

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260830"))
const KNOTINFO_DEFAULT = joinpath(PROJECT_ROOT, "data", "knotinfo-n6.json")
const CAMPAIGN_JSON_DEFAULT = joinpath(PROJECT_ROOT, "data", "knot-n6-campaign.json")
const CAMPAIGN_MD_DEFAULT = joinpath(PROJECT_ROOT, "reports", "knot-n6-campaign.md")
const CAMPAIGN_COMMAND = "julia --startup-file=no --pkgimages=no --project=lib/topology.jl lib/topology.jl/scripts/run_knot_n6_campaign.jl"

"""Parse the bracketed signed integer notation used by the KnotInfo snapshot."""
function parse_gauss_notation(notation)
    notation === nothing && return Int[]
    text = strip(String(notation))
    isempty(text) && return Int[]
    (startswith(text, "[") && endswith(text, "]")) ||
        throw(ArgumentError("Gauss notation must be a bracketed integer list: $text"))
    body = strip(text[2:end-1])
    isempty(body) && return Int[]
    tokens = strip.(split(body, ','))
    all(token -> occursin(r"^[+-]?\d+$", token), tokens) ||
        throw(ArgumentError("invalid Gauss notation token in $text"))
    parse.(Int, tokens)
end

function _snapshot_sha256(path::AbstractString)
    bytes2hex(SHA.sha256(read(path)))
end

function load_official_rows(path::AbstractString=KNOTINFO_DEFAULT)
    snapshot = JSON.parsefile(path)
    haskey(snapshot, "sheets") && !isempty(snapshot["sheets"]) ||
        throw(ArgumentError("KnotInfo snapshot has no sheet"))
    sheet = snapshot["sheets"][1]
    rows = sheet["rows"]
    names = String.(snapshot["included_knot_names"])
    length(rows) == length(names) == snapshot["included_row_count"] ||
        throw(ArgumentError("KnotInfo row-count metadata is inconsistent"))
    result = NamedTuple[]
    for (row, expected_name) in zip(rows, names)
        values = row["values"]
        name = String(values[1])
        name == expected_name || throw(ArgumentError("unexpected KnotInfo row order: $name"))
        crossing_number = Int(row["crossing_number"])
        code = parse_gauss_notation(values[23])
        length(code) == 2crossing_number ||
            throw(ArgumentError("$name has inconsistent Gauss-code length"))
        push!(result, (name=name, crossing_number=crossing_number, gauss=code,
                       source_row_index=Int(row["source_row_index"])))
    end
    result
end

function _orientation_from_index(n::Int, index::Int)
    [isodd(index >>> (label - 1)) ? -1 : 1 for label in 1:n]
end

"""Decorate one official ordinary Gauss code with all 2^n crossing signs."""
function decorate_official_row(row)
    n = row.crossing_number
    decorations = Any[]
    for index in 0:(2^n-1)
        orientation = _orientation_from_index(n, index)
        diagram = KnotGaussDiagram(row.gauss, orientation)
        push!(decorations, Dict{String,Any}(
            "orientation_index" => index,
            "orientation" => orientation,
            "canonical_serialization" => canonical_serialization(diagram),
            "canonical_hash" => canonical_hash(diagram),
        ))
    end
    decorations
end

function _stream_checksum(lines)
    context = SHA.SHA256_CTX()
    visited = 0
    for line in lines
        SHA.update!(context, codeunits(line))
        visited += 1
    end
    visited, bytes2hex(SHA.digest!(context))
end

"""Stream all rooted normalized signed codes for each n in the requested range."""
function stream_census(; n_min::Int=0, n_max::Int=6)
    0 <= n_min <= n_max <= 6 || throw(ArgumentError("census range must lie in 0:6"))
    records = Any[]
    for n in n_min:n_max
        started_ns = time_ns()
        lines = ("n=$(n);gauss=$(join(code, ','))\n"
                 for code in enumerate_rooted_normalized_signed_gc(n))
        visited, checksum = _stream_checksum(lines)
        elapsed_seconds = (time_ns() - started_ns) / 1.0e9
        expected = rooted_normalized_signed_gc_count(n)
        visited == expected || throw(AssertionError("census count mismatch at n=$n"))
        push!(records, Dict{String,Any}(
            "crossing_number" => n,
            "visited_count" => visited,
            "expected_count" => expected,
            "stream_checksum_sha256" => checksum,
            "elapsed_seconds" => elapsed_seconds,
            "codes_per_second" => visited / elapsed_seconds,
            "stream_encoding" => "UTF-8 lines n=<n>;gauss=<comma-separated signed labels>\\n",
        ))
    end
    records
end

function timed_stream_census(; n_min::Int=0, n_max::Int=6)
    total_started_ns = time_ns()
    records = stream_census(n_min=n_min, n_max=n_max)
    total_elapsed_seconds = (time_ns() - total_started_ns) / 1.0e9
    total_visited_count = sum(record["visited_count"] for record in records)
    Dict{String,Any}(
        "records" => records,
        "total_visited_count" => total_visited_count,
        "total_elapsed_seconds" => total_elapsed_seconds,
        "total_codes_per_second" => total_visited_count / total_elapsed_seconds,
        "timing_note" => "elapsed_seconds and codes_per_second are machine-specific and non-deterministic; counts and checksums are reproducible",
    )
end

function _campaign_payload(knotinfo_path::AbstractString)
    started_at_utc = string(Dates.now(Dates.UTC))
    rows = load_official_rows(knotinfo_path)
    snapshot = JSON.parsefile(knotinfo_path)
    stored_hash = String(snapshot["source"]["snapshot_sha256"])
    actual_hash = _snapshot_sha256(knotinfo_path)
    seeds = Any[]
    for row in rows
        decorations = decorate_official_row(row)
        push!(seeds, Dict{String,Any}(
            "name" => row.name,
            "crossing_number" => row.crossing_number,
            "gauss_notation" => isempty(row.gauss) ? "[]" : "[" * join(row.gauss, ", ") * "]",
            "decoration_count" => length(decorations),
            "decorations" => decorations,
        ))
    end
    census = timed_stream_census()
    finished_at_utc = string(Dates.now(Dates.UTC))
    Dict{String,Any}(
        "schema_version" => 1,
        "campaign" => "Area-2 single-component classical knot Gauss-code core",
        "scope" => Dict{String,Any}(
            "included" => ["representation validation", "rotation/relabel canonicalization",
                           "official-seed orientation decorations", "rooted signed GC census"],
            "excluded" => ["Area 2.5", "Area 3", "classical knot/isotopy census",
                           "Reidemeister connectivity run"],
        ),
        "claims" => Dict{String,Any}(
            "canonical_identity" => "rotation/relabel canonical hash is representation identity only",
            "knotinfo_role" => "KnotInfo labels and invariants are oracle/reference data, not isotopy proof",
            "isotopy_certificate" => "only an independently verified Reidemeister move path is an isotopy certificate",
            "move_verifier" => "R1/R2/R3 connectivity is deferred to integration/accepted move verifier",
            "decoration_limit" => "2^n orientation decorations are recorded representations; they are not claimed all realizable or all the same knot",
        ),
        "input" => Dict{String,Any}(
            "path" => relpath(knotinfo_path, pwd()),
            "snapshot_sha256" => actual_hash,
            "source_snapshot_sha256" => stored_hash,
            "included_names" => String.(snapshot["included_knot_names"]),
            "source_url" => String(snapshot["source"]["url"]),
            "source_oracle_status" => String(snapshot["source"]["oracle_status"]),
        ),
        "canonicalization" => Dict{String,Any}(
            "version" => KNOT_GAUSS_CANON_VERSION,
            "serialization" => "canonical_serialization(KnotGaussDiagram)",
            "hash" => "SHA-256 of UTF-8 canonical serialization",
        ),
        "runtime" => Dict{String,Any}(
            "command" => CAMPAIGN_COMMAND,
            "julia_version" => string(VERSION),
            "os" => string(Sys.KERNEL),
            "machine" => string(Sys.MACHINE),
            "run_date" => string(Dates.today()),
            "started_at_utc" => started_at_utc,
            "finished_at_utc" => finished_at_utc,
        ),
        "official_seed_decorations" => seeds,
        "rooted_signed_gc_census" => census,
    )
end

function write_campaign_outputs(; knotinfo_path::AbstractString=KNOTINFO_DEFAULT,
                                json_path::AbstractString=CAMPAIGN_JSON_DEFAULT,
                                markdown_path::AbstractString=CAMPAIGN_MD_DEFAULT)
    payload = _campaign_payload(knotinfo_path)
    mkpath(dirname(json_path)); mkpath(dirname(markdown_path))
    open(json_path, "w") do io
        JSON.print(io, payload, 2)
        println(io)
    end
    open(markdown_path, "w") do io
        _write_markdown(io, payload)
    end
    payload
end

function _write_markdown(io, payload)
    input = payload["input"]
    runtime = payload["runtime"]
    run_date = runtime["run_date"]
    canon_version = payload["canonicalization"]["version"]
    input_path = input["path"]
    input_hash = input["snapshot_sha256"]
    source_snapshot_hash = input["source_snapshot_sha256"]
    source_url = input["source_url"]
    included_names = input["included_names"]
    census = payload["rooted_signed_gc_census"]
    total_count = census["total_visited_count"]
    total_elapsed = census["total_elapsed_seconds"]
    total_throughput = census["total_codes_per_second"]
    julia_version = runtime["julia_version"]
    os = runtime["os"]
    machine = runtime["machine"]
    println(io, "# Knot n<=6 Area-2 campaign record")
    println(io)
    println(io, "Date: $run_date  ")
    println(io, "Agent: Codex  ")
    println(io, "Canonicalization: v$canon_version  ")
    println(io)
    println(io, "## Scope")
    println(io)
    println(io, "This is an Area-2 representation computation: validated one-component signed Gauss codes, rotation/relabel canonical identity, all crossing-orientation decorations of eight local KnotInfo seeds, and a complete rooted normalized signed-GC census for n=0..6.")
    println(io)
    println(io, "Area 2.5 and Area 3 are explicitly excluded. The census is combinatorial and is not a classical knot or isotopy census.")
    println(io)
    println(io, "## Input and provenance")
    println(io)
    println(io, "- Input: `$input_path`")
    println(io, "- Input snapshot SHA-256: `$input_hash`")
    println(io, "- Upstream source snapshot SHA-256 (metadata): `$source_snapshot_hash`")
    println(io, "- Source URL: $source_url")
    println(io, "- Official rows: `$(join(included_names, ", "))`")
    println(io)
    println(io, "The official ordinary Gauss column does not encode crossing signs. Each row therefore has exactly 2^n recorded orientation decorations. These are representations only; no claim is made that all are realizable classical diagrams or that they represent the same knot.")
    println(io)
    println(io, "## Rooted signed-GC census")
    println(io)
    println(io, "| n | visited | expected | elapsed s | codes/s | stream SHA-256 |")
    println(io, "|---:|---:|---:|---:|---:|---|")
    for record in census["records"]
        n = record["crossing_number"]
        visited = record["visited_count"]
        expected = record["expected_count"]
        elapsed = round(record["elapsed_seconds"], digits=6)
        throughput = round(record["codes_per_second"], digits=2)
        checksum = record["stream_checksum_sha256"]
        println(io, "| $n | $visited | $expected | $elapsed | $throughput | `$checksum` |")
    end
    println(io)
    println(io, "Total census: $total_count codes in $(round(total_elapsed, digits=6)) s ($(round(total_throughput, digits=2)) codes/s). Timings are machine-specific and non-deterministic; counts and checksums are reproducible.")
    println(io)
    println(io, "## Identification limits")
    println(io)
    println(io, "The rotation/relabel canonical hash identifies only the chosen Gauss-code representation orbit. KnotInfo names and invariants are oracle/reference data, not proof of isotopy. Only an independently verified R1/R2/R3 path can serve as an isotopy certificate; the R1/R2/R3 connectivity run is deferred to integration with the accepted move verifier.")
    println(io)
    println(io, "## Reproduction")
    println(io)
    println(io, "```text")
    println(io, runtime["command"])
    println(io, "```")
    println(io)
    println(io, "Machine/runtime: Julia $julia_version, $os, $machine. Machine-readable seed decorations and checksums are in `threads/PJ260830/data/knot-n6-campaign.json`.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_campaign_outputs()
end
