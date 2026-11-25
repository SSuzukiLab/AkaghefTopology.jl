# Read-only diagnostics of the generated newcomer entrypoints; not library acceptance.
# cf: lib/topology.jl/docs/specs/api-contracts.md §§4.4–4.4.1, 6 — reader navigation.
empty!(LOAD_PATH)
push!(LOAD_PATH, "@stdlib")
include(joinpath(@__DIR__, "api_index.jl"))
using Markdown

"Reassess the published tree without loading the library or rewriting its output."
function gaps_reading()
    root = APIIndex.ROOT
    out = APIIndex.OUT
    # Snapshot once: other lanes may edit sources while this read-only diagnostic runs.
    pages = Dict(joinpath(dir, f) => read(joinpath(dir, f), String)
        for (dir, _, fs) in walkdir(out) for f in fs if endswith(f, ".md"))
    gaps = pages[joinpath(out, "GAPS.md")]
    sections = Dict{String,Vector{String}}()
    active = ""
    for line in split(gaps, '\n')
        if startswith(line, "## ")
            active = String(strip(line[4:end]))
            sections[active] = String[]
        elseif startswith(line, "- ") && !isempty(active)
            push!(sections[active], line)
        end
    end
    for line in split(first(split(gaps, "\n## "; limit=2)), '\n')
        occursin(r": \d+$", line) && println(line)
    end
    for (heading, rows) in sort!(collect(sections); by=first)
        println("Heading ", heading, ": ", length(rows), " entries")
    end
    issues = String[]
    for (label, heading) in (("No contract", "No contract"), ("Incomplete", "Incomplete"),
            ("Complete", "Complete"), ("External re-exported bindings", "外部パッケージから再exportされた束縛"),
            ("Unassigned public functions", "どの対象にも属さない公開関数"),
            ("Unregistered or multiply registered exported types", "どの対象にも属さない export された型"),
            ("Exports not indexed", "索引に載せない export"),
            ("Items requiring specification clarification", "仕様が置き場所を定めていないもの"),
            ("Unresolved references", "解決できない参照"))
        m = match(Regex("(?m)^" * label * raw": (\d+)$"), gaps)
        m === nothing && (push!(issues, "Missing count: $label"); continue)
        parse(Int, m[1]) == length(get(sections, heading, String[])) || push!(issues, "Count/list mismatch: $label")
    end
    complete = Set{String}()
    for row in sections["Complete"]
        m = match(r"\]\(([^)]+)\)", row)
        m === nothing || push!(complete, normpath(joinpath(out, m[1])))
    end
    evidence = Dict(s => 0 for s in APIIndex.GUARANTEE_STATUSES)
    evidence[""] = 0
    methodrows = 0
    assessed_pages = 0
    for (p, raw) in sort!(collect(pages); by=first)
        basename(p) in ("INDEX.md", "GAPS.md", "SPECS.md", "SOURCES.md") && continue
        chunks = collect(eachmatch(r"(?ms)^## `([^\n]+)`\n\n(.*?)(?=^## `|^## Object indexes|\z)", raw))
        methodrows += length(chunks)
        occursin(r"(?m)^\[external:", raw) && continue
        assessed_pages += 1
        seen = Set{String}()
        for chunk in chunks
            body = String(rstrip(first(split(chunk[2], r"\n### (?:Inherited sources|Definition)\n"; limit=2))))
            body = replace(body, r"^\n+" => "") # Preserve the signature's four-space code indent.
            body == "[no contract]" && continue
            body = replace(body, r"(?m)^\[incomplete: .*\]\s*$" => "")
            body in seen && continue
            push!(seen, body)
            extension = occursin("(...) ->", chunk[1]) || occursin("/ext/", chunk[2])
            c = APIIndex.contract((raw=body,), chunk[1]; extension)
            line = count(==('\n'), SubString(raw, 1, prevind(raw, chunk.offset))) + 1
            locator = "$(APIIndex.slash(relpath(p, root))):$line"
            p in complete && c.status != "complete" && push!(issues,
                "$locator counted complete, now incomplete: $(join(c.missing, ", "))")
            for g in c.guarantees
                evidence[g.status] += 1
                isempty(g.status) && push!(issues,
                    "$locator guarantee $(g.number): missing or invalid status [$(g.annotation)]")
            end
        end
    end
    println("Read library-owned contract pages: ", assessed_pages, "; method rows including external: ", methodrows)
    println("Reassessed guarantees (identical rendered bodies counted once per callable): ", evidence)
    # Repeated rows must agree after resolving their relative contract links.
    repetitions = Dict{Tuple{String,String},Vector{Tuple{String,String,Int}}}()
    for (p, raw) in pages
        basename(p) == "INDEX.md" && dirname(p) != out || continue
        for (line, text) in enumerate(split(raw, '\n'))
            m = match(r"^- `(.+?)` — (.*) \[contract\]\(([^)]+)\)$", text)
            m === nothing && continue
            target = normpath(joinpath(dirname(p), m[3]))
            push!(get!(repetitions, (target, String(m[1])), Tuple{String,String,Int}[]), (String(m[2]), p, line))
        end
    end
    for ((target, signature), rows) in repetitions
        length(unique(first.(rows))) > 1 || continue
        length(unique(r[2] for r in rows)) > 1 || continue
        for (_, p, line) in rows
            push!(issues, "$(APIIndex.slash(relpath(p, root))):$line conflicting repeated row: $signature -> $(APIIndex.slash(relpath(target, root)))")
        end
    end
    foreach(println, sort!(unique(issues)))
    println("GAPS reading diagnostic: ", length(unique(issues)), " findings.")
end

function main()
    ARGS == ["--gaps"] && return gaps_reading()
    root = APIIndex.ROOT
    entrypoints = isempty(ARGS) ? [joinpath(root, "docs/api/SPECS.md"),
        joinpath(root, "docs/api/o-graph/INDEX.md")] : [joinpath(root, p) for p in ARGS]
    pages = copy(entrypoints)
    for p in entrypoints, url in APIIndex.markdown_links(Markdown.parse(read(p, String)))
        occursin(r"^[a-z]+:"i, url) && continue
        target = normpath(joinpath(dirname(p), first(split(url, '#'; limit=2))))
        startswith(APIIndex.slash(relpath(target, root)), "docs/api/") && isfile(target) && push!(pages, target)
    end
    cache = Dict{String,APIIndex.SpecPage}()
    issues = String[]
    for p in sort!(unique(pages))
        raw = read(p, String)
        lines = split(raw, '\n')
        for url in APIIndex.markdown_links(Markdown.parse(raw))
            occursin(r"^[a-z]+:"i, url) && continue
            parts = split(APIIndex.uri_decode(url), '#'; limit=2)
            target = isempty(first(parts)) ? p : normpath(joinpath(dirname(p), first(parts)))
            line = something(findfirst(s -> occursin("](" * url, s), lines), 0)
            reason = if startswith(url, "@ref")
                "unresolved docstring reference"
            elseif !isfile(target)
                "missing file"
            elseif length(parts) == 2 && endswith(target, ".md") && !startswith(parts[2], "^")
                page = get!(cache, target) do
                    APIIndex.spec_page(read(target, String))
                end
                haskey(page.anchors, parts[2]) ? "" : "missing heading/anchor"
            else
                ""
            end
            isempty(reason) || push!(issues, "$(APIIndex.slash(relpath(p, root))):$line — $url — $reason")
        end
    end
    foreach(println, sort!(issues))
    println("Navigation diagnostic: ", length(issues), " issues across ", length(unique(pages)), " generated pages.")
end
main()
