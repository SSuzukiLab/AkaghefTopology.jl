# cf: lib/topology.jl/docs/specs/api-contracts.md §§4.4–4.4.1, 5 — section-scoped reverse lookup.
# cf: tools/paper_db/document_references.py WIKILINK_RE — repository cite-key token syntax.
# Included inside APIIndex. These readers inspect Markdown, never Julia source.

const SPEC_CITATION_RE = r"\[\[([A-Za-z][A-Za-z0-9_-]{1,63})(?:#[^\]|\n]+)?(?:\|[^\]\n]+)?\]\]"

struct SpecSection
    level::Int
    anchor::String
    title::String
    citations::Set{String}
end

struct SpecPage
    sections::Vector{SpecSection}
    # Zero denotes an existing explicit anchor that does not identify a heading.
    anchors::Dict{String,Int}
end

"Walk visible Markdown content, excluding code examples and link destinations."
function markdown_refs!(@nospecialize(node), urls, citations, anchors, inline=false)
    if node isa AbstractString
        text = replace(node, r"<!--.*?-->"s => "")
        for m in eachmatch(SPEC_CITATION_RE, text)
            push!(citations, String(m.captures[1]))
        end
        for m in eachmatch(r"<a\s+(?:id|name)\s*=\s*[\"']([^\"']+)[\"'][^>]*>"i, text)
            push!(anchors, String(m.captures[1]))
        end
    elseif node isa Markdown.Code
        # Several existing specifications spell source wikilinks in inline code.
        # They still cite the key; a fenced/indented example does not.
        if inline
            for m in eachmatch(SPEC_CITATION_RE, node.code)
                push!(citations, String(m.captures[1]))
            end
        end
        return
    elseif node isa Markdown.LaTeX
        return
    elseif node isa Markdown.Link
        push!(urls, node.url)
        # Existing pages also link a registered Work's face note by its cite-key
        # filename. A PDF/source-file link alone does not name a cite key.
        face = match(r"(?:^|/)references/items/[^/]+/[^/]+/([^/#]+)\.md(?:#.*)?$", node.url)
        face === nothing || push!(citations, uri_decode(face.captures[1]))
        markdown_refs!(node.text, urls, citations, anchors, true)
    elseif node isa AbstractVector
        for child in node
            markdown_refs!(child, urls, citations, anchors, inline)
        end
    elseif parentmodule(typeof(node)) === Markdown
        for field in fieldnames(typeof(node))
            markdown_refs!(getfield(node, field), urls, citations, anchors,
                inline || node isa Markdown.Paragraph || node isa Markdown.Header || node isa Markdown.Table)
        end
    end
    nothing
end

function spec_catalog()
    pages = Dict{String,SpecPage}()
    for kind in ("specs", "conventions")
        directory = joinpath(ROOT, "docs", kind)
        for (dir, _, files) in walkdir(directory), name in files
            endswith(name, ".md") || continue
            path = joinpath(dir, name)
            pages[slash(relpath(path, ROOT))] = spec_page(read(path, String))
        end
    end
    pages
end

function markdown_links(md)
    urls = String[]
    markdown_refs!(md, urls, Set{String}(), String[])
    sort!(unique!(urls))
end

"Visible heading text, including inline code but excluding formatting markup."
function heading_text(@nospecialize(node))
    node isa AbstractString && return replace(node, r"<[^>]+>" => "")
    node isa Markdown.Code && return node.code
    node isa Markdown.LaTeX && return node.formula
    node isa AbstractVector && return join(heading_text(x) for x in node)
    hasproperty(node, :text) && return heading_text(node.text)
    hasproperty(node, :content) && return heading_text(node.content)
    ""
end

heading_slug(text) = replace(replace(lowercase(text), r"[^\p{L}\p{N}\p{M}_\-\s]" => ""), r"\s" => "-")

"Collect heading scopes in one pass; a parent includes its descendant sections."
function spec_page(raw)
    sections = SpecSection[]
    anchors = Dict{String,Int}()
    used = Set{String}()
    active = Int[]
    pending = String[]
    for block in Markdown.parse(raw).content
        citations = Set{String}()
        explicit = String[]
        markdown_refs!(block, String[], citations, explicit)
        if block isa Markdown.Header
            level = typeof(block).parameters[1]
            while !isempty(active) && sections[last(active)].level >= level
                pop!(active)
            end
            title = heading_text(block.text)
            base = heading_slug(title)
            slug = base
            suffix = 0
            while slug in used
                suffix += 1
                slug = base * "-" * string(suffix)
            end
            push!(used, slug)
            push!(sections, SpecSection(level, slug, title, Set{String}()))
            n = length(sections)
            anchors[slug] = n
            for alias in vcat(pending, explicit)
                anchors[alias] = n
            end
            empty!(pending)
            push!(active, n)
        else
            # An immediately preceding standalone HTML anchor names the heading.
            # Paragraph anchors have no section and supply no attribution (§4.4).
            for alias in explicit
                anchors[alias] = 0
            end
            rendered = block isa Markdown.Paragraph ? heading_text(block.content) : ""
            standalone = block isa Markdown.Paragraph && !isempty(explicit) && isempty(strip(rendered))
            if standalone
                append!(pending, explicit)
            else
                empty!(pending)
            end
        end
        for n in active
            union!(sections[n].citations, citations)
        end
    end
    SpecPage(sections, anchors)
end

function uri_decode(text)
    bytes = UInt8[]
    input = codeunits(text)
    i = 1
    while i <= length(input)
        if input[i] == UInt8('%') && i + 2 <= length(input)
            byte = tryparse(UInt8, String(input[i+1:i+2]); base=16)
            if byte !== nothing
                push!(bytes, byte)
                i += 3
                continue
            end
        end
        push!(bytes, input[i])
        i += 1
    end
    decoded = String(bytes)
    isvalid(decoded) ? decoded : String(text)
end

"Normalize a contract's docs-relative link to a library-relative page and anchor."
function spec_target(url; base=joinpath(ROOT, "docs", "api"))
    occursin(r"^(?:[A-Za-z][A-Za-z0-9+.-]*:|//|@ref)", url) && return nothing
    parts = split(url, '#'; limit=2)
    path = uri_decode(first(split(first(parts), '?'; limit=2)))
    anchor = length(parts) == 2 ? uri_decode(parts[2]) : ""
    full = startswith(path, "docs/") ? joinpath(ROOT, path) :
        startswith(path, "lib/topology.jl/docs/") ? joinpath(ROOT, "..", "..", path) : joinpath(base, path)
    relative = slash(relpath(normpath(full), ROOT))
    occursin(r"^docs/(?:specs|conventions)/.+\.md$", relative) || return nothing
    (page=relative, anchor=anchor)
end

"Resolve each page once per generation; never borrow citations from another section."
function spec_reference(target, cache)
    path = joinpath(ROOT, target.page)
    isfile(path) || return (error="page does not exist", clarification="", citations=String[])
    page = get!(cache, target.page) do
        spec_page(read(path, String))
    end
    isempty(target.anchor) && return (error="", clarification="", citations=String[])
    haskey(page.anchors, target.anchor) || return (error="anchor does not exist", clarification="", citations=String[])
    n = page.anchors[target.anchor]
    n == 0 && return (error="", clarification="節の無い錨", citations=String[])
    (error="", clarification="", citations=sort!(collect(page.sections[n].citations)))
end

"Readable section title, shared by a heading's canonical anchor and aliases."
function spec_title(target, cache)
    page = cache[target.page]
    if isempty(target.anchor)
        return isempty(page.sections) ? "Page" : first(page.sections).title * " (page)"
    end
    n = page.anchors[target.anchor]
    n == 0 ? "節の無い錨" : page.sections[n].title
end

spec_locator(target) = target.page * (isempty(target.anchor) ? "" : "#" * target.anchor)
spec_link(from, target) = link(from, joinpath(ROOT, target.page)) *
    (isempty(target.anchor) ? "" : "#" * replace(target.anchor, "%" => "%25", " " => "%20"))
