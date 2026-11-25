# Diagnostic checks of the generator's reflection and formatting helpers.
# These fixtures are not library acceptance evidence and never write docs/api.
# cf: lib/topology.jl/docs/specs/api-contracts.md §§2–4.7 — cases exercised below.
empty!(LOAD_PATH)
push!(LOAD_PATH, "@stdlib")
using Markdown, TOML, Test
include(joinpath(@__DIR__, "api_index.jl"))

module ExternalFixture
module Values
"""
    Foreign() -> Foreign

Return a foreign-package value.

See [`foreign_missing`](@ref) and [missing display](@ref nowhere), or
[`Thing`](@ref).

# Guarantees
1. This foreign guarantee is not graded by Topology. [open: external scope]
"""
struct Foreign end
"Return a foreign-package result."
foreign_function(x::Foreign) = x
end
end

module APIFixture
using ..ExternalFixture.Values: Foreign, foreign_function
export Foreign, ForeignAlias, foreign_function
const ForeignAlias = Foreign
export Thing, Other, sample, sample_alias, nodoc, unrelated, ambiguous, majority,
    optional_vote, declaration, documented_declaration, missing_return, long_signature,
    external, union_placement, open_contract, unstated, return_only, mixed,
    source_vote, source_tie, external_untyped, CONSTANT, @fixture
"""A fixture type.

cf:[[Fixture]] §1 — diagnostic type scope
"""
struct Thing end
Base.show(io::IO, ::Thing) = print(io, "Thing()") # Added method, not an exported binding.
"""
    Other() -> Other

Return the other diagnostic object.

# Arguments
No arguments.

# Returns
`Other`.

# Guarantees
1. The result is an `Other`. [proved/agent: by construction]

# Examples
```jldoctest
julia> 1 + 1
2
```

Role: convert
"""
struct Other end
"""
    sample(x::Thing) -> Other

Return the fixture result. This sentence stays out of the index.

# Arguments
`x::Thing` is a fixture.
See [the cited section](../specs/ograph-colour-invariants.md#closed-examples),
[the page](../specs/ograph-colour-invariants.md),
[a paragraph anchor](../specs/ograph-colour-invariants.md#fixture-paragraph),
[a missing section](../specs/ograph-colour-invariants.md#fixture-missing), and
[a missing page](../specs/api-fixture-missing.md#section).

# Returns
`Other`.

# Guarantees
1. The result is an `Other`. [proved/agent: by construction]

# Examples
```jldoctest
julia> 1 + 1
2
```

Role: convert
"""
sample(x::Thing) = Other()
const sample_alias = sample
const CONSTANT = 1
macro fixture(x)
    esc(x)
end
nodoc(x::Thing) = x
unrelated(x) = x
"""
    return_only() -> Other

Return the other registered type without a typed argument.
"""
return_only() = Other()
mixed(x::Thing) = x
mixed(x) = x
source_vote(x::Int) = x
source_tie(x::Int) = x
ambiguous(x::Thing) = x
ambiguous(x::Other) = x
majority(x::Thing) = x
majority(x::Other) = x
majority(x::Other, y::Int) = x
"""
    optional_vote(x::Other, y=1, z=2) -> Other

Return the other object.
"""
optional_vote(x::Other, y=1, z=2) = x
optional_vote(x::Thing) = x
function declaration end
"""
    documented_declaration(x::Thing) -> Thing

Declare a future operation.

# Arguments
`x` is a fixture object.

# Returns
A fixture object.

# Guarantees
1. The extension returns its argument. [unverified: fixture extension has not been loaded]

# Examples
```julia
using FixtureExtension
documented_declaration(Thing())
```

Role: access
"""
function documented_declaration end
"""
    missing_return(x::Thing)

Return the input without documenting its type.
"""
missing_return(x::Thing) = x
"""
    long_signature(a::APIFixture.Thing, b::Tuple{Int,Int}, c, d, e, f; option::Bool=true, another=1) -> APIFixture.Other

Return another fixture object.
"""
long_signature(a::Thing, b, c, d, e, f; option=true, another=1) = Other()
union_placement(x::Union{Thing,Other}) = x
"""
    open_contract(x::Thing, y=1) -> Thing

Return a diagnostic object with unsettled interpretations.

# Arguments
`x` is a diagnostic object.

# Returns
The same object.

# Guarantees
3. An unsettled interpretation. [open: missing identification]
4. An inherited interpretation. [unverified: missing source comparison]
5. A source-backed diagnostic claim. [proved/literature: [[Fixture#^one]] §1]

# Examples
```jldoctest
julia> 1 + 1
2
```

Role: convert
"""
open_contract(x::Thing, y=1) = x
"""
    unstated(x::Thing) -> Thing

Return an object with an unannotated guarantee.

# Arguments
`x` is a diagnostic object.

# Returns
The same object.

# Guarantees
1. The object is unchanged.

# Examples
```jldoctest
julia> 1 + 1
2
```

Role: access
"""
unstated(x::Thing) = x
end
# Synthetic external definition; no file is created at this locator.
Base.include_string(APIFixture, "external(x::Thing) = x", "/external/FixturePackage/src/external.jl")
Base.include_string(APIFixture, "external_untyped(x) = x", "/external/FixturePackage/src/external.jl")
# Reflected definition locations in two existing files; their contents stay intact.
Base.include_string(APIFixture, "source_vote(x::String) = x\nsource_vote(x::Float64) = x\nsource_tie(x::String) = x",
    joinpath(@__DIR__, "api_index.jl"))

@testset "generator helpers (diagnostic only)" begin
    mods = APIIndex.modules(APIFixture)
    @test mods == [APIFixture]
    ds = APIIndex.documents(APIFixture.sample, mods)
    @test length(ds) == 1
    method = only(methods(APIFixture.sample))
    d = APIIndex.choose_doc(method, ds)
    c = APIIndex.contract(d, APIIndex.reflected_signature(method))
    @test c.status == "complete"
    @test c.summary == "Return the fixture result."
    @test c.role == "convert"
    @test c.signature == "sample(x::Thing) -> Other"
    @test only(c.guarantees).status == "proved/agent"
    @test APIIndex.contract(only(APIIndex.documents(APIFixture.Other, mods)), "Other()"; constructor=true).role == "convert"
    @test APIIndex.contract(only(APIIndex.documents(APIFixture.Thing, mods)), "Thing()"; constructor=true).role == "construct"
    @test "Guarantee status" in APIIndex.contract(only(APIIndex.documents(APIFixture.unstated, mods)), "unstated(x)").missing
    @test APIIndex.contract(only(APIIndex.documents(APIFixture.open_contract, mods)), "open_contract(x)").status == "complete"
    @test APIIndex.evidence_annotation("claim [proved/literature: [[Fixture#^one]] §1]") ==
        ("proved/literature", "proved/literature: [[Fixture#^one]] §1")
    @test APIIndex.evidence_annotation("claim [open: missing] extra text") == ("", "")
    @test first(APIIndex.evidence_annotation("claim [invented: missing]")) == ""
    for annotation in ("proved/agent", "proved/agent:", "open:   ")
        @test first(APIIndex.evidence_annotation("claim [" * annotation * "]")) == ""
    end
    for heading in APIIndex.REQUIRED
        raw = replace(d.raw, Regex("(?ms)^# " * heading * "\\n.*?(?=^# |\\z)") => "# " * heading * "\n\n")
        @test heading in APIIndex.contract((raw=raw,), "sample(x)").missing
    end
    empty_subsection = replace(d.raw, "`Other`." => "## Fields\n\n<!-- not content -->")
    @test "Returns" in APIIndex.contract((raw=empty_subsection,), "sample(x)").missing
    # An example elsewhere in the docstring cannot fill an empty Examples section.
    misplaced = replace(d.raw, "# Examples" => "# Sources") * "\n# Examples\n\nRole: convert\n"
    @test "Examples" in APIIndex.contract((raw=misplaced,), "sample(x)").missing
    for annotation in ("proved/agent", "proved/agent: ")
        raw = replace(d.raw, "proved/agent: by construction" => annotation)
        @test "Guarantee status" in APIIndex.contract((raw=raw,), "sample(x)").missing
    end
    decl = only(APIIndex.documents(APIFixture.documented_declaration, mods))
    @test APIIndex.contract(decl, "documented_declaration(...)"; extension=true).status == "complete"
    @test APIIndex.contract(decl, "documented_declaration(...)").status == "incomplete"
    @test APIIndex.contract((raw=replace(decl.raw, "```julia" => "```jldoctest"),),
        "documented_declaration(...)"; extension=true).status == "incomplete"
    @test APIIndex.signature_arguments(APIFixture, "f(x, y::Vector{Thing}=Thing[]; z::Other) -> Other", APIFixture) == Any[Any, Vector{APIFixture.Thing}]
    parsed = APIIndex.guarantees(Markdown.parse("# Guarantees\n\n2. Missing.\n\n   1. Nested [open: not a separate guarantee]\n\n# Examples\n\n1. Not a guarantee.\n"))
    @test length(parsed) == 1
    @test only(parsed).number == 2
    @test isempty(only(parsed).status)
    @test isempty(APIIndex.guarantees(Markdown.parse("# Guarantees\n\n- Not numbered. [open: no]\n")))
    # Section detection must agree with contract()'s accepted Markdown headers.
    @test isempty(only(APIIndex.guarantees(Markdown.parse("## Guarantees\n\n1. Missing.\n"))).status)
    @test isempty(only(APIIndex.guarantees(Markdown.parse("# Guarantees\n\n## Detail\n\n1. Missing evidence. [proved/agent]\n\n# Examples\n"))).status)
    @test APIIndex.return_type(APIFixture, c.returns, APIFixture) === APIFixture.Other
    @test APIIndex.contract(nothing, "nodoc(x::Thing)").status == "no contract"
    @test APIIndex.contract(nothing, "nodoc(x::Thing)").signature == "nodoc(x::Thing) -> ?"
    row = (contract=c, doc=d)
    @test isempty(APIIndex.contract_conflicts([row, row]))
    alternate = merge(row, (doc=merge(d, (raw=replace(d.raw, "Return the fixture result." => "Return an alternate description."),)),))
    @test APIIndex.contract_conflicts([row, alternate]) == [c.signature]
    @test APIIndex.short_signature("f(a::A.B.T, b::Tuple{A.B.T,Int}, c, d, e; kw=1) -> A.B.R") ==
        "f(a::T, b::Tuple{T, Int}, c, d, …; …) -> R"
    @test APIIndex.short_signature("f(a=1, b=2) -> Int") == "f(a = 1, b = 2) -> Int"
    @test APIIndex.short_signature("f(x::T) where {T<:A.B.R} -> T") == "f(x::T) -> T"
    @test APIIndex.short_signature("declaration(...) -> ?") == "declaration(...) -> ?"
    @test APIIndex.short_signature("APIFixture.Thing() -> ?") == "Thing() -> ?"
    @test APIIndex.short_signature("f(x::A.B.T=A.DEFAULT) -> A.R") == "f(x::T = A.DEFAULT) -> R"
    @test !APIIndex.mentions(Any, APIFixture.Thing)
    @test APIIndex.mentions(Vector{APIFixture.Thing}, APIFixture.Thing)
    @test APIIndex.mentions(Union{APIFixture.Thing,Nothing}, APIFixture.Thing)
    @test APIIndex.mentions(Type{APIFixture.Thing}, APIFixture.Thing)
    @test APIIndex.mentions(Number, Int)
    @test !APIIndex.mentions(Number, APIFixture.Thing)
    @test APIIndex.filename("a<b?") == "a%3Cb%3F.md"
    registered = TOML.parsefile(joinpath(@__DIR__, "..", "docs", "objects.toml"))
    @test all(key -> haskey(registered, key), ["casson-walker", "lens-space"])
    registry_gaps = String[]
    APIIndex.source_registry(registered, registry_gaps)
    @test isempty(registry_gaps)
    @test APIIndex.cf_lines("text\ncf:[[Fixture#^one]] §1 — scope\n") == ["cf:[[Fixture#^one]] §1 — scope"]
end

@testset "section reverse lookup (diagnostic only)" begin
    raw = """
    # Page
    [[RootKey]]

    <a id="parent-alias"></a>
    ## Parent `code`
    [[ParentKey#^source]]
    ### Child
    [[ChildKey]]
    A source written as `[[InlineKey#^source]]`.
    A Gauss display `GC [[[-1,-2,+2,+1]];[-,+]]` is not a citation.
    [Face note](../../../../references/items/article/w-fixture/FaceKey.md).
    ```text
    ## Not a heading
    [[ExampleKey]] [example](../specs/example.md)
    ```
    ## Sibling
    [[SiblingKey]]
    ## Parent `code`
    [[RepeatedKey]]
    <a id="inline"></a> An anchor in a paragraph.
    ## 日本語の節
    [[UnicodeKey]]
    """
    page = APIIndex.spec_page(raw)
    cites(anchor) = page.sections[page.anchors[anchor]].citations
    @test cites("parent-code") == Set(["ParentKey", "ChildKey", "InlineKey", "FaceKey"])
    @test cites("parent-alias") == cites("parent-code")
    @test page.sections[page.anchors["parent-alias"]].title == "Parent code"
    @test cites("child") == Set(["ChildKey", "InlineKey", "FaceKey"])
    @test cites("sibling") == Set(["SiblingKey"])
    @test cites("parent-code-1") == Set(["RepeatedKey"])
    @test cites("日本語の節") == Set(["UnicodeKey"])
    @test page.anchors["inline"] == 0
    @test !haskey(page.anchors, "not-a-heading")
    @test !("ExampleKey" in cites("page"))
    @test APIIndex.markdown_links(Markdown.parse(raw)) == ["../../../../references/items/article/w-fixture/FaceKey.md"]
    @test APIIndex.markdown_links(Markdown.parse("[link](../specs/a.md#one) and `[code](../specs/b.md)`")) == ["../specs/a.md#one"]
    @test APIIndex.spec_target("../specs/a.md#one") == (page="docs/specs/a.md", anchor="one")
    @test APIIndex.spec_target("docs/conventions/a.md#one") == (page="docs/conventions/a.md", anchor="one")
    @test APIIndex.spec_target("lib/topology.jl/docs/specs/a.md#one") == (page="docs/specs/a.md", anchor="one")
    @test APIIndex.spec_target("https://example.org/docs/specs/a.md#one") === nothing
    @test APIIndex.spec_target("../src/a.jl") === nothing
    @test APIIndex.uri_decode("%E6%97%A5%E6%9C%AC%E8%AA%9E") == "日本語"
    @test APIIndex.render_doc("[section](docs/specs/a.md#one)",
        joinpath(APIIndex.OUT, "fixture", "sample.md"), Dict{String,Vector{String}}(), String[], "sample") ==
        "[section](../../specs/a.md#one)"
    cache = Dict("docs/specs/ograph-colour-invariants.md" => page)
    target = (page="docs/specs/ograph-colour-invariants.md", anchor="parent-alias")
    @test APIIndex.spec_reference(target, cache).citations == ["ChildKey", "FaceKey", "InlineKey", "ParentKey"]
    @test isempty(APIIndex.spec_reference(merge(target, (anchor="",)), cache).citations)
    @test APIIndex.spec_reference(merge(target, (anchor="absent",)), cache).error == "anchor does not exist"
    @test APIIndex.spec_reference(merge(target, (anchor="inline",)), cache).clarification == "節の無い錨"
    @test isempty(APIIndex.spec_reference(merge(target, (anchor="inline",)), cache).citations)
    @test APIIndex.spec_title(target, cache) == "Parent code"
    @test APIIndex.spec_title(merge(target, (anchor="日本語の節",)), cache) == "日本語の節"
    @test APIIndex.spec_title(merge(target, (anchor="inline",)), cache) == "節の無い錨"
    @test APIIndex.spec_reference((page="docs/specs/api-fixture-missing.md", anchor=""), cache).error == "page does not exist"
end

@testset "source registry (diagnostic only)" begin
    gaps = String[]
    reg = Dict("one" => Dict("sources" => ["scripts/"]))
    owners = APIIndex.source_registry(reg, gaps)
    @test owners[joinpath(@__DIR__, "api_checks.jl")] == "one"
    @test isempty(gaps)
    reg["two"] = Dict("sources" => ["scripts/api_checks.jl"])
    owners = APIIndex.source_registry(reg, gaps)
    @test !haskey(owners, joinpath(@__DIR__, "api_checks.jl"))
    @test owners[joinpath(@__DIR__, "api_index.jl")] == "one"
    @test any(s -> occursin("belongs to both", s), gaps)
    gaps = String[]
    @test isempty(APIIndex.source_registry(Dict("one" => Dict("sources" => ["../escape.jl", "src/api_missing.jl"])), gaps))
    @test length(gaps) == 2
    @test isempty(APIIndex.source_registry(Dict("one" => Dict()), String[]))
end

let
    directory = joinpath(@__DIR__, ".api-check-output")
    registry = joinpath(@__DIR__, ".api-check-objects.toml")
    # Freeze real page inputs for determinism; substitute one in-memory fixture
    # to test aliases and paragraph anchors without editing a specification page.
    pages = APIIndex.spec_catalog()
    pages["docs/specs/ograph-colour-invariants.md"] = APIIndex.spec_page("""
    # Fixture references
    <a id="closed-examples"></a>
    ## Closed examples
    [[BP97]]
    ## Unlinked section
    [[UnlinkedKey]]
    <a id="fixture-paragraph"></a> A paragraph citing [[ParagraphKey]].
    """)
    write(registry, """
    [fixture]
    name = "Diagnostic fixture"
    types = ["Thing", "Foreign"]
    specs = []
    conventions = []
    [other]
    name = "Other diagnostic object"
    types = ["Other"]
    specs = []
    conventions = []
    """)
    APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    @test occursin("Return the fixture result.", read(joinpath(directory, "fixture", "INDEX.md"), String))
    @test !occursin("This sentence", read(joinpath(directory, "fixture", "INDEX.md"), String))
    @test occursin("inherited: type", read(joinpath(directory, "SOURCES.md"), String))
    sources = read(joinpath(directory, "SOURCES.md"), String)
    specs = read(joinpath(directory, "SPECS.md"), String)
    @test occursin("## BP97", sources)
    @test occursin("### 直接の出典", sources)
    @test occursin("### 仕様を介した依存", sources)
    @test occursin("via [docs/specs/ograph-colour-invariants.md#closed-examples]", sources)
    @test !occursin("via [docs/specs/ograph-colour-invariants.md]", sources)
    @test !occursin("ParagraphKey", sources)
    @test !occursin("UnlinkedKey", sources)
    @test occursin("## [docs/specs/ograph-colour-invariants.md]", specs)
    @test occursin("## [docs/specs/ograph-colour-invariants.md#closed-examples]", specs)
    @test occursin("#closed-examples) — Closed examples", specs)
    @test occursin("#fixture-paragraph) — 節の無い錨", specs)
    paragraph_section = split(split(specs, "#fixture-paragraph) — 節の無い錨"; limit=2)[2], "\n## "; limit=2)[1]
    @test occursin("sample(", paragraph_section)
    @test !occursin("fixture-missing", specs)
    @test isfile(joinpath(directory, "unassigned", "unrelated.md"))
    @test isfile(joinpath(directory, "fixture", "ambiguous.md")) # alphabetical tie
    @test isfile(joinpath(directory, "other", "majority.md"))
    @test isfile(joinpath(directory, "other", "optional_vote.md")) # count methods before doc dedup
    @test isfile(joinpath(directory, "unassigned", "declaration.md"))
    @test isfile(joinpath(directory, "fixture", "documented_declaration.md"))
    @test occursin("documented_declaration(...) -> Thing", read(joinpath(directory, "fixture", "INDEX.md"), String))
    @test !occursin("[incomplete", read(joinpath(directory, "fixture", "documented_declaration.md"), String))
    gaps = read(joinpath(directory, "GAPS.md"), String)
    @test occursin("External re-exported bindings: 2", gaps)
    @test !occursin("Base.show", gaps)
    @test occursin("## 外部パッケージから再exportされた束縛", gaps)
    @test occursin("ExternalFixture.Values.Foreign", gaps)
    @test occursin("[external: ExternalFixture]", read(joinpath(directory, "fixture", "INDEX.md"), String))
    foreign_page = read(joinpath(directory, "fixture", "Foreign.md"), String)
    @test occursin("APIFixture.ForeignAlias", foreign_page)
    @test !occursin("[incomplete", foreign_page)
    @test !occursin("[no contract]", foreign_page)
    @test occursin("See `foreign_missing` and `missing display`", foreign_page)
    @test occursin("[`Thing`](Thing.md)", foreign_page)
    @test !occursin("@ref", foreign_page)
    @test occursin("unresolved or ambiguous docstring reference `foreign_missing`", gaps)
    @test occursin("unresolved or ambiguous docstring reference `nowhere`", gaps)
    @test !occursin("external scope", split(gaps, "## Guarantees: open"; limit=2)[2])
    total = parse(Int, match(r"Public functions indexed[^\n]*: (\d+)", gaps)[1])
    assessed = sum(parse(Int, match(Regex("\n" * s * ": (\\d+)"), gaps)[1])
        for s in ("No contract", "Incomplete", "Complete"))
    @test total == assessed + 2
    @test occursin("## 解決できない参照", gaps)
    @test occursin("fixture-missing`: anchor does not exist", gaps)
    @test occursin("api-fixture-missing.md#section`: page does not exist", gaps)
    unresolved_section = split(gaps, "## 解決できない参照"; limit=2)[2]
    @test occursin("#fixture-paragraph`: 節の無い錨", unresolved_section)
    clarification_section = split(split(gaps, "## 仕様が置き場所を定めていないもの"; limit=2)[2], "## 解決できない参照"; limit=2)[1]
    @test !occursin("fixture-paragraph", clarification_section)
    @test !occursin("one method selects multiple primary objects", gaps)
    @test isfile(joinpath(directory, "fixture", "union_placement.md"))
    @test occursin("- proved/agent: 2", gaps)
    @test occursin("- proved/literature: 1", gaps)
    @test occursin("- open: 1", gaps)
    @test occursin("- unverified: 2", gaps)
    @test occursin("- missing or invalid status: 1", gaps)
    @test occursin("Guarantee 3 — [open: missing identification]", gaps)
    @test occursin("Guarantee 4 — [unverified: missing source comparison]", gaps)
    @test occursin("Guarantee status", read(joinpath(directory, "fixture", "unstated.md"), String))
    @test !occursin("no documented return type", gaps)
    @test !occursin("no repository-local source locator", gaps)
    @test occursin("[incomplete: return type", read(joinpath(directory, "fixture", "missing_return.md"), String))
    @test occursin("APIFixture/external.jl:1", read(joinpath(directory, "fixture", "external.md"), String))
    @test !occursin("[external:", read(joinpath(directory, "fixture", "external.md"), String))
    @test !occursin("/external/FixturePackage", read(joinpath(directory, "fixture", "external.md"), String))
    index = read(joinpath(directory, "fixture", "INDEX.md"), String)
    @test !occursin("[incomplete:", index)
    @test occursin("missing_return(x::Thing) -> ?", index)
    @test occursin("long_signature(a::Thing, b::Tuple{Int, Int}, c, d, …; …) -> Other", index)
    @test occursin("another=1", read(joinpath(directory, "fixture", "long_signature.md"), String))
    @test occursin("CONSTANT (constant)", read(joinpath(directory, "GAPS.md"), String))
    @test occursin("@fixture (macro)", read(joinpath(directory, "GAPS.md"), String))
    @test occursin("../fixture/sample.md", read(joinpath(directory, "other", "INDEX.md"), String))
    snapshot() = Dict(relpath(joinpath(dir, f), directory) => read(joinpath(dir, f))
        for (dir, _, fs) in walkdir(directory) for f in fs)
    before = snapshot()
    APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    @test snapshot() == before

    # A source fallback is per function, after all registered argument/return
    # types. No-method declarations and external definitions retain §4.7 placement.
    open(registry, "a") do io
        print(io, """
        sources = ["scripts/api_index.jl"]
        [fallback]
        name = "Source fallback"
        types = []
        specs = []
        conventions = []
        sources = ["scripts/api_checks.jl"]
        """)
    end
    elapsed_first = @elapsed APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    @test isfile(joinpath(directory, "fallback", "unrelated.md"))
    @test !isfile(joinpath(directory, "unassigned", "unrelated.md"))
    @test isfile(joinpath(directory, "other", "return_only.md"))
    @test isfile(joinpath(directory, "fixture", "mixed.md"))
    @test !occursin("mixed(", read(joinpath(directory, "fallback", "INDEX.md"), String))
    @test isfile(joinpath(directory, "other", "source_vote.md"))
    @test isfile(joinpath(directory, "fallback", "source_tie.md"))
    @test isfile(joinpath(directory, "unassigned", "declaration.md"))
    @test isfile(joinpath(directory, "fixture", "documented_declaration.md"))
    @test isfile(joinpath(directory, "unassigned", "external_untyped.md"))
    @test occursin(" — ", only(filter(s -> startswith(s, "## [docs/specs/ograph-colour-invariants.md#closed-examples]"),
        split(read(joinpath(directory, "SPECS.md"), String), '\n'))))
    before = snapshot()
    elapsed_second = @elapsed APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    @test snapshot() == before
    # A loaded extension is not a child of the root module. Admit only methods
    # defined under this library's ext/, retaining the declaration's contract.
    extension = Module(:FixtureExtension)
    Core.eval(extension, :(import Main.APIFixture))
    Base.include_string(extension,
        "APIFixture.documented_declaration(x::APIFixture.Thing) = x",
        joinpath(APIIndex.ROOT, "ext", "FixtureExtension.jl"))
    Base.include_string(extension,
        "APIFixture.documented_declaration(x::APIFixture.Other) = x",
        joinpath(APIIndex.ROOT, "external", "FixtureExtension.jl"))
    @test count(APIIndex.extension_method, methods(APIFixture.documented_declaration)) == 1
    APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    loaded = read(joinpath(directory, "fixture", "documented_declaration.md"), String)
    @test occursin("documented_declaration(x::Thing)", loaded)
    @test !occursin("[incomplete", loaded)
    @test !occursin("FixtureExtension.jl:1", read(joinpath(directory, "other", "INDEX.md"), String))
    @test !occursin("documented_declaration(", read(joinpath(directory, "other", "INDEX.md"), String))
    before = snapshot()
    APIIndex.generate(APIFixture; output=directory, registry_file=registry, spec_pages=pages)
    @test snapshot() == before
    println("Source-fallback diagnostic generations (s): ", round(elapsed_first; digits=3), ", ", round(elapsed_second; digits=3))
    println("Diagnostic output only: scripts/.api-check-output")
end

@testset "source scopes and inherited links (diagnostic only)" begin
    # This transient file is a file-head reader fixture, not Julia discovery.
    path = joinpath(@__DIR__, ".api-check-header.jl")
    write(path, "# cf:[[FileKey]]\n\"\"\"\ncf:[[CallableKey]]\n\"\"\"\nfunction example() end\n")
    @test APIIndex.file_sources(path) == ["cf:[[FileKey]]"]
    write(path, "# cf:[[FileKey]]\n\"\"\"\ncf:[[ModuleKey]]\n\"\"\"\nmodule Example\n\"\"\"\ncf:[[CallableKey]]\n\"\"\"\nfunction example() end\nend\n")
    @test APIIndex.file_sources(path) == ["cf:[[FileKey]]", "cf:[[ModuleKey]]"]
    rendered = APIIndex.render_doc("cf: [PS](../specs/virtual-link-moves.md#125-integral-ps)",
        joinpath(APIIndex.OUT, "o-graph", "example.md"), Dict{String,Vector{String}}(), String[], "fixture")
    @test occursin("](../../specs/virtual-link-moves.md#125-integral-ps)", rendered)
end
