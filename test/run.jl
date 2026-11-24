# Tiered, module-selective test runner. Policy: test/TESTING.md.
#
#   julia --project=. test/run.jl                          # unit tier, all modules
#   julia --project=. test/run.jl heegaard.koda            # unit tests of one module
#   julia --project=. test/run.jl --tier integration heegaard.koda
#   julia --project=. test/run.jl --tier acceptance        # the daily full run
#   julia --project=. test/run.jl --list [--tier T] [modules...]
#   julia --project=. test/run.jl --check                  # manifest covers every test file
#
# A module argument matches a module tag exactly or as a dotted prefix
# ("heegaard" selects "heegaard.koda"), or a file path relative to test/.
# --tier T selects every file whose tier is at or below T (unit < integration <
# acceptance); --only T selects tier T alone. "manual" files run only when named
# by path. Timings are appended to the file named by TOPOLOGY_TEST_TIMINGS
# (default test/timings.tsv, not tracked).

using TOML, Dates

const TESTDIR = @__DIR__
const REPODIR = normpath(joinpath(TESTDIR, "..", "..", ".."))  # the research repository root
const ORDER = Dict("unit" => 1, "integration" => 2, "acceptance" => 3, "manual" => 4)

function parse_args(args)
    tier = get(ENV, "TOPOLOGY_TEST_TIER", "unit"); only = false
    list = false; check = false; selectors = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--tier" || a == "--only"
            only = a == "--only"; tier = args[i + 1]; i += 2; continue
        elseif a == "--list"; list = true
        elseif a == "--check"; check = true
        else push!(selectors, a)
        end
        i += 1
    end
    haskey(ORDER, tier) || error("unknown tier $(tier); use unit, integration, acceptance or manual")
    return (; tier, only, list, check, selectors)
end

module_match(sel, tag) = tag == sel || startswith(tag, sel * ".")

function selected(entry, opts)
    named = any(sel -> normpath(sel) == normpath(entry["file"]), opts.selectors)
    named && return true
    t = entry["tier"]
    t == "manual" && return false
    tier_ok = opts.only ? t == opts.tier : ORDER[t] <= ORDER[opts.tier]
    tier_ok || return false
    isempty(opts.selectors) && return true
    return any(sel -> any(tag -> module_match(sel, tag), entry["modules"]), opts.selectors)
end

function check_manifest(manifest)
    listed = Set(normpath(joinpath(TESTDIR, e["file"])) for e in manifest["test"])
    preludes = Set(normpath(joinpath(TESTDIR, s["prelude"])) for s in values(manifest["suites"]))
    infra = Set(normpath(joinpath(TESTDIR, f)) for f in ("run.jl", "serve.jl", "runtests.jl", "heegaard/runtests.jl"))
    problems = String[]
    for e in manifest["test"]
        isfile(joinpath(TESTDIR, e["file"])) || push!(problems, "listed but missing: $(e["file"])")
    end
    for dir in (TESTDIR, joinpath(TESTDIR, "heegaard"))
        for f in readdir(dir; join=true)
            endswith(f, ".jl") || continue
            p = normpath(f)
            p in listed || p in preludes || p in infra || push!(problems, "not in manifest: $(relpath(p, TESTDIR))")
        end
    end
    return problems
end

function main(args)
    opts = parse_args(args)
    manifest = TOML.parsefile(joinpath(TESTDIR, "manifest.toml"))
    if opts.check
        problems = check_manifest(manifest)
        foreach(println, problems)
        println(isempty(problems) ? "MANIFEST_OK" : "MANIFEST_PROBLEMS $(length(problems))")
        return isempty(problems) ? 0 : 1
    end
    entries = [e for e in manifest["test"] if selected(e, opts)]
    if opts.list
        for e in entries
            println(rpad(e["tier"], 12), rpad(join(e["modules"], ","), 44), e["file"])
        end
        println("SELECTED $(length(entries))")
        return 0
    end
    isempty(entries) && (println("SELECTED 0: nothing to run"); return 1)
    timings = get(ENV, "TOPOLOGY_TEST_TIMINGS", joinpath(TESTDIR, "timings.tsv"))
    suites = Dict{String,Module}()
    loaded = Dict{String,Set{String}}()   # files already included into each suite
    failed = String[]
    for e in entries
        sname = e["suite"]
        mod = get!(suites, sname) do
            m = Module(Symbol("TopologyTests_", sname))
            # Module() has no include/eval of its own; give it the ones a `module` block has.
            Core.eval(m, :(include(p::AbstractString) = Base.include($m, p)))
            Core.eval(m, :(include(f::Function, p::AbstractString) = Base.include(f, $m, p)))
            Core.eval(m, :(eval(x) = Core.eval($m, x)))
            Base.include(m, joinpath(TESTDIR, manifest["suites"][sname]["prelude"]))
            m
        end
        done = get!(loaded, sname, Set{String}())
        # `requires` names files whose definitions this file uses (for example the
        # run_* functions of examples/*.jl); they are included first, once per suite.
        for r in get(e, "requires", String[])
            r in done && continue
            println("REQUIRE ", r)
            cd(TESTDIR) do
                Base.include(mod, joinpath(TESTDIR, r))
            end
            push!(done, r)
        end
        e["file"] in done && continue
        push!(done, e["file"])
        path = joinpath(TESTDIR, e["file"])
        println("RUN ", e["tier"], " ", e["file"])
        t0 = time(); ok = true
        try
            # Pkg.test runs from test/; a few files read repository paths relative to
            # the repository root and say so with cwd = "repo" in the manifest.
            dir = get(e, "cwd", "test") == "repo" ? REPODIR : TESTDIR
            cd(dir) do
                Base.include(mod, path)
            end
        catch err
            ok = false; push!(failed, e["file"])
            showerror(stdout, err); println()
        end
        dt = round(time() - t0; digits=1)
        println("DONE ", ok ? "pass" : "FAIL", " ", dt, "s ", e["file"])
        open(timings, "a") do io
            println(io, join((Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"), e["tier"], e["file"], dt, ok ? "pass" : "fail"), '\t'))
        end
    end
    println(isempty(failed) ? "ALL_PASS $(length(entries))" : "FAILED $(length(failed)): " * join(failed, ", "))
    return isempty(failed) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
