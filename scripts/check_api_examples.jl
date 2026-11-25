# cf: lib/topology.jl/docs/specs/api-contracts.md §5 — runs the `jldoctest`
# blocks of the contracts and compares each printed value with the run.
#
# Run from the repository root through the warm server:
#   sh lib/topology.jl/test/client.sh --exec lib/topology.jl/scripts/check_api_examples.jl [filter...]
# A plain filter selects the bindings whose qualified name contains it; a filter
# written with a leading `!` excludes them. Without filters every contract is run
# (over ten minutes for the whole library: use the long tier).

using Topology

module APIExamples

"Child modules of `root`, `root` first."
function modules(root::Module)
    found = Module[]
    function visit(m)
        m in found && return
        push!(found, m)
        for n in names(m; all=true, imported=false)
            isdefined(m, n) || continue
            v = getfield(m, n)
            v isa Module && v !== m && parentmodule(v) === m && visit(v)
        end
    end
    visit(root)
    found
end

blocks(text) = [String(m.captures[1]) for m in eachmatch(r"```jldoctest[^\n]*\n(.*?)```"s, text)]

"Split a block into (input, expected output) pairs at the `julia> ` prompts."
function steps(block)
    result = Tuple{String,String}[]
    lines = split(block, '\n')
    i = 1
    while i <= length(lines)
        if startswith(lines[i], "julia> ")
            input = String(lines[i][8:end])
            i += 1
            while i <= length(lines) && startswith(lines[i], "       ")
                input *= "\n" * lines[i][8:end]
                i += 1
            end
            expected = String[]
            while i <= length(lines) && !startswith(lines[i], "julia> ")
                push!(expected, lines[i])
                i += 1
            end
            push!(result, (input, String(strip(join(expected, "\n")))))
        else
            i += 1
        end
    end
    result
end

"Run one block in a fresh module; return the mismatches."
function run_block(block, root::Module)
    sandbox = Module(:APIExample)
    Core.eval(sandbox, :(const $(nameof(root)) = $root))
    Core.eval(sandbox, :(using .$(nameof(root))))
    failures = String[]
    for (input, expected) in steps(block)
        actual = try
            value = Core.eval(sandbox, Meta.parseall(input))
            (endswith(rstrip(input), ";") || value === nothing) ? "" :
                sprint(show, MIME"text/plain"(), value; context=:limit => true)
        catch err
            "ERROR: " * sprint(showerror, err)
        end
        actual = String(strip(actual))
        actual == expected || push!(failures,
            "julia> " * input * "\n    written: " * replace(expected, "\n" => "\n             ") *
            "\n    run:     " * replace(actual, "\n" => "\n             "))
    end
    failures
end

"""
Run the examples of the bindings selected by `filters`: a binding is run when its
qualified name contains one of the plain filters (or there is none) and none of
the filters written with a leading `!`.
"""
function check(root::Module; filters::Vector{String}=String[])
    wanted = [f for f in filters if !startswith(f, "!")]
    unwanted = [f[2:end] for f in filters if startswith(f, "!")]
    selected(name) = (isempty(wanted) || any(f -> occursin(f, name), wanted)) &&
        !any(f -> occursin(f, name), unwanted)
    nblocks = 0
    failed = Tuple{String,Vector{String}}[]
    seen = Base.IdSet{Any}()
    for m in sort!(modules(root); by=string)
        meta = Base.Docs.meta(m; autoinit=false)
        meta === nothing && continue
        bindings = sort!([b for b in keys(meta) if b isa Base.Docs.Binding]; by=string)
        for binding in bindings
            name = string(binding)
            selected(name) || continue
            multi = meta[binding]
            multi isa Base.Docs.MultiDoc || continue
            for sig in multi.order
                doc = multi.docs[sig]
                doc in seen && continue
                push!(seen, doc)
                text = sprint(io -> foreach(part -> print(io, part), doc.text))
                for block in blocks(text)
                    nblocks += 1
                    failures = run_block(block, root)
                    isempty(failures) || push!(failed, (name, failures))
                end
            end
        end
    end
    println("contract examples: ", nblocks, " blocks, ", length(failed), " with a mismatch")
    for (name, failures) in failed
        println("\n== ", name)
        foreach(f -> println("  ", f), failures)
    end
    isempty(failed) || error(string(length(failed), " contract example blocks do not match the run"))
    nothing
end

end # module APIExamples

APIExamples.check(Topology; filters=String.(ARGS))
