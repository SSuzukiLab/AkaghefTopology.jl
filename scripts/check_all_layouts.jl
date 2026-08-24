#!/usr/bin/env julia
#
# Run the layout property checks over every o-graph plot instance and emit both
# a human table and a machine-readable record.
#
#     julia --project=julia julia/scripts/check_all_layouts.jl [out.json]
#
# The JSON is the point. `qmd思想.md` requires that generated numbers are not
# transcribed into prose by hand -- documents and status pages should render
# this file rather than quote counts that go stale the moment a router is fixed.
# It is also the seed of the value-record layer in `ATLAS_VISION.md` section 4:
# one row per (object, diagram) with an evidence class attached.

using Topology
using JSON

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "examples", "VLExample.jl"))
include(joinpath(ROOT, "examples", "OgraphPlotSuite.jl"))

"""Every plotted instance, as `(group, name, link, bending)`."""
function instances()
    out = NamedTuple[]
    for (group, calls) in run_ograph_plot_suite(), call in calls
        push!(out, (group=group, name=call.name, link=call.link, bending=call.bending))
    end
    for call in run_vl_plot_calls()
        push!(out, (group="topology/Manifold/VLExample", name=call.name,
                    link=call.link, bending=call.bending))
    end
    out
end

_diag(d) = Dict("rule"=>String(d.rule), "severity"=>String(d.severity),
                "message"=>d.message, "edges"=>d.edges,
                "at"=>d.at === nothing ? nothing : [real(d.at), imag(d.at)])

function inspect(instance)
    record = Dict{String,Any}(
        "group"     => instance.group,
        "name"      => instance.name,
        "realOnly"  => instance.link.real_only,
        "weighted"  => instance.link.is_weighted,
        "crossings" => length(instance.link.orientation),
        "pinned"    => instance.bending !== nothing,
    )
    if isempty(instance.link.orientation)
        record["status"] = "skipped"
        record["reason"] = "no crossings"
        record["diagnostics"] = []
        return record
    end
    local layout
    try
        layout = orthogonal_layout(instance.link; bending_numbers=instance.bending)
    catch err
        record["status"] = "error"
        record["reason"] = sprint(showerror, err)
        record["diagnostics"] = []
        return record
    end
    diagnostics = check_layout(layout)
    failing = [String(r) for r in LAYOUT_RULES
               if any(d -> d.rule === r && d.severity === :error, diagnostics)]
    record["edges"]       = length(layout.edges)
    record["bendTotal"]   = sum(abs, layout.bending_numbers)
    record["status"]      = isempty(failing) ? "pass" : "fail"
    record["failing"]     = failing
    record["diagnostics"] = _diag.(diagnostics)
    record
end

function main(args)
    records = inspect.(instances())
    checked = filter(r -> r["status"] != "skipped", records)
    failing = filter(r -> r["status"] != "pass" && r["status"] != "skipped", checked)

    println("=== o-graph layout property check ===")
    println("checked $(length(checked)) instance(s), $(length(failing)) failing\n")
    if isempty(failing)
        println("all clean")
    else
        width = maximum(length(r["group"] * "/" * r["name"]) for r in failing)
        for r in failing
            label = rpad(r["group"] * "/" * r["name"], width)
            detail = r["status"] == "error" ? "ERROR " * r["reason"] : join(r["failing"], ",")
            println("  ", label, "  ", detail)
        end
    end

    payload = Dict(
        "generated"  => string(round(Int, time())),
        "checked"    => length(checked),
        "failing"    => length(failing),
        "rules"      => [String(r) for r in LAYOUT_RULES],
        "instances"  => records,
    )
    path = isempty(args) ? joinpath(ROOT, "artifacts", "layout_check.json") : args[1]
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, payload, 2)
    end
    println("\nwrote ", path)
    isempty(failing) ? 0 : 1
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && exit(main(ARGS))
