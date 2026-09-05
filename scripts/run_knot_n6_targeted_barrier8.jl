"""Run the unresolved barrier-7 components through targeted barrier 8.

The N<=6 classification and targeted barrier-7 bridge are each reconstructed
once.  The bridge callback supplies only its transient unresolved workspace;
the barrier-8 lane expands those ten roots on demand and never enumerates all
N8 states.  R1/R2 insertion targets above the current barrier are excluded;
all adopted transitions are independently applied and verified.
"""

using Dates
using JSON

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TOPOLOGY_ROOT = joinpath(ROOT, "lib", "topology.jl")
const COMMAND = "julia --startup-file=no --pkgimages=no --project=. " *
                "scripts/run_knot_n6_targeted_barrier8.jl"

include(joinpath(TOPOLOGY_ROOT, "src", "knot_n6_targeted_barrier8.jl"))

function _final_type_counts(classification, bridge7)
    labels = Dict{Int,String}(classification.component_type_labels)
    for (component, label) in bridge7.resolved_components
        haskey(labels, component) && labels[component] != label &&
            error("barrier-7 component $component has conflicting labels")
        labels[component] = label
    end
    canonical = Dict{String,Int}()
    raw = Dict{String,Int}()
    unresolved_raw = 0
    for (index, diagram) in enumerate(classification.states)
        crossings = crossing_count(diagram)
        multiplicity = crossings == 0 ? 1 :
            (2crossings) ÷ rotation_automorphism_count(diagram)
        component = classification.components[index]
        label = get(labels, component, nothing)
        if label === nothing
            unresolved_raw += multiplicity
        else
            canonical[label] = get(canonical, label, 0) + 1
            raw[label] = get(raw, label, 0) + multiplicity
        end
    end
    Dict{String,Any}(
        "final_canonical_knot_type_counts" => canonical,
        "final_raw_rooted_classical_type_counts" => raw,
        "final_raw_rooted_classical_unresolved_count" => unresolved_raw,
        "final_classified_canonical_total" => sum(values(canonical); init=0),
        "final_classified_raw_rooted_total" => sum(values(raw); init=0))
end

function _compact_counts(classification, bridge7)
    Dict{String,Any}(
        "n6_state_counts" => copy(classification.state_counts),
        "n6_total_states" => length(classification.states),
        "n6_component_count" => length(classification.component_roots),
        "n6_unresolved_components_before" => copy(classification.unresolved_components),
        "n6_canonical_knot_type_counts" =>
            copy(classification.canonical_knot_type_counts),
        "n6_raw_rooted_classical_type_counts" =>
            copy(classification.raw_rooted_classical_type_counts),
        _final_type_counts(classification, bridge7)...,
        "barrier6_fingerprint" => classification.fingerprint,
        "barrier6_edge_count" => length(classification.edges),
        "barrier7_unresolved_components" => copy(bridge7.unresolved_components),
        "barrier7_n7_state_count" => bridge7.n7_state_count)
end

function _run(; progress::Bool=true)
    started = time_ns()
    started_at_utc = string(Dates.now(Dates.UTC))
    println(stderr, "[knot-n6-barrier8] reconstructing barrier-6 classification")
    classification = build_knot_n6_classification(
        6; barrier=6, progress=progress)
    workspace_ref = Ref{Any}(nothing)
    println(stderr, "[knot-n6-barrier8] reconstructing targeted barrier-7 bridge")
    bridge7 = run_knot_n6_targeted_bridge(
        classification; barrier=7, progress=progress,
        _workspace_callback=workspace -> (workspace_ref[] = workspace))
    println(stderr, "[knot-n6-barrier8] expanding unresolved barrier-7 roots: " *
        join(bridge7.unresolved_components, ","))
    result = run_knot_n6_targeted_barrier8(
        classification, bridge7, workspace_ref[]; barrier=8,
        progress=progress)
    finished_at_utc = string(Dates.now(Dates.UTC))
    payload = knot_n6_targeted_barrier8_summary(result)
    merge!(payload, _compact_counts(classification, bridge7))
    payload["barrier7_summary"] = knot_n6_targeted_bridge_summary(bridge7)
    payload["runtime"] = Dict{String,Any}(
        "command" => COMMAND,
        "julia_version" => string(VERSION),
        "os" => string(Sys.KERNEL),
        "machine" => string(Sys.MACHINE),
        "started_at_utc" => started_at_utc,
        "finished_at_utc" => finished_at_utc,
        "elapsed_seconds_total" => (time_ns() - started) / 1.0e9)
    payload["scope"] = Dict{String,Any}(
        "input" => "ten unresolved barrier-7 components only",
        "barrier6_reconstruction_count" => 1,
        "barrier7_reconstruction_count" => 1,
        "barrier8" => 8,
        "higher_barrier" => "not implemented",
        "known_boundary_policy" => "record boundary and do not expand known side",
        "insertion_completeness" =>
            "R1/R2 insertions with target <=8 are enumerated; reverse deletion edges supply the <=7 graph",
        "path_policy" => "transient shortest reach only; no full paths persisted")
    println(JSON.json(payload, 2))
    payload
end

if abspath(PROGRAM_FILE) == @__FILE__
    _run()
end
