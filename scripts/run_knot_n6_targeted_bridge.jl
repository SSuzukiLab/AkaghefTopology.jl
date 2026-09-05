"""Run the unresolved barrier-6 components through the barrier-7 bridge.

The default invocation builds the N<=6 classification exactly once in memory,
then runs the targeted closure.  It writes only a deterministic stdout
summary; no N7 census, path archive, or database is produced.
"""

include(joinpath(@__DIR__, "..", "src", "knot_n6_targeted_bridge.jl"))

const COMMAND =
    "julia --startup-file=no --pkgimages=no --project=. " *
    "scripts/run_knot_n6_targeted_bridge.jl"

function _n6tb_script_parse_args(args)
    max_n = 6
    progress = true
    for arg in args
        if startswith(arg, "--max-n=")
            max_n = parse(Int, split(arg, "=", limit=2)[2])
        elseif arg == "--quiet"
            progress = false
        else
            throw(ArgumentError("unsupported argument $arg; use --max-n=N or --quiet"))
        end
    end
    max_n, progress
end

function _n6tb_print_counts(prefix::AbstractString, counts::Dict{String,Int})
    println(prefix * join(("$(label):$(counts[label])" for label in
                           sort!(collect(keys(counts)))), ','))
end

function main(args=ARGS)
    max_n, progress = _n6tb_script_parse_args(args)
    max_n == 6 || throw(ArgumentError(
        "the production targeted bridge requires --max-n=6"))
    started = time_ns()
    classification = build_knot_n6_classification(
        6; barrier=6, progress=progress)
    bridge = run_knot_n6_targeted_bridge(classification; barrier=7,
                                         progress=progress)
    summary = knot_n6_targeted_bridge_summary(bridge)
    println("bridge_version=$(KNOT_N6_TARGETED_BRIDGE_VERSION)")
    println("scope=barrier6_unresolved_only;barrier=7;barrier8_not_run")
    println("n6_input_states=$(length(classification.states));n6_input_components=" *
            "$(length(classification.component_roots));n6_input_fingerprint=" *
            classification.fingerprint)
    println("seed_components=$(join(bridge.seed_components, ','));seed_state_count=" *
            "$(bridge.seed_state_count)")
    println("n7_state_count=$(bridge.n7_state_count);n7_processed_count=" *
            "$(bridge.n7_processed_count);max_frontier=$(bridge.max_frontier)")
    println("candidate_count=$(bridge.candidate_count);adopted_transition_count=" *
            "$(bridge.adopted_transition_count);duplicate_transition_count=" *
            "$(bridge.duplicate_transition_count);reused_n6_candidate_count=" *
            "$(bridge.reused_n6_candidate_count);excluded_over_barrier_count=" *
            "$(bridge.excluded_over_barrier_count)")
    println("transition_verification_failures=" *
            "$(bridge.transition_verification_failures);n6_membership_failures=" *
            "$(bridge.n6_membership_failures);classicality_checks=" *
            "$(bridge.classicality_checks);classicality_failures=" *
            "$(bridge.classicality_failures)")
    _n6tb_print_counts("resolved_component_type_counts=",
                       bridge.resolved_component_type_counts)
    _n6tb_print_counts("resolved_state_type_counts=",
                       bridge.resolved_state_type_counts)
    for component in sort!(collect(keys(bridge.connected_known_labels)))
        labels = bridge.connected_known_labels[component]
        println("component_connection=$(component):$(isempty(labels) ? "unresolved" :
                join(labels, ','))")
    end
    for component in sort!(collect(keys(bridge.known_reach_distance)))
        println("known_reach_distance=$(component):$(bridge.known_reach_distance[component])")
    end
    println("closure_groups=" * join((join(group, ',') for group in
                                      bridge.closure_groups), ';'))
    println("resolved_component_count=$(length(bridge.resolved_components));" *
            "unresolved_component_count=$(length(bridge.unresolved_components));" *
            "unresolved_components=$(isempty(bridge.unresolved_components) ? "none" :
                join(bridge.unresolved_components, ','))")
    println("resolved_components=" * join(("$(component):$(label)" for
                                            (component, label) in
                                            sort!(collect(bridge.resolved_components))), ','))
    println("fingerprint_sha256=$(bridge.fingerprint)")
    println("bridge_elapsed_seconds=$(bridge.elapsed_seconds);total_elapsed_seconds=" *
            string((time_ns() - started) / 1e9))
    println("status=$(isempty(bridge.unresolved_components) ?
            "barrier7_resolved_all_seed_components" :
            "barrier7_unresolved;barrier8_not_run")")
    println("command=$(COMMAND)")
    summary
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

