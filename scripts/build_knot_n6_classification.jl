"""Run the complete barrier-6 classical Gauss classification.

The run is intentionally stdout-summary-only.  Intermediate progress goes to
stderr so the final machine-readable rows remain easy to capture.  It does
not write a database or a census file and does not invoke the N5 insertion
audit.
"""

include(joinpath(@__DIR__, "..", "src", "knot_n6_classification.jl"))

const COMMAND =
    "julia --startup-file=no --pkgimages=no --project=. " *
    "scripts/build_knot_n6_classification.jl"

function _n6_script_parse_args(args)
    max_n = 6
    barrier = 6
    snapshot_path = _n6_default_snapshot()
    progress = true
    for arg in args
        if startswith(arg, "--max-n=")
            max_n = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--barrier=")
            barrier = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--snapshot=")
            snapshot_path = split(arg, "=", limit=2)[2]
        elseif arg == "--quiet"
            progress = false
        else
            throw(ArgumentError("unsupported argument $arg; use --max-n=N, " *
                               "--barrier=N, --snapshot=PATH, or --quiet"))
        end
    end
    max_n, barrier, snapshot_path, progress
end

function main(args=ARGS)
    max_n, barrier, snapshot_path, progress = _n6_script_parse_args(args)
    result = build_knot_n6_classification(max_n; barrier=barrier,
                                          snapshot_path=snapshot_path,
                                          progress=progress)
    summary = knot_n6_summary(result)
    println("classification_version=$(KNOT_N6_CLASSIFICATION_VERSION)")
    println("enumeration_version=$(CLASSICAL_GAUSS_ENUMERATION_VERSION)")
    println("scope=one_component;direct_classical_rooted_signed;barrier=$(result.barrier)")
    println("canonicalization=KnotGaussDiagram v$(KNOT_GAUSS_CANON_VERSION);" *
            "cyclic_root_shift+first_occurrence_relabel_only;mirror_and_reversal_distinct")
    println("known_label_convention=+ local KnotInfo PD representative; - mirror;" *
            " no right/left-handed assertion")
    println("n state_count raw_rooted_classical_count")
    for n in 0:result.max_n
        println(n, ' ', result.state_counts[n + 1], ' ',
                result.raw_rooted_classical_counts_by_n[n + 1])
    end
    println("total_states=$(length(result.states));raw_rooted_total=$(result.raw_rooted_total);" *
            "raw_rooted_classical_total=$(result.raw_rooted_classical_total)")
    println("edge_count=$(length(result.edges));direct_transition_count=$(result.direct_transition_count);" *
            "duplicate_edge_transition_count=$(result.duplicate_edge_transition_count);" *
            "nonclassical_target_count=$(result.nonclassical_target_count)")
    println("component_count=$(length(result.component_roots));known_component_count=" *
            "$(summary["known_component_count"]);known_representative_count=" *
            "$(length(result.known_representatives));unresolved_state_count=" *
            "$(summary["unresolved_state_count"]);unresolved_component_count=" *
            "$(length(result.unresolved_components))")
    println("unresolved_components=$(isempty(result.unresolved_components) ? "none" :
            join(result.unresolved_components, ','))")
    println("canonical_knot_type_counts=" *
            join(("$(label):$(result.canonical_knot_type_counts[label])" for
                  label in sort!(collect(keys(result.canonical_knot_type_counts)))), ','))
    println("raw_rooted_classical_type_counts=" *
            join(("$(label):$(result.raw_rooted_classical_type_counts[label])" for
                  label in sort!(collect(keys(result.raw_rooted_classical_type_counts)))), ','))
    println("raw_rooted_classical_unresolved_count=" *
            string(result.raw_rooted_classical_unresolved_count))
    for (component, label) in sort!(collect(result.component_type_labels))
        println("component_type=$(component):$(label)")
    end
    println("spanning_tree_edge_replay_count=$(result.spanning_tree_edge_replay_count);" *
            "spanning_tree_replay_failures=$(result.spanning_tree_replay_failures)")
    println("edge_invariant_sample_count=$(result.edge_invariant_sample_count);" *
            "edge_invariant_failures=$(length(result.edge_invariant_failures))")
    println("fingerprint_sha256=$(result.fingerprint)")
    println("barrier_status=$(isempty(result.unresolved_components) ?
            "barrier$(result.barrier)_resolved_for_known_components" :
            "barrier$(result.barrier)_unresolved;higher_barrier_not_run")")
    println("elapsed_seconds=$(result.elapsed_seconds)")
    println("command=$(COMMAND)")
    result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
