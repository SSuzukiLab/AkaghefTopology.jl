"""Run the finite barrier-5 classical signed-Gauss classification.

The script writes only a deterministic stdout report.  It does not create or
modify the KnotInfo snapshot or any database.
"""

include(joinpath(@__DIR__, "..", "src", "knot_n5_classification.jl"))

function _n5_script_args(args)
    max_n = 5
    barrier = 5
    snapshot_path = _n5_default_snapshot()
    for argument in args
        if startswith(argument, "--max-n=")
            max_n = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--barrier=")
            barrier = parse(Int, split(argument, "=", limit=2)[2])
        elseif startswith(argument, "--snapshot=")
            snapshot_path = split(argument, "=", limit=2)[2]
        else
            throw(ArgumentError("unsupported argument $argument; use --max-n=N, " *
                               "--barrier=N, or --snapshot=PATH"))
        end
    end
    max_n, barrier, snapshot_path
end

function main(args=ARGS)
    max_n, barrier, snapshot_path = _n5_script_args(args)
    result = build_knot_n5_classification(max_n;
                                          barrier=barrier,
                                          snapshot_path=snapshot_path)
    summary = knot_n5_summary(result)
    println("classification_version=$(KNOT_N5_CLASSIFICATION_VERSION)")
    println("scope=one_component;rooted_normalized_signed;Carter_genus_zero")
    println("canonicalization=KnotGaussDiagram v$(KNOT_GAUSS_CANON_VERSION);" *
            "cyclic_root_shift+first_occurrence_relabel_only;mirror_and_reversal_distinct")
    println("barrier=$(result.barrier);max_n=$(result.max_n)")
    println("state_counts=$(join(result.state_counts, ','));total_states=$(length(result.states))")
    println("canonical_knot_type_counts=" *
            join(("$(type):$(result.canonical_knot_type_counts[type])" for
                  type in KNOT_N5_STRICT_KNOT_TYPES), ","))
    println("raw_rooted_total=$(result.raw_rooted_total);" *
            "raw_rooted_classical_total=$(result.raw_rooted_classical_total)")
    println("raw_rooted_classical_counts=" *
            join(("$(type):$(result.raw_rooted_classical_counts[type])" for
                  type in KNOT_N5_STRICT_KNOT_TYPES), ","))
    println("edge_count=$(length(result.edges));direct_transition_count=$(result.direct_transition_count);" *
            "duplicate_edge_transition_count=$(result.duplicate_edge_transition_count)")
    println("component_count=$(length(result.component_roots));known_representative_count=" *
            "$(length(result.known_representatives));known_component_count=" *
            "$(summary["known_component_count"])")
    println("known_reachable_state_count=$(summary["known_reachable_state_count"])")
    println("unresolved_component_count=$(length(result.unresolved_components))")
    if isempty(result.unresolved_components)
        println("unresolved_components=none")
        if length(result.component_roots) <= 8
            println("barrier_status=barrier$(result.barrier)_resolved_for_enumerated_states")
        else
            println("barrier_status=barrier$(result.barrier)_unresolved;" *
                    "component_count_exceeds_known_low_crossing_bound")
        end
    else
        println("unresolved_components=$(join(result.unresolved_components, ','))")
        println("barrier_status=barrier$(result.barrier)_unresolved;" *
                "do_not_label_as_distinct_knot;higher_barrier_required_for_resolution")
        for component_id in result.unresolved_components
            states_in_component = findall(==(component_id), result.components)
            println("unresolved_component=$(component_id);state_count=" *
                    "$(length(states_in_component));min_crossing=" *
                    "$(minimum(crossing_count(result.states[index]) for index in states_in_component));" *
                    "max_crossing=$(maximum(crossing_count(result.states[index]) for index in states_in_component));" *
                    "higher_barrier_required=true")
        end
    end
    if length(result.component_roots) > 8
        println("barrier_component_excess=true;expected_known_low_crossing_max=8")
    else
        println("barrier_component_excess=false;expected_known_low_crossing_max=8")
    end
    println("path_replay_count=$(result.path_replay_count);path_replay_failures=" *
            "$(result.path_replay_failures)")
    println("edge_invariant_checks=$(result.edge_invariant_checks);edge_invariant_failures=" *
            "$(length(result.edge_invariant_failures))")
    for invariant in result.component_invariants
        println("component=$(invariant.component);representative=$(invariant.representative);" *
                "alexander=$(knot_polynomial_string(invariant.alexander));" *
                "jones=$(knot_polynomial_string(invariant.jones))")
    end
    for representative in result.known_representatives
        println("known=$(representative.label);state=$(representative.state_index);" *
                "component=$(representative.component);alexander=" *
                "$(knot_polynomial_string(representative.alexander));jones=" *
                "$(knot_polynomial_string(representative.jones))")
    end
    println("fingerprint_sha256=$(result.fingerprint)")
    println("elapsed_seconds=$(result.elapsed_seconds)")
    println("command=julia --startup-file=no --pkgimages=no --project=. " *
            "scripts/build_knot_n5_classification.jl")
    result
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
