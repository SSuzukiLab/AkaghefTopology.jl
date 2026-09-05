"""Bounded benchmark for the streaming all-state N=6 enumerator.

The default invocation examines at most one million raw N=6 inputs and at
most five minutes.  It never performs the complete 42,577,920-input N=6
scan.  No DB or census file is written.  A direct-classical canonical set can
be supplied when this script is included and the benchmark function is called
programmatically; the command-line benchmark intentionally leaves planarity
unknown because it must remain a bounded all-state measurement.

Run from `lib/topology.jl`:

    julia --startup-file=no --pkgimages=no --project=. \
        scripts/benchmark_knot_n6_all_state_enumeration.jl
"""

using Printf

if !isdefined(@__MODULE__, :KnotGaussDiagram)
    include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
end
if !isdefined(@__MODULE__, :stream_canonical_knot_gauss)
    include(joinpath(@__DIR__, "..", "src", "knot_n6_all_state_enumeration.jl"))
end

const KNOT_N6_BENCHMARK_RAW_LIMIT = 1_000_000
const KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS = 300.0

function _benchmark_parse_args(args)
    n = 6
    max_raw = KNOT_N6_BENCHMARK_RAW_LIMIT
    max_seconds = KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS
    for arg in args
        if startswith(arg, "--n=")
            n = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--max-raw=")
            max_raw = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--max-seconds=")
            max_seconds = parse(Float64, split(arg, "=", limit=2)[2])
        else
            throw(ArgumentError(
                "unsupported argument $arg; use --n=N --max-raw=M --max-seconds=S"))
        end
    end
    0 <= n <= 6 || throw(ArgumentError("n must lie in 0:6"))
    0 < max_raw <= KNOT_N6_BENCHMARK_RAW_LIMIT ||
        throw(ArgumentError("max_raw must lie in 1:$(KNOT_N6_BENCHMARK_RAW_LIMIT)"))
    isfinite(max_seconds) && 0 < max_seconds <= KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS ||
        throw(ArgumentError(
            "max_seconds must lie in (0,$(KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS)]"))
    n, max_raw, max_seconds
end

"""Run one bounded stream and derive explicitly labelled extrapolations.

`classical_canonical_set` is intentionally an input, never constructed here;
it is expected to be the direct Carter-genus-0 canonical set when planarity is
needed by a caller.  The returned estimates are not exact unless the stream
reaches `raw_total`.
"""
function benchmark_knot_n6_all_state_enumeration(n::Integer=6;
                                                 max_raw::Integer=KNOT_N6_BENCHMARK_RAW_LIMIT,
                                                 max_seconds::Real=KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS,
                                                 classical_canonical_set=nothing)
    n = Int(n)
    max_raw = Int(max_raw)
    max_seconds = Float64(max_seconds)
    0 <= n <= 6 || throw(ArgumentError("n must lie in 0:6"))
    0 < max_raw <= KNOT_N6_BENCHMARK_RAW_LIMIT ||
        throw(ArgumentError("max_raw lies above the benchmark cap"))
    isfinite(max_seconds) && 0 < max_seconds <= KNOT_N6_BENCHMARK_TIME_LIMIT_SECONDS ||
        throw(ArgumentError("max_seconds lies above the benchmark cap"))

    stats = stream_canonical_knot_gauss(n;
        max_raw=max_raw,
        max_seconds=max_seconds,
        classical_canonical_set=classical_canonical_set)
    complete = stats.raw_seen == stats.raw_total
    raw_fraction = stats.raw_seen == 0 ? NaN : stats.raw_seen / stats.raw_total
    raw_per_second = stats.elapsed_seconds <= 0 ? NaN :
        stats.raw_seen / stats.elapsed_seconds
    emit_rate = stats.raw_seen == 0 ? NaN :
        stats.canonical_emitted / stats.raw_seen
    estimated_full_seconds = complete ? stats.elapsed_seconds :
        (stats.raw_seen == 0 ? NaN : stats.elapsed_seconds / raw_fraction)
    estimated_canonical_count = complete ? Float64(stats.canonical_emitted) :
        (stats.raw_seen == 0 ? NaN : stats.canonical_emitted / raw_fraction)
    average_serialization_bytes = stats.canonical_emitted == 0 ? NaN :
        stats.emitted_serialization_bytes / stats.canonical_emitted
    estimated_canonical_payload_bytes =
        isfinite(estimated_canonical_count) && isfinite(average_serialization_bytes) ?
            round(Int, estimated_canonical_count * average_serialization_bytes) : 0
    raw_scalar_lower_bound_bytes = big(stats.raw_total) * (3 * n * sizeof(Int))
    (stats=stats,
     complete=complete,
     raw_fraction=raw_fraction,
     raw_per_second=raw_per_second,
     canonical_emit_rate=emit_rate,
     estimated_full_seconds=estimated_full_seconds,
     estimated_canonical_count=estimated_canonical_count,
     average_serialization_bytes=average_serialization_bytes,
     estimated_canonical_payload_bytes=estimated_canonical_payload_bytes,
     raw_scalar_lower_bound_bytes=raw_scalar_lower_bound_bytes)
end

function _benchmark_print(result, n, max_raw, max_seconds)
    stats = result.stats
    println("knot_n6_all_state_enumeration_version=$(KNOT_N6_ALL_STATE_ENUMERATION_VERSION)")
    println("julia_version=$(VERSION)")
    println("scope=n=$(n);rooted_normalized_signed;orientations=true;canonical_identity_only")
    println("limits=max_raw=$(max_raw);max_seconds=$(max_seconds);full_n6_scan=false")
    println("raw_total=$(stats.raw_total)")
    println("raw_seen=$(stats.raw_seen)")
    println("raw_per_second=$( @sprintf("%.9f", result.raw_per_second) )")
    println("canonical_emitted=$(stats.canonical_emitted)")
    println("canonical_emit_rate=$( @sprintf("%.12g", result.canonical_emit_rate) )")
    println("planar_membership=external_direct_classical_set;command_line_set_supplied=false")
    println("direct_classical_n6_reference_canonical_count=71628;membership_not_constructed_in_bounded_run=true")
    println("planar_emitted=$(stats.planar_emitted)")
    println("nonplanar_emitted=$(stats.nonplanar_emitted)")
    println("planar_unknown=$(stats.planar_unknown)")
    println("complete=$(result.complete)")
    println("stop_reason=$(stats.stop_reason)")
    println("elapsed_seconds=$( @sprintf("%.9f", stats.elapsed_seconds) )")
    println("estimated_full_scan_seconds=$( @sprintf("%.9f", result.estimated_full_seconds) );estimate=$( !result.complete )")
    println("estimated_full_scan_minutes=$( @sprintf("%.9f", result.estimated_full_seconds / 60) );estimate=$( !result.complete )")
    println("estimated_canonical_count=$( @sprintf("%.3f", result.estimated_canonical_count) );estimate=$( !result.complete )")
    if !result.complete && stats.canonical_emitted == 0
        println("estimate_caveat=raw_prefix_has_zero_canonical_emissions;canonical_count_estimate_is_uninformative")
    elseif !result.complete
        println("estimate_caveat=raw_prefix_rate_extrapolation_only;not_an_exact_N6_count")
    end
    println("average_emitted_serialization_bytes=$( @sprintf("%.3f", result.average_serialization_bytes) )")
    println("estimated_canonical_payload_bytes=$(result.estimated_canonical_payload_bytes);estimate=$( !result.complete );utf8_payload_lower_bound=true")
    println("stream_retained_collection_bytes=0;callback_storage_excluded=true")
    println("stream_transient_scalar_lower_bound_bytes=$(3 * n * sizeof(Int));per_diagram_only=true")
    println("all_raw_scalar_lower_bound_bytes=$(result.raw_scalar_lower_bound_bytes);raw_tuple_scalars_only=true")
    println("all_canonical_set_allocator_overhead=not_estimated")
    println("emitted_serialization_bytes=$(stats.emitted_serialization_bytes)")
    println("emitted_fingerprint_sha256=$(stats.emitted_fingerprint)")
    nothing
end

function main(args=ARGS)
    n, max_raw, max_seconds = _benchmark_parse_args(args)
    result = benchmark_knot_n6_all_state_enumeration(n;
        max_raw=max_raw, max_seconds=max_seconds)
    _benchmark_print(result, n, max_raw, max_seconds)
    nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
