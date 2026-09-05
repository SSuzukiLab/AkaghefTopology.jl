using Test
using SHA

# `_all_canonical_states` is the existing repository census helper.  It is
# included only as an independent test oracle; this test does not write the
# builder's DB outputs.
include(joinpath(@__DIR__, "..", "scripts", "build_knot_n5_classification_db.jl"))
include(joinpath(@__DIR__, "..", "src", "knot_n6_all_state_enumeration.jl"))

@testset "streamed canonical all-state enumeration" begin
    states, _, expected_counts, raw_total = _all_canonical_states(5)
    @test expected_counts == collect(KNOT_GAUSS_EXPECTED_CANONICAL_COUNTS)
    @test raw_total == sum(length(enumerate_knot_gauss(n; orientations=true)) for n in 0:5)

    expected = [Set{String}() for _ in 0:5]
    for diagram in states
        push!(expected[crossing_count(diagram) + 1], canonical_serialization(diagram))
    end

    observed = [Set{String}() for _ in 0:5]
    observed_order = [String[] for _ in 0:5]
    streamed_stats = Vector{Any}(undef, 6)
    for n in 0:5
        # Keep the callback deliberately one-argument here so the complete
        # record (including the planarity state) is checked by the test.
        streamed_stats[n + 1] = stream_canonical_knot_gauss(n) do record
            @test record.planar === missing
            push!(observed[n + 1], record.serialization)
            push!(observed_order[n + 1], record.serialization)
            @test record.serialization == canonical_serialization(record.diagram)
            @test record.diagram == canonicalize(record.diagram)
        end
    end

    @test [stats.canonical_emitted for stats in streamed_stats] ==
        collect(KNOT_GAUSS_EXPECTED_CANONICAL_COUNTS)
    @test [length(set) for set in observed] ==
        collect(KNOT_GAUSS_EXPECTED_CANONICAL_COUNTS)
    @test observed == expected
    @test all(stats -> stats.raw_seen == stats.raw_total, streamed_stats)
    @test all(stats -> !stats.stopped && stats.stop_reason === :complete,
              streamed_stats)
    @test all(stats -> stats.planar_unknown == stats.canonical_emitted &&
                      stats.planar_emitted == 0 && stats.nonplanar_emitted == 0,
              streamed_stats)

    # This word has a nontrivial rotational automorphism.  The rooted input
    # is not itself canonical, but its canonical representative is emitted
    # exactly once; no periodic-orbit multiplicity is hidden by the callback.
    periodic = KnotGaussDiagram([1, -1, 2, -2], [1, 1])
    periodic_serialization = canonical_serialization(periodic)
    @test rotation_automorphism_count(periodic) == 2
    @test count(==(periodic_serialization), observed_order[3]) == 1

    # Independently recompute the aggregate stream digest from callback data.
    # The digest is over emission order plus a newline, not the per-state
    # `canonical_hash` convention.
    for n in 0:5
        context = SHA.SHA256_CTX()
        for serialization in observed_order[n + 1]
            SHA.update!(context, codeunits(serialization))
            SHA.update!(context, UInt8[0x0a])
        end
        @test bytes2hex(SHA.digest!(context)) == streamed_stats[n + 1].emitted_fingerprint
    end

    # A supplied direct-classical canonical set is queried by serialization;
    # this fixture checks the plumbing only (it is not a new realizability
    # claim).  The other canonical N=1 state must be reported nonplanar under
    # this intentionally partial injected membership relation.
    one_crossing = sort!(collect(expected[2]))
    injected_planar = Set([one_crossing[1]])
    two_arg_flags = Bool[]
    injected_stats = stream_canonical_knot_gauss(1;
        classical_canonical_set=injected_planar) do diagram, planar
        @test diagram isa KnotGaussDiagram
        @test planar isa Bool
        push!(two_arg_flags, planar)
    end
    @test injected_stats.canonical_emitted == 2
    @test injected_stats.planar_emitted == 1
    @test injected_stats.nonplanar_emitted == 1
    @test injected_stats.planar_unknown == 0
    @test sort(two_arg_flags) == [false, true]

    predicate_stats = stream_canonical_knot_gauss(1;
        planar_membership=serialization -> serialization in injected_planar)
    @test predicate_stats.planar_emitted == injected_stats.planar_emitted
    @test predicate_stats.nonplanar_emitted == injected_stats.nonplanar_emitted
    @test predicate_stats.emitted_fingerprint == injected_stats.emitted_fingerprint

    bounded = stream_canonical_knot_gauss(3; max_raw=17)
    @test bounded.raw_total == 960
    @test bounded.raw_seen == 17
    @test bounded.stopped && bounded.stop_reason === :max_raw
    @test bounded.canonical_emitted <= bounded.raw_seen
end
