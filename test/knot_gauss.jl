using Test
using Random

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))

@testset "KnotGaussDiagram signed Gauss representation validation" begin
    @test collect(KnotGaussDiagram([1, -1], [1]).gauss) == [1, -1]
    @test_throws ArgumentError KnotGaussDiagram([1, -1, 2], [1])
    @test_throws ArgumentError KnotGaussDiagram([1, 1], [1])
    @test_throws ArgumentError KnotGaussDiagram([1, -1], [0])
    @test_throws ArgumentError KnotGaussDiagram([1, -1, 3, -3], [1])
    @test_throws DimensionMismatch reindex(KnotGaussDiagram([1, -1], [1]), [1, 2])

    # The pairing checks are representation checks, not a classical-planarity
    # certificate.  This known non-realizable Gauss word is intentionally
    # accepted here and must be filtered by a separately justified predicate.
    nonrealizable_word = KnotGaussDiagram([1, 2, 3, -1, -3, -2], [1, 1, 1])
    @test crossing_count(nonrealizable_word) == 3
    @test KNOT_GAUSS_CENSUS_KIND == :signed_gauss_representation
    @test !classical_knot_census_available()
    @test_throws ArgumentError require_classical_realizability(nonrealizable_word)

    d = KnotGaussDiagram([1, -2, 2, -1], [1, -1])
    @test crossing_count(d) == 2
    @test collect(cyclic_shift(d, 1).gauss) == [-1, 1, -2, 2]
    @test collect(cyclic_shift(d, -1).gauss) == [-2, 2, -1, 1]
    relabeled = reindex(d, [2, 1])
    @test canonicalize(relabeled) == canonicalize(d)
    @test reindex(d, Dict(1 => 2, 2 => 1)) == relabeled
    @test reindex(d, i -> 3 - i) == relabeled
    original_hash = canonical_hash(d)
    @test_throws MethodError setindex!(d.gauss, -1, 1)
    @test_throws MethodError setindex!(d.orientation, -1, 1)
    @test canonical_hash(d) == original_hash
end

@testset "Knot Gauss canonicalization" begin
    d = KnotGaussDiagram([1, -2, 2, -1], [1, -1])
    @test canonicalize_slow(d) == canonicalize(d)
    @test canonicalize_fast(d) == canonicalize(d)
    @test canonical_hash(d) == canonical_hash(cyclic_shift(d, 3))
    @test canonical_hash(d) == canonical_hash(reindex(d, [2, 1]))
    @test canonical_serialization(d) == canonical_serialize(d)
    @test occursin("canon=v1", canonical_serialization(d))
    @test_throws ArgumentError canonical_serialization(d; canon_version=2)
    @test rotation_automorphism_count(KnotGaussDiagram([1, -1], [1])) == 1

    # These transformations are intentionally outside the equivalence relation.
    asymmetric = KnotGaussDiagram([1, 2, -1, -2], [1, -1])
    @test canonical_hash(asymmetric) != canonical_hash(mirror(asymmetric))
    @test canonical_hash(asymmetric) != canonical_hash(orientation_reverse(asymmetric))
    @test orientation_reverse(asymmetric) == KnotGaussDiagram([-2, -1, 2, 1], [1, -1])
end

@testset "Deterministic rooted enumeration" begin
    expected = [1, 2, 12, 120, 1680, 30240, 665280]
    @test [rooted_normalized_signed_gc_count(n) for n in 0:6] == expected
    @test [length(enumerate_knot_gauss(n)) for n in 0:6] == expected .* [2^n for n in 0:6]
    # n=6 is checked through the lazy iterator interface; do not collect it.
    stream6 = enumerate_knot_gauss(6)
    @test Base.IteratorSize(typeof(stream6)) == Base.HasLength()
    first6, state6 = iterate(stream6)
    @test first6 isa KnotGaussDiagram
    @test state6 == 1
    @test collect(enumerate_rooted_normalized_signed_gc(2)) == [
        [1, -1, 2, -2], [1, -1, -2, 2], [1, 2, -1, -2],
        [1, 2, -2, -1], [1, -2, -1, 2], [1, -2, 2, -1],
        [-1, 1, 2, -2], [-1, 1, -2, 2], [-1, 2, 1, -2],
        [-1, 2, -2, 1], [-1, -2, 1, 2], [-1, -2, 2, 1]]
    @test collect(first(enumerate_knot_gauss(3)).orientation) == [1, 1, 1]
    @test collect(last(collect(enumerate_knot_gauss(2))).orientation) == [-1, -1]
end

@testset "Slow/fast oracle agreement" begin
    for n in 0:3, d in enumerate_knot_gauss(n)
        @test canonicalize_slow(d) == canonicalize(d)
    end
    rng = MersenneTwister(0x210821)
    for n in 1:6
        for _ in 1:25
            code_index = rand(rng, 0:rooted_normalized_signed_gc_count(n)-1)
            code = _unrank_rooted_signed(n, code_index)
            orientation = rand(rng, (-1, 1), n)
            d = KnotGaussDiagram(code, orientation)
            @test canonicalize_slow(d) == canonicalize(d)
        end
    end
end
