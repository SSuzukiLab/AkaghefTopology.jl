using Test

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_realizability.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_enumeration.jl"))

function _assert_rooted_unsigned_word(word, n)
    @test length(word) == 2n
    n == 0 && return
    @test word[1] == 1
    @test sort(unique(word)) == collect(1:n)
    @test all(label -> count(==(label), word) == 2, 1:n)
    first_seen = Int[]
    for label in word
        label in first_seen || push!(first_seen, label)
    end
    @test first_seen == collect(1:n)
end
@testset "rooted unsigned restricted-growth stream" begin
    expected_counts = [1, 1, 3, 15, 105]
    for n in 0:4
        words = enumerate_rooted_unsigned_gauss_words(n)
        streamed = Tuple{Vararg{Int}}[]
        foreach_rooted_unsigned_gauss_word(n) do word
            push!(streamed, Tuple(word))
        end
        @test length(words) == expected_counts[n + 1]
        @test length(streamed) == expected_counts[n + 1]
        @test streamed == Tuple.(words)
        @test length(unique(streamed)) == expected_counts[n + 1]
        for word in words
            _assert_rooted_unsigned_word(word, n)
        end
    end
end

@testset "direct classical generation and exhaustive equivalence through n=4" begin
    expected_spherical = [1, 2, 8, 42, 260]
    expected_raw = [1, 4, 32, 336, 4160]
    expected_canonical = [1, 2, 10, 60, 536]
    direct_results = [direct_classical_gauss_enumeration(n) for n in 0:4]

    for n in 0:4
        result = direct_results[n + 1]
        @test result.unsigned_words == [1, 1, 3, 15, 105][n + 1]
        @test result.rotations_checked == result.unsigned_words * 2^n
        @test result.spherical_immersions == expected_spherical[n + 1]
        @test result.generated_decorations == expected_raw[n + 1]
        @test result.raw_duplicates == 0
        @test length(result.raw_keys) == expected_raw[n + 1]

        canonical = Set{String}()
        for key in result.raw_keys
            code, orientation = classical_gauss_key_to_code(key, n)
            @test _cge_pack_signed(code, orientation) == key
            @test classical_gauss_realizable_fast(code, orientation)
            push!(canonical,
                  canonical_serialization(KnotGaussDiagram(code, orientation)))
        end
        @test length(canonical) == expected_canonical[n + 1]

        exhaustive = Set{UInt64}()
        for diagram in enumerate_knot_gauss(n; orientations=true)
            classical_gauss_realizable_fast(diagram) || continue
            push!(exhaustive, _cge_pack_signed(diagram.gauss, diagram.orientation))
        end
        @test result.raw_keys == exhaustive
        @test isempty(setdiff(result.raw_keys, exhaustive))
        @test isempty(setdiff(exhaustive, result.raw_keys))
    end
end

@testset "packed key diagrams preserve repository convention" begin
    result = direct_classical_gauss_enumeration(3)
    for key in result.raw_keys
        code, orientation = classical_gauss_key_to_code(key, 3)
        diagram = classical_gauss_key_to_diagram(key, 3)
        @test diagram isa KnotGaussDiagram
        @test diagram.gauss == Tuple(code)
        @test diagram.orientation == Tuple(orientation)
    end
end
