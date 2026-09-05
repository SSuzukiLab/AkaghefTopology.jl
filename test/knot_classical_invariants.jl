using Test

const INVARIANT_SRC = joinpath(@__DIR__, "..", "src", "knot_classical_invariants.jl")
include(INVARIANT_SRC)
include(joinpath(@__DIR__, "..", "src", "reidemeister.jl"))

terms(pairs...) = Dict(Int(exponent) => BigInt(coefficient) for
                        (exponent, coefficient) in pairs)

const KNOTINFO_PD = Dict(
    "0_1" => Int[],
    "3_1" => [[1, 5, 2, 4], [3, 1, 4, 6], [5, 3, 6, 2]],
    "4_1" => [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]],
    "5_1" => [[2, 8, 3, 7], [4, 10, 5, 9], [6, 2, 7, 1],
              [8, 4, 9, 3], [10, 6, 1, 5]],
    "5_2" => [[1, 5, 2, 4], [3, 9, 4, 8], [5, 1, 6, 10],
              [7, 3, 8, 2], [9, 7, 10, 6]],
)

const EXPECTED_ALEXANDER = Dict(
    "0_1" => terms(0 => 1),
    "3_1" => terms(0 => 1, 1 => -1, 2 => 1),
    "4_1" => terms(0 => 1, 1 => -3, 2 => 1),
    "5_1" => terms(0 => 1, 1 => -1, 2 => 1, 3 => -1, 4 => 1),
    "5_2" => terms(0 => 2, 1 => -3, 2 => 2),
)

const EXPECTED_JONES = Dict(
    "0_1" => terms(0 => 1),
    "3_1" => terms(1 => 1, 3 => 1, 4 => -1),
    "4_1" => terms(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1),
    "5_1" => terms(2 => 1, 4 => 1, 5 => -1, 6 => 1, 7 => -1),
    "5_2" => terms(1 => 1, 2 => -1, 3 => 2, 4 => -1, 5 => 1, 6 => -1),
)

@testset "classical PD bridge and independent local invariants" begin
    diagrams = Dict{String,KnotGaussDiagram}()
    for name in keys(KNOTINFO_PD)
        diagram = pd_to_gauss(KNOTINFO_PD[name])
        diagrams[name] = diagram
        pd = pd_from_gauss(diagram)
        @test verify_pd(pd)
        @test crossing_count(diagram) == (name == "0_1" ? 0 : parse(Int, name[1]))
        @test polynomial_terms(alexander_polynomial(diagram)) == EXPECTED_ALEXANDER[name]
        @test polynomial_terms(jones_polynomial(diagram)) == EXPECTED_JONES[name]
    end

    trefoil = diagrams["3_1"]
    @testset "cyclic shift and relabel" begin
        for amount in -7:7
            shifted = cyclic_shift(trefoil, amount)
            @test polynomial_terms(alexander_polynomial(shifted)) ==
                polynomial_terms(alexander_polynomial(trefoil))
            @test jones_polynomial(shifted) == jones_polynomial(trefoil)
        end
        for permutation in ([3, 1, 2], [2, 3, 1], [1, 3, 2])
            relabelled = reindex(trefoil, permutation)
            @test polynomial_terms(alexander_polynomial(relabelled)) ==
                polynomial_terms(alexander_polynomial(trefoil))
            @test jones_polynomial(relabelled) == jones_polynomial(trefoil)
        end
    end

    @testset "orientation reversal and mirror" begin
        reversed = orientation_reverse(trefoil)
        @test polynomial_terms(alexander_polynomial(reversed)) ==
            polynomial_terms(alexander_polynomial(trefoil))
        @test jones_polynomial(reversed) == jones_polynomial(trefoil)

        reflected = mirror(trefoil)
        @test polynomial_terms(alexander_polynomial(reflected)) ==
            polynomial_terms(alexander_polynomial(trefoil))
        @test jones_polynomial(reflected) == invert_knot_variable(jones_polynomial(trefoil))
        @test jones_polynomial(reflected) != jones_polynomial(trefoil)
    end

    @testset "Jones A-to-t normalization is strict" begin
        @test _jones_t_terms(Dict(0 => big(1), 4 => big(-1))) ==
            Dict(0 => big(1), -1 => big(-1))
        @test_throws ArgumentError _jones_t_terms(Dict(2 => big(1), 6 => big(-1)))
    end

    @testset "R1/R2/R3 invariance" begin
        empty = KnotGaussDiagram(Int[], Int[])
        one_crossing = apply_reidemeister(empty,
            first(movable(empty, "R1"; param=:insert)))
        @test alexander_polynomial(one_crossing) == alexander_polynomial(empty)
        @test jones_polynomial(one_crossing) == jones_polynomial(empty)

        r2_insert = first(movable(empty, "R2"; param=:insert))
        two_crossings = apply_reidemeister(empty, r2_insert)
        @test alexander_polynomial(two_crossings) == alexander_polynomial(empty)
        @test jones_polynomial(two_crossings) == jones_polynomial(empty)
        r2_delete = first(movable(two_crossings, "R2"; param=:delete))
        @test verify_reidemeister_transition(two_crossings,
            apply_reidemeister(two_crossings, r2_delete), r2_delete)

        r3 = KnotGaussDiagram([-1, -2, -3, 1, 2, 3], [1, -1, 1])
        r3_candidate = first(movable(r3, "R3"))
        r3_after = apply_reidemeister(r3, r3_candidate)
        @test verify_reidemeister_transition(r3, r3_after, r3_candidate)
        @test alexander_polynomial(r3_after) == alexander_polynomial(r3)
        @test jones_polynomial(r3_after) == jones_polynomial(r3)
    end
end

@testset "PD bridge rejects links and malformed rotations" begin
    link_rows = [[1, 5, 2, 4], [3, 1, 4, 6], [5, 3, 6, 2],
                 [7, 11, 8, 10], [9, 7, 10, 12], [11, 9, 12, 8]]
    @test_throws ArgumentError pd_to_gauss(link_rows)
    @test_throws ArgumentError pd_to_gauss([[1, 2, 3, 4]])
end
