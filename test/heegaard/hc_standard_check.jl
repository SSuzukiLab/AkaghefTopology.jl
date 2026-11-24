using Test
using Topology

const H = Topology.Heegaard

function standard_lc(words::Dict{String,Vector{Int}})
    H.LickorishCode(1, words,
        Dict{String,String}(), Dict{String,Vector{String}}(),
        Dict{String,String}(), Dict{String,Int}(), nothing)
end

@testset "same signed word is standard" begin
    bound = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1], "m1" => [1, 2], "b" => [-2])),
        Dict("m1" => :alpha, "b" => :beta))
    judged = H.hc_standard_check(bound)
    @test judged.standard
    @test judged.diagram === bound
    @test judged.correspondence == ["m1" => "m1"]
    @test isempty(judged.failures)
    @test H.is_standard(judged.diagram)

    parallel = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 2], "b" => [-2, -4], "a" => [3, 4])),
        Dict("a" => :alpha, "b" => :beta))
    moved = H.hc_standard_check(parallel)
    @test moved.standard
    @test moved.correspondence == ["a" => "m1"]
    @test moved.diagram.lc.words["l1"] == [-1]
    @test moved.diagram.lc.words["m1"] == [1, 2]
    @test moved.diagram.lc.words["b"] == [-2]
    @test !haskey(moved.diagram.lc.words, "a")
    @test moved.diagram.roles == Dict("m1" => :alpha, "b" => :beta)
    @test H.is_standard(moved.diagram)

    rotated = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 2], "b" => [-2, -4], "a" => [4, 3])),
        Dict("a" => :alpha, "b" => :beta))
    turned = H.hc_standard_check(rotated)
    @test turned.standard
    @test turned.correspondence == ["a" => "m1"]
    @test turned.diagram.lc.words["m1"] == [1, 2]
    @test turned.diagram.lc.words["b"] == [-2]
end

@testset "opposite orientation is the same meridian" begin
    # a = [-4, -3] has class [-1, 0] = -[m1]. Reversing a negates the reversed
    # word, so a becomes [3, 4]. Labels 3 and 4 flip on l1 and b:
    # l1 = [-1, -3], b = [-2, -4]. That is the parallel copy above.
    reversed = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, 3], "m1" => [1, 2], "b" => [-2, 4], "a" => [-4, -3])),
        Dict("a" => :alpha, "b" => :beta))
    judged = H.hc_standard_check(reversed)
    @test judged.standard
    @test judged.correspondence == ["a" => "m1"]
    @test judged.diagram.lc.words["l1"] == [-1]
    @test judged.diagram.lc.words["m1"] == [1, 2]
    @test judged.diagram.lc.words["b"] == [-2]
    @test !haskey(judged.diagram.lc.words, "a")
    @test judged.diagram.roles == Dict("m1" => :alpha, "b" => :beta)
    @test H.is_standard(judged.diagram)
    @test isempty(judged.witness)
    @test isempty(judged.failures)
end

@testset "a different signed word is not standard" begin
    # Label 3: i(a, l1) = -1. Label 5: i(a, m1) = -1. Class [1, -1], meridian [1, 0].
    wound = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 2, -5], "b" => [-2, -4], "a" => [3, 5, 4])),
        Dict("a" => :alpha, "b" => :beta))
    held = H.hc_standard_check(wound)
    @test !held.standard
    @test held.diagram === nothing
    @test held.failures[1].curve == "a"
    @test held.failures[1].reason == :homology_class
    @test held.failures[1].detail["alpha_class"] == [1, -1]
    @test held.failures[1].detail["meridian"] == "m1"
    @test held.failures[1].detail["meridian_class"] == [1, 0]
    @test isempty(held.witness)

    occupied = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 2], "b" => [-2, -4], "a" => [3, 4])),
        Dict("a" => :alpha, "m1" => :beta, "b" => :beta))
    blocked = H.hc_standard_check(occupied)
    @test !blocked.standard
    @test blocked.failures[1].reason == :meridian_occupied

    marked = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 2], "b" => [-2, -4], "a" => [3, 4])),
        Dict("a" => :alpha, "b" => :beta); points=Dict("p" => "a"))
    decorated = H.hc_standard_check(marked)
    @test !decorated.standard
    @test decorated.failures[1].reason == :decoration_blocks_removal
end

@testset "genus 2 pairs each parallel with its meridian" begin
    lc = H.LickorishCode(2, Dict(
        "l1" => [-1, -3, -5, -7],
        "l2" => [-2, 4, -6],
        "m1" => [1],
        "m2" => [2],
        "n1" => [3, -4],
        "a1" => [5],
        "a2" => [6],
        "b" => [7],
    ), Dict{String,String}(), Dict{String,Vector{String}}(),
        Dict{String,String}(), Dict{String,Int}(), nothing)
    hc = H.resolve_hc(lc, Dict("a1" => :alpha, "a2" => :alpha, "b" => :beta))
    judged = H.hc_standard_check(hc)
    @test judged.standard
    @test sort(judged.correspondence) == sort(["a1" => "m1", "a2" => "m2"])
    @test judged.diagram.lc.words["l1"] == [-1, -3, -7]
    @test judged.diagram.lc.words["l2"] == [-2, 4]
    @test judged.diagram.lc.words["m1"] == [1]
    @test judged.diagram.lc.words["m2"] == [2]
    @test judged.diagram.lc.words["n1"] == [3, -4]
    @test judged.diagram.lc.words["b"] == [7]
    @test judged.diagram.roles == Dict("m1" => :alpha, "m2" => :alpha, "b" => :beta)

    both = H.LickorishCode(2, Dict(
        "l1" => [-1, -3, -5, -7],
        "l2" => [-2, 4],
        "m1" => [1],
        "m2" => [2],
        "n1" => [3, -4],
        "a1" => [5],
        "a2" => [7],
    ), Dict{String,String}(), Dict{String,Vector{String}}(),
        Dict{String,String}(), Dict{String,Int}(), nothing)
    crowded = H.hc_standard_check(H.resolve_hc(both, Dict("a1" => :alpha, "a2" => :alpha)))
    @test !crowded.standard
    @test crowded.failures[1].reason == :no_bijection
end

@testset "homology refuses m+l and one bigon moves onto m1" begin
    # +3 on a and -3 on l1 give i(a,l1)=-1. -5 on a and +5 on m1 give i(a,m1)=+1.
    # Coefficients are x=-i(a,l1), y=i(a,m1), so the class is [1, 1].
    summed = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3], "m1" => [1, 5], "a" => [-5, 3])),
        Dict("a" => :alpha))
    off = H.hc_standard_check(summed)
    @test !off.standard
    @test off.diagram === nothing
    @test off.failures[1].curve == "a"
    @test off.failures[1].reason == :homology_class
    @test off.failures[1].detail["alpha_class"] == [1, 1]
    @test off.failures[1].detail["meridian"] == "m1"
    @test off.failures[1].detail["meridian_class"] == [1, 0]

    # a meets l1 at 3 with the same passage as m1 at 1, plus the canceling pair 4, 5.
    # Removing that empty bigon leaves a : [3], l1 : [-1, -3], then the same-word step drops 3.
    bigon = H.resolve_hc(standard_lc(Dict(
        "l1" => [-1, -3, 5, -4], "m1" => [1], "a" => [3, -5, 4])),
        Dict("a" => :alpha))
    moved = H.hc_standard_check(bigon)
    @test moved.standard
    @test moved.correspondence == ["a" => "m1"]
    @test moved.diagram.lc.words["l1"] == [-1]
    @test moved.diagram.lc.words["m1"] == [1]
    @test !haskey(moved.diagram.lc.words, "a")
    @test moved.diagram.roles == Dict("m1" => :alpha)
    @test H.is_standard(moved.diagram)
    @test moved.witness == ["lc.isotopy.digon.remove.00"]
    @test isempty(moved.failures)
end
