using Test, Topology

@testset "Dedekind sums and lens spaces" begin
    @test dedekind_sum(1, 1) == 0
    @test dedekind_sum(1, 3) == 1//18
    @test dedekind_sum(1, 4) == 1//8
    @test dedekind_sum(2, 7) == 1//14
    @test dedekind_sum(2, 3) == -1//18
    @test casson_walker_lens(3, 1) == -1//36
    @test casson_walker_lens(4, 1) == -1//16
    @test casson_walker_lens(7, 2) == -1//28
    @test casson_walker_lens(2, 1) == 0
    @test casson_walker_lens(3, 2) == 1//36
end

@testset "two surgery writings agree" begin
    for p in 1:25, q in 1:25
        gcd(p, q) == 1 || continue
        @test casson_walker_lens(p, q) == casson_walker_knot_surgery(p, q, 0)
        @test casson_walker_knot_surgery(p, q, 0) == casson_walker_surgery_pup(p, q, 0)
        @test casson_walker_knot_surgery(-p, q, 0) == -casson_walker_lens(p, q)
    end
    @test casson_walker_knot_surgery(1, 1, 2) == 1
    @test casson_walker_knot_surgery(-1, 1, 2) == -1
    @test casson_walker_knot_surgery(7, 1, -2) == -9//28
    @test casson_walker_knot_surgery(-7, 1, -2) == 9//28
    @test casson_walker_knot_surgery(8, 3, 0) == -1//32
    @test casson_walker_knot_surgery(-8, 3, 0) == 1//32
    @test casson_walker_surgery_pup(1, 1, 2) == 1
    @test casson_walker_surgery_pup(-7, 1, -2) == 9//28
end

@testset "Brieskorn spheres and the Seifert closed form" begin
    @test casson_walker_brieskorn(2, 3, 5) == -1
    @test casson_walker_brieskorn(2, 3, 7) == -1
    @test casson_walker_brieskorn(2, 3, 11) == -2
    @test casson_walker_brieskorn(2, 5, 7) == -2
    @test casson_walker_brieskorn(2, 3, 7) == -(2^2 - 1) * (3^2 - 1)//24
    @test casson_walker_brieskorn(2, 5, 11) == -(2^2 - 1) * (5^2 - 1)//24
    @test casson_walker_brieskorn(3, 4, 13) == -(3^2 - 1) * (4^2 - 1)//24
    @test casson_walker_seifert(-5, Tuple{Int,Int}[]) == (lambda = -1//10, h1 = 5, E = 5)
    @test casson_walker_seifert(-1, [(2, 1), (3, 1), (5, 1)]) ==
        (lambda = -1, h1 = 1, E = -1//30)
    @test casson_walker_seifert(1, [(2, -1), (3, -1), (7, -1)]) ==
        (lambda = -1, h1 = 1, E = -1//42)
    @test casson_walker_seifert(0, [(2, 1), (3, 1), (5, 1)]).lambda == -29//124
    pos = casson_walker_seifert(0, [(2, 1), (3, 1), (5, 1)])
    neg = casson_walker_seifert(0, [(2, -1), (3, -1), (5, -1)])
    @test pos.lambda + neg.lambda == 0
end
