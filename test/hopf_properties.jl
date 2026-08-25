using Test, Topology

const S3_HOPF_TABLE = [
    1 2 3 4 5 6
    2 3 1 6 4 5
    3 1 2 5 6 4
    4 5 6 1 2 3
    5 6 4 3 1 2
    6 4 5 2 3 1
]

@testset "finite Hopf algebra predicates" begin
    cases=[
        (name="C[S3]",algebra=group_algebra(S3_HOPF_TABLE),involutory=true,unimodular=true),
        (name="Sweedler",algebra=sweedler_algebra(),involutory=false,unimodular=false),
        (name="exterior(2)",algebra=exterior_algebra(2),involutory=true,unimodular=true),
        (name="Kac-Paljutkin",algebra=kac_paljutkin_algebra(),involutory=true,unimodular=true),
        (name="uqsl2(1/2)",algebra=uqsl2_borel_small(1//2),involutory=false,unimodular=false),
        (name="uqsl2(2/5)",algebra=uqsl2_borel_small(2//5),involutory=false,unimodular=false),
    ]
    for case in cases
        @testset "$(case.name)" begin
            @test is_involutory(case.algebra)==case.involutory
            @test is_unimodular(case.algebra)==case.unimodular
        end
    end
end

@testset "unimodularity ignores integral normalization" begin
    algebra=group_algebra(S3_HOPF_TABLE)
    algebra.left_integral.terms[(1,)]=2
    @test is_unimodular(algebra)
    empty!(algebra.left_integral.terms)
    @test !is_unimodular(algebra)
end
