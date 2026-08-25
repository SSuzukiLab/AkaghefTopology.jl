using Test, Topology

const S3_TABLE = [
    1 2 3 4 5 6
    2 3 1 6 4 5
    3 1 2 5 6 4
    4 5 6 1 2 3
    5 6 4 3 1 2
    6 4 5 2 3 1
]

@testset "group algebra from a multiplication table" begin
    algebra=group_algebra(S3_TABLE)
    @test algebra.dimension==6
    @test algebra.multiplication.terms[(2,4,6)]==1
    @test algebra.multiplication.terms[(4,2,5)]==1
    @test !haskey(algebra.multiplication.terms,(2,4,5))
    @test algebra.unit.terms==Dict{Tuple,ComplexF64}((1,)=>1)
    @test algebra.coproduct.terms==Dict{Tuple,ComplexF64}((i,i,i)=>1 for i in 1:6)
    @test algebra.counit.terms==Dict{Tuple,ComplexF64}((i,)=>1 for i in 1:6)
    @test algebra.antipode.terms==Dict{Tuple,ComplexF64}(
        (1,1)=>1,(3,2)=>1,(2,3)=>1,(4,4)=>1,(5,5)=>1,(6,6)=>1)
    @test all(values(verify_hopf(algebra)))
    @test algebra.right_integral.terms==Dict{Tuple,ComplexF64}((1,)=>1)
    @test algebra.left_integral.terms==algebra.right_integral.terms
    @test algebra.right_cointegral.terms==Dict{Tuple,ComplexF64}((i,)=>1 for i in 1:6)
    @test algebra.left_cointegral.terms==algebra.right_cointegral.terms
    @test sum(algebra.right_integral.terms[(i,)]*algebra.right_cointegral.terms[(i,)]
              for i in 1:6 if haskey(algebra.right_integral.terms,(i,)))==1
end

@testset "group algebra from ordered labels and an operation" begin
    canonical=[:e,:r,:r2,:s,:sr,:sr2]
    elements=[:s,:e,:r,:r2,:sr,:sr2]
    multiply(left,right)=canonical[S3_TABLE[findfirst(==(left),canonical),
                                                  findfirst(==(right),canonical)]]
    algebra=group_algebra(elements,multiply)
    @test algebra.dimension==6
    @test algebra.unit.terms==Dict{Tuple,ComplexF64}((2,)=>1)
    @test all(values(verify_hopf(algebra)))
end

@testset "group algebra rejects nongroups" begin
    @test_throws ArgumentError group_algebra(zeros(Int,0,0))
    @test_throws ArgumentError group_algebra([1 2 1; 2 1 2])
    @test_throws ArgumentError group_algebra([1 2; 2 3])
    @test group_algebra([2 1; 1 2]).unit.terms==Dict{Tuple,ComplexF64}((2,)=>1)
    nonassociative=copy(S3_TABLE); nonassociative[2,2]=2
    @test_throws ArgumentError group_algebra(nonassociative)
    @test_throws ArgumentError group_algebra([1 2; 2 2])
    @test_throws ArgumentError group_algebra([:e,:e],(left,right)->:e)
    @test_throws ArgumentError group_algebra([:e,:a],(left,right)->:outside)
end

@testset "cyclic group algebra remains compatible" begin
    legacy=cyclic_group_algebra(5)
    general=group_algebra(collect(0:4),(left,right)->mod(left+right,5))
    @test legacy.dimension==general.dimension
    @test all(field -> getfield(legacy,field).terms==getfield(general,field).terms,
              filter(!=(:dimension),fieldnames(FiniteHopfAlgebra)))
    @test (basis_element(legacy,2)^5).coefficients≈basis_element(legacy,1).coefficients
    @test all(values(verify_hopf(legacy)))
end
