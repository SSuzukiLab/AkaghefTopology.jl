using Test, Topology

@testset "SparseFactor preserves exact rational coefficients" begin
    left=SparseFactor(Int[],Dict{Tuple,Rational{Int}}(
        (1,1)=>1//2,
        (1,2)=>1//3,
        (2,1)=>2//5,
    ))
    right=SparseFactor(Int[],Dict{Tuple,Rational{Int}}(
        (1,1)=>3//7,
        (2,1)=>5//11,
    ))
    result=evaluate_tensor_expression(
        "A{1,2}B{2,3}",Dict("A"=>left,"B"=>right),[1,3])

    @test result isa SparseFactor{Rational{Int}}
    @test result.terms==Dict{Tuple,Rational{Int}}(
        (1,1)=>(1//2)*(3//7)+(1//3)*(5//11),
        (2,1)=>(2//5)*(3//7),
    )
    @test eltype(dense_tensor(result,2))==Rational{Int}
end

@testset "mixed coefficient types promote without floating conversion" begin
    integer_factor=SparseFactor(Int[],Dict{Tuple,Int}((1,1)=>2,(2,2)=>3))
    rational_factor=SparseFactor(Int[],Dict{Tuple,Rational{Int}}(
        (1,1)=>1//3,(2,1)=>2//5))
    result=evaluate_tensor_expression(
        "A{1,2}B{2,3}",Dict("A"=>integer_factor,"B"=>rational_factor),[1,3])
    @test result isa SparseFactor{Rational{Int}}
    @test result.terms==Dict{Tuple,Rational{Int}}((1,1)=>2//3,(2,1)=>6//5)
end

@testset "repeated labels and matrix powers preserve coefficient type" begin
    diagonal_source=SparseFactor(Int[],Dict{Tuple,Rational{Int}}(
        (1,1)=>1//2,(1,2)=>9//10,(2,2)=>2//3))
    diagonal=evaluate_tensor_expression(
        "A{1,1}",Dict("A"=>diagonal_source),[1];index_rule=:network)
    @test diagonal isa SparseFactor{Rational{Int}}
    @test diagonal.terms==Dict{Tuple,Rational{Int}}((1,)=>1//2,(2,)=>2//3)

    matrix=SparseFactor(Int[],Dict{Tuple,Rational{Int}}(
        (1,1)=>1//2,(2,2)=>2//3))
    squared=evaluate_tensor_expression(
        "(A^2){1,2}",Dict("A"=>matrix),[1,2];dimension=2,index_rule=:network)
    @test squared isa SparseFactor{Rational{Int}}
    @test squared.terms==Dict{Tuple,Rational{Int}}((1,1)=>1//4,(2,2)=>4//9)
end
