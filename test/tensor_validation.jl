using Test, Topology

# Real expressions lifted verbatim from examples/C250406TensorNetwork.jl and
# the MATLAB source C250815TensorNetwork.m. These must validate clean.
const VALID_NO_DICT = [
    ("C{3,4,1}M{4,2,5}", [1,2,5,3]),                       # calcTensorExpression worked example
    ("Il{1}Cr{2}(An^-4){1,2}", Int[]),                     # scalar output, ord = []
    ("Tp{0,1,2,3}(An^2){3,0}", [2,1]),                     # index 0 is a normal contraction label
    ("Mu{1,2,3}Cl{3}(An^-2){1,2}", Int[]),
    ("Ep{1}Et{2}", [1,2]),
]

@testset "validator accepts real expressions (no dict, arity off)" begin
    for (expr,ord) in VALID_NO_DICT
        @test isempty(validate_tensor_expression(expr, ord))
    end
end

@testset "validator accepts a full C250406 expression with a tensor dict" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    expr="Mu{1,2,3}Il{3}An{4,2}De{7,5,6}De{6,4,8}(An^-2){8,9}" *
         "Mu{9,12,10}Mu{10,5,11}De{13,14,12}Cl{13}(An){14,15}"
    @test isempty(validate_tensor_expression(expr, [11,1,7,15]; tensors=tensors))
end

@testset "network mode retains a multiply-occurring output label" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    # Labels 1:4 occur twice in the expression and are retained as outputs.
    # This is legal in the Julia factor-network extension, but not under the
    # MATLAB expr+ord exactly-twice rule.
    expr="Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}Tp{16,9,10,15}" *
         "(An^4){7,6}(An^-2){9,8}(An^-4){11,10}(An^2){1,3}" *
         "(An^-2){14,2}(An^2){15,16}(An^-2){13,12}(An^2){4,14}"
    @test !isempty(validate_tensor_expression(expr,1:4;tensors,index_rule=:matlab))
    @test isempty(validate_tensor_expression(expr,1:4;tensors,index_rule=:network))
end

@testset "L2 multiset: every index appears exactly twice" begin
    # index 3 dangles: once in pattern, absent from ord
    v=validate_tensor_expression("Mu{1,2,3}", [1,2])
    @test !isempty(v)
    @test any(occursin("3", m) && occursin("twice", m) for m in v)

    # index 2 over-contracted: three occurrences
    v=validate_tensor_expression("An{1,2}An{2,3}An{2,4}", [1,3,4])
    @test any(occursin("2", m) && occursin("appears 3", m) for m in v)

    # ord index absent from the pattern
    v=validate_tensor_expression("Mu{1,2,3}", [1,2,9])
    @test !isempty(v)
end

@testset "ord must not repeat a label" begin
    v=validate_tensor_expression("Mu{1,2,3}", [3,3])
    @test any(occursin("ord", m) && occursin("duplicate", m) for m in v)
end

@testset "lexical layer rejects residue and malformed legs" begin
    @test !isempty(validate_tensor_expression("Mu{1,2,3} garbage An{4,5}", Int[]))
    @test !isempty(validate_tensor_expression("Mu{1,-2,3}", Int[]))   # negatives are powers, not indices
    @test !isempty(validate_tensor_expression("Mu{1,,3}", Int[]))     # empty leg
    @test !isempty(validate_tensor_expression("nonsense", Int[]))     # no factors at all
end

@testset "power syntax requires exactly two legs" begin
    v=validate_tensor_expression("(Mu^2){1,2,3}", [1,2,3])
    @test any(occursin("power", m) for m in v)
end

@testset "arity and unknown-label checks need a tensor dict" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    # Mu has rank 3; two legs is an arity error
    v=validate_tensor_expression("Mu{1,2}", [1,2]; tensors=tensors)
    @test any(occursin("Mu", m) && occursin("rank", m) for m in v)
    # unknown label
    v=validate_tensor_expression("Foo{1,2}", [1,2]; tensors=tensors)
    @test any(occursin("Foo", m) for m in v)
    @test any(occursin("column 1",m) for m in v)
end

@testset "diagnostics identify the source position" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    v=validate_tensor_expression("Mu{1,2,3} junk",[1,2,3];tensors,index_rule=:network)
    @test any(occursin("column 11",m) for m in v)
    v=validate_tensor_expression("Et{1}Mu{2,3}",[1,2,3];tensors,index_rule=:network)
    @test any(occursin("column 6",m) && occursin("rank",m) for m in v)
end

@testset "evaluator rejects an invalid expression before contracting" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    # residue between factors must be rejected, not skipped
    @test_throws ArgumentError evaluate_tensor_expression("Mu{1,2,3} junk", tensors, [1,2,3]; dimension=4)
    # Strict MATLAB compatibility remains available explicitly.
    @test_throws ArgumentError evaluate_tensor_expression(
        "Mu{1,2,3}",tensors,[1,2];dimension=4,index_rule=:matlab)
end

@testset "evaluator still accepts every real C250406 expression" begin
    tensors=legacy_hopf_tensors(sweedler_algebra())
    @test_nowarn evaluate_tensor_expression("Mu{1,2,3}Cl{3}(An^-2){1,2}", tensors, Int[]; dimension=4)
    @test_nowarn evaluate_tensor_expression("Ep{1}Et{2}", tensors, [1,2]; dimension=4)
end

@testset "full framing path remains valid in network mode" begin
    include(joinpath(@__DIR__,"..","examples","C250406TensorNetwork.jl"))
    @test_nowarn run_C250406_sweedler_framing()
end
