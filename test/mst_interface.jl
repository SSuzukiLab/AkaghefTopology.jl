using Test, Topology
import AkaghefAlgebra
isdefined(@__MODULE__,:vl_example) || include(joinpath(@__DIR__,"../examples/VLExample.jl"))

@testset "MST staged Taft interface" begin
    v=vl_example("s3_abalone"); calculate_weight(v;assign=true)
    before=deepcopy(v.weights)
    H1=TaftAlgebra(9); H2=TaftAlgebra(9;M=2)
    raw=invariant_expression(v,H1)
    @test raw.body==invariant_expression(v,H2).body
    reduced=simplify_invariant(raw;strategy=:qa)
    @test length(raw.remaining)==4
    @test isempty(reduced.remaining)
    @test v.weights==before
    base=evaluate_invariant(reduced,9)
    z1=normalize_root(base,RootSpec(9,1)); z2=normalize_root(base,RootSpec(9,2))
    @test z1 isa Cyclotomic{9,1}
    @test z2 isa Cyclotomic{9,2}
    @test z1==AkaghefAlgebra.TaftExact.zeta(9)^(-64)
    @test z2==AkaghefAlgebra.TaftExact.zeta(9)^(-128)
    @test Z_MST(v,H2)==z2
    @test hopf_invariant(v,2//9)==z2
    @test Z_MST(v,TaftAlgebra(9;M=11)) isa Cyclotomic{9,2}
    @test Z_MST(v,TaftAlgebra(:N)) isa AkaghefAlgebra.TaftExpression
    @test_throws ArgumentError evaluate_invariant(reduced,:N)
    @test_throws ArgumentError TaftAlgebra(9;M=3)
    @test_throws ArgumentError TaftAlgebra(1)
    @test_throws ArgumentError TaftAlgebra(:N;M=2)
    @test_throws ArgumentError Z_MST(v,H1;strategy=:unknown)
    @test_throws ArgumentError Z_MST(v,TaftAlgebra(3);strategy=:none,maxterms=1)
    @test Z_MST(v,TaftAlgebra(3);strategy=:none)==Z_MST(v,TaftAlgebra(3);strategy=:qa)
    generic=sweedler_algebra()
    expr=invariant_expression(v,generic)
    count_before=length(expr.factors)
    @test evaluate_invariant(expr)==evaluate_invariant(expr)
    @test length(expr.factors)==count_before
    @test evaluate_invariant(simplify_invariant(expr;strategy=:qa))==Z_MST(v,generic)
    @test Z_MST(v,generic)==hopf_invariant(v,generic;method=:numeric)
end

@testset "four registered MST values and root reuse" begin
    for name in ("s3_abalone","lens_2_1","lens_3_1_mst2_example23","poincare")
        v=vl_example(name); calculate_weight(v;assign=true)
        raw=invariant_expression(v,TaftAlgebra(:N))
        qa=simplify_invariant(raw;strategy=:qa)
        net=simplify_invariant(raw;strategy=:network)
        for n in (3,5)
            base=evaluate_invariant(qa,n)
            @test base==evaluate_invariant(net,n)
            @test base==AkaghefAlgebra.SymbolicTaft.evaluate_exact(qa.expression,qa.context,Dict(:N=>n))
            for m in (1,2)
                z=normalize_root(base,RootSpec(n,m))
                expected=AkaghefAlgebra.TaftNetworks.galois(base,n,m)
                @test canonical_value(z)==expected
                @test ComplexF64(z)≈ComplexF64(expected)
                @test normalize_root(z,RootSpec(n,m+n))===z
            end
        end
        @test canonical_value(Z_MST(v,TaftAlgebra(3)))==evaluate_invariant(qa,3)
    end
end
