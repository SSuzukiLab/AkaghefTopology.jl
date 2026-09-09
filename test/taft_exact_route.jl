using Test, Topology
import AkaghefAlgebra

@testset "Taft exact public contraction" begin
    v=VirtualLink(); set_data!(v;rgauss=[[1,-1]],orientation=[1])
    set_weight!(v,[[-2,1]])
    q=AkaghefAlgebra.TaftExact.zeta(5)
    @test hopf_invariant(v,1//5)==q^-1
    @test hopf_invariant(v,2//5)==q^-2
    reduced=taft_reduce(v)
    @test isempty(reduced.remaining)
    @test taft_evaluate(reduced,9)==AkaghefAlgebra.TaftExact.zeta(9)^-1
    @test_throws ArgumentError taft_evaluate(reduced,9;M=3)
    @test_throws ArgumentError hopf_invariant(v,1//5;style=:R)
    @test_throws ArgumentError hopf_invariant(v,1//1)
    cut=deepcopy(v); push!(cut.gauss[1],0)
    @test_throws ArgumentError hopf_invariant(cut,1//5)
    fractional=deepcopy(v); set_weight!(fractional,[[0.5,0.0]])
    @test_throws ArgumentError hopf_invariant(fractional,1//5)
    h=uqsl2_borel_small(1//5)
    @test hopf_invariant(v,h)==q^-1
    @test ComplexF64(hopf_invariant(v,h))≈hopf_invariant(v,h;method=:numeric)
    empty!(h.positive_vertex.terms)
    @test hopf_invariant(v,h)==0
    @test_throws ArgumentError hopf_invariant(v,h;method=:exact)
    # Generic algebras and the explicitly requested legacy path remain available.
    @test hopf_invariant(v,sweedler_algebra())≈-1
    h4=uqsl2_borel_small(1//4)
    @test ComplexF64(hopf_invariant(v,h4))≈hopf_invariant(v,h4;method=:numeric)
end
