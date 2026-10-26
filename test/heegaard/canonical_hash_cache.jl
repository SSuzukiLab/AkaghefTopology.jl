module CanonicalHashCacheTests
using Test, Topology, SHA
const H=Topology.Heegaard
@testset "Canonical hash cache preserves values and notices mutations" begin
    code=H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)])
    expected=bytes2hex(SHA.sha256(H._canonical_json(H._canonical_structural_dict(code))))
    @test H.canonical_hash(code)==expected
    @test H.canonical_hash(deepcopy(code))==expected
    # Changing a sign on the same object must recompute the oriented code.
    code.intersections[1]=H.Intersection(1,1,1,-1)
    changed=bytes2hex(SHA.sha256(H._canonical_json(H._canonical_structural_dict(code))))
    @test changed!=expected
    @test H.canonical_hash(code)==changed
    push!(code.alpha[1].visits,1)
    @test_throws ArgumentError H.canonical_hash(code)
end
end
