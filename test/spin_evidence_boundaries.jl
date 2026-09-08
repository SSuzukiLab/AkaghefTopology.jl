using Test
using Topology

@testset "Z2 state cannot enter integer or unweighted evidence" begin
    integer=lens_space_link(2,1)
    set_weight!(integer,[1,-1,-1,0])
    h=only(movable(integer,"H";param=1,vertices=[1]))
    baseline=Topology.verify_weighted_move(integer,h;framed=true)
    @test baseline.valid
    spin=deepcopy(integer)
    set_weight_mode!(spin,:z2)
    snapshot=deepcopy(spin)
    @test_throws ArgumentError Topology.weighted_raw_hash(spin)
    @test Topology.verify_weighted_move(spin,h).reason==:unsupported_weight_mode
    @test Topology.verify_weighted_move(spin,h;framed=true).reason==:unsupported_weight_mode
    cp=MoveCandidate("CP","",Int[],NTuple{2,Int}[],:forward,nothing,nothing,Tuple{Int,Int}[])
    @test Topology.verify_weighted_relative_result(spin,integer,cp).reason==:unsupported_weight_mode
    @test Topology.verify_weighted_relative_result(integer,spin,cp).reason==:unsupported_weight_mode
    @test_throws ArgumentError ograph_state(spin)
    unassigned=deepcopy(spin)
    unassigned.is_weighted=false
    @test_throws ArgumentError ograph_state(unassigned)
    @test spin.weights==snapshot.weights && weight_mode(spin)==:z2
    for family in ("H","MP","PS","02","BMP","B02")
        @test isempty(movable(spin,family))
        @test_throws ArgumentError move!(spin,family)
        @test spin.gauss==snapshot.gauss && spin.weights==snapshot.weights && weight_mode(spin)==:z2
    end
    @test_throws ArgumentError move!(spin,h)

    # These existing schemas certify integer arithmetic only. A domain tag
    # cannot be smuggled into them, even if all numeric entries are 0 or 1.
    evidence=Topology.weighted_move_evidence_data(baseline.evidence)
    tagged=deepcopy(evidence)
    tagged["source_raw"]["weight_mode"]="z2"
    @test !Topology.replay_weighted_move(tagged).valid
    @test_throws ArgumentError Topology.weighted_move_evidence_from_data(tagged)
    @test Topology.replay_weighted_move(evidence).valid

    # Explicit internal forgetting is safe for topological transport only.
    forgotten=Topology._weighted_unweight(spin)
    @test !forgotten.is_weighted && weight_mode(forgotten)==:integer
    @test all(iszero,Iterators.flatten(forgotten.weights))
end
