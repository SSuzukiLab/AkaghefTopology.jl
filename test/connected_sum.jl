using Test
using Topology

@testset "connected-sum integral weight lift" begin
    weighted_sum=VirtualLink()
    set_data!(weighted_sum;rgauss=[[1,-1]],orientation=[1])
    solve_integer_weights(weighted_sum;assign=true)
    weighted_copy=deepcopy(weighted_sum)

    disjoint_sum!(weighted_sum,weighted_copy)
    @test weighted_sum.is_weighted
    connected_sum!(weighted_sum)
    @test weighted_sum.is_weighted
    @test BigInt.(only(weighted_sum.weights))==[-2,1,0,0,1,0,0,1,-2,0,1,-1]
    @test is_framed_weight(weighted_sum)
    @test is_closed(weighted_sum)==(true,(true,true,true))

    anchored_weighted=deepcopy(weighted_copy)
    disjoint_sum!(anchored_weighted,weighted_copy)
    connected_sum!(anchored_weighted,[1 2; -1 -2])
    @test is_framed_weight(anchored_weighted)
    @test BigInt.(only(anchored_weighted.weights))==[1,-2,0,0,1,0,0,-2,1,0,1,-1]

    discarded=deepcopy(weighted_copy)
    disjoint_sum!(discarded,weighted_copy)
    connected_sum!(discarded;transport_weights=false)
    @test !discarded.is_weighted

    lift_tangle=VirtualLink()
    set_data!(lift_tangle;rgauss=[[-2,-1,1,4],[-4,-3,3,2]],
              orientation=[-1,1,-1,1],cut=[true,true])
    for (u,v) in ((0,0),(1,0),(0,1),(-3,5))
        local_weights=[0,0,1,0,0,v,0,1,-1,u]
        set_weight!(lift_tangle,[local_weights[1:5],local_weights[6:10]])
        boundary=boundary_weight_data(lift_tangle)
        @test boundary.closed_residual[[2,4,8]]==[0,0,0]
        @test boundary.measured[[1,5,6]]==[2u,2u,-2u]
        @test boundary.measured[[3,7,9]]==[2v,2v,-2v]
    end

    mutated_weights=[0,0,1,0,0,0,0,1,0,0]
    set_weight!(lift_tangle,[mutated_weights[1:5],mutated_weights[6:10]])
    @test any(!iszero,boundary_weight_data(lift_tangle).closed_residual[[2,4,8]])
end
