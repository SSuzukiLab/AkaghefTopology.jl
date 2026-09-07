using Test, Topology

@testset "Unweighted evidence cannot silently forget a weighted PS split" begin
    plain=VirtualLink(); set_data!(plain;rgauss=[[-1,1]],orientation=[-1])
    candidate=only(movable(plain,"PS";param="2",edge_occurrences=[(1,1),(1,1)]))
    report=verify_ograph_move(plain,weighted_ps_candidate(candidate;k=3))
    @test !report.valid
    @test report.reason==:weighted_split_unsupported
end

@testset "Integral PS source figure and directed partial arcs" begin
    # Original source figure transcription, catalogue-oriented weights:
    # I/II bottom-left=-a, top-left=a+k, both internal weights=0;
    # III top-left=a+k, internal right=1, bottom-left=-a-1;
    # IV top-right=l+1, internal left=-1, bottom-left=-a.
    v=VirtualLink(); set_data!(v;odata=([[-1,1]],[-1],[[-2,1]]))
    @test Topology._ps_check_framed(v)
    expected=(
        (p2=1,a=-1,weights=[0,0,1,0,-3,1],nodes=[6,4],dots=[1,0]),
        (p2=-7,a=3,weights=[-3,0,1,0,0,1],nodes=[6,4,2],dots=[1,0,0]),
        (p2=6,a=3,weights=[0,0,1,1,-4,1],nodes=[1,3,6,4],dots=[0,1,1,0]),
        (p2=0,a=-1,weights=[1,-1,1,0,-3,1],nodes=Int[],dots=Int[]))
    for p in 1:4
        c=only(movable(Topology._ps_unweighted(v),"PS";param=string(p),
                       edge_occurrences=[(1,1),(1,1)]))
        plan=Topology._weighted_ps_plan(v,c)
        oracle=expected[p]
        @test (plan.detail.p2,plan.detail.a)==(oracle.p2,oracle.a)
        @test plan.detail.nodes==oracle.nodes
        @test plan.detail.dots==oracle.dots
        @test only(plan.target.weights)==oracle.weights
        @test plan.detail.k+plan.detail.l==-2
        @test plan.detail.full_rotation2==-2
        @test Topology._ps_check_framed(plan.target)
        inverse=only(movable(Topology._ps_unweighted(plan.target),"PS";
                             param=string(p),vertices=[2,3],strict=true))
        restored=Topology._weighted_ps_plan(plan.target,inverse)
        @test restored.target.gauss==v.gauss
        @test restored.target.weights==v.weights
        # A source internal edge is prescribed by the original figure. It
        # cannot be discarded merely because an unweighted inverse is legal.
        tampered=deepcopy(plan.target)
        component,start=inverse.locations[1]
        tampered.weights[component][mod1(start+1,length(tampered.gauss[component]))]+=1
        @test_throws ArgumentError Topology._weighted_ps_plan(tampered,inverse)
        # k/l are weights of portions, with k+l equal to the one source edge.
        # An inverse must also recognize legal nondefault splits.
        for k in (-5,0,3)
            split_candidate=Topology.weighted_ps_candidate(c;k=k)
            split_plan=Topology._weighted_ps_plan(v,split_candidate)
            @test move(v,split_candidate).weights==split_plan.target.weights
            split_inverse=only(movable(Topology._ps_unweighted(split_plan.target),"PS";
                                       param=string(p),vertices=[2,3],strict=true))
            split_back=Topology._weighted_ps_plan(split_plan.target,split_inverse)
            @test split_back.target.weights==v.weights
            @test split_back.detail.k+split_back.detail.l==-2
            if p in (2,4)
                @test split_back.detail.k==k
            end
            frozen=Topology.weighted_ps_candidate(c;k=split_back.detail.forward_split_k)
            @test move(split_back.target,frozen).weights==split_plan.target.weights
        end
    end
    @test Topology._ps_integral_a("1",-5)==2
    @test Topology._ps_integral_a("2",7)==-4
    @test Topology._ps_integral_a("3",-6)==-3
    @test Topology._ps_integral_a("4",8)==-5
    @test_throws ArgumentError Topology._ps_integral_a("1",0)
    @test_throws ArgumentError Topology._ps_integral_a("3",1)
    @test v.gauss==[[-1,1]] && v.weights==[[-2,1]]
end

@testset "Integral PS public application and edge occurrence coverage" begin
    source=VirtualLink(); set_data!(source;odata=([[-1,1]],[-1],[[-2,1]]))
    seed=move(source,only(movable(source,"PS";param="1",edge_occurrences=[(1,1),(1,1)])))
    counts=zeros(Int,4,2)
    for p in 1:4, left in 1:6, right in 1:6
        geometry=movable(Topology._ps_unweighted(seed),"PS";param=string(p),
                         edge_occurrences=[(1,left),(1,right)])
        for candidate in geometry
            target=move(seed,candidate)
            inverse=only(movable(target,"PS";param=string(p),vertices=[4,5],strict=true))
            restored=move(target,inverse)
            @test restored.gauss==seed.gauss && restored.weights==seed.weights
            counts[p,left==right ? 2 : 1]+=1
        end
    end
    # Finite acceptance counts, not a proof that arbitrary representations
    # admit PS or that a bounded enumerator is universally complete.
    @test counts==[8 2;5 1;17 3;6 6]
    for p in 1:4
        shifted=VirtualLink(); set_data!(shifted;odata=([[1,-1]],[-1],[[1,-2]]))
        candidate=only(movable(shifted,"PS";param=string(p),edge_occurrences=[(1,2),(1,2)]))
        target=move(shifted,candidate)
        inverse=only(movable(target,"PS";param=string(p),vertices=[2,3],strict=true))
        restored=move(target,inverse)
        @test restored.gauss==shifted.gauss && restored.weights==shifted.weights
    end
end
