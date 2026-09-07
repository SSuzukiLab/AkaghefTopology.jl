using Test
using Topology

_transport_raw(v)=(deepcopy(v.gauss),copy(v.orientation),deepcopy(v.weights),v.real_only,
                   v.is_weighted,copy(v.strand_types))
function _transport_cycles(v)
    code,ori=real_gauss_code(v)
    cycles=[isempty(g) ? () : minimum(Tuple(circshift(collect(zip(g,w)),k))
                                      for k in 0:length(g)-1) for (g,w) in zip(code,v.weights)]
    (sort!(cycles),ori)
end
function _transport_link(code,orientation,weights)
    v=VirtualLink(); set_data!(v;odata=(code,orientation,weights)); v
end

@testset "weighted 02 boundary-edge transport" begin
    source=_transport_link([[1,-1]],[1],[[5,7]])
    snapshot=_transport_raw(source)
    # The old selected edge weight goes on the incoming external portion;
    # the five added local portions are exactly [0,0,1,-1,0].
    for (pair,g,w) in (([1,-1],[[1,-2,-3,3,2,-1]],[[5,7,0,1,-1,0]]),
                       ([-1,1],[[1,-1,-2,-3,3,2]],[[0,7,5,0,1,-1]]))
        candidate=only(movable(source,"02";param="dir",vertices=pair))
        target=move(source,candidate)
        @test target.gauss==g
        @test target.weights==w
        @test target.is_weighted
        @test _transport_raw(source)==snapshot
        inverse=only(movable(target,"02";param="inv",vertices=[2,3],strict=true))
        @test _transport_cycles(move(target,inverse))==_transport_cycles(source)
        mutating=deepcopy(source); move!(mutating,candidate)
        @test _transport_raw(mutating)==_transport_raw(target)
        for visit in (-3,3,2)
            broken=deepcopy(target)
            broken.weights[1][findfirst(==(visit),only(broken.gauss))]+=1
            prior=_transport_raw(broken)
            @test isempty(movable(broken,"02";param="inv",vertices=[2,3],strict=true))
            @test_throws ArgumentError move!(broken,inverse)
            @test _transport_raw(broken)==prior
        end
        detail=Topology._weighted_transport_detail(source,candidate)
        @test detail.kind==:fragment_transport
        @test count(f -> f.kind==:exterior,detail.fragments)==1
    end
    patch=_transport_link([[-1,-2,2,1]],[1,-1],[[0,0,1,-1]])
    @test move(patch,"02";param="inv",vertices=[1,2]).weights==[Float64[]]
    loop=_transport_link([[-1,-2,2,1]],[1,-1],[[3,0,1,-1]])
    prior=_transport_raw(loop)
    @test_throws ArgumentError move!(loop,"02";param="inv",vertices=[1,2])
    @test _transport_raw(loop)==prior
end

@testset "MP exterior portions, repeated component and cyclic seams" begin
    data=make_move_data()
    @test length(data["MP-L"])==length(data["MP-R"])==16
    for before in data["MP-L"]
        after=only(filter(p -> p.parameter==before.parameter,data["MP-R"]))
        free=[5,8,11]
        weights=[vcat(free[i]+first(w)+last(w),w[2:end-1]) for (i,w) in enumerate(before.weights)]
        source=_transport_link(before.gauss,before.orientation,weights)
        candidate=only(movable(source,"MP";param=before.parameter,vertices=[1,2],strict=true))
        target=move(source,candidate)
        expected=[vcat(free[i]+first(w)+last(w),w[2:end-1]) for (i,w) in enumerate(after.weights)]
        @test target.gauss==after.gauss
        @test target.weights==expected
        @test _transport_cycles(move(target,"MP";param=before.parameter,vertices=[1,2,3]))==_transport_cycles(source)
        joined=reduce(vcat,before.gauss)
        joined_weights=reduce(vcat,[vcat(free[i]+last(before.weights[mod1(i-1,3)])+first(w),w[2:end-1])
                                   for (i,w) in enumerate(before.weights)])
        for shift in 0:length(joined)-1
            source=_transport_link([circshift(joined,shift)],before.orientation,[circshift(joined_weights,shift)])
            candidate=only(movable(source,"MP";param=before.parameter,vertices=[1,2],strict=true))
            target=move(source,candidate)
            restored=move(target,"MP";param=before.parameter,vertices=[1,2,3])
            @test _transport_cycles(restored)==_transport_cycles(source)
            @test target.is_weighted
        end
    end
end

@testset "H integer parameters and occurrence multiplicity" begin
    source=_transport_link([[1,-1]],[1],[[5,7]])
    candidate=only(movable(source,"H";param=23,vertices=[1]))
    @test move(source,candidate).weights==source.weights
    @test length(Topology._weighted_transport_detail(source,candidate).incidences)==4
    lens=lens_space_link(2,1); set_weight!(lens,[1,-1,-1,0])
    before=_transport_raw(lens)
    move!(lens,"H";param=-7,vertices=[2])
    @test reduce(vcat,lens.weights)==[1,-8,-1,7]
    move!(lens,"H";param=7,vertices=[2])
    @test _transport_raw(lens)==before
    for bad in (1.5,NaN,Inf,big(2)^53,big(2)^100)
        @test_throws ArgumentError move!(lens,"H";param=bad,vertices=[2])
        @test _transport_raw(lens)==before
    end
    limit=Float64(big(2)^53-1)
    large=_transport_link([[-2,-1,1,2]],[1,1],[[0,limit,0,0]])
    prior=_transport_raw(large)
    @test_throws ArgumentError move!(large,"H";param=1,vertices=[2])
    @test _transport_raw(large)==prior
end

@testset "generic capacities retain exterior edge occurrences" begin
    # Mechanical converter fixtures, not additional mathematically admitted
    # move families: exchange the two caps to join and split two loops.
    before=MovePattern("transport",[[1],[-1]],[1],[1,2])
    after=MovePattern("transport",[[1],[-1]],[1],[2,1])
    candidate=MoveCandidate("transport","transport",[1],NTuple{2,Int}[],:forward,
        before,after,[(1,1),(2,1)])
    mapping=Dict(1=>1)
    code,weights,ori=Topology._reconnect_pattern([[1],[-1]],[[5.0],[7.0]],candidate,
        before.gauss,after.gauss,mapping,[1];weighted=true)
    joined=_transport_link(code,ori,weights)
    @test _transport_cycles(joined)==_transport_cycles(_transport_link([[1,-1]],[1],[[5,7]]))
    locations=[(1,findfirst(==(signed),only(code))) for signed in (1,-1)]
    inverse=MoveCandidate("transport","transport",[1],NTuple{2,Int}[],:inverse,
        after,before,locations)
    code,weights,ori=Topology._reconnect_pattern(code,weights,inverse,
        after.gauss,before.gauss,mapping,[1];weighted=true)
    @test _transport_cycles(_transport_link(code,ori,weights))==
          _transport_cycles(_transport_link([[1],[-1]],[1],[[5],[7]]))

    before=MovePattern("transport",[[1,-1]],[1],[1])
    after=MovePattern("transport",[[1,-1],[2,-2]],[1,1],[1,0];weights=[[0,0,0],[2,3,5]])
    candidate=MoveCandidate("transport","transport",[1],NTuple{2,Int}[],:forward,
        before,after,[(1,1)])
    code,weights,ori=Topology._reconnect_pattern([[1,-1]],[[11.0,0.0]],candidate,
        before.gauss,after.gauss,mapping,[1,1];weighted=true)
    @test weights==[[11,0],[7,3]]
    inverse=MoveCandidate("transport","transport",[1,2],NTuple{2,Int}[],:inverse,
        after,before,[(1,1),(2,1)])
    restored=Topology._reconnect_pattern(code,weights,inverse,after.gauss,before.gauss,mapping,[1];weighted=true)
    @test restored[2]==[[11,0]]
    weights[2][1]+=1
    @test_throws ArgumentError Topology._reconnect_pattern(code,weights,inverse,
        after.gauss,before.gauss,mapping,[1];weighted=true)
end

@testset "weighted capability and frozen candidate safety" begin
    source=_transport_link([[1,-1]],[1],[[5,7]])
    candidate=only(movable(source,"02";param="dir",vertices=[1,-1]))
    for type in ("CP","BMP","B02","S","PS")
        prior=_transport_raw(source)
        @test isempty(movable(source,type;vertices=[1,2]))
        @test_throws ArgumentError move!(source,type;vertices=[1,2])
        @test _transport_raw(source)==prior
    end
    for weights in ([[NaN,0]],[[Inf,0]],[[0.5,0]],[[Float64(big(2)^53),0]],[[1]],[[1,2],[]])
        bad=deepcopy(source); bad.weights=weights
        prior=_transport_raw(bad)
        @test_throws Union{ArgumentError,DimensionMismatch} move!(bad,candidate)
        @test isequal(_transport_raw(bad),prior)
    end
    for bad in (MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,
                    candidate.edges,candidate.direction,candidate.before,candidate.after,[(1,2),(1,1)]),
                MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,
                    candidate.edges,:inverse,candidate.before,candidate.after,candidate.locations),
                MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,
                    candidate.edges,candidate.direction,candidate.before,candidate.after,candidate.locations,big(2)))
        prior=_transport_raw(source)
        @test_throws ArgumentError move!(source,bad)
        @test _transport_raw(source)==prior
    end
    malformed=deepcopy(source); malformed.gauss=[[1,1]]
    @test_throws ArgumentError move!(malformed,"H";param=1,vertices=[1])
    cut=deepcopy(source); push!(cut.gauss[1],0)
    prior=_transport_raw(cut)
    @test_throws ArgumentError move!(cut,candidate)
    @test _transport_raw(cut)==prior
end
