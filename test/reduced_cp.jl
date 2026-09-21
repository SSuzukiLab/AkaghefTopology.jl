using Test, Topology

function _rcp_link(words, signs)
    v=VirtualLink()
    set_data!(v;rgauss=words,orientation=signs)
end

@testset "rCP 2-to-5 rooted replacement" begin
    # Source and target are the independently established R6 and R0 patches.
    small=[-1,2,-2,1]
    large=[-1,-2,3,-4,2,4,5,-5,1,-3]
    signs=[-1,-1,1,-1,1]
    source=_rcp_link([small],[-1,1])
    candidate=only(movable(source,"rCP";vertices=[1,2],strict=true))
    @test candidate.direction==:forward
    @test candidate.type=="RCP"
    target=move(source,candidate)
    @test real_gauss_code(target)==([large],signs)
    @test real_gauss_code(source)==([small],[-1,1])
    @test certify_closed_ograph(ograph_state(source)) isa ClosedNormalOGraph
    @test certify_closed_ograph(ograph_state(target)) isa ClosedNormalOGraph
    @test length(disk_table(target).ID)==length(disk_table(source).ID)+3
    inplace=deepcopy(source)
    move!(inplace,"rCP",(v=[1,2],p="dir"))
    @test real_gauss_code(inplace)==real_gauss_code(target)

    # Stored cyclic basepoints must not change applicability or the result.
    for shift in 0:3
        rotated=_rcp_link([circshift(small,shift)],[-1,1])
        result=move(rotated,:rCP;vertices=[1,2])
        @test canonicalize_ograph(ograph_state(result))==canonicalize_ograph(ograph_state(target))
    end

    cut=_rcp_link([vcat(small,0)],[-1,1])
    transported=move_boundary_transport(cut,"rCP";vertices=[1,2],strict=true)
    @test real_gauss_code(transported.target)==([vcat(large,0)],signs)
    @test length(transported.endpoints)==2
    @test all(row.source==row.target for row in transported.endpoints)
    signature(v)=Set((p.start.side,p.start.lane,p.stop.side,p.stop.lane,p.dots)
                     for p in circuit_b_boundary(v))
    @test signature(cut)==signature(transported.target)==Set([
        (:in,1,:out,2,0),(:in,2,:out,1,0),(:out,3,:in,3,4)])

    # A cyclically equivalent closed word is not the same rooted open patch.
    bad=_rcp_link([[1,-1,2,-2,0]],[-1,1])
    @test isempty(movable(bad,"rCP";vertices=[1,2]))
    snapshot=deepcopy(real_gauss_code(bad))
    @test_throws ArgumentError move!(bad,"rCP";vertices=[1,2])
    @test real_gauss_code(bad)==snapshot
    @test_throws ArgumentError move!(bad,candidate)
    @test real_gauss_code(bad)==snapshot
    @test isempty(movable(_rcp_link([small],[1,-1]),"rCP";vertices=[1,2],strict=true))
    @test isempty(movable(_rcp_link([small],[1,1]),"rCP";vertices=[1,2]))
    @test isempty(movable(_rcp_link([[-1,2,0],[-2,1,0]],[-1,1]),"rCP";vertices=[1,2]))
    @test_throws ArgumentError move(target,"rCP";vertices=1:5)
end

@testset "rCP retains exterior and respects selected vertex roles" begin
    # a=3, b=1; vertex 2 and the second cut strand belong to the exterior.
    source=_rcp_link([[2,-3,1,-1,3,-2,0],[4,-4,0]],[1,1,-1,1])
    @test isempty(movable(source,"rCP";vertices=[1,3],strict=true))
    candidate=only(movable(source,"rCP";vertices=[1,3]))
    @test candidate.vertices==[3,1]
    target=move(source,candidate)
    @test real_gauss_code(target)==(
        [[2,-3,-1,5,-6,1,6,7,-7,3,-5,-2,0],[4,-4,0]],
        [-1,1,-1,1,1,-1,1])
    @test strand_cuts(target)==strand_cuts(source)
    @test length(real_gauss_code(target)[2])==length(real_gauss_code(source)[2])+3
end

@testset "rCP has no implicit integer or spin weight transport" begin
    source=_rcp_link([[-1,2,-2,1]],[-1,1])
    candidate=only(movable(source,"rCP";vertices=[1,2],strict=true))
    for mode in (:integer,:z2)
        weighted=deepcopy(source)
        set_weight!(weighted,[0,0,0,0];mode=mode)
        snapshot=(deepcopy(real_gauss_code(weighted)),deepcopy(weighted.weights),weight_mode(weighted))
        @test isempty(movable(weighted,"rCP";vertices=[1,2]))
        @test_throws ArgumentError move!(weighted,"rCP";vertices=[1,2])
        @test_throws ArgumentError move!(weighted,candidate)
        @test (real_gauss_code(weighted),weighted.weights,weight_mode(weighted))==snapshot
    end
end
