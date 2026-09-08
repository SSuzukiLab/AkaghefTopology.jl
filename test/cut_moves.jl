using Test
using Topology

function _cut_fixture()
    v=VirtualLink()
    set_data!(v;rgauss=[[-1,2,-2,3,-3,1,0]],orientation=[1,1,1])
    v
end

@testset "cut move coordinates preserve the open strand" begin
    source=_cut_fixture()
    candidate=only(movable(source,"MP";param="C3",vertices=[1,2],strict=true))
    target=move(source,candidate)
    @test count(identity,Topology.strand_cuts(target))==1
    target_code=first(real_gauss_code(target))
    @test last(only(target_code))==0
    @test all(!=(0),only(target_code)[1:end-1])
end

@testset "cut candidates do not cross the boundary marker" begin
    source=_cut_fixture()
    word=only(gauss_code(source))
    @test isempty(Topology._pattern_occurrences(word,[1,-1]))
    @test Topology._pattern_occurrences(word,[2,-2])==[2]
end

@testset "H transports endpoint portions on an open strand" begin
    source=_cut_fixture()
    source.weights=[[0,1,2,3,4,5,6]]
    source.is_weighted=true
    move!(source,"H";param=1,vertices=[1])
    @test count(identity,Topology.strand_cuts(source))==1
    source_code=first(real_gauss_code(source))
    @test last(only(source_code))==0
    @test length(only(source.weights))==length(only(source_code))
end

@testset "weighted MP transports distinct cut edge portions" begin
    source=_cut_fixture()
    geometric=move(source,"MP";param="C3",vertices=[1,2])
    # C3's second source fragment has outgoing rule 0, replaced by -1.
    # Thus the outgoing cut portion 7 becomes 6; it must not be discarded
    # and restored unchanged after a cyclic weight calculation.
    set_weight!(source,[2,0,3,4,5,6,7])
    target=move(source,"MP";param="C3",vertices=[1,2])
    @test target.gauss==geometric.gauss
    @test first(only(target.weights))==2
    @test last(only(target.weights))==6
    @test length(only(target.weights))==length(real_edge_table(target))
    @test only(source.weights)==[2,0,3,4,5,6,7]
end

@testset "weighted raw evidence retains cut endpoint data" begin
    source=_cut_fixture()
    source.weights=[[0,1,2,3,4,5,6]]
    source.is_weighted=true
    candidate=only(movable(source,"H";param="1",vertices=[1],strict=true))
    targets=zeros(BigInt,length(disk_table(source).ID))
    report=verify_weighted_move(source,candidate;
        boundary=(source=targets,target=targets))
    @test report.valid
    data=weighted_move_evidence_data(report.evidence)
    @test data["source_raw"]["cuts"]==[true]
    @test data["source_raw"]["boundary_weights"]==[[0,6]]
    @test replay_weighted_move(data).valid
end

@testset "weighted cut PS uses explicit measured boundary targets" begin
    source=_cut_fixture()
    source.weights=[[0,0,1,0,1,0,0]]
    source.is_weighted=true
    geometric=Topology._ps_unweighted(source)
    candidate=only(movable(geometric,"PS";param="4",
        edge_occurrences=[(1,1),(1,1)]))
    disks=disk_table(source)
    open_rows=[i for i in eachindex(disks.ID) if !disks.isCycle[i]]
    values=[BigInt(2),BigInt(0),BigInt(0)]
    targets=Dict((disks.boundary_start[row],disks.boundary_end[row])=>values[index]
                 for (index,row) in enumerate(open_rows))
    before=(gauss=deepcopy(source.gauss),weights=deepcopy(source.weights))
    report=verify_weighted_move(source,candidate;framed=true,
        boundary=(source=targets,target=targets))
    @test report.valid
    @test source.gauss==before.gauss && source.weights==before.weights
    if report.valid
        @test replay_weighted_move(report.evidence).valid
        @test report.target.gauss==[[-5,-4,4,5,-1,2,-2,3,-3,1,0]]
        @test report.target.weights==[[1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,0.0,1.0,0.0,0.0]]
        bad_targets=Dict(key=>(value==2 ? BigInt(4) : value)
                         for (key,value) in targets)
        bad=verify_weighted_move(source,candidate;framed=true,
            boundary=(source=targets,target=bad_targets))
        @test !bad.valid
    end
end

@testset "relative assigned results bind the unweighted target topology" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[-1,1,-2,2,-3,3]],orientation=[-1,-1,-1])
    source.weights=[zeros(length(only(gauss_code(source))))]
    source.is_weighted=true
    geometric=Topology._weighted_unweight(source)
    candidate=only(movable(geometric,"CP";vertices=[1,2,3],strict=true))
    expected=deepcopy(geometric)
    Topology._apply_candidate!(expected,candidate)
    expected.weights=[zeros(length(word)) for word in expected.gauss]
    expected.is_weighted=true
    report=Topology.verify_weighted_relative_result(source,expected,candidate)
    @test report.valid
    wrong=Topology.verify_weighted_relative_result(source,source,candidate)
    @test !wrong.valid
    @test wrong.reason==:raw_target_mismatch
end
