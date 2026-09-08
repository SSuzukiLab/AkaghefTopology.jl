using Test
using Topology

@testset "weighted PS on an empty interval retains another cut" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[0],[0]],orientation=Int[])
    set_weight!(source,[[2],[5]])
    measured=boundary_weight_data(source)
    targets=Dict(id=>value for (id,value) in zip(measured.boundary_identities,measured.measured))
    candidate=only(movable(Topology._weighted_unweight(source),"PS";
        param="4",edge_occurrences=[(1,1),(1,1)]))
    result=verify_weighted_move(source,candidate;framed=true,
        boundary=(source=targets,target=targets))
    @test result.valid
    if result.valid
        @test result.target.gauss[2]==[0]
        @test result.target.weights[2]==[5]
        @test replay_weighted_move(result.evidence).valid
        inverse=only(movable(Topology._weighted_unweight(result.target),"PS";param="4",vertices=[1,2]))
        back=move(result.target,inverse;boundary=(source=targets,target=targets))
        @test back.gauss==source.gauss
        @test back.weights==source.weights
    end
end

function _evidence_bmp_source()
    source=VirtualLink()
    set_data!(source;rgauss=[[1,2,0],[-1,0],[-2,0]],orientation=[-1,1])
    set_weight!(source,[[10,20,30],[40,50],[60,70]])
    source
end

function _evidence_mp_source()
    source=VirtualLink()
    set_data!(source;rgauss=[[-1,2,0],[1,0],[-2,0]],orientation=[1,1])
    set_weight!(source,[[11,0,13],[21,22],[31,32]])
    source
end

function _evidence_endpoint_map(data)
    rows=data["transport_detail"]["boundary_transport"]["endpoints"]
    Dict((row["source"]["strand"],row["source"]["side"]) =>
         (row["target"]["strand"],row["target"]["side"]) for row in rows)
end

@testset "weighted multi-cut endpoint lineage" begin
    source=_evidence_bmp_source()
    geometric=Topology._weighted_unweight(source)
    candidate=only(movable(geometric,"BMP";param=0,vertices=[1,2],strict=true))
    transport=move_boundary_transport(geometric,candidate)
    target=transport.target
    set_weight!(target,[zeros(length(word)) for word in target.gauss])

    result=verify_weighted_relative_result(source,target,candidate)
    @test result.valid
    if result.valid
        evidence=weighted_move_evidence_data(result.evidence)
        endpoints=_evidence_endpoint_map(evidence)
        @test endpoints[(2,"out")] == (3,"out")
        @test endpoints[(3,"out")] == (2,"out")
        @test evidence["source_equation"]["boundary_targets"] === nothing
        @test evidence["source_equation"]["measured_boundary_residual"] isa Vector
        @test evidence["target_equation"]["boundary_targets"] === nothing
        @test replay_weighted_relative_result(evidence).valid

        tampered=deepcopy(evidence)
        row=only(filter(item -> item["source"]["strand"]==2 &&
                        item["source"]["side"]=="out",
                     tampered["transport_detail"]["boundary_transport"]["endpoints"]))
        row["target"]["strand"]=1
        @test !replay_weighted_relative_result(tampered).valid
    end
end

@testset "weighted multi-cut executor records endpoint lineage" begin
    source=_evidence_mp_source()
    candidate=only(movable(Topology._weighted_unweight(source),"MP";
        param="C3",vertices=[1,2],strict=true))
    result=verify_weighted_move(source,candidate)
    @test result.valid
    if result.valid
        data=weighted_move_evidence_data(result.evidence)
        boundary_transport=data["transport_detail"]["boundary_transport"]
        @test boundary_transport["kind"]=="endpoint-lineage"
        @test boundary_transport["source_cut_count"]==3
        @test boundary_transport["target_cut_count"]==3
        @test replay_weighted_move(data).valid

        tampered=deepcopy(data)
        row=only(filter(item -> item["source"]["strand"]==2 &&
                        item["source"]["side"]=="out",
                     tampered["transport_detail"]["boundary_transport"]["endpoints"]))
        row["target"]["side"]="in"
        @test !replay_weighted_move(tampered).valid
    end
end
