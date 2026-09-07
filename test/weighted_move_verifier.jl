using Test
using Topology
using JSON

@testset "versioned weighted executor evidence" begin
    # Independent explicit lens H incidence: incoming -1, outgoing +1;
    # the self-edge contributions cancel, leaving [0,-1,0,1].
    source=lens_space_link(2,1)
    set_weight!(source,[1,-1,-1,0])
    snapshot=deepcopy(source)
    candidate=only(movable(source,"H";param=1,vertices=[1]))
    result=Topology.verify_weighted_move(source,candidate;framed=true)
    @test result.valid
    if result.valid
        @test result.target.weights==[[1.0,-2.0,-1.0,1.0]]
        @test source.gauss==snapshot.gauss && source.weights==snapshot.weights
        @test source.is_weighted==snapshot.is_weighted
        evidence=Topology.weighted_move_evidence_data(result.evidence)
        @test evidence["identity"]=="raw-weighted-only"
        @test evidence["source_equation"]["residual"]==["0","0","0"]
        @test evidence["target_equation"]["residual"]==["0","0","0"]
        @test Topology.weighted_raw_hash(source)!=Topology.weighted_raw_hash(result.target)
        decoded=Topology.weighted_move_evidence_from_data(JSON.parse(ograph_canonical_json(evidence)))
        @test Topology.replay_weighted_move(decoded).valid
        inverse=only(movable(result.target,"H";param=-1,vertices=[1]))
        back=Topology.verify_weighted_move(result.target,inverse;framed=true)
        @test back.valid && back.target.weights==source.weights

        for mutate in (
            d -> (d["source_raw"]["weights"][1][1]+=1),
            d -> (d["target_raw"]["weights"][1][2]+=1),
            d -> (d["source_raw_sha256"]="forged"),
            d -> (d["target_raw_sha256"]="forged"),
            d -> (d["candidate"]["parameter"]="2"),
            d -> (d["candidate"]["direction"]="inverse"),
            d -> (d["candidate"]["locations"]=[[1,1]]),
            d -> (d["candidate"]["edges"]=[[1,-1]]),
            d -> (d["candidate"]["ps_split_k"]=1),
            d -> (d["source_equation"]["residual"][1]="5"),
            d -> (d["transport_detail"]["transport"]["incidences"][1]["delta"]=-9),
            d -> (d["source_rule_version"]="forged"),
            d -> (d["schema"]="akaghef.ograph/1"),
        )
            tampered=deepcopy(evidence)
            mutate(tampered)
            @test !Topology.replay_weighted_move(tampered).valid
        end

        # Also detect source edits after recomputing its raw hash: they must
        # reproduce target, equation, and occurrence evidence, not only a hash.
        tampered=deepcopy(evidence)
        tampered["source_raw"]["weights"][1][2]+=1
        changed=Topology._weighted_link_from_data(tampered["source_raw"])
        tampered["source_raw_sha256"]=Topology.weighted_raw_hash(changed)
        @test !Topology.replay_weighted_move(tampered).valid

        # Fresh process reads JSON, reconstructs source and current rules, and
        # runs the public replay entry. Artifacts stay in the approved scratch.
        root=normpath(joinpath(@__DIR__,"..","..",".."))
        scratch=joinpath(root,".research-scratch","PJ260830","260908-weighted-implementation","evidence")
        mkpath(scratch)
        file=joinpath(scratch,"weighted-h-replay.json")
        write(file,ograph_canonical_json(evidence))
        project=normpath(joinpath(@__DIR__,".."))
        script="using Topology, JSON; r=Topology.replay_weighted_move(JSON.parsefile(ARGS[1])); r.valid || error(string(r.reason, \" \", r.error)); println(\"weighted-fresh-process-replay-ok\")"
        output=read(`$(Base.julia_cmd()) --startup-file=no --pkgimages=no --compiled-modules=yes --project=$project -e $script $file`,String)
        @test occursin("weighted-fresh-process-replay-ok",output)
    end

    for family in ("CP","BMP","B02","S")
        unsupported=MoveCandidate(family,"",Int[],NTuple{2,Int}[],:forward,nothing,nothing,Tuple{Int,Int}[])
        @test Topology.verify_weighted_move(source,unsupported).reason==:unsupported_move
        @test source.weights==snapshot.weights && source.gauss==snapshot.gauss
    end
    fractional=deepcopy(source); fractional.weights[1][1]=0.5
    @test !Topology.verify_weighted_move(fractional,candidate).valid
    not_framed=deepcopy(source); not_framed.weights[1][1]=2
    @test !Topology.verify_weighted_move(not_framed,candidate;framed=true).valid
    @test Topology.verify_weighted_move(not_framed,candidate;framed=false).valid
    @test !verify_ograph_move(source,candidate).valid # Old schema still refuses weights.
end

@testset "weighted PS recomputes arc evidence" begin
    source=VirtualLink()
    set_data!(source;odata=([[-1,1]],[-1],[[-2,1]]))
    saved=Any[]
    for parameter in 1:4
        candidate=only(movable(source,"PS";param=string(parameter),
            edge_occurrences=[(1,1),(1,1)]))
        result=Topology.verify_weighted_move(source,candidate;framed=true)
        @test result.valid
        if result.valid
            @test Topology.replay_weighted_move(result.evidence).valid
            push!(saved,Topology.weighted_move_evidence_data(result.evidence))
            details=result.evidence.data["transport_detail"]["ps"]
            @test details["p2"]==[1,-7,6,0][parameter]
            @test details["a"]==[-1,3,3,-1][parameter]
            @test result.evidence.data["topological"]["status"]=="verified-by-existing-unweighted-verifier"
            for key in ("p2","a","left_lane","right_lane","dot_count","partial_edge_sum")
                data=Topology.weighted_move_evidence_data(result.evidence)
                data["transport_detail"]["ps"][key]+=1
                @test !Topology.replay_weighted_move(data).valid
            end
            inverse=only(movable(result.target,"PS";param=string(parameter),vertices=[2,3],strict=true))
            restored=Topology.verify_weighted_move(result.target,inverse;framed=true)
            @test restored.valid
            if restored.valid
                @test Topology.replay_weighted_move(restored.evidence).valid
                push!(saved,Topology.weighted_move_evidence_data(restored.evidence))
                @test restored.target.gauss==source.gauss && restored.target.weights==source.weights
            end
        end
        if parameter in (2,4)
            for k in (-5,0,3)
                split=Topology.weighted_ps_candidate(candidate;k=k)
                split_result=Topology.verify_weighted_move(source,split;framed=true)
                @test split_result.valid
                if split_result.valid
                    @test Topology.replay_weighted_move(split_result.evidence).valid
                    push!(saved,Topology.weighted_move_evidence_data(split_result.evidence))
                    @test split_result.evidence.data["candidate"]["ps_split_k"]==k
                    altered=Topology.weighted_move_evidence_data(split_result.evidence)
                    altered["candidate"]["ps_split_k"]=k+1
                    @test !Topology.replay_weighted_move(altered).valid
                    inverse=only(movable(split_result.target,"PS";param=string(parameter),vertices=[2,3],strict=true))
                    restored=Topology.verify_weighted_move(split_result.target,inverse;framed=true)
                    @test restored.valid
                    if restored.valid
                        @test Topology.replay_weighted_move(restored.evidence).valid
                        push!(saved,Topology.weighted_move_evidence_data(restored.evidence))
                        @test restored.target.weights==source.weights
                        restored_split=Topology.weighted_ps_candidate(candidate;
                            k=restored.evidence.data["transport_detail"]["ps"]["k"])
                        again=Topology.verify_weighted_move(restored.target,restored_split;framed=true)
                        @test again.valid && again.target.weights==split_result.target.weights
                    end
                end
            end
        end
    end
    @test length(saved)==20
    if length(saved)==20
        root=normpath(joinpath(@__DIR__,"..","..",".."))
        scratch=joinpath(root,".research-scratch","PJ260830","260908-weighted-implementation","evidence")
        mkpath(scratch)
        file=joinpath(scratch,"weighted-ps-replay.json")
        write(file,ograph_canonical_json(saved))
        project=normpath(joinpath(@__DIR__,".."))
        script="using Topology, JSON; for d in JSON.parsefile(ARGS[1]); r=Topology.replay_weighted_move(d); r.valid || error(string(r.reason, \" \", r.error)); end; println(\"weighted-ps-fresh-process-replay-ok\")"
        output=read(`$(Base.julia_cmd()) --startup-file=no --pkgimages=no --compiled-modules=yes --project=$project -e $script $file`,String)
        @test occursin("weighted-ps-fresh-process-replay-ok",output)
    end
end

@testset "weighted 02 and MP replay preserves full rules" begin
    source=lens_space_link(2,1)
    set_weight!(source,[1,-1,-1,0])
    direct=only(movable(source,"02";param="dir",vertices=[-2,-1]))
    result=Topology.verify_weighted_move(source,direct;framed=true)
    @test result.valid
    if result.valid
        @test Topology.replay_weighted_move(result.evidence).valid
        data=Topology.weighted_move_evidence_data(result.evidence)
        data["candidate"]["after"]["weights"][1][3]+=1
        @test Topology.replay_weighted_move(data).reason==:candidate_mismatch
        data=Topology.weighted_move_evidence_data(result.evidence)
        data["candidate"]["before"]["capacity"][1]+=1
        @test Topology.replay_weighted_move(data).reason==:candidate_mismatch
        inverse=only(movable(result.target,"02";param="inv",vertices=[3,4]))
        restored=Topology.verify_weighted_move(result.target,inverse;framed=true)
        @test restored.valid
        @test Topology.replay_weighted_move(restored.evidence).valid
        @test restored.target.weights==source.weights
        @test restored.target.gauss==source.gauss
    end

    # A real/virtual, multiple-component raw source is kept exactly in evidence;
    # it cannot silently become a one-circuit OGraphState certificate.
    multi=VirtualLink()
    set_data!(multi;gauss=[[1,-1,2,-3],[-2,3]],orientation=[1,1,0])
    # C3's source edge -1 -> 2 is internal and must have weight zero.
    set_weight!(multi,[[-3,1,0],[-1]])
    candidates=movable(multi,"MP";param="C3",vertices=[1,2],strict=true)
    @test length(candidates)==1
    if length(candidates)==1
        invalid_internal=deepcopy(multi)
        invalid_internal.weights[1][3]=1
        @test Topology.verify_weighted_move(invalid_internal,only(candidates)).reason==:application_failed
        @test invalid_internal.weights==[[-3.0,1.0,1.0],[-1.0]]
        result=Topology.verify_weighted_move(multi,only(candidates))
        @test result.valid
        if result.valid
            @test Topology.replay_weighted_move(result.evidence).valid
            @test result.evidence.data["source_raw"]["orientation"]==[1,1,0]
            @test result.evidence.data["topological"]["status"]=="outside-existing-unweighted-verifier-scope"
            @test result.target.gauss==[[3,1,2,-3],[-1,-2]]
            @test result.target.weights==[[-3.0,1.0,0.0,0.0],[-1.0,0.0]]
        end
    end
end
