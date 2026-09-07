using Test, Topology, JSON

@testset "PS four types: actual insertion, inverse, and global rejection" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[-1,1]],orientation=[-1])
    @test first(is_closed(source))
    for p in 1:4
        candidates=movable(source,"PS";param=string(p),edge_occurrences=[(1,1),(1,1)])
        @test length(candidates)==1
        candidate=only(candidates)
        target=move(source,candidate)
        forward=verify_ograph_move(source,candidate)
        @test forward.valid
        @test replay_ograph_move(ograph_move_witness_data(forward.witness)).valid
        inverse=only(movable(target,"PS";param=string(p),vertices=[2,3],strict=true))
        @test verify_ograph_move(target,inverse).valid
        @test replay_ograph_move(ograph_move_witness_data(verify_ograph_move(target,inverse).witness)).valid
        @test canonicalize_ograph(ograph_state(move(target,inverse)))==canonicalize_ograph(ograph_state(source))
    end
    # Honest local candidate on a closed source, but its specified lanes lie
    # on different Method-B circuits. This is not a corrupted certificate.
    data=make_move_data()
    bad=MoveCandidate("PS","1",Int[],NTuple{2,Int}[],:forward,
        only(data["PS-L"]),data["PS-R"][1],[(1,2),(1,2)])
    condition=ps_global_condition(source,bad)
    @test !condition.valid
    @test condition.reason==:different_method_b_circuits
    @test isempty(movable(source,"PS";param="1",edge_occurrences=[(1,2),(1,2)]))
end

@testset "PS-III actual application rejected on a closed input" begin
    source=VirtualLink()
    set_data!(source;rgauss=[[-1,2,-2,1]],orientation=[-1,1])
    @test is_closed(source)==(true,(true,true,true))
    @test_throws ArgumentError move(source,"PS";param="3",
        edge_occurrences=[(1,1),(1,2)])
    data=make_move_data()
    candidate=MoveCandidate("PS","3",Int[],NTuple{2,Int}[],:forward,
        only(data["PS-L"]),data["PS-R"][3],[(1,1),(1,2)])
    validation=Topology.validate_ps_move(source,candidate)
    @test validation.reason==:target_not_closed
    @test validation.target_conditions==(true,true,false)
    @test real_gauss_code(validation.target)==
        ([[3,4,-1,-4,-3,2,-2,1]],[-1,1,1,-1])
    @test ps_global_condition(source,candidate).reason==:different_method_b_circuits
    @test real_gauss_code(source)==([[-1,2,-2,1]],[-1,1])
    # Existing fixed-store row 193: the legacy vertex-selected deletion entry
    # also rejects a genuine local PS-III occurrence on a closed input.
    inverse_source=VirtualLink()
    set_data!(inverse_source;rgauss=[[-1,-2,2,1,-3,-4,4,3]],
        orientation=[-1,1,1,-1])
    @test first(is_closed(inverse_source))
    @test_throws ArgumentError move(inverse_source,"PS";param="3",
        vertices=[2,1],strict=true)
    inverse=MoveCandidate("PS","3",[2,1],NTuple{2,Int}[],:inverse,
        data["PS-R"][3],only(data["PS-L"]),[(1,3),(1,1)])
    inverse_validation=Topology.validate_ps_move(inverse_source,inverse)
    @test inverse_validation.reason==:target_not_closed
    @test real_gauss_code(inverse_validation.target)==([[-1,-2,2,1]],[1,-1])
    @test ps_global_condition(inverse_source,inverse).reason==:different_method_b_circuits
end

@testset "INS S is a separate zero-to-four move" begin
    source=VirtualLink(); set_data!(source;rgauss=[[-1,1]],orientation=[-1])
    candidate=only(movable(source,"S";edge_occurrences=[(1,1)]))
    target=move(source,candidate)
    report=verify_ograph_move(source,candidate)
    @test report.valid
    @test length(ograph_state(target).vertex_signs)==5
    @test replay_ograph_move(ograph_move_witness_data(report.witness)).valid
    inverse=only(movable(target,"S";vertices=[2,3,4,5],strict=true))
    @test verify_ograph_move(target,inverse).valid
    @test canonicalize_ograph(ograph_state(move(target,inverse)))==canonicalize_ograph(ograph_state(source))
end

@testset "INS S registers only in the Area-3 manifold profile" begin
    @test Topology.MoveGraphStore.DomainTaxonomy.validate_topological_move("S",3)
    @test_throws ArgumentError Topology.MoveGraphStore.DomainTaxonomy.validate_topological_move("S",2)
    store=open_store()
    try
        migrate_ograph_store!(store)
        stamp="2026-09-08T00:00:00Z"
        presentation=OGraphRelationProfile("presentation";moves_complete=true,
            evidence_status="identity-only",created_at_utc=stamp)
        combed=OGraphRelationProfile("combed";parent_profile_id=presentation.profile_id,
            created_at_utc=stamp)
        manifold=OGraphRelationProfile("manifold";profile_version=2,
            parent_profile_id=combed.profile_id,created_at_utc=stamp)
        for profile in (presentation,combed,manifold)
            register_ograph_profile!(store,profile)
        end
        @test "S" in manifold.moves
        @test !("S" in combed.moves)
        @test verify_ograph_store(store).profiles==3
    finally
        close_store!(store)
    end
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__) && !isempty(ARGS)
    @testset "All MP types on small existing closed states" begin
        store=open_store(ARGS[1];migrate=false)
        covered=Set{Tuple{String,Symbol}}()
        try
            for row in Topology.MoveGraphStore._rows(store.db,"SELECT payload FROM states WHERE area=3 AND complexity<=3")
                payload=JSON.parse(String(row.payload))
                state=OGraphState(Int.(payload["gauss"]),Int.(payload["vertex_signs"]))
                source=Topology._link_from_ograph_state(state)
                n=length(state.vertex_signs)
                for arity in (2,3)
                    arity>n && continue
                    for roles in Topology._permutations(collect(1:n))
                        for c in movable(source,"MP";vertices=roles[1:arity],strict=true)
                            key=(c.parameter,c.direction); key in covered && continue
                            r=verify_ograph_move(source,c)
                            r.valid || continue
                            @test replay_ograph_move(ograph_move_witness_data(r.witness)).valid
                            target=move(source,c)
                            inverse_roles=c.direction==:forward ? vcat(c.vertices,n+1) :
                                r.witness.source_to_target_vertices[c.vertices[1:2]]
                            reverse_candidates=movable(target,"MP";param=c.parameter,
                                vertices=inverse_roles,strict=true)
                            @test any(rc -> canonicalize_ograph(ograph_state(move(target,rc)))==
                                canonicalize_ograph(state),reverse_candidates)
                            push!(covered,key)
                        end
                    end
                end
            end
        finally
            close_store!(store)
        end
        println("MP type/direction coverage: ",sort(collect(covered)))
        @test length(covered)==32
    end
    @testset "Existing fixed store remains replayable" begin
        store=open_store(ARGS[1];migrate=false)
        try
            @test verify_ograph_store(store).transitions>=254
        finally
            close_store!(store)
        end
    end
end
