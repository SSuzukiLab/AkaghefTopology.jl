using Test, Topology

@testset "exact path weights after geometric selection" begin
    mktempdir() do directory
        store=open_store(joinpath(directory,"path.sqlite3"))
        try
            migrate_ograph_search_store!(store)
            source=VirtualLink()
            set_data!(source;rgauss=[[-1,-2,3,4,-4,-3,2,1]],orientation=ones(Int,4))
            set_weight!(source,[1,-1,-1,0,1,0,-1,0])
            sid=put_ograph_search_state!(store,source).state_id
            run=begin_ograph_search_run!(store,@__FILE__;label="affine relation smoke")
            outgoing=expand_ograph_search_state!(store,run,source,"02";
                selector=Dict("direction"=>"forward"))
            edge=only(Topology._ograph_search_path_rows(store,[String(first(outgoing).operation_id)]))
            relation=ograph_move_weight_relation(store,edge)
            cached=ograph_move_weight_relation(store,edge)
            @test cached==relation
            replay=replay_weighted_ograph_path(store,source,[edge])
            @test replay.valid
            @test relation.C*Topology._opath_flat(source)==relation.d
            @test relation.T*Topology._opath_flat(source)+relation.b==Topology._opath_flat(replay.target)
            # One nonzero H on the destination gives a concrete reachable
            # target; solve all coefficients and verify through public replay.
            desired=deepcopy(replay.target)
            move!(desired,"H";vertices=[1],param=3)
            solution=solve_ograph_path_weights(store,source,[edge];
                target_weights=Topology._opath_flat(desired))
            @test solution.solution.status==:integer_solution
            actual=record_weighted_ograph_path_trial!(store,run,source,[edge];h_moves=solution.h_moves)
            @test actual.replay.valid
            @test Topology._opath_flat(actual.replay.target)==Topology._opath_flat(desired)
            record_ograph_path_constraint_trial!(store,run,source,[edge],solution)

            # H incidence matrices agree with actual operations, including
            # self-adjacencies in this presentation.
            for v in (source,replay.target),vertex in 1:size(Topology._opath_H(v),2)
                h=deepcopy(v);move!(h,"H";vertices=[vertex],param=1)
                @test Topology._opath_flat(h)-Topology._opath_flat(v)==Topology._opath_H(v)[:,vertex]
            end

            # An actual geometric two-edge inverse loop must be queryable when
            # its endpoint is the EXACT source, not only a geometry-hash alias.
            inverse=expand_ograph_search_state!(store,run,replay.target,"02";
                selector=Dict("direction"=>"inverse"))
            @test any(r->r.target_state_id==sid,inverse)
            loops=query_ograph_search_paths(store,sid;target_state_id=sid,
                kind_pattern=["02","02"],limit=100)
            @test !isempty(loops)
            loop=first(loops)
            impossible=solve_ograph_path_weights(store,source,loop;
                target_weights=[1,-1,-1,2,1,0,-1,-2])
            @test impossible.solution.status!=:integer_solution
            iso=VirtualLink()
            code,signs=real_gauss_code(source)
            set_data!(iso;rgauss=[circshift(only(code),-3)],orientation=signs)
            iso_id=put_ograph_search_state!(store,iso).state_id
            permutation=[6,7,8,1,2,3,4,5]
            @test permutation in ograph_path_endpoint_permutations(iso,source)
            lifted=lift_ograph_search_path!(store,run,[edge],iso_id;
                source_permutation=permutation)
            @test only(lifted).source_state_id==iso_id
            manual=move(iso,"02";param="dir",vertices=Int.(edge.candidate["vertices"]))
            @test only(lifted).target_state_id==ograph_search_state(manual).state_id
            @test_throws ArgumentError lift_ograph_search_path!(store,run,[edge],iso_id;
                source_permutation=collect(1:8))
            finish_ograph_search_run!(store,run;status="completed-found")
            # A terminated writer's incomplete expansion can be resumed with
            # retained operations; a live competing run remains protected.
            Topology.MoveGraphStore._dbexec(store.db,
                "UPDATE ograph_search_expansions SET status='building' WHERE source_state_id=? AND move_kind='02'",
                (sid,))
            resume=begin_ograph_search_run!(store,@__FILE__;label="resume interrupted enumeration")
            resumed=expand_ograph_search_state!(store,resume,source,"02";
                selector=Dict("direction"=>"forward"))
            @test length(resumed)==length(outgoing)
            @test all(r->String(r.outcome)=="success",resumed)
            @test Int(Topology.MoveGraphStore._one(store.db,"PRAGMA busy_timeout").timeout)==10_000
            finish_ograph_search_run!(store,resume;status="completed-found")
        finally
            close_store!(store)
        end
    end
end
