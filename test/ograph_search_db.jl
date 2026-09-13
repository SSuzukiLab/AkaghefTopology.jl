using Test
using Topology

function _search_s2xs1(weight=nothing)
    link = VirtualLink()
    set_data!(link;
        gauss=[[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],
        orientation=[1,1,1,1,0,0])
    weight === nothing || set_weight!(link, weight)
    link
end

@testset "unweighted o-graph search DB" begin
    mktempdir() do directory
        store = open_store(joinpath(directory, "ograph-search.sqlite3"))
        try
            migrate_ograph_search_store!(store)
            source = _search_s2xs1()
            weighted_source = _search_s2xs1([1,-1,-1,0,1,0,-1,0])
            shifted_source = _search_s2xs1([1,-1,-1,2,1,0,-1,-2])

            # Search identity is geometry-only by design. Weight vectors are
            # admitted only by the later replay function.
            source_state = put_ograph_search_state!(store, source)
            @test ograph_search_state(weighted_source).state_id == source_state.state_id
            @test ograph_search_state(shifted_source).state_id == source_state.state_id
            @test !occursin("weights", source_state.serialization)
            @test load_ograph_search_state(store, source_state.state_id).is_weighted == false

            run = begin_ograph_search_run!(store, @__FILE__;
                label="literal-02-smoke",
                scope=Dict("moves" => ["02"], "weights" => "excluded"))
            edges = expand_ograph_search_state!(store, run, source, "02";
                selector=Dict("direction" => "forward"))
            @test !isempty(edges)
            @test all(String(edge.move_kind) == "02" for edge in edges)
            @test all(String(edge.outcome) == "success" for edge in edges)

            target_id = String(first(edges).target_state_id)
            path = query_ograph_search_path(store, source_state.state_id;
                target_state_id=target_id, kinds=("02",), max_depth=1)
            @test path !== nothing
            @test length(path) == 1
            @test only(path).move_kind == "02"
            @test query_ograph_search_paths(store, source_state.state_id;
                target_state_id=target_id, kinds=("02",), max_depth=1,
                limit=10) == [path]

            replay = replay_weighted_ograph_path(store, weighted_source, path)
            @test replay.valid
            @test ograph_search_state(replay.target).state_id == target_id
            trial = record_weighted_ograph_path_trial!(store, run,
                weighted_source, path)
            @test trial.replay.valid
            @test verify_ograph_search_store(store).weighted_trials == 1
            @test verify_ograph_search_store(store).states == 1 + length(edges)

            # The same expansion is reused by a second script run instead of
            # re-enumerating the state neighborhood.
            second_run = begin_ograph_search_run!(store, @__FILE__;
                label="literal-02-reuse", scope=Dict("moves" => ["02"]))
            cached = expand_ograph_search_state!(store, second_run, source, "02";
                selector=Dict("direction" => "forward"))
            @test length(cached) == length(edges)
            second_summary = finish_ograph_search_run!(store, second_run;
                status="completed-none", counts=Dict("reused" => true))
            @test second_summary.expansions == 1

            summary = finish_ograph_search_run!(store, run;
                status="completed-found", counts=Dict("paths" => 1))
            @test summary.candidates == length(edges)
            @test summary.successes == length(edges)

            @test_throws ArgumentError ograph_search_candidates(source, "H")
            @test_throws ArgumentError expand_ograph_search_state!(
                store, run, source, "H")

            # Literal 02 and PS have disjoint move kinds and selectors.
            ps_candidates = ograph_search_candidates(source, "PS";
                selector=Dict("parameters" => ["1"], "direction" => "inverse"))
            @test all(candidate.type == "PS" for candidate in ps_candidates)
            @test all(candidate.type == "02" for candidate in
                ograph_search_candidates(source, "02";
                    selector=Dict("direction" => "forward")))

            verified = verify_ograph_search_store(store)
            @test verified.states >= 2
            @test verified.runs == 2
            @test verified.expansions == 1
            @test verified.operations == length(edges)
        finally
            close_store!(store)
        end
    end
end
