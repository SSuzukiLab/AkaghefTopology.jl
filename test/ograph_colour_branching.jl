using Test
using Topology

include(joinpath(@__DIR__, "..", "src", "ograph_colour_branching.jl"))
include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))

const BRANCHING_ROUND_TRIP_NAMES = (
    "s3_abalone",
    "lens_2_1",
    "lens_3_1_mst2_example23",
    "lens_3_2",
    "poincare",
    "quaternion",
    "s2xs1",
    "s2xs1_2vertex",
    "odata_variant",
)

function _presentation(name)
    ograph_state(vl_example(name))
end

function _canonical_set(states)
    Set(canonicalize_ograph(state) for state in states)
end

@testset "P2 colour of one step" begin
    @test p2_colour(true, true) == 0
    @test p2_colour(false, false) == 0
    @test p2_colour(false, true) == 1
    @test p2_colour(true, false) == -1
end

@testset "round trip through Prop. 3.2.1" begin
    for name in BRANCHING_ROUND_TRIP_NAMES
        state = _presentation(name)
        graph = coloured_ograph_from_normal(state)
        report = ograph_conditions(graph)
        @test conditions_passed(report)
        recovered = recovered_normal_states(graph; convention_id=state.convention_id)
        canonical = _canonical_set(recovered)
        @test canonical == Set([canonicalize_ograph(state)])
        reversed = reversed_normal_presentation(state)
        negated = negate_colours(graph)
        @test _canonical_set(recovered_normal_states(
            negated; convention_id=state.convention_id)) ==
              Set([canonicalize_ograph(reversed)])
        branching = compatible_branching(graph; convention_id=state.convention_id)
        @test branching.discs == :unread_fig25
        @test branching.provisional
        @test _canonical_set(branching.states) == canonical
        @test length(branching.diagrams) == length(branching.states)
    end
end

@testset "registered closed inputs and the two-vertex negative index" begin
    certified = String[]
    refused = String[]
    for input in vl_virtual_gauss_inputs()
        link = vl_example(input.name)
        try
            certify_closed_ograph(link)
            push!(certified, input.name)
        catch error
            error isa ArgumentError || rethrow()
            push!(refused, input.name)
        end
    end
    @test Set(certified) == Set(BRANCHING_ROUND_TRIP_NAMES)
    @test refused == ["koda_fig13"]
    negative = _presentation("s2xs1_2vertex")
    @test -1 in negative.vertex_signs
end

@testset "negated colours and a one-edge shift" begin
    lens = coloured_ograph_from_normal(_presentation("lens_2_1"))
    @test admits_compatible_branching(negate_colours(lens))
    shifted = shift_colour(lens, 2)
    @test conditions_passed(ograph_conditions(shifted))
    @test !admits_compatible_branching(shifted)
    @test any(!=(0), lens.colours)
end

@testset "colour tuples that do not come from a branching" begin
    for name in ("s3_abalone", "lens_2_1", "s2xs1_2vertex")
        graph = coloured_ograph_from_normal(_presentation(name))
        census = branching_colour_census(graph)
        nedges = length(graph.colours)
        @test census.total == 3^nedges
        @test census.branching == 2
        @test census.nonbranching == census.total - 2
        report = ograph_conditions(graph)
        @test report.O0 && report.O1 && report.O2 && report.O3
    end
    abalone = coloured_ograph_from_normal(_presentation("s3_abalone"))
    failures = 0
    values = (-1, 0, 1)
    for index in CartesianIndices((3, 3))
        colours = [values[index[1]], values[index[2]]]
        coloured = BranchColourGraph(abalone.nvertices, abalone.over_slots,
                                      abalone.edge_darts, colours)
        @test conditions_passed(ograph_conditions(coloured))
        admits_compatible_branching(coloured) || (failures += 1)
    end
    @test failures == 7
end

@testset "O0–O3 failures name the condition" begin
    empty = BranchColourGraph(0, NTuple{2,Int}[], NTuple{2,Int}[], Int[])
    empty_report = ograph_conditions(empty)
    @test !empty_report.O1
    @test any(text -> startswith(text, "O1:"), empty_report.failures)

    residue = coloured_ograph_from_normal(OGraphState([1, -1], [1]))
    residue.colours[1] = 2
    colour_report = ograph_conditions(residue)
    @test colour_report.O0 && colour_report.O1 && colour_report.O2 && colour_report.O3
    @test admits_compatible_branching(residue)
    short = BranchColourGraph(1, [(2, 4)], [(2, 1), (3, 4)], Int[0])
    short_report = ograph_conditions(short)
    @test !short_report.O3
    @test any(text -> startswith(text, "O3:"), short_report.failures)

    disconnected = BranchColourGraph(
        2,
        [(2, 4), (2, 4)],
        [(1, 2), (3, 4), (5, 6), (7, 8)],
        [0, 0, 0, 0],
    )
    disconnected_report = ograph_conditions(disconnected)
    @test !disconnected_report.O0
    @test any(text -> occursin("not connected", text), disconnected_report.failures)

    adjacent = coloured_ograph_from_normal(OGraphState([1, -1], [1]))
    adjacent_over = BranchColourGraph(adjacent.nvertices, [(1, 2)],
                                       adjacent.edge_darts, adjacent.colours)
    adjacent_report = ograph_conditions(adjacent_over)
    @test !adjacent_report.O2
    @test any(text -> startswith(text, "O2:"), adjacent_report.failures)

    missing_port = BranchColourGraph(1, [(2, 4)], [(1, 2), (3, 3)], [0, 0])
    missing_report = ograph_conditions(missing_port)
    @test !missing_report.O0
    @test any(text -> startswith(text, "O0:"), missing_report.failures)
end

function _core_from_branch(graph::BranchColourGraph)
    nvertices = graph.nvertices
    vertices = Tuple(Tuple((vertex - 1) * 4 + 1:(vertex - 1) * 4 + 4)
                     for vertex in 1:nvertices)
    mate = zeros(Int, 4 * nvertices)
    for (tail, head) in graph.edge_darts
        mate[tail] = head
        mate[head] = tail
    end
    overs = Tuple(graph.over_slots[vertex] for vertex in 1:nvertices)
    records = Tuple(EdgeColour(tail, colour)
                    for ((tail, _), colour) in zip(graph.edge_darts, graph.colours))
    ColouredOGraph(vertices, mate, overs, records)
end

@testset "adapter to the shared coloured o-graph" begin
    @test isdefined(Topology, :ColouredOGraph)
    @test_throws ArgumentError coloured_ograph("not a coloured o-graph")
    state = _presentation("quaternion")
    branch = coloured_ograph_from_normal(state)
    adapted = coloured_ograph(_core_from_branch(branch))
    @test _canonical_set(recovered_normal_states(adapted; convention_id=state.convention_id)) ==
          Set([canonicalize_ograph(state)])
    invalid = _core_from_branch(branch)
    adjacent = ColouredOGraph(
        invalid.vertex_halfedges,
        invalid.mate,
        Tuple(ntuple(_ -> (1, 2), length(invalid.vertex_halfedges))),
        invalid.edge_colours,
    )
    @test_throws ArgumentError coloured_ograph_from_core(adjacent)
end

@testset "one circuit cannot be all colour zero" begin
    lens = coloured_ograph_from_normal(_presentation("lens_2_1"))
    blank = BranchColourGraph(lens.nvertices, lens.over_slots, lens.edge_darts,
                              fill(0, length(lens.colours)))
    @test !admits_compatible_branching(blank)
    @test sum(colour_mod3(colour) for colour in lens.colours) == 0
end

@testset "two homogeneous circuits" begin
    # Under ports 1↔3 and over ports 2↔4 are separate straight-through circuits.
    graph = BranchColourGraph(
        2,
        [(2, 4), (2, 4)],
        [(1, 5), (3, 7), (2, 6), (4, 8)],
        [0, 0, 0, 0],
    )
    @test conditions_passed(ograph_conditions(graph))
    @test length(compatible_orientations(graph)) == 4
    @test isempty(recovered_normal_states(graph))
    diagrams = recovered_normal_diagrams(graph)
    @test length(diagrams) == 4
    @test all(diagram -> length(diagram.components) == 2, diagrams)
    @test length(unique(diagram.vertex_indices for diagram in diagrams)) > 1
    @test !admits_compatible_branching(shift_colour(graph, 1))
    census = branching_colour_census(graph)
    @test census.total == 81
    @test census.branching == 1
    branching = compatible_branching(graph)
    @test isempty(branching.states)
    @test length(branching.diagrams) == 4
end

@testset "two mixed circuits have four branching colourings" begin
    graph = BranchColourGraph(
        2,
        [(2, 4), (2, 4)],
        [(2, 5), (7, 4), (1, 6), (8, 3)],
        [-1, 1, 1, -1],
    )
    @test conditions_passed(ograph_conditions(graph))
    @test length(compatible_orientations(graph)) == 1
    @test admits_compatible_branching(negate_colours(graph))
    @test !admits_compatible_branching(shift_colour(graph, 1))
    census = branching_colour_census(graph)
    @test census.total == 81
    @test census.branching == 4
end

@testset "Fig. 2.5 port chart stays unread" begin
    @test FIG25_PORT_ASSIGNMENT == :unread
    @test PRO_REVIEW_MOMENT1_STATUS == :read
    @test PRO_REVIEW_MOMENT1 ==
          "threads/PJ261120/sources/pro/branching-20261122-000145-answer.md"
    @test PRO_REVIEW_MOMENT3_STATUS == :read
    @test PRO_REVIEW_MOMENT3 ==
          "threads/PJ261120/sources/pro/branching-20261122-013056.md"
end
