using Test, Topology

include(joinpath(@__DIR__, "..", "src", "ograph_colour_moves.jl"))
include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))

using .OGraphColourMoves

const REGISTERED = (
    "s3_abalone",
    "lens_2_1",
    "lens_3_1_mst2_example23",
    "poincare",
    "quaternion",
    "s2xs1",
    "s2xs1_2vertex",
)

const KNOWN_H1 = Dict(
    "s3_abalone" => (0, BigInt[]),
    "lens_2_1" => (0, BigInt[2]),
    "lens_3_1_mst2_example23" => (0, BigInt[3]),
    "poincare" => (0, BigInt[]),
    "quaternion" => (0, BigInt[2, 2]),
    "s2xs1" => (1, BigInt[]),
    "s2xs1_2vertex" => (1, BigInt[]),
)

function show_sanity(label, report)
    println(label,
            " n=", report.n,
            " edges=", report.edges,
            " halfedges=", report.halfedges,
            " germs=", report.germs,
            " discs=", report.discs,
            " chi=", report.euler,
            " boundary_length=", report.boundary_length,
            " component_chi=", report.component_euler,
            " chi_boundary=", report.boundary_euler,
            " H1=(", report.h1_free, ", ", report.h1_torsion, ")",
            " dual_edges=", report.dual_edges,
            " dual_pairs=", report.dual_pairs,
            " calibration=", report.calibration,
            " ", report.status)
    report
end

function loop_graph()
    ColourMoveGraph(
        [(1, 2, 3, 4)],
        [3, 4, 1, 2],
        [(1, 3)],
        Dict((1, 3) => 0, (2, 4) => 0),
        OGRAPH_COLOUR_CONVENTION_ID,
    )
end

function four_edge_graph()
    ColourMoveGraph(
        [(1, 2, 3, 4), (5, 6, 7, 8)],
        [5, 6, 7, 8, 1, 2, 3, 4],
        [(1, 3), (1, 3)],
        Dict((1, 5) => 0, (2, 6) => 0, (3, 7) => 0, (4, 8) => 0),
        OGRAPH_COLOUR_CONVENTION_ID,
    )
end

# Figure 3.10 panel 1. Slots are counter-clockwise from west: west, south, east, north.
function figure_310()
    ColourMoveGraph(
        [(1, 4, 3, 2), (5, 8, 7, 6)],
        [7, 6, 5, 8, 3, 2, 1, 4],
        [(1, 3), (1, 3)],
        Dict((1, 7) => -1, (2, 6) => -1, (3, 5) => 0, (4, 8) => -1),
        OGRAPH_COLOUR_CONVENTION_ID,
    )
end

function compact_colour(ids, colour)
    order = sort(collect(ids))
    remap = Dict(old => new for (new, old) in enumerate(order))
    Dict(_pair_key(remap[left], remap[right]) => value for ((left, right), value) in colour)
end

_pair_key(left, right) = left < right ? (left, right) : (right, left)

function edge_roles(graph::ColourMoveGraph)
    owner = Dict{Int,Int}()
    for (vertex, cyclic) in enumerate(graph.cyclic)
        for halfedge in cyclic
            owner[halfedge] = vertex
        end
    end
    loops = Dict{Int,Vector{Int}}()
    spans = Dict{NTuple{2,Int},Int}()
    for ((left, right), value) in graph.colour
        a, b = owner[left], owner[right]
        if a == b
            push!(get!(loops, a, Int[]), value)
        else
            spans[_pair_key(a, b)] = value
        end
    end
    loops, spans
end

function c_with_counts(graph, budget)
    found = ColourMoveGraph[]
    function walk(current, left)
        if all(iszero, left)
            push!(found, current)
            return
        end
        for move in applicable_c(current)
            move.direction == :forward || continue
            left[move.vertex] <= 0 && continue
            left[move.vertex] -= 1
            walk(c_move(current, move), left)
            left[move.vertex] += 1
        end
    end
    walk(graph, copy(budget))
    found
end

function normal_h1(link)
    closed = certify_closed_ograph(ograph_state(link))
    ograph_h1(ograph_pi1_presentation(closed))
end

function h1_tuple(report)
    (report.h1_free, report.h1_torsion)
end

@testset "subsegment sum" begin
    @test subsegment_sum([1, 1, 1]) == 0
    @test subsegment_sum([0, -1, -1]) == 1
    @test subsegment_sum(Int[]) == 0
    @test_throws ArgumentError subsegment_sum(["x"])
    try
        subsegment_sum(["x"])
    catch error
        @test occursin("SUB", error.msg)
    end
end

@testset "C changes colours and reverses" begin
    graph = loop_graph()
    forward = filter(move -> move.direction == :forward, applicable_c(graph))
    inverse = filter(move -> move.direction == :inverse, applicable_c(graph))
    @test length(forward) == 2
    @test length(inverse) == 2
    @test all(move -> move.id == "C" && move.delta == 1, forward)
    @test all(move -> move.delta == -1, inverse)
    move = only(filter(item -> item.frame == 1, forward))
    moved = c_move(graph, move)
    @test moved.cyclic[1] == (1, 3, 2, 4)
    @test moved.over[1] == (2, 4)
    @test moved.colour[(1, 3)] == -1
    @test moved.colour[(2, 4)] == -1
    back = only(filter(item -> item.frame == 1 && item.direction == :inverse, applicable_c(moved)))
    @test c_move(moved, back) == graph
    @test_throws ArgumentError c_move(graph, CMove("C", 1, 2, :forward, 1, 4, 3))
    try
        c_move(graph, CMove("C", 1, 2, :forward, 1, 4, 3))
    catch error
        @test occursin("C", error.msg)
    end

    spread = four_edge_graph()
    spread_move = only(filter(item -> item.vertex == 1 && item.frame == 1, applicable_c(spread)))
    spread_moved = c_move(spread, spread_move)
    @test all(value -> value == 1, values(spread_moved.colour))
    back_spread = only(filter(item ->
        item.vertex == 1 && item.frame == 1 && item.direction == :inverse,
        applicable_c(spread_moved)))
    @test c_move(spread_moved, back_spread) == spread
end

@testset "MP on the printed ports and Figure 3.10" begin
    graph = figure_310()
    forward = filter(move -> move isa MPForward, applicable_mp(graph))
    @test length(forward) == 2
    @test isempty(filter(move -> move isa MPInverse, applicable_mp(graph)))
    printed = MPForward("MP", :forward, 1, 2, 1, 1, 0, -1, 1, 0)
    @test printed in forward
    turned = MPForward("MP", :forward, 2, 1, 3, 3, 0, -1, 1, 0)
    @test turned in forward
    moved = mp_move(graph, printed)
    @test length(moved.cyclic) == 3
    loops, spans = edge_roles(moved)
    @test sort(loops[1]) == [0]
    @test sort(loops[2]) == [-1]
    @test sort(loops[3]) == [1]
    @test spans[(1, 2)] == -1
    @test spans[(2, 3)] == 1
    @test spans[(1, 3)] == 0
    # Ports after compaction of surviving ids 1,2,4,6,7,8 and new ids 9:14.
    # Top (A north, new, new, B north), middle (A west, new, B east, new),
    # bottom (A south, B south, new, new). Over is vertical, horizontal, vertical.
    @test moved.cyclic == [(2, 7, 11, 4), (1, 9, 5, 8), (3, 6, 12, 10)]
    @test moved.over == [(2, 4), (1, 3), (2, 4)]
    # Surviving ids 1,2,4,6,7,8 and new ids 9:14, compacted in increasing order.
    # Upper (2,6) -> 0, inner (4,8) -> +1, outer (1,7) -> -1.
    expected_ids = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    raw = Dict(
        (2, 6) => 0, (4, 8) => 1, (1, 7) => -1,
        (9, 10) => -1, (11, 12) => 1, (13, 14) => 0,
    )
    @test moved.colour == compact_colour(expected_ids, raw)
    other = mp_move(graph, turned)
    other_raw = Dict(
        (2, 6) => 1, (4, 8) => 0, (1, 7) => -1,
        (9, 10) => -1, (11, 12) => 1, (13, 14) => 0,
    )
    @test other.colour == compact_colour(expected_ids, other_raw)
    # Panel 1 gives the three outer edges the same colour, so these two
    # witnesses are one coloured o-graph after vertex renumbering and slot
    # rotation. The dicts above are different before that gauge.
    @test coloured_ograph_hash(to_coloured(moved)) == coloured_ograph_hash(to_coloured(other))
    restored = filter(move -> move isa MPInverse, applicable_mp(moved))
    @test length(restored) == 1
    @test coloured_ograph_hash(to_coloured(mp_move(moved, only(restored)))) ==
          coloured_ograph_hash(to_coloured(graph))
    reached = c_with_counts(moved, [2, 1, 0])
    @test any(item -> begin
        item_loops, item_spans = edge_roles(item)
        sort(vcat(values(item_loops)...)) == [1, 1, 1] &&
            sort(collect(values(item_spans))) == [-1, -1, -1]
    end, reached)
    nonzero = ColourMoveGraph(graph.cyclic, graph.mate, graph.over,
                              merge(graph.colour, Dict((3, 5) => 1)), graph.convention_id)
    @test isempty(filter(move -> move isa MPForward, applicable_mp(nonzero)))
    @test_throws ArgumentError mp_move(nonzero, printed)
    bad = MPForward("MP", :forward, 1, 2, 1, 1, 1, -1, 1, 0)
    @test_throws ArgumentError mp_move(graph, bad)
    try
        mp_move(graph, bad)
    catch error
        @test occursin("MP", error.msg)
    end
    looped = MPForward("MP", :forward, 1, 1, 1, 1, 0, -1, 1, 0)
    @test_throws ArgumentError mp_move(graph, looped)
end

@testset "normal MP agreement after Proposition 3.2.1 recolouring" begin
    # Direct comparison on every registered example. One C before MP, and one C
    # after MP, only when the source has at most three vertices. calibration B.
    for name in REGISTERED
        link = vl_example(name)
        n = length(ograph_state(link).vertex_signs)
        source_graph = from_coloured(normal_coloured_ograph(link).graph)
        direct = Set{String}()
        for move in applicable_mp(source_graph)
            move isa MPForward || continue
            push!(direct, coloured_ograph_hash(to_coloured(mp_move(source_graph, move))))
        end
        conjugated = Set{String}()
        if n <= 3
            for cmove in applicable_c(source_graph)
                cmove.direction == :forward || continue
                shifted = c_move(source_graph, cmove)
                for move in applicable_mp(shifted)
                    move isa MPForward || continue
                    push!(conjugated, coloured_ograph_hash(to_coloured(mp_move(shifted, move))))
                end
            end
            for move in applicable_mp(source_graph)
                move isa MPForward || continue
                shifted = mp_move(source_graph, move)
                for cmove in applicable_c(shifted)
                    cmove.direction == :forward || continue
                    push!(conjugated, coloured_ograph_hash(to_coloured(c_move(shifted, cmove))))
                end
            end
        end
        accepted = 0
        direct_hits = 0
        conjugated_hits = 0
        unmatched = 0
        for a in 1:n, b in 1:n
            a == b && continue
            for candidate in movable(link, "MP"; vertices=[a, b], strict=true)
                candidate.direction == :forward || continue
                report = verify_ograph_move(link, candidate)
                report.valid || continue
                accepted += 1
                target = normal_coloured_ograph(move(link, candidate)).graph
                target_hash = coloured_ograph_hash(target)
                if target_hash in direct
                    direct_hits += 1
                elseif target_hash in conjugated
                    conjugated_hits += 1
                else
                    unmatched += 1
                end
            end
        end
        println(name, " accepted_forward_MP=", accepted,
                " direct=", direct_hits,
                " one_C=", n <= 3 ? conjugated_hits : "not_searched",
                " unmatched=", n <= 3 ? unmatched : "not_classified")
        # Locked to the run that printed these counts. A change is a finding.
        expected = Dict(
            "s3_abalone" => (accepted=0, direct=0),
            "lens_2_1" => (accepted=2, direct=2),
            "lens_3_1_mst2_example23" => (accepted=2, direct=0),
            "poincare" => (accepted=10, direct=0),
            "quaternion" => (accepted=6, direct=0),
            "s2xs1" => (accepted=6, direct=2),
            "s2xs1_2vertex" => (accepted=4, direct=2),
        )[name]
        @test (accepted=accepted, direct=direct_hits) == expected
        if n <= 3
            one_c_expected = Dict(
                "s3_abalone" => (one_C=0, unmatched=0),
                "lens_2_1" => (one_C=0, unmatched=0),
                "lens_3_1_mst2_example23" => (one_C=0, unmatched=2),
                "quaternion" => (one_C=0, unmatched=6),
                "s2xs1_2vertex" => (one_C=0, unmatched=2),
            )[name]
            @test (one_C=conjugated_hits, unmatched=unmatched) == one_c_expected
        end
    end
end

@testset "N1 rejects a zero-vertex side" begin
    @test OGRAPH_MOVE_VERIFIER_VERSION == "ograph-move-source-replay/3"
    empty_state = OGraphState(Int[], Int[])
    empty_witness = closedness_witness(empty_state)
    @test empty_witness.n1 == false
    @test empty_witness.criteria == (true, true, true)
    empty_link = VirtualLink()
    set_data!(empty_link; rgauss=[Int[]], orientation=Int[])
    @test is_closed(empty_link) == (false, (true, true, false))
    @test Topology._ograph_move_class(empty_witness, :source_outside_n1, :source_not_closed) ==
          :source_outside_n1
    closed_failing = closedness_witness(ograph_state(vl_example("s3_abalone")))
    @test Topology._ograph_move_class(closed_failing, :source_outside_n1, :source_not_closed) == :ok

    source = VirtualLink()
    set_data!(source; rgauss=[[1, 2, -3, -4, 4, 3, -2, -1]], orientation=[-1, -1, -1, -1])
    @test is_closed(source) == (true, (true, true, true))
    candidate = only(movable(source, "S"; vertices=[1, 2, 3, 4], strict=true))
    @test candidate.direction == :inverse
    moved = move(source, candidate)
    @test length(ograph_state(moved).vertex_signs) == 0
    report = verify_ograph_move(source, candidate)
    @test !report.valid
    @test report.reason == :target_outside_n1

    one = VirtualLink()
    set_data!(one; rgauss=[[-1, 1]], orientation=[-1])
    forward = only(movable(one, "PS"; param="1", edge_occurrences=[(1, 1), (1, 1)]))
    target = move(one, forward)
    inverse = only(movable(target, "PS"; param="1", vertices=[2, 3], strict=true))
    @test verify_ograph_move(target, inverse).valid
    @test length(ograph_state(move(target, inverse)).vertex_signs) == 1
end

function boundary_preserved(graph)
    poly = face_cycles(graph; calibration=:b)
    boundary = thickening_boundary(graph; calibration=:b)
    pairs = dual_face_pairings(graph; calibration=:b)
    homology = cellular_homology(graph; calibration=:b)
    orientable = all(component ->
        component.orientable && component.closed_surface && component.euler <= 2,
        boundary.components)
    reversing = all(pairing -> pairing.orientation_reversing, pairs)
    preserved = orientable && reversing && boundary.euler_sum == 2 * poly.euler &&
                length(pairs) == 2 * poly.vertex_count
    preserved, poly, boundary, homology
end

function assert_same_manifold(graph, before; vertex_delta::Int)
    preserved, poly, boundary, homology = boundary_preserved(graph)
    @test preserved
    @test poly.vertex_count == before.n + vertex_delta
    @test poly.edge_count == before.edges + 2 * vertex_delta
    @test poly.disc_count == before.discs + vertex_delta
    @test poly.euler == before.euler
    @test boundary.euler_sum == 2
    @test [component.euler for component in boundary.components] == [2]
    @test (homology.free_rank, homology.torsion_factors) ==
          (before.h1_free, before.h1_torsion)
    poly, boundary, homology
end

@testset "sanity and H1 on the registered closed examples" begin
    # calibration B. Both forward C seatings and every forward MP are applied.
    # The stored graph and its canonical representative must satisfy the same
    # counts. An earlier inequality was a half-edge numbering defect in
    # thickening_boundary, since fixed.
    @test OGRAPH_FACE_CYCLE_STATUS == "prop321-calibration-b/draft"
    for name in REGISTERED
        link = vl_example(name)
        colouring = normal_coloured_ograph(link)
        before = show_sanity("before C $name", move_sanity(colouring.graph; calibration=:b))
        @test before.discs == before.n + 1
        @test before.euler == 1
        @test before.component_euler == [2]
        @test before.boundary_euler == 2
        pipeline = normal_h1(link)
        @test (pipeline.free_rank, collect(pipeline.torsion_invariant_factors)) == KNOWN_H1[name]
        @test h1_tuple(before) == KNOWN_H1[name]
        graph = from_coloured(colouring.graph)
        for move in applicable_c(graph)
            move.direction == :forward || continue
            after_graph = to_coloured(c_move(graph, move))
            poly, boundary, homology = assert_same_manifold(after_graph, before; vertex_delta=0)
            canonical = canonical_coloured_ograph(after_graph)
            _, cboundary, _ = assert_same_manifold(canonical, before; vertex_delta=0)
            @test cboundary.euler_sum == boundary.euler_sum
            println("after C $name v=", move.vertex, " frame=", move.frame,
                    " n=", poly.vertex_count, " discs=", poly.disc_count,
                    " chi=", poly.euler,
                    " component_chi=", [component.euler for component in boundary.components],
                    " chi_boundary=", boundary.euler_sum,
                    " canonical_chi_boundary=", cboundary.euler_sum,
                    " H1=", (homology.free_rank, homology.torsion_factors))
        end
        forward_mp = filter(item -> item isa MPForward, applicable_mp(graph))
        println(name, " coloured forward MP sites=", length(forward_mp))
        for move in forward_mp
            mp_graph = to_coloured(mp_move(graph, move))
            poly, boundary, homology = assert_same_manifold(mp_graph, before; vertex_delta=1)
            canonical = canonical_coloured_ograph(mp_graph)
            _, cboundary, _ = assert_same_manifold(canonical, before; vertex_delta=1)
            @test cboundary.euler_sum == boundary.euler_sum
            println("after MP $name A=", move.vertex_a, " frameA=", move.frame_a,
                    " B=", move.vertex_b, " frameB=", move.frame_b,
                    " n=", poly.vertex_count, " discs=", poly.disc_count,
                    " chi=", poly.euler,
                    " component_chi=", [component.euler for component in boundary.components],
                    " chi_boundary=", boundary.euler_sum,
                    " canonical_chi_boundary=", cboundary.euler_sum,
                    " H1=", (homology.free_rank, homology.torsion_factors))
        end
    end
    sample = from_coloured(normal_coloured_ograph(vl_example("lens_2_1")).graph)
    swapped = ColourMoveGraph(reverse(sample.cyclic), sample.mate, reverse(sample.over),
                              sample.colour, sample.convention_id)
    @test coloured_ograph_hash(to_coloured(swapped)) == coloured_ograph_hash(to_coloured(sample))
    @test h1_tuple(move_sanity(to_coloured(swapped))) == h1_tuple(move_sanity(to_coloured(sample)))
end
