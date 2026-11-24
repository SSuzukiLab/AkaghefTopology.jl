using Test
using Topology
using Random

isdefined(@__MODULE__, :OGraphColourBoundary) ||
    include(joinpath(@__DIR__, "..", "src", "ograph_colour_boundary.jl"))

# cf:[[BP97#^bp97-sec-0549]] Fig. 2.5 — twelve physical pairings of the locked right/left chart; unverified

const _GOLDEN_PAIRINGS = (
    (true, true, 0, (((1, 1), (1, 1)), ((0, 1), (0, 1)), ((-1, -1), (-1, -1)))),
    (true, true, 1, (((1, 1), (0, 1)), ((0, 1), (-1, -1)), ((-1, -1), (1, 1)))),
    (true, true, -1, (((1, 1), (-1, -1)), ((0, 1), (1, 1)), ((-1, -1), (0, 1)))),
    (true, false, 0, (((-1, -1), (1, 1)), ((0, -1), (0, 1)), ((1, 1), (-1, -1)))),
    (true, false, 1, (((-1, -1), (0, 1)), ((0, -1), (-1, -1)), ((1, 1), (1, 1)))),
    (true, false, -1, (((-1, -1), (-1, -1)), ((0, -1), (1, 1)), ((1, 1), (0, 1)))),
    (false, true, 0, (((1, 1), (-1, -1)), ((0, 1), (0, -1)), ((-1, -1), (1, 1)))),
    (false, true, 1, (((1, 1), (0, -1)), ((0, 1), (1, 1)), ((-1, -1), (-1, -1)))),
    (false, true, -1, (((1, 1), (1, 1)), ((0, 1), (-1, -1)), ((-1, -1), (0, -1)))),
    (false, false, 0, (((-1, -1), (-1, -1)), ((0, -1), (0, -1)), ((1, 1), (1, 1)))),
    (false, false, 1, (((-1, -1), (0, -1)), ((0, -1), (1, 1)), ((1, 1), (-1, -1)))),
    (false, false, -1, (((-1, -1), (1, 1)), ((0, -1), (-1, -1)), ((1, 1), (0, -1)))),
)

function _pairing_key(pairs)
    sort!(collect(pairs))
end

function _read_gauged_colour(left_up, right_up, alpha, right_shift, left_shift)
    pairs = fig25_leg_pairing(left_up, right_up, alpha)
    right_number = Dict(fig25_leg_direction(right_up, leg) => fig25_tau(right_shift, leg) for leg in 1:3)
    left_number = Dict(fig25_leg_direction(left_up, leg) => fig25_tau(left_shift, leg) for leg in 1:3)
    image = zeros(Int, 3)
    for (right_direction, left_direction) in pairs
        image[right_number[right_direction]] = left_number[left_direction]
    end
    colour_from_gluing((image[1], image[2], image[3]))
end

function _checks(graph)
    Dict(check.id => check for check in validate_coloured_ograph(graph))
end

function _one_vertex(; over=(1, 3), colours=(0, 1), rights=(1, 3))
    ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), (over,),
        ((rights[1], colours[1]), (rights[2], colours[2])))
end

@testset "Fig. 2.5 locked leg directions" begin
    @test fig25_leg_direction(true, 1) == (1, 1)
    @test fig25_leg_direction(true, 2) == (0, 1)
    @test fig25_leg_direction(true, 3) == (-1, -1)
    @test fig25_leg_direction(false, 1) == (-1, -1)
    @test fig25_leg_direction(false, 2) == (0, -1)
    @test fig25_leg_direction(false, 3) == (1, 1)
    @test OGRAPH_COLOUR_CONVENTION_ID == "bp97-fig2.5-rl-tau/draft1"
end

@testset "Fig. 2.5 golden pairings" begin
    @test length(_GOLDEN_PAIRINGS) == 12
    for (left_up, right_up, alpha, expected) in _GOLDEN_PAIRINGS
        @test fig25_leg_pairing(left_up, right_up, alpha) == expected
        @test colour_from_gluing(edge_gluing_from_colour(alpha)) == alpha
    end
    @test edge_gluing_from_colour(0) == (1, 2, 3)
    @test edge_gluing_from_colour(1) == (2, 3, 1)
    @test edge_gluing_from_colour(-1) == (3, 1, 2)
    @test edge_gluing_from_colour(2) == edge_gluing_from_colour(-1)
    @test edge_gluing_from_colour(-2) == edge_gluing_from_colour(1)
    @test edge_gluing_from_colour(4) == edge_gluing_from_colour(1)
    @test fig25_leg_pairing(true, true, 2) == fig25_leg_pairing(true, true, -1)
    @test_throws ArgumentError colour_from_gluing((1, 3, 2))
    @test_throws ArgumentError DiscGerm(1, 4)
    @test_throws ArgumentError fig25_frame_reversal(4)
end

@testset "Fig. 2.5 gauge and frame reversal" begin
    for (left_up, right_up, alpha, _) in _GOLDEN_PAIRINGS
        for right_shift in 0:2, left_shift in 0:2
            @test _read_gauged_colour(left_up, right_up, alpha, right_shift, left_shift) ==
                  fig25_gauged_colour(alpha, right_shift, left_shift)
        end
    end
    for alpha in (-1, 0, 1), beta in (-1, 0, 1)
        composed = ntuple(leg -> edge_gluing_from_colour(beta)[edge_gluing_from_colour(alpha)[leg]], 3)
        @test composed == edge_gluing_from_colour(fig25_gauged_colour(alpha, 0, beta))
    end
    for alpha in -1:1, leg in 1:3
        @test fig25_frame_reversal(fig25_tau(-alpha, fig25_frame_reversal(leg))) ==
              fig25_tau(alpha, leg)
    end
    mate = (2, 1, 4, 3)
    for alpha in (-1, 0, 1)
        forward = edge_disc_germs(alpha, 4, mate)
        opposite = edge_disc_germs(alpha, 3, mate)
        transported = (
            (DiscGerm(3, fig25_frame_reversal(left_germ.leg)),
             DiscGerm(4, fig25_frame_reversal(right_germ.leg)))
            for (right_germ, left_germ) in forward)
        @test Set(opposite) == Set(transported)
    end
    # Freezing the printed numbers and only reversing the bijection yields -α.
    # That algebra is not the stored colour: the stored transport is ρ with α unchanged.
    for (left_up, right_up, alpha, _) in _GOLDEN_PAIRINGS
        frozen = fig25_leg_pairing(right_up, left_up, -alpha)
        flipped = ((left_direction, right_direction)
                   for (right_direction, left_direction) in
                   fig25_leg_pairing(left_up, right_up, alpha))
        @test _pairing_key(frozen) == _pairing_key(flipped)
        @test fig25_gauged_colour(alpha, 0, 0) == alpha
    end
end

@testset "disc germs of one drawing" begin
    germs = edge_disc_germs(1, 4, (2, 1, 4, 3))
    @test germs == ((DiscGerm(4, 1), DiscGerm(3, 2)),
                    (DiscGerm(4, 2), DiscGerm(3, 3)),
                    (DiscGerm(4, 3), DiscGerm(3, 1)))
    @test edge_disc_germs(2, 4, (2, 1, 4, 3)) == edge_disc_germs(-1, 4, (2, 1, 4, 3))
    @test_throws ArgumentError edge_disc_germs(1, 5, (2, 1, 4, 3))
end

@testset "O0-O3 validator" begin
    ok = _one_vertex()
    passed = _checks(ok)
    for id in ("O0", "O1", "O2", "O3")
        @test passed[id].passed
        @test passed[id].id == id
        @test occursin(id, passed[id].message)
    end

    loops = _one_vertex(colours=(1, -1))
    @test all(check.passed for check in validate_coloured_ograph(loops))

    multi = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((1, 3), (2, 4)),
        ((1, 0), (2, 1), (3, -1), (4, 0)))
    @test all(check.passed for check in validate_coloured_ograph(multi))

    disconnected = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (2, 1, 4, 3, 6, 5, 8, 7),
        ((1, 3), (1, 3)),
        ((1, 0), (3, 1), (5, -1), (7, 0)))
    split = _checks(disconnected)
    @test !split["O0"].passed
    @test occursin("O0", split["O0"].message)
    @test occursin("connected", split["O0"].message)
    @test split["O1"].passed && split["O2"].passed && split["O3"].passed

    odd = ColouredOGraph(((1, 2, 3),), (2, 1, 3), ((1, 3),), ())
    degree = _checks(odd)
    @test !degree["O0"].passed
    @test occursin("O0", degree["O0"].message)
    @test occursin("quadrivalent", degree["O0"].message)

    missing = _checks(_one_vertex(over=()))
    @test missing["O0"].passed && missing["O1"].passed && missing["O3"].passed
    @test !missing["O2"].passed
    @test occursin("O2", missing["O2"].message)
    @test occursin("over pair", missing["O2"].message)

    adjacent = _checks(_one_vertex(over=(1, 2)))
    @test !adjacent["O2"].passed
    @test occursin("not opposite", adjacent["O2"].message)

    both = _checks(_one_vertex(over=(1, 2, 3, 4)))
    @test !both["O2"].passed
    @test occursin("more than one", both["O2"].message)

    congruent = ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((1, 3),),
        ((1, 2), (2, -1), (3, 0)))
    @test _checks(congruent)["O3"].passed

    two_colours = _checks(ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((1, 3),),
        ((1, 1), (2, 0), (3, 0))))
    @test !two_colours["O3"].passed
    @test occursin("more than one colour", two_colours["O3"].message)

    uncoloured = _checks(ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((1, 3),),
        ((1, 0),)))
    @test !uncoloured["O3"].passed
    @test occursin("no colour", uncoloured["O3"].message)

    written_two = _checks(_one_vertex(colours=(2, 0)))
    @test written_two["O0"].passed && written_two["O1"].passed
    @test written_two["O2"].passed && written_two["O3"].passed

    empty = _checks(ColouredOGraph((), (), (), ()))
    @test !empty["O0"].passed
    @test !empty["O1"].passed
    @test occursin("O1", empty["O1"].message)
    @test_throws ArgumentError ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((1, 3),),
        ((1, 0), (3, 0)); convention_id="  ")
    @test_throws ArgumentError ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), (),
        ((1, 0), (3, 0)))
end

@testset "canonical hash keeps the colour and does not quotient C or MP" begin
    original = _one_vertex(colours=(1, 0))
    rotated = ColouredOGraph(((2, 3, 4, 1),), (2, 1, 4, 3), ((2, 4),),
        ((1, 1), (3, 0)))
    other_end = _one_vertex(colours=(1, 0), rights=(2, 4))
    negated = _one_vertex(colours=(-1, 0), rights=(2, 4))
    @test coloured_ograph_hash(original) == coloured_ograph_hash(rotated)
    @test coloured_ograph_hash(original) == coloured_ograph_hash(other_end)
    @test coloured_ograph_hash(original) != coloured_ograph_hash(negated)
    @test canonical_coloured_ograph(original) == canonical_coloured_ograph(rotated)
    @test canonical_coloured_ograph(original) == canonical_coloured_ograph(other_end)
    @test canonical_coloured_ograph(canonical_coloured_ograph(original)) ==
          canonical_coloured_ograph(original)
    @test coloured_ograph_hash(original) != coloured_ograph_hash(_one_vertex(colours=(1, 1)))
    @test coloured_ograph_hash(_one_vertex(colours=(2, 0))) ==
          coloured_ograph_hash(_one_vertex(colours=(-1, 0)))
    reduced = canonical_coloured_ograph(_one_vertex(colours=(2, 4)))
    @test all(record.colour in (-1, 0, 1) for record in reduced.edge_colours)
    @test any(record.colour == -1 for record in reduced.edge_colours)
    @test any(record.colour == 1 for record in reduced.edge_colours)

    mirror = ColouredOGraph(((1, 4, 3, 2),), (2, 1, 4, 3), ((1, 3),),
        ((1, 1), (3, 0)))
    @test coloured_ograph_hash(original) != coloured_ograph_hash(mirror)

    adjacent_loops = ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((1, 3),),
        ((1, 1), (3, 0)))
    opposite_loops = ColouredOGraph(((1, 2, 3, 4),), (3, 4, 1, 2), ((1, 3),),
        ((1, 1), (2, 0)))
    @test coloured_ograph_hash(adjacent_loops) != coloured_ograph_hash(opposite_loops)

    parallel = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((1, 3), (2, 4)),
        ((1, 0), (2, 1), (3, -1), (4, 0)))
    other_ends = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((1, 3), (2, 4)),
        ((5, 0), (6, 1), (7, -1), (8, 0)))
    reordered = ColouredOGraph(
        ((5, 6, 7, 8), (1, 2, 3, 4)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((2, 4), (1, 3)),
        ((1, 0), (2, 1), (3, -1), (4, 0)))
    @test coloured_ograph_hash(parallel) == coloured_ograph_hash(other_ends)
    @test coloured_ograph_hash(parallel) == coloured_ograph_hash(reordered)
    changed = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((1, 3), (2, 4)),
        ((1, 1), (2, 1), (3, -1), (4, 0)))
    @test coloured_ograph_hash(parallel) != coloured_ograph_hash(changed)
    other_over = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((2, 4), (2, 4)),
        ((1, 0), (2, 1), (3, -1), (4, 0)))
    @test coloured_ograph_hash(parallel) != coloured_ograph_hash(other_over)
    @test_throws ArgumentError coloured_ograph_hash(
        ColouredOGraph(((1, 2, 3, 4), (5, 6, 7, 8)),
            (2, 1, 4, 3, 6, 5, 8, 7), ((1, 3), (1, 3)),
            ((1, 0), (3, 0), (5, 0), (7, 0))))
end

isdefined(@__MODULE__, :vl_example) || include(joinpath(@__DIR__, "../examples/VLExample.jl"))

const _CLOSED_COLOUR_NAMES = (
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

function _direction_set(edges)
    Set((edge.tail, edge.head) for edge in edges)
end

function _replace_colour(graph, halfedge, colour)
    records = map(graph.edge_colours) do record
        record.halfedge == halfedge ? (halfedge, colour) : (record.halfedge, record.colour)
    end
    ColouredOGraph(graph.vertex_halfedges, graph.mate, graph.over_positions, records;
                   convention_id=graph.convention_id)
end

@testset "Prop. 3.2.1 colour of one step" begin
    @test prop321_colour(true, true) == 0
    @test prop321_colour(false, false) == 0
    @test prop321_colour(false, true) == 1
    @test prop321_colour(true, false) == -1
    @test OGRAPH_FACE_CYCLE_STATUS == "prop321-calibration-b/draft"
end

@testset "closed normal colouring of the example registry" begin
    seen = Int[]
    signs = Int[]
    for name in _CLOSED_COLOUR_NAMES
        link = vl_example(name)
        state = ograph_state(link)
        colouring = normal_coloured_ograph(state)
        @test normal_coloured_ograph(link).colours == colouring.colours
        certified = normal_coloured_ograph(certify_closed_ograph(state))
        @test sort(collect(certified.colours)) == sort(collect(colouring.colours))
        @test length(certified.graph.vertex_halfedges) == length(colouring.graph.vertex_halfedges)
        @test all(check.passed for check in validate_coloured_ograph(colouring.graph))
        normal = validate_normal_colouring(colouring)
        @test map(check -> check.id, normal) == ("N1", "N2", "C1", "C2", "C3", "P1", "P2")
        @test all(check.passed for check in normal)
        @test all(occursin(check.id, check.message) for check in normal)
        expected = [prop321_colour(edge.crossing[1] > 0, edge.crossing[2] > 0)
                    for edge in real_edge_table(link)]
        @test collect(colouring.colours) == expected
        @test any(found -> _direction_set(found) == _direction_set(colouring.edges),
                  prop321_orientations(colouring.graph))
        append!(seen, colouring.colours)
        append!(signs, state.vertex_signs)
    end
    @test -1 in seen && 0 in seen && 1 in seen
    @test -1 in signs && 1 in signs

    lens = normal_coloured_ograph(ograph_state(vl_example("lens_2_1")))
    @test lens.colours == (-1, 0, 1, 0)
    @test coloured_gc_display(lens) == "coloured GC [[[-2,-1,+1,+2]];[+,+];[-1,0,1,0]]"
    abalone = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")))
    @test abalone.colours == (1, -1)
    @test coloured_gc_display(abalone) == "coloured GC [[[+1,-1]];[+];[1,-1]]"

    empty_error = try
        normal_coloured_ograph(OGraphState(Int[], Int[]))
        nothing
    catch error
        error
    end
    @test empty_error isa ArgumentError
    @test occursin("N1", empty_error.msg)

    koda_error = try
        normal_coloured_ograph(ograph_state(vl_example("koda_fig13")))
        nothing
    catch error
        error
    end
    @test koda_error isa ArgumentError
    @test occursin("C1-C3", koda_error.msg)
end

@testset "P2 counts on undirected colourings" begin
    several = ColouredOGraph(
        ((1, 2, 3, 4), (5, 6, 7, 8)),
        (5, 6, 7, 8, 1, 2, 3, 4),
        ((2, 4), (2, 4)),
        ((1, 0), (2, 0), (3, 0), (4, 0)))
    @test all(check.passed for check in validate_coloured_ograph(several))
    @test length(normal_orientations(several)) == 4
    @test length(prop321_orientations(several)) == 4

    none = _one_vertex(over=(2, 4), colours=(0, 0))
    @test length(normal_orientations(none)) == 2
    @test isempty(prop321_orientations(none))

    one = _one_vertex(over=(2, 4), colours=(1, -1))
    @test length(normal_orientations(one)) == 2
    @test length(prop321_orientations(one)) == 1

    lens = normal_coloured_ograph(ograph_state(vl_example("lens_2_1")))
    shifted = _replace_colour(lens.graph, lens.edges[1].tail, lens.colours[1] + 1)
    @test all(check.passed for check in validate_coloured_ograph(shifted))
    @test !any(found -> _direction_set(found) == _direction_set(lens.edges),
               prop321_orientations(shifted))
    @test edge_gluing_from_colour(lens.colours[1]) != edge_gluing_from_colour(lens.colours[1] + 1)
    @test !isempty(normal_orientations(shifted))
end

@testset "plot_svg draws Z3 colours" begin
    link = vl_example("lens_2_1")
    colouring = normal_coloured_ograph(ograph_state(link))
    path = joinpath(tempdir(), "ograph-colour-lens.svg")
    plot_svg(link, path; edge_labels=:colour, edge_colours=collect(colouring.colours))
    svg = read(path, String)
    rm(path)
    @test occursin(">-1</text>", svg)
    @test occursin(">0</text>", svg)
    @test occursin(">1</text>", svg)
    @test !occursin("paint-order", svg)
    @test occursin("fill=\"#faf9f7\"", svg)
    @test_throws ArgumentError plot_svg(link, path; edge_labels=:colour)
end

@testset "one false normal criterion at a time" begin
    lens = normal_coloured_ograph(ograph_state(vl_example("lens_2_1")))
    function replaced(; graph=lens.graph, edges=lens.edges, colours=lens.colours,
                      witness=lens.witness)
        ClosedNormalColouring(graph, edges, colours, witness, lens.state)
    end
    empty_checks = validate_normal_colouring(replaced(;
        graph=ColouredOGraph((), (), (), ()), edges=(), colours=()))
    @test map(check -> check.id, empty_checks) == ("N1", "N2", "C1", "C2", "C3", "P1", "P2")
    @test !empty_checks[1].passed
    @test occursin("at least one vertex", empty_checks[1].message)
    @test all(check.passed for check in empty_checks[2:end])

    short = ColouredOGraph(((1, 2, 3),), (2, 1, 3), ((1,),), ((1, 0),))
    short_checks = validate_normal_colouring(replaced(; graph=short, edges=(), colours=()))
    @test !short_checks[1].passed
    @test short_checks[1].id == "N1"
    @test occursin("not a normal crossing", short_checks[1].message)

    first = lens.edges[1]
    reversed = (NormalEdgeOrientation(first.head, first.tail), lens.edges[2:end]...)
    reversed_checks = validate_normal_colouring(replaced(; edges=reversed))
    @test reversed_checks[1].id == "N1" && reversed_checks[1].passed
    @test reversed_checks[2].id == "N2" && !reversed_checks[2].passed
    @test all(check.passed for check in reversed_checks[3:5])
    @test reversed_checks[6].id == "P1" && !reversed_checks[6].passed
    @test reversed_checks[7].id == "P2" && !reversed_checks[7].passed

    old = lens.witness
    for (index, id, criteria) in (
            (3, "C1", (false, true, true)),
            (4, "C2", (true, false, true)),
            (5, "C3", (true, true, false)))
        witness = ClosednessWitness(
            old.canonical_hash, old.verifier_version, criteria,
            old.c1_oriented_circuits, old.c2_trivalent_components,
            old.c3_method_circuits, old.true_vertices, old.method_circuit_lengths)
        checks = validate_normal_colouring(replaced(; witness))
        @test checks[index].id == id
        @test !checks[index].passed
        @test all(checks[j].passed for j in 1:7 if j != index)
    end
end

const _REFUTED_A_COUNTS = Dict(
    "s3_abalone" => (2, 1),
    "lens_2_1" => (3, 1),
    "lens_3_1_mst2_example23" => (3, 1),
    "s2xs1_2vertex" => (5, 3),
    "lens_3_2" => (2, -1),
    "quaternion" => (4, 1),
    "odata_variant" => (5, 1),
    "s2xs1" => (5, 1),
    "poincare" => (2, -3),
)

const _ADOPTED_B_LENGTHS = Dict(
    "s3_abalone" => [1, 5],
    "lens_2_1" => [1, 3, 8],
    "lens_3_1_mst2_example23" => [2, 5, 5],
    "s2xs1_2vertex" => [2, 4, 6],
    "lens_3_2" => [1, 3, 7, 7],
    "quaternion" => [3, 5, 5, 5],
    "odata_variant" => [1, 1, 2, 10, 10],
    "s2xs1" => [1, 1, 6, 8, 8],
    "poincare" => [5, 5, 5, 5, 5, 5],
)

const _NAMED_CLOSED_H1 = Dict(
    "s3_abalone" => (0, ()),
    "lens_2_1" => (0, (2,)),
    "lens_3_1_mst2_example23" => (0, (3,)),
    "lens_3_2" => (0, (3,)),
    "quaternion" => (0, (2, 2)),
    "s2xs1" => (1, ()),
    "s2xs1_2vertex" => (1, ()),
    "poincare" => (0, ()),
)

function _cycle_lengths(poly)
    sort!(collect(length(cycle.germs) for cycle in poly.cycles))
end

function _h1_pair(homology)
    (Int(homology.free_rank), Tuple(Int(factor) for factor in homology.torsion_factors))
end

function _corner_partition(poly)
    sorted_germ_partition((cycle.corners for cycle in poly.cycles))
end

function _rotate_origins(graph, rotation)
    cyclics = Vector{Int}[]
    overs = Vector{Int}[]
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        turn = rotation[vertex]
        rotated = vcat(collect(cyclic)[turn + 1:end], collect(cyclic)[1:turn])
        push!(cyclics, rotated)
        push!(overs, [mod1(slot - turn, 4) for slot in graph.over_positions[vertex]])
    end
    records = [(record.halfedge, record.colour) for record in graph.edge_colours]
    ColouredOGraph(Tuple(Tuple(cyclic) for cyclic in cyclics), graph.mate,
        Tuple(Tuple(positions) for positions in overs), records;
        convention_id=graph.convention_id)
end

function _swap_two_vertices(graph)
    n = length(graph.vertex_halfedges)
    new_of_old = Dict{Int,Int}()
    cyclics = Vector{Int}[]
    overs = Vector{Int}[]
    next_id = 1
    for old in (2, 1)
        rotated = Int[]
        for halfedge in graph.vertex_halfedges[old]
            new_of_old[halfedge] = next_id
            push!(rotated, next_id)
            next_id += 1
        end
        push!(cyclics, rotated)
        push!(overs, collect(Int, graph.over_positions[old]))
    end
    mate = zeros(Int, 4n)
    colour_of = Dict{NTuple{2,Int},Int}()
    for record in graph.edge_colours
        other = Int(graph.mate[record.halfedge])
        colour_of[minmax(record.halfedge, other)] = record.colour
    end
    for (ends, colour) in colour_of
        left = new_of_old[ends[1]]
        right = new_of_old[ends[2]]
        mate[left] = right
        mate[right] = left
    end
    records = [(min(new_of_old[ends[1]], new_of_old[ends[2]]), colour)
               for (ends, colour) in colour_of]
    ColouredOGraph(Tuple(Tuple(cyclic) for cyclic in cyclics), Tuple(mate),
        Tuple(Tuple(positions) for positions in overs), records;
        convention_id=graph.convention_id)
end

function _abalone_both_positive()
    colouring = normal_coloured_ograph(ograph_state(vl_example("s3_abalone")))
    graph = colouring.graph
    for record in colouring.graph.edge_colours
        record.colour == -1 || continue
        graph = _replace_colour(graph, record.halfedge, 1)
    end
    graph
end

function _invariant_key(graph)
    poly = face_cycles(graph)
    boundary = thickening_boundary(graph)
    homology = cellular_homology(graph)
    component_rows = [(component.euler, component.face_count, component.edge_count,
        component.vertex_count, component.orientable, component.closed_surface)
                       for component in boundary.components]
    (poly.disc_count, poly.euler, _cycle_lengths(poly), component_rows,
     boundary.euler_sum, _h1_pair(homology))
end

@testset "face cycles of the adopted calibration" begin
    failures = Tuple{Int,String,Int,Int,Int}[]
    println("EXAMPLE closed face cycles")
    for name in _CLOSED_COLOUR_NAMES
        link = vl_example(name)
        colouring = normal_coloured_ograph(ograph_state(link))
        graph = colouring.graph
        for calibration in (:a, :b)
            poly = face_cycles(graph; calibration)
            boundary = thickening_boundary(graph; calibration)
            pairings = dual_face_pairings(graph; calibration)
            orbits = dual_edge_orbits(graph; calibration)
            homology = cellular_homology(graph; calibration)
            n = poly.vertex_count
            lengths = _cycle_lengths(poly)
            corners = collect(Iterators.flatten(cycle.corners for cycle in poly.cycles))
            pipeline = ograph_h1(ograph_pi1_presentation(certify_closed_ograph(link)))
            pipe_h1 = (Int(pipeline.free_rank), Tuple(Int(factor)
                for factor in pipeline.torsion_invariant_factors))
            ours_h1 = _h1_pair(homology)
            println("EXAMPLE ", name, " ", calibration,
                " n=", n, " e=", poly.edge_count, " half=", poly.halfedge_count,
                " germs=", poly.germ_count, " d=", poly.disc_count, " chi=", poly.euler,
                " blen=", poly.boundary_length, " lengths=", lengths,
                " bcomps=", length(boundary.components),
                " beuler=", [component.euler for component in boundary.components],
                " bsum=", boundary.euler_sum, " h1=", ours_h1, " pipe=", pipe_h1)
            @test poly.edge_count == 2n
            @test poly.halfedge_count == 4n
            @test poly.germ_count == 6n
            @test poly.boundary_length == 6n
            @test poly.euler == poly.disc_count - n
            @test length(corners) == 6n
            @test length(unique(corners)) == 6n
            @test sum(lengths) == 6n
            @test length(pairings) == 2n
            @test all(pairing.orientation_reversing for pairing in pairings)
            @test length(orbits) == poly.disc_count
            @test sorted_germ_partition(orbits) == _corner_partition(poly)
            @test all(component.closed_surface && component.orientable &&
                component.euler <= 2 &&
                component.euler == component.vertex_count - component.edge_count + component.face_count
                for component in boundary.components)
            @test boundary.euler_sum == 2 * poly.euler
            if calibration == :a
                @test (poly.disc_count, poly.euler) == _REFUTED_A_COUNTS[name]
                @test !compare_disc_partition(link, colouring; calibration=:a).equal
                if poly.disc_count != n + 1 || poly.euler != 1
                    push!(failures, (n, name, poly.disc_count, poly.euler, poly.germ_count))
                end
            else
                @test poly.disc_count == n + 1
                @test poly.euler == 1
                @test lengths == _ADOPTED_B_LENGTHS[name]
                @test length(boundary.components) == 1
                component = only(boundary.components)
                @test component.euler == 2
                @test component.face_count == 4n
                @test component.edge_count == 6n
                @test component.vertex_count == 2n + 2
                comparison = compare_disc_partition(link, colouring)
                @test comparison.equal
                @test comparison.germ_count == 6n
                @test comparison.unique_germs == 6n
                @test ours_h1 == pipe_h1
                if haskey(_NAMED_CLOSED_H1, name)
                    @test ours_h1 == _NAMED_CLOSED_H1[name]
                else
                    @test name == "odata_variant"
                    @test ours_h1 == (0, ())
                end
            end
        end
    end
    sort!(failures)
    @test length(failures) == 3
    @test failures[1] == (2, "s2xs1_2vertex", 5, 3, 12)
    @test Set(item[2] for item in failures) == Set(("s2xs1_2vertex", "lens_3_2", "poincare"))

    println("EXAMPLE non-normal face cycles")
    zero_zero = _one_vertex(over=(2, 4), colours=(0, 0))
    one_one = _one_vertex(over=(2, 4), colours=(1, 1))
    plus = _abalone_both_positive()
    for (label, graph) in (("one_vertex_00", zero_zero), ("one_vertex_11", one_one),
                           ("s3_abalone_both_plus", plus))
        poly = face_cycles(graph)
        boundary = thickening_boundary(graph)
        pairings = dual_face_pairings(graph)
        orbits = dual_edge_orbits(graph)
        homology = cellular_homology(graph)
        n = poly.vertex_count
        lengths = _cycle_lengths(poly)
        println("EXAMPLE ", label, " b n=", n, " d=", poly.disc_count, " chi=", poly.euler,
            " lengths=", lengths, " bcomps=", length(boundary.components),
            " bsum=", boundary.euler_sum, " h1=", _h1_pair(homology))
        @test n == 1
        @test poly.edge_count == 2 && poly.halfedge_count == 4 && poly.germ_count == 6
        @test poly.disc_count == 2 && poly.euler == 1 && lengths == [2, 4]
        @test poly.boundary_length == 6
        @test length(boundary.components) == 1
        @test only(boundary.components).euler == 2
        @test boundary.euler_sum == 2
        @test only(boundary.components).orientable && only(boundary.components).closed_surface
        @test length(pairings) == 2 && all(pairing.orientation_reversing for pairing in pairings)
        @test length(orbits) == 2
        @test sorted_germ_partition(orbits) == _corner_partition(poly)
        @test _h1_pair(homology) == (0, (4,))
    end
    @test _corner_partition(face_cycles(zero_zero)) != _corner_partition(face_cycles(one_one))
    @test _corner_partition(face_cycles(one_one)) == _corner_partition(face_cycles(plus))
    plus_a = face_cycles(plus; calibration=:a)
    println("EXAMPLE s3_abalone_both_plus a d=", plus_a.disc_count, " chi=", plus_a.euler,
        " lengths=", _cycle_lengths(plus_a))
    @test plus_a.disc_count == 3 && plus_a.euler == 2
    @test _cycle_lengths(plus_a) == [1, 1, 4]

    println("EXAMPLE gauge")
    for name in ("s3_abalone", "lens_2_1", "s2xs1_2vertex", "lens_3_2", "poincare")
        graph = normal_coloured_ograph(ograph_state(vl_example(name))).graph
        base = _invariant_key(graph)
        n = length(graph.vertex_halfedges)
        for turn in 0:3
            @test _invariant_key(_rotate_origins(graph, fill(turn, n))) == base
        end
        n < 2 && continue
        mixed = _rotate_origins(graph, vcat([1, 0], fill(2, n - 2)))
        @test _invariant_key(mixed) == base
        println("EXAMPLE ", name, " gauge-invariant")
    end
    lens = normal_coloured_ograph(ograph_state(vl_example("lens_2_1"))).graph
    @test _invariant_key(_swap_two_vertices(lens)) == _invariant_key(lens)
    println("EXAMPLE lens_2_1 vertex-swap gauge-invariant")
end

function _apply_presentation_gauge(graph, σ, vertex_map, turns)
    n = length(graph.vertex_halfedges)
    cyclics = Vector{NTuple{4,Int}}(undef, n)
    overs = Vector{NTuple{2,Int}}(undef, n)
    for old in 1:n
        cyclic = collect(Int, graph.vertex_halfedges[old])
        turn = turns[old]
        rotated = vcat(cyclic[turn + 1:end], cyclic[1:turn])
        cyclics[vertex_map[old]] = ntuple(index -> σ[rotated[index]], 4)
        overs[vertex_map[old]] = ntuple(
            index -> mod1(Int(graph.over_positions[old][index]) - turn, 4), 2)
    end
    mate = zeros(Int, length(σ))
    for (halfedge, other) in enumerate(graph.mate)
        mate[σ[halfedge]] = σ[Int(other)]
    end
    records = [(σ[Int(record.halfedge)], Int(record.colour)) for record in graph.edge_colours]
    ColouredOGraph(Tuple(cyclics), Tuple(mate), Tuple(overs), records;
                   convention_id=graph.convention_id)
end

function _transport_corner(corner, vertex_map, turns)
    vertex, left, right = corner
    turn = turns[vertex]
    a = mod1(left - turn, 4)
    b = mod1(right - turn, 4)
    (vertex_map[vertex], min(a, b), max(a, b))
end

function _transport_partition(poly, vertex_map, turns)
    discs = Vector{NTuple{3,Int}}[]
    for cycle in poly.cycles
        push!(discs, NTuple{3,Int}[_transport_corner(corner, vertex_map, turns)
                                    for corner in cycle.corners])
    end
    sorted_germ_partition(discs)
end

function _ordered_pair(left, right, reversing)
    if right < left
        left, right = right, left
    end
    (left, right, reversing)
end

function _transported_flags(pairings, vertex_map, turns)
    keys = Tuple{NTuple{2,Int},NTuple{2,Int},Bool}[]
    for pairing in pairings
        left = (vertex_map[pairing.source_vertex],
                mod1(pairing.source_slot - turns[pairing.source_vertex], 4))
        right = (vertex_map[pairing.target_vertex],
                 mod1(pairing.target_slot - turns[pairing.target_vertex], 4))
        push!(keys, _ordered_pair(left, right, pairing.orientation_reversing))
    end
    sort!(keys)
end

function _boundary_rows(boundary)
    [(component.euler, component.face_count, component.edge_count,
      component.vertex_count, component.orientable, component.closed_surface)
     for component in boundary.components]
end

function _section22_eulers(graph)
    surface = OGraphColourBoundary.boundary_from_coloured(graph)
    sort!(Int[component.euler for component in surface.components])
end

function _assert_surface_counts(graph)
    poly = face_cycles(graph)
    boundary = thickening_boundary(graph)
    pairings = dual_face_pairings(graph)
    orbits = dual_edge_orbits(graph)
    homology = cellular_homology(graph)
    n = poly.vertex_count
    @test poly.edge_count == 2n
    @test poly.halfedge_count == 4n
    @test poly.germ_count == 6n
    @test poly.boundary_length == 6n
    @test poly.euler == poly.disc_count - n
    @test sum(_cycle_lengths(poly)) == 6n
    @test boundary.euler_sum == 2 * poly.euler
    @test all(component.closed_surface && component.orientable && component.euler <= 2 &&
              component.euler == component.vertex_count - component.edge_count + component.face_count
              for component in boundary.components)
    @test length(pairings) == 2n
    @test all(pairing -> pairing.orientation_reversing, pairings)
    @test length(orbits) == poly.disc_count
    @test sorted_germ_partition(orbits) == _corner_partition(poly)
    ours = sort!(Int[component.euler for component in boundary.components])
    @test ours == _section22_eulers(graph)
    poly, boundary, pairings, homology
end

@testset "half-edge renaming leaves the public surface unchanged" begin
    samples = Tuple{String,ColouredOGraph,Bool}[]
    for name in _CLOSED_COLOUR_NAMES
        graph = normal_coloured_ograph(ograph_state(vl_example(name))).graph
        push!(samples, (name, graph, true))
    end
    push!(samples, ("one_vertex_00", _one_vertex(over=(2, 4), colours=(0, 0)), false))
    push!(samples, ("one_vertex_11", _one_vertex(over=(2, 4), colours=(1, 1)), false))
    push!(samples, ("one_vertex_over13", _one_vertex(over=(1, 3), colours=(0, 1)), false))
    lens = normal_coloured_ograph(ograph_state(vl_example("lens_2_1"))).graph
    record = lens.edge_colours[1]
    push!(samples, ("lens_2_1_colour_shift",
                    _replace_colour(lens, record.halfedge, record.colour + 1), false))
    rng = MersenneTwister(20261122)
    println("GAUGE renaming")
    for (label, graph, closed) in samples
        poly, boundary, pairings, homology = _assert_surface_counts(graph)
        n = poly.vertex_count
        base_rows = _boundary_rows(boundary)
        base_h1 = _h1_pair(homology)
        section = _section22_eulers(graph)
        println("GAUGE ", label, " n=", n, " d=", poly.disc_count, " chi=", poly.euler,
            " bd=", boundary.euler_sum, " bd22=", section,
            " comps=", [component.euler for component in boundary.components],
            " rev=", count(pairing -> pairing.orientation_reversing, pairings), "/",
            length(pairings), " h1=", base_h1)
        if closed
            @test poly.disc_count == n + 1
            @test poly.euler == 1
            @test boundary.euler_sum == 2
            @test length(boundary.components) == 1
        end
        for trial in 1:6
            σ = trial == 1 ? collect(1:length(graph.mate)) : randperm(rng, length(graph.mate))
            vertex_map = trial == 1 ? collect(1:n) : randperm(rng, n)
            turns = trial == 1 ? zeros(Int, n) : rand(rng, 0:3, n)
            image = _apply_presentation_gauge(graph, σ, vertex_map, turns)
            image_poly, image_boundary, image_pairs, image_homology =
                _assert_surface_counts(image)
            @test _transport_partition(poly, vertex_map, turns) == _corner_partition(image_poly)
            @test _transported_flags(pairings, vertex_map, turns) ==
                  _transported_flags(image_pairs, collect(1:n), zeros(Int, n))
            @test _boundary_rows(image_boundary) == base_rows
            @test image_boundary.euler_sum == boundary.euler_sum
            @test _h1_pair(image_homology) == base_h1
        end
    end
end
