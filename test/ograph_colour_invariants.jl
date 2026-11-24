using Test, Topology
using AbstractAlgebra

include(joinpath(@__DIR__, "..", "src", "ograph_colour_invariants.jl"))
include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))
const Moves = include(joinpath(@__DIR__, "..", "src", "ograph_colour_moves.jl"))

# Registry names from VLExample. The homology column is the group named for
# that registry entry. quaternion keeps the stored source string beside the
# name. odata_variant matches the pipeline and has no manifold source string.
const CLOSED = (
    (name="s3_abalone", free_rank=0, factors=BigInt[]),
    (name="lens_2_1", free_rank=0, factors=BigInt[2]),
    (name="lens_3_1_mst2_example23", free_rank=0, factors=BigInt[3]),
    (name="lens_3_2", free_rank=0, factors=BigInt[3]),
    (name="odata_variant", free_rank=0, factors=BigInt[]),
    (name="poincare", free_rank=0, factors=BigInt[]),
    (name="quaternion", free_rank=0, factors=BigInt[2, 2]),
    (name="s2xs1", free_rank=1, factors=BigInt[]),
    (name="s2xs1_2vertex", free_rank=1, factors=BigInt[]),
)

const SOURCE_STRINGS = Dict(
    "lens_3_1_mst2_example23" => "MST2 Example 23",
    "quaternion" => "K10 Fig. 14; positive trefoil o-graph",
    "s2xs1_2vertex" => "PJ260830: 261009_koda-bp-ograph-conventions",
)

function h1_of(graph)
    homology = cellular_homology(graph)
    (homology.free_rank, Tuple(homology.torsion_factors))
end

function h1_text(free_rank, factors)
    free_rank == 0 && isempty(factors) && return "0"
    torsion = join(("Z/$(factor)" for factor in factors), " x ")
    free_rank == 0 && return torsion
    isempty(factors) && return free_rank == 1 ? "Z" : "Z^$(free_rank)"
    "Z^$(free_rank) x $(torsion)"
end

function pipeline_h1(name)
    presentation = ograph_pi1_presentation(certify_closed_ograph(vl_example(name)))
    homology = ograph_h1(presentation)
    (homology.free_rank, Tuple(BigInt.(homology.torsion_invariant_factors)))
end

function abelianize(presentation::ColourPi1)
    rows = length(presentation.relators)
    cols = presentation.generator_count
    matrix = zeros(BigInt, rows, cols)
    for (row, relator) in enumerate(presentation.relators)
        for letter in relator
            matrix[row, abs(letter)] += sign(letter)
        end
    end
    diagonal = BigInt[]
    if rows > 0 && cols > 0
        smith = AbstractAlgebra.snf(AbstractAlgebra.matrix(AbstractAlgebra.ZZ, matrix))
        diagonal = [abs(BigInt(smith[index, index])) for index in 1:min(size(smith)...)]
    end
    nonzero = filter(!iszero, diagonal)
    factors = Tuple(value for value in nonzero if value > 1)
    (cols - length(nonzero), factors)
end

"""Test-only disc corners.

Chart A places leg 1 on the counter-clockwise side of an over slot. The
transport is ρ∘τ^{−α}. `docs/conventions/OGraphColours.md` records that this
gluing is the same as calibration B with transport ρ∘τ^α. This function is
not part of the library.
"""
function independent_disc_corners(graph::ColouredOGraph)
    owner = Dict{Int,Tuple{Int,Int}}()
    for (vertex, cyclic) in enumerate(graph.vertex_halfedges)
        for (slot, halfedge) in enumerate(cyclic)
            owner[Int(halfedge)] = (vertex, slot)
        end
    end
    colour_of = Dict{Tuple{Int,Int},Int}()
    for record in graph.edge_colours
        other = Int(graph.mate[record.halfedge])
        reduced = mod(Int(record.colour), 3)
        colour_of[minmax(Int(record.halfedge), other)] = reduced == 2 ? -1 : reduced
    end
    charts = Vector{NTuple{3,Int}}(undef, length(graph.mate))
    for (halfedge, (vertex, slot)) in owner
        cyclic = graph.vertex_halfedges[vertex]
        opposite = Int(cyclic[mod1(slot + 2, 4)])
        ccw = Int(cyclic[mod1(slot + 1, 4)])
        cw = Int(cyclic[mod1(slot - 1, 4)])
        upright = slot in graph.over_positions[vertex]
        charts[halfedge] = upright ? (ccw, opposite, cw) : (cw, opposite, ccw)
    end
    n = length(graph.mate)
    used = falses(n, 3)
    discs = Vector{NTuple{3,Int}}[]
    for start_halfedge in 1:n, start_leg in 1:3
        used[start_halfedge, start_leg] && continue
        corners = NTuple{3,Int}[]
        halfedge, leg = start_halfedge, start_leg
        while true
            used[halfedge, leg] && error("reused sheet end")
            used[halfedge, leg] = true
            other = Int(graph.mate[halfedge])
            alpha = colour_of[minmax(halfedge, other)]
            other_leg = fig25_frame_reversal(fig25_tau(-alpha, leg))
            used[other, other_leg] && error("reused arrival")
            used[other, other_leg] = true
            partner = charts[other][other_leg]
            back = findfirst(==(other), charts[partner])
            vertex, slot_a = owner[other]
            _, slot_b = owner[partner]
            push!(corners, (vertex, min(slot_a, slot_b), max(slot_a, slot_b)))
            if (partner, back) == (start_halfedge, start_leg)
                break
            end
            halfedge, leg = partner, back
        end
        push!(discs, corners)
    end
    discs
end

function corner_partition(graph)
    poly = face_cycles(graph)
    sorted_germ_partition(collect(cycle.corners) for cycle in poly.cycles)
end

function loop_graph(colours)
    ColouredOGraph(((1, 2, 3, 4),), (2, 1, 4, 3), ((2, 4),),
        (EdgeColour(1, colours[1]), EdgeColour(3, colours[2])))
end

function with_edge_colour(graph, edge_index, colour)
    records = EdgeColour[]
    seen = 0
    for (halfedge, other) in enumerate(graph.mate)
        halfedge < Int(other) || continue
        seen += 1
        value = seen == edge_index ? colour :
            only(item.colour for item in graph.edge_colours if item.halfedge == halfedge)
        push!(records, EdgeColour(halfedge, value))
    end
    ColouredOGraph(graph.vertex_halfedges, graph.mate, graph.over_positions, records;
                   convention_id=graph.convention_id)
end

function edge_number(graph, halfedge)
    seen = 0
    for (left, right) in enumerate(graph.mate)
        left < Int(right) || continue
        seen += 1
        if left == halfedge || Int(right) == halfedge
            return seen
        end
    end
    error("missing edge")
end

function rotate_slots(graph, shift)
    vertices = Tuple{Vararg{Int}}[]
    overs = Tuple{Vararg{Int}}[]
    for (cyclic, over) in zip(graph.vertex_halfedges, graph.over_positions)
        push!(vertices, Tuple(circshift(collect(Int.(cyclic)), shift)))
        push!(overs, Tuple(sort(mod1.(Int.(over) .+ shift, 4))))
    end
    ColouredOGraph(Tuple(vertices), graph.mate, Tuple(overs), graph.edge_colours;
                   convention_id=graph.convention_id)
end

function exchange_vertices(graph, left, right)
    vertices = collect(graph.vertex_halfedges)
    overs = collect(graph.over_positions)
    vertices[left], vertices[right] = vertices[right], vertices[left]
    overs[left], overs[right] = overs[right], overs[left]
    ColouredOGraph(Tuple(vertices), graph.mate, Tuple(overs), graph.edge_colours;
                   convention_id=graph.convention_id)
end

function counts_of(graph)
    poly = face_cycles(graph)
    boundary = thickening_boundary(graph)
    homology = h1_of(graph)
    (
        h1 = homology,
        discs = poly.disc_count,
        euler = poly.euler,
        boundary_euler = boundary.euler_sum,
        components = length(boundary.components),
        germs = poly.germ_count,
        pairs = length(dual_face_pairings(graph)),
    )
end

function assert_general(graph)
    poly = face_cycles(graph)
    n = poly.vertex_count
    @test poly.calibration == :b
    @test poly.status == "prop321-calibration-b/draft"
    @test poly.status == OGRAPH_FACE_CYCLE_STATUS
    @test poly.edge_count == 2 * n
    @test poly.halfedge_count == 4 * n
    @test poly.germ_count == 6 * n
    @test poly.boundary_length == 6 * n
    @test sum(length(cycle.germs) for cycle in poly.cycles) == 6 * n
    @test poly.euler == poly.disc_count - n
    boundary = thickening_boundary(graph)
    @test boundary.euler_sum == 2 * poly.euler
    for component in boundary.components
        @test component.closed_surface
        @test component.orientable
        @test component.euler <= 2
    end
    pairings = dual_face_pairings(graph)
    @test length(pairings) == 2 * n
    @test all(pairing -> pairing.orientation_reversing, pairings)
    @test all(pairing -> pairing.source_slot in 1:4 && pairing.target_slot in 1:4, pairings)
    @test length(dual_edge_orbits(graph)) == poly.disc_count
    @test length(dual_edge_orbits(graph)) == poly.vertex_count + poly.euler
    presentation = coloured_pi1(graph)
    @test abelianize(presentation) == h1_of(graph)
    @test sorted_germ_partition(independent_disc_corners(graph)) == corner_partition(graph)
    report = Moves.move_sanity(graph)
    @test report.calibration == :b
    @test report.discs == poly.disc_count
    @test report.boundary_euler == boundary.euler_sum
    poly, boundary
end

function move_label(move)
    if move isa Moves.CMove
        "C v=$(move.vertex) f=$(move.frame) $(move.direction)"
    elseif move isa Moves.MPForward
        "MP+ a=$(move.vertex_a) b=$(move.vertex_b) fa=$(move.frame_a) fb=$(move.frame_b)"
    else
        "MP- t=$(move.vertex_t) m=$(move.vertex_m) d=$(move.vertex_d)"
    end
end

function image_of(graph, move)
    moved = move isa Moves.CMove ? Moves.c_move(graph, move) : Moves.mp_move(graph, move)
    Moves.to_coloured(moved)
end

@testset "calibration B on the registered closed examples" begin
    @test OGRAPH_FACE_CYCLE_STATUS == "prop321-calibration-b/draft"
    println("closed name n d chi edges half germs bd comp h1 pipeline pi1")
    for item in CLOSED
        input = vl_example_input(item.name)
        if haskey(SOURCE_STRINGS, item.name)
            @test input.source == SOURCE_STRINGS[item.name]
        else
            @test !hasproperty(input, :source)
        end
        graph = normal_coloured_ograph(vl_example(item.name)).graph
        poly, boundary = assert_general(graph)
        n = poly.vertex_count
        @test poly.disc_count == n + 1
        @test poly.euler == 1
        @test length(boundary.components) == 1
        @test only(boundary.components).euler == 2
        @test boundary.euler_sum == 2
        homology = h1_of(graph)
        @test homology == (item.free_rank, Tuple(item.factors))
        @test pipeline_h1(item.name) == homology
        println(item.name, " ", n, " ", poly.disc_count, " ", poly.euler, " ",
            poly.edge_count, " ", poly.halfedge_count, " ", poly.germ_count, " ",
            boundary.euler_sum, " ",
            [(component.euler, component.orientable, component.closed_surface)
             for component in boundary.components], " ",
            h1_text(homology...), " ", h1_text(pipeline_h1(item.name)...), " ",
            h1_text(abelianize(coloured_pi1(graph))...))
    end
end

@testset "koda_fig13 is not a closed normal colouring" begin
    @test_throws ArgumentError normal_coloured_ograph(vl_example("koda_fig13"))
end

@testset "non-normal one-vertex colourings" begin
    partitions = Dict{Tuple{Int,Int},Any}()
    for colours in ((0, 0), (1, 1))
        graph = loop_graph(colours)
        @test graph.mate == (2, 1, 4, 3)
        @test graph.over_positions == ((2, 4),)
        poly, boundary = assert_general(graph)
        @test poly.vertex_count == 1
        @test poly.edge_count == 2
        @test poly.disc_count == 2
        @test poly.euler == 1
        @test length(boundary.components) == 1
        @test only(boundary.components).euler == 2
        homology = h1_of(graph)
        @test homology == (0, (BigInt(4),))
        partitions[colours] = corner_partition(graph)
        component = only(boundary.components)
        println("loops ", colours, " n=1 d=2 chi=1 H1=", h1_text(homology...),
            " boundary=(", component.euler, ", orientable=", component.orientable,
            ", closed=", component.closed_surface, ", F=", component.face_count,
            ", E=", component.edge_count, ", V=", component.vertex_count, ")")
    end
    @test partitions[(0, 0)] != partitions[(1, 1)]
end

@testset "negative control on s3_abalone" begin
    graph = normal_coloured_ograph(vl_example("s3_abalone")).graph
    record = only(item for item in graph.edge_colours if item.colour == -1)
    control = with_edge_colour(graph, edge_number(graph, record.halfedge), 1)
    loops = loop_graph((1, 1))
    @test control.vertex_halfedges == loops.vertex_halfedges
    @test control.mate == loops.mate
    @test control.over_positions == loops.over_positions
    @test control.edge_colours == loops.edge_colours
    poly, boundary = assert_general(control)
    homology = h1_of(control)
    @test poly.disc_count == 2
    @test poly.euler == 1
    @test homology == (0, (BigInt(4),))
    @test boundary.euler_sum == 2
    println("negative control H1=", h1_text(homology...),
        " d=", poly.disc_count, " chi=", poly.euler,
        " boundary=", boundary.euler_sum)
end

@testset "gauge" begin
    graphs = [normal_coloured_ograph(vl_example(item.name)).graph for item in CLOSED]
    push!(graphs, loop_graph((0, 0)), loop_graph((1, 1)))
    for graph in graphs
        base = counts_of(graph)
        for shift in 1:3
            turned = rotate_slots(graph, shift)
            @test counts_of(turned) == base
            Moves.move_sanity(turned)
            @test abelianize(coloured_pi1(turned)) == base.h1
        end
        if length(graph.vertex_halfedges) >= 2
            swapped = exchange_vertices(graph, 1, 2)
            @test counts_of(swapped) == base
            Moves.move_sanity(swapped)
        end
        canon = canonical_coloured_ograph(graph)
        @test counts_of(canon) == base
        Moves.move_sanity(canon)
    end
    println("gauge: slot rotations 1, 2, 3, the exchange of vertices 1 and 2, and canonical_coloured_ograph preserve H1, disc count, chi, and boundary euler")
end

@testset "candidate A fails the closed-spine count on s2xs1_2vertex" begin
    graph = normal_coloured_ograph(vl_example("s2xs1_2vertex")).graph
    poly = face_cycles(graph; calibration=:a)
    boundary = thickening_boundary(graph; calibration=:a)
    @test poly.vertex_count == 2
    @test poly.disc_count == 5
    @test poly.euler == 3
    @test poly.germ_count == 12
    @test boundary.euler_sum == 6
    @test length(boundary.components) == 3
    @test all(component -> component.euler == 2 && component.orientable &&
            component.closed_surface, boundary.components)
    println("candidate A s2xs1_2vertex n=2 d=5 chi=3 boundary_components=3 each chi 2 sum 6")
end

@testset "ordered tetrahedra of s3_abalone" begin
    ordered = disk_tetrahedra(vl_example("s3_abalone"))
    @test ordered.face_orientations_verified
    @test ordered.edge_orientations_verified
    @test length(ordered.tetrahedra) == 1
end

@testset "C and MP where the image satisfies the sanity rule" begin
    failures = String[]
    passing_c = 0
    passing_mp = 0
    lens_forward_c = false
    lens_forward_mp = false
    abalone_c = 0
    for item in CLOSED
        coloured = normal_coloured_ograph(vl_example(item.name)).graph
        graph = Moves.from_coloured(coloured)
        before = h1_of(coloured)
        before_counts = counts_of(coloured)
        for move in vcat(Moves.applicable_c(graph), Moves.applicable_mp(graph))
            image = image_of(graph, move)
            label = item.name * " " * move_label(move)
            local report
            try
                report = Moves.move_sanity(image)
            catch error
                push!(failures, label * " " * sprint(showerror, error))
                continue
            end
            got = (report.h1_free, Tuple(report.h1_torsion))
            @test got == before
            if move isa Moves.CMove
                passing_c += 1
                item.name == "s3_abalone" && (abalone_c += 1)
                if item.name == "lens_2_1" && move.direction == :forward
                    lens_forward_c = true
                    @test report.discs == before_counts.discs
                    @test report.euler == before_counts.euler
                    @test report.boundary_euler == 2
                end
            else
                passing_mp += 1
                if move isa Moves.MPForward && item.name == "lens_2_1"
                    lens_forward_mp = true
                    @test report.n == length(coloured.vertex_halfedges) + 1
                    @test report.discs == before_counts.discs + 1
                    @test report.euler == before_counts.euler
                    @test report.boundary_euler == 2
                    @test report.component_euler == [2]
                end
            end
        end
    end
    @test abalone_c == 4
    @test lens_forward_c
    @test lens_forward_mp
    @test passing_c > 0
    @test passing_mp > 0
    println("C images passing sanity: ", passing_c)
    println("MP images passing sanity: ", passing_mp)
    println("C/MP images failing sanity: ", length(failures))
    for failure in failures
        println(failure)
    end
end
