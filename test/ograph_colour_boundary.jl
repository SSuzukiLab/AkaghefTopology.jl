using Test
using Topology

include(joinpath(@__DIR__, "..", "src", "ograph_colour_boundary.jl"))

include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))

const B = OGraphColourBoundary

const BOUNDARY_REPORT_NAMES = (
    "s3_abalone",
    "lens_2_1",
    "lens_3_1_mst2_example23",
    "poincare",
    "quaternion",
    "s2xs1",
    "s2xs1_2vertex",
)

# Closed normal o-graphs. `koda_fig13` is outside this list.
const CLOSED_NORMAL_NAMES = (
    BOUNDARY_REPORT_NAMES...,
    "lens_3_2",
    "odata_variant",
)

# PD slots 1..4. This is the DIR placement in the spec.
const _SLOT_PORTS = Dict(
    1 => (:mi, :po, :mo, :pi),
    -1 => (:mi, :pi, :mo, :po),
)

function _surface_rows(surface)
    rows = NTuple{4,Int}[
        (component.theta_vertices, component.theta_edges,
         component.ribbon_cycles, component.euler)
        for component in surface.components
    ]
    sort!(rows)
end

function _coloured_from_words(words, signs)
    n = length(signs)
    mate = zeros(Int, 4n)
    over = [Int[] for _ in 1:n]
    colours = EdgeColour[]
    for vertex in 1:n
        for (slot, port) in enumerate(_SLOT_PORTS[signs[vertex]])
            (port === :po || port === :pi) && push!(over[vertex], slot)
        end
    end
    for word in words
        m = length(word)
        for index in 1:m
            tail = Int(word[index])
            head = Int(word[mod1(index + 1, m)])
            tail_vertex = abs(tail)
            head_vertex = abs(head)
            tail_slot = findfirst(==(tail > 0 ? :po : :mo), _SLOT_PORTS[signs[tail_vertex]])
            head_slot = findfirst(==(head > 0 ? :pi : :mi), _SLOT_PORTS[signs[head_vertex]])
            tail_halfedge = 4 * (tail_vertex - 1) + tail_slot
            head_halfedge = 4 * (head_vertex - 1) + head_slot
            mate[tail_halfedge] = head_halfedge
            mate[head_halfedge] = tail_halfedge
            push!(colours, EdgeColour(min(tail_halfedge, head_halfedge),
                                       B.prop321_colour(tail > 0, head > 0)))
        end
    end
    cycles = Tuple(ntuple(slot -> 4 * (vertex - 1) + slot, 4) for vertex in 1:n)
    ColouredOGraph(cycles, Tuple(mate), Tuple(Tuple(slots) for slots in over), Tuple(colours))
end

@testset "Proposition 3.2.1 colour of one step" begin
    @test B.prop321_colour(true, true) == 0
    @test B.prop321_colour(false, false) == 0
    @test B.prop321_colour(false, true) == 1
    @test B.prop321_colour(true, false) == -1
end

@testset "one-vertex abalone, including the reversed mixed edge" begin
    words = [[1, -1]]
    signs = [1]
    gauss = B.boundary_from_gauss(words, signs)
    # Slots S, E, N, W are half-edges 1, 2, 3, 4.
    # East–South has colour -1, so the smaller half-edge is the south end.
    # North–West has colour +1, and the smaller half-edge is the north end.
    hand = ColouredOGraph(
        ((1, 2, 3, 4),),
        (2, 1, 4, 3),
        ((2, 4),),
        (EdgeColour(1, -1), EdgeColour(3, 1)),
    )
    @test gauss.reading == :provisional
    @test _surface_rows(gauss) == _surface_rows(B.boundary_from_coloured(hand))
    @test length(gauss.components) == 1
    component = gauss.components[1]
    # §2.2 for one o-graph vertex: V = 4, E = 6. One disc on each unoriented circle.
    @test component.theta_vertices == 4
    @test component.theta_edges == 6
    @test component.ribbon_cycles == 4
    @test component.euler == 2
    @test sum(component.orbit_degrees) == 12
    @test length(component.orbit_degrees) == 4
    println("ABALONE degrees=$(component.orbit_degrees)")
end

@testset "a both-above colour +1 is independent of which end is the tail" begin
    # Slots 2 and 4 are the over pair. East–West has colour +1.
    east_tail = ColouredOGraph(
        ((1, 2, 3, 4),),
        (3, 4, 1, 2),
        ((2, 4),),
        (EdgeColour(1, 0), EdgeColour(2, 1)),
    )
    # Renumber so the west end is the smaller half-edge of the same edge.
    west_tail = ColouredOGraph(
        ((1, 4, 3, 2),),
        (3, 4, 1, 2),
        ((2, 4),),
        (EdgeColour(1, 0), EdgeColour(2, 1)),
    )
    @test _surface_rows(B.boundary_from_coloured(east_tail)) ==
          _surface_rows(B.boundary_from_coloured(west_tail))
end

@testset "registered inputs under the Proposition 3.2.1 colouring and one shift" begin
    seen = String[]
    for input in vl_virtual_gauss_inputs()
        push!(seen, input.name)
        signs = Int.(input.orientation)
        n = length(signs)
        surface = B.boundary_from_gauss(input.rgauss, signs)
        altered = B.boundary_from_gauss(input.rgauss, signs; alter_edge=1)
        coloured = B.boundary_from_coloured(_coloured_from_words(input.rgauss, signs))
        @test surface.reading == :provisional
        @test altered.reading == :provisional
        @test _surface_rows(surface) == _surface_rows(coloured)
        @test sum(component.theta_vertices for component in surface.components) == 4n
        @test sum(component.theta_edges for component in surface.components) == 6n
        @test sum(component.theta_vertices for component in altered.components) == 4n
        @test sum(component.theta_edges for component in altered.components) == 6n
        for component in surface.components
            @test component.euler ==
                  component.theta_vertices - component.theta_edges + component.ribbon_cycles
            @test sum(component.orbit_degrees) == 2 * component.theta_edges
            @test length(component.orbit_degrees) == component.ribbon_cycles
        end
        for component in altered.components
            @test component.euler ==
                  component.theta_vertices - component.theta_edges + component.ribbon_cycles
            @test sum(component.orbit_degrees) == 2 * component.theta_edges
        end
        # Acceptance check after the circle count. The arm placement is not chosen from it.
        # A colour shifted by +1 is printed and is not required to be a sphere.
        if input.name in CLOSED_NORMAL_NAMES
            @test length(surface.components) == 1
            if length(surface.components) == 1
                @test surface.components[1].euler == 2
                @test surface.components[1].ribbon_cycles == 2n + 2
            end
        end
        println("EXAMPLE $(input.name) n=$n canonical=$(_surface_rows(surface)) degrees=$([component.orbit_degrees for component in surface.components]) alter1=$(_surface_rows(altered))")
    end
    for name in BOUNDARY_REPORT_NAMES
        @test name in seen
    end
end

@testset "an all-zero colouring still runs" begin
    signs = [1]
    graph = _coloured_from_words([[1, -1]], signs)
    zeros = ColouredOGraph(
        graph.vertex_halfedges,
        graph.mate,
        graph.over_positions,
        Tuple(EdgeColour(record.halfedge, 0) for record in graph.edge_colours),
    )
    surface = B.boundary_from_coloured(zeros)
    @test surface.reading == :provisional
    @test sum(component.theta_vertices for component in surface.components) == 4
    @test sum(component.theta_edges for component in surface.components) == 6
end

@testset "a non-opposite over pair is rejected" begin
    bad = ColouredOGraph(
        ((1, 2, 3, 4),),
        (2, 1, 4, 3),
        ((1, 2),),
        (EdgeColour(1, 0), EdgeColour(3, 0)),
    )
    @test_throws ArgumentError B.boundary_from_coloured(bad)
end
