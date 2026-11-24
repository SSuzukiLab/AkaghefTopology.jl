using Test
using Topology

isdefined(@__MODULE__, :OGraphColourBoundary) ||
    include(joinpath(@__DIR__, "..", "src", "ograph_colour_boundary.jl"))
include(joinpath(@__DIR__, "..", "examples", "BingHouse.jl"))
isdefined(@__MODULE__, :vl_example) ||
    include(joinpath(@__DIR__, "..", "examples", "VLExample.jl"))

# cf:[[BP97#^bp97-sec-0565]] Corollary 3.1.7 — each incident germ carries its own direction
# cf:[[Cas65#^cas65-sec-0016]] X minus the singular 1-skeleton is three open disks

function _canon_letters(word)
    best = nothing
    seq = collect(Int.(word))
    for shift in 0:length(seq)-1
        rotated = vcat(seq[shift+1:end], seq[1:shift])
        (best === nothing || rotated < best) && (best = rotated)
    end
    rev = Int[-value for value in reverse(seq)]
    for shift in 0:length(rev)-1
        rotated = vcat(rev[shift+1:end], rev[1:shift])
        (best === nothing || rotated < best) && (best = rotated)
    end
    best
end

function _word_key(poly)
    sort!(collect(_canon_letters(cycle.edge_letters) for cycle in poly.cycles))
end

const _BING_WORDS = sort!(collect(_canon_letters(word) for word in (
    [-4],
    [-4, -3, -1, 3, -2, 1, 2, 4, -2, 3],
    [-1],
)))

"""Region orientations of Corollary 3.1.7. A repeated passage is one germ."""
function _branching_failures(poly)
    words = [cycle.edge_letters for cycle in poly.cycles]
    edges = sort(unique(abs(letter) for word in words for letter in word))
    failures = 0
    for mask in 0:(1 << length(words)) - 1
        for edge in edges
            signs = Int[]
            for (region, word) in enumerate(words)
                bit = ((mask >> (region - 1)) & 1) == 1 ? 1 : -1
                for letter in word
                    abs(letter) == edge && push!(signs, bit * sign(letter))
                end
            end
            if length(signs) != 3 || abs(sum(signs)) == 3
                failures += 1
                break
            end
        end
    end
    failures
end

function _strand_lengths(graph, orientation)
    nxt = Dict{Int,Int}(edge.tail => edge.head for edge in orientation)
    opposite = Dict{Int,Int}()
    for cyclic in graph.vertex_halfedges
        for slot in 1:4
            opposite[cyclic[slot]] = cyclic[mod1(slot + 2, 4)]
        end
    end
    seen = Set{Int}()
    lengths = Int[]
    for edge in orientation
        edge.tail in seen && continue
        tail = edge.tail
        count = 0
        while !(tail in seen)
            push!(seen, tail)
            count += 1
            tail = opposite[nxt[tail]]
        end
        push!(lengths, count)
    end
    sort!(lengths)
end

function _rotate_slots(graph, vertex, turn)
    cyclic = collect(graph.vertex_halfedges[vertex])
    rotated = vcat(cyclic[turn+1:end], cyclic[1:turn])
    vertices = ntuple(index -> index == vertex ? Tuple(rotated) : graph.vertex_halfedges[index],
                      length(graph.vertex_halfedges))
    overs = ntuple(length(graph.vertex_halfedges)) do index
        index == vertex || return graph.over_positions[index]
        Tuple(sort(Int[mod1(position - turn, 4) for position in graph.over_positions[index]]))
    end
    ColouredOGraph(vertices, graph.mate, overs, graph.edge_colours)
end

@testset "Bing's house coloured o-graph" begin
    graph = bing_house_coloured_ograph()
    @test all(check -> check.passed, validate_coloured_ograph(graph))
    poly = face_cycles(graph)
    boundary = thickening_boundary(graph)
    pairings = dual_face_pairings(graph)
    orbits = dual_edge_orbits(graph)
    homology = cellular_homology(graph)
    section = OGraphColourBoundary.boundary_from_coloured(graph)
    n = poly.vertex_count
    println("BING n=", n,
        " E=", poly.edge_count,
        " halfedges=", poly.halfedge_count,
        " germs=", poly.germ_count,
        " d=", poly.disc_count,
        " chi=", poly.euler,
        " length_sum=", poly.boundary_length,
        " bd=", [component.euler for component in boundary.components],
        " bd_sum=", boundary.euler_sum,
        " bd22=", [component.euler for component in section.components],
        " pairs=", length(pairings),
        " reversing=", count(pairing -> pairing.orientation_reversing, pairings),
        " dual_edges=", length(orbits),
        " h1=", homology.free_rank, " ", homology.torsion_factors,
        " n2=", length(normal_orientations(graph)),
        " p2=", length(prop321_orientations(graph)))
    @test n == 2
    @test poly.edge_count == 4
    @test poly.halfedge_count == 8
    @test poly.germ_count == 12
    @test poly.disc_count == 3
    @test poly.euler == 1
    @test poly.boundary_length == 12
    @test sum(length(cycle.edge_letters) for cycle in poly.cycles) == 12
    @test _word_key(poly) == _BING_WORDS
    @test boundary.euler_sum == 2
    @test length(boundary.components) == 1
    @test boundary.components[1].euler == 2
    @test boundary.components[1].closed_surface
    @test boundary.components[1].orientable
    @test sort(Int[component.euler for component in section.components]) == [2]
    @test length(pairings) == 4
    @test all(pairing -> pairing.orientation_reversing, pairings)
    @test length(orbits) == 3
    @test homology.free_rank == 0
    @test isempty(homology.torsion_factors)
    @test _branching_failures(poly) == 0
    orientations = prop321_orientations(graph)
    @test length(orientations) == 4
    @test all(orientation -> _strand_lengths(graph, orientation) == [1, 1, 2], orientations)
    canonical = coloured_ograph_serialization(canonical_coloured_ograph(graph))
    for turn in 0:3
        turned = _rotate_slots(graph, 1, turn)
        @test coloured_ograph_serialization(canonical_coloured_ograph(turned)) == canonical
    end
    for input in vl_virtual_gauss_inputs()
        recovered = try
            normal_coloured_ograph(vl_example(input.name)).graph
        catch
            continue
        end
        other = coloured_ograph_serialization(canonical_coloured_ograph(recovered))
        @test other != canonical
    end
end
