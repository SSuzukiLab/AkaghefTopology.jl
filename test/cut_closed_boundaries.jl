using Test
using Topology

@testset "closed APIs reject raw cut boundaries" begin
    closed = lens_space_link(2, 1)
    certified = certify_closed_ograph(closed)
    @test certified.witness.criteria == (true, true, true)

    cut = lens_space_link(2, 1)
    cut.gauss[1] = [cut.gauss[1]..., 0]

    # The raw cut marker must survive the real-Gauss view for boundary-aware
    # callers, while closed-only APIs reject it before normalization.
    code, orientation = real_gauss_code(cut)
    @test last(only(code)) == 0
    @test orientation == [1, 1]

    @test_throws ArgumentError ograph_state(cut)
    @test_throws ArgumentError certify_closed_ograph(cut)
    @test_throws ArgumentError c2_graph(cut)
    @test_throws ArgumentError disk_tetrahedra(cut)
    @test_throws ArgumentError hopf_invariant(cut, exterior_algebra(1))
end

@testset "raw cut edge and disk tables retain open boundaries" begin
    cut = lens_space_link(2, 1)
    cut.gauss[1] = [cut.gauss[1]..., 0]

    edges = real_edge_table(cut)
    @test length(edges) == 5
    @test first(edges).crossing[1] == 0
    @test last(edges).crossing[2] == 0

    disks = disk_table(cut)
    @test !isempty(disks.ID)
    @test any(!, disks.isCycle)
    @test length(disks.paths) == length(disks.ID)
    @test length(disks.edges_s) == length(disks.ID)
end
