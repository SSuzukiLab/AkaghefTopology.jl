using Test
using Topology

@testset "cut core representations" begin
    closed=VirtualLink()
    set_data!(closed;rgauss=[[1,-1]],orientation=[1])
    @test strand_cuts(closed)==[false]
    @test real_gauss_code(closed)==([[1,-1]],[1])
    @test [edge.crossing for edge in real_edge_table(closed)]==[(-1,1),(1,-1)]

    open=deepcopy(closed)
    set_cuts!(open,[true])
    @test strand_cuts(open)==[true]
    @test gauss_code(open)==[[1,-1,0]]
    @test real_gauss_code(open)==([[1,-1,0]],[1])
    @test [edge.crossing for edge in real_edge_table(open)]==[(0,1),(1,-1),(-1,0)]
    set_weight!(open,[[10,20,30]])
    @test [edge.weight for edge in real_edge_table(open)]==[10.0,20.0,30.0]
    @test head_map(open)[1,2]==-1
    @test head_map(open)[1,1]==0
    @test_throws ArgumentError cyclic_shift!(open,1,1)
    @test cyclic_shift!(open,1,0)===open

    repeated=VirtualLink()
    set_data!(repeated;gauss=[[1,2,-1,-2]],orientation=[1,0],cut=[true])
    arcs=real_edge_table(repeated)
    @test [edge.crossing for edge in arcs]==[(0,1),(1,-1),(-1,0)]
    @test arcs[1].arc==[1]
    @test arcs[2].arc==[1,2,-1]
    @test arcs[3].arc==[-1,-2]

    empty_interval=VirtualLink()
    set_data!(empty_interval;rgauss=[[0]],orientation=Int[])
    @test strand_cuts(empty_interval)==[true]
    @test real_gauss_code(empty_interval)==([[0]],Int[])
    @test [edge.crossing for edge in real_edge_table(empty_interval)]==[(0,0)]
    disk=disk_table(empty_interval)
    @test length(disk.boundary_nodes)==2
    @test length(disk.boundary_lanes)==6
    @test length(disk.open_paths)==3
    @test all(path -> !isempty(path),disk.open_paths)

    two=VirtualLink()
    set_data!(two;gauss=[[1,-1]],orientation=[1],cut=[true])
    pd=pd_code(two)
    @test length(unique(pd))>=3
    roundtrip=VirtualLink(); set_data!(roundtrip;pd=pd)
    @test strand_cuts(roundtrip)==[true]
    @test real_gauss_code(roundtrip)==([[1,-1,0]],[1])

    mixed=VirtualLink()
    set_data!(mixed;gauss=[[1,2,0],[-1,-2]],orientation=[1,-1])
    mixed_pd=pd_code(mixed)
    mixed_roundtrip=VirtualLink()
    set_data!(mixed_roundtrip;pd=mixed_pd,orientation=[1,-1])
    @test strand_cuts(mixed_roundtrip)==[true,false]
    @test real_gauss_code(mixed_roundtrip)==([[1,2,0],[-1,-2]],[1,-1])

    mirrored=deepcopy(open)
    # No source-backed framing lift: reject without changing the cut or weights.
    @test_throws ArgumentError mirror_manifold!(mirrored)
    @test mirrored.gauss==open.gauss && mirrored.weights==open.weights
    mirrored=deepcopy(two) # same cut geometry, explicitly unweighted
    mirror_manifold!(mirrored)
    @test strand_cuts(mirrored)==[true]
    sum_link=VirtualLink()
    set_data!(sum_link;gauss=[[1,-1,0]],orientation=[1])
    other=VirtualLink(); set_data!(other;gauss=[[1,-1,0]],orientation=[1])
    disjoint_sum!(sum_link,other)
    @test strand_cuts(sum_link)==[true,true]
    @test_throws ArgumentError connected_sum!(sum_link)
end

@testset "cut marker and flag validation" begin
    v=VirtualLink(); set_data!(v;gauss=[[1,-1]],orientation=[1])
    @test_throws DimensionMismatch set_cuts!(v,Bool[])
    @test_throws ArgumentError set_cuts!(v,[1])
    v.gauss[1]=[1,0,-1]
    @test_throws ArgumentError strand_cuts(v)
end
