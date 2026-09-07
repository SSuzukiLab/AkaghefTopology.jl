using Test, Topology

@testset "disk to ordered tetrahedra" begin
    for (vertex_sign,word,expected) in
        ((1,[1,-1],(1,5,6,3,2,4)),(-1,[-1,1],(1,3,5,2,6,4)))
        v=VirtualLink()
        set_data!(v;rgauss=[word],orientation=[vertex_sign])
        result=disk_tetrahedra(v)
        t=only(result.tetrahedra)
        @test t.local_nodes==expected # source Figure 2.2 + production germ labels
        @test t.sign == -vertex_sign
        @test result.disk_count==2
        @test all(result.closedness.criteria)
        @test result.edge_orientations_verified
        @test result.face_orientations_verified
        @test sort(unique(result.node_to_disk))==collect(1:result.disk_count)
        @test t.edges==Tuple(result.node_to_disk[germ] for germ in expected)
    end
    empty=VirtualLink()
    set_data!(empty;rgauss=[Int[]],orientation=Int[])
    @test_throws ArgumentError disk_tetrahedra(empty)
end
