using Test, Topology

@testset "Lens-space vg2g issue #1" begin
    l09=VirtualLink()
    set_data!(l09;rgauss=[[-9,-7,-5,-3,-1,-8,-6,-4,-2,1,2,3,4,5,6,7,8,9]],
              orientation=ones(Int,9))
    source=deepcopy(l09.gauss)
    gauss,orientation=planarize_virtual_gauss(l09.gauss,l09.orientation;DEBUG_MODE=true,debug_io=devnull)
    @test l09.gauss==source
    @test length(gauss)==1
    @test orientation[1:9]==ones(Int,9)
    @test count(==(0),orientation)==9
    labels=abs.(filter(!=(0),gauss[1]))
    @test sort(unique(labels))==collect(1:18)
    @test all(label -> count(==(label),labels)==2,unique(labels))
    preserved=[entry for entry in gauss[1] if entry!=0 && abs(entry)<=9]
    @test preserved==source[1]

    planar=VirtualLink()
    set_data!(planar;gauss=gauss,orientation=orientation)
    pd=pd_code(planar)
    @test length(regions(pd))==size(pd,1)+2
    @test isempty(check_layout(orthogonal_layout(l09)))

    failures=[(5,2),(7,2),(7,3),(8,3),(9,2),(9,4),(10,3)]
    for (s,t) in failures
        link=lens_space_link(s,t)
        @test_nowarn planarize_virtual_gauss(link.gauss,link.orientation)
    end
end
