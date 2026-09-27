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

@testset "Deterministic framed lens-space examples" begin
    input=framed_lens_space_input(5,2)
    @test input.schema=="framed-lens-o-graph-v2"
    @test input.parameters.r==4
    @test input.parameters.A==4
    @test length(input.orientation)==20
    @test input.weights==[[4,0,0,0,-1,5,-3,-1,5,-3,-1,1,0,0,0,0,0,0,-1,0,
                           0,0,0,0,0,0,0,1,0,0,-1,-1,0,-1,0,0,-1,1,-1,-3]]

    link=framed_lens_space_link(5,2)
    @test real_gauss_code(link)==(input.rgauss,input.orientation)
    @test BigInt.(only(link.weights))==only(input.weights)
    @test is_closed(link)==(true,(true,true,true))
    @test is_framed_weight(link)

    qminus=framed_lens_space_input(4,3)
    @test qminus.construction=="previous_q_minus_one_base"
    @test qminus.rgauss==[[-4,-1,-2,-3,1,2,3,4]]
    @test qminus.weights==[[1,0,0,0,-1,0,0,-1]]

    series=framed_lens_space_series(3)
    @test length(series)==3
    @test [(row.p,row.q) for row in series]==
          [(2,1),(3,1),(3,2)]
    series[1].link.gauss[1][1]=999
    @test series[2].link.gauss[1][1]!=999

    @test_throws ArgumentError framed_lens_space_input(5,0)
    @test_throws ArgumentError framed_lens_space_input(6,2)
    @test_throws ArgumentError framed_lens_space_input(true,1)
    @test_throws ArgumentError framed_lens_space_series(5;minimum_p=6)
end
