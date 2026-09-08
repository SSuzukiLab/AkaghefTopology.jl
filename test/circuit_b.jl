using Test, Topology

@testset "Circuit B endpoint queries" begin
    v=VirtualLink()
    set_data!(v;rgauss=[[0]],orientation=Int[])
    paths=circuit_b_boundary(v)
    @test length(paths)==3
    @test all(p -> length(p.edges_s)==1 && p.dots==0,paths)
    incoming=circuit_b_at(v,1,:in)
    @test [p.endpoint.lane for p in incoming]==[1,2,3]
    @test [p.direction for p in incoming]==[:outgoing,:outgoing,:incoming]
    @test all(p -> p.other.strand==1 && p.other.side==:out,incoming)
    @test [only(p.edges_s) for p in incoming]==[1,1,-1]
    @test length(circuit_b_at(v,1,:out;lane=3))==1
    @test_throws ArgumentError circuit_b_at(v,1,:end)
    @test_throws ArgumentError circuit_b_at(v,1,:in;lane=4)
    closed=lens_space_link(2,1)
    @test isempty(circuit_b_boundary(closed))
    @test_throws ArgumentError circuit_b_at(closed,1,:in)
    mktempdir() do dir
        path=plot_circuit_b_svg(v,joinpath(dir,"interval.svg"))
        svg=read(path,String)
        @test length(collect(eachmatch(r"data-cycle=\"false\"",svg)))==3
        @test occursin("S1 in / lane 3",svg)
        @test occursin("S1 out / lane 3",svg)
    end
end
