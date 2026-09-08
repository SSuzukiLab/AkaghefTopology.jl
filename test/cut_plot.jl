using Test
using Topology

@testset "cut SVG rendering" begin
    cut=VirtualLink()
    set_data!(cut;gauss=[[1,-1,0]],orientation=[1])
    @test_throws ArgumentError orthogonal_layout(cut)

    output=mktempdir()
    path=plot_svg(cut,joinpath(output,"cut.svg");edge_labels=:all)
    svg=read(path,String)
    edges=real_edge_table(cut)
    @test length(edges)==3
    @test count("boundary-endpoint",svg)==2
    @test count("data-side=",svg)==2
    @test count("data-boundary=\"true\"",svg)==2
    @test count("class=\"arc\" data-arc",svg)==length(edges)
    @test count("data-arc=\"$(first(edges).id)\"",svg)>=2
    @test !occursin("Z",svg)

    crossing_positions=Dict(1=>0.0+0.0im)
    replay=[ComplexF64[
        edge.crossing[1]==0 ? -1.0im : crossing_positions[abs(edge.crossing[1])],
        edge.crossing[2]==0 ? 1.0im : crossing_positions[abs(edge.crossing[2])]]
             for edge in edges]
    replay_path=plot_svg(cut,joinpath(output,"cut-replay.svg");
                         replay_positions=replay,edge_labels=:id)
    @test count("boundary-endpoint",read(replay_path,String))==2

    mixed=VirtualLink()
    set_data!(mixed;gauss=[[1,-1,0],[2,-2]],orientation=[1,1])
    mixed_path=plot_svg(mixed,joinpath(output,"mixed.svg");edge_labels=:id)
    mixed_svg=read(mixed_path,String)
    @test count("boundary-endpoint",mixed_svg)==2
    @test count("class=\"arc\"",mixed_svg)==length(real_edge_table(mixed))

    mixed_empty=VirtualLink()
    set_data!(mixed_empty;gauss=[Int[],[0]],orientation=Int[])
    mixed_empty_path=plot_svg(mixed_empty,joinpath(output,"mixed-empty.svg"))
    mixed_empty_svg=read(mixed_empty_path,String)
    @test count("boundary-endpoint",mixed_empty_svg)==2
    @test count("<circle class=\"arc\"",mixed_empty_svg)==1

    mixed_empty_real=VirtualLink()
    set_data!(mixed_empty_real;gauss=[Int[],[1,-1,0]],orientation=[1])
    mixed_empty_real_path=plot_svg(mixed_empty_real,joinpath(output,"mixed-empty-real.svg"))
    mixed_empty_real_svg=read(mixed_empty_real_path,String)
    @test count("boundary-endpoint",mixed_empty_real_svg)==2
    @test count("<circle class=\"arc\"",mixed_empty_real_svg)==1

    pure=VirtualLink()
    set_data!(pure;gauss=[[0]],orientation=Int[])
    pure_path=plot_svg(pure,joinpath(output,"pure.svg"))
    pure_svg=read(pure_path,String)
    @test count("boundary-endpoint",pure_svg)==2
    @test !occursin("title>unknot",pure_svg)
end
