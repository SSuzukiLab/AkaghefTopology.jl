using Topology, JSON3, Test
const H=Topology.Heegaard
viewpath=joinpath(@__DIR__,"../../artifacts/heegaard/full_curves/KodaCommon.json")
@testset "saved view retains the full marked conversion source" begin
    view=JSON3.read(read(viewpath,String),Dict{String,Any})
    source=H.marked_display_source(view)
    @test source==view["source_document"]
    bad=deepcopy(view);pop!(bad["disk"])
    @test_throws ArgumentError H.marked_display_source(bad)
    bad=deepcopy(view);bad["crossings"][1]["sign"]*=-1
    @test_throws ArgumentError H.marked_display_source(bad)
    bad=deepcopy(view);reverse!(bad["curves"][1]["segments"])
    @test_throws ArgumentError H.marked_display_source(bad)
end
