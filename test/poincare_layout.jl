using Test, Topology, JSON

@testset "Poincare MATLAB/Sage routing orientation" begin
    fixture=JSON.parsefile(joinpath(@__DIR__,"fixtures","poincare_matlab_layout_v1.json"))
    gauss=[Int.(component) for component in fixture["gauss"]]
    orientation=Int.(fixture["orientation"])
    source=VirtualLink()
    set_data!(source;gauss=gauss,orientation=orientation)

    pd=pd_code(source)
    bending=minimal_bending_numbers(pd)
    layout=orthogonal_layout(source)

    @test pd==reduce(vcat,(permutedims(Int.(row)) for row in fixture["pd"]))
    @test bending==Int.(fixture["bending_numbers"])
    expected_crossings=ComplexF64[complex(point[1],point[2]) for point in fixture["crossings"]]
    @test layout.crossings==expected_crossings
    @test layout.orientation==orientation
    @test replace(layout.orientation,0=>1)==Int.(fixture["routing_orientation"])
    @test isempty(check_layout(layout))
    @test all(edge -> first(edge.points)==layout.crossings[edge.tail] &&
                      last(edge.points)==layout.crossings[edge.head],layout.edges)
end
