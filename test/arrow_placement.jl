using Test, Topology

@testset "Direction marker avoids polyline vertices" begin
    points=ComplexF64[0,1,1+im]
    arclength_midpoint,_=Topology._arclength_midpoint(points)
    marker_midpoint,marker_tangent=Topology._middle_segment_midpoint(points)
    @test arclength_midpoint==1
    @test marker_midpoint==0.5
    @test marker_tangent==1
    @test all(point -> !isapprox(marker_midpoint,point;atol=1e-9),points)

    whitehead=VirtualLink()
    set_data!(whitehead;gauss=[[1,-4,5,-3],[3,-1,2,-5,4,-2]],
              orientation=[-1,-1,-1,1,1])
    layout=orthogonal_layout(whitehead)
    arcs=Topology._real_arc_layouts(whitehead,layout)
    vertices=[layout.crossings[index] for index in eachindex(layout.crossings)
              if layout.orientation[index]!=0]
    for arc in arcs
        midpoint,_=Topology._middle_segment_midpoint(arc.points)
        @test all(point -> !isapprox(midpoint,point;atol=1e-9),arc.points)
        @test all(vertex -> abs(midpoint-vertex)>Topology.LINK_STYLE.vertex,vertices)
    end

    output=joinpath(mktempdir(),"whitehead.svg")
    svg=read(plot_svg(whitehead,output),String)
    @test count("class=\"arrow\" data-arc=",svg)==length(arcs)
    @test occursin("stroke-width=\"0.07\"",svg)
    @test occursin("stroke-width=\"0.06\"",svg)
end
