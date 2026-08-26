using Test, Topology, JSON

nested_int(value)=[Int.(collect(component)) for component in value]

@testset "MATLAB vg2g exact-output fixtures" begin
    fixture=JSON.parsefile(joinpath(@__DIR__,"fixtures","vg2g_matlab_v1.json"))
    @test fixture["schema_version"]==1
    @test length(fixture["cases"])==5
    for case in fixture["cases"]
        input=nested_int(case["input_rgauss"])
        orientation=Int.(case["input_orientation"])
        source=deepcopy(input)
        gauss,planar_orientation=planarize_virtual_gauss(input,orientation)
        @test input==source
        @test gauss==nested_int(case["gauss"])
        @test planar_orientation==Int.(case["orientation"])

        trace=IOBuffer()
        debug_gauss,debug_orientation=planarize_virtual_gauss(
            input,orientation;DEBUG_MODE=true,debug_io=trace)
        trace_text=String(take!(trace))
        @test debug_gauss==gauss
        @test debug_orientation==planar_orientation
        @test occursin("face incidences",trace_text)
        @test occursin("phase=complete",trace_text)
    end
end

@testset "MATLAB splice-boundary translation" begin
    for (prefix,expected) in [
        (0,[9,1,2,3]),
        (2,[1,2,9,3]),
        (3,[1,2,3,9]),
    ]
        values=[1,2,3]
        Topology._insert_after_prefix!(values,prefix,9;context="boundary regression")
        @test values==expected
    end
    @test_throws ArgumentError Topology._insert_after_prefix!([1,2,3],-1,9)
    @test_throws ArgumentError Topology._insert_after_prefix!([1,2,3],4,9)
end

@testset "DevLinkVG2Gscr input corpus" begin
    cases=[
        (name="unknot",rgauss=[Int[]],orientation=Int[]),
        (name="abalone",rgauss=[[1,-1]],orientation=[1]),
        (name="trefoil",rgauss=[[1,-2,3,-1,2,-3]],orientation=[1,1,1]),
        (name="hopf",rgauss=[[1,-2],[-1,2]],orientation=[1,1]),
        (name="virtual_hopf",rgauss=[[1],[-1]],orientation=[1]),
        (name="mirror_trefoil",rgauss=[[1,-2,3,-1,2,-3]],orientation=[-1,-1,-1]),
        (name="koda",rgauss=[[-1,-2,-8,3,9,1,-4,6,5,4,-6,-5,-7,8,-3,7,2,-9]],
         orientation=[-1,1,1,1,-1,1,1,1,1]),
        (name="koda_virtual",rgauss=[[-1,-2,3,1,-4,5,4,-5,-3,2]],
         orientation=[-1,1,1,1,-1]),
        (name="virtual_trefoil",rgauss=[[-1,2,1,-2]],orientation=[1,1]),
        (name="multi_component",rgauss=[[3,2],[1,-3],[-1,-2]],orientation=[-1,1,1]),
    ]
    for case in cases
        source=deepcopy(case.rgauss)
        trace=IOBuffer()
        gauss,orientation=planarize_virtual_gauss(
            case.rgauss,case.orientation;DEBUG_MODE=true,debug_io=trace)
        @test case.rgauss==source
        @test orientation[1:length(case.orientation)]==case.orientation
        @test all(==(0),orientation[length(case.orientation)+1:end])
        expected=[[entry for entry in component if entry!=0] for component in case.rgauss]
        preserved=[[entry for entry in component
                    if entry!=0 && abs(entry)<=length(case.orientation)] for component in gauss]
        @test preserved==expected
        labels=abs.(filter(!=(0),reduce(vcat,gauss;init=Int[])))
        @test all(label -> count(==(label),labels)==2,unique(labels))
        @test occursin("[vg2g]",String(take!(trace)))
    end
end
