module JoiningDiskTests
using Test, Topology, JSON3
const H=Topology.Heegaard
include(joinpath(@__DIR__,"../../examples/HeegaardAbalone.jl"))

@testset "Joining disk and abalone smoke" begin
    ex=HeegaardAbaloneExample.fixture(); h,d=ex.diagram,ex.disk
    @test H.surface_genus(h)==1
    @test H.intersection_matrix(h)==reshape([1],1,1)
    @test d.boundary_darts==[1,8]
    @test H.validate_joining_disk(h,d)
    @test H.validate_joining_disk(h,H.JoiningDisk(d.labelled_diagram_hash,d.region,[8,1]))
    @test_throws ArgumentError H.make_punctured(h)
    @test_throws ArgumentError H.to_odata(h)
    @test_throws ArgumentError H.joining_disk(h,2) # repeated boundary vertices
    @test_throws ArgumentError H.joining_disk(h,0)
    @test_throws ArgumentError H.validate_joining_disk(h,H.JoiningDisk(d.labelled_diagram_hash,d.region,[1,7]))
    @test_throws ArgumentError H.validate_joining_disk(h,H.JoiningDisk("stale",d.region,d.boundary_darts))
    minimal=H.certify_heegaard(H.decorate(H.heegaard_code([[1]],[[1]],
        [H.Intersection(1,1,1,1)]);cellular=true))
    @test_throws ArgumentError H.joining_disk(minimal,1) # disk closure not embedded
    ph=H.make_punctured(h;joining_disk=d)
    b=H.to_odata(ph)
    @test b.ograph.gauss==[[1,-1]]
    @test b.ograph.orientation==[1]
    @test b.certificate.label_map["smoothed_intersections"]==[1,2]
    @test b.certificate.label_map["remaining_intersections"]==[3]
    @test [e.source_darts for e in b.spine.edges]==[[5,12],[10,3]]
    @test H.to_heegaard(b.spine)===h # retained provenance, not an inverse algorithm
    @test b.comparison["manifold_comparison"]=="not performed"
    @test !haskey(b.ograph.metadata,"h1")
    link=VirtualLink(); set_data!(link;rgauss=b.ograph.gauss,orientation=b.ograph.orientation)
    @test real_gauss_code(link)==([[1,-1]],[1])
    @test is_closed(link)==(true,(true,true,true))
    # The second embedded joining disk gives the same one-vertex presentation.
    other=H.to_odata(h;joining_disk=H.joining_disk(h,3)).ograph
    @test (other.gauss,other.orientation)==(b.ograph.gauss,b.ograph.orientation)
    # Reverse an input parameter orientation; disk-induced directions compensate.
    reversed=H.certify_heegaard(H.decorate(H.heegaard_code([[3,2,1]],[[1,2,3]],
        [H.Intersection(i,1,1,s) for (i,s) in enumerate([-1,1,-1])]);cellular=true))
    r=only([i for (i,cy) in enumerate(H.boundary_cycles(reversed))
        if Set(H.ribbon_graph(reversed).darts[j].origin for j in cy)==Set([1,2]) && length(cy)==2])
    @test H.to_odata(reversed;joining_disk=H.joining_disk(reversed,r)).ograph.orientation==[1]
    @test_throws ArgumentError H.make_punctured(reversed;joining_disk=d)
    # Reflect the oriented surface instead: the vertex sign must change.
    reflected=H.certify_heegaard(H.decorate(H.heegaard_code([[1,2,3]],[[1,2,3]],
        [H.Intersection(i,1,1,s) for (i,s) in enumerate([-1,1,-1])]);cellular=true))
    r=only([i for (i,cy) in enumerate(H.boundary_cycles(reflected))
        if Set(H.ribbon_graph(reflected).darts[j].origin for j in cy)==Set([1,2]) && length(cy)==2])
    @test H.to_odata(reflected;joining_disk=H.joining_disk(reflected,r)).ograph.orientation==[-1]
    mktempdir() do dir
        path=joinpath(dir,"example.punctured.hgd.json")
        H.save_punctured(path,ph); loaded=H.load_punctured(path)
        @test loaded.joining_disk.boundary_darts==d.boundary_darts
        @test H.canonical_hash(loaded.diagram)==H.canonical_hash(h)
        @test H.to_odata(loaded).ograph.gauss==b.ograph.gauss
        @test H.to_odata(loaded).ograph.orientation==b.ograph.orientation
        value=JSON3.read(read(path,String),Dict{String,Any})
        value["joining_disk"]["boundary_darts"]=[1,7]
        write(path,JSON3.write(value))
        @test_throws ArgumentError H.load_punctured(path)
        core=joinpath(dir,"legacy.hgd.json")
        H.save_heegaard(core,h)
        @test H.canonical_hash(H.load_heegaard(core))==H.canonical_hash(h)
    end
    svg=HeegaardAbaloneExample.render_svg(ex)
    @test occursin("data-joining-disk=\"1\"",svg)
    @test length(collect(eachmatch(r"data-intersection-id=",svg)))==3
end
end
