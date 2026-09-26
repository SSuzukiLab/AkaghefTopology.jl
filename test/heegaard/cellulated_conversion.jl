module CellulatedConversionTests
using Topology, Test, JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const ROOT=normpath(joinpath(@__DIR__,"../.."))
@testset "Marked genus-3 to existing Koda conversion" begin
    doc=JSON3.read(read(joinpath(ROOT,"vendor/koda_genus3/data/figure6.heegaard.marked.json"),String),Dict{String,Any})
    original=deepcopy(doc)
    p=H.cellulated_punctured(doc)
    @test doc==original
    @test length(p.provenance["region_map"])==7
    b=H.to_odata(p.object)
    @test (length(b.spine.vertex_signs),length(b.spine.edges),length(b.spine.complex.regions))==(5,10,6)
    @test b.spine.complex.checks["oriented_branching"]
    source=JSON3.read(read(joinpath(ROOT,"vendor/koda_genus3/data/figure6.odata.json"),String))
    @test H._cyclic_equal(only(b.ograph.gauss),reverse(Int.(only(source[1]))))
    @test b.ograph.orientation==Int.(source[2])
    v=VirtualLink();set_data!(v;rgauss=b.ograph.gauss,orientation=b.ograph.orientation)
    @test is_closed(v)==(true,(true,true,true))
    @test size(disk_table(v).delta,1)==6
    @test startswith(b.certificate.source_hash,"labelled-sha256:")
    @test H.compare_representations(p.object.diagram,b.ograph)["source_matches"]
    # The original target differs even after the same convention transport.
    target=JSON3.read(read(joinpath(ROOT,"vendor/koda_genus3/data/target_odata.json"),String))
    @test !H._cyclic_equal(only(b.ograph.gauss),reverse(Int.(only(target[1]))))
    # Source mutation is rejected rather than fitted to the expected output.
    bad=deepcopy(doc);only(x for x in bad["crossings"] if x["id"]=="x2")["sign"]*=-1
    @test_throws M.DiagramError H.cellulated_punctured(bad)
    bad=deepcopy(doc);empty!(bad["region_marks"])
    @test_throws ArgumentError H.cellulated_punctured(bad)
    bad=deepcopy(doc);push!(bad["region_marks"],deepcopy(only(bad["region_marks"])))
    @test_throws Exception H.cellulated_punctured(bad)
    # Existing geometric-only D is not automatically reoriented.
    bad=M.TOML.parsefile(joinpath(ROOT,"vendor/marked_surface_diagrams/examples/torus_three_crossings.toml"))
    bad["region_marks"]=[Dict("id"=>"D","seed"=>["P",68],"claim"=>"embedded_closed_disk")]
    @test_throws ArgumentError H.cellulated_punctured(bad)
    println("Koda genus3 repository o-data: ",(b.ograph.gauss,b.ograph.orientation))
end
end
