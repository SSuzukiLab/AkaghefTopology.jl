module MarkedValidityTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
const ROOT=normpath(joinpath(@__DIR__,"../.."))
loadcase(name)=JSON3.read(read(joinpath(ROOT,"vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})

# Source: marked_topology.jl::validate/marked_topology_report and existing
# marked_topology.jl tests. Each negative probes a different public boundary;
# no diagram is repaired or reoriented to make conversion applicable.
@testset "marked validity and operation preconditions" begin
    torus=loadcase("torus_mer_lon");before=deepcopy(torus)
    @test H.is_valid(torus)
    @test H.validate(torus;level=:surface).valid
    report=H.marked_topology_report(torus)
    @test report["diagram_valid"] && report["closed_heegaard"]
    @test !report["conversion"]["applicable"] # no stored D
    @test "one_stored_disk" in report["conversion"]["reasons"]
    @test torus==before

    for (mutate,code) in (
        (d->(d["surface"]["genus"]=2),:SURFACE_GENUS),
        (d->(d["surface"]["marking_id"]=""),:MARKING_REQUIRED),
        (d->(d["crossings"][1]["sign"]*=-1),:CROSSING_SIGN),
        (d->(empty!(d["crossings"])),:UNDECLARED_INTERSECTION),
        (d->(d["curves"][2]["role"]="alpha"),:SAME_FAMILY_INTERSECTION))
        bad=deepcopy(torus);mutate(bad)
        r=H.validate(bad)
        @test !r.valid && only(r.issues).code==code
        @test !H.marked_topology_report(bad)["diagram_valid"]
    end

    # Actual positive-genus compression is different from malformed embedding.
    contractible=B.grid_document(6)
    B.add_curve!(contractible,"a","alpha",true,[[(1,1),(5,1),(5,5),(1,5),(1,1)]])
    B.add_curve!(contractible,"b","beta",true,[[(2,2),(4,2),(4,4),(2,4),(2,2)]])
    @test H.validate(contractible;level=:surface).valid
    @test !H.is_valid(contractible)
    r=H.marked_topology_report(contractible)
    @test r["diagram_valid"] && !r["closed_heegaard"]
    @test sort(r["compression"].alpha_boundary_genera)==[0,1]
    @test "complete_cut_systems" in r["conversion"]["reasons"]

    # Proper arcs are valid embedded objects, but not closed compression curves.
    arc=B.grid_document(4;surface=:disk)
    B.add_curve!(arc,"a","alpha",false,[[(2,0),(2,4)]])
    @test H.validate(arc;level=:surface).valid
    @test !H.is_valid(arc)
    ar=H.marked_topology_report(arc)
    @test ar["diagram_valid"] && !ar["compression_applicable"]
    @test "closed_curves" in ar["conversion"]["reasons"]

    # A false D claim carries both its kernel code and failed disk conditions.
    bad=deepcopy(contractible)
    bad["region_marks"]=[Dict("id"=>"D","seed"=>["P",0],"claim"=>"embedded_closed_disk")]
    invalid=H.validate(bad;level=:surface)
    @test !invalid.valid && only(invalid.issues).code==:REGION_CLAIM_FAILED
    @test occursin("embedded_closed_disk",only(invalid.issues).message)
    koda=JSON3.read(read(joinpath(ROOT,"vendor/koda_genus3/data/figure6.fullmarked.json"),String),Dict{String,Any})
    @test H.marked_topology_report(koda)["conversion"]["applicable"]
    @test H.is_valid(koda)

    # Source: marked_topology.jl::_marked_states reserved metadata contract.
    # Source metadata otherwise stays opaque, including arrays and decimal
    # values nested under unrelated fields. Invalid geometric metadata must
    # fail before source_manifest/computation or move consumers dereference it.
    for metadata in (Any[],7,nothing)
        bad=deepcopy(torus);bad["metadata"]=metadata
        invalid=H.validate(bad;level=:surface)
        @test !invalid.valid && only(invalid.issues).code==:METADATA_OBJECT
        @test !H.marked_topology_report(bad)["diagram_valid"]
    end
    for display in (Any[],7,nothing)
        bad=deepcopy(torus);bad["metadata"]=Dict("display"=>display)
        invalid=H.validate(bad;level=:surface)
        @test !invalid.valid && only(invalid.issues).code==:DISPLAY_OBJECT
    end
    for infinity in (nothing,Any[],["P"],["P",true],["P",-1],["P",0.5],[7,0])
        bad=deepcopy(torus);bad["metadata"]=Dict("display"=>Dict("infinity_triangle"=>infinity))
        invalid=H.validate(bad;level=:surface)
        @test !invalid.valid && only(invalid.issues).code==:INFINITY_TRIANGLE
    end
    for infinity in (["absent",0],["P",100000])
        bad=deepcopy(torus);bad["metadata"]=Dict("display"=>Dict("infinity_triangle"=>infinity))
        invalid=H.validate(bad;level=:surface)
        @test !invalid.valid && only(invalid.issues).code==:INFINITY_TRIANGLE_MISSING
    end
    opaque=deepcopy(torus)
    opaque["metadata"]=Dict("source"=>Any[0.25,nothing,Dict("custom"=>true)],"display"=>Dict("infinity_triangle"=>["P",0]))
    @test H.validate(opaque;level=:surface).valid
    delete!(opaque["metadata"]["display"],"infinity_triangle")
    @test H.validate(opaque;level=:surface).valid
end
end
