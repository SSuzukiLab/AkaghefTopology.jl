module MarkedSchemaTests
using Test, Topology, JSON3
const H = Topology.Heegaard
const B = H.MarkedDiagramBuilders
const M = H.MarkedSurfaceDiagrams
include(joinpath(@__DIR__, "../../examples/HeegaardInput.jl"))

@testset "lossless marked mathematical JSON" begin
    mktempdir() do directory
        path = joinpath(directory, "source.hgd.json")
        # Existing independent source fixtures include disconnected/zero-crossing
        # components; do not reconstruct them from their intersection projection.
        for name in ("sphere_empty", "nested_zero_crossings", "torus_parallel_zero_crossings", "torus_mer_lon")
            original = JSON3.read(read(joinpath(@__DIR__, "../../vendor/marked_surface_diagrams/examples", name*".json"), String), Dict{String,Any})
            original["metadata"] = Dict("source"=>Dict("figure"=>"α / β", "locator"=>[3, "Fig. 2"]),
                "display"=>Dict("infinity_triangle"=>["P", 0]), "weight"=>0.125, "reviewed"=>false, "pending"=>nothing)
            H.save_heegaard(path, original)
            loaded = H.load_heegaard(path)
            @test loaded isa Dict{String,Any}
            @test loaded == original
            @test M.representation_digest(loaded) == M.representation_digest(original)
            manifest = H.source_manifest(loaded)
            @test manifest == original["metadata"]
            manifest["source"]["figure"] = "changed copy"
            @test loaded["metadata"] == original["metadata"]
        end
        # The existing reference constructor supplies genera independently of
        # alpha/beta counts. Preserve its auxiliary roles and source provenance.
        for genus in 1:4
            original = H.lickorish_reference(genus)
            H.save_heegaard(path, original)
            @test H.load_heegaard(path) == original
        end
        redundant = B.grid_document(6)
        B.add_curve!(redundant, "a1", "alpha", true, [[(0, 1), (6, 1)]])
        B.add_curve!(redundant, "a2", "alpha", true, [[(0, 3), (6, 3)]])
        B.add_curve!(redundant, "b", "beta", true, [[(0, 5), (6, 5)]])
        redundant["patches"][1]["side_parameters"] = Dict(s=>["$(i)/6" for i in 0:6] for s in ("south", "east", "north", "west"))
        H.save_heegaard(path, redundant)
        @test H.load_heegaard(path) == redundant
        # Proper arcs remain valid source input even without a closed splitting.
        arc = B.grid_document(4; surface=:disk)
        B.add_curve!(arc, "a", "alpha", false, [[(0, 2), (4, 2)]])
        H.save_heegaard(path, arc)
        @test H.load_heegaard(path) == arc

        # marked_topology.jl::_marked_states uses alpha U beta for D; the
        # auxiliary cut crosses and splits D only in its display realization.
        split = B.grid_document(8)
        B.add_curve!(split, "a", "alpha", true, [[(2, 2), (6, 2), (6, 6), (2, 6), (2, 2)]])
        B.add_curve!(split, "cut", "aux", true, [[(4, 0), (4, 8)]])
        B.declare_crossing!(split, "p", 4, 2, "a", "cut", 1)
        B.declare_crossing!(split, "q", 4, 6, "a", "cut", -1)
        split["region_marks"] = [Dict("id"=>"D", "seed"=>["P", 54], "claim"=>"embedded_closed_disk")]
        H.save_heegaard(path, split)
        @test H.load_heegaard(path) == split
        @test length(H.prepare_marked_display(H.load_heegaard(path); cut_curves=["cut"])["region_marks"][1]["cut_pieces"]) == 2
        accepted = read(path, String)
        malformed = deepcopy(split); malformed["crossings"][1]["sign"] = -1
        @test_throws M.DiagramError H.save_heegaard(path, malformed)
        @test read(path, String) == accepted
        for (key, value) in (("version", 2), ("unexpected_mathematics", []), ("metadata", []))
            malformed = deepcopy(split); malformed[key] = value
            expected_error = key == "metadata" ? ArgumentError : M.DiagramError
            @test_throws expected_error H.save_heegaard(path, malformed)
            @test read(path, String) == accepted
            write(joinpath(directory, "bad.json"), JSON3.write(malformed))
            @test_throws expected_error H.load_heegaard(joinpath(directory, "bad.json"))
        end
        malformed = deepcopy(split); delete!(malformed["surface"], "marking_id")
        @test_throws M.DiagramError H.save_heegaard(path, malformed)
        @test read(path, String) == accepted
        @test_throws ArgumentError H.save_heegaard(path, Dict("schema"=>"marked-cut-view/1"))
        write(joinpath(directory, "bad.json"), "[]")
        @test_throws ArgumentError H.load_heegaard(joinpath(directory, "bad.json"))

        authored = heegaard_input_example(path)
        @test length(authored["crossings"]) == 1
        @test authored["crossings"][1]["at"] == ["P", "v3_1"]
        @test authored == H.load_heegaard(path)
        # Existing legacy serialization and public return type remain intact.
        legacy = H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)]), [[1]])
        H.save_heegaard(path, legacy)
        @test H.load_heegaard(path) isa H.DecoratedHeegaardCode
        @test H.canonical_hash(H.load_heegaard(path)) == H.canonical_hash(legacy)
    end
end
end
