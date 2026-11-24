module KodaCmr3Tests
using Test, Topology, JSON
const H = Topology.HeegaardCodes
const ROOT = abspath(joinpath(@__DIR__, "..", "..", "..", ".."))
const EXPECTED = joinpath(ROOT, "threads", "PJ260922", "data", "261114-conversion", "cmr3-expected.json")

matrix(rows) = [Int(rows[i][j]) for i in eachindex(rows), j in eachindex(rows[1])]
column(values) = [Int(value) for value in values]

function basis(record)
    H.CellularBasis(String.(record["vertices"]), String.(record["edges"]), String.(record["faces"]))
end

function vertex_map(record)
    Dict{Int,Union{Int,Nothing}}(parse(Int, key) => value === nothing ? nothing : Int(value)
                                 for (key, value) in record)
end

function edge_pairs(record)
    [String(pair[1]) => String(pair[2]) for pair in record]
end

function cellular_witness(document, phi1)
    cell = document["cellular"]
    H.CellularMapWitness(
        basis(cell["source_complex_basis"]),
        basis(cell["spine_complex_basis"]),
        matrix(cell["phi0"]),
        matrix(phi1),
        matrix(cell["phi2"]),
        matrix(cell["source_boundary_1"]),
        matrix(cell["source_boundary_2"]),
        matrix(cell["spine_boundary_1"]),
        matrix(cell["spine_boundary_2"]),
        vertex_map(cell["original_to_smoothed_vertex"]),
        vertex_map(cell["smoothed_to_gauss_vertex"]),
        edge_pairs(cell["ograph_edge_to_spine_edge"]),
        String(cell["conventions"]),
        String(cell["source_hash"]),
        String(cell["target_hash"]))
end

function diagram_of(case)
    link = VirtualLink()
    set_data!(link; gauss=case["gauss"], orientation=column(case["orientation"]))
    link
end

@testset "CMR3 frozen Euler classes and cellular witnesses" begin
    document = JSON.parsefile(EXPECTED)

    combed = H.cmr3_combed()
    @test combed.kind == Symbol(document["combed"]["expected_kind"])
    @test combed.order === nothing
    @test combed.diagnostic === nothing
    held = H.combing_framing_from_cmr3(combed; b1=0)
    @test held.b1 == 0
    @test held.presentation == :cmr3
    @test held.framing_status == :unjudged
    @test held.second_field === nothing
    @test held.spin_c === nothing
    @test held.p1 === nothing
    @test held.theta === nothing
    @test_throws ArgumentError H.combing_framing_from_cmr3(combed; b1=1)

    for item in document["frozen"]
        delta = matrix(item["delta"])
        c_p = column(item["c_p"])
        witness = H.ModularWitness(Int(item["witness"]["modulus"]), column(item["witness"]["u"]))
        reading = H.classify_euler_class(delta, c_p; witness=witness)
        @test reading.state.kind == Symbol(item["expected_kind"])
        @test reading.state.order == Int(item["expected_order"])
        @test reading.state.diagnostic == Symbol(item["expected_diagnostic"])
        @test reading.witness_annihilates
        @test reading.witness_residue == Int(item["witness"]["expected_residue"])
        @test reading.integer_solution === nothing
        mod2 = H.classify_mod2_solution(delta, c_p)
        @test mod2.kind == Symbol(item["mod2_expected_kind"])
        @test mod2.kind != :framing_extendable
        packed = H.combing_framing_from_cmr3(reading, c_p; b1=0, mod2=mod2)
        @test packed.framing_status == :obstructed
        @test packed.second_field === nothing
        @test packed.euler_representative == c_p
        @test packed.combing.order == Int(item["expected_order"])
        @test packed.mod2_status == :mod2_solution
        @test packed.p1 === nothing
        @test packed.theta === nothing
    end

    good = cellular_witness(document, document["cellular"]["good"]["phi1"])
    broken = cellular_witness(document, document["cellular"]["broken"]["phi1"])
    good_chain = H.cellular_map_commutes(good)
    broken_chain = H.cellular_map_commutes(broken)
    @test good_chain.commutes == document["cellular"]["good"]["commutes"]
    @test good_chain.degree1 && good_chain.degree2
    @test broken_chain.commutes == document["cellular"]["broken"]["commutes"]
    @test !(broken_chain.degree1 && broken_chain.degree2)
    @test H.cellular_cochain_natural(good).natural == document["cellular"]["good"]["cochain_natural"]
    @test H.cellular_cochain_natural(broken).natural == document["cellular"]["broken"]["cochain_natural"]
    @test good.phi1[1, 2] == -1
    @test broken.phi1[1, 2] == 1
    @test good.phi1[3, 4] == 2
    @test count(pair -> pair == ("d" => "g"), good.ograph_edge_to_spine_edge) == 1
    @test count(pair -> pair == ("b" => "e"), good.ograph_edge_to_spine_edge) == 1
    @test [H.gauss_vertex_of(good, id) for id in 1:5] == [1, 2, nothing, nothing, 3]
    @test good.source_hash == "hand-hc-cmr3"
    @test good.target_hash == "hand-spine-cmr3"

    for case in document["mst22"]["cases"]
        disks = disk_table(diagram_of(case))
        delta = matrix(case["delta"])
        c_p = column(case["c_p"])
        @test disks.delta == delta
        @test disks.cp2 == column(case["cp2"])
        reading = H.classify_euler_class(delta, c_p)
        @test reading.state.kind == Symbol(case["expected_kind"])
        if case["integer_solution"]
            @test reading.state.kind == :framing_extendable
            @test reading.integer_solution !== nothing
            mod2 = H.classify_mod2_solution(delta, c_p)
            @test mod2.kind == :mod2_solution
            packed = H.combing_framing_from_cmr3(reading, c_p; b1=0, mod2=mod2)
            @test packed.framing_status == :framed
            @test packed.second_field == reading.integer_solution
            @test packed.p1 === nothing
            @test packed.theta === nothing
            @test packed.mod2_status == :mod2_solution
        else
            @test reading.state.kind == :euler_obstructed
            @test reading.state.order == Int(case["expected_order"])
            @test reading.integer_solution === nothing
            packed = H.combing_framing_from_cmr3(reading, c_p; b1=0)
            @test packed.framing_status == :obstructed
            @test packed.second_field === nothing
            @test packed.p1 === nothing
        end
        if haskey(case, "printed_weights")
            weights = case["printed_weights"]
            vectors = weights[1] isa AbstractVector ? weights : [weights]
            for vector in vectors
                @test delta * column(vector) == -c_p
            end
        end
    end
end

end
