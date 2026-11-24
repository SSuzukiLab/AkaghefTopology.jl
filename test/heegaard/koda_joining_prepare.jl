module KodaJoiningPrepareTests
using Test, Topology
const H = Topology.HeegaardCodes
const ROOT = abspath(joinpath(@__DIR__, "..", "..", "..", ".."))

function _case(name)
    H.read_hc(joinpath(ROOT, "threads", "PJ260922", "data", "261114-usecases", name, "diagram.hc"))
end

function _square_match(left, right, modulus)
    a = mod(numerator(left), modulus)
    b = mod(numerator(right), modulus)
    any(u -> gcd(u, modulus) == 1 && mod(u * u * a, modulus) == b, 1:modulus-1)
end

@testset "prepare a joining disk on the three judgment diagrams" begin
    for name in ("g1-D-5-2", "stab-D-5-2", "g2-S3-pairs")
        source = _case(name)
        prepared, witnesses = H.hc_prepare_joining_disk(source)
        @test !isempty(H.enumerate_joining_disk_candidates(H.forget_reference(H.realize_hc(prepared))))
        before = H.hc_first_homology(source)
        after = H.hc_first_homology(prepared)
        @test before.free_rank == after.free_rank
        @test before.invariant_factors == after.invariant_factors
        left = H.hc_linking_form(source)
        right = H.hc_linking_form(prepared)
        @test left.free_rank == right.free_rank
        @test left.torsion == right.torsion
        if name == "stab-D-5-2"
            @test left.matrix[1, 1] == 3//5
            @test right.matrix[1, 1] == 2//5
            @test _square_match(3//5, 2//5, 5)
            @test first(witnesses).row_id == "hc.slide.pants"
            @test last(witnesses).row_id == "lc.isotopy.digon.create.ab"
        else
            @test left.matrix == right.matrix
            @test all(witness.row_id == "lc.isotopy.digon.create.ab" for witness in witnesses)
        end
        @test H.hc_reidemeister_torsion(source).values == H.hc_reidemeister_torsion(prepared).values
    end
end

end
