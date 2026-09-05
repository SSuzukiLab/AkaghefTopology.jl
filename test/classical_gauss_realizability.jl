using Test

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_realizability.jl"))

function _check_ribbon_certificate(c, n)
    @test c.vertices == n
    @test c.edges == 2n
    n == 0 || @test c.faces == length(c.face_cycles)
    @test c.genus == (2 - (c.vertices - c.edges + c.faces)) ÷ 2
    n == 0 && return
    @test sort(reduce(vcat, c.face_cycles)) == collect(1:4n)
    @test all(order -> length(unique(order)) == 4, c.rotation)
end

@testset "input validation and ribbon-map certificates" begin
    @test_throws ArgumentError unsigned_gauss_certificate([1, 2, 1])
    @test_throws ArgumentError unsigned_gauss_certificate([1, 2, 1, 3])
    @test_throws ArgumentError unsigned_gauss_certificate([1, 2, 1, 2, 3, 4])
    @test_throws DimensionMismatch classical_gauss_certificate([1, -1], Int[])
    @test_throws ArgumentError classical_gauss_certificate([1, 1], [1])

    empty_unsigned = unsigned_gauss_certificate(Int[])
    _check_ribbon_certificate(empty_unsigned, 0)
    @test empty_unsigned.realizable
    @test empty_unsigned.faces == 2

    for word in ([1, 1], [1, 2, 2, 1], [1, 2, 3, 1, 2, 3])
        certificate = unsigned_gauss_certificate(word)
        _check_ribbon_certificate(certificate, length(word) ÷ 2)
        @test certificate.realizable
        @test construct_classical_immersion(word).genus == 0
    end
    for word in ([1, 2, 1, 2], [1, 2, 3, 4, 1, 2, 3, 4])
        certificate = unsigned_gauss_certificate(word)
        _check_ribbon_certificate(certificate, length(word) ÷ 2)
        @test !certificate.realizable
        @test_throws ArgumentError construct_classical_immersion(word)
    end
end

@testset "Kurlin Carter-surface signed convention" begin
    trefoil = KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, 1, 1])
    trefoil_certificate = classical_gauss_certificate(trefoil)
    _check_ribbon_certificate(trefoil_certificate, 3)
    @test trefoil_certificate.route == :carter_surface
    @test trefoil_certificate.realizable
    @test trefoil_certificate.faces == 5
    @test classical_gauss_realizable(trefoil)

    # This is a one-component fixture whose cyclic letter sequence is the
    # concatenation [1,-2,-1,2] of Kurlin's two words {12+, 1-2}; it is not
    # the multi-component paragraph itself.  With under-head exponents
    # [-1,+1], this one-component Carter surface is a torus, not a sphere.
    kurlin_torus = KnotGaussDiagram([1, -2, -1, 2], [-1, 1])
    torus_certificate = classical_gauss_certificate(kurlin_torus)
    _check_ribbon_certificate(torus_certificate, 2)
    @test torus_certificate.faces == 2
    @test torus_certificate.genus == 1
    @test !torus_certificate.realizable

    # The unsigned shadow 123123 is realizable, but an arbitrary independent
    # orientation vector need not be compatible with the chosen over/under
    # visits.  This is precisely why the signed Carter test is not replaced by
    # an unsigned-shadow predicate.
    virtual_decoration = KnotGaussDiagram([1, 2, 3, -1, -2, -3], [1, 1, 1])
    @test unsigned_gauss_realizable(abs.(virtual_decoration.gauss))
    @test !classical_gauss_realizable(virtual_decoration)

    # Every over/under choice of one fixed planar immersion is classical when
    # the crossing signs are attached consistently with its local rotation.
    immersion = construct_classical_immersion([1, 2, 3, 1, 2, 3])
    decorations = enumerate_classical_decorations(immersion.word;
                                                   immersion=immersion)
    @test length(decorations) == 2^3
    @test all(classical_gauss_realizable, decorations)
end

@testset "Carter rotation reversal regression" begin
    fixtures = (([1, -2, 3, -1, 2, -3], [1, 1, 1]),
                ([1, -2, -1, 2], [-1, 1]),
                ([1, 2, 3, -1, -2, -3], [1, 1, 1]))
    for (code, orientation) in fixtures
        certificate = classical_gauss_certificate(code, orientation)
        reversed_rotation = [
            (order[1], order[4], order[3], order[2])
            for order in certificate.rotation
        ]
        reversed_faces = _cg_face_cycles(reversed_rotation, certificate.vertices)
        reversed_genus = _cg_genus(certificate.vertices, certificate.edges,
                                   length(reversed_faces))
        @test length(reversed_faces) == certificate.faces
        @test reversed_genus == certificate.genus
        @test (reversed_genus == 0) == certificate.realizable
        # Raw face-cycle orientation/order is convention-dependent and is not
        # a database-persisted identity; only face count/genus are compared.
    end
end

@testset "repository identity convention" begin
    d = KnotGaussDiagram([1, -2, 3, -1, 2, -3], [1, -1, 1])
    @test canonical_serialization(d) == canonical_serialization(cyclic_shift(d, 2))
    @test classical_gauss_certificate(d).realizable ==
          classical_gauss_certificate(cyclic_shift(d, 2)).realizable
    @test unsigned_gauss_certificate([1, 2, 3, 1, 2, 3]).genus ==
          unsigned_gauss_certificate([3, 1, 2, 3, 1, 2]).genus
end

@testset "fast Carter decision agrees with certificates" begin
    # This covers every rooted normalized signed diagram through n=3 while
    # keeping the focused test independent of the N=5 campaign runtime.
    for n in 0:3
        for diagram in enumerate_knot_gauss(n; orientations=true)
            @test classical_gauss_realizable_fast(diagram) ==
                  classical_gauss_realizable(diagram)
        end
    end
end

# Existing virtual planarisation is intentionally used only as an independent
# smoke cross-check: zero inserted virtual crossings is evidence of agreement
# on these examples, not a replacement for the Carter theorem decision.
if !isdefined(Main, :Topology)
    include(joinpath(@__DIR__, "..", "src", "Topology.jl"))
end

@testset "independent virtual-planarisation cross-check" begin
    classical_code = [1, -2, 3, -1, 2, -3]
    classical_output = Topology.planarize_virtual_gauss([classical_code], [1, 1, 1])
    @test count(==(0), classical_output[2]) == 0

    virtual_code = [1, 2, 3, -1, -2, -3]
    virtual_output = Topology.planarize_virtual_gauss([virtual_code], [1, 1, 1])
    @test count(==(0), virtual_output[2]) > 0
end
