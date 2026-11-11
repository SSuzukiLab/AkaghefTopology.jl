using Test
using Topology

const H = Topology.Heegaard

@testset "one HeegaardDiagram from standard or explicit roles" begin
    mktempdir() do dir
        write(joinpath(dir, "figure.lc"),
            "genus: 1\nl1 : -1\nm1 : 1 2\nb : -2\n")
        write(joinpath(dir, "compact.hc"),
            "format : hc/2\nlc : figure.lc\nstandard\nrole b : beta\n")
        write(joinpath(dir, "explicit.hc"),
            "format : hc/2\nlc : figure.lc\nrole m1 : alpha\nrole b : beta\n")
        a = H.read_hc(joinpath(dir, "compact.hc"))
        b = H.read_hc(joinpath(dir, "explicit.hc"))
        @test a isa Topology.HeegaardDiagram
        @test b isa Topology.HeegaardDiagram
        @test typeof(a) === typeof(b) === H.HeegaardLayer
        @test Topology.HeegaardDiagram === H.HeegaardLayer
        @test !(:standard in fieldnames(typeof(a)))
        @test a.roles == b.roles == Dict("m1" => :alpha, "b" => :beta)
        @test a.lc.words == b.lc.words
        @test a.lc.words["m1"] == [1, 2]
        @test !haskey(a.lc.words, "alpha1")
        @test Topology.is_standard(a) && H.is_standard(b)
        auto = H.hc_text(a; lc_filename="figure.lc")
        explicit = H.hc_text(a; lc_filename="figure.lc", encoding=:explicit)
        @test startswith(auto, "format : hc/2\n")
        @test occursin("\nstandard\n", auto)
        @test !occursin("role m1", auto)
        @test occursin("role m1 : alpha\n", explicit)
        @test !occursin("\nstandard\n", explicit)
        round = joinpath(dir, "round.hc")
        H.write_hc(round, a; lc_filename="figure.lc", encoding=:explicit)
        restored = H.read_hc(round)
        @test restored.roles == a.roles
        @test restored.lc.words["m1"] == a.lc.words["m1"]
        @test H.is_standard(restored)
        bare = joinpath(dir, "bare.hc")
        write(bare, "lc : figure.lc\nstandard\nrole b : beta\n")
        @test_throws H.KodaInputError H.read_hc(bare)
        general = joinpath(dir, "general.hc")
        write(general, "lc : figure.lc\nrole b : alpha\n")
        g = H.read_hc(general)
        @test !H.is_standard(g)
        kept = joinpath(dir, "keep.hc")
        write(kept, "do-not-truncate")
        @test_throws ArgumentError H.write_hc(kept, g; lc_filename="figure.lc", encoding=:standard)
        @test read(kept, String) == "do-not-truncate"
    end
end

@testset "legacy HC is the same HeegaardDiagram" begin
    path = normpath(joinpath(@__DIR__, "..", "..", "..", "..",
        "threads", "PJ260922", "data", "261107-lickorish-data", "group-15", "case-01.hc"))
    hc = H.read_hc(path)
    @test hc isa Topology.HeegaardDiagram
    @test hc.roles == Dict("alpha1" => :alpha, "beta1" => :beta)
    @test !haskey(hc.roles, "m1")
    @test !H.is_standard(hc)
end
