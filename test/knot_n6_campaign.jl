using Test
using JSON

include(joinpath(@__DIR__, "..", "scripts", "run_knot_n6_campaign.jl"))

@testset "KnotInfo Gauss parser and official seeds" begin
    @test parse_gauss_notation("[]") == Int[]
    @test parse_gauss_notation("[1, -2, 3, -1, 2, -3]") == [1, -2, 3, -1, 2, -3]
    @test_throws ArgumentError parse_gauss_notation("1, -1")
    @test_throws ArgumentError parse_gauss_notation("[1, nope]")

    rows = load_official_rows()
    @test [row.name for row in rows] == ["0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3"]
    @test [row.crossing_number for row in rows] == [0, 3, 4, 5, 5, 6, 6, 6]
    @test all(row -> KnotGaussDiagram(row.gauss, ones(Int, row.crossing_number)) isa KnotGaussDiagram,
              rows)
    @test sum(length(decorate_official_row(row)) for row in rows) == 281
    @test all(row -> length(decorate_official_row(row)) == 2^row.crossing_number, rows)
end

@testset "Streaming rooted census" begin
    first_run = stream_census(n_min=0, n_max=2)
    second_run = stream_census(n_min=0, n_max=2)
    @test [record["visited_count"] for record in first_run] == [1, 2, 12]
    @test [record["visited_count"] for record in first_run] ==
        [record["visited_count"] for record in second_run]
    @test [record["stream_checksum_sha256"] for record in first_run] ==
        [record["stream_checksum_sha256"] for record in second_run]
    @test all(record -> record["elapsed_seconds"] > 0 && record["codes_per_second"] > 0,
              vcat(first_run, second_run))
    @test [record["stream_checksum_sha256"] for record in first_run] == [
        "bc5fb9ef58674d407bd16f150ad13c0666b331083ce8d9f842f4689be4493b07",
        "360af86de85fa4ff1a795197264f05646ad39df76ac0b780c5bda023c119b188",
        "7d04a0bb275ca7b8ec3e545adfd6974178d5ad8121684737cd09e5fd9cff822a",
    ]
end

@testset "Campaign schema and claims" begin
    path = joinpath(@__DIR__, "..", "..", "..", "threads", "PJ260830",
                    "data", "knot-n6-campaign.json")
    @test isfile(path)
    payload = JSON.parsefile(path)
    @test payload["schema_version"] == 1
    @test payload["canonicalization"]["version"] == 1
    @test payload["input"]["snapshot_sha256"] ==
        "5d4326898fcbe2f59c8717e897a3ac0412818287ba485c76ca0535a81843e7c8"
    @test payload["input"]["source_snapshot_sha256"] ==
        "578be92fc1d37b1351c8273c3beab00654913ea009ca751a71e99bac2209751f"
    @test payload["scope"]["excluded"] == ["Area 2.5", "Area 3", "classical knot/isotopy census", "Reidemeister connectivity run"]
    @test occursin("representation identity only", payload["claims"]["canonical_identity"])
    @test occursin("not isotopy proof", payload["claims"]["knotinfo_role"])
    @test occursin("independently verified Reidemeister", payload["claims"]["isotopy_certificate"])
    @test occursin("deferred", payload["claims"]["move_verifier"])
    @test length(payload["official_seed_decorations"]) == 8
    @test all(seed -> seed["decoration_count"] == 2^seed["crossing_number"],
              payload["official_seed_decorations"])
    census = payload["rooted_signed_gc_census"]
    @test length(census["records"]) == 7
    @test all(record -> record["visited_count"] == record["expected_count"] &&
                     record["elapsed_seconds"] > 0 && record["codes_per_second"] > 0,
              census["records"])
    @test [record["stream_checksum_sha256"] for record in census["records"]] == [
        "bc5fb9ef58674d407bd16f150ad13c0666b331083ce8d9f842f4689be4493b07",
        "360af86de85fa4ff1a795197264f05646ad39df76ac0b780c5bda023c119b188",
        "7d04a0bb275ca7b8ec3e545adfd6974178d5ad8121684737cd09e5fd9cff822a",
        "07d408c8ed571280eefa1bd150a74cfc98c442b41b4ba7bc3f2513050812205b",
        "a0f11fec866baaa7a3618806ee07eae938ffdfce8a1fc119122a8eef9e785cd8",
        "3ad41249bf99f4ca7c4e2da298f7169201d1f661a93c1923ff23251e0868c86c",
        "84b3a985db5b759d7779fcd1448646fbe7f1f491d457863399c690a6f6829771",
    ]
    @test census["total_visited_count"] == 697335
    @test census["total_elapsed_seconds"] > 0 && census["total_codes_per_second"] > 0
    @test occursin("machine-specific", census["timing_note"])
end
