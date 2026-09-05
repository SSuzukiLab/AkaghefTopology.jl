using Test
using SHA

if Base.find_package("JSON") === nothing
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using JSON

const INVARIANT_SRC = joinpath(@__DIR__, "..", "src", "knot_classical_invariants.jl")
include(INVARIANT_SRC)
include(joinpath(@__DIR__, "..", "src", "classical_gauss_realizability.jl"))

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SNAPSHOT_PATH = joinpath(REPO_ROOT, "threads", "PJ260830", "data",
                               "knotinfo-n6.json")
const TARGET_NAMES = ("0_1", "3_1", "4_1", "5_1", "5_2", "6_1", "6_2", "6_3")
const SOURCE_SNAPSHOT_SHA256 =
    "578be92fc1d37b1351c8273c3beab00654913ea009ca751a71e99bac2209751f"
const SNAPSHOT_FILE_SHA256 =
    "2fee6d3ac73e1aeeb8eb7e9b9c9a299623f91749568118818da57cb97a3fdfc7"

terms(pairs...) = Dict(Int(exponent) => BigInt(coefficient)
                        for (exponent, coefficient) in pairs)

const EXPECTED_ALEXANDER = Dict(
    "0_1" => terms(0 => 1),
    "3_1" => terms(0 => 1, 1 => -1, 2 => 1),
    "4_1" => terms(0 => 1, 1 => -3, 2 => 1),
    "5_1" => terms(0 => 1, 1 => -1, 2 => 1, 3 => -1, 4 => 1),
    "5_2" => terms(0 => 2, 1 => -3, 2 => 2),
    "6_1" => terms(0 => 2, 1 => -5, 2 => 2),
    "6_2" => terms(0 => 1, 1 => -3, 2 => 3, 3 => -3, 4 => 1),
    "6_3" => terms(0 => 1, 1 => -3, 2 => 5, 3 => -3, 4 => 1),
)

const EXPECTED_JONES = Dict(
    "0_1" => terms(0 => 1),
    "3_1" => terms(1 => 1, 3 => 1, 4 => -1),
    "4_1" => terms(-2 => 1, -1 => -1, 0 => 1, 1 => -1, 2 => 1),
    "5_1" => terms(2 => 1, 4 => 1, 5 => -1, 6 => 1, 7 => -1),
    "5_2" => terms(1 => 1, 2 => -1, 3 => 2, 4 => -1, 5 => 1, 6 => -1),
    "6_1" => terms(-2 => 1, -1 => -1, 0 => 2, 1 => -2,
                    2 => 1, 3 => -1, 4 => 1),
    "6_2" => terms(-1 => 1, 0 => -1, 1 => 2, 2 => -2,
                    3 => 2, 4 => -2, 5 => 1),
    "6_3" => terms(-3 => -1, -2 => 2, -1 => -2, 0 => 3,
                    1 => -2, 2 => 2, 3 => -1),
)

function _snapshot_terms(text::AbstractString)
    source = replace(strip(String(text)), r"\s+" => "")
    source = replace(source, "*" => "")
    isempty(source) && return Dict{Int,BigInt}()
    startswith(source, "+") || startswith(source, "-") || (source = "+" * source)
    terms = Dict{Int,BigInt}()
    bytes = collect(codeunits(source))
    position = 1
    while position <= length(bytes)
        coefficient_sign = bytes[position] == UInt8('+') ? big(1) : big(-1)
        position += 1
        start = position
        parenthesis_depth = 0
        while position <= length(bytes)
            byte = bytes[position]
            byte == UInt8('(') && (parenthesis_depth += 1)
            byte == UInt8(')') && (parenthesis_depth -= 1)
            parenthesis_depth == 0 &&
                (byte == UInt8('+') || byte == UInt8('-')) && break
            position += 1
        end
        token = replace(String(bytes[start:position - 1]), "(" => "", ")" => "")
        if occursin("t", token)
            pieces = split(token, "^"; limit=2)
            coefficient_text = replace(pieces[1], "t" => "")
            exponent_text = length(pieces) == 2 ? pieces[2] : "1"
            exponent = parse(Int, exponent_text)
            coefficient = isempty(coefficient_text) ? big(1) : parse(BigInt, coefficient_text)
        else
            exponent = 0
            coefficient = parse(BigInt, token)
        end
        terms[exponent] = get(terms, exponent, big(0)) + coefficient_sign * coefficient
        terms[exponent] == 0 && delete!(terms, exponent)
    end
    terms
end

_field(row, headers, name) = row["values"][findfirst(==(name), headers)]
_file_sha256(path) = bytes2hex(SHA.sha256(read(path)))

function _load_snapshot(path)
    payload = JSON.parsefile(path)
    sheet = first(payload["sheets"])
    headers = String.(sheet["headers"])
    name_index = findfirst(==("name"), headers)
    rows = Dict{String,Any}()
    for row in sheet["rows"]
        rows[String(row["values"][name_index])] = row
    end
    payload, sheet, headers, rows
end

@testset "local KnotInfo N<=6 PD signed-Gauss invariant contract" begin
    @test isfile(SNAPSHOT_PATH)
    payload, sheet, headers, rows = _load_snapshot(SNAPSHOT_PATH)
    @test payload["schema_version"] == 1
    @test payload["included_row_count"] == length(TARGET_NAMES)
    @test Tuple(String.(payload["included_knot_names"])) == TARGET_NAMES
    @test String(payload["source"]["snapshot_sha256"]) == SOURCE_SNAPSHOT_SHA256
    @test _file_sha256(SNAPSHOT_PATH) == SNAPSHOT_FILE_SHA256
    @test headers[27] == "pd_notation"
    @test headers[29] == "crossing_number"
    @test headers[63] == "alexander_polynomial"
    @test headers[67] == "jones_polynomial"
    @test headers[127] == "symmetry_type"
    @test headers[129] == "full_symmetry_group"

    @test [Int(rows[name]["source_row_index"]) for name in TARGET_NAMES] == collect(2:9)
    @test sort(collect(keys(rows))) == sort(collect(TARGET_NAMES))

    diagrams = Dict{String,KnotGaussDiagram}()
    for name in TARGET_NAMES
        row = rows[name]
        local_crossings = Int(row["crossing_number"])
        raw_pd = String(_field(row, headers, "pd_notation"))
        pd_rows = isempty(strip(raw_pd)) ? Any[] : JSON.parse(raw_pd)
        diagram = pd_to_gauss(pd_rows)
        diagrams[name] = diagram

        pd = pd_from_gauss(diagram)
        certificate = classical_gauss_certificate(diagram)
        @test length(pd_rows) == local_crossings
        @test crossing_count(diagram) == local_crossings
        @test verify_pd(pd)
        @test certificate.route == :carter_surface
        @test certificate.realizable
        @test certificate.genus == 0
        @test certificate.faces == local_crossings + 2
        @test classical_gauss_realizable_fast(diagram)

        computed_alexander = alexander_polynomial(diagram)
        computed_jones = jones_polynomial(diagram)
        @test polynomial_terms(computed_alexander) == EXPECTED_ALEXANDER[name]
        @test polynomial_terms(computed_jones) == EXPECTED_JONES[name]

        snapshot_alexander = _snapshot_terms(String(_field(row, headers,
                                                            "alexander_polynomial")))
        snapshot_jones = _snapshot_terms(String(_field(row, headers,
                                                       "jones_polynomial")))
        @test snapshot_alexander == EXPECTED_ALEXANDER[name]
        @test snapshot_jones == EXPECTED_JONES[name]
        @test snapshot_alexander == polynomial_terms(computed_alexander)
        @test snapshot_jones == polynomial_terms(computed_jones)

        reflected = mirror(diagram)
        reflected_pd = pd_from_gauss(reflected)
        reflected_certificate = classical_gauss_certificate(reflected)
        @test verify_pd(reflected_pd)
        @test reflected_certificate.realizable
        @test polynomial_terms(alexander_polynomial(reflected)) ==
              polynomial_terms(computed_alexander)
        @test jones_polynomial(reflected) == invert_knot_variable(computed_jones)
        @test polynomial_terms(alexander_polynomial(orientation_reverse(diagram))) ==
              polynomial_terms(computed_alexander)
        @test jones_polynomial(orientation_reverse(diagram)) == computed_jones
    end

    for name in TARGET_NAMES
        base = diagrams[name]
        reflected = mirror(base)
        @test canonical_serialization(base) isa String
        @test canonical_serialization(reflected) isa String
        @test canonical_hash(base) isa String
        @test canonical_hash(reflected) isa String
    end

    # The N=6 Jones comparisons are diagnostic only: equality is necessary for
    # an amphicheiral knot but is not, by itself, an amphicheirality proof.
    @test jones_polynomial(mirror(diagrams["6_1"])) != jones_polynomial(diagrams["6_1"])
    @test jones_polynomial(mirror(diagrams["6_2"])) != jones_polynomial(diagrams["6_2"])
    @test jones_polynomial(mirror(diagrams["6_3"])) == jones_polynomial(diagrams["6_3"])
end
