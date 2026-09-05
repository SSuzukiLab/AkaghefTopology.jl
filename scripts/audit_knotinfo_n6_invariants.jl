"""Audit the independent PD -> signed-Gauss -> invariant bridge for N<=6.

Run from the repository root with:

    julia --startup-file=no --pkgimages=no --project=. \
      lib/topology.jl/scripts/audit_knotinfo_n6_invariants.jl

The local KnotInfo snapshot supplies the eight PD representatives and the
final reference comparison only.  The Alexander and Jones values below are
computed by `knot_classical_invariants.jl`; no KnotInfo value is used as a
decision oracle or to select a crossing-sign decoration.
"""

if Base.find_package("JSON") === nothing
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using JSON
using SHA

include(joinpath(@__DIR__, "..", "src", "knot_classical_invariants.jl"))
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

function _check!(condition::Bool, message::AbstractString, failures::Vector{String})
    condition || push!(failures, String(message))
    nothing
end

payload, sheet, headers, rows = _load_snapshot(SNAPSHOT_PATH)
failures = String[]

_check!(payload["schema_version"] == 1, "schema_version != 1", failures)
_check!(payload["included_row_count"] == length(TARGET_NAMES),
        "included_row_count mismatch", failures)
_check!(Tuple(String.(payload["included_knot_names"])) == TARGET_NAMES,
        "included names/order mismatch", failures)
_check!(String(payload["source"]["snapshot_sha256"]) == SOURCE_SNAPSHOT_SHA256,
        "upstream source snapshot hash mismatch", failures)
_check!(_file_sha256(SNAPSHOT_PATH) == SNAPSHOT_FILE_SHA256,
        "local snapshot file hash mismatch", failures)
for (index, expected_header) in ((23, "gauss_notation"), (27, "pd_notation"),
                                  (29, "crossing_number"), (45, "braid_notation"),
                                  (63, "alexander_polynomial"),
                                  (67, "jones_polynomial"),
                                  (127, "symmetry_type"),
                                  (129, "full_symmetry_group"))
    _check!(headers[index] == expected_header,
            "header[$index] != $expected_header", failures)
end
_check!([Int(rows[name]["source_row_index"]) for name in TARGET_NAMES] == collect(2:9),
        "source row indices mismatch", failures)
_check!(sort(collect(keys(rows))) == sort(collect(TARGET_NAMES)),
        "row-name set mismatch", failures)

println("SNAPSHOT file=$(SNAPSHOT_PATH) file_sha256=$(_file_sha256(SNAPSHOT_PATH))")
source_hash = String(payload["source"]["snapshot_sha256"])
source_rows = join([Int(rows[name]["source_row_index"]) for name in TARGET_NAMES], ",")
println("SNAPSHOT source_sha256=$source_hash source_rows=$source_rows")
println("LABEL_CONVENTION plus=local KnotInfo pd_notation representative; " *
        "minus=mirror(diagram); no right/left-handedness claim")
println("BRIDGE pd_notation -> pd_to_gauss -> signed KnotGaussDiagram; " *
        "braid_notation/ordinary gauss_notation are retained as non-decision metadata")

diagrams = Dict{String,KnotGaussDiagram}()
for name in TARGET_NAMES
    row = rows[name]
    local_crossings = Int(row["crossing_number"])
    raw_pd = String(_field(row, headers, "pd_notation"))
    pd_rows = isempty(strip(raw_pd)) ? Any[] : JSON.parse(raw_pd)
    diagram = pd_to_gauss(pd_rows)
    diagrams[name] = diagram
    certificate = classical_gauss_certificate(diagram)
    computed_alexander = alexander_polynomial(diagram)
    computed_jones = jones_polynomial(diagram)
    reflected = mirror(diagram)
    reflected_certificate = classical_gauss_certificate(reflected)
    reflected_alexander = alexander_polynomial(reflected)
    reflected_jones = jones_polynomial(reflected)

    local_alexander = _snapshot_terms(String(_field(row, headers,
                                                    "alexander_polynomial")))
    local_jones = _snapshot_terms(String(_field(row, headers, "jones_polynomial")))
    _check!(length(pd_rows) == local_crossings, "$name PD row count mismatch", failures)
    _check!(crossing_count(diagram) == local_crossings, "$name crossing count mismatch", failures)
    _check!(verify_pd(pd_from_gauss(diagram)), "$name reconstructed PD failed verification", failures)
    _check!(certificate.route == :carter_surface && certificate.realizable &&
            certificate.genus == 0 && certificate.faces == local_crossings + 2,
            "$name Carter certificate is not planar", failures)
    _check!(reflected_certificate.realizable && reflected_certificate.genus == 0,
            "$name mirror Carter certificate is not planar", failures)
    _check!(polynomial_terms(computed_alexander) == EXPECTED_ALEXANDER[name],
            "$name independent Alexander mismatch", failures)
    _check!(polynomial_terms(computed_jones) == EXPECTED_JONES[name],
            "$name independent Jones mismatch", failures)
    _check!(local_alexander == EXPECTED_ALEXANDER[name] &&
            local_jones == EXPECTED_JONES[name],
            "$name local snapshot polynomial mismatch", failures)
    _check!(local_alexander == polynomial_terms(computed_alexander) &&
            local_jones == polynomial_terms(computed_jones),
            "$name computed/snapshot polynomial mismatch", failures)
    _check!(polynomial_terms(reflected_alexander) ==
            polynomial_terms(computed_alexander), "$name mirror Alexander mismatch", failures)
    _check!(reflected_jones == invert_knot_variable(computed_jones),
            "$name mirror Jones inversion mismatch", failures)

    symmetry = String(_field(row, headers, "symmetry_type"))
    full_symmetry = String(_field(row, headers, "full_symmetry_group"))
    signature = String(_field(row, headers, "signature"))
    source_row_index = Int(row["source_row_index"])
    println("ROW $name source_row_index=$source_row_index n=$local_crossings " *
            "planar=$(certificate.realizable) signs=$(collect(diagram.orientation)) " *
            "Alexander=$(knot_polynomial_string(computed_alexander)) " *
            "Jones=$(knot_polynomial_string(computed_jones)) " *
            "snapshot_match=true metadata_signature=$(repr(signature)) " *
            "symmetry_type=$(repr(symmetry)) full_symmetry_group=$(repr(full_symmetry))")
    println("IDENTITY $name+ canonical_hash=$(canonical_hash(diagram)) " *
            "canonical_serialization=$(canonical_serialization(diagram))")
    println("IDENTITY $name- canonical_hash=$(canonical_hash(reflected)) " *
            "canonical_serialization=$(canonical_serialization(reflected))")
    alexander_equal = polynomial_terms(reflected_alexander) ==
                      polynomial_terms(computed_alexander)
    jones_t_inverse = reflected_jones == invert_knot_variable(computed_jones)
    println("MIRROR $name alexander_equal=$alexander_equal jones_t_inverse=$jones_t_inverse")
end

for name in ("6_1", "6_2", "6_3")
    base_jones = jones_polynomial(diagrams[name])
    mirror_jones = jones_polynomial(mirror(diagrams[name]))
    equal = base_jones == mirror_jones
    println("N6_JONES_SEPARATION $name base_vs_mirror_equal=$equal separated=$(!equal)")
end
println("NOTE Jones equality is not a proof of amphicheirality; strict R-equivalence " *
        "identification remains a verified move-graph question.")

if isempty(failures)
    println("AUDIT PASS: eight local PD representatives and all eight mirrors passed " *
            "Carter-planarity, Alexander, Jones, and mirror checks")
else
    for failure in failures
        println("CONVENTION_MISMATCH: $failure")
    end
    error("KnotInfo N<=6 invariant audit failed")
end
