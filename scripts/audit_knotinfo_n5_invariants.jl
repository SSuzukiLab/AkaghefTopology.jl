"""Audit the independent local PD -> signed-Gauss -> invariant bridge.

Run from the repository root with:

    julia --startup-file=no --pkgimages=no --project=. \
      lib/topology.jl/scripts/audit_knotinfo_n5_invariants.jl

The KnotInfo snapshot is read only for the PD rows and for the final displayed
cross-check.  Neither invariant implementation consults a KnotInfo value.
"""

if Base.find_package("JSON") === nothing
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end
using JSON

include(joinpath(@__DIR__, "..", "src", "knot_classical_invariants.jl"))

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SNAPSHOT_PATH = joinpath(REPO_ROOT, "threads", "PJ260830", "data",
                               "knotinfo-n6.json")
const TARGET_NAMES = ("0_1", "3_1", "4_1", "5_1", "5_2")

const EXPECTED_ALEXANDER = Dict(
    "0_1" => Dict(0 => big(1)),
    "3_1" => Dict(0 => big(1), 1 => big(-1), 2 => big(1)),
    "4_1" => Dict(0 => big(1), 1 => big(-3), 2 => big(1)),
    "5_1" => Dict(0 => big(1), 1 => big(-1), 2 => big(1),
                  3 => big(-1), 4 => big(1)),
    "5_2" => Dict(0 => big(2), 1 => big(-3), 2 => big(2)),
)

const EXPECTED_JONES = Dict(
    "0_1" => Dict(0 => big(1)),
    "3_1" => Dict(1 => big(1), 3 => big(1), 4 => big(-1)),
    "4_1" => Dict(-2 => big(1), -1 => big(-1), 0 => big(1),
                  1 => big(-1), 2 => big(1)),
    "5_1" => Dict(2 => big(1), 4 => big(1), 5 => big(-1),
                  6 => big(1), 7 => big(-1)),
    "5_2" => Dict(1 => big(1), 2 => big(-1), 3 => big(2),
                  4 => big(-1), 5 => big(1), 6 => big(-1)),
)

function _snapshot_terms(text::AbstractString)
    source = replace(replace(strip(String(text)), " " => ""), "*" => "")
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
            parenthesis_depth == 0 && (byte == UInt8('+') || byte == UInt8('-')) && break
            position += 1
        end
        token = String(bytes[start:position - 1])
        token = replace(token, "(" => "", ")" => "")
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

function _pd_rows(text)
    isempty(strip(String(text))) && return Any[]
    JSON.parse(String(text))
end

snapshot = JSON.parsefile(SNAPSHOT_PATH)
sheet = first(snapshot["sheets"])
headers = String.(sheet["headers"])
index(name) = findfirst(==(name), headers)
name_index = index("name")
crossing_index = index("crossing_number")
pd_index = index("pd_notation")
alexander_index = index("alexander_polynomial")
jones_index = index("jones_polynomial")
rows = Dict(String(row["values"][name_index]) => row for row in sheet["rows"])

failures = String[]
for name in TARGET_NAMES
    row = rows[name]
    local_crossings = parse(Int, String(row["values"][crossing_index]))
    diagram = pd_to_gauss(_pd_rows(row["values"][pd_index]))
    pd = pd_from_gauss(diagram)
    verify_pd(pd)
    computed_alexander = alexander_polynomial(diagram)
    computed_jones = jones_polynomial(diagram)
    expected_alexander = EXPECTED_ALEXANDER[name]
    expected_jones = EXPECTED_JONES[name]
    local_alexander = _snapshot_terms(row["values"][alexander_index])
    local_jones = _snapshot_terms(row["values"][jones_index])
    selected = canonical_serialization(diagram)
    println(name, " crossings=", local_crossings,
            " selected=", selected,
            " crossing_signs=", collect(diagram.orientation))
    println("  Alexander: ", knot_polynomial_string(computed_alexander),
            " | snapshot=", row["values"][alexander_index])
    println("  Jones:     ", knot_polynomial_string(computed_jones),
            " | snapshot=", row["values"][jones_index])
    if crossing_count(diagram) != local_crossings
        push!(failures, "$name crossing count mismatch")
    end
    polynomial_terms(computed_alexander) == expected_alexander ||
        push!(failures, "$name Alexander does not match expected independent value")
    polynomial_terms(computed_jones) == expected_jones ||
        push!(failures, "$name Jones does not match expected independent value")
    local_alexander == expected_alexander ||
        push!(failures, "$name local snapshot Alexander parse mismatch")
    local_jones == expected_jones ||
        push!(failures, "$name local snapshot Jones parse mismatch")
end

if isempty(failures)
    println("AUDIT PASS: five local PD rows independently reconstructed and separated")
else
    for failure in failures
        println("CONVENTION_MISMATCH: ", failure)
    end
    error("KnotInfo invariant audit failed")
end
