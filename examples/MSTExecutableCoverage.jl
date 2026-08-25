using Topology
using JSON
using LinearAlgebra

include(joinpath(@__DIR__, "C250811MSTinvUR.jl"))
include(joinpath(@__DIR__, "C250818MSTinv_uqsl2.jl"))
include(joinpath(@__DIR__, "C250823MSTinv_extalg.jl"))
include(joinpath(@__DIR__, "C250825Dsl2.jl"))
include(joinpath(@__DIR__, "C260220LensSpaces.jl"))

_jsonable(x::Complex) = Dict("re" => real(x), "im" => imag(x))
_jsonable(x::AbstractArray) = [_jsonable(v) for v in x]
_jsonable(x::AbstractDict) = Dict(string(k) => _jsonable(v) for (k, v) in x)
_jsonable(x::Tuple) = [_jsonable(v) for v in x]
_jsonable(x::NamedTuple) = Dict(string(k) => _jsonable(v) for (k, v) in pairs(x))
_jsonable(x) = x

"""Execute paper-relevant vertical slices and write bounded JSON/SVG evidence.

This is not a test runner. It records which mathematical objects, intermediate
representations, calculations, moves, and renderings are currently executable.
"""
function run_mst_executable_coverage(; output_dir=joinpath(@__DIR__, "..", "artifacts", "executable_coverage"))
    mkpath(output_dir)

    # MST1: finite Hopf -> Heisenberg double, plus concrete o-graph invariant.
    double_result = run_C250825_dsl2(order=3)
    move_result = run_C250811_mp_moves()
    conversions = run_C250811_conversion_cells()

    # MST2: integral o-graph -> weights -> Hopf values -> layout -> SVG.
    uqsl2 = run_C250818_uqsl2()
    exterior = run_C250823_extalg()
    lens = run_C260220_lens_spaces(maximum_s=10)
    lens_svg = plot_svg(lens.lens21, joinpath(output_dir, "MST_lens_L2_1.svg"))

    # MST2 Example 28: Z(L(2,1),f;u_q(sl2^+)) = q^{-1} Tr(S). Reproduce the RHS
    # independently and compare with the diagram-based invariant value.
    borel = uqsl2_borel_small(2 // 5; style=:L)
    antipode_matrix = zeros(ComplexF64, borel.dimension, borel.dimension)
    for ((o, i), v) in borel.antipode.terms
    antipode_matrix[o, i] = v
    end
    q = exp(4pi * im / 5)
    trace_formula = q^(-1) * tr(antipode_matrix)
    lens_formula_gap = maximum(abs.(uqsl2.lens_values .- trace_formula))

    move_shapes = Dict(name => Dict(
        "gauss" => link.gauss,
        "orientation" => link.orientation,
        "crossing_count" => length(link.orientation),
    ) for (name, link) in move_result.trials)

    report = Dict(
        "MST1" => Dict(
            "heisenberg_double" => Dict(
                "base_hopf_dimension" => double_result.algebra.dimension,
                "double_dimension" => double_result.heisenberg.dimension,
                "multiplication_term_count" => length(double_result.heisenberg.multiplication.terms),
                "pairing" => _jsonable(double_result.pairing),
                "commutator_nonzero_entries" => count(!iszero, double_result.commutator.coefficients),
            ),
            "paper_relevant_o_graph" => Dict(
                "lens_trace" => string(conversions[:lens].trace),
                "available_conversion_cases" => sort!(string.(collect(keys(conversions)))),
            ),
            "boundary" => [
                "Heisenberg double multiplication is executable.",
                "The paper's canonical element, Fock character, Tr(S) examples, and all 16 MP tensor equalities are not exposed as one reproduced workflow.",
            ],
        ),
        "MST2" => Dict(
            "lens_L2_1" => Dict(
                "gauss" => lens.lens21.gauss,
                "orientation" => lens.lens21.orientation,
                "weight_base" => lens.weight.base,
                "weight_kernel" => lens.weight.kernel,
                "hperp" => lens.weight.hperp,
                "disk_edge_sums" => lens.disk_edge_sums,
                "layout_valid" => is_valid_layout(lens.lens21),
                "failing_layout_rules" => string.(failing_rules(lens.lens21)),
                "svg" => basename(lens_svg),
            ),
            "lens_sweep_data_cases" => length(lens.sweep),
            "uqsl2_lens_values" => _jsonable(uqsl2.lens_values),
            "exterior_lens_values" => _jsonable(exterior.lens_values),
            "paper_example_28" => Dict(
                "claim" => "Z(L(2,1),f;u_q(sl2^+)) = q^{-1} Tr(S)",
                "q_inv_trace_S" => _jsonable(trace_formula),
                "diagram_invariant" => _jsonable(uqsl2.lens_values),
                "max_abs_gap" => lens_formula_gap,
                "reproduced" => lens_formula_gap < 1e-9,
            ),
            "boundary" => [
                "Concrete integral weights and finite-Hopf contractions are executable.",
                "Example 28's q^-1 Tr(S) is reproduced end-to-end for L(2,1); the general framing/WRT conjecture (Conjecture 29) is not certified.",
            ],
        ),
        "MST3" => Dict(
            "direct_MP_2_to_3" => move_shapes,
            "chained_A4_then_C2" => Dict(
                "gauss" => move_result.chained.gauss,
                "orientation" => move_result.chained.orientation,
                "crossing_count" => length(move_result.chained.orientation),
            ),
            "boundary" => [
                "Only direct non-weighted 2-to-3 MP moves used by C250811 are executable.",
                "Inverse 3-to-2, weighted transport, PS-I--IV condition checking, pure-sliding sequences, CP moves, and Theorem 3 decomposition are not executable workflows.",
            ],
        ),
        "qmd" => Dict(
            "qmd_files_present" => 0,
            "reading_copies" => ["threads/PJ260825/MST1.md", "threads/PJ260825/MST2.md", "threads/PJ260825/MST3.md"],
            "boundary" => "The current Markdown copies are mechanically converted reading text, not executable qmd documents.",
        ),
    )

    json_path = joinpath(output_dir, "mst_executable_coverage.json")
    open(json_path, "w") do io
        JSON.print(io, _jsonable(report), 2)
    end
    (; report, json_path, lens_svg)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    result = run_mst_executable_coverage()
    println("coverage_json=", result.json_path)
    println("lens_svg=", result.lens_svg)
end
