# Diagnostic direct driver, not a replacement for fresh-process library acceptance.
# cf: lib/topology.jl/docs/specs/api-contracts.md §§2–2.1 — contract completeness.
# cf: .research-scratch/codex-lanes/lib-api/CONTEXT.md — docstring-only ownership.
# Example provenance for turn 9 (paths relative to lib/topology.jl):
# - MP and template/catalogue examples: test/cut_moves.jl, test/topology_core.jl.
# - Integer weights and ranks: test/weight_moves.jl, test/cut_weights.jl,
#   test/integer_weight_solver.jl; the one-by-one parity obstruction also
#   follows directly from its scalar equation.
# - Binary/relative spin equations: test/spin_weights.jl,
#   test/cut_spin_weights.jl; CP case zero: test/spin_cp.jl.
# - H verification/replay/decoding: test/weighted_move_verifier.jl.
#   Rejected mixed coefficient modes: test/spin_evidence_boundaries.jl.
# - Empty planarisation: test/virtual_planarization.jl and the two source
#   empty branches. Empty layouts/diagnostics: the corresponding source branches.
# - Candidate constructors, PS split, envelope copying and hash length:
#   direct field/copy/digest definitions in their source files.
# - Search queries: empty-filter branches in src/heegaard/hc_search_db.jl
#   and src/ograph_search_db.jl; literal HC identity: the genus-one input
#   and CanonicalState construction in the HC source.
# - UR diagnostic: test/topology_core.jl's three-edge endpoint order, with
#   numeric weights; the asserted title is emitted by src/svg_plot.jl.
# These derivations are proposals for the lead's full example checker.
empty!(LOAD_PATH)
push!(LOAD_PATH, "@stdlib")
include(joinpath(@__DIR__, "api_index.jl"))

isdoc(e) = e isa Expr && e.head === :macrocall && length(e.args) >= 4 &&
    (e.args[1] === Symbol("@doc") || e.args[1] isa GlobalRef && e.args[1].name === Symbol("@doc"))
stripdocs(x) = x
function stripdocs(e::Expr)
    isdoc(e) && return stripdocs(e.args[end])
    Expr(e.head, (stripdocs(x) for x in e.args if !(x isa LineNumberNode))...)
end
function docs!(out, e)
    e isa Expr || return out
    e.head in (:error, :incomplete) && error("Julia parse error: $e")
    if isdoc(e)
        e.args[3] isa String && push!(out, e.args[3])
        docs!(out, e.args[end])
    else
        foreach(x -> docs!(out, x), e.args)
    end
    out
end

function main()
root = normpath(joinpath(@__DIR__, ".."))
files = isempty(ARGS) ? ["weighted_move_verifier.jl", "weight_validation.jl",
    "spin_weights.jl", "spin_cp.jl", "integer_weight_solver.jl", "weighted_ps.jl",
    "moves.jl", "helpers.jl", "Topology.jl", "virtual_planarization.jl",
    "layer8_planarization.jl", "layout_invariants.jl", "orthogonal_layout.jl", "svg_plot.jl"] : ARGS
cache = Dict{String,APIIndex.SpecPage}()
failures = String[]
changed = 0
targets = Set(split("search_neighbors search_state MoveCandidate MovePattern boundary_weight_data check_layout failing_rules first_betti_number is_framed_weight is_spin_weight is_valid_layout make_move_data minimal_bending_numbers move_boundary_transport mp_move! orthogonal_layout planarize_virtual_gauss planarize_virtual_gauss_layer8 replay_weighted_move solve_integer_weight_lift solve_integer_weight_system solve_integer_weights solve_spin_weight_system solve_spin_weights spin_boundary_data spin_cp_move spin_cp_transport spin_weight_equation validate_weight verify_integer_weight_solution verify_weighted_move verify_weighted_relative_result weight_residual weighted_move_evidence_data weighted_move_evidence_from_data weighted_ps_candidate weighted_raw_hash plot_urd_structure_svg"))
reached = Set{String}()
target_contracts = 0
for file in files
    relative = "src/" * file
    old = Meta.parseall(read(Cmd(Cmd(["git", "show", "HEAD:" * relative]); dir=root), String))
    current = Meta.parseall(read(joinpath(root, relative), String))
    stripdocs(old) == stripdocs(current) || push!(failures, "$file: code changed")
    previous = Set(docs!(String[], old))
    for raw in docs!(String[], current)
        raw in previous && continue
        changed += 1
        result = APIIndex.contract((raw=raw,), "unknown(...)")
        name = first(split(result.signature, '('))
        if name in targets
            push!(reached, name)
            target_contracts += 1
        end
        label = "$file: $(result.signature)"
        result.status == "complete" || push!(failures, "$label: $(join(result.missing, ", "))")
        for link in result.links
            target = APIIndex.spec_target(link)
            target === nothing && continue
            ref = APIIndex.spec_reference(target, cache)
            isempty(ref.error) || push!(failures, "$label: $link: $(ref.error)")
        end
    end
end
println("Changed docstrings versus HEAD: ", changed)
if isempty(ARGS)
    println("Turn 9 target bindings: ", length(reached), "; contracts: ", target_contracts)
    isempty(setdiff(targets, reached)) || push!(failures, "Unreached: $(sort!(collect(setdiff(targets, reached))))")
else
    println("Selected source files: ", length(files))
end
foreach(println, failures)
isempty(failures) || error("$(length(failures)) contract/code diagnostics failed")
println("Docstring-stripped ASTs agree; changed contracts and specification anchors pass the diagnostic.")
end
main()
