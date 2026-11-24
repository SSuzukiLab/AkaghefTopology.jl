"""O-graphs and links: the `VirtualLink` presentation, its moves, weights, and search database."""
module VirtualLinks

using LinearAlgebra
using JuMP
using GLPK
using Graphs
using ..MoveGraphStore
import ..Algebra: SparseFactor, FiniteHopfAlgebra, TaftKind, _contract, _factor_with_variables
import ..URDiagrams: URDiagram, URAlgebra, URMSTExpression, set_data_cvs!, trace_value
import ..MoveAPI: movable, move, move!,
    search_state, put_search_state!, load_search_state,
    begin_search_run!, finish_search_run!, expand_search_state!,
    search_neighbors, query_search_paths, verify_search_store
import AkaghefAlgebra
import AkaghefAlgebra: simplify_invariant, evaluate_invariant

include("virtual_link.jl")
include("framed_lens_space.jl")
include("circuit_b.jl")
include("integer_weight_solver.jl")
include("ograph_state.jl")
include("ograph_colours.jl")
include("ograph_colour_normal.jl")
include("ograph_face_cycles.jl")
include("weight_validation.jl")
include("spin_weights.jl")
include("layer8_planarization.jl")
include("virtual_planarization.jl")
include("orthogonal_layout.jl")
include("layout_invariants.jl")
include("helpers.jl")
include("moves.jl")
include("ograph_move_verifier.jl")
include("ograph_tetrahedra.jl")
include("weighted_ps.jl")
include("weighted_move_verifier.jl")
include("spin_cp.jl")
include("hopf_link.jl")
include("taft_invariant.jl")
include("ur_mst_link.jl")
include("ograph_gauss_diagram.jl")
include("ograph_move_db.jl")
include("ograph_search_db.jl")
include("ograph_path_constraints.jl")
include("ograph_path_lifting.jl")
include("ograph_invariants.jl")

# VirtualLink methods of the MoveAPI search names. The ograph_search_* functions stay.
search_state(link::VirtualLink) = ograph_search_state(link)
put_search_state!(store::MoveGraphStoreDB, link::VirtualLink) =
    put_ograph_search_state!(store, link)
load_search_state(store::MoveGraphStoreDB, state_id::AbstractString) =
    load_ograph_search_state(store, state_id)
function begin_search_run!(store::MoveGraphStoreDB, script_path::AbstractString;
                           label::AbstractString, scope=Dict{String,Any}())
    begin_ograph_search_run!(store, script_path; label, scope)
end
function finish_search_run!(store::MoveGraphStoreDB, run_id::AbstractString;
                            status::AbstractString, counts=Dict{String,Any}())
    finish_ograph_search_run!(store, run_id; status, counts)
end
function expand_search_state!(store::MoveGraphStoreDB, run_id::AbstractString,
                              source::VirtualLink, move_kind::AbstractString;
                              selector=Dict{String,Any}())
    expand_ograph_search_state!(store, run_id, source, move_kind; selector)
end
function search_neighbors(store::MoveGraphStoreDB, state_id::AbstractString; kinds=("MP", "02"))
    ograph_search_neighbors(store, state_id; kinds)
end
function query_search_paths(store::MoveGraphStoreDB, source_state_id::AbstractString; kwargs...)
    query_ograph_search_paths(store, source_state_id; kwargs...)
end
verify_search_store(store::MoveGraphStoreDB) = verify_ograph_search_store(store)

end
