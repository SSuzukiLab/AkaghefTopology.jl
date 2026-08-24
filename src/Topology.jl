module Topology

using LinearAlgebra
using JuMP
using GLPK

export URDiagram, VirtualLink, SymExpr, sym,
       LaurentPolynomial, LaurentRational, laurent_variable, laurent_variables, evaluate, limit_zero,
       set_data_cvw!, set_data_cvs!, set_data_v!, set_data_cm!, set_data!, set_weight!,
       set_weight_by_crossing!,
       validate_data, check_data, swap!, delete_edge!, add_edges!, dilation!,
       compose_dilation!, put_zero!, reduct3!, reduct4!, reduct5!, reduct6!, reduct7!,
       reduct8!, simplify1!, simplify2!, simplify2_history!, simplify!,
       trace_value, trace1, trace2, trace3, standard_matrix_adjustment, standard_matrix_trace,
       rank, odata, weight_expression, unicode_art, trans_knot_complement,
       gauss_code, real_gauss_code, real_edge_table, RealEdge, pd_code, edge_directions, regions,
       head_map, is_closed,
       hmove_matrix, disk_table, calculate_weight, mirror_manifold!, cyclic_shift!, mp_move!,
       connected_sum!, disjoint_sum!, knot_complement!, vl_to_urdiagram,
       ur_invariant, ur_prime_invariant, EdgeLayout, LinkLayout,
       minimal_bending_numbers, orthogonal_layout, plot_svg, plot_urd_structure_svg, bernstein, bezier,
       LayoutDiagnostic, check_layout, is_valid_layout, failing_rules, LAYOUT_RULES,
       planarize_virtual_gauss,
       MovePattern, make_move_data, SmallObject
export SparseFactor, FiniteHopfAlgebra, uqsl2_borel_small, sweedler_algebra,
       exterior_algebra, hopf_invariant, label_tensor, contract_network,
       dense_tensor, hopf_tensors, evaluate_tensor_expression, verify_hopf
export legacy_hopf_tensors
export FiniteAlgebra, AlgebraElement, basis_element, cyclic_group_algebra,
       kac_paljutkin_algebra, dual_algebra, hopf_pairing, heisenberg_double,
       drinfeld_double

include("symbolic.jl")
include("laurent.jl")
include("ur_diagram.jl")
include("virtual_link.jl")
include("virtual_planarization.jl")
include("orthogonal_layout.jl")
include("layout_invariants.jl")
include("helpers.jl")
include("hopf_invariant.jl")
include("tensor_network.jl")
include("finite_algebra.jl")
include("svg_plot.jl")

end
