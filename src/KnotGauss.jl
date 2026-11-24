"""Classical knot Gauss diagrams, Reidemeister moves, and the bounded censuses."""
module KnotGauss

using ..MoveGraphStore
import ..MoveAPI: movable, move, move!
import ..Invariants: mirror

include("knot_gauss.jl")
include("reidemeister.jl")
include("classical_gauss_realizability.jl")
include("knot_classical_invariants.jl")
include("classical_gauss_enumeration.jl")
include("knot_n6_classification.jl")
include("knot_n6_targeted_bridge.jl")
include("knot_n6_targeted_barrier8.jl")
include("knot_n6_all_state_enumeration.jl")
include("knot_move_db.jl")
include("knot_n5_classification.jl")

end
