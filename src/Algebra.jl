"""Presentation-independent algebra: symbols, Laurent polynomials, finite and Hopf algebras, tensor contraction."""
module Algebra

using LinearAlgebra

include("symbolic.jl")
include("laurent.jl")
include("hopf_invariant.jl")
include("tensor_network.jl")
include("hopf_properties.jl")
include("finite_algebra.jl")

end
