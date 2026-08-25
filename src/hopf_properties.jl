"""Return whether the antipode satisfies `S² = id` in the stored basis."""
function is_involutory(algebra::FiniteHopfAlgebra;atol=1e-10,rtol=1e-10)
    dimension=algebra.dimension
    antipode=zeros(ComplexF64,dimension,dimension)
    for ((output,input),value) in algebra.antipode.terms
        antipode[output,input]=value
    end
    isapprox(antipode*antipode,Matrix{ComplexF64}(I,dimension,dimension);atol,rtol)
end

function _factor_vector(factor::SparseFactor,dimension::Int)
    vector=zeros(ComplexF64,dimension)
    for ((index,),value) in factor.terms
        vector[index]=value
    end
    vector
end

"""Return whether the stored left and right integrals span the same line.

The representatives may use different nonzero normalizations. Empty or
numerically zero representatives are rejected instead of being classified as
unimodular.
"""
function is_unimodular(algebra::FiniteHopfAlgebra;atol=1e-10,rtol=1e-10)
    left=_factor_vector(algebra.left_integral,algebra.dimension)
    right=_factor_vector(algebra.right_integral,algebra.dimension)
    left_norm=norm(left); right_norm=norm(right)
    left_norm>atol && right_norm>atol || return false
    pivot=argmax(abs.(right))
    scale=left[pivot]/right[pivot]
    norm(left-scale*right)<=atol+rtol*max(left_norm,abs(scale)*right_norm)
end