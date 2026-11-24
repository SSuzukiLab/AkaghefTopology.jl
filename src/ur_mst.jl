"""U(R) MST trace convention; independent of URDiagram rendering and rewrites.

The involutory U(R) specialization does not use the integer edge weights.
The current evaluator covers one nonempty closed real-Gauss component and
returns an exact rational number when the defining matrix is nonsingular.
"""
struct URAlgebra end

"The defining U(R) trace matrix, in reversed Gauss traversal order."
struct URMSTExpression
    matrix::Matrix{BigInt}
end

function simplify_invariant(x::URMSTExpression; strategy=:auto)
    strategy in (:auto, :determinant, :none) ||
        throw(ArgumentError("U(R) strategies are :auto, :determinant and :none"))
    URMSTExpression(copy(x.matrix))
end

function evaluate_invariant(x::URMSTExpression)
    determinant = det(Rational{BigInt}.(x.matrix))
    iszero(determinant) && throw(DomainError(determinant,
        "singular U(R) trace matrix; no finite rational value is defined by this evaluator"))
    -inv(determinant)
end

export URAlgebra, URMSTExpression
