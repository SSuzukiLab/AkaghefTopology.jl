"""
    URAlgebra() -> URAlgebra

Select the exact U(R) determinant specialization of the o-graph contraction.

# Arguments
There are no arguments or input checks. Its later application uses U8–U9
of the [U(R) specification](../specs/ur-diagrams.md#numbered-conditions).

# Returns
A fieldless `URAlgebra` selector for the
[exact MST expression](../specs/ur-diagrams.md#exact-mst-expression).

# Guarantees
1. The result selects the exact rational determinant evaluator and its supported one-component domain. [proved/agent: by construction]
2. Construction does not modify any external state. [proved/agent: by construction]

# Conventions
[Exact MST word and matrix](../conventions/URDiagrams.md#exact-mst-word-and-matrix).

# Not guaranteed
This selector is not a finite-dimensional Hopf tensor record. It does not
certify the general identification of its determinant with the literature's
MST manifold invariant; that comparison remains unverified in the specification.

# Sources
cf:lib/topology.jl/docs/specs/ur-diagrams.md Exact MST expression — the determinant target and geometric evidence boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> abs(Z_MST(v, URAlgebra())) == 1//BigInt(1)
true
```

Role: construct
"""
struct URAlgebra end

"""
    URMSTExpression(matrix::Matrix{BigInt}) -> URMSTExpression

Package an integer trace matrix for the exact U(R) determinant evaluation.

# Arguments
`matrix` is a `Matrix{BigInt}` or converts to it; conversion errors propagate.
The raw constructor does not check squareness, U9, or whether the matrix
comes from the U8 conversion ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Those are caller obligations for a finite o-graph value.

# Returns
A `URMSTExpression` with its `matrix` field in the
[exact MST convention](../conventions/URDiagrams.md#exact-mst-word-and-matrix).

# Guarantees
1. The matrix field contains the supplied entries after conversion. An already correctly typed matrix is shared. [proved/agent: by construction]
2. Construction does not mutate the input or evaluate a determinant. [proved/agent: by construction]

# Not guaranteed
A raw matrix record does not certify an o-graph, a finite value or manifold
invariance. Mutating a shared matrix changes subsequent evaluations.

# Sources
cf:lib/topology.jl/docs/specs/ur-diagrams.md Exact MST expression — integer matrix record and nonsingular evaluation domain.

# Examples
```jldoctest
julia> x = URMSTExpression(zeros(BigInt,1,1));

julia> try evaluate_invariant(x); false catch err; err isa DomainError end
true
```

Role: construct
"""
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
