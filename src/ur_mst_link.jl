"""
    invariant_expression(v::VirtualLink, H::URAlgebra) -> URMSTExpression

Return the exact integer trace matrix associated with one closed real-Gauss component.

# Arguments
Requires U8 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Trailing cut markers, component count, emptiness, signed visits and crossing
signs are checked with `ArgumentError` identifying the unsupported input.
The caller supplies coherent VirtualLink data; other conversion/indexing
errors propagate. Full closed-normality is not certified by these guards.
`H` selects the exact determinant convention.

# Returns
A `URMSTExpression` with the matrix of the
[exact MST convention](../conventions/URDiagrams.md#exact-mst-word-and-matrix).

# Guarantees
1. Matrix entries are exact integers in the specified real-Gauss coordinate construction. Assigned edge weights do not enter them. [proved/agent: by construction]
2. Inputs are unchanged and the matrix is newly allocated. [proved/agent: by construction]

# Not guaranteed
The matrix can be singular. No equality with the unreduced legacy trace or
general identification with a framed/combed manifold invariant is established;
the specification records the finite evidence and missing source comparison.

# Sources
cf:lib/topology.jl/docs/specs/ur-diagrams.md Exact MST expression — determinant target and geometric evidence boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> abs(evaluate_invariant(invariant_expression(v,URAlgebra()))) == 1//BigInt(1)
true
```

Role: convert
"""
function invariant_expression(v::VirtualLink, ::URAlgebra)
    _require_closed_ograph_input(v, "U(R) MST")
    components, signs = real_gauss_code(v)
    length(components) == 1 && !isempty(only(components)) ||
        throw(ArgumentError("U(R) MST requires one nonempty real-Gauss component"))
    word = reverse(only(components))
    n = length(signs)
    all(s -> s in (-1, 1), signs) &&
        sort(word) == vcat(collect(-n:-1), collect(1:n)) ||
        throw(ArgumentError("U(R) MST requires one positive and negative visit per vertex"))
    word = circshift(word, 1 - findfirst(>(0), word))
    section = Dict(zip(word, cumsum(word .> 0)))
    matrix = Matrix{BigInt}(I, n, n)
    matrix -= circshift(Matrix{BigInt}(I, n, n), (1, 0))
    for vertex in 1:n
        matrix[section[vertex], section[-vertex]] += signs[vertex]
    end
    URMSTExpression(matrix)
end

"""
    Z_MST(v::VirtualLink, H::URAlgebra; strategy=:auto) -> Rational{BigInt}

Evaluate the exact U(R) determinant expression of one closed real-Gauss component.

# Arguments
Requires U8–U9 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
The guards and errors of `invariant_expression(v,H)` apply. A singular
matrix raises `DomainError` stating that no finite rational value is defined.
`strategy` accepts `:auto`, `:determinant` or `:none`; unsupported names raise
`ArgumentError`. Other keywords are not accepted.

# Returns
The exact `Rational{BigInt}` value in the
[MST determinant convention](../conventions/URDiagrams.md#exact-mst-word-and-matrix).

# Guarantees
1. The result is the specified rational determinant expression and is independent of assigned edge weights and accepted strategy. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
No analytic regularization is provided for a singular matrix. A general
source-invariant identification, equivalence to the legacy regulated trace,
or universal positive sign is not established by the finite fixture checks.

# Sources
cf:lib/topology.jl/docs/specs/ur-diagrams.md Exact MST expression — exact coordinate value and its finite geometric comparisons.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> abs(Z_MST(v,URAlgebra())) == 1//BigInt(1)
true
```

Role: invariant
"""
Z_MST(v::VirtualLink, H::URAlgebra; strategy=:auto) =
    evaluate_invariant(simplify_invariant(invariant_expression(v,H); strategy=strategy))
