"""
    is_involutory(algebra::FiniteHopfAlgebra; atol=1e-10, rtol=1e-10) -> Bool

Numerically test involutivity of the stored antipode on the algebra's coordinate space.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

- `algebra` supplies a nonnegative dimension and an antipode with valid matrix
  coordinates and coefficients convertible to `ComplexF64`. Conversion,
  allocation and indexing errors propagate; no Hopf-axiom check is performed.
- `atol` and `rtol` are finite nonnegative tolerances accepted by matrix
  `isapprox`. The caller supplies meaningful tolerances and finite tensor
  coefficients; no independent validation is performed here.

# Returns
A `Bool` indicating agreement of the squared stored antipode with the identity
map under the requested matrix comparison tolerances.

# Guarantees
1. The result is the numerical involutivity check on the supplied coordinate
   map after conversion to `ComplexF64`. Its intended algebraic property is
   involutivity of that Hopf algebra, conditional on valid Hopf data.
   [proved/agent: by construction]
2. The input algebra and its tensor containers are not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
`true` is not an exact proof of involutivity or of the Hopf axioms; `false`
can reflect floating-point error. The call does not test unimodularity or
establish topological invariance of another operation.

# Sources
cf:[[Kup91#^kup91-blk-0055]] involutory Hopf algebras — the mathematical target of this numerical predicate.

# Examples
```jldoctest
julia> (is_involutory(sweedler_algebra()), is_involutory(exterior_algebra(2)))
(false, true)
```

See also [`is_unimodular`](@ref), [`verify_hopf`](@ref).

Role: condition
"""
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

"""
    is_unimodular(algebra::FiniteHopfAlgebra; atol=1e-10, rtol=1e-10) -> Bool

Numerically compare the lines represented by the stored left and right integral vectors.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

- `algebra` supplies its dimension and integral vectors with valid coordinate
  indices and coefficients convertible to `ComplexF64`. Conversion, indexing
  and allocation errors propagate; the integral identities are not checked.
- `atol` and `rtol` must be finite nonnegative tolerances with meaningful
  arithmetic for the supplied finite coefficients. This requirement is the
  caller's responsibility and is not validated separately.

# Returns
A `Bool` from the numerical proportionality test of the two stored vectors,
rejecting either vector when its norm does not exceed `atol`.

# Guarantees
1. With valid finite inputs, `true` means that the stored vectors are nonzero
   above the absolute threshold and pass the proportionality comparison with
   the requested tolerances. This compares lines in the supplied coordinate
   space. [proved/agent: by construction]
2. Empty or numerically zero representatives return `false` for nonnegative
   tolerances. [proved/agent: by construction]
3. The algebra and its tensor containers are not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
This does not verify that the stored fields are actual left and right
integrals, identify their primal/dual convention, or prove exact
unimodularity. Very small rescalings can change the answer through the
absolute threshold. It does not certify counimodularity or the hypotheses
of a topological-invariance theorem.

# Sources
cf:[[MST25#^mst25-blk-0027]] unimodularity and counimodularity are separate hypotheses in the literature invariant; this predicate checks only the supplied integral-vector data.
cf:lib/topology.jl/src/hopf_properties.jl, is_unimodular — numerical coordinate-line comparison and its zero-vector boundary.

# Examples
```jldoctest
julia> (is_unimodular(sweedler_algebra()), is_unimodular(exterior_algebra(2)))
(false, true)
```

See also [`is_involutory`](@ref), [`verify_hopf`](@ref).

Role: condition
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
