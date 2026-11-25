"""
    hopf_invariant(v::VirtualLink, algebra::FiniteHopfAlgebra; method=:auto, strategy=:auto, kwargs...) -> Any

Evaluate the repository's closed o-graph tensor contraction with the supplied Hopf-algebra data.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

- For the closed-normal geometric interpretation, the caller supplies N1,
  N2 and C1–C3 of [the terminology specification §11](../specs/ograph-terminology.md#11-closed-normal-conditions-n1-n2-c1-c2-c3).
  This entry point rejects cut markers but does not certify those conditions.
- Assigned weights must satisfy the integral-weight checks of
  [`validate_weight`](@ref); the contraction checks integer mode and finite
  integer real-arc weights, with `ArgumentError` identifying a failed check.
  It does not check the framing equation. Unweighted input uses zero edge weights.
- `algebra` must supply coherent structural and crossing tensors in the
  repository's basis convention. The Hopf axioms and tensor shapes are not
  validated. Singular inverse antipode powers and arithmetic failures propagate.
- `method` is `:auto` or `:numeric`; otherwise `ArgumentError` identifies
  the supported numeric path. Extra keywords are rejected as unsupported.
  With `:auto`, `strategy` is passed to [`Z_MST`](@ref); with `:numeric` it
  does not select a reduction strategy.

# Returns
The scalar contraction value; its type follows the supplied tensor arithmetic
and selected backend, so the return type is `Any`.

# Guarantees
1. On supported coherent data, the result evaluates the repository's
   crossing-tensor network with its integer edge-weight insertions.
   [proved/agent: by construction]
2. The supplied diagram and algebra tensors are not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
This call does not certify that the raw network matches every convention of
the literature invariant or that its input represents a framed manifold with
the required hypotheses. Consequently the result alone does not establish
invariance under changes of diagram, combing or framing. Floating-point paths
provide no exactness or error bound, and nonempty, well-formed contraction
data remain the caller's responsibility.

# Sources
cf:[[MST25#^mst25-blk-0025]] positive crossings and Heisenberg-double tensors — mathematical target of the network construction, not a certification of the repository's port dictionary.
cf:[[MST25#^mst25-blk-0027]] framed, combed and topological invariant hypotheses — additional geometric meanings not automatically certified by this evaluator.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]);

julia> hopf_invariant(v, sweedler_algebra()) ≈ -1
true
```

See also [`weight_equation`](@ref), [`verify_hopf`](@ref).

Role: invariant
"""
function hopf_invariant(v::VirtualLink,algebra::FiniteHopfAlgebra;method=:auto,strategy=:auto,kwargs...)
    method in (:auto,:numeric) || throw(ArgumentError("generic algebras support method=:numeric"))
    isempty(kwargs) || throw(ArgumentError("unsupported numeric contraction keywords"))
    method==:numeric ? _hopf_invariant_numeric(v,algebra) : Z_MST(v,algebra;strategy=strategy)
end

function _hopf_invariant_factors(v::VirtualLink,algebra::FiniteHopfAlgebra)
    _require_closed_ograph_input(v, "hopf_invariant")
    edges=real_edge_table(v)
    integer_weights=v.is_weighted ? _integer_edge_weights(v) : Int[]
    orientation=last(real_gauss_code(v))
    factors=SparseFactor[]
    zero_map=Dict{Int,Int}()
    for (edge_index,edge) in enumerate(edges)
        labels=4 .* abs.([edge.crossing[1],edge.crossing[2]]) .+
               ([edge.crossing[1] > 0,edge.crossing[2] > 0]) .- [1,3]
        exponent=v.is_weighted ? 2integer_weights[edge_index] : 0
        if exponent == 0
            zero_map[labels[1]]=labels[2]
        else
            matrix=zeros(ComplexF64,algebra.dimension,algebra.dimension)
            for ((output,input),value) in algebra.antipode.terms; matrix[output,input]=value; end
            powered=matrix^exponent
            terms=Dict((output,input)=>powered[output,input]
                       for output in 1:algebra.dimension,input in 1:algebra.dimension
                       if abs(powered[output,input]) > 1e-12)
            push!(factors,SparseFactor(labels,terms))
        end
    end
    resolve(label)=haskey(zero_map,label) ? resolve(zero_map[label]) : label
    for vertex in eachindex(orientation)
        labels=resolve.(4(vertex-1) .+ (1:4))
        tensor=orientation[vertex] < 0 ? algebra.negative_vertex : algebra.positive_vertex
        push!(factors,_factor_with_variables(tensor,labels))
    end
    factors
end

_hopf_invariant_numeric(v::VirtualLink,algebra::FiniteHopfAlgebra) =
    _contract(_hopf_invariant_factors(v,algebra))
