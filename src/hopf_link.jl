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
