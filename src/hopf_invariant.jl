struct SparseFactor
    variables::Vector{Int}
    terms::Dict{Tuple,ComplexF64}
end

struct FiniteHopfAlgebra
    dimension::Int
    multiplication::SparseFactor
    unit::SparseFactor
    coproduct::SparseFactor
    counit::SparseFactor
    positive_vertex::SparseFactor
    negative_vertex::SparseFactor
    antipode::SparseFactor
    right_integral::SparseFactor
    right_cointegral::SparseFactor
    left_integral::SparseFactor
    left_cointegral::SparseFactor
end

_sparse_data(terms)=SparseFactor(Int[],Dict{Tuple,ComplexF64}(
    Tuple(key)=>ComplexF64(value) for (key,value) in terms))

_qbinomial(q,n,k)=begin
    k < 0 || k > n ? 0.0+0im :
    prod((1-q^(n-j))/(1-q^(j+1)) for j in 0:k-1;init=1.0+0im)
end

const _UQSL2_BOREL_CACHE=Dict{Tuple{Any,Symbol},FiniteHopfAlgebra}()

"""Build one canonical small Borel algebra instance (uncached internal path)."""
function _build_uqsl2_borel_small(ratio::Rational;style=:L)
    numerator(ratio) > 0 && numerator(ratio) < denominator(ratio) ||
        throw(ArgumentError("ratio must lie in (0,1)"))
    M,N=numerator(ratio),denominator(ratio)
    q=N == 2 || N == 4 ? im^(M*4÷N) : exp(2pi*im*M/N)
    dimension=N^2
    basis(i,j)=i+N*j+1

    multiplication=Dict{Tuple,ComplexF64}()
    for i1 in 0:N-1,i2 in 0:N-1,j1 in 0:N-1,j2 in 0:N-1
        j1+j2 >= N && continue
        key=(basis(i1,j1),basis(i2,j2),basis(mod(i1+i2,N),j1+j2))
        multiplication[key]=q^(-j1*i2)
    end
    coproduct=Dict{Tuple,ComplexF64}()
    for i in 0:N-1,j in 0:N-1,k in 0:j
        key=(basis(i,j),basis(i,k),basis(mod(i+k,N),j-k))
        coproduct[key]=_qbinomial(q,j,k)*q^(-k*(j-k))
    end
    antipode=Dict{Tuple,ComplexF64}()
    for i in 0:N-1,j in 0:N-1
        antipode[(basis(mod(-i-j,N),j),basis(i,j))]=(-1)^j*q^(i*j+j*(j+1)÷2)
    end
    if style == :L
        coproduct=Dict((input,right,left)=>value for ((input,left,right),value) in coproduct)
        matrix=zeros(ComplexF64,dimension,dimension)
        for ((output,input),value) in antipode; matrix[output,input]=value; end
        inverse_matrix=inv(matrix)
        antipode=Dict((output,input)=>inverse_matrix[output,input]
                      for output in 1:dimension,input in 1:dimension
                      if abs(inverse_matrix[output,input]) > 1e-12)
    elseif style != :R
        throw(ArgumentError("style must be :L or :R"))
    end

    positive=Dict{Tuple,ComplexF64}()
    for ((delta_input,left,right),delta_value) in coproduct
        for ((mu_input,middle,output),mu_value) in multiplication
            middle == left || continue
            key=(output,right,mu_input,delta_input)
            positive[key]=get(positive,key,0)+mu_value*delta_value
        end
    end
    negative=Dict{Tuple,ComplexF64}()
    antipode_by_input=Dict{Int,Vector{Tuple{Int,ComplexF64}}}()
    for ((output,input),value) in antipode
        push!(get!(antipode_by_input,input,Tuple{Int,ComplexF64}[]),(output,value))
    end
    for ((delta_input,left,right),delta_value) in coproduct
        for (antipode_output,antipode_value) in get(antipode_by_input,left,Tuple{Int,ComplexF64}[])
            for ((mu_input,middle,output),mu_value) in multiplication
                middle == antipode_output || continue
                key=(output,right,mu_input,delta_input)
                negative[key]=get(negative,key,0)+mu_value*antipode_value*delta_value
            end
        end
    end
    unit=Dict{Tuple,ComplexF64}((1,)=>1)
    counit=Dict{Tuple,ComplexF64}((basis(i,0),)=>1 for i in 0:N-1)
    right_integral=Dict{Tuple,ComplexF64}((basis(0,N-1),)=>q)
    right_cointegral=Dict{Tuple,ComplexF64}(
        (basis(i,N-1),)=>q^(i-1) for i in 0:N-1)
    antipode_matrix=zeros(ComplexF64,dimension,dimension)
    for ((output,input),value) in antipode; antipode_matrix[output,input]=value; end
    inverse_antipode=inv(antipode_matrix)
    right_cointegral_vector=zeros(ComplexF64,dimension)
    for ((index,),value) in right_cointegral; right_cointegral_vector[index]=value; end
    right_integral_vector=zeros(ComplexF64,dimension)
    for ((index,),value) in right_integral; right_integral_vector[index]=value; end
    left_cointegral_vector=inverse_antipode*right_cointegral_vector
    left_integral_vector=transpose(antipode_matrix)*right_integral_vector
    left_cointegral=Dict{Tuple,ComplexF64}((index,)=>value
        for (index,value) in enumerate(left_cointegral_vector) if abs(value)>1e-12)
    left_integral=Dict{Tuple,ComplexF64}((index,)=>value
        for (index,value) in enumerate(left_integral_vector) if abs(value)>1e-12)
    FiniteHopfAlgebra(dimension,_sparse_data(multiplication),_sparse_data(unit),
                      _sparse_data(coproduct),_sparse_data(counit),
                      _sparse_data(positive),_sparse_data(negative),_sparse_data(antipode),
                      _sparse_data(right_integral),_sparse_data(right_cointegral),
                      _sparse_data(left_integral),_sparse_data(left_cointegral))
end

"""Return an isolated small Borel algebra while caching its expensive build.

The algebra contains mutable sparse dictionaries.  The cached object is never
returned directly, so callers retain the old ownership semantics even if they
mutate a factor's terms.
"""
function uqsl2_borel_small(ratio::Rational;style=:L)
    key=(ratio,style)
    prototype=get!(_UQSL2_BOREL_CACHE,key) do
        _build_uqsl2_borel_small(ratio;style=style)
    end
    deepcopy(prototype)
end

function sweedler_algebra()
    multiplication=Dict{Tuple,ComplexF64}(
        (1,1,1)=>1,(1,2,2)=>1,(2,1,2)=>1,(1,3,3)=>1,(3,1,3)=>1,
        (1,4,4)=>1,(4,1,4)=>1,(2,2,1)=>1,(2,3,4)=>-1,(3,2,4)=>1,
        (2,4,3)=>-1,(4,2,3)=>1)
    coproduct=Dict{Tuple,ComplexF64}(
        (1,1,1)=>1,(2,2,2)=>1,(3,1,3)=>1,(3,3,2)=>1,
        (4,2,4)=>1,(4,4,1)=>1)
    antipode=Dict{Tuple,ComplexF64}((1,1)=>1,(2,2)=>1,(4,3)=>-1,(3,4)=>1)
    positive=Dict{Tuple,ComplexF64}()
    negative=Dict{Tuple,ComplexF64}()
    for ((delta_input,left,right),delta_value) in coproduct
        for ((mu_input,middle,output),mu_value) in multiplication
            middle == left || continue
            key=(output,right,mu_input,delta_input)
            positive[key]=get(positive,key,0)+mu_value*delta_value
        end
        for ((antipode_output,antipode_input),antipode_value) in antipode
            antipode_input == left || continue
            for ((mu_input,middle,output),mu_value) in multiplication
                middle == antipode_output || continue
                key=(output,right,mu_input,delta_input)
                negative[key]=get(negative,key,0)+mu_value*antipode_value*delta_value
            end
        end
    end
    unit=Dict{Tuple,ComplexF64}((1,)=>1)
    counit=Dict{Tuple,ComplexF64}((1,)=>1,(2,)=>1)
    right_integral=Dict{Tuple,ComplexF64}((4,)=>1)
    right_cointegral=Dict{Tuple,ComplexF64}((3,)=>1,(4,)=>1)
    left_integral=Dict{Tuple,ComplexF64}((3,)=>-1)
    left_cointegral=Dict{Tuple,ComplexF64}((3,)=>-1,(4,)=>1)
    FiniteHopfAlgebra(4,_sparse_data(multiplication),_sparse_data(unit),
                      _sparse_data(coproduct),_sparse_data(counit),
                      _sparse_data(positive),_sparse_data(negative),_sparse_data(antipode),
                      _sparse_data(right_integral),_sparse_data(right_cointegral),
                      _sparse_data(left_integral),_sparse_data(left_cointegral))
end

function exterior_algebra(generator_count::Integer)
    generator_count >= 0 || throw(ArgumentError("generator count must be nonnegative"))
    dimension=2^generator_count
    multiplication=Dict{Tuple,ComplexF64}()
    for left in 0:dimension-1,right in 0:dimension-1
        left & right == 0 || continue
        inversions=sum((1 for a in 0:generator_count-1,b in 0:generator_count-1
                        if (left >> a)&1 == 1 && (right >> b)&1 == 1 && a > b);init=0)
        multiplication[(left+1,right+1,(left|right)+1)]=(-1)^inversions
    end
    coproduct=Dict{Tuple,ComplexF64}()
    for source in 0:dimension-1
        left=source
        while true
            right=xor(source,left)
            inversions=sum((1 for a in 0:generator_count-1,b in 0:generator_count-1
                            if (left >> a)&1 == 1 && (right >> b)&1 == 1 && a > b);init=0)
            coproduct[(source+1,left+1,right+1)]=(-1)^inversions
            left == 0 && break
            left=(left-1)&source
        end
    end
    antipode=Dict{Tuple,ComplexF64}((mask+1,mask+1)=>(-1)^count_ones(mask)
                                    for mask in 0:dimension-1)
    positive=Dict{Tuple,ComplexF64}(); negative=Dict{Tuple,ComplexF64}()
    for ((delta_input,left,right),delta_value) in coproduct
        for ((mu_input,middle,output),mu_value) in multiplication
            middle == left || continue
            key=(output,right,mu_input,delta_input)
            positive[key]=get(positive,key,0)+mu_value*delta_value
        end
        for ((antipode_output,antipode_input),antipode_value) in antipode
            antipode_input == left || continue
            for ((mu_input,middle,output),mu_value) in multiplication
                middle == antipode_output || continue
                key=(output,right,mu_input,delta_input)
                negative[key]=get(negative,key,0)+mu_value*antipode_value*delta_value
            end
        end
    end
    unit=Dict{Tuple,ComplexF64}((1,)=>1)
    counit=Dict{Tuple,ComplexF64}((1,)=>1)
    top=(dimension,)=>1
    integral=Dict{Tuple,ComplexF64}(top)
    FiniteHopfAlgebra(dimension,_sparse_data(multiplication),_sparse_data(unit),
                      _sparse_data(coproduct),_sparse_data(counit),
                      _sparse_data(positive),_sparse_data(negative),_sparse_data(antipode),
                      _sparse_data(integral),_sparse_data(integral),
                      _sparse_data(integral),_sparse_data(integral))
end

function _factor_with_variables(factor::SparseFactor,variables)
    source_variables=Int.(variables)
    unique_variables=unique(source_variables)
    positions=[[index for (index,variable) in enumerate(source_variables) if variable==target]
               for target in unique_variables]
    terms=Dict{Tuple,ComplexF64}()
    for (key,value) in factor.terms
        all(all(key[index]==key[first(group)] for index in group) for group in positions) || continue
        reduced=Tuple(key[first(group)] for group in positions)
        terms[reduced]=get(terms,reduced,0)+value
    end
    filter!(pair -> abs(last(pair))>1e-12,terms)
    SparseFactor(unique_variables,terms)
end

function _multiply_factors(left::SparseFactor,right::SparseFactor)
    variables=unique(vcat(left.variables,right.variables))
    left_positions=[findfirst(==(variable),variables) for variable in left.variables]
    right_positions=[findfirst(==(variable),variables) for variable in right.variables]
    shared=intersect(left.variables,right.variables)
    left_shared=[findfirst(==(variable),left.variables) for variable in shared]
    right_shared=[findfirst(==(variable),right.variables) for variable in shared]
    terms=Dict{Tuple,ComplexF64}()
    right_index=Dict{Tuple,Vector{Tuple{Tuple,ComplexF64}}}()
    for (right_key,right_value) in right.terms
        join_key=Tuple(right_key[position] for position in right_shared)
        push!(get!(right_index,join_key,Tuple{Tuple,ComplexF64}[]),(right_key,right_value))
    end
    for (left_key,left_value) in left.terms
        join_key=Tuple(left_key[position] for position in left_shared)
        for (right_key,right_value) in get(right_index,join_key,Tuple{Tuple,ComplexF64}[])
        key=fill(0,length(variables))
        for (position,value) in zip(left_positions,left_key); key[position]=value; end
        for (position,value) in zip(right_positions,right_key); key[position]=value; end
        tuple=Tuple(key)
        terms[tuple]=get(terms,tuple,0)+left_value*right_value
        end
    end
    filter!(pair -> abs(last(pair))>1e-12,terms)
    SparseFactor(variables,terms)
end

function _sum_variable(factor::SparseFactor,variable::Int)
    position=findfirst(==(variable),factor.variables)
    position === nothing && return factor
    variables=factor.variables[setdiff(eachindex(factor.variables),position)]
    terms=Dict{Tuple,ComplexF64}()
    for (key,value) in factor.terms
        reduced=Tuple(key[index] for index in eachindex(key) if index != position)
        terms[reduced]=get(terms,reduced,0)+value
    end
    SparseFactor(variables,terms)
end

function _contract(factors::Vector{SparseFactor})
    while any(!isempty(factor.variables) for factor in factors)
        variables=unique(vcat((factor.variables for factor in factors)...))
        costs=[prod(max(length(factor.terms),1)
                    for factor in factors if variable in factor.variables) for variable in variables]
        variable=variables[argmin(costs)]
        selected=findall(factor -> variable in factor.variables,factors)
        selected_factors=sort(factors[selected];by=factor -> length(factor.terms))
        combined=reduce(_multiply_factors,selected_factors)
        for index in reverse(selected); deleteat!(factors,index); end
        push!(factors,_sum_variable(combined,variable))
    end
    final=reduce(_multiply_factors,factors)
    get(final.terms,(),0.0+0im)
end

function hopf_invariant(v::VirtualLink,algebra::FiniteHopfAlgebra)
    edges=real_edge_table(v)
    orientation=last(real_gauss_code(v))
    factors=SparseFactor[]
    zero_map=Dict{Int,Int}()
    for edge in edges
        labels=4 .* abs.([edge.crossing[1],edge.crossing[2]]) .+
               ([edge.crossing[1] > 0,edge.crossing[2] > 0]) .- [1,3]
        exponent=v.is_weighted ? round(Int,2edge.weight) : 0
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
    _contract(factors)
end
