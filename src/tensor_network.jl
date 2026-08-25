"""Attach Einstein-index labels to a stored sparse tensor."""
label_tensor(tensor::SparseFactor,labels)=_coalesce_repeated(
    SparseFactor(Int.(collect(labels)),tensor.terms))

_coefficient_iszero(value::AbstractFloat)=abs(value)<=1e-12
_coefficient_iszero(value::Complex{<:AbstractFloat})=abs(value)<=1e-12
_coefficient_iszero(value)=iszero(value)

function _coalesce_repeated(factor::SparseFactor{T}) where {T}
    variables=unique(factor.variables)
    length(variables)==length(factor.variables) && return factor
    positions=[[index for (index,value) in enumerate(factor.variables) if value==variable]
               for variable in variables]
    terms=Dict{Tuple,T}()
    for (key,value) in factor.terms
        all(all(key[index]==key[first(group)] for index in group) for group in positions) || continue
        reduced=Tuple(key[first(group)] for group in positions)
        terms[reduced]=get(terms,reduced,zero(T))+value
    end
    SparseFactor(variables,terms)
end

function _contract_to(factors::AbstractVector{<:SparseFactor},outputs::Vector{Int})
    isempty(factors) && return SparseFactor(copy(outputs),Dict{Tuple,ComplexF64}(() => 1))
    factors=_coalesce_repeated.(factors)
    present=unique(vcat((factor.variables for factor in factors)...))
    all(output in present for output in outputs) ||
        throw(ArgumentError("every output index must occur in the tensor network"))
    eliminate=setdiff(present,outputs)
    while !isempty(eliminate)
        costs=[prod(max(length(factor.terms),1)
                    for factor in factors if variable in factor.variables)
               for variable in eliminate]
        variable=eliminate[argmin(costs)]
        selected=findall(factor -> variable in factor.variables,factors)
        combined=reduce(_multiply_factors,factors[selected])
        for index in reverse(selected); deleteat!(factors,index); end
        push!(factors,_sum_variable(combined,variable))
        deleteat!(eliminate,findfirst(==(variable),eliminate))
    end
    combined=reduce(_multiply_factors,factors)
    extra=setdiff(combined.variables,outputs)
    for variable in extra; combined=_sum_variable(combined,variable); end
    permutation=[findfirst(==(variable),combined.variables) for variable in outputs]
    T=eltype(values(combined.terms))
    terms=Dict{Tuple,T}()
    for (key,value) in combined.terms
        ordered=Tuple(key[index] for index in permutation)
        terms[ordered]=get(terms,ordered,zero(T))+value
    end
    SparseFactor(copy(outputs),terms)
end

"""Contract labeled sparse tensors and retain `outputs` in the requested order."""
contract_network(factors,outputs=Int[])=_contract_to(collect(factors),Int.(collect(outputs)))

function dense_tensor(factor::SparseFactor{T},dimension::Integer) where {T}
    dimensions=ntuple(_->Int(dimension),length(factor.variables))
    result=zeros(T,dimensions)
    for (key,value) in factor.terms; result[key...]=value; end
    result
end

function _matrix_power_factor(factor::SparseFactor,power::Integer,dimension::Int)
    matrix=dense_tensor(label_tensor(factor,[1,2]),dimension)
    powered=matrix^power
    T=eltype(powered)
    terms=Dict{Tuple,T}()
    for output in 1:dimension,input in 1:dimension
        !_coefficient_iszero(powered[output,input]) && (terms[(output,input)]=powered[output,input])
    end
    SparseFactor(Int[],terms)
end

function hopf_tensors(algebra::FiniteHopfAlgebra)
    Dict{String,SparseFactor}(
        "Mu"=>algebra.multiplication,"De"=>algebra.coproduct,
        "Et"=>algebra.unit,"Ep"=>algebra.counit,
        "An"=>algebra.antipode,"Ir"=>algebra.right_integral,
        "Il"=>algebra.left_integral,"Cr"=>algebra.right_cointegral,
        "Cl"=>algebra.left_cointegral,"Tp"=>algebra.positive_vertex,
        "Tm"=>algebra.negative_vertex,
        "Un"=>_matrix_power_factor(algebra.antipode,0,algebra.dimension))
end

"""Tensor convention used by `C250406TensorNetwork.mlx` before the 2025-08-16 antipode transpose."""
function legacy_hopf_tensors(algebra::FiniteHopfAlgebra)
    tensors=hopf_tensors(algebra)
    tensors["An"]=_sparse_data(Dict((input,output)=>value
        for ((output,input),value) in algebra.antipode.terms))
    tensors["Tp"]=evaluate_tensor_expression(
        "De{1,4,0}Mu{3,4,2}",tensors,[1,0,3,2];dimension=algebra.dimension)
    tensors["Tm"]=evaluate_tensor_expression(
        "De{1,4,0}An{4,5}Mu{3,5,2}",tensors,[0,1,2,3];dimension=algebra.dimension)
    tensors["Un"]=_matrix_power_factor(tensors["An"],0,algebra.dimension)
    tensors
end

_source_column(text::AbstractString,index::Integer)=
    index<=firstindex(text) ? 1 : length(SubString(text,firstindex(text),prevind(text,index)))+1

function _first_nonspace(text::AbstractString,start::Integer,stop::Integer)
    index=start
    while index<=stop && isspace(text[index])
        index=nextind(text,index)
    end
    index
end

"""Validate MATLAB-compatible tensor notation and return all diagnostics."""
function validate_tensor_expression(expression::AbstractString,outputs=Int[];
                                    tensors=nothing,index_rule::Symbol=:matlab)
    index_rule in (:matlab,:network) ||
        throw(ArgumentError("index_rule must be :matlab or :network"))
    token=r"(?:\(([A-Za-z_]\w*)(?:\^(-?\d+))?\)|([A-Za-z_]\w*))\{([0-9]+(?:,[0-9]+)*)?\}"
    matches=collect(eachmatch(token,expression))
    problems=String[]
    isempty(matches) && push!(problems,"no tensor factors found")

    cursor=1
    labels=Int[]
    for match in matches
        if match.offset>cursor
            residue=SubString(expression,cursor,prevind(expression,match.offset))
            if !isempty(strip(residue))
                bad=_first_nonspace(expression,cursor,prevind(expression,match.offset))
                push!(problems,"invalid syntax at column $(_source_column(expression,bad)) near `$(strip(residue))`")
            end
        end
        cursor=match.offset+ncodeunits(match.match)

        name=match.captures[1]===nothing ? match.captures[3] : match.captures[1]
        power=match.captures[2]
        legs=match.captures[4]===nothing ? Int[] : parse.(Int,split(match.captures[4],','))
        append!(labels,legs)
        power!==nothing && length(legs)!=2 &&
            push!(problems,"tensor power $name^$power at column $(_source_column(expression,match.offset)) must have exactly 2 legs")

        if tensors!==nothing
            if !haskey(tensors,name)
                push!(problems,"unknown tensor label $name at column $(_source_column(expression,match.offset))")
            else
                tensor=tensors[name]
                tensor_rank=isempty(tensor.terms) ? length(tensor.variables) :
                    length(first(keys(tensor.terms)))
                tensor_rank==length(legs) ||
                    push!(problems,"tensor $name at column $(_source_column(expression,match.offset)) has rank $tensor_rank, got $(length(legs)) legs")
            end
        end
    end
    if cursor<=ncodeunits(expression)
        residue=SubString(expression,cursor,lastindex(expression))
        if !isempty(strip(residue))
            bad=_first_nonspace(expression,cursor,lastindex(expression))
            push!(problems,"invalid syntax at column $(_source_column(expression,bad)) near `$(strip(residue))`")
        end
    end

    output_labels=Int.(collect(outputs))
    length(unique(output_labels))==length(output_labels) ||
        push!(problems,"ord contains duplicate index labels")
    if index_rule==:matlab
        counts=Dict{Int,Int}()
        for label in Iterators.flatten((labels,output_labels))
            counts[label]=get(counts,label,0)+1
        end
        for label in sort!(collect(keys(counts)))
            count=counts[label]
            count==2 || push!(problems,"index $label appears $count times; every index must appear exactly twice across the expression and ord")
        end
    else
        for label in output_labels
            label in labels || push!(problems,"output index $label does not occur in the expression")
        end
    end
    problems
end

"""Evaluate MATLAB `makeTensorExpression` notation as a sparse Einstein network."""
function evaluate_tensor_expression(expression::AbstractString,
                                    tensors::AbstractDict,
                                    outputs=Int[];
                                    dimension::Union{Nothing,Integer}=nothing,
                                    index_rule::Symbol=:network)
    problems=validate_tensor_expression(expression,outputs;tensors,index_rule)
    isempty(problems) || throw(ArgumentError(
        "invalid tensor expression:\n- " * join(problems,"\n- ")))
    pattern=r"(?:\(([A-Za-z_]\w*)(?:\^(-?\d+))?\)|([A-Za-z_]\w*))\{([^}]*)\}"
    factors=SparseFactor[]
    for match in eachmatch(pattern,expression)
        name=match.captures[1]===nothing ? match.captures[3] : match.captures[1]
        haskey(tensors,name) || throw(KeyError(name))
        tensor=tensors[name]
        if match.captures[2]!==nothing
            dimension===nothing && throw(ArgumentError("dimension is required for tensor powers"))
            tensor=_matrix_power_factor(tensor,parse(Int,match.captures[2]),Int(dimension))
        end
        labels=isempty(strip(match.captures[4])) ? Int[] :
            parse.(Int,split(match.captures[4],','))
        length(labels)==length(first(keys(tensor.terms))) ||
            throw(DimensionMismatch("$name has $(length(first(keys(tensor.terms)))) indices, got $(length(labels))"))
        push!(factors,label_tensor(tensor,labels))
    end
    isempty(factors) && throw(ArgumentError("no tensor factors found"))
    contract_network(factors,outputs)
end

function _factors_approx(left::SparseFactor,right::SparseFactor;atol=1e-10)
    left.variables==right.variables || return false
    keys_union=union(keys(left.terms),keys(right.terms))
    all(abs(get(left.terms,key,0)-get(right.terms,key,0))<=atol for key in keys_union)
end

"""Verify the unit, counit, associativity, coassociativity, bialgebra and antipode axioms."""
function verify_hopf(algebra::FiniteHopfAlgebra;atol=1e-10,degrees=nothing)
    t=hopf_tensors(algebra)
    checks=Dict{String,Bool}()
    identity=_matrix_power_factor(algebra.antipode,0,algebra.dimension)
    checks["left_unit"]=_factors_approx(
        evaluate_tensor_expression("Et{1}Mu{1,2,3}",t,[2,3]),label_tensor(identity,[2,3]);atol)
    checks["right_unit"]=_factors_approx(
        evaluate_tensor_expression("Et{2}Mu{1,2,3}",t,[1,3]),label_tensor(identity,[1,3]);atol)
    checks["left_counit"]=_factors_approx(
        evaluate_tensor_expression("Ep{2}De{1,2,3}",t,[1,3]),label_tensor(identity,[1,3]);atol)
    checks["right_counit"]=_factors_approx(
        evaluate_tensor_expression("Ep{3}De{1,2,3}",t,[1,2]),label_tensor(identity,[1,2]);atol)
    checks["associativity"]=_factors_approx(
        evaluate_tensor_expression("Mu{1,2,4}Mu{4,3,5}",t,[1,2,3,5]),
        evaluate_tensor_expression("Mu{2,3,4}Mu{1,4,5}",t,[1,2,3,5]);atol)
    checks["coassociativity"]=_factors_approx(
        evaluate_tensor_expression("De{1,2,4}De{4,3,5}",t,[1,2,3,5]),
        evaluate_tensor_expression("De{1,4,5}De{4,2,3}",t,[1,2,3,5]);atol)
    unit_counit=evaluate_tensor_expression("Ep{1}Et{5}",t,[1,5])
    checks["left_antipode"]=_factors_approx(
        evaluate_tensor_expression("De{1,2,3}An{4,2}Mu{4,3,5}",t,[1,5]),unit_counit;atol)
    checks["right_antipode"]=_factors_approx(
        evaluate_tensor_expression("De{1,2,3}An{4,3}Mu{2,4,5}",t,[1,5]),unit_counit;atol)
    bialgebra_right=if degrees===nothing
        evaluate_tensor_expression("De{1,6,7}De{2,8,9}Mu{6,8,4}Mu{7,9,5}",t,[1,2,4,5])
    else
        length(degrees)==algebra.dimension || throw(DimensionMismatch("one degree per basis element required"))
        koszul=Dict{Tuple,ComplexF64}((left,right)=>(-1)^(degrees[left]*degrees[right])
                                     for left in 1:algebra.dimension,right in 1:algebra.dimension)
        graded=copy(t); graded["Koszul"]=_sparse_data(koszul)
        evaluate_tensor_expression(
            "De{1,6,7}De{2,8,9}Koszul{7,8}Mu{6,8,4}Mu{7,9,5}",graded,[1,2,4,5])
    end
    checks["bialgebra"]=_factors_approx(
        evaluate_tensor_expression("Mu{1,2,3}De{3,4,5}",t,[1,2,4,5]),
        bialgebra_right;atol)
    checks
end
