"""
    SparseFactor(variables, terms::AbstractDict) -> SparseFactor

Store coefficients and index labels of a sparse tensor.

# Arguments
`variables` supplies integer labels convertible to `Vector{Int}`; `terms`
maps tuple keys to coefficients of type `T`. Failed conversions propagate.
Construction does not check key lengths, coordinate ranges, repeated labels
or zero coefficients. The caller supplies coherent tensor data for the
operations that consume the result.

# Returns
`SparseFactor{T}` with `variables::Vector{Int}` and `terms::Dict{Tuple,T}`.
These are the label and coefficient containers of the
[Algebra tensor API](../specs/algebra.md#definitions).

# Guarantees
1. Coefficient type and values are preserved; construction does not convert
   exact coefficients to floating point or discard stored zeros. [proved/agent: by construction]
2. With a consistent choice of basis and index labels, the coefficient data
   represent the supplied tensor in that basis. [proved/agent: by construction]
3. Supplying exactly the declared field types retains those containers;
   the general dictionary method collects labels and entries into new containers. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
No tensor identities or contraction rules are checked. In particular,
unlabelled structure tensors must acquire appropriate labels before being
used as network factors; construction does not infer tensor variance.

# Sources
cf:[[Kup91#^kup91-blk-0024]] basis coordinates of tensors; sparse storage and integer network labels are repository representations.

# Examples
```jldoctest
julia> a = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>1//2, (1,2)=>1//3, (2,1)=>2//5));

julia> b = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>3//7, (2,1)=>5//11));

julia> r = evaluate_tensor_expression("A{1,2}B{2,3}", Dict("A"=>a, "B"=>b), [1,3]);

julia> r isa SparseFactor{Rational{Int}}
true
```

See also [`label_tensor`](@ref), [`dense_tensor`](@ref).

Role: construct
"""
struct SparseFactor{T}
    variables::Vector{Int}
    terms::Dict{Tuple,T}
end

SparseFactor(variables,terms::AbstractDict{<:Tuple,T}) where {T}=
    SparseFactor{T}(Int.(collect(variables)),
                    Dict{Tuple,T}(Tuple(key)=>value for (key,value) in terms))

"""
    FiniteHopfAlgebra(dimension, multiplication, unit, coproduct, counit, positive_vertex, negative_vertex, antipode, right_integral, right_cointegral, left_integral, left_cointegral) -> FiniteHopfAlgebra

Store finite-dimensional Hopf-algebra tensor data without validating the axioms.

# Arguments
`dimension` must be convertible to `Int`; the remaining eleven arguments
must be `SparseFactor` values. Failed field conversions and unsupported
argument lists propagate their Julia exceptions. The caller is responsible
for a consistent dimension, basis, tensor shapes and algebraic identities.

# Returns
`FiniteHopfAlgebra{Nothing}` for the unparameterized constructor.
`dimension` is the supplied vector-space dimension; `multiplication`, `unit`,
`coproduct`, `counit`, and `antipode` store its structural maps.
`positive_vertex` and `negative_vertex` store crossing tensors;
`right_integral`, `right_cointegral`, `left_integral`, and `left_cointegral`
store the supplied distinguished tensors. Their API belongs to
[Algebra](../specs/algebra.md#definitions).

# Guarantees
1. The fields retain the supplied values after conversion; tensor containers
   are shared rather than copied. [proved/agent: by construction]
2. When the supplied structural maps satisfy the Hopf axioms in their
   declared basis, the record represents that finite-dimensional Hopf algebra.
   An instance of the type alone is not evidence for those hypotheses. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Construction checks neither the Hopf axioms nor integral normalization,
involutivity, unimodularity or compatibility of the crossing tensors.
The coefficient types may differ between fields. The record does not assert
that a contraction gives a topological invariant.

# Sources
cf:[[Kup91#^kup91-sec-0038]] Hopf algebras and their structural tensors; the storage constructor also permits noninvolutory and unchecked input data.

# Examples
```jldoctest
julia> h = sweedler_algebra();

julia> hcopy = FiniteHopfAlgebra(h.dimension, h.multiplication, h.unit, h.coproduct, h.counit, h.positive_vertex, h.negative_vertex, h.antipode, h.right_integral, h.right_cointegral, h.left_integral, h.left_cointegral);

julia> all(values(verify_hopf(hcopy)))
true
```

See also [`verify_hopf`](@ref), [`hopf_tensors`](@ref).

Role: convert
"""
struct FiniteHopfAlgebra{Kind}
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

# Untagged user-built algebras retain the generic tensor contraction path.
FiniteHopfAlgebra(args...) = FiniteHopfAlgebra{Nothing}(args...)
struct TaftKind{N,M,Style} end
struct GroupAlgebraKind end

_sparse_data(terms)=SparseFactor(Int[],terms)

_qbinomial(q,n,k)=begin
    k < 0 || k > n ? 0.0+0im :
    prod((1-q^(n-j))/(1-q^(j+1)) for j in 0:k-1;init=1.0+0im)
end

const _UQSL2_BOREL_CACHE=Dict{Tuple{Any,Symbol},FiniteHopfAlgebra}()

"""Build one canonical small Borel algebra instance (uncached internal path)."""
function _build_uqsl2_borel_small(ratio::Rational;style=:L)
    numerator(ratio) > 0 && numerator(ratio) < denominator(ratio) ||
        throw(ArgumentError("ratio must lie in (0,1)"))
    M,N=Int(numerator(ratio)),Int(denominator(ratio))
    q=N == 2 || N == 4 ? ComplexF64(im)^(M*4÷N) : exp(2pi*im*M/N)
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
    FiniteHopfAlgebra{TaftKind{N,M,style}}(dimension,_sparse_data(multiplication),_sparse_data(unit),
                      _sparse_data(coproduct),_sparse_data(counit),
                      _sparse_data(positive),_sparse_data(negative),_sparse_data(antipode),
                      _sparse_data(right_integral),_sparse_data(right_cointegral),
                      _sparse_data(left_integral),_sparse_data(left_cointegral))
end

"""
    uqsl2_borel_small(ratio::Rational; style=:L) -> FiniteHopfAlgebra

Return the library's small Borel tensor model at the rational root parameter.

# Arguments
- `ratio` must lie strictly between zero and one; otherwise `ArgumentError`
  reports that interval. Its numerator and denominator must fit `Int`, and
  the resulting tensor dimensions must fit available index and memory ranges.
  Conversion, allocation and numerical linear-algebra errors propagate.
- `style` must be `:L` or `:R`; otherwise `ArgumentError` identifies the
  accepted styles. These are the repository's two small-Borel presentations
  in the [Algebra API](../specs/algebra.md#definitions).

# Returns
A `FiniteHopfAlgebra` with dimension `denominator(ratio)^2`, `ComplexF64`
tensor coefficients, and a type tag recording the root parameters and style.
Fields have the roles described by [`FiniteHopfAlgebra`](@ref).

# Guarantees
1. Within the representable range, the result supplies the fixed small-Borel
   structural and crossing tensors for the requested parameter and style on
   that complex coordinate space. [proved/agent: by construction]
2. The result's mutable tensor dictionaries are isolated from results of
   other calls. Mutating a returned tensor does not change later results
   for the same arguments. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The coefficients are numerical, not exact cyclotomic values or a proof of
Hopf identities. The parameter and basis naming do not certify a translation
to every literature convention for a small quantum group. The type tag alone
does not guarantee eligibility for an exact invariant route after tensors
have been modified; consult [`hopf_invariant`](@ref).

# Sources
cf:lib/topology.jl/src/hopf_invariant.jl, _build_uqsl2_borel_small — the repository parameterized tensor model and its two styles.

# Examples
```jldoctest
julia> h = uqsl2_borel_small(1//2; style=:L); other = uqsl2_borel_small(1//2; style=:L); other.antipode.terms[(1,1)] = 999;

julia> (h !== other, h.antipode.terms[(1,1)] != 999)
(true, true)
```

See also [`verify_hopf`](@ref), [`hopf_invariant`](@ref).

Role: construct
"""
function uqsl2_borel_small(ratio::Rational;style=:L)
    key=(ratio,style)
    prototype=get!(_UQSL2_BOREL_CACHE,key) do
        _build_uqsl2_borel_small(ratio;style=style)
    end
    deepcopy(prototype)
end

"""
    sweedler_algebra() -> FiniteHopfAlgebra

Return the library's four-dimensional Sweedler Hopf-algebra model.

# Arguments
No arguments are required.

# Returns
`FiniteHopfAlgebra{Nothing}` of dimension four with `ComplexF64` tensor
coefficients and the stored structural, crossing and distinguished tensors.
Its field roles belong to [`FiniteHopfAlgebra`](@ref).

# Guarantees
1. The result is the fixed finite-dimensional Hopf model on a complex
   four-dimensional coordinate space used by the
   [Algebra API](../specs/algebra.md#definitions). [proved/agent: by construction]
2. Its antipode is not involutory. Each call returns independently allocated
   sparse tensor data, so mutations do not alter another call's result. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The current specification does not fix a translation from these basis and
integral-field names to a particular literature presentation of Sweedler's
algebra. The model is not unimodular or an admissible involutory input merely
because it has type `FiniteHopfAlgebra`. No manifold invariant is computed.

# Sources
cf:[[Kup91#^kup91-sec-0038]] finite-dimensional Hopf structural tensors; the named noninvolutory model and its basis require a separate convention identification.

# Examples
```jldoctest
julia> h = sweedler_algebra();

julia> (all(values(verify_hopf(h))), is_involutory(h), is_unimodular(h))
(true, false, false)
```

See also [`verify_hopf`](@ref), [`is_involutory`](@ref).

Role: construct
"""
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

"""
    exterior_algebra(generator_count::Integer) -> FiniteHopfAlgebra

Return the library's finite exterior-algebra model with its graded Hopf tensors.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

`generator_count` must be nonnegative; otherwise `ArgumentError` states that
requirement. The dimension `2^generator_count` must be representable and
small enough to allocate the tensor data. No explicit overflow or resource
bound is enforced; arithmetic and allocation failures propagate.

# Returns
`FiniteHopfAlgebra{Nothing}` with dimension `2^generator_count` and
`ComplexF64` coefficients. The fields have the roles described by
[`FiniteHopfAlgebra`](@ref), with exterior grading for the bialgebra identity.

# Guarantees
1. Within the representable range, the coordinate dimension is the dimension
   of the exterior algebra on the requested number of generators. [proved/agent: by construction]
2. The stored multiplication, coproduct and antipode give the repository's
   exterior model; the antipode is involutory. [proved/agent: by construction]
3. Each call returns independently allocated sparse tensors. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The ordinary, ungraded bialgebra check is not the graded Hopf check. No
translation to another literature basis or sign convention is certified.
No ordinary Hopf or topological-invariance claim is inferred by ignoring
the grading.

# Sources
cf:[[Kup91#^kup91-sec-0038]] structural tensors for Hopf algebras; the exterior model uses an additional grading not specified by this ordinary Hopf-algebra source.

# Examples
```jldoctest
julia> h = exterior_algebra(2);

julia> (h.dimension, all(values(verify_hopf(h; degrees=count_ones.(0:3)))))
(4, true)
```

See also [`verify_hopf`](@ref), [`is_involutory`](@ref).

Role: construct
"""
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

function _factor_with_variables(factor::SparseFactor{T},variables) where {T}
    source_variables=Int.(variables)
    unique_variables=unique(source_variables)
    positions=[[index for (index,variable) in enumerate(source_variables) if variable==target]
               for target in unique_variables]
    terms=Dict{Tuple,T}()
    for (key,value) in factor.terms
        all(all(key[index]==key[first(group)] for index in group) for group in positions) || continue
        reduced=Tuple(key[first(group)] for group in positions)
        terms[reduced]=get(terms,reduced,zero(T))+value
    end
    filter!(pair -> !_coefficient_iszero(last(pair)),terms)
    SparseFactor(unique_variables,terms)
end

function _multiply_factors(left::SparseFactor{L},right::SparseFactor{R}) where {L,R}
    T=promote_type(L,R)
    variables=unique(vcat(left.variables,right.variables))
    left_positions=[findfirst(==(variable),variables) for variable in left.variables]
    right_positions=[findfirst(==(variable),variables) for variable in right.variables]
    shared=intersect(left.variables,right.variables)
    left_shared=[findfirst(==(variable),left.variables) for variable in shared]
    right_shared=[findfirst(==(variable),right.variables) for variable in shared]
    terms=Dict{Tuple,T}()
    right_index=Dict{Tuple,Vector{Tuple{Tuple,R}}}()
    for (right_key,right_value) in right.terms
        join_key=Tuple(right_key[position] for position in right_shared)
        push!(get!(right_index,join_key,Tuple{Tuple,R}[]),(right_key,right_value))
    end
    for (left_key,left_value) in left.terms
        join_key=Tuple(left_key[position] for position in left_shared)
        for (right_key,right_value) in get(right_index,join_key,Tuple{Tuple,R}[])
        key=fill(0,length(variables))
        for (position,value) in zip(left_positions,left_key); key[position]=value; end
        for (position,value) in zip(right_positions,right_key); key[position]=value; end
        tuple=Tuple(key)
        terms[tuple]=get(terms,tuple,zero(T))+left_value*right_value
        end
    end
    filter!(pair -> !_coefficient_iszero(last(pair)),terms)
    SparseFactor(variables,terms)
end

function _sum_variable(factor::SparseFactor{T},variable::Int) where {T}
    position=findfirst(==(variable),factor.variables)
    position === nothing && return factor
    variables=factor.variables[setdiff(eachindex(factor.variables),position)]
    terms=Dict{Tuple,T}()
    for (key,value) in factor.terms
        reduced=Tuple(key[index] for index in eachindex(key) if index != position)
        terms[reduced]=get(terms,reduced,zero(T))+value
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

