"""
    FiniteAlgebra(dimension, multiplication, unit, name) -> FiniteAlgebra

Store multiplication and unit data for a finite-dimensional algebra in a chosen basis.

# Arguments
`dimension::Int`, `multiplication::SparseFactor`, `unit::SparseFactor`, and
`name::String`, or values convertible to those field types. Conversion
failures propagate. No dimension, tensor-shape, unit or associativity check
is performed; the caller supplies coherent algebra data.

# Returns
`FiniteAlgebra` with the supplied dimension, multiplication tensor, unit
tensor and descriptive name. The representation-independent algebra API is
owned by [Algebra](../specs/algebra.md#definitions).

# Guarantees
1. Fields retain the supplied values after conversion; supplied sparse
   factors are shared with the result. [proved/agent: by construction]
2. If the supplied tensors are the multiplication and unit of an algebra,
   the result stores that algebra in the supplied basis. Construction does
   not certify the hypothesis. [proved/agent: by construction]
3. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The name does not identify an isomorphism class. Invalid dimensions and
inconsistent structure constants are not rejected. No coalgebra, antipode
or integral data are added.

# Sources
cf:[[Kup91#^kup91-blk-0039]] multiplication tensor of a finite-dimensional algebra; the container and its validation boundary are repository-specific.

# Examples
```jldoctest
julia> h = uqsl2_borel_small(1//3; style=:L); d = heisenberg_double(h);

julia> a = FiniteAlgebra(d.dimension, d.multiplication, d.unit, d.name);

julia> a.dimension
81
```

See also [`basis_element`](@ref), [`FiniteHopfAlgebra`](@ref).

Role: convert
"""
struct FiniteAlgebra
    dimension::Int
    multiplication::SparseFactor
    unit::SparseFactor
    name::String
end

"""
    AlgebraElement(algebra::FiniteAlgebra, coefficients_or_index) -> AlgebraElement

Return a coordinate element of the supplied finite-dimensional algebra.

# Arguments
- `algebra` supplies the basis and algebra operations.
- A coefficient vector must be convertible to `Vector{ComplexF64}` and
  should have `algebra.dimension` entries. Conversion failures propagate;
  the vector's length is not checked by this constructor.
- An integer `index` selects a basis element and must lie in
  `1:algebra.dimension`; otherwise a `BoundsError` contains that range and
  the offending index.

# Returns
`AlgebraElement` with `algebra` the supplied algebra and
`coefficients::Vector{ComplexF64}` its coordinates. This is part of the
[Algebra API](../specs/algebra.md#definitions).

# Guarantees
1. The vector form preserves the supplied coordinates after conversion.
   A supplied `Vector{ComplexF64}` is retained without copying. [proved/agent: by construction]
2. The index form returns the indicated basis vector: its indexed
   coefficient is one and all other coefficients are zero, in a fresh vector. [proved/agent: by construction]
3. For coherent algebra data and a correctly sized vector, the result is
   the corresponding element of that algebra, not an element of its dual. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The algebra axioms and the vector length are not certified. Coefficients
are stored in floating-point complex arithmetic, even when supplied exactly.

# Sources
cf:[[Kup91#^kup91-blk-0024]] coordinates of a vector in a given finite-dimensional basis; complex floating-point storage is a repository choice.

# Examples
```jldoctest
julia> h = cyclic_group_algebra(5); b = basis_element(h, 2);

julia> a = AlgebraElement(b.algebra, 2);

julia> (a^5).coefficients ≈ basis_element(h, 1).coefficients
true
```

See also [`basis_element`](@ref), [`hopf_pairing`](@ref).

Role: convert
"""
struct AlgebraElement
    algebra::FiniteAlgebra
    coefficients::Vector{ComplexF64}
end

function AlgebraElement(algebra::FiniteAlgebra,index::Integer)
    1<=index<=algebra.dimension || throw(BoundsError(1:algebra.dimension,index))
    coefficients=zeros(ComplexF64,algebra.dimension); coefficients[index]=1
    AlgebraElement(algebra,coefficients)
end

function _finite_algebra(algebra::FiniteHopfAlgebra,name="Hopf algebra")
    FiniteAlgebra(algebra.dimension,algebra.multiplication,algebra.unit,name)
end

"""
    basis_element(algebra::FiniteAlgebra, index::Integer) -> AlgebraElement

Return the indicated basis element of a finite algebra or the algebra underlying a Hopf algebra.

# Arguments
`index` must lie in `1:algebra.dimension`; otherwise `BoundsError` reports
the admissible range and index. The supplied algebra must have coherent
dimension and multiplication data for subsequent algebra operations.

# Returns
An `AlgebraElement` with a fresh `Vector{ComplexF64}` of coordinates.
For a Hopf input, its `algebra` field is a `FiniteAlgebra` with the same
dimension, multiplication and unit, named `"Hopf algebra"`.

# Guarantees
1. The coordinate vector has `algebra.dimension` entries, one at `index`
   and zero elsewhere. It represents that basis vector of the input algebra. [proved/agent: by construction]
2. A Hopf input's multiplication and unit factors are shared with the
   underlying algebra in the result. The input is unchanged. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
This does not validate the algebra or choose a canonical basis. The result
stores no coproduct or antipode. Basis order belongs to the particular
algebra constructor; see the [Algebra scope](../specs/algebra.md#definitions).

# Sources
cf:[[Kup91#^kup91-blk-0024]] coordinates in a chosen basis; extracting an underlying algebra is the repository API boundary.

# Examples
```jldoctest
julia> h = cyclic_group_algebra(5);

julia> (basis_element(h, 2)^5).coefficients ≈ basis_element(h, 1).coefficients
true
```

See also [`AlgebraElement`](@ref).

Role: access
"""
basis_element(algebra::FiniteAlgebra,index::Integer)=AlgebraElement(algebra,index)
"""
    basis_element(algebra::FiniteHopfAlgebra, index::Integer) -> AlgebraElement

Return the indicated basis element of the algebra underlying a Hopf algebra.

# Arguments
`index` must lie in `1:algebra.dimension`; otherwise `BoundsError` reports
the admissible range and index. The input must have coherent dimension,
multiplication and unit data for subsequent algebra operations.

# Returns
An `AlgebraElement` with a fresh `Vector{ComplexF64}` of coordinates and a
`FiniteAlgebra` carrying the input's dimension, multiplication and unit,
named `"Hopf algebra"`.

# Guarantees
1. The vector has `algebra.dimension` entries, one at `index` and zero
   elsewhere. It represents that basis vector of the input algebra. [proved/agent: by construction]
2. Multiplication and unit factors are shared with the input. The input is
   unchanged; no coproduct or antipode is stored in the result. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
No algebra identities are checked and no canonical basis is chosen.
The basis belongs to the input constructor in the
[Algebra API](../specs/algebra.md#definitions).

# Sources
cf:[[Kup91#^kup91-blk-0024]] coordinates in a chosen finite-dimensional basis.

# Examples
```jldoctest
julia> h = cyclic_group_algebra(5);

julia> (basis_element(h, 2)^5).coefficients ≈ basis_element(h, 1).coefficients
true
```

See also [`AlgebraElement`](@ref).

Role: access
"""
basis_element(algebra::FiniteHopfAlgebra,index::Integer)=basis_element(_finite_algebra(algebra),index)

_same_algebra(left::FiniteAlgebra,right::FiniteAlgebra)=left===right ||
    (left.dimension==right.dimension && left.name==right.name &&
     left.multiplication.terms==right.multiplication.terms && left.unit.terms==right.unit.terms)

function Base.:*(left::AlgebraElement,right::AlgebraElement)
    _same_algebra(left.algebra,right.algebra) || throw(ArgumentError("elements belong to different algebras"))
    result=zeros(ComplexF64,left.algebra.dimension)
    for ((i,j,k),value) in left.algebra.multiplication.terms
        result[k]+=value*left.coefficients[i]*right.coefficients[j]
    end
    AlgebraElement(left.algebra,result)
end
Base.:+(left::AlgebraElement,right::AlgebraElement)=begin
    _same_algebra(left.algebra,right.algebra) || throw(ArgumentError("elements belong to different algebras"))
    AlgebraElement(left.algebra,left.coefficients+right.coefficients)
end
Base.:-(left::AlgebraElement,right::AlgebraElement)=begin
    _same_algebra(left.algebra,right.algebra) || throw(ArgumentError("elements belong to different algebras"))
    AlgebraElement(left.algebra,left.coefficients-right.coefficients)
end
Base.:*(scalar::Number,element::AlgebraElement)=AlgebraElement(element.algebra,scalar.*element.coefficients)
Base.:*(element::AlgebraElement,scalar::Number)=scalar*element
Base.:^(element::AlgebraElement,power::Integer)=begin
    power>=0 || throw(ArgumentError("negative algebra-element powers are unsupported"))
    unit=zeros(ComplexF64,element.algebra.dimension)
    for ((index,),value) in element.algebra.unit.terms; unit[index]=value; end
    result=AlgebraElement(element.algebra,unit); base=element; exponent=power
    while exponent>0
        isodd(exponent) && (result=result*base)
        exponent>>=1; exponent>0 && (base=base*base)
    end
    result
end
Base.isapprox(left::AlgebraElement,right::AlgebraElement;kwargs...)=
    _same_algebra(left.algebra,right.algebra) && isapprox(left.coefficients,right.coefficients;kwargs...)

function _vertex_constants(multiplication,coproduct,antipode)
    positive=Dict{Tuple,ComplexF64}(); negative=Dict{Tuple,ComplexF64}()
    antipode_by_input=Dict{Int,Vector{Tuple{Int,ComplexF64}}}()
    for ((output,input),value) in antipode
        push!(get!(antipode_by_input,input,Tuple{Int,ComplexF64}[]),(output,value))
    end
    for ((delta_input,left,right),delta_value) in coproduct
        for ((mu_input,middle,output),mu_value) in multiplication
            middle==left || continue
            key=(output,right,mu_input,delta_input)
            positive[key]=get(positive,key,0)+mu_value*delta_value
        end
        for (antipode_output,antipode_value) in get(antipode_by_input,left,Tuple{Int,ComplexF64}[])
            for ((mu_input,middle,output),mu_value) in multiplication
                middle==antipode_output || continue
                key=(output,right,mu_input,delta_input)
                negative[key]=get(negative,key,0)+mu_value*antipode_value*delta_value
            end
        end
    end
    positive,negative
end

function _derived_integrals(multiplication,coproduct,antipode,dimension)
    tensors=Dict(
        "Mu"=>_sparse_data(multiplication),
        "De"=>_sparse_data(coproduct),
        "An"=>_sparse_data(antipode),
    )
    projection=dense_tensor(evaluate_tensor_expression(
        "Mu{1,2,3}An{3,4}De{4,5,6}An{5,1}",tensors,[2,6];dimension=dimension),
        dimension)
    decomposition=svd(projection)
    right_integral=decomposition.U[:,1]
    right_cointegral=decomposition.S[1].*decomposition.V[:,1]
    right_cointegral ./= transpose(right_integral)*right_cointegral
    antipode_matrix=zeros(ComplexF64,dimension,dimension)
    for ((output,input),value) in antipode
        antipode_matrix[output,input]=value
    end
    left_cointegral=inv(antipode_matrix)*right_cointegral
    left_integral=transpose(antipode_matrix)*right_integral
    sparse(vector)=Dict{Tuple,ComplexF64}((index,)=>value
        for (index,value) in enumerate(vector) if abs(value)>1e-10)
    sparse(right_integral),sparse(right_cointegral),sparse(left_integral),sparse(left_cointegral)
end

function _make_hopf_algebra(dimension,multiplication,unit,coproduct,counit,antipode;
                            integrals=nothing)
    positive,negative=_vertex_constants(multiplication,coproduct,antipode)
    right_integral,right_cointegral,left_integral,left_cointegral =
        integrals === nothing ? _derived_integrals(multiplication,coproduct,antipode,dimension) : integrals
    FiniteHopfAlgebra(dimension,_sparse_data(multiplication),_sparse_data(unit),
        _sparse_data(coproduct),_sparse_data(counit),_sparse_data(positive),
        _sparse_data(negative),_sparse_data(antipode),_sparse_data(right_integral),
        _sparse_data(right_cointegral),_sparse_data(left_integral),
        _sparse_data(left_cointegral))
end

"""
    group_algebra(table::AbstractMatrix{<:Integer}) -> FiniteHopfAlgebra

Return exact structure tensors for the Hopf algebra of a finite group supplied by its multiplication table.

# Arguments
`table` is a nonempty square group multiplication table with entries
convertible to `Int`. Shape, closure, a two-sided identity, associativity and
unique two-sided inverses are checked. Violations raise `ArgumentError`
naming the failed group-table requirement; conversion and allocation errors
propagate. The table's index convention belongs to the
[Algebra API](../specs/algebra.md#definitions).

# Returns
A `FiniteHopfAlgebra` of dimension `size(table,1)` with `BigInt` structure
coefficients. Fields have the API roles described by [`FiniteHopfAlgebra`](@ref).

# Guarantees
1. Multiplication, unit, coproduct, counit and antipode are the group-algebra
   structural maps for the validated table in its supplied group basis.
   [proved/agent: by construction]
2. The crossing tensors and distinguished tensors use the same basis and
   agree with the library's cyclic specialization when its table is supplied.
   [proved/agent: by construction]
3. The input table is not modified. Result tensors are newly allocated.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The `integral` and `cointegral` field names do not certify agreement with a
literature naming convention. Exact structural coefficients do not make every
downstream operation exact or certify a topological invariant.

# Sources
cf:[[Kup91#^kup91-blk-0138]] the group-basis specialization of Hopf tensors; storage and distinguished-field names follow the repository API.

# Examples
```jldoctest
julia> h = group_algebra([2 1; 1 2]);

julia> h.unit.terms == Dict{Tuple,ComplexF64}((2,)=>1)
true
```

See also [`cyclic_group_algebra`](@ref), [`verify_hopf`](@ref).

Role: construct
"""
function group_algebra(table::AbstractMatrix{<:Integer})
    rows,columns=size(table)
    rows>0 && rows==columns || throw(ArgumentError("group table must be nonempty and square"))
    order=rows
    values=Matrix{Int}(table)
    all(value -> 1<=value<=order,values) ||
        throw(ArgumentError("group table entries must lie in 1:$order"))
    identity=findfirst(candidate ->
        all(i -> values[candidate,i]==i && values[i,candidate]==i,1:order),1:order)
    identity===nothing && throw(ArgumentError("group table has no two-sided identity"))
    for i in 1:order,j in 1:order,k in 1:order
        values[values[i,j],k]==values[i,values[j,k]] ||
            throw(ArgumentError("group table is not associative at ($i,$j,$k)"))
    end
    inverses=Vector{Int}(undef,order)
    for i in 1:order
        candidates=[j for j in 1:order if values[i,j]==identity && values[j,i]==identity]
        length(candidates)==1 || throw(ArgumentError("group element $i has no unique two-sided inverse"))
        inverses[i]=only(candidates)
    end
    multiplication=Dict{Tuple,BigInt}()
    coproduct=Dict{Tuple,BigInt}(); antipode=Dict{Tuple,BigInt}()
    for i in 1:order,j in 1:order
        multiplication[(i,j,values[i,j])]=1
    end
    for i in 1:order
        coproduct[(i,i,i)]=1
        antipode[(inverses[i],i)]=1
    end
    unit=Dict{Tuple,BigInt}((identity,)=>1)
    counit=Dict{Tuple,BigInt}((i,)=>1 for i in 1:order)
    positive=Dict{Tuple,BigInt}((values[a,b],b,a,b)=>1 for a in 1:order,b in 1:order)
    negative=Dict{Tuple,BigInt}((values[a,inverses[b]],b,a,b)=>1 for a in 1:order,b in 1:order)
    integral=Dict{Tuple,BigInt}((identity,)=>1)
    cointegral=Dict{Tuple,BigInt}((i,)=>1 for i in 1:order)
    FiniteHopfAlgebra{GroupAlgebraKind}(order,_sparse_data(multiplication),_sparse_data(unit),
        _sparse_data(coproduct),_sparse_data(counit),_sparse_data(positive),
        _sparse_data(negative),_sparse_data(antipode),_sparse_data(integral),
        _sparse_data(cointegral),_sparse_data(integral),_sparse_data(cointegral))
end

"""
    group_algebra(elements::AbstractVector, multiply) -> FiniteHopfAlgebra

Return exact group Hopf tensors for a finite group supplied by ordered elements and a product.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

- `elements` must be nonempty with pairwise distinct labels under `isequal`.
  Violations raise `ArgumentError` naming nonemptiness or uniqueness.
- `multiply` must accept two labels and return a label equal under `isequal`
  to exactly one supplied element; failure raises `ArgumentError` reporting
  lack of closure. Callback exceptions propagate. The resulting group-table
  requirements are checked as in [`group_algebra(::AbstractMatrix)`](@ref group_algebra).
- Stable equality and a deterministic, nonmutating product callback are the
  caller's responsibility. A stateful callback can invalidate the intended
  group interpretation and can modify data passed to it.

# Returns
A `FiniteHopfAlgebra` with `BigInt` tensor coefficients, of dimension
`length(elements)`, in the supplied element order.

# Guarantees
1. For a stable product on the supplied labels, the result is the group Hopf
   algebra in that ordered basis, with the structural interpretation of the
   table overload. [proved/agent: by construction]
2. The function itself does not modify the element vector or its labels;
   side effects of `multiply` are not controlled. Tensor containers are new.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The returned record does not retain the label vector or decide isomorphism
between different group presentations. Distinguished-field naming and
downstream exactness have the limitations of the table overload.

# Sources
cf:[[Kup91#^kup91-blk-0138]] group-algebra tensors in a group basis; the ordered labels choose the repository coordinate basis.

# Examples
```jldoctest
julia> h = group_algebra(collect(0:4), (a,b)->mod(a+b,5));

julia> (h.dimension, all(values(verify_hopf(h))))
(5, true)
```

See also [`cyclic_group_algebra`](@ref).

Role: construct
"""
function group_algebra(elements::AbstractVector,multiply)
    isempty(elements) && throw(ArgumentError("group elements must be nonempty"))
    labels=collect(elements)
    for i in eachindex(labels),j in firstindex(labels):i-1
        isequal(labels[i],labels[j]) && throw(ArgumentError("group elements must be unique"))
    end
    table=Matrix{Int}(undef,length(labels),length(labels))
    for i in eachindex(labels),j in eachindex(labels)
        product=multiply(labels[i],labels[j])
        matches=findall(label -> isequal(label,product),labels)
        length(matches)==1 || throw(ArgumentError("group operation is not closed on the supplied elements"))
        table[i,j]=only(matches)
    end
    group_algebra(table)
end

"""
    cyclic_group_algebra(order::Integer) -> FiniteHopfAlgebra

Return exact structure tensors for the group Hopf algebra of a cyclic group.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

`order` must be positive; otherwise `ArgumentError` states that the group
order must be positive. The dimension must fit the machine's array indices
and available storage. Conversion or allocation failures propagate.

# Returns
`FiniteHopfAlgebra` of dimension `order`, with `BigInt` coefficients in its
sparse tensors. Its fields have the roles described by
[`FiniteHopfAlgebra`](@ref).

# Guarantees
1. The multiplication, unit, coproduct, counit and antipode are those of
   the cyclic group's group algebra in the ordered group basis. [proved/agent: by construction]
2. The result has the same tensor values as
   `group_algebra(collect(0:order-1), (a,b) -> mod(a+b,order))`. [proved/agent: by construction]
3. Each call returns independently allocated tensor data. [proved/agent: by construction]
4. The supplied arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The integer tensor coefficients do not make all downstream operations exact:
`basis_element`, duals and doubles use complex floating-point coordinates.
No translation to another basis or integral normalization is certified, and
no manifold invariant is computed by this constructor.

# Sources
cf:[[Kup91#^kup91-blk-0138]] group-algebra weights in a group basis; the cyclic specialization and integer coefficient representation are reconstructed.

# Examples
```jldoctest
julia> h = cyclic_group_algebra(5);

julia> all(values(verify_hopf(h)))
true

julia> (basis_element(h, 2)^5).coefficients ≈ basis_element(h, 1).coefficients
true
```

See also [`group_algebra`](@ref), [`basis_element`](@ref).

Role: construct
"""
function cyclic_group_algebra(order::Integer)
    order>0 || throw(ArgumentError("group order must be positive"))
    table=[mod(i+j-2,order)+1 for i in 1:order,j in 1:order]
    group_algebra(table)
end

_kp_index(x,y,z)=1+x+2y+4z
_kp_bits(index)=((index-1)&1,(index-1)>>1&1,(index-1)>>2&1)

"""Multiply basis elements `x^a*y^b*z^c` of the Kac–Paljutkin algebra H8."""
function _kp_basis_product(left::Int,right::Int,::Type{T}=ComplexF64) where T
    a,b,c=_kp_bits(left); d,e,f=_kp_bits(right)
    # z*x=y*z and z*y=x*z, so a left z swaps the x/y powers on its right.
    p=xor(a,c==1 ? e : d); q=xor(b,c==1 ? d : e)
    if c+f<2
        return Dict(_kp_index(p,q,c+f)=>one(T))
    end
    # 2z^2=1+x+y-xy.
    Dict(
        _kp_index(p,q,0)=>one(T)/2,
        _kp_index(xor(p,1),q,0)=>one(T)/2,
        _kp_index(p,xor(q,1),0)=>one(T)/2,
        _kp_index(xor(p,1),xor(q,1),0)=>-one(T)/2,
    )
end

function _kp_tensor_product(left::Dict{Tuple{Int,Int},T},
                            right::Dict{Tuple{Int,Int},T}) where T
    result=Dict{Tuple{Int,Int},T}()
    for ((l1,r1),v1) in left,((l2,r2),v2) in right
        for (lo,lv) in _kp_basis_product(l1,l2,T),(ro,rv) in _kp_basis_product(r1,r2,T)
            key=(lo,ro); result[key]=get(result,key,0)+v1*v2*lv*rv
        end
    end
    approximate=T<:Union{AbstractFloat,Complex{<:AbstractFloat}}
    filter!(pair -> approximate ? abs(last(pair))>1e-12 : !iszero(last(pair)),result)
    result
end

"""Local structure constants of H8 in `[1,x,y,xy,z,xz,yz,xyz]` order.

cf:src/finite_algebra.jl @ b224700, kac_paljutkin_algebra — same local definition, now parameterized independently of numerical integrals.
"""
function _kp_structure_constants(::Type{T}=ComplexF64) where T
    dimension=8
    multiplication=Dict{Tuple,T}()
    for left in 1:dimension,right in 1:dimension,(output,value) in _kp_basis_product(left,right,T)
        multiplication[(left,right,output)]=value
    end

    delta_x=Dict((_kp_index(1,0,0),_kp_index(1,0,0))=>one(T))
    delta_y=Dict((_kp_index(0,1,0),_kp_index(0,1,0))=>one(T))
    delta_z=Dict(
        (_kp_index(0,0,1),_kp_index(0,0,1))=>one(T)/2,
        (_kp_index(1,0,1),_kp_index(0,0,1))=>one(T)/2,
        (_kp_index(0,0,1),_kp_index(0,1,1))=>one(T)/2,
        (_kp_index(1,0,1),_kp_index(0,1,1))=>-one(T)/2,
    )
    coproduct=Dict{Tuple,T}()
    for source in 1:dimension
        a,b,c=_kp_bits(source)
        terms=Dict((1,1)=>one(T))
        a==1 && (terms=_kp_tensor_product(terms,delta_x))
        b==1 && (terms=_kp_tensor_product(terms,delta_y))
        c==1 && (terms=_kp_tensor_product(terms,delta_z))
        for ((left,right),value) in terms
            coproduct[(source,left,right)]=value
        end
    end

    antipode=Dict{Tuple,T}()
    for source in 1:dimension
        a,b,c=_kp_bits(source)
        terms=Dict(1=>one(T))
        for factor in vcat(c==1 ? [5] : Int[],b==1 ? [3] : Int[],a==1 ? [2] : Int[])
            next=Dict{Int,T}()
            for (left,left_value) in terms,(output,value) in _kp_basis_product(left,factor,T)
                next[output]=get(next,output,0)+left_value*value
            end
            terms=next
        end
        for (output,value) in terms
            antipode[(output,source)]=value
        end
    end
    unit=Dict{Tuple,T}((1,)=>one(T))
    counit=Dict{Tuple,T}((index,)=>one(T) for index in 1:dimension)
    (;dimension,multiplication,unit,coproduct,counit,antipode)
end

# cf:_kp_structure_constants — same local definition, retaining the existing approximate integral construction.
"""
    kac_paljutkin_algebra() -> FiniteHopfAlgebra

Return the library's eight-dimensional Kac–Paljutkin tensor model with numerical distinguished tensors.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

No arguments are required. Numerical linear-algebra and allocation errors
can propagate when constructing the model.

# Returns
`FiniteHopfAlgebra{Nothing}` of dimension eight with `ComplexF64` coefficients
in the fields described by [`FiniteHopfAlgebra`](@ref).

# Guarantees
1. The result is the repository's fixed eight-dimensional algebra and
   coalgebra model, with crossing, antipode and distinguished tensors in
   the same coordinate space. [proved/agent: by construction]
2. Each call returns newly allocated tensor data; there are no input
   arguments to modify. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The distinguished tensors are numerical, not exact algebraic certificates.
Their normalization and basis identification with another literature
presentation are not certified. Construction does not return a proof of
Hopf identities or topological invariance.

# Sources
cf:lib/topology.jl/src/finite_algebra.jl, _kp_structure_constants — the fixed local eight-dimensional model; no literature basis identification is asserted.

# Examples
```jldoctest
julia> h = kac_paljutkin_algebra(); x = basis_element(h,2); y = basis_element(h,3); z = basis_element(h,5);

julia> (z*x ≈ y*z, z*y ≈ x*z)
(true, true)
```

See also [`verify_hopf`](@ref), [`is_involutory`](@ref).

Role: construct
"""
function kac_paljutkin_algebra()
    h=_kp_structure_constants()
    _make_hopf_algebra(h.dimension,h.multiplication,h.unit,h.coproduct,h.counit,h.antipode)
end

"""
    dual_algebra(algebra::FiniteHopfAlgebra) -> FiniteHopfAlgebra

Return tensor data for the linear dual of a finite-dimensional Hopf algebra.

# Arguments
`algebra` must supply coherent finite-dimensional Hopf structure tensors
in a fixed basis. Multiplication and coproduct keys must have three entries
and antipode keys two entries; their coefficients must be convertible to
`ComplexF64`. Invalid keys and conversion failures propagate their Julia
exceptions. No Hopf-axiom or dimension validation is performed.

# Returns
`FiniteHopfAlgebra{Nothing}` of the same dimension. Its multiplication and
unit correspond to the input coproduct and counit; its coproduct and counit
correspond to the input multiplication and unit; its antipode is the dual
map. The distinguished left/right tensors exchange the integral and
cointegral field roles. See the [Algebra API](../specs/algebra.md#definitions).

# Guarantees
1. For valid input Hopf data, the structural maps are those of its linear
   dual in the dual basis, subject to conversion of multiplication,
   coproduct and antipode coefficients to `ComplexF64`. [proved/agent: by construction]
2. The input is unchanged. The four distinguished integral/cointegral
   factors are shared with the result; the other tensor dictionaries are new. [proved/agent: by construction]
3. The returned dimension equals the input dimension. No involutivity
   assumption is checked or imposed. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Floating-point conversion need not preserve exact coefficients. Applying
duality twice need not recover the input's type tag or bitwise tensor data.
The stored integral-field names do not certify a normalization or resolve
conventions absent from the algebra specification.

# Sources
cf:[[Kup91#^kup91-sec-0038]] algebra structures on a finite-dimensional Hopf algebra and its dual; dual antipode and distinguished tensors use the repository field convention.

# Examples
```jldoctest
julia> d = dual_algebra(uqsl2_borel_small(1//3; style=:L));

julia> (d.dimension, all(values(verify_hopf(d))))
(9, true)
```

See also [`hopf_pairing`](@ref), [`verify_hopf`](@ref).

Role: convert
"""
function dual_algebra(algebra::FiniteHopfAlgebra)
    multiplication=Dict{Tuple,ComplexF64}((left,right,output)=>value
        for ((output,left,right),value) in algebra.coproduct.terms)
    coproduct=Dict{Tuple,ComplexF64}((output,left,right)=>value
        for ((left,right,output),value) in algebra.multiplication.terms)
    unit=copy(algebra.counit.terms); counit=copy(algebra.unit.terms)
    antipode=Dict{Tuple,ComplexF64}((input,output)=>value
        for ((output,input),value) in algebra.antipode.terms)
    positive,negative=_vertex_constants(multiplication,coproduct,antipode)
    FiniteHopfAlgebra(algebra.dimension,_sparse_data(multiplication),_sparse_data(unit),
        _sparse_data(coproduct),_sparse_data(counit),_sparse_data(positive),
        _sparse_data(negative),_sparse_data(antipode),
        algebra.right_cointegral,algebra.right_integral,
        algebra.left_cointegral,algebra.left_integral)
end

"""
    hopf_pairing(dual::AlgebraElement, primal::AlgebraElement) -> ComplexF64

Return the bilinear coordinate pairing of a dual element with a primal element.

# Arguments
The two stored algebra dimensions must agree; otherwise `DimensionMismatch`
is thrown without a custom message. For the evaluation interpretation,
`dual` and `primal` must be expressed in corresponding dual bases and have
coordinate vectors of the declared dimension. That correspondence and those
lengths are not checked; incompatible broadcast lengths can throw
`DimensionMismatch`.

# Returns
A `ComplexF64` scalar. This operation belongs to the
[Algebra API](../specs/algebra.md#definitions).

# Guarantees
1. For corresponding bases and correctly sized vectors, the result is the
   value of the supplied linear functional on the supplied vector, subject
   to complex floating-point arithmetic. [proved/agent: by construction]
2. The pairing is bilinear: the dual coefficients are not conjugated. [proved/agent: by construction]
3. Both inputs are unchanged. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Equal dimensions do not establish that the two stored algebras are dual.
This is not a Hermitian inner product or an exact pairing for arbitrary
exact mathematical coefficients.

# Sources
cf:[[Kup91#^kup91-blk-0023]] evaluation of a dual vector on a vector by contraction in corresponding bases.

# Examples
```jldoctest
julia> h = uqsl2_borel_small(1//3; style=:L); d = dual_algebra(h);

julia> K = basis_element(h, 2); E = basis_element(h, 4); pE = basis_element(d, 4);

julia> hopf_pairing(pE, E*K) ≈ 0
true
```

See also [`dual_algebra`](@ref), [`basis_element`](@ref).

Role: invariant
"""
hopf_pairing(dual::AlgebraElement,primal::AlgebraElement)=begin
    dual.algebra.dimension==primal.algebra.dimension || throw(DimensionMismatch())
    sum(dual.coefficients.*primal.coefficients)
end

_combined_index(first,second,dimension)=first+dimension*(second-1)

function _double_algebra(algebra::FiniteHopfAlgebra,kind::Symbol)
    dimension=algebra.dimension
    tensors=hopf_tensors(algebra)
    if kind==:heisenberg
        factor=evaluate_tensor_expression(
            "De{5,1,7}De{2,9,8}Mu{7,9,3}Mu{8,4,6}",tensors,1:6;
            dimension=dimension)
        name="Heisenberg double"
    elseif kind==:drinfeld
        tensors["Si"]=_matrix_power_factor(algebra.antipode,-1,dimension)
        factor=evaluate_tensor_expression(
            "Mu{7,8,3}Mu{9,10,8}De{2,10,11}De{11,12,13}" *
            "Si{7,13}Mu{12,4,6}De{5,1,9}",tensors,1:6;dimension=dimension)
        name="Drinfeld double"
    else
        throw(ArgumentError("unknown double kind"))
    end
    multiplication=Dict{Tuple,ComplexF64}()
    for ((f1,a1,f2,a2,fo,ao),value) in factor.terms
        key=(_combined_index(f1,a1,dimension),_combined_index(f2,a2,dimension),
             _combined_index(fo,ao,dimension))
        multiplication[key]=get(multiplication,key,0)+value
    end
    unit=Dict{Tuple,ComplexF64}()
    for ((functional,),epsilon_value) in algebra.counit.terms
        for ((element,),unit_value) in algebra.unit.terms
            unit[(_combined_index(functional,element,dimension),)]=epsilon_value*unit_value
        end
    end
    FiniteAlgebra(dimension^2,_sparse_data(multiplication),_sparse_data(unit),name)
end

"""
    heisenberg_double(algebra::FiniteHopfAlgebra) -> FiniteAlgebra

Return the repository's Heisenberg-double product data on the dual/primal tensor-product space.

# Arguments
`algebra` must contain coherent dimension and structure-tensor data with
valid coordinate indices and coefficients convertible to complex floating
point. Malformed tensors propagate indexing, conversion or tensor-expression
errors; the Hopf axioms are not checked here.

# Returns
`FiniteAlgebra` of dimension `algebra.dimension^2`, named
`"Heisenberg double"`, with complex floating-point multiplication data and
the tensor-product unit. The algebra operations belong to
[Algebra](../specs/algebra.md#definitions).

# Guarantees
1. The underlying coordinate space is the tensor product of the input's
   dual and primal vector spaces. Multiplication and unit are the stored
   double-model tensors on this space. [proved/agent: by construction]
2. The returned dimension is the square of the input dimension, when that
   value is representable. The input is unchanged and result tensors are new. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The current specification does not establish the identification of this
product with a literature Heisenberg double under a fixed action and
opposite-algebra convention. Construction does not validate associativity.
No Hopf structure or universal R-matrix is returned, and exact input
coefficients do not imply exact output arithmetic.

# Sources
cf:[[Kup91#^kup91-sec-0018]] tensor products and contractions in finite-dimensional coordinates; this supports the tensor-space interpretation, not the identification of the repository product with a literature double.

# Examples
```jldoctest
julia> d = heisenberg_double(uqsl2_borel_small(1//3; style=:L));

julia> (d.dimension, !isempty(d.multiplication.terms))
(81, true)
```

See also [`drinfeld_double`](@ref), [`basis_element`](@ref).

Role: convert
"""
heisenberg_double(algebra::FiniteHopfAlgebra)=_double_algebra(algebra,:heisenberg)

"""
    drinfeld_double(algebra::FiniteHopfAlgebra) -> FiniteAlgebra

Return the repository's Drinfeld-double product data on the dual/primal tensor-product space.

# Arguments
`algebra` must contain coherent finite-dimensional Hopf data with valid
coordinate indices and coefficients convertible to complex floating point.
Its antipode matrix must be invertible; singular inversion propagates
`LinearAlgebra.SingularException`. Invalid tensor data can propagate
indexing, conversion or tensor-expression errors. The Hopf axioms are not
checked here.

# Returns
`FiniteAlgebra` of dimension `algebra.dimension^2`, named
`"Drinfeld double"`, with complex floating-point multiplication data and
the tensor-product unit. This belongs to the
[Algebra API](../specs/algebra.md#definitions).

# Guarantees
1. The underlying coordinate space is the tensor product of the input's
   dual and primal vector spaces. Multiplication and unit are the stored
   double-model tensors on this space. [proved/agent: by construction]
2. The returned dimension is the square of the input dimension, when that
   value is representable. The input is unchanged and result tensors are new. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The current specification does not establish which opposite/copposite and
basis conventions identify this product with a literature Drinfeld double.
Construction does not validate associativity or return coproduct, antipode
or universal R-matrix data for the double. Floating-point arithmetic and
tensor zero filtering do not certify exact identities.

# Sources
cf:[[Kup91#^kup91-sec-0018]] tensor products and contractions in finite-dimensional coordinates; the double-product identification requires additional convention evidence.

# Examples
```jldoctest
julia> d = drinfeld_double(uqsl2_borel_small(1//3; style=:L));

julia> (d.dimension, !isempty(d.multiplication.terms))
(81, true)
```

See also [`heisenberg_double`](@ref), [`basis_element`](@ref).

Role: convert
"""
drinfeld_double(algebra::FiniteHopfAlgebra)=_double_algebra(algebra,:drinfeld)
