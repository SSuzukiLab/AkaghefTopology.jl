struct FiniteAlgebra
    dimension::Int
    multiplication::SparseFactor
    unit::SparseFactor
    name::String
end

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

basis_element(algebra::FiniteAlgebra,index::Integer)=AlgebraElement(algebra,index)
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

function cyclic_group_algebra(order::Integer)
    order>0 || throw(ArgumentError("group order must be positive"))
    multiplication=Dict{Tuple,ComplexF64}()
    coproduct=Dict{Tuple,ComplexF64}(); antipode=Dict{Tuple,ComplexF64}()
    for i in 0:order-1,j in 0:order-1
        multiplication[(i+1,j+1,mod(i+j,order)+1)]=1
    end
    for i in 0:order-1
        coproduct[(i+1,i+1,i+1)]=1
        antipode[(mod(-i,order)+1,i+1)]=1
    end
    unit=Dict{Tuple,ComplexF64}((1,)=>1)
    counit=Dict{Tuple,ComplexF64}((i,)=>1 for i in 1:order)
    positive,negative=_vertex_constants(multiplication,coproduct,antipode)
    integral=Dict{Tuple,ComplexF64}((1,)=>1)
    cointegral=Dict{Tuple,ComplexF64}((i,)=>1 for i in 1:order)
    FiniteHopfAlgebra(order,_sparse_data(multiplication),_sparse_data(unit),
        _sparse_data(coproduct),_sparse_data(counit),_sparse_data(positive),
        _sparse_data(negative),_sparse_data(antipode),_sparse_data(integral),
        _sparse_data(cointegral),_sparse_data(integral),_sparse_data(cointegral))
end

_kp_index(x,y,z)=1+x+2y+4z
_kp_bits(index)=((index-1)&1,(index-1)>>1&1,(index-1)>>2&1)

"""Multiply basis elements `x^a*y^b*z^c` of the Kac–Paljutkin algebra H8."""
function _kp_basis_product(left::Int,right::Int)
    a,b,c=_kp_bits(left); d,e,f=_kp_bits(right)
    # z*x=y*z and z*y=x*z, so a left z swaps the x/y powers on its right.
    p=xor(a,c==1 ? e : d); q=xor(b,c==1 ? d : e)
    if c+f<2
        return Dict(_kp_index(p,q,c+f)=>1.0+0im)
    end
    # 2z^2=1+x+y-xy.
    Dict(
        _kp_index(p,q,0)=>0.5+0im,
        _kp_index(xor(p,1),q,0)=>0.5+0im,
        _kp_index(p,xor(q,1),0)=>0.5+0im,
        _kp_index(xor(p,1),xor(q,1),0)=>-0.5+0im,
    )
end

function _kp_tensor_product(left::Dict{Tuple{Int,Int},ComplexF64},
                            right::Dict{Tuple{Int,Int},ComplexF64})
    result=Dict{Tuple{Int,Int},ComplexF64}()
    for ((l1,r1),v1) in left,((l2,r2),v2) in right
        for (lo,lv) in _kp_basis_product(l1,l2),(ro,rv) in _kp_basis_product(r1,r2)
            key=(lo,ro); result[key]=get(result,key,0)+v1*v2*lv*rv
        end
    end
    filter!(pair -> abs(last(pair))>1e-12,result)
    result
end

"""Eight-dimensional Kac–Paljutkin Hopf algebra in MATLAB's `[1,x,y,xy,z,xz,yz,xyz]` basis."""
function kac_paljutkin_algebra()
    dimension=8
    multiplication=Dict{Tuple,ComplexF64}()
    for left in 1:dimension,right in 1:dimension,(output,value) in _kp_basis_product(left,right)
        multiplication[(left,right,output)]=value
    end

    delta_x=Dict((_kp_index(1,0,0),_kp_index(1,0,0))=>1.0+0im)
    delta_y=Dict((_kp_index(0,1,0),_kp_index(0,1,0))=>1.0+0im)
    delta_z=Dict(
        (_kp_index(0,0,1),_kp_index(0,0,1))=>0.5+0im,
        (_kp_index(1,0,1),_kp_index(0,0,1))=>0.5+0im,
        (_kp_index(0,0,1),_kp_index(0,1,1))=>0.5+0im,
        (_kp_index(1,0,1),_kp_index(0,1,1))=>-0.5+0im,
    )
    coproduct=Dict{Tuple,ComplexF64}()
    for source in 1:dimension
        a,b,c=_kp_bits(source)
        terms=Dict((1,1)=>1.0+0im)
        a==1 && (terms=_kp_tensor_product(terms,delta_x))
        b==1 && (terms=_kp_tensor_product(terms,delta_y))
        c==1 && (terms=_kp_tensor_product(terms,delta_z))
        for ((left,right),value) in terms
            coproduct[(source,left,right)]=value
        end
    end

    antipode=Dict{Tuple,ComplexF64}()
    for source in 1:dimension
        a,b,c=_kp_bits(source)
        terms=Dict(1=>1.0+0im)
        for factor in vcat(c==1 ? [5] : Int[],b==1 ? [3] : Int[],a==1 ? [2] : Int[])
            next=Dict{Int,ComplexF64}()
            for (left,left_value) in terms,(output,value) in _kp_basis_product(left,factor)
                next[output]=get(next,output,0)+left_value*value
            end
            terms=next
        end
        for (output,value) in terms
            antipode[(output,source)]=value
        end
    end
    unit=Dict{Tuple,ComplexF64}((1,)=>1)
    counit=Dict{Tuple,ComplexF64}((index,)=>1 for index in 1:dimension)
    _make_hopf_algebra(dimension,multiplication,unit,coproduct,counit,antipode)
end

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

heisenberg_double(algebra::FiniteHopfAlgebra)=_double_algebra(algebra,:heisenberg)
drinfeld_double(algebra::FiniteHopfAlgebra)=_double_algebra(algebra,:drinfeld)
