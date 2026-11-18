# cf:threads/PJ261117/code/local_hopf_adapter.jl @ 23ff0cbff — local maps and shared MST wiring.
# cf:src/hopf_invariant.jl @ b326299, _hopf_invariant_factors — edge port order.
using AkaghefAlgebra.FiniteFactorEvaluation: TableKernel, LocalFactor, Network

function local_hopf(data,::Type{T}) where {T<:Union{BigInt,Rational{BigInt}}}
    d=data.dimension
    (;dimension=d,scalar=one(T),
      mu=TableKernel((d,d,d),data.multiplication,T),
      delta=TableKernel((d,d,d),data.coproduct,T),
      antipode=TableKernel((d,d),data.antipode,T),
      unit=TableKernel((d,),data.unit,T),counit=TableKernel((d,),data.counit,T))
end
"Exact Kac–Paljutkin local maps in the constructor's basis order."
function kp_adapter()
    owner=parentmodule(FiniteHopfAlgebra)
    local_hopf(owner._kp_structure_constants(Rational{BigInt}),Rational{BigInt})
end
"Read the local maps of a finite Hopf algebra with integer coefficients."
function group_adapter(h::FiniteHopfAlgebra)
    data=(;dimension=h.dimension,multiplication=h.multiplication.terms,
        coproduct=h.coproduct.terms,antipode=h.antipode.terms,
        unit=h.unit.terms,counit=h.counit.terms)
    local_hopf(data,BigInt)
end

"Construct the closed MST factor network from exact local Hopf maps."
function mst_network(v::VirtualLink,a)
    _require_closed_ograph_input(v,"mst_network")
    edges=real_edge_table(v); signs=last(real_gauss_code(v))
    weights=v.is_weighted ? _integer_edge_weights(v) : zeros(Int,length(edges))
    d=a.dimension; n=4length(signs); factors=LocalFactor[]
    matrix=zeros(typeof(a.scalar),d,d)
    for ((o,i),c) in a.antipode.terms; matrix[o,i]=c; end
    identity=Matrix{typeof(a.scalar)}(I,d,d)
    involutory=matrix*matrix==identity
    for (e,w) in zip(edges,weights)
        labels=4 .* abs.([e.crossing[1],e.crossing[2]]) .+
            ([e.crossing[1]>0,e.crossing[2]>0]) .- [1,3]
        exponent=involutory ? 0 : Base.Checked.checked_mul(2,w)
        if !involutory && exponent<0 && a.scalar isa BigInt
            throw(ArgumentError("negative powers of a non-involutory antipode require rational coefficients"))
        end
        powered=involutory ? identity : matrix^exponent
        rows=Dict((o,i)=>powered[o,i] for o in 1:d,i in 1:d if !iszero(powered[o,i]))
        push!(factors,LocalFactor(TableKernel((d,d),rows,typeof(a.scalar)),labels))
    end
    for (j,sign) in enumerate(signs)
        o,r,left,input=4(j-1).+(1:4)
        n+=1; t=n
        push!(factors,LocalFactor(a.mu,[left,t,o]))
        if sign<0
            n+=1; u=n
            push!(factors,LocalFactor(a.antipode,[t,u]))
            push!(factors,LocalFactor(a.delta,[input,u,r]))
        else
            push!(factors,LocalFactor(a.delta,[input,t,r]))
        end
    end
    Network([1:d for _ in 1:n],factors;scalar=a.scalar)
end
