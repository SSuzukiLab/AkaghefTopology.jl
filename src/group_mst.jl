# cf:lib/topology.jl/src/hopf_invariant.jl @ af90b62, _hopf_invariant_factors — identical four-port/edge ordering; working-tree split: src/hopf_link.jl.
# cf:lib/topology.jl/src/finite_algebra.jl group_algebra — group-basis W(a,b)=(ab,b), inverse vertex (ab⁻¹,b).
# cf:Graphs.connected_components (version in Manifest.toml) — identify equal indices and split independent sums.
const _GroupAlgebraKind=parentmodule(FiniteHopfAlgebra).GroupAlgebraKind

"Finite sum of group-valued indices, subject to ternary product constraints."
struct GroupConstraintExpression
    multiplication::Matrix{Int}
    inverse::Vector{Int}
    constraints::Vector{NTuple{4,Int}} # (a,b,c,sign): c=a*b^sign
    domains::Vector{BitVector}
    consistent::Bool
end

function _group_data(H)
    d=H.dimension
    length(H.multiplication.terms)==d^2 || return nothing
    table=zeros(Int,d,d)
    for (key,value) in H.multiplication.terms
        length(key)==3 && value==1 || return nothing
        a,b,c=key
        all(k->1<=k<=d,key) && table[a,b]==0 || return nothing
        table[a,b]=c
    end
    inv=zeros(Int,d)
    length(H.antipode.terms)==d || return nothing
    for (key,value) in H.antipode.terms
        length(key)==2 && value==1 || return nothing
        a,b=key
        1<=a<=d && 1<=b<=d && inv[b]==0 || return nothing
        inv[b]=a
    end
    all(k->inv[k]>0 && inv[inv[k]]==k,1:d) || return nothing
    length(H.coproduct.terms)==d && all(get(H.coproduct.terms,(i,i,i),0)==1 for i in 1:d) || return nothing
    for (tensor,sign) in ((H.positive_vertex,1),(H.negative_vertex,-1))
        length(tensor.terms)==d^2 || return nothing
        for a in 1:d,b in 1:d
            get(tensor.terms,(table[a,sign==1 ? b : inv[b]],b,a,b),0)==1 || return nothing
        end
    end
    (table,inv)
end

"Existing sparse network, with its integer coefficients retained for comparison."
function _group_sparse_reference(v,H)
    factors=SparseFactor[]
    for f in _hopf_invariant_factors(v,H)
        terms=Dict{Tuple,BigInt}()
        for (key,value) in f.terms
            iszero(imag(value)) && isinteger(real(value)) || error("noninteger group reference coefficient")
            terms[key]=BigInt(real(value))
        end
        push!(factors,SparseFactor(copy(f.variables),terms))
    end
    value=_contract(factors)
    iszero(imag(value)) && isinteger(real(value)) || error("noninteger group reference value")
    BigInt(real(value))
end

function invariant_expression(v::VirtualLink,H::FiniteHopfAlgebra{_GroupAlgebraKind})
    _require_closed_ograph_input(v,"group invariant_expression")
    data=_group_data(H)
    data===nothing && return MSTTensorExpression(_hopf_invariant_factors(v,H))
    v.is_weighted && _integer_edge_weights(v)
    signs=last(real_gauss_code(v))
    isempty(signs) && throw(ArgumentError("group invariant requires a real vertex"))
    equalities=Graphs.SimpleGraph(4length(signs))
    for edge in real_edge_table(v)
        labels=4 .* abs.([edge.crossing[1],edge.crossing[2]]) .+
            ([edge.crossing[1]>0,edge.crossing[2]>0]) .- [1,3]
        labels[1]!=labels[2] && Graphs.add_edge!(equalities,labels...)
    end
    # Delta(b)=b⊗b and S^(2w)=id, including negative integer weights.
    for j in eachindex(signs); Graphs.add_edge!(equalities,4j-2,4j); end
    components=Graphs.connected_components(equalities)
    labels=zeros(Int,4length(signs))
    for (i,component) in enumerate(components); labels[component].=i; end
    constraints=[(labels[4j-1],labels[4j],labels[4j-3],Int(signs[j])) for j in eachindex(signs)]
    GroupConstraintExpression(data...,constraints,[trues(H.dimension) for _ in components],true)
end

mutable struct GroupSumStats
    nodes::Int
    branches::Int
    support_checks::Int
    max_depth::Int
    components::Int
end

function _group_propagate!(domains,constraints,table,inv,stats)
    all(any,domains) || return false
    d=size(table,1); changed=true
    while changed
        changed=false
        for (ia,ib,ic,sign) in constraints
            sa,sb,sc=falses(d),falses(d),falses(d)
            for a in 1:d
                domains[ia][a] || continue
                for b in 1:d
                    domains[ib][b] || continue
                    stats.support_checks+=1
                    ia==ib && a!=b && continue
                    c=table[a,sign==1 ? b : inv[b]]
                    domains[ic][c] || continue
                    ia==ic && a!=c && continue
                    ib==ic && b!=c && continue
                    sa[a]=true; sb[b]=true; sc[c]=true
                end
            end
            for (i,support) in ((ia,sa),(ib,sb),(ic,sc))
                old=count(domains[i]); domains[i] .&= support
                any(domains[i]) || return false
                changed |= old!=count(domains[i])
            end
        end
    end
    true
end

function _group_tautology(constraint,domains,table,inv)
    ia,ib,ic,sign=constraint
    total=prod((BigInt(count(domains[i])) for i in unique((ia,ib,ic)));init=BigInt(1))
    valid=BigInt(0)
    for a in findall(domains[ia]),b in findall(domains[ib])
        ia==ib && a!=b && continue
        c=table[a,sign==1 ? b : inv[b]]
        domains[ic][c] && (ia!=ic || a==c) && (ib!=ic || b==c) && (valid+=1)
    end
    valid==total
end

function simplify_invariant(x::GroupConstraintExpression;strategy=:auto)
    strategy in (:auto,:qa,:network,:none,:numeric) || throw(ArgumentError("unknown MST strategy"))
    domains=deepcopy(x.domains)
    stats=GroupSumStats(0,0,0,0,0)
    ok=x.consistent && (strategy==:none || _group_propagate!(domains,x.constraints,x.multiplication,x.inverse,stats))
    constraints=strategy==:none || !ok ? copy(x.constraints) :
        filter(c->!_group_tautology(c,domains,x.multiplication,x.inverse),x.constraints)
    GroupConstraintExpression(x.multiplication,x.inverse,constraints,domains,ok)
end

function _group_search(domains,constraints,variables,table,inv,stats,depth)
    stats.nodes+=1; stats.max_depth=max(stats.max_depth,depth)
    _group_propagate!(domains,constraints,table,inv,stats) || return BigInt(0)
    constraints=filter(c->!_group_tautology(c,domains,table,inv),constraints)
    live=filter(i->count(domains[i])>1,variables)
    isempty(live) && return BigInt(1)
    if isempty(constraints)
        return prod((BigInt(count(domains[i])) for i in live);init=BigInt(1))
    end
    degree(i)=count(c->i in c[1:3],constraints)
    chosen=first(sort(live;by=i->(count(domains[i]),-degree(i),i)))
    result=BigInt(0)
    for value in findall(domains[chosen])
        stats.branches+=1
        next=deepcopy(domains); next[chosen].=false; next[chosen][value]=true
        result+=_group_search(next,constraints,variables,table,inv,stats,depth+1)
    end
    result
end

"Count compatible assignments exactly, with no intermediate assignment tables."
function evaluate_invariant(x::GroupConstraintExpression;details=false)
    stats=GroupSumStats(0,0,0,0,0); domains=deepcopy(x.domains)
    ok=x.consistent && _group_propagate!(domains,x.constraints,x.multiplication,x.inverse,stats)
    result=BigInt(0)
    if ok
        n=length(domains); dependencies=Graphs.SimpleGraph(n)
        live=filter(i->count(domains[i])>1,1:n)
        for c in x.constraints
            active=unique(filter(i->count(domains[i])>1,c[1:3]))
            for i in active
                i!=first(active) && Graphs.add_edge!(dependencies,first(active),i)
            end
        end
        components=filter(!isempty,[filter(i->i in live,c) for c in Graphs.connected_components(dependencies)])
        stats.components=length(components); result=BigInt(1)
        for variables in components
            constraints=filter(c->any(i->i in variables,c[1:3]),x.constraints)
            result*=_group_search(deepcopy(domains),constraints,variables,x.multiplication,x.inverse,stats,0)
            iszero(result) && break
        end
    end
    details ? (value=result,variables=length(domains),free_indices=count(d->count(d)>1,domains),
        domain_cells=sum(length,domains;init=0),nodes=stats.nodes,branches=stats.branches,
        support_checks=stats.support_checks,max_depth=stats.max_depth,components=stats.components) : result
end
