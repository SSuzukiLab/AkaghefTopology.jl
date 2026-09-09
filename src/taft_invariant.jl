const _TaftST = AkaghefAlgebra.SymbolicTaft
const _TaftTE = AkaghefAlgebra.TaftExact
const _TaftNetworks = AkaghefAlgebra.TaftNetworks
const _TaftEvaluation = AkaghefAlgebra.TaftEvaluation
import AkaghefAlgebra: Cyclotomic, RootSpec, normalize_root, root_spec,
    specified_root, canonical_value, with_root, simplify_invariant, evaluate_invariant

"""Taft parameters. Integer N fixes the order; :N requests a symbolic result."""
struct TaftAlgebra{T}
    N::T
    M::Int
    style::Symbol
    function TaftAlgebra(N;M::Integer=1,style=:L)
        style in (:L,:R) || throw(ArgumentError("style must be :L or :R"))
        if N isa Integer
            root=RootSpec(N,M)
            return new{Int}(root.N,root.M,style)
        end
        n=N isa Symbol ? _TaftST.integer(N) : N
        n isa _TaftST.Term && n.op==:var || throw(ArgumentError("N must be an integer or an integer symbol"))
        M==1 || throw(ArgumentError("symbolic N supports M=1; specify another M after fixing N"))
        new{typeof(n)}(n,1,style)
    end
end
root_spec(H::TaftAlgebra{Int})=RootSpec(H.N,H.M)
specified_root(H::TaftAlgebra{Int})=_TaftTE.zeta(H.N)^H.M

"""Build the defining Taft network, without reduction or a choice of M."""
function invariant_expression(v::VirtualLink,H::TaftAlgebra)
    _require_closed_ograph_input(v,"invariant_expression")
    H.style==:L || throw(ArgumentError("symbolic Taft contraction currently supports style=:L"))
    edges=real_edge_table(v)
    weights=v.is_weighted ? _integer_edge_weights(v) : zeros(Int,length(edges))
    signs=last(real_gauss_code(v))
    ports=[4 .* abs.([e.crossing[1],e.crossing[2]]) .+
           ([e.crossing[1]>0,e.crossing[2]>0]) .- [1,3] for e in edges]
    ctx,vars,body=_TaftNetworks.network(ports,signs,weights)
    _TaftNetworks.TaftExpression(body,vars,ctx;ports=ports,weights=weights)
end

"The existing tensor expression for algebras without a symbolic reduction backend."
struct MSTTensorExpression
    factors::Vector{SparseFactor}
end
invariant_expression(v::VirtualLink,H::FiniteHopfAlgebra)=
    MSTTensorExpression(_hopf_invariant_factors(v,H))
function simplify_invariant(x::MSTTensorExpression;strategy=:auto)
    strategy in (:auto,:qa,:network,:none,:numeric) || throw(ArgumentError("unknown MST strategy"))
    # No algebra-specific rewrite is available here; preserve the definition.
    MSTTensorExpression(deepcopy(x.factors))
end
evaluate_invariant(x::MSTTensorExpression)=_contract(copy(x.factors))

function _canonical_taft(H::FiniteHopfAlgebra{TaftKind{N,M,S}}) where {N,M,S}
    prototype=get(_UQSL2_BOREL_CACHE,(M//N,S),nothing)
    prototype!==nothing && H.dimension==prototype.dimension && all(name->begin
        a=getfield(H,name); b=getfield(prototype,name)
        a.variables==b.variables && a.terms==b.terms
    end,fieldnames(typeof(H))[2:end])
end
function invariant_expression(v::VirtualLink,H::FiniteHopfAlgebra{TaftKind{N,M,S}}) where {N,M,S}
    S==:L && _canonical_taft(H) && return invariant_expression(v,TaftAlgebra(N;style=S))
    MSTTensorExpression(_hopf_invariant_factors(v,H))
end

"""Compute the MST invariant. Taft outputs carry N,M and canonical ζ_N coefficients.

Taft strategies :qa (default), :network and :none preserve the final value.
Symbolic N returns a TaftExpression; a concrete evaluation failure throws.
"""
function Z_MST(v::VirtualLink,H::TaftAlgebra;strategy=:qa,details=false,maxterms=3_000_000)
    reduced=simplify_invariant(invariant_expression(v,H);strategy=strategy)
    H.N isa Integer || return reduced
    run=@timed evaluate_invariant(reduced,H.N;maxterms=maxterms)
    value=normalize_root(run.value,root_spec(H))
    details ? (value=value,base_value=run.value,reduced=reduced,evaluation_seconds=run.time) : value
end
function Z_MST(v::VirtualLink,H::FiniteHopfAlgebra;strategy=:auto,details=false,maxterms=3_000_000)
    details && throw(ArgumentError("details is supported by the exact Taft backend"))
    evaluate_invariant(simplify_invariant(invariant_expression(v,H);strategy=strategy))
end
function Z_MST(v::VirtualLink,H::FiniteHopfAlgebra{TaftKind{N,M,S}};
               strategy=:qa,details=false,maxterms=3_000_000) where {N,M,S}
    if S==:L && _canonical_taft(H)
        return Z_MST(v,TaftAlgebra(N;M=M,style=S);strategy=strategy,details=details,maxterms=maxterms)
    end
    details && throw(ArgumentError("details requires unmodified style-L Taft data"))
    evaluate_invariant(simplify_invariant(MSTTensorExpression(_hopf_invariant_factors(v,H));strategy=strategy))
end

# Compatibility entrances use the same stages as Z_MST.
function taft_reduce(v::VirtualLink;style=:L,strategy=:qa)
    simplify_invariant(invariant_expression(v,TaftAlgebra(:N;style=style));strategy=strategy)
end
function taft_evaluate(reduced,N::Integer;M::Integer=1,maxterms=3_000_000)
    root=RootSpec(N,M)
    normalize_root(evaluate_invariant(reduced,root.N;maxterms=maxterms),root)
end
function hopf_invariant(v::VirtualLink,ratio::Rational;style=:L,strategy=:qa,details=false,maxterms=3_000_000)
    0<ratio<1 || throw(ArgumentError("ratio must lie in (0,1)"))
    Z_MST(v,TaftAlgebra(denominator(ratio);M=numerator(ratio),style=style);
          strategy=strategy,details=details,maxterms=maxterms)
end
hopf_invariant(v::VirtualLink,H::TaftAlgebra;kwargs...)=Z_MST(v,H;kwargs...)
function hopf_invariant(v::VirtualLink,H::FiniteHopfAlgebra{TaftKind{N,M,S}};
                        method=:auto,strategy=:qa,details=false,maxterms=3_000_000) where {N,M,S}
    method in (:auto,:exact,:numeric) || throw(ArgumentError("method must be :auto, :exact or :numeric"))
    if method==:numeric
        details && throw(ArgumentError("details is available for exact Taft evaluation"))
        return _hopf_invariant_numeric(v,H)
    end
    method==:exact && !(S==:L && _canonical_taft(H)) &&
        throw(ArgumentError("exact path requires unmodified style-L Taft tensors"))
    Z_MST(v,H;strategy=strategy,details=details,maxterms=maxterms)
end

export Z_MST, TaftAlgebra, Cyclotomic, RootSpec, invariant_expression,
       simplify_invariant, evaluate_invariant, normalize_root, root_spec,
       specified_root, canonical_value, with_root, taft_reduce, taft_evaluate
