const _TaftST = AkaghefAlgebra.SymbolicTaft
const _TaftTE = AkaghefAlgebra.TaftExact
const _TaftNetworks = AkaghefAlgebra.TaftNetworks
const _TaftEvaluation = AkaghefAlgebra.TaftEvaluation

"""Reduce a closed, integer-weighted style-L Taft network without fixing N.

Unweighted input uses zero edge weights, as in generic `hopf_invariant`.
The returned expression can be reused at different N without rebuilding the graph.
"""
function taft_reduce(v::VirtualLink;style=:L)
    _require_closed_ograph_input(v,"taft_reduce")
    style==:L || throw(ArgumentError("symbolic Taft contraction currently supports style=:L"))
    edges=real_edge_table(v)
    weights=v.is_weighted ? _integer_edge_weights(v) : zeros(Int,length(edges))
    signs=last(real_gauss_code(v))
    ports=[4 .* abs.([e.crossing[1],e.crossing[2]]) .+
           ([e.crossing[1]>0,e.crossing[2]>0]) .- [1,3] for e in edges]
    ctx,vars,body=_TaftNetworks.network(ports,signs,weights)
    run=@timed _TaftNetworks.reduce_with_qrules(body,vars,ctx)
    reduced,remaining,trace=run.value
    (expression=_TaftNetworks.close_sums(reduced,remaining,ctx),body=reduced,
     context=ctx,remaining=remaining,trace=trace,seconds=run.time,
     ports=ports,weights=weights)
end

"""Evaluate a reduced Taft expression exactly in Q(ζ_N), with q mapped to ζ_N^M."""
function taft_evaluate(reduced,N::Integer;M::Integer=1,maxterms=3_000_000)
    N>=2 || throw(ArgumentError("Taft order must be at least 2"))
    gcd(N,M)==1 || throw(ArgumentError("q must be a primitive N-th root"))
    value=_TaftEvaluation.evaluate_fast(reduced.expression,reduced.context,
                                        Dict(:N=>N);maxterms=maxterms)
    value=convert(_TaftTE.Cyclotomic{Int(N)},value)
    _TaftNetworks.galois(value,Int(N),mod(M,N))
end

"""`hopf_invariant(v, M//N; style=:L, details=false)` computes the exact Taft MST value.

This entry avoids constructing numerical structure tensors. With `details=true`,
return the value, residual expression, and separate reduction/evaluation timings.
"""
function hopf_invariant(v::VirtualLink,ratio::Rational;style=:L,details=false,maxterms=3_000_000)
    0<ratio<1 || throw(ArgumentError("ratio must lie in (0,1)"))
    reduced=taft_reduce(v;style=style)
    run=@timed taft_evaluate(reduced,denominator(ratio);M=numerator(ratio),maxterms=maxterms)
    details ? (value=run.value,reduced=reduced,evaluation_seconds=run.time) : run.value
end

function hopf_invariant(v::VirtualLink,algebra::FiniteHopfAlgebra{TaftKind{N,M,S}};
                        method=:auto,details=false,maxterms=3_000_000) where {N,M,S}
    method in (:auto,:exact,:numeric) || throw(ArgumentError("method must be :auto, :exact or :numeric"))
    # Public tensors are mutable. Modified algebras must retain their actual
    # tensor semantics rather than silently evaluating the original Taft model.
    prototype=get(_UQSL2_BOREL_CACHE,(M//N,S),nothing)
    canonical=prototype!==nothing && algebra.dimension==prototype.dimension && all(name->begin
        a=getfield(algebra,name); b=getfield(prototype,name)
        a.variables==b.variables && a.terms==b.terms
    end,fieldnames(typeof(algebra))[2:end])
    if method==:numeric || S!=:L || !canonical
        method==:exact && throw(ArgumentError("exact path requires unmodified style-L Taft tensors"))
        details && throw(ArgumentError("details is available for exact Taft evaluation"))
        return _hopf_invariant_numeric(v,algebra)
    end
    hopf_invariant(v,M//N;style=S,details=details,maxterms=maxterms)
end

export taft_reduce, taft_evaluate
