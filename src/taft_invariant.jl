const _TaftST = AkaghefAlgebra.SymbolicTaft
const _TaftTE = AkaghefAlgebra.TaftExact
const _TaftNetworks = AkaghefAlgebra.TaftNetworks
const _TaftEvaluation = AkaghefAlgebra.TaftEvaluation
import AkaghefAlgebra: Cyclotomic, RootSpec, normalize_root, root_spec,
    specified_root, canonical_value, with_root, simplify_invariant, evaluate_invariant

"""
    TaftAlgebra(N; M::Integer=1, style=:L) -> TaftAlgebra

Select concrete or symbolic Taft parameters for the o-graph contraction.

# Arguments
Requires H7 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
The order, primitive exponent and style checks raise `ArgumentError` naming
the invalid order, nonprimitive root, unsupported style or symbolic exponent.
Other symbolic/type conversion errors propagate. Style `:R` is accepted by
this constructor but is not supported by its network evaluator.

# Returns
A `TaftAlgebra` with `N`, normalized `M` and `style` in the
[exact Taft coordinate domain](../specs/hopf-invariants.md#exact-taft-route).

# Guarantees
1. A concrete selector has the root data of `RootSpec(N,M)`; a symbolic selector retains its integer-variable term. [proved/agent: by construction]
2. Construction does not mutate inputs or evaluate an invariant. [proved/agent: by construction]

# Not guaranteed
The selector does not contain the full structure tensors or certify the
geometric hypotheses H2–H3. A symbolic order does not produce a scalar value.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Exact Taft route — root/selector data and supported coordinate model.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> Z_MST(v,uqsl2_borel_small(1//5)) == Z_MST(v,TaftAlgebra(5))
true
```

Role: construct
"""
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

"""
    invariant_expression(v::VirtualLink, H::TaftAlgebra) -> AkaghefAlgebra.TaftExpression

Return the symbolic Taft coordinate contraction associated with the o-graph data.

# Arguments
Requires H2–H4 and H7 ([conditions](../specs/hopf-invariants.md#numbered-conditions))
for the intended geometric interpretation. Checked: trailing cut markers,
assigned integer-weight mode and finite integral real-edge weights, style
`:L`, valid crossing signs and complete port incidence. Violations raise
`ArgumentError` naming the failure. Full closed-normality and the framing
equation are caller obligations. Conversion/indexing failures propagate.

# Returns
An `AkaghefAlgebra.TaftExpression` with body, closed expression, residual sum
indices, context, trace/timing metadata, ports and weights, in the
[exact Taft route](../specs/hopf-invariants.md#exact-taft-route).

# Guarantees
1. The expression is the defining symbolic coordinate network for the supplied real-Gauss data and edge weights; no concrete order or root exponent is substituted. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
The expression is not a scalar or a certificate of H2–H3. The general port,
normalization and source-invariant identification remains unverified.

# Sources
cf:[[MST25#^mst25-blk-0025]] the Heisenberg-double contraction target; the repository coordinate comparison has the specification's evidence boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> x = simplify_invariant(invariant_expression(v,TaftAlgebra(:N)); strategy=:qa);

julia> isempty(x.remaining)
true
```

Role: convert
"""
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
"""
    invariant_expression(v::VirtualLink, H::FiniteHopfAlgebra) -> VirtualLinks.MSTTensorExpression

Return the finite structural-tensor network associated with the o-graph data.

# Arguments
Requires H1–H3 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
The input and exception boundary is that of `hopf_invariant`'s numeric route:
cut markers and invalid assigned integer weights raise `ArgumentError`.
Coherent Hopf tensors and geometric/framing hypotheses are caller obligations.
Singular inverse antipode powers, indexing and arithmetic errors propagate.

# Returns
A `VirtualLinks.MSTTensorExpression` with `factors` representing the
[structural-tensor contraction](../specs/hopf-invariants.md#structural-tensor-contraction).

# Guarantees
1. The factors define the supported numeric coordinate network of the supplied diagram and algebra. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
No Hopf-axiom or general manifold-invariance certification is performed.
Floating arithmetic and tolerance-based coefficient omission do not supply
an exact algebraic certificate or a certified error bound.

# Sources
cf:[[MST25#^mst25-blk-0025]] structural-tensor target; full source-port correspondence remains unverified.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> evaluate_invariant(invariant_expression(v,sweedler_algebra())) ≈ -1
true
```

Role: convert
"""
invariant_expression(v::VirtualLink,H::FiniteHopfAlgebra)=
    MSTTensorExpression(_hopf_invariant_factors(v,H))
function simplify_invariant(x::MSTTensorExpression;strategy=:auto)
    strategy in (:auto,:qa,:network,:none,:numeric) || throw(ArgumentError("unknown MST strategy"))
    # No algebra-specific rewrite is available here; preserve the definition.
    MSTTensorExpression(deepcopy(x.factors))
end
evaluate_invariant(x::MSTTensorExpression)=_contract(copy(x.factors))

function _canonical_taft(H::FiniteHopfAlgebra{TaftKind{N,M,S}}) where {N,M,S}
    # cf:hopf_invariant.jl _UQSL2_BOREL_CACHE — owned by the algebra's defining module, including split-module layouts.
    cache=getfield(parentmodule(FiniteHopfAlgebra),:_UQSL2_BOREL_CACHE)
    prototype=get(cache,(M//N,S),nothing)
    prototype!==nothing && H.dimension==prototype.dimension && all(name->begin
        a=getfield(H,name); b=getfield(prototype,name)
        a.variables==b.variables && a.terms==b.terms
    end,fieldnames(typeof(H))[2:end])
end
"""
    invariant_expression(v::VirtualLink, H::FiniteHopfAlgebra{TaftKind{N,M,S}}) -> Union{AkaghefAlgebra.TaftExpression,VirtualLinks.MSTTensorExpression}

Return the supported exact or numeric Taft-tagged coordinate network of an o-graph.

# Arguments
Requires H1–H4 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
Input checks and errors are those of the corresponding Taft selector or
generic finite-Hopf `invariant_expression` contract. The raw kind tag alone
does not validate tensor contents; geometric hypotheses are caller obligations.

# Returns
The symbolic expression for the supported unmodified style-L prototype, or
the finite tensor expression otherwise, as defined by the
[exact Taft bridge](../specs/hopf-invariants.md#exact-taft-route).

# Guarantees
1. The result follows that bridge and reflects the supplied tensor contents. Its symbolic branch leaves order and root exponent for later evaluation. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
A Taft tag does not guarantee use of exact arithmetic. Mutated tensors or
other styles can select the numeric path; neither path certifies H2–H3.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Exact Taft route — supported prototype and numeric fallback.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> x = simplify_invariant(invariant_expression(v,uqsl2_borel_small(1//5)); strategy=:qa);

julia> isempty(x.remaining)
true
```

Role: convert
"""
function invariant_expression(v::VirtualLink,H::FiniteHopfAlgebra{TaftKind{N,M,S}}) where {N,M,S}
    S==:L && _canonical_taft(H) && return invariant_expression(v,TaftAlgebra(N;style=S))
    MSTTensorExpression(_hopf_invariant_factors(v,H))
end

"""
    Z_MST(v::VirtualLink, H::TaftAlgebra; strategy=:qa, details=false, maxterms=3000000) -> Any

Evaluate or symbolically reduce the o-graph's Taft coordinate contraction.

# Arguments
Requires H2–H4 and H7 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
The domain checks and failures of `invariant_expression(v,H)` apply.
`strategy` is `:qa`, `:network`, `:none` or `:auto`; another value raises
`ArgumentError`. The concrete evaluator checks its nonnegative `Int` term
budget and raises `ArgumentError` on invalid/exceeded budgets. Symbolic
evaluation does not consume `maxterms`. Unsupported expression operations
and arithmetic errors propagate. `details` is used as a Boolean on the
concrete branch.

# Returns
For concrete `N`, a rooted `Cyclotomic` value, or with `details=true` a record
`(value, base_value, reduced, evaluation_seconds)`: rooted value, canonical
root value, reduced expression, and elapsed evaluation time. A symbolic order
returns the reduced `AkaghefAlgebra.TaftExpression` regardless of `details`.

# Guarantees
1. The concrete value is the exact coordinate contraction evaluated at the selected primitive root. The symbolic result retains the residual sums. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
Exact arithmetic does not establish the unresolved geometric and source-port
identifications H2–H3. The budget is not a wall-clock limit and timing is not
reproducible metadata. Symbolic reduction need not eliminate every sum.

# Sources
cf:[[MST25#^mst25-blk-0025]] and [[MST25#^mst25-blk-0027]] — literature contraction and additional manifold hypotheses, subject to the coordinate comparison in the Hopf specification.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> Z_MST(v,TaftAlgebra(5; M=2)) isa Cyclotomic{5,2}
true
```

Role: invariant
"""
function Z_MST(v::VirtualLink,H::TaftAlgebra;strategy=:qa,details=false,maxterms=3_000_000)
    reduced=simplify_invariant(invariant_expression(v,H);strategy=strategy)
    H.N isa Integer || return reduced
    run=@timed evaluate_invariant(reduced,H.N;maxterms=maxterms)
    value=normalize_root(run.value,root_spec(H))
    details ? (value=value,base_value=run.value,reduced=reduced,evaluation_seconds=run.time) : value
end
"""
    Z_MST(v::VirtualLink, H::FiniteHopfAlgebra; strategy=:auto, details=false, maxterms=3000000) -> Any

Evaluate the finite structural-tensor contraction of an o-graph.

# Arguments
Requires H1–H3 ([conditions](../specs/hopf-invariants.md#numbered-conditions))
and the domain of `invariant_expression(v,H)`. `details=true` raises
`ArgumentError` identifying the exact Taft backend. `strategy` accepts
`:auto`, `:qa`, `:network`, `:none`, `:numeric`; unsupported names raise
`ArgumentError`. `maxterms` is unused by this numeric method. Tensor,
coefficient and allocation failures propagate.

# Returns
The scalar contraction in the
[structural-tensor model](../specs/hopf-invariants.md#structural-tensor-contraction),
with type determined by the tensor arithmetic.

# Guarantees
1. The value is the evaluation of the finite expression for these inputs; accepted strategies preserve that expression's value. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
There is no exactness, floating error bound or effective term budget here.
General manifold invariance requires the unverified source-coordinate and
geometric comparisons in H2–H3.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Structural-tensor contraction and MST geometric target.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> Z_MST(v,sweedler_algebra()) ≈ -1
true
```

Role: invariant
"""
function Z_MST(v::VirtualLink,H::FiniteHopfAlgebra;strategy=:auto,details=false,maxterms=3_000_000)
    details && throw(ArgumentError("details is supported by the exact Taft backend"))
    evaluate_invariant(simplify_invariant(invariant_expression(v,H);strategy=strategy))
end
"""
    Z_MST(v::VirtualLink, H::FiniteHopfAlgebra{TaftKind{N,M,S}}; strategy=:qa, details=false, maxterms=3000000) -> Any

Evaluate the exact or numeric coordinate contraction selected by the supplied Taft-tagged tensors.

# Arguments
Requires H1–H4 ([conditions](../specs/hopf-invariants.md#numbered-conditions)).
The supported unmodified style-L prototype uses the `TaftAlgebra` method's
keyword checks and exact term budget. Other tensors use the numeric method's
checks, ignore `maxterms` and reject `details=true` with `ArgumentError`
requiring unmodified style-L data. Input, conversion and evaluation errors
propagate as in those methods.

# Returns
The rooted exact scalar or details record on the exact branch; otherwise the
numeric contraction scalar, following the
[Taft bridge](../specs/hopf-invariants.md#exact-taft-route).

# Guarantees
1. The selected branch uses the supplied tensors; the exact branch evaluates the kind tag's normalized primitive root. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]

# Not guaranteed
A Taft kind tag does not certify exact dispatch or the Hopf axioms. Numeric
fallback has no certified error bound, and general geometric invariant
identification remains unverified for both branches.

# Sources
cf:lib/topology.jl/docs/specs/hopf-invariants.md Exact Taft route — prototype comparison, selected root and fallback.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v,[[-2,1]]);

julia> Z_MST(v,uqsl2_borel_small(2//5)) isa Cyclotomic{5,2}
true
```

Role: invariant
"""
function Z_MST(v::VirtualLink,H::FiniteHopfAlgebra{TaftKind{N,M,S}};
               strategy=:qa,details=false,maxterms=3_000_000) where {N,M,S}
    if S==:L && _canonical_taft(H)
        return Z_MST(v,TaftAlgebra(N;M=M,style=S);strategy=strategy,details=details,maxterms=maxterms)
    end
    details && throw(ArgumentError("details requires unmodified style-L Taft data"))
    evaluate_invariant(simplify_invariant(MSTTensorExpression(_hopf_invariant_factors(v,H));strategy=strategy))
end

# Compatibility entrances use the same stages as Z_MST.
"""
    taft_reduce(v::VirtualLink; style=:L, strategy=:qa) -> AkaghefAlgebra.TaftExpression

Return a symbolic finite-sum expression for the closed o-graph's Taft contraction.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

`v` has the closed-normal and integral-weight caller obligations of
[`hopf_invariant`](@ref). Cut markers and invalid assigned weights are
rejected with `ArgumentError`, but full closedness and the framing equation
are not certified. The current route requires `style=:L`; unsupported styles
raise `ArgumentError`. `strategy` is `:qa`, `:network`, `:none`, or `:auto`;
an unknown strategy raises `ArgumentError`. Unsupported symbolic operations
and inconsistent network data propagate their exceptions.

# Returns
`AkaghefAlgebra.TaftExpression` with symbolic expression/body, remaining sum
indices, algebraic context, reduction trace and timing, and the input network's
port and weight data. The order remains symbolic and no root parameter is fixed.

# Guarantees
1. The result represents the repository's Taft contraction before a concrete
   order and primitive root are chosen. [proved/agent: by construction]
2. The source diagram is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
All sums need not be eliminated. This is not a numerical value or a
certificate of diagram-independent, framed-manifold invariance; geometric
hypotheses and the full source-convention dictionary are not checked.

# Sources
cf:[[MST25#^mst25-blk-0025]] tensor-network target of the symbolic calculation.
cf:lib/algebra.jl/src/Hopf/TaftNetworks.jl TaftExpression and simplify_invariant — symbolic network and residual finite sums.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]);

julia> isempty(taft_reduce(v).remaining)
true
```

See also [`taft_evaluate`](@ref).

Role: convert
"""
function taft_reduce(v::VirtualLink;style=:L,strategy=:qa)
    simplify_invariant(invariant_expression(v,TaftAlgebra(:N;style=style));strategy=strategy)
end
"""
    taft_evaluate(reduced, N::Integer; M::Integer=1, maxterms=3000000) -> Cyclotomic

Evaluate a symbolic Taft contraction exactly at a specified primitive root.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

`reduced` is an expression supported by `evaluate_invariant(reduced,N)`,
normally returned by [`taft_reduce`](@ref). `N` and `M` must be accepted by
[`RootSpec`](@ref); invalid order or nonprimitive root data raise
`ArgumentError`. `maxterms` is a nonnegative `Int`; a negative budget or
exhausted summand budget raises `ArgumentError`. Unsupported expression
types, unspecialized symbols and arithmetic failures propagate their errors.

# Returns
An exact `Cyclotomic` value, with the normalized order/root metadata and
coefficients in the canonical cyclotomic field.

# Guarantees
1. On supported expressions, the result evaluates the expression at `N`
   and applies the primitive-root specialization selected by `M` once.
   [proved/agent: by construction]
2. The input expression is not modified by the supported Taft evaluator.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Exact arithmetic does not validate the geometric origin of an expression
or establish topological invariance. The summand budget does not bound all
symbolic work, memory, or elapsed time.

# Sources
cf:lib/algebra.jl/src/Hopf/TaftEvaluation.jl evaluate_invariant — exact residual finite-sum evaluation.
cf:lib/algebra.jl/src/CyclotomicValues.jl RootSpec and normalize_root — exact selected-root interpretation.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]); reduced = taft_reduce(v);

julia> taft_evaluate(reduced, 9) == specified_root(RootSpec(9))^-1
true
```

Role: invariant
"""
function taft_evaluate(reduced,N::Integer;M::Integer=1,maxterms=3_000_000)
    root=RootSpec(N,M)
    normalize_root(evaluate_invariant(reduced,root.N;maxterms=maxterms),root)
end
"""
    hopf_invariant(v::VirtualLink, ratio::Rational; style=:L, strategy=:qa, details=false, maxterms=3000000) -> Union{Cyclotomic,NamedTuple}

Evaluate the closed o-graph Taft contraction at the primitive root specified by a rational parameter.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

- `ratio` must lie strictly between zero and one; otherwise `ArgumentError`
  reports that interval. Its denominator and numerator select the order
  and root parameter accepted by [`RootSpec`](@ref); order conversion errors
  are reported as `ArgumentError`.
- `v` has the closed-normal and integral-weight caller obligations of
  [`hopf_invariant(::VirtualLink, ::FiniteHopfAlgebra)`](@ref hopf_invariant).
  Cut markers and nonintegral assigned weights are rejected; full closedness
  and the framing equation are not checked.
- The current symbolic Taft route supports `style=:L`; other styles raise
  `ArgumentError` identifying the unsupported style. `strategy` is `:qa`,
  `:network`, `:none`, or `:auto`; an unknown strategy raises `ArgumentError`.
- `details` is Boolean. `maxterms` is a nonnegative `Int` finite-sum budget;
  a negative budget or exceeded budget raises `ArgumentError`, and arithmetic
  or unsupported symbolic-evaluation failures propagate.

# Returns
An exact `Cyclotomic` value carrying root metadata. With `details=true`,
a `NamedTuple` contains `value`, the canonical-root `base_value`, the reduced
expression `reduced`, and `evaluation_seconds` for the final evaluation.

# Guarantees
1. On successful supported evaluation, the result is the repository's Taft
   network value at the selected primitive root, with exact cyclotomic
   coefficients. [proved/agent: by construction]
2. The input diagram is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
An exact algebraic value is not by itself a proof of framed-manifold
invariance. The source's framing, first-Betti-number and normalization
hypotheses are not certified here, nor is the entire port/basis translation
to that source. The term budget is not a bound on runtime or all reductions.

# Sources
cf:[[MST25#^mst25-blk-0025]] Heisenberg-double network target.
cf:[[MST25#^mst25-blk-0027]] geometric scope of the literature invariant; exact numerical evaluation alone does not check its hypotheses.
cf:lib/algebra.jl/src/CyclotomicValues.jl normalize_root — exact primitive-root specialization of the calculated value.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]);

julia> hopf_invariant(v, 1//5) == specified_root(TaftAlgebra(5))^-1
true
```

See also [`taft_reduce`](@ref), [`taft_evaluate`](@ref).

Role: invariant
"""
function hopf_invariant(v::VirtualLink,ratio::Rational;style=:L,strategy=:qa,details=false,maxterms=3_000_000)
    0<ratio<1 || throw(ArgumentError("ratio must lie in (0,1)"))
    Z_MST(v,TaftAlgebra(denominator(ratio);M=numerator(ratio),style=style);
          strategy=strategy,details=details,maxterms=maxterms)
end
"""
    hopf_invariant(v::VirtualLink, H::TaftAlgebra; strategy=:qa, details=false, maxterms=3000000) -> Union{Cyclotomic,NamedTuple,AkaghefAlgebra.TaftExpression}

Return the exact or symbolic Taft contraction of a closed o-graph.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

`v` satisfies the closed-normal and weight obligations of
[`hopf_invariant(::VirtualLink, ::FiniteHopfAlgebra)`](@ref hopf_invariant),
and `H` is a `TaftAlgebra`. Cut markers, invalid assigned integer weights,
unsupported style and unknown strategy raise `ArgumentError` as in the
rational-parameter overload. Keywords and failures are those of
[`Z_MST(::VirtualLink, ::TaftAlgebra)`](@ref Z_MST).
Concrete evaluation requires a nonnegative `Int` term budget; exceeded
budgets and symbolic-evaluation errors propagate.

# Returns
For integer order, an exact `Cyclotomic`, or the `value`, `base_value`,
`reduced`, `evaluation_seconds` tuple when `details=true`. For symbolic order,
an `AkaghefAlgebra.TaftExpression` is returned even when `details=true`.

# Guarantees
1. The result has the algebraic semantics of `Z_MST(v,H; kwargs...)`:
   an exact selected-root value or a symbolic finite-sum expression for
   the repository's network. [proved/agent: by construction]
2. The input diagram is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
A symbolic result need not have all sums eliminated. Geometric input
hypotheses and the full dictionary to the literature invariant are not
certified; no unconditional knot, combed-manifold or framed-manifold
invariance follows from successful evaluation alone.

# Sources
cf:[[MST25#^mst25-blk-0025]] network interpretation targeted by the algebraic calculation.
cf:[[MST25#^mst25-blk-0027]] additional geometric hypotheses for the literature invariant.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]);

julia> hopf_invariant(v, TaftAlgebra(5)) == specified_root(TaftAlgebra(5))^-1
true
```

Role: invariant
"""
hopf_invariant(v::VirtualLink,H::TaftAlgebra;kwargs...)=Z_MST(v,H;kwargs...)

"""
    hopf_invariant(v::VirtualLink, H::FiniteHopfAlgebra{TaftKind{N,M,S}}; method=:auto, strategy=:qa, details=false, maxterms=3000000) -> Any

Evaluate a closed o-graph with tagged finite Taft data through the applicable exact or numeric route.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/hopf-invariants.md).

`v` and coherent finite tensors satisfy the caller obligations of
[`hopf_invariant(::VirtualLink, ::FiniteHopfAlgebra)`](@ref hopf_invariant).
`method` must be `:auto`, `:exact`, or `:numeric`; invalid choices raise
`ArgumentError`. `:exact` additionally requires unmodified style-L Taft
tensors matching the cached constructor prototype, or raises `ArgumentError`.
`details=true` is rejected on the numeric route or when the exact route is
unavailable. On the exact route `strategy` and `maxterms` have the Taft
overload's checks and failures; forced numeric evaluation does not use them.

# Returns
An exact rooted cyclotomic value, its optional detail tuple, or the scalar
contraction of the supplied finite tensors on a fallback/numeric route.
The return type is `Any` because the fallback follows the supplied tensors.

# Guarantees
1. `method=:numeric` evaluates the supplied finite tensor data. Eligible
   unmodified style-L inputs on the exact route use their tagged root order
   and parameter. [proved/agent: by construction]
2. With `:auto`, modified or ineligible tensors retain their actual supplied
   contraction semantics instead of acquiring the canonical exact value.
   [proved/agent: by construction]
3. The supplied diagram and tensor data are not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The type tag alone does not prove the tensors are canonical or valid Hopf
data. Numeric results have no certified error bound. Neither route checks
all geometric and framing hypotheses or proves the repository's network
dictionary matches the literature invariant in complete generality.

# Sources
cf:[[MST25#^mst25-blk-0025]] Heisenberg-double tensor-network target.
cf:[[MST25#^mst25-blk-0027]] literature invariant hypotheses beyond scalar contraction.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]); set_weight!(v, [[-2,1]]); h = uqsl2_borel_small(1//5);

julia> ComplexF64(hopf_invariant(v, h)) ≈ hopf_invariant(v, h; method=:numeric)
true
```

Role: invariant
"""
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

include("group_mst.jl")
include("local_hopf_adapter.jl")
