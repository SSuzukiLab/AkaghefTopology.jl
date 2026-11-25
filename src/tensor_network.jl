"""
    label_tensor(tensor::SparseFactor, labels) -> SparseFactor

Return a labelled coordinate tensor with repeated labels restricted to their common diagonal.

# Arguments
- `tensor` supplies coherent coefficient keys and arithmetic; `labels` must
  collect to integer labels convertible to `Int`, one per coefficient axis.
  Conversion errors propagate. Key arity is not validated; a mismatch can
  produce indexing errors or data without the intended tensor interpretation.
- Labels and basis interpretations are supplied by the caller in the
  [Algebra tensor API](../specs/algebra.md#definitions); no variance or
  dimension consistency check is performed.

# Returns
A `SparseFactor` with the input coefficient type, and distinct labels in
their first-occurrence order. The original factor's labels are replaced.

# Guarantees
1. For coherent input, distinct labels preserve the coefficient tensor;
   repeated labels select its corresponding diagonal coordinates without
   summing over the remaining label. [proved/agent: by construction]
2. Inputs are not modified. The label vector is new. With distinct labels
   the coefficient dictionary can be shared with the input; repeated labels
   produce a new dictionary. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Labelling is not contraction over all labels, tensor-identity validation,
or automatic inference of a primal/dual convention. Coefficients that are
mutable objects are not deep-copied.

# Sources
cf:[[Kup91#^kup91-sec-0018]] coordinate indices for tensors; repeated-label diagonal extraction is the repository network convention.

# Examples
```jldoctest
julia> a = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>1//2, (1,2)=>9//10, (2,2)=>2//3));

julia> r = label_tensor(a,[1,1]); (r.terms[(1,)], r.terms[(2,)])
(1//2, 2//3)
```

See also [`contract_network`](@ref), [`dense_tensor`](@ref).

Role: convert
"""
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

"""
    contract_network(factors, outputs=Int[]) -> SparseFactor

Return the contraction of a labelled sparse tensor network with the requested free indices.

# Arguments
- `factors` must collect to an `AbstractVector{<:SparseFactor}`. Each factor
  must have coherent labels, coefficient keys and arithmetic; these structural
  preconditions are the caller's responsibility. Malformed keys or unsupported
  coefficient operations can propagate indexing or arithmetic exceptions.
- `outputs` supplies distinct free labels convertible to `Int`, in the desired
  axis order. For a nonempty network, absent output labels raise `ArgumentError`
  stating that every output index must occur in the network. Distinctness is
  not checked by this entry point.
- Compatible basis and dual-basis interpretations are the caller's
  responsibility; integer labels alone do not certify tensor variance. See
  the [Algebra tensor API](../specs/algebra.md#definitions).

# Returns
A `SparseFactor` whose `variables` are the requested output labels and whose
`terms` are the contracted coefficients. With no free indices this is a
rank-zero tensor, not a scalar value.

# Guarantees
1. For coherent labels and supported exact arithmetic without overflow, the
   result is the coordinate contraction defined in
   [coordinate tensors](../specs/algebra.md#coordinate-tensors), with the
   requested outputs. [proved/agent: by construction]
2. Coefficient types are promoted by their Julia arithmetic without an
   unconditional floating-point conversion. [proved/agent: by construction]
3. Input factors, their containers and the output-label collection are not
   modified. [proved/agent: by construction]
4. An empty typed factor collection with empty outputs returns the unit
   rank-zero tensor with a `ComplexF64` coefficient. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Floating coefficients can lose products of magnitude at most `1e-12`;
there is no error bound. This call does not enforce the MATLAB index rule,
check dimensions, or establish a topological interpretation. Empty networks
with nonempty outputs and duplicate output labels are not validated and
do not have the ordinary free-tensor interpretation asserted above.

# Sources
cf:[[Kup91#^kup91-sec-0018]] finite-dimensional tensor products and contractions; sparse labels are the repository representation.

# Examples
```jldoctest
julia> a = SparseFactor([1,2], Dict{Tuple,Rational{Int}}((1,1)=>1//2, (1,2)=>1//3, (2,1)=>2//5));

julia> b = SparseFactor([2,3], Dict{Tuple,Rational{Int}}((1,1)=>3//7, (2,1)=>5//11));

julia> r = contract_network([a,b], [1,3]); (r.terms[(1,1)], r.terms[(2,1)])
(169//462, 6//35)
```

See also [`label_tensor`](@ref), [`evaluate_tensor_expression`](@ref).

Role: convert
"""
contract_network(factors,outputs=Int[])=_contract_to(collect(factors),Int.(collect(outputs)))

"""
    dense_tensor(factor::SparseFactor{T}, dimension::Integer) -> Array

Return a dense coefficient array for a labelled sparse tensor.

# Arguments
`dimension` must be a nonnegative integer convertible to `Int`. Each tensor
coordinate must be a valid array index, and keys should have one entry per
label in `factor.variables`. Invalid dimensions or conversions propagate
allocation/conversion exceptions; invalid array indices propagate
`BoundsError`. The input's shape is not separately validated.

# Returns
An `Array{T}` with one axis of length `dimension` per label. A tensor with
no labels has a zero-dimensional array result, rather than a scalar.
The tensor API belongs to [Algebra](../specs/algebra.md#definitions).

# Guarantees
1. For well-formed keys, each stored coefficient occupies its corresponding
   array entry and all unstored entries are `zero(T)`. [proved/agent: by construction]
2. The result is the same tensor in dense coordinates, with axis order
   inherited from the input labels and no summation or relabelling. [proved/agent: by construction]
3. The result is newly allocated, the input is unchanged, and the coefficient
   type is preserved, including exact rational coefficient types. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Equal labels do not trigger contraction or diagonal extraction. An unlabelled
Hopf structure tensor does not acquire its axes automatically; use
[`label_tensor`](@ref) first. No tensor variance is inferred.

# Sources
cf:[[Kup91#^kup91-blk-0024]] tensor entries in a fixed basis; dense array storage preserves the supplied coefficient type.

# Examples
```jldoctest
julia> a = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>1//2, (1,2)=>1//3, (2,1)=>2//5));

julia> b = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>3//7, (2,1)=>5//11));

julia> r = evaluate_tensor_expression("A{1,2}B{2,3}", Dict("A"=>a, "B"=>b), [1,3]);

julia> eltype(dense_tensor(r, 2)) == Rational{Int}
true
```

See also [`SparseFactor`](@ref), [`label_tensor`](@ref).

Role: convert
"""
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

"""
    hopf_tensors(algebra::FiniteHopfAlgebra) -> Dict{String,SparseFactor}

Return the structure-tensor dictionary used for contractions of a finite Hopf algebra.

# Arguments
`algebra` must have coherent structural tensors and a nonnegative dimension
that fits array indices. Its antipode must have valid matrix coordinates and
a coefficient type supporting matrix zeroth powers. Malformed data can
propagate indexing, conversion or matrix-operation errors; Hopf axioms are
not checked.

# Returns
A new `Dict{String,SparseFactor}` containing the eleven stored structural,
crossing and distinguished tensors, together with an identity-map tensor.
Its consumers are [`evaluate_tensor_expression`](@ref) and
[`verify_hopf`](@ref), in the [Algebra API](../specs/algebra.md#definitions).

# Guarantees
1. For coherent input, the dictionary supplies tensor coordinates for
   structural maps on that input's vector space and its dual, and an identity
   map of dimension `algebra.dimension`. [proved/agent: by construction]
2. The eleven stored tensor values are shared with the input; the dictionary
   and identity tensor are new. Replacing a dictionary entry does not replace
   a field in the input, but mutating a shared tensor changes that tensor. [proved/agent: by construction]
3. The call leaves the input unchanged and does not impose involutivity. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
This is not the legacy antipode-transpose convention of
[`legacy_hopf_tensors`](@ref). Returning the dictionary certifies neither
Hopf identities nor topological invariance of an arbitrary contraction.

# Sources
cf:[[Kup91#^kup91-sec-0038]] structural tensors of a finite-dimensional Hopf algebra; string keys and crossing tensors are repository interfaces.

# Examples
```jldoctest
julia> tensors = hopf_tensors(sweedler_algebra());

julia> any(occursin("Mu", message) && occursin("rank", message) for message in validate_tensor_expression("Mu{1,2}", [1,2]; tensors=tensors))
true
```

See also [`legacy_hopf_tensors`](@ref), [`evaluate_tensor_expression`](@ref).

Role: access
"""
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

"""
    legacy_hopf_tensors(algebra::FiniteHopfAlgebra) -> Dict{String,SparseFactor}

Return Hopf-coordinate tensors for the legacy C250406 tensor-expression interface.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

`algebra` must meet [`hopf_tensors`](@ref)'s requirements. Multiplication,
coproduct and antipode data must also be nonempty, coherently indexed factors
accepted by [`evaluate_tensor_expression`](@ref). The caller supplies valid
Hopf data; tensor-validation, indexing, conversion and arithmetic errors
propagate without a separate Hopf-axiom check.

# Returns
A new `Dict{String,SparseFactor}` with the same keys as `hopf_tensors`,
including the legacy antipode, crossing tensors and identity tensor.

# Guarantees
1. For coherent input, the dictionary supplies the repository's legacy
   coordinate presentation of maps on the input vector space and its dual
   for use by the C250406 expression interface. [proved/agent: by construction]
2. The algebra is not modified. Multiplication, coproduct, unit, counit and
   distinguished tensor values are shared with it; the dictionary and the
   legacy antipode, crossing and identity tensors are new.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
The dictionary is not interchangeable with `hopf_tensors` in expressions
whose port conventions differ. No literature basis/port identification,
Hopf-axiom certificate or topological-invariance theorem is returned.
Floating-point contractions have the limitations of `contract_network`.

# Sources
cf:lib/topology.jl/examples/C250406TensorNetwork.jl, run_C250406_tensor_network — the legacy coordinate-expression interface consuming these tensors.

# Examples
```jldoctest
julia> tensors = legacy_hopf_tensors(sweedler_algebra());

julia> expr = "Mu{1,2,3}Il{3}An{4,2}De{7,5,6}De{6,4,8}(An^-2){8,9}Mu{9,12,10}Mu{10,5,11}De{13,14,12}Cl{13}(An){14,15}";

julia> isempty(validate_tensor_expression(expr,[11,1,7,15]; tensors=tensors))
true
```

See also [`hopf_tensors`](@ref), [`evaluate_tensor_expression`](@ref).

Role: access
"""
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

"""
    validate_tensor_expression(expression::AbstractString, outputs=Int[]; tensors=nothing, index_rule::Symbol=:matlab) -> Vector{String}

Return syntax, index and optional tensor-rank diagnostics for a coordinate tensor expression.

# Arguments
- `expression` is a string; invalid expression syntax is reported in the
  returned diagnostics rather than evaluated as Julia code.
- `outputs` is a collectable sequence of labels convertible to `Int`.
  Conversion and integer-parsing errors propagate.
- `index_rule` is `:matlab` or `:network`; other values raise `ArgumentError`
  identifying these choices. These select the legacy expression or extended
  factor-network conditions A5–A6 of the
  [Algebra specification](../specs/algebra.md#numbered-conditions).
- `tensors` is `nothing` or a name-to-`SparseFactor` dictionary. When supplied,
  it enables name and rank checks. Values must have the tensor fields expected
  by this API; malformed values can cause field-access errors. Coordinate
  validity and coefficient arithmetic are not checked.

# Returns
A `Vector{String}` of diagnostics. An empty vector means the expression
passes the checks enabled by these arguments. Messages identify syntax
locations, unknown names, tensor ranks, powered-factor arity, repeated outputs
or index-rule violations as applicable.

# Guarantees
1. The result describes acceptance by the selected textual coordinate-network
   interface, not evaluation of a tensor or an algebraic identity.
   [proved/agent: by construction]
2. Without `tensors`, unknown names and stored tensor ranks are not diagnosed.
   With `tensors`, referenced names and ranks are checked.
   [proved/agent: by construction]
3. The expression, output-label collection, dictionary and tensor data are
   not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
An empty result does not guarantee successful evaluation: an empty coefficient
dictionary, invalid coordinates, missing matrix-power dimension, singular
negative power or unsupported coefficient arithmetic can still fail later.
This does not check variance, Hopf axioms or topological invariance. It is not
a general Julia parser or a mathematical identity checker.

# Sources
cf:[[Kup91#^kup91-sec-0018]] finite-dimensional tensor notation — the mathematical setting of the coordinate expressions, not the repository's textual grammar.
cf:lib/topology.jl/src/tensor_network.jl, validate_tensor_expression — the enabled syntax, name, rank and label checks.

# Examples
```jldoctest
julia> isempty(validate_tensor_expression("C{3,4,1}M{4,2,5}", [1,2,5,3]))
true

julia> tensors = hopf_tensors(sweedler_algebra());

julia> any(occursin("Mu", message) && occursin("rank", message) for message in validate_tensor_expression("Mu{1,2}", [1,2]; tensors=tensors))
true
```

See also [`evaluate_tensor_expression`](@ref), [`label_tensor`](@ref).

Role: condition
"""
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

"""
    evaluate_tensor_expression(expression::AbstractString, tensors::AbstractDict, outputs=Int[]; dimension::Union{Nothing,Integer}=nothing, index_rule::Symbol=:network) -> SparseFactor

Evaluate a named tensor expression as a sparse coordinate contraction.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

- `expression`, the name-to-`SparseFactor` dictionary `tensors`, and `outputs`
  must pass [`validate_tensor_expression`](@ref) with the selected
  `index_rule` (`:network` or `:matlab`). Invalid rules raise `ArgumentError`;
  validation failures raise `ArgumentError` beginning `invalid tensor
  expression` and listing syntax, name, rank or label diagnostics.
- Each referenced tensor must have a nonempty coefficient dictionary with
  coherent key lengths and compatible arithmetic. Empty dictionaries are
  not handled as zero tensors by this evaluator; key-access errors can result.
  Tensor semantics and coefficient validity remain the caller's responsibility.
- Integer matrix powers require `dimension`; omitting it raises
  `ArgumentError` stating that dimension is required. It must fit `Int` and
  the matrix coordinates. Conversion, indexing, singular-matrix and arithmetic
  errors propagate. Negative powers require an invertible matrix.

# Returns
A `SparseFactor` with free labels in `outputs` order and the corresponding
coefficient dictionary, as in [`contract_network`](@ref).

# Guarantees
1. With valid tensor data and supported exact arithmetic without overflow,
   the result is the coordinate value of the named tensor expression, using
   ordinary matrix powers for powered factors. [proved/agent: by construction]
2. Exact rational coefficient types are preserved where the requested
   arithmetic supports them; output is not unconditionally complex floating
   point. [proved/agent: by construction]
3. The expression, tensor dictionary, stored tensors and output labels are
   not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Validation does not certify tensor variance, Hopf identities or topological
invariance. Floating-point contractions and powers can discard coefficients
of magnitude at most `1e-12` and have no certified error bound. The syntax is
not an arbitrary Julia expression or a general symbolic algebra language.

# Sources
cf:[[Kup91#^kup91-sec-0018]] coordinate tensor contraction; the textual interface belongs to the repository's Algebra API.

# Examples
```jldoctest
julia> a = SparseFactor(Int[], Dict{Tuple,Int}((1,1)=>2, (2,2)=>3));

julia> b = SparseFactor(Int[], Dict{Tuple,Rational{Int}}((1,1)=>1//3, (2,1)=>2//5));

julia> r = evaluate_tensor_expression("A{1,2}B{2,3}", Dict("A"=>a,"B"=>b), [1,3]); (r.terms[(1,1)], r.terms[(2,1)])
(2//3, 6//5)
```

See also [`validate_tensor_expression`](@ref), [`dense_tensor`](@ref).

Role: convert
"""
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

"""
    verify_hopf(algebra::FiniteHopfAlgebra; atol=1e-10, degrees=nothing) -> Dict{String,Bool}

Return coordinate checks of the named Hopf identities for the supplied structural tensors.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/algebra.md).

- `algebra` must have coherent tensor ranks, valid coordinate indices and
  supported coefficient arithmetic accepted by [`hopf_tensors`](@ref) and
  [`evaluate_tensor_expression`](@ref). Malformed tensor data can raise
  exceptions rather than returning failed checks.
- `atol` is a finite nonnegative absolute coefficient tolerance. The caller
  supplies a meaningful tolerance; it is not validated independently.
- `degrees` is `nothing` for the ordinary bialgebra comparison, or an indexed
  collection of basis degrees for its graded comparison. A length mismatch
  raises `DimensionMismatch` stating that one degree per basis element is
  required. The caller supplies a valid homogeneous grading; parity arithmetic
  and coefficient-conversion errors propagate.

# Returns
A `Dict{String,Bool}` with checks named `left_unit`, `right_unit`,
`left_counit`, `right_counit`, `associativity`, `coassociativity`,
`left_antipode`, `right_antipode`, and `bialgebra`. Each value reports the
corresponding tensor comparison within `atol`.

# Guarantees
1. For coherent tensors, the named comparisons test the corresponding
   structural identities on the input coordinate space, using the graded
   bialgebra comparison when `degrees` is supplied.
   [proved/agent: by construction]
2. Checks use absolute coefficient differences; missing sparse coefficients
   are treated as zero. There is no relative tolerance parameter.
   [proved/agent: by construction]
3. The algebra, its tensor containers and the supplied degrees are not
   modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/Algebra.md).

# Not guaranteed
Passing all comparisons is not an exact proof of every Hopf axiom. The result
contains only the named checks; it does not separately certify unit/coproduct
or counit/product compatibility, the grading, distinguished integral
identities, crossing tensors, involutivity or unimodularity. Floating-point
zero filtering and numerical matrix operations have no certified error bound.
This is not a topological-invariance test.

# Sources
cf:[[Kup91#^kup91-sec-0038]] Hopf structural identities — the mathematical target; the graded comparison is an additional repository option.

# Examples
```jldoctest
julia> all(values(verify_hopf(sweedler_algebra())))
true

julia> all(values(verify_hopf(exterior_algebra(2); degrees=count_ones.(0:3))))
true
```

See also [`is_involutory`](@ref), [`is_unimodular`](@ref).

Role: condition
"""
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
