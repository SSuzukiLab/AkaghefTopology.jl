"""
    URDiagram(args...) -> URDiagram

Construct a mutable record for a weighted U(R) diagram.

# Arguments
The zero-argument constructor supplies an empty diagram. The four-field
constructor accepts a coefficient `C`, vertex word `V::Vector{Int}`, edge
table `E::Matrix{Int}`, and weights `W::Vector{Any}`, or values convertible
to those field types. It performs no consistency check; conversion errors
propagate. Use [`set_data_cvw!`](@ref) to build the edge table from a word and
[`validate_data`](@ref) to check the supported structural conditions.

# Returns
`URDiagram` with fields `C`, `V`, `E` and `W` in the
[URDiagrams API](../specs/ur-diagrams.md#definitions).

# Guarantees
1. The zero-argument result has coefficient one and empty vertex, edge and
   weight data, with independently allocated containers.
   [proved/agent: by construction]
2. The field constructor retains values after conversion; already correctly
   typed arrays are shared. Inputs are not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The constructor does not certify a valid operator diagram, a finite trace,
or a represented knot or manifold. Mutable aliases can invalidate its data.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m properties and constructor — source weighted-diagram container; no mathematical validity theorem is asserted.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> validate_data(d)
true
```

Role: construct
"""
mutable struct URDiagram
    C::Any
    V::Vector{Int}
    E::Matrix{Int}
    W::Vector{Any}
end

URDiagram() = URDiagram(1, Int[], zeros(Int, 5, 0), Any[])
LinearAlgebra.rank(d::URDiagram) = count(==(0), d.V) + 1
next_index(d::URDiagram) = isempty(d.V) ? 1 : maximum(abs.(d.V)) + 1

"""
    set_data_cvw!(d::URDiagram, C, V, W) -> URDiagram

Replace a U(R) diagram's coefficient, vertex word and weights and rebuild its edge table.

# Arguments
`V` must be collectable and convertible entrywise to `Int`. `W` is an array
of weights or one scalar weight; the weight count is checked against the
distinct positive labels of `V`. A mismatch raises `ArgumentError` with
`one weight is required for each positive edge`. Conversion errors propagate.
`C` is arbitrary. The caller supplies paired, consistent diagram data;
this function does not perform [`validate_data`](@ref).

# Returns
The same `URDiagram` with fields replaced by the supplied data and a new
edge table, in the [URDiagrams API](../specs/ur-diagrams.md#definitions).

# Guarantees
1. The coefficient is retained, the vertex sequence is converted to a new
   integer vector, and weights occupy a new `Vector{Any}`.
   [proved/agent: by construction]
2. On success the rebuilt table has one column for each distinct positive
   label; auxiliary table entries are reset. [proved/agent: by construction]
3. `d` is modified. Input containers are not modified, but coefficient and
   weight objects are not deep-copied. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Duplicate or unpaired signed endpoints are not rejected here. Construction
does not certify an operator interpretation, finite trace or represented
space. Rebuilding discards previous auxiliary edge-table metadata.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m setDataFromCVW — source data-setting operation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> validate_data(d)
true
```

Role: move
"""
function set_data_cvw!(d::URDiagram, C, V, W)
    vertices = Int.(collect(V))
    edges = sort(unique(filter(>(0), vertices)))
    weights = W isa AbstractArray ? Any[W...] : Any[W]
    length(weights) == length(edges) || throw(ArgumentError("one weight is required for each positive edge"))
    d.C, d.V, d.E, d.W = C, vertices, zeros(Int, 5, length(edges)), weights
    d.E[1, :] = edges
    d
end

"""
    set_data_cvs!(d::URDiagram, C, V, S; epsilon=nothing) -> URDiagram

Replace a U(R) diagram's data with weights perturbed by supplied or formal epsilon values.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`C`, `V` and `S` satisfy the input requirements of [`set_data_cvw!`](@ref),
with `S` supplying base weights. Omitted `epsilon` supplies formal symbolic
perturbations; a scalar is used for every weight, and an array must match
the weight count. A length mismatch raises `ArgumentError` naming the
epsilon/edge-count mismatch. Weight addition must be defined; conversion
and arithmetic errors propagate. No smallness or sign check is performed.

# Returns
The same `URDiagram`, with the base weights replaced by their sums with
the corresponding perturbations.

# Guarantees
1. On success the vertex and edge data are those of `set_data_cvw!` and each
   weight equals its base value plus its selected perturbation.
   [proved/agent: by construction]
2. The diagram is modified; input containers are not modified. An epsilon
   length or arithmetic failure can leave the base diagram already installed.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The perturbation is not a limit operation and need not remove singularities.
Symbolic epsilon names do not imply convergence or analytic regularization.
No knot or manifold invariant is certified by this constructor.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m setDataFromCVS — source perturbed-weight construction.

# Examples
```jldoctest
julia> code = [-1,2,-4,1,-3,4,-2,3]; d = URDiagram();

julia> set_data_cvs!(d, 1, trans_knot_complement(code), repeat([1,-1,-1,1], outer=4); epsilon=1e-5); simplify2!(d);

julia> isapprox(trace2(d), -4.526591341323682e-5; rtol=1e-12)
true
```

Role: move
"""
function set_data_cvs!(d::URDiagram, C, V, S; epsilon=nothing)
    set_data_cvw!(d, C, V, S)
    eps = epsilon === nothing ? Any[sym("eps$(e)") for e in d.E[1, :]] :
          epsilon isa AbstractArray ? Any[epsilon...] : fill(epsilon, length(d.W))
    length(eps) == length(d.W) || throw(ArgumentError("epsilon length does not match edge count"))
    d.W = Any[d.W[i] + eps[i] for i in eachindex(d.W)]
    d
end

"""
    set_data_v!(d::URDiagram, V) -> URDiagram

Set a U(R) diagram's vertex word with unspecified numeric weights.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`V` must support positive-label filtering and collection into integer data.
Conversion/iteration errors propagate. Pairing and structural consistency
are the caller's responsibility and can be checked by [`validate_data`](@ref).

# Returns
The same `URDiagram` with coefficient one, the supplied vertex word, a
rebuilt edge table and `NaN` weights.

# Guarantees
1. Each distinct positive edge label receives one unspecified `NaN` weight.
   [proved/agent: by construction]
2. The diagram's fields are replaced; the supplied word is not modified
   or retained as its mutable backing vector. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The result is structural data, not a finite scalar-weight evaluation.
Subsequent arithmetic can propagate `NaN`; valid pairing is not checked here.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m setDataFromV — source word-only construction.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_v!(d, [1,2,-3,3,4,-1,-2,0,-4]); reduct3!(d);

julia> Tuple(d.V)
(1, -3, 3, 4, -1, 0, -4)
```

Role: move
"""
set_data_v!(d::URDiagram, V) = set_data_cvw!(d, 1, V, fill(NaN, length(unique(filter(>(0), V)))))

"""
    set_data_cm!(d::URDiagram, C, matrix::AbstractMatrix) -> URDiagram

Replace a U(R) diagram by the weighted endpoint data of a standard coefficient matrix.

# Arguments
Requires U11 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Squareness is checked; `ArgumentError` states that the coefficient matrix
must be square. Conversion/allocation errors propagate. Edge labels and
dimensions must fit `Int`; no operator-domain check is performed.

# Returns
The same `d`, holding the [standard matrix data](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. The output is the standard endpoint representation of the supplied matrix, with coefficient `C`. [proved/agent: by construction]
2. The diagram fields are replaced. The matrix is not modified; its entries are retained by value or reference according to their own types. [proved/agent: by construction]

# Not guaranteed
No finite trace, trace-class operator interpretation or equivalence with
another matrix convention is certified.

# Sources
cf:lib/topology.jl/docs/specs/ur-diagrams.md Standard matrix data — coefficient matrix endpoint representation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cm!(d, 1, [2.0 1.0; 3.0 4.0]);

julia> Tuple(d.V)
(1, 3, -1, -2, 0, 2, 4, -3, -4)
```

Role: move
"""
function set_data_cm!(d::URDiagram,C,matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    positive=[Int[] for _ in 1:rows]
    negative=[Int[] for _ in 1:rows]
    edge=0
    for column in 1:columns, row in 1:rows
        edge += 1
        push!(positive[row],edge)
        push!(negative[column],-edge)
    end
    vertices=Int[]
    for section in 1:rows
        append!(vertices,positive[section]); append!(vertices,negative[section])
        section < rows && push!(vertices,0)
    end
    set_data_cvw!(d,C,vertices,vec(matrix))
end

"""
    standard_matrix_adjustment(matrix::AbstractMatrix) -> Matrix

Return the cyclic correction matrix for the legacy standard-matrix trace.

# Arguments
Requires U11 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Squareness is checked with `ArgumentError` stating that the coefficient
matrix must be square. The first entry must support the required zero/one
and subtraction operations when nonempty; arithmetic/type errors propagate.

# Returns
A fresh square matrix with the input dimensions and the
[legacy trace correction](../conventions/URDiagrams.md#legacy-trace-formulas).

# Guarantees
1. The output has the specified cyclic correction entries, including the empty and one-dimensional cases. [proved/agent: by construction]
2. The input matrix is unchanged. [proved/agent: by construction]

# Not guaranteed
This correction alone is not a trace value or a certificate that a matrix
comes from an operator or an o-graph.

# Sources
cf:lib/topology.jl/docs/conventions/URDiagrams.md Legacy trace formulas — cyclic correction and small-size cases.

# Examples
```jldoctest
julia> Tuple(vec(standard_matrix_adjustment(fill(1//1,2,2))))
(1//1, -1//1, -1//1, 1//1)
```

Role: construct
"""
function standard_matrix_adjustment(matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    rows==0 && return zeros(eltype(matrix),0,0)
    prototype=matrix[1]
    adjustment=fill(zero(prototype),rows,rows)
    rows==1 && return adjustment
    for index in 1:rows
        adjustment[index,index]=one(prototype)
        adjustment[index==1 ? rows : index-1,index]-=one(prototype)
    end
    adjustment
end

"""
    standard_matrix_trace(C, matrix::AbstractMatrix) -> Any

Evaluate the legacy standard-matrix trace expression at the supplied coefficients.

# Arguments
Requires U3 and U11 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Squareness is checked with `ArgumentError`. Coefficient arithmetic and a
nonsingular denominator are caller obligations; arithmetic errors or
nonfinite values follow the coefficient types, without a regularization step.

# Returns
The scalar expression in the [legacy trace convention](../conventions/URDiagrams.md#legacy-trace-formulas),
with type determined by the coefficient arithmetic.

# Guarantees
1. On its finite arithmetic domain the value is the convention's standard-matrix expression. [proved/agent: by construction]
2. Arguments are not modified. [proved/agent: by construction]

# Not guaranteed
No analytic trace-classness, regulator limit, certified floating error bound
or general manifold-invariant interpretation is established.

# Sources
cf:lib/topology.jl/docs/conventions/URDiagrams.md Legacy trace formulas — algebraic determinant expression.

# Examples
```jldoctest
julia> isapprox(standard_matrix_trace(1, [2.0 1.0; 3.0 4.0]), -1/15)
true
```

Role: invariant
"""
function standard_matrix_trace(C,matrix::AbstractMatrix)
    rows,columns=size(matrix)
    rows == columns || throw(ArgumentError("the coefficient matrix must be square"))
    -C/_det(matrix+standard_matrix_adjustment(matrix))
end

function edge_index(d::URDiagram, e::Integer)
    i = findfirst(==(abs(e)), vec(d.E[1, :]))
    i === nothing ? 0 : i
end
vertex_index(d::URDiagram, v::Integer) = something(findfirst(==(v), d.V), 0)
get_weight(d::URDiagram, e::Integer) = d.W[edge_index(d, e)]

"""
    validate_data(d::URDiagram) -> Bool

Check the supported structural consistency conditions of a weighted U(R) diagram.

# Arguments
`d` is a `URDiagram`. Conditions U1–U2 of the
[specification](../specs/ur-diagrams.md#numbered-conditions) are enforced. A failed check raises
`ArgumentError` naming duplicate labels, unpaired vertices, table mismatch
or weight-count mismatch; an unusable edge-table shape can raise `BoundsError`.
These checks belong to [URDiagrams](../specs/ur-diagrams.md#definitions).

# Returns
`true` on success; invalid data raise an exception rather than returning
`false`.

# Guarantees
1. Success means U1–U2 hold for the current diagram.
   [proved/agent: by construction]
2. The diagram is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Coefficient and weight arithmetic, auxiliary edge-table rows, admissibility
of a subsequent rewrite, operator interpretation and finite trace are not
checked. This is not a validity test for an associated knot or manifold.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m properties — paired endpoint and weight data; the Julia validator checks its supported structural subset.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> validate_data(d)
true
```

Role: condition
"""
function validate_data(d::URDiagram)
    positives = sort(filter(>(0), d.V))
    negatives = sort(-filter(<(0), d.V))
    positives == unique(positives) || throw(ArgumentError("positive edge labels must be unique"))
    negatives == unique(negatives) || throw(ArgumentError("negative edge labels must be unique"))
    mismatch = union(setdiff(Set(positives), Set(negatives)), setdiff(Set(negatives), Set(positives)))
    positives == negatives || throw(ArgumentError("vertices are not paired: $mismatch"))
    positives == vec(d.E[1, :]) || throw(ArgumentError("edge table does not match vertices"))
    length(d.W) == length(positives) || throw(ArgumentError("weight count does not match edges"))
    true
end

"""
    check_data(d::URDiagram, C, V, W) -> Bool

Check a U(R) diagram's coefficient, vertex word and weights against expected data.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`C` is an expected coefficient, `V` is collectable expected vertex data,
and `W` is iterable expected weight data. Equality must yield Boolean
decisions. A mismatch raises `AssertionError` naming `C`, `V`, or `W`;
collection, iteration and equality errors propagate.

# Returns
`true` when all three comparisons succeed, rather than `false` on mismatch.

# Guarantees
1. Success means equality of those three fields with the supplied expected
   data, using their Julia equality operations. [proved/agent: by construction]
2. The diagram and arguments are not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The edge table is not compared. This is not [`validate_data`](@ref), a
symbolic identity prover, or a test of equivalence of represented spaces.

# Sources
cf:lib/topology.jl/src/ur_diagram.jl, check_data — structural data checking and diagnostic display.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]); swap!(d, [3,-3]);

julia> check_data(d, 1/(1-sym(:z)), [1,-2,3,-3,-1,2], [sym(:x),sym(:y),sym(:z)/(1-sym(:z))])
true
```

Role: condition
"""
function check_data(d::URDiagram, C, V, W)
    d.C == C || throw(AssertionError("C differs"))
    d.V == collect(V) || throw(AssertionError("V differs"))
    d.W == Any[W...] || throw(AssertionError("W differs"))
    true
end

"""
    swap!(d::URDiagram, pair; skip=false) -> URDiagram

Apply the selected endpoint interchange to a weighted U(R) diagram.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must satisfy caller-supplied U1–U3 of the
[specification](../specs/ur-diagrams.md#numbered-conditions). `pair` supplies two
integer-convertible signed endpoints. Unless `skip=true`, absent endpoints
or failure of U4 raises `ArgumentError` with `vertices must exist` or
`vertices must be adjacent`. Overall structural validity and nonsingular
arithmetic are the caller's responsibility. Conversion/indexing/arithmetic
failures propagate; skipped placement checks provide no invalid-input guarantee.

# Returns
The same `URDiagram` with the local endpoint order and associated weights
updated.

# Guarantees
1. On valid nonsingular input, the selected endpoint order is exchanged.
   Same-sign swaps preserve coefficient and weights; a same-edge swap can
   rescale its weight and coefficient, and an opposite-sign different-edge
   swap introduces one additional edge. [proved/agent: by construction]
2. The diagram is modified; the pair argument is not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
This does not certify operator equality, trace preservation or a move on a
represented space. Singular arithmetic can yield infinities or exceptions;
the endpoint order may already have changed when arithmetic fails.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m swap — source endpoint interchange and coefficient transformations.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]); swap!(d, [3,-3]);

julia> (Tuple(d.V), d.C == 1/(1-sym(:z)), d.W[3] == sym(:z)/(1-sym(:z)))
((1, -2, 3, -3, -1, 2), true, true)
```

Role: move
"""
function swap!(d::URDiagram, pair; skip=false)
    a, b = Int.(pair)
    p1, p2 = vertex_index(d, a), vertex_index(d, b)
    !skip && (p1 == 0 || p2 == 0) && throw(ArgumentError("vertices must exist"))
    !skip && abs(p1 - p2) != 1 && throw(ArgumentError("vertices must be adjacent"))
    s1, s2 = sign(d.V[p1]), sign(d.V[p2])
    d.V[p1], d.V[p2] = d.V[p2], d.V[p1]
    if a == -b
        wi = edge_index(d, a); w = d.W[wi]
        if !xor(p1 < p2, s1 < 0)
            d.W[wi] = w / (1 - w); d.C = d.C / (1 - w)
        else
            d.W[wi] = w / (1 + w); d.C = d.C / (1 + w)
        end
    elseif s1 != s2
        w1, w2 = get_weight(d, a), get_weight(d, b)
        q1, q2 = vertex_index(d, -a), vertex_index(d, -b)
        nw = xor(s1 > 0, p2 < p1) ? -w1 * w2 : w1 * w2
        n = next_index(d)
        lo, hi = min(q1, q2), max(q1, q2)
        vlo = q1 < q2 ? -s1*n : -s2*n
        vhi = q1 < q2 ? -s2*n : -s1*n
        insert!(d.V, lo, vlo); insert!(d.V, hi + 1, vhi)
        d.E = hcat(d.E, zeros(Int, 5)); d.E[1, end] = n; push!(d.W, nw)
    end
    d
end

"""
    delete_edge!(d::URDiagram, edge::Integer) -> URDiagram

Remove one edge and its two endpoint labels from a weighted U(R) diagram.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have consistent vertex, edge and weight data as checked by
[`validate_data`](@ref); this operation does not run that check. `edge`
selects an edge by absolute label. A missing edge raises `ArgumentError`
with `edge does not exist`. Malformed data may cause indexing failures or
an inconsistent result, and are the caller's responsibility.

# Returns
The same mutated `URDiagram`.

# Guarantees
1. On consistent input, the selected edge's table entry, weight and signed
   endpoints are absent from the result. Remaining entries retain their
   order and coefficient `C` is unchanged. [proved/agent: by construction]
2. The diagram is modified in place; shared mutable containers may expose
   part of that modification to aliases. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Deletion alone is not certified to preserve an operator, trace, knot or
manifold. Failure is not transactional.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m del — source edge deletion operation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]);

julia> dilation!(d, 3, -2); add_edges!(d, (2,3)); dilation!(d, 2, -1); delete_edge!(d, 1);

julia> Tuple(d.V)
(-2, 2)
```

Role: move
"""
function delete_edge!(d::URDiagram, edge::Integer)
    i = edge_index(d, edge); i == 0 && throw(ArgumentError("edge does not exist"))
    d.V = filter(v -> abs(v) != abs(edge), d.V)
    d.E = d.E[:, setdiff(axes(d.E, 2), i)]; deleteat!(d.W, i); d
end

"""
    add_edges!(d::URDiagram, pair; skip=false) -> URDiagram

Merge a selected pair of U(R) diagram edges by adding their weights.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must satisfy U1–U3; those checks are not performed here.
`pair` supplies two distinct integer-convertible
edge labels, interpreted by absolute value. Unless `skip=true`, the endpoint
placement check U5 runs; failure raises `ArgumentError`
with `both signed pairs must be adjacent`. `skip` bypasses these checks and
makes valid placement the caller's responsibility. Weight addition must be
defined; conversion, indexing and arithmetic errors propagate.

# Returns
The same mutated `URDiagram`.

# Guarantees
1. On consistent input passing the placement checks, the first selected
   edge remains with the sum of the two weights and the second edge is
   removed. Other weights and `C` are unchanged. [proved/agent: by construction]
2. The diagram is modified; the pair argument is not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
No general theorem identifying this rewrite with an equality of operators
or preserving a represented space is certified. Invalid unchecked placement
can corrupt data, and errors do not roll back earlier mutations.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m add — source weighted-edge merge, without a general trace-preservation proof.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]);

julia> dilation!(d, 3, -2); add_edges!(d, (2,3)); dilation!(d, 2, -1); delete_edge!(d, 1);

julia> Tuple(d.V)
(-2, 2)
```

Role: move
"""
function add_edges!(d::URDiagram, pair; skip=false)
    a, b = abs.(Int.(pair))
    pp = vertex_index.(Ref(d), (a, b)); pn = vertex_index.(Ref(d), (-a, -b))
    !skip && (any(==(0), pp) || any(==(0), pn) || abs(pp[1]-pp[2]) != 1 || abs(pn[1]-pn[2]) != 1) &&
        throw(ArgumentError("both signed pairs must be adjacent"))
    for p in sort([pp[2], pn[2]], rev=true); deleteat!(d.V, p); end
    ia, ib = edge_index(d, a), edge_index(d, b)
    d.W[ia] = d.W[ia] + d.W[ib]
    d.E = d.E[:, setdiff(axes(d.E, 2), ib)]; deleteat!(d.W, ib); d
end

"""
    dilation!(d::URDiagram, dil::Integer, arg::Integer; skip=false) -> URDiagram

Move a selected signed endpoint across a U(R) dilation pair and update its edge weight.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must satisfy U1–U3. `dil` selects an edge
by absolute label; `arg` selects the signed endpoint of a different edge.
With `skip=false`, failure of the dilation placement condition U6 raises
`ArgumentError` with `dilation vertices and target must be adjacent`.
The caller supplies valid data and nonsingular arithmetic; these are not
fully checked. `skip=true` additionally bypasses the placement check.
Indexing, conversion and arithmetic failures propagate.

# Returns
The same mutated `URDiagram` after the local dilation rewrite.

# Guarantees
1. On valid nonsingular input, the endpoint multiset, edge labels and `C`
   are unchanged; only the selected local order and the target edge's weight
   can change. [proved/agent: by construction]
2. The diagram is modified in place. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
No analytic operator identity, trace preservation, or equivalence of
represented spaces is certified. Singular arithmetic can produce infinities
or exceptions. Failure need not leave the diagram unchanged.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m dilation — source local rewrite; operator-domain hypotheses are not established here.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]);

julia> dilation!(d, 3, -2); add_edges!(d, (2,3)); dilation!(d, 2, -1); delete_edge!(d, 1);

julia> trace3(d) == -1
true
```

Role: move
"""
function dilation!(d::URDiagram, dil::Integer, arg::Integer; skip=false)
    dil = abs(dil); pos = [vertex_index(d, dil), vertex_index(d, -dil), vertex_index(d, arg)]
    !skip && (any(==(0), pos) || sort(pos) != collect(minimum(pos):maximum(pos))) &&
        throw(ArgumentError("dilation vertices and target must be adjacent"))
    e1, e2, e3 = pos[2]-pos[1], sign(arg), sign(pos[3]-pos[1])
    w = get_weight(d, dil); mul = e1 < 0 ? w/(1-w) : w+1
    exponent=e2*e3
    powered=exponent<0 && mul isa Integer ? 1/(mul^(-exponent)) : mul^exponent
    d.W[edge_index(d, arg)] = get_weight(d, arg) * powered
    old = copy(d.V[pos]); shift = e3 .* [1, 1, -2]
    d.V[pos .+ shift] = old; d
end

"""
    compose_dilation!(d::URDiagram, pair; skip=false) -> URDiagram

Compose a selected pair of adjacent U(R) dilation edges into one edge.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have consistent data as checked by [`validate_data`](@ref), and
`pair` supplies two distinct integer-convertible edge labels. Unless `skip`,
failed endpoint existence or block adjacency raises `ArgumentError` with
`dilations must form one adjacent block`; an additional required endpoint
swap can raise the placement errors of [`swap!`](@ref). The caller supplies
admissible dilation data and nonsingular scalar arithmetic. `skip=true`
bypasses placement checks, with no correctness promise for invalid input.

# Returns
The same mutated `URDiagram` with the first selected label retained.

# Guarantees
1. On admissible nonsingular input, the second selected edge is removed
   and the retained edge receives the composed dilation weight. No new
   edge label is introduced. [proved/agent: by construction]
2. The diagram is modified, including a possible coefficient change when
   a required swap changes its normalization; `pair` is not modified.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The call does not certify an operator identity or trace preservation for
arbitrary weights. Arithmetic errors and failed intermediate placement
checks can leave a partially modified diagram.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m composeDilation — source local composition operation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [-5,1,-4,4,3,2,-3,5,-2,-1], [-1,-1,-1,1,1]);

julia> swap!(d,(2,3)); dilation!(d,3,5); swap!(d,(2,5)); add_edges!(d,(2,3)); swap!(d,(-5,1)); swap!(d,(5,6)); add_edges!(d,(4,5)); add_edges!(d,(2,6)); dilation!(d,4,1); add_edges!(d,(1,2)); swap!(d,(-1,1)); compose_dilation!(d,(1,4));

julia> Tuple(d.V)
(-1, 1)
```

Role: move
"""
function compose_dilation!(d::URDiagram, pair; skip=false)
    a, b = abs.(Int.(pair)); pos = collect(vertex_index.(Ref(d), (a, b, -a, -b)))
    !skip && (any(==(0), pos) || sort(pos) != collect(minimum(pos):maximum(pos))) &&
        throw(ArgumentError("dilations must form one adjacent block"))
    # MATLAB `composeDilation`: first normalize the two possible crossed
    # orderings, then choose the weight law from the endpoint directions.
    direction=(pos[1]>pos[3],pos[2]>pos[4])
    if direction[1] && !direction[2]
        swap!(d,(d.V[pos[1]],d.V[pos[3]]);skip=skip)
    elseif direction[2] && !direction[1]
        swap!(d,(d.V[pos[2]],d.V[pos[4]]);skip=skip)
    end
    ia, ib = edge_index(d, a), edge_index(d, b); w1, w2 = d.W[ia], d.W[ib]
    d.W[ia] = all(direction) ? 1-(1-w1)*(1-w2) : (1+w1)*(1+w2)-1
    for p in sort([pos[2], pos[4]], rev=true); deleteat!(d.V, p); end
    d.E = d.E[:, setdiff(axes(d.E, 2), ib)]; deleteat!(d.W, ib); d
end

"""
    put_zero!(d::URDiagram, positive_position::Integer, negative_position::Integer) -> URDiagram

Insert the paired endpoints of a zero-weight dilation edge into a U(R) diagram.

# Arguments
Requires U10 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Both gap positions are checked against the original word bounds; an invalid
position raises `BoundsError`. Structural consistency, label arithmetic and
available storage are caller obligations. Allocation and arithmetic failures
can leave partial mutation.

# Returns
The same `URDiagram`, with one additional edge in its table and weight vector.

# Guarantees
1. The new edge has weight zero and the endpoint positions of the [insertion convention](../conventions/URDiagrams.md#convention-dictionary). Existing endpoints keep their relative order. [proved/agent: by construction]
2. Coefficient and number of word separators are unchanged; the diagram is mutated in place. [proved/agent: by construction]

# Not guaranteed
No zero separator is inserted. General operator or trace preservation is
unverified, particularly across the unequal-gap source-convention difference.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m put0 — zero-weight edge; the unequal-gap sign difference is recorded in the convention dictionary.

# Examples
```jldoctest
julia> t = laurent_variable(:t); d = URDiagram(); set_data_cvw!(d,1,[1,2,-3,3,-1,-2],[t,t^2,t^3]);

julia> put_zero!(d,length(d.V),length(d.V)); reduct4!(d,4);

julia> Tuple(d.V)
(1, 2, -1, -2, 4, -4)
```

Role: move
"""
function put_zero!(d::URDiagram, positive_position::Integer, negative_position::Integer)
    0 <= positive_position <= length(d.V) || throw(BoundsError(d.V,positive_position))
    0 <= negative_position <= length(d.V) || throw(BoundsError(d.V,negative_position))
    new_edge=next_index(d)
    if positive_position <= negative_position
        insert!(d.V,positive_position+1,new_edge)
        insert!(d.V,negative_position+2,-new_edge)
    else
        insert!(d.V,negative_position+1,-new_edge)
        insert!(d.V,positive_position+2,new_edge)
    end
    d.E=hcat(d.E,zeros(Int,5)); d.E[1,end]=new_edge; push!(d.W,0)
    d
end

function _reduct_add!(d::URDiagram)
    changed = false
    for i in 1:max(length(d.V)-1, 0)
        a, b = d.V[i], d.V[i+1]
        if a > 0 && b > 0
            na, nb = vertex_index(d, -a), vertex_index(d, -b)
            if abs(na-nb) == 1
                add_edges!(d, (a, b)); changed = true; break
            end
        end
    end
    changed
end

function _groups(d::URDiagram)
    boundaries=vcat(0,findall(!=(0),diff(sign.(d.V))))
    group(position)=count(boundary -> position > boundary,boundaries)
    ([group(vertex_index(d,edge)) for edge in vec(d.E[1,:])],
     [group(vertex_index(d,-edge)) for edge in vec(d.E[1,:])])
end

"""
    reduct3!(d::URDiagram) -> URDiagram

Apply the same-sign-group edge-merge relation to a weighted U(R) diagram.

# Arguments
Requires U1–U3 and U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
No structural validation occurs at entry. The caller supplies coherent
weights and labels supporting the [local reduction relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).
Invalid intermediate data or unsupported arithmetic can raise indexing,
conversion or arithmetic errors after partial mutation.

# Returns
The same `URDiagram` after the supported merges.

# Guarantees
1. On successful return, the result is the weighted endpoint diagram obtained by the specified same-group merges; every performed merge adds the two edge weights. [proved/agent: by construction]
2. The diagram is mutated in place; its coefficient is unchanged. [proved/agent: by construction]

# Not guaranteed
No complete normal form or general trace/operator preservation theorem is
established.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct3 — local same-group relation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_v!(d,[1,2,-3,3,4,-1,-2,0,-4]); reduct3!(d);

julia> Tuple(d.V)
(1, -3, 3, 4, -1, 0, -4)
```

Role: move
"""
function reduct3!(d::URDiagram)
    positive_groups,negative_groups=_groups(d)
    deleted=Set{Int}()
    pairs=Tuple{Int,Int}[]
    edges=vec(d.E[1,:])
    for left in eachindex(edges), right in left+1:length(edges)
        right in deleted && continue
        if positive_groups[left] == positive_groups[right] &&
           negative_groups[left] == negative_groups[right]
            push!(pairs,(edges[left],edges[right])); push!(deleted,right)
        end
    end
    for pair in pairs
        add_edges!(d,pair;skip=true)
    end
    d
end

"""
    reduct4!(d::URDiagram, dilation_edges) -> URDiagram

Apply the adjacent-group chord relation relative to supplied U(R) section dilations.

# Arguments
Requires U12 and the referenced U1–U6 conditions
([specification](../specs/ur-diagrams.md#numbered-conditions)).
`dilation_edges` is an iterable of integer-convertible edge labels in section
order; a single integer selects one dilation. No global preflight is made.
A missing section dilation raises `ArgumentError` naming its rank block;
local swap/dilation checks raise `ArgumentError` about absent or nonadjacent
vertices. Conversion, indexing and coefficient errors propagate.

# Returns
The same diagram after the [local relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. Successful return gives the weighted endpoint diagram of the specified adjacent-group relation with the supplied dilations. [proved/agent: by construction]
2. The diagram is modified; the supplied iterable is not modified. Failure may retain partial diagram changes. [proved/agent: by construction]

# Not guaranteed
The operation does not verify operator domains, finite traces or general
trace preservation, and does not certify a canonical reduction.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct4 — adjacent-group dilation relation.

# Examples
```jldoctest
julia> t = laurent_variable(:t); d = URDiagram(); set_data_cvw!(d,1,[1,2,-3,3,-1,-2],[t,t^2,t^3]);

julia> put_zero!(d,length(d.V),length(d.V)); reduct4!(d,[4]);

julia> Tuple(d.V)
(1, 2, -1, -2, 4, -4)
```

Role: move
"""
function reduct4!(d::URDiagram, dilation_edges)
    dilations=Int.(collect(dilation_edges))
    positive_groups,negative_groups=_groups(d)
    edge_labels=vec(d.E[1,:])
    dilation_indices=Set(edge_index.(Ref(d),abs.(dilations)))
    # MATLAB `reduct4` computes both the candidates and their directions once,
    # before any swap.  Reclassifying later edges after an earlier swap skips
    # source-selected reductions and changes the resulting trace.
    candidates=[index for index in eachindex(edge_labels)
                if index ∉ dilation_indices && abs(positive_groups[index]-negative_groups[index])==1]
    composed=Tuple{Int,Int}[]
    for index in candidates
        edge=edge_labels[index]
        direction=positive_groups[index]-negative_groups[index]
        positive,negative=vertex_index(d,edge),vertex_index(d,-edge)
        between=direction==1 ? d.V[negative+1:positive-1] : d.V[positive+1:negative-1]
        positive>negative && (between=reverse(between))
        for other in filter(>(0),between)
            swap!(d,(edge,other))
        end
        for other in reverse(filter(<(0),between))
            swap!(d,(-edge,other))
        end
        direction==1 && swap!(d,(edge,-edge))
        section=count(==(0),d.V[1:vertex_index(d,edge)])+1
        section<=length(dilations) || throw(ArgumentError("missing separator dilation for rank block $section"))
        dilation=dilations[section]
        negative_position,separator_position=vertex_index(d,-edge),vertex_index(d,dilation)
        for other in d.V[negative_position+1:separator_position-1]
            dilation!(d,edge,other)
        end
        push!(composed,(dilation,edge))
    end
    for pair in reverse(composed)
        compose_dilation!(d,pair)
    end
    d
end

"""
    reduct4!(d::URDiagram, dilation_edge::Integer) -> URDiagram

Apply the adjacent-group chord relation relative to supplied U(R) section dilations.

# Arguments
Requires U12 and the referenced U1–U6 conditions
([specification](../specs/ur-diagrams.md#numbered-conditions)).
`dilation_edges` is an iterable of integer-convertible edge labels in section
order; a single integer selects one dilation. No global preflight is made.
A missing section dilation raises `ArgumentError` naming its rank block;
local swap/dilation checks raise `ArgumentError` about absent or nonadjacent
vertices. Conversion, indexing and coefficient errors propagate.

# Returns
The same diagram after the [local relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. Successful return gives the weighted endpoint diagram of the specified adjacent-group relation with the supplied dilations. [proved/agent: by construction]
2. The diagram is modified; the supplied iterable is not modified. Failure may retain partial diagram changes. [proved/agent: by construction]

# Not guaranteed
The operation does not verify operator domains, finite traces or general
trace preservation, and does not certify a canonical reduction.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct4 — adjacent-group dilation relation.

# Examples
```jldoctest
julia> t = laurent_variable(:t); d = URDiagram(); set_data_cvw!(d,1,[1,2,-3,3,-1,-2],[t,t^2,t^3]);

julia> put_zero!(d,length(d.V),length(d.V)); reduct4!(d,4);

julia> Tuple(d.V)
(1, 2, -1, -2, 4, -4)
```

Role: move
"""
reduct4!(d::URDiagram,dilation_edge::Integer)=reduct4!(d,[dilation_edge])

"""
    reduct5!(d::URDiagram) -> URDiagram

Apply the rank-one zero-weight-dilation relation to a U(R) diagram.

# Arguments
Requires U7 and U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
Rank and consistency are caller obligations, not checked at entry. Local
placement checks may raise `ArgumentError`; conversion, indexing, allocation
and coefficient errors propagate after any earlier mutation.

# Returns
The same weighted endpoint diagram after the
[rank-one relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. A successful result is the rank-one dilation reduction in the linked specification. [proved/agent: by construction]
2. The diagram is mutated in place. [proved/agent: by construction]

# Not guaranteed
One call need not finish a reduction or decrease the word length. General
operator and trace preservation remains open in the specification.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct5 — rank-one dilation relation.

# Examples
```jldoctest
julia> t = laurent_variable(:t); d = URDiagram(); set_data_cvw!(d,1,[1,2,-3,3,-1,-2],[t,t^2,t^3]); reduct5!(d);

julia> Tuple(d.V)
(1, 2, -1, -2, 4, -4)

julia> reduct5!(d); Tuple(d.V)
(5, -5)
```

Role: move
"""
function reduct5!(d::URDiagram)
    put_zero!(d,length(d.V),length(d.V))
    reduct4!(d,d.V[end-1])
end

"""
    reduct6!(d::URDiagram) -> URDiagram

Apply the first sign-inversion exchange to a weighted U(R) endpoint word.

# Arguments
Requires U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
The caller supplies coherent data and admissible local arithmetic. There
is no global validity check; the delegated swap's `ArgumentError` for
invalid placement and arithmetic/indexing errors propagate.

# Returns
The same diagram after the [first-inversion relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records),
or unchanged when no negative entry precedes a positive entry.

# Guarantees
1. At most one local exchange is performed. A word with no such inversion is unchanged. [proved/agent: by construction]
2. Any change is in place and a later failure does not undo it. [proved/agent: by construction]

# Not guaranteed
The relation is not a complete sort or reduction, and general operator or
trace preservation remains open.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct6 — first-inversion exchange.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d,1,[1,-2,3,-1,2,-3],[-1,-1,-1]); simplify1!(d); reduct6!(d);

julia> Tuple(d.V)
(4, -4)
```

Role: move
"""
function reduct6!(d::URDiagram)
    seen_negative=false
    for index in eachindex(d.V)
        if d.V[index] < 0
            seen_negative=true
        elseif d.V[index] > 0 && seen_negative
            swap!(d,d.V[[index,index-1]])
            break
        end
    end
    d
end

"""
    reduct7!(d::URDiagram) -> Tuple{URDiagram,Vector{Bool}}

Apply a bounded sign-inversion pass and return its interval flags.

# Arguments
Requires U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
No global check is made; local placement `ArgumentError` and indexing or
arithmetic errors propagate, potentially after partial mutation.

# Returns
`(d, complete)`, with the same diagram and the Boolean visited-interval
record defined in [local reduction records](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. Each returned flag records whether an exchange was absent in its visited interval, with the specified bounded-pass meaning. [proved/agent: by construction]
2. The diagram is modified in place; the flag vector is newly allocated. [proved/agent: by construction]

# Not guaranteed
The flags are not a certificate that all sections of the final word are
reduced. General trace preservation and operator equality remain open.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct7 — bounded sign-inversion pass.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d,1,[1,-2,3,-1,2,-3],[-1,-1,-1]); simplify1!(d); reduct7!(d);

julia> Tuple(d.V)
(4, -4)
```

Role: move
"""
function reduct7!(d::URDiagram)
    seen_negative=false; changed=false; complete=Bool[]
    # MATLAB uses `for idx=1:length(obj.V)`: the upper bound is captured
    # before a swap inserts its auxiliary edge pair.  Traversing a growing
    # Julia vector in the same pass performs extra inversions and makes the
    # multi-rank reduction diverge from the source implementation.
    initial_length=length(d.V)
    for index in 1:initial_length
        vertex=d.V[index]
        if vertex < 0
            seen_negative=true
        elseif vertex > 0 && seen_negative && !changed
            swap!(d,d.V[[index,index-1]])
            changed=true
        elseif vertex == 0
            push!(complete,!changed); seen_negative=false; changed=false
        end
    end
    push!(complete,!changed)
    d,complete
end

"""
    reduct8!(d::URDiagram) -> URDiagram

Apply the sectionwise zero-weight-dilation relation to a U(R) diagram.

# Arguments
Requires U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)).
No global validity or operator check occurs. Local placement and missing
dilation checks raise `ArgumentError`; indexing, allocation, conversion
and coefficient errors propagate, possibly after partial mutation.

# Returns
The same weighted diagram after the
[sectionwise relation](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. On success the output is the weighted endpoint diagram of that sectionwise dilation relation. [proved/agent: by construction]
2. The diagram is mutated in place. [proved/agent: by construction]

# Not guaranteed
No complete reduction, trace preservation or general operator identity is
certified. Word length need not decrease.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m reduct8 — sectionwise dilation relation.

# Examples
```jldoctest
julia> t = laurent_variable(:t); d = URDiagram(); set_data_cvw!(d,1,[1,2,-3,3,-1,-2],[t,t^2,t^3]); reduct5!(d); reduct8!(d);

julia> Tuple(d.V)
(5, -5)
```

Role: move
"""
function reduct8!(d::URDiagram)
    positions=vcat(findall(==(0),d.V).-1,length(d.V))
    for index in eachindex(positions)
        positions[index]+=2*(index-1)
        put_zero!(d,positions[index],positions[index])
    end
    dilations=Int[d.V[position+1] for position in positions]
    reduct4!(d,dilations)
    previous=vcat(-2,vertex_index.(Ref(d),dilations))
    snapshot=copy(d.V)
    for index in eachindex(dilations)
        start=previous[index]+3
        stop=previous[index+1]-1
        start>stop && continue
        for vertex in reverse(snapshot[start:stop])
            vertex>0 && break
            dilation!(d,dilations[index],vertex)
        end
    end
    d
end

"""
    simplify1!(d::URDiagram; maxiter=50) -> URDiagram

Apply the rank-one U(R) diagram reduction routine with an iteration bound.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

The caller supplies rank-one, structurally consistent data and scalar
arithmetic suitable for the reductions. Neither rank nor overall validity
is checked at entry; failed local checks or arithmetic propagate their
exceptions. `maxiter` should be a nonnegative integer; it is not validated.
For other values no iteration-bound interpretation is promised.

# Returns
The same `URDiagram`, possibly only partially reduced.

# Guarantees
1. For a nonnegative integer bound, at most `maxiter` reduction rounds are
   applied; a word already of length at most two is unchanged.
   [proved/agent: by construction]
2. The diagram can be modified in place. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Successful return does not certify complete reduction, a canonical diagram,
operator equality, trace preservation, or an invariant of a represented
space. Errors do not roll back reductions already applied.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m simplify1 — source rank-one reduction routine; no general preservation theorem is asserted.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,3,-1,2,-3], [-1,-1,-1]); simplify1!(d);

julia> Tuple(d.V)
(4, -4)
```

Role: move
"""
function simplify1!(d::URDiagram;maxiter=50)
    count=0
    while length(d.V)>2 && count<maxiter
        count+=1
        reduct3!(d); reduct5!(d); reduct6!(d)
    end
    d
end

function _simplify2!(d::URDiagram;maxiter=50,history=nothing)
    count=0
    while count<maxiter
        count+=1
        reduct3!(d); reduct8!(d)
        _,complete=reduct7!(d)
        history === nothing || push!(history,(iteration=count,vertices=copy(d.V),complete=copy(complete)))
        all(complete) && break
    end
    reduct3!(d)
    d
end

"""
    simplify2!(d::URDiagram; maxiter=50) -> URDiagram

Apply the multi-rank U(R) diagram reduction routine with an iteration bound.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

The caller supplies structurally consistent data and scalar arithmetic
admitting the selected local reductions. These are not validated at entry;
local placement and arithmetic exceptions propagate. `maxiter` should be
a nonnegative integer; this is not checked.

# Returns
The same `URDiagram`, possibly only partially reduced.

# Guarantees
1. For a nonnegative integer bound, at most `maxiter` main reduction rounds
   are applied. A final merge pass can still change the diagram when the
   bound is zero. [proved/agent: by construction]
2. The diagram is modified in place and no completion flag is returned.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The bound is not a bound on elementary moves or runtime. Successful return
does not certify complete or canonical reduction, general operator equality,
trace preservation or a represented-space invariant. Failure may leave a
partially changed diagram.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m simplify2 — source multi-rank reduction routine.

# Examples
```jldoctest
julia> code = [-1,2,-4,1,-3,4,-2,3]; d = URDiagram(); set_data_cvw!(d, 1, trans_knot_complement(code), Float64.(repeat([1,-1,-1,1], outer=4)).+1e-5); simplify2!(d);

julia> isapprox(trace2(d), -4.526591341323682e-5; rtol=1e-12)
true
```

Role: move
"""
function simplify2!(d::URDiagram;maxiter=50)
    _simplify2!(d;maxiter=maxiter)
end

"""
    simplify2_history!(d::URDiagram; maxiter=50) -> Tuple{URDiagram,Vector{NamedTuple}}

Reduce a U(R) diagram and retain the endpoint and interval-flag history of its rounds.

# Arguments
Requires U12 ([conditions](../specs/ur-diagrams.md#numbered-conditions)) with
the same caller obligations and local failures as `simplify2!`.
`maxiter` is a nonnegative integer supplied by the caller; it is not checked.

# Returns
`(d, history)`: the same mutated diagram and the
[reduction-history records](../specs/ur-diagrams.md#standard-matrix-data-and-local-reduction-records).

# Guarantees
1. The final diagram has the same reduction semantics as `simplify2!` at the same bound. At most `maxiter` history entries are produced. [proved/agent: by construction]
2. Each entry owns copied endpoint and flag vectors, independent of later diagram mutation. [proved/agent: by construction]
3. The diagram is mutated; a final merge can change it even when the history is empty. [proved/agent: by construction]

# Not guaranteed
History omits coefficients, weights and elementary witnesses, and its last
entry can precede further final changes. It is not a replay certificate or
proof of complete reduction or trace preservation.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m simplify2 — bounded reduction rounds.

# Examples
```jldoctest
julia> code = [-1,2,-4,1,-3,4,-2,3]; d = URDiagram(); set_data_cvw!(d,1,trans_knot_complement(code),Float64.(repeat([1,-1,-1,1],outer=4)).+1e-5); simplify2_history!(d);

julia> isapprox(trace2(d), -4.526591341323682e-5; rtol=1e-12)
true
```

Role: move
"""
function simplify2_history!(d::URDiagram;maxiter=50)
    history=NamedTuple[]
    _simplify2!(d;maxiter=maxiter,history=history)
    d,history
end

"""
    simplify!(d::URDiagram; maxiter=50) -> URDiagram

Reduce a weighted U(R) diagram using the multi-rank reduction contract.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` and `maxiter` have the caller obligations and exception behavior of
[`simplify2!`](@ref): consistent data, admissible scalar arithmetic, and
a nonnegative integer round bound. These are not globally checked at entry.

# Returns
The same mutated `URDiagram`.

# Guarantees
1. The result and bound semantics are those of `simplify2!(d; maxiter)`.
   [proved/agent: by construction]
2. The input diagram is modified, including any changes made before a later
   error. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
No complete reduction, canonical representative, operator identity, or
invariant of a represented space is certified by a successful return.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m simplify2 — the source reduction routine used by this entry point.

# Examples
```jldoctest
julia> code = [-1,2,-4,1,-3,4,-2,3]; d = URDiagram(); set_data_cvw!(d, 1, trans_knot_complement(code), Float64.(repeat([1,-1,-1,1], outer=4)).+1e-5); simplify!(d);

julia> isapprox(trace2(d), -4.526591341323682e-5; rtol=1e-12)
true
```

Role: move
"""
function simplify!(d::URDiagram; maxiter=50)
    simplify2!(d;maxiter=maxiter)
end

function _det(A::AbstractMatrix)
    n = size(A, 1); n == 0 && return 1; n == 1 && return A[1,1]
    acc = 0
    for j in 1:n
        minor = A[2:end, setdiff(1:n, j)]
        acc = acc + (isodd(j) ? 1 : -1) * A[1,j] * _det(minor)
    end
    acc
end

function _matrix_display(d::URDiagram)
    n = rank(d); M = Matrix{Any}(undef, n, n); fill!(M, 0)
    section(v) = count(==(0), d.V[1:vertex_index(d, v)]) + 1
    for (i, e) in enumerate(vec(d.E[1, :]))
        M[section(e), section(-e)] = d.W[i]
    end
    M
end

"""Rank-one block matrix used by MATLAB `URDiagram.trace3` and `simplify3`."""
function _block_matrix(d::URDiagram)
    rank(d)==1 || throw(ArgumentError("block matrix is defined only for rank-one diagrams"))
    block=1; active=true; blocks=Int[]
    for vertex in d.V
        !active && vertex>0 && (active=true; block+=1)
        vertex<0 && (active=false)
        push!(blocks,block)
    end
    matrix=Matrix{Any}(undef,block,block); fill!(matrix,0)
    for (index,edge) in enumerate(vec(d.E[1,:]))
        matrix[blocks[vertex_index(d,edge)],blocks[vertex_index(d,-edge)]]=d.W[index]
    end
    matrix
end

"""
    trace2(d::URDiagram) -> Any

Evaluate the unregularized matrix-trace expression assigned to a U(R) diagram.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have consistent endpoint/edge data and scalar weights supporting
determinant arithmetic and division. Its section matrix must be an
appropriate trace representation for the intended application; neither
structural validity nor that mathematical hypothesis is checked here.
For a finite value the denominator must be nonzero. Indexing and arithmetic
exceptions propagate, and floating-point singularities can yield infinities
or `NaN`.

# Returns
A scalar or symbolic expression in the arithmetic of `d.C` and its weights;
the declared return type is `Any` because those fields are unrestricted.

# Guarantees
1. The result is the section-matrix expression in
   [legacy trace formulas](../conventions/URDiagrams.md#legacy-trace-formulas).
   [proved/agent: by construction]
2. The diagram is not modified and no symbolic epsilon limit is taken.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Multiple edges occupying the same matrix position are not summed, and this
function does not check that the diagram has been sufficiently reduced.
No analytic operator trace, general sign normalization, trace-class property,
or knot/manifold invariant is established by the returned expression.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m getMatrixDisplay and trace2 — source matrix expression; the source marks its sign as provisional, and Julia omits epsilon limits.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]); dilation!(d,3,-2); add_edges!(d,(2,3)); dilation!(d,2,-1); delete_edge!(d,1);

julia> trace2(d) == -1
true
```

See also [`simplify2!`](@ref), [`trace_value`](@ref).

Role: invariant
"""
trace2(d::URDiagram) = d.C / _det(-_matrix_display(d))

"""
    trace1(d::URDiagram) -> Any

Evaluate the unregularized single-edge U(R) trace expression.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have exactly one weight; otherwise `ArgumentError` reports
`trace1 requires one edge`. Its coefficient and weight must support negation
and division. A nonzero weight is the caller's responsibility for a finite
value; singular arithmetic can return an infinity or raise an exception.
The vertex word, rank and edge-table consistency are not checked.

# Returns
The scalar or symbolic single-edge expression in
[legacy trace formulas](../conventions/URDiagrams.md#legacy-trace-formulas),
with type determined by the input arithmetic and therefore declared as `Any`.

# Guarantees
1. On supported nonsingular scalar data, the result is the indicated
   single-edge expression. [proved/agent: by construction]
2. The diagram is not modified and no epsilon limit is taken.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
One stored weight alone does not certify an admissible trace diagram.
No analytic trace or knot/manifold invariant is established. The epsilon
limits of the MATLAB `trace1` routine are not part of this result.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m trace1 — source single-edge expression before epsilon limits.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]); dilation!(d,3,-2); add_edges!(d,(2,3)); dilation!(d,2,-1); delete_edge!(d,1);

julia> trace1(d) == -1
true
```

Role: invariant
"""
trace1(d::URDiagram) = length(d.W) == 1 ? -d.C/d.W[1] : throw(ArgumentError("trace1 requires one edge"))

"""
    trace3(d::URDiagram) -> Any

Evaluate the rank-one block-matrix U(R) trace expression.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have rank one; otherwise `ArgumentError` reports
`trace3 requires rank one`. The caller supplies consistent diagram data,
supported coefficient/weight arithmetic and a nonsingular denominator for
a finite value. Structural validity and the operator interpretation of the
block matrix are not checked. Indexing and arithmetic failures propagate.

# Returns
A scalar or symbolic expression in the arithmetic of the diagram's
coefficient and weights. Its type is unrestricted (`Any`).

# Guarantees
1. The result is the rank-one block-matrix expression in
   [legacy trace formulas](../conventions/URDiagrams.md#legacy-trace-formulas).
   [proved/agent: by construction]
2. For a valid single-edge rank-one diagram at a nonsingular weight, it
   equals `trace1(d)`. [proved/agent: by construction]
3. The diagram is not modified and no epsilon limit is taken.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The call does not certify an analytic operator trace, equality with every
other trace routine on unreduced data, or an invariant of a knot or manifold.
Weights sharing a block-matrix position are not summed. Singular arithmetic
is not regularized.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m getBlockMatrix and trace3 — source block-matrix expression; analytic hypotheses and a general topological identification are not asserted.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]);

julia> trace3(d) == -1/(1+(-1)*(sym(:x)+(-1))*(sym(:z)+(sym(:y)+(-1))))
true
```

Role: invariant
"""
trace3(d::URDiagram) = rank(d) == 1 ? -d.C/_det(_block_matrix(d)+standard_matrix_adjustment(_block_matrix(d))) :
    throw(ArgumentError("trace3 requires rank one"))
"""
    trace_value(d::URDiagram) -> Any

Return the U(R) matrix-trace expression after reducing the supplied diagram in place.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must satisfy the caller obligations of [`simplify!`](@ref) and, after
reduction, [`trace2`](@ref). Reduction uses the default iteration bound.
Local placement, indexing and arithmetic exceptions propagate; no global
validity or finite-trace check is performed.

# Returns
The scalar or symbolic `trace2` value of the reduced diagram, with type
determined by its coefficient arithmetic (`Any`).

# Guarantees
1. Successful return equals `trace2(d)` for the diagram's state on return.
   [proved/agent: by construction]
2. The supplied diagram is modified by reduction, including when a later
   trace evaluation fails. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
This is not a nonmutating accessor. Complete reduction, epsilon limits,
analytic trace existence, general preservation of the original expression,
and knot/manifold invariance are not certified.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m simplify2 and trace2 — source operations; Julia uses this route for every rank and does not take epsilon limits.

# Examples
```jldoctest
julia> code = [-1,2,-4,1,-3,4,-2,3]; d = URDiagram(); set_data_cvw!(d, 1, trans_knot_complement(code), Float64.(repeat([1,-1,-1,1], outer=4)).+1e-5);

julia> isapprox(trace_value(d), -4.526591341323682e-5; rtol=1e-12)
true
```

Role: invariant
"""
trace_value(d::URDiagram) = (simplify!(d); trace2(d))

"""
    odata(d::URDiagram) -> String

Render a unit-weight U(R) diagram as its source-format O-data word.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

Every weight must be numeric with absolute value one; otherwise
`ArgumentError` reports `O-data requires weights ±1`. The weights must
also support comparison with zero. The caller supplies consistent endpoint
and edge data, which are not validated here. Indexing or unsupported
comparison errors propagate.

# Returns
A `String` encoding the diagram's signed endpoint word and weight signs in
the O-data rendering format.

# Guarantees
1. The result contains the source-format rendering of each endpoint and
   separator, without changing the diagram. [proved/agent: by construction]
2. The coefficient `C` is not encoded. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The string is not a certified branched spine or manifold description, nor
a lossless serialization of the entire diagram. Unlike the source MATLAB
routine, this call takes no epsilon limits before checking weights.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m getOdata — source text rendering, before any independent geometric interpretation.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]); dilation!(d,3,-2); add_edges!(d,(2,3)); dilation!(d,2,-1); delete_edge!(d,1); odata(d);

julia> trace3(d) == -1
true
```

Role: convert
"""
function odata(d::URDiagram)
    all(w -> w isa Number && abs(w) == 1, d.W) || throw(ArgumentError("O-data requires weights ±1"))
    parts = String[]
    for v in reverse(d.V)
        v == 0 ? push!(parts, "⊗") : push!(parts, "v_{$(abs(v))}^{" * (v > 0 ? "+" : "-") * (get_weight(d,v)>0 ? "r}" : "l}"))
    end
    join(parts)
end

"""
    weight_expression(d::URDiagram) -> String

Render the U(R) diagram's endpoint word as a formal word of weighted-operator labels.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` must have consistent endpoint/edge data and weights that support a
Boolean comparison with zero. These are the caller's responsibility;
indexing or comparison errors propagate. No unit-weight or structural
validity check is performed.

# Returns
A `String` in the repository's weighted-operator label notation.

# Guarantees
1. The result is a textual word reflecting the endpoint order and the
   positive/nonpositive branch of each weight comparison.
   [proved/agent: by construction]
2. The diagram is not modified. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
Weight magnitudes and coefficient `C` are not encoded. The string is not
an evaluated operator or a lossless diagram serialization. No epsilon
limits or knot/manifold-invariance checks are performed.

# Sources
cf:references/library/articles/MST-Rx/sources/selected-main/04_trace_computation/SRC-262_computations_URDiagram.m getWExpression — source operator-word rendering; Julia uses its own text notation and omits the source's epsilon limits and unit-weight check.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-3,3,-2,2,-1], [1,1,1]); dilation!(d,3,-2); add_edges!(d,(2,3)); dilation!(d,2,-1); delete_edge!(d,1); weight_expression(d);

julia> trace3(d) == -1
true
```

Role: convert
"""
function weight_expression(d::URDiagram)
    parts = String[]
    for v in d.V
        v == 0 ? push!(parts, "⊗") : push!(parts, (get_weight(d,v)>0 ? "" : "bar(") * "WO_$(abs(v))$(v>0 ? 2 : 1)" * (get_weight(d,v)>0 ? "" : ")"))
    end
    join(parts)
end

"""
    unicode_art(d::URDiagram) -> String

Return a textual diagnostic of a weighted U(R) diagram.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`d` is a `URDiagram` whose coefficient and weights support text display.
No validity check is performed; malformed tables or display methods can
raise indexing or display errors.

# Returns
A newline-terminated `String` showing coefficient, rank, the signed endpoint
word and available edge labels with their weights.

# Guarantees
1. The displayed data are taken from the supplied diagram without changing
   it. [proved/agent: by construction]
2. The text includes a rank and coefficient line even for empty data.
   [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
This is not a planar embedding, a complete edge-table dump, a canonical
serialization, or an invariant of a represented operator or space. Invalid
data can be displayed; output depends on coefficient display methods.

# Sources
cf:lib/topology.jl/src/ur_diagram.jl, unicode_art — structural data checking and diagnostic display.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, [1,-2,-3,3,-1,2], [sym(:x),sym(:y),sym(:z)]); unicode_art(d);

julia> validate_data(d)
true
```

Role: draw
"""
function unicode_art(d::URDiagram)
    io = IOBuffer()
    println(io, "C=$(d.C), rank=$(rank(d))")
    println(io, join((v == 0 ? "⊗" : v > 0 ? "+$(v)" : "$(v)" for v in d.V), "─"))
    for (e,w) in zip(vec(d.E[1,:]), d.W); println(io, "e$(e): $(w)"); end
    String(take!(io))
end

"""
    trans_knot_complement(code) -> Vector{Int}

Return the U(R) endpoint word assigned by the repository's knot-complement conversion recipe.

# Arguments
The data definitions and numbered conditions are in the
[specification](../specs/ur-diagrams.md).

`code` must be collectable and convertible entrywise to `Int`. The caller
supplies signed knot-code data in the source recipe's convention; no knot-code
validity, pairing, crossing-sign or planarity check is made. Conversion errors
propagate. Use labels small enough to avoid machine-integer overflow; this
is not checked.

# Returns
A `Vector{Int}` of endpoint data for subsequent use with [`set_data_cvw!`](@ref),
of length `4*length(code)+1` when collection and arithmetic are supported.
No weights or coefficient are returned.

# Guarantees
1. The recipe returns a new endpoint-word vector with one inserted section
   separator. The input is not modified. [proved/agent: by construction]
2. For paired nonzero input labels without arithmetic overflow, the output
   contains paired signed endpoint labels. [proved/agent: by construction]

# Conventions
[Coordinate conventions](../conventions/URDiagrams.md).

# Not guaranteed
The function alone does not construct a weighted diagram, identify the
homeomorphism type of a knot complement, or establish an associated trace
as a complement invariant. Invalid or overflowing input is not certified.

# Sources
cf:lib/topology.jl/examples/URDforKnotCmpl.jl run_URDforKnotCmpl — use of the endpoint recipe for the registered unknot, trefoil and figure-eight computations; no general complement-identification theorem is asserted.

# Examples
```jldoctest
julia> d = URDiagram(); set_data_cvw!(d, 1, trans_knot_complement([1,-1]), [1,-1,-1,1]);

julia> validate_data(d)
true
```

Role: convert
"""
function trans_knot_complement(code)
    a = Int.(collect(code)); pos = Int[]; neg = Int[]
    for v in a
        push!(pos, 4v - (v>0 ? 2 : 0) + (v<0 ? 3 : 0))
        push!(pos, 4v - (v>0 ? 3 : 0) + (v<0 ? 1 : 0))
        push!(neg, 4v - (v>0 ? 0 : 0) + (v<0 ? 2 : 0))
        push!(neg, 4v - (v>0 ? 1 : 0) + (v<0 ? 0 : 0))
    end
    vcat(pos, 0, reverse(neg))
end
