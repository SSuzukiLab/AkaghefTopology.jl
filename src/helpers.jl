"""
    SmallObject(value) -> SmallObject

Wrap a value in a mutable typed cell.

# Arguments
`value` can have any type. `SmallObject{T}(value)` additionally converts to
`T`, propagating conversion errors. No mathematical validity check is made.

# Returns
`SmallObject{typeof(value)}` for the inferred constructor, or `SmallObject{T}`
for the explicitly typed constructor, with one mutable field `value`.

# Guarantees
1. The inferred constructor retains the supplied value; the typed constructor
   retains its converted value. [proved/agent: by construction]
2. Wrapping does not modify or copy the supplied object. Assigning the field
   replaces the cell's value subject to its field type.
   [proved/agent: by construction]

# Not guaranteed
This cell adds no mathematical structure, ownership isolation, synchronization,
or validity certificate to the wrapped value.

# Sources
cf:lib/topology.jl/src/helpers.jl, SmallObject — the current cell and polynomial-coordinate operations.

# Examples
```jldoctest
julia> s = SmallObject(3); s.value = 4;

julia> s.value
4
```

Role: convert
"""
mutable struct SmallObject{T}
    value::T
end

"""
    bernstein(n::Integer, t) -> Matrix

Evaluate the degree-n Bernstein basis at one or more parameters.

# Arguments
`n` must satisfy B1 of the [curve specification](../specs/curves.md#numbered-conditions);
otherwise `ArgumentError` reports `degree must be positive`.
`t` is a number or a collectable one-dimensional sequence satisfying B2.
Violation of B2 raises
`ArgumentError` naming that interval. Comparison, binomial, power and
arithmetic errors propagate. Exact arithmetic claims require coefficient
types and degree sizes that do not overflow.

# Returns
A matrix with `n+1` basis rows in the
[sample-matrix convention](../specs/curves.md#definitions).
Element type follows the parameter arithmetic.

# Guarantees
1. On supported arithmetic, column entries are the Bernstein basis values
   at the corresponding parameter, in increasing basis-index order.
   [proved/agent: by construction]
2. The result is newly allocated and inputs are not modified.
   [proved/agent: by construction]

# Not guaranteed
Degree zero is not supported. Floating-point partition-of-unity identities
need not hold bit-for-bit, and no overflow or conditioning bound is provided.

# Sources
cf:lib/topology.jl/src/helpers.jl, bernstein — the current cell and polynomial-coordinate operations.

# Examples
```jldoctest
julia> b = bernstein(3, [0.0,0.5,1.0]);

julia> (size(b), all(isapprox.(sum(b; dims=1), 1)))
((4, 3), true)
```

See also [`bezier`](@ref).

Role: invariant
"""
function bernstein(n::Integer, t)
    n >= 1 || throw(ArgumentError("degree must be positive"))
    ts = t isa Number ? [t] : collect(t)
    all(x -> 0 <= x <= 1, ts) || throw(ArgumentError("parameters must lie in [0,1]"))
    [binomial(n,i) * (1-x)^(n-i) * x^i for i in 0:n, x in ts]
end

"""
    bezier(t, control_points::AbstractMatrix) -> AbstractMatrix

Evaluate the polynomial Bezier curve with the supplied control-point columns.

# Arguments
`control_points` must satisfy B3 of the
[curve specification](../specs/curves.md#numbered-conditions); a failed
column-count check raises `ArgumentError` with
`at least two control points are required`.
`t` satisfies B2, checked by [`bernstein`](@ref). Coordinate and
parameter types must support the required scalar arithmetic and matrix
multiplication; arithmetic and conversion failures propagate.

# Returns
A matrix in the [sample-matrix convention](../specs/curves.md#definitions),
with the same coordinate-row count as `control_points`.

# Guarantees
1. For supported arithmetic, each result column is the Bezier point at its
   parameter for the supplied ordered control points.
   [proved/agent: by construction]
2. The function does not modify the control points or parameters.
   [proved/agent: by construction]

# Not guaranteed
This does not interpolate all control points, validate a link embedding,
or bound floating-point error. Out-of-interval extrapolation is rejected
by the Bernstein evaluator.

# Sources
cf:lib/topology.jl/src/helpers.jl, bezier — the current cell and polynomial-coordinate operations.

# Examples
```jldoctest
julia> p = [0.0 1 2; 0.0 1 0];

julia> bezier(0.5, p)[:,1] ≈ [1.0,0.5]
true
```

Role: convert
"""
function bezier(t, control_points::AbstractMatrix)
    size(control_points,2) >= 2 || throw(ArgumentError("at least two control points are required"))
    control_points * bernstein(size(control_points,2)-1,t)
end

"""
    MovePattern(parameter::AbstractString, gauss, orientation, capacity=Int[]; weights=nothing) -> MovePattern

Construct a local move template with its boundary and weight-rule metadata.

# Arguments
Use the fields and default conventions of [MovePattern](../specs/virtual-link-moves.md#51-movepattern).
Checked: the capacity and weight-rule collections must have one entry per
fragment; `DimensionMismatch` identifies the mismatched collection. Numeric
and string conversion errors propagate.
Caller's responsibility: local labels, orientations, capacities and rule
lengths must satisfy the linked template contract for later move execution.
The full eight-field constructor only performs field conversions.

# Returns
A `MovePattern` with the field meanings in the linked specification. The
convenience constructor derives endpoint and vertex metadata from the fragments.

# Guarantees
1. The convenience constructor supplies the documented default metadata and
   weight rules when omitted. [proved/agent: by construction]
2. Inputs are not modified. The full field constructor may retain arrays
   supplied by the caller. [proved/agent: by construction]

# Not guaranteed
Construction is not an applicability or geometric-preservation check.
An unspecified weight rule is not permission to use a weighted move.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §§5.1–5.2 — template fields and weight-rule coordinates.

# Examples
```jldoctest
julia> p = MovePattern("transport", [[1],[-1]], [1], [1,2]);

julia> (p.vertex_count, p.gauss_first, length.(p.weights))
(1, [1, -1], [2, 2])
```

Role: construct
"""
struct MovePattern
    parameter::String
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    capacity::Vector{Int}
    gauss_first::Vector{Int}
    gauss_last::Vector{Int}
    vertex_count::Int
    weights::Vector{Vector{Float64}}
end

function MovePattern(parameter::AbstractString, gauss, orientation, capacity=Int[]; weights=nothing)
    fragments=[Int.(collect(fragment)) for fragment in gauss]
    caps=isempty(capacity) ? collect(1:length(fragments)) : Int.(collect(capacity))
    length(caps)==length(fragments) || throw(DimensionMismatch("one capacity label per Gauss fragment is required"))
    firsts=[isempty(fragment) ? 0 : first(fragment) for fragment in fragments]
    lasts=[isempty(fragment) ? 0 : last(fragment) for fragment in fragments]
    labels=abs.(reduce(vcat,fragments;init=Int[]))
    vertex_count=isempty(labels) ? 0 : maximum(labels)
    fragment_weights=weights===nothing ?
        [zeros(Float64,length(fragment)+1) for fragment in fragments] :
        [Float64.(collect(weight)) for weight in weights]
    length(fragment_weights)==length(fragments) ||
        throw(DimensionMismatch("one weight rule per Gauss fragment is required"))
    MovePattern(String(parameter),fragments,Int.(collect(orientation)),caps,
                firsts,lasts,vertex_count,fragment_weights)
end

"""Apply MATLAB `VirtualLink.calcReverseStrand` to a local move template."""
function _reverse_template(gauss,orientation,flags)
    reversed=[Int.(copy(fragment)) for fragment in gauss]
    changed=Int[]
    for index in eachindex(reversed)
        flags[index] || continue
        reversed[index]=reverse(reversed[index])
        append!(changed,reversed[index])
    end
    result=Int.(copy(orientation))
    positive=filter(>(0),changed)
    paired=filter(vertex -> -vertex in changed,positive)
    affected=setdiff(unique(abs.(changed)),unique(paired))
    virtual=filter(vertex -> 1<=vertex<=length(result) && result[vertex]==0,affected)
    for vertex in affected
        result[vertex]*=-1
    end
    # MATLAB keeps a virtual orientation at zero and swaps the two signed visits
    # when only one side of that virtual crossing belongs to a reversed strand.
    for vertex in virtual
        for fragment in reversed, index in eachindex(fragment)
            if fragment[index]==vertex
                fragment[index]=-vertex
            elseif fragment[index]==-vertex
                fragment[index]=vertex
            end
        end
    end
    reversed,result
end

function _mp_patterns(side::Symbol)
    labels=["A","B","C","D"]
    gauss=side===:left ? [
        [[-1],[1,2],[-2]],
        [[1],[-1,-2],[2]],
        [[-1,2],[1],[-2]],
        [[2,-1],[1],[-2]],
    ] : [
        [[1,-3,2],[3],[-1,-2]],
        [[-1,3,-2],[-3],[1,2]],
        [[2,-3],[1,3],[-1,-2]],
        [[-3,1],[3,2],[-1,-2]],
    ]
    orientations=side===:left ? [[1,1],[-1,-1],[-1,1],[1,-1]] :
                                 [[1,-1,1],[-1,1,-1],[1,1,-1],[-1,-1,1]]
    flags=[[false,false,false],[false,false,true],[false,true,false],[false,true,true]]
    patterns=MovePattern[]
    for (label,fragments,orientation) in zip(labels,gauss,orientations)
        for (number,reverse_flags) in enumerate(flags)
            transformed,transformed_orientation=_reverse_template(fragments,orientation,reverse_flags)
            parameter=label*string(number)
            weights=[zeros(Float64,length(fragment)+1) for fragment in transformed]
            # Integral MP weights: MST25 Fig. 3.4 (Z_framed_MP_move.pdf),
            # with the individually labelled diagrams in Taguchi's thesis
            # source, MP_list_AB_int.jpg and MP_list_CD_int.jpg.  Weights
            # run incoming boundary, internal edges, outgoing boundary.
            # In A4/B4 fragment 1 runs lower-left -> lower vertex 1 ->
            # right vertex 3 -> upper vertex 2 -> upper-left.  The source
            # figure places +1 on its incoming leg and -1 on 1 -> 3;
            # In A2/B2 the -1 is on 3 -> 2 and +1 on the outgoing leg.
            # The former fragment-3 rule [1,-1,0] is H-related at vertex 2,
            # but is not the weights on this signed-port source diagram.
            if side===:right && parameter in ("A4","B4")
                weights[1]=Float64[1,-1,0,0]
            elseif side===:right && parameter in ("A2","B2")
                weights[1]=Float64[0,0,-1,1]
            elseif side===:right && parameter=="C3"
                weights[2]=Float64[0,1,-1]
            elseif side===:right && parameter=="D3"
                weights[2]=Float64[-1,1,0]
            end
            push!(patterns,MovePattern(parameter,transformed,transformed_orientation;weights=weights))
        end
    end
    patterns
end

"""
    make_move_data() -> Dict{String,Vector{MovePattern}}

Construct the library's catalogue of local virtual-link move patterns.

# Arguments
Checked: no arguments or input conditions.
Caller: select and match patterns using the
[candidate rules](../specs/virtual-link-moves.md#6-candidate-contract).

# Returns
A dictionary of pattern vectors, keyed by the
[catalogue names](../specs/virtual-link-moves.md#53-catalogue-keys).

# Guarantees
1. The values are the local diagram templates and weight rules of the
   [move catalogue](../specs/virtual-link-moves.md#5-move-data-contract),
   including the sixteen MP variants on each side. [proved/agent: by construction]
2. Each call constructs fresh catalogue data and writes no external file.
   [proved/agent: by construction]

# Not guaranteed
A template alone does not certify applicability, weighted admissibility, or
preservation of a represented manifold on an arbitrary input.

# Sources
cf: [virtual-link moves §5](../specs/virtual-link-moves.md#5-move-data-contract) — catalogue and local template semantics.

# Examples
```jldoctest
julia> d = make_move_data(); (length(d), length(d["MP-L"]), length(d["MP-R"]))
(16, 16, 16)
```

Role: construct
"""
function make_move_data()
    Dict(
        "MP-L" => _mp_patterns(:left),
        "MP-R" => _mp_patterns(:right),
        "02-L" => [MovePattern("dir",[Int[],Int[]],Int[];weights=[Float64[0],Float64[0]])],
        # MST25 Fig. 2.9 + Proposition A.1: switched = Γ× ∘ standard ∘ Γ×.
        # Only both crossings' over/under data change, never strand directions.
        # [[MST25#^mst25-fig-0066]], [[MST25#^mst25-blk-0188]].
        # No framed lift of this conjugation has been specified.
        "02-R" => [MovePattern("inv",[[-1,-2,2,1]],[1,-1];weights=[Float64[0,0,1,-1,0]]),
                   MovePattern("inv-switched",[[1,2,-2,-1]],[-1,1];weights=[fill(NaN,5)])],
        "CP-L" => [MovePattern("dir",[[1,-2,2,-3,3,-1]],[-1,-1,-1])],
        "CP-R" => [MovePattern("inv",[[1,2,-5,3,-2,-3,-4,4,-1,5]],[1,1,1,-1,-1])],
        # Rooted reduced CP: R6 -> R0 (2 -> 5), not a cyclically rerooted patch.
        # Research definition: PJ260830 / 260922_reduced_cp_move.
        "RCP-L" => [MovePattern("dir",[[-1,2,-2,1]],[-1,1])],
        "RCP-R" => [MovePattern("dir",[[-1,-2,3,-4,2,4,5,-5,1,-3]],[-1,-1,1,-1,1])],
        # INS24 Fig.16, independently transcribed in the source visual record.
        # S is not an alias for BP's 3-to-5 CP: its local patch is 0-to-4.
        "S-L" => [MovePattern("INS",[Int[]],Int[])],
        "S-R" => [MovePattern("INS",[[1,2,-3,-4,4,3,-2,-1]],[-1,-1,-1,-1])],
        "PS-L" => [MovePattern("",[Int[],Int[]],Int[])],
        "PS-R" => [
            MovePattern("1",[[1,2],[-1,-2]],[1,-1]),
            MovePattern("2",[[1,2],[-1,-2]],[-1,1]),
            MovePattern("3",[[1,2],[-2,-1]],[1,-1]),
            MovePattern("4",[[1,2],[-2,-1]],[-1,1]),
        ],
        "BMP-L" => [
            MovePattern("E1",[[1,2],[-1],[-2]],[-1,1],[1,2,3]),
            MovePattern("F1",[[-1,-2],[1],[2]],[1,-1],[1,2,3]),
            MovePattern("E2",[[1,2],[-1],[-2]],[-1,1],[1,2,3]),
            MovePattern("F2",[[-1,-2],[1],[2]],[1,-1],[1,2,3]),
        ],
        "BMP-R" => [
            MovePattern("E1",[[3],[2],[1],[-1,-2,-3]],[-1,1,1],[1,3,2,0]),
            MovePattern("F1",[[-3],[-2],[-1],[1,2,3]],[1,-1,-1],[1,3,2,0]),
            MovePattern("E2",[[3],[2],[1],[-3,-2,-1]],[1,-1,-1],[1,3,2,0]),
            MovePattern("F2",[[-3],[-2],[-1],[3,2,1]],[-1,1,1],[1,3,2,0]),
        ],
        "B02-L" => [MovePattern("",[Int[],Int[]],[1,-1],[1,2])],
        "B02-R" => [
            MovePattern("1",[[1],[2],[-1,-2]],[1,-1],[2,1,3]),
            MovePattern("2",[[1],[2],[-2,-1]],[-1,1],[2,1,3]),
        ],
    )
end
