"""Small, explicit taxonomy for the common move graph.

The numeric values are intentionally user-facing: Area 2, Area 2.5, and Area
3 are different mathematical domains, not storage buckets.  A representation
change is never a topological move.
"""
module DomainTaxonomy

export AREA_2, AREA_2_5, AREA_3, VALID_AREAS, normalize_area,
       normalize_representation, is_topological_move, validate_topological_move,
       move_area

const AREA_2 = 2.0
const AREA_2_5 = 2.5
const AREA_3 = 3.0
const VALID_AREAS = (AREA_2, AREA_2_5, AREA_3)

"""
    normalize_area(area) -> Float64

Return the registered domain number represented by `area`.

# Arguments
`area` must convert to `Float64` and match an entry of `VALID_AREAS` exactly.
Conversion errors propagate; an unregistered value raises `ArgumentError`
listing the accepted areas. The taxonomy belongs to the
[MoveGraphStore module](../specs/move-graph-store.md#definitions).

# Returns
The accepted `Float64` domain number.

# Guarantees
1. The result is in `VALID_AREAS` and equals `Float64(area)`.
   [proved/agent: by construction]
2. No supplied object is modified. [proved/agent: by construction]

# Not guaranteed
Acceptance of a number does not establish the domain of a mathematical object.

# Sources
cf:lib/topology.jl/src/domain_taxonomy.jl, normalize_area — the common-store domain and representation vocabulary; name recognition does not certify diagram moves.

# Examples
```jldoctest
julia> Topology.MoveGraphStore.DomainTaxonomy.normalize_area(Topology.AREA_3) == 3.0
true
```

Role: convert
"""
function normalize_area(area)
    x = Float64(area)
    any(a -> isapprox(x, a; atol=0, rtol=0) , VALID_AREAS) ||
        throw(ArgumentError("area must be 2, 2.5, or 3"))
    x
end

"""
    normalize_representation(x) -> String

Return the normalized textual label of a stored representation.

# Arguments
`x` must support conversion with `String`; conversion errors propagate.
No registry membership or nonemptiness check is performed. Representation
labels belong to [MoveGraphStore](../specs/move-graph-store.md#definitions).

# Returns
A `String` suitable for the representation field of [`CanonicalState`](@ref).

# Guarantees
1. The result equals the representation label produced by the serialization
   constructor of `CanonicalState` from `x`. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not convert a diagram or verify that a serialization has the named
representation. Distinct inputs may have the same normalized label.

# Sources
cf:lib/topology.jl/src/domain_taxonomy.jl, normalize_representation — the common-store domain and representation vocabulary; name recognition does not certify diagram moves.

# Examples
```jldoctest
julia> s = CanonicalState("state-a"; area=Topology.AREA_2, representation="KnotGauss");

julia> Topology.MoveGraphStore.DomainTaxonomy.normalize_representation("KnotGauss") == s.representation
true
```

Role: convert
"""
normalize_representation(x) = lowercase(strip(String(x)))

const _AREA2 = Set(["r1", "r2", "r3", "reidemeister1", "reidemeister2", "reidemeister3"])
const _AREA25 = Set(["mixed", "mixed_move", "mixed_link", "lo", "lo_graph", "lo-graph",
                     "mixed/lo-graph", "mixed/lo-graph link", "lo-graph link"])
const _AREA3 = Set(["mp", "ps", "cp", "s", "0-2", "02", "0_2",
                    "mp_inverse", "ps_inverse", "cp_inverse", "s_inverse",
                    "bmp", "b02", "h",
                    "3-manifold mp", "3-manifold ps", "3-manifold cp", "3-manifold s"])

"""
    is_topological_move(move) -> Bool

Recognize a name in the repository's topological-move vocabulary.

# Arguments
`move` must convert to `String`; conversion errors propagate. Recognition
uses the [domain vocabulary](../specs/move-graph-store.md#domain-vocabulary).
No diagram or witness is supplied or validated.

# Returns
`Bool`: whether the normalized name belongs to that vocabulary.

# Guarantees
1. Recognized names can be passed to `move_area` and `validate_topological_move` with that returned area. [proved/agent: by construction]
2. The argument is not modified. [proved/agent: by construction]

# Not guaranteed
Name recognition does not establish that any transition is a valid move
or that a represented space is preserved.

# Sources
cf:lib/topology.jl/docs/specs/move-graph-store.md Domain vocabulary — permitted move names and areas.

# Examples
```jldoctest
julia> Topology.MoveGraphStore.DomainTaxonomy.is_topological_move("MP")
true

julia> Topology.MoveGraphStore.DomainTaxonomy.is_topological_move("representation_transform")
false
```

Role: condition
"""
function is_topological_move(move)
    s = lowercase(strip(String(move)))
    isempty(s) && return false
    occursin("transform", s) && return false
    occursin("representation", s) && return false
    occursin("serialize", s) && return false
    s in _AREA2 || s in _AREA25 || s in _AREA3
end

"""
    move_area(move) -> Float64

Return the domain assigned to a recognized topological move name.

# Arguments
`move` must convert to `String` and name a move in the
[MoveGraphStore taxonomy](../specs/move-graph-store.md#definitions).
The lookup accepts surrounding whitespace and ignores letter case.
An unknown name raises `ArgumentError` identifying the unknown or
non-topological move; string conversion errors propagate.

# Returns
The move name's domain as one of the exported area constants.

# Guarantees
1. A successful result agrees with `validate_topological_move(move, result)`.
   [proved/agent: by construction]
2. The argument is not modified. [proved/agent: by construction]

# Not guaranteed
Recognizing a name neither supplies a move on a diagram nor verifies a
transition between two mathematical objects.

# Sources
cf:lib/topology.jl/src/domain_taxonomy.jl, move_area — the common-store domain and representation vocabulary; name recognition does not certify diagram moves.

# Examples
```jldoctest
julia> Topology.MoveGraphStore.DomainTaxonomy.move_area("MP") == Topology.AREA_3
true
```

Role: access
"""
function move_area(move)
    s = lowercase(strip(String(move)))
    s in _AREA2 && return AREA_2
    s in _AREA25 && return AREA_2_5
    s in _AREA3 && return AREA_3
    throw(ArgumentError("unknown or non-topological move type: $move"))
end

"""
    validate_topological_move(move, area) -> Bool

Check that a recognized move name belongs to the supplied domain.

# Arguments
`move` is convertible to `String`; `area` is accepted by [`normalize_area`](@ref).
Checks the name and domain using the
[MoveGraphStore taxonomy](../specs/move-graph-store.md#definitions).
`ArgumentError` identifies a forbidden or unknown move, an invalid area,
or a mismatch between the move's assigned area and the supplied area.
Conversion failures propagate.

# Returns
`true` on success; an invalid pair raises an exception rather than returning
`false`.

# Guarantees
1. Success means `move_area(move) == normalize_area(area)` and the move name
   is recognized as topological. [proved/agent: by construction]
2. Neither argument is modified. [proved/agent: by construction]

# Not guaranteed
This name-level check does not establish that a proposed diagram transition
is an instance of the move or preserves a represented space.

# Sources
cf:lib/topology.jl/src/domain_taxonomy.jl, validate_topological_move — the common-store domain and representation vocabulary; name recognition does not certify diagram moves.

# Examples
```jldoctest
julia> Topology.MoveGraphStore.DomainTaxonomy.validate_topological_move("MP", Topology.AREA_3)
true
```

Role: condition
"""
function validate_topological_move(move, area)
    is_topological_move(move) ||
        throw(ArgumentError("representation changes and non-topological moves cannot enter transitions"))
    normalize_area(area) == move_area(move) ||
        throw(ArgumentError("move $(move) belongs to Area $(move_area(move)), not Area $(area)"))
    true
end

end
