# External states of one fixed combing. The integer-weight solver keeps
# `:integer_solution`, `:rational_only`, and `:inconsistent`. This file does
# not change a combing and does not search for a CP move.

import AbstractAlgebra

"""
    CMR3State(kind, order, diagnostic) -> CMR3State

Store an external classification state for one fixed combing.

# Arguments
`kind::Symbol`, `order::Union{Nothing,Int}` and
`diagnostic::Union{Nothing,Symbol}`, or convertible values. Conversion errors
propagate. Status consistency and the asserted obstruction order are not
checked; use the classification functions for computed records.

# Returns
A state tag, an optional finite obstruction order, and an optional solver
diagnostic, as supplied.

# Guarantees
1. Fields preserve the supplied values after conversion; no classification
   or change of combing is performed. [proved/agent: by construction]
2. Already-typed arguments are not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
The constructor does not certify an Euler class or a framing. In particular,
a mod-two solution does not certify an integer second field. The [CMR3 specification](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the status and obstruction interpretation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §10.3 — combing claims require a compatible geometric correspondence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> s = CMR3State(:combed, nothing, nothing); (s.kind, s.order, s.diagnostic)
(:combed, nothing, nothing)
```

See also [`classify_euler_class`](@ref), [`classify_mod2_solution`](@ref).

Role: construct
"""
struct CMR3State
    kind::Symbol
    order::Union{Nothing,Int}
    diagnostic::Union{Nothing,Symbol}
end

"""
    ModularWitness(modulus, u) -> ModularWitness

Store a candidate modular witness for a nonzero Euler-obstruction class.

# Arguments
`modulus::Int` and `u::Vector{Int}`, or convertible values. Conversion errors
propagate. The constructor does not check the modulus, vector dimensions or
annihilation/residue conditions; use `modular_witness_check` for those checks.

# Returns
A `ModularWitness` with the supplied modulus and witness vector.

# Guarantees
1. Fields preserve supplied values after conversion; a typed witness vector
   is retained without mutation. [proved/agent: by construction]
2. The record is candidate algebraic evidence, with cellular coordinates
   governed by [the cellular-basis convention](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
   [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Construction does not prove a nonzero cokernel class or any geometric
framing obstruction. The [CMR3 specification](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the modular-obstruction definition.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> w = ModularWitness(2, [1]); (w.modulus, Tuple(w.u))
(2, (1,))
```

See also [`modular_witness_check`](@ref), [`classify_euler_class`](@ref).

Role: construct
"""
struct ModularWitness
    modulus::Int
    u::Vector{Int}
end

"""
    CMR3Reading(state, witness_annihilates, witness_residue, integer_solution) -> CMR3Reading

Store an Euler-obstruction classification and its optional algebraic witnesses.

# Arguments
`state::CMR3State`, `witness_annihilates::Union{Nothing,Bool}`,
`witness_residue::Union{Nothing,Int}` and
`integer_solution::Union{Nothing,Vector{BigInt}}`, or convertible values.
Conversion errors propagate. Compatibility of these fields is the caller's
responsibility; this constructor does not run the classification checks.

# Returns
A `CMR3Reading` with the state, optional modular-witness results and optional
integer solution for the producer's fixed-spine equation.

# Guarantees
1. Fields retain the supplied values after conversion. A typed solution
   vector is shared and not modified. [proved/agent: by construction]
2. The record carries the supplied algebraic reading; it supplies no
   geometric modification of the spine. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
The presence of an integer vector or a Boolean witness result does not prove
the asserted equation. The [CMR3 specification](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the obstruction convention; obtain computed records from `classify_euler_class`.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.5, 10.3 — cellular cochain coordinates and geometric scope of combing claims.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> r = CMR3Reading(CMR3State(:combed, nothing, nothing), nothing, nothing, nothing); (r.state.kind, r.integer_solution)
(:combed, nothing)
```

See also [`classify_euler_class`](@ref).

Role: convert
"""
struct CMR3Reading
    state::CMR3State
    witness_annihilates::Union{Nothing,Bool}
    witness_residue::Union{Nothing,Int}
    integer_solution::Union{Nothing,Vector{BigInt}}
end

"""
    CellularBasis(vertices, edges, faces) -> CellularBasis

Store ordered cell labels for the three degrees of a cellular-map witness.

# Arguments
`vertices::Vector{String}`, `edges::Vector{String}`, `faces::Vector{String}`,
or values convertible to those types. Conversion errors propagate. Empty
lists and repeated labels are accepted.

# Returns
`CellularBasis` with the supplied vertex, edge and face label lists, in the
supplied order. It can be used for either side of a `CellularMapWitness`;
orientations and coordinates follow [cellular basis conventions](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).

# Guarantees
1. Fields preserve the supplied values after conversion; typed vectors are
   retained without copying. Already-typed arguments are not modified. [proved/agent: by construction]
2. The record supplies cell labels for the correspondence between the
   source complex and spine in
   [the conversion specification](../specs/heegaard-koda-conversion.md#64-返す対応).
   It supplies neither attaching maps nor boundary matrices. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Construction does not establish a cell complex, uniqueness of labels, a
chain map or a homology equivalence. Cell orientation and matrix-coordinate
compatibility must be supplied with the witness.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.4 — source-to-spine cell correspondence; this record only stores labels.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> path = joinpath(pkgdir(Topology), "..", "..", "threads", "PJ260922", "data", "261114-conversion", "cmr3-expected.json");

julia> c = H.JSON3.read(read(path, String))["cellular"];

julia> basis(r) = H.CellularBasis(String.(r["vertices"]), String.(r["edges"]), String.(r["faces"]));

julia> mat(r) = [Int(r[i][j]) for i in eachindex(r), j in eachindex(r[1])];

julia> vm(r) = Dict{Int,Union{Int,Nothing}}(parse(Int, String(k)) => v for (k,v) in pairs(r));

julia> w = H.CellularMapWitness(basis(c["source_complex_basis"]), basis(c["spine_complex_basis"]), mat(c["phi0"]), mat(c["good"]["phi1"]), mat(c["phi2"]), mat(c["source_boundary_1"]), mat(c["source_boundary_2"]), mat(c["spine_boundary_1"]), mat(c["spine_boundary_2"]), vm(c["original_to_smoothed_vertex"]), vm(c["smoothed_to_gauss_vertex"]), [String(p[1]) => String(p[2]) for p in c["ograph_edge_to_spine_edge"]], String(c["conventions"]), String(c["source_hash"]), String(c["target_hash"]));

julia> (H.cellular_map_commutes(w).commutes, H.cellular_cochain_natural(w).natural)
(true, true)
```

See also [`CellularMapWitness`](@ref), [`cellular_map_commutes`](@ref).

Role: construct
"""
struct CellularBasis
    vertices::Vector{String}
    edges::Vector{String}
    faces::Vector{String}
end

"""
    CellularMapWitness(source_complex_basis, spine_complex_basis, phi0, phi1, phi2, source_boundary_1, source_boundary_2, spine_boundary_1, spine_boundary_2, original_to_smoothed_vertex, smoothed_to_gauss_vertex, ograph_edge_to_spine_edge, conventions, source_hash, target_hash) -> CellularMapWitness

Store a candidate cellular map from a Heegaard source complex to its spine complex.

# Arguments
The two bases are `CellularBasis` records. Map and boundary matrices must
have the dimensions of [the cellular coordinate convention](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
Shapes are checked with `ArgumentError` naming the matrix and its size.
`conventions`, `source_hash` and `target_hash` must be nonempty, checked with
`ArgumentError`. Matrix entries convert to `Int`; vertex dictionaries,
edge-pair vectors and strings must convert to the declared field types.
Conversion errors propagate. Chain-map equations are not checked here.

# Returns
A witness holding both bases, the three map matrices, both pairs of boundary
matrices, two vertex correspondences, edge correspondence, convention text
and source/target identity strings.

# Guarantees
1. Matrix dimensions agree with the supplied bases. The record describes a
   candidate map in the linked cellular coordinates.
   [proved/agent: by construction]
2. Already-typed inputs are retained without copying or modification.
   [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Construction does not establish boundary-square or chain-map identities,
cellular realizability, or a homology equivalence. Identity strings and
convention text are not independently verified.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.4–6.5 — source-to-spine cellular correspondence and matrix coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> b = CellularBasis(["v"], String[], String[]); w = CellularMapWitness(b, b, ones(Int,1,1), zeros(Int,0,0), zeros(Int,0,0), zeros(Int,1,0), zeros(Int,0,0), zeros(Int,1,0), zeros(Int,0,0), Dict{Int,Union{Int,Nothing}}(), Dict{Int,Union{Int,Nothing}}(), Pair{String,String}[], "identity", "source", "target");

julia> (cellular_map_commutes(w).commutes, cellular_cochain_natural(w).natural)
(true, true)
```

See also [`cellular_map_commutes`](@ref), [`cellular_cochain_natural`](@ref).

Role: convert
"""
struct CellularMapWitness
    source_complex_basis::CellularBasis
    spine_complex_basis::CellularBasis
    phi0::Matrix{Int}
    phi1::Matrix{Int}
    phi2::Matrix{Int}
    source_boundary_1::Matrix{Int}
    source_boundary_2::Matrix{Int}
    spine_boundary_1::Matrix{Int}
    spine_boundary_2::Matrix{Int}
    original_to_smoothed_vertex::Dict{Int,Union{Int,Nothing}}
    smoothed_to_gauss_vertex::Dict{Int,Union{Int,Nothing}}
    ograph_edge_to_spine_edge::Vector{Pair{String,String}}
    conventions::String
    source_hash::String
    target_hash::String
    function CellularMapWitness(
            source_complex_basis, spine_complex_basis,
            phi0, phi1, phi2,
            source_boundary_1, source_boundary_2,
            spine_boundary_1, spine_boundary_2,
            original_to_smoothed_vertex, smoothed_to_gauss_vertex,
            ograph_edge_to_spine_edge, conventions, source_hash, target_hash)
        sv, se, sf = length(source_complex_basis.vertices), length(source_complex_basis.edges), length(source_complex_basis.faces)
        pv, pe, pf = length(spine_complex_basis.vertices), length(spine_complex_basis.edges), length(spine_complex_basis.faces)
        size(phi0) == (pv, sv) || throw(ArgumentError("phi0 has shape $(size(phi0)), expected $((pv, sv))"))
        size(phi1) == (pe, se) || throw(ArgumentError("phi1 has shape $(size(phi1)), expected $((pe, se))"))
        size(phi2) == (pf, sf) || throw(ArgumentError("phi2 has shape $(size(phi2)), expected $((pf, sf))"))
        size(source_boundary_1) == (sv, se) || throw(ArgumentError("source boundary 1 has shape $(size(source_boundary_1))"))
        size(source_boundary_2) == (se, sf) || throw(ArgumentError("source boundary 2 has shape $(size(source_boundary_2))"))
        size(spine_boundary_1) == (pv, pe) || throw(ArgumentError("spine boundary 1 has shape $(size(spine_boundary_1))"))
        size(spine_boundary_2) == (pe, pf) || throw(ArgumentError("spine boundary 2 has shape $(size(spine_boundary_2))"))
        isempty(conventions) && throw(ArgumentError("conventions are required"))
        isempty(source_hash) && throw(ArgumentError("source hash is required"))
        isempty(target_hash) && throw(ArgumentError("target hash is required"))
        new(source_complex_basis, spine_complex_basis,
            phi0, phi1, phi2,
            source_boundary_1, source_boundary_2,
            spine_boundary_1, spine_boundary_2,
            original_to_smoothed_vertex, smoothed_to_gauss_vertex,
            ograph_edge_to_spine_edge, conventions, source_hash, target_hash)
    end
end

"""
    cmr3_combed() -> CMR3State

Return the unclassified branching-combing state without a second-field equation result.

# Arguments
No arguments or geometric validation.

# Returns
`CMR3State` with kind `:combed` and absent order and diagnostic.

# Guarantees
1. The result records the initial combing state and supplies no obstruction
   classification or integer solution. [proved/agent: by construction]
2. No existing state or spine is modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
This state alone does not identify a vector field on a particular manifold
or certify that it extends to a framing. The [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the meaning of this status.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §10.3 — geometric correspondence required for additional-structure claims.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> s = cmr3_combed(); (s.kind, s.order, s.diagnostic)
(:combed, nothing, nothing)
```

See also [`classify_euler_class`](@ref).

Role: construct
"""
cmr3_combed() = CMR3State(:combed, nothing, nothing)

function _cmr3_snf(A::AbstractMatrix{<:Integer})
    rows, cols = size(A)
    exact = AbstractAlgebra.matrix(AbstractAlgebra.ZZ, [BigInt(A[i, j]) for i in 1:rows, j in 1:cols])
    smith, left, right = AbstractAlgebra.snf_with_transform(exact)
    S = [BigInt(smith[i, j]) for i in 1:rows, j in 1:cols]
    U = [BigInt(left[i, j]) for i in 1:rows, j in 1:rows]
    V = [BigInt(right[i, j]) for i in 1:cols, j in 1:cols]
    for i in 1:min(rows, cols)
        if S[i, i] < 0
            S[i, :] .*= -1
            U[i, :] .*= -1
        end
    end
    S, U, V
end

function _cmr3_class_order(S::AbstractMatrix{<:Integer}, U::AbstractMatrix{<:Integer}, c::AbstractVector{<:Integer})
    transformed = U * BigInt.(c)
    rows, cols = size(S)
    rank = count(i -> S[i, i] != 0, 1:min(rows, cols))
    any(!iszero, transformed[rank + 1:rows]) && return nothing
    order = 1
    for i in 1:rank
        divisor = S[i, i]
        order = lcm(order, Int(div(divisor, gcd(divisor, transformed[i]))))
    end
    order
end

function _cmr3_particular(S, U, V, target::AbstractVector{<:Integer})
    rows, cols = size(S)
    transformed = U * BigInt.(target)
    rank = count(i -> S[i, i] != 0, 1:min(rows, cols))
    any(!iszero, transformed[rank + 1:rows]) && return nothing
    any(i -> mod(transformed[i], S[i, i]) != 0, 1:rank) && return nothing
    coefficients = zeros(BigInt, cols)
    for i in 1:rank
        coefficients[i] = div(transformed[i], S[i, i])
    end
    V * coefficients
end

"""
    modular_witness_check(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}, witness::ModularWitness) -> Tuple{Bool,Int}

Check whether a modular functional detects the supplied integer cokernel class.

# Arguments
The witness modulus must exceed one and its vector length must equal the row count of `delta` and the length of `c_p`; these checks raise `ArgumentError` naming the mismatch. The caller supplies compatible [cellular coordinates](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates) for geometric interpretation.

# Returns
`(annihilates, residue)`, where the flag records annihilation of the matrix image modulo the witness modulus and the residue evaluates the supplied class in that quotient.

# Guarantees
1. If `annihilates` is true and `residue` is nonzero, `c_p` is not in the integer image of `delta`. [proved/agent: by construction]
2. The arguments are not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
A zero residue does not prove that the integer cokernel class vanishes. The routine does not certify that `c_p` is a geometric Euler representative or inspect other stored state.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md#cellular-basis-coordinates — cellular matrix domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> H.modular_witness_check(reshape([2],1,1), [1], H.ModularWitness(2,[1]))
(true, 1)
```

Role: condition
"""
function modular_witness_check(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}, witness::ModularWitness)
    modulus = witness.modulus
    modulus > 1 || throw(ArgumentError("modular witness modulus must exceed 1"))
    u = witness.u
    length(u) == size(delta, 1) == length(c_p) || throw(ArgumentError("witness length must match the rows of delta and c_P"))
    image = [mod(sum(BigInt(u[i]) * delta[i, j] for i in eachindex(u)), modulus) for j in 1:size(delta, 2)]
    residue = mod(sum(BigInt(u[i]) * c_p[i] for i in eachindex(u)), modulus)
    (all(iszero, image), Int(residue))
end

"""
    classify_euler_class(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}; witness=nothing) -> CMR3Reading

Classify the integer cokernel class of a supplied Euler vector for a fixed cellular matrix.

# Arguments
`delta` is an integer matrix and `c_p` has one entry per row, checked with
`ArgumentError`. Optional `witness` is a `ModularWitness`; its dimension and
modulus checks propagate. For geometric use the caller supplies the actual
cellular data and Euler representative in compatible
[cellular coordinates](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
This operation does not verify that geometric identification.

# Returns
A `CMR3Reading` with the class status, optional finite obstruction order,
optional modular checks and an integer solution when the class vanishes.

# Guarantees
1. A returned integer solution satisfies `delta * x == -c_p` exactly;
   otherwise the class is reported as nonzero, with its finite order when
   finite. [proved/agent: by construction]
2. Supplying a modular witness records its checks without changing the
   cokernel classification. [proved/agent: by construction]
3. Inputs are not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
This arithmetic classification does not identify the supplied vector with
a geometric Euler obstruction or construct a framing. The [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the compatibility obligations for that interpretation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.5, 10.3 — cellular coordinates and limits on combing interpretation.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> r = classify_euler_class(reshape([2],1,1), [1]); (r.state.kind, r.state.order, r.integer_solution)
(:euler_obstructed, 2, nothing)
```

See also [`classify_mod2_solution`](@ref), [`combing_framing_from_cmr3`](@ref).

Role: invariant
"""
function classify_euler_class(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}; witness::Union{Nothing,ModularWitness}=nothing)
    size(delta, 1) == length(c_p) || throw(ArgumentError("one Euler value per row of delta"))
    S, U, V = _cmr3_snf(delta)
    order = _cmr3_class_order(S, U, c_p)
    annihilates = nothing
    residue = nothing
    if witness !== nothing
        annihilates, residue = modular_witness_check(delta, c_p, witness)
    end
    if order === nothing
        return CMR3Reading(CMR3State(:euler_obstructed, nothing, :not_torsion), annihilates, residue, nothing)
    end
    if order == 1
        solution = _cmr3_particular(S, U, V, -BigInt.(c_p))
        solution === nothing && throw(ArgumentError("order 1 did not produce an integer solution of delta*x = -c_P"))
        BigInt.(delta) * solution == -BigInt.(c_p) || throw(ArgumentError("integer solution does not satisfy delta*x = -c_P"))
        return CMR3Reading(CMR3State(:framing_extendable, nothing, nothing), annihilates, residue, solution)
    end
    CMR3Reading(CMR3State(:euler_obstructed, order, :rational_only), annihilates, residue, nothing)
end

function _mod2_solvable(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer})
    rows, cols = size(delta)
    matrix = [mod(Int(delta[i, j]), 2) for i in 1:rows, j in 1:cols]
    target = [mod(Int(c_p[i]), 2) for i in 1:rows]
    row = 1
    for col in 1:cols
        row > rows && break
        pivot = findfirst(i -> matrix[i, col] == 1, row:rows)
        pivot === nothing && continue
        pivot += row - 1
        matrix[row, :], matrix[pivot, :] = matrix[pivot, :], matrix[row, :]
        target[row], target[pivot] = target[pivot], target[row]
        for i in 1:rows
            i == row && continue
            matrix[i, col] == 1 || continue
            matrix[i, :] .= mod.(matrix[i, :] .+ matrix[row, :], 2)
            target[i] = mod(target[i] + target[row], 2)
        end
        row += 1
    end
    all(iszero, target[row:rows])
end

"""
    classify_mod2_solution(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}) -> CMR3State

Classify solvability of the supplied cellular equation modulo two.

# Arguments
`delta` is an integer matrix and `c_p` has one entry per row, checked with
`ArgumentError`. Entries must convert to `Int`; overflow or conversion errors
propagate. The caller supplies any claimed geometric interpretation in
[compatible cellular coordinates](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).

# Returns
A `CMR3State` reporting modular solvability or its obstruction, without
an integer solution vector or a framing judgment.

# Guarantees
1. The returned kind is `:mod2_solution` exactly when the equation is
   solvable modulo two on the accepted integer domain.
   [proved/agent: by construction]
2. The inputs are not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Mod-two solvability does not establish an integer solution, a framing, or
the geometric meaning of the supplied Euler vector. The [CMR3 interpretation](../specs/heegaard-koda-conversion.md#cmr3-euler-framing)
states the obstruction convention.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.5, 10.3 — compatible cochain coordinates and geometric evidence scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> (classify_mod2_solution(reshape([1],1,1), [1]).kind, classify_mod2_solution(reshape([2],1,1), [1]).diagnostic)
(:mod2_solution, :mod2_unsolvable)
```

See also [`classify_euler_class`](@ref).

Role: condition
"""
function classify_mod2_solution(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer})
    size(delta, 1) == length(c_p) || throw(ArgumentError("one Euler value per row of delta"))
    if _mod2_solvable(delta, c_p)
        return CMR3State(:mod2_solution, nothing, nothing)
    end
    CMR3State(:euler_obstructed, nothing, :mod2_unsolvable)
end

"""
    cellular_map_commutes(witness::CellularMapWitness) -> NamedTuple

Check the degree-one and degree-two boundary identities of a candidate cellular map.

# Arguments
`witness` uses [the cellular coordinate convention](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
Its constructor checks dimensions; later container mutations remain the
caller's responsibility. Matrix dimension errors propagate. Integer entries
must be small enough to avoid `Int` overflow during products; this is not checked.

# Returns
`(degree1, degree2, commutes)` with one Boolean per identity and their conjunction.

# Guarantees
1. On the stated exact-arithmetic domain, the flags give the boundary-map
   equalities in degrees one and two of the linked convention.
   [proved/agent: by construction]
2. The witness is not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
Commutation alone does not prove that either input is a cell complex, that
the matrices are geometrically realized, or that the map is a homology equivalence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.5 — chain-map coordinates and boundary identities.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> b = CellularBasis(["v"], String[], String[]); w = CellularMapWitness(b,b,ones(Int,1,1),zeros(Int,0,0),zeros(Int,0,0),zeros(Int,1,0),zeros(Int,0,0),zeros(Int,1,0),zeros(Int,0,0),Dict{Int,Union{Int,Nothing}}(),Dict{Int,Union{Int,Nothing}}(),Pair{String,String}[],"identity","source","target"); cellular_map_commutes(w)
(degree1 = true, degree2 = true, commutes = true)
```

See also [`cellular_cochain_natural`](@ref).

Role: condition
"""
function cellular_map_commutes(witness::CellularMapWitness)
    degree1 = witness.spine_boundary_1 * witness.phi1 == witness.phi0 * witness.source_boundary_1
    degree2 = witness.spine_boundary_2 * witness.phi2 == witness.phi1 * witness.source_boundary_2
    (degree1=degree1, degree2=degree2, commutes=degree1 && degree2)
end

"""
    cellular_cochain_natural(witness::CellularMapWitness) -> NamedTuple

Check degree-two coboundary naturality on each spine-edge basis cochain.

# Arguments
`witness` uses [the cellular coordinate convention](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
Its dimensions must remain valid after construction; incompatible matrices
raise dimension errors. The caller bounds integer products to avoid unchecked
`Int` overflow.

# Returns
`(rows, natural)`, with one Boolean per spine-edge basis cochain and their
conjunction. An empty basis yields an empty row vector and `natural=true`.

# Guarantees
1. On the stated arithmetic domain, the flags test the degree-two cochain
   naturality identity in the linked coordinates.
   [proved/agent: by construction]
2. The witness is not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
This does not test degree-one chain-map commutation, geometric realizability
or a homology equivalence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.5 — transpose pullback and cochain naturality coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> b = CellularBasis(["v"], String[], String[]); w = CellularMapWitness(b,b,ones(Int,1,1),zeros(Int,0,0),zeros(Int,0,0),zeros(Int,1,0),zeros(Int,0,0),zeros(Int,1,0),zeros(Int,0,0),Dict{Int,Union{Int,Nothing}}(),Dict{Int,Union{Int,Nothing}}(),Pair{String,String}[],"identity","source","target"); r = cellular_cochain_natural(w); (isempty(r.rows), r.natural)
(true, true)
```

See also [`cellular_map_commutes`](@ref).

Role: condition
"""
function cellular_cochain_natural(witness::CellularMapWitness)
    left = witness.phi1 * witness.source_boundary_2
    right = witness.spine_boundary_2 * witness.phi2
    rows = [left[i, :] == right[i, :] for i in 1:size(left, 1)]
    (rows=rows, natural=all(rows))
end

"""
    gauss_vertex_of(witness::CellularMapWitness, original::Integer) -> Union{Int,Nothing}

Return an original Heegaard intersection's Gauss-vertex label from a cellular witness.

# Arguments
`witness` supplies the two vertex correspondences of the
[cellular conversion](../specs/heegaard-koda-conversion.md#cellular-basis-coordinates).
`original` must convert to `Int`; conversion errors propagate. Membership and
geometric correctness of the stored dictionaries are not checked.

# Returns
The corresponding Gauss label, or `nothing` when either correspondence is
absent or records a deleted vertex.

# Guarantees
1. The result is the correspondence recorded by the witness, not a fresh
   inference of a geometric vertex. [proved/agent: by construction]
2. The witness is not modified. [proved/agent: by construction]

# Conventions
CMR3 cellular and obstruction coordinates: [convention entry](../conventions/HeegaardCodes.md#computed-coordinate-conventions), [definition](../specs/heegaard-koda-conversion.md#cmr3-euler-framing).

# Not guaranteed
A stored label does not establish that the witness is a chain map or that
it identifies the original and smoothed manifolds.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §6.4 — original-to-spine correspondences.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> b = H.CellularBasis(["v"], String[], String[]); z = zeros(Int,0,0); w = H.CellularMapWitness(b,b,ones(Int,1,1),z,z,zeros(Int,1,0),z,zeros(Int,1,0),z,Dict(1=>2),Dict(2=>7),Pair{String,String}[],"test","s","t");

julia> (H.gauss_vertex_of(w,1), H.gauss_vertex_of(w,3))
(7, nothing)
```

Role: access
"""
function gauss_vertex_of(witness::CellularMapWitness, original::Integer)
    smoothed = get(witness.original_to_smoothed_vertex, Int(original), nothing)
    smoothed === nothing && return nothing
    get(witness.smoothed_to_gauss_vertex, smoothed, nothing)
end
