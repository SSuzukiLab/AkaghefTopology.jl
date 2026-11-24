# External states of one fixed combing. The integer-weight solver keeps
# `:integer_solution`, `:rational_only`, and `:inconsistent`. This file does
# not change a combing and does not search for a CP move.

import AbstractAlgebra

"""One external state of a fixed combing.

`kind` is `:combed`, `:framing_extendable`, `:euler_obstructed`, or
`:mod2_solution`. `order` is the order of `[c_P]` for `:euler_obstructed`
when that order is finite. `diagnostic` may hold the solver's internal name
`:rational_only`. A `:mod2_solution` is not a framing.
"""
struct CMR3State
    kind::Symbol
    order::Union{Nothing,Int}
    diagnostic::Union{Nothing,Symbol}
end

"""Modular witness that `[c_P]` is nonzero in the cokernel of `δ` modulo `modulus`."""
struct ModularWitness
    modulus::Int
    u::Vector{Int}
end

"""Classification of `δx = -c_P` on one fixed spine."""
struct CMR3Reading
    state::CMR3State
    witness_annihilates::Union{Nothing,Bool}
    witness_residue::Union{Nothing,Int}
    integer_solution::Union{Nothing,Vector{BigInt}}
end

struct CellularBasis
    vertices::Vector{String}
    edges::Vector{String}
    faces::Vector{String}
end

"""Chain map `Φ: C_*(X_HC) → C_*(P)`, with the two-stage vertex renumbering.

`∂(u→v) = v-u`. Commutativity is not required at construction, so a map with
one broken coefficient can be stored and rejected by `cellular_map_commutes`.
A missing vertex value is a deleted or unmapped cell, not a zero label.
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

"""Branching combing with no weight equation."""
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

"""Check `u^T δ ≡ 0` and return `(annihilates, residue)` with `residue = u^T c_P mod n`."""
function modular_witness_check(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer}, witness::ModularWitness)
    modulus = witness.modulus
    modulus > 1 || throw(ArgumentError("modular witness modulus must exceed 1"))
    u = witness.u
    length(u) == size(delta, 1) == length(c_p) || throw(ArgumentError("witness length must match the rows of delta and c_P"))
    image = [mod(sum(BigInt(u[i]) * delta[i, j] for i in eachindex(u)), modulus) for j in 1:size(delta, 2)]
    residue = mod(sum(BigInt(u[i]) * c_p[i] for i in eachindex(u)), modulus)
    (all(iszero, image), Int(residue))
end

"""Classify `[c_P]` for the fixed matrix `δ`.

An integer solution of `δx = -c_P` gives `:framing_extendable`. A nonzero class
of finite order `n` gives `:euler_obstructed` with that order and the internal
diagnostic `:rational_only`. Passing `witness` stores the two modular checks
and does not change the state.
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

"""Read the same equation over `Z/2`.

A solution is `:mod2_solution`. This function never returns `:framing_extendable`.
"""
function classify_mod2_solution(delta::AbstractMatrix{<:Integer}, c_p::AbstractVector{<:Integer})
    size(delta, 1) == length(c_p) || throw(ArgumentError("one Euler value per row of delta"))
    if _mod2_solvable(delta, c_p)
        return CMR3State(:mod2_solution, nothing, nothing)
    end
    CMR3State(:euler_obstructed, nothing, :mod2_unsolvable)
end

"""Return whether `∂^P Φ = Φ ∂^HC` in degrees 1 and 2."""
function cellular_map_commutes(witness::CellularMapWitness)
    degree1 = witness.spine_boundary_1 * witness.phi1 == witness.phi0 * witness.source_boundary_1
    degree2 = witness.spine_boundary_2 * witness.phi2 == witness.phi1 * witness.source_boundary_2
    (degree1=degree1, degree2=degree2, commutes=degree1 && degree2)
end

"""Check `δ_HC Φ1* e* = Φ2* δ_P e*` on each spine-edge basis cochain."""
function cellular_cochain_natural(witness::CellularMapWitness)
    left = witness.phi1 * witness.source_boundary_2
    right = witness.spine_boundary_2 * witness.phi2
    rows = [left[i, :] == right[i, :] for i in 1:size(left, 1)]
    (rows=rows, natural=all(rows))
end

"""Compose original intersection → smoothed vertex → Gauss vertex.

A deleted corner or a missing second-stage label returns `nothing`.
"""
function gauss_vertex_of(witness::CellularMapWitness, original::Integer)
    smoothed = get(witness.original_to_smoothed_vertex, Int(original), nothing)
    smoothed === nothing && return nothing
    get(witness.smoothed_to_gauss_vertex, smoothed, nothing)
end
