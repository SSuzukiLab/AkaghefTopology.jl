"""Exact presentation, intersection, and first-homology calculations."""

function _presentation_code(object)
    _undecorated(object)
end

function relator(object, beta_id::Integer)
    code = _presentation_code(object)
    curve = findfirst(curve -> curve.id == beta_id, code.beta)
    curve === nothing && throw(ArgumentError("unknown beta curve id $beta_id"))
    points = _intersection_map(code)
    [(points[point].alpha, points[point].sign) for point in code.beta[curve].visits]
end

function fundamental_group_presentation(object)
    code = _presentation_code(object)
    GroupPresentation(["x$(curve.id)" for curve in code.alpha],
        [relator(code, curve.id) for curve in code.beta],
        [curve.id for curve in code.beta],
        Dict{String,Any}(
            "generator_convention" => "alpha curve x_i",
            "relator_convention" => "beta cyclic word with intersection signs",
            "source_hash" => canonical_hash(object),
        ))
end

function intersection_matrix(object)
    code = _presentation_code(object)
    alpha_ids = [curve.id for curve in code.alpha]
    beta_ids = [curve.id for curve in code.beta]
    row = Dict(id => index for (index, id) in enumerate(alpha_ids))
    col = Dict(id => index for (index, id) in enumerate(beta_ids))
    result = zeros(Int, length(alpha_ids), length(beta_ids))
    for point in code.intersections
        haskey(row, point.alpha) && haskey(col, point.beta) || continue
        result[row[point.alpha], col[point.beta]] += point.sign
    end
    result
end

function abelianization_matrix(presentation::GroupPresentation)
    result = zeros(Int, length(presentation.relators), length(presentation.generators))
    for (row, relator_word) in enumerate(presentation.relators)
        for (generator, exponent) in relator_word
            1 <= generator <= size(result, 2) || throw(ArgumentError("relator references unknown generator $generator"))
            result[row, generator] += exponent
        end
    end
    result
end

function _swap_rows!(A, i, j)
    i == j && return A
    left, right = copy(A[i, :]), copy(A[j, :])
    A[i, :] = right
    A[j, :] = left
    A
end

function _swap_cols!(A, i, j)
    i == j && return A
    left, right = copy(A[:, i]), copy(A[:, j])
    A[:, i] = right
    A[:, j] = left
    A
end

"""Compute an integer Smith normal form using row and column operations."""
function smith_normal_form(matrix::AbstractMatrix{<:Integer})
    A = Int.(matrix)
    rows, cols = size(A)
    k = 1
    while k <= min(rows, cols)
        location = nothing
        for i in k:rows, j in k:cols
            if A[i, j] != 0
                location = (i, j)
                break
            end
        end
        location === nothing && break
        _swap_rows!(A, k, location[1])
        _swap_cols!(A, k, location[2])

        while true
            reduced = false
            for i in k+1:rows
                A[i, k] == 0 && continue
                q = div(A[i, k], A[k, k])
                A[i, :] .-= q .* A[k, :]
                if A[i, k] != 0 && abs(A[i, k]) < abs(A[k, k])
                    _swap_rows!(A, i, k)
                end
                reduced = true
                break
            end
            reduced && continue
            for j in k+1:cols
                A[k, j] == 0 && continue
                q = div(A[k, j], A[k, k])
                A[:, j] .-= q .* A[:, k]
                if A[k, j] != 0 && abs(A[k, j]) < abs(A[k, k])
                    _swap_cols!(A, j, k)
                end
                reduced = true
                break
            end
            reduced && continue

            bad = nothing
            for i in k+1:rows, j in k+1:cols
                if mod(A[i, j], A[k, k]) != 0
                    bad = (i, j)
                    break
                end
            end
            bad === nothing && break
            # The pivot divides the first row and column, so adding this row
            # introduces the non-divisible entry into the pivot row.
            A[k, :] .+= A[bad[1], :]
        end
        sign(A[k, k]) == -1 && (A[k, :] .*= -1)
        k += 1
    end
    A
end

function first_homology(object)
    matrix = abelianization_matrix(fundamental_group_presentation(object))
    smith = smith_normal_form(matrix)
    diagonal = [abs(smith[i, i]) for i in 1:min(size(smith)...)]
    nonzero = filter(!=(0), diagonal)
    torsion = sort([value for value in nonzero if value > 1])
    free_rank = size(matrix, 2) - length(nonzero)
    (; free_rank, torsion, smith_normal_form=smith,
       exponent_matrix=matrix,
       invariant_factors=copy(nonzero))
end

function is_homology_sphere(object)
    h = first_homology(object)
    h.free_rank == 0 && isempty(h.torsion)
end

"""Optional external group check; core operations never depend on GAP."""
function group_probe(object; backend=:gap)
    backend == :gap || throw(ArgumentError("unsupported group probe backend: $backend"))
    (; status=:unavailable,
       backend=:gap,
       reason="GAP.jl is optional and is not loaded in this core environment",
       presentation=fundamental_group_presentation(object))
end

function manifold_fingerprint(object)
    h = first_homology(object)
    (; h1_free_rank=h.free_rank,
       h1_torsion=h.torsion,
       surface_genus=surface_genus(object),
       intersection_matrix=intersection_matrix(object),
       presentation=fundamental_group_presentation(object))
end
