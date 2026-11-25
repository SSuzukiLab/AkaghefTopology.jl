"""Exact presentation, intersection, and first-homology calculations."""

function _presentation_code(object)
    _undecorated(object)
end

"""
    relator(object, beta_id::Integer) -> Vector{Tuple{Int,Int}}

Return the signed alpha-incidence word of one legacy beta curve.

# Arguments
`object` supplies a legacy undecorated code through `CombinatorialHeegaardCode`,
`DecoratedHeegaardCode` or `CellularHeegaardDiagram`; the caller supplies the
valid incidence of [legacy records](../specs/heegaard-terminology.md#legacy-cellular-records).
`beta_id` must identify a stored beta curve. Checked: an absent ID raises
`ArgumentError` naming it. Invalid incidence can raise lookup errors.

# Returns
A vector of alpha-ID and sign pairs, in the selected beta curve's stored
visit order, using [presentation-word encoding](../specs/heegaard-terminology.md#presentation-word-encoding).

# Guarantees
1. The result records exactly the signed alpha intersections of that beta
   curve, without word reduction. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This word alone does not certify a manifold presentation or a complete cut
system. Alpha IDs are not automatically replaced by dense generator indices.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §群表示の語 — legacy incidence and word coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); Tuple(H.relator(c, 1))
((1, 1),)
```

See also [`fundamental_group_presentation`](@ref).

Role: access
"""
function relator(object, beta_id::Integer)
    code = _presentation_code(object)
    curve = findfirst(curve -> curve.id == beta_id, code.beta)
    curve === nothing && throw(ArgumentError("unknown beta curve id $beta_id"))
    points = _intersection_map(code)
    [(points[point].alpha, points[point].sign) for point in code.beta[curve].visits]
end

"""
    fundamental_group_presentation(object) -> GroupPresentation

Return the alpha-generator and beta-relator presentation of a decorated legacy Heegaard diagram.

# Arguments
`object` is a `DecoratedHeegaardCode` or `CellularHeegaardDiagram` with the
valid incidence and decoration of [legacy records](../specs/heegaard-terminology.md#legacy-cellular-records).
Alpha completeness is checked by `compression_boundary_report`; failure
throws `ArgumentError` requiring a complete alpha cut system. Surface,
incidence, canonicalization and serialization errors propagate. The caller
supplies the remaining [Heegaard hypotheses](../specs/heegaard-koda-conversion.md#51-曲線系と多様体);
the operation does not fully validate them before returning a presentation.

# Returns
`GroupPresentation` with alpha generator names, beta relation words,
original beta IDs and convention/source metadata. Word coordinates and
IDs follow [presentation encoding](../specs/heegaard-terminology.md#presentation-word-encoding).

# Guarantees
1. The returned relations encode the signed alpha incidences along the
   stored beta visits in the legacy presentation convention.
   [proved/agent: by construction]
2. The intended geometric interpretation is the fundamental group of
   `Mhat(H)` for a valid closed Heegaard diagram with the required complete
   alpha system, as specified in
   [the manifold domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
   [unverified: general identification of this legacy presentation with the stated manifold group]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The general manifold-group identification in guarantee 2 is not certified
by this operation or by its stored source hash. Curve counts alone are not
complete-cut-system evidence. No word simplification or group classification
is performed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §群表示の語 and lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — legacy presentation coordinates and intended manifold Mhat(H).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.abelianization_matrix(H.fundamental_group_presentation(d)) == H.intersection_matrix(d)'
true
```

See also [`abelianization_matrix`](@ref), [`Topology.HeegaardCodes.first_homology`](@ref).

Role: invariant
"""
function fundamental_group_presentation(object)
    code = _presentation_code(object)
    compression_boundary_report(object).alpha_complete || throw(ArgumentError(
        "alpha generators require a complete cut system; use the full marked document for redundant systems"))
    generator_index=Dict(c.id=>i for (i,c) in enumerate(code.alpha))
    GroupPresentation(["x$(curve.id)" for curve in code.alpha],
        [[(generator_index[id],sign) for (id,sign) in relator(code,curve.id)] for curve in code.beta],
        [curve.id for curve in code.beta],
        Dict{String,Any}(
            "generator_convention" => "alpha curve x_i",
            "relator_convention" => "beta cyclic word with intersection signs",
            "source_hash" => canonical_hash(object),
        ))
end

"""
    intersection_matrix(object) -> Matrix{Int}

Return the ordered alpha–beta intersection matrix of a legacy Heegaard record.

# Arguments
`object` is a `CombinatorialHeegaardCode`, `DecoratedHeegaardCode` or
`CellularHeegaardDiagram` in the
[legacy cellular format](../specs/heegaard-terminology.md#legacy-cellular-records).
The caller guarantees consistent unique curve IDs and valid signed incidence;
these conditions are not checked here. Unknown curve references in intersection
records are omitted from the matrix. Integer sums must fit in `Int`.

# Returns
A fresh integer matrix whose rows and columns follow the stored alpha and
beta curve orders, respectively.

# Guarantees
1. For valid oriented input without integer overflow, the entries are the
   algebraic intersections of the corresponding oriented curves on the
   decorated surface. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Conventions
[HC versus ordered-intersection signs](../specs/heegaard-code.md#intersection-sign-translation)
and [the convention dictionary](../conventions/HeegaardCodes.md#comparison-conventions).

# Not guaranteed
This is not `hc_intersection` on an LC-bound HC. The two conventions
require the linked sign and coordinate translation. A returned matrix alone
does not certify a closed Heegaard diagram.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md, legacy cellular records; lib/topology.jl/docs/specs/heegaard-code.md, intersection-sign translation — ordered curve pairing and coordinate conversion.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.intersection_matrix(d) == [1;;]
true
```

See also [`hc_intersection`](@ref), [`abelianization_matrix`](@ref).

Role: invariant
"""
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

"""
    abelianization_matrix(presentation::GroupPresentation) -> Matrix{Int}

Return the integer relation matrix of the abelianization of a finite group presentation.

# Arguments
`presentation` uses the [presentation-word encoding](../specs/heegaard-terminology.md#presentation-word-encoding).
Generator references are checked: an out-of-range index throws `ArgumentError`
identifying the unknown generator. The caller must ensure that exponent sums
fit in `Int`; overflow is not checked.

# Returns
A fresh `Matrix{Int}` in the relation and generator coordinates of the
[presentation conventions](../specs/heegaard-terminology.md#presentation-word-encoding),
including zero rows for empty relations.

# Guarantees
1. With no integer overflow, the rows generate the relation subgroup in the
   free abelian group on the presentation's generators, so their quotient
   is the presented group's abelianization. [proved/agent: by construction]
2. Neither the presentation nor its word vectors are modified.
   [proved/agent: by construction]

# Not guaranteed
This operation does not simplify the group presentation or compute Smith
factors. A manifold interpretation requires an independently justified
identification of the input presentation with that manifold's group.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §群表示の語 — integer relation matrix and presentation coordinates.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> p = GroupPresentation(["x1"], [[(1,1)]], [1], Dict{String,Any}());

julia> (size(abelianization_matrix(p)), abelianization_matrix(p)[1,1])
((1, 1), 1)
```

See also [`fundamental_group_presentation`](@ref), [`smith_normal_form`](@ref), [`Topology.HeegaardCodes.first_homology`](@ref).

Role: invariant
"""
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

"""
    smith_normal_form(matrix::AbstractMatrix{<:Integer}) -> Matrix{Int}

Return the diagonal Smith form of an integer matrix using machine-integer storage.

# Arguments
Entries must convert to `Int`; out-of-range conversions raise `InexactError`. The caller must ensure intermediate integer arithmetic does not overflow. Without that bound the exact Smith-form guarantee does not apply.

# Returns
A newly allocated matrix of the same dimensions, with nonnegative Smith diagonal entries and zeros elsewhere.

# Guarantees
1. When arithmetic is exact, the returned matrix is integrally equivalent to the input, and its nonzero diagonal entries divide their successors. It presents the same integer cokernel. [proved/agent: integral elementary operations and divisibility termination]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The unimodular transformations are not returned. Arbitrary-precision arithmetic and overflow detection are not provided by this legacy helper.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 — presentation abelianization and integer homology.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> H.smith_normal_form([2 4; 6 8]) == [2 0; 0 4]
true
```

Role: other
"""
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

"""
    first_homology(object) -> NamedTuple

Return the abelianized legacy Heegaard presentation as free rank and torsion factors.

# Arguments
`object` must be in the domain of `fundamental_group_presentation`, including
the [complete alpha system](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).
Its domain and source-serialization errors propagate. The caller must ensure
integer calculations fit in `Int`; overflow is not checked.

# Returns
A named tuple with `free_rank`, sorted nontrivial `torsion` factors,
`smith_normal_form`, the relation `exponent_matrix`, and nonzero
`invariant_factors` including units.

# Guarantees
1. Subject to the integer range, the fields describe the abelianization of
   `fundamental_group_presentation(object)`.
   [proved/agent: by construction]
2. Under an established identification of that presentation with the group
   of `Mhat(H)`, this is `H1(Mhat(H); Z)`.
   [unverified: the general manifold identification of the legacy presentation]
3. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The general identification in guarantee 2 is not established by the returned
numbers. They neither classify the manifold nor prove that a homology sphere
is a sphere. This function is the legacy `HeegaardCodes` binding; resolved
HC uses `hc_first_homology`.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 and lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — first homology as abelianization and the intended manifold Mhat(H).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> h = H.first_homology(d); (h.free_rank, Tuple(h.torsion))
(0, ())
```

See also [`fundamental_group_presentation`](@ref), [`hc_first_homology`](@ref).

Role: invariant
"""
first_homology(object)=first_homology(fundamental_group_presentation(object))

"""
    first_homology(presentation::GroupPresentation) -> NamedTuple

Return the free rank and torsion factors of a finitely presented group's abelianization.

# Arguments
`presentation` uses [the word encoding](../specs/heegaard-terminology.md#presentation-word-encoding).
Invalid generator references throw `ArgumentError` naming the reference.
The caller must ensure that exponent sums and Smith calculations fit in
`Int`; overflow is not checked.

# Returns
A named tuple with `free_rank`, sorted nontrivial `torsion` factors,
`smith_normal_form`, the relation `exponent_matrix`, and nonzero
`invariant_factors` including units.

# Guarantees
1. Subject to the integer range, the presented group's abelianization is
   the direct sum of the free group of rank `free_rank` and cyclic groups
   of orders in `torsion`. [proved/agent: by construction]
2. The returned matrices have the input relation/generator dimensions,
   and the presentation is not modified. [proved/agent: by construction]

# Not guaranteed
This does not identify a supplied presentation with any manifold group,
and does not determine the nonabelian group. Arbitrary-precision safety is
not provided by the machine-integer output.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §群表示の語 — relation matrix and abelianization.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> p = H.GroupPresentation(["x1"], [fill((1,1), 5)], [1], Dict{String,Any}());

julia> h = H.first_homology(p); (h.free_rank, Tuple(h.torsion))
(0, (5,))
```

See also [`abelianization_matrix`](@ref), [`smith_normal_form`](@ref).

Role: invariant
"""
function first_homology(presentation::GroupPresentation)
    matrix = abelianization_matrix(presentation)
    smith = smith_normal_form(matrix)
    diagonal = [abs(smith[i, i]) for i in 1:min(size(smith)...)]
    nonzero = filter(!=(0), diagonal)
    torsion = sort([value for value in nonzero if value > 1])
    free_rank = size(matrix, 2) - length(nonzero)
    (; free_rank, torsion, smith_normal_form=smith,
       exponent_matrix=matrix,
       invariant_factors=copy(nonzero))
end

"""
    is_homology_sphere(object) -> Bool

Test whether a legacy presentation has trivial abelianization.

# Arguments
`object` is accepted by `Topology.HeegaardCodes.first_homology`, including
a `GroupPresentation` or a legacy Heegaard record. Its domain checks and
errors propagate; the caller supplies consistent presentation data and sums
within `Int`. For a manifold interpretation the caller must separately
establish the closed connected orientable 3-manifold hypotheses of
[§5.1](../specs/heegaard-koda-conversion.md#51-曲線系と多様体).

# Returns
`Bool`: whether the computed first-homology summary has free rank zero and
no nontrivial torsion factors.

# Guarantees
1. Within the algebraic domain of `first_homology`, the result tests whether
   the presented group's abelianization is trivial.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The function does not validate a closed manifold or identify the presented
group with its fundamental group. Without that identification, `true` is
only a statement about the presentation. It does not recognize the 3-sphere.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1; lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — homology interpretation requires the diagram's manifold domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.is_homology_sphere(d)
true
```

See also [`Topology.HeegaardCodes.first_homology`](@ref), [`is_closed_heegaard`](@ref).

Role: condition
"""
function is_homology_sphere(object)
    h = first_homology(object)
    h.free_rank == 0 && isempty(h.torsion)
end

"""
    group_probe(object; backend=:gap) -> NamedTuple

Return the legacy fundamental-group presentation with an unavailable external-probe status.

# Arguments
`object` must be accepted by `fundamental_group_presentation`, with the
[legacy incidence and region preconditions](../specs/heegaard-terminology.md#legacy-cellular-records).
Those errors propagate. Checked: a backend other than `:gap` raises
`ArgumentError` naming the unsupported backend.

# Returns
A named tuple with `status=:unavailable`, `backend=:gap`, a human-readable
`reason`, and `presentation::GroupPresentation`.

# Guarantees
1. The presentation is the legacy presentation associated with the supplied
   diagram in the scope of its presentation contract. [proved/agent: by construction]
2. The input is not modified and no external GAP calculation is performed.
   [proved/agent: by construction]

# Not guaranteed
The result does not test group isomorphism, finiteness or triviality. Backend
availability is not probed dynamically; this entry reports unavailable.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md#legacy-cellular-records — scope of legacy presentation words.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); p = H.group_probe(H.decorate(c, [[1]]));

julia> (p.status, p.backend)
(:unavailable, :gap)
```

Role: other
"""
function group_probe(object; backend=:gap)
    backend == :gap || throw(ArgumentError("unsupported group probe backend: $backend"))
    (; status=:unavailable,
       backend=:gap,
       reason="GAP.jl is optional and is not loaded in this core environment",
       presentation=fundamental_group_presentation(object))
end

"""
    manifold_fingerprint(object) -> NamedTuple

Collect the legacy presentation's homology and the decorated surface's algebraic data.

# Arguments
`object` must be a decorated or certified legacy Heegaard record accepted
by `first_homology`, `surface_genus`, `intersection_matrix` and
`fundamental_group_presentation`. Their checks and exceptions propagate.
The caller supplies valid incidence and planar region data in
[the legacy format](../specs/heegaard-terminology.md#legacy-cellular-records),
and machine-integer computations without overflow.

# Returns
A named tuple with `h1_free_rank`, `h1_torsion`, `surface_genus`,
`intersection_matrix` and `presentation`, carrying the corresponding
component results without additional comparison or identification.

# Guarantees
1. The homology fields describe the abelianization of the returned
   presentation on the stated algebraic domain.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This tuple is not a complete manifold invariant: surface genus and the
chosen presentation depend on the diagram. Equality does not prove
homeomorphism, and the general identification of the legacy presentation
with a manifold group is not established by collecting these fields.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9.1 and legacy cellular records — presentation, abelianization and decorated surface data.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> f = H.manifold_fingerprint(d); (f.surface_genus, f.h1_free_rank, Tuple(f.h1_torsion))
(1, 0, ())
```

See also [`Topology.HeegaardCodes.first_homology`](@ref), [`fundamental_group_presentation`](@ref).

Role: invariant
"""
function manifold_fingerprint(object)
    h = first_homology(object)
    (; h1_free_rank=h.free_rank,
       h1_torsion=h.torsion,
       surface_genus=surface_genus(object),
       intersection_matrix=intersection_matrix(object),
       presentation=fundamental_group_presentation(object))
end
