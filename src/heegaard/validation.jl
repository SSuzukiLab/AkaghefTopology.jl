"""Layered, non-mutating validation and certification."""

const _VALIDATION_LEVELS = (:syntax, :incidence, :ribbon, :decoration, :surface, :closed_heegaard)

function _validation_level(level)
    level == :full && return length(_VALIDATION_LEVELS) + 1
    level in _VALIDATION_LEVELS || throw(ArgumentError("unknown validation level: $level"))
    findfirst(==(level), _VALIDATION_LEVELS)
end

function _add_issue!(issues, level, code, message, path="")
    push!(issues, ValidationIssue(level, code, message, path))
end

function _occurrences(code::CombinatorialHeegaardCode)
    alpha = Dict{Int,Vector{Tuple{Int,Int}}}()
    beta = Dict{Int,Vector{Tuple{Int,Int}}}()
    for curve in code.alpha
        for (position, point) in enumerate(curve.visits)
            push!(get!(alpha, point, Tuple{Int,Int}[]), (curve.id, position))
        end
    end
    for curve in code.beta
        for (position, point) in enumerate(curve.visits)
            push!(get!(beta, point, Tuple{Int,Int}[]), (curve.id, position))
        end
    end
    alpha, beta
end

function _syntax_issues(code::CombinatorialHeegaardCode)
    issues = ValidationIssue[]
    for family in (:alpha, :beta)
        curves = family == :alpha ? code.alpha : code.beta
        ids = [curve.id for curve in curves]
        length(ids) == length(unique(ids)) || _add_issue!(issues, :syntax, :duplicate_curve_id,
            "duplicate $family curve id", "curves.$family")
        for (index, curve) in enumerate(curves)
            for (position, point) in enumerate(curve.visits)
                point > 0 || _add_issue!(issues, :syntax, :invalid_visit,
                    "intersection visit must be a positive id", "curves.$family[$index].visits[$position]")
            end
        end
    end
    point_ids = [point.id for point in code.intersections]
    length(point_ids) == length(unique(point_ids)) || _add_issue!(issues, :syntax, :duplicate_intersection_id,
        "duplicate intersection id", "intersections")
    alpha_ids = Set(curve.id for curve in code.alpha)
    beta_ids = Set(curve.id for curve in code.beta)
    for (index, point) in enumerate(code.intersections)
        point.sign in (-1, 1) || _add_issue!(issues, :syntax, :invalid_sign,
            "intersection sign must be +1 or -1", "intersections[$index].sign")
        point.alpha in alpha_ids || _add_issue!(issues, :syntax, :dangling_alpha,
            "intersection references unknown alpha curve", "intersections[$index].alpha")
        point.beta in beta_ids || _add_issue!(issues, :syntax, :dangling_beta,
            "intersection references unknown beta curve", "intersections[$index].beta")
    end
    issues
end

function _incidence_issues(code::CombinatorialHeegaardCode)
    issues = ValidationIssue[]
    alpha, beta = _occurrences(code)
    point_map = _intersection_map(code)
    for point in code.intersections
        length(get(alpha, point.id, Tuple{Int,Int}[])) == 1 ||
            _add_issue!(issues, :incidence, :alpha_occurrence,
                "intersection $(point.id) must occur exactly once on alpha curve $(point.alpha)", "intersections.$(point.id)")
        length(get(beta, point.id, Tuple{Int,Int}[])) == 1 ||
            _add_issue!(issues, :incidence, :beta_occurrence,
                "intersection $(point.id) must occur exactly once on beta curve $(point.beta)", "intersections.$(point.id)")
        for occurrence in get(alpha, point.id, Tuple{Int,Int}[])
            occurrence[1] == point.alpha || _add_issue!(issues, :incidence, :alpha_owner,
                "intersection $(point.id) is listed on the wrong alpha curve", "curves.alpha")
        end
        for occurrence in get(beta, point.id, Tuple{Int,Int}[])
            occurrence[1] == point.beta || _add_issue!(issues, :incidence, :beta_owner,
                "intersection $(point.id) is listed on the wrong beta curve", "curves.beta")
        end
    end
    for point in setdiff(union(keys(alpha), keys(beta)), keys(point_map))
        _add_issue!(issues, :incidence, :dangling_visit,
            "curve visits unknown intersection $point", "curves")
    end
    issues
end

function _decoration_issues(object::DecoratedHeegaardCode, cycles)
    issues = ValidationIssue[]
    blocks = object.regions.boundary_cycle_blocks
    flattened = reduce(vcat, blocks; init=Int[])
    expected = collect(eachindex(cycles))
    sort(flattened) == expected || _add_issue!(issues, :decoration, :incomplete_partition,
        "region partition must contain each boundary cycle exactly once", "regions.boundary_cycle_blocks")
    length(flattened) == length(unique(flattened)) || _add_issue!(issues, :decoration, :overlapping_partition,
        "region partition blocks must be pairwise disjoint", "regions.boundary_cycle_blocks")
    for (index, block) in enumerate(blocks)
        isempty(block) && _add_issue!(issues, :decoration, :empty_region,
            "region block cannot be empty", "regions.boundary_cycle_blocks[$index]")
        for cycle in block
            1 <= cycle <= length(cycles) || _add_issue!(issues, :decoration, :dangling_cycle,
                "region references unknown boundary cycle $cycle", "regions.boundary_cycle_blocks[$index]")
        end
    end
    issues
end

function _validate_code(object)
    object isa HeegaardDiagram && return object.code
    object
end

"""
    validate(object; level=:full) -> ValidationReport

Report the requested consistency and surface checks of a legacy Heegaard record.

# Arguments
`object` is a legacy combinatorial code, decorated code or cellular diagram. `level` is `:syntax`, `:incidence`, `:ribbon`, `:decoration`, `:surface`, `:closed_heegaard` or `:full`; unknown levels raise `ArgumentError`. Checks use [legacy records](../specs/heegaard-terminology.md#legacy-cellular-records) and [the compression domain](../specs/heegaard-koda-conversion.md#51-曲線系と多様体). Unsupported objects and malformed field values can throw.

# Returns
A `ValidationReport` containing the requested level, overall `valid` flag, issue records and per-stage check flags. Failed prerequisites can produce blocked-stage issues.

# Guarantees
1. `valid` means no issue was recorded among the requested checks; full validation includes decoration, reconstructed surface and spherical compression-boundary checks. [proved/agent: by construction]
2. The object is not modified. [proved/agent: by construction]

# Not guaranteed
A partial-level success does not imply full validity. Even a full success is not a general manifold-equivalence proof, a joining-disk certificate or a Koda-conversion guarantee.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — surface and compression conditions.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); (H.validate(c;level=:syntax).valid, H.validate(c).valid)
(true, false)
```

Role: condition
"""
function validate(object; level=:full)
    requested = _validation_level(level)
    source = _validate_code(object)
    code = _undecorated(source)
    issues = ValidationIssue[]
    checks = Dict{Symbol,Bool}()

    syntax = _syntax_issues(code)
    append!(issues, syntax)
    checks[:syntax] = isempty(syntax)
    requested == 1 && return ValidationReport(level, isempty(issues), issues, checks)

    incidence = checks[:syntax] ? _incidence_issues(code) : ValidationIssue[]
    append!(issues, incidence)
    checks[:incidence] = isempty(incidence) && checks[:syntax]
    requested == 2 && return ValidationReport(level, isempty(issues), issues, checks)

    graph = nothing
    cycles = Vector{Vector{Int}}()
    ribbon = ValidationIssue[]
    if checks[:incidence]
        try
            graph = ribbon_graph(code)
            cycles = boundary_cycles(graph)
            length(unique(vcat(cycles...))) == length(graph.darts) ||
                _add_issue!(ribbon, :ribbon, :dart_coverage, "boundary cycles do not cover every dart", "ribbon")
        catch err
            _add_issue!(ribbon, :ribbon, :construction_failure, sprint(showerror, err), "ribbon")
        end
    else
        _add_issue!(ribbon, :ribbon, :blocked, "ribbon validation blocked by earlier failures", "ribbon")
    end
    append!(issues, ribbon)
    checks[:ribbon] = isempty(ribbon)
    requested == 3 && return ValidationReport(level, isempty(issues), issues, checks)

    decoration = ValidationIssue[]
    if source isa DecoratedHeegaardCode && checks[:ribbon]
        append!(decoration, _decoration_issues(source, cycles))
    elseif !(source isa DecoratedHeegaardCode)
        _add_issue!(decoration, :decoration, :missing_decoration,
            "full certification requires explicit region decoration", "regions")
    else
        _add_issue!(decoration, :decoration, :blocked, "decoration validation blocked by ribbon failures", "regions")
    end
    append!(issues, decoration)
    checks[:decoration] = isempty(decoration)
    requested == 4 && return ValidationReport(level, isempty(issues), issues, checks)

    surface = ValidationIssue[]
    if checks[:ribbon] && checks[:decoration]
        try
            genus = surface_genus(source)
            genus >= 0 || _add_issue!(surface, :surface, :negative_genus, "reconstructed genus is negative", "surface")
            isinteger(genus) || _add_issue!(surface, :surface, :nonintegral_genus, "reconstructed genus is not integral", "surface")
        catch err
            _add_issue!(surface, :surface, :construction_failure, sprint(showerror, err), "surface")
        end
    else
        _add_issue!(surface, :surface, :blocked, "surface validation blocked by earlier failures", "surface")
    end
    append!(issues, surface)
    checks[:surface] = isempty(surface)
    requested == 5 && return ValidationReport(level, isempty(issues), issues, checks)

    closed = ValidationIssue[]
    if checks[:surface]
        try
            report = compression_boundary_report(source)
            report.closed || _add_issue!(closed, :closed_heegaard, :not_closed,
                "alpha and beta systems do not define genus-zero compression boundaries", "closedness")
        catch err
            _add_issue!(closed,:closed_heegaard,:compression_failure,sprint(showerror,err),"closedness")
        end
    else
        _add_issue!(closed, :closed_heegaard, :blocked, "closedness validation blocked by surface failures", "closedness")
    end
    append!(issues, closed)
    checks[:closed_heegaard] = isempty(closed)
    ValidationReport(level, isempty(issues), issues, checks)
end

"""
    certify_heegaard(object) -> CellularHeegaardDiagram

Wrap a decorated legacy diagram with its successful full validation report.

# Arguments
`object` must be a `DecoratedHeegaardCode`; other values raise `ArgumentError`
requiring explicit decoration. Full validation checks the applicable
[legacy surface and Heegaard conditions](../specs/heegaard-terminology.md#legacy-cellular-records).
A failing report raises `ValidationError` carrying its issues; malformed
low-level inputs may propagate their validation errors.

# Returns
A `CellularHeegaardDiagram` retaining the supplied decorated code and the
successful `ValidationReport`.

# Guarantees
1. At return time the input passed `validate(object; level=:full)`.
   [proved/agent: by construction]
2. The input is not modified and is retained by reference.
   [proved/agent: by construction]

# Not guaranteed
Shared nested containers remain mutable; later edits can invalidate the
stored report. Certification does not identify the represented 3-manifold
or validate a later Koda conversion.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §legacy-cellular-records and lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — legacy validation and the Heegaard domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]); H.certify_heegaard(d) isa H.CellularHeegaardDiagram
true
```

See also [`validate`](@ref), [`joining_disk`](@ref).

Role: convert
"""
function certify_heegaard(object)
    object isa DecoratedHeegaardCode || throw(ArgumentError("certify_heegaard requires DecoratedHeegaardCode; decorate explicitly"))
    report = validate(object; level=:full)
    report.valid || throw(ValidationError(report))
    HeegaardDiagram(object, report, Val(:certified))
end

"""
    is_valid(object) -> Bool

Test full validity of a legacy decorated Heegaard record, including closedness.

# Arguments
`object` is a legacy combinatorial, decorated or certified Heegaard record
of [the cellular domain](../specs/heegaard-terminology.md#legacy-cellular-records).
Unsupported object types raise `ArgumentError` in `validate`. Recorded
validation failures return `false`; full validation requires explicit region
decoration, so an undecorated code returns `false`.

# Returns
The `valid` flag of `validate(object; level=:full)`.

# Guarantees
1. `true` means the implemented syntax, incidence, ribbon, decoration,
   surface and compression-boundary checks all passed on this record.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not provide an independent embedding certificate or require a
minimal complete cutting system. Use `compression_boundary_report` for its
separate completeness flags; closed compression boundaries alone do not imply
those flags.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1; lib/topology.jl/docs/specs/heegaard-terminology.md, legacy cellular records — closedness and the supported decorated surface domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.is_valid(d)
true
```

See also [`validate`](@ref), [`is_closed_heegaard`](@ref), [`compression_boundary_report`](@ref).

Role: condition
"""
is_valid(object) = validate(object; level=:full).valid

"""
    is_closed_heegaard(object) -> Bool

Test full legacy Heegaard validity with spherical compression boundaries.

# Arguments
`object` is a legacy record accepted by `validate`, in the
[cellular domain](../specs/heegaard-terminology.md#legacy-cellular-records).
Unsupported types raise `ArgumentError`. Failure of any full-validation
check, including missing region decoration, returns `false`.

# Returns
`Bool`, equal to `is_valid(object)` for this legacy API.

# Guarantees
1. `true` means full validation found no failure, including the test of
   spherical alpha and beta compression boundaries.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is the legacy full-validity predicate, not a test on arbitrary LC-bound
HC data. It does not identify a closed 3-manifold, certify minimal cut systems
or replace the separate completeness checks needed by Koda conversion.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §5.1 — closed compression bodies and the separate complete-system domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> H.is_closed_heegaard(d)
true
```

See also [`validate`](@ref), [`is_valid`](@ref), [`compression_boundary_report`](@ref).

Role: condition
"""
is_closed_heegaard(object) = validate(object; level=:full).valid

"""
    repair_suggestions(report::ValidationReport) -> Vector{String}

Return advisory repair messages for recognized legacy validation issues.

# Arguments
`report` is a validation report. No independent verification of its issue
codes or validity flag is performed; the caller supplies a meaningful report.

# Returns
A vector of distinct messages, in first-occurrence order, for recognized
decoration, incidence and boundary-partition problems.

# Guarantees
1. Each returned message addresses an issue code present in the supplied
   report; unsupported issue codes contribute no message.
   [proved/agent: by construction]
2. The report and the diagram it describes are not modified.
   [proved/agent: by construction]

# Not guaranteed
Suggestions are not repairs or a certificate that a repair exists. An empty
result does not imply a valid diagram.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> r = H.ValidationReport(:full, false, H.ValidationIssue[], Dict{Symbol,Bool}()); isempty(H.repair_suggestions(r))
true
```

See also [`validate`](@ref).

Role: access
"""
function repair_suggestions(report::ValidationReport)
    suggestions = String[]
    for issue in report.issues
        if issue.code == :missing_decoration
            push!(suggestions, "call decorate(code, regions) or pass cellular=true after checking the hypothesis")
        elseif issue.code in (:alpha_occurrence, :beta_occurrence, :dangling_visit)
            push!(suggestions, "reconcile intersection table and cyclic visits before certification")
        elseif issue.code == :incomplete_partition
            push!(suggestions, "partition all computed boundary cycle ids exactly once")
        end
    end
    unique(suggestions)
end
