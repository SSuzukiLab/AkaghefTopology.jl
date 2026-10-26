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

function _occurrences(code::HeegaardCode)
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

function _syntax_issues(code::HeegaardCode)
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

function _incidence_issues(code::HeegaardCode)
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

"""Return a structured validation report without modifying `object`."""
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

"""Certify a decorated code as an immutable `HeegaardDiagram`."""
function certify_heegaard(object)
    object isa DecoratedHeegaardCode || throw(ArgumentError("certify_heegaard requires DecoratedHeegaardCode; decorate explicitly"))
    report = validate(object; level=:full)
    report.valid || throw(ValidationError(report))
    HeegaardDiagram(object, report, Val(:certified))
end

is_valid(object) = validate(object; level=:full).valid
is_closed_heegaard(object) = validate(object; level=:full).valid

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
