"""Representation-only canonicalization and SHA-256 identity."""

function _canonical_json(value)
    if value === nothing
        return "null"
    elseif value === missing
        return "null"
    elseif value isa AbstractDict
        keys_sorted = sort!(String[string(key) for key in keys(value)])
        parts = String[]
        for key in keys_sorted
            original = haskey(value, key) ? key : Symbol(key)
            push!(parts, JSON3.write(key) * ":" * _canonical_json(value[original]))
        end
        return "{" * join(parts, ",") * "}"
    elseif value isa NamedTuple
        return _canonical_json(Dict(string(name) => getfield(value, name) for name in propertynames(value)))
    elseif value isa Tuple || value isa AbstractVector
        return "[" * join((_canonical_json(item) for item in value), ",") * "]"
    elseif value isa Symbol
        return JSON3.write(String(value))
    else
        return JSON3.write(value)
    end
end

function _permutations(values::Vector{Int})
    result = Vector{Vector{Int}}()
    function visit!(prefix, remaining)
        isempty(remaining) && (push!(result, copy(prefix)); return)
        for (position, value) in enumerate(remaining)
            push!(prefix, value)
            visit!(prefix, [remaining[1:position-1]; remaining[position+1:end]])
            pop!(prefix)
        end
    end
    visit!(Int[], values)
    result
end

_rotations(values::Vector{Int}) = isempty(values) ? [Int[]] :
    [vcat(values[index:end], values[1:index-1]) for index in eachindex(values)]

function _ordered_variants(curves::Vector{Curve}, order::Vector{Int})
    variants = Vector{Vector{Vector{Int}}}()
    function visit!(position, current)
        if position > length(order)
            push!(variants, deepcopy(current))
            return
        end
        for rotation in _rotations(curves[order[position]].visits)
            push!(current, rotation)
            visit!(position + 1, current)
            pop!(current)
        end
    end
    visit!(1, Vector{Vector{Int}}())
    variants
end

function _region_signatures(object, intersection_map, curve_maps)
    object isa DecoratedHeegaardCode || return Vector{Vector{Any}}()
    graph = ribbon_graph(object.code)
    cycles = boundary_cycles(graph)
    dart_map = _dart_map(graph)
    signatures = Any[]
    for cycle in cycles
        raw = [(intersection_map[dart_map[id].origin],
                dart_map[id].origin_port,
                (dart_map[id].edge, dart_map[id].edge),
                intersection_map[dart_map[id].target],
                dart_map[id].target_port) for id in cycle]
        # Replace the edge id with a stable (family, curve) descriptor.
        for (index, id) in enumerate(cycle)
            edge = graph.edges[dart_map[id].edge]
            raw[index] = (intersection_map[dart_map[id].origin],
                          dart_map[id].origin_port,
                          (edge.family, get(curve_maps, edge.family, Dict{Int,Int}())[edge.curve_id]),
                          intersection_map[dart_map[id].target],
                          dart_map[id].target_port)
        end
        isempty(raw) || (raw = _lexmin([vcat(raw[i:end], raw[1:i-1]) for i in eachindex(raw)]))
        push!(signatures, raw)
    end
    signatures
end

function _candidate(code::HeegaardCode, alpha_order, alpha_visits, beta_order, beta_visits,
                    decorated::Bool, original_partition=nothing)
    alpha_map = Dict(code.alpha[old].id => new for (new, old) in enumerate(alpha_order))
    beta_map = Dict(code.beta[old].id => new for (new, old) in enumerate(beta_order))
    intersection_map = Dict{Int,Int}()
    next_id = 1
    remap(visits) = begin
        result = Int[]
        for point in visits
            if !haskey(intersection_map, point)
                intersection_map[point] = next_id
                next_id += 1
            end
            push!(result, intersection_map[point])
        end
        result
    end
    new_alpha = [remap(visits) for visits in alpha_visits]
    new_beta = [remap(visits) for visits in beta_visits]
    points = sort([begin
        original = findfirst(point -> point.id == old_id, code.intersections)
        point = code.intersections[original]
        (intersection_map[old_id], alpha_map[point.alpha], beta_map[point.beta], point.sign)
    end for (old_id, _) in intersection_map]; by=first)
    regions = Vector{Vector{Any}}()
    if decorated
        curve_maps = Dict(:alpha => alpha_map, :beta => beta_map)
        # Encode the actual region blocks by transformed cycle signatures,
        # making region labels and cycle numbering immaterial.
        graph = ribbon_graph(code)
        cycles = boundary_cycles(graph)
        dart_map = _dart_map(graph)
        cycle_signatures = Any[]
        for cycle in cycles
            raw = [(intersection_map[dart_map[id].origin], dart_map[id].origin_port,
                    (graph.edges[dart_map[id].edge].family,
                     curve_maps[graph.edges[dart_map[id].edge].family][graph.edges[dart_map[id].edge].curve_id]),
                    intersection_map[dart_map[id].target], dart_map[id].target_port) for id in cycle]
            push!(cycle_signatures, _lexmin([vcat(raw[i:end], raw[1:i-1]) for i in eachindex(raw)]))
        end
        original_partition === nothing && throw(ArgumentError("decorated canonicalization requires a region partition"))
        regions = sort([sort([cycle_signatures[cycle] for cycle in block]; by=repr)
                        for block in original_partition.boundary_cycle_blocks]; by=repr)
    end
    key = Dict{String,Any}(
        "alpha" => new_alpha,
        "beta" => new_beta,
        "intersections" => points,
        "regions" => regions,
    )
    (; key, alpha=new_alpha, beta=new_beta, intersections=points, regions)
end

"""Canonicalize curve/intersection labels and cyclic starting points.

For small codes this searches all family relabelings and cyclic rotations,
which makes the representation equivalence exact.  Larger codes use the
stable existing order to avoid putting canonicalization on normal compute
paths; the returned hash remains deterministic but is explicitly not a
complete graph-isomorphism solver for that case.
"""
function canonicalize_representation(object)
    code = _undecorated(object)
    decorated = object isa Union{DecoratedHeegaardCode,HeegaardDiagram}
    isempty(_syntax_issues(code)) || throw(ArgumentError("cannot canonicalize syntactically invalid code"))
    isempty(_incidence_issues(code)) || throw(ArgumentError("cannot canonicalize code with incidence errors"))
    alpha_orders = length(code.alpha) <= 6 ? _permutations(collect(eachindex(code.alpha))) : [collect(eachindex(code.alpha))]
    beta_orders = length(code.beta) <= 6 ? _permutations(collect(eachindex(code.beta))) : [collect(eachindex(code.beta))]
    estimate = sum(prod(length(_rotations(code.alpha[i].visits)) for i in order) for order in alpha_orders) *
               sum(prod(length(_rotations(code.beta[i].visits)) for i in order) for order in beta_orders)
    if estimate > 200_000
        alpha_orders = [sort(collect(eachindex(code.alpha)); by=i -> (length(code.alpha[i].visits), code.alpha[i].id))]
        beta_orders = [sort(collect(eachindex(code.beta)); by=i -> (length(code.beta[i].visits), code.beta[i].id))]
    end
    best = nothing
    best_key = nothing
    for alpha_order in alpha_orders, alpha_visits in _ordered_variants(code.alpha, alpha_order),
        beta_order in beta_orders, beta_visits in _ordered_variants(code.beta, beta_order)
        partition = object isa HeegaardDiagram ? object.code.regions :
                    object isa DecoratedHeegaardCode ? object.regions : nothing
        candidate = _candidate(code, alpha_order, alpha_visits, beta_order, beta_visits,
                              decorated, partition)
        key = _canonical_json(candidate.key)
        if best_key === nothing || key < best_key
            best_key, best = key, candidate
        end
    end
    alpha = [Curve(:alpha, i, visits) for (i, visits) in enumerate(best.alpha)]
    beta = [Curve(:beta, i, visits) for (i, visits) in enumerate(best.beta)]
    points = [Intersection(point[1], point[2], point[3], point[4]) for point in best.intersections]
    canonical_code = HeegaardCode(alpha, beta, points, deepcopy(code.metadata))
    decorated || return canonical_code

    graph = ribbon_graph(canonical_code)
    cycles = boundary_cycles(graph)
    # Recompute cycle signatures of the canonical graph and use the saved
    # signature blocks to restore the exact non-cellular partition.
    dart_map = _dart_map(graph)
    curve_maps = Dict(:alpha => Dict(i => i for i in eachindex(alpha)),
                      :beta => Dict(i => i for i in eachindex(beta)))
    signature_to_id = Dict{Any,Int}()
    for (cycle_id, cycle) in enumerate(cycles)
        raw = [(dart_map[id].origin, dart_map[id].origin_port,
                 (graph.edges[dart_map[id].edge].family,
                  curve_maps[graph.edges[dart_map[id].edge].family][graph.edges[dart_map[id].edge].curve_id]),
                 dart_map[id].target, dart_map[id].target_port) for id in cycle]
        signature_to_id[_lexmin([vcat(raw[i:end], raw[1:i-1]) for i in eachindex(raw)])] = cycle_id
    end
    blocks = [[signature_to_id[signature] for signature in block] for block in best.regions]
    DecoratedHeegaardCode(canonical_code, RegionPartition(blocks), deepcopy(code.metadata))
end

function _canonical_structural_dict(object)
    canonical = canonicalize_representation(object)
    code = _undecorated(canonical)
    result = Dict{String,Any}(
        "schema" => HGD_SCHEMA,
        "curves" => Dict("alpha" => [copy(curve.visits) for curve in code.alpha],
                         "beta" => [copy(curve.visits) for curve in code.beta]),
        "intersections" => [[point.alpha, point.beta, point.sign] for point in code.intersections],
    )
    canonical isa DecoratedHeegaardCode &&
        (result["regions"] = [sort(copy(block)) for block in canonical.regions.boundary_cycle_blocks])
    result
end

canonical_hash(object) = bytes2hex(SHA.sha256(_canonical_json(_canonical_structural_dict(object))))
