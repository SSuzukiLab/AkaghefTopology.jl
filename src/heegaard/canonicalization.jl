"""Representation-only canonicalization and SHA-256 identity."""

function _write_canonical_json(io::IO, value)
    if value === nothing
        write(io, "null")
    elseif value === missing
        write(io, "null")
    elseif value isa AbstractDict
        keys_sorted = sort!(String[string(key) for key in keys(value)])
        write(io, '{')
        for (index, key) in enumerate(keys_sorted)
            index > 1 && write(io, ',')
            original = haskey(value, key) ? key : Symbol(key)
            write(io, JSON3.write(key), ':')
            _write_canonical_json(io, value[original])
        end
        write(io, '}')
    elseif value isa NamedTuple
        names = sort!(collect(propertynames(value)); by=string)
        write(io, '{')
        for (index, name) in enumerate(names)
            index > 1 && write(io, ',')
            write(io, JSON3.write(string(name)), ':')
            _write_canonical_json(io, getfield(value, name))
        end
        write(io, '}')
    elseif value isa Tuple || value isa AbstractVector
        write(io, '[')
        for (index, item) in enumerate(value)
            index > 1 && write(io, ',')
            _write_canonical_json(io, item)
        end
        write(io, ']')
    elseif value isa Symbol
        write(io, JSON3.write(String(value)))
    else
        write(io, JSON3.write(value))
    end
    io
end

function _canonical_json(value)
    io = IOBuffer()
    _write_canonical_json(io, value)
    String(take!(io))
end

function _cyclic_lexmin(values::AbstractVector)
    isempty(values) && return collect(values)
    keys = repr.(values)
    best = 1
    for candidate in 2:length(values)
        for offset in 0:length(values)-1
            candidate_key = keys[mod1(candidate + offset, length(values))]
            best_key = keys[mod1(best + offset, length(values))]
            candidate_key == best_key && continue
            candidate_key < best_key && (best = candidate)
            break
        end
    end
    vcat(values[best:end], values[1:best-1])
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
        isempty(raw) || (raw = _cyclic_lexmin(raw))
        push!(signatures, raw)
    end
    signatures
end

function _candidate(code::CombinatorialHeegaardCode, alpha_order, alpha_visits, beta_order, beta_visits,
                    decorated::Bool, original_partition=nothing, ribbon_data=nothing)
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
        ribbon_data === nothing && throw(ArgumentError("decorated canonicalization requires ribbon data"))
        graph, cycles, dart_map = ribbon_data
        cycle_signatures = Any[]
        for cycle in cycles
            raw = [(intersection_map[dart_map[id].origin], dart_map[id].origin_port,
                    (graph.edges[dart_map[id].edge].family,
                     curve_maps[graph.edges[dart_map[id].edge].family][graph.edges[dart_map[id].edge].curve_id]),
                    intersection_map[dart_map[id].target], dart_map[id].target_port) for id in cycle]
            push!(cycle_signatures, _cyclic_lexmin(raw))
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

"""
    canonicalize_representation(object) -> Union{CombinatorialHeegaardCode,DecoratedHeegaardCode}

Return a deterministic relabelled legacy representation with normalized cyclic starts.

# Arguments
`object` is a legacy combinatorial, decorated or certified diagram under
[legacy records](../specs/heegaard-terminology.md#legacy-cellular-records).
Syntax and incidence are checked; failures raise `ArgumentError` describing
the invalid code. Decorated inputs require a consistent boundary partition;
ribbon and correspondence errors propagate. Full Heegaard validity is not checked.

# Returns
A `CombinatorialHeegaardCode` for an undecorated input, otherwise a
`DecoratedHeegaardCode`, with relabelled curves/intersections and corresponding
region partition. A certified input does not retain its certification wrapper.

# Guarantees
1. On consistent input, the returned representation differs by family-preserving
   relabelling and cyclic starts, with the supplied region correspondence
   transported. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is not a complete representation-isomorphism solver on every size and
does not decide 3-manifold homeomorphism. Large inputs use bounded choices;
curve labels/tags and the certified wrapper are not preserved. The complete
canonical-order and cutoff convention needs an owning specification entry.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9 — representation identity versus manifold equivalence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); n = H.canonicalize_representation(c); (length(n.alpha), length(n.beta), length(n.intersections))
(1, 1, 1)
```

See also [`Topology.HeegaardCodes.canonical_hash`](@ref).

Role: move
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
    ribbon_data = if decorated
        graph = ribbon_graph(code)
        cycles = boundary_cycles(graph)
        (graph, cycles, _dart_map(graph))
    else
        nothing
    end
    for alpha_order in alpha_orders, alpha_visits in _ordered_variants(code.alpha, alpha_order),
        beta_order in beta_orders, beta_visits in _ordered_variants(code.beta, beta_order)
        partition = object isa HeegaardDiagram ? object.code.regions :
                    object isa DecoratedHeegaardCode ? object.regions : nothing
        candidate = _candidate(code, alpha_order, alpha_visits, beta_order, beta_visits,
                              decorated, partition, ribbon_data)
        key = _canonical_json(candidate.key)
        if best_key === nothing || key < best_key
            best_key, best = key, candidate
        end
    end
    alpha = [Curve(:alpha, i, visits) for (i, visits) in enumerate(best.alpha)]
    beta = [Curve(:beta, i, visits) for (i, visits) in enumerate(best.beta)]
    points = [Intersection(point[1], point[2], point[3], point[4]) for point in best.intersections]
    canonical_code = CombinatorialHeegaardCode(alpha, beta, points, deepcopy(code.metadata))
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
        signature_to_id[_cyclic_lexmin(raw)] = cycle_id
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

# Rendering/validation repeatedly asks for the same canonical identity. Key the
# bounded process-local cache by serialized input CONTENT, never object identity:
# curves/regions are mutable, and an edited input must not reuse an old hash.
# This preserves `_canonical_structural_dict` and every existing hash value.
const _canonical_hash_cache = Dict{String,String}()
const _canonical_hash_lock = ReentrantLock()
"""
    canonical_hash(object) -> String

Return a structural fingerprint of a legacy cellular Heegaard representation.

# Arguments
`object` is a `CombinatorialHeegaardCode`, `DecoratedHeegaardCode` or
`CellularHeegaardDiagram` with consistent
[legacy records](../specs/heegaard-terminology.md#legacy-cellular-records).
Syntax and incidence failures throw `ArgumentError` identifying the failed
class of checks. Decorations must reference existing boundary cycles;
invalid references can raise indexing or lookup errors. Metadata must be
serializable by the JSON interface, whose errors propagate. A resolved HC
is not accepted by this legacy operation.

# Returns
A lowercase hexadecimal SHA-256 string for the normalized structural data.

# Guarantees
1. For unchanged data, the result is deterministic and ignores display
   labels, source paths and metadata in the fingerprinted structure.
   [proved/agent: by construction]
2. The registered relabeled examples and reordered region blocks have equal
   fingerprints. [verified-computationally: test/heegaard/core.jl, Representation canonicalization, two-curve signed examples]
3. The input and its containers are not modified.
   [proved/agent: by construction]

# Not guaranteed
This is not a manifold invariant or a collision-free identity certificate.
For large codes, representation normalization limits the considered curve
orders; equivalent representations can therefore have different hashes.
Neither equal nor different hashes alone decide homeomorphism. Full surface
and Heegaard validity is not checked.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-terminology.md §9 — representation identity is distinct from moves and manifold comparison.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1], [2]], [[1], [2]], [H.Intersection(1,1,1,1), H.Intersection(2,2,2,-1)]);

julia> H.canonical_hash(H.decorate(c, [[2], [1]])) == H.canonical_hash(H.decorate(c, [[1], [2]]))
true
```

See also [`canonicalize_representation`](@ref), [`hc_canonical`](@ref).

Role: invariant
"""
function canonical_hash(object)
    key = _canonical_json(_json_dict(object))
    lock(_canonical_hash_lock) do
        haskey(_canonical_hash_cache,key) && return _canonical_hash_cache[key]
        result = bytes2hex(SHA.sha256(_canonical_json(_canonical_structural_dict(object))))
        length(_canonical_hash_cache) >= 32 && empty!(_canonical_hash_cache)
        _canonical_hash_cache[key] = result
        result
    end
end
