"""Pure Heegaard moves with explicit certificates."""

function _invariant_snapshot(object)
    h = first_homology(object)
    Dict{String,Any}(
        "surface_genus" => surface_genus(object),
        "h1_free_rank" => h.free_rank,
        "h1_torsion" => copy(h.torsion),
        "intersection_matrix" => intersection_matrix(object),
    )
end

function _resolve_curve(code::HeegaardCode, selector)
    family = _spec_field(selector, :family, nothing)
    id = _spec_field(selector, :id, selector)
    if selector isa Tuple && length(selector) == 2
        family, id = selector
    end
    family !== nothing && (family = Symbol(family))
    if family === nothing
        alpha = findall(curve -> curve.id == id, code.alpha)
        beta = findall(curve -> curve.id == id, code.beta)
        length(alpha) + length(beta) == 1 || throw(ArgumentError("curve selector is ambiguous; specify family and id"))
        return isempty(alpha) ? (:beta, only(beta)) : (:alpha, only(alpha))
    end
    family in (:alpha, :beta) || throw(ArgumentError("curve family must be alpha or beta"))
    curves = family == :alpha ? code.alpha : code.beta
    index = findfirst(curve -> curve.id == id, curves)
    index === nothing && throw(ArgumentError("unknown $family curve id $id"))
    family, index
end

function _replace_curves(code::HeegaardCode, alpha, beta, points)
    HeegaardCode(copy(alpha), copy(beta), copy(points), deepcopy(code.metadata))
end

function _new_partition(object, new_code; region=1, append_cycles=true)
    old_partition = object.code.regions
    old_cycles = length(boundary_cycles(object))
    new_cycles = length(boundary_cycles(new_code))
    blocks = [copy(block) for block in old_partition.boundary_cycle_blocks]
    if isempty(blocks) && new_cycles > 0
        blocks = [Int[]]
    end
    if append_cycles && new_cycles > old_cycles
        1 <= region <= length(blocks) || throw(ArgumentError("region index is out of range"))
        append!(blocks[region], old_cycles+1:new_cycles)
    end
    if new_cycles < old_cycles
        for block in blocks
            filter!(cycle -> cycle <= new_cycles, block)
        end
        filter!(!isempty, blocks)
    end
    if new_cycles > 0 && isempty(blocks)
        blocks = [[i] for i in 1:new_cycles]
    end
    decorate(new_code, blocks; metadata=object.code.metadata)
end

function _finish_move(H::HeegaardDiagram, decorated, move_type, attachment, inverse)
    result = certify_heegaard(decorated)
    certificate = MoveCertificate(move_type, _copy_meta(attachment), _copy_meta(Dict("identity" => true)),
        canonical_hash(H), canonical_hash(result), _copy_meta(inverse), _invariant_snapshot(H))
    result, certificate
end

"""Reverse a curve orientation and all incident intersection signs."""
function reverse_curve_orientation(H::HeegaardDiagram, curve)
    code = H.code.code
    family, index = _resolve_curve(code, curve)
    selected = family == :alpha ? code.alpha[index] : code.beta[index]
    curves_alpha = copy(code.alpha)
    curves_beta = copy(code.beta)
    reversed = Curve(family, selected.id, reverse(selected.visits);
                     label=selected.label, tags=selected.tags)
    family == :alpha ? (curves_alpha[index] = reversed) : (curves_beta[index] = reversed)
    points = [Intersection(point.id, point.alpha, point.beta,
                           (point.alpha == selected.id || point.beta == selected.id) ? -point.sign : point.sign;
                           source_label=point.source_label) for point in code.intersections]
    new_code = _replace_curves(code, curves_alpha, curves_beta, points)
    decorated = DecoratedHeegaardCode(new_code, H.code.regions, deepcopy(H.code.metadata))
    _finish_move(H, decorated, :epsilon,
        Dict("family" => String(family), "curve_id" => selected.id),
        Dict("move" => "reverse_curve_orientation", "family" => String(family), "curve_id" => selected.id))
end

function _insert_after(values::Vector{Int}, position::Integer, value::Int)
    isempty(values) && (position in (0, 1) || throw(ArgumentError("empty curve insertion position must be 0 or 1")))
    !isempty(values) && (1 <= position <= length(values) || throw(ArgumentError("curve insertion position is out of range")))
    isempty(values) ? [value] : [values[1:position]; value; values[position+1:end]]
end

"""Add a signed canceling pair at explicit alpha/beta segment positions."""
function add_canceling_pair(H::HeegaardDiagram, location)
    code = H.code.code
    alpha_id = _spec_field(location, :alpha_curve, _spec_field(location, :alpha, 0))
    beta_id = _spec_field(location, :beta_curve, _spec_field(location, :beta, 0))
    alpha_id > 0 && beta_id > 0 || throw(ArgumentError("R move requires alpha_curve and beta_curve"))
    ai = findfirst(curve -> curve.id == alpha_id, code.alpha)
    bi = findfirst(curve -> curve.id == beta_id, code.beta)
    ai === nothing && throw(ArgumentError("unknown alpha curve id $alpha_id"))
    bi === nothing && throw(ArgumentError("unknown beta curve id $beta_id"))
    apos = _spec_field(location, :alpha_position, length(code.alpha[ai].visits))
    bpos = _spec_field(location, :beta_position, length(code.beta[bi].visits))
    sign_pair = _spec_field(location, :signs, (1, -1))
    length(sign_pair) == 2 || throw(ArgumentError("R move signs must have length two"))
    sign_pair[1] == -sign_pair[2] || throw(ArgumentError("R move requires opposite signs"))
    p = isempty(code.intersections) ? 1 : maximum(point.id for point in code.intersections) + 1
    q = p + 1
    alpha = copy(code.alpha); beta = copy(code.beta)
    alpha[ai] = Curve(:alpha, alpha[ai].id, _insert_after(alpha[ai].visits, apos, p);
                      label=alpha[ai].label, tags=alpha[ai].tags)
    alpha[ai] = Curve(:alpha, alpha[ai].id, _insert_after(alpha[ai].visits, apos + 1, q);
                      label=alpha[ai].label, tags=alpha[ai].tags)
    beta[bi] = Curve(:beta, beta[bi].id, _insert_after(beta[bi].visits, bpos, p);
                     label=beta[bi].label, tags=beta[bi].tags)
    beta[bi] = Curve(:beta, beta[bi].id, _insert_after(beta[bi].visits, bpos + 1, q);
                     label=beta[bi].label, tags=beta[bi].tags)
    points = [copy(code.intersections);
              Intersection(p, alpha_id, beta_id, sign_pair[1]);
              Intersection(q, alpha_id, beta_id, sign_pair[2])]
    new_code = _replace_curves(code, alpha, beta, points)
    old_cycles = length(boundary_cycles(H))
    new_cycles = length(boundary_cycles(new_code))
    blocks = [copy(block) for block in H.code.regions.boundary_cycle_blocks]
    for cycle in old_cycles+1:new_cycles
        push!(blocks, [cycle])
    end
    decorated = decorate(new_code, blocks; metadata=H.code.metadata)
    _finish_move(H, decorated, :R,
        Dict("alpha_curve" => alpha_id, "beta_curve" => beta_id, "intersection_ids" => [p, q]),
        Dict("move" => "remove_canceling_pair", "pair" => [p, q]))
end

function remove_canceling_pair(H::HeegaardDiagram, pair)
    ids = Int.(collect(pair))
    length(ids) == 2 && ids[1] != ids[2] || throw(ArgumentError("R inverse requires two distinct intersection ids"))
    points = _intersection_map(H.code.code)
    all(haskey(points, id) for id in ids) || throw(ArgumentError("R inverse references an unknown intersection"))
    first_point, second_point = points[ids[1]], points[ids[2]]
    (first_point.alpha == second_point.alpha && first_point.beta == second_point.beta &&
     first_point.sign == -second_point.sign) || throw(ArgumentError("pair is not a signed canceling pair"))
    code = H.code.code
    alpha = [Curve(curve.family, curve.id, filter(point -> !(point in ids), curve.visits);
                   label=curve.label, tags=curve.tags) for curve in code.alpha]
    beta = [Curve(curve.family, curve.id, filter(point -> !(point in ids), curve.visits);
                  label=curve.label, tags=curve.tags) for curve in code.beta]
    new_points = [point for point in code.intersections if !(point.id in ids)]
    new_code = _replace_curves(code, alpha, beta, new_points)
    decorated = _new_partition(H, new_code; append_cycles=false)
    _finish_move(H, decorated, :R_inverse,
        Dict("intersection_ids" => ids),
        Dict("move" => "add_canceling_pair", "pair" => ids))
end

"""Apply a fully explicit combinatorial band-sum replacement.

The specification deliberately does not infer a band from a picture.  The
caller must provide `path`, `replacement_visits`, and the post-slide region
partition.  This keeps the operation auditable even though the geometry of a
band is not represented by the core code.
"""
function handle_slide(H::HeegaardDiagram, slide_spec)
    path = _spec_field(slide_spec, :path, nothing)
    replacement = _spec_field(slide_spec, :replacement_visits, nothing)
    regions = _spec_field(slide_spec, :regions, nothing)
    path === nothing && throw(ArgumentError("handle_slide requires explicit band/path data"))
    replacement === nothing && throw(ArgumentError("handle_slide requires replacement_visits"))
    regions === nothing && throw(ArgumentError("handle_slide requires explicit post-slide region partition"))
    source = _spec_field(slide_spec, :source_curve, nothing)
    target = _spec_field(slide_spec, :target_curve, nothing)
    source === nothing && throw(ArgumentError("handle_slide requires source_curve and target_curve"))
    target === nothing && throw(ArgumentError("handle_slide requires source_curve and target_curve"))
    family, source_index = _resolve_curve(H.code.code, source)
    target_family, target_index = _resolve_curve(H.code.code, target)
    family == target_family || throw(ArgumentError("handle_slide source and target must be in one curve family"))
    source_index != target_index || throw(ArgumentError("handle_slide source and target must differ"))
    alpha = copy(H.code.code.alpha); beta = copy(H.code.code.beta)
    source_curve = family == :alpha ? alpha[source_index] : beta[source_index]
    changed = Curve(family, source_curve.id, replacement;
                    label=source_curve.label, tags=source_curve.tags)
    family == :alpha ? (alpha[source_index] = changed) : (beta[source_index] = changed)
    new_code = _replace_curves(H.code.code, alpha, beta, H.code.code.intersections)
    decorated = decorate(new_code, regions; metadata=H.code.metadata)
    _finish_move(H, decorated, :H,
        Dict("source_curve" => source, "target_curve" => target, "path" => deepcopy(path)),
        Dict("move" => "handle_slide", "slide_spec" => deepcopy(slide_spec)))
end

function stabilize(H::HeegaardDiagram, specification=NamedTuple())
    code = H.code.code
    region = _spec_field(specification, :region, 1)
    sign_value = _spec_field(specification, :sign, 1)
    sign_value in (-1, 1) || throw(ArgumentError("stabilization sign must be +1 or -1"))
    alpha_id = isempty(code.alpha) ? 1 : maximum(curve.id for curve in code.alpha) + 1
    beta_id = isempty(code.beta) ? 1 : maximum(curve.id for curve in code.beta) + 1
    point_id = isempty(code.intersections) ? 1 : maximum(point.id for point in code.intersections) + 1
    alpha = [copy(code.alpha); Curve(:alpha, alpha_id, [point_id])]
    beta = [copy(code.beta); Curve(:beta, beta_id, [point_id])]
    points = [copy(code.intersections); Intersection(point_id, alpha_id, beta_id, sign_value)]
    new_code = _replace_curves(code, alpha, beta, points)
    decorated = _new_partition(H, new_code; region=region, append_cycles=true)
    _finish_move(H, decorated, :S,
        Dict("alpha_curve" => alpha_id, "beta_curve" => beta_id, "intersection_id" => point_id),
        Dict("move" => "destabilize", "pair" => [alpha_id, beta_id, point_id]))
end

function destabilize(H::HeegaardDiagram, pair)
    ids = Int.(collect(pair))
    length(ids) == 3 || throw(ArgumentError("destabilization requires [alpha_curve, beta_curve, intersection]"))
    alpha_id, beta_id, point_id = ids
    code = H.code.code
    ai = findfirst(curve -> curve.id == alpha_id, code.alpha)
    bi = findfirst(curve -> curve.id == beta_id, code.beta)
    ai === nothing && throw(ArgumentError("unknown alpha curve id $alpha_id"))
    bi === nothing && throw(ArgumentError("unknown beta curve id $beta_id"))
    code.alpha[ai].visits == [point_id] || throw(ArgumentError("alpha curve is not an isolated stabilization pair"))
    code.beta[bi].visits == [point_id] || throw(ArgumentError("beta curve is not an isolated stabilization pair"))
    alpha = [curve for (index, curve) in enumerate(code.alpha) if index != ai]
    beta = [curve for (index, curve) in enumerate(code.beta) if index != bi]
    points = [point for point in code.intersections if point.id != point_id]
    new_code = _replace_curves(code, alpha, beta, points)
    decorated = _new_partition(H, new_code; append_cycles=false)
    _finish_move(H, decorated, :S_inverse,
        Dict("alpha_curve" => alpha_id, "beta_curve" => beta_id, "intersection_id" => point_id),
        Dict("move" => "stabilize", "pair" => ids))
end

function apply_move(H::HeegaardDiagram, move)
    kind = Symbol(_spec_field(move, :kind, _spec_field(move, :move, "")))
    kind == :epsilon && return reverse_curve_orientation(H, _spec_field(move, :curve, move))
    kind == :R && return add_canceling_pair(H, move)
    kind == :R_inverse && return remove_canceling_pair(H, _spec_field(move, :pair, Int[]))
    kind == :H && return handle_slide(H, move)
    kind == :S && return stabilize(H, move)
    kind == :S_inverse && return destabilize(H, _spec_field(move, :pair, Int[]))
    throw(ArgumentError("unknown Heegaard move $kind"))
end

inverse_move(certificate::MoveCertificate) = deepcopy(certificate.inverse_data)

function enumerate_local_moves(H::HeegaardDiagram; kind=:all)
    result = Any[]
    kind in (:all, :epsilon) && append!(result, [(kind=:epsilon, curve=(curve.family, curve.id))
                                                  for curve in vcat(H.code.code.alpha, H.code.code.beta)])
    if kind in (:all, :destabilize)
        for alpha in H.code.code.alpha, beta in H.code.code.beta
            if length(alpha.visits) == 1 && alpha.visits == beta.visits
                push!(result, (kind=:S_inverse, pair=[alpha.id, beta.id, only(alpha.visits)]))
            end
        end
    end
    result
end
