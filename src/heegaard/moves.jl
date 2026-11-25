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

function _resolve_curve(code::CombinatorialHeegaardCode, selector)
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

function _replace_curves(code::CombinatorialHeegaardCode, alpha, beta, points)
    CombinatorialHeegaardCode(copy(alpha), copy(beta), copy(points), deepcopy(code.metadata))
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

"""
    reverse_curve_orientation(H::CellularHeegaardDiagram, curve) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Return the legacy curve-reversal result and its move certificate.

# Arguments
`curve` selects a family and curve ID, for example `(:alpha,1)`; a bare ID must be unambiguous across families. Ambiguous, invalid-family or missing selectors raise `ArgumentError`. The resulting decorated record must pass `certify_heegaard`; validation and hashing errors propagate.

# Returns
The recertified diagram and a `MoveCertificate` recording the selection, inverse command and source/target identities.

# Guarantees
1. The selected curve's visit order is reversed and its label and tags are retained. [proved/agent: by construction]
2. The input is not modified; some unchanged curve and partition objects remain shared with the result. [proved/agent: by construction]

# Not guaranteed
The legacy sign update compares the selected numeric ID against both families. If an unrelated curve in the other family shares that ID, extra crossing signs can change. Thus general geometric orientation reversal and manifold preservation are not guaranteed on that input; see [the legacy reversal limitation](../specs/heegaard-moves.md#legacy-reversal-id-collision).

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md — curve moves and validation of their witnesses.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> h = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); cellular=true)); r,c = H.reverse_curve_orientation(h, (:alpha,1)); r.code.code.intersections[1].sign
-1
```

Role: move
"""
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

"""
    add_canceling_pair(H::CellularHeegaardDiagram, location) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Insert two opposite-sign intersections at specified legacy alpha and beta positions.

# Arguments
`H` is a certified legacy diagram. `location` supplies alpha and beta curve
IDs, optional insertion positions and opposite signs. Missing/unknown curves,
invalid positions or nonopposite signs raise `ArgumentError` with the failed
requirement; result certification can raise `ValidationError`. The caller
supplies the geometric support required by [the move domain](../specs/heegaard-moves.md#comparison-and-move-witnesses);
the routine does not check a supporting disk between the chosen segments.

# Returns
The certified resulting diagram and a `MoveCertificate` with new
intersection IDs, hashes, inverse annotation and input invariants.

# Guarantees
1. The output contains the two new intersections in the selected visits,
   with opposite signs, and passes current cellular certification.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]
3. The registered genus-one insertion/removal round trip recovers the source
   canonical hash. [verified-computationally: test/heegaard/core.jl, Moves and certificates]

# Not guaranteed
Certification alone does not prove that arbitrary selected positions define
a local manifold-preserving move. The witness boundary in the linked section
applies.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §comparison-and-move-witnesses — legacy combinatorial move records and geometric support.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]])); a, _ = H.add_canceling_pair(d, (alpha_curve=1,beta_curve=1,alpha_position=1,beta_position=1)); length(a.code.code.intersections)
3
```

See also [`remove_canceling_pair`](@ref).

Role: move
"""
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

"""
    remove_canceling_pair(H::CellularHeegaardDiagram, pair) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Remove a specified opposite-sign intersection pair from a legacy cellular diagram.

# Arguments
`H` is a certified diagram in the [legacy move domain](../specs/heegaard-moves.md#comparison-and-move-witnesses).
`pair` must collect to two distinct integer IDs of existing intersections on
the same alpha and beta curves with opposite signs. Checked failures raise
`ArgumentError` identifying the malformed or noncanceling pair; conversion
errors propagate. Result certification can raise `ValidationError`.
The caller supplies geometric cancellation data when a manifold-preserving
interpretation is required; adjacency and a local disk are not checked here.

# Returns
The certified resulting diagram and its `MoveCertificate`, including the
removed IDs, hashes, inverse-operation annotation and input invariant snapshot.

# Guarantees
1. On success the selected intersections are absent from the intersection
   table and all visits, and the output passes cellular certification.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]
3. Removal after the tested genus-one pair insertion recovers the original
   canonical diagram. [verified-computationally: test/heegaard/core.jl, Moves and certificates, S3 fixture]

# Not guaranteed
Opposite signs alone do not certify a local geometric cancellation. The
certificate's inverse annotation does not store all insertion positions.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §comparison-and-move-witnesses — distinction between combinatorial records and geometric witnesses.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]])); a, _ = H.add_canceling_pair(d, (alpha_curve=1, beta_curve=1, alpha_position=1, beta_position=1));

julia> b, _ = H.remove_canceling_pair(a, [2,3]); H.canonical_hash(b) == H.canonical_hash(d)
true
```

See also [`add_canceling_pair`](@ref).

Role: move
"""
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

"""
    handle_slide(H::CellularHeegaardDiagram, slide_spec) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Replace one legacy curve by an explicitly supplied candidate handle slide.

# Arguments
`H` is a certified legacy diagram. `slide_spec` supplies `path`,
`replacement_visits`, `regions`, `source_curve` and `target_curve`.
Missing fields, unresolved curves, different families or equal source and
target raise `ArgumentError`; result certification errors propagate.
The caller supplies a genuine band and compatible region partition in the
[move-witness domain](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Returns
The resulting certified diagram and a `MoveCertificate` recording source,
target, path, hashes and inverse-description data.

# Guarantees
1. The source curve has the supplied replacement visits; the other curve
   words and intersection records remain present. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The supplied path is not geometrically checked against the replacement.
Certification alone does not prove that this is a handle slide preserving
the represented 3-manifold, and the stored inverse description alone is not
a validated inverse band.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md#comparison-and-move-witnesses — checked records versus geometric move witnesses.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> c = H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]); h = H.certify_heegaard(H.decorate(c, [[1]]));

julia> try H.handle_slide(h, (;)) catch e; e isa ArgumentError end
true
```

Role: move
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

"""
    stabilize(H::CellularHeegaardDiagram, specification=NamedTuple()) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Add an isolated one-intersection alpha–beta pair to a legacy cellular diagram.

# Arguments
`H` is a certified diagram in the [legacy move domain](../specs/heegaard-moves.md#comparison-and-move-witnesses).
`specification` supplies optional `region` (default 1) and `sign` (default 1).
Checked: sign must be `1` or `-1`; an invalid sign or selected region raises
`ArgumentError`. Result certification and invariant errors propagate.
Geometric placement in the chosen region is the caller's responsibility
when interpreting the result as a stabilization of a 3-manifold diagram.

# Returns
The certified diagram and a `MoveCertificate` with new curve/intersection
IDs, source and target hashes, destabilization data and input invariants.

# Guarantees
1. The new alpha and beta curves meet only each other once with the selected
   sign; the original curve words and intersections remain present.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]
3. For the tested genus-one fixture, the resulting surface has genus two
   and destabilization recovers the source canonical hash.
   [verified-computationally: test/heegaard/core.jl, Moves and certificates]

# Not guaranteed
The certificate alone does not prove manifold preservation for arbitrary
hand-built region data; see the linked move-witness boundary.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §§4.1, 6.4 — stabilization and its witness domain.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]])); s, _ = H.stabilize(d);

julia> H.surface_genus(s)
2
```

See also [`destabilize`](@ref).

Role: move
"""
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

"""
    destabilize(H::CellularHeegaardDiagram, pair) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Remove a specified isolated stabilization pair from a legacy cellular diagram.

# Arguments
`H` is a certified cellular diagram. `pair` must collect and convert to
three integers identifying the alpha curve, beta curve and intersection.
Both curve visits must consist of that single intersection. Violations throw
`ArgumentError` identifying the malformed tuple, absent curve or non-isolated
pair; conversion errors propagate. Failed output certification throws
`ValidationError`. Geometric interpretation follows
[the move-witness boundary](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Returns
The certified resulting diagram and a `MoveCertificate` with the selected
IDs, input/output hashes, inverse data and input invariant snapshot.

# Guarantees
1. The returned code omits the selected two curves and intersection and
   passes current cellular validation. [proved/agent: by construction]
2. Inputs are not modified. [proved/agent: by construction]
3. The tested stabilization followed by this operation recovers the original
   legacy canonical hash. [verified-computationally: test/heegaard/core.jl, Moves and certificates, genus-one stabilization]

# Not guaranteed
The inverse record does not retain every choice of the original stabilization.
The certificate alone does not prove manifold preservation for arbitrary
hand-built input; its evidence scope is the linked witness boundary.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §§4.1, 6.4 — stabilization and the boundary of executable move records.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> s, _ = H.stabilize(d); back, _ = H.destabilize(s, [2,2,2]);

julia> H.canonical_hash(back) == H.canonical_hash(d)
true
```

See also [`stabilize`](@ref), [`apply_move`](@ref).

Role: move
"""
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

"""
    apply_move(H::CellularHeegaardDiagram, move) -> Tuple{CellularHeegaardDiagram,MoveCertificate}

Apply a selected legacy cellular Heegaard operation and return its certified output and record.

# Arguments
`H` is a certified legacy diagram. `move` is a named tuple, dictionary with
symbol or string keys, or property-bearing object. Its `kind` (or fallback
`move`) selects `:epsilon`, `:R`, `:R_inverse`, `:H`, `:S` or `:S_inverse`.
These call `reverse_curve_orientation`, `add_canceling_pair`,
`remove_canceling_pair`, `handle_slide`, `stabilize` and `destabilize`,
respectively, with their argument requirements. `:epsilon` uses `curve`;
the inverse moves use `pair`; other operations take the whole specification.
An unknown kind throws `ArgumentError` naming it. Selected-operation argument
errors propagate, and failed output certification throws `ValidationError`.
Geometric support is the caller's responsibility within
[the witness boundary](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Returns
The output `CellularHeegaardDiagram` and a `MoveCertificate` recording the
operation, attachment, correspondence, source and target hashes, inverse
data and input invariant snapshot.

# Guarantees
1. A returned output passes the current full cellular validation checks.
   [proved/agent: by construction]
2. The input diagram and move specification are not modified.
   [proved/agent: by construction]
3. The stabilization example raises surface genus from one to two and retains
   trivial first-homology torsion. [verified-computationally: test/heegaard/core.jl, Moves and certificates, genus-one fixture]

# Not guaranteed
Cellular certification alone does not certify that a caller's band is
geometrically valid or that every supported operation preserves a manifold.
The stored inverse data are descriptive and are not necessarily accepted
unchanged by this dispatcher. This is not the resolved-HC catalogue API.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §4.1 — executable records and the boundary of move certification.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> out, cert = apply_move(d, (kind=:S,));

julia> (surface_genus(out), isempty(Topology.HeegaardCodes.first_homology(out).torsion))
(2, true)
```

See also [`stabilize`](@ref), [`handle_slide`](@ref), [`Topology.HeegaardCodes.move`](@ref).

Role: move
"""
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

"""
    inverse_move(certificate::MoveCertificate) -> Dict{String,Any}

Return an independent copy of a legacy move certificate's inverse metadata.

# Arguments
`certificate` is a `MoveCertificate`. The stored inverse metadata is accepted
without checking whether it describes an executable inverse.

# Returns
A deep copy of `certificate.inverse_data`.

# Guarantees
1. The result preserves the stored inverse metadata and its nested values
   under Julia's deep-copy semantics. [proved/agent: by construction]
2. The certificate is not modified. [proved/agent: by construction]

# Not guaranteed
This does not apply or verify an inverse move. In particular, the metadata's
`"move"` names need not be tags accepted by `apply_move`. The distinction
between recorded assertions and checked witnesses is specified in
[the move specification](../specs/heegaard-moves.md#comparison-and-move-witnesses).

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> H = Topology.HeegaardCodes; d = H.certify_heegaard(H.decorate(H.heegaard_code([[1]], [[1]], [H.Intersection(1,1,1,1)]), [[1]]));

julia> _, c = H.reverse_curve_orientation(d, (:alpha, 1)); H.inverse_move(c)["move"]
"reverse_curve_orientation"
```

See also [`apply_move`](@ref), [`MoveCertificate`](@ref).

Role: access
"""
inverse_move(certificate::MoveCertificate) = deepcopy(certificate.inverse_data)

"""
    enumerate_local_moves(H::CellularHeegaardDiagram; kind=:all) -> Vector{Any}

Return legacy orientation-reversal specifications and isolated destabilization candidates.

# Arguments
`H` is a certified cellular diagram. `kind=:epsilon` selects orientation
reversals, `:destabilize` selects isolated-pair candidates, and `:all` selects
both. Other values return an empty vector without error. No move is applied
or output diagram certified during enumeration.

# Returns
A vector of specifications accepted in form by `apply_move`, carrying a
kind and curve selector or pair of curve IDs with an intersection ID.

# Guarantees
1. The selected reversal specifications cover the stored curves; selected
   destabilization specifications refer to alpha/beta curves whose visit
   lists consist of the same single intersection.
   [proved/agent: by construction]
2. These are candidates in the scope of
   [legacy move records](../specs/heegaard-moves.md#comparison-and-move-witnesses).
   The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is not an exhaustive move search: it lists no handle slides, pair
insertions or stabilizations. A candidate can fail when applied and certified.
It is not the resolved-HC catalogue's `movable` operation.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-moves.md §4.1 — candidate records and executable validation are distinct.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> d = certify_heegaard(decorate(heegaard_code([[1]], [[1]], [Intersection(1,1,1,1)]), [[1]]));

julia> !isempty(enumerate_local_moves(d; kind=:epsilon))
true
```

See also [`apply_move`](@ref), [`destabilize`](@ref), [`Topology.HeegaardCodes.movable`](@ref).

Role: access
"""
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
