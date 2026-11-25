"""Raw integer-weighted move evidence, separate from the unweighted OGraph schema.

This envelope certifies deterministic executor replay and records equation checks.
It is not an independent mathematical proof of framing transport. Its hashes bind
the full raw presentation and integer weights; they are NOT canonical state IDs.
"""
const WEIGHTED_MOVE_EVIDENCE_SCHEMA = "akaghef.weighted-move-evidence/2"
const WEIGHTED_MOVE_VERIFIER_VERSION = "weighted-executor-replay/3"
const WEIGHTED_MOVE_SOURCE_RULE_VERSION = "integer-mp-h-02-ps-runtime-cut/2"
const WEIGHTED_RELATIVE_RESULT_SCHEMA = "akaghef.weighted-relative-result/2"
const WEIGHTED_BOUNDARY_TRANSPORT_VERSION = 1

"""
    WeightedMoveEvidence(data) -> WeightedMoveEvidence

Wrap the data of an integer-weighted move evidence packet.

# Arguments
`data` must be or convert to `Dict{String,Any}`; conversion errors propagate.
This constructor does not validate the schema or the packet's contents.
The caller supplies evidence conforming to
[weighted evidence §12.7](../specs/virtual-link-moves.md#127-weighted-evidence)
when validated replay is intended.

# Returns
`WeightedMoveEvidence` with a `data::Dict{String,Any}` field containing the
supplied packet data.

# Guarantees
1. The dictionary is retained after conversion, without replaying the packet.
   [proved/agent: by construction]
2. Input data are not modified; a correctly typed dictionary is shared.
   [proved/agent: by construction]

# Not guaranteed
The wrapper alone is not validated evidence. The executor-replay boundary
of §12.7 applies: no independent framing-preservation proof or canonical
manifold identity is conferred by the packet or its raw hashes.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.7 — weighted evidence and the executor-replay boundary.

# Examples
```jldoctest
julia> source = lens_space_link(2,1); set_weight!(source, [1,-1,-1,0]); candidate = only(movable(source, "H"; param=1, vertices=[1])); r = verify_weighted_move(source, candidate; framed=true);

julia> packet = WeightedMoveEvidence(weighted_move_evidence_data(r.evidence));

julia> replay_weighted_move(packet).valid
true
```

Role: construct
"""
struct WeightedMoveEvidence
    data::Dict{String,Any}
end

"""
    WeightedMoveVerification(valid, reason, evidence, target, error) -> WeightedMoveVerification

Store the outcome and optional artifacts of a weighted move verification.

# Arguments
Fields must have or convert to `Bool`, `Symbol`, optional
`WeightedMoveEvidence`, optional `VirtualLink`, and optional `String`,
respectively. Conversion errors propagate. No consistency check between
the flag, reason, packet and target is performed by this constructor.
Verification semantics are owned by
[weighted evidence §12.7](../specs/virtual-link-moves.md#127-weighted-evidence).

# Returns
`WeightedMoveVerification` with a validity flag, reason, optional evidence,
optional target diagram and optional explanatory error text.

# Guarantees
1. Fields retain the supplied values after conversion. Evidence and target
   objects are shared, not independently verified or copied.
   [proved/agent: by construction]
2. The constructor does not modify any input. [proved/agent: by construction]

# Not guaranteed
Manually setting `valid=true` does not establish successful replay or a
mathematical theorem. Even a verified executor result has the framing and
canonical-identity limits of §12.7.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.7 — results of weighted executor verification.

# Examples
```jldoctest
julia> source = lens_space_link(2,1); set_weight!(source, [1,-1,-1,0]); candidate = only(movable(source, "H"; param=1, vertices=[1])); r = verify_weighted_move(source, candidate; framed=true);

julia> stored = WeightedMoveVerification(r.valid, r.reason, r.evidence, r.target, r.error);

julia> stored.valid && stored.target.weights == [[1.0,-2.0,-1.0,1.0]]
true
```

Role: convert
"""
struct WeightedMoveVerification
    valid::Bool
    reason::Symbol
    evidence::Union{Nothing,WeightedMoveEvidence}
    target::Union{Nothing,VirtualLink}
    error::Union{Nothing,String}
end

_weighted_invalid(reason; error=nothing) = WeightedMoveVerification(
    false, reason, nothing, nothing, error===nothing ? nothing : sprint(showerror,error))
_weighted_json_copy(value) = JSON.parse(ograph_canonical_json(value); dicttype=Dict{String,Any})
_weighted_digest(value) = bytes2hex(SHA.sha256(codeunits(ograph_canonical_json(value))))
_weighted_decode_int(value) = value isa AbstractString ? parse(BigInt,value) : BigInt(value)
_weighted_decode_targets(value) = value===nothing ? nothing : _weighted_decode_int.(value)

function _weighted_exact_int(value)
    value isa Real && !(value isa Bool) && isfinite(value) && isinteger(value) ||
        throw(ArgumentError("weighted evidence requires finite integer weights"))
    abs(big(value)) <= big(2)^53-1 || throw(ArgumentError(
        "weighted evidence exceeds the consecutive exact Float64 integer range"))
    Int(value)
end

function _weighted_raw_data(v::VirtualLink)
    weight_mode(v)==:integer || throw(ArgumentError(
        "integer weighted evidence does not support Z2 mode; no implicit coefficient conversion is allowed"))
    v.is_weighted || throw(ArgumentError("weighted evidence requires assigned weights"))
    _validate_gauss(v.gauss,v.orientation;real_only=v.real_only)
    all(sign -> sign in (-1,0,1),v.orientation) ||
        throw(ArgumentError("invalid crossing orientations"))
    length(v.weights)==length(v.gauss)==length(v.strand_types) ||
        throw(DimensionMismatch("raw component metadata does not match Gauss components"))
    code,_=real_gauss_code(v)
    for (word,weights) in zip(code,v.weights)
        length(weights)==length(word) ||
            throw(DimensionMismatch("raw weights must follow real crossing occurrences"))
    end
    (gauss=deepcopy(v.gauss), orientation=copy(v.orientation),
     weights=[_weighted_exact_int.(w) for w in v.weights], real_only=v.real_only,
     strand_types=string.(v.strand_types), is_weighted=true,
     cuts=[_is_cut_component(word) for word in v.gauss],
     boundary_weights=[_is_cut_component(word) ?
         [_weighted_exact_int(first(v.weights[index])),_weighted_exact_int(last(v.weights[index]))] : nothing
         for (index,word) in enumerate(v.gauss)])
end

"""
    weighted_raw_hash(v::VirtualLink) -> String

Hash the complete raw integer-weighted presentation in the weighted evidence schema.

# Arguments
Checked: E1 of [weighted evidence](../specs/virtual-link-moves.md#127-weighted-evidence),
including assigned integer mode, Gauss/orientation validity, component metadata
and safely representable integer weights. Invalid values raise `ArgumentError`,
misalignment `DimensionMismatch`; conversion or encoding errors propagate.
Caller: retain the exact raw presentation when using the digest for replay.

# Returns
A lowercase hexadecimal SHA-256 digest string of the versioned raw packet.

# Guarantees
1. The same validated raw packet under this schema has the same digest,
   including its weight and cut metadata. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is not a cyclic/relabeling canonical ID, an injective encoding, or a
manifold invariant. Hash equality alone is not authenticated evidence.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — raw weighted identity scope.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); length(weighted_raw_hash(v))
64
```

Role: invariant
"""
weighted_raw_hash(v::VirtualLink) = _weighted_digest((
    schema="akaghef.raw-integer-weighted-virtual-link/2",state=_weighted_raw_data(v)))

function _weighted_link_from_data(value)
    _move_data_object(value,("gauss","orientation","weights","real_only","strand_types","is_weighted",
                             "cuts","boundary_weights"))
    field(key)=_move_data_field(value,key)
    _move_data_bool(field("is_weighted"),"is_weighted") ||
        throw(ArgumentError("unweighted states do not belong to this evidence schema"))
    code=_move_data_vector(field("gauss"),nothing,_move_data_int_vector,"gauss")
    signs=_move_data_int_vector(field("orientation"),"orientation")
    weights=_move_data_vector(field("weights"),length(code),_move_data_int_vector,"weights")
    foreach(w -> foreach(_weighted_exact_int,w),weights)
    real_only=_move_data_bool(field("real_only"),"real_only")
    strands=_move_data_vector(field("strand_types"),length(code),_move_data_string,"strand_types")
    all(s -> length(s)==1,strands) || throw(ArgumentError("strand types must be single characters"))
    cuts=_move_data_vector(field("cuts"),length(code),_move_data_bool,"cuts")
    boundary=_move_data_vector(field("boundary_weights"),length(code),
        (entry,label)->begin
            entry===nothing && return nothing
            _move_data_vector(entry,2,_move_data_int,label)
        end,"boundary_weights")
    all(cuts[i] == _is_cut_component(code[i]) for i in eachindex(code)) ||
        throw(ArgumentError("raw cut flags do not match trailing markers"))
    link=VirtualLink(code,signs,[Float64.(w) for w in weights],real_only,true,first.(strands))
    _weighted_raw_data(link)
    link
end

function _weighted_pattern_data(pattern)
    pattern===nothing && return nothing
    (parameter=pattern.parameter,gauss=deepcopy(pattern.gauss),
     orientation=copy(pattern.orientation),capacity=copy(pattern.capacity),
     gauss_first=copy(pattern.gauss_first),gauss_last=copy(pattern.gauss_last),
     vertex_count=pattern.vertex_count,
     weights=[[isnan(x) ? "NaN" : _weighted_exact_int(x) for x in w] for w in pattern.weights])
end

function _weighted_pattern_from_data(value)
    value===nothing && return nothing
    _move_data_object(value,("parameter","gauss","orientation","capacity","gauss_first",
                             "gauss_last","vertex_count","weights"))
    field(key)=_move_data_field(value,key)
    weight(x,label)=x=="NaN" ? NaN : Float64(_weighted_exact_int(_move_data_int(x,label)))
    weights=_move_data_vector(field("weights"),nothing,
        (w,label)->_move_data_vector(w,nothing,weight,label),"weights")
    MovePattern(_move_data_string(field("parameter"),"parameter"),
        _move_data_vector(field("gauss"),nothing,_move_data_int_vector,"gauss"),
        _move_data_int_vector(field("orientation"),"orientation"),
        _move_data_int_vector(field("capacity"),"capacity"),
        _move_data_int_vector(field("gauss_first"),"gauss_first"),
        _move_data_int_vector(field("gauss_last"),"gauss_last"),
        _move_data_int(field("vertex_count"),"vertex_count"),weights)
end

_weighted_candidate_data(c::MoveCandidate) = (
    type=c.type,parameter=c.parameter,vertices=copy(c.vertices),edges=collect.(c.edges),
    direction=String(c.direction),before=_weighted_pattern_data(c.before),
    after=_weighted_pattern_data(c.after),locations=collect.(c.locations),
    ps_split_k=c.ps_split_k===nothing ? nothing : _weighted_exact_int(c.ps_split_k))

function _weighted_candidate_from_data(value)
    _move_data_object(value,("type","parameter","vertices","edges","direction","before","after","locations","ps_split_k"))
    field(key)=_move_data_field(value,key)
    pairs(key)=Tuple.(_move_data_vector(field(key),nothing,
        (v,label)->_move_data_int_vector(v,label;length_=2),key))
    direction=_move_data_string(field("direction"),"direction")
    direction in ("forward","inverse","weight") || throw(ArgumentError("invalid move direction"))
    MoveCandidate(_move_data_string(field("type"),"type"),
        _move_data_string(field("parameter"),"parameter"),
        _move_data_int_vector(field("vertices"),"vertices"),pairs("edges"),Symbol(direction),
        _weighted_pattern_from_data(field("before")),_weighted_pattern_from_data(field("after")),pairs("locations"),
        field("ps_split_k")===nothing ? nothing : BigInt(_weighted_exact_int(_move_data_int(field("ps_split_k"),"ps_split_k"))))
end

function _weighted_catalogue_candidate(source,candidate)
    # Match the frozen site AND all template fields against freshly generated
    # source rules. A valid move at another site never rescues a stale candidate.
    geometry_source=_weighted_unweight(source)
    choices = if candidate.type=="PS" && candidate.direction==:forward
        _movable_geometry(geometry_source,"PS";param=candidate.parameter,
            edge_occurrences=candidate.locations,strict=true)
    else
        _movable_geometry(geometry_source,candidate.type;param=candidate.parameter,
            vertices=candidate.vertices,edges=candidate.edges,strict=true)
    end
    geometric=MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,candidate.edges,
        candidate.direction,candidate.before,candidate.after,candidate.locations)
    identity=ograph_canonical_json(_weighted_candidate_data(geometric))
    index=findfirst(c -> ograph_canonical_json(_weighted_candidate_data(c))==identity,choices)
    index===nothing && return nothing
    canonical=deepcopy(choices[index])
    MoveCandidate(canonical.type,canonical.parameter,canonical.vertices,canonical.edges,
        canonical.direction,canonical.before,canonical.after,canonical.locations,candidate.ps_split_k)
end

# Core transport traces use Julia symbols and Pair vertex mappings. Their JSON
# representation is deliberately mechanical; replay recalculates every field.
_weighted_detail_data(value::Symbol)=String(value)
_weighted_detail_data(value::Pair)=[_weighted_detail_data(first(value)),_weighted_detail_data(last(value))]
_weighted_detail_data(value::NamedTuple)=Dict(String(k)=>_weighted_detail_data(getfield(value,k)) for k in keys(value))
_weighted_detail_data(value::AbstractDict)=Dict(String(k)=>_weighted_detail_data(v) for (k,v) in value)
_weighted_detail_data(value::Union{Tuple,AbstractVector})=[_weighted_detail_data(v) for v in value]
_weighted_detail_data(value)=value

function _weighted_equation_data(v;framed,boundary=nothing)
    # Accumulate in BigInt independently of the machine-Int public residual.
    measurement_only=boundary===nothing && any(_is_cut_component, v.gauss) && !framed
    equation=if measurement_only
        data=boundary_weight_data(v)
        (edge_weights=data.edge_weights, residual=data.measured,
         boundary_targets=nothing,
         boundary_start=data.boundary_start,boundary_end=data.boundary_end,
         delta=data.delta,cp2=data.cp2)
    else
        weight_equation(v;boundary=boundary)
    end
    measured=measurement_only ? equation.residual : boundary_weight_data(v).measured
    disks=disk_table(v)
    residual=BigInt.(equation.residual)
    closed=all(!=(0),v.orientation) && all(w->all(!=(0),w),v.gauss) && is_closed(v)[1]
    framed && !closed && boundary===nothing && throw(ArgumentError("open framed evidence requires explicit boundary targets"))
    framed && any(!=(0),residual) && throw(ArgumentError("nonzero framing residual"))
    (framed_requested=framed,closed_normal=closed,residual=string.(residual),
     real_edge_weights=Int.(equation.edge_weights),
     measured_boundary_residual=string.(measured),
     boundary_targets=equation.boundary_targets===nothing ? nothing : string.(equation.boundary_targets),
     boundary_targets_supplied=boundary!==nothing,
     boundary_start=_weighted_detail_data(equation.boundary_start),
     boundary_end=_weighted_detail_data(equation.boundary_end))
end

function _weighted_unweight(v)
    result=deepcopy(v)
    result.is_weighted=false
    result.weights=[zeros(length(w)) for w in v.weights]
    result.weight_mode=:integer
    result
end

"""Return the raw presentation fields that a relative weighted result must share."""
function _weighted_topology_raw(v::VirtualLink)
    _validate_gauss(v.gauss,v.orientation;real_only=v.real_only)
    (gauss=deepcopy(v.gauss),orientation=copy(v.orientation),real_only=v.real_only,
     strand_types=string.(v.strand_types),
     cuts=[_is_cut_component(word) for word in v.gauss])
end

"""Replay a relative candidate on the source after discarding assigned weights."""
function _weighted_relative_topology_target(source::VirtualLink,candidate::MoveCandidate)
    replay_source=_weighted_unweight(source)
    _apply_candidate!(replay_source,candidate)
    replay_source
end

function _weighted_topological_data(source,target,candidate)
    applicable=candidate.type in ("MP","PS") && length(source.gauss)==1 &&
        all(!=(0),source.orientation) && all(!=(0),only(source.gauss)) && is_closed(source)[1]
    applicable || return (status="outside-existing-unweighted-verifier-scope",witness=nothing)
    topcandidate=MoveCandidate(candidate.type,candidate.parameter,copy(candidate.vertices),
        NTuple{2,Int}[],candidate.direction,candidate.before,candidate.after,copy(candidate.locations))
    result=verify_ograph_move(_weighted_unweight(source),topcandidate;
        expected_raw_target=_weighted_unweight(target))
    result.valid || throw(ArgumentError("unweighted topological verification failed: $(result.reason)"))
    (status="verified-by-existing-unweighted-verifier",witness=ograph_move_witness_data(result.witness))
end

"""Recompute and encode the actual open-end lineage of a frozen move.

The transport executor works on the unweighted presentation because CP/BMP/B02
have no assigned integer transport rule.  The returned target is nevertheless
checked against the weighted executor target by the callers.  Thus this record
binds endpoint identity to the executed geometry, while leaving measured and
prescribed boundary weights in their separate equation fields.
"""
function _weighted_boundary_transport_data(source::VirtualLink,
                                           candidate::MoveCandidate)
    geometric_source=candidate.type=="H" ? source : _weighted_unweight(source)
    geometric_candidate=candidate.ps_split_k===nothing ? candidate :
        MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,candidate.edges,
            candidate.direction,candidate.before,candidate.after,candidate.locations)
    transport=move_boundary_transport(geometric_source,geometric_candidate)
    (version=WEIGHTED_BOUNDARY_TRANSPORT_VERSION,
     kind="endpoint-lineage",
     source_cut_count=count(_is_cut_component,geometric_source.gauss),
     target_cut_count=count(_is_cut_component,transport.target.gauss),
     endpoints=transport.endpoints,
     target=_weighted_topology_raw(transport.target))
end

function _weighted_checked_boundary_transport_data(source::VirtualLink,
                                                   target::VirtualLink,
                                                   candidate::MoveCandidate)
    detail=_weighted_boundary_transport_data(source,candidate)
    expected=_weighted_topology_raw(target)
    actual=_weighted_detail_data(detail.target)
    ograph_canonical_json(actual)==ograph_canonical_json(_weighted_detail_data(expected)) ||
        throw(ArgumentError("boundary endpoint trace and executed move disagree"))
    detail
end

function _weighted_boundary_pair(boundary)
    boundary===nothing && return (nothing,nothing)
    if boundary isa NamedTuple && haskey(boundary,:source) && haskey(boundary,:target)
        return (boundary.source,boundary.target)
    elseif boundary isa AbstractDict && haskey(boundary,"source") && haskey(boundary,"target")
        return (boundary["source"],boundary["target"])
    elseif boundary isa Tuple && length(boundary)==2 &&
           !(boundary[1] isa Integer && boundary[2] isa Integer)
        return (boundary[1],boundary[2])
    end
    (boundary,boundary)
end

"""
    verify_weighted_move(source::VirtualLink, candidate::MoveCandidate; framed=false, expected_raw_target=nothing, boundary=nothing) -> WeightedMoveVerification

Return a verification report with the executed integer-weighted move's target, raw-state evidence and equation checks.

# Arguments
Checked: E1–E4 of [weighted evidence](../specs/virtual-link-moves.md#127-weighted-evidence),
including integer mode, supported candidate, raw storage, Boolean `framed`,
and source/target equations when requested. Unsupported or mismatched data
return an invalid report; ordinary execution or encoding errors become its
reason/error fields. Interruption propagates.
Caller: choose `boundary` as the source/target prescription of §12.7, and
supply any geometric interpretation beyond executor evidence.
`expected_raw_target` is optional raw data or a `VirtualLink` to compare.

# Returns
A verification record; success has `valid=true`, `reason=:verified`,
the executed target and a versioned evidence packet. Failure has no target
or evidence and may contain error text.

# Guarantees
1. Success binds the current catalogue candidate, raw source and target,
   measured equations, boundary lineage and any requested target comparison.
   [proved/agent: by construction]
2. With `framed=true`, both source and target have zero supported equation
   residuals. Existing unweighted verification is included only within its
   documented scope. [proved/agent: by construction]
3. The source, candidate and expected target are not modified.
   [proved/agent: by construction]

# Not guaranteed
The report does not independently prove preservation of a framing. In
particular, a supplied boundary target does not turn a nonclosed presentation
into a closed-normal one; the recorded closed-normal flag must be read.
Raw hashes are not manifold identities. Measurement-only cut evidence
does not assert prescribed relative targets.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — executor verification and its geometric boundary.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); c = only(movable(v, "H"; param=1, vertices=[1]));

julia> r = verify_weighted_move(v, c; framed=true); (r.valid, r.reason)
(true, :verified)
```

Role: condition
"""
function verify_weighted_move(source::VirtualLink,candidate::MoveCandidate;
                              framed=false,expected_raw_target=nothing,boundary=nothing)
    weight_mode(source)==:integer || return _weighted_invalid(:unsupported_weight_mode)
    candidate.type in ("MP","H","02","PS") || return _weighted_invalid(:unsupported_move)
    # Current supported weighted sites use vertices or explicit occurrences.
    # An ignored endpoint selector must not become unauthenticated evidence.
    isempty(candidate.edges) || return _weighted_invalid(:candidate_mismatch)
    candidate.ps_split_k===nothing ||
        (candidate.type=="PS" && candidate.direction==:forward) ||
        return _weighted_invalid(:candidate_mismatch)
    stage=:invalid_source
    try
        framed isa Bool || throw(ArgumentError("framed must be Boolean"))
        source_boundary,target_boundary=_weighted_boundary_pair(boundary)
        raw_source=_weighted_raw_data(source)
        source_equation=_weighted_equation_data(source;framed=framed,boundary=source_boundary)
        stage=:candidate_mismatch
        frozen=_weighted_catalogue_candidate(source,candidate)
        frozen===nothing && return _weighted_invalid(:candidate_mismatch)
        stage=:application_failed
        details = if frozen.type=="PS"
            plan=_weighted_ps_plan(source,frozen;
                boundary=(source=source_boundary,target=target_boundary))
            (coordinate_system="real-gauss-components/head-occurrences",
             ps=plan.detail,transport=_weighted_transport_detail(source,plan.candidate))
        else
            (coordinate_system="real-gauss-components/head-occurrences",
             ps=nothing,transport=_weighted_transport_detail(source,frozen))
        end
        target=move(source,frozen;
            boundary=(source=source_boundary,target=target_boundary))
        stage=:invalid_target
        raw_target=_weighted_raw_data(target)
        target_equation=_weighted_equation_data(target;framed=framed,boundary=target_boundary)
        if expected_raw_target!==nothing
            expected=expected_raw_target isa VirtualLink ? _weighted_raw_data(expected_raw_target) : expected_raw_target
            ograph_canonical_json(expected)==ograph_canonical_json(raw_target) ||
                return _weighted_invalid(:raw_target_mismatch)
        end
        stage=:topological_verification_failed
        topological=_weighted_topological_data(source,target,frozen)
        boundary_transport=_weighted_checked_boundary_transport_data(source,target,frozen)
        stage=:evidence_encoding_failed
        data=_weighted_json_copy((schema=WEIGHTED_MOVE_EVIDENCE_SCHEMA,
            verifier_version=WEIGHTED_MOVE_VERIFIER_VERSION,
            source_rule_version=WEIGHTED_MOVE_SOURCE_RULE_VERSION,
            evidence_kind="executor-replay-with-equation-checks",identity="raw-weighted-only",
            source_raw=raw_source,target_raw=raw_target,
            source_raw_sha256=weighted_raw_hash(source),target_raw_sha256=weighted_raw_hash(target),
            candidate=_weighted_candidate_data(frozen),source_equation=source_equation,
            target_equation=target_equation,
            boundary=(source=source_equation.boundary_targets_supplied ? source_equation.boundary_targets : nothing,
                      target=target_equation.boundary_targets_supplied ? target_equation.boundary_targets : nothing),
            topological=topological,
            transport_detail=_weighted_detail_data(merge(details,
                (boundary_transport=boundary_transport,)))))
        WeightedMoveVerification(true,:verified,WeightedMoveEvidence(data),target,nothing)
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(stage;error=error)
    end
end

"""
    verify_weighted_relative_result(source::VirtualLink, target::VirtualLink, candidate::MoveCandidate; boundary=nothing) -> WeightedMoveVerification

Return a verification report for an externally assigned weighted CP, BMP or B02 target and its boundary data.

# Arguments
Checked: E1, E3 and E5 of [weighted evidence](../specs/virtual-link-moves.md#127-weighted-evidence),
including assigned integer-mode source and target, supported move family,
raw data and equality to the unweighted geometric target.
Ordinary validation/execution failures return invalid reports; interruption
propagates. Caller: supply a genuine current geometric candidate, independently
chosen target weights and, if desired, source/target `boundary` prescriptions;
this interface does not choose them or freshly authenticate the catalogue templates.

# Returns
A verification record with `reason=:assigned_result`, evidence and a copied
target on success; otherwise an invalid report without target or evidence.

# Guarantees
1. Success binds the actual geometric target, assigned raw weights, directed
   endpoint lineage and equation measurements. [proved/agent: by construction]
2. Source, target, candidate and boundary data are not modified; the successful
   returned target is copied. [proved/agent: by construction]

# Not guaranteed
No generic integer-weight transport rule for these move families is asserted.
Raw replay of supplied candidate data does not authenticate a hand-built template.
Prescribed targets are recorded but residuals need not vanish; success is not
a framing-preservation or geometric boundary-extension theorem.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — assigned-result evidence.

# Examples
```jldoctest
julia> spin = VirtualLink(; weight_mode=:z2); integer = VirtualLink(); c = MoveCandidate("CP", "dir", Int[], NTuple{2,Int}[], :forward, nothing, nothing, Tuple{Int,Int}[]);

julia> verify_weighted_relative_result(spin, integer, c).reason
:unsupported_weight_mode
```

Role: condition
"""
function verify_weighted_relative_result(source::VirtualLink,target::VirtualLink,
                                         candidate::MoveCandidate;boundary=nothing)
    weight_mode(source)==weight_mode(target)==:integer ||
        return _weighted_invalid(:unsupported_weight_mode)
    candidate.type in ("CP","BMP","B02") || return _weighted_invalid(:unsupported_move)
    source.is_weighted && target.is_weighted || return _weighted_invalid(:invalid_source)
    try
        source_boundary,target_boundary=_weighted_boundary_pair(boundary)
        source_raw=_weighted_raw_data(source); target_raw=_weighted_raw_data(target)
        expected_topology=_weighted_relative_topology_target(source,candidate)
        _weighted_digest(_weighted_topology_raw(expected_topology)) ==
            _weighted_digest(_weighted_topology_raw(target)) ||
            return _weighted_invalid(:raw_target_mismatch)
        boundary_transport=_weighted_checked_boundary_transport_data(source,target,candidate)
        source_equation=_weighted_equation_data(source;framed=false,boundary=source_boundary)
        target_equation=_weighted_equation_data(target;framed=false,boundary=target_boundary)
        data=_weighted_json_copy((schema=WEIGHTED_RELATIVE_RESULT_SCHEMA,
            evidence_kind="assigned-result-with-boundary-evidence",
            source_raw=source_raw,target_raw=target_raw,
            source_raw_sha256=weighted_raw_hash(source),target_raw_sha256=weighted_raw_hash(target),
            candidate=_weighted_candidate_data(candidate),
            boundary=(source=source_equation.boundary_targets_supplied ? source_equation.boundary_targets : nothing,
                      target=target_equation.boundary_targets_supplied ? target_equation.boundary_targets : nothing),
            source_equation=source_equation,target_equation=target_equation,
            transport_detail=_weighted_detail_data((coordinate_system="real-gauss-components/head-occurrences",
                boundary_transport=boundary_transport))))
        WeightedMoveVerification(true,:assigned_result,WeightedMoveEvidence(data),deepcopy(target),nothing)
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(:assigned_result_invalid;error=error)
    end
end

"""
    replay_weighted_relative_result(value) -> WeightedMoveVerification

Return a verification report comparing an assigned weighted move result with its saved boundary evidence.

# Arguments
`value` is a `WeightedMoveEvidence` or dictionary-like packet in the
assigned-result schema of
[weighted evidence §12.7](../specs/virtual-link-moves.md#127-weighted-evidence).
Schema, raw hashes, candidate, target, boundary data and reconstructed envelope
are checked. Ordinary malformed-input and verification failures return an
invalid result; `InterruptException` propagates.

# Returns
`WeightedMoveVerification` with `valid`, `reason`, optional recomputed
evidence, optional reconstructed target and optional error text. A successful
result has the assigned-result semantics of [`verify_weighted_relative_result`](@ref).

# Guarantees
1. Success means the raw hashes and complete saved envelope agree with
   the current assigned-result verification for the reconstructed source,
   target and candidate. [proved/agent: by construction]
2. The supplied packet is not modified. [proved/agent: by construction]

# Not guaranteed
This validates the assigned target and boundary record, not an integer-weight
transport law for CP/BMP/B02. It is not an independent proof of framing
preservation or a canonical manifold identity. An invalid result can mean
malformed or obsolete evidence rather than mathematical inequivalence.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-moves.md §12.7 — assigned-result replay and its distinct evidence scope.

# Examples
```jldoctest
julia> source = VirtualLink(); set_data!(source; rgauss=[[1,2,0],[-1,0],[-2,0]], orientation=[-1,1]); set_weight!(source, [[10,20,30],[40,50],[60,70]]);

julia> geometric = Topology._weighted_unweight(source); candidate = only(movable(geometric, "BMP"; param=0, vertices=[1,2], strict=true)); target = move_boundary_transport(geometric, candidate).target; set_weight!(target, [zeros(length(word)) for word in target.gauss]);

julia> result = verify_weighted_relative_result(source, target, candidate);

julia> replay_weighted_relative_result(weighted_move_evidence_data(result.evidence)).valid
true
```

Role: condition
"""
function replay_weighted_relative_result(value)
    try
        data=value isa WeightedMoveEvidence ? value.data : value
        data["schema"]==WEIGHTED_RELATIVE_RESULT_SCHEMA || return _weighted_invalid(:malformed_evidence)
        _move_data_object(data,("schema","evidence_kind","source_raw","target_raw",
            "source_raw_sha256","target_raw_sha256","candidate","boundary",
            "source_equation","target_equation","transport_detail"))
        source=_weighted_link_from_data(data["source_raw"])
        target=_weighted_link_from_data(data["target_raw"])
        weighted_raw_hash(source)==data["source_raw_sha256"] || return _weighted_invalid(:source_hash_mismatch)
        weighted_raw_hash(target)==data["target_raw_sha256"] || return _weighted_invalid(:target_hash_mismatch)
        candidate=_weighted_candidate_from_data(data["candidate"])
        decode_boundary(value)=_weighted_decode_targets(value)
        replay=verify_weighted_relative_result(source,target,candidate;
            boundary=(source=decode_boundary(data["boundary"]["source"]),
                      target=decode_boundary(data["boundary"]["target"])))
        replay.valid || return replay
        ograph_canonical_json(replay.evidence.data)==ograph_canonical_json(data) ||
            return _weighted_invalid(:evidence_mismatch)
        replay
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(:malformed_evidence;error=error)
    end
end

"""
    weighted_move_evidence_data(evidence::WeightedMoveEvidence) -> Dict{String,Any}

Return a detached copy of a weighted evidence envelope.

# Arguments
Checked: only the `WeightedMoveEvidence` argument type; no schema or validity
check is performed. Caller: a valid
[weighted packet](../specs/virtual-link-moves.md#127-weighted-evidence)
is required for the returned data to have evidence semantics.

# Returns
A deep copy of the wrapper's data dictionary.

# Guarantees
1. Stored contents are retained without adding a verification claim.
   [proved/agent: by construction]
2. The input is unchanged and ordinary nested packet containers are detached.
   [proved/agent: by construction]

# Not guaranteed
An arbitrary manually wrapped dictionary need not be valid or JSON-compatible.
This function does not replay the evidence.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — packet copying versus verification.

# Examples
```jldoctest
julia> d = Dict{String,Any}("sample" => [1]); e = WeightedMoveEvidence(d); c = weighted_move_evidence_data(e); c["sample"][1] = 2;

julia> d["sample"] == [1]
true
```

Role: convert
"""
weighted_move_evidence_data(evidence::WeightedMoveEvidence)=deepcopy(evidence.data)

"""
    weighted_move_evidence_from_data(value) -> WeightedMoveEvidence

Decode a current integer-weighted executor evidence envelope without replaying it.

# Arguments
Checked: required packet fields, schema and executor/rule versions, raw state
and candidate decoding under E1–E4 of
[weighted evidence](../specs/virtual-link-moves.md#127-weighted-evidence).
Unsupported versions or invalid fields raise `ArgumentError`; malformed
containers and encoding failures may propagate other exceptions.
Caller: provide dictionary-like packet data, not an assigned-result envelope.

# Returns
A `WeightedMoveEvidence` containing detached JSON-compatible packet data.

# Guarantees
1. Successful decoding establishes the current envelope shape and decodability
   of its raw states and candidate. [proved/agent: by construction]
2. The supplied value is not modified. [proved/agent: by construction]

# Not guaranteed
Hashes, equations and saved outputs are not authenticated by decoding;
call `replay_weighted_move` before treating the packet as verified evidence.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — versioned decoding and replay.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); c = only(movable(v, "H"; param=1, vertices=[1])); r = verify_weighted_move(v, c; framed=true);

julia> e = weighted_move_evidence_from_data(weighted_move_evidence_data(r.evidence)); replay_weighted_move(e).valid
true
```

Role: construct
"""
function weighted_move_evidence_from_data(value)
    _move_data_object(value,("schema","verifier_version","source_rule_version","evidence_kind",
        "identity","source_raw","target_raw","source_raw_sha256","target_raw_sha256",
        "candidate","source_equation","target_equation","boundary","topological","transport_detail"))
    field(key)=_move_data_field(value,key)
    field("schema")==WEIGHTED_MOVE_EVIDENCE_SCHEMA || throw(ArgumentError("unsupported weighted evidence schema"))
    field("verifier_version")==WEIGHTED_MOVE_VERIFIER_VERSION || throw(ArgumentError("unsupported weighted verifier version"))
    field("source_rule_version")==WEIGHTED_MOVE_SOURCE_RULE_VERSION || throw(ArgumentError("unsupported weighted source rule version"))
    _weighted_link_from_data(field("source_raw"))
    _weighted_link_from_data(field("target_raw"))
    _weighted_candidate_from_data(field("candidate"))
    WeightedMoveEvidence(_weighted_json_copy(value))
end

"""
    replay_weighted_move(value) -> WeightedMoveVerification

Return a verification report from replaying a weighted move packet's candidate from its saved raw source.

# Arguments
Checked: current schema and versions, decoded raw data, hashes, candidate,
requested equations, target and full envelope agreement. Ordinary failures
return an invalid result; interruption propagates.
Caller: `value` is a `WeightedMoveEvidence` or packet dictionary in the
[executor schema](../specs/virtual-link-moves.md#127-weighted-evidence).

# Returns
A verification record with validity, reason, and optional reconstructed target,
evidence and error text; success has reason `:verified`.

# Guarantees
1. Success means the current executor's result agrees with the entire packet,
   including raw states and requested equation checks.
   [proved/agent: by construction and envelope comparison]
2. The supplied packet is not modified. [proved/agent: by construction]

# Not guaranteed
Replay is not an independent proof of framing preservation or canonical
manifold identity. Invalid evidence may be malformed or obsolete rather than
mathematically incorrect. Assigned-result packets use `replay_weighted_relative_result`.

# Sources
cf: [virtual-link moves §12.7](../specs/virtual-link-moves.md#127-weighted-evidence) — executor evidence and version scope.

# Examples
```jldoctest
julia> v = lens_space_link(2,1); set_weight!(v, [1,-1,-1,0]); c = only(movable(v, "H"; param=1, vertices=[1])); r = verify_weighted_move(v, c; framed=true);

julia> replay_weighted_move(r.evidence).valid
true
```

Role: condition
"""
function replay_weighted_move(value)
    try
        evidence=weighted_move_evidence_from_data(value isa WeightedMoveEvidence ? value.data : value)
        data=evidence.data
        source=_weighted_link_from_data(data["source_raw"])
        target=_weighted_link_from_data(data["target_raw"])
        weighted_raw_hash(source)==data["source_raw_sha256"] || return _weighted_invalid(:source_hash_mismatch)
        weighted_raw_hash(target)==data["target_raw_sha256"] || return _weighted_invalid(:target_hash_mismatch)
        equation=data["source_equation"]
        framed=_move_data_bool(_move_data_field(equation,"framed_requested"),"framed_requested")
        boundary=data["boundary"]
        decode_boundary(value)=_weighted_decode_targets(value)
        source_boundary=decode_boundary(_move_data_field(boundary,"source"))
        target_boundary=decode_boundary(_move_data_field(boundary,"target"))
        replay=verify_weighted_move(source,_weighted_candidate_from_data(data["candidate"]);
            framed=framed,expected_raw_target=target,
            boundary=(source=source_boundary,target=target_boundary))
        replay.valid || return replay
        ograph_canonical_json(replay.evidence.data)==ograph_canonical_json(data) ||
            return _weighted_invalid(:evidence_mismatch)
        replay
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(:malformed_evidence;error=error)
    end
end
