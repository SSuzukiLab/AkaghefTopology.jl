"""Raw integer-weighted move evidence, separate from the unweighted OGraph schema.

This envelope certifies deterministic executor replay and records equation checks.
It is not an independent mathematical proof of framing transport. Its hashes bind
the full raw presentation and integer weights; they are NOT canonical state IDs.
"""
const WEIGHTED_MOVE_EVIDENCE_SCHEMA = "akaghef.weighted-move-evidence/1"
const WEIGHTED_MOVE_VERIFIER_VERSION = "weighted-executor-replay/1"
const WEIGHTED_MOVE_SOURCE_RULE_VERSION = "integer-mp-h-02-ps-runtime/1"

struct WeightedMoveEvidence
    data::Dict{String,Any}
end

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

function _weighted_exact_int(value)
    value isa Real && !(value isa Bool) && isfinite(value) && isinteger(value) ||
        throw(ArgumentError("weighted evidence requires finite integer weights"))
    abs(big(value)) <= big(2)^53-1 || throw(ArgumentError(
        "weighted evidence exceeds the consecutive exact Float64 integer range"))
    Int(value)
end

function _weighted_raw_data(v::VirtualLink)
    v.is_weighted || throw(ArgumentError("weighted evidence requires assigned weights"))
    _validate_gauss(v.gauss,v.orientation;real_only=v.real_only)
    all(sign -> sign in (-1,0,1),v.orientation) ||
        throw(ArgumentError("invalid crossing orientations"))
    length(v.weights)==length(v.gauss)==length(v.strand_types) ||
        throw(DimensionMismatch("raw component metadata does not match Gauss components"))
    for (word,weights) in zip(v.gauss,v.weights)
        length(weights)==count(x -> x!=0 && v.orientation[abs(x)]!=0,word) ||
            throw(DimensionMismatch("raw weights must follow real crossing occurrences"))
    end
    (gauss=deepcopy(v.gauss), orientation=copy(v.orientation),
     weights=[_weighted_exact_int.(w) for w in v.weights], real_only=v.real_only,
     strand_types=string.(v.strand_types), is_weighted=true)
end

"Hash a raw weighted presentation, without cyclic/relabeling canonicalization."
weighted_raw_hash(v::VirtualLink) = _weighted_digest((
    schema="akaghef.raw-integer-weighted-virtual-link/1",state=_weighted_raw_data(v)))

function _weighted_link_from_data(value)
    _move_data_object(value,("gauss","orientation","weights","real_only","strand_types","is_weighted"))
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

function _weighted_equation_data(v;framed)
    # Accumulate in BigInt independently of the machine-Int public residual.
    weights=BigInt[_weighted_exact_int(e.weight) for e in real_edge_table(v)]
    disks=disk_table(v)
    residual=2*BigInt.(disks.delta)*weights+BigInt.(disks.cp2)
    closed=all(!=(0),v.orientation) && all(w->all(!=(0),w),v.gauss) && is_closed(v)[1]
    framed && !closed && throw(ArgumentError("framed evidence requires a closed normal o-graph"))
    framed && any(!=(0),residual) && throw(ArgumentError("nonzero framing residual"))
    # Strings retain arbitrary-size exact residuals in JSON without overflow.
    (framed_requested=framed,closed_normal=closed,residual=string.(residual),
     real_edge_weights=Int.(weights))
end

function _weighted_unweight(v)
    result=deepcopy(v)
    result.is_weighted=false
    result.weights=[zeros(length(w)) for w in v.weights]
    result
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

"""Apply a frozen weighted candidate on a copy, check it, and produce raw evidence.

`framed=true` checks closed-normal source and target plus the exact framing
equation. Residual zero alone does not prove preservation of the same framing.
Unsupported/unmatched moves return an invalid report without modifying source.
"""
function verify_weighted_move(source::VirtualLink,candidate::MoveCandidate;
                              framed=false,expected_raw_target=nothing)
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
        raw_source=_weighted_raw_data(source)
        source_equation=_weighted_equation_data(source;framed=framed)
        stage=:candidate_mismatch
        frozen=_weighted_catalogue_candidate(source,candidate)
        frozen===nothing && return _weighted_invalid(:candidate_mismatch)
        stage=:application_failed
        details = if frozen.type=="PS"
            plan=_weighted_ps_plan(source,frozen)
            (coordinate_system="real-gauss-components/head-occurrences",
             ps=plan.detail,transport=_weighted_transport_detail(source,plan.candidate))
        else
            (coordinate_system="real-gauss-components/head-occurrences",
             ps=nothing,transport=_weighted_transport_detail(source,frozen))
        end
        target=move(source,frozen)
        stage=:invalid_target
        raw_target=_weighted_raw_data(target)
        target_equation=_weighted_equation_data(target;framed=framed)
        if expected_raw_target!==nothing
            expected=expected_raw_target isa VirtualLink ? _weighted_raw_data(expected_raw_target) : expected_raw_target
            ograph_canonical_json(expected)==ograph_canonical_json(raw_target) ||
                return _weighted_invalid(:raw_target_mismatch)
        end
        stage=:topological_verification_failed
        topological=_weighted_topological_data(source,target,frozen)
        stage=:evidence_encoding_failed
        data=_weighted_json_copy((schema=WEIGHTED_MOVE_EVIDENCE_SCHEMA,
            verifier_version=WEIGHTED_MOVE_VERIFIER_VERSION,
            source_rule_version=WEIGHTED_MOVE_SOURCE_RULE_VERSION,
            evidence_kind="executor-replay-with-equation-checks",identity="raw-weighted-only",
            source_raw=raw_source,target_raw=raw_target,
            source_raw_sha256=weighted_raw_hash(source),target_raw_sha256=weighted_raw_hash(target),
            candidate=_weighted_candidate_data(frozen),source_equation=source_equation,
            target_equation=target_equation,topological=topological,
            transport_detail=_weighted_detail_data(details)))
        WeightedMoveVerification(true,:verified,WeightedMoveEvidence(data),target,nothing)
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(stage;error=error)
    end
end

"Return a detached JSON-compatible envelope; no unweighted certificate is changed."
weighted_move_evidence_data(evidence::WeightedMoveEvidence)=deepcopy(evidence.data)

"Decode the separate versioned envelope. Replay is required before trusting it."
function weighted_move_evidence_from_data(value)
    _move_data_object(value,("schema","verifier_version","source_rule_version","evidence_kind",
        "identity","source_raw","target_raw","source_raw_sha256","target_raw_sha256",
        "candidate","source_equation","target_equation","topological","transport_detail"))
    field(key)=_move_data_field(value,key)
    field("schema")==WEIGHTED_MOVE_EVIDENCE_SCHEMA || throw(ArgumentError("unsupported weighted evidence schema"))
    field("verifier_version")==WEIGHTED_MOVE_VERIFIER_VERSION || throw(ArgumentError("unsupported weighted verifier version"))
    field("source_rule_version")==WEIGHTED_MOVE_SOURCE_RULE_VERSION || throw(ArgumentError("unsupported weighted source rule version"))
    _weighted_link_from_data(field("source_raw"))
    _weighted_link_from_data(field("target_raw"))
    _weighted_candidate_from_data(field("candidate"))
    WeightedMoveEvidence(_weighted_json_copy(value))
end

"""Recompute all evidence from its raw source and compare the entire envelope.

Saved target, residuals, PS totals, template rules, and topological witnesses are
never used as trusted computed inputs. Old unweighted witnesses keep their own
parser and replay path. This is executor replay, not a second transport proof.
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
        replay=verify_weighted_move(source,_weighted_candidate_from_data(data["candidate"]);
            framed=framed,expected_raw_target=target)
        replay.valid || return replay
        ograph_canonical_json(replay.evidence.data)==ograph_canonical_json(data) ||
            return _weighted_invalid(:evidence_mismatch)
        replay
    catch error
        error isa InterruptException && rethrow()
        _weighted_invalid(:malformed_evidence;error=error)
    end
end
