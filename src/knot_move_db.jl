"""Area-2 adapter between signed Gauss diagrams and the verified move store."""

using .MoveGraphStore

const KNOT_MOVE_ADAPTER_VERSION = 1

function _knot_payload(d::KnotGaussDiagram)
    Dict("gauss" => collect(d.gauss), "orientation" => collect(d.orientation))
end

"Construct the Area-2 representation state used by the persistent move graph." 
function knot_canonical_state(d::KnotGaussDiagram)
    c = canonicalize(d)
    CanonicalState(canonical_serialization(c); area=AREA_2, representation="KnotGauss",
                   payload=_knot_payload(c), complexity=Float64(crossing_count(c)))
end

function _candidate_metadata(before, after, candidate)
    Dict("adapter_version" => KNOT_MOVE_ADAPTER_VERSION,
         "theorem" => "Reidemeister local witness with independent source-bound verifier",
         "move_type" => candidate.type,
         "direction" => String(candidate.direction),
         "variant" => String(candidate.variant),
         "vertices" => collect(candidate.vertices),
         "edge_positions" => [collect(edge) for edge in candidate.edge_positions],
         "local_before" => collect(candidate.before),
         "local_after" => collect(candidate.after),
         "source_serialization" => canonical_serialization(before),
         "target_serialization" => canonical_serialization(after))
end

"Apply and independently verify one source-bound candidate, then persist it." 
function persist_reidemeister_transition!(store::MoveGraphStore.MoveGraphStoreDB,
                                          before::KnotGaussDiagram,
                                          candidate::ReidemeisterCandidate)
    candidate.source == before || throw(ArgumentError("candidate is stale for before diagram"))
    after = apply_reidemeister(before, candidate)
    verify_reidemeister_transition(before, after, candidate) ||
        throw(ArgumentError("independent Reidemeister verification failed"))
    source = knot_canonical_state(before)
    target = knot_canonical_state(after)
    metadata = _candidate_metadata(before, after, candidate)
    verifier = (s, t, m) -> verify_reidemeister_transition(before, after, candidate)
    direction = candidate.direction == :insert ? :forward :
                candidate.direction == :delete ? :inverse : :forward
    MoveGraphStore.insert_verified_transition!(store, source, target, candidate.type;
        direction=direction,
        location=metadata["edge_positions"], parameters=metadata,
        engine_version="Topology/$(KNOT_MOVE_ADAPTER_VERSION)",
        verifier_version="Reidemeister/$(KNOT_MOVE_ADAPTER_VERSION)",
        verifier=verifier, verified=true,
        certificate_paths=[metadata])
end

"Persist at most `store_cap` verified one-step neighbors of `before`." 
function expand_reidemeister_neighbors!(store::MoveGraphStore.MoveGraphStoreDB,
                                        before::KnotGaussDiagram; store_cap::Integer,
                                        max_crossings::Integer=typemax(Int))
    store_cap >= 0 || throw(ArgumentError("store_cap must be nonnegative"))
    max_crossings >= 0 || throw(ArgumentError("max_crossings must be nonnegative"))
    attempted = persisted = rejected = 0
    seen = Set{String}()
    for typ in ("R1", "R2", "R3")
        for candidate in movable(before, typ)
            attempted += 1
            after = try
                apply_reidemeister(before, candidate)
            catch
                rejected += 1
                continue
            end
            crossing_count(after) <= max_crossings || continue
            sid = knot_canonical_state(after).state_id
            sid in seen && continue
            push!(seen, sid)
            persisted >= store_cap && return (attempted=attempted, persisted=persisted,
                                              rejected=rejected, capped=true)
            try
                persist_reidemeister_transition!(store, before, candidate)
                persisted += 1
            catch
                rejected += 1
            end
        end
    end
    (attempted=attempted, persisted=persisted, rejected=rejected, capped=false)
end
