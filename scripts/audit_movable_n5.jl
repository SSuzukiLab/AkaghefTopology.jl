"""Stream an exhaustive signed-Gauss Reidemeister audit through N=5.

Deletion and R3 witnesses are enumerated for every rooted normalized signed
Gauss diagram and every crossing-orientation decoration.  Insertion witnesses
are enumerated directly only through `--direct-insert-nmax` (default 3).  For
larger N their count is obtained by an exact gap/type factorization:

* with `g` directed rooted gaps, R1 has `4g` witnesses;
* R2 has `8g(g-1)` witnesses on distinct ordered gaps and `16g` witnesses
  on one gap, hence `8g(g+1)` in total.

The factorization is independent of the source code contents; it is not a
claim that a naive hundreds-of-millions insertion loop was run.  Candidate
keys and transition results are streamed into SHA-256 contexts, while only
one diagram's small candidate lists are retained at a time.
"""

using Dates
using JSON
using SHA

const AUDIT_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(AUDIT_ROOT, "src", "knot_gauss.jl"))
include(joinpath(AUDIT_ROOT, "src", "reidemeister.jl"))
include(joinpath(AUDIT_ROOT, "test", "support", "reidemeister_oracle.jl"))

function audit_int_arg(name::AbstractString, default::Int)
    prefix = "--$(name)="
    for argument in ARGS
        startswith(argument, prefix) || continue
        return parse(Int, argument[length(prefix) + 1:end])
    end
    default
end

const AUDIT_NMIN = audit_int_arg("nmin", 0)
const AUDIT_NMAX = audit_int_arg("nmax", 5)
const AUDIT_DIRECT_INSERT_NMAX = audit_int_arg("direct-insert-nmax", 3)
const AUDIT_PROGRESS = "--progress" in ARGS

function audit_production_as_oracle(candidate::ReidemeisterCandidate)
    OracleReidemeisterWitness(candidate.type, candidate.direction,
                              candidate.vertices, candidate.edge_positions,
                              candidate.before, candidate.after,
                              candidate.variant, candidate.source)
end

function audit_production_candidates(d::KnotGaussDiagram,
                                     move_type::AbstractString,
                                     direction)
    direction === nothing ? movable(d, move_type) :
        movable(d, move_type; param=direction)
end

function audit_diagram_text(d::KnotGaussDiagram)
    "gauss=$(join(d.gauss, ','));orientation=$(join(d.orientation, ','))"
end

function audit_candidate_line(move_type::AbstractString,
                              candidate_texts::Vector{String})
    "$(move_type)|$(join(candidate_texts, ";"))\n"
end

function audit_candidate_texts(candidates)
    sort!([oracle_witness_text(candidate) for candidate in candidates])
end

function audit_target_text(d::KnotGaussDiagram)
    "$(join(d.gauss, ','))|$(join(d.orientation, ','))"
end

function audit_increment!(table::Dict{String,Int}, key::AbstractString,
                          amount::Int=1)
    table[String(key)] = get(table, String(key), 0) + amount
end

function audit_first_difference(left::Dict{String,Int},
                                right::Dict{String,Int})
    keys_union = sort!(collect(union(keys(left), keys(right))))
    for key in keys_union
        get(left, key, 0) == get(right, key, 0) ||
            return Dict{String,Any}("key" => key,
                                    "production" => get(left, key, 0),
                                    "oracle" => get(right, key, 0))
    end
    nothing
end

function audit_record_failure!(failures::Vector{Any}, kind, n, d,
                               move_type, details)
    length(failures) >= 20 && return
    push!(failures, Dict{String,Any}(
        "kind" => String(kind),
        "crossing_number" => n,
        "diagram" => audit_diagram_text(d),
        "move_type" => String(move_type),
        "details" => details,
    ))
end

function audit_sha256(context::SHA.SHA256_CTX)
    bytes2hex(SHA.digest!(context))
end

function audit_rule_record(move_type::AbstractString, direction)
    Dict{String,Any}(
        "type" => String(move_type),
        "direction" => direction === nothing ? "self_inverse" : String(direction),
        "diagrams_with_candidates" => 0,
        "candidate_count" => 0,
        "application_checks" => 0,
        "inverse_checks" => 0,
        "oracle_inverse_checks" => 0,
        "boundary_edge_count" => 0,
        "variant_counts" => Dict{String,Int}(),
        "three_sign_counts" => Dict{String,Int}(),
    )
end

function audit_update_rule!(record::Dict{String,Any}, d::KnotGaussDiagram,
                           candidates, oracle_candidates_for_diagram)
    isempty(candidates) || (record["diagrams_with_candidates"] += 1)
    record["candidate_count"] += length(candidates)
    variants = record["variant_counts"]::Dict{String,Int}
    for candidate in candidates
        audit_increment!(variants, String(candidate.variant))
        l = length(d.gauss)
        any(edge -> l > 0 && (edge[1] == l || edge[2] == 1),
            candidate.edge_positions) && (record["boundary_edge_count"] += 1)
        if uppercase(candidate.type) == "R3"
            oracle_data = oracle_r3_data(d, candidate.edge_positions)
            oracle_data === nothing || begin
                three_signs = record["three_sign_counts"]::Dict{String,Int}
                audit_increment!(three_signs, string(first(oracle_data.three_signs)))
            end
        end
    end
end

function audit_pair_candidates(d::KnotGaussDiagram, move_type, direction,
                               production, oracle, failures, n)
    production_oracle = audit_production_as_oracle.(production)
    production_keys = oracle_multiset(production_oracle)
    oracle_keys = oracle_multiset(oracle)
    if production_keys != oracle_keys
        difference = audit_first_difference(production_keys, oracle_keys)
        audit_record_failure!(failures, "candidate_multiset", n, d,
                              move_type, difference)
        return false, production_keys, oracle_keys
    end
    true, production_keys, oracle_keys
end

function audit_transition_checks!(record::Dict{String,Any},
                                  transition_context::SHA.SHA256_CTX,
                                  d::KnotGaussDiagram, move_type, direction,
                                  production, oracle, failures, n;
                                  inverse_enabled::Bool=false)
    by_key = Dict{String,Vector{ReidemeisterCandidate}}()
    for candidate in production
        key = oracle_witness_text(audit_production_as_oracle(candidate))
        push!(get!(by_key, key, ReidemeisterCandidate[]), candidate)
    end
    for oracle_candidate in oracle
        key = oracle_witness_text(oracle_candidate)
        choices = get(by_key, key, ReidemeisterCandidate[])
        isempty(choices) && begin
            audit_record_failure!(failures, "candidate_binding", n, d,
                                  move_type, Dict("key" => key))
            continue
        end
        production_candidate = popfirst!(choices)
        production_target = nothing
        oracle_target = nothing
        try
            production_target = apply_reidemeister(d, production_candidate)
            oracle_target = oracle_apply(d, oracle_candidate)
        catch error
            audit_record_failure!(failures, "application_exception", n, d,
                                  move_type,
                                  Dict("key" => key,
                                       "exception" => sprint(showerror, error)))
            continue
        end
        record["application_checks"] += 1
        transition_line = "$(key)|target=$(audit_target_text(production_target))\n"
        SHA.update!(transition_context, codeunits(transition_line))
        production_target == oracle_target ||
            audit_record_failure!(failures, "application_target", n, d,
                                  move_type,
                                  Dict("key" => key,
                                       "production" => audit_target_text(production_target),
                                       "oracle" => audit_target_text(oracle_target)))
        verify_reidemeister_transition(d, production_target,
                                       production_candidate) ||
            audit_record_failure!(failures, "transition_verifier", n, d,
                                  move_type, Dict("key" => key))

        if inverse_enabled
            inverse_direction = production_candidate.direction == :insert ? :delete :
                                production_candidate.direction == :delete ? :insert :
                                nothing
            inverse_candidates = audit_production_candidates(
                production_target, move_type, inverse_direction)
            oracle_inverse_candidates = oracle_candidates(
                production_target, move_type; direction=inverse_direction)
            inverse_keys = oracle_multiset(
                audit_production_as_oracle.(inverse_candidates))
            oracle_inverse_keys = oracle_multiset(oracle_inverse_candidates)
            inverse_keys == oracle_inverse_keys || audit_record_failure!(
                failures, "inverse_candidate_multiset", n, d, move_type,
                Dict("key" => key,
                     "production" => audit_first_difference(inverse_keys,
                                                            oracle_inverse_keys)))
            inverse_ok = any(inverse_candidate -> begin
                try
                    canonicalize(apply_reidemeister(
                        production_target, inverse_candidate)) == canonicalize(d)
                catch
                    false
                end
            end, inverse_candidates)
            oracle_inverse_ok = any(inverse_candidate -> begin
                try
                    canonicalize(oracle_apply(
                        production_target, inverse_candidate)) == canonicalize(d)
                catch
                    false
                end
            end, oracle_inverse_candidates)
            record["inverse_checks"] += 1
            record["oracle_inverse_checks"] += 1
            inverse_ok || audit_record_failure!(
                failures, "inverse_closure", n, d, move_type,
                Dict("key" => key, "inverse_candidate_count" => length(inverse_candidates)))
            oracle_inverse_ok || audit_record_failure!(
                failures, "oracle_inverse_closure", n, d, move_type,
                Dict("key" => key,
                     "inverse_candidate_count" => length(oracle_inverse_candidates)))
        end
    end
    for (key, choices) in by_key
        isempty(choices) || audit_record_failure!(
            failures, "unmatched_production_candidate", n, d, move_type,
            Dict("key" => key, "remaining" => length(choices)))
    end
end

function audit_direct_insert_record(move_type::AbstractString)
    Dict{String,Any}(
        "type" => String(move_type),
        "mode" => "direct_enumeration",
        "diagram_count" => 0,
        "candidate_count" => 0,
        "application_checks" => 0,
        "inverse_checks" => 0,
        "oracle_inverse_checks" => 0,
        "production_stream_sha256" => nothing,
        "oracle_stream_sha256" => nothing,
        "transition_stream_sha256" => nothing,
        "factorization" => nothing,
    )
end

function audit_factorized_insert_record(move_type::AbstractString)
    Dict{String,Any}(
        "type" => String(move_type),
        "mode" => "factorized_proved",
        "diagram_count" => 0,
        "candidate_count" => 0,
        "application_checks" => nothing,
        "inverse_checks" => nothing,
        "oracle_inverse_checks" => nothing,
        "production_stream_sha256" => nothing,
        "oracle_stream_sha256" => nothing,
        "transition_stream_sha256" => nothing,
        "factorization" => nothing,
    )
end

function audit_insert_factorization(move_type::AbstractString, n::Int,
                                    diagram_count::Int)
    gaps = n == 0 ? 1 : 2n
    if move_type == "R1"
        per_diagram = 4gaps
        formula = "2 endpoint orders * 2 crossing signs * g gaps"
    else
        distinct = 8gaps * (gaps - 1)
        same_gap = 16gaps
        per_diagram = distinct + same_gap
        formula = "8 choices per ordered distinct gap pair + 16 choices per same gap (both HH-TT and TT-HH concatenations)"
    end
    Dict{String,Any}(
        "directed_gap_count" => gaps,
        "per_diagram_count" => per_diagram,
        "diagram_count" => diagram_count,
        "total_candidate_count" => diagram_count * per_diagram,
        "formula" => formula,
        "proof_status" => "proved combinatorial factorization; source-independent",
    )
end

function audit_direct_insertions!(records::Dict{String,Any}, n::Int,
                                  diagrams, failures)
    contexts = Dict{String,SHA.SHA256_CTX}(move_type => SHA.SHA256_CTX()
                                           for move_type in ("R1", "R2"))
    oracle_contexts = Dict{String,SHA.SHA256_CTX}(move_type => SHA.SHA256_CTX()
                                                 for move_type in ("R1", "R2"))
    counts = Dict{String,Int}(move_type => 0 for move_type in ("R1", "R2"))
    transition_contexts = Dict{String,SHA.SHA256_CTX}(
        move_type => SHA.SHA256_CTX() for move_type in ("R1", "R2"))
    diagram_count = 0
    for d in diagrams
        diagram_count += 1
        for move_type in ("R1", "R2")
            production = movable(d, move_type; param=:insert)
            oracle = oracle_candidates(d, move_type; direction=:insert)
            production_keys = audit_candidate_texts(
                audit_production_as_oracle.(production))
            oracle_keys = audit_candidate_texts(oracle)
            production_counts = oracle_multiset(audit_production_as_oracle.(production))
            oracle_counts = oracle_multiset(oracle)
            production_counts == oracle_counts || begin
                audit_record_failure!(failures, "insertion_multiset", n, d,
                                      move_type,
                                      audit_first_difference(production_counts,
                                                             oracle_counts))
            end
            # Direct insertion is the bounded application benchmark as well as
            # a candidate-list comparison.  Every paired witness is applied by
            # both implementations, then its inverse is sought modulo the
            # diagram's cyclic/root and crossing-label equivalence.
            audit_transition_checks!(
                records[move_type], transition_contexts[move_type], d,
                move_type, :insert,
                production, oracle, failures, n; inverse_enabled=true)
            SHA.update!(contexts[move_type], codeunits(
                "$(audit_diagram_text(d))|$(audit_candidate_line(move_type, production_keys))"))
            SHA.update!(oracle_contexts[move_type], codeunits(
                "$(audit_diagram_text(d))|$(audit_candidate_line(move_type, oracle_keys))"))
            counts[move_type] += length(production)
        end
    end
    for move_type in ("R1", "R2")
        record = records[move_type]
        record["diagram_count"] = diagram_count
        record["candidate_count"] = counts[move_type]
        record["production_stream_sha256"] = audit_sha256(contexts[move_type])
        record["oracle_stream_sha256"] = audit_sha256(oracle_contexts[move_type])
        record["transition_stream_sha256"] =
            audit_sha256(transition_contexts[move_type])
        record["factorization"] = audit_insert_factorization(move_type, n,
                                                              diagram_count)
    end
end

function audit_factorized_insertions!(records::Dict{String,Any}, n::Int,
                                      diagram_count::Int)
    for move_type in ("R1", "R2")
        record = records[move_type]
        record["diagram_count"] = diagram_count
        factorization = audit_insert_factorization(move_type, n, diagram_count)
        record["candidate_count"] = factorization["total_candidate_count"]
        record["factorization"] = factorization
    end
end

function audit_payload()
    0 <= AUDIT_NMIN <= AUDIT_NMAX <= 5 ||
        throw(ArgumentError("require 0 <= --nmin <= --nmax <= 5"))
    0 <= AUDIT_DIRECT_INSERT_NMAX <= 5 ||
        throw(ArgumentError("require 0 <= --direct-insert-nmax <= 5"))

    failures = Any[]
    oracle_source = read(joinpath(AUDIT_ROOT, "test", "support",
                                  "reidemeister_oracle.jl"), String)
    forbidden_oracle_tokens = ["_rm_", "movable(", "apply_reidemeister",
                               "verify_reidemeister_candidate"]
    for token in forbidden_oracle_tokens
        occursin(token, oracle_source) && push!(failures, Dict{String,Any}(
            "kind" => "oracle_dependency", "token" => token))
    end

    records = Any[]
    totals = Dict{String,Int}(
        "diagram_count" => 0,
        "expected_diagram_count" => 0,
        "deletion_candidate_count" => 0,
        "r3_candidate_count" => 0,
        "insertion_candidate_count" => 0,
    )
    started_at_utc = Dates.now(Dates.UTC)
    started_ns = time_ns()

    for n in AUDIT_NMIN:AUDIT_NMAX
        diagrams = enumerate_knot_gauss(n)
        expected_diagram_count = rooted_normalized_signed_gc_count(n) * 2^n
        rule_records = Dict{String,Any}(
            "R1" => audit_rule_record("R1", :delete),
            "R2" => audit_rule_record("R2", :delete),
            "R3" => audit_rule_record("R3", nothing),
        )
        insertion_records = Dict{String,Any}(
            "R1" => AUDIT_DIRECT_INSERT_NMAX >= n ?
                audit_direct_insert_record("R1") : audit_factorized_insert_record("R1"),
            "R2" => AUDIT_DIRECT_INSERT_NMAX >= n ?
                audit_direct_insert_record("R2") : audit_factorized_insert_record("R2"),
        )
        rule_contexts = Dict{String,SHA.SHA256_CTX}(move_type => SHA.SHA256_CTX()
                                                   for move_type in ("R1", "R2", "R3"))
        oracle_contexts = Dict{String,SHA.SHA256_CTX}(move_type => SHA.SHA256_CTX()
                                                     for move_type in ("R1", "R2", "R3"))
        transition_context = SHA.SHA256_CTX()
        diagram_count = 0
        for d in diagrams
            diagram_count += 1
            for (move_type, direction) in (("R1", :delete),
                                           ("R2", :delete),
                                           ("R3", nothing))
                production = audit_production_candidates(d, move_type, direction)
                oracle = oracle_candidates(d, move_type; direction=direction)
                matched, production_keys, oracle_keys = audit_pair_candidates(
                    d, move_type, direction, production, oracle, failures, n)
                audit_update_rule!(rule_records[move_type], d, production, oracle)
                SHA.update!(rule_contexts[move_type], codeunits(
                    "$(audit_diagram_text(d))|$(audit_candidate_line(move_type,
                                                                        audit_candidate_texts(audit_production_as_oracle.(production))))"))
                SHA.update!(oracle_contexts[move_type], codeunits(
                    "$(audit_diagram_text(d))|$(audit_candidate_line(move_type,
                                                                        audit_candidate_texts(oracle)))"))
                matched && audit_transition_checks!(
                    rule_records[move_type], transition_context, d, move_type,
                    direction, production, oracle, failures, n;
                    # Insertion enumeration is avoided beyond the direct
                    # benchmark.  Deletion inverses are therefore checked
                    # directly only in the same bounded direct phase; R3 has
                    # no combinatorial insertion blow-up.
                    inverse_enabled=(move_type == "R3" ||
                                     (n <= AUDIT_DIRECT_INSERT_NMAX &&
                                      move_type in ("R1", "R2"))))
            end
            AUDIT_PROGRESS && diagram_count % 100_000 == 0 &&
                println(stderr, "[movable-audit] n=$n diagrams=$diagram_count")
        end

        if n <= AUDIT_DIRECT_INSERT_NMAX
            # A fresh iterator is required: the first iterator was consumed by
            # the deletion/R3 stream above.
            audit_direct_insertions!(insertion_records, n,
                                     enumerate_knot_gauss(n), failures)
        else
            audit_factorized_insertions!(insertion_records, n,
                                          diagram_count)
        end

        for move_type in ("R1", "R2", "R3")
            rule_records[move_type]["production_stream_sha256"] =
                audit_sha256(rule_contexts[move_type])
            rule_records[move_type]["oracle_stream_sha256"] =
                audit_sha256(oracle_contexts[move_type])
        end
        expected_diagram_count == diagram_count || push!(failures,
            Dict{String,Any}("kind" => "enumeration_count", "n" => n,
                             "expected" => expected_diagram_count,
                             "visited" => diagram_count))
        totals["diagram_count"] += diagram_count
        totals["expected_diagram_count"] += expected_diagram_count
        totals["deletion_candidate_count"] +=
            sum(rule_records[move_type]["candidate_count"] for move_type in ("R1", "R2"))
        totals["r3_candidate_count"] += rule_records["R3"]["candidate_count"]
        totals["insertion_candidate_count"] +=
            sum(insertion_records[move_type]["candidate_count"]
                for move_type in ("R1", "R2"))
        push!(records, Dict{String,Any}(
            "crossing_number" => n,
            "expected_diagram_count" => expected_diagram_count,
            "visited_diagram_count" => diagram_count,
            "deletion_and_r3" => rule_records,
            "insertion" => insertion_records,
        ))
    end

    elapsed_seconds = (time_ns() - started_ns) / 1e9
    Dict{String,Any}(
        "schema_version" => 1,
        "status" => isempty(failures) ? "pass" : "fail",
        "campaign" => "signed Gauss Reidemeister movable oracle audit",
        "scope" => Dict{String,Any}(
            "included" => ["every rooted normalized signed Gauss diagram",
                           "every crossing-orientation decoration",
                           "R1/R2 deletion", "R3 self-inverse",
                           "direct R1/R2 insertion benchmark through configured bound"],
            "excluded" => ["classical realizability", "isotopy classification",
                           "naive large-N insertion loop"],
            "crossing_number_range" => [AUDIT_NMIN, AUDIT_NMAX],
            "direct_insertion_nmax" => AUDIT_DIRECT_INSERT_NMAX,
        ),
        "conventions" => Dict{String,Any}(
            "positive_visit" => "head/over endpoint",
            "negative_visit" => "tail/under endpoint",
            "gap" => "directed rooted edge (i,next(i)); (0,0) for empty code",
            "R1" => "one adjacent head/tail pair",
            "R2" => "adjacent HH and TT pairs with opposite crossing orientations",
            "R3" => "matched HH/TT/HT arcs with equal sign*parity*direction",
            "R3_application" => "swap endpoints on each selected arc",
        ),
        "enumeration" => Dict{String,Any}(
            "diagram_count_formula" => "(2n)!/n! * 2^n",
            "rooted_code_count_formula" => "(2n)!/n!",
            "streaming" => true,
            "retained_candidate_scope" => "one diagram at a time",
        ),
        "insertion_factorization" => Dict{String,Any}(
            "R1" => "4g",
            "R2" => "8g(g-1)+16g = 8g(g+1)",
            "same_gap_note" => "The extra 8 same-gap witnesses are TT-HH concatenations; HH-TT contributes the first 8.",
            "independence" => "depends only on the directed gap count and local type/order/sign choices, not source helper functions",
            "large_N_statement" => "No naive 800M insertion run is claimed; n>direct-insert-nmax rows use this exact factorization.",
        ),
        "records" => records,
        "totals" => totals,
        "coverage" => Dict{String,Any}(
            "near_miss_fixtures" => ["unmatched R3 kinds", "unequal R3 3-signs",
                                     "all HT arcs", "matched shape with wrong 3-signs"],
            "equivariance" => "cyclic shift and label relabel tests are in test/movable_oracle.jl",
            "source_binding" => "stale and forged witness tests are in test/movable_oracle.jl",
            "inverse_closure" => "direct audit for R3 and direct-insertion bound for R1/R2; representative full suite in test/movable_oracle.jl",
        ),
        "runtime" => Dict{String,Any}(
            "command" => "julia --startup-file=no --pkgimages=no --project=. scripts/audit_movable_n5.jl --nmin=$(AUDIT_NMIN) --nmax=$(AUDIT_NMAX) --direct-insert-nmax=$(AUDIT_DIRECT_INSERT_NMAX)",
            "julia_version" => string(VERSION),
            "kernel" => string(Sys.KERNEL),
            "machine" => string(Sys.MACHINE),
            "started_at_utc" => string(started_at_utc),
            "elapsed_seconds" => elapsed_seconds,
        ),
        "failures" => failures,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    payload = audit_payload()
    JSON.print(stdout, payload, 2)
    println()
    payload["status"] == "pass" || exit(1)
end
