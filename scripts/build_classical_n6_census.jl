"""Build the direct classical rooted signed Gauss census through `n == 6`.

The N=6 lane never iterates the repository's 42,577,920 arbitrary signed
inputs.  It visits each rooted unsigned double-occurrence word once, checks
all local alternating rotations, and decorates only spherical immersions.
Progress is written to stderr; machine-readable counts and fingerprints are
written to stdout.

The per-n raw fingerprint hashes sorted lines
`n=<n>|raw|gauss=<signed visits>;orientation=<crossing signs>`, and the
per-n canonical fingerprint hashes sorted lines
`n=<n>|canonical|<KnotGaussDiagram canonical_serialization>`; each line ends
with a newline.  The `_all_n` fingerprints concatenate those same line streams
in increasing n order.

Optional `--audit-exhaustive-through=5` performs the expensive independent
comparison against `classical_gauss_realizable_fast` only through the requested
boundary.  There is deliberately no exhaustive N=6 filter in this script.
"""

using SHA

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_realizability.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_enumeration.jl"))

function _cge_script_update!(context, line::AbstractString)
    # SHA.update! handles CodeUnits linearly; avoid SHA.sha256(CodeUnits),
    # whose Julia 1.12 generic one-shot path is pathological on large Strings.
    SHA.update!(context, codeunits(string(line, '\n')))
    nothing
end

function _cge_script_raw_serialization(key::UInt64, n::Int)
    code, orientation = classical_gauss_key_to_code(key, n)
    "gauss=$(join(code, ','));orientation=$(join(orientation, ','))"
end

"""Hash and canonicalize one direct raw set, without retaining a census file."""
function _cge_script_summarize(result::ClassicalGaussEnumerationResult,
                               raw_all_context, canonical_all_context)
    n = result.n
    raw_keys = sort!(collect(result.raw_keys))
    raw_context = SHA.SHA256_CTX()
    canonical_set = Set{String}()

    for key in raw_keys
        raw_serialization = _cge_script_raw_serialization(key, n)
        raw_line = "n=$(n)|raw|$(raw_serialization)"
        _cge_script_update!(raw_context, raw_line)
        _cge_script_update!(raw_all_context, raw_line)
        diagram = classical_gauss_key_to_diagram(key, n)
        push!(canonical_set, canonical_serialization(diagram))
    end

    canonical_context = SHA.SHA256_CTX()
    canonical_serializations = sort!(collect(canonical_set))
    for serialization in canonical_serializations
        canonical_line = "n=$(n)|canonical|$(serialization)"
        _cge_script_update!(canonical_context, canonical_line)
        _cge_script_update!(canonical_all_context, canonical_line)
    end

    (raw_classical=length(raw_keys),
     canonical_classical=length(canonical_serializations),
     canonical_state=canonical_set,
     raw_fingerprint=bytes2hex(SHA.digest!(raw_context)),
     canonical_fingerprint=bytes2hex(SHA.digest!(canonical_context)))
end

function _cge_script_audit_one(n::Int, direct_keys::Set{UInt64})
    exhaustive_keys = Set{UInt64}()
    rooted_total = 0
    for diagram in enumerate_knot_gauss(n; orientations=true)
        rooted_total += 1
        classical_gauss_realizable_fast(diagram) || continue
        push!(exhaustive_keys,
              _cge_pack_signed(diagram.gauss, diagram.orientation))
    end

    direct_only = setdiff(direct_keys, exhaustive_keys)
    exhaustive_only = setdiff(exhaustive_keys, direct_keys)
    symmetric_difference = length(direct_only) + length(exhaustive_only)
    symmetric_difference == 0 ||
        throw(AssertionError(
            "direct/exhaustive disagreement at n=$n: direct_only=$(length(direct_only)) " *
            "exhaustive_only=$(length(exhaustive_only))"))
    (n=n, rooted_total=rooted_total,
     direct_raw=length(direct_keys), exhaustive_expected=length(exhaustive_keys),
     symmetric_difference=symmetric_difference)
end

function _cge_script_parse_args(args)
    max_n = 6
    audit_through = nothing
    for arg in args
        if startswith(arg, "--max-n=")
            max_n = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--audit-exhaustive-through=")
            audit_through = parse(Int, split(arg, "=", limit=2)[2])
        else
            throw(ArgumentError(
                "unsupported argument $arg; use --max-n=N and/or " *
                "--audit-exhaustive-through=N"))
        end
    end
    0 <= max_n <= 6 || throw(ArgumentError("max_n must lie in 0:6"))
    if audit_through !== nothing
        0 <= audit_through <= 5 ||
            throw(ArgumentError("exhaustive audit boundary must lie in 0:5"))
        audit_through <= max_n ||
            throw(ArgumentError("audit boundary cannot exceed max_n"))
    end
    max_n, audit_through
end

"""Run direct generation and optional N<=5 exhaustive equivalence audit."""
function build_classical_n6_census(max_n::Integer=6;
                                   audit_exhaustive_through=nothing)
    maximum_n = _cge_validate_n(max_n)
    audit = audit_exhaustive_through === nothing ? nothing :
        _cge_validate_n(audit_exhaustive_through)
    audit !== nothing && audit > 5 &&
        throw(ArgumentError("exhaustive audit is restricted to n<=5"))
    audit !== nothing && audit > maximum_n &&
        throw(ArgumentError("audit boundary cannot exceed max_n"))

    rows = NamedTuple[]
    audit_rows = NamedTuple[]
    canonical_sets = Set{String}[]
    raw_all_context = SHA.SHA256_CTX()
    canonical_all_context = SHA.SHA256_CTX()
    cumulative_words = 0
    cumulative_spheres = 0
    cumulative_generated = 0
    cumulative_duplicates = 0
    cumulative_raw = 0
    cumulative_canonical = 0
    started_all = time_ns()

    for n in 0:maximum_n
        started_n = time_ns()
        expected_words = rooted_unsigned_gauss_count(n)
        progress = function (words, spheres, generated)
            words % 1000 == 0 || return nothing
            println(stderr,
                    "progress n=$(n) unsigned_words=$(words)/$(expected_words) " *
                    "spherical_immersions=$(spheres) generated=$(generated)")
            flush(stderr)
            nothing
        end
        result = direct_classical_gauss_enumeration(n; progress=progress)
        result.unsigned_words == expected_words ||
            throw(AssertionError("rooted unsigned word count mismatch at n=$n"))
        generated_seconds = (time_ns() - started_n) / 1e9

        if audit !== nothing && n <= audit
            push!(audit_rows, _cge_script_audit_one(n, result.raw_keys))
        end

        summary = _cge_script_summarize(result, raw_all_context,
                                        canonical_all_context)
        push!(canonical_sets, summary.canonical_state)
        cumulative_words += result.unsigned_words
        cumulative_spheres += result.spherical_immersions
        cumulative_generated += result.generated_decorations
        cumulative_duplicates += result.raw_duplicates
        cumulative_raw += summary.raw_classical
        cumulative_canonical += summary.canonical_classical
        elapsed_seconds = (time_ns() - started_n) / 1e9
        push!(rows,
              (n=n,
               unsigned_words=result.unsigned_words,
               rotations_checked=result.rotations_checked,
               spherical_immersions=result.spherical_immersions,
               generated_decorations=result.generated_decorations,
               raw_duplicates=result.raw_duplicates,
               raw_classical=summary.raw_classical,
               canonical_classical=summary.canonical_classical,
               generation_seconds=generated_seconds,
               elapsed_seconds=elapsed_seconds,
               raw_fingerprint=summary.raw_fingerprint,
               canonical_fingerprint=summary.canonical_fingerprint,
               cumulative_words=cumulative_words,
               cumulative_spherical_immersions=cumulative_spheres,
               cumulative_generated_decorations=cumulative_generated,
               cumulative_raw_duplicates=cumulative_duplicates,
               cumulative_raw_classical=cumulative_raw,
               cumulative_canonical_classical=cumulative_canonical))
    end

    (rows=rows,
     audit_rows=audit_rows,
     canonical_sets=canonical_sets,
     raw_fingerprint=bytes2hex(SHA.digest!(raw_all_context)),
     canonical_fingerprint=bytes2hex(SHA.digest!(canonical_all_context)),
     elapsed_seconds=(time_ns() - started_all) / 1e9)
end

function main(args=ARGS)
    max_n, audit = _cge_script_parse_args(args)
    result = build_classical_n6_census(max_n;
                                      audit_exhaustive_through=audit)
    println("classical_gauss_enumeration_version=$(CLASSICAL_GAUSS_ENUMERATION_VERSION)")
    println("julia_version=$(VERSION)")
    println("scope=max_n=$(max_n);direct_planar_shadow_generation;one_component")
    println("convention=+i_over_tail;-i_under_head;orientation=undercrossing_exponent")
    audit_text = audit === nothing ? "none" : string(audit)
    println("boolean_validation=every_decoration_inherits_an_accepted_spherical_rotation;" *
            "optional_fast_audit_through=$(audit_text)")
    println("signed_theorem=Kurlin_Definition_3.1;Lemma_3.2;$(KURLIN_GAUSS_PARAGRAPH_URL)")
    println("unsigned_theorem=Marx_Theorem_2;$(MARX_GAUSS_REALIZABILITY_URL)")
    println("identity=KnotGaussDiagram_canonicalization_v$(KNOT_GAUSS_CANON_VERSION);rotations+first-occurrence-relabel_only")
    println("n unsigned_words rotations_checked spherical_immersions generated_decorations " *
            "raw_duplicates raw_classical canonical_classical generation_seconds " *
            "elapsed_seconds cumulative_words cumulative_spherical_immersions " *
            "cumulative_generated_decorations cumulative_raw_duplicates " *
            "cumulative_raw_classical cumulative_canonical_classical raw_fingerprint " *
            "canonical_fingerprint")
    for row in result.rows
        println(row.n, ' ', row.unsigned_words, ' ', row.rotations_checked, ' ',
                row.spherical_immersions, ' ', row.generated_decorations, ' ',
                row.raw_duplicates, ' ', row.raw_classical, ' ',
                row.canonical_classical, ' ', row.generation_seconds, ' ',
                row.elapsed_seconds, ' ', row.cumulative_words, ' ',
                row.cumulative_spherical_immersions, ' ',
                row.cumulative_generated_decorations, ' ',
                row.cumulative_raw_duplicates, ' ', row.cumulative_raw_classical,
                ' ', row.cumulative_canonical_classical, ' ',
                row.raw_fingerprint, ' ', row.canonical_fingerprint)
    end
    for row in result.audit_rows
        println("audit_n=$(row.n);rooted_total=$(row.rooted_total);direct_raw=$(row.direct_raw);" *
                "exhaustive_expected=$(row.exhaustive_expected);" *
                "symmetric_difference=$(row.symmetric_difference)")
    end
    println("raw_fingerprint_sha256_all_n=$(result.raw_fingerprint)")
    println("canonical_fingerprint_sha256_all_n=$(result.canonical_fingerprint)")
    println("elapsed_seconds=$(result.elapsed_seconds)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
