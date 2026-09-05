"""Build the bounded N<=5 full signed-Gauss classical census.

Run from `lib/topology.jl` with:

    julia --startup-file=no --pkgimages=no --project=. scripts/build_classical_n5_census.jl

The loop uses the repository's `enumerate_knot_gauss` iterator, including all
`2^n` crossing-orientation assignments.  Canonical keys below are a direct
allocation-light transcription of `canonicalize`/`canonical_serialization`:
cyclic shifts and first-occurrence label normalization are the only quotient.
Mirror and traversal reversal are deliberately not identified.
"""

using SHA

include(joinpath(@__DIR__, "..", "src", "knot_gauss.jl"))
include(joinpath(@__DIR__, "..", "src", "classical_gauss_realizability.jl"))

"""Return the repository identity key without constructing a temporary object."""
function _repository_identity_key(code, orientation)
    raw_code = Int[x for x in code]
    signs = Int[x for x in orientation]
    length(raw_code) == 0 && return ()
    length(raw_code) % 2 == 0 || error("invalid even-length code invariant")
    n = length(signs)
    length(raw_code) == 2n || error("invalid orientation-length invariant")
    best = nothing
    for start in 1:length(raw_code)
        labels = Dict{Int,Int}()
        canonical_code = Vector{Int}(undef, length(raw_code))
        canonical_orientation = Vector{Int}(undef, n)
        next_label = 1
        for offset in 0:length(raw_code)-1
            position = mod1(start + offset, length(raw_code))
            value = raw_code[position]
            old_label = abs(value)
            if !haskey(labels, old_label)
                labels[old_label] = next_label
                canonical_orientation[next_label] = signs[old_label]
                next_label += 1
            end
            canonical_code[offset + 1] = sign(value) * labels[old_label]
        end
        candidate = Tuple(vcat(canonical_code, 0, canonical_orientation))
        (best === nothing || candidate < best) && (best = candidate)
    end
    best
end

function _serialize_identity_key(key, n::Int)
    n == 0 && return "KnotGaussDiagram;canon=v$(KNOT_GAUSS_CANON_VERSION);n=0;gauss=;orientation="
    code = key[1:2n]
    orientation = key[(2n + 2):length(key)]
    "KnotGaussDiagram;canon=v$(KNOT_GAUSS_CANON_VERSION);n=$(n);" *
        "gauss=$(join(code, ','));orientation=$(join(orientation, ','))"
end

function _max_n_from_args(args)
    max_n = 5
    for arg in args
        startswith(arg, "--max-n=") ||
            throw(ArgumentError("unsupported argument $arg; use --max-n=N"))
        max_n = parse(Int, split(arg, "=", limit=2)[2])
    end
    0 <= max_n <= 5 || throw(ArgumentError("max_n must lie in 0:5"))
    max_n
end

"""Enumerate and classify all repository-rooted signed diagrams through `max_n`."""
function build_classical_n5_census(max_n::Integer=5)
    max_n = Int(max_n)
    0 <= max_n <= 5 || throw(ArgumentError("max_n must lie in 0:5"))
    rows = NamedTuple[]
    record_lines = String[]
    elapsed_ns = UInt128(0)
    progress = get(ENV, "CLASSICAL_GAUSS_PROGRESS", "") == "1"

    for n in 0:max_n
        started = time_ns()
        rooted_total = 0
        rooted_classical = 0
        all_keys = Set{Tuple}()
        classical_keys = Set{Tuple}()
        # `enumerate_knot_gauss` is the authoritative repository iterator.  It
        # already supplies one normalized signed Gauss code for each rooted
        # visit word and all crossing orientations in deterministic order.
        for diagram in enumerate_knot_gauss(n; orientations=true)
            rooted_total += 1
            key = _repository_identity_key(diagram.gauss, diagram.orientation)
            push!(all_keys, key)
            if classical_gauss_realizable_fast(diagram)
                rooted_classical += 1
                push!(classical_keys, key)
            end
        end
        all_serialized = sort!([_serialize_identity_key(key, n) for key in all_keys])
        classical_serialized =
            sort!([_serialize_identity_key(key, n) for key in classical_keys])
        virtual_serialized = setdiff(all_serialized, classical_serialized)
        append!(record_lines, ["n=$(n)|classical|$(serial)" for serial in classical_serialized])
        append!(record_lines, ["n=$(n)|virtual|$(serial)" for serial in virtual_serialized])
        elapsed = time_ns() - started
        elapsed_ns += UInt128(elapsed)
        push!(rows, (n=n,
                     rooted_total=rooted_total,
                     rooted_classical=rooted_classical,
                     rooted_virtual=rooted_total - rooted_classical,
                     canonical_total=length(all_serialized),
                     canonical_classical=length(classical_serialized),
                     canonical_virtual=length(virtual_serialized),
                     elapsed_ns=UInt128(elapsed)))
        progress && (println(stderr, "completed n=$(n) in $(elapsed / 1e9)s"); flush(stderr))
    end

    # `record_lines` is already deterministic: rows are appended in increasing
    # n, with each class's canonical serializations sorted before appending.
    # Avoid a second global string sort over the N=5 record set.
    progress && (println(stderr, "assembled records=$(length(record_lines))"); flush(stderr))
    payload = isempty(record_lines) ? "" : join(record_lines, "\n") * "\n"
    progress && (println(stderr, "joined payload bytes=$(sizeof(payload))"); flush(stderr))
    # Materialize the UTF-8 bytes: on Julia 1.12 the SHA stdlib's generic
    # CodeUnits path can spend minutes traversing a large String, whereas the
    # concrete byte vector is linear in the 9.8 MB N=5 payload.
    fingerprint = bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(payload))))
    progress && (println(stderr, "hashed payload"); flush(stderr))
    (rows=rows,
     records=record_lines,
     payload_bytes=sizeof(payload),
     fingerprint=fingerprint,
     elapsed_ns=elapsed_ns)
end

function main(args=ARGS)
    max_n = _max_n_from_args(args)
    started = time_ns()
    result = build_classical_n5_census(max_n)
    total_elapsed = time_ns() - started

    println("classical_gauss_census_version=$(CLASSICAL_GAUSS_REALIZABILITY_VERSION)")
    println("scope=max_n=$(max_n);one_component;rooted_normalized_signed;all_crossing_orientations")
    println("identity=KnotGaussDiagram canonicalization v$(KNOT_GAUSS_CANON_VERSION);rotations+first-occurrence-relabel only")
    println("signed_theorem=Kurlin Definition 3.1;Lemma 3.2;Theorem 3.6;$(KURLIN_GAUSS_PARAGRAPH_URL)")
    println("unsigned_theorem=Marx Theorem 2;$(MARX_GAUSS_REALIZABILITY_URL)")
    println("n rooted_total rooted_classical rooted_virtual canonical_total canonical_classical canonical_virtual elapsed_seconds")
    for row in result.rows
        println(row.n, ' ', row.rooted_total, ' ', row.rooted_classical, ' ',
                row.rooted_virtual, ' ', row.canonical_total, ' ',
                row.canonical_classical, ' ', row.canonical_virtual, ' ',
                row.elapsed_ns / 1e9)
    end
    println("records=$(length(result.records));payload_bytes=$(result.payload_bytes)")
    println("fingerprint_sha256=$(result.fingerprint)")
    println("elapsed_seconds=$(total_elapsed / 1e9)")
    nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
