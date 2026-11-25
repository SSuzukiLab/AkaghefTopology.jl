# Source: docs/specs/lickorish-code.md §3–4;
# lib/topology.jl/scripts/heegaard/lc_input.py::parse.
"""Lickorish Code (LC): words and placement metadata without Heegaard roles."""

const LC_META = (:place, :disk, :point, :weight)

"""
    LickorishCode(genus, words, places, disks, points, weights, source_path) -> LickorishCode

Store a curve arrangement's Lickorish words and placement metadata.

# Arguments
`genus::Int`, `words::Dict{String,Vector{Int}}`,
`places::Dict{String,String}`, `disks::Dict{String,Vector{String}}`,
`points::Dict{String,String}`, `weights::Dict{String,Int}`, and
`source_path::Union{Nothing,String}`, or convertible values. Conversion
errors propagate. This constructor does not check
[LC validity](../specs/lickorish-code.md) or the value of the genus.

# Returns
`LickorishCode` with genus, named cyclic words, placement references, disk,
point and weight dictionaries, and optional source path. The mathematical
arrangement and reference surface are defined by
[the LC specification](../specs/lickorish-code.md).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed containers
   are retained and not modified. [proved/agent: by construction]
2. The value stores the supplied LC data without assigning Heegaard roles
   or realizing a surface. Construction does not read `source_path`.
   [proved/agent: by construction]

# Conventions
Use [LC passage signs and reference-surface conventions](../conventions/HeegaardCodes.md).

# Not guaranteed
Arbitrary words and places need not realize a valid curve arrangement on
the stated surface. The constructor accepts legacy annotation dictionaries;
it does not transport or validate them. HC roles, disk and mark semantics
belong to [HC annotations](../specs/heegaard-code.md#3-hc-に置く追加情報).
Shared dictionaries and word vectors remain mutable.

# Sources
cf: lib/topology.jl/docs/specs/lickorish-code.md §§3–4 — Lickorish words, placement and their represented surface arrangement.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> d = H.LickorishCode(1, Dict("l1" => [-1,2,-3], "m1" => [1], "c" => [-2,3]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> length(H.movable(d, "lc.isotopy.digon.remove.00"))
1
```

See also [`read_lc`](@ref), [`validate_lc_input`](@ref), [`realize_lc`](@ref), [`resolve_hc`](@ref).

Role: construct
"""
struct LickorishCode
    genus::Int
    words::Dict{String, Vector{Int}}
    places::Dict{String, String}
    disks::Dict{String, Vector{String}}
    points::Dict{String, String}
    weights::Dict{String, Int}
    source_path::Union{Nothing, String}
end

function _lc_k_names(g::Integer)
    names = String["l$i" for i in 1:g]
    append!(names, String["m$i" for i in 1:g])
    g > 1 && append!(names, String["n$i" for i in 1:(g - 1)])
    names
end

function _lc_parse_text(text::AbstractString)
    genus = nothing
    words = Dict{String, Vector{Int}}()
    places = Dict{String, String}()
    disks = Dict{String, Vector{String}}()
    points = Dict{String, String}()
    weights = Dict{String, Int}()
    for raw in split(text, '\n')
        line = strip(split(raw, '#', limit=2)[1])
        isempty(line) && continue
        if startswith(line, "genus:")
            genus = parse(Int, strip(split(line, ':', limit=2)[2]))
            continue
        end
        parts_colon = split(line, ':', limit=2)
        length(parts_colon) != 2 && throw(KodaInputError(:read_lc, "not a component line", Dict("line" => raw)))
        left, rest = parts_colon[1], parts_colon[2]
        head = split(strip(left))[1]
        if Symbol(head) in LC_META
            kind = Symbol(head)
            parts = split(strip(left))
            length(parts) == 2 || throw(KodaInputError(:read_lc, "bad $kind line", Dict("line" => raw)))
            name = parts[2]
            if kind == :place
                ref = strip(rest)
                isempty(ref) && throw(KodaInputError(:read_lc, "place without region", Dict("line" => raw)))
                haskey(places, name) && throw(KodaInputError(:read_lc, "duplicate place", Dict("name" => name)))
                places[name] = ref
            elseif kind == :disk
                refs = split(strip(rest))
                isempty(refs) && throw(KodaInputError(:read_lc, "disk without regions", Dict("line" => raw)))
                haskey(disks, name) && throw(KodaInputError(:read_lc, "duplicate disk", Dict("name" => name)))
                disks[name] = collect(refs)
            elseif kind == :point
                target = strip(rest)
                isempty(target) && throw(KodaInputError(:read_lc, "point without target", Dict("line" => raw)))
                haskey(points, name) && throw(KodaInputError(:read_lc, "duplicate point", Dict("name" => name)))
                points[name] = target
            elseif kind == :weight
                val = parse(Int, strip(rest))
                haskey(weights, name) && throw(KodaInputError(:read_lc, "duplicate weight", Dict("name" => name)))
                weights[name] = val
            end
            continue
        end
        name = strip(left)
        name in string.(LC_META) && throw(KodaInputError(:read_lc, "reserved component name", Dict("name" => name)))
        haskey(words, name) && throw(KodaInputError(:read_lc, "duplicate component", Dict("name" => name)))
        words[name] = parse.(Int, split(strip(rest)))
    end
    genus === nothing && throw(KodaInputError(:read_lc, "missing genus line", Dict()))
    LickorishCode(genus, words, places, disks, points, weights, nothing)
end

"""
    read_lc(path::AbstractString) -> LickorishCode

Read a Lickorish-code record from a text file.

# Arguments
`path` names a readable file in the [LC grammar](../specs/lickorish-code.md#3-lc-のデータ形式gauss-code-形式).
File errors propagate. Checked syntax errors raise `KodaInputError` naming
missing or duplicate declarations; invalid integer tokens raise parse errors.
Geometric validity is the caller's responsibility until `validate_lc_input`
and the relevant realization checks are run.

# Returns
A `LickorishCode` with the parsed genus, words and annotation dictionaries,
and `source_path` equal to the absolute input path.

# Guarantees
1. The record represents the declarations read from the file, in the linked
   LC format. [proved/agent: by construction]
2. Reading does not modify the file. [proved/agent: by construction]

# Not guaranteed
Successful parsing alone does not certify an embedded curve arrangement or
the surface's genus.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> mktemp() do path, io; write(io, "genus: 1\\nl1 : -1\\nm1 : 1\\n"); close(io); lc = read_lc(path); (lc.genus, lc.words["m1"]) end
(1, [1])
```

See also [`validate_lc_input`](@ref), [`realize_lc`](@ref), [`write_lc`](@ref).

Role: construct
"""
function read_lc(path::AbstractString)
    text = read(path, String)
    lc = _lc_parse_text(text)
    LickorishCode(lc.genus, lc.words, lc.places, lc.disks, lc.points, lc.weights, abspath(path))
end

"""
    validate_lc_input(lc::LickorishCode) -> Bool

Check the preliminary genus, reference-name and signed-incidence requirements of an LC record.

# Arguments
`lc` is an LC record. The implemented preliminary checks concern
[LC input validity](../specs/lickorish-code.md).
Failures raise `KodaInputError` with stage `:validate_lc` and detail naming
the negative genus, missing K component or conflicting label incidence.

# Returns
`true` if these preliminary checks succeed; invalid input raises an exception.

# Guarantees
1. Successful return establishes the checked genus, K-name and signed-label
   consistency of the input. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This is not the full LC validator: it does not check the complete K words,
the realized surface, placements, disks or marks. Use the applicable
realization and diagram checks before a geometric interpretation.

# Sources
cf: lib/topology.jl/docs/specs/lickorish-code.md §4 — domain of LC consistency checks; this entry covers the preliminary subset.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); validate_lc_input(lc)
true
```

See also [`realize_lc`](@ref), [`validate_diagram`](@ref).

Role: condition
"""
function validate_lc_input(lc::LickorishCode)
    g = lc.genus
    g < 0 && throw(KodaInputError(:validate_lc, "genus must be nonnegative", Dict("genus" => g)))
    if g > 0
        knames = _lc_k_names(g)
        for k in knames
            haskey(lc.words, k) || throw(KodaInputError(:validate_lc, "missing K component", Dict("component" => k)))
        end
    end
    labels = Int[]
    for w in values(lc.words)
        append!(labels, abs.(w))
    end
    for k in labels
        count(==(k), labels) == 2 || throw(KodaInputError(:validate_lc, "label must appear twice", Dict("label" => k)))
    end
    owners = Dict{Tuple{Int, Symbol}, String}()
    for (name, word) in lc.words
        for x in word
            key = (abs(x), x > 0 ? :+ : :-)
            haskey(owners, key) && throw(KodaInputError(:validate_lc, "label owner collision", Dict("label" => abs(x))))
            owners[key] = name
        end
    end
    true
end

"""
    HeegaardCode(lc, roles, disks, points, weights, hc_path) -> HeegaardCode

Store a Lickorish-based Heegaard diagram's curve roles and annotations.

# Arguments
`lc::LickorishCode`, `roles::Dict{String,Symbol}`,
`disks::Dict{String,Vector{String}}`, `points::Dict{String,String}`,
`weights::Dict{String,Int}`, and `hc_path::Union{Nothing,String}`, or
convertible values. Conversion errors propagate. This constructor does not
check [HC input conditions](../specs/heegaard-code.md#4-hc-の入力検査) or
[the public diagram domain](../specs/heegaard-diagram.md).

# Returns
`HeegaardCode`, also exposed as `Topology.HeegaardDiagram`, holding the LC,
curve-role map, named disks, points, region weights and optional source path.
The meanings of the annotations belong to
[the HC specification](../specs/heegaard-code.md#3-hc-に置く追加情報).

# Guarantees
1. Fields preserve supplied values after conversion; already-typed objects
   and dictionaries are retained and not modified.
   [proved/agent: by construction]
2. The value carries HC data over the supplied LC. It makes no automatic
   standard-role binding; use [`resolve_hc`](@ref) for that operation.
   [proved/agent: by construction]

# Not guaranteed
Construction does not realize a surface or certify a Heegaard splitting.
The name `HeegaardDiagram` is not evidence that arbitrary fields satisfy
the diagram conditions. Later mutation of shared containers can invalidate
previously checked data. A source path is recorded without reading it.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §§1–4 and lib/topology.jl/docs/specs/heegaard-diagram.md — public diagram data, role binding and validation scope.

# Examples
```jldoctest
julia> H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1" => [-1], "m1" => [1,2], "b" => [-2]), Dict{String,String}(), Dict{String,Vector{String}}(), Dict{String,String}(), Dict{String,Int}(), nothing);

julia> a = H.resolve_hc(lc, Dict("b" => :beta); standard=true); b = H.HeegaardCode(a.lc, a.roles, a.disks, a.points, a.weights, a.hc_path);

julia> (b.roles == Dict("m1" => :alpha, "b" => :beta), H.is_standard(b), b isa Topology.HeegaardDiagram)
(true, true, true)
```

See also [`resolve_hc`](@ref), [`validate_hc_input`](@ref), [`realize_hc`](@ref).

Role: convert
"""
struct HeegaardCode
    lc::LickorishCode
    roles::Dict{String, Symbol}
    disks::Dict{String, Vector{String}}
    points::Dict{String, String}
    weights::Dict{String, Int}
    hc_path::Union{Nothing, String}
end

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/overlay/lib/topology.jl/docs/specs/heegaard-code.md §2.3–§2.4 — `standard` は各 m_i への α 束縛であり、順序付き α 支持が [m1,…,mg] のときだけ standard
function hc_component_names(lc::LickorishCode)
    reference = lc.genus > 0 ? _lc_k_names(lc.genus) : String[]
    vcat(reference, sort([name for name in keys(lc.words) if !(name in reference)]))
end

"""
    resolve_hc(lc::LickorishCode, roles::AbstractDict; standard=false, disks=Dict{String,Vector{String}}(), points=Dict{String,String}(), weights=Dict{String,Int}(), hc_path=nothing) -> HeegaardCode

Bind role and mark data to an existing LC.

# Arguments
`roles` uses the [HC role and standard conventions](../specs/heegaard-code.md). Standard syntax conflicts, missing reference meridians, missing or invalid role assignments and nonexistent role targets raise `KodaInputError`. `validate_hc_input` errors and dictionary conversions propagate. The caller separately establishes geometric validity of disk, point and weight supports.

# Returns
A `HeegaardCode` retaining the supplied LC and `hc_path`, with newly constructed roles and copied mark dictionaries; disk vectors are deep-copied.

# Guarantees
1. Non-K components receive a declared role and standard mode supplies the prescribed reference-meridian binding. [proved/agent: by construction]
2. Inputs are not modified; the LC is retained by reference, so later LC mutation is visible through the result. [proved/agent: by construction]

# Not guaranteed
This does not create duplicate meridian curves, certify full cut systems or validate all mark geometry. The returned HC does not freeze its shared LC.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §§2–4 — role assignment, standard binding and additional data.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = H.resolve_hc(lc, Dict("l1"=>:beta, "m1"=>:alpha));

julia> (h.lc === lc, H.is_standard(h))
(true, true)
```

Role: convert
"""
function resolve_hc(lc::LickorishCode, roles::AbstractDict;
                    standard::Bool=false, disks=Dict{String,Vector{String}}(),
                    points=Dict{String,String}(), weights=Dict{String,Int}(), hc_path=nothing)
    r = Dict{String,Symbol}(String(k) => Symbol(v) for (k, v) in roles)
    if standard
        any(==(:alpha), values(r)) && throw(KodaInputError(:resolve_hc,
            "standard and explicit alpha roles are mutually exclusive", Dict{String,Any}()))
        for i in 1:lc.genus
            name = "m$i"
            haskey(lc.words, name) || throw(KodaInputError(:resolve_hc,
                "missing reference meridian", Dict{String,Any}("curve" => name)))
            haskey(r, name) && throw(KodaInputError(:resolve_hc,
                "meridian already has a conflicting role", Dict{String,Any}("curve" => name)))
            r[name] = :alpha
        end
    end
    knames = Set(_lc_k_names(lc.genus))
    for name in keys(lc.words)
        name in knames && continue
        haskey(r, name) || throw(KodaInputError(:resolve_hc,
            "non-K curve has no role", Dict{String,Any}("curve" => String(name))))
    end
    for (name, role) in r
        haskey(lc.words, name) || throw(KodaInputError(:resolve_hc,
            "role target does not exist", Dict{String,Any}("curve" => name)))
        role in (:alpha, :beta) || throw(KodaInputError(:resolve_hc,
            "invalid role", Dict{String,Any}("curve" => name, "role" => String(role))))
    end
    hc = HeegaardCode(lc, r,
        deepcopy(Dict{String,Vector{String}}(disks)),
        copy(Dict{String,String}(points)),
        copy(Dict{String,Int}(weights)),
        hc_path)
    validate_hc_input(hc)
    hc
end

hc_curve_names(hc::HeegaardCode, family::Symbol) =
    String[n for n in hc_component_names(hc.lc) if get(hc.roles, n, :none) == family]

"""
    is_standard(hc::HeegaardCode) -> Bool

Test whether an HC has the standard alpha binding to its fixed reference curves.

# Arguments
`hc` must satisfy [HC input requirements](../specs/heegaard-code.md#4-hc-の入力検査).
`validate_hc_input` checks them and its `KodaInputError` diagnostics propagate;
an invalid record is not converted to a negative answer.

# Returns
`Bool` for the standard binding defined in
[the HC format](../specs/heegaard-code.md).

# Guarantees
1. The result tests the existing standard binding of the HC independently
   of whether it was read using compact or explicit text encoding.
   [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
The predicate does not change roles, introduce reference curves, certify
complete cutting systems or identify the represented 3-manifold.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1" => [-1], "m1" => [1,2], "b" => [-2]), Dict(), Dict(), Dict(), Dict(), nothing);

julia> h = HeegaardCode(lc, Dict("m1" => :alpha, "b" => :beta), Dict(), Dict(), Dict(), nothing);

julia> is_standard(h)
true
```

See also [`resolve_hc`](@ref), [`hc_text`](@ref).

Role: condition
"""
function is_standard(hc::HeegaardCode)
    validate_hc_input(hc)
    hc_curve_names(hc, :alpha) == String["m$i" for i in 1:hc.lc.genus]
end

function _hc_parse_text(text::AbstractString, lc_path_hint=nothing)
    roles = Dict{String,Symbol}()
    disks = Dict{String,Vector{String}}()
    points = Dict{String,String}()
    weights = Dict{String,Int}()
    lc_ref = lc_path_hint
    saw_lc = lc_path_hint !== nothing
    format_version = nothing
    saw_standard = false
    for (lineno, raw) in enumerate(split(replace(text, "\ufeff" => ""), '\n'))
        line = strip(split(raw, '#', limit=2)[1])
        isempty(line) && continue
        detail = Dict{String,Any}("line" => lineno, "text" => String(raw))
        if line == "standard"
            saw_standard && throw(KodaInputError(:read_hc, "duplicate standard directive", detail))
            saw_standard = true
            continue
        end
        legacy = match(r"^role\s+(\S+)\s+(alpha|beta)$"i, line)
        if legacy !== nothing && !occursin(':', line)
            name, family = String.(legacy.captures)
            haskey(roles, name) && throw(KodaInputError(:read_hc, "duplicate role", detail))
            roles[name] = Symbol(lowercase(family))
            continue
        end
        colon = match(r"^([^:]+)\s*:\s*(.*)$", line)
        colon === nothing && throw(KodaInputError(:read_hc, "unknown HC directive", detail))
        lhs, rhs = strip.(String.(colon.captures))
        isempty(rhs) && throw(KodaInputError(:read_hc, "empty directive value", detail))
        if lhs == "format"
            format_version === nothing || throw(KodaInputError(:read_hc, "duplicate format directive", detail))
            rhs == "hc/2" || throw(KodaInputError(:read_hc, "unsupported HC format version", detail))
            format_version = 2
            continue
        end
        if lhs == "lc"
            saw_lc && throw(KodaInputError(:read_hc, "duplicate lc directive", detail))
            saw_lc = true
            lc_ref = rhs
            continue
        end
        head = split(lhs)
        length(head) == 2 || throw(KodaInputError(:read_hc, "unknown HC directive", detail))
        kind, name = String.(head)
        if kind == "role"
            lowercase(rhs) in ("alpha", "beta") || throw(KodaInputError(:read_hc, "role must be alpha or beta", detail))
            haskey(roles, name) && throw(KodaInputError(:read_hc, "duplicate role", detail))
            roles[name] = Symbol(lowercase(rhs))
        elseif kind == "disk"
            haskey(disks, name) && throw(KodaInputError(:read_hc, "duplicate disk", detail))
            disks[name] = String.(split(rhs))
        elseif kind == "point"
            haskey(points, name) && throw(KodaInputError(:read_hc, "duplicate point", detail))
            points[name] = rhs
        elseif kind == "weight"
            haskey(weights, name) && throw(KodaInputError(:read_hc, "duplicate weight", detail))
            value = tryparse(Int, rhs)
            value === nothing && throw(KodaInputError(:read_hc, "weight must be an integer", detail))
            weights[name] = value
        else
            throw(KodaInputError(:read_hc, "unknown HC directive", detail))
        end
    end
    lc_ref === nothing && throw(KodaInputError(:read_hc, "missing lc reference", Dict{String,Any}()))
    (; lc_ref=String(lc_ref), roles, disks, points, weights, standard=saw_standard, format_version)
end

"""
    read_hc(path::AbstractString) -> HeegaardCode

Read an HC file and resolve its referenced LC with the declared roles and marks.

# Arguments
`path` must name a readable file in [the HC format](../specs/heegaard-code.md). The LC reference must be relative; standard or reference-role syntax requires format `hc/2`. Parser and role violations raise `KodaInputError` with the relevant stage and directive. LC reading, I/O and `resolve_hc` errors propagate.

# Returns
A resolved `HeegaardCode` whose `hc_path` is absolute and whose LC is read relative to the HC file.

# Guarantees
1. Roles and supported disk, point and weight declarations are retained in the resolved HC; the standard syntax has `resolve_hc`'s meaning. [proved/agent: by construction]
2. Input files are not modified. [proved/agent: by construction]

# Not guaranteed
Reading does not certify complete cut systems, joining disks, a successful drawing or a successful Koda conversion.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §§2–5 — HC file format and validation boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> mktempdir() do dir
           write(joinpath(dir,"d.lc"), "genus: 1\\nl1 : -1\\nm1 : 1\\n")
           write(joinpath(dir,"d.hc"), "format : hc/2\\nlc : d.lc\\nrole m1 : alpha\\nrole l1 : beta\\n")
           h = H.read_hc(joinpath(dir,"d.hc"))
           (h.lc.genus, H.is_standard(h))
       end
(1, true)
```

Role: construct
"""
function read_hc(path::AbstractString)
    absolute = abspath(path)
    parsed = _hc_parse_text(read(absolute, String))
    isabspath(parsed.lc_ref) && throw(KodaInputError(:read_hc,
        "lc reference must be relative", Dict{String,Any}("path" => parsed.lc_ref)))
    lc = read_lc(normpath(joinpath(dirname(absolute), parsed.lc_ref)))
    if parsed.standard || any(name -> haskey(parsed.roles, name), _lc_k_names(lc.genus))
        parsed.format_version == 2 || throw(KodaInputError(:read_hc,
            "standard/reference-role input requires 'format : hc/2' to protect against legacy readers",
            Dict{String,Any}()))
    end
    resolve_hc(lc, parsed.roles; standard=parsed.standard, disks=parsed.disks,
        points=parsed.points, weights=parsed.weights, hc_path=absolute)
end

"""
    validate_hc_input(hc::HeegaardCode) -> Bool

Check preliminary LC consistency and selected role requirements of an HC record.

# Arguments
`hc` is an HC record in the terminology of [HC input checks](../specs/heegaard-code.md#4-hc-の入力検査).
LC input errors propagate. Invalid non-K role values or the absence of every
alpha and beta role raise `KodaInputError` with stage `:validate_hc`.

# Returns
`true` if the implemented preliminary checks succeed.

# Guarantees
1. Successful return means `validate_lc_input` succeeds and the implemented
   role checks accept the record. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
This does not enforce every HC requirement: it does not check that every role
target exists, that every non-K curve has a role, or that role values on K
curves are valid. It does not certify cut systems, disk data or a 3-manifold.
Use `resolve_hc` and the applicable diagram checks for those input workflows.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-code.md §4 — HC input domain; this preliminary routine checks only a subset.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> lc = LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = HeegaardCode(lc, Dict("m1"=>:alpha, "l1"=>:beta), Dict(), Dict(), Dict(), nothing); validate_hc_input(h)
true
```

See also [`resolve_hc`](@ref), [`validate_diagram`](@ref).

Role: condition
"""
function validate_hc_input(hc::HeegaardCode)
    validate_lc_input(hc.lc)
    g = hc.lc.genus
    for (name, _) in hc.lc.words
        haskey(hc.roles, name) || continue
    end
    knames = Set(_lc_k_names(g))
    for (name, role) in hc.roles
        name in knames && continue
        role in (:alpha, :beta) || throw(KodaInputError(:validate_hc, "invalid role", Dict("curve" => name)))
    end
    alpha = [k for (k, r) in hc.roles if r == :alpha]
    beta = [k for (k, r) in hc.roles if r == :beta]
    isempty(alpha) && isempty(beta) &&
        throw(KodaInputError(:validate_hc, "no alpha or beta roles assigned", Dict()))
    true
end
