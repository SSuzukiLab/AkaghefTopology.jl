# Source: docs/specs/lickorish-code.md §3–4;
# lib/topology.jl/scripts/heegaard/lc_input.py::parse.
"""Lickorish Code (LC): words and placement metadata without Heegaard roles."""

const LC_META = (:place, :disk, :point, :weight)

struct LickorishConfiguration
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
    LickorishConfiguration(genus, words, places, disks, points, weights, nothing)
end

function read_lc(path::AbstractString)
    text = read(path, String)
    lc = _lc_parse_text(text)
    LickorishConfiguration(lc.genus, lc.words, lc.places, lc.disks, lc.points, lc.weights, abspath(path))
end

function validate_lc_input(lc::LickorishConfiguration)
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

struct HeegaardLayer
    lc::LickorishConfiguration
    roles::Dict{String, Symbol}
    disks::Dict{String, Vector{String}}
    points::Dict{String, String}
    weights::Dict{String, Int}
    hc_path::Union{Nothing, String}
end

# cf: .research-scratch/PJ260922/261111-standard-revision/PJ260922_standard_bundle/overlay/lib/topology.jl/docs/specs/heegaard-code.md §2.3–§2.4 — `standard` は各 m_i への α 束縛であり、順序付き α 支持が [m1,…,mg] のときだけ standard
function hc_component_names(lc::LickorishConfiguration)
    reference = lc.genus > 0 ? _lc_k_names(lc.genus) : String[]
    vcat(reference, sort([name for name in keys(lc.words) if !(name in reference)]))
end

"""Bind `standard` onto the reference meridians. The LC words are not copied or rewritten."""
function resolve_hc(lc::LickorishConfiguration, roles::AbstractDict;
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
    hc = HeegaardLayer(lc, r,
        deepcopy(Dict{String,Vector{String}}(disks)),
        copy(Dict{String,String}(points)),
        copy(Dict{String,Int}(weights)),
        hc_path)
    validate_hc_input(hc)
    hc
end

hc_curve_names(hc::HeegaardLayer, family::Symbol) =
    String[n for n in hc_component_names(hc.lc) if get(hc.roles, n, :none) == family]

function is_standard(hc::HeegaardLayer)
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

"""Read an `.hc` file as one resolved HC.

`standard` assigns α to the existing meridians `m1,…,mg`. The value is a
`HeegaardLayer`, the same struct `Topology` exports as `HeegaardDiagram`.
No second curve and no stored standard flag are created.
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

function validate_hc_input(hc::HeegaardLayer)
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
