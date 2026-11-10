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

function _hc_parse_text(text::AbstractString, lc_path_hint=nothing)
    roles = Dict{String, Symbol}()
    disks = Dict{String, Vector{String}}()
    points = Dict{String, String}()
    weights = Dict{String, Int}()
    lc_ref = lc_path_hint
    for raw in split(text, '\n')
        line = strip(split(raw, '#', limit=2)[1])
        isempty(line) && continue
        if startswith(line, "lc")
            parts = split(strip(line))
            length(parts) >= 2 || throw(KodaInputError(:read_hc, "bad lc reference line", Dict()))
            lc_ref = strip(parts[2])
            continue
        end
        if startswith(line, "role")
            parts = split(strip(line))
            length(parts) == 3 || throw(KodaInputError(:read_hc, "bad role line", Dict("line" => raw)))
            fam = Symbol(lowercase(parts[3]))
            fam in (:alpha, :beta) || throw(KodaInputError(:read_hc, "role must be alpha or beta", Dict()))
            roles[parts[2]] = fam
            continue
        end
        parts_colon = split(line, ':', limit=2)
        length(parts_colon) != 2 && continue
        left, rest = parts_colon[1], parts_colon[2]
        head = split(strip(left))[1]
        name = split(strip(left))[2]
        if head == "disk"
            refs = split(strip(rest))
            disks[name] = refs
        elseif head == "point"
            points[name] = strip(rest)
        elseif head == "weight"
            weights[name] = parse(Int, strip(rest))
        end
    end
    lc_ref === nothing && throw(KodaInputError(:read_hc, "HC file must reference an lc file", Dict()))
    throw(KodaImplementationError(:_hc_parse_text, "parser output must be completed by read_hc", Dict()))
end

function read_hc(path::AbstractString)
    dir = dirname(abspath(path))
    text = read(path, String)
    lc_ref = nothing
    roles = Dict{String, Symbol}()
    disks = Dict{String, Vector{String}}()
    points = Dict{String, String}()
    weights = Dict{String, Int}()
    for raw in split(text, '\n')
        line = strip(split(raw, '#', limit=2)[1])
        isempty(line) && continue
        if occursin(r"^lc\s*:", line)
            lc_ref = joinpath(dir, strip(split(line, ':', limit=2)[2]))
            continue
        end
        if startswith(line, "role")
            parts = split(strip(line))
            length(parts) >= 3 || throw(KodaInputError(:read_hc, "bad role line", Dict("line" => raw)))
            if occursin(':', line)
                left, right = split(line, ':', limit=2)
                name = strip(split(strip(left))[2])
                fam = strip(right)
            else
                name = parts[2]
                fam = parts[3]
            end
            roles[name] = Symbol(lowercase(fam))
            continue
        end
        parts_colon = split(line, ':', limit=2)
        length(parts_colon) != 2 && continue
        left, rest = parts_colon[1], parts_colon[2]
        head, name = split(strip(left))[1:2]
        head == "disk" && (disks[name] = split(strip(rest)))
        head == "point" && (points[name] = strip(rest))
        head == "weight" && (weights[name] = parse(Int, strip(rest)))
    end
    lc_ref === nothing && throw(KodaInputError(:read_hc, "missing lc reference", Dict("hc" => path)))
    lc = read_lc(lc_ref)
    HeegaardLayer(lc, roles, disks, points, weights, abspath(path))
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
