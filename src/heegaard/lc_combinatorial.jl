# Source: docs/specs/heegaard-koda-conversion.md §3.1; docs/specs/lickorish-code.md §4;
# lib/topology.jl/scripts/heegaard/lc_input.py (rotation_occ / face_names / surface).
"""Pure combinatorial LC rotation system (no marked mesh)."""

const LCCorner = Tuple{Int, String, String}

struct LCFace
    id::String
    corners::Vector{LCCorner}
end

struct LCRealization
    genus::Int
    lc::LickorishConfiguration
    faces::Vector{LCFace}
    corner_face::Dict{LCCorner, String}
    occurrence::Dict{Int, Dict{String, Tuple{String, Int}}}
    surface_chi::Int
end

_lc_sign_char(x::Integer) = x > 0 ? "+" : "-"
const _LC_ROT_SEQ = [("-", "out"), ("+", "out"), ("-", "in"), ("+", "in")]

function _lc_rot_next(sign::AbstractString, way::AbstractString)
    i = findfirst(t -> t == (sign, way), _LC_ROT_SEQ)
    _LC_ROT_SEQ[mod1(i + 1, 4)]
end

function _lc_rot_prev(sign::AbstractString, way::AbstractString)
    i = findfirst(t -> t == (sign, way), _LC_ROT_SEQ)
    _LC_ROT_SEQ[mod1(i - 1, 4)]
end

function _lc_nonempty_words(words::AbstractDict)
    Dict(k => copy(v) for (k, v) in words if !isempty(v))
end

function _lc_rotation_occ(words::AbstractDict)
    occ = Dict{Int, Dict{String, Tuple{String, Int}}}()
    for (name, word) in words
        for (pos, x) in enumerate(word)
            lab = abs(x)
            sg = _lc_sign_char(x)
            haskey(occ, lab) || (occ[lab] = Dict{String, Tuple{String, Int}}())
            occ[lab][sg] = (name, pos - 1)
        end
    end
    occ
end

function _lc_along(words::AbstractDict, occ, label::Int, sign::AbstractString, way::AbstractString)
    name, pos = occ[label][sign]
    word = words[name]
    n = length(word)
    idx = mod(pos + (way == "out" ? 1 : -1), n)
    x = word[idx + 1]
    abs(x), _lc_sign_char(x), way == "out" ? "in" : "out"
end

function _lc_tau(words::AbstractDict, occ, c::LCCorner)
    _lc_along(words, occ, c...)
end

function _lc_faces(words::AbstractDict)
    occ = _lc_rotation_occ(words)
    isempty(occ) && return LCFace[], Dict{LCCorner, String}()
    seen = Set{LCCorner}()
    faces = LCFace[]
    corner_face = Dict{LCCorner, String}()
    for label in sort(collect(keys(occ)))
        for (sign, way) in _LC_ROT_SEQ
            start = (label, sign, way)
            start in seen && continue
            corners = LCCorner[]
            z = start
            while !(z in seen)
                push!(seen, z)
                push!(corners, z)
                other, s2, w2 = _lc_tau(words, occ, z)
                s3, w3 = _lc_rot_next(s2, w2)
                z = (other, s3, w3)
            end
            fid = "F$(length(faces) + 1)"
            push!(faces, LCFace(fid, corners))
            for c in corners
                corner_face[c] = fid
            end
        end
    end
    faces, corner_face
end

function _lc_surface_stats(words::AbstractDict)
    occ = _lc_rotation_occ(words)
    isempty(occ) && return (0, 0, 0)
    v = length(occ)
    e = 2v
    faces, = _lc_faces(words)
    f = length(faces)
    v, e, f
end

function _lc_occurrence(words::AbstractDict)
    occ = Dict{Int, Dict{String, Tuple{String, Int}}}()
    for (name, word) in words
        for (pos, x) in enumerate(word)
            lab = abs(x)
            sg = _lc_sign_char(x)
            haskey(occ, lab) || (occ[lab] = Dict{String, Tuple{String, Int}}())
            occ[lab][sg] = (name, pos)
        end
    end
    occ
end

function _lc_band_edges(words::AbstractDict)
    edges = Tuple{String, Int, Int, Int, Int, String}[]
    for (name, word) in words
        n = length(word)
        n < 1 && continue
        n == 1 && continue
        for i in 1:n
            j = mod1(i + 1, n)
            push!(edges, (name, i, abs(word[i]), abs(word[j]), word[i], _lc_sign_char(word[i])))
        end
    end
    edges
end

function _lc_corner_at_occ(word::Vector{Int}, pos::Int, dir::AbstractString)
    x = word[pos]
    lab = abs(x)
    sg = _lc_sign_char(x)
    way = dir == "out" ? "out" : "in"
    (lab, sg, way)
end

function _lc_faces_incident_to_band(words::AbstractDict, corner_face::Dict{LCCorner, String}, comp::AbstractString, i::Int)
    word = words[comp]
    n = length(word)
    n >= 2 || return String[], String[]
    c_out = _lc_corner_at_occ(word, i, "out")
    c_in = _lc_corner_at_occ(word, mod1(i + 1, n), "in")
    f1 = get(corner_face, c_out, nothing)
    f2 = get(corner_face, c_in, nothing)
    f1 === nothing && f2 === nothing && return String[], String[]
    unique(filter(!isnothing, [f1, f2]))
end

function _lc_validate_surface!(lc::LickorishConfiguration)
    g = lc.genus
    active = _lc_nonempty_words(lc.words)
    if isempty(active)
        g == 0 || throw(KodaInputError(:realize_lc, "no curve data", Dict("genus" => g)))
        return
    end
    v, e, f = _lc_surface_stats(active)
    chi = v - e + f
    expected = 2 - 2g
    chi == expected ||
        throw(KodaInputError(
            :realize_lc,
            "Euler characteristic does not match genus",
            Dict("chi" => chi, "expected" => expected, "genus" => g, "vertices" => v, "edges" => e, "faces" => f),
        ))
end

function _lc_realize_combinatorial(lc::LickorishConfiguration)
    _lc_validate_surface!(lc)
    active = isempty(_lc_nonempty_words(lc.words)) ? lc.words : _lc_nonempty_words(lc.words)
    faces, corner_face = _lc_faces(active)
    v, e, f = _lc_surface_stats(active)
    LCRealization(lc.genus, lc, faces, corner_face, _lc_rotation_occ(active), v - e + f)
end
