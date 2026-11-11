# cf: lib/topology.jl/docs/specs/heegaard-moves.md §4–§6 — embedded support witnesses
# cf: .research-scratch/PJ260922/261111-gpt-pro-review.md §4 — left-face convention; reconstructed

function _move_left_boundary(words, occ, face)
    [(component=occ[c[1]][c[2]][1], position=occ[c[1]][c[2]][2] + 1,
      forward=c[3] == "out", corner=c)
     for c in (_lc_tau(words, occ, z) for z in reverse(face.corners))]
end

function _move_corner(value)
    length(value) == 3 || throw(LCMoveError(:support, "a half-edge is (label, passage sign, in/out)", Dict()))
    c = (Int(value[1]), String(value[2]), String(value[3]))
    c[2] in ("+", "-") && c[3] in ("in", "out") ||
        throw(LCMoveError(:support, "invalid half-edge direction", Dict("corner" => value)))
    c
end

"""Validate a finite embedded path across distinct disk faces.

`start`, `finish`, and each entry of `crossings` are LEFT-face half-edge
occurrences (label, passage sign, in/out). The path enters the face left of
start, crosses each specified edge to its opposite side, and ends on finish.
Distinct visited faces make the interior simple; repeated-face paths need
additional ordering data and are explicitly unsupported. No path is chosen.
"""
function lc_move_path(lc::LickorishConfiguration; start, crossings=Any[], finish, avoid=String[])
    validate_lc_input(lc)
    active = _lc_nonempty_words(lc.words)
    faces, right = _lc_faces(active)
    occ = _lc_rotation_occ(active)
    left = Dict(_lc_tau(active, occ, c) => id for (c,id) in right)
    a, b = _move_corner(start), _move_corner(finish)
    all(haskey(left,c) for c in (a,b)) || throw(LCMoveError(:support, "endpoint is not a boundary occurrence", Dict()))
    current = left[a]
    visited = String[current]
    steps = Dict{String,Any}[]
    for raw in crossings
        c = _move_corner(raw)
        get(left,c,nothing) == current || throw(LCMoveError(:disconnected_path,
            "crossing is not on the current face", Dict("face" => current, "corner" => raw)))
        name, pos0 = occ[c[1]][c[2]]
        name in avoid && throw(LCMoveError(:path_meets_family,
            "path crosses an excluded curve", Dict("component" => name)))
        opposite = _lc_tau(active,occ,c)
        next = left[opposite]
        next in visited && throw(LCMoveError(:path_order_required,
            "a repeated face needs explicit disjoint arc ordering", Dict("face" => next)))
        gap = c[3] == "out" ? mod1(pos0 + 2,length(active[name])) : pos0 + 1
        push!(steps,Dict("from"=>current,"to"=>next,"corner"=>collect(c),
            "component"=>name,"gap"=>gap,"forward"=>c[3]=="out"))
        push!(visited,next)
        current = next
    end
    left[b] == current || throw(LCMoveError(:disconnected_path,
        "finish is not on the final face", Dict("face"=>current,"finish"=>collect(b))))
    a != b || throw(LCMoveError(:support,"path endpoints must be distinct occurrences",Dict()))
    any(id in values(lc.places) for id in visited) && throw(LCMoveError(:path_order_required,
        "place-only circles require a path in the punctured face",Dict("faces"=>visited)))
    Dict{String,Any}("kind"=>"disk_face_path","faces"=>visited,
        "start"=>collect(a),"finish"=>collect(b),"crossings"=>steps,
        "range"=>"paths visiting each disk face at most once")
end

function _move_require_separated(input, components)
    input isa HeegaardLayer || return
    shared = [n for n in components if _is_k(input.lc,n) && haskey(input.roles,n)]
    isempty(shared) || throw(LCMoveError(:shared_support_requires_separation,
        "moving a role without moving K requires separating its shared support first",
        Dict("components"=>shared)))
end

export lc_move_path
