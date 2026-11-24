# cf: threads/PJ260922/ACTIVE.md — Owner decisions 2026-11-12; free mark placement.
# cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-4-1-full.md — automatic reference transport; reconstructed.

function _move_untouched_cells(lc,witness)
    faces,right=_lc_faces(_lc_nonempty_words(lc.words));occ=_lc_rotation_occ(lc.words)
    touched=Set{String}();edges=Set{Tuple{String,Int}}()
    function edge!(name,i)
        word=get(lc.words,name,Int[]);isempty(word) && return
        i=mod1(i,length(word));push!(edges,(name,i));x=word[i]
        d=(abs(x),x>0 ? "+" : "-","out")
        push!(touched,right[d]);push!(touched,right[_lc_tau(lc.words,occ,d)])
    end
    for (j,(ci,pos)) in enumerate(witness.locations)
        name=witness.names[ci];word=get(lc.words,name,Int[]);isempty(word) && continue
        # The splice can alter both endpoint germs, so exclude the neighboring edges too.
        for i in pos-1:pos+length(witness.before_frags[j])-1;edge!(name,i);end
    end
    for run in get(witness.extra,"runs",Any[])
        run.kind==:curve || continue
        for i in eachindex(lc.words[run.component]);edge!(run.component,i);end
    end
    Set(f.id for f in faces if !(f.id in touched)),edges
end

function _move_ref_key(s::AbstractString)
    m = match(r"^F([0-9]+)$", s)
    m === nothing ? (1, 0, String(s)) :
                    (0, parse(Int, m.captures[1]), String(s))
end

function _move_cyclic_equal(a, b)
    n = length(a)
    n == length(b) || return false
    n == 0 && return true
    any(k -> all(a[i] == b[mod1(i + k, n)] for i in 1:n), 0:n-1)
end

function _move_exact_face_map(oldlc, newlc, label_map)
    oldfaces, _ = _lc_faces(_lc_nonempty_words(oldlc.words))
    newfaces, _ = _lc_faces(_lc_nonempty_words(newlc.words))
    fmap = Dict{String,String}()

    for f in oldfaces
        all(haskey(label_map, c[1]) for c in f.corners) || continue
        cs = [(label_map[c[1]], c[2], c[3]) for c in f.corners]
        hits = sort!(
            String[g.id for g in newfaces
                   if _move_cyclic_equal(cs, g.corners)];
            by=_move_ref_key,
        )
        length(hits) == 1 && (fmap[f.id] = only(hits))
    end
    fmap
end

_move_edge_ref(name, x::Int) =
    string(name, "/", x > 0 ? "+" : "", x)

function _move_edge_image(ref::String, oldlc, newlc, label_map)
    p = split(ref, '/'; limit=2)
    length(p) == 2 || return nothing
    name, token = p
    x = tryparse(Int, token)
    x === nothing && return nothing

    haskey(oldlc.words, name) && haskey(newlc.words, name) ||
        return nothing
    oldword, newword = oldlc.words[name], newlc.words[name]
    (isempty(oldword) || isempty(newword)) && return nothing

    oldhits = findall(==(x), oldword)
    length(oldhits) == 1 || return nothing
    i = only(oldhits)
    y = oldword[mod1(i + 1, length(oldword))]

    haskey(label_map, abs(x)) && haskey(label_map, abs(y)) ||
        return nothing
    a = sign(x) * label_map[abs(x)]
    b = sign(y) * label_map[abs(y)]

    newhits = [
        j for j in eachindex(newword)
        if newword[j] == a &&
           newword[mod1(j + 1, length(newword))] == b
    ]
    length(newhits) == 1 ? _move_edge_ref(name, a) : nothing
end

function _move_mark_cells(lc)
    fs, _ = _lc_faces(_lc_nonempty_words(lc.words))
    faces = sort!(String[f.id for f in fs]; by=_move_ref_key)
    edges, empties, sides = String[], String[], String[]

    for name in sort!(collect(keys(lc.words)))
        word = lc.words[name]
        if isempty(word) && haskey(lc.places, name)
            push!(empties, name)
            append!(sides, (name * "+", name * "-"))
        else
            append!(edges, (_move_edge_ref(name, x) for x in word))
        end
    end

    (
        faces=faces,
        edges=sort!(unique(edges); by=_move_ref_key),
        empties=empties,
        weight_sides=sides,
    )
end

function _move_point_image(ref, oldlc, newlc, label_map, fmap, cells)
    haskey(fmap, ref) && return fmap[ref], :exact_face

    edge = _move_edge_image(ref, oldlc, newlc, label_map)
    edge !== nothing && return edge, :exact_edge

    if haskey(oldlc.words, ref) &&
       isempty(oldlc.words[ref]) &&
       ref in cells.empties
        return ref, :exact_empty
    end

    !isempty(cells.faces) &&
        return first(cells.faces), :reassigned_face
    !isempty(cells.edges) &&
        return first(cells.edges), :reassigned_edge
    !isempty(cells.empties) &&
        return first(cells.empties), :reassigned_empty

    nothing, :deferred
end

function _move_weight_image(ref, oldlc, fmap, cells)
    haskey(fmap, ref) && return fmap[ref], :exact_face

    if length(ref) > 1 &&
       (endswith(ref, "+") || endswith(ref, "-"))
        name = chop(ref; tail=1)
        if haskey(oldlc.words, name) &&
           isempty(oldlc.words[name]) &&
           ref in cells.weight_sides
            return ref, :exact_empty_side
        end
    end

    !isempty(cells.faces) &&
        return first(cells.faces), :reassigned_face
    !isempty(cells.weight_sides) &&
        return first(cells.weight_sides), :reassigned_empty_side

    nothing, :deferred
end

_move_add_weight!(d::Dict{String,Int}, ref::String, value::Int) =
    (d[ref] = get(d, ref, 0) + value)

_move_rebuild_marks(x::HeegaardCode, disks, points, weights) =
    HeegaardCode(
        x.lc, x.roles, disks, points, weights, x.hc_path,
    )

_move_rebuild_marks(x::LickorishCode, disks, points, weights) =
    LickorishCode(
        x.genus, x.words, x.places,
        disks, points, weights, x.source_path,
    )

function _move_transport_marks(input, output, witness, label_map)
    _ = witness  # geometric move has already been completed
    isempty(_missing_marks(input)) && return output, (faces=NamedTuple[],refs=NamedTuple[],weight_collisions=NamedTuple[])

    oldlc = _lc_of(input)
    newlc = _lc_of(output)
    marks = _marks(input)
    fmap = _move_exact_face_map(oldlc, newlc, label_map)
    cells = _move_mark_cells(newlc)
    refs = NamedTuple[]

    disks = Dict{String,Vector{String}}()
    for name in sort!(collect(keys(marks.disks)))
        oldrefs = marks.disks[name]
        mapped = String[get(fmap, ref, "") for ref in oldrefs]
        exact =
            !isempty(oldrefs) &&
            all(ref -> !isempty(ref), mapped) &&
            length(unique(mapped)) == length(mapped)

        if exact
            disks[name] = mapped
        else
            if isempty(cells.faces)
                for oldref in oldrefs
                    push!(refs,(kind=:disk,mark=name,old=oldref,new=nothing,mode=:deferred,value=nothing))
                end
                continue
            end
            mapped = fill(first(cells.faces), length(oldrefs))
            disks[name] = [first(cells.faces)]
        end

        mode = exact ? :exact_face : :reassigned_face
        for (oldref, newref) in zip(oldrefs, mapped)
            push!(refs, (
                kind=:disk, mark=name, old=oldref, new=newref,
                mode=mode, value=nothing,
            ))
        end
    end

    points = Dict{String,String}()
    for name in sort!(collect(keys(marks.points)))
        oldref = marks.points[name]
        newref, mode = _move_point_image(
            oldref, oldlc, newlc, label_map, fmap, cells,
        )
        newref===nothing || (points[name] = newref)
        push!(refs, (
            kind=:point, mark=name, old=oldref, new=newref,
            mode=mode, value=nothing,
        ))
    end

    weights = Dict{String,Int}()
    sources = Dict{String,Vector{String}}()
    for oldref in sort!(
        collect(keys(marks.weights)); by=_move_ref_key
    )
        newref, mode =
            _move_weight_image(oldref, oldlc, fmap, cells)
        value = marks.weights[oldref]
        if newref===nothing
            push!(refs,(kind=:weight,mark="",old=oldref,new=nothing,mode=:deferred,value=value))
            continue
        end

        # The assignment inserts newref even when value or the sum is zero.
        _move_add_weight!(weights, newref, value)
        push!(get!(sources, newref, String[]), oldref)
        push!(refs, (
            kind=:weight, mark="", old=oldref, new=newref,
            mode=mode, value=value,
        ))
    end

    collisions = [
        (new=newref, old=copy(sources[newref]), value=weights[newref])
        for newref in sort!(
            collect(keys(sources)); by=_move_ref_key
        )
        if length(sources[newref]) > 1
    ]

    receipt = (
        faces=[
            (old=oldref, new=fmap[oldref])
            for oldref in sort!(
                collect(keys(fmap)); by=_move_ref_key
            )
        ],
        refs=refs,
        weight_collisions=collisions,
    )

    _move_rebuild_marks(output, disks, points, weights), receipt
end

# The geometric inverse recomputes an annulus from the current chart. That
# chart no longer knows which face held a mark before the forward move, so a
# failed exact face map leaves the mark on the first face. Remember the input
# marks and put them back when the inverse restores the same words and places.
function _recorded_marks(input)
    marks = _marks(input)
    Dict{String,Any}(
        "disks" => deepcopy(marks.disks),
        "points" => deepcopy(marks.points),
        "weights" => deepcopy(marks.weights),
        "lc" => lc_text(_lc_of(input)),
    )
end

function _restore_recorded_marks(result::LCMoveResult, witness::LCMoveWitness)
    snap = get(witness.extra, "mark_before", nothing)
    snap === nothing && return result
    lc_text(_lc_of(result.output)) == snap["lc"] || return result
    output = _move_rebuild_marks(
        result.output,
        deepcopy(snap["disks"]),
        deepcopy(snap["points"]),
        deepcopy(snap["weights"]),
    )
    LCMoveResult(output, result.inverse, result.correspondence)
end
