function _hc_rank!(signatures)
    order = sort!(unique(signatures))
    index = Dict(signature => i for (i, signature) in enumerate(order))
    [index[signature] for signature in signatures]
end

function _hc_colors(words, roles, genus)
    kn = _hc_k_set(genus)
    occs = NamedTuple{(:comp, :pos, :sign, :label, :init), Tuple{String,Int,Int,Int,Tuple}}[]
    for name in sort!(collect(String.(keys(words))))
        word = words[name]
        if name in kn
            init_role = (:k, name)
        elseif get(roles, name, nothing) == :alpha || get(roles, name, nothing) == :beta
            init_role = (roles[name],)
        else
            init_role = (:unassigned, name)
        end
        for (pos, x) in enumerate(word)
            push!(occs, (comp=name, pos=pos, sign=sign(x), label=abs(x), init=(:occ, init_role, sign(x))))
        end
    end
    isempty(occs) && return Dict{String,Vector{Int}}(), Int[]
    labels = sort!(unique(o.label for o in occs))
    lab_index = Dict(lab => i for (i, lab) in enumerate(labels))
    n_occ = length(occs)
    n = n_occ + length(labels)
    neighbors = [Tuple{Int,Int}[] for _ in 1:n]
    by_comp = Dict{String,Vector{Int}}()
    for (i, occ) in enumerate(occs)
        push!(get!(() -> Int[], by_comp, occ.comp), i)
        push!(neighbors[i], (3, n_occ + lab_index[occ.label]))
        push!(neighbors[n_occ + lab_index[occ.label]], (occ.sign > 0 ? 4 : 5, i))
    end
    for idxs in values(by_comp)
        sort!(idxs, by=i -> occs[i].pos)
        m = length(idxs)
        for t in 1:m
            push!(neighbors[idxs[t]], (1, idxs[mod1(t + 1, m)]))
            push!(neighbors[idxs[t]], (2, idxs[mod1(t - 1, m)]))
        end
    end
    colors = _hc_rank!(vcat([occ.init for occ in occs], [(:lab,) for _ in labels]))
    for _ in 1:n + 2
        signatures = Vector{Any}(undef, n)
        for i in 1:n
            bag = Tuple{Int,Int}[]
            for (edge, j) in neighbors[i]
                push!(bag, (edge, colors[j]))
            end
            sort!(bag)
            signatures[i] = (colors[i], Tuple(bag))
        end
        new_colors = _hc_rank!(signatures)
        new_colors == colors && break
        colors = new_colors
    end
    comp_colors = Dict{String,Vector{Int}}()
    for (name, idxs) in by_comp
        sort!(idxs, by=i -> occs[i].pos)
        comp_colors[name] = [colors[i] for i in idxs]
    end
    comp_colors, colors
end

function _hc_comp_key(cycle::Vector{Int})
    isempty(cycle) && return Tuple{Int}()
    best = cycle
    for start in 2:length(cycle)
        rotated = vcat(cycle[start:end], cycle[1:start-1])
        rotated < best && (best = rotated)
    end
    Tuple(best)
end

function _hc_tie_groups(names::Vector{String}, comp_colors)
    keyed = Pair{Tuple,String}[]
    for name in names
        cycle = get(comp_colors, name, Int[])
        push!(keyed, _hc_comp_key(cycle) => name)
    end
    sort!(keyed, by=pair -> pair.first)
    groups = Pair{Tuple,Vector{String}}[]
    for (key, name) in keyed
        if isempty(groups) || groups[end].first != key
            push!(groups, key => String[name])
        else
            push!(groups[end].second, name)
        end
    end
    groups
end

function _hc_split_empty(names::Vector{String}, words)
    present = String[]
    empty = String[]
    for name in names
        isempty(words[name]) ? push!(empty, name) : push!(present, name)
    end
    present, empty
end

function _hc_transport_ref(ref::AbstractString, name_map, face_map, label_map)
    text = String(strip(ref))
    text == "outer" && return "outer"
    if occursin('/', text)
        curve, signed = split(text, '/', limit=2)
        lab = parse(Int, signed)
        haskey(label_map, abs(lab)) || throw(KodaInputError(
            :hc_canonical, "point names an unknown label", Dict("ref" => text)))
        new_lab = sign(lab) * label_map[abs(lab)]
        new_curve = get(name_map, String(curve), String(curve))
        return string(new_curve, "/", new_lab > 0 ? "+" : "", new_lab)
    end
    if length(text) >= 2 && text[1] == 'F' && all(isdigit, text[2:end])
        return get(face_map, text, text)
    end
    if !isempty(text) && (endswith(text, "+") || endswith(text, "-"))
        side = text[end:end]
        curve = text[1:end-1]
        return string(get(name_map, curve, curve), side)
    end
    get(name_map, text, text)
end

function _hc_face_map(words, label_map)
    active = Dict{String,Vector{Int}}()
    for (name, word) in words
        isempty(word) && continue
        active[String(name)] = Int[sign(x) * label_map[abs(x)] for x in word]
    end
    isempty(active) && return Dict{String,String}()
    original_faces, _ = _lc_faces(words)
    _, new_corner_face = _lc_faces(active)
    out = Dict{String,String}()
    for face in original_faces
        isempty(face.corners) && continue
        corner = face.corners[1]
        haskey(label_map, corner[1]) || throw(KodaImplementationError(
            :hc_canonical, "face corner was not labeled", Dict("label" => corner[1])))
        transported = (label_map[corner[1]], corner[2], corner[3])
        haskey(new_corner_face, transported) || throw(KodaImplementationError(
            :hc_canonical, "transported face corner is missing", Dict("corner" => string(transported))))
        out[face.id] = new_corner_face[transported]
    end
    out
end

function _hc_disk_pairs(dict, transport_ref)
    pairs = Pair{String,Vector{String}}[]
    for name in sort!(collect(String.(keys(dict))))
        push!(pairs, name => String[transport_ref(ref) for ref in dict[name]])
    end
    pairs
end

function _hc_point_pairs(dict, transport_ref)
    pairs = Pair{String,String}[]
    for name in sort!(collect(String.(keys(dict))))
        push!(pairs, name => transport_ref(dict[name]))
    end
    pairs
end

function _hc_weight_pairs(dict, transport_ref)
    pairs = Pair{String,Int}[]
    for name in sort!(collect(String.(keys(dict))))
        push!(pairs, transport_ref(name) => Int(dict[name]))
    end
    sort!(pairs, by=pair -> pair.first)
    pairs
end

function _hc_annotation(hc::HeegaardLayer, name_map, face_map, label_map)
    transport_ref = ref -> _hc_transport_ref(ref, name_map, face_map, label_map)
    places = Pair{String,String}[]
    for name in sort!(collect(String.(keys(hc.lc.places))))
        canonical = get(name_map, name, name)
        push!(places, canonical => transport_ref(hc.lc.places[name]))
    end
    sort!(places, by=pair -> pair.first)
    places,
    _hc_disk_pairs(hc.lc.disks, transport_ref), _hc_disk_pairs(hc.disks, transport_ref),
    _hc_point_pairs(hc.lc.points, transport_ref), _hc_point_pairs(hc.points, transport_ref),
    _hc_weight_pairs(hc.lc.weights, transport_ref), _hc_weight_pairs(hc.weights, transport_ref)
end

function _hc_has_annotations(hc::HeegaardLayer)
    !isempty(hc.lc.places) || !isempty(hc.lc.disks) || !isempty(hc.disks) ||
        !isempty(hc.lc.points) || !isempty(hc.points) ||
        !isempty(hc.lc.weights) || !isempty(hc.weights)
end

function _hc_name_map(k_order, alpha_order, beta_order, other_order)
    name_map = Dict{String,String}()
    for name in k_order
        name_map[name] = name
    end
    for (i, name) in enumerate(alpha_order)
        name_map[name] = "a$i"
    end
    for (i, name) in enumerate(beta_order)
        name_map[name] = "b$i"
    end
    for name in other_order
        name_map[name] = name
    end
    name_map
end

function _hc_best_empty_order(names::Vector{String}, hc, base_map, face_map, label_map, family_prefix::String)
    isempty(names) && return String[]
    length(names) > 8 && throw(HCNotComputed(
        :empty_curve_order, Dict{String,Any}("n" => length(names), "family" => family_prefix)))
    best = nothing
    best_order = names
    for perm in (length(names) == 1 ? [names] : _hc_permutations(names))
        trial = copy(base_map)
        for (i, name) in enumerate(perm)
            trial[name] = string(family_prefix, "e", i)
        end
        record = String[]
        for name in perm
            place = get(hc.lc.places, name, "")
            push!(record, place == "" ? "" : _hc_transport_ref(place, trial, face_map, label_map))
        end
        if best === nothing || record < best
            best = record
            best_order = copy(perm)
        end
    end
    best_order
end
