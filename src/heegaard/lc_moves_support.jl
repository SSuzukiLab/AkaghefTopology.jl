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

"""Certify a route using the common right-face engine.
Legacy corner triples denote the LEFT face. Explicit dictionaries contain a
RIGHT-face `dart` and its rational `u`, and permit ordered repeated-edge visits.
The returned `route` can be supplied as a handle-slide band.
"""
function lc_move_path(lc::LickorishCode; start, crossings=Any[], finish, avoid=String[])
    validate_lc_input(lc)
    active=_lc_nonempty_words(lc.words);occ=_lc_rotation_occ(active)
    point(x)=x isa AbstractDict ? x : Dict("dart"=>collect(_lc_tau(active,occ,_move_corner(x))),"u"=>"1/2")
    band=Dict("start"=>point(start),"finish"=>point(finish),"crossings"=>[point(x) for x in crossings])
    route=_move_certify_route(lc,band;avoid=avoid)
    steps=[Dict("from"=>route.faces[i],"to"=>route.faces[i+1],"component"=>p.component,
        "gap"=>p.edge+1,"dart"=>collect(p.dart),"u"=>string(p.u)) for (i,p) in enumerate(route.crossings)]
    Dict{String,Any}("kind"=>"disk_face_path","faces"=>route.faces,"crossings"=>steps,
        "start"=>band["start"],"finish"=>band["finish"],"route"=>band,
        "range"=>"ordered routes with nonalternating chords in disk faces")
end

function _move_require_separated(input, components)
    input isa HeegaardCode || return
    shared = [n for n in components if _is_k(input.lc,n) && haskey(input.roles,n)]
    isempty(shared) || throw(LCMoveError(:shared_support_requires_separation,
        "moving a role without moving K requires separating its shared support first",
        Dict("components"=>shared)))
end

# A face not split by the support has a unique image via its surviving corners.
# Ambiguity is an absent transport choice, never a reason to reuse its F number.
function _move_transport_places(lc, words, places, label_map;
        digon_create_untouched_corners=Dict{String,Vector{Tuple{Int,String,String}}}())
    isempty(lc.places) && return places
    before,_ = _lc_faces(_lc_nonempty_words(lc.words))
    _,after = _lc_faces(_lc_nonempty_words(words))
    out=copy(places)
    for (name,ref) in lc.places
        isempty(get(words,name,Int[0])) || continue
        startswith(ref,"F") || continue
        fs=[f for f in before if f.id==ref]
        length(fs)==1 || throw(LCMoveError(:missing_place_transport,"place face is not realized",Dict("component"=>name,"face"=>ref)))
        # cf: .research-scratch/PJ260922/astra-lanes/hc-moves/answer-4-3-harvest.md —
        # roots and their nested descendants can be parked in an untouched old
        # corner sector. The source constructor uses open-edge ports u=1/2 and
        # a closed band avoiding every old vertex; a surviving label alone is insufficient.
        safe=get(digon_create_untouched_corners,ref,nothing)
        if safe!==nothing
            !isempty(safe) && all(c->c in only(fs).corners,safe) || throw(LCMoveError(
                :invalid_creation_parking,"parking corner is absent from its source face",Dict("face"=>ref)))
            mapped=[get(after,(get(label_map,c[1],0),c[2],c[3]),nothing) for c in safe]
            all(x->x!==nothing,mapped) || throw(LCMoveError(:invalid_creation_parking,
                "an untouched old corner did not survive",Dict("face"=>ref)))
            out[name]=first(sort!(String.(mapped);by=_move_ref_key))
            continue
        end
        images=Set{String}()
        for (v,s,d) in only(fs).corners
            haskey(label_map,v) || continue
            image=get(after,(label_map[v],s,d),nothing)
            image===nothing || push!(images,image)
        end
        length(images)==1 || throw(LCMoveError(:missing_place_transport,
            "place needs a choice in a split or removed face",Dict("component"=>name,"face"=>ref,"images"=>collect(images))))
        out[name]=only(images)
    end
    out
end

export lc_move_path
# cf: docs/specs/lickorish-code.md §3 — moves require a connected cellular input and placed empty curves.
function _move_validate_input(lc)
    _lc_validate_surface!(lc)
    for (name,word) in lc.words
        isempty(word) && !haskey(lc.places,name) && throw(LCMoveError(
            :invalid_input_place,"empty component has no place",Dict("component"=>name)))
        !isempty(word) && haskey(lc.places,name) && throw(LCMoveError(
            :invalid_input_place,"nonempty component has a place",Dict("component"=>name)))
    end
    active=_lc_nonempty_words(lc.words)
    isempty(active) && return
    remaining=Set(keys(active))
    reached=Set([first(remaining)])
    labels=Set(abs.(active[first(reached)]))
    changed=true
    while changed
        changed=false
        for name in remaining
            name in reached && continue
            any(abs(x) in labels for x in active[name]) || continue
            push!(reached,name); union!(labels,abs.(active[name])); changed=true
        end
    end
    length(reached)==length(remaining) || throw(LCMoveError(
        :invalid_input_disconnected,"curve graph is disconnected",Dict{String,Any}()))
end
