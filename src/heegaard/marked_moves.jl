"""Constructive moves on the native whole marked document.

Source: vendor/marked_surface_diagrams/src/MarkedSurfaceDiagrams.jl,
`DiagramState`, `parse_word`, `regions`, `move_vertex`, `certify_layout`.
The vendor is unchanged. These adapters preserve auxiliary-safe D semantics
through marked_topology.jl::_marked_states; no intersection-only replacement
or caller-supplied post-move diagram is used. Current combinatorial isotopies
are empty interior-triangle moves, not a complete isotopy decision procedure.
"""
_mm_hash(x) = bytes2hex(sha256(_canonical_json(x)))
_mm_need(ok, message) = ok ? nothing : throw(ArgumentError(message))
_mm_cell(x) = (String(x[1]), Int(x[2]))
_mm_cells(region) = Set(_mm_cell(x) for x in region["triangles"])
_mm_triangle(t) = minimum((Tuple(t), Tuple(circshift(t, 1)), Tuple(circshift(t, 2))))

function _mm_bare(doc)
    bare = deepcopy(doc)
    delete!(bare, "region_marks")
    bare
end

function _mm_inverse_seeds(doc)
    Dict{String,Any}(mark["id"] => deepcopy(mark["seed"]) for mark in get(doc, "region_marks", []))
end

function _mm_infinity(doc)
    get(get(get(doc, "metadata", Dict()), "display", Dict()), "infinity_triangle", nothing)
end

# Source: vendor `regions` returns zero-based triangle occurrences. A triangle
# isotopy identifies every region by cells outside its swept support. This is
# a geometric transport test; an unchanged seed or a hash is not a substitute.
function _mm_region_map(before, after, image; omit=Set{Tuple{String,Int}}())
    target = [_mm_cells(r) for r in after]
    out = Any[]
    for (i, region) in enumerate(before)
        cells = Set(image(c) for c in _mm_cells(region) if !(c in omit))
        _mm_need(!isempty(cells), "move support consumes an entire region; explicit transport is required")
        matches = findall(s -> !isdisjoint(cells, s), target)
        _mm_need(length(matches) == 1 && issubset(cells, target[only(matches)]),
                 "region transport is ambiguous or splits an existing region")
        j = only(matches)
        push!(out, Dict("source_region" => i, "target_region" => j,
            "source_triangles" => [collect(c) for c in sort(collect(_mm_cells(region)))],
            "target_triangles" => [collect(c) for c in sort(collect(target[j]))]))
    end
    _mm_need(length(out) == length(after) && length(unique(x["target_region"] for x in out)) == length(after),
             "move does not biject complementary regions")
    out
end

function _mm_transport_regions!(doc, before, image; omit=Set{Tuple{String,Int}}(), specification=Dict())
    old = MarkedSurfaceDiagrams.regions(before.mathematics)
    next = _marked_states(_mm_bare(doc))
    new = MarkedSurfaceDiagrams.regions(next.mathematics)
    correspondence = _mm_region_map(old, new, image; omit=omit)
    function seed_image(seed, source_regions, target_regions, mapping, preferred=nothing)
        source = _mm_cell(seed)
        i = findall(r -> source in _mm_cells(r), source_regions)
        _mm_need(length(i) == 1, "selected source seed is not a unique mathematical region")
        j = mapping[only(i)]["target_region"]
        target = _mm_cells(target_regions[j])
        choice = preferred === nothing ? image(source) : _mm_cell(preferred)
        if preferred === nothing && (source in omit || !(choice in target))
            stable = sort([image(c) for c in _mm_cells(source_regions[only(i)]) if !(c in omit) && image(c) in target])
            _mm_need(!isempty(stable), "selected region has no stable seed outside support")
            choice = first(stable)
        end
        _mm_need(choice in target, "requested inverse seed is outside the transported region")
        collect(choice)
    end
    preferred = get(specification, "region_seeds", Dict())
    for mark in get(doc, "region_marks", [])
        mark["seed"] = seed_image(mark["seed"], old, new, correspondence, get(preferred, mark["id"], nothing))
    end
    infinity = _mm_infinity(before.document)
    if infinity !== nothing
        # Infinity belongs to the full displayed geometry, including aux.
        # D deliberately belongs to alpha union beta above. These transports
        # must differ when an auxiliary loop splits a mathematical region.
        old_geometry=MarkedSurfaceDiagrams.regions(before.geometry)
        new_geometry=MarkedSurfaceDiagrams.regions(next.geometry)
        geometry_map=_mm_region_map(old_geometry,new_geometry,image;omit=omit)
        doc["metadata"]["display"]["infinity_triangle"] =
            seed_image(infinity, old_geometry, new_geometry, geometry_map, get(specification, "infinity_triangle", nothing))
    end
    _marked_states(doc) # D claims are checked on alpha union beta, not aux.
    correspondence
end

function _mm_identity_cells(st)
    Dict("cell_semantics" => "identity of the unchanged background triangulation; not cellwise ambient isotopy transport",
         "vertices" => [Dict("source" => [f, v], "target" => [f, v])
            for (f,v) in sort(collect(keys(st.vertex_dsu.p)))],
         "edges" => [Dict("source" => [f, e...], "target" => [f, e...])
            for f in sort(collect(keys(st.patches))) for e in sort(collect(keys(st.incidence[f])))],
         "triangles" => [Dict("source" => [f, k-1], "target" => [f, k-1])
            for (f, k) in sort(collect(keys(st.triangles)))])
end

function _mm_witness(kind, before, after, specification, inverse_specification, correspondence)
    Dict{String,Any}("schema" => "marked-move/1", "kind" => kind,
        "source_document" => _mm_hash(before), "target_document" => _mm_hash(after),
        "source_marking" => MarkedSurfaceDiagrams.marking_digest(before),
        "target_marking" => MarkedSurfaceDiagrams.marking_digest(after),
        "specification" => deepcopy(specification), "inverse_specification" => deepcopy(inverse_specification),
        "correspondence" => correspondence)
end

"""Expand/contract a curve through one empty interior triangle.

The specification is a JSON-compatible Dict with `curve`, `face`, zero-based
`triangle`, and oriented `path=[u,v]` (expand) or `[u,w,v]` (contract).
The disk-supported isotopy fixes the other curves and marked patch boundaries.
Source: vendor `Dart` and disk triangulation/link checks; the empty-triangle
isotopy construction and region-overlap transport are added here by the agent.
No crossing births/deaths, seam moves, or general isotopy recognition is claimed.
"""
function marked_triangle_isotopy(input::AbstractDict, specification::AbstractDict)
    before = _marked_states(_marked_document(input)); st = before.geometry
    doc = deepcopy(before.document); spec = _plain_value(specification)
    M = MarkedSurfaceDiagrams
    face = String(spec["face"]); k = Int(spec["triangle"]); cid = String(spec["curve"])
    path = String.(spec["path"])
    _mm_need(length(path) in (2, 3) && length(unique(path)) == length(path), "path must have two or three distinct vertices")
    _mm_need(haskey(st.triangles, (face, k+1)), "unknown support triangle")
    triangle = st.triangles[(face, k+1)]
    _mm_need(all(v -> !(v in st.boundaries[face]), triangle), "triangle isotopy support must be interior to a marked patch")
    _mm_need(issubset(Set(path), Set(triangle)), "path is not on the selected triangle")
    _mm_need(haskey(st.curves, cid), "unknown curve")
    curve, darts = st.curves[cid]
    source = [M.Dart(face, a, b) for (a, b) in zip(path[1:end-1], path[2:end])]
    where = findall(j -> j+length(source)-1 <= length(darts) && darts[j:j+length(source)-1] == source, eachindex(darts))
    _mm_need(length(where) == 1, "selected path is not a unique contiguous curve subpath (word-start wrap is not supported)")
    j = only(where)
    middle = only(setdiff(Set(triangle), Set((path[1], path[end]))))
    alternative = length(path) == 2 ? [path[1], middle, path[end]] : [path[1], path[end]]
    replacement = [M.Dart(face, a, b) for (a,b) in zip(alternative[1:end-1], alternative[2:end])]
    # All occupied support edges must be precisely the selected subpath.
    selected_edges = Set(M.edge(d.u,d.v) for d in source)
    for (a,b) in M.tedges(triangle)
        owners = get(st.used, (face, M.edge(a,b)), String[])
        _mm_need(isempty(owners) || (M.edge(a,b) in selected_edges && owners == [cid]), "support triangle contains another occupied edge")
    end
    for (other, (_, ds)) in st.curves
        touched = any(d -> d.face == face && (d.u == middle || d.v == middle), ds)
        _mm_need(!touched || (length(path) == 3 && other == cid), "alternate vertex is occupied")
    end
    changed = vcat(darts[1:j-1], replacement, darts[j+length(source):end])
    word = M.print_word(changed)
    if haskey(spec,"word_format")
        _mm_need(M.parse_word(spec["word_format"]) == changed, "inverse word formatting changes the constructed route")
        word = spec["word_format"]
    end
    only(c for c in doc["curves"] if c["id"] == cid)["word"] = word
    support = Set([(face,k)])
    region_map = _mm_transport_regions!(doc,before,identity;omit=support,specification=spec)
    after = _marked_states(doc)
    _mm_need(doc["crossings"] == before.document["crossings"], "triangle isotopy changed crossing incidence")
    inverse = Dict{String,Any}("curve"=>cid,"face"=>face,"triangle"=>k,"path"=>alternative,
        "region_seeds"=>_mm_inverse_seeds(before.document),"word_format"=>curve["word"])
    _mm_infinity(before.document) !== nothing && (inverse["infinity_triangle"] = deepcopy(_mm_infinity(before.document)))
    correspondence = _mm_identity_cells(st)
    correspondence["regions"] = region_map
    correspondence["geometric_regions"] = _mm_region_map(M.regions(st),M.regions(after.geometry),identity;omit=support)
    correspondence["curves"] = [Dict("source"=>c["id"],"target"=>c["id"],
        "source_word"=>c["word"],"target_word"=>only(x for x in doc["curves"] if x["id"]==c["id"])["word"]) for c in before.document["curves"]]
    correspondence["crossings"] = [Dict("source"=>x["id"],"target"=>x["id"]) for x in doc["crossings"]]
    correspondence["isotopy_support"] = Dict("triangle"=>[face,k],"source_path"=>path,"target_path"=>alternative)
    (;document=doc,witness=_mm_witness("triangle_isotopy",before.document,doc,spec,inverse,correspondence))
end

_mm_rational(x) = x isa Rational ? BigInt(numerator(x))//BigInt(denominator(x)) : MarkedSurfaceDiagrams.qnum(x)
_mm_positions(pos) = Dict(String(f)=>Dict(String(v)=>Tuple(_mm_rational.(p)) for (v,p) in vs) for (f,vs) in pos)
_mm_position_json(pos) = Dict(f=>Dict(v=>MarkedSurfaceDiagrams.qstr.(collect(p)) for (v,p) in vs) for (f,vs) in pos)

"""Exact one-interior-vertex ambient isotopy with replayable rational coordinates.

This is a realization operation: its mathematical document remains the same.
Source: unchanged vendor `move_vertex` and cellulated vendor docs/PROOFS_JA.md,
section "一頂点の線形移動は同位を与える". Only one vertex moves at a time.
"""
function marked_vertex_isotopy(input::AbstractDict, positions::AbstractDict, specification::AbstractDict)
    before = _marked_states(_marked_document(input)); spec = _plain_value(specification)
    M = MarkedSurfaceDiagrams; pos = _mm_positions(positions)
    face=String(spec["face"]); vertex=String(spec["vertex"])
    out, local_witness = M.move_vertex(before.geometry,pos,face,vertex,spec["target"])
    inverse=Dict("face"=>face,"vertex"=>vertex,"target"=>M.qstr.(collect(pos[face][vertex])))
    correspondence=_mm_identity_cells(before.geometry)
    correspondence["realization_isotopy"]=local_witness
    witness=_mm_witness("vertex_isotopy",before.document,before.document,spec,inverse,correspondence)
    witness["source_positions"]=_mm_hash(_mm_position_json(pos))
    witness["target_positions"]=_mm_hash(_mm_position_json(out))
    (;document=deepcopy(before.document),positions=out,witness)
end

"""Apply an explicitly verified oriented automorphism of the fixed patch mesh.

`specification["vertices"][face][vertex]` is a local vertex permutation, with
one full permutation per existing patch. It must preserve oriented triangles,
named-side incidence and paired seams; physical boundary points stay fixed.
Curves, crossing points, region seeds and infinity are constructed as images.
The target triangulation/marking stays fixed. This supports finite simplicial
actions, not a Dehn-twist constructor or arbitrary mapping-class recognition.
Source: vendor patch_info/_validate gluing, triangle orientation and Dart rules.
"""
function marked_simplicial_action(input::AbstractDict, specification::AbstractDict)
    before=_marked_states(_marked_document(input));st=before.geometry;M=MarkedSurfaceDiagrams
    doc=deepcopy(before.document);spec=_plain_value(specification);maps=spec["vertices"]
    _mm_need(Set(keys(maps))==Set(keys(st.patches)),"simplicial action requires every existing patch")
    side_image=Dict();triangle_image=Dict{Tuple{String,Int},Tuple{String,Int}}();inverse_maps=Dict{String,Any}()
    for (f,p) in st.patches
        vm=maps[f];vertices=Set(v for t in p["triangles"] for v in t)
        _mm_need(Set(keys(vm))==vertices && Set(values(vm))==vertices,"vertex map must be a local bijection")
        inverse_maps[f]=Dict(v=>u for (u,v) in vm)
        index=Dict(_mm_triangle(t)=>i-1 for (i,t) in enumerate(p["triangles"]))
        for (i,t) in enumerate(p["triangles"])
            image=_mm_triangle([vm[v] for v in t])
            _mm_need(haskey(index,image),"vertex map does not preserve oriented triangles")
            triangle_image[(f,i-1)]=(f,index[image])
        end
        for s in p["side_order"]
            path=[vm[v] for v in p["sides"][s]]
            matches=[t for t in p["side_order"] if p["sides"][t]==path]
            _mm_need(length(matches)==1,"vertex map does not preserve oriented named-side paths")
            target=only(matches)
            _mm_need(M.side_parameters(p,s)==M.side_parameters(p,target),"side parameters are not preserved")
            side_image[(f,s)]=(f,target)
        end
    end
    seam=Dict{Tuple{String,String},Tuple{String,String}}()
    for pair in doc["surface"]["gluings"]
        a,b=Tuple(pair[1]),Tuple(pair[2]);seam[a]=b;seam[b]=a
    end
    for (a,b) in seam
        _mm_need(get(seam,side_image[a],nothing)==side_image[b],"surface action does not descend across paired seams")
        pa=st.patches[a[1]]["sides"][a[2]];pb=st.patches[b[1]]["sides"][b[2]]
        for (u,v) in zip(pa,reverse(pb))
            _mm_need(M.root!(st.vertex_dsu,(a[1],maps[a[1]][u]))==M.root!(st.vertex_dsu,(b[1],maps[b[1]][v])),"mapped seam vertices disagree")
        end
    end
    for boundary in doc["surface"]["boundary"]
        f,s=boundary["side"]
        _mm_need(all(maps[f][v]==v for v in st.patches[f]["sides"][s]),"physical boundary must be pointwise fixed")
    end
    for curve in doc["curves"]
        ds=M.parse_word(curve["word"])
        curve["word"]=M.print_word([M.Dart(d.face,maps[d.face][d.u],maps[d.face][d.v]) for d in ds])
    end
    for x in doc["crossings"]
        f,v=x["at"];x["at"]=[f,maps[f][v]]
    end
    regions=_mm_transport_regions!(doc,before,c->triangle_image[c];specification=spec)
    inverse=Dict{String,Any}("vertices"=>inverse_maps,"region_seeds"=>_mm_inverse_seeds(before.document))
    _mm_infinity(before.document)!==nothing && (inverse["infinity_triangle"]=deepcopy(_mm_infinity(before.document)))
    # Lexical formatting is supplementary inverse information, accepted only
    # after the inverse map has constructed exactly the same directed darts.
    if haskey(spec,"word_formats")
        for c in doc["curves"]
            word=spec["word_formats"][c["id"]]
            _mm_need(M.parse_word(word)==M.parse_word(c["word"]),"inverse word format changes the mapped route")
            c["word"]=word
        end
    end
    inverse["word_formats"]=Dict(c["id"]=>c["word"] for c in before.document["curves"])
    correspondence=Dict{String,Any}(
        "vertices"=>[Dict("source"=>[f,v],"target"=>[f,maps[f][v]]) for f in sort(collect(keys(maps))) for v in sort(collect(keys(maps[f])))],
        "edges"=>[Dict("source"=>[f,e...],"target"=>[f,maps[f][e[1]],maps[f][e[2]]]) for f in sort(collect(keys(st.patches))) for e in sort(collect(keys(st.incidence[f])))],
        "triangles"=>[Dict("source"=>collect(c),"target"=>collect(triangle_image[c])) for c in sort(collect(keys(triangle_image)))],
        "regions"=>regions,"curves"=>[Dict("source"=>c["id"],"target"=>c["id"],"source_word"=>c["word"],"target_word"=>only(x for x in doc["curves"] if x["id"]==c["id"])["word"]) for c in before.document["curves"]],
        "crossings"=>[Dict("source"=>x["id"],"target"=>x["id"],"source_at"=>x["at"],"target_at"=>only(y for y in doc["crossings"] if y["id"]==x["id"])["at"]) for x in before.document["crossings"]])
    _marked_states(doc)
    (;document=doc,witness=_mm_witness("simplicial_action",before.document,doc,spec,inverse,correspondence))
end

"""Return inverse construction instructions, never a retained source snapshot."""
function inverse_marked_move(witness::AbstractDict)
    w=_plain_value(witness)
    _mm_need(get(w,"schema","")=="marked-move/1","unknown marked move witness schema")
    result=deepcopy(w)
    if haskey(w,"inverse_kind")
        result["kind"],result["inverse_kind"]=w["inverse_kind"],w["kind"]
    end
    for (a,b) in (("source_document","target_document"),("source_marking","target_marking"),("specification","inverse_specification"))
        result[a],result[b]=deepcopy(w[b]),deepcopy(w[a])
    end
    if haskey(w,"source_positions")
        result["source_positions"],result["target_positions"]=w["target_positions"],w["source_positions"]
    end
    delete!(result,"correspondence") # inverse replay reconstructs its own map.
    result
end

"""Reconstruct a move and verify both source and target, including metadata.

For exact realization motion pass `positions`; only their identities are in
the witness. JSON serialization round trips instructions without live objects.
"""
function _mm_execute(kind,doc,spec,positions)
    if kind=="triangle_isotopy"
        marked_triangle_isotopy(doc,spec)
    elseif kind=="simplicial_action"
        marked_simplicial_action(doc,spec)
    elseif kind=="vertex_isotopy"
        _mm_need(positions!==nothing,"vertex isotopy replay requires source positions")
        marked_vertex_isotopy(doc,positions,spec)
    elseif kind=="barycentric_subdivision"
        marked_subdivide(doc,spec)
    elseif kind=="barycentric_collapse"
        marked_collapse_subdivision(doc,spec)
    elseif kind=="dual_star"
        marked_dual_star(doc,spec)
    elseif kind=="dual_star_collapse"
        marked_collapse_dual_star(doc,spec)
    elseif kind=="parallel_copy"
        marked_parallel_copy(doc,spec)
    elseif kind=="parallel_remove"
        marked_remove_parallel(doc,spec)
    elseif kind=="handle_slide"
        marked_handle_slide(doc,spec)
    elseif kind=="handle_slide_inverse"
        marked_inverse_handle_slide(doc,spec)
    elseif kind=="stabilization"
        marked_stabilize(doc,spec)
    elseif kind=="generated_destabilization"
        marked_destabilize(doc,spec)
    elseif kind=="dehn_twist"
        marked_dehn_twist(doc,spec)
    elseif kind=="generated_twist_inverse"
        marked_inverse_dehn_twist(doc,spec)
    else
        throw(ArgumentError("unknown marked move kind $kind"))
    end
end

function replay_marked_move(input::AbstractDict,witness::AbstractDict;positions=nothing)
    w=_plain_value(witness);doc=_marked_document(input)
    _mm_need(get(w,"schema","")=="marked-move/1","unknown marked move witness schema")
    _mm_need(_mm_hash(doc)==w["source_document"],"marked move source document does not match")
    spec=w["specification"];kind=w["kind"]
    if kind=="vertex_isotopy"
        _mm_need(positions!==nothing,"vertex isotopy replay requires source positions")
        _mm_need(_mm_hash(_mm_position_json(_mm_positions(positions)))==w["source_positions"],"marked move source positions do not match")
    end
    result=_mm_execute(kind,doc,spec,positions)
    _mm_need(_mm_hash(result.document)==w["target_document"],"reconstructed marked move target does not match")
    for key in ("source_marking","target_marking")
        _mm_need(result.witness[key]==w[key],"reconstructed marked move $key does not match")
    end
    haskey(w,"correspondence") && _mm_need(_canonical_json(w["correspondence"])==_canonical_json(result.witness["correspondence"]),"supplied correspondence differs from the reconstructed move")
    kind=="vertex_isotopy" && _mm_need(result.witness["target_positions"]==w["target_positions"],"reconstructed target positions do not match")
    restored=_mm_execute(get(w,"inverse_kind",kind),result.document,w["inverse_specification"],kind=="vertex_isotopy" ? result.positions : nothing)
    _mm_need(_mm_hash(restored.document)==w["source_document"],"supplied inverse instructions do not restore the source")
    kind=="vertex_isotopy" && _mm_need(restored.witness["target_positions"]==w["source_positions"],"supplied inverse does not restore source positions")
    result
end
