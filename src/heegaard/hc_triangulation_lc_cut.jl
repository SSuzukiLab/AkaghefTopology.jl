# cf: docs/specs/heegaard-snappy.md — genus-one cut-disk transport.
# cf: docs/specs/lickorish-code.md — passage signs and fixed reference words.
# Algorithm: Sol Pro answer-4.md, §§genus-one construction, passage sign, face correspondence.
# The attachment was unavailable; this is a local implementation of the answer's explicit algorithm.

_hc_cut_det(a,b) = a[1]*b[2]-a[2]*b[1]
_hc_cut_sub(a,b) = (a[1]-b[1],a[2]-b[2])
_hc_cut_lerp(a,b,t) = (a[1]+t*(b[1]-a[1]),a[2]+t*(b[2]-a[2]))

function _hc_cut_intersection(a,b,c,d)
    u,v,w=_hc_cut_sub(b,a),_hc_cut_sub(d,c),_hc_cut_sub(c,a)
    det=_hc_cut_det(u,v)
    det==0 && return nothing
    t,s=_hc_cut_det(w,v)/det,_hc_cut_det(w,u)/det
    (;t,s,sign=Int(sign(det)))
end

include("hc_triangulation_lc_cut_overlay.jl")

"""Transport a checked genus-one or standard genus-two Koda surface into a resolved LC/HC.

The old ribbon graph is retained. Only the reference curves are inserted;
all old crossings and the specified joining region have checked ancestry.
Other polygon words and higher genus require additional transport and are unsupported.
"""
function hc_cut_to_resolved(diagram, cut=hc_surface_cut_system(diagram); joining_region::Integer)
    cut==hc_surface_cut_system(diagram) ||
        throw(HCTriangulationError(:surface_cut_mismatch,"cut certificate differs from the source diagram"))
    g=Int(cut["genus"])
    g in (1,2) || throw(HCTriangulationError(:general_polygon_slides_pending,"genus $g requires arc-carrying polygon slides"))
    graph=ribbon_graph(diagram); code=_undecorated(diagram)
    boundary=Int.(cut["disk_boundary_darts"]); cotree=Int.(cut["dual_tree_edges"])
    word=Int.(cut["polygon_word"])
    length(word)==4g && all(word[4i-1]==-word[4i-3] && word[4i]==-word[4i-2] for i in 1:g) ||
        throw(HCTriangulationError(:polygon_not_commutator,string(word)))
    sourcefaces=boundary_cycles(graph)
    1<=joining_region<=length(sourcefaces) || throw(ArgumentError("joining region out of range"))
    joining_disk(diagram,joining_region) # Validate the given disk, never choose another.
    outgoing=Dict((d.origin,d.origin_port)=>d.id for d in graph.darts)
    # Blow up each retained edge side and each intervening vertex sector.
    lanes=Dict{Int,Vector{Int}}(); side=Dict{Int,Tuple{Int,Int}}(); port=Dict{Int,Int}()
    n=0
    for (i,did) in enumerate(boundary)
        d=graph.darts[did]; next=graph.darts[boundary[mod1(i+1,length(boundary))]]
        d.target==next.origin || throw(HCTriangulationError(:boundary_vertex_mismatch,string(did)))
        side[did]=(n+1,n+2); n+=2
        rot=graph.rotation[d.target]; j=findfirst(==(d.target_port),rot)
        lane=Int[]
        for step in 1:4
            p=rot[mod1(j-step,4)]
            p==next.origin_port && break
            half=outgoing[(d.target,p)]
            graph.darts[half].edge in cotree ||
                throw(HCTriangulationError(:boundary_lane_mismatch,"non-cotree half-dart $half"))
            push!(lane,half);n+=1;port[half]=n
        end
        lanes[did]=lane
    end
    Set(keys(port))==Set(vcat([[2*e-1,2*e] for e in cotree]...)) ||
        throw(HCTriangulationError(:cotree_ports_missing,"each cotree half-dart must occur exactly once"))
    points=[(BigInt(i)//1,BigInt(i)^2//1) for i in 1:n]
    chords=Dict(e=>(points[port[2*e-1]],points[port[2*e]]) for e in cotree)
    between(x,a,b)=mod(x-a,n)<mod(b-a,n)
    for (i,e) in enumerate(cotree), f in cotree[i+1:end]
        xor(between(port[2*f-1],port[2*e-1],port[2*e]),between(port[2*f],port[2*e-1],port[2*e])) &&
            throw(HCTriangulationError(:cotree_chords_cross,"alternating cotree endpoints"))
    end
    endpoint(s,u)=begin
        dart=s>0 ? 2s-1 : -2s
        i,j=side[dart];_hc_cut_lerp(points[i],points[j],u)
    end
    bound=8(length(cotree)+2)^2+17
    layout=_hc_cut_overlay_layout(graph,cut,side,port,points,chords,endpoint;actual_geometry=true)
    labels,events,words,chosen,trials=layout.labels,layout.events,layout.words,layout.chosen,layout.trials
    eps=Dict(v.id=>v.sign for v in code.intersections)
    roles=Dict{String,Symbol}();subdivision=Dict{Int,Vector{Int}}()
    sourcecurves=Dict{String,Vector{Int}}()
    for curve in _curve_vector(code)
        name=string(curve.family,curve.id); w=Int[];roles[name]=curve.family
        sourcecurves[name]=copy(curve.visits)
        for v in curve.visits
            e=graph.darts[outgoing[(v,_port(curve.family,:out))]].edge
            token=(curve.family==:alpha ? -eps[v] : eps[v])*labels[v]
            segment=vcat([token],[x.token for x in sort(events[e];by=x->x.t)])
            append!(w,segment);subdivision[e]=abs.(vcat(segment,[labels[graph.edges[e].to]]))
        end
        words[name]=w
    end
    for name in _lc_k_names(g)
        [x for x in words[name] if abs(x)<=3g-2]==_lc_k_constants(g)[name] ||
            throw(HCTriangulationError(:reference_word_mismatch,name))
    end
    lc=LickorishCode(g,words,Dict{String,String}(),Dict{String,Vector{String}}(),
                              Dict{String,String}(),Dict{String,Int}(),nothing)
    validate_lc_input(lc);realization=realize_lc(lc).value["realization"]
    # Glue the two sides of every reference segment. Components must be exactly
    # the source faces, witnessed by every original crossing corner.
    parent=Dict(f.id=>f.id for f in realization.faces)
    for name in _lc_k_names(g), i in eachindex(words[name])
        j=mod1(i+1,length(words[name]))
        a=_lc_corner_at_occ(words[name],i,"out");b=_lc_corner_at_occ(words[name],j,"in")
        _uf_union!(parent,realization.corner_face[a],realization.corner_face[b])
    end
    face_of=Dict(d=>i for (i,face) in enumerate(sourcefaces) for d in face)
    root_source=Dict{String,Int}();cornerchecks=Any[];rotationchecks=Any[]
    for v in graph.vertices
        label=labels[v];lcrot=Symbol[]
        for (sg,way) in _LC_ROT_SEQ
            family=(sg==_lc_sign_char(-eps[v])) ? :alpha : :beta
            push!(lcrot,_port(family,Symbol(way)))
            out=outgoing[(v,_port(family,Symbol(way)))]
            # Source walks the face on its left; LC walks the opposite side.
            opposite=isodd(out) ? out+1 : out-1
            sf=face_of[opposite]
            root=_uf_find!(parent,realization.corner_face[(label,sg,way)])
            haskey(root_source,root) && root_source[root]!=sf &&
                throw(HCTriangulationError(:face_transport_conflict,"source faces merged at vertex $v"))
            root_source[root]=sf
            push!(cornerchecks,Dict("vertex"=>v,"sign"=>sg,"way"=>way,"source_face"=>sf))
        end
        _cyclic_equal(lcrot,graph.rotation[v]) || throw(HCTriangulationError(:rotation_transport_mismatch,string(v)))
        push!(rotationchecks,Dict("vertex"=>v,"preserved"=>true))
    end
    roots=Set(_uf_find!(parent,f.id) for f in realization.faces)
    Set(keys(root_source))==roots && Set(values(root_source))==Set(eachindex(sourcefaces)) &&
        length(root_source)==length(sourcefaces) ||
        throw(HCTriangulationError(:face_transport_not_bijective,"forgotten LC faces do not match source faces"))
    face_map=Dict(i=>[f.id for f in realization.faces if root_source[_uf_find!(parent,f.id)]==i] for i in eachindex(sourcefaces))
    joining_faces=face_map[Int(joining_region)]
    hc=resolve_hc(lc,roles;disks=Dict("joining"=>joining_faces))
    forgotten=forget_reference(realize_hc(hc))
    matched=[r for r in forgotten.regions if Set(r.face_roots)==Set(joining_faces)]
    length(matched)==1 || throw(HCTriangulationError(:joining_face_transport_mismatch,
        "transported face union does not identify one forgotten region"))
    joining_region_id=only(matched).id
    certificate=Dict{String,Any}("genus"=>g,"polygon_word"=>word,"polygon_slides"=>Any[],
        "standard_boundary_word"=>(g==1 ? [word[2],-word[1],-word[2],word[1]] : word),
        "boundary_lanes"=>lanes,"old_vertex_to_lc_label"=>labels,"old_edge_subdivision"=>subdivision,
        "vertex_rotation_checks"=>rotationchecks,"source_curves"=>sourcecurves,
        "source_face_to_lc_faces"=>face_map,"face_corner_checks"=>cornerchecks,
        "joining_region_source"=>joining_region,"joining_lc_faces"=>joining_faces,
        "joining_region_id"=>joining_region_id,
        "grid_bound"=>bound,"grid_trials"=>trials,"endpoint_parameters"=>string.(chosen.us),
        "coordinate_model"=>"convex polygon (i,i^2), doubled retained sides and cotree ports")
    if g==2
        certificate["overlay_grid_rule"]="x=1/(r+1), r=1:bound"
        certificate["even_handle_orientation_reversed"]=true
        certificate["reference_arcs"]=Dict(name=>[Dict("from_side"=>p.from.side,"from_u"=>string(p.from.u),
            "to_side"=>p.to.side,"to_u"=>string(p.to.u)) for p in pieces] for (name,pieces) in layout.overlay.reference_arcs)
    end
    (;lc,hc,joining_disk_name="joining",joining_faces,joining_region_id,certificate)
end

function hc_triangulation_to_resolved(input;orders=nothing,max_nodes::Int=100_000)
    flow=hc_flow_certificate(input;orders,max_nodes)
    hc_cut_to_resolved(flow.diagram;joining_region=flow.region)
end

export hc_cut_to_resolved, hc_triangulation_to_resolved
