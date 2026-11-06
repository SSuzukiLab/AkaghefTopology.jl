"""Free reduction of a signed-edge word; no commutativity is assumed."""
function _marked_reduce(word)
    out=Int[]
    for x in word
        !isempty(out) && out[end]==-x ? pop!(out) : push!(out,x)
    end
    out
end
_marked_inverse(word)=-reverse(word)

"""Cellular surface presentation and saved curve words from the SAME embedding.

Reuse MarkedSurfaceDiagrams' fine edge involution and vertex quotient. Collapse a
primal spanning tree and eliminate a disjoint dual spanning tree by triangle
relations. Each elimination solves a generator occurring once, so it is an
explicit presentation transformation, not a homology-only approximation. For a
closed surface one face relation remains; with boundary the exterior is the dual
root and no face relation remains. This is an implementation-derived cellular
construction; it does not certify any source's mapping-class-group relations.
"""
function _marked_surface_words(st)
    M=MarkedSurfaceDiagrams;M.check_state(st)
    faces=sort(collect(keys(st.triangles)));fi=Dict(f=>i for (i,f) in enumerate(faces))
    occurrences=[(f,i,u,v) for (f,i) in faces for (u,v) in M.tedges(st.triangles[(f,i)])]
    pairs=M.fine_edge_occurrences(st)
    representative(o)=min(o,get(pairs,o,o))
    reps=sort(unique(representative(o) for o in occurrences));ei=Dict(o=>i for (i,o) in enumerate(reps))
    signed(o)=(o==representative(o) ? 1 : -1)*ei[representative(o)]
    vertices=sort(unique(M.root!(st.vertex_dsu,(f,v)) for (f,k) in faces for v in st.triangles[(f,k)]))
    vi=Dict(v=>i for (i,v) in enumerate(vertices))
    ends=[(vi[M.root!(st.vertex_dsu,(o[1],o[3]))],vi[M.root!(st.vertex_dsu,(o[1],o[4]))]) for o in reps]
    adj=[Tuple{Int,Int}[] for _ in vertices]
    for (e,(a,b)) in enumerate(ends)
        push!(adj[a],(e,b));push!(adj[b],(e,a))
    end
    tree=Set{Int}();visited=Set([1]);queue=[1];tree_paths=Dict(1=>Int[])
    for a in queue, (e,b) in sort(adj[a])
        b in visited && continue
        push!(visited,b);push!(queue,b);push!(tree,e)
        tree_paths[b]=vcat(tree_paths[a],ends[e][1]==a ? e : -e)
    end
    length(visited)==length(vertices) || error("Disconnected surface edge complex")
    with_boundary=st.surface_report["boundary_components"]>0
    root=with_boundary ? 0 : 1
    dual=Dict(i=>Tuple{Int,Int}[] for i in (with_boundary ? (0:length(faces)) : (1:length(faces))))
    for (e,o) in enumerate(reps)
        e in tree && continue
        a=fi[o[1:2]];p=get(pairs,o,nothing);b=p===nothing ? 0 : fi[p[1:2]]
        a==b && continue
        push!(dual[a],(e,b));push!(dual[b],(e,a))
    end
    parent=Dict{Int,Tuple{Int,Int}}();seen=Set([root]);order=[root]
    for a in order, (e,b) in sort(dual[a])
        b in seen && continue
        push!(seen,b);push!(order,b);parent[b]=(a,e)
    end
    length(seen)==length(dual) || error("Tree complement has disconnected dual")
    dual_edges=Set(p[2] for p in values(parent))
    generators=sort(collect(setdiff(Set(eachindex(reps)),union(tree,dual_edges))))
    gi=Dict(e=>i for (i,e) in enumerate(generators))
    solutions=Dict{Int,Vector{Int}}();eliminations=Any[]
    function rewrite(word)
        out=Int[]
        for x in word
            abs(x) in tree && continue
            w=get(solutions,abs(x),[abs(x)])
            append!(out,x>0 ? w : _marked_inverse(w))
        end
        _marked_reduce(out)
    end
    faceword(k)=[signed((faces[k][1],faces[k][2],u,v)) for (u,v) in M.tedges(st.triangles[faces[k]])]
    for k in reverse(order[2:end])
        e=parent[k][2];w=rewrite(faceword(k));where=findall(x->abs(x)==e,w)
        length(where)==1 || error("Dual leaf relation does not contain its parent edge exactly once")
        j=only(where);rest=_marked_reduce(vcat(w[j+1:end],w[1:j-1]))
        solutions[e]=w[j]>0 ? _marked_inverse(rest) : rest
        all(x->haskey(gi,abs(x)),solutions[e]) || error("Uneliminated dual variable")
        push!(eliminations,Dict("source_triangle"=>[faces[k][1],faces[k][2]-1],
            "edge"=>e,"relation"=>w,"replacement"=>solutions[e]))
    end
    numbered(w)=[(gi[abs(x)],sign(x)) for x in rewrite(w)]
    relations=with_boundary ? Vector{Tuple{Int,Int}}[] : [numbered(faceword(root))]
    # Source: fine_edge_occurrences/tedges above and the SAME rewrite/numbered
    # elimination, with no second reducer. The image of u->v means t_u e t_v^-1,
    # where t_u is the actual saved primal-tree path from the base vertex. An
    # open edge remains a groupoid path; its image is a transported loop in this
    # presentation, not a claim that its two source endpoints coincide.
    edge_words=Dict{Tuple{String,Int,String,String},Vector{Tuple{Int,Int}}}()
    edge_correspondence=Any[]
    for o in occurrences
        w=numbered([signed(o)])
        for (u,v,image,direction) in ((o[3],o[4],w,1),
                (o[4],o[3],[(g,-s) for (g,s) in reverse(w)],-1))
            key=(o[1],o[2]-1,u,v)
            edge_words[key]=image
            push!(edge_correspondence,Dict("occurrence"=>collect(key),"word"=>image,
                "kind"=>"edge_groupoid_path","cellular_edge"=>abs(signed(o)),
                "cellular_direction"=>direction*sign(signed(o)),
                "start_vertex"=>collect(M.root!(st.vertex_dsu,(o[1],u))),
                "end_vertex"=>collect(M.root!(st.vertex_dsu,(o[1],v)))))
        end
    end
    # Every old triangle bounds a disk, hence its boundary is identity IN THE
    # SURFACE PRESENTATION. On a closed surface the root boundary is retained
    # as a relator, not generally the empty free word. Consumers retracting new
    # overlay edges through old cells must keep this distinction.
    face_words=Dict((f,i-1)=>numbered(faceword(k)) for (k,(f,i)) in enumerate(faces))
    face_correspondence=[Dict("source_triangle"=>[f,i-1],"word"=>face_words[(f,i-1)],
        "retained_relator"=>(!with_boundary && k==root ? 1 : nothing),
        "identity_scope"=>"identity modulo surface relators, not a canonical free-word equality")
        for (k,(f,i)) in enumerate(faces)]
    words=Dict{String,Vector{Tuple{Int,Int}}}()
    for (cid,(_,darts)) in st.curves
        word=Int[]
        for d in darts
            i,u,v=first(st.incidence[d.face][M.edge(d.u,d.v)])
            push!(word,(d.u==u && d.v==v ? 1 : -1)*signed((d.face,i,u,v)))
        end
        words[cid]=numbered(word)
    end
    expected=2st.surface_report["genus"]+max(st.surface_report["boundary_components"]-1,0)
    length(generators)==expected || error("Surface generator count disagrees with Euler characteristic")
    meta=Dict{String,Any}("space"=>"marked surface","source_hash"=>st.fingerprint,
        "base_vertex"=>collect(vertices[1]),"primal_tree_edges"=>sort(collect(tree)),
        "vertex_tree_paths"=>[Dict("vertex"=>collect(v),"signed_cellular_edges"=>tree_paths[i])
            for (i,v) in enumerate(vertices)],
        "dual_eliminations"=>eliminations,
        "cellular_edge_sources"=>[Dict("edge"=>e,"source_triangle"=>[o[1],o[2]-1],
            "oriented_vertices"=>[o[3],o[4]],
            "start_vertex"=>collect(vertices[ends[e][1]]),"end_vertex"=>collect(vertices[ends[e][2]]))
            for (e,o) in enumerate(reps)],
        "generator_edges"=>[Dict("generator"=>i,"edge"=>e,"occurrence"=>collect(reps[e]),
            "source_triangle"=>[reps[e][1],reps[e][2]-1],
            "oriented_vertices"=>[reps[e][3],reps[e][4]]) for (i,e) in enumerate(generators)],
        "relator_sources"=>(with_boundary ? Any[] : [Dict("relator"=>1,"kind"=>"surface_triangle",
            "source_triangle"=>[faces[root][1],faces[root][2]-1])]),
        "index_convention"=>"generators, relators and cellular edges are one-based; source_triangle is zero-based; legacy occurrence triangle is one-based",
        "edge_word_convention"=>"ordered source-generator basis; occurrence keys (face,zero-based triangle,u,v); word represents t_u edge t_v^-1 via vertex_tree_paths; source endpoints retained as groupoid data",
        "face_word_convention"=>"triangle boundary identity modulo the saved surface presentation; the one closed-surface root boundary is a retained relator and may be nonempty as a free word",
        "curve_word_convention"=>"saved orientation; spanning-tree transport to the base vertex; only closed curves are based loops, proper paths retain distinct endpoints")
    p=GroupPresentation(["x$i" for i in eachindex(generators)],relations,Int[],meta)
    (;presentation=p,curve_words=words,edge_words,edge_correspondence,face_words,face_correspondence)
end

# Source: schema.jl::source_manifest and vendor/marked_surface_diagrams/src/
# MarkedSurfaceDiagrams.jl::endpoint_anchors/representation_digest/marking_digest.
# The working geometry omits region claims, but every public result is bound to
# the FULL saved document (including auxiliary curves and joining-disk claims).
function _marked_bound_words(states)
    M=MarkedSurfaceDiagrams;st=states.geometry
    data=_marked_surface_words(st);p=data.presentation
    correspondence=Dict{String,Any}()
    for cid in sort(collect(keys(st.curves)))
        c,ds=st.curves[cid]
        correspondence[cid]=Dict("curve_id"=>cid,"role"=>c["role"],"closed"=>c["closed"],
            "kind"=>(c["closed"] ? "based_loop" : "unbased_proper_path"),
            "source_word"=>c["word"],"word"=>data.curve_words[cid],
            "endpoint_anchors"=>M.endpoint_anchors(st,cid),
            "start_vertex"=>[first(ds).face,first(ds).u],
            "end_vertex"=>[last(ds).face,last(ds).v],
            "word_semantics"=>(c["closed"] ? "saved loop transported to base_vertex via primal tree" :
                "path coordinates after spanning-tree collapse; not a based loop of the original path"))
    end
    meta=merge(p.metadata,Dict("working_cellulation_sha256"=>st.fingerprint,
        "source_hash"=>M.representation_digest(states.document),
        "source_representation_sha256"=>M.representation_digest(states.document),
        "source_marking_sha256"=>M.marking_digest(states.document),
        "source_manifest"=>source_manifest(states.document),"curve_correspondence"=>correspondence))
    (;presentation=GroupPresentation(p.generators,p.relators,p.beta_ids,meta),curve_words=data.curve_words,
      edge_words=data.edge_words,edge_correspondence=data.edge_correspondence,
      face_words=data.face_words,face_correspondence=data.face_correspondence)
end

surface_group_presentation(doc::AbstractDict)=_marked_bound_words(_marked_states(doc)).presentation
"""Saved curve and oriented fine-edge images in one fixed source surface basis.

`edge_words[(face, triangle, u, v)]` uses a ZERO-based triangle index, and
contains both directions of every fine-edge occurrence. Words use the returned
presentation's one-based generator indices. `edge_correspondence` and
`face_correspondence` are JSON-compatible source records; `face_words` uses
`(face, triangle)` keys. Metadata supplies the base vertex, actual signed
primal-tree paths, source generator edges and original document/marking hashes.
These are transported edge-groupoid paths modulo the saved surface relators,
not canonical free words. A consumer of a changed mesh must explicitly map its
new edge chains to these OLD source cells instead of recomputing a new basis.
"""
marked_curve_words(doc::AbstractDict)=_marked_bound_words(_marked_states(doc))

# Source: the cellular surface construction above and marked_topology.jl::
# _marked_compression. Attach precisely the mathematical alpha/beta curves;
# auxiliary marking cuts never become manifold relations.
function _marked_manifold_presentation(states,data)
    st=states.mathematics
    _marked_compression(st).closed || throw(ArgumentError("The marked systems do not compress to spherical boundaries"))
    p=data.presentation;curves=sort(collect(keys(st.curves)))
    relations=vcat(p.relators,[data.curve_words[c] for c in curves])
    origins=vcat(p.metadata["relator_sources"],[Dict("relator"=>length(p.relators)+i,
        "kind"=>"attached_curve","curve_id"=>c,"role"=>st.curves[c][1]["role"])
        for (i,c) in enumerate(curves)])
    meta=merge(p.metadata,Dict("space"=>"closed Heegaard 3-manifold",
        "attached_curve_relations"=>curves,"relator_sources"=>origins))
    GroupPresentation(p.generators,relations,Int[],meta)
end

function fundamental_group_presentation(doc::AbstractDict)
    states=_marked_states(doc)
    _marked_manifold_presentation(states,_marked_bound_words(states))
end

# Source: groups.jl::intersection_matrix and the vendor legacy_projection sign
# convention. Preserve the bare matrix API; expose its ordering and each saved
# crossing contribution through the public report below.
function _marked_intersections(st)
    a=sort([c for (c,(d,_)) in st.curves if d["role"]=="alpha"])
    b=sort([c for (c,(d,_)) in st.curves if d["role"]=="beta"])
    ai=Dict(c=>i for (i,c) in enumerate(a));bi=Dict(c=>i for (i,c) in enumerate(b))
    matrix=zeros(Int,length(a),length(b))
    contributions=Any[]
    for x in sort(collect(values(st.crossing_at));by=x->x["id"])
        c,d=x["curves"]
        alpha,beta,sign=haskey(ai,c) ? (c,d,x["sign"]) : (d,c,-x["sign"])
        matrix[ai[alpha],bi[beta]]+=sign
        push!(contributions,Dict("crossing_id"=>x["id"],"alpha"=>alpha,"beta"=>beta,
            "row"=>ai[alpha],"column"=>bi[beta],"contribution"=>sign,
            "source_crossing"=>deepcopy(x)))
    end
    Dict("matrix"=>matrix,"matrix_shape"=>collect(size(matrix)),
        "matrix_rows"=>[collect(matrix[i,:]) for i in axes(matrix,1)],
        "row_curves"=>a,"column_curves"=>b,"contributions"=>contributions,
        "convention"=>"one-based rows alpha, columns beta; signed transverse crossings in saved orientation; proper paths retain fixed endpoints, not an absolute H1 pairing")
end
intersection_matrix(doc::AbstractDict)=_marked_intersections(_marked_states(doc).mathematics)["matrix"]

# Source: groups.jl::first_homology/abelianization_matrix. Smith reduction returns
# invariant factors only, not change-of-basis matrices; keep the correspondence
# at the cellular presentation level rather than inventing H1 basis coordinates.
function _marked_group_report(p)
    Dict("generators"=>p.generators,"relators"=>p.relators,"metadata"=>p.metadata)
end
function _marked_homology_report(p)
    h=first_homology(p)
    Dict("free_rank"=>h.free_rank,"torsion"=>h.torsion,"invariant_factors"=>h.invariant_factors,
        "exponent_matrix"=>h.exponent_matrix,"exponent_matrix_shape"=>collect(size(h.exponent_matrix)),
        "exponent_matrix_rows"=>[collect(h.exponent_matrix[i,:]) for i in axes(h.exponent_matrix,1)],
        "smith_normal_form"=>h.smith_normal_form,
        "row_relators"=>p.metadata["relator_sources"],"column_generators"=>p.metadata["generator_edges"],
        "element_coordinates_available"=>false,
        "scope"=>"H1 invariants and cellular presentation matrix; Smith change-of-basis and H1 element coordinates are not computed")
end

"""Computations and source-element correspondence for one native marked object.

Source: marked_topology.jl::_marked_states and marked_topology_report, and the
cellular/group functions in this file. Regions and intersections use alpha U
beta; surface words also retain auxiliary curves. Surface pi1/H1 apply even to
physical-boundary or incomplete diagrams. Closed-manifold pi1/H1 apply only to
closed alpha/beta families whose actual compressions are spherical. Proper
paths expose endpoint anchors and are never identified with original loops.
No isotopy, group-isomorphism decision or manifold recognition is asserted.
"""
function marked_computation_report(doc::AbstractDict)
    validity=validate(doc;level=:surface)
    if !validity.valid
        issues=[Dict("code"=>string(i.code),"message"=>i.message,"path"=>i.path) for i in validity.issues]
        return Dict("schema"=>"marked-computation-report/1","diagram_valid"=>false,"issues"=>issues,
            "surface"=>Dict("applicable"=>false,"reasons"=>["invalid_marked_diagram"]),
            "manifold"=>Dict("applicable"=>false,"reasons"=>["invalid_marked_diagram"]))
    end
    M=MarkedSurfaceDiagrams;states=_marked_states(doc);st=states.mathematics
    data=_marked_bound_words(states);p=data.presentation
    closed=all(c["closed"] for (c,_) in values(st.curves))
    compression=closed ? _marked_compression(st) : nothing
    reasons=!closed ? ["closed-manifold attachment requires closed alpha/beta curves, not proper paths"] :
        compression.closed ? String[] : ["actual compression retains positive genus or original physical boundary"]
    manifold=Dict{String,Any}("applicable"=>isempty(reasons),"reasons"=>reasons,"compression"=>compression)
    if isempty(reasons)
        q=_marked_manifold_presentation(states,data)
        manifold["presentation"]=_marked_group_report(q);manifold["first_homology"]=_marked_homology_report(q)
    end
    Dict("schema"=>"marked-computation-report/1","diagram_valid"=>true,
        "source"=>Dict("representation_sha256"=>M.representation_digest(states.document),
            "marking_sha256"=>M.marking_digest(states.document),"manifest"=>source_manifest(states.document),
            "curve_ids"=>sort(collect(keys(states.geometry.curves))),
            "mathematical_curve_ids"=>sort(collect(keys(st.curves)))),
        "surface"=>Dict("applicable"=>true,"reasons"=>String[],"topology"=>states.geometry.surface_report,
            "presentation"=>_marked_group_report(p),"first_homology"=>_marked_homology_report(p)),
        "curves"=>p.metadata["curve_correspondence"],"regions"=>M.regions(st),
        "oriented_edge_words"=>data.edge_correspondence,"face_boundary_words"=>data.face_correspondence,
        "region_convention"=>"complement of mathematical alpha U beta; auxiliary cuts are not barriers; source triangle indices are zero-based",
        "intersection_cycles"=>M.intersection_cycles(st),"intersections"=>_marked_intersections(st),
        "manifold"=>manifold)
end
