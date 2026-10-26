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
    tree=Set{Int}();visited=Set([1]);queue=[1]
    for a in queue, (e,b) in sort(adj[a])
        b in visited && continue
        push!(visited,b);push!(queue,b);push!(tree,e)
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
        "dual_eliminations"=>eliminations,
        "generator_edges"=>[Dict("generator"=>i,"edge"=>e,"occurrence"=>collect(reps[e])) for (i,e) in enumerate(generators)],
        "curve_word_convention"=>"saved orientation, based via the primal spanning tree; closed curve start is the saved word start")
    p=GroupPresentation(["x$i" for i in eachindex(generators)],relations,Int[],meta)
    (;presentation=p,curve_words=words)
end

surface_group_presentation(doc::AbstractDict)=_marked_surface_words(_marked_states(doc).geometry).presentation
marked_curve_words(doc::AbstractDict)=_marked_surface_words(_marked_states(doc).geometry)

function fundamental_group_presentation(doc::AbstractDict)
    states=_marked_states(doc);st=states.mathematics
    _marked_compression(st).closed || throw(ArgumentError("The marked systems do not compress to spherical boundaries"))
    data=_marked_surface_words(st);p=data.presentation
    curves=sort(collect(keys(st.curves)))
    relations=vcat(p.relators,[data.curve_words[c] for c in curves])
    meta=merge(p.metadata,Dict("space"=>"closed Heegaard 3-manifold",
        "attached_curve_relations"=>curves,"source_hash"=>MarkedSurfaceDiagrams.representation_digest(states.document)))
    GroupPresentation(p.generators,relations,Int[],meta)
end

function intersection_matrix(doc::AbstractDict)
    st=_marked_states(doc).mathematics
    a=sort([c for (c,(d,_)) in st.curves if d["role"]=="alpha"])
    b=sort([c for (c,(d,_)) in st.curves if d["role"]=="beta"])
    ai=Dict(c=>i for (i,c) in enumerate(a));bi=Dict(c=>i for (i,c) in enumerate(b))
    matrix=zeros(Int,length(a),length(b))
    for x in values(st.crossing_at)
        c,d=x["curves"]
        if haskey(ai,c);matrix[ai[c],bi[d]]+=x["sign"]
        else;matrix[ai[d],bi[c]]-=x["sign"];end
    end
    matrix
end
