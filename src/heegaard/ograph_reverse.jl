"""
Independently construct a Heegaard surface from a normal o-graph's dual
tetrahedra. Source metadata is never read. This is a new splitting of the
oriented manifold, not recovery of a source marking, D, or branching field.

The local tetrahedron recipe is reused from the research root's Topology.jl
commit 5a372b64f3a11519d9c7787b67eefcb4db6613df:
src/ograph_move_verifier.jl::_cp_recipe_germs/_cp_recipe_face_map and
_CP_GKS_FACE_OMITTED_VERTEX; src/ograph_tetrahedra.jl::disk_tetrahedra;
docs/specs/ograph-tetrahedra.md and ograph-move-verifier.md (GKS26 Fig.2.2).
Only the pure local recipe is imported here; newer database/move dependencies
are not required. The flag-level splitting argument below is agent-derived.
"""
const _REVERSE_RECIPE_SOURCE = "Topology.jl@5a372b64f3a11519d9c7787b67eefcb4db6613df/ograph_tetrahedra"
const _REVERSE_RHO = Dict(1=>(2,0,3,1), -1=>(3,1,2,0))
_reverse_perms() = [(a,b,c,d) for a in 0:3 for b in 0:3 for c in 0:3 for d in 0:3 if length(Set((a,b,c,d)))==4]
_reverse_parity(p) = iseven(sum(p[i]>p[j] for i in 1:4 for j in i+1:4)) ? 1 : -1

function _reverse_recipe_germs(s)
    # Exact production disk dictionary, shared by the existing local disk_table.
    dp=zeros(Int,12);dm=zeros(Int,12)
    dp[[1,2,4,5,9,12,7,6,8,11,10,3]]=repeat(1:6,2)
    dm[[1,2,4,5,9,12,7,10,3,11,6,8]]=repeat(1:6,2)
    x=s==1 ? dp : dm
    ntuple(f->Tuple(x[(3*f-2):(3*f)]),4)
end
function _reverse_face_map(s,f,t,g)
    left,right=_reverse_recipe_germs(s),_reverse_recipe_germs(t)
    edge(table,x)=Tuple(i for i in 0:3 if !(x in table[i+1]))
    candidates=[p for p in _reverse_perms() if p[f+1]==g &&
        all(Set(p[i+1] for i in edge(left,a))==Set(edge(right,b)) for (a,b) in zip(left[f+1],right[g+1]))]
    length(candidates)==1 || error("dual face recipe is not unique")
    raw=only(candidates);a,b=_REVERSE_RHO[s],_REVERSE_RHO[t]
    ordered=ntuple(i->b[raw[findfirst(==(i-1),a)]+1],4)
    # Same ordered-edge and opposite induced-face checks as disk_tetrahedra.
    source=setdiff(0:3,[a[f+1]]);target=setdiff(0:3,[b[g+1]])
    [ordered[x+1] for x in source]==target || error("dual face map reverses an ordered edge")
    (-s)*(-1)^a[f+1] == -(-t)*(-1)^b[g+1] || error("dual face orientation mismatch")
    ordered
end

function _reverse_tetrahedra(G::OData)
    T=parentmodule(@__MODULE__);M=MarkedSurfaceDiagrams
    v=T.VirtualLink();T.set_data!(v;rgauss=deepcopy(G.gauss),orientation=copy(G.orientation))
    closed,conditions=T.is_closed(v)
    closed || throw(ArgumentError("reverse construction requires closed normal o-data (C1,C2,C3=$conditions)"))
    code,signs=T.real_gauss_code(v);n=length(signs)
    n>0 || throw(ArgumentError("zero-vertex o-data has no dual tetrahedra"))
    pair=Dict();records=[]
    for edge in T.real_edge_table(v)
        a,b=edge.crossing;i,j=abs(a),abs(b);f,g=2+Int(a>0),Int(b>0)
        rho,sigma=_REVERSE_RHO[signs[i]],_REVERSE_RHO[signs[j]]
        omitted,target=rho[f+1],sigma[g+1];p=_reverse_face_map(signs[i],f,signs[j],g)
        haskey(pair,(i,omitted)) && error("tetrahedron face reused")
        haskey(pair,(j,target)) && error("tetrahedron face reused")
        inverse=ntuple(k->findfirst(==(k-1),p)-1,4)
        pair[(i,omitted)]=(j,target,p);pair[(j,target)]=(i,omitted,inverse)
        push!(records,Dict("ograph_edge"=>edge.id,"source"=>[i,omitted],"target"=>[j,target],"vertex_map"=>collect(p)))
    end
    length(pair)==4n || error("unpaired dual face")
    vertex=M.UnionFind([(i,a) for i in 1:n for a in 0:3])
    edges=M.UnionFind([(i,a,b) for i in 1:n for a in 0:3 for b in a+1:3])
    # Link triangle corners (tetrahedron, its vertex, second endpoint).
    link=M.UnionFind([(i,a,b) for i in 1:n for a in 0:3 for b in 0:3 if a!=b])
    for ((i,f),(j,g,p)) in pair, a in setdiff(0:3,[f])
        M.unite!(vertex,(i,a),(j,p[a+1]))
        for b in setdiff(0:3,[f,a])
            x,y=sort([a,b]);u,w=sort([p[a+1],p[b+1]])
            M.unite!(edges,(i,x,y),(j,u,w));M.unite!(link,(i,a,b),(j,p[a+1],p[b+1]))
        end
    end
    # Explicit vertex-link 2-manifold test: paired triangle edges, reversed
    # induced orientations, connected triangles, and circular stars at EVERY
    # link vertex. This avoids replacing the sphere hypothesis by chi alone.
    triangleadj=Dict();corneradj=Dict();link_pairings=[]
    oriented_link(i,a)=(-signs[i])*(-1)^a>0 ? setdiff(0:3,[a]) : reverse(setdiff(0:3,[a]))
    oriented_side(q,f)=only((q[k],q[mod1(k+1,3)]) for k in 1:3 if q[k]!=f && q[mod1(k+1,3)]!=f)
    for ((i,f),(j,g,p)) in pair, a in setdiff(0:3,[f])
        b=p[a+1];push!(get!(triangleadj,(i,a),[]),(j,b))
        u,w=oriented_side(oriented_link(i,a),f);x,y=oriented_side(oriented_link(j,b),g)
        (p[u+1],p[w+1])==(y,x) || error("vertex-link edge orientations do not oppose")
        for c in setdiff(0:3,[f,a]);push!(get!(corneradj,(i,a,c),[]),(j,b,p[c+1]));end
        isless((i,f),(j,g)) && push!(link_pairings,Dict("source"=>[i,a,f],"target"=>[j,b,g],"oriented_edge"=>[u,w],"mapped_edge"=>[p[u+1],p[w+1]]))
    end
    stars=[]
    for r in sort(unique(M.root!(link,x) for x in keys(link.p)))
        cs=Set(x for x in keys(link.p) if M.root!(link,x)==r)
        all(length(get(corneradj,x,[]))==2 for x in cs) || error("vertex-link star is not closed")
        length(M.components(cs,corneradj))==1 || error("vertex-link star is disconnected")
        push!(stars,Dict("link_vertex"=>collect(r),"corners"=>[collect(x) for x in sort(collect(cs))],"circle"=>true))
    end
    roots=sort(unique(M.root!(vertex,x) for x in keys(vertex.p)));links=[]
    for r in roots
        triangles=[x for x in keys(vertex.p) if M.root!(vertex,x)==r]
        all(length(get(triangleadj,x,[]))==3 for x in triangles) || error("vertex link has unpaired triangle edges")
        length(M.components(Set(triangles),triangleadj))==1 || error("tetrahedron vertex link is disconnected")
        points=Set(M.root!(link,(i,a,b)) for (i,a) in triangles for b in 0:3 if a!=b)
        chi=length(points)-3length(triangles)÷2+length(triangles)
        chi==2 || throw(ArgumentError("dual vertex link is not a sphere: chi=$chi"))
        push!(links,Dict("vertex"=>collect(r),"triangles"=>[collect(x) for x in sort(triangles)],"chi"=>chi,
            "connected"=>true,"closed_oriented_surface"=>true,"all_vertex_stars_circular"=>true))
    end
    eroots=sort(unique(M.root!(edges,x) for x in keys(edges.p)))
    # Tree edges specify the meridian disks omitted on each side.
    primal=M.UnionFind(roots);ptree=Set()
    for e in eroots
        i,a,b=e;x=M.root!(vertex,(i,a));y=M.root!(vertex,(i,b))
        if M.root!(primal,x)!=M.root!(primal,y);push!(ptree,e);M.unite!(primal,x,y);end
    end
    dual=M.UnionFind(collect(1:n));dtree=Set()
    for r in records
        i,f=r["source"];j,g=r["target"]
        if M.root!(dual,i)!=M.root!(dual,j)
            push!(dtree,min((i,f),(j,g)));M.unite!(dual,i,j)
        end
    end
    length(ptree)==length(roots)-1 && length(dtree)==n-1 || error("dual triangulation is disconnected")
    (signs=signs,pair=pair,records=records,vertex=vertex,edges=edges,links=links,link_pairings=link_pairings,link_stars=stars,
     primal_tree=ptree,dual_tree=dtree,genus=length(eroots)-length(roots)+1,conditions=conditions)
end

"""Flag quads are the level 1/2 boundary of the primal one-skeleton neighborhood.

On each ordered barycentric flag (v,e,f,t), low vertices are v,e and high
vertices f,t. Quad order (vf,vt,et,ef) points low→high for a positive flag.
Actual face pairings, rather than Gauss-word drawing, glue the quads. Alpha
uses et–ef; beta uses ef–vf. Every local flag and side remains in the witness.
"""
function _reverse_flag_overlay(data)
    M=MarkedSurfaceDiagrams;n=length(data.signs);perms=_reverse_perms()
    flags=[(i,p) for i in 1:n for p in perms];index=Dict(x=>j for (j,x) in enumerate(flags))
    corners=M.UnionFind([(j,k) for j in eachindex(flags) for k in 1:4])
    sides=M.UnionFind([(j,k) for j in eachindex(flags) for k in 1:4]);adj=Dict()
    point(p,k)= k==1 ? ((p[1],),Tuple(sort(collect(p[1:3])))) :
        k==2 ? ((p[1],),(0,1,2,3)) : k==3 ? (Tuple(sort(collect(p[1:2]))),(0,1,2,3)) :
        (Tuple(sort(collect(p[1:2]))),Tuple(sort(collect(p[1:3]))))
    for (j,(i,p)) in enumerate(flags),s in 1:4
        q=collect(p);target=i;map=Tuple(0:3)
        if s==1;q[2],q[3]=q[3],q[2]
        elseif s==2;q[3],q[4]=q[4],q[3]
        elseif s==3;q[1],q[2]=q[2],q[1]
        else
            target,_,map=data.pair[(i,p[4])];q=[map[x+1] for x in p]
        end
        h=index[(target,Tuple(q))];adj[(j,s)]=(h,s);M.unite!(sides,(j,s),(h,s))
        for k in (s,mod1(s+1,4))
            a,b=point(p,k);mapped=(Tuple(sort([map[x+1] for x in a])),Tuple(sort([map[x+1] for x in b])))
            l=only(x for x in (s,mod1(s+1,4)) if point(q,x)==mapped)
            M.unite!(corners,(j,k),(h,l))
        end
    end
    all(get(adj,b,nothing)==a for (a,b) in adj) || error("flag gluing is not involutive")
    # Build oriented vertex links and retain the local quadrilateral identities.
    links=Dict();edgeends=Dict();kinds=Dict();quads=[];sectors=Dict()
    for (j,(i,p)) in enumerate(flags)
        order=(-data.signs[i])*_reverse_parity(p)>0 ? [1,2,3,4] : [4,3,2,1]
        for (offset,k) in enumerate(order)
            before=order[mod1(offset-1,4)];after=order[mod1(offset+1,4)]
            side(a,b)=mod1(a+1,4)==b ? a : b
            incoming=M.root!(sides,(j,side(before,k)));outgoing=M.root!(sides,(j,side(k,after)))
            r=M.root!(corners,(j,k));link=get!(links,r,Dict())
            haskey(link,outgoing) && link[outgoing]!=incoming && error("nonmanifold quad vertex link")
            link[outgoing]=incoming
            k==4 && (sectors[(r,outgoing)]=(i,p[1]))
        end
        for s in 1:4
            e=M.root!(sides,(j,s));a=M.root!(corners,(j,s));b=M.root!(corners,(j,mod1(s+1,4)))
            edgeends[e]=(a,b);kinds[e]=s==3 ? "alpha" : s==4 ? "beta" : "none"
        end
        push!(quads,Dict("tetrahedron"=>i,"flag"=>collect(p),"orientation"=>order,
            "corners"=>[collect(M.root!(corners,(j,k))) for k in 1:4],
            "sides"=>[collect(M.root!(sides,(j,k))) for k in 1:4]))
    end
    # Recheck the actual induced orientations after corner identifications.
    for ((j,s),(h,_)) in adj
        i,p=flags[j];k,q=flags[h]
        a,b=(s,mod1(s+1,4));c,d=(s,mod1(s+1,4))
        (-data.signs[i])*_reverse_parity(p)<0 && ((a,b)=(b,a))
        (-data.signs[k])*_reverse_parity(q)<0 && ((c,d)=(d,c))
        M.root!(corners,(j,a))==M.root!(corners,(h,d)) &&
            M.root!(corners,(j,b))==M.root!(corners,(h,c)) || error("quad pairing does not reverse induced boundary orientation")
    end
    for (r,link) in links
        start=minimum(keys(link));seen=Set();x=start
        while !(x in seen);push!(seen,x);x=get(link,x,nothing);x===nothing && error("open vertex link");end
        x==start && length(seen)==length(link) || error("disconnected vertex link")
    end
    curves=Dict();visits=Dict();curve_sources=Dict();roles=Dict()
    for family in ("alpha","beta")
        edges=sort([e for e in keys(edgeends) if kinds[e]==family]);incident=Dict()
        for e in edges, r in edgeends[e];push!(get!(incident,r,[]),e);end
        all(length(es)==2 for es in values(incident)) || error("meridian is not an embedded circle")
        unseen=Set(edges);number=0
        while !isempty(unseen)
            number+=1;cid="$(family)_$number";firstedge=minimum(unseen)
            start=minimum(edgeends[firstedge]);x=start;e=firstedge;cycle=[];walk=[]
            while true
                e in unseen || error("meridian repeats before closing");delete!(unseen,e)
                a,b=edgeends[e];y=x==a ? b : a
                push!(walk,(x,e,y));push!(cycle,x)
                x=y;x==start && break;e=only(setdiff(incident[x],[e]))
            end
            crossing_nodes=[r for r in cycle if r[2]==4]
            isempty(crossing_nodes) && error("meridian lacks full-overlay crossings")
            curves[cid]=crossing_nodes
            for k in eachindex(walk)
                r,out,_=walk[k];r[2]==4 || continue
                incoming=walk[mod1(k-1,length(walk))][2];visits[(cid,r)]=(incoming,out)
            end
            j,_=first(crossing_nodes);i,p=flags[j]
            src=if family=="alpha"
                a,b=sort(collect(p[1:2]));M.root!(data.edges,(i,a,b))
            else
                h,f,_=data.pair[(i,p[4])];min((i,p[4]),(h,f))
            end
            omitted=src in (family=="alpha" ? data.primal_tree : data.dual_tree)
            roles[cid]=omitted ? "aux" : family
            curve_sources[cid]=Dict("kind"=>family,"cell"=>collect(src),"omitted_tree_disk"=>omitted,
                "quad_edges"=>[collect(z[2]) for z in walk])
        end
    end
    nodecurves=Dict()
    for (c,rs) in curves,r in rs;push!(get!(nodecurves,r,[]),c);end
    labels=Dict(r=>"x$j" for (j,r) in enumerate(sort(collect(keys(nodecurves)))))
    pairs=Dict();signs=Dict();crossing_sources=Dict();dart_sectors=Dict()
    for r in sort(collect(keys(nodecurves)))
        cs=nodecurves[r];length(cs)==2 || error("full overlay crossing multiplicity")
        a=only(c for c in cs if startswith(c,"alpha"));b=only(c for c in cs if startswith(c,"beta"))
        ai,ao=visits[(a,r)];bi,bo=visits[(b,r)];link=links[r]
        order=[];e=ao
        while isempty(order) || e!=ao;push!(order,e);e=link[e];end
        order==[ao,bo,ai,bi] || order==[ao,bi,ai,bo] || error("nontransverse flag meridians")
        name=labels[r];pairs[name]=(a,b);signs[name]=order==[ao,bo,ai,bi] ? 1 : -1
        for (cid,incoming,outgoing) in ((a,ai,ao),(b,bi,bo))
            dart_sectors[(cid,name,1)]=sectors[(r,outgoing)]
            dart_sectors[(cid,name,-1)]=sectors[(r,incoming)]
        end
        j,k=r;i,p=flags[j]
        crossing_sources[name]=Dict("tetrahedron"=>i,"edge"=>sort(collect(p[1:2])),"face"=>sort(collect(p[1:3])),"quad_corner"=>collect(r))
    end
    (cycles=Dict(c=>[labels[r] for r in rs] for (c,rs) in curves),pairs=pairs,signs=signs,
     roles=roles,curve_sources=curve_sources,crossing_sources=crossing_sources,dart_sectors=dart_sectors,quads=quads,
     quad_vertex_count=length(links),quad_edge_count=length(edgeends))
end

"""Reuse reference_surface.jl's vertex-disks/ribbon-gluings/capped-faces model.
The rotations now come from flag vertex links, not a chosen example. Subdivide
each gap into three edges so monogons/digons are genuine triangulated disks.
Tree meridians are omitted as mathematical curves; their chart scaffold stays.
"""
function _reverse_marked_overlay(overlay,genus;display_extension=nothing)
    rotations=Dict();sigma=Dict();twin=Dict()
    for v in sort(collect(keys(overlay.pairs)))
        a,b=overlay.pairs[v];s=overlay.signs[v]
        r=[(a,v,1),(b,v,s),(a,v,-1),(b,v,-s)];rotations[v]=r
        for j in 1:4;sigma[r[j]]=r[mod1(j+1,4)];end
    end
    for c in sort(collect(keys(overlay.cycles)))
        vs=overlay.cycles[c]
        for (x,y) in zip(vs,circshift(vs,-1));twin[(c,x,1)]=(c,y,-1);twin[(c,y,-1)]=(c,x,1);end
    end
    inverse=Dict(v=>k for (k,v) in sigma);unseen=Set(keys(twin));faces=[]
    while !isempty(unseen)
        start=minimum(unseen);z=start;face=[]
        while true
            z in unseen || error("overlay face repeats");delete!(unseen,z);push!(face,z)
            z=inverse[twin[z]];z==start && break
        end
        push!(faces,face)
    end
    # Every complementary disk must correspond to precisely one original
    # tetrahedron vertex. Counts or homology alone do not establish this bridge.
    face_sources=[];strips=[]
    for (j,face) in enumerate(faces)
        sources=unique(overlay.dart_sectors[d] for d in face if haskey(overlay.dart_sectors,d))
        if isempty(sources) && display_extension!==nothing
            alpha=[d for d in face if d[1] in display_extension["source_curves"] && d[3]==1]
            parallel=[d for d in face if d[1] in display_extension["cut_curves"] && d[3]==-1]
            length(face)==4 && length(alpha)==1 && length(parallel)==1 || error("parallel-copy complementary piece is not a rectangle")
            a=only(alpha);p=only(parallel)
            display_extension["parallel_of"][p[1]]==a[1] || error("parallel rectangle joins different meridians")
            push!(face_sources,(0,j));push!(strips,Dict("face"=>j,"source_alpha_dart"=>collect(a),"parallel_dart"=>collect(p),"boundary_darts"=>[collect(d) for d in face]))
        else
            push!(face_sources,only(sources))
        end
    end
    length(unique(face_sources))==length(faces) || error("overlay face-to-tetrahedron-vertex map is not injective")
    if display_extension!==nothing
        length(strips)==length(display_extension["crossing_copies"]) || error("parallel annulus rectangle coverage mismatch")
        display_extension["rectangles"]=strips
    end
    patches=[];gluings=[];ports=Dict();gaps=Dict()
    for v in sort(collect(keys(rotations)))
        sides=Dict();order=[];bd=[]
        for (j,d) in enumerate(rotations[v])
            i=j-1;a="a$i";g="g$i"
            sides[a]=["b$(2i)","p$i","b$(2i+1)"]
            sides[g]=["b$(2i+1)","q$(i)a","q$(i)b","b$(mod(2i+2,8))"]
            append!(order,[a,g]);append!(bd,sides[a][1:end-1]);append!(bd,sides[g][1:end-1])
            ports[d]=["V_$v",a,"p$i"];gaps[d]=["V_$v",g]
        end
        push!(patches,Dict("face"=>"V_$v","side_order"=>order,"sides"=>sides,
            "triangles"=>[["z",bd[j],bd[mod1(j+1,length(bd))]] for j in eachindex(bd)]))
    end
    for d in sort(collect(keys(twin)))
        e=twin[d];isless(d,e) && push!(gluings,[ports[d][1:2],ports[e][1:2]])
    end
    for (j,face) in enumerate(faces)
        f="F$j";n=length(face);bd=["t$k" for k in 1:3n];order=["e$k" for k in 1:n]
        push!(patches,Dict("face"=>f,"side_order"=>order,
            "sides"=>Dict("e$k"=>[bd[mod1(3k-2+h,3n)] for h in 0:3] for k in 1:n),
            "triangles"=>[["z",bd[k],bd[mod1(k+1,3n)]] for k in 1:3n]))
        for (k,d) in enumerate(face);push!(gluings,[[f,"e$k"],gaps[d]]);end
    end
    curves=[Dict("id"=>c,"role"=>overlay.roles[c],"closed"=>true,
        "word"=>join(["V_$v{"*ports[(c,v,-1)][3]*" z "*ports[(c,v,1)][3]*"}" for v in overlay.cycles[c]]," | "))
        for c in sort(collect(keys(overlay.cycles)))]
    doc=Dict{String,Any}("version"=>1,"surface"=>Dict("marking_id"=>"dual-tetrahedron-flag-splitting-v1",
        "genus"=>genus,"gluings"=>gluings,"boundary"=>[]),"patches"=>patches,"curves"=>curves,
        "crossings"=>[Dict("id"=>v,"at"=>["V_$v","z"],"curves"=>collect(overlay.pairs[v]),"sign"=>overlay.signs[v]) for v in sort(collect(keys(overlay.pairs)))])
    chartmap=Dict("crossing_charts"=>Dict(v=>["V_$v","z"] for v in keys(overlay.pairs)),
        "curve_chart_words"=>Dict(c["id"]=>c["word"] for c in curves),
        "meridian_ports"=>[Dict("dart"=>collect(d),"chart_side_gate"=>ports[d]) for d in sort(collect(keys(ports)))],
        "cap_charts"=>[Dict("chart"=>"F$j","tetrahedron_vertex"=>collect(face_sources[j]),"boundary_darts"=>[collect(d) for d in faces[j]]) for j in eachindex(faces)])
    doc,faces,face_sources,chartmap
end

"""Construct a native marked Heegaard diagram independently from o-data.

All original tetrahedral vertex links must be spheres. Both emitted meridian
systems are validated by actual cut/cap topology. The witness gives the dual
face maps, 24 oriented quads per tetrahedron, meridian cells and crossing cells.
The only metadata consumed is newly constructed here, never G.metadata.
"""
function _ograph_to_marked_heegaard(G::OData;display_cuts::Bool=false)
    data=_reverse_tetrahedra(G);overlay=_reverse_flag_overlay(data)
    doc,faces,face_sources,chartmap=_reverse_marked_overlay(overlay,data.genus)
    length(faces)==4length(data.signs) || error("flag overlay does not have one disk per tetrahedron vertex")
    overlay.quad_vertex_count-overlay.quad_edge_count+length(overlay.quads)==2-2data.genus || error("flag surface Euler mismatch")
    display_extension=nothing;infinity=nothing
    if display_cuts
        extended,display_extension=_reverse_parallel_display(overlay)
        doc,faces,face_sources,chartmap=_reverse_marked_overlay(extended,data.genus;display_extension=display_extension)
        Set(x for x in face_sources if x[1]!=0)==Set((i,a) for i in 1:length(data.signs) for a in 0:3) || error("parallel extension lost a source complementary disk")
        infinity=["F$(only(findall(==((1,0)),face_sources)))",0]
    end
    merges=_reverse_merge_charts!(doc)
    if infinity!==nothing
        for m in merges
            infinity[1]==m["absorbed_chart"] || continue
            infinity=[m["kept_chart"],infinity[2]+m["absorbed_triangle_index_offset"]]
        end
    end
    doc["metadata"]=Dict("construction"=>"dual tetrahedron primal-dual handle splitting",
        "source_recipe"=>_REVERSE_RECIPE_SOURCE,"source_marking_recovered"=>false,"joining_disk_recovered"=>false,
        "branching_field_recovered"=>false,"retained_source_used"=>false,
        "correspondence"=>Dict("tetrahedron_faces"=>data.records,"tetrahedron_vertex_links"=>data.links,
            "vertex_link_edge_pairings"=>data.link_pairings,"vertex_link_circular_stars"=>data.link_stars,
            "flag_quadrilaterals"=>overlay.quads,"curves"=>overlay.curve_sources,"crossings"=>overlay.crossing_sources,
            "primal_tree_edges"=>[collect(x) for x in sort(collect(data.primal_tree))],
            "dual_tree_faces"=>[collect(x) for x in sort(collect(data.dual_tree))],
            "overlay_face_darts"=>[[collect(d) for d in face] for face in faces],
            "overlay_face_tetrahedron_vertices"=>[collect(x) for x in face_sources],"initial_chart_correspondence"=>chartmap,
            "chart_vertex_prefix_policy"=>"initial chart local vertex v becomes chart_v before merge identifications",
            "chart_merges"=>merges),
        "scope"=>"new marked splitting of the oriented dual-triangulated manifold; no original marked diagram, D or branching recovery")
    if display_extension!==nothing
        doc["metadata"]["display_cuts"]=display_extension["cut_curves"]
        doc["metadata"]["display"]=Dict("cut_curves"=>copy(display_extension["cut_curves"]),"infinity_triangle"=>infinity,
            "infinity_choice"=>"new choice in the cap corresponding to tetrahedron 1 vertex 0; not recovered from a source",
            "auxiliary_construction"=>display_extension)
    end
    st=MarkedSurfaceDiagrams.validate(doc)
    result=_marked_compression(_marked_states(doc).mathematics)
    result.alpha_complete && result.beta_complete || error("constructed meridians fail actual complete-cut verification")
    if display_extension!==nothing
        cuts=Set(display_extension["cut_curves"])
        regions=MarkedSurfaceDiagrams.regions(st;barriers=Set(e for (e,cs) in st.used if !isdisjoint(Set(cs),cuts)))
        length(regions)==1 && only(regions)["genus"]==0 && only(regions)["boundary_components"]==2data.genus || error("parallel auxiliary cuts are not a complete display cut system")
    end
    doc
end

"""Explicit alpha-left parallel copies in the oriented meridian ribbon surface.

This is a local annular construction, not a general neighborhood operation on
arbitrary marked meshes. At a crossing of sign epsilon the beta orientation
meets the left parallel AFTER alpha for epsilon=+1, BEFORE for epsilon=-1.
The same rule applies when the beta meridian is an omitted-tree auxiliary.
Source alpha/beta rotations are unchanged. The chart builder checks the added
complementary pieces are exactly the rectangles of these annuli.
"""
function _reverse_parallel_display(overlay)
    cycles=deepcopy(overlay.cycles);pairs=deepcopy(overlay.pairs);signs=deepcopy(overlay.signs);roles=deepcopy(overlay.roles)
    selected=sort([c for c in keys(cycles) if roles[c]=="alpha"])
    copies=Dict();parallel_of=Dict();cut_curves=String[];markers=Dict()
    for a in selected
        c="display_"*a;parallel_of[c]=a;push!(cut_curves,c);roles[c]="aux"
        cycles[c]=["parallel_"*x for x in overlay.cycles[a]]
        for x in overlay.cycles[a]
            y="parallel_"*x;alpha,beta=overlay.pairs[x];alpha==a || error("parallel source family mismatch")
            pairs[y]=(c,beta);signs[y]=signs[x];copies[x]=y
        end
    end
    for (b,vs) in overlay.cycles
        startswith(b,"beta") || continue
        cycles[b]=reduce(vcat,[haskey(copies,x) ? (signs[x]==1 ? [x,copies[x]] : [copies[x],x]) : [x] for x in vs];init=String[])
    end
    # The right side of alpha and left side of its parallel belong to the old
    # cap disks; the other sides bound the new annular rectangles. These known
    # sectors identify every old complementary disk after adding the copies.
    for (x,(a,b)) in overlay.pairs
        markers[(a,x,-1)]=overlay.dart_sectors[(a,x,-1)]
        if haskey(copies,x)
            markers[("display_"*a,copies[x],1)]=overlay.dart_sectors[(a,x,1)]
        else
            markers[(a,x,1)]=overlay.dart_sectors[(a,x,1)]
        end
    end
    for (c,vs) in overlay.cycles
        filter(x->!startswith(x,"parallel_"),cycles[c])==vs || error("parallel extension changed the original meridian order")
    end
    extra=Dict{String,Any}("kind"=>"explicit alpha-left annular push-off","source_curves"=>selected,
        "cut_curves"=>cut_curves,"parallel_of"=>parallel_of,"crossing_copies"=>copies,
        "original_alpha_beta_preserved"=>true,"scope"=>"construction-time local ribbon annuli; not a generic neighborhood API")
    merge(overlay,(cycles=cycles,pairs=pairs,signs=signs,roles=roles,dart_sectors=markers)),extra
end

"""Merge distinct chart disks along one paired side, retaining every other seam.

This is an explicit disk-gluing operation, not curve rerouting. It only adapts
the exact overlay to the immutable kernel's 64-patch/256-vertex resource bounds.
Curves and crossing vertices use the same identification as the triangles.
`maximum_refined_vertices` optionally bounds the EXACT V+E+F vertex count of
one subsequent barycentric subdivision, both before and after each merge.
The optional `correspondence` dictionary receives original-to-final vertices,
zero-based triangles and side identities; removed marking sides map to nothing.
Region/infinity seeds follow those same triangle maps. Opaque metadata is not
rewritten. Default merge choices and the returned merge-record vector are kept.
"""
function _reverse_merge_charts!(doc;maximum_vertices::Int=224,maximum_patches::Int=64,
                                maximum_refined_vertices::Union{Nothing,Int}=nothing,
                                correspondence::Union{Nothing,AbstractDict}=nothing)
    M=MarkedSurfaceDiagrams;patches=Dict(p["face"]=>p for p in doc["patches"])
    3<=maximum_vertices<=256 || throw(ArgumentError("merged chart vertex budget must lie in 3:256"))
    1<=maximum_patches<=64 || throw(ArgumentError("chart count budget must lie in 1:64"))
    maximum_refined_vertices===nothing || 3<=maximum_refined_vertices<=256 ||
        throw(ArgumentError("one-subdivision vertex budget must lie in 3:256"))
    # Reuse the immutable kernel's disk triangulation check, not a new disk
    # recognition heuristic. For two disks glued along a side with L vertices,
    # (V,E,F) becomes (Va+Vb-L,Ea+Eb-(L-1),Fa+Fb). Barycentric V'=V+E+F.
    counts=Dict()
    vertex_map=Dict();triangle_map=Dict();side_map=Dict()
    for (f,p) in patches
        _,inc,vs=M.patch_info(p)
        counts[f]=(length(vs),length(inc),length(p["triangles"]))
        maximum_refined_vertices===nothing || sum(counts[f])<=maximum_refined_vertices ||
            throw(ArgumentError("chart $f needs $(sum(counts[f])) vertices after one subdivision; budget is $maximum_refined_vertices"))
        for v in vs;vertex_map[(f,v)]=(f,f*"_"*v);end
        for i in eachindex(p["triangles"]);triangle_map[(f,i-1)]=(f,i-1);end
        for s in p["side_order"];side_map[(f,s)]=(f,s);end
    end
    darts=Dict(c["id"]=>M.parse_word(c["word"]) for c in doc["curves"])
    # Unique names before chart identifications, so unrelated local vertices
    # never merge accidentally when two disks both use a local name such as z.
    for (f,p) in patches
        for t in p["triangles"];t .= f.*"_".*t;end
        for path in values(p["sides"]);path .= f.*"_".*path;end
    end
    for (c,ds) in darts;darts[c]=[M.Dart(d.face,d.face*"_"*d.u,d.face*"_"*d.v) for d in ds];end
    for x in doc["crossings"];x["at"][2]=x["at"][1]*"_"*x["at"][2];end
    record=[];gluings=doc["surface"]["gluings"]
    while length(patches)>maximum_patches
        selected=nothing
        for (i,(a,b)) in enumerate(gluings)
            a[1]==b[1] && continue
            pa,pb=patches[a[1]],patches[b[1]]
            va,ea,ta=counts[a[1]];vb,eb,tb=counts[b[1]];shared=length(pa["sides"][a[2]])
            merged=(va+vb-shared,ea+eb-shared+1,ta+tb)
            merged[1]<=maximum_vertices || continue
            maximum_refined_vertices===nothing || sum(merged)<=maximum_refined_vertices || continue
            selected=i;break
        end
        selected===nothing && throw(ArgumentError("reverse marked output exceeds native chart resource bounds"))
        a,b=gluings[selected];fa,fb=a[1],b[1];pa,pb=patches[fa],patches[fb]
        triangle_offset=length(pa["triangles"])
        va,vb=pa["sides"][a[2]],pb["sides"][b[2]]
        length(va)==length(vb) || error("chart merger encountered unequal side subdivisions")
        M.side_parameters(pa,a[2])==[1-t for t in reverse(M.side_parameters(pb,b[2]))] ||
            error("chart merger encountered unequal oriented side parameters")
        replace=Dict(zip(vb,reverse(va)));rv(v)=get(replace,v,v)
        # The two cyclic boundary lists splice after deleting the glued arcs.
        oa,ob=pa["side_order"],pb["side_order"];ia=findfirst(==(a[2]),oa);ib=findfirst(==(b[2]),ob)
        tail(o,i)=[o[mod1(i+k,length(o))] for k in 1:length(o)-1]
        rename=Dict(s=>fb*"_"*s for s in ob if s!=b[2])
        pa["side_order"]=vcat(tail(oa,ia),[rename[s] for s in tail(ob,ib)])
        delete!(pa["sides"],a[2]);for (s,path) in pb["sides"];s==b[2] || (pa["sides"][rename[s]]=rv.(path));end
        # Preserve explicit nonuniform seam parameters when either disk has them.
        if haskey(pa,"side_parameters") || haskey(pb,"side_parameters")
            params=get!(pa,"side_parameters",Dict());delete!(params,a[2])
            for (s,values) in get(pb,"side_parameters",Dict());s==b[2] || (params[rename[s]]=deepcopy(values));end
        end
        append!(pa["triangles"],[[rv(v) for v in t] for t in pb["triangles"]])
        ca,cb=counts[fa],counts[fb];shared=length(va)
        counts[fa]=(ca[1]+cb[1]-shared,ca[2]+cb[2]-shared+1,ca[3]+cb[3]);delete!(counts,fb)
        deleteat!(gluings,selected)
        for pair in gluings,endpoint in pair
            endpoint[1]==fb || continue;endpoint[1]=fa;endpoint[2]=rename[endpoint[2]]
        end
        for (c,ds) in darts;darts[c]=[d.face==fb ? M.Dart(fa,rv(d.u),rv(d.v)) : d for d in ds];end
        for x in doc["crossings"]
            x["at"][1]==fb || continue;x["at"]=[fa,rv(x["at"][2])]
        end
        for (source,(f,v)) in vertex_map;f==fb && (vertex_map[source]=(fa,rv(v)));end
        for (source,(f,i)) in triangle_map;f==fb && (triangle_map[source]=(fa,i+triangle_offset));end
        for (source,target) in side_map
            target===nothing && continue
            if target==(fa,a[2]) || target==(fb,b[2]);side_map[source]=nothing
            elseif target[1]==fb;side_map[source]=(fa,rename[target[2]])
            end
        end
        push!(record,Dict("kept_chart"=>fa,"absorbed_chart"=>fb,"identified_vertices"=>replace,
            "absorbed_triangle_index_offset"=>triangle_offset,"absorbed_triangle_count"=>length(pb["triangles"]),
            "removed_sides"=>[a[2],b[2]],"absorbed_side_names"=>rename,
            "result_cell_counts"=>collect(counts[fa]),"result_one_subdivision_vertices"=>sum(counts[fa])));delete!(patches,fb)
    end
    doc["patches"]=[patches[f] for f in sort(collect(keys(patches)))]
    for c in doc["curves"];c["word"]=M.print_word(darts[c["id"]]);end
    for mark in get(doc,"region_marks",[]);mark["seed"]=collect(triangle_map[Tuple(mark["seed"])]);end
    display=get(get(doc,"metadata",Dict()),"display",Dict())
    if haskey(display,"infinity_triangle");display["infinity_triangle"]=collect(triangle_map[Tuple(display["infinity_triangle"])]);end
    if correspondence!==nothing
        correspondence["vertices"]=[Dict("source"=>collect(s),"target"=>collect(vertex_map[s])) for s in sort(collect(keys(vertex_map)))]
        correspondence["triangles"]=[Dict("source"=>collect(s),"target"=>collect(triangle_map[s])) for s in sort(collect(keys(triangle_map)))]
        correspondence["sides"]=[Dict("source"=>collect(s),"target"=>side_map[s]===nothing ? nothing : collect(side_map[s])) for s in sort(collect(keys(side_map)))]
        correspondence["cell_counts"]=[Dict("chart"=>f,"vertices"=>counts[f][1],"edges"=>counts[f][2],"triangles"=>counts[f][3],"one_subdivision_vertices"=>sum(counts[f])) for f in sort(collect(keys(counts)))]
        correspondence["vertex_prefix_policy"]="original local vertex v is prefixed chart_v before explicit seam identifications"
    end
    record
end
