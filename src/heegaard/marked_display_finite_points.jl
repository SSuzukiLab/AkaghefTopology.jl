# Source: marked_display_source_map.jl::_mdsm_map_point/_mdsm_collar_layer;
# marked_display_source_map.py::build_rendered_mesh, build_graph_free_collars,
# build_fixed_quotient; finite point composition extracted from the earlier
# FixedReferencePaths prototype. This maintained helper has NO scratch imports.
# Include inside Heegaard after marked_display_source_map.jl.
module MarkedDisplayFinitePoints
const H=parentmodule(@__MODULE__)
const M=H.MarkedSurfaceDiagrams
const Q=M.Rat
need(ok,message)=H._mm_need(ok,"finite display point: "*message)
point(x)=H._mdsm_point(x)
cell(x)=H._mm_cell(x)
xy(w)=(w[2],w[3])
sub(a,b)=(a[1]-b[1],a[2]-b[2])
blend(a,b,t)=ntuple(k->a[k]+t*(b[k]-a[k]),2)
function mesh_atoms(P,R,T;label)
    [begin
        source=[P[Int(v)] for v in t];target=[R[Int(v)] for v in t]
        need(M.orient(source...)>0 && M.orient(target...)>0,"folded PL component")
        (;source,target,label=(label,i))
    end for (i,t) in enumerate(T)]
end

# Source: existing exact component mapper owns barycentric location/tie images.
# Local collar complement alone is identity; a missing global cell is an error.
function point_stage(p,stage)
    hits=Tuple{Q,Q}[]
    for atom in stage.atoms
        a,b,c=atom.source
        all(>=(0),(M.orient(p,b,c),M.orient(a,p,c),M.orient(a,b,p))) || continue
        q=H._mdsm_map_point(p,Dict(enumerate(atom.source)),Dict(enumerate(atom.target)),[[1,2,3]])
        push!(hits,q.point)
    end
    if isempty(hits)
        need(stage.identity_elsewhere,"point leaves a complete finite PL domain")
        return p
    end
    need(all(==(first(hits)),hits),"component-edge images disagree")
    first(hits)
end

# Source: production build_rendered_mesh exactly, including all refined
# constraint vertices. These are actual SVG units; no second rounding is added.
function rendered_stages(recipe)
    pm=recipe["port_map"];U=H._mdsm_points(pm["source_positions"])
    V=H._mdsm_points(pm["target_positions"]);T=pm["triangles"]
    port=(;atoms=mesh_atoms(U,V,T;label="port"),identity_elsewhere=false)
    S=Dict(i=>ntuple(k->600*q[k],2) for (i,q) in V)
    R=copy(S);ids=pm["original_vertex_ids"]
    for i in 0:length(ids)-1
        R[i]=ntuple(k->H._mdsm_q(Float64(S[i][k]/600)*600),2)
    end
    idx=Dict(Int(v)=>i-1 for (i,v) in enumerate(ids));assigned=Dict()
    for (edge,marker) in pm["refined_constraint_edges"]
        a,b=Int.(pm["constraints"][Int(marker)]);u,v=idx[a],idx[b]
        axis=abs(S[v][1]-S[u][1])>=abs(S[v][2]-S[u][2]) ? 1 : 2
        for w in Int.(edge)
            t=(S[w][axis]-S[u][axis])/(S[v][axis]-S[u][axis])
            need(0<=t<=1 && blend(S[u],S[v],t)==S[w],"refined constraint left original edge")
            q=blend(R[u],R[v],t)
            need(!haskey(assigned,w) || assigned[w]==q,"rounded constraint ties disagree")
            assigned[w]=q
        end
    end
    merge!(R,assigned);A=copy(R)
    for (i,q) in recipe["exact_circle_anchors"]["corrected_boundary_positions"]
        A[parse(Int,string(i))]=ntuple(k->H._mdsm_q(q[k])*600,2)
    end
    # Scaling is itself an exact affine stage, so all components share one API.
    scaled=Dict(i=>ntuple(k->600*q[k],2) for (i,q) in V)
    scale=(;atoms=mesh_atoms(V,scaled,T;label="SVG scale"),identity_elsewhere=false)
    rounded=(;atoms=mesh_atoms(S,R,T;label="SVG rounding"),identity_elsewhere=false)
    anchored=(;atoms=mesh_atoms(R,A,T;label="exact anchor"),identity_elsewhere=false)
    [port,scale,rounded,anchored]
end

function compile_finite_base(recipe,document;rendered_svg=false)
    need(recipe["source_document_sha256"]==H._mm_hash(document),"wrong native source")
    chart=recipe["source_chart"];P=H._mdsm_points(chart["positions"])
    native=Dict{Tuple{String,Int},Vector{Any}}()
    for (label,weights,t) in zip(chart["reduction_source_triangle_labels"],
            chart["reduction_source_triangle_barycentric"],chart["reduction_source_triangles"])
        label===nothing && continue
        source=[xy(H._mdsm_q.(w)) for w in weights];target=[P[Int(v)] for v in t]
        need(M.orient(source...)>0 && M.orient(target...)>0,"native chart cell reverses")
        push!(get!(native,cell(label),Any[]),(;source,target,label=copy(label)))
    end
    layers=recipe["collars"]
    need(length(layers)==2 && layers[1]["layer"]=="native_to_clipped_source" &&
        layers[2]["layer"]=="final_mesh_polygon_to_circle_polygon","collar order changed")
    function collar_stage(layer)
        atoms=Any[]
        for a in H._mdsm_collar_layer(layer)
            append!(atoms,mesh_atoms(a.source,a.target,a.triangles;label=(a.cut,a.side)))
        end
        (;atoms,identity_elsewhere=true)
    end
    fixed=recipe["fixed_map"];S=H._mdsm_points(fixed["source"]["source_positions"])
    R=H._mdsm_points(fixed["target_positions"])
    middle=(;atoms=mesh_atoms(S,R,fixed["source"]["triangles"];label="fixed"),identity_elsewhere=false)
    stages=Any[collar_stage(layers[1]),middle,collar_stage(layers[2])]
    rendered_svg && append!(stages,rendered_stages(recipe))
    (;native,stages,recipe,rendered_svg)
end

function old_point(base,parent,p)
    atoms=get(base.native,cell(parent),Any[])
    stage=(;atoms,identity_elsewhere=false)
    if !any(all(>=(0),(M.orient(p,a.source[2],a.source[3]),M.orient(a.source[1],p,a.source[3]),M.orient(a.source[1],a.source[2],p))) for a in atoms)
        need(cell(parent)==cell(base.recipe["infinity_extension"]["source_triangle"]),"missing native chart point outside cap")
        return (;finite=false,old_parent=cell(parent),old_point=p)
    end
    q=point_stage(p,stage);pre_rendered_point=nothing
    for (i,stage) in enumerate(base.stages)
        q=point_stage(q,stage)
        i==3 && (pre_rendered_point=q)
    end
    (;finite=true,old_parent=cell(parent),old_point=p,point=q,pre_rendered_point)
end

# Source: build_fixed_quotient's F1=K1 composed with R1 composed with f1
# composed with B, and final_graph. Outside the FINAL graph-free outer polygon
# both R1 and K1 are identity. K0 is absent from the ordinary native endpoint.
# This does not evaluate points inside support or the native infinity cap.
function require_final_identity_attachment(base,q)
    p=base.rendered_svg ? ntuple(k->q[k]/600,2) : q
    for c in base.recipe["graph_free_collars"]["collars"]
        z=sub(p,point(c["center"]))
        need(any(sum(z[k]*H._mdsm_q(n[k]) for k in 1:2)>H._mdsm_q(h)
            for (n,h) in zip(c["halfspace_normals"],c["halfspace_offsets"])),
            "point inside final attachment support; explicit source replacement required")
    end
    q
end

# Source: old-reference C0=K0 composed with R0 composed with B additionally
# requires initial attachment evidence. Preserve this stricter existing gate
# for old-reference consumers; ordinary current-endpoint annotations use above.
function require_identity_attachment(base,q;initial_point=nothing)
    require_final_identity_attachment(base,q)
    p=base.rendered_svg ? ntuple(k->q[k]/600,2) : q
    need(!base.rendered_svg || initial_point!==nothing,"initial B point required for separate K0 support check")
    initial=initial_point===nothing ? p : initial_point
    for cut in base.recipe["fixed_quotient"]["cuts"],b in cut["boundaries"]
        z=sub(initial,point(b["center"]));r=H._mdsm_q(b["initial_support_radius"])
        need(sum(x*x for x in z)>r*r,"point inside initial quotient support")
    end
    q
end
end
