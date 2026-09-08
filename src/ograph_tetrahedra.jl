"""
    disk_tetrahedra(v::VirtualLink)

Return the ordered tetrahedra dual to a nonempty, operationally certified
closed normal o-graph. `tetrahedra[i].edges` lists disk IDs in the order
`(01,02,03,12,13,23)`, oriented from the smaller local vertex to the larger.
`sign` is the tetrahedron orientation, equal to minus the repository vertex
sign under the reviewed GKS Figure 2.2 convention.

The six-node dictionary is shared with `disk_table`; the face-to-ordered-vertex
map is shared with the existing CP recipe audit. Order-preserving face gluing
preserves all these edge directions, so no inverse-edge flag is needed.
This is combinatorial input data, not a state sum or a proof of MST equivalence.
The zero-vertex presentation is rejected because it has no tetrahedral data.
"""
function disk_tetrahedra(v::VirtualLink)
    _require_closed_ograph_input(v, "disk_tetrahedra")
    _,vertex_signs=real_gauss_code(v)
    n=length(vertex_signs)
    n>0 || throw(ArgumentError("zero-vertex o-graphs have no disk-local tetrahedra"))
    certified=certify_closed_ograph(v)
    disk=disk_table(v)
    all(disk.isCycle) || throw(ArgumentError("disk circuits must be closed"))
    node_to_disk=zeros(Int,6n)
    for (id,path) in zip(disk.ID,disk.paths), node in path
        1<=node<=6n && node_to_disk[node]==0 ||
            throw(ArgumentError("disk table does not partition local nodes"))
        node_to_disk[node]=id
    end
    all(>(0),node_to_disk) || throw(ArgumentError("disk table misses local nodes"))
    ordered_edges=((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
    local_nodes=Dict{Int,NTuple{6,Int}}()
    germ_edges_by_sign=Dict{Int,NTuple{6,Tuple{Int,Int}}}()
    for vertex_sign in (-1,1)
        faces=_cp_recipe_germs(vertex_sign)
        rho=_CP_GKS_FACE_OMITTED_VERTEX[vertex_sign]
        germ_edges=ntuple(germ -> Tuple(sort([rho[face+1] for face in 0:3
                                              if !(germ in faces[face+1])])),6)
        length(unique(germ_edges))==6 || error("invalid tetrahedral germ recipe")
        germ_edges_by_sign[vertex_sign]=germ_edges
        local_nodes[vertex_sign]=ntuple(i -> only(germ for germ in 1:6
                                                if germ_edges[germ]==ordered_edges[i]),6)
    end
    # Validate directions, not only unordered edge incidence. Each generator
    # of the disk equivalence relation preserves the local ascending direction;
    # therefore its transitive closure supplies a consistent global direction.
    for edge in real_edge_table(v)
        tail,head=edge.crossing
        s,t=vertex_signs[abs(tail)],vertex_signs[abs(head)]
        source_face=2+Int(tail>0); target_face=Int(head>0)
        sigma=_cp_recipe_face_map(s,source_face,t,target_face)
        sigma===nothing && error("face recipe has no unique simplicial map")
        source_rho=_CP_GKS_FACE_OMITTED_VERTEX[s]
        target_rho=_CP_GKS_FACE_OMITTED_VERTEX[t]
        (-s)*(-1)^source_rho[source_face+1] == -(-t)*(-1)^target_rho[target_face+1] ||
            throw(ArgumentError("paired tetrahedral faces have incompatible orientations"))
        source_inverse=[only(i-1 for i in 1:4 if source_rho[i]==a) for a in 0:3]
        for (a,b) in zip(_cp_recipe_germs(s)[source_face+1],
                         _cp_recipe_germs(t)[target_face+1])
            mapped=Tuple(target_rho[sigma[source_inverse[x+1]+1]+1]
                         for x in germ_edges_by_sign[s][a])
            mapped==germ_edges_by_sign[t][b] ||
                throw(ArgumentError("face gluing reverses a local ordered edge"))
            node_to_disk[6(abs(tail)-1)+a]==node_to_disk[6(abs(head)-1)+b] ||
                error("face pairing disagrees with disk partition")
        end
    end
    tetrahedra=[(vertex=i,vertex_sign=vertex_signs[i],sign=-vertex_signs[i],
                 local_nodes=local_nodes[vertex_signs[i]],
                 edges=ntuple(j -> node_to_disk[6(i-1)+local_nodes[vertex_signs[i]][j]],6))
                for i in 1:n]
    (;tetrahedra,disk_count=length(disk.ID),node_to_disk,
      closedness=certified.witness,ordered_edges,edge_orientations_verified=true,
      face_orientations_verified=true)
end
