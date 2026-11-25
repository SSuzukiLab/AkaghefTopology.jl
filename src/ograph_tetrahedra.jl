"""
    disk_tetrahedra(v::VirtualLink) -> NamedTuple

Return the disk-labelled ordered tetrahedra dual to a certified closed normal o-graph.

# Arguments
- `v`: a nonempty unweighted closed normal presentation satisfying C1–C3
  in [the terminology specification](../specs/ograph-terminology.md).
  Checked: `ArgumentError` rejects open, weighted, Z2, empty or uncertifiable
  input, nonclosed disk circuits, invalid disk-node coverage, or incompatible
  face/edge orientations. See [`certify_closed_ograph`](@ref).

# Returns
A `NamedTuple` containing `tetrahedra` (one record per true vertex, with
`vertex`, `vertex_sign`, `sign`, `local_nodes` and `edges`), `disk_count`,
`node_to_disk` (the local-node to disk map), `closedness` (the replayed witness),
`ordered_edges`, `edge_orientations_verified`, `face_orientations_verified`,
and `face_gluings`. Each gluing records `tet`, `omitted`, `to_tet`,
`to_omitted`, `permutation`, `ograph_edge`, `ograph_strand`, `source_crossing`
and `target_crossing`; their coordinates are fixed by the linked face-map spec.

# Guarantees
1. The result describes the ordered dual tetrahedral data of the presentation under the closedness and duality interpretation of the linked specification. [proved/agent: by construction]
2. Each local node belongs to exactly one returned disk, and paired local edges have equal disk IDs and compatible directions. [proved/agent: by construction]
3. Both orientation verification flags are `true` on successful return. [proved/agent: by construction]
4. `v` is not modified. [proved/agent: by construction]

# Conventions
Use [face maps and orientations](../specs/heegaard-snappy.md#面写像と向き)
for tetrahedral numbering, local ordering, orientation signs and permutations.

# Not guaranteed
The operational C1–C3 witness is not a general proof of the reconstruction
theorem. No state sum, manifold recognition, or identification with the
separate coloured slot model is returned. Inconsistent internal recipes can
raise `ErrorException` instead of a domain error.

# Sources
cf:docs/specs/heegaard-snappy.md — domain, dual closed triangulation, face maps and orientations.
cf:docs/specs/ograph-move-verifier.md — reviewed ordered tetrahedral face recipe and its source comparison.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[1,-1]], orientation=[1]);

julia> t = disk_tetrahedra(v);

julia> (only(t.tetrahedra).local_nodes, only(t.tetrahedra).sign, t.disk_count)
((1, 5, 6, 3, 2, 4), -1, 2)
```

Role: convert
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
    face_gluings = NamedTuple[]
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
        # cf: docs/specs/heegaard-snappy.md — expose the checked simplicial map.
        permutation = ntuple(i -> target_rho[sigma[source_inverse[i]+1]+1], 4)
        push!(face_gluings, (tet=abs(tail)-1,
            omitted=source_rho[source_face+1], to_tet=abs(head)-1,
            to_omitted=target_rho[target_face+1], permutation=permutation,
            ograph_edge=edge.id, ograph_strand=edge.strand,
            source_crossing=tail, target_crossing=head))
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
      face_orientations_verified=true,face_gluings)
end
