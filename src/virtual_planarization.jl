mutable struct _FaceGraph
    node_count::Int
    edges::Vector{NTuple{3,Int}}
end
_FaceGraph() = _FaceGraph(0,NTuple{3,Int}[])

function _add_edges!(graph::_FaceGraph,source::Integer,targets,weights)
    target_values = targets isa Integer ? [targets] : collect(targets)
    weight_values = weights isa Integer ? [weights] : collect(weights)
    length(target_values) == length(weight_values) || throw(DimensionMismatch())
    graph.node_count=max(graph.node_count,source,maximum(target_values;init=0))
    append!(graph.edges,[(Int(source),Int(target),Int(weight))
                         for (target,weight) in zip(target_values,weight_values)])
    graph
end
_find_edges(graph::_FaceGraph,a::Integer,b::Integer) =
    findall(edge -> (edge[1] == a && edge[2] == b) || (edge[1] == b && edge[2] == a),graph.edges)
_neighbors(graph::_FaceGraph,node::Integer) = unique(vcat(
    [edge[2] for edge in graph.edges if edge[1] == node],
    [edge[1] for edge in graph.edges if edge[2] == node]))

function _remove_edges!(graph::_FaceGraph,indices)
    deleteat!(graph.edges,sort(unique(collect(indices)))); graph
end
function _remove_node!(graph::_FaceGraph,node::Integer)
    graph.edges=[(a-(a>node),b-(b>node),weight) for (a,b,weight) in graph.edges
                 if a != node && b != node]
    graph.node_count-=1
    graph
end
function _shortest_path(graph::_FaceGraph,source::Integer,target::Integer)
    source == target && return [Int(source)]
    queue=[Int(source)]; parent=Dict(Int(source)=>0)
    while !isempty(queue)
        node=popfirst!(queue)
        for neighbor in _neighbors(graph,node)
            haskey(parent,neighbor) && continue
            parent[neighbor]=node
            if neighbor == target
                path=Int[target]
                while last(path) != source
                    push!(path,parent[last(path)])
                end
                return reverse(path)
            end
            push!(queue,neighbor)
        end
    end
    Int[]
end

function _find_nested(key,arrays)
    locations=[(array_index,findfirst(==(key),array)) for (array_index,array) in enumerate(arrays)
               if key in array]
    length(locations)==1 || throw(ArgumentError("expected one occurrence of $key, found $(length(locations))"))
    locations[1]
end

"""Translate MATLAB `[a(1:k), value, a(k+1:end)]` without an index shift.

`prefix_length` may be zero. The original port passed this quantity directly
to Julia's `insert!`, whose argument is instead the one-based insertion
position.
"""
function _insert_after_prefix!(values::Vector{T},prefix_length::Integer,value::T;
                               context::AbstractString="") where T
    0 <= prefix_length <= length(values) || throw(ArgumentError(
        "invalid splice prefix $prefix_length for length $(length(values))"*
        (isempty(context) ? "" : " ($context)")))
    insert!(values,Int(prefix_length)+1,value)
end

function _debug_face_connections(graph::_FaceGraph,region_vertices,region_edges)
    [(faces=(source,target),edge=weight,incidences=[
        (face=face,vertices=(region_vertices[face][mod1(position-1,length(region_vertices[face]))],
                             region_vertices[face][position]))
        for face in eachindex(region_edges)
        for position in findall(==(abs(weight)),abs.(region_edges[face]))
    ]) for (source,target,weight) in graph.edges]
end

_debug_incidence_faces(weight,region_edges)=[
    face for face in eachindex(region_edges)
    for _ in findall(==(abs(weight)),abs.(region_edges[face]))
]

function _debug_planarization_checkpoint(graph::_FaceGraph,region_vertices,region_edges,
        strand_edges,gauss,orientation,source_gauss,source_orientation;
        phase, strand_index=0, vertex_position=0, current_face=0, io::IO=stderr)
    println(io,"[vg2g] phase=$phase strand=$strand_index visit=$vertex_position face=$current_face")
    println(io,"  region vertices = ",region_vertices)
    println(io,"  region edges    = ",region_edges)
    println(io,"  strand edges    = ",strand_edges)
    println(io,"  face graph edges= ",graph.edges)
    println(io,"  gauss           = ",gauss)
    println(io,"  orientation     = ",orientation)

    face_count=length(region_vertices)
    length(region_edges)==face_count || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: region vertex/edge face counts differ"))
    graph.node_count==face_count || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: graph has $(graph.node_count) nodes for $face_count faces"))
    all(length(region_vertices[i])==length(region_edges[i]) for i in 1:face_count) ||
        throw(ArgumentError("vg2g debug invariant failed at $phase: a face has mismatched vertex/edge boundary lengths"))
    all(1<=a<=face_count && 1<=b<=face_count for (a,b,_) in graph.edges) ||
        throw(ArgumentError("vg2g debug invariant failed at $phase: face graph endpoint is out of range"))
    println(io,"  face incidences = ",_debug_face_connections(graph,region_vertices,region_edges))

    graph_weights=Set(abs(weight) for (_,_,weight) in graph.edges if weight!=0)
    boundary_weights=Set(abs(edge) for boundary in region_edges for edge in boundary if edge!=0)
    graph_weights==boundary_weights || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: graph weights $graph_weights != boundary weights $boundary_weights"))
    graph_weight_counts=Dict(weight=>count(edge -> abs(edge[3])==weight,graph.edges)
                             for weight in graph_weights)
    boundary_weight_counts=Dict(weight=>count(==(weight),
        abs.(filter(!=(0),reduce(vcat,region_edges;init=Int[])))) for weight in boundary_weights)
    all(==(1),values(graph_weight_counts)) || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: a nonzero edge ID occurs more than once in the face graph"))
    all(==(2),values(boundary_weight_counts)) || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: a nonzero edge ID does not occur twice on region boundaries"))
    for (source,target,weight) in graph.edges
        weight==0 && continue
        incidence_faces=sort(_debug_incidence_faces(weight,region_edges))
        incidence_faces==sort([source,target]) || throw(ArgumentError(
            "vg2g debug invariant failed at $phase: edge $weight joins faces "*
            "$(sort([source,target])) in the graph but occurs on faces $incidence_faces"))
    end

    real_count=length(source_orientation)
    orientation[1:real_count]==source_orientation || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: real orientations changed"))
    all(==(0),orientation[real_count+1:end]) || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: an inserted crossing is not virtual"))
    preserved=[[entry for entry in component if entry!=0 && abs(entry)<=real_count]
               for component in gauss]
    expected=[[entry for entry in component if entry!=0] for component in source_gauss]
    preserved==expected || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: ordered real crossing visits changed"))
    labels=abs.(filter(!=(0),reduce(vcat,gauss;init=Int[])))
    all(label -> count(==(label),labels)==2,unique(labels)) || throw(ArgumentError(
        "vg2g debug invariant failed at $phase: a crossing does not occur exactly twice"))
    nothing
end

function _update_face!(graph::_FaceGraph,boundary,face1,face2,new_face,new_edge)
    # Update the face graph with the new boundary. The MATLAB source retains
    # `bnd=unique(bnd)` as a commented alternative; duplicates instead encode
    # a self-adjacent boundary and are kept explicitly here.
    counts=Dict{Int,Int}()
    for edge in boundary, counts[edge]=get(counts,edge,0)+1; end
    self_edges=sort([edge for (edge,count) in counts if count>1])
    ordinary=sort([edge for (edge,count) in counts if count==1])
    indices=[findfirst(item -> item[3] == edge,graph.edges) for edge in ordinary]
    any(isnothing,indices) && throw(ArgumentError("face boundary is absent"))
    other_faces=Int[]
    for index in Int.(indices)
        a,b,_=graph.edges[index]
        push!(other_faces,a == face1 ? b : a)
    end
    _remove_edges!(graph,Int.(indices))
    _add_edges!(graph,new_face,other_faces,ordinary)
    _add_edges!(graph,face2,face2,new_edge)
    if !isempty(self_edges)
        self_indices=[findfirst(item -> item[3] == edge,graph.edges) for edge in self_edges]
        _remove_edges!(graph,Int.(self_indices))
        _add_edges!(graph,new_face,fill(new_face,length(self_edges)),self_edges)
    end
    graph
end

"""Insert virtual crossings so a real-Gauss/o-data diagram becomes planar.

With `DEBUG_MODE=true`, emit the MATLAB-style region/edge/face state at stable
checkpoints and validate the face-graph and Gauss-code invariants. Plotting is
not a hot path, so this route favors diagnostic clarity over minimal overhead.
"""
function planarize_virtual_gauss(real_gauss,real_orientation;
                                 DEBUG_MODE::Bool=false,debug_io::IO=stderr)
    input=[Int.(component) for component in real_gauss]
    orientation=Int.(real_orientation)
    if isempty(orientation)
        all(isempty,input) || throw(ArgumentError(
            "empty orientation requires crossing-free Gauss components"))
        DEBUG_MODE && println(debug_io,"[vg2g] crossing-free input; planarisation is a pass-through")
        return deepcopy(input),orientation
    end

    # Store information about regions using vertex/edge/face arrays. The face
    # graph represents adjacency of regions and its weights are the edge IDs
    # connecting them. The object is retained deliberately as an independently
    # inspectable validation model when DEBUG_MODE is enabled.
    region_vertices=Vector{Vector{Int}}()
    region_edges=Vector{Vector{Int}}()
    gauss=deepcopy(input)
    strand_edges=[Int[] for _ in gauss]
    used_vertices=Set([0])
    edge_count=0
    vertex_count=maximum(abs.(reduce(vcat,input;init=Int[])))+1
    face_graph=_FaceGraph()

    for strand_index in eachindex(input)
        # Loop through each strand.
        strand=copy(input[strand_index])
        (isempty(strand) || strand == [0]) && continue
        cut=last(strand)==0
        strand=vcat(strand[1:end-cut],[0]) # Delete the cut suffix before sweeping.
        length(strand)==1 && (input[strand_index]=Int[]; continue)
        edge_count+=1
        push!(strand_edges[strand_index],edge_count)
        first_vertex=abs(strand[1])
        face1=0
        face2=0
        if first_vertex in used_vertices
            # Start strand sweeping from a used vertex.
            face2,vertex_index2=_find_nested(-first_vertex,region_vertices)
            region_vertices[face2]=vcat(vertex_count+1,region_vertices[face2][vcat(vertex_index2:length(region_vertices[face2]),1:vertex_index2)])
            region_edges[face2]=vcat(-edge_count,region_edges[face2][vertex_index2:end],region_edges[face2][1:vertex_index2-1],edge_count)
            _add_edges!(face_graph,face2,[face2,face2],[0,edge_count])
            face1,vertex_index=_find_nested(first_vertex,region_vertices)
            region_vertices[face1]=vcat(region_vertices[face1][1:vertex_index-1],-first_vertex,0,-first_vertex,region_vertices[face1][vertex_index+1:end])
            region_edges[face1]=vcat(region_edges[face1][1:vertex_index-1],[0,0],region_edges[face1][vertex_index:end])
            face1=face2
        else
            # Start strand sweeping from a new vertex.
            face1=length(region_vertices)+1
            sign_value=orientation[first_vertex]*strand[1]
            push!(region_vertices,[vertex_count+1,-sign_value,0,sign_value])
            push!(region_edges,[-edge_count,0,0,edge_count])
            _add_edges!(face_graph,face1,[face1,face1],[0,edge_count])
            push!(used_vertices,first_vertex)
        end

        for vertex_position in 2:length(strand)
            # Loop through each vertex in the strand. This checkpoint matches
            # MATLAB's `fd()` debug view before the next edge is attached.
            DEBUG_MODE && _debug_planarization_checkpoint(face_graph,region_vertices,region_edges,
                strand_edges,gauss,orientation,input,real_orientation;phase=:visit_entry,
                strand_index=strand_index,vertex_position=vertex_position,
                current_face=face1,io=debug_io)
            edge_count+=1
            push!(strand_edges[strand_index],edge_count)
            vertex=abs(strand[vertex_position])
            if vertex in used_vertices
                # Patterns 1, 2, 3, 4, and 6 meet an already-used vertex.
                face2,vertex_index2=_find_nested(vertex,region_vertices)
                new_face,_=_find_nested(-vertex,region_vertices)
                # Find the shortest path from face1 to face2 (and new_face) in
                # the graph object. It records which boundaries must be crossed.
                path=vcat(_shortest_path(face_graph,face1,face2),new_face)
                path_length=length(path)
                if path_length==1
                    # Pattern 1: connect to another component (f1⋯>f2->fnew).
                    new_face,vertex_index3=_find_nested(-vertex,region_vertices)
                    region_vertices[face2]=vcat(region_vertices[face2][1:vertex_index2-1],-vertex,
                        region_vertices[face1][2:end],-vertex,region_vertices[face2][vertex_index2+1:end])
                    region_vertices[new_face]=vcat(vertex_count+1,-vertex,
                        region_vertices[new_face][vcat(vertex_index3+1:end,1:vertex_index3-1)],-vertex)
                    deleteat!(region_vertices,face1)
                    region_edges[face2]=vcat(region_edges[face2][1:vertex_index2-1],region_edges[face1],region_edges[face2][vertex_index2:end])
                    region_edges[new_face]=vcat(-edge_count,region_edges[new_face][vcat(vertex_index3:end,1:vertex_index3-1)],edge_count)
                    deleteat!(region_edges,face1)
                    for neighbor in filter(!=(face1),_neighbors(face_graph,face1))
                        indices=_find_edges(face_graph,face1,neighbor)
                        _add_edges!(face_graph,face2,fill(neighbor,length(indices)),[face_graph.edges[i][3] for i in indices])
                    end
                    indices=_find_edges(face_graph,face1,face1)
                    _add_edges!(face_graph,face2,fill(face2,length(indices)),[face_graph.edges[i][3] for i in indices])
                    _add_edges!(face_graph,new_face,new_face,edge_count)
                    _remove_node!(face_graph,face1)
                    face1=new_face-(face1<new_face)
                elseif path[1]==path[2] && vertex_position<length(strand)
                    # Patterns 2 and 3: connect to a vertex in the same face.
                    vertex_index=vertex_index2
                    face2,vertex_index2=_find_nested(-vertex,region_vertices)
                    new_face=length(region_vertices)+1
                    if vertex_index2<vertex_index
                        # Pattern 2.
                        push!(region_vertices,region_vertices[face1][vcat(vertex_index2,vertex_index+1:length(region_vertices[face1]))])
                        region_vertices[face1]=region_vertices[face1][vcat(1,vertex_index2:vertex_index-1,vertex_index2,2:vertex_index2)]
                        push!(region_edges,region_edges[face1][vertex_index:end])
                        region_edges[face1]=vcat(-edge_count,region_edges[face1][vcat(vertex_index2:vertex_index-1,1:vertex_index2-1)],edge_count)
                    elseif vertex_index<vertex_index2
                        # Pattern 3.
                        push!(region_vertices,region_vertices[face1][vcat(vertex_index2,2:vertex_index-1)])
                        region_vertices[face1]=region_vertices[face1][vcat(1,vertex_index2:length(region_vertices[face1]),vertex_index2,vertex_index+1:vertex_index2)]
                        push!(region_edges,region_edges[face1][1:vertex_index-1])
                        region_edges[face1]=vcat(-edge_count,region_edges[face1][vcat(vertex_index2:length(region_edges[face1]),vertex_index:vertex_index2-1)],edge_count)
                    else
                        throw(ArgumentError("vg2g pattern 2/3 is ambiguous: boundary indices coincide"))
                    end
                    _update_face!(face_graph,abs.(region_edges[new_face]),face1,face2,new_face,edge_count)
                    face1=face2
                else
                    # Patterns 4 and 6: traverse the face path. Intermediate
                    # boundaries receive virtual crossings.
                    pop!(strand_edges[strand_index]); edge_count-=1
                    for path_index in 1:path_length-1
                        face1=path[path_index]; face2=path[path_index+1]
                        if path_index<path_length-1
                            # Add a virtual crossing if we have to cut through
                            # an edge between adjacent faces.
                            edge_count+=1
                            boundary_indices=_find_edges(face_graph,face1,face2)
                            isempty(boundary_indices) && throw(ArgumentError(
                                "vg2g face path has no boundary between faces $face1 and $face2"))
                            # MATLAB: `bnd=bnds(1)` — ここが選択の余地あり.
                            boundary=face_graph.edges[first(boundary_indices)][3]
                            vertex_index=findfirst(==(boundary),abs.(region_edges[face1]))+1
                            edge_orientation=region_edges[face1][vertex_index-1]>0
                            edge_sign=2edge_orientation-1
                            # Put the new vertex in the middle of
                            # region_vertices[face1][vertex_index-[0,1]].
                            _insert_after_prefix!(region_vertices[face1],vertex_index-1,vertex_count;
                                context="pattern 4 face1 vertices")
                            _insert_after_prefix!(region_edges[face1],vertex_index-1-(!edge_orientation),
                                edge_sign*edge_count;context="pattern 4 face1 edges")
                            vertex_index2=findfirst(==(boundary),abs.(region_edges[face2]))+1
                            _insert_after_prefix!(region_vertices[face2],vertex_index2-1,-vertex_count;
                                context="pattern 4 face2 vertices")
                            _insert_after_prefix!(region_edges[face2],vertex_index2-1-edge_orientation,
                                -edge_sign*edge_count;context="pattern 4 face2 edges")
                            _add_edges!(face_graph,face1,face2,edge_count)
                            virtual_sign=vertex_count
                            gauss_index=findfirst(==(strand[vertex_position]),vcat(gauss[strand_index],[0]))
                            _insert_after_prefix!(gauss[strand_index],gauss_index-1,
                                -edge_sign*virtual_sign;context="current strand Gauss insertion")
                            previous_strand,previous_edge=_find_nested(boundary,strand_edges)
                            _insert_after_prefix!(strand_edges[previous_strand],previous_edge,
                                edge_count;context="previous strand edge insertion")
                            _insert_after_prefix!(gauss[previous_strand],previous_edge,
                                edge_sign*virtual_sign;context="previous strand Gauss insertion")
                            push!(orientation,0)
                            vertex=vertex_count
                            vertex_count+=1
                        elseif vertex_position==length(strand)
                            # Pattern 6: connect ends.
                            break
                        else
                            vertex=abs(strand[vertex_position])
                            vertex_index=findfirst(==(vertex),region_vertices[face1])
                            face2,vertex_index2=_find_nested(-vertex,region_vertices)
                        end
                        edge_count+=1
                        new_face=length(region_vertices)+1
                        # Pattern 4: cut through edges and add a virtual crossing.
                        push!(region_vertices,vcat(-vertex,region_vertices[face1][vertex_index+1:end]))
                        # MATLAB source: `ここの2はなぜ？ 1にしたらどうなる？`
                        region_vertices[face1]=vcat(-vertex,region_vertices[face1][2:vertex_index-1])
                        region_vertices[face2]=vcat(vertex_count+1,region_vertices[face2][vcat(vertex_index2:end,1:vertex_index2)])
                        push!(region_edges,region_edges[face1][vertex_index:end])
                        region_edges[face1]=region_edges[face1][1:vertex_index-1]
                        region_edges[face2]=vcat(-edge_count,region_edges[face2][vcat(vertex_index2:end,1:vertex_index2-1)],edge_count)
                        _update_face!(face_graph,abs.(region_edges[new_face]),face1,face2,new_face,edge_count)
                        face1=face2
                        push!(strand_edges[strand_index],edge_count)
                    end
                end
            else
                # Pattern 5: add a new crossing to the temporary graph.
                push!(used_vertices,vertex)
                signed_vertex=orientation[vertex]*strand[vertex_position]
                region_vertices[face1]=vcat(vertex_count+1,-signed_vertex,region_vertices[face1][2:end],signed_vertex)
                region_edges[face1]=vcat(-edge_count,region_edges[face1],edge_count)
                _add_edges!(face_graph,face1,face1,edge_count)
            end
            DEBUG_MODE && _debug_planarization_checkpoint(face_graph,region_vertices,region_edges,
                strand_edges,gauss,orientation,input,real_orientation;phase=:visit_exit,
                strand_index=strand_index,vertex_position=vertex_position,
                current_face=face1,io=debug_io)
        end
        # Pattern 6: connect ends of the strand.
        zero_index=findfirst(==(0),region_vertices[face1])
        new_face=length(region_vertices)+1
        push!(region_vertices,region_vertices[face1][zero_index+1:end])
        push!(region_edges,region_edges[face1][zero_index+1:end])
        region_vertices[face1]=region_vertices[face1][2:zero_index-1]
        region_edges[face1]=region_edges[face1][vcat(2:zero_index-2,1)]
        _update_face!(face_graph,abs.(region_edges[new_face]),face1,face2,new_face,edge_count)
        self_edges=_find_edges(face_graph,face2,face2)
        # MATLAB `rmedge(f2,f2)` removes every parallel self-edge.
        !isempty(self_edges) && _remove_edges!(face_graph,self_edges)
        DEBUG_MODE && _debug_planarization_checkpoint(face_graph,region_vertices,region_edges,
            strand_edges,gauss,orientation,input,real_orientation;phase=:strand_closed,
            strand_index=strand_index,current_face=face1,io=debug_io)
    end
    DEBUG_MODE && _debug_planarization_checkpoint(face_graph,region_vertices,region_edges,
        strand_edges,gauss,orientation,input,real_orientation;phase=:complete,io=debug_io)
    gauss,orientation
end
