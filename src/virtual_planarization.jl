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

function _update_face!(graph::_FaceGraph,boundary,face1,face2,new_face,new_edge)
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

"""Insert virtual crossings so a real-Gauss/o-data diagram becomes planar."""
function planarize_virtual_gauss(real_gauss,real_orientation)
    input=[Int.(component) for component in real_gauss]
    orientation=Int.(real_orientation)
    isempty(orientation) && return deepcopy(input),orientation

    region_vertices=Vector{Vector{Int}}()
    region_edges=Vector{Vector{Int}}()
    gauss=deepcopy(input)
    strand_edges=[Int[] for _ in gauss]
    used_vertices=Set([0])
    edge_count=0
    vertex_count=maximum(abs.(reduce(vcat,input;init=Int[])))+1
    face_graph=_FaceGraph()

    for strand_index in eachindex(input)
        strand=copy(input[strand_index])
        (isempty(strand) || strand == [0]) && continue
        cut=last(strand)==0
        strand=vcat(strand[1:end-cut],[0])
        length(strand)==1 && (input[strand_index]=Int[]; continue)
        edge_count+=1
        push!(strand_edges[strand_index],edge_count)
        first_vertex=abs(strand[1])
        face1=0
        face2=0
        if first_vertex in used_vertices
            face2,vertex_index2=_find_nested(-first_vertex,region_vertices)
            region_vertices[face2]=vcat(vertex_count+1,region_vertices[face2][vcat(vertex_index2:length(region_vertices[face2]),1:vertex_index2)])
            region_edges[face2]=vcat(-edge_count,region_edges[face2][vertex_index2:end],region_edges[face2][1:vertex_index2-1],edge_count)
            _add_edges!(face_graph,face2,[face2,face2],[0,edge_count])
            face1,vertex_index=_find_nested(first_vertex,region_vertices)
            region_vertices[face1]=vcat(region_vertices[face1][1:vertex_index-1],-first_vertex,0,-first_vertex,region_vertices[face1][vertex_index+1:end])
            region_edges[face1]=vcat(region_edges[face1][1:vertex_index-1],[0,0],region_edges[face1][vertex_index:end])
            face1=face2
        else
            face1=length(region_vertices)+1
            sign_value=orientation[first_vertex]*strand[1]
            push!(region_vertices,[vertex_count+1,-sign_value,0,sign_value])
            push!(region_edges,[-edge_count,0,0,edge_count])
            _add_edges!(face_graph,face1,[face1,face1],[0,edge_count])
            push!(used_vertices,first_vertex)
        end

        for vertex_position in 2:length(strand)
            edge_count+=1
            push!(strand_edges[strand_index],edge_count)
            vertex=abs(strand[vertex_position])
            if vertex in used_vertices
                face2,vertex_index2=_find_nested(vertex,region_vertices)
                new_face,_=_find_nested(-vertex,region_vertices)
                path=vcat(_shortest_path(face_graph,face1,face2),new_face)
                path_length=length(path)
                if path_length==1
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
                    vertex_index=vertex_index2
                    face2,vertex_index2=_find_nested(-vertex,region_vertices)
                    new_face=length(region_vertices)+1
                    if vertex_index2<vertex_index
                        push!(region_vertices,region_vertices[face1][vcat(vertex_index2,vertex_index+1:length(region_vertices[face1]))])
                        region_vertices[face1]=region_vertices[face1][vcat(1,vertex_index2:vertex_index-1,vertex_index2,2:vertex_index2)]
                        push!(region_edges,region_edges[face1][vertex_index:end])
                        region_edges[face1]=vcat(-edge_count,region_edges[face1][vcat(vertex_index2:vertex_index-1,1:vertex_index2-1)],edge_count)
                    else
                        push!(region_vertices,region_vertices[face1][vcat(vertex_index2,2:vertex_index-1)])
                        region_vertices[face1]=region_vertices[face1][vcat(1,vertex_index2:length(region_vertices[face1]),vertex_index2,vertex_index+1:vertex_index2)]
                        push!(region_edges,region_edges[face1][1:vertex_index-1])
                        region_edges[face1]=vcat(-edge_count,region_edges[face1][vcat(vertex_index2:length(region_edges[face1]),vertex_index:vertex_index2-1)],edge_count)
                    end
                    _update_face!(face_graph,abs.(region_edges[new_face]),face1,face2,new_face,edge_count)
                    face1=face2
                else
                    pop!(strand_edges[strand_index]); edge_count-=1
                    for path_index in 1:path_length-1
                        face1=path[path_index]; face2=path[path_index+1]
                        if path_index<path_length-1
                            edge_count+=1
                            boundary_indices=_find_edges(face_graph,face1,face2)
                            boundary=face_graph.edges[first(boundary_indices)][3]
                            vertex_index=findfirst(==(boundary),abs.(region_edges[face1]))+1
                            edge_orientation=region_edges[face1][vertex_index-1]>0
                            edge_sign=2edge_orientation-1
                            insert!(region_vertices[face1],vertex_index,vertex_count)
                            insert!(region_edges[face1],vertex_index-1-(!edge_orientation),edge_sign*edge_count)
                            vertex_index2=findfirst(==(boundary),abs.(region_edges[face2]))+1
                            insert!(region_vertices[face2],vertex_index2,-vertex_count)
                            insert!(region_edges[face2],vertex_index2-1-edge_orientation,-edge_sign*edge_count)
                            _add_edges!(face_graph,face1,face2,edge_count)
                            virtual_sign=vertex_count
                            gauss_index=findfirst(==(strand[vertex_position]),vcat(gauss[strand_index],[0]))
                            insert!(gauss[strand_index],gauss_index,-edge_sign*virtual_sign)
                            previous_strand,previous_edge=_find_nested(boundary,strand_edges)
                            insert!(strand_edges[previous_strand],previous_edge+1,edge_count)
                            insert!(gauss[previous_strand],previous_edge+1,edge_sign*virtual_sign)
                            push!(orientation,0)
                            vertex=vertex_count
                            vertex_count+=1
                        elseif vertex_position==length(strand)
                            break
                        else
                            vertex=abs(strand[vertex_position])
                            vertex_index=findfirst(==(vertex),region_vertices[face1])
                            face2,vertex_index2=_find_nested(-vertex,region_vertices)
                        end
                        edge_count+=1
                        new_face=length(region_vertices)+1
                        push!(region_vertices,vcat(-vertex,region_vertices[face1][vertex_index+1:end]))
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
                push!(used_vertices,vertex)
                signed_vertex=orientation[vertex]*strand[vertex_position]
                region_vertices[face1]=vcat(vertex_count+1,-signed_vertex,region_vertices[face1][2:end],signed_vertex)
                region_edges[face1]=vcat(-edge_count,region_edges[face1],edge_count)
                _add_edges!(face_graph,face1,face1,edge_count)
            end
        end
        zero_index=findfirst(==(0),region_vertices[face1])
        new_face=length(region_vertices)+1
        push!(region_vertices,region_vertices[face1][zero_index+1:end])
        push!(region_edges,region_edges[face1][zero_index+1:end])
        region_vertices[face1]=region_vertices[face1][2:zero_index-1]
        region_edges[face1]=region_edges[face1][vcat(2:zero_index-2,1)]
        _update_face!(face_graph,abs.(region_edges[new_face]),face1,face2,new_face,edge_count)
        self_edges=_find_edges(face_graph,face2,face2)
        !isempty(self_edges) && _remove_edges!(face_graph,[first(self_edges)])
    end
    gauss,orientation
end
