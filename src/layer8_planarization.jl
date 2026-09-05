"""Layer8/Braid fringe translation for real Gauss diagrams.

This file translates the bounded conversion core in Andrew Bartholomew's
Braid 31.3 `gauss-to-peer.cpp`: Gauss-code edge data, PD data, fringe
selection/reordering, and virtual-crossing insertion.  The Julia adapter keeps
only the virtual-Gauss witness needed by `VirtualLink`; invariant calculations
and the rest of the Braid programme are deliberately out of scope.

The source algorithm is not claimed to be virtual-crossing optimal.  Its
published source comments explicitly describe greedy tie choices.  This route
therefore establishes an executable Layer8-derived planarisation path, not a
minimum-virtual-crossing theorem.
"""

struct Layer8FringeResult
    pd::Matrix{Int}
    virtual_crossings::Vector{NTuple{2,Int}}
    arc_directions::Vector{Int}
    virtual_counts::Vector{Int}
    initial_crossing::Int
end

# The source uses 0-based edge labels and the following direction enum.
const _L8_UNKNOWN = 0
const _L8_UPWARDS = 1
const _L8_DOWNWARDS = 2
const _L8_SIDEWAYS = 3

_l8_mod(x,n) = mod(x,n)

function _l8_remove_fringe_edge!(fringe::Vector{Int},edge::Integer)
    kept = [value for value in fringe if value != edge]
    length(kept) == length(fringe)-2 || throw(ArgumentError(
        "Layer8 fringe edge $edge did not occur twice"))
    empty!(fringe)
    append!(fringe,kept)
    fringe
end

function _l8_add_virtual_crossing!(location::Integer,fringe::Vector{Int},
                                   counts::Vector{Int},
                                   crossings::Vector{NTuple{2,Int}})
    0 <= location < length(fringe)-1 || throw(ArgumentError(
        "Layer8 virtual-crossing location $location is outside the fringe"))
    first_edge,second_edge=fringe[Int(location)+1],fringe[Int(location)+2]
    push!(crossings,(first_edge,second_edge))
    counts[first_edge+1]+=1
    counts[second_edge+1]+=1
    fringe[Int(location)+1],fringe[Int(location)+2]=second_edge,first_edge
    fringe
end

function _l8_assign_gauss_arc_direction!(arc1::Integer,arc2::Integer,
                                         table, crossing::Integer,
                                         directions::Vector{Int})
    odd_term=table.odd_terminating[Int(crossing)+1]
    even_term=table.even_terminating[Int(crossing)+1]
    even_origin=table.even_originating[Int(crossing)+1]
    odd_origin=table.odd_originating[Int(crossing)+1]
    if (odd_term == arc1 && even_origin == arc2) ||
       (even_term == arc1 && odd_origin == arc2)
        directions[Int(arc1)+1]=_L8_UPWARDS
        directions[Int(arc2)+1]=_L8_DOWNWARDS
    else
        directions[Int(arc1)+1]=_L8_DOWNWARDS
        directions[Int(arc2)+1]=_L8_UPWARDS
    end
    directions
end

_l8_matching(row,edge) = any(==(edge),row)

function _l8_rotate_left(values::Vector{Int},amount::Integer)
    isempty(values) && return Int[]
    amount=mod(Int(amount),length(values))
    amount == 0 && return copy(values)
    vcat(values[amount+1:end],values[1:amount])
end

function _l8_rotate_right(values::Vector{Int},amount::Integer)
    isempty(values) && return Int[]
    amount=mod(Int(amount),length(values))
    amount == 0 && return copy(values)
    vcat(values[end-amount+1:end],values[1:end-amount])
end

"""Build the Layer8 generic-code rows needed by `gauss_to_peer_code`.

The signed visits in `VirtualLink` use `+` for the over visit and `-` for the
under visit.  Layer8's Gauss reader stores the second visit's over/under label,
so the conversion below agrees with the source `read_gauss_code` convention.
"""
function _l8_gauss_table(real_gauss,real_orientation)
    components=[Int.(component) for component in real_gauss]
    orientation=Int.(real_orientation)
    n=length(orientation)
    length(components)==0 && (components=[Int[]])
    all(component -> all(!=(0),component),components) || throw(ArgumentError(
        "Layer8 planarisation currently requires closed Gauss components (no cut marker 0)"))

    occurrences=[Tuple{Int,Int}[] for _ in 1:n]
    first_edge_on_component=Int[]
    edge_component=Int[]
    edge_position=Int[]
    edge=0
    for (component_index,component) in enumerate(components)
        push!(first_edge_on_component,edge)
        for (position,value) in enumerate(component)
            crossing=abs(value)-1
            0 <= crossing < n || throw(ArgumentError(
                "Gauss label $value is outside the Layer8 orientation table"))
            push!(occurrences[crossing+1],(edge,value))
            push!(edge_component,component_index)
            push!(edge_position,position)
            edge+=1
        end
    end
    num_edges=edge
    num_edges==2n || throw(ArgumentError(
        "Layer8 Gauss input has $num_edges visits for $n crossings"))
    all(length(occurrence)==2 for occurrence in occurrences) || throw(ArgumentError(
        "each Layer8 Gauss crossing must occur exactly twice"))
    num_component_edges=length.(components)

    next_edge=zeros(Int,num_edges)
    for e in 0:num_edges-1
        component=edge_component[e+1]
        first_edge=first_edge_on_component[component]
        next_edge[e+1]=first_edge+mod(edge_position[e+1],num_component_edges[component])
    end

    odd_terminating=zeros(Int,n)
    even_terminating=zeros(Int,n)
    even_originating=zeros(Int,n)
    odd_originating=zeros(Int,n)
    labels=zeros(Int,n)
    types=zeros(Int,n)
    for crossing in 1:n
        first_visit,second_visit=sort(occurrences[crossing];by=item -> item[1])
        odd_terminating[crossing]=first_visit[1]
        even_terminating[crossing]=second_visit[1]
        even_originating[crossing]=next_edge[first_visit[1]+1]
        odd_originating[crossing]=next_edge[second_visit[1]+1]
        # The second visit determines Layer8's peer-code label.
        labels[crossing]=sign(second_visit[2])
        types[crossing]=orientation[crossing]==labels[crossing] ? 1 : -1
    end

    pd=zeros(Int,n,4)
    for crossing in 1:n
        if labels[crossing]==1
            if types[crossing]==-1
                pd[crossing,:]=[odd_terminating[crossing],even_terminating[crossing],
                                 even_originating[crossing],odd_originating[crossing]]
            else
                pd[crossing,:]=[odd_terminating[crossing],odd_originating[crossing],
                                 even_originating[crossing],even_terminating[crossing]]
            end
        else
            if types[crossing]==-1
                pd[crossing,:]=[even_terminating[crossing],even_originating[crossing],
                                 odd_originating[crossing],odd_terminating[crossing]]
            else
                pd[crossing,:]=[even_terminating[crossing],odd_terminating[crossing],
                                 odd_originating[crossing],even_originating[crossing]]
            end
        end
    end
    (;components,orientation,occurrences,first_edge_on_component,num_component_edges,
      next_edge,odd_terminating,even_terminating,even_originating,odd_originating,
      labels,types,pd)
end

function _l8_candidate(real_gauss,table,initial_crossing::Int)
    pd=table.pd
    n=size(pd,1)
    num_edges=2n
    directions=fill(_L8_UNKNOWN,num_edges)
    counts=zeros(Int,num_edges)
    virtual_crossings=NTuple{2,Int}[]
    fringe=Int[pd[initial_crossing+1,j] for j in 1:4]

    for i in 0:1
        _l8_assign_gauss_arc_direction!(fringe[i+1],fringe[i+3],table,
                                        initial_crossing,directions)
    end

    # Initial Reidemeister-I/single-component handling from the source.
    crossing_loop_detected=false
    for i in 0:2
        crossing_loop_detected && break
        for j in i+1:3
            if fringe[i+1]==fringe[j+1]
                directions[fringe[i+1]+1]=_L8_SIDEWAYS
                if mod(i-j+4,4)==2
                    _l8_add_virtual_crossing!(i,fringe,counts,virtual_crossings)
                    edge=fringe[j+1]
                    _l8_remove_fringe_edge!(fringe,edge)
                else
                    edge=fringe[i+1]
                    _l8_remove_fringe_edge!(fringe,edge)
                end
                crossing_loop_detected=true
                break
            end
        end
    end

    if get(ENV,"L8_TRACE","0")=="1"
        println("L8TRACE initial=",initial_crossing," fringe=",fringe)
    end
    considered=falses(n)
    considered[initial_crossing+1]=true
    while true
        complete=true
        best_crossing=0
        best_match_count=0
        best_precedence=5
        best_last_fringe_index=0
        for crossing in 0:n-1
            considered[crossing+1] && continue
            complete=false
            row=vec(pd[crossing+1,:])
            match_count=0
            last_fringe_index=-1
            last_pd_index=0
            for pd_index in 0:3
                position=findfirst(==(row[pd_index+1]),fringe)
                position===nothing && continue
                match_count+=1
                fringe_index=position-1
                if fringe_index>last_fringe_index
                    last_fringe_index=fringe_index
                    last_pd_index=pd_index
                end
            end
            match_count==0 && continue

            precedence=0
            for j in 0:match_count-1
                row[mod(last_pd_index+j,4)+1] == fringe[last_fringe_index-j+1] ||
                    (precedence=1)
            end
            if precedence==1
                if last_fringe_index==length(fringe)-1
                    scan=0
                    while scan<length(fringe) && _l8_matching(row,fringe[scan+1])
                        last_pd_index=findfirst(==(fringe[scan+1]),row)-1
                        scan+=1
                    end
                    if scan != 0
                        last_fringe_index=scan-1
                        for j in 0:match_count-1
                            row[mod(last_pd_index+j,4)+1] ==
                                fringe[mod(last_fringe_index-j,length(fringe))+1] ||
                                (precedence=2)
                        end
                    else
                        last_fringe_index=length(fringe)-1
                        precedence=2
                    end
                else
                    precedence=2
                end
            end
            if precedence==2
                for j in 0:match_count-1
                    index=last_fringe_index-j
                    found=index>=0 && any(==(fringe[index+1]),row)
                    found || (precedence=3;break)
                end
            end
            if precedence<best_precedence ||
               (precedence==best_precedence && match_count>best_match_count)
                best_crossing=crossing
                best_precedence=precedence
                best_match_count=match_count
                best_last_fringe_index=last_fringe_index
            end
        end
        complete && break
        best_match_count>0 || throw(ArgumentError(
            "Layer8 fringe stalled before all crossings were attached"))
        last_fringe_index=best_last_fringe_index
        if get(ENV,"L8_TRACE","0")=="1"
            println("L8TRACE choose=",best_crossing," match=",best_match_count,
                    " prec=",best_precedence," last=",last_fringe_index,
                    " fringe=",fringe," events=",virtual_crossings)
        end
        row=vec(pd[best_crossing+1,:])
        matching_indices=Int[]
        for index in eachindex(fringe)
            _l8_matching(row,fringe[index]) && push!(matching_indices,index-1)
        end
        length(matching_indices)==best_match_count || throw(ArgumentError(
            "Layer8 fringe matching count is inconsistent"))

        if best_precedence==1
            left_count=0
            for index in eachindex(fringe)
                _l8_matching(row,fringe[index]) || break
                left_count+=1
            end
            right_count=0
            for index in length(fringe):-1:1
                _l8_matching(row,fringe[index]) || break
                right_count+=1
            end
            if right_count<left_count
                fringe=_l8_rotate_right(fringe,right_count)
            else
                fringe=_l8_rotate_left(fringe,left_count)
            end
        elseif best_precedence==3
            m=best_match_count
            anchor=0
            min_virtual_count=-1
            best_wrap=falses(m)
            wrap=falses(m)
            for i in 0:m-1
                fill!(wrap,false)
                count=0
                pdi=[findfirst(==(fringe[matching_indices[k]+1]),row)-1 for k in 1:m]
                for j in 0:m-1
                    direct=abs(matching_indices[j+1]-matching_indices[i+1])
                    wrapped=length(fringe)-direct
                    count+=min(direct,wrapped)
                    if wrapped==direct
                        if (matching_indices[i+1]<matching_indices[j+1] &&
                            pdi[i+1]==mod(pdi[j+1]+1,4)) ||
                           (matching_indices[j+1]<matching_indices[i+1] &&
                            pdi[j+1]==mod(pdi[i+1]+1,4))
                            wrap[j+1]=false
                        else
                            wrap[j+1]=true
                        end
                    else
                        wrap[j+1]=wrapped<direct
                    end
                end
                if min_virtual_count<0 || count<min_virtual_count
                    min_virtual_count=count
                    anchor=i
                    best_wrap=copy(wrap)
                end
            end

            left_count=0
            for i in 0:anchor-1
                best_wrap[i+1] || break
                left_count+=1
                while matching_indices[i+1]>0 &&
                      (i>0 ? matching_indices[i+1] != matching_indices[i] + 1 : true)
                    _l8_add_virtual_crossing!(matching_indices[i+1]-1,fringe,counts,virtual_crossings)
                    matching_indices[i+1]-=1
                end
            end
            if left_count>0
                fringe=_l8_rotate_left(fringe,left_count)
                matching_indices=[mod(matching_indices[mod(left_count+i,m)+1]-left_count,
                                      length(fringe)) for i in 0:m-1]
                anchor-=left_count
            else
                right_count=0
                for i in m-1:-1:anchor+1
                    best_wrap[i+1] || break
                    right_count+=1
                    while matching_indices[i+1]<length(fringe)-1 &&
                          (i<m-1 ? matching_indices[i+1] != matching_indices[i+2]-1 : true)
                        _l8_add_virtual_crossing!(matching_indices[i+1],fringe,counts,virtual_crossings)
                        matching_indices[i+1]+=1
                    end
                end
                if right_count>0
                    fringe=_l8_rotate_right(fringe,right_count)
                    new_indices=fill(0,m)
                    for i in 0:m-1
                        position=mod(right_count+i,m)+1
                        new_indices[position]=mod(matching_indices[i+1]+right_count,
                                                   length(fringe))
                    end
                    matching_indices=new_indices
                    anchor=mod(anchor+right_count,m)
                end
            end
            for i in anchor+1:m-1
                while matching_indices[i+1] != matching_indices[i]+1
                    _l8_add_virtual_crossing!(matching_indices[i+1]-1,fringe,counts,virtual_crossings)
                    matching_indices[i+1]-=1
                end
            end
            for i in anchor-1:-1:0
                while matching_indices[i+1] != matching_indices[i+2]-1
                    _l8_add_virtual_crossing!(matching_indices[i+1],fringe,counts,virtual_crossings)
                    matching_indices[i+1]+=1
                end
            end
            first_pd_index=findfirst(==(fringe[matching_indices[1]+1]),row)-1
            for i in 0:m-1
                if fringe[matching_indices[i+1]+1] != row[mod(first_pd_index-i,4)+1]
                    best_precedence=2
                    last_fringe_index=matching_indices[end]
                    break
                end
            end
        end

        diametric=false
        crossing_edge_c=0
        crossing_edge_d=0
        if best_precedence==2
            m=best_match_count
            first_pd_index=-1
            if m<4
                block_start=last_fringe_index-m+1
                pd_indices=[findfirst(==(fringe[block_start+i+1]),row)-1 for i in 0:m-1]
                if m==2 && abs(pd_indices[2]-pd_indices[1])==2
                    diametric=true
                    push!(virtual_crossings,(row[mod(pd_indices[1]+1,4)+1],
                                             fringe[last_fringe_index]))
                    counts[end] # keep this expression type-stable; counts are updated below
                    pair=last(virtual_crossings)
                    counts[pair[1]+1]+=1
                    counts[pair[2]+1]+=1
                    crossing_edge_c=row[mod(pd_indices[1]+1,4)+1]
                    crossing_edge_d=row[mod(pd_indices[2]+1,4)+1]
                else
                    nonmatching=findfirst(index -> !(index in pd_indices),0:3)-1
                    for i in 0:3
                        candidate=mod(nonmatching-i+4,4)
                        candidate in pd_indices && (first_pd_index=candidate;break)
                    end
                end
            else
                first_fringe_index=last_fringe_index-m+1
                least=10
                for start in 0:3
                    virtual_count=0
                    for j in 0:3
                        for k in 0:3
                            if row[mod(start-j,4)+1] == fringe[first_fringe_index+k+1]
                                virtual_count=max(virtual_count,abs(k-j))
                                break
                            end
                        end
                    end
                    if virtual_count<least
                        least=virtual_count
                        first_pd_index=start
                    end
                end
            end
            if !diametric
                first_fringe_index=last_fringe_index-m+1
                for i in 0:m-1
                    for j in 0:m-1
                        if fringe[first_fringe_index+j+1] == row[mod(first_pd_index-i,4)+1]
                            start=first_fringe_index+j-1
                            stop=first_fringe_index+i
                            if start>=stop
                                for k in start:-1:stop
                                    _l8_add_virtual_crossing!(k,fringe,counts,virtual_crossings)
                                end
                            end
                        end
                    end
                end
            end
        end

        first_fringe_index=length(fringe)
        first_pd_index=0
        for j in 0:3
            position=findfirst(==(row[j+1]),fringe)
            position===nothing && continue
            if position-1<first_fringe_index
                first_fringe_index=position-1
                first_pd_index=j
            end
        end
        first_fringe_index<length(fringe) || throw(ArgumentError(
            "Layer8 selected crossing has no fringe incidence"))

        if diametric
            first_fringe_index+1<length(fringe) || throw(ArgumentError(
                "Layer8 diametric replacement exceeds the fringe"))
            fringe[first_fringe_index+1]=crossing_edge_d
            fringe[first_fringe_index+2]=crossing_edge_c
            if crossing_edge_c==crossing_edge_d
                directions[crossing_edge_c+1]=_L8_SIDEWAYS
                _l8_remove_fringe_edge!(fringe,crossing_edge_c)
            else
                _l8_assign_gauss_arc_direction!(crossing_edge_d,crossing_edge_c,
                                                table,best_crossing,directions)
            end
        else
            match=best_match_count
            new_edges=Int[row[mod(first_pd_index+1+i,4)+1] for i in 0:3-match]
            prefix=first_fringe_index==0 ? Int[] : fringe[1:first_fringe_index]
            suffix=first_fringe_index+match==length(fringe) ? Int[] :
                   fringe[first_fringe_index+match+1:end]
            old_matching=fringe[first_fringe_index+1:first_fringe_index+match]
            new_fringe=vcat(prefix,new_edges,suffix)
            if match==1
                arc1=row[mod(first_pd_index+1,4)+1]
                arc2=row[mod(first_pd_index+2,4)+1]
                arc3=row[mod(first_pd_index+3,4)+1]
                if arc2==arc1
                    directions[arc2+1]=_L8_SIDEWAYS
                    directions[arc3+1]=directions[old_matching[1]+1]
                    _l8_remove_fringe_edge!(new_fringe,arc2)
                elseif arc2==arc3
                    directions[arc1+1]=directions[old_matching[1]+1]
                    directions[arc2+1]=_L8_SIDEWAYS
                    _l8_remove_fringe_edge!(new_fringe,arc2)
                elseif arc1==arc3
                    directions[arc1+1]=_L8_SIDEWAYS
                    directions[arc2+1]=directions[old_matching[1]+1]
                    pair=(arc1,arc2)
                    push!(virtual_crossings,pair)
                    counts[arc1+1]+=1; counts[arc2+1]+=1
                    _l8_remove_fringe_edge!(new_fringe,arc1)
                else
                    _l8_assign_gauss_arc_direction!(arc1,arc3,table,best_crossing,directions)
                    directions[arc2+1]=directions[old_matching[1]+1]
                end
            elseif match==2
                edge1=row[mod(first_pd_index+1,4)+1]
                edge2=row[mod(first_pd_index+2,4)+1]
                direction1=directions[old_matching[1]+1]
                direction2=directions[old_matching[2]+1]
                if direction1==_L8_UPWARDS
                    if direction2==_L8_UPWARDS
                        directions[edge1+1]=_L8_UPWARDS; directions[edge2+1]=_L8_UPWARDS
                    elseif edge1==edge2
                        directions[edge1+1]=_L8_SIDEWAYS
                        _l8_remove_fringe_edge!(new_fringe,edge1)
                    else
                        directions[edge1+1]=_L8_DOWNWARDS; directions[edge2+1]=_L8_UPWARDS
                    end
                elseif direction2==_L8_UPWARDS
                    if edge1==edge2
                        directions[edge1+1]=_L8_SIDEWAYS
                        _l8_remove_fringe_edge!(new_fringe,edge1)
                    else
                        directions[edge1+1]=_L8_UPWARDS; directions[edge2+1]=_L8_DOWNWARDS
                    end
                else
                    directions[edge1+1]=_L8_DOWNWARDS; directions[edge2+1]=_L8_DOWNWARDS
                end
            elseif match==3
                up_count=count(==( _L8_UPWARDS),directions[old_matching .+ 1])
                edge1=row[mod(first_pd_index+1,4)+1]
                directions[edge1+1]=up_count==2 ? _L8_UPWARDS : _L8_DOWNWARDS
            end
            fringe=new_fringe
        end
        considered[best_crossing+1]=true
    end

    Layer8FringeResult(copy(pd),copy(virtual_crossings),directions,counts,initial_crossing)
end

function _l8_best_fringe(real_gauss,table)
    n=size(table.pd,1)
    best=nothing
    for initial_crossing in 0:n-1
        candidate=_l8_candidate(real_gauss,table,initial_crossing)
        if best===nothing || length(candidate.virtual_crossings)<length(best.virtual_crossings)
            best=candidate
        end
    end
    best===nothing && throw(ArgumentError("Layer8 could not choose an initial crossing"))
    best
end

function _l8_cycle_crossing_peers!(crossing_peers::Matrix{Int},first_edge::Int,
                                   component_edges::Int)
    last_edge=first_edge+component_edges-1
    cycled=copy(crossing_peers)
    for edge in first_edge:last_edge
        for row in axes(crossing_peers,1)
            if crossing_peers[row,1]==edge
                cycled[row,1]=edge==last_edge ? first_edge : edge+1
                break
            elseif crossing_peers[row,2]==edge
                cycled[row,2]=edge==last_edge ? first_edge : edge+1
                break
            end
        end
    end
    crossing_peers .= cycled
    crossing_peers
end

function _l8_peer_gauss(real_gauss,real_orientation,table,result::Layer8FringeResult)
    num_gauss_crossings=length(real_orientation)
    num_virtual_crossings=length(result.virtual_crossings)
    num_immersion_crossings=num_gauss_crossings+num_virtual_crossings
    num_gauss_arcs=2num_gauss_crossings
    immersion_peers=zeros(Int,num_immersion_crossings,2)
    first_immersion_edge=fill(0,length(table.components))
    next_immersion_edge=0
    component=1

    # Translate the source's Gauss-arc-to-immersion-edge allocation.  A
    # downwards arc consumes the virtual-crossing list forward; an upwards or
    # sideways arc consumes it in reverse.
    for arc in 0:num_gauss_arcs-1
        if component<length(table.components) &&
           arc==table.first_edge_on_component[component]+table.num_component_edges[component]
            component+=1
            first_immersion_edge[component]=next_immersion_edge
        end
        if result.virtual_counts[arc+1] != 0
            order=result.arc_directions[arc+1]==_L8_DOWNWARDS ?
                collect(1:num_virtual_crossings) :
                collect(num_virtual_crossings:-1:1)
            for event in order
                pair=result.virtual_crossings[event]
                row=num_gauss_crossings+event
                if pair[1]==arc
                    immersion_peers[row,1]=next_immersion_edge
                    next_immersion_edge+=1
                end
                if pair[2]==arc
                    immersion_peers[row,2]=next_immersion_edge
                    next_immersion_edge+=1
                end
            end
        end
        for crossing in 1:num_gauss_crossings
            if table.even_terminating[crossing]==arc
                immersion_peers[crossing,1]=next_immersion_edge
                next_immersion_edge+=1
                break
            elseif table.odd_terminating[crossing]==arc
                immersion_peers[crossing,2]=next_immersion_edge
                next_immersion_edge+=1
                break
            end
        end
    end
    next_immersion_edge==2num_immersion_crossings || throw(ArgumentError(
        "Layer8 immersion-edge allocation ended at $next_immersion_edge, expected $(2num_immersion_crossings)"))

    num_immersion_component_edges=Int[]
    for index in eachindex(first_immersion_edge)
        next_first=index<length(first_immersion_edge) ?
            first_immersion_edge[index+1] : next_immersion_edge
        push!(num_immersion_component_edges,next_first-first_immersion_edge[index])
    end

    # The source cycles a component when a peer pair has equal parity.  This
    # restores the even/odd terminating-edge convention before peer output.
    if length(table.components)>1
        queue=Int[1]
        considered=falses(length(table.components))
        while !isempty(queue)
            current=popfirst!(queue)
            considered[current] && continue
            considered[current]=true
            first_edge=first_immersion_edge[current]
            last_edge=first_edge+num_immersion_component_edges[current]-1
            for edge in first_edge:last_edge
                row=findfirst(r -> immersion_peers[r,1]==edge ||
                                      immersion_peers[r,2]==edge,
                              axes(immersion_peers,1))
                row===nothing && throw(ArgumentError("Layer8 peer edge $edge is unpaired"))
                peer=immersion_peers[row,1]==edge ? immersion_peers[row,2] : immersion_peers[row,1]
                peer_component=1
                for candidate in 2:length(first_immersion_edge)
                    peer >= first_immersion_edge[candidate] && (peer_component=candidate)
                end
                push!(queue,peer_component)
                if iseven(peer)==iseven(edge)
                    _l8_cycle_crossing_peers!(immersion_peers,
                        first_immersion_edge[peer_component],
                        num_immersion_component_edges[peer_component])
                end
            end
        end
    end

    peer_odd=fill(-1,num_immersion_crossings)
    peer_even=fill(-1,num_immersion_crossings)
    for row in axes(immersion_peers,1)
        a,b=immersion_peers[row,1],immersion_peers[row,2]
        even_edge=iseven(a) ? a : b
        odd_edge=iseven(a) ? b : a
        iseven(even_edge) && isodd(odd_edge) || throw(ArgumentError(
            "Layer8 peer crossing $row does not have one even and one odd edge"))
        peer_odd[even_edge÷2+1]=odd_edge
        peer_even[(odd_edge-1)÷2+1]=even_edge
    end

    component_of_edge(edge)=begin
        value=1
        for index in 2:length(first_immersion_edge)
            edge >= first_immersion_edge[index] && (value=index)
        end
        value
    end
    peer_successor(edge)=begin
        owner=component_of_edge(edge)
        first_edge=first_immersion_edge[owner]
        first_edge+mod(edge-first_edge+1,num_immersion_component_edges[owner])
    end

    peer_labels=zeros(Int,num_immersion_crossings)
    peer_types=fill(1,num_immersion_crossings)
    # First assign the real crossings.  The peer crossing index is the naming
    # (even) edge divided by two, as in gauss_to_peer_code's output table.
    for crossing in 1:num_gauss_crossings
        a,b=immersion_peers[crossing,1],immersion_peers[crossing,2]
        even_in_place_zero=iseven(a)
        even_edge=even_in_place_zero ? a : b
        peer_crossing=even_edge÷2+1
        peer_types[peer_crossing]=even_in_place_zero ? table.types[crossing] : -table.types[crossing]
        peer_labels[peer_crossing]=even_in_place_zero ? table.labels[crossing] : -table.labels[crossing]
    end
    for event in 1:num_virtual_crossings
        row=num_gauss_crossings+event
        a,b=immersion_peers[row,1],immersion_peers[row,2]
        even_in_place_zero=iseven(a)
        even_edge=even_in_place_zero ? a : b
        peer_crossing=even_edge÷2+1
        first_direction=result.arc_directions[result.virtual_crossings[event][1]+1]
        second_direction=result.arc_directions[result.virtual_crossings[event][2]+1]
        same_direction=(first_direction in (_L8_UPWARDS,_L8_DOWNWARDS) &&
                        first_direction==second_direction)
        # This is the non-diametric branch of gauss_to_peer_code's virtual
        # TYPE assignment.  Diametric events are retained as a deferred
        # extension until the peer witness carries that source flag.
        peer_types[peer_crossing]=same_direction ?
            (even_in_place_zero ? -1 : 1) :
            (even_in_place_zero ? 1 : -1)
    end

    if get(ENV,"L8_DEBUG_PEER","0")=="1"
        println("L8DEBUG immersion_peers=",immersion_peers)
        println("L8DEBUG first=",first_immersion_edge," counts=",num_immersion_component_edges)
        println("L8DEBUG peer_odd=",peer_odd," peer_even=",peer_even)
        println("L8DEBUG labels=",peer_labels," types=",peer_types)
    end

    if get(ENV,"L8_DEBUG_PEER","0")=="1"
        target=zeros(Int,num_immersion_crossings,4)
        for crossing in 1:num_immersion_crossings
            odd_term=peer_odd[crossing]
            even_term=2(crossing-1)
            even_origin=peer_successor(odd_term)
            odd_origin=2(crossing-1)+1
            if peer_labels[crossing]==1
                target[crossing,:]=peer_types[crossing]==-1 ?
                    [odd_term,even_term,even_origin,odd_origin] :
                    [odd_term,odd_origin,even_origin,even_term]
            else
                target[crossing,:]=peer_types[crossing]==-1 ?
                    [even_term,even_origin,odd_origin,odd_term] :
                    [even_term,odd_term,odd_origin,even_origin]
            end
        end
        println("L8DEBUG target_pd=",target)
    end

    # Serialise the peer traversal back to signed Gauss visits.  The
    # component is traced through each crossing, not by incrementing the edge
    # label: at a crossing the outgoing edge is selected from the parity of
    # the terminating edge exactly as in `write_gauss_code`.
    term_crossing=fill(-1,2num_immersion_crossings)
    odd_origin=collect(1:2:2num_immersion_crossings-1)
    even_origin=zeros(Int,num_immersion_crossings)
    for crossing in 1:num_immersion_crossings
        even_edge=2(crossing-1)
        odd_edge=peer_odd[crossing]
        term_crossing[even_edge+1]=crossing
        term_crossing[odd_edge+1]=crossing
        odd_origin[crossing]=2(crossing-1)+1
        even_origin[crossing]=peer_successor(odd_edge)
    end
    gauss=[Int[] for _ in table.components]
    for component in eachindex(table.components)
        start=first_immersion_edge[component]
        edge=start
        while true
            crossing=term_crossing[edge+1]
            crossing>0 || throw(ArgumentError("Layer8 edge $edge has no terminating crossing"))
            label=peer_labels[crossing]==0 ? 1 : peer_labels[crossing]
            under=(isodd(edge) && label==1) || (iseven(edge) && label==-1)
            sign_value=under ? -1 : 1
            push!(gauss[component],sign_value*crossing)
            edge=isodd(edge) ? even_origin[crossing] : odd_origin[crossing]
            edge==start && break
        end
    end
    peer_orientation=[peer_labels[index]==0 ? 0 :
                      (peer_types[index]==peer_labels[index] ? 1 : -1)
                      for index in eachindex(peer_labels)]
    gauss,peer_orientation
end

"""Translate Layer8 fringe events into a virtual Gauss code.

The original source ultimately serializes a labelled peer code.  For the
Julia plotting boundary we reconstruct the equivalent signed Gauss traversal:
each fringe arc receives its Layer8 event order, and the two visits of each
inserted virtual crossing get opposite temporary signs.  Since its stored
orientation is `0`, that sign pair is only the local frame needed by
`pd_code`; it does not make the crossing real.
"""
function _l8_virtual_gauss(real_gauss,real_orientation,result::Layer8FringeResult)
    components=[Int.(component) for component in real_gauss]
    real_count=length(real_orientation)
    event_count=length(result.virtual_crossings)
    arc_occurrences=[Tuple{Int,Int}[] for _ in 1:2real_count]
    for (event,pair) in enumerate(result.virtual_crossings)
        push!(arc_occurrences[pair[1]+1],(event,1))
        push!(arc_occurrences[pair[2]+1],(event,2))
    end
    offsets=Int[]
    edge=0
    for component in components
        push!(offsets,edge)
        edge+=length(component)
    end
    length(arc_occurrences)==edge || throw(ArgumentError(
        "Layer8 Gauss arc count does not match the input traversal"))
    seen=Dict{Int,Int}()
    output=Vector{Vector{Int}}()
    for (component_index,component) in enumerate(components)
        values=Int[]
        offset=offsets[component_index]
        for (position,value) in enumerate(component)
            push!(values,value)
            arc=offset+position-1
            occurrences=arc_occurrences[arc+1]
            if result.arc_directions[arc+1] != _L8_DOWNWARDS
                occurrences=reverse(occurrences)
            end
            for (event,_) in occurrences
                label=real_count+event
                sign_value=get(seen,event,0)==0 ? 1 : -1
                seen[event]=get(seen,event,0)+1
                push!(values,sign_value*label)
            end
        end
        push!(output,values)
    end
    all(get(seen,event,0)==2 for event in 1:event_count) || throw(ArgumentError(
        "Layer8 virtual crossing reconstruction lost an event visit"))
    output,vcat(Int.(real_orientation),zeros(Int,event_count))
end

"""Layer8/Braid fringe planarisation for a real Gauss presentation."""
function planarize_virtual_gauss_layer8(real_gauss,real_orientation;
                                        DEBUG_MODE::Bool=false,
                                        debug_io::IO=stderr)
    input=[Int.(component) for component in real_gauss]
    orientation=Int.(real_orientation)
    isempty(orientation) && return deepcopy(input),orientation
    table=_l8_gauss_table(input,orientation)
    result=_l8_best_fringe(input,table)
    gauss,planar_orientation=_l8_peer_gauss(input,orientation,table,result)
    DEBUG_MODE && begin
        println(debug_io,"[layer8] initial_crossing=",result.initial_crossing,
                " virtual_crossings=",length(result.virtual_crossings))
        println(debug_io,"[layer8] virtual_pairs=",result.virtual_crossings)
        println(debug_io,"[layer8] gauss=",gauss)
        println(debug_io,"[layer8] orientation=",planar_orientation)
    end
    gauss,planar_orientation
end
