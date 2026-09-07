"""A concrete applicability witness returned by `movable`.

The candidate freezes the selected ordering, direction, templates, and one
location per source fragment.  `move!` consumes this witness instead of
searching a second time with potentially different choices.
"""
struct MoveCandidate
    type::String
    parameter::String
    vertices::Vector{Int}
    edges::Vector{NTuple{2,Int}}
    direction::Symbol
    before::Union{Nothing,MovePattern}
    after::Union{Nothing,MovePattern}
    locations::Vector{Tuple{Int,Int}}
end

struct _MoveSegment
    start_port::Int
    end_port::Int
    word::Vector{Int}
    weights::Vector{Float64}
    leading_weight::Float64
end

_normalize_move_type(type)=uppercase(String(type))
_clean_component(component)=filter(!=(0),Int.(component))
_is_cut_component(component)=!isempty(component) && last(component)==0

function _permutations(values::Vector{Int})
    length(values)<=1 && return [copy(values)]
    result=Vector{Vector{Int}}()
    for index in eachindex(values)
        rest=copy(values); value=splice!(rest,index)
        for suffix in _permutations(rest)
            push!(result,vcat(value,suffix))
        end
    end
    result
end

function _cyclic_occurrences(component::Vector{Int},pattern::Vector{Int})
    isempty(pattern) && return Int[]
    length(pattern)<=length(component) || return Int[]
    [start for start in eachindex(component)
        if all(component[mod1(start+offset,length(component))]==pattern[offset+1]
               for offset in 0:length(pattern)-1)]
end

_interval_positions(start,length_,component_length)=
    Set(mod1(start+offset,component_length) for offset in 0:length_-1)

function _locations_compatible(locations,fragments,code)
    occupied=Set{Tuple{Int,Int}}()
    for (fragment,(component,start)) in zip(fragments,locations)
        n=length(code[component])
        for position in _interval_positions(start,length(fragment),n)
            key=(component,position)
            key in occupied && return false
            push!(occupied,key)
        end
    end
    true
end

function _location_products(options::Vector{Vector{Tuple{Int,Int}}},fragments,code)
    result=Vector{Vector{Tuple{Int,Int}}}()
    current=Tuple{Int,Int}[]
    function visit(index)
        if index>length(options)
            _locations_compatible(current,fragments,code) && push!(result,copy(current))
            return
        end
        for location in options[index]
            push!(current,location); visit(index+1); pop!(current)
        end
    end
    visit(1)
    result
end

function _map_fragment(fragment,vertices)
    [entry==0 ? 0 : sign(entry)*vertices[abs(entry)] for entry in fragment]
end

function _pattern_pairs(type::String,count::Int,data)
    if type=="MP"
        count==2 && return data["MP-L"],data["MP-R"],:forward
        count==3 && return data["MP-R"],data["MP-L"],:inverse
    elseif type=="CP"
        count==3 && return data["CP-L"],data["CP-R"],:forward
        count==5 && return data["CP-R"],data["CP-L"],:inverse
    elseif type=="PS"
        count==2 && return data["PS-R"],data["PS-L"],:inverse
    elseif type=="02"
        count==2 && return data["02-R"],data["02-L"],:inverse
    elseif type=="BMP"
        count==2 && return data["BMP-L"],data["BMP-R"],:forward
        count==3 && return data["BMP-R"],data["BMP-L"],:inverse
    elseif type=="B02"
        count==2 && return data["B02-R"],data["B02-L"],:inverse
    end
    throw(ArgumentError("invalid number of vertices for $type move"))
end

function _after_patterns(before::MovePattern,after_table,param)
    same=filter(pattern -> pattern.parameter==before.parameter,after_table)
    !isempty(same) && return same
    length(after_table)==1 && return after_table
    if isempty(before.parameter) && param!==nothing
        return filter(pattern -> pattern.parameter==string(param),after_table)
    end
    MovePattern[]
end

function _coerce_edges(edges)
    edges===nothing && return NTuple{2,Int}[]
    edges isa AbstractMatrix && return [Tuple(Int.(row)) for row in eachrow(edges)]
    result=NTuple{2,Int}[]
    for edge in edges
        values=Int.(collect(edge)); length(values)==2 || throw(DimensionMismatch("edges must have two endpoints"))
        push!(result,(values[1],values[2]))
    end
    result
end

"""Return all data-table matches for a MATLAB-style VirtualLink move request."""
function movable(v::VirtualLink,type;param=nothing,vertices=Int[],edges=nothing,strict=false)
    move_type=_normalize_move_type(type)
    selected=Int.(collect(vertices))
    edge_list=_coerce_edges(edges)
    code,orientation=real_gauss_code(v)
    code=[_clean_component(component) for component in code]

    if move_type in ("R1","R2","R3")
        return MoveCandidate[]
    elseif move_type=="H"
        length(selected)==1 || throw(ArgumentError("H move requires one vertex"))
        1<=abs(selected[1])<=length(orientation) || return MoveCandidate[]
        return [MoveCandidate(move_type,param===nothing ? "" : string(param),selected,edge_list,
                              :weight,nothing,nothing,Tuple{Int,Int}[])]
    elseif move_type=="02" && string(param)=="dir"
        length(selected)==2 || throw(ArgumentError("direct 0-2 move requires two signed adjacent visits"))
        locations=Tuple{Int,Int}[]
        for vertex in selected
            choices=[(component,position) for (component,word) in enumerate(code)
                     for position in findall(==(vertex),word)]
            length(choices)==1 || return MoveCandidate[]
            push!(locations,only(choices))
        end
        first_location,second_location=locations
        same_component=first_location[1]==second_location[1]
        adjacent=same_component && mod1(first_location[2]+1,length(code[first_location[1]]))==second_location[2]
        adjacent || return MoveCandidate[]
        data=make_move_data()
        return [MoveCandidate(move_type,"dir",selected,edge_list,:forward,
                              only(data["02-L"]),only(data["02-R"]),locations)]
    elseif move_type=="BMP" && length(selected)==2
        selected[1]!=selected[2] || throw(ArgumentError("BMP vertices must be distinct"))
        values=abs.(selected)
        all(vertex -> vertex>0 && vertex<=length(orientation),values) || return MoveCandidate[]
        positions=Dict{Int,Tuple{Int,Int}}()
        for vertex in (values[1],values[2],-values[1],-values[2])
            occurrences=[(component,position) for (component,word) in enumerate(code)
                         for position in findall(==(vertex),word)]
            length(occurrences)==1 || return MoveCandidate[]
            positions[vertex]=only(occurrences)
        end
        v1,v2=values
        p1,p2,p3,p4=positions[v1],positions[v2],positions[-v1],positions[-v2]
        acceptable=Bool[
            p1[1]==p2[1] && p1[2]+1==p2[2] && orientation[values]==[-1,1],
            p3[1]==p4[1] && p3[2]+1==p4[2] && orientation[values]==[1,-1],
            p2[1]==p1[1] && p2[2]+1==p1[2] && orientation[values]==[1,-1],
            p4[1]==p3[1] && p4[2]+1==p3[2] && orientation[values]==[-1,1],
        ]
        if strict
            preferred=selected[1]>0 ? 1 : 2
            acceptable=[acceptable[preferred]]
            !acceptable[1] && return MoveCandidate[]
            values=selected[1]>0 ? values : -values
        else
            any(acceptable) || return MoveCandidate[]
        end
        return [MoveCandidate(move_type,"",Int.(values),edge_list,:forward,
                              nothing,nothing,Tuple{Int,Int}[])]
    elseif move_type=="B02" && length(selected)==2
        all(selected .> 0) || throw(ArgumentError("B02 vertices must be positive"))
        selected[1]!=selected[2] || throw(ArgumentError("B02 vertices must be distinct"))
        values=copy(selected)
        all(vertex -> vertex<=length(orientation),values) || return MoveCandidate[]
        orientation[values[1]]*orientation[values[2]]==-1 || return MoveCandidate[]
        orientation[values[1]]==-1 && (values=reverse(values))
        matches=any(component -> component==values || component==-values ||
                         component==reverse(values) || component==-reverse(values),code)
        matches || return MoveCandidate[]
        return [MoveCandidate(move_type,"",copy(selected),edge_list,:inverse,
                              nothing,nothing,Tuple{Int,Int}[])]
    end

    isempty(selected) && throw(ArgumentError("$move_type move requires selected vertices"))
    all(vertex -> vertex>0 && vertex<=length(orientation),selected) || return MoveCandidate[]
    data=make_move_data()
    before_table,after_table,direction=_pattern_pairs(move_type,length(selected),data)
    orders=strict ? [selected] : _permutations(selected)
    candidates=MoveCandidate[]
    seen=Set{Any}()
    for ordered in orders
        for before in before_table
            param!==nothing && !(string(param) in (before.parameter,"")) && continue
            length(before.orientation)==length(ordered) || continue
            orientation[ordered]==before.orientation || continue
            mapped_before=[_map_fragment(fragment,ordered) for fragment in before.gauss]
            any(isempty,mapped_before) && continue
            options=Vector{Vector{Tuple{Int,Int}}}()
            valid=true
            for fragment in mapped_before
                locations=Tuple{Int,Int}[]
                for (component_index,component) in enumerate(code)
                    append!(locations,[(component_index,start)
                                       for start in _cyclic_occurrences(component,fragment)])
                end
                isempty(locations) && (valid=false; break)
                push!(options,locations)
            end
            valid || continue
            for after in _after_patterns(before,after_table,param)
                for locations in _location_products(options,mapped_before,code)
                    parameter=isempty(before.parameter) ? after.parameter : before.parameter
                    key=(parameter,Tuple(ordered),Tuple(locations),direction)
                    key in seen && continue
                    push!(seen,key)
                    push!(candidates,MoveCandidate(move_type,parameter,copy(ordered),edge_list,
                                                   direction,before,after,locations))
                end
            end
        end
    end
    # For PS, local template matching is only a candidate generator.  The
    # closed-normal contract makes the explicit validator the single filter;
    # all other move families retain their historical behavior unchanged.
    move_type=="PS" ? [candidate for candidate in candidates
                        if validate_ps_move(v,candidate).valid] : candidates
end

function movable(v::VirtualLink,type,arg::NamedTuple)
    param=haskey(arg,:param) ? arg.param : haskey(arg,:p) ? arg.p : nothing
    vertices=haskey(arg,:vertices) ? arg.vertices : haskey(arg,:v) ? arg.v : Int[]
    edges=haskey(arg,:edges) ? arg.edges : haskey(arg,:e) ? arg.e : nothing
    strict=haskey(arg,:strict) ? arg.strict : false
    movable(v,type;param=param,vertices=vertices,edges=edges,strict=strict)
end

function _vertex_maps(orientation,selected,before,after;role_ordered=false)
    old_count=length(orientation)
    delta=after.vertex_count-before.vertex_count
    if delta>=0
        old_to_new=Dict(index=>index for index in 1:old_count)
        after_vertices=vcat(selected,collect(old_count+1:old_count+delta))
        new_orientation=vcat(copy(orientation),zeros(Int,delta))
    else
        remove_count=-delta
        ordered_selected=role_ordered ? selected : sort(selected)
        removed=Set(ordered_selected[end-remove_count+1:end])
        remaining=[vertex for vertex in 1:old_count if !(vertex in removed)]
        old_to_new=Dict(vertex=>index for (index,vertex) in enumerate(remaining))
        kept=ordered_selected[1:after.vertex_count]
        after_vertices=[old_to_new[vertex] for vertex in kept]
        new_orientation=orientation[remaining]
    end
    length(after_vertices)==after.vertex_count || throw(ArgumentError("move-data vertex counts are inconsistent"))
    for local_vertex in eachindex(after.orientation)
        new_orientation[after_vertices[local_vertex]]=after.orientation[local_vertex]
    end
    old_to_new,after_vertices,new_orientation
end

_relabel_entry(entry,mapping)=entry==0 ? 0 : sign(entry)*get(mapping,abs(entry),abs(entry))
_relabel_word(word,mapping)=[_relabel_entry(entry,mapping) for entry in word]

function _cyclic_gap_positions(start,length_,next_start,n)
    positions=Int[]
    position=mod1(start+length_,n)
    while position!=next_start
        push!(positions,position)
        position=mod1(position+1,n)
        length(positions)>n && throw(ArgumentError("overlapping cyclic move fragments"))
    end
    positions
end

function _check_weight_rule(before_rule,after_rule)
    (any(isnan,before_rule) || any(isnan,after_rule)) &&
        throw(ArgumentError("weight transport is not specified for this move parameter"))
end

function _simple_pattern_rewrite(code,weights,candidate,mapped_before,mapped_after,
                                 old_to_new,after_orientation;weighted=false)
    # Relabel retained source entries before inserting target fragments.  Target
    # fragments already use the post-move numbering and must not be mapped twice.
    result=[_relabel_word(component,old_to_new) for component in code]
    result_weights=deepcopy(weights)
    by_component=Dict{Int,Vector{Int}}()
    for (index,(component,_)) in enumerate(candidate.locations)
        push!(get!(by_component,component,Int[]),index)
    end
    for (component,indices) in by_component
        n=length(result[component])
        all(index -> candidate.locations[index][2]+length(mapped_before[index])-1<=n,indices) ||
            return nothing
        sort!(indices;by=index -> candidate.locations[index][2],rev=true)
        for index in indices
            start=candidate.locations[index][2]
            consumed=length(mapped_before[index])
            replacement=mapped_after[index]
            before_rule=candidate.before.weights[index]
            after_rule=candidate.after.weights[index]
            weighted && _check_weight_rule(before_rule,after_rule)
            old_component=result[component]
            old_weights=result_weights[component]
            incoming=weighted ? old_weights[start] : 0.0
            replacement_weights=zeros(Float64,length(replacement))
            if weighted && !isempty(replacement)
                replacement_weights[1]=incoming+first(after_rule)-first(before_rule)
                for position in 2:length(replacement)
                    replacement_weights[position]=after_rule[position]
                end
            end
            result[component]=vcat(old_component[1:start-1],replacement,old_component[start+consumed:end])
            result_weights[component]=vcat(old_weights[1:start-1],replacement_weights,old_weights[start+consumed:end])
            if weighted && !isempty(result_weights[component])
                next_index=start+length(replacement)
                next_index>length(result_weights[component]) && (next_index=1)
                result_weights[component][next_index]+=last(after_rule)-last(before_rule)
            end
        end
    end
    result,result_weights,after_orientation
end

function _append_segment!(word,weights,segment,pending)
    next_pending=pending+segment.leading_weight
    if !isempty(segment.word)
        segment_weights=copy(segment.weights)
        isempty(segment_weights) && (segment_weights=zeros(Float64,length(segment.word)))
        segment_weights[1]+=next_pending
        append!(word,segment.word); append!(weights,segment_weights)
        next_pending=0.0
    end
    next_pending
end

function _reconnect_pattern(code,weights,candidate,mapped_before,mapped_after,
                            old_to_new,new_orientation;weighted=false)
    before=candidate.before; after=candidate.after
    common=min(length(before.gauss),length(after.gauss))
    extra_before=collect(common+1:length(before.gauss))
    removed_components=Set{Int}()
    for index in extra_before
        component,start=candidate.locations[index]
        fragment=mapped_before[index]
        length(code[component])==length(fragment) ||
            throw(ArgumentError("a source-only move fragment must occupy its whole component"))
        _cyclic_occurrences(code[component],fragment)==[start] ||
            throw(ArgumentError("cannot remove a partially matched component"))
        push!(removed_components,component)
    end

    segments=_MoveSegment[]
    affected=Set(component for (index,(component,_)) in enumerate(candidate.locations) if index<=common)
    for component in affected
        _is_cut_component(code[component]) && throw(ArgumentError("cap reconnection on cut components is not supported by MATLAB move data"))
        indices=[index for index in 1:common if candidate.locations[index][1]==component]
        sort!(indices;by=index -> candidate.locations[index][2])
        for (position,index) in enumerate(indices)
            next_index=indices[mod1(position+1,length(indices))]
            start=candidate.locations[index][2]
            next_start=candidate.locations[next_index][2]
            gap_positions=_cyclic_gap_positions(start,length(mapped_before[index]),next_start,length(code[component]))
            gap=_relabel_word(code[component][gap_positions],old_to_new)
            gap_weights=weighted ? copy(weights[component][gap_positions]) : zeros(Float64,length(gap))
            leading=0.0
            if weighted
                before_rule=before.weights[index]; after_rule=after.weights[index]
                _check_weight_rule(before_rule,after_rule)
                delta=last(after_rule)-last(before_rule)
                isempty(gap_weights) ? (leading=delta) : (gap_weights[1]+=delta)
            end
            push!(segments,_MoveSegment(before.capacity[index],-next_index,gap,gap_weights,leading))
        end
    end

    for index in 1:common
        word=mapped_after[index]
        fragment_weights=zeros(Float64,length(word))
        leading=0.0
        if weighted
            before_rule=before.weights[index]; after_rule=after.weights[index]
            _check_weight_rule(before_rule,after_rule)
            component,start=candidate.locations[index]
            incoming=weights[component][start]
            if isempty(word)
                leading=incoming+first(after_rule)-first(before_rule)
            else
                fragment_weights[1]=incoming+first(after_rule)-first(before_rule)
                for position in 2:length(word)
                    fragment_weights[position]=after_rule[position]
                end
            end
        end
        push!(segments,_MoveSegment(-index,after.capacity[index],word,fragment_weights,leading))
    end

    result=Vector{Vector{Int}}()
    result_weights=Vector{Vector{Float64}}()
    for component in eachindex(code)
        (component in affected || component in removed_components) && continue
        push!(result,_relabel_word(code[component],old_to_new))
        push!(result_weights,weighted ? copy(weights[component]) : zeros(Float64,length(code[component])))
    end

    starts=Dict{Int,Int}()
    for (index,segment) in enumerate(segments)
        haskey(starts,segment.start_port) && throw(ArgumentError("ambiguous move capacity $(segment.start_port)"))
        starts[segment.start_port]=index
    end
    used=falses(length(segments))
    while !all(used)
        first_index=findfirst(!,used)
        current=first_index
        word=Int[]; component_weights=Float64[]; pending=0.0
        open_component=false
        while true
            used[current] && current!=first_index && throw(ArgumentError("move capacity graph joins an already used strand"))
            used[current]=true
            segment=segments[current]
            pending=_append_segment!(word,component_weights,segment,pending)
            if segment.end_port==0
                open_component=true
                break
            end
            haskey(starts,segment.end_port) || throw(ArgumentError("unmatched move capacity $(segment.end_port)"))
            next=starts[segment.end_port]
            next==first_index && break
            current=next
        end
        if weighted && pending!=0 && !isempty(component_weights)
            component_weights[1]+=pending
        end
        open_component && push!(word,0)
        push!(result,word)
        push!(result_weights,component_weights)
    end

    for index in common+1:length(after.gauss)
        word=mapped_after[index]
        push!(result,copy(word))
        if weighted
            rule=after.weights[index]; _check_weight_rule(rule,rule)
            component_weights=zeros(Float64,length(word))
            for position in eachindex(component_weights)
                component_weights[position]=rule[min(position+1,length(rule))]
            end
            !isempty(component_weights) && (component_weights[1]+=first(rule)+last(rule))
            push!(result_weights,component_weights)
        else
            push!(result_weights,zeros(Float64,length(word)))
        end
    end
    result,result_weights,new_orientation
end

function _apply_pattern_move!(v::VirtualLink,candidate::MoveCandidate)
    before=something(candidate.before); after=something(candidate.after)
    code,orientation=real_gauss_code(v)
    code=[_clean_component(component) for component in code]
    weights=v.is_weighted ? deepcopy(v.weights) : [zeros(Float64,length(component)) for component in code]
    v.is_weighted && candidate.type!="MP" &&
        throw(ArgumentError("weighted $(candidate.type) move is not implemented in the MATLAB source"))
    old_to_new,after_vertices,new_orientation=_vertex_maps(
        orientation,candidate.vertices,before,after;
        role_ordered=candidate.type=="CP")
    mapped_before=[_map_fragment(fragment,candidate.vertices) for fragment in before.gauss]
    mapped_after=[_map_fragment(fragment,after_vertices) for fragment in after.gauss]

    simple=length(before.gauss)==length(after.gauss) && before.capacity==after.capacity
    converted=simple ? _simple_pattern_rewrite(code,weights,candidate,mapped_before,mapped_after,
                                                old_to_new,new_orientation;weighted=v.is_weighted) : nothing
    converted===nothing &&
        (converted=_reconnect_pattern(code,weights,candidate,mapped_before,mapped_after,
                                      old_to_new,new_orientation;weighted=v.is_weighted))
    result,result_weights,result_orientation=converted
    if v.is_weighted
        set_data!(v;odata=(result,result_orientation,result_weights))
    else
        set_data!(v;rgauss=result,orientation=result_orientation)
    end
    v
end

"""Validate one frozen PS witness in the closed-normal Julia subcategory.

Callers must use `Topology.validate_ps_move` explicitly; `movable` remains the
public applicability gate.  This validator checks that the frozen candidate is
an authentic current-source PS-R witness before applying it, then requires the
source and tentative target to be closed.  It never calls `movable`, and never
mutates the source.  Its closed-normal criterion is intentionally not a
boundary or weighted-PS criterion: those cases retain explicit rejection.

The returned named tuple is deliberately inspectable.  `source_conditions` and
`target_conditions` are `(C1,C2,C3)`, `target` is an independently constructed
 tentative link when application reaches that stage, and `error` is `nothing`
except for an expected tentative-application input failure.
"""
function _same_move_pattern(left::MovePattern,right::MovePattern)
    left.parameter==right.parameter && left.gauss==right.gauss &&
    left.orientation==right.orientation && left.capacity==right.capacity &&
    left.gauss_first==right.gauss_first && left.gauss_last==right.gauss_last &&
    left.vertex_count==right.vertex_count && left.weights==right.weights
end

function validate_ps_move(v::VirtualLink,candidate::MoveCandidate)
    source_closed,source_conditions=is_closed(v)
    invalid(reason,target=nothing,target_closed=false,
            target_conditions=(false,false,false);error=nothing)=
        (valid=false,reason=reason,source_closed=source_closed,
         source_conditions=source_conditions,target_closed=target_closed,
         target_conditions=target_conditions,target=target,error=error)
    _normalize_move_type(candidate.type)=="PS" || return invalid(:not_ps)
    v.is_weighted && return invalid(:weighted_unsupported)
    isempty(candidate.edges) || return invalid(:edge_selected_unsupported)
    (candidate.before===nothing || candidate.after===nothing) &&
        return invalid(:malformed_candidate)
    candidate.direction===:inverse || return invalid(:source_mismatch)
    parameter=try parse(Int,candidate.parameter) catch; return invalid(:source_mismatch) end
    1<=parameter<=4 || return invalid(:source_mismatch)
    data=make_move_data()
    canonical_before=data["PS-R"][parameter]
    canonical_after=only(data["PS-L"])
    _same_move_pattern(candidate.before,canonical_before) &&
        _same_move_pattern(candidate.after,canonical_after) ||
        return invalid(:source_mismatch)
    source_closed || return invalid(:source_not_closed)
    code,orientation=real_gauss_code(v)
    code=[_clean_component(component) for component in code]
    length(candidate.vertices)==2 || return invalid(:source_mismatch)
    vertices=candidate.vertices
    all(vertex -> vertex>0 && vertex<=length(orientation),vertices) ||
        return invalid(:source_mismatch)
    vertices[1]!=vertices[2] || return invalid(:source_mismatch)
    candidate.before.orientation==orientation[vertices] || return invalid(:source_mismatch)
    length(candidate.locations)==length(candidate.before.gauss) ||
        return invalid(:source_mismatch)
    mapped_before=[_map_fragment(fragment,vertices) for fragment in candidate.before.gauss]
    for (fragment,(component,start)) in zip(mapped_before,candidate.locations)
        1<=component<=length(code) || return invalid(:source_mismatch)
        n=length(code[component]); 1<=start<=n || return invalid(:source_mismatch)
        start in _cyclic_occurrences(code[component],fragment) ||
            return invalid(:source_mismatch)
    end
    _locations_compatible(candidate.locations,mapped_before,code) ||
        return invalid(:source_mismatch)
    target=deepcopy(v)
    try
        _apply_pattern_move!(target,candidate)
    catch err
        err isa InterruptException && rethrow()
        err isa Union{ArgumentError,BoundsError,DimensionMismatch} || rethrow()
        return invalid(:tentative_application_failed,error=err)
    end
    target_closed,target_conditions=is_closed(target)
    target_closed || return invalid(:target_not_closed,target,target_closed,target_conditions)
    (valid=true,reason=:closed_local_ps,source_closed=source_closed,
     source_conditions=source_conditions,target_closed=target_closed,
     target_conditions=target_conditions,target=target,error=nothing)
end

function _apply_direct_02!(v::VirtualLink,candidate::MoveCandidate)
    code,orientation=real_gauss_code(v)
    v.is_weighted && throw(ArgumentError("weighted direct 0-2 transport follows MATLAB's unfinished ad-hoc branch"))
    component,start=candidate.locations[1]
    new_vertices=length(orientation).+(1:2)
    fragment=[-new_vertices[1],-new_vertices[2],new_vertices[2],new_vertices[1]]
    word=code[component]
    code[component]=vcat(word[1:start],fragment,word[start+1:end])
    set_data!(v;rgauss=code,orientation=vcat(orientation,[1,-1]))
end

function _apply_bmp_two!(v::VirtualLink,candidate::MoveCandidate,parameter)
    parameter in (0,1) || throw(ArgumentError("BMP circle orientation must be 0 or 1"))
    code,orientation=real_gauss_code(v)
    values=Int.(candidate.vertices)
    S=sign(values[1])
    v1,v2=abs.(values)
    locate(value)=begin
        locations=[(component,position) for (component,word) in enumerate(code)
                   for position in findall(==(value),word)]
        length(locations)==1 || throw(ArgumentError("BMP vertex occurrence is ambiguous"))
        only(locations)
    end
    component1,position1=locate(-values[1])
    word=code[component1]
    code[component1]=vcat(word[position1+1:end],word[1:position1-1])
    component2,_=locate(-values[2])
    component0,position0=locate(values[1])
    new_vertex=S*(maximum(reduce(vcat,code;init=Int[]))+1)
    code[component0][position0]=new_vertex
    position0=findfirst(==(values[2]),code[component0])
    position0===nothing && throw(ArgumentError("BMP second vertex disappeared during rotation"))
    deleteat!(code[component0],position0)
    component2,position2=locate(-values[2])
    circle=parameter==1 ? [-values[1],-values[2],-new_vertex] :
                           [-values[2],-values[1],-new_vertex]
    new_orientation=vcat(orientation,0)
    circle_orientation=parameter==1 ? -[-1,1,1]*S : [-1,1,1]*S
    for (vertex,value) in zip((v1,v2,abs(new_vertex)),circle_orientation)
        new_orientation[vertex]=value
    end
    if component1==component2
        word=code[component2]
        code=vcat(code,
                  [[values[1];word[1:position2-1]],
                   [values[2];word[position2+1:end]]])
        code[component2]=circle
    else
        word=code[component2]
        code[component2]=vcat([values[2]],word[position2+1:end],word[1:position2-1],
                              [values[1]],code[component1])
        code[component1]=circle
    end
    set_data!(v;rgauss=code,orientation=new_orientation)
end

function _apply_b02_two!(v::VirtualLink,candidate::MoveCandidate)
    code,orientation=real_gauss_code(v)
    values=copy(candidate.vertices)
    nv=length(orientation)
    mapping=Dict(old=>new for (new,old) in enumerate(setdiff(1:nv,values)))
    is_positive=any(component -> component==values || component==-values,code)
    is_positive && (values=-values)
    locate(value)=begin
        locations=[(component,position) for (component,word) in enumerate(code)
                   for position in findall(==(value),word)]
        length(locations)==1 || throw(ArgumentError("B02 vertex occurrence is ambiguous"))
        only(locations)
    end
    component0,_=locate(-values[1])
    component2,position2=locate(values[2])
    word=code[component2]
    code[component2]=vcat(word[position2+1:end],word[1:position2-1])
    component1,position1=locate(values[1])
    if component1==component2
        word=code[component1]
        code[component0]=word[1:position1-1]
        code[component1]=word[position1+1:end]
    else
        word1=code[component1]; word2=code[component2]
        code[component0]=vcat(word2,word1[position1+1:end],word1[1:position1-1])
        deleteat!(code,sort([component1,component2]))
    end
    result=[_relabel_word(word,mapping) for word in code]
    new_orientation=copy(orientation)
    deleteat!(new_orientation,sort(abs.(values)))
    set_data!(v;rgauss=result,orientation=new_orientation)
end

function _apply_hmove!(v::VirtualLink,vertex::Int,amount)
    v.is_weighted || throw(ArgumentError("H move requires a weighted VirtualLink"))
    code,_=real_gauss_code(v)
    pairs=NTuple{2,Int}[]
    for signed_vertex in (vertex,-vertex)
        locations=[(component,position) for (component,word) in enumerate(code)
                   for position in findall(==(signed_vertex),word)]
        length(locations)==1 || throw(ArgumentError("H-move vertex occurrence is ambiguous"))
        component,position=only(locations); word=code[component]
        push!(pairs,(word[mod1(position-1,length(word))],signed_vertex))
        push!(pairs,(signed_vertex,word[mod1(position+1,length(word))]))
    end
    adjustments=Dict{NTuple{2,Int},Float64}()
    for (pair,delta) in zip(pairs,(-1,1,-1,1))
        adjustments[pair]=get(adjustments,pair,0.0)+delta*Float64(amount)
    end
    table=real_edge_table(v)
    by_strand=[Float64[] for _ in code]
    for edge in table
        push!(by_strand[edge.strand],edge.weight+get(adjustments,edge.crossing,0.0))
    end
    set_weight!(v,by_strand)
end

"""Apply one MATLAB-style move after strict applicability selection."""
function move!(v::VirtualLink,type;param=nothing,vertices=Int[],edges=nothing,strict=true)
    move_type=_normalize_move_type(type)
    if move_type in ("R1","R2","R3")
        throw(ArgumentError("$move_type move is not implemented in the MATLAB source"))
    elseif move_type=="H"
        param===nothing && throw(ArgumentError("H move requires a weight parameter"))
        candidates=movable(v,move_type;param=param,vertices=vertices,edges=edges,strict=strict)
        isempty(candidates) && throw(ArgumentError("No valid moves found."))
        return _apply_hmove!(v,abs(only(candidates).vertices[1]),param)
    end
    candidates=movable(v,move_type;param=param,vertices=vertices,edges=edges,strict=strict)
    isempty(candidates) && throw(ArgumentError("No valid moves found."))
    length(candidates)>1 && @warn "moves may be ambiguous; applying the first strict candidate" count=length(candidates)
    candidate=first(candidates)
    if move_type=="02" && candidate.parameter=="dir"
        _apply_direct_02!(v,candidate)
    elseif move_type=="BMP" && length(candidate.vertices)==2
        _apply_bmp_two!(v,candidate,param)
    elseif move_type=="B02" && length(candidate.vertices)==2 && candidate.before===nothing
        _apply_b02_two!(v,candidate)
    else
        _apply_pattern_move!(v,candidate)
    end
end

function move!(v::VirtualLink,type,arg::NamedTuple)
    param=haskey(arg,:param) ? arg.param : haskey(arg,:p) ? arg.p : nothing
    vertices=haskey(arg,:vertices) ? arg.vertices : haskey(arg,:v) ? arg.v : Int[]
    edges=haskey(arg,:edges) ? arg.edges : haskey(arg,:e) ? arg.e : nothing
    strict=haskey(arg,:strict) ? arg.strict : true
    move!(v,type;param=param,vertices=vertices,edges=edges,strict=strict)
end

"""Pure-value wrapper around MATLAB-compatible `move!`."""
function move(v::VirtualLink,type;kwargs...)
    result=deepcopy(v); move!(result,type;kwargs...); result
end
move(v::VirtualLink,type,arg::NamedTuple)=move!(deepcopy(v),type,arg)

"""Compatibility wrapper for the first Julia migration API."""
mp_move!(v::VirtualLink,parameter::AbstractString,vertices)=
    move!(v,"MP";param=parameter,vertices=vertices,strict=true)
