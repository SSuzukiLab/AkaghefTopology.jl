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
    ps_split_k::Union{Nothing,BigInt}
end

# Preserve all existing geometric candidate callers and unweighted evidence.
MoveCandidate(type,parameter,vertices,edges,direction,before,after,locations)=
    MoveCandidate(type,parameter,vertices,edges,direction,before,after,locations,nothing)

function _check_candidate_split(v::VirtualLink,candidate::MoveCandidate)
    candidate.ps_split_k===nothing && return
    v.is_weighted && candidate.type=="PS" && candidate.direction==:forward ||
        throw(ArgumentError("a PS split parameter requires a weighted forward PS candidate"))
    _move_integer(candidate.ps_split_k)
    nothing
end

struct _MoveSegment
    start_port::Int
    end_port::Int
    word::Vector{Int}
    weights::Vector{Float64}
    leading_weight::Float64
end

# Capacity labels describe interior seams; boundary labels describe actual
# source endpoints and must never identify two different cut ends.
_move_boundary_port(common,strand,side)=-(common+2strand-Int(side==:in))
function _move_boundaries(code,affected,common)
    Dict(_move_boundary_port(common,strand,side)=>(strand=strand,side=side)
         for strand in affected if _is_cut_component(code[strand]) for side in (:in,:out))
end

function _move_segment_paths(segments,boundaries)
    starts=Dict{Int,Int}(); ends=Dict{Int,Int}()
    for (index,segment) in enumerate(segments)
        haskey(starts,segment.start_port) && throw(ArgumentError("duplicate fragment start marker"))
        haskey(ends,segment.end_port) && throw(ArgumentError("duplicate fragment end marker"))
        starts[segment.start_port]=index; ends[segment.end_port]=index
    end
    for (port,index) in starts
        if haskey(boundaries,port)
            boundaries[port].side==:in || throw(ArgumentError("outgoing boundary used as fragment start"))
        else
            haskey(ends,port) || throw(ArgumentError("unmatched fragment start marker $port"))
        end
    end
    for port in keys(ends)
        if haskey(boundaries,port)
            boundaries[port].side==:out || throw(ArgumentError("incoming boundary used as fragment end"))
        else
            haskey(starts,port) || throw(ArgumentError("unmatched fragment end marker $port"))
        end
    end
    roots=sort!([index for (port,index) in starts if haskey(boundaries,port)];
                by=index->boundaries[segments[index].start_port].strand)
    append!(roots,eachindex(segments))
    used=falses(length(segments)); paths=NamedTuple[]
    for root in roots
        used[root] && continue
        incoming=get(boundaries,segments[root].start_port,nothing)
        indices=Int[]; current=root; outgoing=nothing
        while true
            used[current] && throw(ArgumentError("fragment path enters an already consumed segment"))
            used[current]=true; push!(indices,current)
            port=segments[current].end_port
            if haskey(boundaries,port)
                outgoing=boundaries[port]
                incoming===nothing && throw(ArgumentError("open path has no incoming boundary"))
                break
            end
            current=starts[port]
            current==root && break
        end
        (incoming===nothing)==(outgoing===nothing) || throw(ArgumentError("boundary path closed before its endpoint"))
        push!(paths,(indices=indices,incoming=incoming,outgoing=outgoing))
    end
    paths
end

function _install_move_paths(code,weights,affected,removed,mapping,paths,words,values;
                             boundary_trace=nothing)
    result=Vector{Vector{Int}}(); result_weights=Vector{Vector{Float64}}()
    open_by_source=Dict(path.incoming.strand=>index for (index,path) in enumerate(paths)
                        if path.incoming!==nothing)
    function record_endpoint(source,side)
        boundary_trace===nothing || push!(boundary_trace,
            (source=source,target=(strand=length(result),side=side)))
    end
    for strand in eachindex(code)
        if haskey(open_by_source,strand)
            index=open_by_source[strand]; path=paths[index]
            push!(result,words[index]); push!(result_weights,values[index])
            record_endpoint(path.incoming,:in); record_endpoint(path.outgoing,:out)
        elseif !(strand in affected || strand in removed)
            push!(result,_relabel_word(code[strand],mapping)); push!(result_weights,copy(weights[strand]))
            if _is_cut_component(code[strand])
                record_endpoint((strand=strand,side=:in),:in)
                record_endpoint((strand=strand,side=:out),:out)
            end
        end
    end
    for (index,path) in enumerate(paths)
        path.incoming===nothing || continue
        push!(result,words[index]); push!(result_weights,values[index])
    end
    result,result_weights
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

"""Contiguous occurrences on the linear prefix of an open component."""
function _linear_occurrences(component::Vector{Int},pattern::Vector{Int})
    isempty(pattern) && return Int[]
    n=_is_cut_component(component) ? length(component)-1 : length(component)
    length(pattern)<=n || return Int[]
    [start for start in 1:(n-length(pattern)+1)
        if component[start:start+length(pattern)-1]==pattern]
end

_pattern_occurrences(component,pattern)=
    _is_cut_component(component) ? _linear_occurrences(component,pattern) :
                                    _cyclic_occurrences(component,pattern)

_interval_positions(start,length_,component_length)=
    Set(mod1(start+offset,component_length) for offset in 0:length_-1)

function _locations_compatible(locations,fragments,code)
    occupied=Set{Tuple{Int,Int}}()
    for (fragment,(component,start)) in zip(fragments,locations)
        n=_is_cut_component(code[component]) ? length(code[component])-1 : length(code[component])
        positions=if _is_cut_component(code[component])
            start>=1 && start+length(fragment)-1<=n ? collect(start:start+length(fragment)-1) : Int[]
        else
            collect(_interval_positions(start,length(fragment),n))
        end
        isempty(fragment) || isempty(positions) && return false
        for position in positions
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
        count==0 && return data["PS-L"],data["PS-R"],:forward
        count==2 && return data["PS-R"],data["PS-L"],:inverse
    elseif type=="S"
        count==0 && return data["S-L"],data["S-R"],:forward
        count==4 && return data["S-R"],data["S-L"],:inverse
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

"""Enumerate geometric data-table matches before integer transport checks."""
function _movable_geometry(v::VirtualLink,type;param=nothing,vertices=Int[],edges=nothing,
                 edge_occurrences=nothing,strict=false)
    move_type=_normalize_move_type(type)
    selected=Int.(collect(vertices))
    edge_list=_coerce_edges(edges)
    code,orientation=real_gauss_code(v)

    if move_type in ("PS","S") && (edge_occurrences!==nothing || !isempty(edge_list))
        isempty(selected) || (edge_occurrences===nothing ? (return MoveCandidate[]) :
            throw(ArgumentError("insertion selects edges, not vertices")))
        # Explicit occurrence is (component, head position): insert before head.
        locations = if edge_occurrences!==nothing
            _coerce_edges(edge_occurrences)
        else
            mapped=Tuple{Int,Int}[]
            for (tail,head) in edge_list
                matches=[(c,p) for (c,w) in enumerate(code)
                         for p in 1:length(w)
                         if w[p]==head &&
                            ((p>1 && w[p-1]==tail) ||
                             (p==1 && w[end]==tail))]
                length(matches)==1 || return MoveCandidate[]
                push!(mapped,only(matches))
            end
            mapped
        end
        required=move_type=="PS" ? 2 : 1
        length(locations)==required || throw(ArgumentError("$move_type insertion requires $required edge occurrences"))
        all(cp -> 1<=cp[1]<=length(code) &&
                  1<=cp[2]<=length(code[cp[1]]),locations) ||
            return MoveCandidate[]
        data=make_move_data()
        return [candidate for after in data[move_type*"-R"]
                if param===nothing || string(param)==after.parameter
                for candidate in [MoveCandidate(move_type,after.parameter,Int[],NTuple{2,Int}[],
                    :forward,only(data[move_type*"-L"]),after,copy(locations))]
                if move_type!="PS" || validate_ps_move(v,candidate).valid]
    end

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
        adjacent=same_component &&
            (first_location[2]+1==second_location[2] ||
             (!_is_cut_component(code[first_location[1]]) &&
              mod1(first_location[2]+1,length(code[first_location[1]]))==second_location[2]))
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
        return [MoveCandidate(move_type,param===nothing ? "" : string(param),Int.(values),edge_list,:forward,
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

    numeric_bmp_inverse=move_type=="BMP" && length(selected)==3 && string(param) in ("0","1")
    bmp_sign=numeric_bmp_inverse && 1<=abs(first(selected))<=length(orientation) ?
        (string(param)=="0" ? -orientation[abs(first(selected))] : orientation[abs(first(selected))]) : 1
    numeric_bmp_inverse && (selected=abs.(selected))
    isempty(selected) && throw(ArgumentError("$move_type move requires selected vertices"))
    all(vertex -> vertex>0 && vertex<=length(orientation),selected) || return MoveCandidate[]
    data=make_move_data()
    before_table,after_table,direction=if numeric_bmp_inverse
        p=parse(Int,string(param)); s=bmp_sign
        circle=p==0 ? [-2,-1,-3] : [-1,-2,-3]
        ori=p==0 ? [-1,1,1].*s : [1,-1,-1].*s
        b=MovePattern(string(p),[s.*word for word in [[3],[2],[1],circle]],ori,[1,3,2,0])
        a=MovePattern(string(p),[s.*word for word in [[1,2],[-1],[-2]]],[-s,s],[1,2,3])
        ([b],[a],:inverse)
    else
        _pattern_pairs(move_type,length(selected),data)
    end
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
                                       for start in _pattern_occurrences(component,fragment)])
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

"""Enumerate applicable moves in the selected coefficient mode.

Weighted candidates satisfy their transport rules. Z2 mode supports spin CP
on closed and cut presentations; other weighted families require integer mode.
"""
function movable(v::VirtualLink,type;boundary=nothing,kwargs...)
    if weight_mode(v)==:z2
        _normalize_move_type(type)=="CP" || return MoveCandidate[]
        !any(_is_cut_component,v.gauss) && boundary!==nothing &&
            throw(ArgumentError("spin CP boundary targets require a cut source"))
        _spin_edge_weights(v)
        candidates=_movable_geometry(_weighted_unweight(v),"CP";kwargs...)
        return filter(candidates) do candidate
            try
                spin_cp_move(v,candidate;boundary=boundary)
                true
            catch error
                error isa Union{ArgumentError,DimensionMismatch,BoundsError} || rethrow()
                false
            end
        end
    end
    v.is_weighted || return _movable_geometry(v,type;kwargs...)
    _check_move_weights(v)
    move_type=_normalize_move_type(type)
    move_type in ("MP","H","02","PS") || return MoveCandidate[]
    candidates=_movable_geometry(move_type=="PS" ? _ps_unweighted(v) : v,move_type;kwargs...)
    filter(candidates) do candidate
        try
            _apply_candidate!(deepcopy(v),candidate;boundary=boundary)
            true
        catch error
            error isa Union{ArgumentError,DimensionMismatch,BoundsError} || rethrow()
            false
        end
    end
end

function movable(v::VirtualLink,type,arg::NamedTuple)
    param=haskey(arg,:param) ? arg.param : haskey(arg,:p) ? arg.p : nothing
    vertices=haskey(arg,:vertices) ? arg.vertices : haskey(arg,:v) ? arg.v : Int[]
    edges=haskey(arg,:edges) ? arg.edges : haskey(arg,:e) ? arg.e : nothing
    occurrences=haskey(arg,:edge_occurrences) ? arg.edge_occurrences : nothing
    strict=haskey(arg,:strict) ? arg.strict : false
    boundary=haskey(arg,:boundary) ? arg.boundary : nothing
    movable(v,type;param=param,vertices=vertices,edges=edges,edge_occurrences=occurrences,
            strict=strict,boundary=boundary)
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
    if weighted
        reference=_simple_pattern_rewrite(code,weights,candidate,mapped_before,mapped_after,
                                           old_to_new,after_orientation;weighted=false)
        reference===nothing && return nothing
        transported,values,_=_weighted_pattern_rewrite(code,weights,candidate,mapped_before,
            mapped_after,old_to_new,after_orientation)
        return reference[1],_align_move_weights(reference[1],transported,values),after_orientation
    end
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
                            old_to_new,new_orientation;weighted=false,boundary_trace=nothing)
    if weighted
        reference=_reconnect_pattern(code,weights,candidate,mapped_before,mapped_after,
                                     old_to_new,new_orientation;weighted=false)
        transported,values,_=_weighted_pattern_rewrite(code,weights,candidate,mapped_before,
            mapped_after,old_to_new,new_orientation)
        return reference[1],_align_move_weights(reference[1],transported,values),new_orientation
    end
    before=candidate.before; after=candidate.after
    common=min(length(before.gauss),length(after.gauss))
    extra_before=collect(common+1:length(before.gauss))
    removed_components=Set{Int}()
    for index in extra_before
        component,start=candidate.locations[index]
        fragment=mapped_before[index]
        _is_cut_component(code[component]) &&
            throw(ArgumentError("a source-only fragment cannot remove cut endpoints"))
        (_is_cut_component(code[component]) ? length(code[component])-1 : length(code[component]))==length(fragment) ||
            throw(ArgumentError("a source-only move fragment must occupy its whole component"))
        _pattern_occurrences(code[component],fragment)==[start] ||
            throw(ArgumentError("cannot remove a partially matched component"))
        push!(removed_components,component)
    end

    segments=_MoveSegment[]
    affected=sort!(unique([component for (index,(component,_)) in enumerate(candidate.locations) if index<=common]))
    boundaries=_move_boundaries(code,affected,common)
    for component in affected
        indices=[index for index in 1:common if candidate.locations[index][1]==component]
        sort!(indices;by=index -> (candidate.locations[index][2],-index))
        if _is_cut_component(code[component])
            n=length(code[component])-1
            first_index=first(indices)
            first_start=candidate.locations[first_index][2]
            prefix_positions=first_start>1 ? collect(1:first_start-1) : Int[]
            prefix=_relabel_word(code[component][prefix_positions],old_to_new)
            prefix_weights=weighted ? copy(weights[component][prefix_positions]) : zeros(Float64,length(prefix))
            push!(segments,_MoveSegment(_move_boundary_port(common,component,:in),-first_index,prefix,prefix_weights,0.0))
        end
        for (position,index) in enumerate(indices)
            start=candidate.locations[index][2]
            is_last=_is_cut_component(code[component]) && position==length(indices)
            next_index=is_last ? 0 : indices[mod1(position+1,length(indices))]
            next_start=is_last ? 0 : candidate.locations[next_index][2]
            gap_positions=if is_last
                n_start=start+length(mapped_before[index])
                n_start<=n ? collect(n_start:n) : Int[]
            elseif _is_cut_component(code[component])
                n_start=start+length(mapped_before[index])
                n_start<next_start ? collect(n_start:next_start-1) : Int[]
            else
                n_closed=length(code[component])
                gap_count=next_start+(position==length(indices) ? n_closed : 0)-start-length(mapped_before[index])
                gap_count<0 && throw(ArgumentError("overlapping cyclic move fragments"))
                [mod1(start+length(mapped_before[index])+offset,n_closed) for offset in 0:gap_count-1]
            end
            gap=_relabel_word(code[component][gap_positions],old_to_new)
            gap_weights=weighted ? copy(weights[component][gap_positions]) : zeros(Float64,length(gap))
            leading=0.0
            if weighted
                before_rule=before.weights[index]; after_rule=after.weights[index]
                _check_weight_rule(before_rule,after_rule)
                delta=last(after_rule)-last(before_rule)
                isempty(gap_weights) ? (leading=delta) : (gap_weights[1]+=delta)
            end
            push!(segments,_MoveSegment(before.capacity[index],is_last ? _move_boundary_port(common,component,:out) : -next_index,
                                        gap,gap_weights,leading))
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

    paths=_move_segment_paths(segments,boundaries)
    words=Vector{Vector{Int}}(); values=Vector{Vector{Float64}}()
    for path in paths
        word=reduce(vcat,(segments[index].word for index in path.indices);init=Int[])
        path.incoming===nothing || push!(word,0)
        push!(words,word); push!(values,zeros(Float64,length(word)))
    end
    result,result_weights=_install_move_paths(code,weights,affected,removed_components,
        old_to_new,paths,words,values;boundary_trace=boundary_trace)

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

# Use exact accumulation even though the public storage is Float64.  Exclude
# the first ambiguous integer boundary (2^53), not merely nonfinite values.
const _MOVE_MAX_INTEGER=big(2)^53-1
function _move_integer(value)
    value isa Real && isfinite(value) && isinteger(value) ||
        throw(ArgumentError("weighted moves require finite integer values"))
    integer=BigInt(value)
    abs(integer)<=_MOVE_MAX_INTEGER ||
        throw(ArgumentError("weight exceeds the exact Float64 integer range"))
    integer
end

function _check_move_weights(v::VirtualLink)
    weight_mode(v)==:integer || throw(ArgumentError("integer move transport does not support Z2 mode"))
    v.is_weighted || throw(ArgumentError("move requires assigned integer weights"))
    raw_code,orientation=real_gauss_code(v)
    code=[_clean_component(component) for component in raw_code]
    all(sign -> sign in (-1,1),orientation) || throw(ArgumentError("invalid real crossing orientation"))
    visits=reduce(vcat,code;init=Int[])
    all(vertex -> count(==(vertex),visits)==1 && count(==(-vertex),visits)==1,
        eachindex(orientation)) || throw(ArgumentError("one visit of each sign is required"))
    length(v.weights)==length(raw_code) &&
        all(length(w)==length(c) for (w,c) in zip(v.weights,raw_code)) ||
        throw(DimensionMismatch("one weight per real edge occurrence is required"))
    foreach(_move_integer,Iterators.flatten(v.weights))
    nothing
end

function _commit_move!(v::VirtualLink,target::VirtualLink)
    for field in fieldnames(VirtualLink)
        setfield!(v,field,getfield(target,field))
    end
    v
end

function _normalized_02(candidate::MoveCandidate)
    candidate.type=="02" || return candidate
    empty=MovePattern("dir",[Int[]],Int[];weights=[[0]])
    if candidate.direction==:forward
        # The legacy public request selects two adjacent visits.  Both delimit
        # ONE insertion edge, unlike the two distinct slots of a PS move.
        MoveCandidate("02",candidate.parameter,Int[],candidate.edges,:forward,
                      empty,candidate.after,[candidate.locations[2]])
    else
        MoveCandidate("02",candidate.parameter,candidate.vertices,candidate.edges,:inverse,
                      candidate.before,empty,candidate.locations)
    end
end

struct _WeightedMoveSegment
    start_port::Int
    end_port::Int
    word::Vector{Int}
    # n vertices have n+1 edge portions; concatenation adds the shared portions.
    edges::Vector{BigInt}
end

function _weighted_pattern_rewrite(code,weights,candidate,mapped_before,mapped_after,
                                   old_to_new,new_orientation;trace=nothing)
    before=candidate.before; after=candidate.after
    rules0=[[_move_integer(x) for x in w] for w in before.weights]
    rules1=[[_move_integer(x) for x in w] for w in after.weights]
    for (pattern,rules) in ((before,rules0),(after,rules1))
        all(length(w)==length(g)+1 for (w,g) in zip(rules,pattern.gauss)) ||
            throw(DimensionMismatch("weight rules require both boundary edge portions"))
    end
    common=min(length(before.gauss),length(after.gauss))
    removed=Set{Int}()
    for index in eachindex(mapped_before)
        component,start=candidate.locations[index]
        fragment=mapped_before[index]; cut=_is_cut_component(code[component])
        n=length(code[component])-Int(cut)
        cut && start+length(fragment)-1>n &&
            throw(ArgumentError("weighted source fragment crosses a cut boundary"))
        for offset in 1:length(fragment)-1
            position=cut ? start+offset : mod1(start+offset,n)
            weights[component][position]==rules0[index][offset+1] ||
                throw(ArgumentError("source internal weight mismatch at fragment $index, edge $(offset+1)"))
        end
        if index>common
            cut && throw(ArgumentError("a source-only fragment cannot remove cut endpoints"))
            length(fragment)==n || throw(ArgumentError("source-only fragment must occupy a whole component"))
            weights[component][start]==first(rules0[index])+last(rules0[index]) ||
                throw(ArgumentError("source-only component closing weight mismatch"))
            push!(removed,component)
        end
    end
    affected=sort!(unique([candidate.locations[i][1] for i in 1:common]))
    boundaries=_move_boundaries(code,affected,common)
    segments=_WeightedMoveSegment[]
    for component in affected
        component in removed && throw(ArgumentError("source-only component shares an active fragment"))
        cut=_is_cut_component(code[component])
        n=length(code[component])-Int(cut)
        indices=[i for i in 1:common if candidate.locations[i][1]==component]
        sort!(indices;by=i -> (candidate.locations[i][2],-i))
        if cut
            first_index=first(indices); first_start=candidate.locations[first_index][2]
            positions=collect(1:first_start-1)
            edge_positions=collect(1:first_start)
            exterior=BigInt[_move_integer(weights[component][p]) for p in edge_positions]
            exterior[end]-=first(rules0[first_index])
            word=_relabel_word(code[component][positions],old_to_new)
            push!(segments,_WeightedMoveSegment(_move_boundary_port(common,component,:in),-first_index,word,exterior))
            trace===nothing || push!(trace,(kind=:boundary_prefix,component=component,
                head_positions=edge_positions,word=word,portions=copy(exterior),
                start_port=_move_boundary_port(common,component,:in),end_port=-first_index))
        end
        for (j,index) in enumerate(indices)
            terminal=cut && j==length(indices)
            next_index=terminal ? 0 : indices[mod1(j+1,length(indices))]
            start=candidate.locations[index][2]; consumed=length(mapped_before[index])
            next_start=terminal ? n+1 : candidate.locations[next_index][2]
            # Distinguish a full turn at a single empty slot from the zero path
            # between two empty slots at the same occurrence.
            gap_count=next_start+(!cut && j==length(indices) ? n : 0)-start-consumed
            if cut
                gap_count>=0 || throw(ArgumentError("overlapping weighted cut fragments"))
            else
                gap_count<0 && (gap_count+=n)
            end
            positions=[cut ? start+consumed+k : mod1(start+consumed+k,n) for k in 0:gap_count-1]
            edge_positions=vcat(positions,next_start)
            exterior=BigInt[_move_integer(weights[component][p]) for p in edge_positions]
            # An empty source fragment has no consumed incoming edge.  The
            # edge belongs to its preceding exterior path exactly once.
            if isempty(mapped_before[index])
                exterior[1]=0
            end
            exterior[1]-=last(rules0[index])
            incoming=terminal ? big(0) : first(rules0[next_index])
            exterior[end]-=incoming
            trace===nothing || push!(trace,(kind=:exterior,source_fragment=index,
                next_fragment=next_index,component=component,
                head_positions=edge_positions,
                source_weights=BigInt[_move_integer(weights[component][p]) for p in edge_positions],
                first_portion_omitted=isempty(mapped_before[index]),
                source_outgoing=last(rules0[index]),source_incoming=incoming,
                start_port=before.capacity[index],end_port=terminal ? _move_boundary_port(common,component,:out) : -next_index,
                word=_relabel_word(code[component][positions],old_to_new),portions=copy(exterior)))
            push!(segments,_WeightedMoveSegment(before.capacity[index],terminal ? _move_boundary_port(common,component,:out) : -next_index,
                  _relabel_word(code[component][positions],old_to_new),exterior))
        end
    end
    for index in 1:common
        trace===nothing || push!(trace,(kind=:target_fragment,index=index,
            source_location=candidate.locations[index],source_word=mapped_before[index],
            source_rule=rules0[index],word=mapped_after[index],portions=rules1[index],
            start_port=-index,end_port=after.capacity[index]))
        push!(segments,_WeightedMoveSegment(-index,after.capacity[index],mapped_after[index],rules1[index]))
    end
    paths=_move_segment_paths(segments,boundaries)
    words=Vector{Vector{Int}}(); output_weights=Vector{Vector{Float64}}()
    for path in paths
        word=Int[]; values=BigInt[]; pending=big(0)
        for current in path.indices
            segment=segments[current]
            if isempty(segment.word)
                pending+=only(segment.edges)
            else
                append!(word,segment.word)
                local_values=copy(segment.edges[1:end-1]); local_values[1]+=pending
                append!(values,local_values); pending=last(segment.edges)
            end
        end
        if path.incoming!==nothing
            push!(word,0)
            push!(values,pending)
        elseif isempty(values)
            pending==0 || throw(ArgumentError("nonzero weight on a vertex-free component cannot be represented"))
        else
            values[1]+=pending
        end
        push!(words,word)
        push!(output_weights,Float64.(_move_integer.(values)))
    end
    endpoint_trace=trace===nothing ? nothing : NamedTuple[]
    result,result_weights=_install_move_paths(code,weights,affected,removed,old_to_new,
        paths,words,output_weights;boundary_trace=endpoint_trace)
    trace===nothing || push!(trace,(kind=:boundary_transport,endpoints=endpoint_trace))
    for index in common+1:length(after.gauss)
        values=copy(rules1[index][1:end-1])
        isempty(values) ? (only(rules1[index])==0 || throw(ArgumentError("nonzero empty component weight"))) :
                         (values[1]+=last(rules1[index]))
        push!(result,copy(mapped_after[index]))
        push!(result_weights,Float64.(_move_integer.(values)))
    end
    result,result_weights,new_orientation
end

# Retain the unweighted converter's component order and basepoints.  Only a
# closed component may be shifted cyclically; a cut word has an anchored root.
function _align_move_weights(reference,transported,weights)
    remaining=collect(eachindex(transported)); result=Vector{Vector{Float64}}()
    for word in reference
        found=false
        for (j,index) in enumerate(remaining)
            length(word)==length(transported[index]) || continue
            if _is_cut_component(word) || _is_cut_component(transported[index])
                word==transported[index] || continue
                push!(result,copy(weights[index]))
                deleteat!(remaining,j); found=true; break
            end
            starts=isempty(word) ? [1] : _cyclic_occurrences(transported[index],word)
            isempty(starts) && continue
            push!(result,isempty(word) ? Float64[] : circshift(weights[index],1-first(starts)))
            deleteat!(remaining,j); found=true; break
        end
        found || throw(ArgumentError("weighted transport and Gauss reconnection disagree"))
    end
    isempty(remaining) || throw(ArgumentError("weighted transport has unmatched components"))
    result
end

function _apply_pattern_move!(v::VirtualLink,candidate::MoveCandidate;legacy_mp_labels=false,
                              prepared_weights=false,boundary_trace=nothing)
    _check_candidate_split(v,candidate)
    v.is_weighted && _check_move_weights(v)
    direct_02=candidate.type=="02" && candidate.direction==:forward ? first(candidate.locations) : nothing
    candidate=_normalized_02(candidate)
    before=something(candidate.before); after=something(candidate.after)
    raw_code,orientation=real_gauss_code(v)
    cut_flags=[_is_cut_component(component) for component in raw_code]
    code=raw_code
    weights=v.is_weighted ? deepcopy(v.weights) :
        [zeros(Float64,length(component)) for component in code]
    v.is_weighted && !(candidate.type in ("MP","02") || (candidate.type=="PS" && prepared_weights)) &&
        throw(ArgumentError("integer weight transport is not defined for $(candidate.type)"))
    old_to_new,after_vertices,new_orientation=_vertex_maps(
        orientation,candidate.vertices,before,after;
        role_ordered=candidate.type=="CP" || (candidate.type=="MP" && !legacy_mp_labels) ||
                     (candidate.type=="BMP" && candidate.parameter in ("0","1")))
    mapped_before=[_map_fragment(fragment,candidate.vertices) for fragment in before.gauss]
    mapped_after=[_map_fragment(fragment,after_vertices) for fragment in after.gauss]

    simple=!any(cut_flags) && length(before.gauss)==length(after.gauss) && before.capacity==after.capacity
    converted=simple ? _simple_pattern_rewrite(code,weights,candidate,mapped_before,mapped_after,
                                                old_to_new,new_orientation;weighted=false) : nothing
    converted===nothing &&
        (converted=_reconnect_pattern(code,weights,candidate,mapped_before,mapped_after,
                                      old_to_new,new_orientation;weighted=false,boundary_trace=boundary_trace))
    result,result_weights,result_orientation=converted
    if direct_02!==nothing
        component,start=direct_02
        result=deepcopy(code)
        result[component]=vcat(code[component][1:start],only(mapped_after),code[component][start+1:end])
    end
    target=VirtualLink()
    if v.is_weighted
        transported,transported_weights,_=_weighted_pattern_rewrite(code,weights,candidate,
            mapped_before,mapped_after,old_to_new,new_orientation)
        result_weights=_align_move_weights(result,transported,transported_weights)
        set_data!(target;odata=(result,result_orientation,result_weights))
        _check_move_weights(target)
    else
        set_data!(target;rgauss=result,orientation=result_orientation)
    end
    _commit_move!(v,target)
end

"Internal recomputable edge-occurrence evidence for integer move transport."
function _weighted_transport_detail(v::VirtualLink,candidate::MoveCandidate)
    _check_candidate_split(v,candidate)
    _check_move_weights(v)
    raw_code,orientation=real_gauss_code(v)
    code=raw_code
    if candidate.type=="H"
        amount=_move_integer(parse(BigFloat,candidate.parameter))
        vertex=abs(only(candidate.vertices)); incidences=Any[]
        for signed in (vertex,-vertex)
            for (component,word) in enumerate(code), head in findall(==(signed),word)
                push!(incidences,(component=component,head_position=head,delta=-amount))
                push!(incidences,(component=component,head_position=mod1(head+1,length(word)),delta=amount))
            end
        end
        return (version=1,kind=:H,incidences=incidences)
    end
    candidate=_normalized_02(candidate)
    before=candidate.before; after=candidate.after
    mapping,vertices,new_orientation=_vertex_maps(orientation,candidate.vertices,before,after;
        role_ordered=candidate.type in ("MP","CP"))
    trace=Any[]
    _weighted_pattern_rewrite(code,v.weights,candidate,
        [_map_fragment(g,candidate.vertices) for g in before.gauss],
        [_map_fragment(g,vertices) for g in after.gauss],mapping,new_orientation;trace=trace)
    (version=1,kind=:fragment_transport,vertex_mapping=sort!(collect(mapping);by=first),
     fragments=trace)
end

"""Validate one frozen PS witness in the closed-normal Julia subcategory.

Callers must use `Topology.validate_ps_move` explicitly; `movable` remains the
public applicability gate.  This validator checks that the frozen candidate is
an authentic current-source PS-R witness before applying it, then requires the
source and tentative target to be closed.  It never calls `movable`, and never
mutates the source. Integer-weighted input additionally uses the source-specific
PS rotation plan and exact source/target framing checks. Cut components still
require a separate boundary contract.

The returned named tuple is deliberately inspectable.  `source_conditions` and
`target_conditions` are `(C1,C2,C3)`, `target` is an independently constructed
 tentative link when application reaches that stage, and `error` is `nothing`
except for an expected tentative-application input failure.
"""
function _same_move_pattern(left::MovePattern,right::MovePattern)
    left.parameter==right.parameter && left.gauss==right.gauss &&
    left.orientation==right.orientation && left.capacity==right.capacity &&
    left.gauss_first==right.gauss_first && left.gauss_last==right.gauss_last &&
    left.vertex_count==right.vertex_count && isequal(left.weights,right.weights)
end

function validate_ps_move(v::VirtualLink,candidate::MoveCandidate;boundary=nothing)
    source_closed,source_conditions=is_closed(v)
    cut_source=any(_is_cut_component, v.gauss)
    invalid(reason,target=nothing,target_closed=false,
            target_conditions=(false,false,false);error=nothing)=
        (valid=false,reason=reason,source_closed=source_closed,
         source_conditions=source_conditions,target_closed=target_closed,
         target_conditions=target_conditions,target=target,error=error)
    _normalize_move_type(candidate.type)=="PS" || return invalid(:not_ps)
    if v.is_weighted
        try
            _check_candidate_split(v,candidate)
            target=_weighted_ps_target(v,candidate;boundary=boundary)
            target_closed,target_conditions=is_closed(target)
            return (valid=true,reason=:weighted_ps,source_closed=source_closed,
                source_conditions=source_conditions,target_closed=target_closed,
                target_conditions=target_conditions,target=target,error=nothing)
        catch error
            error isa Union{ArgumentError,DimensionMismatch,BoundsError} || rethrow()
            return invalid(:weighted_ps_inapplicable,error=error)
        end
    end
    isempty(candidate.edges) || return invalid(:edge_selected_unsupported)
    (candidate.before===nothing || candidate.after===nothing) &&
        return invalid(:malformed_candidate)
    candidate.direction in (:forward,:inverse) || return invalid(:source_mismatch)
    parameter=try parse(Int,candidate.parameter) catch; return invalid(:source_mismatch) end
    1<=parameter<=4 || return invalid(:source_mismatch)
    data=make_move_data()
    canonical_before=candidate.direction==:inverse ? data["PS-R"][parameter] : only(data["PS-L"])
    canonical_after=candidate.direction==:inverse ? only(data["PS-L"]) : data["PS-R"][parameter]
    _same_move_pattern(candidate.before,canonical_before) &&
        _same_move_pattern(candidate.after,canonical_after) ||
        return invalid(:source_mismatch)
    code,orientation=real_gauss_code(v)
    length(candidate.vertices)==canonical_before.vertex_count || return invalid(:source_mismatch)
    vertices=candidate.vertices
    all(vertex -> vertex>0 && vertex<=length(orientation),vertices) ||
        return invalid(:source_mismatch)
    length(unique(vertices))==length(vertices) || return invalid(:source_mismatch)
    candidate.before.orientation==orientation[vertices] || return invalid(:source_mismatch)
    length(candidate.locations)==length(candidate.before.gauss) ||
        return invalid(:source_mismatch)
    mapped_before=[_map_fragment(fragment,vertices) for fragment in candidate.before.gauss]
    for (fragment,(component,start)) in zip(mapped_before,candidate.locations)
        1<=component<=length(code) || return invalid(:source_mismatch)
        n=length(code[component]); 1<=start<=n || return invalid(:source_mismatch)
        (isempty(fragment) || start in _pattern_occurrences(code[component],fragment)) ||
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
    !source_closed && !cut_source && return invalid(:source_not_closed,target,target_closed,target_conditions)
    source_closed && !target_closed && return invalid(:target_not_closed,target,target_closed,target_conditions)
    source_closed && !ps_global_condition(v,candidate).valid &&
        return invalid(:ps_method_b_condition_failed,target,target_closed,target_conditions)
    (valid=true,reason=source_closed ? :closed_local_ps : :open_local_ps,
     source_closed=source_closed,source_conditions=source_conditions,target_closed=target_closed,
     target_conditions=target_conditions,target=target,error=nothing)
end

function _apply_direct_02!(v::VirtualLink,candidate::MoveCandidate)
    v.is_weighted && return _apply_pattern_move!(v,candidate)
    code,orientation=real_gauss_code(v)
    component,start=candidate.locations[1]
    new_vertices=length(orientation).+(1:2)
    fragment=[-new_vertices[1],-new_vertices[2],new_vertices[2],new_vertices[1]]
    word=code[component]
    code[component]=vcat(word[1:start],fragment,word[start+1:end])
    set_data!(v;rgauss=code,orientation=vcat(orientation,[1,-1]))
end

function _apply_bmp_two!(v::VirtualLink,candidate::MoveCandidate,parameter;boundary_trace=nothing)
    v.is_weighted && throw(ArgumentError("integer weight transport is not defined for BMP"))
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
    new_vertex=S*(length(orientation)+1)
    circle=parameter==1 ? [-values[1],-values[2],-new_vertex] :
                           [-values[2],-values[1],-new_vertex]
    new_orientation=vcat(orientation,0)
    circle_orientation=parameter==1 ? -[-1,1,1]*S : [-1,1,1]*S
    for (vertex,value) in zip((v1,v2,abs(new_vertex)),circle_orientation)
        new_orientation[vertex]=value
    end
    # Replace the adjacent passage, interchange the opposite passage ports,
    # and add the new closed circle.  The same capacity graph handles circles,
    # intervals, and any number of boundary endpoints.
    before=MovePattern("BMP",[[1,2],[-1],[-2]],orientation[[v1,v2]],[1,2,3])
    after=MovePattern("BMP",[[3],[2],[1],[-1,-2,-3]],circle_orientation,[1,3,2,4])
    locations=[locate(values[1]),locate(-values[1]),locate(-values[2])]
    local_candidate=MoveCandidate("BMP",string(parameter),values,candidate.edges,:forward,
                                  before,after,locations)
    mapped_before=[copy(values),[-values[1]],[-values[2]]]
    mapped_after=[[new_vertex],[values[2]],[values[1]],circle]
    weights=[zeros(Float64,length(word)) for word in code]
    mapping=Dict(index=>index for index in eachindex(orientation))
    result,_,_=_reconnect_pattern(code,weights,local_candidate,mapped_before,mapped_after,
                                  mapping,new_orientation;boundary_trace=boundary_trace)
    set_data!(v;rgauss=result,orientation=new_orientation)
end

function _apply_b02_two!(v::VirtualLink,candidate::MoveCandidate;boundary_trace=nothing)
    v.is_weighted && throw(ArgumentError("integer weight transport is not defined for B02"))
    code,orientation=real_gauss_code(v)
    values=copy(candidate.vertices)
    nv=length(orientation)
    mapping=Dict(old=>new for (new,old) in enumerate(setdiff(1:nv,values)))
    circle_index=findfirst(word -> length(word)==2 &&
        Set(abs.(word))==Set(values) && sign(word[1])==sign(word[2]),code)
    circle_index===nothing && throw(ArgumentError("B02 requires its actual closed two-vertex circle"))
    circle=code[circle_index]
    values=-circle
    locate(value)=begin
        locations=[(component,position) for (component,word) in enumerate(code)
                   for position in findall(==(value),word)]
        length(locations)==1 || throw(ArgumentError("B02 vertex occurrence is ambiguous"))
        only(locations)
    end
    new_orientation=copy(orientation)
    deleteat!(new_orientation,sort(abs.(values)))
    before=MovePattern("B02",[[1],[2],[-1,-2]],orientation[abs.(values)],[2,1,3])
    after=MovePattern("B02",[Int[],Int[]],Int[],[1,2])
    locations=[locate(values[1]),locate(values[2]),(circle_index,1)]
    local_candidate=MoveCandidate("B02",candidate.parameter,candidate.vertices,candidate.edges,
                                  :inverse,before,after,locations)
    weights=[zeros(Float64,length(word)) for word in code]
    result,_,_=_reconnect_pattern(code,weights,local_candidate,[[values[1]],[values[2]],circle],
        [Int[],Int[]],mapping,new_orientation;boundary_trace=boundary_trace)
    set_data!(v;rgauss=result,orientation=new_orientation)
end

function _apply_hmove!(v::VirtualLink,vertex::Int,amount)
    v.is_weighted || throw(ArgumentError("H move requires a weighted VirtualLink"))
    _check_move_weights(v)
    amount=_move_integer(amount)
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
    adjustments=Dict{NTuple{2,Int},BigInt}()
    for (pair,delta) in zip(pairs,(-1,1,-1,1))
        adjustments[pair]=get(adjustments,pair,big(0))+delta*amount
    end
    table=real_edge_table(v)
    by_strand=[Float64[] for _ in code]
    for edge in table
        push!(by_strand[edge.strand],Float64(_move_integer(BigInt(edge.weight)+get(adjustments,edge.crossing,big(0)))))
    end
    set_weight!(v,by_strand)
end

"""Apply one MATLAB-style move after strict applicability selection."""
function _move_unchecked!(v::VirtualLink,type;param=nothing,vertices=Int[],edges=nothing,
               edge_occurrences=nothing,strict=true,boundary=nothing)
    move_type=_normalize_move_type(type)
    if weight_mode(v)==:z2
        move_type=="CP" || throw(ArgumentError("Z2 move transport is only implemented for CP"))
        _spin_edge_weights(v)
    elseif v.is_weighted
        _check_move_weights(v)
        move_type in ("MP","H","02","PS") ||
            throw(ArgumentError("integer weight transport is not defined for $move_type"))
    end
    if move_type in ("R1","R2","R3")
        throw(ArgumentError("$move_type move is not implemented in the MATLAB source"))
    elseif move_type=="H"
        param===nothing && throw(ArgumentError("H move requires a weight parameter"))
        candidates=movable(v,move_type;param=param,vertices=vertices,edges=edges,strict=strict,
                           boundary=boundary)
        isempty(candidates) && throw(ArgumentError("No valid moves found."))
        return _apply_hmove!(v,abs(only(candidates).vertices[1]),param)
    end
    candidates=movable(v,move_type;param=param,vertices=vertices,edges=edges,
                       edge_occurrences=edge_occurrences,strict=strict,boundary=boundary)
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
        _apply_candidate!(v,candidate;boundary=boundary)
    end
end

function move!(v::VirtualLink,type;kwargs...)
    target=deepcopy(v)
    _move_unchecked!(target,type;kwargs...)
    _commit_move!(v,target)
end

function move!(v::VirtualLink,type,arg::NamedTuple)
    param=haskey(arg,:param) ? arg.param : haskey(arg,:p) ? arg.p : nothing
    vertices=haskey(arg,:vertices) ? arg.vertices : haskey(arg,:v) ? arg.v : Int[]
    edges=haskey(arg,:edges) ? arg.edges : haskey(arg,:e) ? arg.e : nothing
    strict=haskey(arg,:strict) ? arg.strict : true
    occurrences=haskey(arg,:edge_occurrences) ? arg.edge_occurrences : nothing
    boundary=haskey(arg,:boundary) ? arg.boundary : nothing
    move!(v,type;param=param,vertices=vertices,edges=edges,edge_occurrences=occurrences,
          strict=strict,boundary=boundary)
end

"Apply this exact candidate, without selecting its site again or canonicalizing."
function _apply_candidate!(v::VirtualLink,candidate::MoveCandidate;boundary=nothing)
    _check_candidate_split(v,candidate)
    if weight_mode(v)==:z2
        candidate.type=="CP" || throw(ArgumentError("Z2 move transport is only implemented for CP"))
        return _commit_move!(v,spin_cp_move(v,candidate;boundary=boundary))
    elseif v.is_weighted && candidate.type=="PS"
        return _commit_move!(v,_weighted_ps_target(v,candidate;boundary=boundary))
    elseif candidate.type=="H"
        amount=tryparse(BigInt,candidate.parameter)
        amount===nothing && (amount=tryparse(BigFloat,candidate.parameter))
        amount===nothing && throw(ArgumentError("H move requires an integer parameter"))
        return _apply_hmove!(v,abs(only(candidate.vertices)),amount)
    elseif candidate.type=="02" && candidate.direction==:forward
        return _apply_direct_02!(v,candidate)
    elseif candidate.type=="BMP" && candidate.before===nothing
        parameter=tryparse(Int,candidate.parameter)
        return _apply_bmp_two!(v,candidate,parameter)
    elseif candidate.type=="B02" && candidate.before===nothing
        return _apply_b02_two!(v,candidate)
    end
    candidate.type=="PS" && !validate_ps_move(v,candidate).valid &&
        throw(ArgumentError("PS candidate fails closed-normal applicability"))
    _apply_pattern_move!(v,candidate)
end

function _same_candidate(left::MoveCandidate,right::MoveCandidate)
    left.type==right.type && left.parameter==right.parameter &&
    left.ps_split_k==right.ps_split_k &&
    left.vertices==right.vertices && left.edges==right.edges &&
    left.direction==right.direction && left.locations==right.locations &&
    ((left.before===nothing && right.before===nothing) ||
     (left.before!==nothing && right.before!==nothing && _same_move_pattern(left.before,right.before))) &&
    ((left.after===nothing && right.after===nothing) ||
     (left.after!==nothing && right.after!==nothing && _same_move_pattern(left.after,right.after)))
end

function move!(v::VirtualLink,candidate::MoveCandidate;boundary=nothing)
    _check_candidate_split(v,candidate)
    if weight_mode(v)==:z2
        candidate.type=="CP" || throw(ArgumentError("Z2 move transport is only implemented for CP"))
        _spin_edge_weights(v)
    else
        v.is_weighted && _check_move_weights(v)
    end
    # Frozen candidates are untrusted input: authenticate templates, orientation,
    # vertex roles and exact occurrences against the current geometric source.
    geometric_source=v.is_weighted && candidate.type=="PS" ? _ps_unweighted(v) : v
    options=if candidate.type in ("PS","S") && candidate.direction==:forward
        _movable_geometry(geometric_source,candidate.type;param=candidate.parameter,
            edge_occurrences=candidate.locations,strict=true)
    else
        _movable_geometry(geometric_source,candidate.type;param=candidate.parameter,
            vertices=candidate.vertices,edges=candidate.edges,strict=true)
    end
    # Geometric enumeration does not choose the optional weighted split. Its
    # validity is checked by the source-specific PS plan after authentication.
    geometric_candidate=candidate.ps_split_k===nothing ? candidate :
        MoveCandidate(candidate.type,candidate.parameter,candidate.vertices,candidate.edges,
            candidate.direction,candidate.before,candidate.after,candidate.locations)
    any(option -> _same_candidate(option,geometric_candidate),options) ||
        throw(ArgumentError("move candidate does not match the current source"))
    target=deepcopy(v)
    _apply_candidate!(target,candidate;boundary=boundary)
    _commit_move!(v,target)
end
move(v::VirtualLink,candidate::MoveCandidate;boundary=nothing)=
    move!(deepcopy(v),candidate;boundary=boundary)

"""Apply a frozen move on a copy and expose its actual boundary endpoint map.

Each row of `endpoints` maps one source `(strand,side)` to one target
`(strand,side)`.  In particular, a target out-end need not descend from the
same source strand as its in-end.  Endpoint lanes retain their lane number.
"""
function move_boundary_transport(source::VirtualLink,candidate::MoveCandidate;boundary=nothing)
    target=move(source,candidate;boundary=boundary)
    trace=NamedTuple[]
    probe=VirtualLink()
    code,orientation=real_gauss_code(source)
    set_data!(probe;rgauss=deepcopy(code),orientation=copy(orientation))
    if candidate.type=="H" || (candidate.type=="02" && candidate.direction==:forward)
        for strand in eachindex(code), side in (:in,:out)
            _is_cut_component(code[strand]) || continue
            endpoint=(strand=strand,side=side)
            push!(trace,(source=endpoint,target=endpoint))
        end
    elseif candidate.type=="BMP" && candidate.before===nothing
        _apply_bmp_two!(probe,candidate,parse(Int,candidate.parameter);boundary_trace=trace)
    elseif candidate.type=="B02" && candidate.before===nothing
        _apply_b02_two!(probe,candidate;boundary_trace=trace)
    else
        _apply_pattern_move!(probe,candidate;boundary_trace=trace)
    end
    if !(candidate.type=="H" || (candidate.type=="02" && candidate.direction==:forward))
        probe.gauss==target.gauss && probe.orientation==target.orientation ||
            throw(ArgumentError("boundary trace and executed move disagree"))
    end
    length(trace)==2count(_is_cut_component,code) || throw(ArgumentError("move lost a source boundary endpoint"))
    length(unique(row.source for row in trace))==length(trace) || throw(ArgumentError("duplicate source boundary endpoint"))
    length(unique(row.target for row in trace))==length(trace) || throw(ArgumentError("duplicate target boundary endpoint"))
    length(trace)==2count(_is_cut_component,target.gauss) || throw(ArgumentError("move changed the boundary endpoint count"))
    (target=target,endpoints=trace)
end

function move_boundary_transport(source::VirtualLink,type;boundary=nothing,kwargs...)
    candidates=movable(source,type;boundary=boundary,kwargs...)
    length(candidates)==1 || throw(ArgumentError("boundary transport requires one explicit move candidate"))
    move_boundary_transport(source,only(candidates);boundary=boundary)
end

"""Pure-value wrapper around MATLAB-compatible `move!`."""
function move(v::VirtualLink,type;kwargs...)
    result=deepcopy(v); move!(result,type;kwargs...); result
end
move(v::VirtualLink,type,arg::NamedTuple)=move!(deepcopy(v),type,arg)

"""Compatibility wrapper for the first Julia migration API."""
mp_move!(v::VirtualLink,parameter::AbstractString,vertices)=
    move!(v,"MP";param=parameter,vertices=vertices,strict=true)
