"""Gauss-code representation of a virtual link or o-graph.

`orientation[i]` has two distinct roles that must not be conflated:

- `+1` / `-1`: crossing `i` is real, with the corresponding signed local
  crossing convention;
- `0`: crossing `i` is virtual. It has no over/under topological information,
  is removed by `real_gauss_code`, and is not drawn as a real vertex.

A virtual crossing still needs a cyclic order when converted to PD slots and
routed geometrically. For that representation-only local frame, this library
uses the same slot convention as `+1`, matching Sage/MATLAB. This does **not**
turn the virtual crossing into a positive real crossing; the stored value
remains `0` and continues to control topology and rendering.
"""
mutable struct VirtualLink
    gauss::Vector{Vector{Int}}
    orientation::Vector{Int}
    weights::Vector{Vector{Float64}}
    real_only::Bool
    is_weighted::Bool
    strand_types::Vector{Char}
    weight_mode::Symbol
end

struct RealEdge
    id::Int
    crossing::NTuple{2,Int}
    weight::Float64
    strand::Int
    arc::Vector{Int}
end

function _check_weight_mode(mode)
    mode in (:integer,:z2) || throw(ArgumentError("weight mode must be :integer or :z2"))
    mode
end
VirtualLink(;weight_mode=:integer) = VirtualLink(Vector{Vector{Int}}(), Int[], Vector{Vector{Float64}}(), false, false, Char[], _check_weight_mode(weight_mode))
# Preserve the historical positional constructor as integral mode.
VirtualLink(g,o,w,r,i,s) = VirtualLink(g,o,w,r,i,s,:integer)
weight_mode(v::VirtualLink)=v.weight_mode

function set_weight_mode!(v::VirtualLink, mode::Symbol)
    _check_weight_mode(mode)
    mode==v.weight_mode && return v
    v.is_weighted && mode==:integer && throw(ArgumentError("Z2 weights have no canonical integer lift; assign explicit integer weights with set_weight!(...;mode=:integer)"))
    v.is_weighted ? set_weight!(v,v.weights;mode=mode) : (v.weight_mode=mode; v)
end

"""Return the real-Gauss input for the lens space `L(s,t)`.

This is the `RGauss` family used in MATLAB `VLExample.mlx` and
`C260220LensSpaces.mlx`. Crossing labels are reduced modulo `s`; `s` and `t`
must be coprime.
"""
function lens_space_rgauss(s::Integer,t::Integer)
    s>=2 || throw(ArgumentError("s must be at least 2"))
    gcd(s,t)==1 || throw(ArgumentError("s and t must be coprime"))
    code=vcat(-mod.((s:-1:1).*t.-1,s).-1,1:s)
    (rgauss=[Int.(code)],orientation=ones(Int,s))
end

"""Construct the real-Gauss `VirtualLink` presentation of `L(s,t)`."""
function lens_space_link(s::Integer,t::Integer)
    input=lens_space_rgauss(s,t)
    link=VirtualLink()
    set_data!(link;rgauss=input.rgauss,orientation=input.orientation)
end

function _normalize_gauss(code)
    code isa AbstractVector{<:Integer} && return [Int.(code)]
    [Int.(c) for c in code]
end

_is_cut_component(component::AbstractVector{<:Integer}) =
    !isempty(component) && component[end] == 0

function _component_visits(component::AbstractVector{<:Integer})
    _is_cut_component(component) ? Int.(component[begin:end-1]) : Int.(component)
end

function _validate_cut_markers(code::Vector{Vector{Int}})
    for (strand,component) in enumerate(code)
        length(component) > 1 && any(==(0),@view component[1:end-1]) &&
            throw(ArgumentError("cut marker 0 is allowed only at the end of strand $strand"))
    end
    true
end

function _validate_cut_flags(flags, strand_count::Integer)
    flags isa AbstractVector || throw(ArgumentError("cut flags must be a vector of Bool"))
    all(flag -> flag isa Bool, flags) ||
        throw(ArgumentError("cut flags must contain actual Bool values"))
    length(flags) == strand_count ||
        throw(DimensionMismatch("one cut flag per Gauss strand is required"))
    Bool[flags...]
end

"""Return one cut flag per strand, derived from its trailing zero marker."""
function strand_cuts(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    Bool[_is_cut_component(component) for component in v.gauss]
end

function _strand_arc_count(v::VirtualLink, strand::Integer)
    visits = _component_visits(v.gauss[strand])
    real_count = count(entry -> entry != 0 &&
                              abs(entry) <= length(v.orientation) &&
                              v.orientation[abs(entry)] != 0, visits)
    real_count + (_is_cut_component(v.gauss[strand]) ? 1 : 0)
end

"""Set trailing-zero cut markers on an unweighted link.

Existing arc-vector prefixes are retained while the derived vector shape is
resized. A weighted link must be cleared and assigned explicit new arc weights
before its cut status is changed.
"""
function set_cuts!(v::VirtualLink, flags)
    new_flags = _validate_cut_flags(flags, length(v.gauss))
    _validate_cut_markers(v.gauss)
    old_flags = Bool[_is_cut_component(component) for component in v.gauss]
    old_flags == new_flags && return v
    v.is_weighted &&
        throw(ArgumentError("changing cuts on a weighted link requires clearing weights and assigning explicit open-arc weights"))
    old_weights = deepcopy(v.weights)
    for (strand,cut) in enumerate(new_flags)
        component = v.gauss[strand]
        currently_cut = _is_cut_component(component)
        if cut && !currently_cut
            push!(component, 0)
        elseif !cut && currently_cut
            pop!(component)
        end
    end
    v.weights = [begin
        count_ = _strand_arc_count(v, strand)
        values = zeros(Float64, count_)
        old = strand <= length(old_weights) ? old_weights[strand] : Float64[]
        copy_count = min(length(old), count_)
        copy_count > 0 && (values[1:copy_count] = Float64.(old[1:copy_count]))
        values
    end for strand in eachindex(v.gauss)]
    v
end

function _validate_gauss(code::Vector{Vector{Int}}, ori::Vector{Int}; real_only=false)
    _validate_cut_markers(code)
    labels = filter(!=(0), reduce(vcat, code; init=Int[]))
    isempty(labels) && return true
    maximum(abs.(labels)) <= length(ori) || throw(ArgumentError("orientation is missing a crossing"))
    for i in 1:length(ori)
        count(==(i), abs.(labels)) == 2 || throw(ArgumentError("crossing $i must occur exactly twice"))
        real_only && ori[i] == 0 && throw(ArgumentError("real Gauss code cannot contain virtual orientation"))
    end
    true
end

"""Set Gauss, real-Gauss, O-data, PD-code, or tabular-knot input.

A `pd` matrix does not round-trip: it is rebuilt through `_pd_to_gauss` and
then re-derived by `pd_code`, which may relabel every edge. Per-edge data
carried alongside a PD matrix -- bending numbers, replay positions -- must be
indexed to what `pd_code` returns. See `PLOT_PIPELINE.md` IS4.
"""
function set_data!(v::VirtualLink; gauss=nothing, rgauss=nothing, orientation=nothing,
                   odata=nothing, pd=nothing, table=nothing, cut=Bool[], strand_types=nothing)
    supplied_weights = odata !== nothing && length(odata) == 3
    if odata !== nothing
        gauss, orientation = odata[1], odata[2]
        length(odata) == 3 && (v.weights = [Float64.(w) for w in odata[3]]; v.is_weighted = true)
        rgauss = gauss; gauss = nothing
    elseif table !== nothing
        tuple(table...) == (3,1) || throw(ArgumentError("only the trefoil table example [3,1] is defined without SageMath"))
        pd = [1 2 4 3; 3 4 6 5; 5 6 2 1]
        orientation = [-1,-1,-1]
    end
    if pd !== nothing
        if isempty(pd)
            v.gauss = [Int[]]; v.orientation = Int[]; v.real_only = false
        else
            v.orientation = orientation === nothing ? fill(1, size(pd,1)) : Int.(orientation)
            length(v.orientation) == size(pd,1) ||
                throw(DimensionMismatch("one orientation value per PD crossing is required"))
            all(sign -> sign in (-1,0,1), v.orientation) ||
                throw(ArgumentError("PD orientation values must be -1, 0, or 1"))
            v.gauss = _pd_to_gauss(Int.(pd); orientation=v.orientation)
            v.real_only = false
        end
    elseif gauss !== nothing || rgauss !== nothing
        code = _normalize_gauss(gauss === nothing ? rgauss : gauss)
        ori = orientation === nothing ? throw(ArgumentError("orientation is required")) : Int.(orientation)
        _validate_gauss(code, ori; real_only=rgauss !== nothing)
        v.gauss, v.orientation, v.real_only = code, ori, rgauss !== nothing
    else
        throw(ArgumentError("one input representation is required"))
    end
    v.is_weighted = false
    isempty(cut) || set_cuts!(v, cut)
    v.strand_types = strand_types === nothing ? fill('o', length(v.gauss)) : collect(strand_types)
    if supplied_weights
        v.is_weighted = false
        set_weight!(v,odata[3])
    else
        v.is_weighted = false
        v.weights = [zeros(Float64,_strand_arc_count(v,strand))
                     for strand in eachindex(v.gauss)]
    end
    v
end

function _pd_to_gauss(pd::Matrix{Int}; orientation=nothing)
    # PD edges form oriented components. Crossing sign is encoded by under/over slots.
    signs = orientation === nothing ? fill(1, size(pd,1)) : Int.(orientation)
    length(signs) == size(pd,1) ||
        throw(DimensionMismatch("one orientation value per PD crossing is required"))
    next = Dict{Int,Tuple{Int,Int}}()
    for (i,row) in enumerate(eachrow(pd))
        next[row[1]] = (row[3], -i)
        # The over-strand transition follows the local PD frame.  Positive
        # (and virtual) rows enter at slot 4 and leave at slot 2; negative
        # rows enter at slot 2 and leave at slot 4.  Keeping this directed
        # incidence avoids self-crossing edge IDs overwriting the transition.
        if signs[i] == -1
            next[row[2]] = (row[4], i)
        else
            next[row[4]] = (row[2], i)
        end
    end
    result = Vector{Vector{Int}}(); unseen = Set(keys(next))
    outgoing = Set(first(value) for value in values(next))
    starts = sort!(collect(setdiff(Set(keys(next)), outgoing)))
    function walk!(start, open)
        edge = start; code = Int[]
        while edge in unseen
            delete!(unseen, edge)
            edge, crossing = next[edge]
            push!(code, crossing)
        end
        open && push!(code, 0)
        push!(result, code)
    end
    for start in starts
        walk!(start, true)
    end
    while !isempty(unseen)
        walk!(minimum(unseen), false)
    end
    result
end

function gauss_code(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    deepcopy(v.gauss)
end
function real_gauss_code(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    real = findall(!=(0), v.orientation); mapping = Dict(x=>i for (i,x) in enumerate(real))
    code = Vector{Vector{Int}}()
    for component in v.gauss
        visits = _component_visits(component)
        filtered = [sign(x)*mapping[abs(x)] for x in visits
                    if x != 0 && haskey(mapping, abs(x))]
        _is_cut_component(component) && push!(filtered, 0)
        push!(code, filtered)
    end
    code, v.orientation[real]
end

"""Return the stable fields of MATLAB `REdgeTable` in traversal order."""
function real_edge_table(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    real_indices=findall(!=(0),v.orientation)
    mapping=Dict(index=>mapped for (mapped,index) in enumerate(real_indices))
    result=RealEdge[]
    edge_id=0
    for (strand,component0) in enumerate(v.gauss)
        component=_component_visits(component0)
        cut=_is_cut_component(component0)
        real_positions=findall(position -> haskey(mapping,abs(component[position])),eachindex(component))
        weights=strand <= length(v.weights) ? v.weights[strand] : Float64[]
        convert_crossing(x)=x == 0 ? 0 : sign(x)*mapping[abs(x)]
        push_edge!(crossing,arc,weight_index)=begin
            edge_id += 1
            weight=weight_index <= length(weights) ? weights[weight_index] : NaN
            push!(result,RealEdge(edge_id,crossing,weight,strand,Int.(arc)))
        end
        if cut
            if isempty(real_positions)
                push_edge!((0,0),component,1)
            else
                first_position=first(real_positions)
                push_edge!((0,convert_crossing(component[first_position])),
                           component[1:first_position],1)
                for weight_index in 2:length(real_positions)
                    previous_position=real_positions[weight_index-1]
                    target_position=real_positions[weight_index]
                    push_edge!((convert_crossing(component[previous_position]),
                                convert_crossing(component[target_position])),
                               component[previous_position:target_position],weight_index)
                end
                last_position=last(real_positions)
                push_edge!((convert_crossing(component[last_position]),0),
                           component[last_position:end],length(real_positions)+1)
            end
        elseif !isempty(real_positions)
            for (weight_index,target_position) in enumerate(real_positions)
                previous_position=real_positions[mod1(weight_index-1,length(real_positions))]
                arc=Int[component[previous_position]]
                cursor=previous_position
                while cursor != target_position
                    cursor=mod1(cursor+1,length(component))
                    push!(arc,component[cursor])
                end
                crossing=(convert_crossing(first(arc)),convert_crossing(last(arc)))
                push_edge!(crossing,arc,weight_index)
            end
        end
    end
    result
end

function head_map(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    n = length(v.orientation); hm = fill(0, n, 2)
    for component in v.gauss
        clean = filter(!=(0), component); isempty(clean) && continue
        closed = isempty(component) || last(component) != 0
        stop = closed ? length(clean) : length(clean)-1
        for i in 1:max(stop,0)
            current = clean[i]; successor = clean[mod1(i+1, length(clean))]
            hm[abs(current), current > 0 ? 2 : 1] = successor
        end
    end
    hm
end

function pd_code(v::VirtualLink)
    # Sage 10.1+ convention: start at the entering under-strand and walk
    # counter-clockwise. This is the convention used by the MATLAB oracle.
    _validate_cut_markers(v.gauss)
    code = Vector{Vector{Int}}()
    cuts = Bool[]
    for component in v.gauss
        visits = _component_visits(component)
        isempty(visits) && continue
        push!(code,visits); push!(cuts,_is_cut_component(component))
    end
    isempty(code) && return zeros(Int,0,4)

    edge_at_visit = Dict{Int,NTuple{2,Int}}()
    edge = 0
    for (component,cut) in zip(code,cuts)
        first_edge = edge+1
        cut && (edge = first_edge)
        for crossing in component
            edge += 1
            edge_at_visit[crossing] = (edge,edge+1)
        end
        if !cut
            last_crossing = last(component)
            incoming,outgoing = edge_at_visit[last_crossing]
            edge_at_visit[last_crossing] = (incoming,first_edge)
        else
            edge += 1 # reserve the outgoing boundary edge
        end
    end

    pd = zeros(Int,length(v.orientation),4)
    for crossing in eachindex(v.orientation)
        under = edge_at_visit[-crossing]
        over = edge_at_visit[crossing]
        if v.orientation[crossing] == -1
            pd[crossing,:] = [under[1],over[1],under[2],over[2]]
        else
            # `0` is topologically virtual, but a PD row still needs a cyclic
            # slot order. Sage/MATLAB use the +1 local frame for that purely
            # combinatorial choice; virtual status remains in `orientation`.
            pd[crossing,:] = [under[1],over[2],under[2],over[1]]
        end
    end
    pd
end

"""Return edge => crossing-index maps for tails and heads of a PD code."""
function _pd_has_boundary(pd::AbstractMatrix{<:Integer})
    isempty(pd) && return false
    values=vec(Int.(pd))
    any(count(==(edge),values)==1 for edge in unique(values))
end

function edge_directions(pd::AbstractMatrix{<:Integer})
    _pd_has_boundary(pd) &&
        throw(ArgumentError("edge directions require a closed PD code; boundary edge IDs are present"))
    tails = Dict{Int,Int}()
    heads = Dict{Int,Int}()
    rows = [Int.(collect(row)) for row in eachrow(pd)]
    for crossing_index in eachindex(rows)
        crossing = rows[crossing_index]
        tails[crossing[3]] = crossing_index
        edge = crossing[3]
        current = crossing_index
        while !haskey(heads,edge)
            next_crossings = [i for i in eachindex(rows) if i != current && edge in rows[i]]
            if isempty(next_crossings)
                heads[edge] = current
                tails[edge] = current
                current_row = rows[current]
                edge = current_row[1] == edge ? current_row[3] :
                       current_row[4] == edge ? current_row[2] : current_row[4]
            else
                next_crossing = first(next_crossings)
                heads[edge] = next_crossing
                tails[edge] = current
                current = next_crossing
                current_row = rows[current]
                slot = findfirst(==(edge),current_row)
                edge = current_row[mod1(slot+2,4)]
            end
        end
    end

    unassigned = setdiff(Set(vec(Int.(pd))),Set(keys(tails)))
    while !isempty(unassigned)
        edge = pop!(unassigned)
        current = findfirst(row -> edge in row,rows)
        while !haskey(heads,edge)
            tails[edge] = current
            next_crossing = findfirst(i -> i != current && edge in rows[i],eachindex(rows))
            next_crossing === nothing && throw(ArgumentError("invalid PD edge $edge"))
            heads[edge] = next_crossing
            current = next_crossing
            current_row = rows[current]
            slot = findfirst(==(edge),current_row)
            edge = current_row[mod1(slot+2,4)]
            delete!(unassigned,edge)
        end
    end
    tails,heads
end

"""Return signed, left-turning region boundaries of a PD code.

Output and ordering match Sage `Link.regions()`, which the minimum-bend MILP
consumes after a stable sort by length. Verified on the hopf link; see
`PLOT_PIPELINE.md` section 4.
"""
function regions(pd::AbstractMatrix{<:Integer})
    _pd_has_boundary(pd) &&
        throw(ArgumentError("regions require a closed PD code; boundary edge IDs are present"))
    rows = [Int.(collect(row)) for row in eachrow(pd)]
    isempty(rows) && return Vector{Vector{Int}}()
    if length(rows) == 1
        crossing = rows[1]
        return crossing[1] == crossing[4] ?
            [[-crossing[3]],[crossing[1]],[crossing[3],-crossing[1]]] :
            [[crossing[3]],[-crossing[1]],[-crossing[3],crossing[1]]]
    end

    tails,heads = edge_directions(pd)
    unsigned = Set(vec(Int.(pd)))
    loops = sort([edge for edge in unsigned if heads[edge] == tails[edge]])
    available = union(unsigned,Set(-edge for edge in unsigned))
    result = Vector{Vector{Int}}()
    for edge in loops
        crossing = rows[heads[edge]]
        push!(result,crossing[4] == edge ? [edge] : [-edge])
        delete!(available,edge)
        delete!(available,-edge)
    end

    while !isempty(available)
        edge = maximum(available)
        delete!(available,edge)
        region = Int[]
        while !(edge in region)
            push!(region,edge)
            crossing_index = edge > 0 ? heads[edge] : tails[-edge]
            crossing = rows[crossing_index]
            slot = findfirst(==(abs(edge)),crossing)
            next_edge = crossing[mod1(slot-1,4)]
            if [next_edge] in result
                push!(region,-next_edge)
                next_edge = crossing[mod1(slot+1,4)]
            elseif [-next_edge] in result
                push!(region,next_edge)
                next_edge = crossing[mod1(slot+1,4)]
            end
            edge = tails[next_edge] == crossing_index ? next_edge : -next_edge
            delete!(available,edge)
        end
        push!(result,region)
    end
    result
end

function hmove_matrix(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    head=zeros(Int,vertex_count,2)
    for (strand,component0) in enumerate(code)
        component=_component_visits(component0)
        isempty(component) && continue
        cut=_is_cut_component(v.gauss[strand])
        for index in eachindex(component)
            vertex=component[index]
            successor = cut && index == lastindex(component) ? 0 :
                component[mod1(index+1,length(component))]
            head[abs(vertex),vertex > 0 ? 2 : 1]=successor
        end
    end
    predecessor(vertex)=begin
        negative=findfirst(==(vertex),view(head,:,1))
        positive=findfirst(==(vertex),view(head,:,2))
        something(positive,0)-something(negative,0)
    end
    vertices=collect(1:vertex_count)
    signed=hcat(-vertices,vertices)
    previous=predecessor.(signed)
    pairs=vcat(hcat(vec(previous),vcat(-vertices,vertices)),
               hcat(vcat(-vertices,vertices),vec(head)))
    crossings=[edge.crossing for edge in real_edge_table(v)]
    indices=[something(findfirst(==(Tuple(pair)),crossings),0) for pair in eachrow(pairs)]
    any(==(0),indices) && throw(ArgumentError("real edge table is incomplete"))
    indices=reshape(indices,vertex_count,4)
    matrix=zeros(Int,length(crossings),vertex_count)
    for vertex in vertices
        matrix[indices[vertex,1],vertex]=-1
        matrix[indices[vertex,2],vertex]=-1
        matrix[indices[vertex,3],vertex]+=1
        matrix[indices[vertex,4],vertex]+=1
    end
    matrix
end

function _disk_component_path(nodes::Vector{Int},
                              edges::Vector{Tuple{Int,Int,Int}},
                              is_cycle::Bool)
    isempty(nodes) && return (Int[],Int[])
    outgoing=Dict{Int,Vector{Tuple{Int,Int}}}()
    incoming=Dict(node=>0 for node in nodes)
    for (source,target,coefficient) in edges
        push!(get!(outgoing,source,Tuple{Int,Int}[]),(target,coefficient))
        incoming[target]=get(incoming,target,0)+1
    end
    all(length(get(outgoing,node,Tuple{Int,Int}[])) <= 1 for node in nodes) ||
        throw(ArgumentError("disk component has branching outgoing edges"))
    all(get(incoming,node,0) <= 1 for node in nodes) ||
        throw(ArgumentError("disk component has branching incoming edges"))

    if is_cycle
        all(length(get(outgoing,node,Tuple{Int,Int}[])) == 1 for node in nodes) ||
            throw(ArgumentError("disk cycle has a missing outgoing edge"))
        all(get(incoming,node,0) == 1 for node in nodes) ||
            throw(ArgumentError("disk cycle has a missing incoming edge"))
        start=minimum(nodes)
    else
        starts=sort([node for node in nodes if get(incoming,node,0) == 0])
        length(starts) == 1 ||
            throw(ArgumentError("open disk component must have one initial node"))
        start=only(starts)
    end

    path=Int[]; edge_signs=Int[]; current=start
    while true
        push!(path,current)
        next_edges=get(outgoing,current,Tuple{Int,Int}[])
        isempty(next_edges) && break
        target,coefficient=only(next_edges)
        push!(edge_signs,coefficient)
        current=target
        is_cycle && current == start && break
        current in path && throw(ArgumentError("disk component repeats before closing"))
    end
    length(path) == length(nodes) && length(edge_signs) == length(edges) ||
        throw(ArgumentError("disk component path does not cover its graph"))
    path,edge_signs
end

"""Return the MATLAB `getDiskTable` data for a virtual link.

The MATLAB columns `ID`, `paths`, `edges_s`, `Ndots`, `isCycle`, `delta`, and
`cp` are retained verbatim as NamedTuple fields.  `cp2 == 2cp` is the exact
integer-scaled form used by Julia's weight solver.  For open paths it uses the
implementation convention `2*isCycle - Ndots`; boundary lane identities are
reported separately and are not folded into a closed framing claim.
"""
function _disk_node_dictionary()
    dp=zeros(Int,12); dm=zeros(Int,12)
    dp[[1,2,4,5,9,12,7,6,8,11,10,3]]=repeat(1:6,2)
    dm[[1,2,4,5,9,12,7,10,3,11,6,8]]=repeat(1:6,2)
    Dict(1=>dp,-1=>dm)
end

function _empty_disk_table()
    (ID=Int[],paths=Vector{Vector{Int}}(),edges_s=Vector{Vector{Int}}(),
     Ndots=Int[],isCycle=Bool[],delta=zeros(Int,0,0),
     cp=Rational{Int}[],cp2=Int[],path_lengths=Int[],
     boundary_start=Union{Nothing,NamedTuple}[],
     boundary_end=Union{Nothing,NamedTuple}[],
     boundary_nodes=NamedTuple[],boundary_lanes=NamedTuple[],
     open_paths=Vector{Vector{Int}}(),
     open_path_rows=Int[])
end

function disk_table(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    edges=real_edge_table(v)
    edge_count=length(edges)
    vertex_count == 0 && isempty(edges) && return _empty_disk_table()
    dictionary=_disk_node_dictionary()

    endpoint_records=NamedTuple[]
    endpoint_by_edge=Dict{Tuple{Int,Symbol},Int}()
    next_boundary_node=6vertex_count+1
    for edge in edges
        for (side,endpoint) in ((:start,edge.crossing[1]),(:end,edge.crossing[2]))
            endpoint == 0 || continue
            public_side=side == :start ? :in : :out
            ids=(next_boundary_node,next_boundary_node+1,next_boundary_node+2)
            next_boundary_node += 3
            record=(strand=edge.strand,side=public_side,lane=0,node_ids=ids,edge_id=edge.id)
            record_index=length(endpoint_records)+1
            push!(endpoint_records,record)
            endpoint_by_edge[(edge.id,side)]=record_index
        end
    end
    node_count=next_boundary_node-1
    parent=collect(1:node_count)
    find_root(node)=begin
        while parent[node] != node
            parent[node]=parent[parent[node]]
            node=parent[node]
        end
        node
    end
    unite(a,b)=begin
        ra,rb=find_root(a),find_root(b)
        ra != rb && (parent[rb]=ra)
    end
    graph_edges=Tuple{Int,Int,Int}[]
    for edge in edges
        source,target=edge.crossing
        source_nodes = if source == 0
            endpoint_records[endpoint_by_edge[(edge.id,:start)]].node_ids
        else
            source_slots=(1:3) .+ 6 .+ 3*(source > 0)
            Tuple(6*(abs(source)-1) .+ dictionary[orientation[abs(source)]][source_slots])
        end
        target_nodes = if target == 0
            endpoint_records[endpoint_by_edge[(edge.id,:end)]].node_ids
        else
            target_slots=(1:3) .+ 3*(target > 0)
            Tuple(6*(abs(target)-1) .+ dictionary[orientation[abs(target)]][target_slots])
        end
        for index in 1:3
            coefficient=index == 3 ? -edge.id : edge.id
            a,b=index == 3 ?
                (target_nodes[index],source_nodes[index]) :
                (source_nodes[index],target_nodes[index])
            push!(graph_edges,(a,b,coefficient)); unite(a,b)
        end
    end
    roots=unique(find_root.(1:node_count))
    root_index=Dict(root=>index for (index,root) in enumerate(roots))
    delta=zeros(Int,length(roots),edge_count)
    for (a,_,coefficient) in graph_edges
        row=root_index[find_root(a)]
        delta[row,abs(coefficient)] += sign(coefficient)
    end
    component_nodes=[[node for node in 1:node_count if find_root(node) == root]
                     for root in roots]
    component_edges=[[(a,b,coefficient) for (a,b,coefficient) in graph_edges
                      if find_root(a) == root] for root in roots]
    is_cycle=[length(edges_) == length(nodes_)
              for (nodes_,edges_) in zip(component_nodes,component_edges)]
    paths=Vector{Vector{Int}}(undef,length(roots))
    edge_signs=Vector{Vector{Int}}(undef,length(roots))
    for row in eachindex(roots)
        paths[row],edge_signs[row]=_disk_component_path(
            component_nodes[row],component_edges[row],is_cycle[row])
    end

    positive_dots=[0,1,0,0,1,0]
    negative_dots=[0,0,1,0,0,1]
    node_dots=zeros(Int,node_count)
    for vertex in 1:vertex_count, local_node in 1:6
        dots=orientation[vertex] > 0 ? positive_dots : negative_dots
        node_dots[6*(vertex-1)+local_node]=dots[local_node]
    end
    ndots=[sum(node_dots[nodes_]) for nodes_ in component_nodes]
    cp2=2 .* Int.(is_cycle) .- ndots
    cp=cp2 .// 2
    path_lengths=length.(paths)
    ids=collect(eachindex(roots))
    endpoint_at(node) = begin
        for record in endpoint_records
            lane=findfirst(==(node),record.node_ids)
            lane === nothing || return (strand=record.strand,side=record.side,lane=lane)
        end
        nothing
    end
    boundary_start=Union{Nothing,NamedTuple}[
        is_cycle[row] ? nothing : endpoint_at(first(paths[row])) for row in ids]
    boundary_end=Union{Nothing,NamedTuple}[
        is_cycle[row] ? nothing : endpoint_at(last(paths[row])) for row in ids]
    boundary_nodes=[(strand=record.strand,side=record.side,
                     node_ids=record.node_ids,edge_id=record.edge_id)
                    for record in endpoint_records]
    boundary_lanes=[(strand=record.strand,side=record.side,lane=lane,
                     node_id=record.node_ids[lane],node_ids=record.node_ids,
                     edge_id=record.edge_id)
                    for record in endpoint_records for lane in 1:3]
    open_path_rows=findall(!,is_cycle)
    open_paths=paths[open_path_rows]
    (ID=ids,paths=paths,edges_s=edge_signs,Ndots=ndots,isCycle=is_cycle,
     delta=delta,cp=cp,cp2=cp2,path_lengths=path_lengths,
     boundary_start=boundary_start,boundary_end=boundary_end,
     boundary_nodes=boundary_nodes,boundary_lanes=boundary_lanes,
     open_paths=open_paths,
     open_path_rows=open_path_rows)
end

function _nullspace_integer(A::Matrix{Int})
    rows,columns=size(A)
    reduced=Rational{BigInt}.(A)
    pivots=Int[]
    pivot_row=1
    for column in 1:columns
        source=findfirst(row -> reduced[row,column] != 0,pivot_row:rows)
        source === nothing && continue
        source=source+pivot_row-1
        reduced[pivot_row,:],reduced[source,:]=copy(reduced[source,:]),copy(reduced[pivot_row,:])
        reduced[pivot_row,:] ./= reduced[pivot_row,column]
        for row in 1:rows
            row == pivot_row && continue
            reduced[row,:] .-= reduced[row,column].*reduced[pivot_row,:]
        end
        push!(pivots,column); pivot_row += 1
        pivot_row > rows && break
    end
    free=setdiff(collect(1:columns),pivots)
    isempty(free) && return zeros(Int,0,columns)
    basis=Vector{Vector{Int}}()
    for free_column in free
        vector=zeros(Rational{BigInt},columns); vector[free_column]=1
        for (row,pivot) in enumerate(pivots)
            vector[pivot]=-reduced[row,free_column]
        end
        common_denominator=foldl(lcm,(value.den for value in vector);init=BigInt(1))
        integers=BigInt[value.num for value in vector.*common_denominator]
        divisor=foldl(gcd,abs.(integers);init=BigInt(0))
        divisor != 0 && (integers .÷= divisor)
        push!(basis,Int.(integers))
    end
    reduce(vcat,(permutedims(vector) for vector in basis))
end

function _rational_particular(A::Matrix{Int},b::Vector{Int})
    rows,columns=size(A)
    augmented=hcat(Rational{BigInt}.(A),Rational{BigInt}.(b))
    pivot_row=1
    pivots=Int[]
    for column in 1:columns
        source=findfirst(row -> augmented[row,column] != 0,pivot_row:rows)
        source === nothing && continue
        source += pivot_row-1
        augmented[pivot_row,:],augmented[source,:]=copy(augmented[source,:]),copy(augmented[pivot_row,:])
        augmented[pivot_row,:] ./= augmented[pivot_row,column]
        for row in 1:rows
            row == pivot_row && continue
            augmented[row,:] .-= augmented[row,column].*augmented[pivot_row,:]
        end
        push!(pivots,column); pivot_row += 1
        pivot_row > rows && break
    end
    any(row -> all(==(0),augmented[row,1:columns]) && augmented[row,end] != 0,1:rows) &&
        throw(ArgumentError("weight equations are inconsistent"))
    result=zeros(Rational{BigInt},columns)
    for (row,pivot) in enumerate(pivots)
        result[pivot]=augmented[row,end]
    end
    result
end

function _integer_particular(A::Matrix{Int},b::Vector{Int})
    rows,columns=size(A)
    model=Model(GLPK.Optimizer); set_silent(model)
    @variable(model,weight[1:columns],Int)
    for row in 1:rows
        @constraint(model,sum(A[row,column]*weight[column] for column in 1:columns)==b[row])
    end
    @objective(model,Min,0)
    optimize!(model)
    status=termination_status(model)
    status == MOI.INFEASIBLE && return nothing
    status == MOI.OPTIMAL ||
        throw(ArgumentError("integer weight feasibility did not terminate: $status"))
    candidate=round.(Int,value.(weight))
    all(A*candidate .== b) ||
        throw(ArgumentError("integer weight solver returned a non-solution"))
    candidate
end

"""MATLAB `minimizeWeightSolution`: nearest lattice representative in a kernel coset.

The particular solution may be rational.  MATLAB's `Matrix(ZZ).solve_right`
returns that rational solution when the Euler cochain is not in the integral
image; only the homogeneous kernel and its H-move quotient are integral.
"""
function _minimize_weight_solution(weight::AbstractVector{<:Number},kernel::Matrix{Int})
    weight=collect(weight)
    size(kernel,1)==0 && return weight
    # MATLAB's `-W / H`: minimize ||xH+W|| over real row vectors x.
    coefficients=vec(round.(Int,(-permutedims(Float64.(weight))) / Float64.(kernel)))
    current=weight+vec(transpose(coefficients)*kernel)
    value=sum(abs2,current)
    improved=true
    while improved
        improved=false
        for index in eachindex(coefficients)
            for delta in (1,-1)
                candidate_coefficients=copy(coefficients); candidate_coefficients[index]+=delta
                candidate=weight+vec(transpose(candidate_coefficients)*kernel)
                candidate_value=sum(abs2,candidate)
                if candidate_value<value
                    coefficients=candidate_coefficients; current=candidate; value=candidate_value
                    improved=true
                    break
                end
            end
            improved && break
        end
    end
    current
end

"""Legacy MATLAB-compatible `(base, kernel, hperp)` diagnostic tuple.

`base` may be rational when `assign=false`. `kernel` contains integer-valued
rows spanning the rational nullspace, NOT necessarily the full integer lattice;
`hperp` is a heuristic rational complement, not an integer quotient certificate.
Use `solve_integer_weights` for a complete certified integer solution lattice.
`assign=true` refuses rational-only systems and inexact/out-of-range storage.
The historical small representative is retained for existing invariant callers.
"""
function calculate_weight(v::VirtualLink; assign=false)
    v.weight_mode==:integer || throw(ArgumentError("calculate_weight requires integer mode; use solve_spin_weights for Z2"))
    any(strand_cuts(v)) &&
        throw(ArgumentError("calculate_weight currently requires closed strands; use disk_table boundary metadata for open data"))
    if assign
        exact=solve_integer_weights(v)
        exact.status==:integer_solution ||
            throw(ArgumentError("no integral weight assignment ($(exact.status)); use solve_integer_weights for the obstruction"))
    end
    disks=disk_table(v)
    edge_count=size(disks.delta,2)
    base=_rational_particular(2disks.delta,-disks.cp2)
    if any(value -> value.den != 1,base)
        integer_base=_integer_particular(2disks.delta,-disks.cp2)
        integer_base !== nothing && (base=integer_base)
    else
        base=Int[value.num for value in base]
    end
    kernel=_nullspace_integer(disks.delta)
    base=_minimize_weight_solution(base,kernel)
    hmoves=permutedims(hmove_matrix(v))
    selected_rows=Vector{Vector{Int}}()
    span=Float64.(hmoves)
    old_rank=rank(span)
    for row in eachrow(kernel)
        candidate=vcat(span,permutedims(Float64.(row)))
        new_rank=rank(candidate)
        if new_rank > old_rank
            push!(selected_rows,collect(row)); span=candidate; old_rank=new_rank
        end
    end
    selected=isempty(selected_rows) ? zeros(Int,0,edge_count) :
        reduce(vcat,(permutedims(row) for row in selected_rows))
    for selected_index in axes(selected,1)
        current=collect(selected[selected_index,:])
        improved=true
        while improved
            improved=false
            for move in eachrow(hmoves),direction in (-1,1)
                candidate=current+direction.*move
                if sum(abs2,candidate) < sum(abs2,current) ||
                   (sum(abs2,candidate) == sum(abs2,current) && Tuple(candidate) > Tuple(current))
                    current=collect(candidate); improved=true
                end
            end
        end
        selected[selected_index,:]=current
    end
    if assign
        all(x -> isinteger(x) && abs(x)<=big(2)^53-1,base) ||
            throw(ArgumentError("legacy representative is not a safely storable integer; use solve_integer_weights"))
        all(iszero,2 .* BigInt.(disks.delta)*BigInt.(base)+BigInt.(disks.cp2)) ||
            throw(ArgumentError("legacy representative fails the exact integer weight equation"))
        weights = Vector{Vector{Float64}}(); offset = 1
        code,_=real_gauss_code(v)
        for g in code
            n=length(g)
            push!(weights,Float64.(base[offset:offset+n-1])); offset += n
        end
        set_weight!(v, weights)
    end
    base,kernel,selected
end

function set_weight!(v::VirtualLink, weight; mode=v.weight_mode)
    _check_weight_mode(mode)
    real_counts=[_strand_arc_count(v,strand) for strand in eachindex(v.gauss)]
    convert_value = mode==:z2 ? (x -> begin
        x isa Real && isfinite(x) && isinteger(x) || throw(ArgumentError("Z2 weights require finite integer-valued inputs"))
        Float64(mod(BigInt(x),2))
    end) : Float64
    weights=Vector{Vector{Float64}}()
    if weight isa AbstractMatrix || weight isa AbstractVector{<:Number}
        flat=vec(weight)
        length(flat)==sum(real_counts) || throw(DimensionMismatch("one weight per real arc is required"))
        idx=1
        for n in real_counts
            push!(weights,convert_value.(flat[idx:idx+n-1])); idx+=n
        end
    else
        length(weight)==length(real_counts) || throw(DimensionMismatch("one weight vector per Gauss component is required"))
        for (values,count_) in zip(weight,real_counts)
            length(values)==count_ || throw(DimensionMismatch("one weight per real arc is required"))
            push!(weights,convert_value.(values))
        end
    end
    v.weights=weights; v.weight_mode=mode; v.is_weighted=true; v
end

"""Assign selected MATLAB `REdgeTable.Crossing` weights without relying on row order."""
function set_weight_by_crossing!(v::VirtualLink,crossings,weights)
    pairs=[Tuple(Int.(pair)) for pair in eachrow(Int.(crossings))]
    length(pairs) == length(weights) || throw(DimensionMismatch("one weight per crossing pair is required"))
    convert_value = v.weight_mode==:z2 ? (x -> begin
        x isa Real && isfinite(x) && isinteger(x) || throw(ArgumentError("Z2 weights require finite integer-valued inputs"))
        Float64(mod(BigInt(x),2))
    end) : Float64
    assignments=Dict(pair=>convert_value(weight) for (pair,weight) in zip(pairs,weights))
    table=real_edge_table(v)
    by_strand=[Float64[] for _ in v.gauss]
    for edge in table
        push!(by_strand[edge.strand],get(assignments,edge.crossing,edge.weight))
    end
    set_weight!(v,by_strand)
end

function _graph_a_connected(v::VirtualLink)
    all(!=(0), v.orientation) || return false
    any(word -> any(==(0), word), v.gauss) && return false
    code, orientation = real_gauss_code(v)
    isempty(orientation) && return length(code) <= 1
    any(isempty, code) && return false
    length(connected_components(ograph_graph_a(v))) == 1
end

function is_closed(v::VirtualLink)
    code,orientation=real_gauss_code(v)
    c1 = length(code) == 1
    unsupported = any(==(0), v.orientation) ||
                  any(word -> any(==(0), word), v.gauss) ||
                  (!isempty(orientation) && any(isempty, code))
    unsupported && return false, (c1, false, false)
    c2 = _graph_a_connected(v)
    c3 = size(disk_table(v).delta,1) == length(orientation)+1
    all((c1,c2,c3)), (c1,c2,c3)
end

function mirror_manifold!(v::VirtualLink)
    _validate_cut_markers(v.gauss)
    v.gauss = [-g for g in v.gauss]; v.orientation .*= -1; v
end
function cyclic_shift!(v::VirtualLink, strand::Integer, amount::Integer)
    _validate_cut_markers(v.gauss)
    _is_cut_component(v.gauss[strand]) && amount != 0 &&
        throw(ArgumentError("cyclic_shift! on a cut strand is undefined; only amount=0 is allowed"))
    v.gauss[strand] = circshift(v.gauss[strand], amount); v
end

function disjoint_sum!(a::VirtualLink, b::VirtualLink)
    _validate_cut_markers(a.gauss); _validate_cut_markers(b.gauss)
    offset=length(a.orientation)
    append!(a.gauss, [[sign(x)*offset+x for x in g] for g in b.gauss])
    append!(a.orientation,b.orientation)
    a.real_only = a.real_only && b.real_only
    a.is_weighted = false
    a.weights = [zeros(Float64,_strand_arc_count(a,strand))
                 for strand in eachindex(a.gauss)]
    a
end

function _connected_sum_selection(selection)
    if selection isa AbstractMatrix
        size(selection,2) == 2 && size(selection,1) in (1,2) ||
            throw(ArgumentError("connected-sum selection must have shape 1×2 or 2×2"))
        all(value -> value isa Integer,selection) ||
            throw(ArgumentError("connected-sum selection entries must be integers"))
        strands=(Int(selection[1,1]),Int(selection[1,2]))
        anchors=size(selection,1) == 2 ?
            (Int(selection[2,1]),Int(selection[2,2])) : nothing
        return strands,anchors
    end

    applicable(iterate,selection) ||
        throw(ArgumentError("connected-sum selection must contain two strand indices"))
    values=collect(selection)
    length(values) == 2 ||
        throw(ArgumentError("connected-sum selection must contain two strand indices"))
    all(value -> value isa Integer,values) ||
        throw(ArgumentError("connected-sum strand indices must be integers"))
    (Int(values[1]),Int(values[2])),nothing
end

function _rotate_connected_sum_component(component::Vector{Int},anchor::Int,strand::Int)
    positions=findall(==(anchor),component)
    length(positions) == 1 ||
        throw(ArgumentError("connected-sum anchor $anchor must occur exactly once on strand $strand"))
    index=only(positions)
    vcat(component[index:end],component[begin:index-1])
end

"""Join two closed real-Gauss components by the normal-o-graph connected sum.

`selection` may be a two-element strand-index iterable, a `1×2` matrix of
strand indices, or a `2×2` matrix whose second row gives the signed crossing
at which each component traversal starts.  Anchor selection is performed on the
real-Gauss code.  The four-vertex splice and labels are specified in
`docs/specs/ograph-connected-sum.md`.
"""
function connected_sum!(v::VirtualLink, selection=(1,2))
    any(strand_cuts(v)) &&
        throw(ArgumentError("connected_sum! currently requires closed strands; open boundaries would be discarded"))
    code,orientation=real_gauss_code(v)
    (i,j),anchors=_connected_sum_selection(selection)
    i==j && throw(ArgumentError("distinct strands required"))
    1 <= i <= length(code) && 1 <= j <= length(code) || throw(BoundsError(code,(i,j)))
    first_component=code[i]
    second_component=code[j]
    if anchors !== nothing
        first_component=_rotate_connected_sum_component(first_component,anchors[1],i)
        second_component=_rotate_connected_sum_component(second_component,anchors[2],j)
    end
    vertex_count=length(orientation)
    a=vertex_count+1
    b=vertex_count+2
    c=vertex_count+3
    d=vertex_count+4
    bridge1=[-b,-a,a,d]
    bridge2=[-d,-c,c,b]
    code[i]=vcat(first_component,bridge1,second_component,bridge2)
    deleteat!(code,j)
    v.gauss=code
    v.orientation=vcat(orientation,[-1,1,-1,1])
    v.real_only=true
    v.is_weighted=false
    v.weights=[zeros(Float64,length(g)) for g in code]
    v
end

"""Convert a knot diagram to the version-2 knot-complement o-graph."""
function knot_complement!(v::VirtualLink)
    any(strand_cuts(v)) &&
        throw(ArgumentError("knot_complement! currently requires closed strands"))
    code,orientation=real_gauss_code(v)
    vertex_count=length(orientation)
    expand(vertex,inverted)=begin
        base=4(abs(vertex)-1)+1
        offset=(vertex < 0 ? 1+2*(orientation[abs(vertex)] < 0) : 0)+2inverted
        sign(vertex).*(base .+ mod.([0,1].+offset,4))
    end
    forward=[reduce(vcat,(expand(vertex,false) for vertex in component);init=Int[]) for component in code]
    inverse=[reduce(vcat,(expand(vertex,true) for vertex in reverse(component));init=Int[]) for component in code]
    v.gauss=vcat(forward,inverse)
    v.orientation=repeat([-1,1,-1,1],vertex_count)
    v.real_only=true
    v.is_weighted=false
    v.weights=[zeros(Float64,length(g)) for g in v.gauss]
    v.strand_types=fill('o',length(v.gauss))
    v
end

"""Convert a VirtualLink to the URDiagram ordering used by MATLAB `VL2URD`."""
function vl_to_urdiagram(v::VirtualLink)
    v.weight_mode==:integer || throw(ArgumentError("UR invariants require integer mode; Z2 weights are not integer lifts"))
    any(strand_cuts(v)) &&
        throw(ArgumentError("vl_to_urdiagram currently requires closed strands"))
    code,orientation=real_gauss_code(v)
    vertices=Int[]
    for (index,component) in enumerate(code)
        index > 1 && push!(vertices,0)
        append!(vertices,reverse(component))
    end
    diagram=URDiagram()
    set_data_cvs!(diagram,1,vertices,orientation)
end

function ur_invariant(v::VirtualLink)
    trace_value(vl_to_urdiagram(v))
end

"""Compatibility entry for the exact U(R) MST backend (corrected word order)."""
ur_prime_invariant(v::VirtualLink) = Z_MST(v, URAlgebra())
