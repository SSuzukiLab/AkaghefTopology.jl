# Source: Mih23, managed Work w-f2194518-7624-4b07-9eb5-3e04e0a7d036,
# spin CP figure ^m23-fig-0094 and H figure ^m23-fig-0093.
# Weights are attached to incoming edge occurrences, not unordered vertex pairs.
# RHS roles in the published drawing: 1 bottom-right, 2 bottom-left,
# 3 top-left, 4 top-right (loop), 5 left-middle.

function _spin_cp_patterns(t::Integer)
    t in (0,1) || throw(ArgumentError("spin CP case must be 0 or 1"))
    left=MovePattern("dir",[[1,-2,2,-3,3,-1]],[-1,-1,-1];
        weights=[[0,0,1,t,1,0,0]])
    right=MovePattern("inv",[[1,2,-5,3,-2,-3,-4,4,-1,5]],[1,1,1,-1,-1];
        weights=[[0,0,1-t,0,1,0,0,1,1,0,0]])
    left,right
end

function _spin_cp_hgauge(v::VirtualLink,vertices,mask::Integer)
    result=deepcopy(v)
    code,_=real_gauss_code(result)
    toggled=Int[]
    for (bit,vertex) in enumerate(vertices)
        iszero(mask & (1 << (bit-1))) && continue
        push!(toggled,vertex)
        for (component,word) in enumerate(code), position in eachindex(word)
            abs(word[position])==vertex || continue
            result.weights[component][position]=1-result.weights[component][position]
            outgoing=_is_cut_component(word) ? position+1 : mod1(position+1,length(word))
            result.weights[component][outgoing]=1-result.weights[component][outgoing]
        end
    end
    result,toggled
end

function _spin_cp_matches(v,candidate,pattern)
    component,start=only(candidate.locations)
    word=only(pattern.gauss); values=only(pattern.weights)
    code,_=real_gauss_code(v)
    cut=_is_cut_component(code[component])
    n=cut ? length(code[component])-1 : length(code[component])
    cut && start+length(word)-2>n && return false
    all(v.weights[component][cut ? start+offset : mod1(start+offset,n)]==values[offset+1]
        for offset in 1:length(word)-1)
end

"""Capture and validate relative Z4 targets for a cut CP source."""
function _spin_cp_boundary_data(v::VirtualLink;boundary=nothing)
    any(_is_cut_component,v.gauss) ||
        return (source=nothing,target=nothing,measured=nothing)
    data=spin_boundary_data(v)
    open_rows=findall(!,Bool.(data.isCycle))
    isempty(open_rows) && throw(ArgumentError("a cut CP source has no open disk paths"))
    source_input=if boundary===nothing
        Dict((data.boundary_start[row],data.boundary_end[row])=>data.measured4[row]
             for row in open_rows)
    elseif boundary isa NamedTuple && haskey(boundary,:source)
        boundary.source
    elseif boundary isa AbstractDict || boundary isa AbstractVector
        boundary
    else
        throw(ArgumentError("spin CP boundary must be a directed endpoint map or (source=..., target=...)"))
    end
    # Normalize row-aligned and row-key spellings to the directed endpoint
    # dictionary used by endpoint transport below.
    targets=_spin_boundary_targets(disk_table(v),source_input)
    source_boundary=Dict((data.boundary_start[row],data.boundary_end[row])=>targets[row]
                         for row in open_rows)
    validate_spin_weight(v;boundary=source_boundary)
    (source=source_boundary,
     target=(boundary isa NamedTuple && haskey(boundary,:target) ? boundary.target : nothing),
     measured=data.measured4)
end

function _spin_cp_transport_boundary(source,target,endpoint_trace,source_boundary;
                                     explicit_target=nothing)
    source_boundary===nothing && return (source=nothing,target=nothing)
    sdata=spin_boundary_data(source); tdata=spin_boundary_data(target)
    source_rows=findall(!,Bool.(sdata.isCycle)); target_rows=findall(!,Bool.(tdata.isCycle))
    length(endpoint_trace)==2count(_is_cut_component,source.gauss) ||
        throw(ArgumentError("spin CP boundary transport lost a source endpoint"))
    endpoint_map=Dict((row.source.strand,row.source.side)=>
                      (row.target.strand,row.target.side) for row in endpoint_trace)
    map_endpoint(endpoint)=begin
        mapped=get(endpoint_map,(endpoint.strand,endpoint.side),nothing)
        mapped===nothing && throw(ArgumentError("spin CP boundary transport has no endpoint lineage"))
        (strand=mapped[1],side=mapped[2],lane=endpoint.lane)
    end
    mapped_source=Dict{Tuple,BigInt}()
    for row in source_rows
        key=(sdata.boundary_start[row],sdata.boundary_end[row])
        value=get(source_boundary,key,nothing)
        value===nothing && throw(ArgumentError("spin CP source boundary map is incomplete"))
        mapped_source[(map_endpoint(key[1]),map_endpoint(key[2]))]=mod(BigInt(value),4)
    end
    # The transported source beta is the preservation contract.  An explicit
    # target map, when supplied, is an additional assertion and never replaces
    # this check.
    validate_spin_weight(target;boundary=mapped_source)
    explicit_target===nothing || validate_spin_weight(target;boundary=explicit_target)
    actual=Dict((tdata.boundary_start[row],tdata.boundary_end[row])=>tdata.measured4[row]
                for row in target_rows)
    Set(keys(actual))==Set(keys(mapped_source)) ||
        throw(ArgumentError("spin CP target boundary identities do not match transported endpoints"))
    (source=source_boundary,target=mapped_source,
     explicit_target=explicit_target)
end

"""Transport a closed-normal Z2-weighted graph by the two published spin CP rules.

The source is Mih23, managed Work w-f2194518-7624-4b07-9eb5-3e04e0a7d036,
anchors ^m23-fig-0094 (CP) and ^m23-fig-0093 (H).
`normalize=false` accepts exactly the displayed internal edge weights; the
external edge weights remain arbitrary and unchanged. `normalize=true` first
tries explicit local H moves (in deterministic bit-mask order). This is H*CP,
not an arbitrary solution of the target spin equations. The returned record
names every H vertex; it is not a durable independent replay certificate.

The exact displayed rules are mutually inverse after aligning cyclic edge
occurrences (the existing move engine may change the cycle's starting point).
`target_vertices` gives the output CP roles, including newly allocated labels.
With H normalization a reverse
CP returns an H-equivalent representative, not necessarily the original weights.
For an uncut source, both source and target must pass closed-normal C1-C3
and the spin equation.
For a cut source, `boundary=nothing` captures the measured open-path residues
as directed endpoint targets; an explicit map is validated before transport.
"""
function spin_cp_transport(v::VirtualLink,candidate::MoveCandidate;
                           normalize::Bool=true,boundary=nothing)
    candidate.type=="CP" || throw(ArgumentError("spin CP requires a CP candidate"))
    _spin_edge_weights(v)
    cut_source=any(_is_cut_component,v.gauss)
    !cut_source && boundary!==nothing &&
        throw(ArgumentError("spin CP boundary targets require a cut source"))
    boundary_record=_spin_cp_boundary_data(v;boundary=boundary)
    geometry=_weighted_unweight(v)
    if cut_source
        _validate_gauss(geometry.gauss,geometry.orientation)
    else
        certify_closed_ograph(geometry)
        validate_spin_weight(v)
    end
    options=_movable_geometry(geometry,"CP";param=candidate.parameter,
        vertices=candidate.vertices,edges=candidate.edges,strict=true)
    any(option -> _same_candidate(option,candidate),options) ||
        throw(ArgumentError("spin CP candidate does not match the current source"))
    length(candidate.locations)==1 || throw(ArgumentError("spin CP requires one source interval"))
    vertices=abs.(candidate.vertices)
    for mask in 0:(normalize ? (1 << length(vertices))-1 : 0)
        prepared,toggled=_spin_cp_hgauge(v,vertices,mask)
        if cut_source
            try
                validate_spin_weight(prepared;boundary=boundary_record.source)
            catch error
                error isa Union{ArgumentError,DimensionMismatch,BoundsError} || rethrow()
                continue
            end
        end
        for t in 0:1
            left,right=_spin_cp_patterns(t)
            before,after=candidate.direction==:forward ? (left,right) : (right,left)
            _spin_cp_matches(prepared,candidate,before) || continue
            frozen=MoveCandidate("CP",candidate.parameter,candidate.vertices,candidate.edges,
                candidate.direction,before,after,candidate.locations)
            code,orientation=real_gauss_code(prepared)
            mapping,after_vertices,new_orientation=_vertex_maps(orientation,
                candidate.vertices,before,after;role_ordered=true)
            words,weights,signs=_weighted_pattern_rewrite(code,prepared.weights,frozen,
                [_map_fragment(g,candidate.vertices) for g in before.gauss],
                [_map_fragment(g,after_vertices) for g in after.gauss],mapping,new_orientation)
            # Execute geometry once through the generic endpoint-aware
            # transport.  This keeps cut marker placement and endpoint
            # lineage independent of weighted arithmetic.
            geometry_transport=cut_source ? move_boundary_transport(geometry,candidate) : nothing
            target=if cut_source
                deepcopy(geometry_transport.target)
            else
                VirtualLink(;weight_mode=:z2)
            end
            if cut_source
                aligned=_align_move_weights(target.gauss,words,weights)
                set_weight!(target,aligned;mode=:z2)
                boundary_record=_spin_cp_transport_boundary(v,target,
                    geometry_transport.endpoints,boundary_record.source;
                    explicit_target=boundary_record.target)
            else
                set_data!(target;rgauss=words,orientation=signs)
                set_weight!(target,reduce(vcat,weights;init=Float64[]);mode=:z2)
                certify_closed_ograph(_weighted_unweight(target))
                validate_spin_weight(target)
            end
            target_word=only([_map_fragment(g,after_vertices) for g in after.gauss])
            target_locations=[(i,s) for (i,g) in enumerate(real_gauss_code(target)[1])
                for s in _pattern_occurrences(g,target_word)]
            return (target=target,case=t,direction=candidate.direction,
                target_vertices=after_vertices,target_locations=target_locations,
                source_h_vertices=toggled,normalized_source=prepared,
                source_figure="w-f2194518-7624-4b07-9eb5-3e04e0a7d036#^m23-fig-0094",
                source_already_normalized=isempty(toggled),modulus=2,
                boundary_modulus=cut_source ? 4 : nothing,
                boundary=boundary_record,
                source_boundary=boundary_record.source,
                target_boundary=boundary_record.target,
                boundary_transport=cut_source ? geometry_transport.endpoints : nothing)
        end
    end
    throw(ArgumentError("spin CP source cannot match either published case" *
        (normalize ? " by local H moves" : "; enable normalize for explicit local H moves")))
end

"Return the target of `spin_cp_transport`; see it for normalization and inverse semantics."
spin_cp_move(v::VirtualLink,candidate::MoveCandidate;
             normalize::Bool=true,boundary=nothing)=
    spin_cp_transport(v,candidate;normalize=normalize,boundary=boundary).target
