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
            outgoing=mod1(position+1,length(word))
            result.weights[component][outgoing]=1-result.weights[component][outgoing]
        end
    end
    result,toggled
end

function _spin_cp_matches(v,candidate,pattern)
    component,start=only(candidate.locations)
    word=only(pattern.gauss); values=only(pattern.weights)
    n=length(v.weights[component])
    all(v.weights[component][mod1(start+offset,n)]==values[offset+1]
        for offset in 1:length(word)-1)
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
Both source and target must pass closed-normal C1-C3 and the spin equation.
Relative/cut spin boundary conventions are deliberately unsupported.
"""
function spin_cp_transport(v::VirtualLink,candidate::MoveCandidate;normalize::Bool=true)
    candidate.type=="CP" || throw(ArgumentError("spin CP requires a CP candidate"))
    _spin_edge_weights(v)
    any(_is_cut_component,v.gauss) && throw(ArgumentError("spin CP currently requires a closed source"))
    geometry=_weighted_unweight(v)
    certify_closed_ograph(geometry)
    validate_spin_weight(v)
    options=_movable_geometry(geometry,"CP";param=candidate.parameter,
        vertices=candidate.vertices,edges=candidate.edges,strict=true)
    any(option -> _same_candidate(option,candidate),options) ||
        throw(ArgumentError("spin CP candidate does not match the current source"))
    length(candidate.locations)==1 || throw(ArgumentError("spin CP requires one source interval"))
    vertices=abs.(candidate.vertices)
    for mask in 0:(normalize ? (1 << length(vertices))-1 : 0)
        prepared,toggled=_spin_cp_hgauge(v,vertices,mask)
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
            target=VirtualLink(;weight_mode=:z2)
            set_data!(target;rgauss=words,orientation=signs)
            set_weight!(target,reduce(vcat,weights;init=Float64[]);mode=:z2)
            certify_closed_ograph(_weighted_unweight(target))
            validate_spin_weight(target)
            target_word=only([_map_fragment(g,after_vertices) for g in after.gauss])
            target_locations=[(i,s) for (i,g) in enumerate(real_gauss_code(target)[1])
                for s in _cyclic_occurrences(g,target_word)]
            return (target=target,case=t,direction=candidate.direction,
                target_vertices=after_vertices,target_locations=target_locations,
                source_h_vertices=toggled,normalized_source=prepared,
                source_figure="w-f2194518-7624-4b07-9eb5-3e04e0a7d036#^m23-fig-0094",
                source_already_normalized=isempty(toggled),modulus=2)
        end
    end
    throw(ArgumentError("spin CP source cannot match either published case" *
        (normalize ? " by local H moves" : "; enable normalize for explicit local H moves")))
end

"Return the target of `spin_cp_transport`; see it for normalization and inverse semantics."
spin_cp_move(v::VirtualLink,candidate::MoveCandidate;normalize::Bool=true)=
    spin_cp_transport(v,candidate;normalize=normalize).target
