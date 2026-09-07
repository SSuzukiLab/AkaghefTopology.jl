# Integral PS rules: MST25 submit.tex, Section;IPS, Figs IPSmove2/Ci_ver3;
# Taguchi23 framing.tex definition IPS-move, original IPS.jpg and Ci.png.
# Internal API: candidate remains the frozen unweighted catalogue witness.

"""Freeze a choice of the paper's edge-portion weight `k` for a forward PS.

When both selected portions belong to one real edge of weight `W`, the other
portion has weight `l=W-k`. Applicability and exact integer bounds are checked
again against the actual weighted source on application.
"""
function weighted_ps_candidate(candidate::MoveCandidate;k)
    candidate.type=="PS" && candidate.direction==:forward ||
        throw(ArgumentError("a PS split belongs only to a forward PS candidate"))
    MoveCandidate(candidate.type,candidate.parameter,copy(candidate.vertices),
        copy(candidate.edges),candidate.direction,candidate.before,candidate.after,
        copy(candidate.locations),BigInt(_move_integer(k)))
end

function _ps_geometry_candidate(candidate)
    MoveCandidate(candidate.type,candidate.parameter,copy(candidate.vertices),
        copy(candidate.edges),candidate.direction,candidate.before,candidate.after,
        copy(candidate.locations))
end

function _ps_unweighted(v::VirtualLink)
    result=VirtualLink()
    set_data!(result;rgauss=deepcopy(first(real_gauss_code(v))),
              orientation=copy(last(real_gauss_code(v))))
end

function _ps_pattern_weights(pattern,weights)
    weights=[[_move_integer(value) for value in rule] for rule in weights]
    MovePattern(pattern.parameter,pattern.gauss,pattern.orientation,pattern.capacity;
                weights=weights)
end

function _ps_check_framed(v)
    weights=BigInt.(_integer_edge_weights(v))
    disks=disk_table(v)
    residual=2BigInt.(disks.delta)*weights+BigInt.(disks.cp2)
    all(iszero,residual) || throw(ArgumentError("PS weights do not satisfy the framed integral equation"))
    true
end

function _ps_candidate_weights(candidate,before_weights,after_weights)
    MoveCandidate(candidate.type,candidate.parameter,copy(candidate.vertices),
        copy(candidate.edges),candidate.direction,
        _ps_pattern_weights(candidate.before,before_weights),
        _ps_pattern_weights(candidate.after,after_weights),copy(candidate.locations))
end

_ps_node_dot(orientation,node)=begin
    vertex=div(node-1,6)+1; local_node=mod1(node,6)
    orientation[vertex]>0 ? Int(local_node in (2,5)) : Int(local_node in (3,6))
end

function _ps_integral_a(parameter,p2::Integer)
    numerator=parameter in ("1","2") ? -p2-1 : parameter=="3" ? p2 :
              parameter=="4" ? -p2-2 : throw(ArgumentError("unknown PS type"))
    iseven(numerator) || throw(ArgumentError("PS rotation has wrong parity for integral a"))
    div(numerator,2)
end

function _ps_forward_rotation(v,candidate;split_k=nothing)
    condition=ps_global_condition(_ps_unweighted(v),_ps_geometry_candidate(candidate))
    condition.valid || throw(ArgumentError("PS Method-B condition failed: $(condition.reason)"))
    witness=condition.witness
    trace=witness.circuit_trace
    left=first(trace)
    right_index=findfirst(s -> s.edge_id==witness.right_edge_id &&
                              s.lane==witness.right_lane,trace)
    right_index===nothing && throw(ArgumentError("PS right boundary is absent from circuit"))
    right=trace[right_index]
    edge_weights=BigInt.(_integer_edge_weights(v))
    same_edge=left.edge_id==right.edge_id
    split_k!==nothing && !same_edge && throw(ArgumentError("PS split parameter requires one shared real edge"))
    k=split_k===nothing ? edge_weights[left.edge_id] : BigInt(_move_integer(split_k))
    l=same_edge ? edge_weights[left.edge_id]-k : edge_weights[right.edge_id]
    # Coincident empty slots follow converter order fragment 2 then fragment 1.
    # By default a shared edge's weight is placed wholly in paper left cut k;
    # an explicit candidate may instead freeze another k/l partition.
    # At a second occurrence of this edge every full lane still contributes k+l.
    fragment_left=candidate.parameter=="2" ? 2 : 1
    fragment_right=3-fragment_left
    no_vertex=same_edge && left.lane==right.lane
    middle=no_vertex ? eltype(trace)[] : trace[2:right_index-1]
    nodes=no_vertex ? Int[] : vcat(left.to_node,[s.to_node for s in middle])
    contributions=[(edge_id=s.edge_id,lane=s.lane,
        component=s.component,tail_occurrence=s.tail_occurrence,
        head_occurrence=s.head_occurrence,from_node=s.from_node,to_node=s.to_node,
        coefficient=s.lane==3 ? -1 : 1,weight=edge_weights[s.edge_id]) for s in middle]
    partial=big(0)
    partial_contributions=NamedTuple{(:boundary,:cut_fragment,:coefficient,:weight),Tuple{String,Int,Int,BigInt}}[]
    if same_edge && !no_vertex
        # On each endpoint edge count other selected cuts met between the
        # bottom cut and its attached vertex. Order 2<1 is directed edge order.
        cut_order=Dict(2=>1,1=>2)
        cut_weights=Dict(fragment_left=>k,fragment_right=>l)
        for (boundary,role,is_start) in ((left,fragment_left,true),(right,fragment_right,false))
            direction=boundary.lane==3 ? -1 : 1
            for other in (1,2)
                other==role && continue
                after=cut_order[other]>cut_order[role]
                selected=is_start ? (direction==1 ? after : !after) :
                                    (direction==1 ? !after : after)
                if selected
                    partial+=direction*cut_weights[other]
                    push!(partial_contributions,(boundary=is_start ? "left" : "right",
                        cut_fragment=other,coefficient=direction,weight=cut_weights[other]))
                end
            end
        end
    end
    dots=[_ps_node_dot(v.orientation,node) for node in nodes]
    integral_sum=sum(c.coefficient*c.weight for c in contributions;init=big(0))+partial
    p2=2integral_sum-sum(dots)
    a=_ps_integral_a(candidate.parameter,p2)
    full_rotation2=2sum((s.lane==3 ? -1 : 1)*edge_weights[s.edge_id] for s in trace)-
        sum(_ps_node_dot(v.orientation,s.from_node) for s in trace)
    full_rotation2==-2 || throw(ArgumentError("PS selected disk has nonframed full rotation"))
    (version=1,source="MST25 Section;IPS / Taguchi23 IPS-move",parameter=candidate.parameter,
     left_edge_id=left.edge_id,right_edge_id=right.edge_id,
     left_lane=left.lane,right_lane=right.lane,k=k,l=l,p2=p2,a=a,
     no_vertex=no_vertex,cut_order=[2,1],nodes=nodes,dots=dots,
     forward_locations=[collect(location) for location in candidate.locations],
     forward_split_k=same_edge ? k : nothing,
     contributions=contributions,partial_contributions=partial_contributions,partial_edge_sum=partial,
     integer_arc_sum=integral_sum,dot_count=sum(dots),full_rotation2=full_rotation2)
end

function _ps_forward_plan(v,candidate;split_k=candidate.ps_split_k)
    detail=_ps_forward_rotation(v,candidate;split_k=split_k)
    a,k,l=detail.a,detail.k,detail.l
    # Read in oriented catalogue fragment order, NOT left-to-right page order.
    rules=candidate.parameter=="1" ? [[-a,0,a+k],[0,0,l]] :
          candidate.parameter=="2" ? [[0,0,l],[-a,0,a+k]] :
          candidate.parameter=="3" ? [[a+k,1,-a-1],[0,0,l]] :
                                    [[-a,0,a+k],[l+1,-1,0]]
    # The generic converter retains each incoming old edge once. Subtract it
    # at its first cut, leaving the exact six edge portions of the source figure.
    for i in eachindex(candidate.locations)
        earlier=any(j -> j>i && candidate.locations[j]==candidate.locations[i],1:2)
        earlier && continue
        component,head=candidate.locations[i]
        rules[i][1]-=BigInt(v.weights[component][head])
    end
    prepared=_ps_candidate_weights(candidate,[[0],[0]],rules)
    target=deepcopy(v)
    _apply_pattern_move!(target,prepared;prepared_weights=true)
    _ps_check_framed(target)
    (candidate=prepared,detail=detail,target=target)
end

function _ps_inverse_plan(v,candidate)
    internal=candidate.parameter=="3" ? [1,0] : candidate.parameter=="4" ? [0,-1] : [0,0]
    # Contract by exact sums, including nonzero internal 1/-1. This recovers
    # the LHS cochain before p is calculated; no weight equation is re-solved.
    contraction=_ps_candidate_weights(candidate,[[0,x,0] for x in internal],[[x] for x in internal])
    lhs=deepcopy(v)
    _apply_pattern_move!(lhs,contraction;prepared_weights=true)
    _ps_check_framed(lhs)
    condition=ps_global_condition(_ps_unweighted(v),candidate)
    condition.valid || throw(ArgumentError("inverse PS Method-B condition failed"))
    boundaries=condition.witness.gap_boundaries
    data=make_move_data()
    forward=MoveCandidate("PS",candidate.parameter,Int[],NTuple{2,Int}[],:forward,
        only(data["PS-L"]),data["PS-R"][parse(Int,candidate.parameter)],
        [(b.target_component,b.target_head_occurrence) for b in boundaries])
    plan=_ps_forward_plan(lhs,forward)
    if plan.detail.left_edge_id==plan.detail.right_edge_id && candidate.parameter in ("2","4")
        # These two types can retain the chosen k/l split as different edge
        # cochains after insertion. Infer it from the source figure instead of
        # silently rejecting every split except the forward default k=W,l=0.
        if candidate.parameter=="2"
            component,start=candidate.locations[2] # paper left lower edge=-a
            recovered_k=BigInt(v.weights[component][start])+plan.detail.a+plan.detail.k
        else
            component,start=candidate.locations[1] # paper left upper edge=a+k
            recovered_k=BigInt(v.weights[component][mod1(start+2,length(v.gauss[component]))])-plan.detail.a
        end
        plan=_ps_forward_plan(lhs,forward;split_k=recovered_k)
    end
    # Undo the inverse deletion's label map, including nontrailing selected IDs.
    _,orientation=real_gauss_code(v)
    mapping,_,_=_vertex_maps(orientation,candidate.vertices,candidate.before,candidate.after)
    restoration=Dict(new=>old for (old,new) in mapping)
    n=length(lhs.orientation)
    for role in 1:2
        restoration[n+role]=candidate.vertices[role]
    end
    restored=[_relabel_word(word,restoration) for word in plan.target.gauss]
    restored_weights=_align_move_weights(v.gauss,restored,plan.target.weights)
    restored_weights==v.weights ||
        throw(ArgumentError("inverse PS source boundary weights do not reproduce the integral source rule"))
    (candidate=contraction,detail=merge(plan.detail,(direction="inverse",)),target=lhs)
end

"Recompute the source-specific integral PS plan, with no persistent candidate cache."
function _weighted_ps_plan(v::VirtualLink,candidate::MoveCandidate)
    v.is_weighted || throw(ArgumentError("integral PS requires assigned weights"))
    _check_move_weights(v)
    candidate.ps_split_k!==nothing && candidate.direction!=:forward &&
        throw(ArgumentError("inverse PS infers its split from the actual source"))
    validation=validate_ps_move(_ps_unweighted(v),_ps_geometry_candidate(candidate))
    validation.valid || throw(ArgumentError("PS applicability failed: $(validation.reason)"))
    _ps_check_framed(v)
    candidate.direction==:forward ? _ps_forward_plan(v,candidate) : _ps_inverse_plan(v,candidate)
end

_weighted_ps_target(v,candidate)=_weighted_ps_plan(v,candidate).target
