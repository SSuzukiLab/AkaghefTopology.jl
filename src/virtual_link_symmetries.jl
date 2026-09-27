# Source: MST25 Appendix A, definitions and Proposition A.1,
# [[MST25#^mst25-blk-0186]], [[MST25#^mst25-blk-0188]].
# Encoding and scope: docs/specs/ograph-symmetries.md.
# No integer/Z2 framing transport is asserted by these operations.
function _check_diagram_symmetry(v::VirtualLink; reflection=false)
    _validate_cut_markers(v.gauss)
    v.is_weighted && throw(ArgumentError("symmetry transport of framing weights is not specified"))
    reflection && any(iszero,v.orientation) &&
        throw(ArgumentError("plane reflection of stored virtual-crossing cyclic orders is not supported"))
end

"""Exchange over/under at every real crossing; do not reverse strand directions.
MST25 Γ×: negate real passage signs and vertex signs. Virtual visits are unchanged.
This is a diagram operation, not a certificate that the input is a normal o-graph.
"""
function switch_crossings!(v::VirtualLink)
    _check_diagram_symmetry(v)
    v.gauss=[[x==0 || v.orientation[abs(x)]==0 ? x : -x for x in w] for w in v.gauss]
    v.orientation .*= -1
    v
end
switch_crossings(v::VirtualLink)=switch_crossings!(deepcopy(v))

"""Reflect the plane of a real-crossing diagram (MST25 Γ*).
Passage signs and traversal stay fixed; vertex signs change. This is NOT over/under exchange.
Explicit virtual crossings and weighted inputs are unsupported.
"""
function reflect_diagram!(v::VirtualLink)
    _check_diagram_symmetry(v;reflection=true)
    v.orientation .*= -1
    v
end
reflect_diagram(v::VirtualLink)=reflect_diagram!(deepcopy(v))

"""Reverse every strand's traversal (MST25 →Γ), keeping crossing data fixed.
Each cut marker stays last: the input/output endpoints exchange roles.
This primitive does not introduce a direction-reversed 0–2 move.
"""
function reverse_strands!(v::VirtualLink)
    _check_diagram_symmetry(v)
    v.gauss=[_is_cut_component(w) ? vcat(reverse(w[1:end-1]),0) : reverse(w) for w in v.gauss]
    v
end
reverse_strands(v::VirtualLink)=reverse_strands!(deepcopy(v))

function _check_mst_symmetry(v::VirtualLink)
    _check_diagram_symmetry(v;reflection=true)
    first(is_closed(v)) || throw(ArgumentError("MST semantic methods require a closed normal o-graph"))
end

"""MST25 Proposition A.1: Γ× represents (M^op, -v^op).
Requires a closed normal unweighted o-graph; does not define a framed lift.
Source: [[MST25#^mst25-blk-0188]] (not merely a manifold-only mirror).
"""
function mst_sign_change!(v::VirtualLink)
    _check_mst_symmetry(v)
    switch_crossings!(v)
end
mst_sign_change(v::VirtualLink)=mst_sign_change!(deepcopy(v))

"""MST25 Proposition A.1: →(Γ*) represents (M^op, v^op).
Requires a closed normal unweighted o-graph; does not define a framed lift.
Source: [[MST25#^mst25-blk-0188]].
"""
function mst_reflect_reverse!(v::VirtualLink)
    _check_mst_symmetry(v)
    reflect_diagram!(v)
    reverse_strands!(v)
end
mst_reflect_reverse(v::VirtualLink)=mst_reflect_reverse!(deepcopy(v))

"""Legacy diagram-level alias of `switch_crossings!`, NOT plane reflection.
For a checked combed interpretation use `mst_sign_change!` (MST25 Proposition A.1).
Weighted inputs now fail closed rather than silently retaining their weights.
"""
mirror_manifold!(v::VirtualLink)=switch_crossings!(v)
