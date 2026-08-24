"""Property checks for a computed orthogonal layout.

Plot acceptance is property-based, not a comparison against MATLAB: the
minimum-bend MILP is degenerate and any optimum is a valid answer, so a stored
reference embedding proves nothing. See `PLOT_PIPELINE.md` section 1 and
`LIVESCRIPT_PORTS.md` for the criterion these checks replace.

The rules below are what "the algorithm delivered what it specifies" means
operationally. `P4v` is the cheapest and the strongest: it reads coordinates
only, never `EdgeLayout.tail`/`head`, so a router that mis-routes an edge cannot
hide behind its own bookkeeping.
"""

struct LayoutDiagnostic
    rule::Symbol
    severity::Symbol
    message::String
    edges::Vector{Int}
    at::Union{Nothing,ComplexF64}
end

LayoutDiagnostic(rule,severity,message;edges=Int[],at=nothing) =
    LayoutDiagnostic(rule,severity,message,edges,at)

function Base.show(io::IO,d::LayoutDiagnostic)
    print(io,"[",d.rule,"/",d.severity,"] ",d.message)
    isempty(d.edges) || print(io," (edge ",join(d.edges,", "),")")
    d.at === nothing || print(io," at ",_fmt(d.at))
end

_fmt(z::ComplexF64) = "($(real(z)), $(imag(z)))"

const _TOL = 1e-9
_same(a::ComplexF64,b::ComplexF64) = abs(a-b) < _TOL
_key(z::ComplexF64) = (round(real(z);digits=9),round(imag(z);digits=9))

_segments(layout::LinkLayout) =
    [(edge.id,index,edge.points[index],edge.points[index+1])
     for edge in layout.edges for index in 1:length(edge.points)-1]

"""Meeting of two closed segments: `(:none,nothing)`, `(:point,z)` or `(:overlap,nothing)`."""
function _meet(a1::ComplexF64,a2::ComplexF64,b1::ComplexF64,b2::ComplexF64)
    min(real(a1),real(a2)) > max(real(b1),real(b2))+_TOL && return (:none,nothing)
    max(real(a1),real(a2)) < min(real(b1),real(b2))-_TOL && return (:none,nothing)
    min(imag(a1),imag(a2)) > max(imag(b1),imag(b2))+_TOL && return (:none,nothing)
    max(imag(a1),imag(a2)) < min(imag(b1),imag(b2))-_TOL && return (:none,nothing)

    da=a2-a1; db=b2-b1; offset=b1-a1
    cross=real(da)*imag(db)-imag(da)*real(db)
    if abs(cross) < 1e-12                       # parallel
        abs(real(da)*imag(offset)-imag(da)*real(offset)) > 1e-12 && return (:none,nothing)
        low  = max(min(real(a1),real(a2)),min(real(b1),real(b2)))
        high = min(max(real(a1),real(a2)),max(real(b1),real(b2)))
        low2 = max(min(imag(a1),imag(a2)),min(imag(b1),imag(b2)))
        high2= min(max(imag(a1),imag(a2)),max(imag(b1),imag(b2)))
        (high-low > _TOL || high2-low2 > _TOL) && return (:overlap,nothing)
        return (:point,complex(low,low2))
    end
    t=(real(offset)*imag(db)-imag(offset)*real(db))/cross
    u=(real(offset)*imag(da)-imag(offset)*real(da))/cross
    (-_TOL <= t <= 1+_TOL && -_TOL <= u <= 1+_TOL) || return (:none,nothing)
    (:point,a1+t*da)
end

"""P1: every segment runs along an axis."""
function _check_p1(layout::LinkLayout)
    out=LayoutDiagnostic[]
    for (id,_,a,b) in _segments(layout)
        (real(a)==real(b) || imag(a)==imag(b)) && continue
        push!(out,LayoutDiagnostic(:P1,:error,
            "segment is neither horizontal nor vertical";edges=[id],at=a))
    end
    out
end

"""P3: each edge starts on its tail crossing and ends on its head crossing."""
function _check_p3(layout::LinkLayout)
    out=LayoutDiagnostic[]
    for edge in layout.edges
        checkbounds(Bool,layout.crossings,edge.tail) || continue
        checkbounds(Bool,layout.crossings,edge.head) || continue
        _same(first(edge.points),layout.crossings[edge.tail]) ||
            push!(out,LayoutDiagnostic(:P3,:error,
                "edge does not start on crossing $(edge.tail); wanted $(_fmt(layout.crossings[edge.tail]))";
                edges=[edge.id],at=first(edge.points)))
        _same(last(edge.points),layout.crossings[edge.head]) ||
            push!(out,LayoutDiagnostic(:P3,:error,
                "edge does not end on crossing $(edge.head); wanted $(_fmt(layout.crossings[edge.head]))";
                edges=[edge.id],at=last(edge.points)))
    end
    out
end

"""P4v: every crossing is the endpoint of exactly four edge-ends.

Coordinates only -- `tail`/`head` are deliberately unused, so this survives a
router that records the wrong incidence. A self-loop contributes two ends at
the same crossing, which keeps the total at four. Two crossings sharing a
position show up as degree eight; a lost edge as a deficit.
"""
function _check_p4v(layout::LinkLayout)
    out=LayoutDiagnostic[]
    isempty(layout.edges) && return out
    degree=Dict{Tuple{Float64,Float64},Int}()
    owners=Dict{Tuple{Float64,Float64},Vector{Int}}()
    for edge in layout.edges, z in (first(edge.points),last(edge.points))
        k=_key(z)
        degree[k]=get(degree,k,0)+1
        push!(get!(owners,k,Int[]),edge.id)
    end
    crossings=Set(_key(c) for c in layout.crossings)
    for (k,d) in sort!(collect(degree))
        at=complex(k[1],k[2])
        if !(k in crossings)
            push!(out,LayoutDiagnostic(:P4v,:error,
                "$d edge-end(s) terminate where there is no crossing";
                edges=sort(owners[k]),at=at))
        elseif d != 4
            push!(out,LayoutDiagnostic(:P4v,:error,
                "crossing has degree $d, expected 4";edges=sort(owners[k]),at=at))
        end
    end
    for c in layout.crossings
        haskey(degree,_key(c)) && continue
        push!(out,LayoutDiagnostic(:P4v,:error,
            "crossing has no incident edge-end";at=c))
    end
    out
end

"""P5: crossing positions are pairwise distinct."""
function _check_p5(layout::LinkLayout)
    out=LayoutDiagnostic[]
    seen=Dict{Tuple{Float64,Float64},Int}()
    for (index,c) in enumerate(layout.crossings)
        k=_key(c)
        if haskey(seen,k)
            push!(out,LayoutDiagnostic(:P5,:error,
                "crossings $(seen[k]) and $index share a position";at=c))
        else
            seen[k]=index
        end
    end
    out
end

"""P6: two segments may meet only at a crossing position."""
function _check_p6(layout::LinkLayout)
    out=LayoutDiagnostic[]
    segments=_segments(layout)
    crossings=[c for c in layout.crossings]
    at_crossing(z)=any(c->_same(c,z),crossings)
    for i in 1:length(segments), j in i+1:length(segments)
        (id1,pos1,a1,a2)=segments[i]
        (id2,pos2,b1,b2)=segments[j]
        id1==id2 && abs(pos1-pos2)<=1 && continue      # same edge, consecutive
        kind,point=_meet(a1,a2,b1,b2)
        kind === :none && continue
        if kind === :overlap
            push!(out,LayoutDiagnostic(:P6,:error,
                "segments overlap along a shared span";edges=sort([id1,id2])))
        elseif !at_crossing(point)
            push!(out,LayoutDiagnostic(:P6,:error,
                "segments meet away from any crossing";edges=sort([id1,id2]),at=point))
        end
    end
    out
end

"""P7: an edge does not run over a crossing it is not incident to."""
function _check_p7(layout::LinkLayout)
    out=LayoutDiagnostic[]
    crossings=[c for c in layout.crossings]
    for edge in layout.edges
        length(edge.points) < 3 && continue
        for z in edge.points[2:end-1]
            any(c->_same(c,z),crossings) || continue
            push!(out,LayoutDiagnostic(:P7,:error,
                "edge passes through a crossing";edges=[edge.id],at=z))
        end
    end
    out
end

const LAYOUT_RULES = (:P1,:P3,:P4v,:P5,:P6,:P7)

"""Run every layout property and return the diagnostics, most structural first.

`rules` selects a subset. An empty result means the layout is a valid
bend-minimal orthogonal embedding as far as these rules can tell; it says
nothing about whether the bend count is optimal, which is checked separately
against the MILP objective.
"""
function check_layout(layout::LinkLayout;rules=LAYOUT_RULES)
    isempty(layout.edges) && return LayoutDiagnostic[]
    out=LayoutDiagnostic[]
    :P1  in rules && append!(out,_check_p1(layout))
    :P5  in rules && append!(out,_check_p5(layout))
    :P4v in rules && append!(out,_check_p4v(layout))
    :P3  in rules && append!(out,_check_p3(layout))
    :P7  in rules && append!(out,_check_p7(layout))
    :P6  in rules && append!(out,_check_p6(layout))
    out
end

check_layout(v::VirtualLink;bending_numbers=nothing,rules=LAYOUT_RULES) =
    check_layout(orthogonal_layout(v;bending_numbers=bending_numbers);rules=rules)

"""True when no rule reports an error."""
is_valid_layout(x;kwargs...) =
    !any(d->d.severity === :error,check_layout(x;kwargs...))

"""Rules that reported at least one error, in `LAYOUT_RULES` order."""
failing_rules(x;kwargs...) =
    [r for r in LAYOUT_RULES if any(d->d.rule===r && d.severity===:error,check_layout(x;kwargs...))]
