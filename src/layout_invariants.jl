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

"""
    LayoutDiagnostic(rule, severity, message; edges=Int[], at=nothing) -> LayoutDiagnostic

Store one diagnostic about a routed link diagram's geometry.

# Arguments
`rule` and `severity` convert to `Symbol`, `message` to `String`, `edges` to
`Vector{Int}`, and `at` is absent or convertible to `ComplexF64`. The
five-positional-field constructor is also available. Conversion errors
propagate; rule names, severity names, edge membership and location are not
validated. Geometric rules are owned by
[plot specification §11.2](../specs/virtual-link-plot.md#112-mandatory-geometric-properties).

# Returns
`LayoutDiagnostic` with the named rule, severity, explanatory message,
related edge identifiers and optional layout-space location.

# Guarantees
1. Fields retain the supplied diagnostic data after conversion.
   [proved/agent: by construction]
2. Inputs are not modified; an already typed edge vector is shared.
   [proved/agent: by construction]

# Not guaranteed
Constructing a diagnostic does not run a geometric check or establish that
the reported violation occurred. Arbitrary rule and severity symbols are
accepted.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §11.2 — geometric properties reported by the layout checker.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; gauss=[[1,-2],[-1,2]], orientation=[1,1]); layout = orthogonal_layout(v); e = first(layout.edges);

julia> note = LayoutDiagnostic(:P3, :info, "endpoint"; edges=[e.id], at=first(e.points));

julia> note.at == layout.crossings[e.tail]
true
```

Role: construct
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

"""
    check_layout(layout::LinkLayout; rules=Topology.VirtualLinks.LAYOUT_RULES) -> Vector{LayoutDiagnostic}

Return diagnostics for selected geometric properties of a routed link diagram.

# Arguments
The caller supplies L4 [LinkLayout data](../specs/virtual-link-plot.md#6-layout-stage-contracts);
malformed records may
raise indexing or arithmetic errors and have no certified geometric meaning.
`rules` supports membership tests for the property symbols of
[§11.2](../specs/virtual-link-plot.md#112-mandatory-geometric-properties).
Unknown symbols are ignored. No complete input-shape validation is performed.

# Returns
A vector of `LayoutDiagnostic` records, each with the rule, severity, message,
related edges and optional location. Empty-edge layouts return no diagnostics.

# Guarantees
1. The implemented selected predicates of §11.2 are evaluated and their
   reported failures returned. [proved/agent: by construction]
2. The layout is not modified. [proved/agent: by construction]

# Not guaranteed
An empty result does not prove bend optimality, completeness of the geometric
predicates, planar isotopy, or correspondence to a specified source diagram.
In particular, empty-edge layouts bypass coordinate checks.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §11.2 — geometric property diagnostics.

# Examples
```jldoctest
julia> l = LinkLayout(zeros(Int,0,4), ComplexF64[], EdgeLayout[], Int[], Int[]);

julia> isempty(check_layout(l))
true
```

Role: condition
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

"""
    check_layout(v::VirtualLink; bending_numbers=nothing, rules=Topology.VirtualLinks.LAYOUT_RULES) -> Vector{LayoutDiagnostic}

Return geometric diagnostics for an orthogonal layout of a virtual-link diagram.

# Arguments
Checked: the input and optional bending data must satisfy
[`orthogonal_layout`](@ref); its input, optimization and routing errors propagate.
The caller supplies coherent diagram data. `rules` selects the predicates of
[§11.2](../specs/virtual-link-plot.md#112-mandatory-geometric-properties).

# Returns
The diagnostics obtained by `check_layout` on the computed layout.

# Guarantees
1. The result concerns the computed diagram geometry, with the limitations
   of the layout overload of `check_layout`. [proved/agent: by construction]
2. The input diagram is not modified. [proved/agent: by construction]

# Not guaranteed
The layout is not returned. A later layout computation can choose another
optimum. Successful diagnostics alone do not certify bend optimality.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §§6 and 11.2 — layout and its geometric diagnostics.

# Examples
```jldoctest
julia> isempty(check_layout(VirtualLink()))
true
```

Role: condition
"""
check_layout(v::VirtualLink;bending_numbers=nothing,rules=LAYOUT_RULES) =
    check_layout(orthogonal_layout(v;bending_numbers=bending_numbers);rules=rules)

"""
    is_valid_layout(x; kwargs...) -> Bool

Test whether the selected layout predicates report no error diagnostic.

# Arguments
`x` is a `LinkLayout` or `VirtualLink`; keyword arguments and checked/caller
preconditions are those of [`check_layout`](@ref). Its exceptions propagate.

# Returns
A `Bool`, excluding only diagnostics whose severity is `:error`.

# Guarantees
1. The result is true exactly when `check_layout(x; kwargs...)` contains no
   error diagnostic in the [property model](../specs/virtual-link-plot.md#112-mandatory-geometric-properties). [proved/agent: by construction]
2. The supplied layout or diagram is not modified. [proved/agent: by construction]

# Not guaranteed
This does not certify unselected properties, bend optimality or topology.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §11 — levels of geometric acceptance.

# Examples
```jldoctest
julia> is_valid_layout(LinkLayout(zeros(Int,0,4), ComplexF64[], EdgeLayout[], Int[], Int[]))
true
```

Role: condition
"""
is_valid_layout(x;kwargs...) =
    !any(d->d.severity === :error,check_layout(x;kwargs...))

"""
    failing_rules(x; kwargs...) -> Vector{Symbol}

List the selected geometric property names that report layout errors.

# Arguments
`x` and keyword arguments obey the checked and caller's preconditions of
[`check_layout`](@ref). Layout, optimization and indexing errors propagate.
For a `VirtualLink`, repeated layout construction may select different optima.

# Returns
Property symbols in the declared rule order of
[§11.2](../specs/virtual-link-plot.md#112-mandatory-geometric-properties).

# Guarantees
1. Each returned property reported an error in its check. Warning-only
   diagnostics do not cause inclusion. [proved/agent: by construction]
2. Supplied layout or diagram data are not modified. [proved/agent: by construction]

# Not guaranteed
For diagram input, the returned names need not describe one common layout:
the current entry point may construct a layout separately for each property.
For a fixed layout record, the result does concern that single record.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §11.2 — property names and diagnostic meaning.

# Examples
```jldoctest
julia> isempty(failing_rules(LinkLayout(zeros(Int,0,4), ComplexF64[], EdgeLayout[], Int[], Int[])))
true
```

Role: condition
"""
failing_rules(x;kwargs...) =
    [r for r in LAYOUT_RULES if any(d->d.rule===r && d.severity===:error,check_layout(x;kwargs...))]
