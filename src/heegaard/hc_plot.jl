# cf: lib/topology.jl/docs/specs/heegaard-plot.md — Julia drawing entry and acceptance boundary
# Final crossing-angle criterion (heegaard-plot.md); it does not constrain relaxation steps.
const _HP_MIN_CROSSING_ANGLE = 20.0
include("hc_plot_geometry.jl")
include("hc_plot_empty.jl")
include("hc_plot_placement.jl")
include("hc_plot_audit.jl")
include("hc_plot_certificates.jl")
include("hc_plot_svg.jl")
include("hc_plot_forces.jl")
include("hc_plot_gate.jl")
include("hc_plot_relax.jl")
include("hc_plot_reference_chart.jl")
include("hc_plot_reference.jl")
include("hc_plot_marks.jl")

"""
    plot(hc::HeegaardCode; crossing_labels=false, show=:both, DEBUG_MODE=false, rounds=24, local_steps=20, log=devnull, reference_alpha=:auto, mark_style=:filled, infinity_face=nothing) -> NamedTuple

Produce an audited planar SVG rendering of a resolved Heegaard diagram.

# Arguments
`hc` must pass `validate_diagram`. Options use the [plot domain and acceptance conditions](../specs/heegaard-plot.md). Invalid `show`, nonpositive smoothing counts or invalid `mark_style` raises `ArgumentError`; chart, mark, placement and audit errors propagate. `crossing_labels` requests alpha–beta crossing labels.

# Returns
A named tuple containing `svg`, `placement`, `audit`, `success` and path-dependent layout data. Labelled output also contains `crossing_labels`; marked and reference-chart outputs retain their corresponding mark and contact records.

# Guarantees
1. The SVG and placement depict the supplied diagram in the selected planar chart; `success` records the drawing path's acceptance checks. It must be inspected before treating the rendering as accepted. [proved/agent: by construction]
2. Input HC data are not modified. [proved/agent: by construction]

# Not guaranteed
Success is not a manifold classification or a guarantee that `drawing(hc)` succeeds; the [single-crossing defect](../specs/heegaard-plot.md#single-crossing-drawing-failure) remains. This function returns data and does not save a file or open a viewer.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-plot.md — planar rendering and acceptance.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); h = H.resolve_hc(lc, Dict("l1"=>:beta, "m1"=>:alpha));

julia> try H.plot(h; rounds=0) catch e; e isa ArgumentError end
true
```

Role: draw
"""
function plot(hc::HeegaardCode;crossing_labels::Bool=false,kwargs...)
    result=_hp_plot(hc;kwargs...)
    crossing_labels || return result
    rows=_hp_alpha_beta_labels(hc,result)
    merge(result,(svg=_hp_crossing_labels_svg(result.svg,result.placement,rows),crossing_labels=rows))
end
# cf: lib/topology.jl/docs/specs/heegaard-solid-view.md §7.4 item 4 — the same labels on the original figure and on the drawing
function _hp_alpha_beta_labels(hc,result)
    p=result.placement;rows=Tuple{Int,HCPlotPoint}[]
    for (i,c) in enumerate(result.audit.crossings)
        Set((p.runs[c.run_a].role,p.runs[c.run_b].role))==Set((:alpha,:beta)) || continue
        label=get(result.audit.labels,i,nothing);isnothing(label) || push!(rows,(label,c.point))
    end
    family=get(result,:reference_alpha,nothing)
    if !isnothing(family)
        contacts=result.reference_contacts
        labels=Set(abs(x) for i in 1:hc.lc.genus for x in hc.lc.words[string(family)*string(i)])
        for r in p.runs
            r.role==:beta && !r.closed || continue
            for side in (:start,:end)
                label=get(contacts,(r.name,side),0);label in labels || continue
                push!(rows,(label,side==:start ? first(r.points) : last(r.points)))
            end
        end
    end
    sort!(rows)
end
function _hp_crossing_labels_svg(svg,p,rows)
    io=IOBuffer()
    for (label,q) in rows
        at="x=\"$(q[1]+7)\" y=\"$(q[2]-7)\" font-size=\"12\" font-family=\"sans-serif\""
        println(io,"<text $at fill=\"none\" stroke=\"#fcfbf7\" stroke-width=\"3.5\" stroke-linejoin=\"round\">V$label</text>")
        println(io,"<text data-crossing=\"$label\" $at fill=\"#222222\">V$label</text>")
    end
    replace(svg,"</svg>"=>String(take!(io))*"</svg>")
end
function _hp_plot(hc::HeegaardCode;show::Symbol=:both,DEBUG_MODE::Bool=false,
              rounds::Int=24,local_steps::Int=20,log::IO=devnull,reference_alpha::Symbol=:auto,
              mark_style::Symbol=:filled,infinity_face=nothing)
    show in (:both,:alpha,:beta) || throw(ArgumentError("show must be :both, :alpha or :beta"))
    rounds>0 && local_steps>0 || throw(ArgumentError("zero smoothing steps cannot be reported as success"))
    mark_style in (:filled,:outline) || throw(ArgumentError("mark_style must be :filled or :outline"))
    validate_diagram(hc)
    family=_hp_reference_family(hc,reference_alpha)
    if !isnothing(family)
        result=_hp_plot_reference(hc,family;show,DEBUG_MODE,rounds,local_steps,log,infinity_face)
        result=merge(result,(contacts=result.reference_contacts,))
        isempty(hc.points) && isempty(hc.disks) && return result
        marks=_hp_marks(hc,result.placement,result.reference_contacts;geometry=result.audit)
        return merge(result,(marks=marks,svg=_hp_marks_svg(result.svg,marks,result.placement;style=mark_style)))
    end
    chart = _hp_disk_chart(hc)
    if isnothing(chart) && (!isempty(hc.points) || !isempty(hc.disks))
        ctx=_hp_place(hc.lc,hc.roles;reference=true,infinity_face)
        result=_hp_relax(ctx.placement;rounds,local_steps,DEBUG_MODE,log,fixed=ctx.fixed,check_final_angle=false,
            placement_check=_hp_empty_guard_check(hc,ctx.placement))
        marks=_hp_marks(hc,result.placement,ctx.contacts;geometry=result.audit)
        visible=_hp_rebuild(result.placement,[r for r in result.placement.runs if r.role!=:none])
        angles=[result.audit.stencil_angles[i] for (i,c) in enumerate(result.audit.crossings)
            if Set((result.placement.runs[c.run_a].role,result.placement.runs[c.run_b].role))==Set((:alpha,:beta))]
        angle_ok=isempty(angles) || minimum(angles) >= _HP_MIN_CROSSING_ANGLE
        audit=merge(result.audit,(angle_ok=angle_ok,ok=result.audit.ok && angle_ok))
        svg=_hp_marks_svg(_hp_svg(visible;show,title="Heegaard diagram"),marks,result.placement;style=mark_style)
        return merge(result,(audit=audit,marks=marks,svg=svg,success=result.steps>0 && audit.ok,
            mark_contacts=ctx.contacts,contacts=ctx.contacts,visual_acceptance=:pending))
    end
    placed,contacts=_hp_place(hc.lc,hc.roles;with_contacts=true,infinity_face)
    result=_hp_relax(placed;rounds,local_steps,DEBUG_MODE,log,placement_check=_hp_empty_guard_check(hc,placed))
    ordinary = merge(result,(svg=_hp_svg(result.placement;show,title=DEBUG_MODE ? "Heegaard diagram · DEBUG_MODE" : "Heegaard diagram"),
        success=result.steps>0 && result.audit.ok,contacts=contacts,visual_acceptance=:pending))
    isnothing(chart) && return ordinary
    marks = _hp_marks(chart.chart,result.placement,contacts;geometry=result.audit)
    disks = [merge(d,(support=hc.disks[d.name],)) for d in marks.disks]
    polygons = Dict(only(ids) => marks.face_polygons[id] for (id,ids) in chart.supports if length(ids)==1)
    marks = merge(marks,(disks=disks,face_polygons=polygons))
    merge(ordinary,(marks=marks,svg=_hp_marks_svg(ordinary.svg,marks,result.placement;style=mark_style)))
end

include("hc_drawing.jl")
