include("solid3d_smoothing.jl")

# cf: docs/specs/heegaard-plot.md and solid3d.md — reuse the full reference chart.
"""Place a role-free LC in the planar cut chart, retaining all reference runs."""
function plot(lc::LickorishCode; infinity_face=nothing)
    validate_lc_input(lc)
    _validate_diagram_passages(lc.words)
    _lc_validate_surface!(lc)
    ctx = _hp_place(lc; reference=true, infinity_face)
    placement = ctx.placement
    audit = _hp_audit(placement)
    if !audit.structure
        # cf: hc_plot_relax.jl::_hp_relax — initial fine sampling separates
        # interior crossings from coincident segment endpoints.
        runs = [r.name in ctx.fixed ? r :
            _hp_run(r,_hp_sample(r.points,fill(10.0,length(r.points)),r.closed))
            for r in placement.runs]
        placement = _hp_rebuild(placement,runs)
        audit = _hp_audit(placement)
    end
    merge(ctx, (; placement, audit, success=audit.ok, svg=_hp_svg(placement; title="Lickorish Code")))
end

# Use the plotted full K chart. General HC smoothing can omit m,n (and move
# curves across them). In that case select the public full LC chart as a unit.
function _solid_chart(hc::HeegaardCode)
    validate_diagram(hc)
    family = _hp_reference_family(hc,:auto)
    marked = !isempty(hc.points) || !isempty(hc.disks)
    if !isnothing(family) || (marked && isnothing(_hp_disk_chart(hc)))
        return plot(hc)
    end
    plot(hc.lc)
end
_solid_chart(lc::LickorishCode) = plot(lc)

# cf: hc_plot_marks.jl::_hp_mark_regions — oriented edges of the actual plotted faces.
function _solid_data(input, crossing_labels)
    hc = input isa HeegaardCode
    lc = hc ? input.lc : input
    lc.genus >= 1 || throw(ArgumentError("solid chart requires positive genus"))
    result = _solid_chart(input)
    result.success || throw(ArgumentError("planar drawing failed its geometry checks"))
    result = _solid_smooth_chart(result,lc)
    p = result.placement
    regions = _hp_mark_regions((;lc),p,result.contacts; geometry=result.audit)
    paths = NamedTuple[]
    for face in sort!(collect(keys(regions.faces))), arc in regions.faces[face]
        copy_on_cut = arc.orientation == -1 && !isnothing(arc.cut)
        arc.orientation == 1 || copy_on_cut || continue
        name, label = arc.edge
        points = copy_on_cut ? reverse(arc.points) : arc.points
        rim = !isnothing(arc.cut) && first(arc.cut) == :hole
        push!(paths,(;name, label, points, closed=false, copy_on_cut, rim))
    end
    for run in p.runs
        isempty(lc.words[run.component]) || continue
        push!(paths,(name=run.component,label=0,points=run.points,closed=run.closed,
                     copy_on_cut=false,rim=false))
    end
    holes = [(x=c[1], y=c[2], r=r) for (_, (c,r)) in p.holes]
    roles = hc ? input.roles : Dict{String,Symbol}()
    kinds = Dict(n=>get(roles,n,n in _lc_k_names(lc.genus) ? :reference : :curve) for n in keys(lc.words))
    (;width=p.width,height=p.height,square=p.square,holes,paths,kinds,crossing_labels)
end

# Embed the chart as inert SVG. XML escaping is shared with the planar writer.
function _solid_svg(data)
    io = IOBuffer()
    println(io, """<svg id="chart" width="$(data.width)" height="$(data.height)" data-square="$(data.square)" data-labels="$(data.crossing_labels)">""")
    for h in data.holes
        println(io, """<circle cx="$(h.x)" cy="$(h.y)" r="$(h.r)"/>""")
    end
    for p in data.paths
        points = join(("$(q[1]),$(q[2])" for q in p.points), " ")
        name = _hp_xml(p.name)
        println(io, """<polyline data-name="$name" data-kind="$(data.kinds[p.name])" data-label="$(p.label)" data-copy="$(p.copy_on_cut)" data-rim="$(p.rim)" data-closed="$(p.closed)" points="$points"/>""")
    end
    print(io, "</svg>")
    String(take!(io))
end

"""Write an interactive, rotatable HTML view of an LC or HC and return its absolute path."""
function draw_heegaard_3d(input::Union{LickorishCode,HeegaardCode}, path::AbstractString;
                         crossing_labels::Bool=false)
    data = _solid_data(input,crossing_labels)
    template = read(joinpath(@__DIR__,"solid3d.html"),String)
    target = abspath(path)
    mkpath(dirname(target))
    write(target,replace(template,"__SOLID_CHART__"=>_solid_svg(data)))
    target
end
