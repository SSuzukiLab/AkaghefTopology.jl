include("solid3d_smoothing.jl")

# cf: docs/specs/heegaard-plot.md and solid3d.md — reuse the full reference chart.
"""
    plot(lc::LickorishCode; infinity_face=nothing) -> NamedTuple

Produce a planar SVG of a role-free LC while retaining its reference curves.

# Arguments
The LC must pass `validate_lc_input`, passage consistency and surface validation; their exceptions propagate. Optional `infinity_face` is interpreted by [the planar chart](../specs/heegaard-plot.md). Chart and placement failures propagate.

# Returns
A named tuple with chart data, `placement`, `audit`, `success` and `svg`. The success flag equals `audit.ok`.

# Guarantees
1. The output retains the LC reference runs in its planar chart and reports the geometric audit of that placement. [proved/agent: by construction]
2. The LC is not modified. [proved/agent: by construction]

# Not guaranteed
This overload does not assign alpha/beta roles or certify a Heegaard splitting. A returned record can have `success=false`; it does not open a viewer or write a file.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-plot.md — reference-chart geometry and acceptance.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1], "m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing); r = H.plot(lc); r.success == r.audit.ok
true
```

Role: draw
"""
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

"""
    draw_heegaard_3d(input::Union{LickorishCode,HeegaardCode}, path::AbstractString; crossing_labels=false) -> String

Write an interactive solid-chart display of an LC or resolved Heegaard diagram.

# Arguments
`input` belongs to the [solid-chart domain](../specs/solid3d.md).
Checked: nonpositive genus raises `ArgumentError` requiring positive genus;
failed planar geometry raises `ArgumentError`.
`path` must be writable. `crossing_labels::Bool` controls crossing labels.
Diagram validation, template-reading and filesystem errors propagate.

# Returns
The absolute path to the written HTML file, with its embedded chart.

# Guarantees
1. On success the HTML contains the full reference chart associated with the
   input in the linked solid-display convention. [proved/agent: by construction]
2. The input is not modified. Parent directories are created as needed and
   an existing destination is overwritten. [proved/agent: by construction]

# Not guaranteed
The displayed folding is not a triangulation or a certificate of the represented
3-manifold. Successful drawing is not guaranteed for every valid input.

# Sources
cf: lib/topology.jl/docs/specs/solid3d.md — full reference chart and interactive HTML display.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.draw_heegaard_3d(H.LickorishCode(0, Dict{String,Vector{Int}}(), Dict(), Dict(), Dict(), Dict(), nothing), "unused.html") catch e; e isa ArgumentError end
true
```

Role: draw
"""
function draw_heegaard_3d(input::Union{LickorishCode,HeegaardCode}, path::AbstractString;
                         crossing_labels::Bool=false)
    data = _solid_data(input,crossing_labels)
    template = read(joinpath(@__DIR__,"solid3d.html"),String)
    target = abspath(path)
    mkpath(dirname(target))
    write(target,replace(template,"__SOLID_CHART__"=>_solid_svg(data)))
    target
end
