"""
    circuit_b_boundary(v::VirtualLink) -> Vector

Return the directed open paths of the presentation's Method-B circuit model.

# Arguments
`v` has the coherent real-edge and cut-marker data required by `disk_table`
in the [data model](../specs/virtual-link-data-model.md#21-open-strands-and-cut-markers).
That query checks malformed markers and inconsistent path incidences with
`ArgumentError`; invalid indices or signs can propagate `BoundsError` or
`KeyError`. The caller supplies the geometric interpretation of the data.

# Returns
One record per open row: `row`, directed endpoints `start` and `stop`, copied
node path `nodes`, signed occurrences `edges_s`, dot count `dots` and copied
incidence vector `coefficients`. Their meanings are the
[Method-B model](../specs/ograph-move-verifier.md#ps-gate) and
[relative endpoints](../specs/spin-weights.md#relative-boundary-data).

# Guarantees
1. The records are exactly the noncycle rows of `disk_table(v)`, in its row order. [proved/agent: by construction]
2. The input is not modified. Returned path, occurrence and coefficient vectors are copies. [proved/agent: by construction]

# Not guaranteed
These paths are presentation coordinates, not a boundary-framing or spin
extension certificate. A closed-normal manifold interpretation is not checked.

# Sources
cf:lib/topology.jl/docs/specs/ograph-move-verifier.md PS gate — the Method-B local-node model and source comparison.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[0]], orientation=Int[]);

julia> length(circuit_b_boundary(v))
3
```

Role: access
"""
function circuit_b_boundary(v::VirtualLink)
    disks=disk_table(v)
    [(row=row, start=disks.boundary_start[row], stop=disks.boundary_end[row],
      nodes=copy(disks.paths[row]), edges_s=copy(disks.edges_s[row]),
      dots=disks.Ndots[row], coefficients=copy(disks.delta[row,:]))
     for row in eachindex(disks.ID) if !disks.isCycle[row]]
end

"""
    circuit_b_at(v::VirtualLink, strand::Integer, side::Symbol; lane=nothing) -> Vector{NamedTuple}

Return the directed Method-B paths incident to selected lanes of a cut strand endpoint.

# Arguments
Caller: E1 of [circuit endpoint selection](../specs/virtual-link-data-model.md#26-circuit-endpoint-selection).
Checked: E2–E4; strand bounds raise `BoundsError`, invalid selectors or an
uncut strand raise `ArgumentError`, and nonunique directed incidence raises
`ErrorException`. Selector errors identify the accepted selection; the
incidence error states that exactly one directed B path is required.
Other `disk_table` errors propagate.

# Returns
The selected `circuit_b_boundary` records extended by queried `endpoint`,
`direction` and the `other` endpoint, in increasing lane order.

# Guarantees
1. Each returned path meets its queried endpoint. Direction is `:outgoing` when that endpoint is the path start and `:incoming` otherwise. [proved/agent: by construction]
2. Omitting the lane returns three records; selecting one returns one record. [proved/agent: by construction]
3. The input is unchanged. [proved/agent: by construction]

# Not guaranteed
This does not cut a closed strand, select a boundary trivialization or certify
a relative framing; it queries the fixed presentation's local circuit model.

# Sources
cf:lib/topology.jl/docs/specs/ograph-move-verifier.md PS gate — local Method-B paths and endpoint lanes.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[0]], orientation=Int[]);

julia> Tuple(p.direction for p in circuit_b_at(v,1,:in))
(:outgoing, :outgoing, :incoming)
```

Role: access
"""
function circuit_b_at(v::VirtualLink,strand::Integer,side::Symbol;lane=nothing)
    1<=strand<=length(v.gauss) || throw(BoundsError(v.gauss,strand))
    side in (:in,:out) || throw(ArgumentError("side must be :in or :out"))
    strand_cuts(v)[strand] || throw(ArgumentError("the selected strand has no cut endpoint"))
    lane===nothing || (lane isa Integer && lane in 1:3) ||
        throw(ArgumentError("lane must be 1, 2, 3 or nothing"))
    paths=circuit_b_boundary(v)
    result=NamedTuple[]
    for selected_lane in (lane===nothing ? (1:3) : (lane:lane))
        endpoint=(strand=Int(strand),side=side,lane=Int(selected_lane))
        matches=[path for path in paths if path.start==endpoint || path.stop==endpoint]
        length(matches)==1 || error("boundary lane must belong to exactly one directed B path")
        path=only(matches)
        push!(result,merge(path,(endpoint=endpoint,
            direction=path.start==endpoint ? :outgoing : :incoming,
            other=path.start==endpoint ? path.stop : path.start)))
    end
    result
end

_circuit_endpoint_label(e)=e===nothing ? "closed" : "S$(e.strand) $(e.side) / lane $(e.lane)"

"""
    plot_circuit_b_svg(v::VirtualLink, path::AbstractString) -> AbstractString

Write a row drawing of the directed Method-B paths and cycles as SVG and return the supplied path.

# Arguments
`v` satisfies the `disk_table` domain linked by `circuit_b_boundary`.
`path` names a writable destination whose parent exists; no directory is
created. Path opening, writing and disk-table errors propagate. An existing
file is overwritten, and a write failure can leave a partial file.

# Returns
The supplied `path`, after writing the
[circuit-B view](../specs/virtual-link-plot.md#7-real-arc-aggregation).

# Guarantees
1. The drawing has one row per disk-table path/cycle and retains its ordered signed edge occurrences and dot count. Open endpoints carry their stored lane labels. [proved/agent: by construction]
2. The output file changes; the input presentation is not modified. [proved/agent: by construction]

# Not guaranteed
The row drawing is not a planar o-graph embedding, a geometric disk diagram
or a relative-framing certificate.

# Sources
cf:lib/topology.jl/docs/specs/virtual-link-plot.md §7 — the combinatorial circuit-B view.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[0]], orientation=Int[]);

julia> mktempdir() do dir
           path = plot_circuit_b_svg(v,joinpath(dir,"interval.svg"))
           occursin("S1 in / lane 3",read(path,String))
       end
true
```

Role: draw
"""
function plot_circuit_b_svg(v::VirtualLink,path::AbstractString)
    d=disk_table(v)
    n=length(d.ID)
    max_edges=maximum(length.(d.edges_s);init=0)
    width=max(760,280+80max_edges)
    height=max(130,80+100n)
    open(path,"w") do io
        println(io,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io,"<title>Circuit B: directed paths and boundary lanes</title><rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
        println(io,"<defs><marker id=\"b-arrow\" viewBox=\"0 0 10 10\" refX=\"9\" refY=\"5\" markerWidth=\"6\" markerHeight=\"6\" orient=\"auto-start-reverse\"><path d=\"M0 0 L10 5 L0 10 Z\" fill=\"#244a70\"/></marker></defs>")
        println(io,"<g font-family=\"sans-serif\" font-size=\"13\" fill=\"#18334d\"><text x=\"24\" y=\"28\" font-size=\"18\">Circuit B — ordered signed edge occurrences</text>")
        for row in eachindex(d.ID)
            y=70+100(row-1)
            closed=d.isCycle[row]
            println(io,"<g data-b-row=\"$row\" data-cycle=\"$closed\"><text x=\"24\" y=\"$y\">B$row · dots $(d.Ndots[row]) · $(closed ? "cycle" : "path")</text>")
            x1=190; x2=width-200
            println(io,"<path d=\"M$x1 $(y+24) H$x2\" stroke=\"#244a70\" stroke-width=\"2\" fill=\"none\" marker-end=\"url(#b-arrow)\"/>")
            if closed
                println(io,"<path d=\"M$x2 $(y+24) V$(y+47) H$x1 V$(y+24)\" stroke=\"#244a70\" fill=\"none\" marker-end=\"url(#b-arrow)\"/>")
            else
                println(io,"<circle cx=\"$x1\" cy=\"$(y+24)\" r=\"5\" fill=\"white\" stroke=\"#b54d28\"/><circle cx=\"$x2\" cy=\"$(y+24)\" r=\"5\" fill=\"white\" stroke=\"#b54d28\"/>")
                println(io,"<text x=\"$x1\" y=\"$(y+62)\">$(_circuit_endpoint_label(d.boundary_start[row]))</text><text x=\"$x2\" y=\"$(y+62)\" text-anchor=\"end\">$(_circuit_endpoint_label(d.boundary_end[row]))</text>")
            end
            occurrences=d.edges_s[row]
            for (i,e) in enumerate(occurrences)
                x=x1+(x2-x1)*(i-0.5)/max(1,length(occurrences))
                println(io,"<text data-edge=\"$(abs(e))\" x=\"$x\" y=\"$(y+15)\" text-anchor=\"middle\">$(e>0 ? "+" : "−")E$(abs(e))</text>")
            end
            println(io,"</g>")
        end
        println(io,"</g></svg>")
    end
    path
end
