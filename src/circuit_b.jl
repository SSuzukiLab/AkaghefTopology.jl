"""
    circuit_b_boundary(v)

Return each directed open circuit-B path with its endpoint lane identities,
ordered signed edge occurrences, dots and incidence row. Boundary identity is
`(strand, side, lane)`; row numbers are local to this exact presentation.
This is a combinatorial query, not a relative-framing certificate.
"""
function circuit_b_boundary(v::VirtualLink)
    disks=disk_table(v)
    [(row=row, start=disks.boundary_start[row], stop=disks.boundary_end[row],
      nodes=copy(disks.paths[row]), edges_s=copy(disks.edges_s[row]),
      dots=disks.Ndots[row], coefficients=copy(disks.delta[row,:]))
     for row in eachindex(disks.ID) if !disks.isCycle[row]]
end

"""
    circuit_b_at(v, strand, side; lane=nothing)

Query circuit B at a particular strand endpoint (`side=:in` or `:out`).
With no lane, return the three lane incidences in lane order. Each result
records whether the directed B path starts or ends at that endpoint, its other
endpoint, and its full ordered edge occurrences. A closed strand has no such
endpoint and is rejected rather than silently cut.
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
    plot_circuit_b_svg(v, path)

Draw the actual ordered circuit-B paths as rows. Open endpoints and closed
cycles are distinguished explicitly. Edge labels keep every signed occurrence
(including repeated edges); this is a combinatorial diagram, not a planar
embedding of the o-graph.
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
