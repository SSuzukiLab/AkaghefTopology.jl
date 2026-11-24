# cf: lib/topology.jl/docs/specs/drawing-objects.md §6 — pass placed parts to Makie, unchanged
module TopologyMakieExt

using Makie
import Topology.HeegaardCodes: drawingplot, drawingplot!, HeegaardDrawing, resolved_style,
    element_at, drawing_inspector_label, drawing_label

_emit(prim) = begin
    pts = Point2f[Point2f(p[1], p[2]) for p in prim.points]
    prim.closed && !isempty(pts) && push!(pts, pts[1])
    pts
end

_anchor(el, prim) = begin
    p = el.drawn.labelpoint === nothing ? prim.points[1] : el.drawn.labelpoint
    Point2f(p[1] + 7, p[2] + 7)
end

function _label(d, pid, npoints)
    (self, index, position) -> begin
        idx = index isa Integer ? Int(index) : Int(first(index))
        drawing_inspector_label(d, element_at(d, pid, clamp(idx, 1, npoints)))
    end
end

function drawingplot(d::HeegaardDrawing)
    fig = Figure(; backgroundcolor="#fcfbf7", size=(720, 720))
    ax = Axis(fig[1, 1]; aspect=DataAspect(), backgroundcolor="#fcfbf7")
    hidedecorations!(ax)
    hidespines!(ax)
    drawingplot!(ax, d)
    DataInspector(fig)
    fig
end

function drawingplot!(ax, d::HeegaardDrawing)
    ax.yreversed[] = true
    order = sort!(collect(eachindex(d.elements)); by=i -> resolved_style(d, d.elements[i])[:zorder])
    styles = Dict{Int,Observable}()
    for i in order
        el = d.elements[i]
        st = styles[i] = Observable(resolved_style(d, el))
        vis = lift(s -> s[:visible], st)
        col = lift(s -> s[:color], st)
        width = lift(s -> Float32(s[:linewidth]), st)
        fill = lift(s -> s[:fill], st)
        marker = lift(s -> s[:marker], st)
        msize = lift(s -> s[:markersize], st)
        linestyle = lift(s -> something(s[:linestyle], :solid), st)
        for pid in el.primitives
            prim = d.primitives[pid]
            pts = Observable(_emit(prim))
            label = _label(d, pid, length(prim.points))
            if prim.shape == :polygon
                poly!(ax, pts; color=fill, strokecolor=col, strokewidth=width, visible=vis,
                    inspectable=true, inspector_label=label)
            elseif prim.shape == :point
                scatter!(ax, pts; color=col, marker=marker, markersize=msize, visible=vis,
                    inspectable=true, inspector_label=label)
            else
                lines!(ax, pts; color=col, linewidth=width, linestyle=linestyle, visible=vis,
                    inspectable=true, inspector_label=label)
            end
        end
        shown = el.kind == :crossing ? [d.primitives[pid] for pid in el.primitives] : [d.primitives[first(el.primitives)]]
        for prim in shown
            pos = Observable([_anchor(el, prim)])
            words = lift(s -> something(drawing_label(el, s[:label]), ""), st)
            text!(ax, pos; text=words, fontsize=lift(s -> s[:labelsize], st), color="#222222",
                visible=lift(s -> s[:visible] && drawing_label(el, s[:label]) !== nothing, st),
                inspectable=false)
        end
    end
    push!(d.style.listeners, () -> begin
        for (i, st) in styles
            st[] = resolved_style(d, d.elements[i])
        end
    end)
    ax
end

end
