# cf: docs/specs/solid3d.md — arclength averaging on the public planar chart.
# The primitive integrates each straight segment exactly. Symmetric averaging
# has continuous tangent (P(s+r)-P(s-r))/(2r), even at a vertex of P.
function _solid_smooth_run(run, width)
    points = HCPlotPoint[first(run.points)]
    for q in run.points
        _hp_dist(q, last(points)) > 1e-10 && push!(points, q)
    end
    if run.closed
        _hp_dist(first(points), last(points)) <= 1e-10 && pop!(points)
        push!(points, first(points))
    end
    length(points) <= 2 && return run
    if !run.closed && all(q -> _hp_nearest(q, first(points), last(points))[1] < 1e-9, points)
        return run
    end
    stations = Float64[0]
    integrals = HCPlotPoint[(0., 0.)]
    for i in 2:length(points)
        d = _hp_dist(points[i-1], points[i])
        push!(stations, last(stations) + d)
        area = _hp_mul(_hp_add(points[i-1], points[i]), d/2)
        push!(integrals, _hp_add(last(integrals), area))
    end
    len = last(stations)
    radius = min(width, len / 8)
    function primitive(t)
        i = clamp(searchsortedlast(stations, t), 1, length(points)-1)
        h = t-stations[i]
        q = _hp_lerp(points[i], points[i+1], h/(stations[i+1]-stations[i]))
        _hp_add(integrals[i], _hp_mul(_hp_add(points[i],q), h/2))
    end
    function extended(t)
        if run.closed
            return _hp_add(primitive(mod(t,len)), _hp_mul(last(integrals),floor(t/len)))
        elseif t < 0
            return _hp_add(_hp_mul(first(points),2t), primitive(-t))
        elseif t > len
            return _hp_add(_hp_mul(last(points),2(t-len)), primitive(2len-t))
        end
        primitive(t)
    end
    average(t) = _hp_mul(_hp_sub(extended(t+radius), extended(t-radius)), 1/(2radius))
    n = max(12, ceil(Int, len/min(3., radius/5)))
    samples = [average(len*i/n) for i in 0:(run.closed ? n-1 : n)]
    if !run.closed
        samples[1], samples[end] = first(run.points), last(run.points)
    end
    _hp_run(run, samples)
end

function _solid_smooth_chart(result, lc)
    original = result.placement
    width = isempty(original.holes) ? min(original.width,original.height)/12 :
            minimum(r for (_,(_,r)) in original.holes)*0.8
    has_empty = any(isempty, values(lc.words))
    function checked(runs)
        placement = _hp_rebuild(original, runs)
        audit = _hp_audit(placement)
        audit.ok || return nothing
        has_empty && !_hp_empty_guard(original,placement,lc) && return nothing
        (;placement, audit)
    end
    for reduction in 0:5
        radius = width / 2^reduction
        runs = [_solid_smooth_run(run, radius) for run in original.runs]
        chart = checked(runs)
        isnothing(chart) && continue
        # A single tight crossing should not limit the radius on every other arc.
        for i in eachindex(runs)
            runs[i] === original.runs[i] && continue
            for k in 0:reduction-1
                candidate = copy(runs)
                candidate[i] = _solid_smooth_run(original.runs[i], width / 2^k)
                trial = checked(candidate)
                isnothing(trial) && continue
                runs, chart = candidate, trial
                break
            end
        end
        return merge(result, chart)
    end
    throw(ArgumentError("could not smooth the solid chart while preserving its signed crossings"))
end
