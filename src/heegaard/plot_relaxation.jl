"""Joint curve layout. Fixed ports, continuous collision checks, no geometry in the core."""

function _short_endpoint_guards(points, size)
    length(points) < 2 && return copy(points)
    a, b = first(points), last(points)
    start = _vadd(a, _vscale(min(size, _vnorm(_vsub(points[2], a))/2),
                            _unit(_vsub(points[2], a))))
    finish = _vadd(b, _vscale(min(size, _vnorm(_vsub(points[end-1], b))/2),
                             _unit(_vsub(points[end-1], b))))
    vcat([a, start], points[2:end-1], [finish, b])
end

# Open a disjoint collar around each puncture by a monotone radial map.
# All strands move together: changing only an incident tangent would cut
# through strands that already run along the boundary. The ports stay fixed;
# radial spokes connect them to the mapped exterior. Supports are disjoint
# and d(R)/d(r) >= 1-2*height/width > 0, preserving the exterior isotopy class.
function _sampled_normal_exit_paths(layout, paths; tolerance)
    disks=layout.punctures
    collars=[begin
        gap=minimum((_vnorm(_vsub(d.center,e.center))-d.radius-e.radius
                     for e in disks if e.id!=d.id);init=Inf)
        width=min(0.9d.radius,0.45gap)
        width > 0 || throw(ArgumentError("normal exits require disjoint puncture disks"))
        (disk=d,width=width,height=min(0.08d.radius,0.4width))
    end for d in disks]
    function warp(p)
        for c in collars
            v=_vsub(p,c.disk.center); r=_vnorm(v)
            s=max(0.0,r-c.disk.radius)
            s >= c.width && continue
            return _vadd(p,_vscale(c.height*(1-s/c.width)^2,_unit(v)))
        end
        p
    end
    function edge!(out,a,b,depth=0)
        wa,wb=warp(a),warp(b)
        midpoint=_lerp(a,b,0.5)
        # Quarter points detect a segment entering and leaving a support.
        error=maximum(_vnorm(_vsub(warp(_lerp(a,b,t)),_lerp(wa,wb,t)))
                      for t in (0.25,0.5,0.75))
        maxstep=minimum(c.width for c in collars)/3
        if (error > tolerance || _vnorm(_vsub(b,a))>maxstep) && depth<18
            edge!(out,a,midpoint,depth+1); edge!(out,midpoint,b,depth+1)
        else
            push!(out,wb)
        end
    end
    transformed=[begin
        out=[first(p),warp(first(p))]
        for i in 1:length(p)-1; edge!(out,p[i],p[i+1]); end
        push!(out,last(p))
        out
    end for p in paths]
    transformed,minimum(c.height for c in collars)
end

# Carry every strand through the SAME piecewise-affine homeomorphism. Exact
# segment splitting is essential for tightly packed combinatorial seeds:
# separately approximating a smooth radial map can introduce new crossings.
function _normal_exit_paths(layout, paths; tolerance)
    disks=layout.punctures
    collars=[begin
        euclidean_gap=minimum((_vnorm(_vsub(d.center,e.center))-d.radius-e.radius
            for e in disks if e.id!=d.id);init=Inf)
        square_gap=minimum((max(abs(d.center[1]-e.center[1]),abs(d.center[2]-e.center[2]))-
            e.radius-d.radius for e in disks if e.id!=d.id);init=Inf)
        width=min(0.9d.radius,0.45euclidean_gap,0.45square_gap)
        (disk=d,width=width,height=min(0.08d.radius,0.4width))
    end for d in disks]
    # Disjoint circles need not have disjoint bounding squares. The checked
    # caller refines the radial approximation for that less common placement.
    any(c->c.width<=0,collars) && return _sampled_normal_exit_paths(layout,paths;tolerance)
    mapped=[Route(r.family,r.curve_id,r.segment_id,r.start_port,r.end_port,paths[i])
        for (i,r) in enumerate(layout.routes)]
    for c in collars
        point_map,lines,_=_pl_radial_square_map(c.disk.center,c.disk.radius,
            c.disk.radius+c.height,c.disk.radius+c.width)
        inside=p->max(abs(p[1]-c.disk.center[1]),abs(p[2]-c.disk.center[2]))<=c.disk.radius+c.width+1e-9
        mapped=_map_routes_on_line_arrangement(mapped,point_map,lines;split_predicate=inside)
    end
    [vcat(paths[i][1:1],r.points,paths[i][end:end]) for (i,r) in enumerate(mapped)],
        minimum(c.height for c in collars)
end

function _exit_angle_degrees(layout)
    maximum((begin
        port=only(p for p in layout.ports if p.id==portid)
        disk=only(d for d in layout.punctures if d.id==port.boundary_id)
        n=_unit(_vsub(port.point,disk.center)); t=_unit(tangent)
        acosd(clamp(n[1]*t[1]+n[2]*t[2],-1.0,1.0))
    end for route in layout.routes for (portid,tangent) in (
        (route.start_port,_vsub(first(route.curve.segments).p1,first(route.curve.segments).p0)),
        (route.end_port,_vsub(last(route.curve.segments).p2,last(route.curve.segments).p3))));init=0.0)
end

# Quadrature masses stay fixed during a solve. Consequently discretization
# density does not act as extra charge, and energy and gradient agree exactly.
function _plot_masses(paths)
    [Float64[(i == 1 ? 0.0 : _vnorm(_vsub(p[i],p[i-1]))/2) +
             (i == length(p) ? 0.0 : _vnorm(_vsub(p[i+1],p[i]))/2)
             for i in eachindex(p)] for p in paths]
end

function _plot_energy_gradient(paths, disks, masses; spacing, spread, tautness, clearance, bending=0.012)
    forces = [fill((0.0,0.0), length(p)) for p in paths]
    energy = 0.0
    # Length, not a spring with the initial detour as its rest length.
    for (r,p) in enumerate(paths), i in 1:length(p)-1
        d = _vsub(p[i+1],p[i]); len = _vnorm(d)
        energy += tautness*len
        f = _vscale(tautness, _unit(d))
        forces[r][i] = _vadd(forces[r][i],f)
        forces[r][i+1] = _vsub(forces[r][i+1],f)
    end
    # Elastic bending: neighboring unit tangents turn gradually, including
    # immediately after the fixed normal exit. This is part of the same
    # objective and time evolution, not a final cosmetic spline operation.
    for (r,p) in enumerate(paths), i in 2:length(p)-1
        a=_vsub(p[i],p[i-1]); b=_vsub(p[i+1],p[i])
        la,lb=_vnorm(a),_vnorm(b)
        min(la,lb) <= 1e-10 && continue
        u,v=_vscale(1/la,a),_vscale(1/lb,b)
        delta=_vsub(u,v); weight=bending/max(masses[r][i],1e-6)
        energy += weight*(delta[1]^2+delta[2]^2)/2
        g=_vscale(weight,delta)
        ga=_vscale(1/la,_vsub(g,_vscale(g[1]*u[1]+g[2]*u[2],u)))
        gb=_vscale(-1/lb,_vsub(g,_vscale(g[1]*v[1]+g[2]*v[2],v)))
        forces[r][i-1]=_vadd(forces[r][i-1],ga)
        forces[r][i]=_vadd(forces[r][i],_vsub(gb,ga))
        forces[r][i+1]=_vsub(forces[r][i+1],gb)
    end
    # Two scales: readable adjacent lanes and a whole-corridor interaction.
    # Both endpoints of every interaction receive the opposite derivative.
    for r in eachindex(paths), s in r:length(paths)
        p,q = paths[r],paths[s]
        for i in eachindex(p), j in eachindex(q)
            r == s && j <= i+4 && continue
            d = _vsub(p[i],q[j]); d2=d[1]^2+d[2]^2
            weight = masses[r][i]*masses[s][j]
            for (sigma, amplitude) in ((0.32spread,spacing),(spread,0.18spacing))
                d2 > 36sigma^2 && continue
                e = amplitude*weight*exp(-d2/(2sigma^2))/sigma
                energy += e
                f = _vscale(e/sigma^2,d)
                forces[r][i] = _vadd(forces[r][i],f)
                forces[s][j] = _vsub(forces[s][j],f)
            end
        end
    end
    # Short-range segment barriers protect intervals between samples as well.
    for r in eachindex(paths), s in r:length(paths)
        p,q=paths[r],paths[s]
        for i in 1:length(p)-1, j in 1:length(q)-1
            r == s && j <= i+3 && continue
            target=0.14spread
            any(max(p[i][k],p[i+1][k])+target < min(q[j][k],q[j+1][k]) ||
                max(q[j][k],q[j+1][k])+target < min(p[i][k],p[i+1][k]) for k in 1:2) && continue
            distance,a,b,u,v = _closest_segment_pair(p[i],p[i+1],q[j],q[j+1])
            distance >= target && continue
            k=20spacing*min(masses[r][i]+masses[r][i+1],masses[s][j]+masses[s][j+1])
            energy += k*(target-distance)^2/2
            f=_vscale(k*(target-distance),_unit(_vsub(a,b)))
            forces[r][i]=_vadd(forces[r][i],_vscale(1-u,f))
            forces[r][i+1]=_vadd(forces[r][i+1],_vscale(u,f))
            forces[s][j]=_vsub(forces[s][j],_vscale(1-v,f))
            forces[s][j+1]=_vsub(forces[s][j+1],_vscale(v,f))
        end
    end
    for (r,p) in enumerate(paths), i in 2:length(p)-2, disk in disks
        closest,t=_closest_on_segment_with_parameter(disk.center,p[i],p[i+1])
        radial=_vsub(closest,disk.center); margin=_vnorm(radial)-disk.radius
        margin >= clearance && continue
        k=35spacing*(masses[r][i]+masses[r][i+1])
        energy += k*(clearance-margin)^2/2
        f=_vscale(k*(clearance-margin),_unit(radial))
        forces[r][i]=_vadd(forces[r][i],_vscale(1-t,f))
        forces[r][i+1]=_vadd(forces[r][i+1],_vscale(t,f))
    end
    energy,forces
end

_lerp(a,b,t)=_vadd(a,_vscale(t,_vsub(b,a)))

# A Lipschitz bound certifies the entire linear motion, not just its endpoint.
# Ambiguous intervals are bisected; exhaustion rejects the step.
function _motion_interval_clear(distance_at, speed, required; lo=0.0, hi=1.0, depth=12)
    mid=(lo+hi)/2
    distance=distance_at(mid)
    distance > required+speed*(hi-lo)/2+1e-10 && return true
    distance <= required+1e-10 && return false
    depth == 0 && return false
    _motion_interval_clear(distance_at,speed,required;lo,hi=mid,depth=depth-1) &&
        _motion_interval_clear(distance_at,speed,required;lo=mid,hi,depth=depth-1)
end

function _plot_motion_clear(old,new,disks)
    velocities=[[ _vnorm(_vsub(a,b)) for (a,b) in zip(p,q)] for (p,q) in zip(old,new)]
    for r in eachindex(old), i in 1:length(old[r])-1
        speed=max(velocities[r][i],velocities[r][i+1])
        if speed > 0
            for disk in disks
                distance_at=t->sqrt(_segment_distance_squared(disk.center,
                    _lerp(old[r][i],new[r][i],t),_lerp(old[r][i+1],new[r][i+1],t)))
                _motion_interval_clear(distance_at,speed,disk.radius+1e-8) || return false
            end
        end
        for s in r:length(old), j in (s==r ? i+2 : 1):length(old[s])-1
            other_speed=max(velocities[s][j],velocities[s][j+1])
            speed+other_speed == 0 && continue
            a,b=old[r][i],old[r][i+1]; c,d=old[s][j],old[s][j+1]
            # Bounding boxes of the whole swept segments.
            ap,bp=new[r][i],new[r][i+1]; cp,dp=new[s][j],new[s][j+1]
            any(max(a[k],b[k],ap[k],bp[k])+1e-7 < min(c[k],d[k],cp[k],dp[k]) ||
                max(c[k],d[k],cp[k],dp[k])+1e-7 < min(a[k],b[k],ap[k],bp[k]) for k in 1:2) && continue
            distance_at=t->_segment_distance(_lerp(a,ap,t),_lerp(b,bp,t),_lerp(c,cp,t),_lerp(d,dp,t))
            _motion_interval_clear(distance_at,speed+other_speed,1e-7) || return false
        end
    end
    true
end

function _project_plot_directions!(directions,paths,disks,scale)
    mobile(r,i)=3 <= i <= length(paths[r])-2
    # Project normal velocities at active contacts. A contact is allowed to
    # slide tangentially; it must not freeze unrelated points in the diagram.
    for pass in 1:6
        for r in eachindex(paths), i in 2:length(paths[r])-2, disk in disks
            a,t=_closest_on_segment_with_parameter(disk.center,paths[r][i],paths[r][i+1])
            radial=_vsub(a,disk.center)
            _vnorm(radial)-disk.radius > 0.025scale && continue
            n=_unit(radial); u=mobile(r,i) ? 1-t : 0.0; v=mobile(r,i+1) ? t : 0.0
            velocity=_vadd(_vscale(u,directions[r][i]),_vscale(v,directions[r][i+1]))
            inward=velocity[1]*n[1]+velocity[2]*n[2]
            inward >= 0 && continue
            denominator=u^2+v^2
            denominator < 1e-12 && continue
            correction=_vscale(-inward/denominator,n)
            directions[r][i]=_vadd(directions[r][i],_vscale(u,correction))
            directions[r][i+1]=_vadd(directions[r][i+1],_vscale(v,correction))
        end
        for r in eachindex(paths), s in r:length(paths)
            p,q=paths[r],paths[s]
            for i in 1:length(p)-1, j in (s==r ? i+3 : 1):length(q)-1
                target=0.018scale
                any(max(p[i][k],p[i+1][k])+target < min(q[j][k],q[j+1][k]) ||
                    max(q[j][k],q[j+1][k])+target < min(p[i][k],p[i+1][k]) for k in 1:2) && continue
                distance,a,b,u,v=_closest_segment_pair(p[i],p[i+1],q[j],q[j+1])
                distance > target && continue
                n=_unit(_vsub(a,b))
                ids=((r,i),(r,i+1),(s,j),(s,j+1))
                weights=(mobile(r,i) ? 1-u : 0.0,mobile(r,i+1) ? u : 0.0,
                         mobile(s,j) ? v-1 : 0.0,mobile(s,j+1) ? -v : 0.0)
                closing=sum(w*(directions[r0][i0][1]*n[1]+directions[r0][i0][2]*n[2])
                            for ((r0,i0),w) in zip(ids,weights))
                closing >= 0 && continue
                denominator=sum(w^2 for w in weights)
                denominator < 1e-12 && continue
                for ((r0,i0),w) in zip(ids,weights)
                    directions[r0][i0]=_vadd(directions[r0][i0],_vscale(-w*closing/denominator,n))
                end
            end
        end
    end
    directions
end

function _vertex_motion_clear(paths,r,i,destination,disks)
    origin=paths[r][i]; speed=_vnorm(_vsub(destination,origin))
    for segment in (i-1,i)
        a,b=paths[r][segment],paths[r][segment+1]
        ap=segment==i ? destination : a
        bp=segment+1==i ? destination : b
        for disk in disks
            at=t->sqrt(_segment_distance_squared(disk.center,_lerp(a,ap,t),_lerp(b,bp,t)))
            _motion_interval_clear(at,speed,disk.radius+1e-9) || return false
        end
        for s in eachindex(paths), j in 1:length(paths[s])-1
            s==r && abs(j-segment)<=1 && continue
            c,d=paths[s][j],paths[s][j+1]
            any(max(a[k],b[k],ap[k],bp[k])+1e-7 < min(c[k],d[k]) ||
                max(c[k],d[k])+1e-7 < min(a[k],b[k],ap[k],bp[k]) for k in 1:2) && continue
            at=t->_segment_distance(_lerp(a,ap,t),_lerp(b,bp,t),c,d)
            _motion_interval_clear(at,speed,1e-7) || return false
        end
    end
    true
end

function _local_plot_sweep(paths,forces,disks,scale,factor)
    proposal=deepcopy(paths)
    for r in eachindex(paths), i in 3:length(paths[r])-2
        force=forces[r][i]; magnitude=_vnorm(force)
        magnitude < 1e-8 && continue
        step=factor*min(0.08,0.018scale/magnitude)
        for trial in 1:14
            destination=_vadd(paths[r][i],_vscale(step,force))
            if _vertex_motion_clear(proposal,r,i,destination,disks)
                proposal[r][i]=destination
                break
            end
            step/=2
        end
    end
    proposal
end

function _elastic_plot_metric(points,masses;tautness,bending)
    n=length(points)-4
    matrix=Matrix{Float64}(I,n,n)
    function addterm(ids,weights,strength)
        for (i,a) in zip(ids,weights), (j,b) in zip(ids,weights)
            3<=i<=length(points)-2 && 3<=j<=length(points)-2 || continue
            matrix[i-2,j-2]+=strength*a*b
        end
    end
    for i in 1:length(points)-1
        len=max(_vnorm(_vsub(points[i+1],points[i])),1e-6)
        addterm((i,i+1),(-1.,1.),tautness/len)
    end
    for i in 2:length(points)-1
        a=max(_vnorm(_vsub(points[i],points[i-1])),1e-6)
        b=max(_vnorm(_vsub(points[i+1],points[i])),1e-6)
        addterm((i-1,i,i+1),(-1/a,1/a+1/b,-1/b),bending/max(masses[i],1e-6))
    end
    cholesky(Symmetric(matrix))
end

function _joint_plot_relaxation(H,layout; iterations=500)
    radius=minimum(p.radius for p in layout.punctures;init=0.28)
    scale=radius/0.28
    initial=[_short_endpoint_guards(_route_points(r),0.012scale) for r in layout.routes]
    @info "Heegaard plot: shortening feasible routes"
    shortened,shortcuts=_global_shortcut_polylines(H,layout,initial)
    tightened,geodesics=_geodesically_tighten_polylines(H,layout,shortened;
        clearance=0.12scale,route_clearance=0.002scale)
    # Normal exits are boundary conditions BEFORE the mechanical solve.
    # Never attach a straight stub or radially warp a converged result.
    normal_paths,collar_length=_normal_exit_paths(layout,tightened;tolerance=1e-4scale)
    if _validated_polyline_candidate(H,layout,normal_paths)===nothing
        separation=_route_pair_clearance(tightened)
        tolerance=min(1e-5scale,0.01separation)
        accepted_normal=false
        for _ in 1:8
            normal_paths,collar_length=_sampled_normal_exit_paths(layout,tightened;tolerance)
            if _validated_polyline_candidate(H,layout,normal_paths)!==nothing
                accepted_normal=true;break
            end
            tolerance/=4
        end
        accepted_normal || throw(ArgumentError("normal-exit construction failed geometric validation before relaxation"))
    end
    normal_paths,_=_simplify_polylines(H,layout,normal_paths;tolerance=0.004scale)
    current=[_densify_route(p;max_spacing=0.09scale) for p in normal_paths]
    @info "Heegaard plot: joint relaxation" points=sum(length,current)
    masses=_plot_masses(current)
    params=(spacing=2.2,spread=0.36scale,tautness=4.0,clearance=0.14scale,bending=0.04scale^2)
    energy,forces=_plot_energy_gradient(current,layout.punctures,masses;params...)
    history=[energy]; accepted=0; stop="iteration-limit"
    # A positive-definite elastic metric couples the whole route and accounts
    # for short elements' bending stiffness. It prevents one tiny endpoint
    # element from forcing microscopic steps everywhere else.
    for iteration in 1:iterations
        directions=[fill((0.0,0.0),length(p)) for p in current]
        for r in eachindex(current)
            ids=3:length(current[r])-2
            isempty(ids) && continue
            metric=_elastic_plot_metric(current[r],masses[r];tautness=params.tautness,bending=params.bending)
            fx=metric\[forces[r][i][1] for i in ids]
            fy=metric\[forces[r][i][2] for i in ids]
            for (k,i) in enumerate(ids); directions[r][i]=(fx[k],fy[k]); end
        end
        _project_plot_directions!(directions,current,layout.punctures,scale)
        maxforce=maximum((_vnorm(f) for fs in directions for f in fs);init=0.0)
        maxforce < 1e-5 && (stop="stationary";break)
        step=min(0.5,0.045scale/maxforce)
        moved=false
        for attempt in 1:18
            proposal=[[_vadd(p[i],_vscale(step,directions[r][i])) for i in eachindex(p)]
                      for (r,p) in enumerate(current)]
            nextenergy,nextforces=_plot_energy_gradient(proposal,layout.punctures,masses;params...)
            if nextenergy < energy-1e-10 && _plot_motion_clear(current,proposal,layout.punctures)
                current=proposal; energy=nextenergy; forces=nextforces
                push!(history,energy); accepted+=1; moved=true
                break
            end
            step/=2
        end
        if !moved || step*maxforce < 2e-6scale || iteration%10==0
            localmoved=false
            for retry in 0:5
                proposal=_local_plot_sweep(current,directions,layout.punctures,scale,0.5^retry)
                nextenergy,nextforces=_plot_energy_gradient(proposal,layout.punctures,masses;params...)
                if nextenergy < energy-1e-8
                    current=proposal; energy=nextenergy; forces=nextforces
                    push!(history,energy); accepted+=1; localmoved=true
                    break
                end
            end
            !moved && !localmoved && (stop="line-search-stalled";break)
        end
        iteration%100==0 && @info "Heegaard curve relaxation" iteration energy
    end
    # Do not run another length-only shortest path after establishing spacing.
    simplified,removed=_simplify_polylines(H,layout,current;tolerance=0.001scale)
    separation=_route_pair_clearance(simplified)
    requested=min(0.025scale,0.85separation)
    result,handles,certificate=_fit_bezier_layout(H,layout,simplified;
        handle_scale=1/3,flatness=1e-4scale,max_depth=20,minimum_route_separation=requested)
    exit_angle=_exit_angle_degrees(result)
    exit_angle <= 15 || throw(ArgumentError("normal exit angle exceeds 15 degrees"))
    trace=deepcopy(result.backend_trace)
    trace["smoothing"]=Dict{String,Any}(
        "method"=>"elastic-normal-boundary",
        "accepted_steps"=>accepted,"stop_reason"=>stop,
        "converged"=>stop=="stationary", "energy_history"=>history,
        "initial_energy"=>first(history),"final_energy"=>last(history),
        "length_before"=>sum(_polyline_length,initial),
        "length_after"=>sum(r->_polyline_length(r.points),result.routes),
        "minimum_flattened_separation"=>isfinite(certificate.minimum_separation) ? certificate.minimum_separation : nothing,
        "minimum_flattened_puncture_margin"=>certificate.minimum_puncture_margin,
        "flatness_bound"=>certificate.tolerance,"endpoint_tangents_certified"=>true,
        "parameters"=>Dict(string(k)=>v for (k,v) in pairs(params)),
        "shortcut_count"=>shortcuts.accepted,"geodesic_route_count"=>geodesics.accepted,
        "bezier_segments"=>sum(length(r.curve.segments) for r in result.routes),
        "continuous_relaxation_motion_checked"=>true,
        "exit_angle_limit_degrees"=>15.0,"maximum_exit_angle_degrees"=>exit_angle,
        "minimum_normal_exit_length"=>collar_length)
    trace["plot"]=Dict("profile"=>"fixed-100","review_slider_values"=>Dict(
        "spacing"=>100,"spread_range"=>100,"tautness"=>100,"puncture_clearance"=>100),
        "implementation"=>"joint curve solver; not the retired browser preview")
    HeegaardLayout(result.diagram_hash,result.view,result.chart,result.punctures,result.ports,
        result.routes,result.crossings,trace,result.report)
end
