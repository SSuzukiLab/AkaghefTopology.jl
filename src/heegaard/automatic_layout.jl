"""Move terminal rays to exact polar ports in disjoint annular collars.

At each radius the lifted angles interpolate two strictly ordered sequences.
The PL approximation is subsequently checked by validate_layout, fail-closed.
"""
function _polar_terminal_routes(cut, punctures, ports, routes; samples=64)
    polar=place_ports(cut,punctures)
    collars=Dict{Int,Vector{NTuple{2,Float64}}}()
    for disk in punctures
        local_ports=sort([p for p in ports if p.boundary_id==disk.id];by=p->p.parameter)
        isempty(local_ports) && continue
        angles=[atan(p.point[2]-disk.center[2],p.point[1]-disk.center[1]) for p in local_ports]
        for i in 2:length(angles)
            while angles[i]<=angles[i-1]; angles[i]+=2π; end
        end
        angles[end]-angles[1]<2π || throw(ArgumentError("terminal rays reverse polar rotation"))
        target=2π .* [p.parameter for p in local_ports]
        shift=2π*round(sum(target.-angles)/(2π*length(angles)))
        angles .+= shift
        for (p,a,b) in zip(local_ports,target,angles)
            # First segment is radial, preserving strict exteriority at r=R.
            points=[(disk.center[1]+disk.radius*cos(a),disk.center[2]+disk.radius*sin(a))]
            for i in 0:samples
                t=i/samples; r=disk.radius*(1.1+1.9t); theta=(1-t)*a+t*b
                push!(points,(disk.center[1]+r*cos(theta),disk.center[2]+r*sin(theta)))
            end
            exact=only(q.point for q in polar if q.id==p.id)
            points[1]=exact
            collars[p.id]=points
        end
    end
    fitted=[Route(r.family,r.curve_id,r.segment_id,r.start_port,r.end_port,
        vcat(collars[r.start_port],r.points[2:end-1],reverse(collars[r.end_port]))) for r in routes]
    polar,fitted
end

"""Fixture-free combinatorial planar realization; geometry failures are explicit."""
function _compact_atlas_exterior(punctures,routes)
    bounds=[(minimum(p.center[k]-p.radius for p in punctures)-0.2,
             maximum(p.center[k]+p.radius for p in punctures)+0.2) for k in 1:2]
    extents=[extrema(p[k] for r in routes for p in r.points) for k in 1:2]
    scales=[(min(1.0,1.0/max(bounds[k][1]-extents[k][1],1e-12)),
             min(1.0,1.0/max(extents[k][2]-bounds[k][2],1e-12))) for k in 1:2]
    coordinate(x,k)=x<bounds[k][1] ? bounds[k][1]+scales[k][1]*(x-bounds[k][1]) :
        x>bounds[k][2] ? bounds[k][2]+scales[k][2]*(x-bounds[k][2]) : x
    warp=p->(coordinate(p[1],1),coordinate(p[2],2))
    lines=[(1.0,0.0,-bounds[1][1]),(1.0,0.0,-bounds[1][2]),
           (0.0,1.0,-bounds[2][1]),(0.0,1.0,-bounds[2][2])]
    fitted=_map_routes_on_line_arrangement(routes,warp,lines)
    fitted,Dict("map"=>"coordinatewise-positive-PL; identity near every puncture",
        "fixed_bounds"=>bounds,"exterior_scales"=>scales)
end

function _automatic_layout(H; seed=nothing)
    hash=canonical_hash(H); cut=cut_along_alpha(H;diagram_hash=hash)
    base=place_ports(cut,default_puncture_placement(cut))
    terminal=_terminal_ribbon_report(H,cut,base,build_wiring_problem(H,cut,base))
    terminal.planar || throw(ArgumentError("terminal ribbon map is not planar"))
    if terminal.components>1
        # Side-by-side components must share the prescribed outer region.
        # Do not silently replace a nested region decoration by another one.
        blocks=H.code.regions.boundary_cycle_blocks
        shared=[b for b in blocks if length(b)>1]
        length(shared)==1 && length(only(shared))==terminal.components ||
            throw(ArgumentError("nested disconnected terminal regions require an explicit placement; automatic nesting is not supported"))
        outer=Dict{Int,Int}()
        for component in terminal.component_vertices
            candidates=[f for (f,cycle) in enumerate(terminal.face_cycles)
                if terminal.origin_vertex[first(cycle)] in component &&
                    terminal.face_to_boundary_cycle[f] in only(shared)]
            length(candidates)==1 || throw(ArgumentError("shared outer region does not select one face per terminal component"))
            outer[minimum(component)]=only(candidates)
        end
        terminal=merge(terminal,(outer_faces=outer,))
    end
    punctures,ports,routes,trace=_realize_terminal_map(H,cut,base,terminal)
    # Remove arbitrary translation/rotation/scale of the Tutte outer triangle
    # before moving individual disks. This global similarity preserves the
    # entire embedding and avoids unnecessary braiding in the transport.
    target=default_puncture_placement(cut)
    z=ComplexF64[complex(p.center...) for p in punctures]
    w=ComplexF64[complex(only(q.center for q in target if q.id==p.id)...) for p in punctures]
    zm=sum(z)/length(z); wm=sum(w)/length(w)
    denominator=sum(abs2, z.-zm)
    a=denominator>1e-12 ? sum((w.-wm).*conj.(z.-zm))/denominator : 1.0+0im
    abs(a)>1e-8 || (a=1.0+0im)
    similarity=p->begin q=a*(complex(p...)-zm)+wm; (real(q),imag(q)) end
    punctures,ports,routes=_map_layout_geometry(punctures,ports,routes,similarity;
        radius_updates=Dict(p.id=>abs(a)*p.radius for p in punctures))
    trace["placement_alignment"]=Dict("real"=>real(a),"imag"=>imag(a),"scale"=>abs(a))
    trace["fixture_used"]=false
    for samples in (32,64,128,256)
        polar,fitted=_polar_terminal_routes(cut,punctures,ports,routes;samples=samples)
        backend=Dict{String,Any}("backend"=>"pure_julia","seed"=>seed,
            "post_router"=>"terminal-planar-realization","terminal_embedding"=>trace,
            "polar_collar_samples"=>samples,"fixture_used"=>false)
        candidate=HeegaardLayout(hash,:alpha_cut,Dict{String,Any}("chart"=>"sphere-at-infinity","outer_halo"=>true),
            punctures,polar,fitted,_route_crossings(fitted),backend,
            LayoutValidationReport(false,ValidationIssue[],Dict{Symbol,Bool}()))
        report=validate_layout(H,candidate;expected_hash=hash)
        if report.valid
            # Keep the ordinary readable 2×g atlas placement; transport the
            # whole plane, not selected curves or their endpoint identities.
            moved,polar,fitted,transport=_transport_to_placement(H,hash,punctures,polar,fitted,
                default_puncture_placement(cut))
            backend["transport"]=transport
            fitted,compact=_compact_atlas_exterior(moved,fitted)
            backend["exterior_compaction"]=compact
            report=_layout_candidate_report(H,hash,moved,polar,fitted)
            report.valid || throw(ArgumentError("exterior compaction failed geometric validation"))
            return HeegaardLayout(hash,candidate.view,candidate.chart,moved,polar,fitted,
                _route_crossings(fitted),backend,report)
        end
        samples==256 && throw(ArgumentError("terminal realization failed geometric validation: "*
            join([string(i.code) for i in report.issues],", ")))
    end
end
