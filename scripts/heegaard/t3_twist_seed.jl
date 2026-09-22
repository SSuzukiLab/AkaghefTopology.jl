# OSS Example 2.1.13 / Figure 2.5. Explicit octahedron embedding followed by
# an analytic ambient twist. Source-specific feasible seed, not a general router.
function octahedron_seed(diagram)
    cut=H.cut_along_alpha(diagram)
    cut.genus==3 || error("octahedron seed requires the T3 fixture")
    target=H.default_puncture_placement(cut)
    centers=[(1.,1.),(3.,-1.),(2.,1.),(2.,-1.),(3.,1.),(1.,-1.)]
    disks=[H.Puncture(p.id,p.pair,p.side,centers[p.id],.28,p.orientation) for p in target]
    ports=H.place_ports(cut,disks);pm=Dict(p.id=>p for p in ports)
    wires=H.build_wiring_problem(diagram,cut,ports).wires
    edges=Dict{Tuple{Int,Int},Vector{NTuple{2,Float64}}}()
    for a in 1:6,b in a+1:6
        cld(a,2)==cld(b,2) && continue
        edges[(a,b)]=[centers[a],centers[b]]
    end
    edges[(1,5)]=[centers[1],(1.,1.8),(3.,1.8),centers[5]]
    edges[(2,6)]=[centers[2],(2.6,-.5),(1.4,-.5),centers[6]]
    edges[(1,4)]=[centers[1],(.3,1.),(.3,-1.8),(1.7,-1.8),centers[4]]
    edges[(4,5)]=[centers[4],(2.3,-1.8),(3.7,-1.8),(3.7,1.),centers[5]]
    edge(a,b)=a<b ? edges[(a,b)] : reverse(edges[(b,a)])
    partner=Dict{Int,Int}()
    for w in wires;partner[w.from_port]=w.to_port;partner[w.to_port]=w.from_port;end
    collars=Dict{Int,Vector{NTuple{2,Float64}}}()
    for disk in disks
        ps=sort(filter(p->p.boundary_id==disk.id,ports);by=p->p.parameter)
        angles=[begin
            v=H._vsub(edge(disk.id,pm[partner[p.id]].boundary_id)[2],disk.center)
            atan(v[2],v[1])
        end for p in ps]
        for i in 2:length(angles)
            while angles[i]<=angles[i-1];angles[i]+=2pi;end
        end
        angles[end]-angles[1]<2pi || error("octahedron rotation differs from source")
        starts=[2pi*p.parameter for p in ps]
        angles .-= 2pi*round(sum(angles-starts)/(2pi*length(ps)))
        for (i,p) in enumerate(ps)
            collars[p.id]=[begin
                r=.28+.14t;theta=(1-t)*starts[i]+t*angles[i]
                H._vadd(disk.center,(r*cos(theta),r*sin(theta)))
            end for t in range(0.,1.;length=81)]
            collars[p.id][1]=p.point
        end
    end
    routes=H.Route[]
    for c in diagram.code.code.beta, (i,w) in enumerate(filter(w->w.beta_id==c.id,wires))
        middle=edge(pm[w.from_port].boundary_id,pm[w.to_port].boundary_id)[2:end-1]
        push!(routes,H.Route(:beta,c.id,i,w.from_port,w.to_port,
            vcat(collars[w.from_port],middle,reverse(collars[w.to_port]))))
    end
    trace=Dict{String,Any}("seed"=>"source-octahedron-ambient-twist",
        "source"=>"OSS Example 2.1.13 / Figure 2.5",
        "scope"=>"constructive seed for the transcribed T3 octahedron only")
    function checked(disks,ports,routes)
        layout=H.HeegaardLayout(cut.diagram_hash,:alpha_cut,
            Dict{String,Any}("kind"=>"paired-puncture"),disks,ports,routes,H._route_crossings(routes),trace,
            H.LayoutValidationReport(false,H.ValidationIssue[],Dict{Symbol,Bool}()))
        report=H.validate_layout(diagram,layout;expected_hash=cut.diagram_hash)
        report.valid || error("octahedron seed: "*join([string(i.code)*": "*i.message for i in report.issues],"; "))
        H.HeegaardLayout(layout.diagram_hash,layout.view,layout.chart,disks,ports,routes,
            layout.crossings,trace,report)
    end
    checked(disks,ports,routes)
    @info "Octahedron equatorial seed valid"
    smooth(t)=begin u=clamp(t,0.,1.);u*u*(3-2u) end
    function rotate(p,center,theta)
        x,y=H._vsub(p,center)
        H._vadd(center,(cos(theta)*x-sin(theta)*y,sin(theta)*x+cos(theta)*y))
    end
    # (r,theta) -> (r,theta+f(r)) has inverse (r,theta-f(r)). The plateau
    # rigidly swaps the lower outside disks, fixing the middle and upper disks.
    # Two inverse local disk spins restore the prescribed boundary port angles.
    function transform(p)
        r=H._vnorm(H._vsub(p,(2.,-1.)))
        angle=pi*smooth((r-.4)/.25)*(1-smooth((r-1.35)/.3))
        q=rotate(p,(2.,-1.),angle)
        for center in ((1.,-1.),(3.,-1.))
            r=H._vnorm(H._vsub(q,center))
            q=rotate(q,center,-pi*(1-smooth((r-.28)/.16)))
        end
        q
    end
    # Resolve the analytic map before simplifying. Validate the actual exported
    # PL representation independently of the analytic homeomorphism argument.
    function mapped_segment(a,b,depth=0)
        mid=H._lerp(a,b,.5);fa,fb,fm=transform(a),transform(b),transform(mid)
        if H._vnorm(H._vsub(a,b))<=.01 && H._segment_distance_squared(fm,fa,fb)<1e-10
            return [fa,fb]
        end
        depth<24 || error("octahedron twist sampling limit")
        vcat(mapped_segment(a,mid,depth+1)[1:end-1],mapped_segment(mid,b,depth+1))
    end
    ports=H.place_ports(cut,target);pm=Dict(p.id=>p for p in ports)
    mapped=[begin
        points=NTuple{2,Float64}[]
        for i in 1:length(r.points)-1
            append!(points,mapped_segment(r.points[i],r.points[i+1])[1:end-1])
        end
        push!(points,transform(last(r.points)))
        points[1]=pm[r.start_port].point;points[end]=pm[r.end_port].point
        H.Route(r.family,r.curve_id,r.segment_id,r.start_port,r.end_port,points)
    end for r in routes]
    current=checked(target,ports,mapped)
    compact,_=H._simplify_polylines(diagram,current,[r.points for r in current.routes];tolerance=.001)
    current=H._validated_polyline_candidate(diagram,current,compact)
    trace["ambient_map"]="radial annular half-twist, followed by two inverse disk spins"
    trace["pipeline"]=Dict("puncture_placement"=>"fixed-2-by-3","feasible_router"=>"source-octahedron","smoothing"=>"not-run")
    current
end
