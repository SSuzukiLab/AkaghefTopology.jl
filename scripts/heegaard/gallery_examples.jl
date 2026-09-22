# Catalog constructors, not a replacement for the production layout engine.
# Every example retains its code, explicit decoration, feasible seed and plot.
function lens(p,q)
    p>0 && 0<q<=p && gcd(p,q)==1 || throw(ArgumentError("coprime lens parameters required"))
    code=H.heegaard_code([collect(1:p)],[[mod1(1+(i-1)*q,p) for i in 1:p]],
        [H.Intersection(i,1,1,1) for i in 1:p])
    H.certify_heegaard(H.decorate(code,[[i] for i in 1:p]))
end

function lens_sum(factors)
    alpha=Vector{Int}[]; beta=Vector{Int}[]; points=H.Intersection[]
    offset=0
    for (g,(p,q)) in enumerate(factors)
        push!(alpha,collect(offset+1:offset+p))
        push!(beta,[offset+mod1(1+(i-1)*q,p) for i in 1:p])
        append!(points,[H.Intersection(offset+i,g,g,1) for i in 1:p])
        offset+=p
    end
    code=H.heegaard_code(alpha,beta,points)
    graph=H.ribbon_graph(code); cycles=H.boundary_cycles(graph)
    owner=Dict(p.id=>p.alpha for p in points)
    # Remove a disk from one complementary disk in each summand and join.
    joined=[findfirst(c->owner[graph.darts[first(c)].origin]==g,cycles) for g in eachindex(factors)]
    regions=vcat([Int.(joined)],[[i] for i in eachindex(cycles) if !(i in joined)])
    H.certify_heegaard(H.decorate(code,regions))
end

function gallery_catalog(;extended=false)
    result=NamedTuple[]
    push!(result,(id="T3",title="T³ · three-torus",category="Source fixture",factors=Tuple{Int,Int}[],
        construction="Ozsváth–Stipsicz–Szabó, Heegaard Floer homology, Example 2.1.13 / Figure 2.5; octahedron and reflection gluing"))
    push!(result,(id="Poincare",title="Poincaré homology sphere",category="Source fixture",factors=Tuple{Int,Int}[],construction="Klebanov Fig. 2.15; region decoration is a reconstruction candidate"))
    for p in 1:16, q in 1:max(1,p-1)
        gcd(p,q)==1 || continue
        push!(result,(id=p==1 ? "S3" : "L$(p)_$(q)",title=p==1 ? "S³" : "L($p,$q)",
            category=p==1 ? "Stabilization" : "Lens space",factors=[(p,q)],
            construction="Genus-one cyclic-word construction; q is the successor step modulo p"))
    end
    for g in 2:6
        push!(result,(id="S3_g$g",title="S³ · genus $g",category="Stabilization",factors=fill((1,1),g),
            construction="$(g-1) stabilizations of the genus-one S³ diagram"))
    end
    for (p,q) in ((2,1),(3,1),(4,1),(5,1),(5,2),(7,1),(7,2)), g in 2:3
        push!(result,(id="L$(p)_$(q)_g$g",title="L($p,$q) · genus $g",category="Stabilization",
            factors=vcat([(p,q)],fill((1,1),g-1)),construction="$(g-1) stabilizations of the genus-one L($p,$q) diagram"))
    end
    sums=[[(2,1),(2,1)],[(2,1),(3,1)],[(3,1),(3,1)],[(2,1),(5,2)],
          [(3,1),(5,2)],[(5,1),(5,2)],[(5,2),(7,2)],[(7,1),(7,2)],
          [(2,1),(2,1),(2,1)],[(2,1),(3,1),(5,2)],[(3,1),(3,1),(3,1)],
          [(2,1),(2,1),(2,1),(2,1)]]
    for factors in sums
        push!(result,(id="sum_"*join(["$(p)_$q" for (p,q) in factors],"__"),
            title=join(["L($p,$q)" for (p,q) in factors]," # "),category="Connected sum",factors=factors,
            construction="Connected sum along complementary disks; explicit merged region block"))
    end
    extended && return result
    selection=JSON3.read(read(joinpath(@__DIR__,"..","..","web","heegaard","gallery-selection.json"),String))
    ids=Set(String(e.id) for e in vcat(collect(selection.representatives),collect(selection.comparisons)))
    filter(e->e.id in ids,result)
end

# Map a two-disk exterior conformally to a round annulus. The p strands have
# ordered angular coordinates at every radius. A common middle circle places
# the point at infinity strictly between strands, not on a strand. This gives
# a constructive seed for this family without expensive general route search.
function annulus_seed(diagram; narrow=false)
    cut=H.cut_along_alpha(diagram); disks=H.default_puncture_placement(cut)
    ports=H.place_ports(cut,disks); wires=H.build_wiring_problem(diagram,cut,ports).wires
    portmap=Dict(p.id=>p for p in ports)
    a=sqrt(1.0-0.28^2)
    to_annulus(p,column)=begin z=complex(p[1]-column,p[2]); (z-im*a)/(z+im*a) end
    function to_plane(w,column)
        z=im*a*(1+w)/(1-w); x=real(z)
        # A monotone map fixing the disks puts independent summands in
        # disjoint strips. The production relaxation can later use free space.
        width=narrow ? 0.16 : 2.0
        abs(x)>.30 && (x=sign(x)*(.30+width*tanh((abs(x)-.30)/width)))
        (column+x,imag(z))
    end
    routes=H.Route[]
    for curve in diagram.code.code.beta
        localwires=filter(w->w.beta_id==curve.id,wires)
        p=length(localwires); column=Float64(curve.id)
        angles=[mod(angle(to_annulus(portmap[w.from_port].point,column)),2π) for w in localwires]
        order=sortperm(angles); starts=angles[order]
        ends=[mod(angle(to_annulus(portmap[localwires[i].to_port].point,column)),2π) for i in order]
        for i in 2:p
            while ends[i]<=ends[i-1]; ends[i]+=2π; end
        end
        maximum(ends)-minimum(ends)<2π || error("annulus matching is not cyclically ordered")
        ends .-= 2π*round(sum(ends-starts)/(2π*p))
        midpoint=[2π*(i-.5)/p for i in 1:p]
        smooth(t)=t*t*(3-2t)
        for (rank,index) in enumerate(order)
            wire=localwires[index]; start=portmap[wire.from_port].point; finish=portmap[wire.to_port].point
            lo=log(abs(to_annulus(start,column))); hi=log(abs(to_annulus(finish,column)))
            function point(t)
                t==0 && return start
                t==1 && return finish
                theta=t<=.5 ? (1-smooth(2t))*starts[rank]+smooth(2t)*midpoint[rank] :
                    (1-smooth(2t-1))*midpoint[rank]+smooth(2t-1)*ends[rank]
                to_plane(exp(complex((1-t)*lo+t*hi,theta)),column)
            end
            # Fine constructive seed; regular plot simplification and the full
            # validator check the actual piecewise-linear representation.
            path=[point(t) for t in range(0.,1.;length=max(161,16p+1))]
            push!(routes,H.Route(:beta,curve.id,index,wire.from_port,wire.to_port,path))
        end
    end
    hash=H.canonical_hash(diagram)
    provisional=H.HeegaardLayout(hash,:alpha_cut,Dict{String,Any}("kind"=>"paired-puncture"),disks,ports,routes,
        H._route_crossings(routes),Dict{String,Any}("seed"=>"analytic-annulus-family","narrow_summands"=>narrow),
        H.LayoutValidationReport(false,H.ValidationIssue[],Dict{Symbol,Bool}()))
    report=H.validate_layout(diagram,provisional)
    report.valid || error("annulus seed failed: "*join([i.message for i in report.issues],"; "))
    verified=H.HeegaardLayout(hash,:alpha_cut,provisional.chart,disks,ports,routes,provisional.crossings,provisional.backend_trace,report)
    compact,_=H._simplify_polylines(diagram,verified,[r.points for r in routes];tolerance=0.004)
    candidate=H._validated_polyline_candidate(diagram,verified,compact)
    candidate===nothing ? verified : candidate
end
