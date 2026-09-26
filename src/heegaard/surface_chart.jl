"""Derived 2D atlas: a punctured sphere and one annulus per alpha curve.

The base contains the original beta geometry. Each annulus has an actual closed
alpha core and transverse beta segments. Boundary maps, not projection lines,
join the charts. Mathematical code remains outside this layout object.
"""
struct HeegaardSurfaceChart
    diagram_hash::String
    genus::Int
    base::HeegaardLayout
    annuli::Vector{NamedTuple}
    intersections::Vector{NamedTuple}
    alpha_cycles::Vector{NamedTuple}
    beta_cycles::Vector{NamedTuple}
    provenance::Dict{String,Any}
    joining_disk::Union{Nothing,NamedTuple}
end

_surface_radial(radius, t) = (radius*cospi(-2t), radius*sinpi(-2t))

"""Construct a full-curve atlas using an existing, valid alpha-cut plot.

The base is S² minus 2g open disks. Attach annulus i along B_i+ (outer)
and B_i- (inner). If t is the alpha cyclic parameter, annulus angle is
-2πt, base B_i+ angle is +2πt, and base B_i- angle is -2πt. The annulus
has the usual xy orientation and the base has the opposite xy orientation;
these attachment maps give a consistently oriented closed genus-g surface.
Alpha runs clockwise on the annulus core. Positive beta traverses inner to
outer, so det(alpha tangent, beta tangent) has the prescribed sign.

Scope: the existing software-certified diagrams with g curves per family,
valid alpha-cut geometry, and at least one visit on every beta curve.
Certification here does not independently prove meridian-system completeness.
Raw/decorated inputs, unmatched
polar ports, and zero-intersection beta components are rejected explicitly.
No global router completeness or new manifold certification is asserted.
"""
function surface_chart(H::HeegaardDiagram, layout::HeegaardLayout; joining_disk=nothing)
    layout.view == :alpha_cut || throw(ArgumentError("surface_chart requires an alpha-cut layout"))
    validate(H).valid || throw(ArgumentError("the certified diagram has been modified and is no longer valid"))
    hash = canonical_hash(H)
    code = H.code.code
    genus = surface_genus(H)
    length(code.alpha) == length(code.beta) == genus ||
        throw(ArgumentError("surface_chart requires genus-g alpha and beta curve counts"))
    empty_beta = [c.id for c in code.beta if isempty(c.visits)]
    isempty(empty_beta) || throw(ArgumentError(
        "surface_chart has no geometric representative for zero-intersection beta curves: $(join(empty_beta, ", "))"))
    layout.diagram_hash==hash ||
        throw(ArgumentError("surface_chart layout hash does not match the diagram"))
    validate_layout(H, layout; expected_hash=hash).valid ||
        throw(ArgumentError("surface_chart requires a valid layout of this diagram"))
    if !isempty(layout.routes)
        all(r->r.curve isa BezierChain,layout.routes) && _exit_angle_degrees(layout)<=1e-3 ||
            throw(ArgumentError("surface_chart requires smooth normal exits; first call plot(diagram, feasible)"))
    end

    cut = cut_along_alpha(H; diagram_hash=hash)
    expected_ports = Dict(p.id=>p for p in place_ports(cut,layout.punctures))
    for p in layout.ports
        _vnorm(_vsub(p.point,expected_ports[p.id].point)) <= 1e-8 ||
            throw(ArgumentError("port $(p.id) does not realize its polar gluing parameter"))
    end
    for p in layout.punctures
        isfinite(p.radius) && p.radius>0 && all(isfinite,p.center) ||
            throw(ArgumentError("annulus attachments require finite positive boundary circles"))
    end
    for i in eachindex(layout.punctures), j in i+1:length(layout.punctures)
        p,q=layout.punctures[i],layout.punctures[j]
        _vnorm(_vsub(p.center,q.center))>p.radius+q.radius ||
            throw(ArgumentError("annulus attachment boundaries overlap"))
    end
    boundaries=Dict(p.id=>p for p in layout.punctures)
    ports=Dict((p.intersection_id,boundaries[p.boundary_id].side)=>p for p in layout.ports)
    alpha_positions=Dict(x=>(c.id,k) for c in code.alpha for (k,x) in enumerate(c.visits))
    beta_positions=Dict(x=>(c.id,k) for c in code.beta for (k,x) in enumerate(c.visits))
    annuli=NamedTuple[]
    for (pair,c) in enumerate(code.alpha)
        push!(annuli,(alpha_id=c.id,pair=pair,inner_radius=0.43,core_radius=0.70,
                     outer_radius=1.0,plus_boundary=2pair-1,minus_boundary=2pair,
                     visits=copy(c.visits),orientation="dx wedge dy",
                     outer_map="base theta = -annulus theta",
                     inner_map="base theta = annulus theta"))
    end
    intersections=NamedTuple[]
    for x in sort(code.intersections;by=x->x.id)
        plus,minus=ports[(x.id,:plus)],ports[(x.id,:minus)]
        abs(plus.parameter+minus.parameter-1)<=1e-10 ||
            throw(ArgumentError("intersection $(x.id) has incompatible gluing parameters"))
        t=plus.parameter
        inner,outer=_surface_radial(0.43,t),_surface_radial(1.0,t)
        start,finish=x.sign==1 ? (inner,outer) : (outer,inner)
        point=_surface_radial(0.70,t)
        alpha_tangent=(point[2],-point[1])
        beta_tangent=_vsub(finish,start)
        determinant=alpha_tangent[1]*beta_tangent[2]-alpha_tangent[2]*beta_tangent[1]
        sign(determinant)==x.sign || error("annulus orientation disagrees with intersection $(x.id)")
        label=get(code.metadata,"display_source_labels",false) && x.source_label!==nothing ? x.source_label : string(x.id)
        push!(intersections,(id=x.id,label=label,alpha_id=x.alpha,beta_id=x.beta,sign=x.sign,
            parameter=t,point=point,beta_start=start,beta_end=finish,
            alpha_tangent=alpha_tangent,beta_tangent=beta_tangent,
            alpha_position=alpha_positions[x.id][2],beta_position=beta_positions[x.id][2],
            plus_port=plus.id,minus_port=minus.id))
    end
    alpha_cycles=NamedTuple[(id=c.id,label=c.label,visits=copy(c.visits)) for c in code.alpha]
    beta_cycles=NamedTuple[(id=c.id,label=c.label,visits=copy(c.visits)) for c in code.beta]
    chart=HeegaardSurfaceChart(hash,genus,deepcopy(layout),annuli,intersections,
                         alpha_cycles,beta_cycles,source_manifest(H),nothing)
    joining_disk === nothing && return chart
    overlay=_joining_disk_layout(H,joining_disk,chart)
    HeegaardSurfaceChart(hash,genus,chart.base,annuli,intersections,
                         alpha_cycles,beta_cycles,chart.provenance,overlay)
end

surface_chart(PH::PuncturedHeegaardDiagram,layout::HeegaardLayout)=
    surface_chart(PH.diagram,layout;joining_disk=PH.joining_disk)

"""Use the common full-curve plot pipeline and retain the joining disk."""
function plot(PH::PuncturedHeegaardDiagram,layout::HeegaardLayout;view=:surface)
    view == :surface || throw(ArgumentError("punctured plot requires view=:surface to retain D"))
    validate_joining_disk(PH.diagram,PH.joining_disk)
    surface_chart(PH,plot(PH.diagram,layout))
end
plot(PH::PuncturedHeegaardDiagram;view=:surface,router=:planar)=
    plot(PH,layout_heegaard(PH.diagram;backend=:pure_julia,router);view)

surface_chart(::Union{HeegaardCode,DecoratedHeegaardCode}, ::HeegaardLayout) =
    throw(ArgumentError("surface_chart requires certified HeegaardDiagram; certify decorated code explicitly"))

plot(chart::HeegaardSurfaceChart)=chart
