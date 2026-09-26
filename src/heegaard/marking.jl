"""Relative path classes in the fixed standard alpha-cut surface.

The reference plane has disks (i, +/-1), radius .28, labelled ports, and
infinity fixed. Seam 2i-1 goes upward from (i+.04,1); seam 2i goes downward
from (i+.04,-1). The parts inside their disks are immaterial. Cutting along
these disjoint rays leaves a simply connected domain. Signed, freely reduced
crossing words therefore specify relative path classes, not just homology.
The annulus attachments are the existing zero-twist radial attachments.
See docs/heegaard/marking.md for the precise scope (and non-claims).
"""
const CUT_MARKING_FRAME = "standard-alpha-cut-infinity/1"

struct MarkedRoute
    beta_id::Int
    segment_id::Int
    start_port::Int
    end_port::Int
    word::Vector{Int}
end

struct CutSurfaceMarking
    frame::String
    source_hash::String
    routes::Vector{MarkedRoute}
end

struct MarkedHeegaardDiagram{T<:Union{HeegaardDiagram,PuncturedHeegaardDiagram}}
    object::T
    marking::CutSurfaceMarking
end

_base_diagram(h::HeegaardDiagram) = h
_base_diagram(h::PuncturedHeegaardDiagram) = h.diagram
_base_diagram(h::MarkedHeegaardDiagram) = _base_diagram(h.object)

function _reduce_crossing_word(word)
    result=Int[]
    for letter in word
        letter==0 && throw(ArgumentError("a seam letter cannot be zero"))
        !isempty(result) && last(result)==-letter ? pop!(result) : push!(result,letter)
    end
    result
end

_marking_key(r::MarkedRoute) = (r.beta_id,r.segment_id,r.start_port,r.end_port)
_marking_key(r::Route) = (r.curve_id,r.segment_id,r.start_port,r.end_port)
_exact_point(p) = (Rational{BigInt}(p[1]),Rational{BigInt}(p[2]))
_seam_origin(id) = (BigInt((id+1)÷2)//1 + BigInt(1)//25, isodd(id) ? BigInt(1)//1 : -BigInt(1)//1)

function _marking_frame_check(h,layout)
    cut=cut_along_alpha(h;diagram_hash=layout.diagram_hash)
    expected=default_puncture_placement(cut)
    # Existing transport artifacts can differ by a few rounding ulps. The
    # bounded near-identity frame identification is not arbitrary placement.
    actual=Dict(p.id=>p for p in layout.punctures)
    length(actual)==length(expected) && all(q->haskey(actual,q.id) && begin
        p=actual[q.id]
        (p.pair,p.side,p.orientation)==(q.pair,q.side,q.orientation) &&
        _vnorm(_vsub(p.center,q.center))<=1e-10 && abs(p.radius-q.radius)<=1e-10
    end,expected) ||
        throw(ArgumentError("marked layout requires the fixed standard alpha-cut frame; moving holes requires an explicit marking transport"))
    ports=place_ports(cut,expected)
    actual_ports=Dict(p.id=>p for p in layout.ports)
    length(actual_ports)==length(ports) && all(p->haskey(actual_ports,p.id) &&
        _vnorm(_vsub(p.point,actual_ports[p.id].point))<=1e-10,ports) ||
        throw(ArgumentError("marked layout requires the fixed polar ports"))
    all(c->!isempty(c.visits),vcat(h.code.code.alpha,h.code.code.beta)) ||
        throw(ArgumentError("marking v1 does not encode zero-intersection components"))
    surface_genus(h)==length(h.code.code.alpha) ||
        throw(ArgumentError("marking v1 requires the standard complete alpha-cut presentation"))
end

# Every cubic and its replacement chord lie in this exact rational convex
# hull. A bounding box excluding each puncture point certifies their homotopy
# in the punctured plane. Subdivide rather than trusting sampled flatness.
function _marking_curve_points(curve::BezierChain,n)
    result=Tuple{Rational{BigInt},Rational{BigInt}}[]
    origins=[_seam_origin(i) for i in 1:n]
    midpoint(a,b)=((a[1]+b[1])/2,(a[2]+b[2])/2)
    function append_cubic(ps,depth)
        safe=all(q->any(k->q[k]<minimum(p[k] for p in ps) ||
                          q[k]>maximum(p[k] for p in ps),1:2),origins)
        if safe
            push!(result,last(ps)); return
        end
        depth<40 || throw(ArgumentError("cannot certify Bezier path away from marking puncture points"))
        a,b,c,d=ps; ab=midpoint(a,b); bc=midpoint(b,c); cd=midpoint(c,d)
        abc=midpoint(ab,bc); bcd=midpoint(bc,cd); m=midpoint(abc,bcd)
        append_cubic((a,ab,abc,m),depth+1); append_cubic((m,bcd,cd,d),depth+1)
    end
    push!(result,_exact_point(first(curve.segments).p0))
    for s in curve.segments
        append_cubic(_exact_point.((s.p0,s.p1,s.p2,s.p3)),0)
    end
    result
end
_marking_curve_points(curve::PolylineCurve,n) = _exact_point.(curve.points)

"""Compute a reduced noncommutative seam word (exact predicates on stored coordinates)."""
function _seam_word(points,n)
    points=_exact_point.(points)
    for p in (first(points),last(points)), id in 1:n
        q=_seam_origin(id); direction=isodd(id) ? 1 : -1
        p[1]==q[1] && direction*(p[2]-q[2])>=0 &&
            throw(ArgumentError("route endpoint lies on a reference seam"))
    end
    word=Int[]
    for k in 1:length(points)-1
        a,b=points[k],points[k+1]
        events=Tuple{Rational{BigInt},Int}[]
        for id in 1:n
            q=_seam_origin(id); direction=isodd(id) ? 1 : -1
            # Half-open convention handles vertices on a ray and tangencies.
            ((a[1]<=q[1]<b[1]) || (b[1]<=q[1]<a[1])) || continue
            t=(q[1]-a[1])/(b[1]-a[1]); y=a[2]+t*(b[2]-a[2])
            y==q[2] && throw(ArgumentError("path hits a marking puncture point"))
            direction*(y-q[2])>0 || continue
            # Positive means left-to-right in the reference xy plane.
            push!(events,(t,b[1]>a[1] ? id : -id))
        end
        sort!(events;by=first)
        append!(word,last.(events))
    end
    _reduce_crossing_word(word)
end

function _layout_marking_routes(h,layout)
    n=length(layout.punctures)
    canonical_ports=Dict(p.id=>p.point for p in place_ports(cut_along_alpha(h;diagram_hash=layout.diagram_hash)))
    sort([MarkedRoute(r.curve_id,r.segment_id,r.start_port,r.end_port,
        _seam_word(vcat([_exact_point(canonical_ports[r.start_port])],
            _marking_curve_points(r.curve,n),[_exact_point(canonical_ports[r.end_port])]),n))
        for r in layout.routes];by=_marking_key)
end

"""Explicitly select the displayed representative; never infer a historic marking."""
function mark_heegaard(object::Union{HeegaardDiagram,PuncturedHeegaardDiagram},layout::HeegaardLayout)
    h=_base_diagram(object)
    validate_layout(h,layout).valid || throw(ArgumentError("marking requires a valid geometric witness"))
    _marking_frame_check(h,layout)
    object isa PuncturedHeegaardDiagram && validate_joining_disk(h,object.joining_disk)
    MarkedHeegaardDiagram(deepcopy(object),CutSurfaceMarking(CUT_MARKING_FRAME,
        _labelled_hash(object),_layout_marking_routes(h,layout)))
end

function _check_marking_data(m::MarkedHeegaardDiagram)
    h=_base_diagram(m); marking=m.marking
    validate(h).valid || throw(ArgumentError("marked input has an invalid underlying diagram"))
    m.object isa PuncturedHeegaardDiagram && validate_joining_disk(h,m.object.joining_disk)
    surface_genus(h)==length(h.code.code.alpha) ||
        throw(ArgumentError("marking v1 requires the standard complete alpha-cut presentation"))
    marking.frame==CUT_MARKING_FRAME || throw(ArgumentError("unsupported surface marking frame"))
    marking.source_hash==_labelled_hash(m.object) || throw(ArgumentError("marking belongs to a different labelled diagram or joining disk"))
    cut=cut_along_alpha(h;diagram_hash="marking-incidence")
    ports=place_ports(cut)
    expected=Set((c.id,k) for c in h.code.code.beta for k in eachindex(c.visits))
    length(marking.routes)==length(expected) && Set((r.beta_id,r.segment_id) for r in marking.routes)==expected ||
        throw(ArgumentError("marking must specify every beta segment exactly once"))
    halfedges=_halfedge_port_map(cut,ports)
    for r in marking.routes
        c=only(c for c in h.code.code.beta if c.id==r.beta_id)
        r.start_port==halfedges[(c.visits[r.segment_id],:out)] &&
        r.end_port==halfedges[(c.visits[mod1(r.segment_id+1,length(c.visits))],:in)] ||
            throw(ArgumentError("marked route endpoint incidence mismatch"))
        all(x->1<=abs(x)<=length(cut.boundaries),r.word) || throw(ArgumentError("unknown seam letter"))
        r.word==_reduce_crossing_word(r.word) || throw(ArgumentError("marked words must be freely reduced"))
    end
    all(c->!isempty(c.visits),vcat(h.code.code.alpha,h.code.code.beta)) ||
        throw(ArgumentError("marking v1 does not encode zero-intersection components"))
    true
end

"""Validate a witness against the marking, not merely the unmarked ribbon hash."""
function validate_marking(m::MarkedHeegaardDiagram,layout::HeegaardLayout)
    _check_marking_data(m)
    h=_base_diagram(m)
    validate_layout(h,layout).valid || throw(ArgumentError("invalid marked geometric witness"))
    _marking_frame_check(h,layout)
    actual=_layout_marking_routes(h,layout)
    target=sort(m.marking.routes;by=_marking_key)
    for (a,b) in zip(actual,target)
        _marking_key(a)==_marking_key(b) && a.word==b.word ||
            throw(ArgumentError("marking mismatch on beta $(b.beta_id) segment $(b.segment_id): expected $(b.word), got $(a.word); refusing a change of relative isotopy class"))
    end
    length(actual)==length(target) || throw(ArgumentError("marked route count mismatch"))
    true
end

function _marking_json(m::CutSurfaceMarking)
    Dict("frame"=>m.frame,"source_hash"=>m.source_hash,"annulus_attachment"=>"radial-zero/1",
         "routes"=>[Dict("beta_id"=>r.beta_id,"segment_id"=>r.segment_id,
             "start_port"=>r.start_port,"end_port"=>r.end_port,"word"=>r.word)
             for r in sort(m.routes;by=_marking_key)])
end
_json_dict(m::MarkedHeegaardDiagram)=Dict("schema"=>"akaghef.heegaard.marked/1",
    "diagram"=>_json_dict(m.object),"marking"=>_marking_json(m.marking))

"""Exact labelled marked-data identity, not a canonical unmarked diagram hash."""
marked_hash(m::MarkedHeegaardDiagram) = (_check_marking_data(m); _labelled_hash(m))

function save_marked(path::AbstractString,m::MarkedHeegaardDiagram)
    _check_marking_data(m)
    open(path,"w") do io; println(io,_canonical_json(_json_dict(m))); end
    path
end

function _marked_from_json(value)
    _spec_field(value,:schema)=="akaghef.heegaard.marked/1" || throw(ArgumentError("unsupported marked schema"))
    object=_spec_field(value.diagram,:schema)=="akaghef.heegaard.punctured/1" ?
        _punctured_from_json(value.diagram) : certify_heegaard(_heegaard_from_json(value.diagram))
    spec=value.marking
    _spec_field(spec,:annulus_attachment)=="radial-zero/1" || throw(ArgumentError("unsupported annulus attachment marking"))
    routes=[MarkedRoute(Int(r.beta_id),Int(r.segment_id),Int(r.start_port),Int(r.end_port),Int.(r.word)) for r in spec.routes]
    m=MarkedHeegaardDiagram(object,CutSurfaceMarking(String(spec.frame),String(spec.source_hash),routes))
    _check_marking_data(m); m
end
load_marked(path::AbstractString) = _marked_from_json(JSON3.read(read(path,String)))

function plot(m::MarkedHeegaardDiagram,witness::HeegaardLayout;view=:alpha_cut)
    validate_marking(m,witness)
    result=plot(_base_diagram(m),witness)
    validate_marking(m,result)
    view==:surface ? surface_chart(m,result) :
        view==:alpha_cut ? result : throw(ArgumentError("unsupported plot view"))
end

function layout_heegaard(m::MarkedHeegaardDiagram; witness=nothing,kwargs...)
    witness===nothing && throw(ArgumentError("marked input requires a geometric witness; unmarked automatic routing cannot choose its isotopy class"))
    isempty(kwargs) || throw(ArgumentError("marked layout options require explicit marking transport"))
    validate_marking(m,witness); deepcopy(witness)
end
plot(m::MarkedHeegaardDiagram;view=:surface,witness=nothing) = plot(m,layout_heegaard(m;witness);view)
function surface_chart(m::MarkedHeegaardDiagram,layout::HeegaardLayout)
    validate_marking(m,layout)
    chart=surface_chart(m.object,layout)
    provenance=deepcopy(chart.provenance)
    provenance["marking"]=Dict("hash"=>marked_hash(m),"data"=>_marking_json(m.marking))
    HeegaardSurfaceChart(chart.diagram_hash,chart.genus,chart.base,chart.annuli,
        chart.intersections,chart.alpha_cycles,chart.beta_cycles,provenance,chart.joining_disk)
end

function load_layout(path::AbstractString,m::MarkedHeegaardDiagram)
    layout=load_layout(path,_base_diagram(m))
    validate_marking(m,layout); layout
end
