# cf: lib/topology.jl/docs/specs/heegaard-code.md §3 — HC point/disk supports.
# cf: .research-scratch/PJ260922/astra-lanes/hc-plot/answer-6a.md — oriented final intervals; adapted to HCPlotRun.component and the existing chart convention [reconstructed]
function _hp_mark_cut_atlas(lc,p)
    owner=Dict{Any,String}();forward=Dict{Any,Int}();kk=Dict{Any,Any}()
    if p.square
        for (side,c,d) in ((:L,"m1",-1),(:R,"m1",-1),(:B,"l1",1),(:T,"l1",1))
            owner[(:side,side)]=c;forward[(:side,side)]=d
        end
        # The four copies of the reference l1/m1 crossing are fixed by
        # l1 pointing right and m1 pointing down in the source square.
        for c in ("l1","m1"),way in ("out","in")
            sg=_lc_sign_char(only(x for x in lc.words[c] if abs(x)==1))
            q=c=="l1" ? (way=="out" ? (0.,0.) : (p.width,p.height)) :
                (way=="out" ? (0.,p.height) : (p.width,0.))
            z=(1,sg,way)
            kk[(z,"l1")]=(q,(:side,q[2]==0. ? :B : :T))
            kk[(z,"m1")]=(q,(:side,q[1]==0. ? :L : :R))
        end
    else
        for (name,_) in p.holes
            owner[(:hole,name)]="l"*replace(name,r"^(top|bot)"=>"")
            forward[(:hole,name)]=startswith(name,"top") ? 1 : -1
        end
    end
    (owner=owner,forward=forward,kk=kk)
end

function _hp_mark_regions(hc, p, contacts;geometry=nothing)
    lc=hc.lc; allwords=lc.words; words=_lc_nonempty_words(allwords)
    bad(s)=throw(ArgumentError("_hp_mark_regions: "*s))
    cut=_hp_mark_cut_atlas(hc.lc,p)
    # cut.owner[boundary], cut.forward[boundary], cut.kk[(corner,support)]
    # are created by the chart constructor; forward is relative to CCW (hole)
    # or +x/+y (square side), and kk contains the four genus-one K-K corners.
    sval(c,l)=begin a=[x for x in words[c] if abs(x)==l];
        length(a)==1 || bad("label $l is not unique on $c"); only(a) end
    runs=p.runs; drawn=Set(r.component for r in runs if haskey(words,r.component))
    seams=Set(setdiff(collect(keys(words)),collect(drawn)))
    if isnothing(geometry)
        xs,ov=_hp_crossings(runs); isempty(ov)||bad("overlap in final polylines")
        ok,labs=_hp_match_words(runs,xs); ok||bad("final crossings do not realize run.word")
    else
        geometry.structure && isempty(geometry.overlaps) || bad("invalid final geometry certificate")
        xs=geometry.crossings;labs=geometry.labels
    end
    events=[Any[] for _ in runs]
    for (i,x) in enumerate(xs), (ri,s) in ((x.run_a,x.station_a),(x.run_b,x.station_b))
        c=runs[ri].component; haskey(words,c)||bad("crossing on unsupported run $(runs[ri].name)")
        push!(events[ri],(station=s,value=sval(c,labs[i]),point=x.point))
    end
    function boundary(q)
        if p.square
            d=(abs(q[1]),abs(q[1]-p.width),abs(q[2]),abs(q[2]-p.height))
            k=argmin(d); d[k]<1e-5||bad("contact is off square boundary")
            return (:side,(:L,:R,:B,:T)[k])
        end
        z=[(abs(_hp_dist(q,c)-r),n) for (n,(c,r)) in p.holes]
        isempty(z)&&bad("no opening for contact"); d,n=z[argmin(first.(z))]
        d<1e-5||bad("contact is off every opening"); (:hole,n)
    end
    ports=Dict{Tuple{String,Int,Symbol},Any}()
    for (ri,r) in enumerate(runs)
        r.closed&&continue; m=length(r.points)-1
        rows=((:start,:out,0.0,first(r.points),_hp_sub(r.points[2],r.points[1])),
              (:end,:in,Float64(m),last(r.points),_hp_sub(r.points[end],r.points[end-1])))
        for (side,way,s,q,t) in rows
            haskey(contacts,(r.name,side))||bad("missing contact $(r.name)/$side")
            l=abs(contacts[(r.name,side)]); n=_hp_norm(t); n>1e-10||bad("zero endpoint tangent")
            push!(events[ri],(station=s,value=sval(r.component,l),point=q))
            key=(r.component,l,way); haskey(ports,key)&&bad("duplicate half-port $key")
            ports[key]=(point=q,boundary=boundary(q),tangent=_hp_mul(t,1/n))
        end
    end
    function piece(r,a,b)
        m=r.closed ? length(r.points) : length(r.points)-1
        r.closed&&b<=a+1e-10&&(b+=m)
        at(s)=begin
            !r.closed&&s>=m-1e-10&&return last(r.points)
            u=mod(s,m); i=floor(Int,u)+1
            _hp_lerp(r.points[i],r.points[mod1(i+1,length(r.points))],u-floor(u))
        end
        q=HCPlotPoint[at(a)]
        for k in floor(Int,a)+1:ceil(Int,b)-1
            _hp_dist(last(q),at(k))>1e-10&&push!(q,at(k))
        end
        _hp_dist(last(q),at(b))>1e-10&&push!(q,at(b)); q
    end
    intervals=Dict{Tuple{String,Int},Vector{HCPlotPoint}}()
    for (ri,r) in enumerate(runs)
        haskey(words,r.component)||continue; r.component in seams&&bad("cut curve has a run")
        e=sort(events[ri];by=x->x.station); length(e)>=(r.closed ? 1 : 2)||bad("too few events on $(r.name)")
        for j in 1:(r.closed ? length(e) : length(e)-1)
            k=mod1(j+1,length(e)); a,b=e[j],e[k]; w=words[r.component]
            w[mod1(findfirst(==(a.value),w)+1,length(w))]==b.value||
                bad("nonconsecutive events on $(r.name)")
            key=(r.component,a.value); haskey(intervals,key)&&bad("duplicate interval $key")
            intervals[key]=piece(r,a.station,b.station)
        end
    end
    for (c,w) in words, x in w
        c in seams||haskey(intervals,(c,x))||bad("missing interval ($c,$x)")
    end
    occ=_lc_rotation_occ(words)
    bverts=Dict{Any,Vector{HCPlotPoint}}()
    for u in values(ports); push!(get!(bverts,u.boundary,HCPlotPoint[]),u.point); end
    for (q,b) in values(cut.kk); push!(get!(bverts,b,HCPlotPoint[]),q); end
    function corner(z,support)
        l,sg,way=z; c=occ[l][sg][1]
        if c in seams
            sg,way=_lc_rot_prev(sg,way); c=occ[l][sg][1]
        end
        if !(c in seams)
            u=get(ports,(c,l,Symbol(way)),nothing)
            u===nothing&&bad("missing half-port for $z")
            cut.owner[u.boundary]==support||bad("wrong cut owner at $z")
            return u.point,u.boundary
        end
        get(cut.kk,(z,support),nothing)===nothing&&bad("unrecorded K-K corner $z")
        cut.kk[(z,support)]
    end
    function trace_cut(s,way,a,b,bid)
        d=(way=="out" ? 1 : -1)*cut.forward[bid]
        if bid[1]==:side
            t=bid[2] in (:L,:R) ? (0.0,1.0) : (1.0,0.0)
            total=d*_hp_dot(_hp_sub(b,a),t); total>=-1e-7||bad("wrong side direction")
            any(1e-7<d*_hp_dot(_hp_sub(q,a),t)<total-1e-7 for q in get(bverts,bid,HCPlotPoint[]))&&
                bad("square cut arc skips a vertex")
            return HCPlotPoint[a,b]
        end
        c,r=only(v for (n,v) in p.holes if n==bid[2])
        θ(q)=atan(q[2]-c[2],q[1]-c[1]); aa=θ(a); bb=θ(b)
        total=mod(d*(bb-aa),2pi); total<1e-10&&(total=2pi)
        any(1e-7<mod(d*(θ(q)-aa),2pi)<total-1e-7 for q in get(bverts,bid,HCPlotPoint[]))&&
            bad("circular cut arc skips a vertex")
        Δ=d*total; n=max(2,ceil(Int,total*r/10))
        q=[(c[1]+r*cos(aa+Δ*k/n),c[2]+r*sin(aa+Δ*k/n)) for k in 0:n]
        q[1]=a; q[end]=b; q
    end
    faces,_=_lc_faces(words); out=Dict{String,Vector{Any}}(); copies=Dict{Any,Vector{Any}}()
    for f in faces
        arcs=Any[]
        for z in f.corners
            l,sg,way=z; c,pos=occ[l][sg]; w=words[c]
            key=way=="out" ? (c,w[pos+1]) : (c,w[mod1(pos,length(w))]); dir=way=="out" ? 1 : -1
            if c in seams
                t=_lc_tau(words,occ,z); s2,w2=_lc_rot_next(t[2],t[3]); z2=(t[1],s2,w2)
                a,ba=corner(z,c); b,bb=corner(z2,c); ba==bb||bad("face changes cut copy")
                pts=trace_cut(c,way,a,b,ba); push!(get!(copies,key,Any[]),(cut=ba,points=pts))
                push!(arcs,(edge=key,orientation=dir,curve=c,cut=ba,points=pts))
            else
                pts=dir==1 ? intervals[key] : reverse(intervals[key])
                push!(arcs,(edge=key,orientation=dir,curve=c,cut=nothing,points=pts))
            end
        end
        all(_hp_dist(last(arcs[i].points),first(arcs[mod1(i+1,length(arcs))].points))<1e-5
            for i in eachindex(arcs))||bad("nonclosing boundary $(f.id)")
        out[f.id]=arcs
    end
    incidence=Dict{Any,Vector{Int}}()
    for arcs in values(out), a in arcs; push!(get!(incidence,a.edge,Int[]),a.orientation); end
    all(length(v)==2 && sum(v)==0 for v in values(incidence))||
        bad("final face darts do not pair two-by-two")
    sides=Dict{String,Any}(); splitfaces=Set{String}()
    pending=[c for (c,w) in allwords if isempty(w)]
    while !isempty(pending)
        progress=false
        for c in copy(pending)
            ref=get(lc.places,c,""); host=(ref=="outer" && isempty(words)) ? (face="outer",tests=Any[]) : haskey(out,ref) ? (face=ref,tests=Any[]) :
                get(sides,ref,nothing); host===nothing&&continue
            rr=[r for r in runs if r.component==c]
            length(rr)==1&&rr[1].closed||bad("empty curve $c has no unique closed run")
            q=vcat(rr[1].points,[first(rr[1].points)]); push!(splitfaces,host.face)
            sides[c*"+"]=(face=host.face,tests=vcat(host.tests,[(c,+1)]),boundary=q)
            sides[c*"-"]=(face=host.face,tests=vcat(host.tests,[(c,-1)]),boundary=reverse(q))
            deleteat!(pending,findfirst(==(c),pending)); progress=true
        end
        progress||bad("unsupported/ambiguous empty placements: $pending")
    end
    (faces=out,intervals=intervals,cutcopies=copies,empty_sides=sides,split_faces=splitfaces)
end

# cf: hc_plot_empty.jl::free_disk — final-boundary clearance search; no curve is displaced.
_hp_mark_area(q)=sum(_hp_cross(_hp_sub(q[i],q[1]),_hp_sub(q[mod1(i+1,length(q))],q[1])) for i in eachindex(q))/2
function _hp_mark_inside(q,p)
    inside=false
    for i in eachindex(q)
        a,b=q[i],q[mod1(i+1,length(q))]
        (a[2]>p[2])!=(b[2]>p[2]) && p[1]<(b[1]-a[1])*(p[2]-a[2])/(b[2]-a[2])+a[1] && (inside=!inside)
    end
    inside
end
function _hp_mark_canvas(p)
    pts=[q for r in p.runs for q in r.points]
    x0=min(0.,minimum(first,pts;init=0.))-55.;x1=max(p.width,maximum(first,pts;init=p.width))+55.
    y0=min(0.,minimum(last,pts;init=0.))-55.;y1=max(p.height,maximum(last,pts;init=p.height))+55.
    HCPlotPoint[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]
end
function _hp_mark_polygon(arcs)
    q=HCPlotPoint[]
    for arc in arcs,pt in arc.points
        (isempty(q) || _hp_dist(q[end],pt)>1e-9) && push!(q,pt)
    end
    length(q)>1 && _hp_dist(q[1],q[end])<1e-9 && pop!(q)
    q
end
function _hp_mark_anchor(p,inside,paths)
    canvas=_hp_mark_canvas(p);segs=_hp_segments(p.runs)
    boundary=[(q[i],q[mod1(i+1,length(q))]) for q in paths for i in eachindex(q)]
    function clearance(pt)
        _hp_mark_inside(canvas,pt) || return -Inf
        inside(pt) || return -Inf
        any(_hp_dist(pt,c)<=r+1e-5 for (_,(c,r)) in p.holes) && return -Inf
        p.square && !(1e-5<pt[1]<p.width-1e-5 && 1e-5<pt[2]<p.height-1e-5) && return -Inf
        minimum(vcat([_hp_nearest(pt,a,b)[1] for (_,_,a,b) in segs],
            [_hp_nearest(pt,a,b)[1] for (a,b) in boundary]);init=Inf)
    end
    best=-Inf;point=(0.,0.)
    function consider(pt)
        d=clearance(pt)
        if d>best;best=d;point=pt;end
    end
    # Interior candidates scale with the actual region, including subpixel
    # bigons. Membership, rather than convexity, decides which are retained.
    for q in paths
        consider((sum(first,q)/length(q),sum(last,q)/length(q)))
        for i in 2:length(q)-1
            consider(((q[1][1]+q[i][1]+q[i+1][1])/3,(q[1][2]+q[i][2]+q[i+1][2])/3))
        end
    end
    # Edge-side probes retain narrow nonconvex faces missed by the coarse grid.
    for (a,b) in boundary
        v=_hp_sub(b,a);len=_hp_norm(v);len<1e-9 && continue
        mid=_hp_lerp(a,b,.5)
        for sign in (-1.,1.),d in (.02,.2,1.,3.,8.)
            consider(_hp_add(mid,(-sign*v[2]*d/len,sign*v[1]*d/len)))
        end
    end
    lo,hi=canvas[1],canvas[3]
    for i in 1:32,j in 1:32;consider((lo[1]+(hi[1]-lo[1])*i/33,lo[2]+(hi[2]-lo[2])*j/33));end
    best>1e-7 || throw(ArgumentError("HC mark has no resolved interior in the final chart"))
    step=max(hi[1]-lo[1],hi[2]-lo[2])/32
    for _ in 1:14
        origin=point
        for dx in (-1,0,1),dy in (-1,0,1);consider(_hp_add(origin,(dx*step,dy*step)));end
        step/=2
    end
    point,best
end
# cf: heegaard-code.md §3 and heegaard-solid-view.md §7 — join fine faces
# across undisplayed guides, keeping the cut circles as chart boundaries.
function _hp_disk_chart(hc)
    isempty(hc.disks) && return nothing
    isempty(hc.points) || return nothing
    lc = hc.lc
    retained = Set(keys(hc.roles))
    union!(retained, ["l$i" for i in 1:lc.genus])
    lc.genus == 1 && push!(retained, "m1")
    words = _lc_nonempty_words(lc.words)
    any(isempty, values(lc.words)) && return nothing
    occ = _lc_rotation_occ(words)
    chartwords = Dict(c => filter(x -> all(v -> first(v) in retained,
        values(occ[abs(x)])), w) for (c,w) in words if c in retained)
    any(isempty, values(chartwords)) && return nothing
    faces, corners = _lc_faces(words)
    parent = Dict(f.id => f.id for f in faces)
    for (corner, id) in corners
        label, sign, _ = corner
        first(occ[label][sign]) in retained && continue
        _uf_union!(parent, id, corners[_lc_tau(words, occ, corner)])
    end
    chartfaces, _ = _lc_faces(chartwords)
    supports = Dict{String,Set{String}}()
    for f in chartfaces
        roots = unique(_uf_find!(parent,corners[z]) for z in f.corners)
        length(roots) == 1 || return nothing
        supports[f.id] = Set(id for id in keys(parent) if _uf_find!(parent,id)==only(roots))
    end
    disks = Dict{String,Vector{String}}()
    for (name,ids) in hc.disks
        wanted = Set(ids)
        selected = sort!([id for (id,support) in supports if issubset(support,wanted)])
        isempty(selected) && return nothing
        union((supports[id] for id in selected)...) == wanted || return nothing
        disks[name] = selected
    end
    reduced = LickorishCode(lc.genus,chartwords,lc.places,lc.disks,lc.points,lc.weights,lc.source_path)
    chart = HeegaardCode(reduced,hc.roles,disks,hc.points,hc.weights,hc.hc_path)
    (; chart, supports)
end

function _hp_marks(hc,p,contacts;geometry=nothing)
    data=_hp_mark_regions(hc,p,contacts;geometry);canvas=_hp_mark_canvas(p)
    polygons=Dict(k=>_hp_mark_polygon(v) for (k,v) in data.faces)
    isempty(polygons) && (polygons["outer"]=canvas)
    empty_runs=Dict(r.component=>r.points for r in p.runs if r.closed && isempty(hc.lc.words[r.component]))
    function faceinside(ref,pt)
        q=polygons[ref]
        (_hp_mark_area(q)<0 ? !_hp_mark_inside(q,pt) : _hp_mark_inside(q,pt))
    end
    function region(ref)
        if haskey(polygons,ref)
            ref in data.split_faces && throw(ArgumentError("HC mark on split face $ref needs an explicit empty-curve side"))
            q=polygons[ref];paths=_hp_mark_area(q)<0 ? [canvas,q] : [q]
            return (inside=pt->faceinside(ref,pt),paths=paths,clipped=_hp_mark_area(q)<0)
        end
        haskey(data.empty_sides,ref) || throw(ArgumentError("unsupported HC mark region $ref"))
        side=data.empty_sides[ref]
        tests=side.tests
        contains(pt)=faceinside(side.face,pt) && all(((_hp_mark_inside(empty_runs[c],pt))==(sg>0)) for (c,sg) in tests)
        # Every empty circle is wholly inside its host region. Intersecting
        # nested positive sides picks the innermost circle; negative sides
        # remove only circles contained in that base, not the whole canvas.
        positives=[c for (c,sg) in tests if sg>0]
        base=isempty(positives) ? polygons[side.face] : empty_runs[last(positives)]
        clipped=isempty(positives) && _hp_mark_area(base)<0
        paths=clipped ? [canvas,base] : [base]
        for (c,sg) in tests
            sg<0 || continue
            q=empty_runs[c]
            (isempty(positives) || _hp_mark_inside(base,first(q))) && push!(paths,q)
        end
        (inside=contains,paths=paths,clipped=clipped)
    end
    points=NamedTuple[];disks=NamedTuple[]
    for (name,ref) in sort!(collect(hc.points);by=first)
        if haskey(polygons,ref) || haskey(data.empty_sides,ref)
            r=region(ref);point,clearance=_hp_mark_anchor(p,r.inside,r.paths)
            push!(points,(name=name,support=ref,point=point,clearance=clearance,kind=:region,inside=true))
        else
            match=Base.match(r"^(.+)/([+-]?\d+)$",ref)
            paths=Vector{HCPlotPoint}[]
            if !isnothing(match)
                c=String(match.captures[1]);x=parse(Int,match.captures[2]);key=(c,x)
                haskey(data.intervals,key) && push!(paths,data.intervals[key])
                haskey(data.cutcopies,key) && append!(paths,[v.points for v in data.cutcopies[key]])
            elseif haskey(empty_runs,ref)
                push!(paths,vcat(empty_runs[ref],[first(empty_runs[ref])]))
            end
            isempty(paths) && throw(ArgumentError("HC point $name has no final curve interval for $ref"))
            path=first(paths);total=sum(_hp_dist(path[i],path[i+1]) for i in 1:length(path)-1)
            point=first(_hp_path_slice(path,total/2,total/2))
            # Cut-circle copies are represented by sampled circular paths.
            if !isnothing(match) && match.captures[1] in _lc_k_names(hc.lc.genus) && startswith(match.captures[1],"l") && !p.square
                for (_,(center,radius)) in p.holes
                    abs(_hp_dist(point,center)-radius)<1. && (point=_hp_add(center,_hp_mul(_hp_sub(point,center),radius/_hp_dist(point,center)));break)
                end
            end
            push!(points,(name=name,support=ref,point=point,clearance=0.,kind=:curve,inside=true))
        end
    end
    for (name,refs) in sort!(collect(hc.disks);by=first)
        regions=[region(ref) for ref in refs]
        paths=[q for r in regions for q in r.paths]
        # Surface edge identities cancel even when two cut copies are far apart.
        incidences=Dict{Any,Vector{Any}}()
        for ref in refs
            haskey(data.faces,ref) || continue
            for arc in data.faces[ref];push!(get!(incidences,arc.edge,Any[]),arc);end
        end
        boundary=Vector{HCPlotPoint}[]
        for arcs in values(incidences)
            if length(arcs)==2 && sum(a.orientation for a in arcs)==0;continue;end
            append!(boundary,[a.points for a in arcs])
        end
        for ref in refs
            haskey(data.empty_sides,ref) && push!(boundary,data.empty_sides[ref].boundary)
        end
        point,clearance=_hp_mark_anchor(p,pt->any(r.inside(pt) for r in regions),paths)
        push!(disks,(name=name,support=refs,paths=paths,boundary=boundary,labelpoint=point,
            clearance=clearance,clipped=any(r.clipped for r in regions)))
    end
    (points=points,disks=disks,face_polygons=polygons,
        audit=(point_supports=all(x.inside for x in points),face_boundaries=true,oriented_edge_pairing=true))
end

function _hp_marks_svg(svg,marks,p;style=:filled)
    style in (:filled,:outline) || throw(ArgumentError("mark_style must be :filled or :outline"))
    background=IOBuffer();foreground=IOBuffer()
    path(q;closed=true)="M"*join(("$(a[1]),$(a[2])" for a in q)," L")*(closed ? " Z" : "")
    for disk in marks.disks
        name=_hp_xml(disk.name);support=_hp_xml(join(disk.support," "))
        if style==:filled
            # Each fine face fills separately; shared surface edges have
            # already been cancelled from the visible outline.
            data=join((path(q) for q in disk.paths)," ")
            println(background,"<path data-hc-disk=\"$name\" data-support=\"$support\" d=\"$data\" fill=\"#eee2b8\" fill-rule=\"evenodd\" stroke=\"none\"/>")
        end
        data=join((path(q;closed=false) for q in disk.boundary)," ")
        dash=style==:outline ? " stroke-dasharray=\"5 4\"" : ""
        println(background,"<path data-hc-disk-boundary=\"$name\" d=\"$data\" fill=\"none\" stroke=\"#9a813c\" stroke-width=\"1.7\"$dash/>")
        a=disk.labelpoint
        println(foreground,"<text data-hc-disk-label=\"$name\" x=\"$(a[1]+8)\" y=\"$(a[2]+18)\" fill=\"#806925\" font-size=\"14\">$name</text>")
    end
    for point in marks.points
        a=point.point;name=_hp_xml(point.name);support=_hp_xml(point.support)
        fill=style==:filled ? "#252525" : "#fcfbf7"
        println(foreground,"<circle data-hc-point=\"$name\" data-support=\"$support\" cx=\"$(a[1])\" cy=\"$(a[2])\" r=\"3.5\" fill=\"$fill\" stroke=\"#252525\" stroke-width=\"1.5\"/>")
        println(foreground,"<text x=\"$(a[1]+9)\" y=\"$(a[2]-9)\" fill=\"#252525\" font-size=\"15\">$name</text>")
    end
    stop=findfirst("</text>",svg)
    isnothing(stop) && error("HC SVG title missing")
    at=last(stop)
    svg=svg[1:at]*String(take!(background))*svg[nextind(svg,at):end]
    replace(svg,"</svg>"=>String(take!(foreground))*"</svg>")
end

# cf: .research-scratch/PJ260922/astra-lanes/hc-plot/answer-6c.md — fixed-source signed area and complete containment matrix; adapted argument order [reconstructed]
# `before` is the fixed placement at entry to _hp_relax.
function _hp_empty_guard(before, after, lc)
    fail(s) = throw(ArgumentError("_hp_empty_guard: " * s))
    names = sort([c for (c,w) in lc.words if isempty(w)])
    isempty(names) && return true

    all(c -> haskey(lc.places,c), names) ||
        fail("empty component without place")

    function ring(p,c)
        rs = [r for r in p.runs if r.component == c]
        length(rs) == 1 || return nothing
        r = only(rs)
        r.closed && isempty(r.word) || return nothing

        P = copy(r.points)
        length(P) > 1 && first(P) == last(P) && pop!(P)
        length(P) >= 3 || return nothing
        all(q -> isfinite(q[1]) && isfinite(q[2]), P) ||
            return nothing

        push!(P,first(P))             # satisfy either closed-ring convention
        P
    end

    function signed_area(P)
        o = P[1]
        s = sum(_hp_cross(_hp_sub(P[i],o),_hp_sub(P[i+1],o))
                for i in 1:length(P)-1)
        s == 0 && return 0.0
        copysign(abs(_hp_mark_area(P)),s)
    end

    src = Dict(c => ring(before,c) for c in names)
    dst = Dict(c => ring(after,c) for c in names)
    min_ratio = sqrt(eps(Float64))

    for c in names
        src[c] === nothing &&
            fail("source $c is not one finite closed polygon")
        dst[c] === nothing && return false

        a0 = signed_area(src[c])
        a1 = signed_area(dst[c])
        isfinite(a0) && a0 != 0 ||
            fail("source $c has zero or nonfinite area")
        isfinite(a1) &&
            sign(a1) == sign(a0) &&
            abs(a1) >= min_ratio*abs(a0) || return false
    end

    for host in names, child in names
        host == child && continue
        H0,C0 = src[host],src[child]
        H,C   = dst[host],dst[child]

        wanted = _hp_mark_inside(H0,first(C0))
        all(q -> _hp_mark_inside(H0,q) == wanted,C0) ||
            fail("source $child crosses $host")
        all(q -> _hp_mark_inside(H,q) == wanted,C) ||
            return false
    end

    true
end

function _hp_empty_guard_check(hc,placement)
    any(isempty,values(hc.lc.words)) || return nothing
    candidate->_hp_empty_guard(placement,candidate,hc.lc)
end
