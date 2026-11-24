# cf: lib/topology.jl/docs/specs/heegaard-plot.md — fixed reference alpha, full K-relative certificate and quotient crossings.

function _hp_reference_family(hc,requested)
    requested in (:auto,:l,:m) || throw(ArgumentError("reference_alpha must be :auto, :l or :m"))
    alpha=Set(hc_curve_names(hc,:alpha));g=hc.lc.genus
    family=g>0 && alpha==Set("l$i" for i in 1:g) ? :l :
           g>0 && alpha==Set("m$i" for i in 1:g) ? :m : nothing
    reference_roles=any(k->haskey(hc.roles,k),_lc_k_names(g))
    if isnothing(family)
        (reference_roles || requested!=:auto) && throw(ArgumentError("reference drawing requires alpha exactly l1..lg or m1..mg, with no free alpha"))
        return nothing
    end
    all(get(hc.roles,k,:none) in (:none,:alpha) for k in _lc_k_names(g)) ||
        throw(ArgumentError("reference drawing does not support beta on K"))
    requested in (:auto,family) || throw(ArgumentError("reference_alpha is an assertion: input alpha is bound to $family, not $requested"))
    family
end

_hp_cyclic_equal(a,b)=length(a)==length(b) && (isempty(a) || any(circshift(a,k)==b for k in 0:length(a)-1))

function _hp_reference_positions(p,contacts)
    Dict((r.name,side)=>(side==:start ? first(r.points) : last(r.points))
        for r in p.runs for side in (:start,:end) if haskey(contacts,(r.name,side)))
end

function _hp_square_endpoint_audit(before,after;pairs=nothing)
    old=Dict(r.name=>(first(r.points),last(r.points)) for r in before.runs if !r.closed)
    new=Dict(r.name=>(first(r.points),last(r.points)) for r in after.runs if !r.closed)
    error=Set(keys(old))==Set(keys(new)) ? maximum((_hp_dist(old[k][j],new[k][j]) for k in keys(old) for j in 1:2);init=0.) : Inf
    # cf: docs/specs/lickorish-code.md — the square identifies opposite edges
    # by translation. Fixed endpoints preserve edge order but the source
    # endpoints must first be on an edge and have their same-component mate.
    endpoints=[(r.component,q) for r in before.runs if !r.closed for q in (first(r.points),last(r.points))]
    boundary=true;paired=true
    for (component,q) in endpoints
        x,y=q;W,H=before.width,before.height
        edge=abs(x)<1e-8 ? :left : abs(x-W)<1e-8 ? :right : abs(y)<1e-8 ? :bottom : abs(y-H)<1e-8 ? :top : :none
        boundary &= edge!=:none && -1e-8<=x<=W+1e-8 && -1e-8<=y<=H+1e-8
        mate=edge==:left ? (W,y) : edge==:right ? (0.,y) : edge==:bottom ? (x,H) : (x,0.)
        paired &= count(t->t[1]==component && _hp_dist(t[2],mate)<1e-8,endpoints)==1
    end
    (cyclic_order=error<1e-8,offcircle=boundary ? 0 : 1,mirror_error=error,minimum_gap_degrees=nothing,
        minimum_spacing=true,existing_small_gaps_not_worsened=true,boundary=boundary,paired=paired,
        ok=before.square && after.square && boundary && paired && error<1e-8)
end

function _hp_reference_source_audit(hc,ctx;geometry=nothing)
    p=ctx.placement;contacts=ctx.contacts;words=hc.lc.words
    beta=hc_curve_names(hc,:beta);order=true
    for component in beta
        actual=Int[]
        for r in p.runs
            r.component==component || continue
            if haskey(contacts,(r.name,:start))
                label=contacts[(r.name,:start)]
                push!(actual,only(x for x in words[component] if abs(x)==label))
            end
            append!(actual,r.word)
        end
        order &= _hp_cyclic_equal(actual,words[component])
    end
    ports=true;seam_signs=true;positions=_hp_reference_positions(p,contacts)
    for r in p.runs
        r.closed && continue
        for side in (:start,:end)
            key=(r.name,side);haskey(contacts,key) || continue
            label=contacts[key];pt=positions[key]
            v=side==:start ? _hp_sub(r.points[2],pt) : _hp_sub(pt,r.points[end-1])
            if p.square
                k=only(k for k in ("l1","m1") if label in abs.(words[k]))
                tangent=k=="l1" ? (1.,0.) : (0.,-1.)
            else
                hi=argmin([abs(_hp_dist(pt,c)-rad) for (_,(c,rad)) in p.holes])
                name,(center,_)=p.holes[hi];d=_hp_sub(pt,center)
                tangent=startswith(name,"top") ? (-d[2],d[1]) : (d[2],-d[1])
            end
            observed=_hp_cross(v,tangent)
            expected=sign(only(x for x in words[r.component] if abs(x)==label))
            seam_signs &= abs(observed)>1e-10*max(_hp_norm(v)*_hp_norm(tangent),1.) && sign(observed)==expected
        end
    end
    if p.square
        for k in ("l1","m1")
            w=words[k];at=findfirst(x->abs(x)==1,w);rot=vcat(w[at:end],w[1:at-1])
            for j in 2:length(rot)
                label=abs(rot[j]);t=_hp_square_port_fraction(j-1,length(rot))
                expected=k=="l1" ? [(p.width*t,0.),(p.width*t,p.height)] :
                    [(0.,p.height*(1-t)),(p.width,p.height*(1-t))]
                actual=[pt for (key,pt) in positions if contacts[key]==label]
                ports &= length(actual)==2 && all(minimum(_hp_dist(q,r) for q in actual)<1e-7 for r in expected)
            end
        end
    else
        bylabel=Dict{Int,Vector{HCPlotPoint}}()
        for (key,pt) in positions;push!(get!(bylabel,contacts[key],HCPlotPoint[]),pt);end
        expected_labels=Set(abs(x) for i in 1:hc.lc.genus for x in words["l$i"])
        ports &= Set(keys(bylabel))==expected_labels
        for points in values(bylabel)
            ports &= length(points)==2 && _hp_dist((points[1][1],p.height-points[1][2]),points[2])<1e-7
        end
        for (hole,(center,radius)) in p.holes
            rows=Tuple{Float64,Int}[]
            for (key,pt) in positions
                abs(_hp_dist(pt,center)-radius)<1e-6 || continue
                push!(rows,(mod(atan(pt[2]-center[2],pt[1]-center[1]),2pi),contacts[key]))
            end
            sort!(rows;by=first)
            index=parse(Int,replace(hole,r"^(top|bot)"=>""))
            expected=abs.(words["l$index"])
            startswith(hole,"bot") && reverse!(expected)
            ports &= _hp_cyclic_equal(last.(rows),expected)
        end
    end
    isnothing(geometry) && (geometry=_hp_audit(p))
    (ok=order && ports && seam_signs && geometry.ok,beta_word_order=order,l_port_order=ports,seam_signs=seam_signs,geometry=geometry)
end

function _hp_reference_angles(hc,ctx,p,family,audit)
    rows=NamedTuple[];beta=Set(hc_curve_names(hc,:beta))
    if !p.square && family==:m
        for (i,c) in enumerate(audit.crossings)
            a,b=p.runs[c.run_a],p.runs[c.run_b]
            (a.role==:alpha || b.role==:alpha) || continue
            push!(rows,(label=audit.labels[i],degrees=c.angle,stencil_degrees=audit.stencil_angles[i],kind=:interior))
        end
    else
        visits=Dict{Int,Vector{Tuple{Float64,Float64}}}()
        for r in p.runs
            r.component in beta && !r.closed || continue
            n=length(r.points)
            for (side,pt,v) in ((:start,r.points[1],_hp_sub(r.points[min(n,4)],r.points[1])),
                                (:end,r.points[end],_hp_sub(r.points[end],r.points[max(1,n-3)])))
                label=ctx.contacts[(r.name,side)]
                any(x->abs(x)==label,hc.lc.words[string(family)*"1"]) || p.square || family==:l || continue
                if p.square
                    k=only(k for k in ("l1","m1") if label in abs.(hc.lc.words[k]))
                    startswith(k,string(family)) || continue
                    tangent=family==:l ? (1.,0.) : (0.,1.)
                else
                    hi=argmin([abs(_hp_dist(pt,c)-rad) for (_,(c,rad)) in p.holes])
                    center,_=p.holes[hi].second;d=_hp_sub(pt,center);tangent=(-d[2],d[1])
                end
                angle=acosd(clamp(abs(_hp_dot(v,tangent))/max(_hp_norm(v)*_hp_norm(tangent),1e-12),0.,1.))
                neighbour=side==:start ? r.points[2] : r.points[end-1]
                u=_hp_sub(neighbour,pt)
                local_angle=acosd(clamp(abs(_hp_dot(u,tangent))/max(_hp_norm(u)*_hp_norm(tangent),1e-12),0.,1.))
                push!(get!(visits,label,Tuple{Float64,Float64}[]),(local_angle,angle))
            end
        end
        for label in sort!(collect(keys(visits)))
            v=visits[label]
            length(v)==2 || throw(ArgumentError("quotient crossing $label does not have two cut copies"))
            push!(rows,(label=label,degrees=minimum(first.(v)),stencil_degrees=minimum(last.(v)),kind=:seam))
        end
    end
    rows
end

function _hp_reference_svg(hc,p,family;show=:both,title="Heegaard diagram · reference alpha",contacts=nothing)
    # Hidden K guides remain in the certificate, not in the visible curve set.
    visible=_hp_rebuild(p,[r for r in p.runs if r.role==:beta || (family==:m && !p.square && r.role==:alpha)])
    svg=_hp_svg(visible;show,title)
    io=IOBuffer()
    if show!=:alpha && !isnothing(contacts)
        labels=Set(abs(x) for i in 1:hc.lc.genus for x in hc.lc.words[string(family)*string(i)])
        for r in p.runs
            r.role==:beta && !r.closed || continue
            for side in (:start,:end)
                get(contacts,(r.name,side),0) in labels || continue
                pt=side==:start ? first(r.points) : last(r.points)
                println(io,"<circle cx=\"$(pt[1])\" cy=\"$(pt[2])\" r=\"2.8\" fill=\"#fcfbf7\" stroke=\"#255dba\" stroke-width=\"1.5\"/>")
            end
        end
    end
    if show==:beta
        return replace(svg,"</svg>"=>String(take!(io))*"</svg>")
    elseif p.square
        for side in (0.,family==:m ? p.width : p.height)
            points=family==:m ? [(side,0.),(side,p.height)] : [(0.,side),(p.width,side)]
            data="M"*join(("$(q[1]),$(q[2])" for q in points)," L")
            println(io,"<path data-component=\"$(family)1\" d=\"$data\" stroke=\"#ba433f\" stroke-width=\"3\" fill=\"none\"/>")
        end
    elseif family==:l
        # One alpha curve l_i is cut into the two openings l_i⁺ and l_i⁻. Drawing both rims in the alpha
        # colour reads as two curves, so only the positive opening carries it.
        for (name,(c,r)) in p.holes
            startswith(name,"top") || continue
            index=replace(name,r"^(top|bot)"=>"")
            println(io,"<circle data-component=\"l$index\" cx=\"$(c[1])\" cy=\"$(c[2])\" r=\"$r\" fill=\"none\" stroke=\"#ba433f\" stroke-width=\"3\"/>")
        end
    end
    replace(svg,"</svg>"=>String(take!(io))*"</svg>")
end

function _hp_plot_reference(hc,family;show,DEBUG_MODE,rounds,local_steps,log,infinity_face=nothing)
    started=time_ns();ctx=_hp_place(hc.lc,hc.roles;reference=true,infinity_face)
    placed_ns=time_ns()
    initial=_hp_reference_source_audit(hc,ctx)
    initial.ok || throw(ArgumentError("reference placement disagrees with input beta/K certificate: beta=$(initial.beta_word_order), ports=$(initial.l_port_order), geometry=$(initial.geometry.ok)"))
    initial_ns=time_ns();ctx=_hp_reference_straighten(ctx)
    chart_ns=time_ns()
    chart=_hp_reference_source_audit(hc,ctx)
    chart.ok || throw(ArgumentError("straight reference chart disagrees with the source beta/K certificate"))
    chart_audit_ns=time_ns();port_check=ctx.placement.square ? _hp_square_endpoint_audit : _hp_port_audit
    result=_hp_relax(ctx.placement;rounds,local_steps,DEBUG_MODE,log,fixed=ctx.fixed,port_check,check_final_angle=false,
        placement_check=_hp_empty_guard_check(hc,ctx.placement))
    relaxed_ns=time_ns()
    p=result.placement
    finalctx=(placement=p,contacts=ctx.contacts,fixed=ctx.fixed)
    certificate=_hp_reference_source_audit(hc,finalctx;geometry=result.audit)
    angles=_hp_reference_angles(hc,ctx,p,family,result.audit)
    minimum_angle=minimum((r.stencil_degrees for r in angles);init=90.)
    angle_ok=minimum_angle >= _HP_MIN_CROSSING_ANGLE
    immutable=all(r.points==only(s for s in ctx.placement.runs if s.name==r.name).points for r in p.runs if r.name in ctx.fixed)
    beta_K_crossings=sum(length(intersect(abs.(hc.lc.words[k]),abs.(hc.lc.words[b]))) for k in _lc_k_names(hc.lc.genus) for b in hc_curve_names(hc,:beta))
    audit=merge(result.audit,(reference=(family=family,certificate=certificate.ok,initial_certificate=initial.ok,chart_certificate=chart.ok,
        beta_word_order=certificate.beta_word_order,l_port_order=certificate.l_port_order,seam_signs=certificate.seam_signs,fixed_geometry=immutable,
        beta_K_crossings=beta_K_crossings,alpha_beta_crossings=length(angles),angles=angles,minimum_angle=isempty(angles) ? nothing : minimum_angle),
        angle_ok=angle_ok,ok=result.audit.ok && certificate.ok && immutable && angle_ok))
    checks=merge(result.checks,(reference_placement_gauss=2,total_global_gauss=result.checks.global_gauss+2))
    pipeline=(placement=(placed_ns-started)/1e9,initial_audit=(initial_ns-placed_ns)/1e9,
        chart=(chart_ns-initial_ns)/1e9,chart_audit=(chart_audit_ns-chart_ns)/1e9,
        relaxation=(relaxed_ns-chart_audit_ns)/1e9,final_certificate=(time_ns()-relaxed_ns)/1e9)
    merge(result,(audit=audit,checks=checks,pipeline_timings=pipeline,svg=_hp_reference_svg(hc,p,family;show,contacts=ctx.contacts),success=result.steps>0 && audit.ok,
        reference_alpha=family,reference_contacts=ctx.contacts,initial_placement=ctx.placement,seconds=(time_ns()-started)/1e9,visual_acceptance=:pending))
end
