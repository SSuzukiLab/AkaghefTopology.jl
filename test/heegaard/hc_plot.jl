# cf: lib/topology.jl/docs/specs/heegaard-plot.md — identity, cyclic starts, roles and the shared production/debug loop
using Test, Topology

# cf: lib/topology.jl/docs/conventions/LickorishCodes.md — cut-chart orientation in the SVG and at every reference crossing
@testset "LC display convention" begin
    H = Topology.Heegaard
    for g in 1:3
        @testset "genus $g" begin
            words = H._lc_k_constants(g)
            lc = H.LickorishCode(g, words, Dict{String,String}(), Dict{String,Vector{String}}(),
                Dict{String,String}(), Dict{String,Int}(), nothing)
            ctx = H._hp_place(lc; reference=true)
            p = ctx.placement
            svg = H._hp_svg(p)
            labels = Dict(m[2] => parse(Float64,m[1]) for m in
                eachmatch(r"y=\"([^\"]+)\"[^>]*>[lℓ](\d+[⁺⁻])</text>", svg))
            for i in 1:g
                @test labels["$(i)⁺"] > labels["$(i)⁻"]
            end
            if g == 1
                @test p.square
                continue
            end
            for i in 1:g, (row, clockwise) in (("top",true), ("bot",false))
                center, radius = only(v for (name,v) in p.holes if name=="$row$i")
                @test (center[2]>p.height/2) == clockwise
                angles = Dict{Int,Float64}()
                for run in p.runs, side in (:start,:end)
                    points = run.points
                    q, neighbor = side==:start ? (points[1],points[2]) : (points[end],points[end-1])
                    v = H._hp_sub(q,center)
                    abs(H._hp_norm(v)-radius)<1e-7 || continue
                    k = ctx.contacts[(run.name,side)]
                    signed_k = only(x for x in words["l$i"] if abs(x)==k)
                    angles[k] = atan(v[2],v[1])
                    t_l = clockwise ? (-v[2],v[1]) : (v[2],-v[1])
                    t_other = side==:start ? H._hp_sub(neighbor,q) : H._hp_sub(q,neighbor)
                    t_b, t_a = signed_k<0 ? (t_l,t_other) : (t_other,t_l)
                    @test -H._hp_cross(t_b,t_a)>0
                end
                a = [angles[abs(x)] for x in words["l$i"]]
                direction = clockwise ? 1 : -1
                @test sum(mod(direction*(a[mod1(j+1,length(a))]-a[j]),2pi) for j in eachindex(a)) ≈ 2pi
            end
        end
    end
end

@testset "HC plot signed identities" begin
    H=Topology.Heegaard
    @test isdefined(H,:plot)
    @test :plot in names(H)
    @test H._hp_outline(.125,520.,520.)==(260.,0.)
    @test H._hp_outline(.375,520.,520.)==(520.,260.)
    @test H._hp_outline(.625,520.,520.)==(260.,520.)
    @test H._hp_outline(.875,520.,520.)==(0.,260.)
    run(name,role,points,word;closed=false)=H.HCPlotRun(name,name,role,closed,H.HCPlotPoint[points...],word)
    a=run("x",:alpha,[(-10.,0.),(10.,0.)],[1])
    b=run("y",:beta,[(0.,-10.),(0.,10.)],[-1])
    hits,overlaps=H._hp_crossings([a,b])
    @test length(hits)==1
    @test isempty(overlaps)
    @test H._hp_match_words([a,b],hits)[1]
    @test hits[1].angle≈90
    @test H._hp_stencil_angle([a,b],hits[1])≈90
    fold=run("fold",:alpha,[(0.,0.),(10.,0.),(0.,0.)],Int[])
    @test !isempty(H._hp_crossings([fold])[2])
    degenerate=H.HCPlotCrossing(1,2,.5,.5,1,90.,(0.,0.))
    @test H._hp_stencil_angle([fold,b],degenerate)==0.0
    bad=run("y",:beta,b.points,[-2])
    @test !H._hp_match_words([a,bad],hits)[1]
    # Same partner/sign sequences, but mismatched crossing identities.
    a2=run("x",:alpha,a.points,[1,2])
    b2=run("y",:beta,b.points,[-1,-2])
    fake=[H.HCPlotCrossing(1,2,.2,.8,1,60.,(0.,0.)),H.HCPlotCrossing(1,2,.8,.2,1,60.,(1.,0.))]
    @test !H._hp_match_words([a2,b2],fake)[1]
    bc=run("y",:beta,b.points,[-1,-2];closed=true)
    @test H._hp_match_words([a2,bc],fake)[1]
    p=H.HCPlotPlacement(520.,520.,Pair{String,Tuple{H.HCPlotPoint,Float64}}[],[a,b],Dict{String,Symbol}(),nothing,0.,false)
    @test occursin("#ba433f",H._hp_svg(p;show=:alpha))
    @test !occursin("#255dba",H._hp_svg(p;show=:alpha))
    @test occursin("#255dba",H._hp_svg(p;show=:beta))
    @test_throws ArgumentError H._hp_svg(p;show=:unknown)
    lc=H.LickorishCode(0,Dict("x"=>Int[],"y"=>Int[]),Dict("x"=>"outer","y"=>"x+"),Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
    hc=H.HeegaardCode(lc,Dict("x"=>:alpha,"y"=>:beta),Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
    normal_result=H.plot(H.resolve_hc(lc,hc.roles);rounds=1,local_steps=2)
    @test normal_result.success
    @test normal_result.mode==:production
    @test normal_result.checks.global_gauss==1
    alpha_result=H.plot(H.resolve_hc(lc,hc.roles);show=:alpha,rounds=1,local_steps=2)
    @test occursin("#ba433f",alpha_result.svg)
    @test !occursin("#255dba",alpha_result.svg)
    @test_throws ArgumentError H.plot(hc;DEBUG_MODE=true,rounds=0)
    @test_throws ArgumentError H.draw_heegaard(hc,"unused.svg";png=true)
    placed=H._hp_place(lc,hc.roles)
    @test length(placed.runs)==2
    @test all(r.closed for r in placed.runs)
    closed_result=H.plot(H.resolve_hc(lc,hc.roles);DEBUG_MODE=true,rounds=1,local_steps=2)
    @test closed_result.steps==2
    @test closed_result.audit.ok
    @test closed_result.checks.global_gauss==closed_result.steps+closed_result.remeshed+1
    @test [r.points for r in normal_result.placement.runs]==[r.points for r in closed_result.placement.runs]
    @test all(r.closed for r in closed_result.placement.runs)
    holes=["top1"=>((50.,75.),10.),"bot1"=>((50.,25.),10.)]
    ports=[begin
        pt=(50.0+10cos(t),75.0+10sin(t))
        run(string(i),:alpha,[pt,(pt[1],100-pt[2])],Int[])
    end for (i,t) in enumerate((0.,2.,4.))]
    pp=H.HCPlotPlacement(100.,100.,holes,ports,Dict{String,Symbol}(),nothing,0.,false)
    @test H._hp_port_audit(pp,pp).ok
    @test H._hp_port_audit(pp,pp).minimum_spacing
    swapped=deepcopy(pp)
    swapped.runs[1].points[1],swapped.runs[2].points[1]=swapped.runs[2].points[1],swapped.runs[1].points[1]
    @test !H._hp_port_audit(pp,swapped).cyclic_order
    @test H._hp_port_audit(pp,swapped).mirror_error>1
    drift=deepcopy(pp);drift.runs[1].points[1]=(61.,75.)
    @test H._hp_port_audit(pp,drift).offcircle==1
    motion=[fill((1.0,0.0),length(r.points)) for r in pp.runs]
    portmotion=[fill((0.0,0.0),length(r.points)) for r in pp.runs]
    obstacle=(hole_entries=[("1","top1",1)],frame_entries=Tuple{String,Int}[])
    @test H._hp_hold_obstacles!(motion,portmotion,pp,obstacle,H._hp_ports(pp))
    @test all(==((0.0,0.0)),motion[1])
    @test all(==((1.0,0.0)),motion[2])
    @test H._hp_displacement(pp,deepcopy(pp))==0.0
    @test H._hp_displacement(pp,drift)>0.0
    emptyplot=H.HCPlotPlacement(520.,520.,Pair{String,Tuple{H.HCPlotPoint,Float64}}[],H.HCPlotRun[],Dict{String,Symbol}(),nothing,0.,false)
    @test H._hp_relax(emptyplot;rounds=1,local_steps=1).steps==0
end

# cf: threads/PJ260922/data/261112-original-compare/DRAWING-CHANGE-PROPOSAL.md §1–2 — square boundary regression.
@testset "HC general square boundary" begin
    H=Topology.Heegaard
    lc=H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\na : -3 +4 -5 +6\nb : -4 +5 -6 -2\n")
    hc=H.resolve_hc(lc,Dict("a"=>:alpha,"b"=>:beta))
    before=H._hp_place(lc,hc.roles)
    @test H._hp_port_audit(before,before).ok
    moved=copy(before.runs);r=moved[1];points=copy(r.points);points[1]=H._hp_add(points[1],(1.,1.))
    moved[1]=H._hp_run(r,points)
    @test !H._hp_port_audit(before,H._hp_rebuild(before,moved)).ok
    @test !H._hp_port_audit(H._hp_rebuild(before,moved),H._hp_rebuild(before,moved)).ok
    result=H.plot(hc)
    @test result.success && result.audit.structure && result.audit.ports.paired
    @test result.audit.ports.offcircle==0
    @test occursin("data-chart=\"square\"",result.svg)
    @test occursin("ℓ1⁺",result.svg) && occursin("m1⁻",result.svg)
end

@testset "HC plot local proposal gate" begin
    H=Topology.Heegaard
    run(name,role,points)=H.HCPlotRun(name,name,role,false,H.HCPlotPoint[points...],Int[])
    a=run("x",:alpha,[(-10.,0.),(0.,0.),(10.,0.)])
    b=run("y",:beta,[(0.,-10.),(0.,10.)])
    p=H.HCPlotPlacement(520.,520.,Pair{String,Tuple{H.HCPlotPoint,Float64}}[],[a,b],Dict{String,Symbol}(),nothing,0.,false)
    ws=H.HCPlotGateWorkspace(p);g=H._hp_gate_geometry(p,ws)
    @test g.ok
    @test length(g.crossings)==1
    state=H.HCPlotGateState(p.runs,g.crossings)
    shifted=H._hp_rebuild(p,[a,run("y",:beta,[(1.,-10.),(1.,10.)])])
    @test H._hp_gate_trial(shifted,H._hp_gate_geometry(shifted,ws),state).ok
    flipped=H._hp_rebuild(p,[a,run("y",:beta,[(1.,10.),(1.,-10.)])])
    @test !H._hp_gate_trial(flipped,H._hp_gate_geometry(flipped,ws),state).ok
    lost=H._hp_rebuild(p,[a,run("y",:beta,[(20.,-10.),(20.,10.)])])
    @test !H._hp_gate_trial(lost,H._hp_gate_geometry(lost,ws),state).ok
    folded=H._hp_rebuild(p,[run("fold",:alpha,[(0.,0.),(10.,0.),(0.,0.)])])
    @test !H._hp_gate_geometry(folded,ws).ok
    @test H._hp_same_event_order([1,2,3],[3,1,2],true)
    @test !H._hp_same_event_order([1,2,3],[3,1,2],false)
    left=run("left",:beta,[(-1.,-10.),(-1.,10.)]);right=run("right",:beta,[(1.,-10.),(1.,10.)])
    q=H._hp_rebuild(p,[a,left,right]);qg=H._hp_gate_geometry(q,ws)
    @test qg.ok
    qs=H.HCPlotGateState(q.runs,qg.crossings)
    exchanged=H._hp_rebuild(q,[a,run("left",:beta,[(2.,-10.),(2.,10.)]),right])
    @test !H._hp_gate_trial(exchanged,H._hp_gate_geometry(exchanged,ws),qs).ok
    ring=H.HCPlotRun("ring","ring",:alpha,true,H.HCPlotPoint[(-10.,0.),(0.,-10.),(10.,0.),(0.,10.)],Int[])
    c1=run("c1",:beta,[(-12.,-2.),(-6.,-2.)])
    c2=run("c2",:beta,[(8.,-8.),(8.,0.)])
    c3=run("c3",:beta,[(0.,8.),(4.,8.)])
    circle=H._hp_rebuild(p,[ring,c1,c2,c3]);cg=H._hp_gate_geometry(circle,ws)
    @test cg.ok && length(cg.crossings)==3
    cs=H.HCPlotGateState(circle.runs,cg.crossings)
    rotated=H._hp_rebuild(circle,[H._hp_run(ring,circshift(ring.points,1)),c1,c2,c3])
    @test H._hp_gate_trial(rotated,H._hp_gate_geometry(rotated,ws),cs).ok
end

@testset "HC plot local sufficient certificates" begin
    H=Topology.Heegaard
    A=[(-2.,0.),(0.,0.),(2.,0.)]
    B=[(-.2,-2.),(-.2,2.)];C=[(.2,-2.),(.2,2.)]
    chart=H.HPCrossingChart((0.,0.),(1.,0.,0.,1.),.25,.25,.75,.75,1)
    @test H._hp_cert_chart_linear_motion(chart,A,B,A,C)
    @test H._hp_cert_chart_linear_motion(chart,A,B,A,[(0.,-2.),(0.,2.)])
    @test !H._hp_cert_chart_linear_motion(chart,A,B,[(-2.,0.),(3.,0.),(2.,0.)],C)
    @test H._hp_cert_swept_pair((-1.,0.),(1.,0.),(-1.,1.),(1.,1.),.3,.3,.3,.3)
    @test !H._hp_cert_swept_pair((-1.,0.),(1.,0.),(-1.,1.),(1.,1.),.51,.51,.51,.51)
    @test !H._hp_cert_swept_pair(B[1],B[2],C[1],C[2],.4,.4,.4,.4)
    @test H._hp_closed_segment_contact((0.,0.),(1.,1e-15),(0.,1e-15),(1.,0.))
    @test !H._hp_closed_segment_contact((0.,0.),(1.,0.),(0.,1e-15),(1.,1e-15))
    @test H._hp_orientation_sign((0.,0.),(1.,1.),(nextfloat(1.),1.))==-1
    for (angle,expected) in ((29.,false),(31.,true))
        s,c=sincosd(angle)
        ch=H.HPCrossingChart((0.,0.),(1.,-c/s,0.,1/s),0.,0.,.5,.5,1)
        @test H._hp_chart_angle_certificate(ch)==expected
    end
    a,b,c=(0.,0.),(1.,1.),(2.,0.)
    @test !H._hp_cert_graph_delete(a,b,c,[((1.,.2),(1.,.8))])
    @test H._hp_cert_graph_delete(a,b,c,Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}[])
    @test H._hp_cert_graph_delete(a,b,c,[((-1.,0.),a),(c,(3.,0.))])
    @test !H._hp_cert_graph_delete(a,b,c,[(a,(1.,.2))])
    @test !H._hp_cert_graph_delete(a,b,c,[];holes=[((1.,.4),.1)])
    @test H._hp_cert_graph_delete(a,b,c,[];holes=[((1.,-1.),.1)])
    @test H._hp_port_collar_certificate((1.,0.),(2.,0.),(0.,0.),1.,.01,.01)
    @test !H._hp_port_collar_certificate((1.,0.),(.5,0.),(0.,0.),1.,0.,0.)
    @test !H._hp_port_collar_certificate((1.,0.),(1.,1.),(0.,0.),1.,.01,0.)
    @test !H._hp_port_collar_certificate((1.001,0.),(2.,0.),(0.,0.),1.,.01,0.)
    @test !H._hp_port_collar_certificate((1.,0.),(2.,0.),(0.,0.),1.,.01,-1.)
    @test H._hp_orientation_sign((0//1,0//1),(1//1,1//1),(1//1+1//big(10)^30,1//1))==-1
end

@testset "HC plot whole-circle port spread" begin
    H=Topology.Heegaard
    @test H._hp_port_spread(2pi/3,2pi/3)==0.
    @test H._hp_port_spread(.8,3.)>0 # both gaps exceed the old 35-degree support
    angles=[0.,.8,3.8]
    energy(t)=sum((mod(t[mod1(i+1,length(t))]-t[i],2pi)-2pi/length(t))^2 for i in eachindex(t))
    direction=[H._hp_port_spread(mod(angles[i]-angles[mod1(i-1,3)],2pi),mod(angles[mod1(i+1,3)]-angles[i],2pi)) for i in 1:3]
    @test energy(angles .+ 0.01 .* direction)<energy(angles)
end

@testset "HC reference roles and closed station" begin
    H=Topology.Heegaard
    lc=H._lc_parse_text("genus: 1\nl1 : -1\nm1 : 1 2\nb : -2\n")
    hm=H.resolve_hc(lc,Dict("b"=>:beta);standard=true)
    hl=H.resolve_hc(lc,Dict("b"=>:beta,"l1"=>:alpha))
    @test H._hp_reference_family(hm,:auto)==:m
    @test H._hp_reference_family(hl,:l)==:l
    @test_throws ArgumentError H._hp_reference_family(hm,:l)
    @test_throws ArgumentError H._hp_reference_family(hl,:m)
    @test_throws ArgumentError H._hp_reference_family(hm,:unknown)
    @test H.is_standard(hm)
    @test !H.is_standard(hl)
    ctx=H._hp_place(lc,hm.roles;reference=true)
    @test H._hp_reference_source_audit(hm,ctx).ok
    square_svg=H._hp_reference_svg(hm,ctx.placement,:m;show=:beta,contacts=ctx.contacts)
    @test occursin("data-chart=\"square\"",square_svg)
    @test occursin("m1⁺",square_svg) && occursin("ℓ1⁻",square_svg)
    @test !occursin("#ba433f",square_svg)
    r=only(ctx.placement.runs);P=r.points
    tangent=H._hp_run(r,H.HCPlotPoint[P[1],(P[1][1],P[1][2]+5.),P[end]])
    @test !H._hp_reference_source_audit(hm,(placement=H._hp_rebuild(ctx.placement,[tangent]),contacts=ctx.contacts,fixed=ctx.fixed)).seam_signs
    bad=copy(ctx.contacts);bad[first(keys(bad))]=999
    @test_throws ArgumentError H._hp_reference_source_audit(hm,(placement=ctx.placement,contacts=bad,fixed=ctx.fixed))
    ring=H.HCPlotRun("ring","ring",:beta,true,H.HCPlotPoint[(0.,0.),(10.,0.),(10.,10.),(0.,10.)],[1])
    line=H.HCPlotRun("line","line",:alpha,false,H.HCPlotPoint[(-10.,0.),(10.,0.)],[-1])
    p=H.HCPlotPlacement(520.,520.,Pair{String,Tuple{H.HCPlotPoint,Float64}}[],[ring,line],Dict{String,Symbol}(),nothing,0.,false)
    ws=H.HCPlotForceWorkspace(p)
    crossing=H.HCPlotCrossing(1,2,4.,.5,1,90.,(0.,0.))
    H._hp_force_prepare!(ws,p,[crossing])
    @test only(ws.events[1,2])[1]==0.
    # The circular tangent point must remain fixed from both sides of its
    # extremal height; a disappearing anchor used to make the chart jump.
    d=sqrt(40.0^2-20.0^2)
    meridian=H.HCPlotRun("@K:m1:1","m1",:alpha,false,H.HCPlotPoint[(30.,450.0-d),(0.,300.),(0.,220.),(30.,70.0+d)],Int[])
    probe=H.HCPlotRun("probe","probe",:beta,false,H.HCPlotPoint[(50.,410.0-1e-6),(50.,410.0+1e-6)],Int[])
    holes=["top1"=>((50.,450.),40.),"bot1"=>((50.,70.),40.)]
    chartp=H.HCPlotPlacement(520.,520.,holes,[probe,meridian],Dict{String,Symbol}(),nothing,0.,false)
    chart=H._hp_reference_straighten((placement=chartp,contacts=Dict{Tuple{String,Symbol},Int}(),fixed=Set([meridian.name])))
    @test all(abs(q[1]-50.)<1e-10 for q in chart.placement.runs[1].points)
    @test all(q[1]==30. for q in chart.placement.runs[2].points)
end

# cf: heegaard-code.md §3; original-compare/group-15/case-01 — final LC support, not a visual acceptance claim.
@testset "HC marks on final regions and curves" begin
    H=Topology.Heegaard
    lc=H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\na : -3 +4 -5 +6\nb : -4 +5 -6 -2\n")
    hc=H.resolve_hc(lc,Dict("a"=>:beta,"b"=>:alpha);points=Dict("z"=>"F1"))
    result=H.plot(hc)
    @test result.success
    @test occursin("data-hc-point=\"z\"",result.svg)
    mark=only(result.marks.points);polygons=result.marks.face_polygons
    @test mark.support=="F1" && mark.clearance>0
    @test H._hp_mark_inside(polygons["F1"],mark.point)
    @test !H._hp_mark_inside(polygons["F5"],mark.point)
    @test !H._hp_mark_inside(polygons["F6"],mark.point)
    data=H._hp_mark_regions(hc,result.placement,result.mark_contacts;geometry=result.audit)
    @test length(data.faces)==6
    @test length(data.faces["F5"])==length(data.faces["F6"])==2
    incidence=Dict{Any,Vector{Int}}()
    for arcs in values(data.faces),arc in arcs;push!(get!(incidence,arc.edge,Int[]),arc.orientation);end
    @test all(length(v)==2 && sum(v)==0 for v in values(incidence))
    @test_throws ArgumentError H._hp_mark_regions(hc,result.placement,Dict();geometry=result.audit)
    outline=H.plot(hc;mark_style=:outline)
    @test [r.points for r in outline.placement.runs]==[r.points for r in result.placement.runs]
    @test occursin("fill=\"#fcfbf7\" stroke=\"#252525\"",outline.svg)
    @test_throws ArgumentError H.plot(hc;mark_style=:unknown)
    lc2=H._lc_parse_text("genus: 2\nl1 : -1 -3\nl2 : -2 +4\nm1 : +1 +5\nm2 : +2 +7\nn1 : +3 +6 -4 -8\nb1 : -5 -6\nb2 : -7 +8\n")
    for family in ("l","m")
        roles=merge(Dict("b1"=>:beta,"b2"=>:beta),Dict("$family$i"=>:alpha for i in 1:2))
        marked=H.resolve_hc(lc2,roles;points=Dict("z"=>"F1","q"=>"l1/-1"),disks=Dict("D"=>["F1","F2"]))
        r=H.plot(marked)
        @test r.success && r.marks.audit.face_boundaries
        @test length(r.marks.points)==2 && length(r.marks.disks)==1
        @test occursin("data-hc-disk=\"D\"",r.svg)
        @test occursin("fill-rule=\"evenodd\"",r.svg)
    end
    sphere=H._lc_parse_text("genus: 0\na :\nb :\nplace a : outer\nplace b : a+\n")
    marked=H.resolve_hc(sphere,Dict("a"=>:alpha,"b"=>:beta);points=Dict("p"=>"b"),disks=Dict("D"=>["b+"]))
    r=H.plot(marked)
    @test r.success && only(r.marks.points).kind==:curve
    ring=only(s for s in r.placement.runs if s.component=="b");point=only(r.marks.points).point
    @test minimum(H._hp_nearest(point,ring.points[i],ring.points[mod1(i+1,length(ring.points))])[1] for i in eachindex(ring.points))<1e-8
end

# cf: answer-6c.md — fixed-source signed area and the complete containment matrix, including siblings.
@testset "HC empty curve placement guard" begin
    H=Topology.Heegaard
    lc=H._lc_parse_text("genus: 0\na :\nb :\nplace a : outer\nplace b : a+\n")
    roles=Dict("a"=>:alpha,"b"=>:beta);p=H._hp_place(lc,roles)
    @test H._hp_empty_guard(p,p,lc)
    translated=H._hp_rebuild(p,[H._hp_run(r,[H._hp_add(q,(3.,-2.)) for q in r.points]) for r in p.runs])
    @test H._hp_empty_guard(p,translated,lc)
    escaped=H._hp_rebuild(p,[r.component=="b" ? H._hp_run(r,[H._hp_add(q,(600.,0.)) for q in r.points]) : r for r in p.runs])
    @test !H._hp_empty_guard(p,escaped,lc)
    reversed=H._hp_rebuild(p,[H._hp_run(r,reverse(r.points)) for r in p.runs])
    @test !H._hp_empty_guard(p,reversed,lc)
    collapsed=H._hp_rebuild(p,[H._hp_run(r,[H._hp_add((260.,260.),H._hp_mul(H._hp_sub(q,(260.,260.)),1e-9)) for q in r.points]) for r in p.runs])
    @test !H._hp_empty_guard(p,collapsed,lc)
    siblings=H._lc_parse_text("genus: 0\na :\nb :\nplace a : outer\nplace b : outer\n")
    s=H._hp_place(siblings,roles);a=only(r for r in s.runs if r.component=="a")
    nested=H._hp_rebuild(s,[r.component=="b" ? H._hp_run(r,[H._hp_lerp((260.,260.),q,.25) for q in a.points]) : r for r in s.runs])
    @test !H._hp_empty_guard(s,nested,siblings)
    result=H.plot(H.resolve_hc(lc,roles))
    @test result.success && result.steps>0 && result.audit.empty_places
    @test result.checks.placement_checks>0
end
# cf: lib/topology.jl/docs/specs/heegaard-plot.md 「二角形」 — the corners of a bigon stay apart
@testset "HC plot bigon corners" begin
    H=Topology.Heegaard
    # alpha1 and beta1 meet in three consecutive crossings: two bigons that share a corner.
    lc=H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\nalpha1 : -3 +4 -5 +6\nbeta1 : -4 +5 -6 -2\n")
    result=H.plot(H.resolve_hc(lc,Dict("alpha1"=>:alpha,"beta1"=>:beta)))
    @test result.success
    crossings=result.audit.crossings;pairs=H._hp_bigon_pairs(result.placement.runs,crossings)
    @test length(pairs)==2
    # 30 is the arc-length window of an ordinary crossing; closer corners read as one crossing.
    @test minimum(H._hp_dist(crossings[a].point,crossings[b].point) for (a,b) in pairs)>=30
    # A single transverse crossing has no bigon.
    single=H._lc_parse_text("genus: 1\nl1 : -1 -2\nm1 : +1 +3\nc1 : -4 -3\nc2 : +4 +2\n")
    plain=H.plot(H.resolve_hc(single,Dict("c1"=>:alpha,"c2"=>:beta)))
    @test plain.success
    @test isempty(H._hp_bigon_pairs(plain.placement.runs,plain.audit.crossings))
end
# cf: lib/topology.jl/docs/specs/heegaard-plot.md 「正方形の縁の端点と周期的な写し」「K の同じ辺に続けて交わる弧」「交点ラベル」
@testset "HC plot torus lines, K bigons, labels and single alpha rim" begin
    H=Topology.HeegaardCodes
    @test H._hp_square_port_fraction(0,3)==0. && H._hp_square_port_fraction(3,3)==1.
    @test H._hp_square_port_fraction(1,3)==.25 && H._hp_square_port_fraction(2,3)==.75
    # L(7,2): b crosses m1 seven times and l1 twice. Every run is a straight segment of one common direction.
    l72=H._lc_parse_text("genus: 1\nl1 : -1 -9 -10\nm1 : +1 +2 +3 +4 +5 +6 +7 +8\nb : -5 -7 +10 -2 -4 -6 -8 +9 -3\n")
    line=H.plot(H.resolve_hc(l72,Dict("b"=>:beta);standard=true);crossing_labels=true)
    @test line.success
    directions=HCPlotDirections=Tuple{Float64,Float64}[]
    for r in line.placement.runs
        r.role==:beta || continue
        v=H._hp_sub(last(r.points),first(r.points));push!(directions,H._hp_mul(v,1/H._hp_norm(v)))
    end
    @test length(directions)==9
    @test all(abs(H._hp_cross(d,directions[1]))<2e-3 for d in directions)
    # alpha = m1 is the left and right edge: each of the seven crossings is labelled at both copies.
    @test sort(first.(line.crossing_labels))==sort(vcat(2:8,2:8))
    @test occursin(">V8</text>",line.svg)
    @test !haskey(H.plot(H.resolve_hc(l72,Dict("b"=>:beta);standard=true)),:crossing_labels)
    # Every path returns the LC label at both ends of every open run.
    @test all(haskey(line.contacts,(r.name,side)) for r in line.placement.runs if !r.closed for side in (:start,:end))
    ordinary=H._lc_parse_text("genus: 1\nl1 : -1 +2\nm1 : +1 +3\nalpha1 : -3 +4 -5 +6\nbeta1 : -4 +5 -6 -2\n")
    plain=H.plot(H.resolve_hc(ordinary,Dict("alpha1"=>:alpha,"beta1"=>:beta)))
    @test plain.success
    @test all(haskey(plain.contacts,(r.name,side)) for r in plain.placement.runs if !r.closed for side in (:start,:end))
    @test Set(values(plain.contacts))⊆Set(abs(x) for k in ("l1","m1") for x in ordinary.words[k])
    # A curve that meets the K curve n1 three times in a row has arcs with both ends on one side of a face.
    kbigon=H._lc_parse_text("genus: 2\nl1 : -1 -3\nl2 : -2 +4\nm1 : +1 +5\nm2 : +2 +7\nn1 : +3 +6 -9 +10 -4 -8\nb1 : -5 -6 +9 -10\nb2 : -7 +8\n")
    @test H.plot(H.resolve_hc(kbigon,Dict("b1"=>:beta,"b2"=>:beta);standard=true)).success
    # alpha = l_i: one rim per alpha curve, on the upper opening.
    rims=H._lc_parse_text("genus: 2\nl1 : -1 -3\nl2 : -2 +4\nm1 : +1 +5\nm2 : +2 +7\nn1 : +3 +6 -4 -8\nb1 : -5 -6\nb2 : -7 +8\n")
    svg=H.plot(H.resolve_hc(rims,Dict("l1"=>:alpha,"l2"=>:alpha,"b1"=>:beta,"b2"=>:beta))).svg
    @test count("data-component=\"l",svg)==2
end

# cf: Koda Figure 6, revised genus-three transcription, D bounded by V9,V17,V18,V13,V12,V10.
@testset "fill a Heegaard region on the ordinary smoothed curves" begin
    H = Topology.Heegaard
    lc = H._lc_parse_text("""
genus: 3
l1 : -1 -4 +19 +20 -21 -22
l2 : -2 -6 +5 -23
l3 : -3 -24 +25 +7 +26
m1 : +27 -28 +29 -30 -31 +1 +32
m2 : +33 -34 +35 -36 +2 +37
m3 : -38 +39 -40 +41 +3
n1 : -5 -42 +43 -44 -45 +4 +46
n2 : -7 +47 -48 -49 +50 +6 +51
alpha1 : -27 -10 +9 -8 +31 +22 +11
alpha2 : -26 -47 -15 +49 +14 +36 -13 +12 -33
alpha3 : -19 +45 +16 -17 +18 +23 -46
beta1 : -29 -9 +17 -43 +21 -11 -32 -37 -51
beta2 : +42 -18 +13 -35 -39 +15 +48 -41 +24
beta3 : -20 +44 -16 +8 +30 -25 +40 -14 -50 +38 +34 -12 +10 +28
""")
    roles = Dict("$(family)$i" => family for family in (:alpha,:beta) for i in 1:3)
    plain = H.plot(H.resolve_hc(lc,roles); crossing_labels=true)
    marked = H.plot(H.resolve_hc(lc,roles; disks=Dict("D"=>["F24"])); crossing_labels=true)
    @test plain.success && marked.success
    @test [(r.name,r.points) for r in marked.placement.runs] == [(r.name,r.points) for r in plain.placement.runs]
    @test marked.crossing_labels == plain.crossing_labels
    disk = only(marked.marks.disks)
    @test disk.support == ["F24"]
    @test length(disk.boundary) == 6
    ends = [first(edge) for edge in disk.boundary]
    @test Set(label for (label,q) in marked.crossing_labels if any(p->H._hp_dist(p,q)<1e-6,ends)) == Set([9,17,18,13,12,10])
    @test H._hp_mark_inside(marked.marks.face_polygons["F24"],disk.labelpoint)
end
