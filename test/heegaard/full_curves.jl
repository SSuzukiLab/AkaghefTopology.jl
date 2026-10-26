using Test, Topology, JSON3
const H=Topology.Heegaard
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const SOURCE=joinpath(ROOT,"artifacts","heegaard","app")

@testset "Full curves: geometry, oriented incidence, and baseline preservation" begin
    entries=JSON3.read(read(joinpath(SOURCE,"manifest.json"),String))
    for entry in entries
        d=H.certify_heegaard(H.load_heegaard(joinpath(SOURCE,"$(entry.id).hgd.json")))
        l=H.load_layout(joinpath(SOURCE,"$(entry.id).hgd.layout.json"),d)
        before=H.render(l); core=H.canonical_hash(d)
        chart=H.surface_chart(d,l)
        unchanged=H.render(l)==before
        baseline_equal=before==replace(read(joinpath(SOURCE,"$(entry.id).svg"),String),"\r\n"=>"\n")
        @test unchanged
        @test baseline_equal # Git's Windows checkout converts tracked SVG line endings.
        @test H.canonical_hash(d)==chart.diagram_hash==core
        @test H.plot(chart)===chart
        @test repr(MIME"image/svg+xml"(),chart)==H.render(chart)
        @test chart.genus==entry.genus
        @test length(chart.annuli)==length(d.code.code.alpha)
        @test length(chart.intersections)==length(d.code.code.intersections)
        @test chart.base !== l
        @test Set(x.id for x in chart.intersections)==Set(x.id for x in d.code.code.intersections)
        portmap=Dict(p.id=>p for p in l.ports)
        # Follow each complete beta through annulus -> base -> next annulus.
        for beta in d.code.code.beta
            points=Dict(x.id=>x for x in chart.intersections)
            for (k,id) in enumerate(beta.visits)
                p=points[id]; next=points[beta.visits[mod1(k+1,length(beta.visits))]]
                route=only(r for r in l.routes if r.curve_id==beta.id && r.segment_id==k)
                @test route.start_port==(p.sign==1 ? p.plus_port : p.minus_port)
                @test route.end_port==(next.sign==1 ? next.minus_port : next.plus_port)
                @test portmap[route.end_port].intersection_id==next.id
            end
        end
        # Compute sign independently from geometric tangent directions; check
        # one transverse crossing of the alpha circle per radial beta interval.
        for p in chart.intersections
            s,f=p.beta_start,p.beta_end
            @test min(hypot(s...),hypot(f...))<0.70<max(hypot(s...),hypot(f...))
            @test hypot(p.point...)≈0.70
            t=(p.point[2],-p.point[1])
            v=(f[1]-s[1],f[2]-s[2])
            @test sign(t[1]*v[2]-t[2]*v[1])==only(x.sign for x in d.code.code.intersections if x.id==p.id)
        end
        for a in d.code.code.alpha
            actual=sort(filter(x->x.alpha_id==a.id,chart.intersections);by=x->x.parameter)
            @test [x.id for x in actual]==a.visits
        end
        svg=H.render(chart)
        @test length(collect(eachmatch(r"data-intersection-id=",svg)))==entry.intersections
        @test occursin("data-curve-family=\"alpha\"",svg)
        @test occursin("data-curve-family=\"beta\"",svg)
        @test occursin("No projection crossings",svg)
        @test_throws ArgumentError H.surface_chart(d.code,l)
        @test_throws ArgumentError H.surface_chart(d.code.code,l)
    end
end

@testset "Explicit boundaries and actual plot entrance" begin
    d=H.certify_heegaard(H.load_heegaard(joinpath(SOURCE,"S3.hgd.json")))
    seed=H.load_layout(joinpath(ROOT,"artifacts","heegaard","review_alpha_cut","S3__alpha_cut.hgd.layout.json"),d)
    @test_throws ArgumentError H.surface_chart(d,seed)
    @test_throws ArgumentError H.plot(d,seed;view=:unknown)
    atlas=H.plot(d,seed;view=:surface)
    @test atlas isa H.HeegaardSurfaceChart
    @test atlas.base.backend_trace["plot"]["profile"]=="fixed-100"
    @test H._exit_angle_degrees(atlas.base)<=1e-3
    @test H.validate_layout(d,atlas.base).valid

    raw=H.heegaard_code([Int[]],[Int[]],H.Intersection[])
    @test_throws H.ValidationError H.certify_heegaard(H.decorate(raw,nothing;cellular=true))
    wrong=H.certify_heegaard(H.load_heegaard(joinpath(SOURCE,"L5_2.hgd.json")))
    @test_throws ArgumentError H.surface_chart(wrong,atlas.base)

    # Curve ids are not required to coincide with positional puncture pairs.
    custom=H.certify_heegaard(H.decorate(H.heegaard_code(
        [H.Curve(:alpha,9,[71])],[H.Curve(:beta,4,[71])],
        [H.Intersection(71,9,4,-1)]),[[1]]))
    custom_atlas=H.plot(custom;view=:surface)
    @test only(custom_atlas.annuli).alpha_id==9
    @test only(custom_atlas.intersections).beta_id==4
    @test only(custom_atlas.intersections).sign==-1
end
