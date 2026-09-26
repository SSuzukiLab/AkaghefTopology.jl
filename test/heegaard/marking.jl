module HeegaardMarkingTests
using Test, Topology, JSON3
const H=Topology.Heegaard

function test_layout(h,ports,punctures,routes)
    hash=H.canonical_hash(h)
    l=H.HeegaardLayout(hash,:alpha_cut,Dict{String,Any}(),punctures,ports,routes,
        H._route_crossings(routes),Dict{String,Any}(),H.LayoutValidationReport(false,H.ValidationIssue[],Dict{Symbol,Bool}()))
    report=H.validate_layout(h,l)
    @test report.valid
    H.HeegaardLayout(hash,l.view,l.chart,punctures,ports,routes,l.crossings,l.backend_trace,report)
end

@testset "Marked relative path identity and torus Dehn twist" begin
    h=H.certify_heegaard(H.decorate(H.heegaard_code([[1]],[[1]],[H.Intersection(1,1,1,1)]);cellular=true))
    cut=H.cut_along_alpha(h); punctures=H.default_puncture_placement(cut); ports=H.place_ports(cut)
    a,b=ports[1].point,ports[2].point
    points=[a,(a[1]-.65,a[2]),(b[1]-.65,b[2]),b]
    l=test_layout(h,ports,punctures,[H.Route(:beta,1,1,1,2,points)])
    m=H.mark_heegaard(h,l)
    @test H.validate_marking(m,l)
    @test only(m.marking.routes).word==Int[]
    @test_throws ArgumentError H.layout_heegaard(m)
    @test H.layout_heegaard(m;witness=l).routes[1].points==points
    # An orientation-preserving twist fixes both boundary circles and every
    # labelled port. The unmarked intersection code is unchanged.
    dense=NTuple{2,Float64}[]
    for k in 1:length(points)-1
        a,b=points[k],points[k+1]
        for t in range(0,1;length=151)[1:end-1]
            push!(dense,((1-t)*a[1]+t*b[1],(1-t)*a[2]+t*b[2]))
        end
    end
    push!(dense,last(points))
    twisted=map(dense) do p
        x,y=p[1]-1,p[2]-1; r=hypot(x,y)
        theta=2pi*clamp((.65-r)/(.65-.28),0,1)
        (1+cos(theta)*x-sin(theta)*y,1+sin(theta)*x+cos(theta)*y)
    end
    twisted[1]=first(points); twisted[end]=last(points)
    t=test_layout(h,ports,punctures,[H.Route(:beta,1,1,1,2,twisted)])
    mt=H.mark_heegaard(h,t)
    @test l.diagram_hash==t.diagram_hash
    @test length(only(mt.marking.routes).word)==1
    @test H.marked_hash(m)!=H.marked_hash(mt)
    @test_throws ArgumentError H.validate_marking(m,t)
    @test_throws ArgumentError H.plot(m,t)
    @test H.validate_marking(mt,t)
    # Same geometry with extra collinear points has the same marked identity.
    subdivided=test_layout(h,ports,punctures,[H.Route(:beta,1,1,1,2,dense)])
    @test H.marked_hash(H.mark_heegaard(h,subdivided))==H.marked_hash(m)
    @test H._reduce_crossing_word([1,2,-2,-1,2])==[2]
    @test H._reduce_crossing_word([1,2,-1,-2])==[1,2,-1,-2]
    # Excursions touching a seam are cancelled; crossing at a polyline vertex
    # is counted once. These checks use exact predicates, not an epsilon.
    seamx=BigInt(26)//25
    @test H._seam_word([(0//1,2//1),(seamx,2//1),(0//1,3//1)],2)==Int[]
    @test H._seam_word([(0//1,2//1),(seamx,2//1),(2//1,2//1)],2)==[1]
    mktempdir() do dir
        file=joinpath(dir,"marked.json"); H.save_marked(file,mt)
        @test H.marked_hash(H.load_marked(file))==H.marked_hash(mt)
        data=H._json_dict(mt)
        data["marking"]["routes"][1]["word"]=[99]
        @test_throws ArgumentError H._marked_from_json(JSON3.read(JSON3.write(data)))
    end
    # The common smoothing/render path must preserve the selected class.
    plotted=H.plot(m,l)
    @test H.validate_marking(m,plotted)
    @test H.surface_chart(m,plotted).genus==1
end

@testset "Existing surface witnesses retain D and all segment classes" begin
    root=normpath(joinpath(@__DIR__,"../../artifacts/heegaard/full_curves"))
    for (name,layoutname) in (("Abalone","AbaloneAuto"),("PoincareD","PoincareDAuto"))
        ph=H.load_punctured(joinpath(root,name=="Abalone" ? "abalone/S3.punctured.hgd.json" : name*".punctured.hgd.json"))
        layout=H.load_layout(joinpath(root,layoutname*".hgd.layout.json"),ph.diagram)
        marked=H.mark_heegaard(ph,layout)
        @test H.validate_marking(marked,layout)
        @test length(marked.marking.routes)==length(ph.diagram.code.code.intersections)
        @test marked.object.joining_disk.boundary_darts==ph.joining_disk.boundary_darts
        encoded=H._json_dict(marked)
        @test !haskey(encoded["marking"],"layout")
        loaded=H._marked_from_json(JSON3.read(JSON3.write(encoded)))
        @test H.validate_marking(loaded,layout)
        @test H.marked_hash(loaded)==H.marked_hash(marked)
        damaged=deepcopy(marked)
        push!(damaged.marking.routes[1].word,2length(ph.diagram.code.code.alpha))
        @test_throws ArgumentError H.validate_marking(damaged,layout)
    end
end

@testset "Genus three marking, disconnected terminal components" begin
    code=H.heegaard_code([[1],[2],[3]],[[1],[2],[3]],[H.Intersection(i,i,i,1) for i in 1:3])
    h=H.certify_heegaard(H.decorate(code,[[1,2,3]]))
    layout=H.layout_heegaard(h;backend=:pure_julia,router=:planar)
    m=H.mark_heegaard(h,layout)
    @test H.validate_marking(m,layout)
    @test length(m.marking.routes)==3
end
end
