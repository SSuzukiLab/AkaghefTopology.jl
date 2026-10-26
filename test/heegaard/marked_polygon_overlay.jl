module MarkedPolygonOverlayTests
using Test,Topology
const H=Topology.Heegaard
const Q=Rational{BigInt}
p(x,y)=(Q(x),Q(y))
triangle(id,points)=(;id,points)

# Source: exact elementary triangle areas and the two opposite diagonal
# triangulations of a unit square. Expected cells/areas are independent of
# Polyhedra's conversion output, not a second call into the same clipping code.
@testset "Exact labelled triangle overlay and contact retention" begin
    large=triangle(["source",0],[p(0,0),p(2,0),p(0,2)])
    small=triangle(["target",0],[p(0,0),p(1,0),p(0,1)])
    snapshot=deepcopy((large,small))
    result=H._marked_triangle_overlay([large],[small];registry_points=[p(1//2,0),p(0,1//2)])
    @test (large,small)==snapshot
    @test length(result.cells)==1
    cell=only(result.cells)
    @test cell.source_id==["source",0] && cell.target_id==["target",0]
    @test cell.area==1//2
    @test cell.centroid==p(1//3,1//3)
    @test Set(cell.vertices)==Set([p(0,0),p(1//2,0),p(1,0),p(0,1),p(0,1//2)])
    @test only(result.coverage).area_matches
    @test p(2,0) in result.registry # retained even though outside the intersection
    @test isempty(result.contacts)

    sources=[triangle("s1",[p(0,0),p(1,0),p(1,1)]),triangle("s2",[p(0,0),p(1,1),p(0,1)])]
    targets=[triangle("t1",[p(0,0),p(1,0),p(0,1)]),triangle("t2",[p(1,0),p(1,1),p(0,1)])]
    result=H._marked_triangle_overlay(sources,targets;registry_points=[p(1//2,0)])
    @test length(result.cells)==4
    @test all(c->c.area==1//4,result.cells)
    @test all(c->p(1//2,1//2) in c.vertices,result.cells)
    @test sum(c.area for c in result.cells)==1
    @test all(c->c.target_area==c.covered_area==1//2 && c.area_matches,result.coverage)
    @test all(pnt->pnt isa NTuple{2,Q},result.registry)
    # Global endpoint insertion must agree on every occurrence of a shared edge.
    @test all(c->p(1//2,0) in c.vertices,filter(c->p(0,0) in c.vertices && p(1,0) in c.vertices,result.cells))

    below=triangle("below",[p(0,0),p(0,-1),p(1,0)])
    edge=H._marked_triangle_overlay([small],[below];registry_points=[p(1//3,0)])
    @test isempty(edge.cells)
    @test only(edge.contacts).dimension==1
    @test only(edge.contacts).vertices==[p(0,0),p(1//3,0),p(1,0)]
    @test !only(edge.coverage).area_matches && only(edge.coverage).covered_area==0
    touching=triangle("point",[p(1,0),p(2,0),p(1,1)])
    point=H._marked_triangle_overlay([small],[touching])
    @test isempty(point.cells) && only(point.contacts).dimension==0
    @test only(point.contacts).vertices==[p(1,0)]
    far=triangle("far",[p(3,0),p(4,0),p(3,1)])
    disjoint=H._marked_triangle_overlay([small],[far])
    @test isempty(disjoint.cells) && isempty(disjoint.contacts)

    # Non-dyadic intersection: x+y=3 meets y=2-x/4 at (4/3,5/3).
    a=triangle("a",[p(0,0),p(3,0),p(0,3)])
    b=triangle("b",[p(0,1),p(4,1),p(0,2)])
    fractional=only(H._marked_triangle_overlay([a],[b]).cells)
    @test p(4//3,5//3) in fractional.vertices
    @test fractional.area==4//3
    @test fractional.centroid==p(5//6,17//12)

    epsilon=Q(1,big(10)^40)
    tiny=triangle("tiny",[p(0,0),p(epsilon,0),p(0,epsilon)])
    @test only(H._marked_triangle_overlay([tiny],[tiny]).cells).area==epsilon^2/2
    empty=H._marked_triangle_overlay([],targets)
    @test isempty(empty.cells) && all(c->!c.area_matches,empty.coverage)
    # Area equality alone intentionally does not certify arbitrary input tiling.
    duplicate=H._marked_triangle_overlay([small,triangle("second",small.points)],[small])
    @test only(duplicate.coverage).covered_area==1 && !only(duplicate.coverage).area_matches
    @test_throws ArgumentError H._marked_triangle_overlay([small,small],[small])
    @test_throws ArgumentError H._marked_triangle_overlay([triangle("clockwise",reverse(small.points))],[small])
    @test_throws ArgumentError H._marked_triangle_overlay([triangle("flat",[p(0,0),p(1,0),p(2,0)])],[small])
    @test_throws ArgumentError H._marked_triangle_overlay([triangle("float",[(0.0,0.0),(1.0,0.0),(0.0,1.0)])],[small])
    @test_throws ArgumentError H._marked_triangle_overlay([small],[small];registry_points=[(0.5,0.0)])
end

# Source: exact coordinate separation and a unit right-triangle overlap of
# leg epsilon. A zero gap retains a point contact; a positive gap is empty.
@testset "Strict rational bounding-box rejection preserves contacts" begin
    unit=triangle("unit",[p(0,0),p(1,0),p(0,1)])
    epsilon=Q(1,big(10)^40)
    for shift in (p(1+epsilon,0),p(0,1+epsilon))
        far=triangle("far",[(x+shift[1],y+shift[2]) for (x,y) in unit.points])
        for (source,target) in ((unit,far),(far,unit))
            result=H._marked_triangle_overlay([source],[target])
            @test isempty(result.cells) && isempty(result.contacts)
            @test Set(result.registry)==Set([source.points;target.points])
            @test only(result.coverage).covered_area==0 && !only(result.coverage).area_matches
        end
    end
    touching=triangle("touching",[p(1,0),p(2,0),p(1,1)])
    contact=H._marked_triangle_overlay([unit],[touching])
    @test only(contact.contacts).dimension==0 && only(contact.contacts).vertices==[p(1,0)]
    overlap=triangle("overlap",[p(1-epsilon,0),p(2-epsilon,0),p(1-epsilon,1)])
    @test only(H._marked_triangle_overlay([unit],[overlap]).cells).area==epsilon^2/2
    # Strictly disjoint geometry must not bypass input validation.
    far_clockwise=triangle("bad",[p(3,0),p(3,1),p(4,0)])
    @test_throws ArgumentError H._marked_triangle_overlay([unit],[far_clockwise])
end

# Source: the closed-segment affine parametrization (1-t)*a+t*b, 0<=t<=1.
# Endpoints and interior collinear points survive; tiny extensions do not.
@testset "Exact closed segment bounds preserve reinsertion and order" begin
    epsilon=Q(1,big(10)^40)
    cases=((p(0,0),p(2,0),p(epsilon,0)),
        (p(0,2),p(0,-2),p(0,epsilon)),
        (p(-1,-1),p(1,1),p(0,0)),
        (p(-1,1),p(1,-1),p(0,0)))
    for (a,b,middle) in cases
        before=(a[1]-epsilon*(b[1]-a[1]),a[2]-epsilon*(b[2]-a[2]))
        after=(b[1]+epsilon*(b[1]-a[1]),b[2]+epsilon*(b[2]-a[2]))
        off=(middle[1]+epsilon,middle[2]+2epsilon)
        registry=sort!(unique([b,off,before,middle,after,a,p(-3,4),p(4,-3)]))
        snapshot=copy(registry)
        @test H._mpo_segment_points(a,b,registry)==[a,middle,b]
        @test H._mpo_segment_points(b,a,registry)==[b,middle,a]
        @test registry==snapshot
    end
    @test_throws ArgumentError H._mpo_segment_points(p(0,0),p(0,0),[p(0,0)])
end
end
