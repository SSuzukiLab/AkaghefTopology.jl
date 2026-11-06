@testset "Lescop's decorated Heegaard formula" begin
    function diagram(alpha, beta, intersections)
        code = H.heegaard_code(alpha, beta, intersections)
        H.certify_heegaard(H.decorate(code; cellular=true))
    end

    # Lescop, https://arxiv.org/abs/1402.2261, Figures 1 and 3,
    # Examples 1.2 and 1.3. The complementary region containing w1
    # has oriented boundary turn -1 after omitting the outer boundary.
    # The two oriented beta arcs turn 1/2 and 3/2 in Fig. 3, so this
    # region follows the c-to-d beta edge forward and the d-to-c edge
    # backward, with two clockwise alpha half-turns and four corners.
    D1 = diagram([[1,2]], [[1,2]], [H.Intersection(i,1,1,1) for i in 1:2])
    @test H.lescop_intersection_terms(D1, [1]) == (ell2=0, sl=0)
    @test H.lescop_theta(D1, [1], 2) == (ell2=0, sl=0, e=-1//4, theta=1//4)
    @test H.lescop_theta(D1, [1], 1) == (ell2=0, sl=0, e=1//4, theta=-1//4)
    @test_throws ArgumentError H.lescop_theta(D1, [1], 3)
    @test_throws ArgumentError H.lescop_intersection_terms(D1, [1,2])

    # Fig. 3 D2: crossings c,d,e,f correspond to IDs 1,2,3,4.
    # The beta turns are 1/2, 0, -1/2 on beta1 and 0 on beta2.
    # With these turns the cycle beginning with alpha1:c->d has total
    # turn -1, identifying the region at infinity (w2).
    D2 = diagram([[1,2], [3,4]], [[1,2,4], [3]],
                 [H.Intersection(1,1,1,1), H.Intersection(2,1,1,1),
                  H.Intersection(3,2,2,1), H.Intersection(4,2,1,1)])
    @test H.lescop_theta(D2, [1,3], 1) == (ell2=0, sl=0, e=1//4, theta=-1//4)
end

@testset "Lescop rotation on the Koda lens-family ribbon" begin
    for p in 2:12
        x, y = p + 1, p + 2
        signs = [fill(-1, p); 1; -1]
        code = H.heegaard_code([[x; collect(p:-1:1); y]], [[y; collect(1:p); x]],
                               [H.Intersection(i,1,1,s) for (i,s) in enumerate(signs)])
        diagram = H.certify_heegaard(H.decorate(code; cellular=true))
        # The joining bigon follows alpha y->x and beta x->y.
        graph = H.ribbon_graph(diagram)
        cycles = H.boundary_cycles(graph)
        dart = only(d for d in graph.darts if
            graph.edges[d.edge].family == :alpha && d.origin == y && d.target == x && isodd(d.id))
        region = only(r for (r, cycle) in enumerate(cycles) if dart.id in cycle)
        @test length(cycles[region]) == 2
        # Koda (2008), Fig. 25 and Theorem 5.14. Compare the full cyclic
        # word and vertex signs to the registered all-positive lens family.
        ograph = H.to_odata(diagram; joining_disk=H.joining_disk(diagram, region)).ograph
        expected = [-p:-1; 1:p]
        @test length(ograph.gauss) == 1
        @test any(circshift(only(ograph.gauss), k) == expected for k in eachindex(expected))
        @test ograph.orientation == ones(Int, p)
        @test H.lescop_theta(diagram, [x], region).theta == (p + 1)//(2p)
        @test H.lescop_theta(diagram, [x], region).e == 0
        @test H.lescop_theta(diagram, [y], region).theta == (1 - p)//(2p)
        @test H.lescop_theta(diagram, [y], region).e == 1//2
    end
end
