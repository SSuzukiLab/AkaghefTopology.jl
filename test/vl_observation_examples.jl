using Test, Topology
isdefined(@__MODULE__, :vl_example) ||
    include(joinpath(@__DIR__, "../examples/VLExample.jl"))

@testset "First-observation VL examples" begin
    expected = ((:lens_2_1, 2, 0, (2,)),
                (:lens_3_1_mst2_example23, 2, 0, (3,)),
                (:lens_3_2, 3, 0, (3,)),
                (:poincare, 5, 0, ()),
                (:quaternion, 3, 0, (2,2)),
                (:s2xs1_2vertex, 2, 1, ()),
                (:s2xs1, 4, 1, ()))
    for (name, vertices, rank, torsion) in expected
        @test String(name) in vl_example_names()
        v = vl_example(name)
        g = certify_closed_ograph(v)
        @test g.witness.true_vertices == vertices
        h = ograph_h1(ograph_pi1_presentation(g))
        @test (h.free_rank, h.torsion_invariant_factors) == (rank, torsion)
        @test !v.is_weighted
    end
    v = vl_example(:s2xs1_2vertex)
    @test real_gauss_code(v) == ([[-1,-2,2,1]], [-1,1])
    @test disk_table(v).edges_s == [[2,3,3,4,1,1],[-1,2,-3,4],[-4,-2]]
    # Opposite vertex signs are the nonclosed BP example, not an alias.
    opposite = VirtualLink()
    set_data!(opposite; rgauss=[[-1,-2,2,1]], orientation=[1,-1])
    @test !first(is_closed(opposite))
    set_weight!(v, [[0,0,0,0]])
    @test !vl_example(:s2xs1_2vertex).is_weighted
    @test length(run_vl_plot_calls()) == 19 # historical MATLAB replay stays fixed
end
