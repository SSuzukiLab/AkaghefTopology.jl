# Inline tests moved out of test/heegaard/runtests.jl (2026-11-17). Content unchanged.

@testset "Heegaard schema and certification" begin
    code = H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)])
    @test !H.validate(code).valid
    @test any(issue -> issue.code == :missing_decoration, H.validate(code).issues)
    H.validate(code; level=:syntax)
    H.validate(code; level=:incidence)
    diagram = s3_fixture()
    @test H.is_valid(diagram)
    @test H.is_closed_heegaard(diagram)
    @test H.surface_genus(diagram) == 1
    @test H.euler_characteristic(diagram) == 0
    @test H.is_cellular(diagram)
    @test_throws ArgumentError H.certify_heegaard(code)
    malformed = H.heegaard_code([[2]], [[2]], [H.Intersection(1, 1, 1, 1)])
    @test !H.validate(malformed).valid

    temporary = mktempdir()
    path = joinpath(temporary, "s3.hgd.json")
    H.save_heegaard(path, diagram)
    loaded = H.load_heegaard(path)
    @test loaded isa H.DecoratedHeegaardCode
    @test H.canonical_hash(loaded) == H.canonical_hash(diagram)
    @test H.validate(loaded).valid
end

@testset "Ribbon, surface, and exact algebra" begin
    diagram = s3_fixture()
    graph = H.ribbon_graph(diagram)
    cycles = H.boundary_cycles(graph)
    @test length(graph.darts) == 2length(graph.edges)
    @test sort(vcat(cycles...)) == collect(1:length(graph.darts))
    @test H.complement_regions(diagram)[1].boundary_cycles == [1]
    @test H.intersection_matrix(diagram) == [1;;]
    presentation = H.fundamental_group_presentation(diagram)
    @test H.abelianization_matrix(presentation) == H.intersection_matrix(diagram)'
    @test H.first_homology(diagram).free_rank == 0
    @test isempty(H.first_homology(diagram).torsion)
    @test H.is_homology_sphere(diagram)
    @test H.smith_normal_form([2 4; 6 8]) == [2 0; 0 4]

    lens = lens_fixture(5)
    @test H.first_homology(lens).torsion == [5]
    @test H.intersection_matrix(lens) == [5;;]
    @test !H.is_homology_sphere(lens)

    s2xs1_code = H.heegaard_code([Int[]], [Int[]], H.Intersection[])
    @test !H.validate(H.decorate(s2xs1_code,nothing;cellular=true)).valid
end

@testset "Representation canonicalization" begin
    first_code = H.heegaard_code([[1], [2]], [[1], [2]],
        [H.Intersection(1, 1, 1, 1), H.Intersection(2, 2, 2, -1)])
    second_code = H.heegaard_code(
        [H.Curve(:alpha, 9, [20]), H.Curve(:alpha, 3, [10])],
        [H.Curve(:beta, 8, [10]), H.Curve(:beta, 4, [20])],
        [H.Intersection(20, 9, 4, 1), H.Intersection(10, 3, 8, -1)])
    first = H.decorate(first_code, [[1, 2]])
    second = H.decorate(second_code, [[1, 2]])
    @test H.canonical_hash(first) == H.canonical_hash(second)
    @test H.canonical_hash(H.decorate(first_code, [[2], [1]])) == H.canonical_hash(H.decorate(first_code, [[1], [2]]))
    @test H.canonical_hash(H.canonicalize_representation(second)) == H.canonical_hash(second)
end

@testset "Moves and certificates" begin
    diagram = s3_fixture()
    reversed, epsilon_certificate = H.reverse_curve_orientation(diagram, (:alpha, 1))
    @test H.first_homology(reversed).free_rank == 0
    @test H.inverse_move(epsilon_certificate)["move"] == "reverse_curve_orientation"

    added, add_certificate = H.add_canceling_pair(diagram,
        (alpha_curve=1, beta_curve=1, alpha_position=1, beta_position=1))
    @test H.is_valid(added)
    @test H.first_homology(added).free_rank == 0
    removed, _ = H.remove_canceling_pair(added, [2, 3])
    @test H.canonical_hash(removed) == H.canonical_hash(diagram)
    @test H.inverse_move(add_certificate)["move"] == "remove_canceling_pair"

    stabilized, _ = H.stabilize(diagram)
    @test H.surface_genus(stabilized) == 2
    @test H.first_homology(stabilized).torsion == Int[]
    recovered, _ = H.destabilize(stabilized, [2, 2, 2])
    @test H.canonical_hash(recovered) == H.canonical_hash(diagram)
    @test_throws ArgumentError H.handle_slide(diagram, (source_curve=(:alpha, 1), target_curve=(:alpha, 1)))
    @test !isempty(H.enumerate_local_moves(diagram; kind=:epsilon))
end

@testset "Alpha-cut surface and ports" begin
    diagram = s3_fixture()
    cut = H.cut_along_alpha(diagram)
    placement = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, placement)
    @test length(placement) == 2
    @test length(ports) == 2
    wiring = H.build_wiring_problem(cut, ports)
    @test wiring.wires == [(beta_id=1, from_port=1, to_port=2)]
end

@testset "Poincare homology-sphere source fixture" begin
    fixture_path = joinpath(@__DIR__, "..", "fixtures", "poincare_homology_sphere_v1.hgd.json")
    fixture = H.load_heegaard(fixture_path)
    diagram = H.certify_heegaard(fixture)
    @test H.validate(diagram).valid
    @test H.surface_genus(diagram) == 2
    @test H.intersection_matrix(diagram) == [-3 2; 2 -1]
    presentation = H.fundamental_group_presentation(diagram)
    @test presentation.relators == [[(1, -1), (1, -1), (1, -1), (1, -1), (2, 1), (1, 1), (2, 1)],
                                     [(1, 1), (2, -1), (2, -1), (1, 1), (2, 1)]]
    @test H.abelianization_matrix(presentation) == [-3 2; 2 -1]
    @test H.first_homology(diagram).free_rank == 0
    @test isempty(H.first_homology(diagram).torsion)
    @test H.is_homology_sphere(diagram)
    @test H.source_manifest(diagram)["region_decoration_status"] == "candidate_not_source_transcribed"
    diagram_hash = H.canonical_hash(diagram)
    @test diagram_hash == "691e24337be7458d8b8309229d22bdb1985c137c628312a26622a5f58625dcc0"
    cut = H.cut_along_alpha(diagram; diagram_hash=diagram_hash)
    ports = H.place_ports(cut, H.default_puncture_placement(cut))
    @test length(H.build_wiring_problem(diagram, cut, ports).wires) > 0

    temporary = mktempdir()
    roundtrip_path = joinpath(temporary, "poincare.hgd.json")
    H.save_heegaard(roundtrip_path, diagram)
    roundtrip = H.load_heegaard(roundtrip_path)
    @test H.canonical_hash(roundtrip) == H.canonical_hash(diagram)
end

