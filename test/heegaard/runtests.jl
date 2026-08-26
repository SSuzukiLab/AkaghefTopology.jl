using Test
using Topology

const H = Topology.Heegaard

function s3_fixture()
    code = H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)])
    H.certify_heegaard(H.decorate(code, [[1]]; metadata=Dict("fixture"=>"s3")))
end

function lens_fixture(p)
    code = H.heegaard_code([collect(1:p)], [collect(1:p)],
        [H.Intersection(i, 1, 1, 1) for i in 1:p])
    H.certify_heegaard(H.decorate(code, [[i] for i in 1:p]; metadata=Dict("fixture"=>"lens")))
end

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
    s2xs1 = H.certify_heegaard(H.decorate(s2xs1_code, nothing; cellular=true))
    @test H.surface_genus(s2xs1) == 1
    @test H.first_homology(s2xs1).free_rank == 1
    @test H.first_homology(s2xs1).torsion == Int[]
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

@testset "Alpha-cut layout and semantic SVG" begin
    diagram = s3_fixture()
    cut = H.cut_along_alpha(diagram)
    placement = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, placement)
    @test length(placement) == 2
    @test length(ports) == 1
    wiring = H.build_wiring_problem(cut, ports)
    @test wiring.wires == [(beta_id=1, from_port=1, to_port=1)]
    layout = H.layout_heegaard(diagram)
    @test layout.view == :alpha_cut
    @test layout.report.valid
    svg = H.render(layout)
    @test occursin("data-intersection-id=\"1\"", svg)
    @test occursin("data-curve-family=\"beta\"", svg)
    @test isempty(layout.crossings)
    @test H.validate_layout(diagram, layout).valid
    @test H.gauss_layout(diagram).view == :gauss

    temporary = mktempdir()
    path = joinpath(temporary, "s3.hgd.layout.json")
    H.save_layout(path, layout)
    loaded = H.load_layout(path, diagram)
    @test loaded.report.valid
    @test_throws ArgumentError H.load_layout(path, lens_fixture(5))
end

@testset "Typed spine conversion" begin
    diagram = s3_fixture()
    punctured = H.make_punctured(diagram)
    spine = H.to_flow_spine(punctured)
    ograph = H.to_ograph(spine)
    bundle = H.to_odata(diagram)
    @test punctured isa H.PuncturedHeegaardDiagram
    @test spine isa H.FlowSpine
    @test ograph isa H.OData
    @test bundle.ograph isa H.OData
    @test bundle.certificate.evidence_status == :unverified
    @test H.to_heegaard(spine) === diagram
    @test bundle.comparison["canonical_conversion"] == false
end
