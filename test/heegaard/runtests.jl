using Test
using Topology
using JSON

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
    @test length(ports) == 2
    wiring = H.build_wiring_problem(cut, ports)
    @test wiring.wires == [(beta_id=1, from_port=1, to_port=2)]
    layout = H.layout_heegaard(diagram)
    @test layout.view == :alpha_cut
    @test layout.report.valid
    @test layout.backend_trace["backend"] == "catlab"
    @test layout.backend_trace["catlab_used"] == true
    @test layout.backend_trace["catlab_layout"]["ordered_port_count"] == 2
    @test layout.backend_trace["pipeline"]["puncture_placement"] == "fixed-2-by-g"
    @test layout.backend_trace["pipeline"]["feasible_router"] ==
          "bounded-fixed-endpoint-feasible"
    @test layout.backend_trace["post_router"] == "bounded-fixed-endpoint-feasible"
    svg = H.render(layout)
    @test occursin("data-intersection-id=\"1\"", svg)
    @test occursin("data-curve-family=\"beta\"", svg)
    @test occursin("fill-opacity=\"0.14\"", svg)
    @test occursin("vector-effect=\"non-scaling-stroke\"", svg)
    @test isempty(layout.crossings)
    @test H.validate_layout(diagram, layout).valid
    @test H.gauss_layout(diagram).view == :gauss

    pure_layout = H.layout_heegaard(diagram; backend=:pure_julia)
    @test pure_layout.backend_trace["catlab_used"] == false
    @test pure_layout.backend_trace["backend"] == "pure_julia"

    temporary = mktempdir()
    path = joinpath(temporary, "s3.hgd.layout.json")
    H.save_layout(path, layout)
    loaded = H.load_layout(path, diagram)
    @test loaded.report.valid
    @test_throws ArgumentError H.load_layout(path, lens_fixture(5))
end

@testset "Alpha-cut validator recomputes semantic geometry" begin
    diagram = lens_fixture(2)
    cut = H.cut_along_alpha(diagram)
    punctures = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, punctures)
    routes = H._alpha_cut_routes(diagram, cut, punctures, ports)
    diagram_hash = H.canonical_hash(diagram)
    empty_report = H.LayoutValidationReport(false, H.ValidationIssue[], Dict{Symbol,Bool}())

    wrong_successors = [
        H.Route(:beta, 1, 1, 1, 4, [ports[1].point, ports[4].point]),
        H.Route(:beta, 1, 2, 2, 3, [ports[2].point, ports[3].point]),
    ]
    wrong_layout = H.HeegaardLayout(diagram_hash, :alpha_cut, Dict{String,Any}(),
                                    punctures, ports, wrong_successors,
                                    Tuple{Int,Int}[], Dict{String,Any}(), empty_report)
    wrong_report = H.validate_layout(diagram, wrong_layout; expected_hash=diagram_hash)
    @test any(issue -> issue.code == :beta_successor, wrong_report.issues)

    s3 = s3_fixture()
    s3_cut = H.cut_along_alpha(s3)
    s3_punctures = H.default_puncture_placement(s3_cut)
    s3_ports = H.place_ports(s3_cut, s3_punctures)
    start, finish = s3_ports[1].point, s3_ports[2].point
    self_crossing = H.Route(:beta, 1, 1, 1, 2,
        [start, (0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0), finish])
    self_layout = H.HeegaardLayout(H.canonical_hash(s3), :alpha_cut, Dict{String,Any}(),
                                   s3_punctures, s3_ports, [self_crossing],
                                   Tuple{Int,Int}[], Dict{String,Any}(), empty_report)
    self_report = H.validate_layout(s3, self_layout; expected_hash=H.canonical_hash(s3))
    @test any(issue -> issue.code == :self_crossing, self_report.issues)
    @test any(issue -> issue.code == :puncture_intrusion, self_report.issues)

    @test !isempty(H._route_crossings(routes))
    crossing_layout = H.HeegaardLayout(diagram_hash, :alpha_cut, Dict{String,Any}(),
                                       punctures, ports, routes,
                                       Tuple{Int,Int}[], Dict{String,Any}(), empty_report)
    crossing_report = H.validate_layout(diagram, crossing_layout; expected_hash=diagram_hash)
    @test any(issue -> issue.code == :crossing_mismatch, crossing_report.issues)
end

@testset "Poincare homology-sphere source fixture" begin
    fixture_path = joinpath(@__DIR__, "..", "fixtures", "poincare_homology_sphere_v1.hgd.json")
    layout_fixture_path = joinpath(@__DIR__, "..", "fixtures",
                                   "poincare_homology_sphere_v1.hgd.layout.json")
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
    placement = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, placement)
    terminal = H._terminal_ribbon_report(diagram, cut, ports,
        H.build_wiring_problem(diagram, cut, ports))
    @test terminal.planar
    @test (terminal.vertices, terminal.edges, terminal.faces, terminal.euler, terminal.genus) ==
          (4, 12, 10, 2, 0)
    @test terminal.face_lengths == [2, 2, 2, 2, 2, 2, 3, 3, 3, 3]
    @test terminal.boundary_cycle_bijection
    @test sort(terminal.face_to_boundary_cycle) == collect(1:10)

    temporary = mktempdir()
    roundtrip_path = joinpath(temporary, "poincare.hgd.json")
    H.save_heegaard(roundtrip_path, diagram)
    roundtrip = H.load_heegaard(roundtrip_path)
    @test H.canonical_hash(roundtrip) == H.canonical_hash(diagram)

    generated = H.layout_heegaard(diagram; backend=:pure_julia)
    @test generated.report.valid
    @test generated.backend_trace["post_router"] == "bounded-fixed-endpoint-feasible"
    @test isempty(generated.crossings)
    layout = H.layout_heegaard(diagram; backend=:catlab, fallback=:error,
                               route_certificate=layout_fixture_path)
    @test layout.report.valid
    @test layout.backend_trace["catlab_used"] == true
    @test layout.backend_trace["post_router"] == "certified-fixed-boundary-pl"
    @test layout.backend_trace["route_certificate"]["file"] ==
          "poincare_homology_sphere_v1.hgd.layout.json"
    @test layout.backend_trace["terminal_embedding"]["planar"] == true
    @test isempty(layout.crossings)
    @test sum(length(route.points) for route in layout.routes) == 346
    replay = H.load_layout(layout_fixture_path, diagram)
    @test replay.report.valid
    @test replay.diagram_hash == diagram_hash
    @test all((first.family, first.curve_id, first.segment_id,
               first.start_port, first.end_port, first.points) ==
              (second.family, second.curve_id, second.segment_id,
               second.start_port, second.end_port, second.points)
              for (first, second) in zip(replay.routes, layout.routes))
    svg = H.render(layout; backend=:svg)
    @test occursin("data-heegaard-view=\"alpha_cut\"", svg)
    @test occursin("data-diagram-hash=\"$diagram_hash\"", svg)
    @test occursin("candidate_not_source_transcribed", svg)
    @test occursin("data-curve-family=\"beta\"", svg)
end

include(joinpath(@__DIR__, "feasible_router.jl"))

@testset "Mechanical smoothing preserves the certified Poincare layout" begin
    fixture_path = joinpath(@__DIR__, "..", "fixtures", "poincare_homology_sphere_v1.hgd.json")
    layout_fixture_path = joinpath(@__DIR__, "..", "fixtures",
                                   "poincare_homology_sphere_v1.hgd.layout.json")
    diagram = H.certify_heegaard(H.load_heegaard(fixture_path))
    feasible = H.layout_heegaard(diagram; backend=:catlab, fallback=:error,
                                 route_certificate=layout_fixture_path)
    smoothed = H.smooth_layout(diagram, feasible; step_fraction=0.14, flatness=1e-4)

    @test smoothed.report.valid
    @test smoothed.crossings == Tuple{Int,Int}[]
    @test smoothed.report.checks[:curve_representation]
    @test all(route.curve isa H.BezierChain for route in smoothed.routes)
    @test all(first(route.points) == H._point_for_port(smoothed.ports, route.start_port) &&
              last(route.points) == H._point_for_port(smoothed.ports, route.end_port)
              for route in smoothed.routes)
    @test any(route -> length(route.points) > 346 / length(smoothed.routes), smoothed.routes)
    @test occursin(" C ", H.render(smoothed; backend=:svg))

    temporary = mktempdir()
    path = joinpath(temporary, "poincare-smoothed.hgd.layout.json")
    H.save_layout(path, smoothed)
    encoded = JSON.parsefile(path)
    @test encoded["routes"][1]["curve"]["type"] == "bezier_chain"
    loaded = H.load_layout(path, diagram)
    @test loaded.report.valid
    @test all(route.curve isa H.BezierChain for route in loaded.routes)
    @test all(first == second for (first, second) in zip(
        [route.points for route in smoothed.routes], [route.points for route in loaded.routes]))

    tampered_routes = deepcopy(feasible.routes)
    tampered_routes[1].points[2] = (tampered_routes[1].points[2][1] + 0.01,
                                    tampered_routes[1].points[2][2])
    tampered = H.HeegaardLayout(feasible.diagram_hash, feasible.view, feasible.chart,
                                feasible.punctures, feasible.ports, tampered_routes,
                                feasible.crossings, feasible.backend_trace, feasible.report)
    tampered_report = H.validate_layout(diagram, tampered)
    @test any(issue -> issue.code == :curve_points_mismatch, tampered_report.issues)

    malformed_path = joinpath(temporary, "malformed.hgd.layout.json")
    malformed = deepcopy(encoded)
    malformed["routes"][1]["curve"]["type"] = "unknown"
    open(malformed_path, "w") do io
        JSON.print(io, malformed)
    end
    @test_throws ArgumentError H.load_layout(malformed_path, diagram)
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
