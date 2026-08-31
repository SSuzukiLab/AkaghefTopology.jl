using Test
using Topology

if !isdefined(@__MODULE__, :H)
    const H = Topology.Heegaard
end

function _router_s3_fixture()
    code = H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)])
    H.certify_heegaard(H.decorate(code, [[1]]; metadata=Dict("fixture" => "s3")))
end

function _router_lens_fixture(p, q=1)
    gcd(p, q) == 1 || throw(ArgumentError("lens-space parameters must be coprime"))
    beta_word = [mod1(1 + (index - 1) * q, p) for index in 1:p]
    code = H.heegaard_code([collect(1:p)], [beta_word],
        [H.Intersection(i, 1, 1, 1) for i in 1:p])
    H.certify_heegaard(H.decorate(code, [[i] for i in 1:p]; metadata=Dict("fixture" => "lens")))
end

function _verify_router_routes(diagram; anchors=nothing)
    cut = H.cut_along_alpha(diagram)
    punctures = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, punctures)
    wiring = H.build_wiring_problem(diagram, cut, ports)
    routes = H.feasible_fixed_endpoint_routes(diagram, punctures, ports, wiring;
                                               anchors=anchors, max_attempts=2)
    candidate = H.HeegaardLayout(H.canonical_hash(diagram), :alpha_cut,
                                 Dict{String,Any}(), punctures, ports, routes,
                                 H._route_crossings(routes), Dict{String,Any}(),
                                 H.LayoutValidationReport(false, H.ValidationIssue[],
                                                          Dict{Symbol,Bool}()))
    report = H.validate_layout(diagram, candidate)
    @test report.valid
    @test isempty(candidate.crossings)
    @test all(route.family == :beta for route in routes)
    @test length(routes) == length(wiring.wires)
    routes
end

@testset "Bounded fixed-endpoint feasible router" begin
    poincare = H.certify_heegaard(H.load_heegaard(joinpath(@__DIR__, "..", "fixtures",
                                                            "poincare_homology_sphere_v1.hgd.json")))
    poincare_routes = _verify_router_routes(poincare)
    @test sum(length(route.points) for route in poincare_routes) > 24

    s3 = _router_s3_fixture()
    s3_layout = H.layout_heegaard(s3; backend=:pure_julia, router=:feasible,
                                  anchors=[(0.0, 0.0)], max_attempts=2)
    @test s3_layout.report.valid
    @test s3_layout.backend_trace["post_router"] == "bounded-fixed-endpoint-feasible"
    @test s3_layout.backend_trace["feasible_router"]
    @test s3_layout.backend_trace["router_anchor_count"] == 1

    for (p, q) in ((2, 1), (5, 1), (5, 2), (7, 1), (7, 2))
        lens_routes = _verify_router_routes(_router_lens_fixture(p, q))
        @test length(lens_routes) == p
    end
    cut = H.cut_along_alpha(s3)
    punctures = H.default_puncture_placement(cut)
    ports = H.place_ports(cut, punctures)
    wiring = H.build_wiring_problem(cut, ports)
    @test_throws ArgumentError H.feasible_fixed_endpoint_routes(
        s3, punctures, ports, wiring; max_attempts=1, max_expansions=0)
end
