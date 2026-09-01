using Topology
using JSON3

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIXTURE = joinpath(ROOT, "test", "fixtures", "poincare_homology_sphere_v1.hgd.json")
const LAYOUT = joinpath(ROOT, "test", "fixtures", "poincare_homology_sphere_v1.hgd.layout.json")
const OUTPUT = joinpath(ROOT, "artifacts", "heegaard", "poincare_smoothing_candidates")

diagram = H.certify_heegaard(H.load_heegaard(FIXTURE))
feasible = H.load_layout(LAYOUT, diagram)
mkpath(OUTPUT)
H.render(feasible; backend=:svg, path=joinpath(OUTPUT, "poincare-baseline.svg"))

candidate = H.smooth_layout(diagram, feasible;
    iterations=48, max_spacing=0.15, simplify_tolerance=0.04, flatness=1e-4)
stem = joinpath(OUTPUT, "poincare-relaxed")
H.save_layout(stem * ".hgd.layout.json", candidate)
H.render(candidate; backend=:svg, path=stem * ".svg")
smoothing = candidate.backend_trace["smoothing"]
summary = [Dict{String,Any}(
    "name" => "relaxed",
    "routes" => length(candidate.routes),
    "dense_points" => smoothing["dense_points"],
    "knots" => smoothing["simplified_points"],
    "bezier_segments" => smoothing["bezier_segments"],
    "control_points" => smoothing["control_points"],
    "accepted_steps" => smoothing["accepted_steps"],
    "initial_energy" => smoothing["initial_energy"],
    "final_energy" => smoothing["final_energy"],
    "minimum_route_separation" => smoothing["minimum_flattened_separation"],
    "configured_minimum_route_separation" =>
        smoothing["parameters"]["minimum_route_separation"],
    "crossings" => length(candidate.crossings),
    "valid" => candidate.report.valid,
    "layout_file" => basename(stem * ".hgd.layout.json"),
    "svg_file" => basename(stem * ".svg"))]

open(joinpath(OUTPUT, "summary.json"), "w") do io
    JSON3.write(io, Dict("source_layout" => basename(LAYOUT), "candidates" => summary))
    write(io, '\n')
end
println(JSON3.write(Dict("output" => OUTPUT, "candidates" => summary)))
