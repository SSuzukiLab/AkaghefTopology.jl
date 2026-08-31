using Topology
using JSON3

const H = Topology.Heegaard
const ROOT = normpath(joinpath(@__DIR__, ".."))
const FIXTURE = joinpath(ROOT, "test", "fixtures", "poincare_homology_sphere_v1.hgd.json")
const LAYOUT = joinpath(ROOT, "test", "fixtures", "poincare_homology_sphere_v1.hgd.layout.json")
const OUTPUT = joinpath(ROOT, "artifacts", "heegaard", "poincare_smoothing_candidates")

diagram = H.certify_heegaard(H.load_heegaard(FIXTURE))
feasible = H.layout_heegaard(diagram; backend=:catlab, fallback=:error,
                             route_certificate=LAYOUT)
mkpath(OUTPUT)
H.render(feasible; backend=:svg, path=joinpath(OUTPUT, "poincare-baseline.svg"))

specifications = [("tight", 0.08), ("balanced", 0.14), ("wide", 0.22)]
summary = Any[]
for (label, step_fraction) in specifications
    candidate = H.smooth_layout(diagram, feasible; step_fraction, flatness=1e-4)
    stem = joinpath(OUTPUT, "poincare-$(label)")
    H.save_layout(stem * ".hgd.layout.json", candidate)
    H.render(candidate; backend=:svg, path=stem * ".svg")
    smoothing = candidate.backend_trace["smoothing"]
    metrics = smoothing["route_metrics"]
    push!(summary, Dict{String,Any}(
        "name" => label,
        "step_fraction" => step_fraction,
        "routes" => length(candidate.routes),
        "flattened_points" => sum(route -> length(route.points), candidate.routes),
        "smoothed_corners" => sum(item -> item["smoothed_corners"], metrics),
        "max_step" => maximum(item["max_step"] for item in metrics),
        "crossings" => length(candidate.crossings),
        "self_crossings" => count(issue -> issue.code == :self_crossing, candidate.report.issues),
        "puncture_intrusions" => count(issue -> issue.code == :puncture_intrusion, candidate.report.issues),
        "valid" => candidate.report.valid,
        "layout_file" => basename(stem * ".hgd.layout.json"),
        "svg_file" => basename(stem * ".svg")))
end

open(joinpath(OUTPUT, "summary.json"), "w") do io
    JSON3.write(io, Dict("source_layout" => basename(LAYOUT), "candidates" => summary))
    write(io, '\n')
end
println(JSON3.write(Dict("output" => OUTPUT, "candidates" => summary)))
