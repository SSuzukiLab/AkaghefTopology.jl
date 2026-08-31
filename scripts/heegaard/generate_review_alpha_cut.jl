using Topology

const H = Topology.Heegaard

function lens_fixture(p, q)
    gcd(p, q) == 1 || throw(ArgumentError("lens-space parameters must be coprime"))
    beta_word = [mod1(1 + (index - 1) * q, p) for index in 1:p]
    code = H.heegaard_code([collect(1:p)], [beta_word],
        [H.Intersection(index, 1, 1, 1) for index in 1:p])
    H.certify_heegaard(H.decorate(code, [[index] for index in 1:p];
        metadata=Dict("fixture" => "L($p,$q)", "status" => "generated-formula")))
end

function s3_fixture()
    code = H.heegaard_code([[1]], [[1]], [H.Intersection(1, 1, 1, 1)])
    H.certify_heegaard(H.decorate(code, [[1]];
        metadata=Dict("fixture" => "S3", "status" => "generated-formula")))
end

function write_layout(output, label, diagram)
    layout = H.layout_heegaard(diagram; backend=:pure_julia, router=:feasible)
    layout.report.valid || error("$label did not produce a valid layout")
    isempty(layout.crossings) || error("$label produced route crossings")
    H.render(layout; backend=:svg, path=joinpath(output, "$(label)__alpha_cut.svg"))
    H.save_layout(joinpath(output, "$(label)__alpha_cut.hgd.layout.json"), layout)
    relaxed = H.smooth_layout(diagram, layout;
        iterations=12, max_spacing=0.15, simplify_tolerance=0.04, flatness=1e-4)
    H.render(relaxed; backend=:svg,
             path=joinpath(output, "$(label)__alpha_cut__relaxed.svg"))
    H.save_layout(joinpath(output, "$(label)__alpha_cut__relaxed.hgd.layout.json"), relaxed)
    metrics = relaxed.backend_trace["smoothing"]
    println(label, " routes=", length(layout.routes),
            " points=", sum(length(route.points) for route in layout.routes),
            " router=", layout.backend_trace["post_router"],
            " dense=", metrics["dense_points"],
            " knots=", metrics["simplified_points"],
            " controls=", metrics["control_points"])
end

function main()
    length(ARGS) == 1 || error("usage: generate_review_alpha_cut.jl OUTPUT_DIR")
    output = abspath(ARGS[1])
    mkpath(output)
    write_layout(output, "S3", s3_fixture())
    for (p, q) in ((2, 1), (5, 1), (5, 2), (7, 1), (7, 2))
        write_layout(output, "L$(p)_$(q)", lens_fixture(p, q))
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
