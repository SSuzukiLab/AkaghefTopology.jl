using Topology

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "artifacts", "plots")
const JULIA_PLOTS = joinpath(ROOT, "artifacts", "julia_plots")
const CONTACT_SHEETS = joinpath(ROOT, "artifacts", "julia_contact_sheets")

function render_acceptance_artifacts(; expected_count::Int=46)
    rsvg = Sys.which("rsvg-convert")
    magick = Sys.which("magick")
    rsvg === nothing && error("rsvg-convert is required to render the SVG acceptance artifacts")
    magick === nothing && error("ImageMagick is required to build the visual-review contact sheets")

    svg_files = sort!([joinpath(dir, file) for (dir, _, files) in walkdir(JULIA_PLOTS)
                       for file in files if endswith(file, ".svg")])
    length(svg_files) == expected_count ||
        error("expected $expected_count dynamic Julia SVGs, found $(length(svg_files))")

    for svg in svg_files
        png = splitext(svg)[1] * ".png"
        run(Cmd([rsvg, "-o", png, svg]))
        filesize(png) > 0 || error("empty rendered PNG: $png")
    end

    mkpath(CONTACT_SHEETS)
    plot_directories = sort!(unique(dirname.(svg_files)))
    for directory in plot_directories
        png_files = sort!(filter(file -> endswith(file, ".png"), readdir(directory; join=true)))
        relative = relpath(directory, JULIA_PLOTS)
        output = joinpath(CONTACT_SHEETS, replace(relative, Base.Filesystem.path_separator => '_') * ".png")
        command = [magick, "montage"]
        append!(command, png_files)
        append!(command, ["-thumbnail", "700x500", "-tile", "4x", "-geometry", "+12+12",
                          "-background", "white", output])
        run(Cmd(command))
        filesize(output) > 0 || error("empty contact sheet: $output")
    end
    return length(svg_files), length(plot_directories)
end

for file in ["VLExample.jl","OgraphPlotSuite.jl","C240610URD.jl","C240623URD.jl","URDforKnotCmpl.jl","dev_URDiagram.jl"]
    include(joinpath(ROOT,"examples",file))
end
run_vl_plot_calls(output_dir=joinpath(JULIA_PLOTS,"topology","Manifold","VLExample"))
run_ograph_plot_suite(output_root=JULIA_PLOTS)
run_C240610(output_dir=OUT)
run_C240623(output_dir=OUT)
run_URDforKnotCmpl(output_dir=OUT)
run_dev_urdiagram(output_dir=OUT)
plot_count, sheet_count = render_acceptance_artifacts()
println("generated and rasterized $plot_count dynamic Julia plots; built $sheet_count contact sheets")
println("generated $(length(readdir(OUT))) supplemental workbench files in $OUT")
