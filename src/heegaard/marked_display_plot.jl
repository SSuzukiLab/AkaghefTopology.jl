"""
    plot(view::AbstractDict; minimum_separation=12, puncture_clearance=12,
         crossing_radius=20, maxiter=1800)

Smooth an explicit `marked-cut-view/1` embedding on a two-by-genus circle grid.
The original oriented curve intervals and joining-disk boundary are retained.
Boundary vertex pairs are mandatory: independent endpoint routing is forbidden.
The SciPy backend is shared by all genera and curve families. Its numerical
proposal is accepted only after independent source, Bezier and clearance checks.
This overload does not reconstruct a missing marking or a missing feasible view.
See scripts/heegaard/marked_display_pipeline.py::plot and check_smooth_marked.py.
"""
function plot(view::AbstractDict; minimum_separation=12.,puncture_clearance=12.,
              crossing_radius=20.,maxiter=1800)
    get(view,"schema",nothing)=="marked-cut-view/1" || throw(ArgumentError("Expected an explicit marked-cut-view/1 embedding"))
    all(isfinite(x) && x>0 for x in (minimum_separation,puncture_clearance,crossing_radius)) ||
        throw(ArgumentError("Drawing clearances and crossing radius must be positive finite numbers"))
    maxiter isa Integer && maxiter>0 || throw(ArgumentError("maxiter must be a positive integer"))
    uv=Sys.which("uv");uv===nothing && error("uv is required for the shared SciPy drawing backend")
    root=normpath(joinpath(@__DIR__,"../.."));script=joinpath(root,"scripts/heegaard/marked_display_pipeline.py")
    mktempdir() do dir
        input=joinpath(dir,"input.json");output=joinpath(dir,"output.json")
        write(input,JSON3.write(view))
        command=`$uv run --no-project --offline --python 3.12 --with numpy --with scipy --with shapely --with networkx python $script $input $output --minimum-separation $minimum_separation --puncture-clearance $puncture_clearance --crossing-radius $crossing_radius --maxiter $maxiter`
        diagnostic=IOBuffer()
        process=run(pipeline(ignorestatus(addenv(command,"OPENBLAS_NUM_THREADS"=>"1")),stdout=stdout,stderr=diagnostic))
        success(process) || error("Marked drawing rejected: "*last(String(take!(diagnostic)),3000))
        JSON3.read(read(output,String),Dict{String,Any})
    end
end

"""Render the checked output of the marked-display plot overload."""
function render(view::AbstractDict; backend=:svg)
    get(view,"schema",nothing)=="marked-cut-view/1" || throw(ArgumentError("Expected a marked display"))
    backend==:svg || throw(ArgumentError("Marked display currently supports SVG output"))
    haskey(view,"smooth_check") && haskey(view,"svg") || throw(ArgumentError("Call plot on the explicit embedded layout first"))
    view["svg"]
end
