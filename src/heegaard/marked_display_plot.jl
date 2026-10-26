"""
    plot(view::AbstractDict; minimum_separation=12, puncture_clearance=12,
         crossing_radius=20, maxiter=1800)

Smooth an explicit `marked-cut-view/1` embedding on a two-by-genus circle grid.
The original oriented curve intervals and joining-disk boundary are retained.
Boundary vertex pairs are mandatory: independent endpoint routing is forbidden.
Split D pieces share those curve intervals and oriented cut-circle arcs; the
source infinity triangle is a separate retained display choice.
The SciPy backend is shared by all genera and curve families. Its numerical
proposal is accepted only after independent source, Bezier and clearance checks.
This overload does not reconstruct a missing marking or a missing feasible view.
See scripts/heegaard/marked_display_pipeline.py::plot and check_smooth_marked.py.
"""
function plot(view::AbstractDict; minimum_separation=12.,puncture_clearance=12.,
              crossing_radius=20.,maxiter=1800,cut_curves=nothing,
              outer_triangle=nothing,visible_curves=nothing,_source_map_geometry=false)
    # Native input has no accepted embedded display to preserve. Its ordinary
    # public plot path must construct the full source map; legacy embedded
    # views retain their existing default unless they carry a saved source map.
    fullmap=_source_map_geometry || haskey(view,"source_map") || get(view,"schema",nothing)!="marked-cut-view/1"
    if get(view,"schema",nothing)!="marked-cut-view/1"
        prepared=layout_marked(view;cut_curves=cut_curves,outer_triangle=outer_triangle,visible_curves=visible_curves,_source_map_geometry=fullmap)
        return plot(prepared;minimum_separation=minimum_separation,puncture_clearance=puncture_clearance,
            crossing_radius=crossing_radius,maxiter=maxiter,_source_map_geometry=fullmap)
    end
    all(isfinite(x) && x>0 for x in (minimum_separation,puncture_clearance,crossing_radius)) ||
        throw(ArgumentError("Drawing clearances and crossing radius must be positive finite numbers"))
    maxiter isa Integer && maxiter>0 || throw(ArgumentError("maxiter must be a positive integer"))
    # Source: marked_display.jl::marked_display_source re-prepares the same
    # native object. Infinity is checked separately from its representation hash.
    haskey(view,"source_document") && marked_display_source(view)
    uv=Sys.which("uv");uv===nothing && error("uv is required for the shared SciPy drawing backend")
    root=normpath(joinpath(@__DIR__,"../.."));script=joinpath(root,"scripts/heegaard/marked_display_pipeline.py")
    mktempdir() do dir
      try
        input=joinpath(dir,"input.json");output=joinpath(dir,"output.json")
        write(input,JSON3.write(view))
        command=`$uv run --no-project --offline --python 3.12 --with numpy --with scipy --with shapely --with networkx python $script $input $output --minimum-separation $minimum_separation --puncture-clearance $puncture_clearance --crossing-radius $crossing_radius --maxiter $maxiter`
        proof=joinpath(dir,"source-map.json")
        if fullmap
            # Source: the existing prepare_marked_display boundary/cellulation
            # producer and bounded Python recipe decoder. Re-prepare the actual
            # source and effective infinity; saved proof booleans are ignored.
            p=prepare_marked_display(view["source_document"];cut_curves=view["cuts"],visible_curves=view["visible_curves"])
            p["source_document"]=view["source_document"]
            p["outer_triangle"]=view["display_request"]["infinity_triangle"]
            p["display_request"]=view["display_request"]
            prepared_file=joinpath(dir,"prepared.json");write(prepared_file,JSON3.write(p))
            command=`$uv run --no-project --offline --python 3.12 --with numpy==2.5.4 --with scipy --with shapely==2.1.2 --with networkx --with triangle==20250106 python $script $input $output --minimum-separation $minimum_separation --puncture-clearance $puncture_clearance --crossing-radius $crossing_radius --maxiter $maxiter --prepared $prepared_file --source-map-geometry --source-map-proof-output $proof`
        end
        diagnostic=IOBuffer()
        process=run(pipeline(ignorestatus(addenv(command,"OPENBLAS_NUM_THREADS"=>"1")),stdout=stdout,stderr=diagnostic))
        success(process) || error("Marked drawing rejected: "*last(String(take!(diagnostic)),3000))
        result=JSON3.read(read(output,String),Dict{String,Any})
        if fullmap
            # The helper must be included by the module owner before this
            # private branch is enabled. A transport envelope is not proof.
            0<filesize(proof)<=32*1024*1024 || error("Invalid source-map decoded size")
            raw=read(proof);envelope=result["source_map"]
            length(raw)==envelope["uncompressed_bytes"] && bytes2hex(SHA.sha256(raw))==envelope["sha256"] || error("Source-map proof bytes differ from the output envelope")
            recipe=JSON3.read(String(raw),Dict{String,Any})
            recipe["display_request"]==result["display_request"]==view["display_request"] || error("Source-map effective marking changed")
            recipe["source_representation_sha256"]==result["preparation_representation_sha256"] || error("Source-map representation changed")
            result["source_document"]==view["source_document"] || error("Source-map native document changed")
            bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(result["svg"]))))==recipe["output_svg_sha256"] || error("Source-map SVG differs from the validated output")
            result["source_map_cell_check"]=_mdsm_verify_recipe_cells(recipe)
            delete!(result,"source_map_verification_pending")
            # Source: the existing /api/render saved-view envelope. Count the
            # actual compact final object AFTER adding Julia's cell result;
            # bounding only compressed proof data misses the rest of the view.
            saved_request_bytes=ncodeunits(JSON3.write(Dict("diagram"=>result)))
            0<saved_request_bytes<=5_000_000 || error("Saved marked-view request exceeds the existing 5000000-byte limit: "*string(saved_request_bytes))
        end
        result
      catch err
        if fullmap
            # Source: Base.mktempdir(cleanup=false)/cp. Preserve the caller's
            # actual completed stage files before the scoped temporary input
            # directory is cleaned up; a late assertion must not erase evidence.
            retained=mktempdir(;prefix="marked-source-map-failure-",cleanup=false)
            cp(dir,joinpath(retained,"stages"))
            error(sprint(showerror,err)*"; source-map stage files retained at "*retained)
        end
        rethrow()
      end
    end
end

"""Construct an initial 2×g view from the whole marked cellulation.

The auxiliary cut curves and the source triangle containing infinity must be
specified, either as keywords or in metadata.display. No region is selected for
readability. This bounded numerical realization returns only checked meshes;
failure to find one is explicit, not a claim of topological impossibility.
Explicit keywords override metadata.display. The unmodified source_document is
retained; display_request and infinity_source_triangle record the effective
selection and must be used by consumers replaying operations on this view.
See scripts/heegaard/layout_marked_document.py and marked_cut_layout.py.
"""
function layout_marked(input::AbstractDict;cut_curves=nothing,outer_triangle=nothing,visible_curves=nothing,_source_map_geometry=false)
    doc=_plain_value(input)
    options=get(get(doc,"metadata",Dict()),"display",Dict())
    cuts=cut_curves===nothing ? get(options,"cut_curves",nothing) : cut_curves
    outer=outer_triangle===nothing ? get(options,"infinity_triangle",nothing) : outer_triangle
    visible=visible_curves===nothing ? get(options,"visible_curves",nothing) : visible_curves
    cuts===nothing && throw(ArgumentError("Specify the auxiliary cut_curves; alpha is not inferred as the display cut"))
    outer===nothing && throw(ArgumentError("Specify infinity_triangle in the source cellulation; the exterior region is not selected by the renderer"))
    outer isa AbstractVector && length(outer)==2 && outer[1] isa AbstractString && outer[2] isa Integer && !(outer[2] isa Bool) ||
        throw(ArgumentError("infinity_triangle must be [patch name, zero-based triangle index]"))
    p=prepare_marked_display(doc;cut_curves=cuts,visible_curves=visible)
    any(t["source"]==outer for t in p["triangles"]) || throw(ArgumentError("Unknown infinity_triangle in source cellulation"))
    p["source_document"]=doc;p["outer_triangle"]=collect(outer)
    p["display_request"]=Dict("cut_curves"=>p["cut_curves"],"visible_curves"=>p["visible_curves"],
        "infinity_triangle"=>collect(outer))
    uv=Sys.which("uv");uv===nothing && error("uv is required for the shared SciPy drawing backend")
    root=normpath(joinpath(@__DIR__,"../.."));script=joinpath(root,"scripts/heegaard/layout_marked_document.py")
    mktempdir() do dir
        request=joinpath(dir,"prepared.json");response=joinpath(dir,"layout.json")
        write(request,JSON3.write(p));diagnostic=IOBuffer()
        command=`$uv run --no-project --offline --python 3.12 --with numpy --with scipy --with shapely --with networkx python $script $request $response`
        if _source_map_geometry
            command=`$uv run --no-project --offline --python 3.12 --with numpy==2.5.4 --with scipy --with shapely==2.1.2 --with networkx --with triangle==20250106 python $script $request $response --source-map-geometry`
        end
        process=run(pipeline(ignorestatus(addenv(command,"OPENBLAS_NUM_THREADS"=>"1")),stdout=stdout,stderr=diagnostic))
        success(process) || error("Initial marked layout rejected: "*last(String(take!(diagnostic)),3000))
        JSON3.read(read(response,String),Dict{String,Any})
    end
end

"""Render the checked output of the marked-display plot overload."""
function render(view::AbstractDict; backend=:svg)
    get(view,"schema",nothing)=="marked-cut-view/1" || throw(ArgumentError("Expected a marked display"))
    backend==:svg || throw(ArgumentError("Marked display currently supports SVG output"))
    haskey(view,"smooth_check") && haskey(view,"svg") || throw(ArgumentError("Call plot on the explicit embedded layout first"))
    get(view,"source_map_verification_pending",false) && throw(ArgumentError("The full source-map validation has not completed"))
    view["svg"]
end
