using Topology, JSON

const ROOT=normpath(joinpath(@__DIR__,".."))
include(joinpath(ROOT,"examples","VLExample.jl"))

matlab_string(path)=replace(abspath(path),"'"=>"''")
nested_int(value)=[Int.(collect(component)) for component in value]

function compare_vg2g(;workdir=mktempdir(),show_debug=false)
    mkpath(workdir)
    input_path=joinpath(workdir,"vg2g-inputs.json")
    matlab_path=joinpath(workdir,"vg2g-matlab.json")
    summary_path=joinpath(workdir,"vg2g-parity.json")

    inputs=vl_virtual_gauss_inputs()
    input_names=[input.name for input in inputs]
    length(unique(input_names))==length(input_names) || error("duplicate vg2g input names")
    open(input_path,"w") do io
        JSON.print(io,[Dict(
            "name"=>input.name,
            "rgauss"=>input.rgauss,
            "orientation"=>input.orientation,
        ) for input in inputs],2)
    end

    matlab=get(ENV,"MATLAB_EXECUTABLE",something(Sys.which("matlab"),""))
    isempty(matlab) && error("MATLAB executable not found; set MATLAB_EXECUTABLE")
    matlab_source=normpath(joinpath(ROOT,"..","topology"))
    matlab_scripts=joinpath(ROOT,"scripts","matlab")
    command="addpath(genpath('$(matlab_string(matlab_source))')); "*
            "addpath('$(matlab_string(matlab_scripts))'); "*
            "export_vg2g_json('$(matlab_string(input_path))','$(matlab_string(matlab_path))');"
    isfile(matlab_path) && rm(matlab_path)
    run(`$matlab -batch $command`)
    isfile(matlab_path) || error("MATLAB did not create $matlab_path")

    payload=JSON.parsefile(matlab_path)
    payload isa AbstractDict || error("MATLAB vg2g output must be a JSON object")
    get(payload,"schema_version",nothing)==1 || error("unsupported MATLAB vg2g schema")
    get(payload,"matlab_version",nothing) isa AbstractString || error("MATLAB version is missing")
    raw_results=get(payload,"results",nothing)
    raw_results isa AbstractVector || error("MATLAB vg2g results must be an array")
    all(item -> item isa AbstractDict && all(haskey(item,key) for key in ("name","gauss","orientation")),
        raw_results) || error("MATLAB vg2g result schema is invalid")
    result_names=[item["name"] for item in raw_results]
    length(unique(result_names))==length(result_names) || error("duplicate MATLAB vg2g result names")
    length(result_names)==length(input_names) && Set(result_names)==Set(input_names) ||
        error("MATLAB result names do not match inputs: $(Set(result_names)) != $(Set(input_names))")
    matlab_results=Dict(item["name"]=>item for item in raw_results)
    summary=Any[]
    for input in inputs
        debug_buffer=IOBuffer()
        debug_io=show_debug ? stderr : debug_buffer
        julia_gauss,julia_orientation=planarize_virtual_gauss(
            input.rgauss,input.orientation;DEBUG_MODE=true,debug_io=debug_io)
        reference=matlab_results[input.name]
        matlab_gauss=nested_int(reference["gauss"])
        matlab_orientation=Int.(reference["orientation"])
        gauss_equal=julia_gauss==matlab_gauss
        orientation_equal=julia_orientation==matlab_orientation
        push!(summary,Dict(
            "name"=>input.name,
            "gauss_equal"=>gauss_equal,
            "orientation_equal"=>orientation_equal,
            "inserted_crossings"=>length(julia_orientation)-length(input.orientation),
            "julia_gauss"=>julia_gauss,
            "matlab_gauss"=>matlab_gauss,
            "julia_orientation"=>julia_orientation,
            "matlab_orientation"=>matlab_orientation,
            "debug_trace_bytes"=>show_debug ? 0 : sizeof(take!(debug_buffer)),
        ))
        println(input.name,": gauss=",gauss_equal ? "MATCH" : "MISMATCH",
                " orientation=",orientation_equal ? "MATCH" : "MISMATCH")
    end
    artifact=Dict(
        "schema_version"=>1,
        "matlab_version"=>payload["matlab_version"],
        "julia_version"=>string(VERSION),
        "results"=>summary,
    )
    open(summary_path,"w") do io
        JSON.print(io,artifact,2)
    end
    all(item["gauss_equal"] && item["orientation_equal"] for item in summary) ||
        error("MATLAB/Julia vg2g mismatch; inspect $summary_path")
    println("parity summary: ",abspath(summary_path))
    summary_path
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    workdir=isempty(ARGS) ? mktempdir() : ARGS[1]
    compare_vg2g(workdir=workdir,show_debug=get(ENV,"VG2G_SHOW_DEBUG","0")=="1")
end
