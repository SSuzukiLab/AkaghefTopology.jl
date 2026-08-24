using Topology
using JSON

"""Decode one MATLAB `REdgeTable.Position` polyline exported as N-by-2 rows."""
function _matlab_position_polyline(rows)
    rows isa AbstractVector || throw(ArgumentError("MATLAB Position must be an array of [real, imag] rows"))
    points=ComplexF64[]
    for row in rows
        row isa AbstractVector && length(row)==2 || throw(ArgumentError("every MATLAB Position row must contain [real, imag]"))
        all(value -> value isa Number && isfinite(value),row) || throw(ArgumentError("MATLAB Position contains a non-finite coordinate"))
        push!(points,complex(Float64(row[1]),Float64(row[2])))
    end
    length(points)>=2 || throw(ArgumentError("MATLAB Position polyline requires at least two points"))
    points
end

function _c251020_fixture_link(call_id)
    result=run_C251020_knot_dsl2b()
    call_id==1 && return result.borromean,"borromean"
    call_id==2 && return result.whitehead,"whitehead"
    throw(ArgumentError("C251020 replay only defines plot calls 1 and 2"))
end

function _fixture_crossings(value)
    value isa AbstractVector || throw(ArgumentError("MATLAB REdgeTable.Crossing must be an array"))
    [begin
        row isa AbstractVector && length(row)==2 || throw(ArgumentError("every MATLAB Crossing row must contain two signed crossings"))
        Tuple(Int.(row))
     end for row in value]
end

"""Load the two C251020 MATLAB geometry fixtures from an exported replay JSON.

The loader refuses mismatched arc order or absent/non-finite coordinates so a
topologically compatible but geometrically unrelated replay cannot be drawn.
"""
function load_c251020_replay_fixtures(path::AbstractString)
    records=JSON.parsefile(path)
    records isa AbstractVector || throw(ArgumentError("replay JSON root must be an array"))
    fixtures=Dict{Int,NamedTuple}()
    for record in records
        relative=String(get(record,"relativePath",""))
        occursin("C251020KnotDsl2bfromUR/call_",relative) || continue
        call_id=Int(record["callId"])
        link,name=_c251020_fixture_link(call_id)
        properties=record["objects"]["vl"]["properties"]
        redge=properties["REdgeTable"]
        haskey(redge,"Position") || throw(ArgumentError("C251020 call $call_id lacks REdgeTable.Position"))
        expected=[edge.crossing for edge in real_edge_table(link)]
        actual=_fixture_crossings(redge["Crossing"])
        actual==expected || throw(ArgumentError("C251020 call $call_id REdgeTable traversal differs from Julia link"))
        positions=[_matlab_position_polyline(rows) for rows in redge["Position"]]
        length(positions)==length(expected) || throw(ArgumentError("C251020 call $call_id Position count differs from REdgeTable"))
        fixtures[call_id]=(name=name,link=link,positions=positions,relative_path=relative)
    end
    Set(keys(fixtures))==Set([1,2]) || throw(ArgumentError("replay must contain exactly C251020 calls 1 and 2"))
    fixtures
end

"""Render C251020 exactly along canonical MATLAB position polylines."""
function render_c251020_replay_fixtures(json_path::AbstractString,output_root::AbstractString)
    fixtures=load_c251020_replay_fixtures(json_path)
    mkpath(output_root)
    Dict(call_id => plot_svg(fixture.link,joinpath(output_root,
        "figure_$(lpad(call_id,3,'0'))_$(fixture.name).svg");replay_positions=fixture.positions)
         for (call_id,fixture) in fixtures)
end
