import JSON3

"""A genuine closed sphere retaining the native nested disk as one hemisphere.

Source: vendor/marked_surface_diagrams/examples/nested_zero_crossings.json;
P and its curves are copied unchanged. Q is an oriented fan disk on the same
24-point boundary. Opposite cyclic side correspondence plus reversal gives
the sphere gluing, not a relabelling of the original physical-boundary disk.
Used by drawing and marked-operation tests; this is a small shared fixture.
"""
function capped_nested_sphere()
    path=joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples/nested_zero_crossings.json")
    doc=JSON3.read(read(path,String),Dict{String,Any})
    source=only(doc["patches"])
    boundary=String[v for side in source["side_order"] for v in source["sides"][side][1:end-1]]
    cap=Dict("face"=>"Q","side_order"=>deepcopy(source["side_order"]),"sides"=>deepcopy(source["sides"]),
        "triangles"=>[["cap_center",boundary[i],boundary[mod1(i+1,length(boundary))]] for i in eachindex(boundary)])
    push!(doc["patches"],cap)
    doc["surface"]=Dict("marking_id"=>"nested-sphere-v1","genus"=>0,"boundary"=>[],
        "gluings"=>[[["P","south"],["Q","south"]],[["P","east"],["Q","west"]],
                     [["P","north"],["Q","north"]],[["P","west"],["Q","east"]]])
    doc
end
