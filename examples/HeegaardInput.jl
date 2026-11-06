#= Short marked authoring, correction and lossless JSON save/reload.

Source: vendor/marked_surface_diagrams/src/builders.jl::grid_document/add_curve!/
declare_crossing!, schema.jl::save_heegaard/load_heegaard. Grid indices denote
explicit triangulation edges, never drawing coordinates or inferred routes.
Run this existing Topology package example; the returned path can be opened in
the app's Data Open dialog. Display layout is a separate derived artifact.
=#
using Topology
const H = Topology.Heegaard
const B = H.MarkedDiagramBuilders

function heegaard_input_example(path::AbstractString)
    doc = B.grid_document(4; surface=:torus)
    B.add_curve!(doc, "a", "alpha", true, [[(0, 1), (4, 1)]])
    B.add_curve!(doc, "b", "beta", true, [[(2, 0), (2, 4)]])
    B.declare_crossing!(doc, "x", 2, 1, "a", "b", 1)
    doc["metadata"] = Dict("authoring"=>"explicit grid edges", "review"=>"agent example")
    H.save_heegaard(path, doc)
    corrected = H.load_heegaard(path)
    # Correct the explicit beta route from column 2 to column 3, together with
    # its crossing incidence. This is a source edit, not a certified isotopy.
    corrected["curves"][2]["word"] = B.waypoint_word([[(3, 0), (3, 4)]])
    corrected["crossings"][1]["at"] = ["P", "v3_1"]
    corrected["metadata"]["review"] = "ready for visual comparison"
    H.save_heegaard(path, corrected)
    reloaded = H.load_heegaard(path)
    @assert reloaded == corrected
    reloaded
end
