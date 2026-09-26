using Topology
const M = Topology.Heegaard.MarkedSurfaceDiagrams
const B = Topology.Heegaard.MarkedDiagramBuilders
# Authoring grid indices, not coordinates or an implicit route search.
doc = B.grid_document(4; surface=:annulus)
B.add_curve!(doc,"gamma","aux",false,[[(1,0),(4,3)],[(0,3),(1,4)]])
state = M.validate(doc)
positions = M.realize(state)
out = joinpath(@__DIR__,"../artifacts/heegaard/cellulated")
mkpath(out)
M.write_json(joinpath(out,"authored-annulus.json"),state.doc)
M.write_svg(state,positions,joinpath(out,"authored-annulus.svg"))
println("Open authored-annulus.json in the existing app's データを開く dialog")
