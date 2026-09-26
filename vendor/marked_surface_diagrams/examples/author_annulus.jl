include(joinpath(@__DIR__,"..","src","MarkedSurfaceDiagrams.jl"))
include(joinpath(@__DIR__,"..","src","builders.jl"))
using .MarkedSurfaceDiagrams, .MarkedDiagramBuilders
using TOML
# Same genuine boundary endpoints; one positive passage across the seam.
d=grid_document(4;surface=:annulus)
add_curve!(d,"gamma","aux",false,[[(1,0),(4,3)],[(0,3),(1,4)]])
st=validate(d)               # Must succeed before saving or drawing.
@assert only(seam_word(st,"gamma"))["direction"]==1
out=joinpath(@__DIR__,"..","artifacts")
open(joinpath(out,"authored_annulus.toml"),"w") do io;TOML.print(io,d);end
write_svg(st,realize(st),joinpath(out,"authored_annulus.svg"))
