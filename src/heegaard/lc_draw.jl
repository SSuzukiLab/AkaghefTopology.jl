"""Write the SVG produced by the Julia `plot` entry.

cf: docs/specs/heegaard-plot.md — explicit HC roles, one placement/smoothing implementation
"""
function draw_heegaard(lc::LickorishCode, svg_path::AbstractString;
        roles::AbstractDict{<:AbstractString,<:Any}=Dict{String,Symbol}(),
        png::Bool=false, kwargs...)
    draw_heegaard(resolve_hc(lc,roles),svg_path;png,kwargs...)
end

function draw_heegaard(hc::HeegaardCode, svg_path::AbstractString; png::Bool=false, kwargs...)
    png && throw(ArgumentError("the Julia drawing entry outputs SVG; PNG conversion is external"))
    result=plot(hc;kwargs...)
    result.success || throw(ArgumentError("drawing did not make a verified nonzero smoothing step"))
    path=abspath(svg_path);mkpath(dirname(path));write(path,result.svg)
    path
end
