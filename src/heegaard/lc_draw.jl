"""
    draw_heegaard(lc::LickorishCode, svg_path::AbstractString; roles=Dict(), png=false, kwargs...) -> String

Write the planar SVG of a resolved Heegaard diagram.

# Arguments
The input must satisfy the [drawing domain](../specs/heegaard-plot.md).
Drawing validation errors propagate. `svg_path` is a writable destination;
parent directories are created and an existing file is overwritten.
Checked: `png=true` raises `ArgumentError` explaining that PNG conversion
is external. Other keywords are passed to `plot`.
`roles` supplies the explicit HC roles required by `resolve_hc`; resolution
errors occur before drawing. No role is inferred from the destination.

# Returns
The absolute output path as a `String`. The file contains the SVG returned
by the planar drawing entry.

# Guarantees
1. A successful return writes a drawing whose geometry audit reports success.
   [proved/agent: by construction]
2. The input is not modified; the destination file and missing parent
   directories are written. [proved/agent: by construction]

# Not guaranteed
This entry does not guarantee successful placement of every valid diagram,
produce PNG, or certify a 3-manifold invariant. Drawing and filesystem errors
propagate; output replacement is not transactional.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-plot.md — planar chart and geometry acceptance.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> lc = H.LickorishCode(1, Dict("l1"=>[-1,-2], "m1"=>[1,3], "c1"=>[-4,-3], "c2"=>[4,2]), Dict(), Dict(), Dict(), Dict(), nothing); hc = H.resolve_hc(lc, Dict("c1"=>:alpha, "c2"=>:beta));

julia> try H.draw_heegaard(hc, "unused.svg"; png=true) catch e; e isa ArgumentError end
true
```

Role: draw
"""
function draw_heegaard(lc::LickorishCode, svg_path::AbstractString;
        roles::AbstractDict{<:AbstractString,<:Any}=Dict{String,Symbol}(),
        png::Bool=false, kwargs...)
    draw_heegaard(resolve_hc(lc,roles),svg_path;png,kwargs...)
end

"""
    draw_heegaard(hc::HeegaardCode, svg_path::AbstractString; png=false, kwargs...) -> String

Write the planar SVG of a resolved Heegaard diagram.

# Arguments
`hc` satisfies the [drawing domain](../specs/heegaard-plot.md).
`svg_path` must be writable. Checked: `png=true` raises `ArgumentError`
explaining that PNG conversion is external. Other keywords pass to `plot`;
validation and drawing errors propagate.

# Returns
The absolute output path as a `String`, containing the planar entry's SVG.

# Guarantees
1. On success, the drawing's geometry audit reports success.
   [proved/agent: by construction]
2. The input is not modified; the destination is overwritten and missing
   parent directories are created. [proved/agent: by construction]

# Not guaranteed
Successful placement of every valid diagram, PNG output and transactional
file replacement are not guaranteed.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-plot.md — planar chart and geometry acceptance.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> hc = H.resolve_hc(H.LickorishCode(1, Dict("l1"=>[-1],"m1"=>[1]), Dict(), Dict(), Dict(), Dict(), nothing), Dict("l1"=>:beta,"m1"=>:alpha));

julia> try H.draw_heegaard(hc, "unused.svg"; png=true) catch e; e isa ArgumentError end
true
```

Role: draw
"""
function draw_heegaard(hc::HeegaardCode, svg_path::AbstractString; png::Bool=false, kwargs...)
    png && throw(ArgumentError("the Julia drawing entry outputs SVG; PNG conversion is external"))
    result=plot(hc;kwargs...)
    result.success || throw(ArgumentError("drawing did not make a verified nonzero smoothing step"))
    path=abspath(svg_path);mkpath(dirname(path));write(path,result.svg)
    path
end
