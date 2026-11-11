"""Draw the cut surface of a Lickorish configuration or resolved HC.

cf: lib/topology.jl/scripts/heegaard/lc_diagram_draw.py draw_diagram — `.lc` テキストを渡し、配置器の SVG を書く
"""
function draw_heegaard(lc::LickorishConfiguration, svg_path::AbstractString;
        roles::AbstractDict{<:AbstractString,<:Any}=Dict{String,Symbol}(),
        png::Bool=false)
    svg_path = abspath(svg_path)
    mkpath(dirname(svg_path))
    script = normpath(joinpath(@__DIR__, "..", "..", "scripts", "heegaard", "lc_diagram_draw.py"))
    isfile(script) || throw(ArgumentError("placement entry is missing: $script"))
    mktempdir() do dir
        lc_path = joinpath(dir, "diagram.lc")
        write_lc(lc_path, lc)
        roles_path = joinpath(dir, "roles.json")
        open(roles_path, "w") do io
            write(io, JSON3.write(Dict(string(key) => string(value) for (key, value) in pairs(roles))))
        end
        command = `python $script $lc_path $svg_path --roles $roles_path`
        png && (command = `$command --png`)
        run(command)
    end
    svg_path
end

function draw_heegaard(hc::HeegaardLayer, svg_path::AbstractString; png::Bool=false)
    draw_heegaard(hc.lc, svg_path; roles=hc.roles, png=png)
end
