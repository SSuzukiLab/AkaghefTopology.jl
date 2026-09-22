using Topology
using JSON3
const H=Topology.Heegaard
length(ARGS)==2 || error("usage: render_request.jl REQUEST.json RESPONSE.json")
request=JSON3.read(read(ARGS[1],String))
mktempdir() do dir
    corepath=joinpath(dir,"input.hgd.json")
    write(corepath,JSON3.write(request.diagram))
    diagram=H.certify_heegaard(H.load_heegaard(corepath))
    before=if haskey(request,:layout) && request.layout !== nothing
        path=joinpath(dir,"input.hgd.layout.json")
        write(path,JSON3.write(request.layout)); H.load_layout(path,diagram)
    else
        H.layout_heegaard(diagram;backend=:pure_julia)
    end
    cut=H.cut_along_alpha(diagram)
    expected=H.default_puncture_placement(cut)
    all(any(p.id==q.id && p.center==q.center && p.radius==q.radius for q in expected)
        for p in before.punctures) || error("Please supply the standard 2×g paired-puncture layout.")
    after=H.plot(diagram,before)
    result=Dict("id"=>"imported","title"=>"Imported Heegaard diagram",
        "genus"=>length(diagram.code.code.alpha),"intersections"=>length(diagram.code.code.intersections),
        "beta"=>length(diagram.code.code.beta),"source"=>"User-supplied diagram",
        "seed"=>haskey(request,:layout) && request.layout !== nothing ? "User-supplied feasible layout" : "Bounded fixed-endpoint router",
        "diagram"=>request.diagram,"layout"=>H._layout_json(after),"before"=>H._layout_json(before),
        "svg"=>H.render(after;backend=:svg),"before_svg"=>H.render(before;backend=:svg))
    write(ARGS[2],JSON3.write(result))
end

