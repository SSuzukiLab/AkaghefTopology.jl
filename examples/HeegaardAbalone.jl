"""Small S3 Heegaard diagram on R²/Z², with one cancelling pair and a joining disk."""
module HeegaardAbaloneExample
using Topology
const H = Topology.Heegaard

function fixture()
    # Coordinates use dx ∧ dy. Opposite sides of the unit square are identified.
    # beta winds (0,1); the upper/lower lobes meet alpha only at their endpoints.
    beta = [
        ((.20,0.),(.20,.20),(.24,.35),(.30,.50)),
        ((.30,.50),(.36,.75),(.44,.75),(.50,.50)),
        ((.50,.50),(.56,.25),(.64,.25),(.70,.50)),
        ((.70,.50),(.80,.75),(.90,.78),(.90,.86)),
        ((.90,.86),(.90,.96),(.20,.84),(.20,1.)),
    ]
    points = [beta[i][4] for i in 1:3]
    signs = [Int(sign(beta[i][4][2]-beta[i][3][2])) for i in 1:3]
    alpha_order = sortperm(points;by=first)
    code = H.heegaard_code([alpha_order],[collect(1:3)],
        [H.Intersection(i,1,1,signs[i]) for i in 1:3];metadata=Dict(
            "example"=>"S3 with abalone joining disk",
            "construction"=>"torus meridian/longitude with one cancelling bigon pair",
            "surface_orientation"=>"dx wedge dy",
            "alpha_homology"=>[1,0],"beta_homology"=>[0,1]))
    diagram = H.certify_heegaard(H.decorate(code;cellular=true))
    graph = H.ribbon_graph(diagram)
    # Select by the actual boundary 1 --alpha--> 2 --beta reversed--> 1.
    region = only([r for (r,block) in enumerate(diagram.code.regions.boundary_cycle_blocks)
        if length(block)==1 && [(graph.darts[d].origin,graph.darts[d].origin_port)
            for d in H.boundary_cycles(graph)[only(block)]] == [(1,:alpha_out),(2,:beta_in)]])
    disk = H.joining_disk(diagram,region)
    (;diagram,disk,beta,points)
end

"""Exact example chart; topology is unchanged by the SVG coordinate transform."""
function render_svg(example)
    (;diagram,disk,beta,points) = example
    H.validate_joining_disk(diagram,disk)
    H.canonical_hash(diagram) == H.canonical_hash(fixture().diagram) ||
        throw(ArgumentError("this renderer is specific to the abalone smoke fixture"))
    xy(p) = "$(round(65+500p[1];digits=3)),$(round(575-500p[2];digits=3))"
    cubic(b) = "C $(xy(b[2])) $(xy(b[3])) $(xy(b[4]))"
    beta_path = "M $(xy(beta[1][1])) " * join(cubic.(beta)," ")
    upper = beta[2]
    disk_path = "M $(xy(points[1])) L $(xy(points[2])) C $(xy(upper[3])) $(xy(upper[2])) $(xy(upper[1])) Z"
    io=IOBuffer()
    print(io,"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 640 680" data-heegaard="abalone" data-source-hash="$(H.canonical_hash(diagram))">
    <style>text{font-family:Segoe UI,sans-serif;fill:#172033;font-size:15px}.label{font-size:19px;font-style:italic}.small{font-size:13px;fill:#526078}</style>
    <rect width="640" height="680" fill="white"/>
    <text x="65" y="32" font-size="21">S³ Heegaard diagram + joining disk D</text>
    <text x="65" y="55" class="small">T² = [0,1]² / opposite edges; surface orientation dx ∧ dy</text>
    <rect x="65" y="75" width="500" height="500" fill="#f8fafc" stroke="#94a3b8" stroke-width="1.5"/>
    <path data-joining-disk="$(disk.region)" data-boundary-darts="$(join(disk.boundary_darts,","))" d="$disk_path" fill="#cbd5e1" stroke="none"/>
    <path data-curve-family="alpha" d="M 65,325 L 565,325" fill="none" stroke="#1d4ed8" stroke-width="3"/>
    <path data-curve-family="beta" d="$beta_path" fill="none" stroke="#b45309" stroke-width="3"/>
    <path d="M 471,319 L 481,325 L 471,331" fill="none" stroke="#1d4ed8" stroke-width="2"/>
    <path d="M 159,514 L 165,504 L 171,514" fill="none" stroke="#b45309" stroke-width="2"/>
    <text x="500" y="314" class="label" style="fill:#1d4ed8">α</text>
    <text x="138" y="491" class="label" style="fill:#b45309">β</text>
    <text x="256" y="273" class="label">D</text>
    <text x="280" y="68">a →</text><text x="280" y="599">a →</text>
    <text x="31" y="190">b ↑</text><text x="575" y="190">b ↑</text>
    """)
    for (i,p) in enumerate(points)
        x,y = 65+500p[1],575-500p[2]
        print(io,"<circle data-intersection-id=\"$i\" cx=\"$x\" cy=\"$y\" r=\"4\" fill=\"#172033\"/>")
        print(io,"<text x=\"$(x-9)\" y=\"$(y+24)\">$i</text>")
    end
    print(io,"""<text x="65" y="626">α: (1 2 3)　β: (1 2 3)　intersection signs: (+ − +)</text>
    <text x="65" y="651" class="small">D boundary: 1 → 2 on α, 2 → 1 on β. No over/under on this surface.</text></svg>""")
    String(take!(io))
end
end
