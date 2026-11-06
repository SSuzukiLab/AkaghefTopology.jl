module MarkedAnnulusProductTests
using Test,Topology
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const B=H.MarkedDiagramBuilders

# Source: established native torus meridian/longitude fixture. The selected
# support comes from the existing independently normalized regular-annulus
# constructor, not from an expected coordinate list made by this test.
@testset "Exact product chart on actual derived annulus flags" begin
    doc=B.grid_document(2)
    B.add_curve!(doc,"a","alpha",true,[[(0,1),(2,1)]])
    B.add_curve!(doc,"b","beta",true,[[(1,0),(1,2)]])
    B.declare_crossing!(doc,"x",1,1,"a","b",1)
    charts=Any[]
    for reverse_core in (false,true)
        input=deepcopy(doc)
        if reverse_core
            curve=only(c for c in input["curves"] if c["id"]=="a")
            ds=M.parse_word(curve["word"])
            curve["word"]=M.print_word([M.Dart(d.face,d.v,d.u) for d in reverse(ds)])
            for x in input["crossings"];"a" in x["curves"] && (x["sign"]=-x["sign"]);end
        end
        original=deepcopy(input)
        parallel=H.marked_parallel_copy(input,Dict("curve"=>"a","id"=>"parallel","side"=>"left"))
        coarse=H.marked_dual_star(input,Dict("merge_charts"=>true,"subdivision_budget"=>256)).document
        cst=H._marked_states(coarse).geometry;fst=H._marked_states(parallel.document).geometry
        selected=Set(Tuple(c) for c in parallel.witness["correspondence"]["annulus_triangles"])
        chart=H._marked_annulus_product(cst,fst,"a",selected)
        push!(charts,chart)
        @test input==original
        @test chart.period==1 && chart.support_cells==selected
        @test Set(Tuple(t.id) for t in chart.triangles)==selected
        @test all(H._ann_det(t.points)>0 for t in chart.triangles)
        @test chart.certificate["area"]=="2" && chart.certificate["core_winding"]=="1"
        @test Set(f["side"] for f in chart.certificate["frontiers"])==Set([-1,1])
        @test any(s["u_shift"]!="0" for s in chart.seam_lifts)
        @test_throws ArgumentError H._marked_annulus_product(cst,fst,"a",setdiff(selected,Set([first(selected)])))
        bad=deepcopy(fst);firstcell=first(sort(collect(keys(bad.triangles))))
        bad.triangles[firstcell]=reverse(bad.triangles[firstcell])
        @test_throws ArgumentError H._marked_annulus_product(cst,bad,"a",selected)
    end
    # Reversing the geometric core reverses both longitudinal and transverse
    # coordinates on the same actual mesh. This independently detects a
    # left/right or core-start convention error while preserving orientation.
    positive=Dict(Tuple(t.id)=>t.points for t in charts[1].triangles)
    negative=Dict(Tuple(t.id)=>t.points for t in charts[2].triangles)
    @test keys(positive)==keys(negative)
    @test all(negative[c][j][2]==-positive[c][j][2] for c in keys(positive) for j in 1:3)
    @test all(denominator(negative[c][j][1]+positive[c][j][1])==1 for c in keys(positive) for j in 1:3)
end
end
