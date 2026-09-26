module KodaComplexTests
using Test, Topology
const H=Topology.Heegaard
include("../../examples/HeegaardAbalone.jl")

@testset "Koda sectors and independent hand-computed word" begin
    ab=HeegaardAbaloneExample.fixture()
    a=H.to_odata(ab.diagram;joining_disk=ab.disk)
    @test (length(a.spine.vertex_signs),length(a.spine.edges),length(a.spine.complex.regions))==(1,2,2)
    @test a.ograph.metadata["o_block_count"]==1
    ph=H.load_punctured(joinpath(@__DIR__,"../../artifacts/heegaard/full_curves/PoincareD.punctured.hgd.json"))
    b=H.to_odata(ph)
    @test (length(b.spine.vertex_signs),length(b.spine.edges),length(b.spine.complex.regions))==(8,16,9)
    @test b.ograph.metadata["o_block_count"]==2
    # Independent transcript, not generated from the converter's successors.
    labels=["x3","x2","x1","x7","x5","y5","y3","y2"]
    hand=[-1,-2,-3,-4,-5,6,4,7,-7,-8,-6,5,8,3,2,1]
    types=[:r,:r,:r,:l,:l,:l,:r,:l]
    source=Dict(p.id=>p.source_label for p in ph.diagram.code.code.intersections)
    original=sort(collect(keys(b.spine.vertex_signs)))
    rename=[only(findall(==(source[p]),labels)) for p in original]
    actual=[sign(v)*rename[abs(v)] for v in only(b.ograph.gauss)]
    @test H._cyclic_equal(actual,hand)
    # Transcript epsilon_K=det(oriented alpha,beta); repo sigma=det(beta,alpha).
    @test [b.ograph.orientation[only(findall(==(i),rename))] for i in 1:8]==
        [t==:r ? -1 : 1 for t in types]
    @test sort([r["source_intersections"] for r in b.ograph.metadata["runs"]];by=repr)==
        sort([[3,2,1,7,6],[9,7,11],[11,10,9],[6,10,1,2,3]];by=repr)
    for bundle in (a,b)
        c=bundle.spine.complex
        @test c.checks["euler_characteristic"]==1
        @test all(length(inc)==3 && sort(last.(inc))==[-1,1,1] for inc in c.edge_incidence)
        @test length(c.gluings)==length(bundle.punctured.joining_disk.boundary_darts)
        @test sort([e["spine_edge"] for e in bundle.ograph.metadata["edge_map"]])==collect(eachindex(bundle.spine.edges))
        # Replay every region's attaching map against its retained source darts.
        twin(d)=isodd(d) ? d+1 : d-1
        for r in c.regions
            expanded=vcat([e>0 ? bundle.spine.edges[e].source_darts :
                twin.(reverse(bundle.spine.edges[-e].source_darts)) for e in r.boundary]...)
            @test expanded==r.source_darts
        end
        link=VirtualLink(); set_data!(link;rgauss=bundle.ograph.gauss,orientation=bundle.ograph.orientation)
        @test is_closed(link)==(true,(true,true,true))
    end
end
end
