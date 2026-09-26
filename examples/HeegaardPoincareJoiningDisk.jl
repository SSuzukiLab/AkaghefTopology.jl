"""The Poincare joining disk requested in the HeegaardDiagramProgram conversation.

No diagram-specific coordinates: retain the baseline labelled incidences and
use its validated layout. Source transcription remains owner-reviewable.
"""
module HeegaardPoincareJoiningDiskExample
using Topology
const H=Topology.Heegaard
const LABELS=["x1","x2","x3","x4","x6","x5","x7","y1","y5","y2","y3","y4"]
const SOURCE_CONVERSATION="chatgpt-conversation://6a8ece48-b94c-83ee-9018-e1181ed9134b"

function fixture()
    path=joinpath(@__DIR__,"..","artifacts","heegaard","app","Poincare.hgd.json")
    baseline=H.certify_heegaard(H.load_heegaard(path))
    source=baseline.code.code
    # Compare all four words and every signed incidence with the referenced
    # conversation, not merely the intersection count or homology.
    alpha=[["x1","x2","x3","x4","x6","x5","x7"],["y1","y5","y2","y3","y4"]]
    beta=[["x1","x2","x3","x4","y1","x5","y2"],["x7","y3","y4","x6","y5"]]
    @assert [LABELS[c.visits] for c in source.alpha]==alpha
    @assert [LABELS[c.visits] for c in source.beta]==beta
    signs=Dict(zip(vcat(["x$i" for i in 1:7],["y$i" for i in 1:5]),
                   [-1,-1,-1,-1,1,1,1,1,1,-1,-1,1]))
    @assert all(p.sign==signs[LABELS[p.id]] for p in source.intersections)
    metadata=merge(deepcopy(baseline.code.metadata),Dict{String,Any}(
        "reference_conversation"=>SOURCE_CONVERSATION,
        "reference_message"=>"7c91383a-4724-460a-b098-3107bc03d47d",
        "display_source_labels"=>true,
        "baseline_label_map"=>Dict(string(p.id)=>Dict("baseline"=>p.source_label,"conversation"=>LABELS[p.id]) for p in source.intersections),
        "joining_disk_note"=>"参照会話のD₀：x4 → α1 → x6 → β2⁻¹ → y4 → α2 → y1 → β1⁻¹ → x4。下の境界表示は、同じ四辺をアプリの『Dを左に見る』向きで逆にたどります。",
        "transcription_status"=>"matches referenced conversation after explicit relabelling; original-source review not upgraded"))
    code=H.heegaard_code(source.alpha,source.beta,
        [H.Intersection(p.id,p.alpha,p.beta,p.sign;source_label=LABELS[p.id]) for p in source.intersections];metadata)
    diagram=H.certify_heegaard(H.decorate(code,baseline.code.regions.boundary_cycle_blocks))
    graph=H.ribbon_graph(diagram)
    # These are the reversed oriented source edges, not just an equal corner set.
    expected=[(4,:beta_out,8),(8,:alpha_in,12),(12,:beta_out,5),(5,:alpha_in,4)]
    target=Set(expected)
    region=only([r for (r,block) in enumerate(diagram.code.regions.boundary_cycle_blocks)
        if length(block)==1 && Set((graph.darts[d].origin,graph.darts[d].origin_port,graph.darts[d].target)
            for d in H.boundary_cycles(graph)[only(block)])==target])
    disk=H.joining_disk(diagram,region)
    offset=findfirst(d->graph.darts[d].origin==4,disk.boundary_darts)
    ordered=circshift(disk.boundary_darts,1-offset)
    disk=H.JoiningDisk(disk.labelled_diagram_hash,disk.region,ordered)
    @assert [(graph.darts[d].origin,graph.darts[d].origin_port,graph.darts[d].target) for d in ordered]==expected
    H.make_punctured(diagram;joining_disk=disk)
end
end
