"""One app payload for registered examples and the import/render endpoint."""
function surface_payload(object,layout;id="imported",title="Imported Heegaard diagram",
                         source="User-supplied diagram",seed="Common Heegaard plot")
    diagram=H._base_diagram(object)
    chart=H.surface_chart(object,layout)
    Dict("id"=>id,"title"=>title,"genus"=>chart.genus,
        "intersections"=>length(chart.intersections),"beta"=>length(chart.beta_cycles),
        "source"=>source,"seed"=>seed,"svg"=>H.render(chart),
        "base_svg"=>H.render(layout),"chart"=>H._surface_chart_json(chart),
        "diagram"=>H._json_dict(object),"layout"=>H._layout_json(layout),
        "joining_disk_note"=>get(diagram.code.metadata,"joining_disk_note",""),
        "source_status"=>"software checks only; visual review by owner pending",
        "marking"=>object isa H.MarkedHeegaardDiagram ? Dict(
            "status"=>"verified-witness","hash"=>H.marked_hash(object),
            "frame"=>object.marking.frame,"routes"=>H._marking_json(object.marking)["routes"]) :
            Dict("status"=>"unmarked"))
end
