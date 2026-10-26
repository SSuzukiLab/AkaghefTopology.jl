"""One app payload for registered examples and the import/render endpoint."""
function surface_payload(object,layout;id="imported",title="Imported Heegaard diagram",
                         source="User-supplied diagram",seed="Common Heegaard plot",view=:alpha_cut)
    diagram=H._base_diagram(object)
    chart=H.surface_chart(object,layout;view)
    Dict("id"=>id,"title"=>title,"genus"=>chart.genus,
        "intersections"=>length(chart.intersections),"beta"=>length(chart.beta_cycles),
        "source"=>source,"seed"=>seed,"view"=>String(view),"svg"=>H.render(chart;view),
        "alpha_cut_svg"=>H.render(chart;view=:alpha_cut),"surface_svg"=>H.render(chart;view=:surface),
        "base_svg"=>H.render(layout),"chart"=>H._surface_chart_json(chart),
        "diagram"=>H._json_dict(object),"layout"=>H._layout_json(layout),
        "joining_disk_note"=>get(diagram.code.metadata,"joining_disk_note",""),
        "source_status"=>"software checks only; visual review by owner pending")
end
