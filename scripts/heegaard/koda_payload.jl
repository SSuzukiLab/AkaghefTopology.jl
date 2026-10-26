function koda_payload(object)
    bundle=H.to_odata(object)
    data=H._conversion_bundle_json(bundle)
    link=VirtualLink()
    set_data!(link;rgauss=bundle.ograph.gauss,orientation=bundle.ograph.orientation)
    closed,criteria=is_closed(link)
    data["closedness"]=Dict("closed"=>closed,"criteria"=>collect(criteria))
    data["svg"]=mktempdir() do dir
        path=joinpath(dir,"ograph.svg")
        plot_svg(link,path;edge_labels=:none)
        read(path,String)
    end
    data
end
