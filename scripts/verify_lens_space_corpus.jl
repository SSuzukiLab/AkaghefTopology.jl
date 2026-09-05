using Topology, JSON

function lens_space_corpus(;output_dir=mktempdir())
    mkpath(output_dir)
    records=Any[]
    for s in 2:10, t in 1:s-1
        gcd(s,t)==1 || continue
        name="L_$(s)_$(t)"
        input=lens_space_rgauss(s,t)
        status="ok"
        error_text=""
        crossing_count=0
        virtual_count=0
        region_count=0
        euler_characteristic=0
        genus=0
        diagnostics=0
        svg_path=joinpath(output_dir,"$(name).svg")
        try
            source=VirtualLink()
            set_data!(source;rgauss=input.rgauss,orientation=input.orientation)
            gauss,orientation=planarize_virtual_gauss(source.gauss,source.orientation;DEBUG_MODE=true,debug_io=devnull)
            labels=abs.(filter(!=(0),reduce(vcat,gauss;init=Int[])))
            real_count=length(input.orientation)
            preserved=[[entry for entry in component if entry!=0 && abs(entry)<=real_count]
                       for component in gauss]
            expected=[[entry for entry in component if entry!=0] for component in input.rgauss]
            preserved==expected || throw(ArgumentError("real crossing order changed"))
            all(label -> count(==(label),labels)==2,unique(labels)) ||
                throw(ArgumentError("crossing label does not occur exactly twice"))
            orientation[1:real_count]==input.orientation || throw(ArgumentError("real orientation changed"))
            all(==(0),orientation[real_count+1:end]) || throw(ArgumentError("inserted orientation is not zero"))

            planar=VirtualLink()
            set_data!(planar;gauss=gauss,orientation=orientation)
            pd=pd_code(planar)
            regions_data=regions(pd)
            crossing_count=size(pd,1)
            virtual_count=count(==(0),orientation)
            region_count=length(regions_data)
            euler_characteristic=crossing_count-2crossing_count+region_count
            iseven(round(Int,euler_characteristic)) || throw(ArgumentError("Euler characteristic is nonintegral"))
            genus=round(Int,1-euler_characteristic/2)
            genus==0 || throw(ArgumentError("Julia region genus is $genus"))

            layout=orthogonal_layout(source)
            diagnostics=length(check_layout(layout))
            diagnostics==0 || throw(ArgumentError("layout diagnostics=$diagnostics"))
            plot_svg(source,svg_path)
            isfile(svg_path) && filesize(svg_path)>0 || throw(ArgumentError("SVG was not written"))
        catch error
            status="error"
            error_text=sprint(showerror,error)
        end
        record=Dict(
            "name"=>name,"s"=>s,"t"=>t,"status"=>status,
            "crossings"=>crossing_count,"virtual_crossings"=>virtual_count,
            "regions"=>region_count,"euler_characteristic"=>euler_characteristic,
            "genus"=>genus,"layout_diagnostics"=>diagnostics,"error"=>error_text,
        )
        push!(records,record)
        println(name,": ",status,
                status=="ok" ? " crossings=$(crossing_count) virtual=$(virtual_count) chi=$(euler_characteristic)" :
                " error=$(error_text)")
    end
    summary=Dict("count"=>length(records),"ok"=>count(x->x["status"]=="ok",records),
                 "errors"=>count(x->x["status"]!="ok",records),"records"=>records,
                 "output_dir"=>abspath(output_dir))
    open(joinpath(output_dir,"lens-corpus.json"),"w") do io
        JSON.print(io,summary,2)
    end
    summary
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    output=isempty(ARGS) ? mktempdir() : ARGS[1]
    summary=lens_space_corpus(output_dir=output)
    println("summary=",summary["ok"],"/",summary["count"]," ok; output=",summary["output_dir"])
    summary["errors"]==0 || error("Lens corpus has failures")
end
