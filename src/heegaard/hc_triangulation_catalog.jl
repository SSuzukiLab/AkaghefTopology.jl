# cf: lib/topology.jl/docs/specs/heegaard-snappy.md — full-input catalog protocol
# Standalone verification driver, deliberately not included by Heegaard.jl.
using Topology
const CAT_H = Topology.Heegaard
const CAT_JSON = CAT_H.JSON3

function catalog_error(e)
    result=Dict{String,Any}("status"=>"error", "type"=>string(typeof(e)), "message"=>sprint(showerror,e))
    for field in (:reason,:condition,:detail,:witness,:stage)
        hasproperty(e,field) && (result[String(field)]=getproperty(e,field))
    end
    result
end

function catalog_attempt(f)
    try
        Dict{String,Any}("status"=>"ok", "value"=>f())
    catch e
        catalog_error(e)
    end
end

function run_hc_triangulation_catalog(root; runner=nothing, backend_available=true,
                                      mode="package", output=joinpath(root,"threads","PJ260922","data","261111-snappy-catalog"))
    manifest = CAT_JSON.read(read(joinpath(output,"inputs.json"),String),Dict{String,Any})
    rows = Dict{String,Any}[]
    comparisons = Dict{String,Any}[]
    mkpath(joinpath(output,"triangulations"))
    function emit(row)
        push!(rows,row)
        open(joinpath(output,"rows.jsonl"),"a") do io
            println(io,CAT_JSON.write(row)); flush(io)
        end
    end
    # Refuse to overwrite a previous experiment; caller archives its own run explicitly.
    isfile(joinpath(output,"rows.jsonl")) && error("rows.jsonl already exists; use a fresh output directory or archive the prior run")
    for (index,input) in enumerate(manifest["inputs"])
        label = input["input"]
        base = Dict{String,Any}("input"=>label, "source"=>input["path"], "origin"=>input["origin"],
            "mode"=>mode, "joining_disk"=>nothing, "triangulation"=>nothing,
            "snappy"=>nothing, "h1_match"=>nothing, "hom_matches"=>Dict{String,Any}(g=>nothing for g in ("Z2","Z3","Z5","S3")))
        println("INPUT ",index,"/",length(manifest["inputs"])," ",label); flush(stdout)
        hc = try
            CAT_H.read_hc(joinpath(root,input["path"]))
        catch e
            base["conversion"] = catalog_error(e); emit(base); continue
        end
        base["hc_homology"] = catalog_attempt() do
            h=CAT_H.hc_first_homology(hc)
            Dict("free_rank"=>h.free_rank,"torsion"=>h.torsion,"condition"=>String(h.condition))
        end
        base["hc_presentation"] = catalog_attempt() do
            p=CAT_H.hc_presentation(hc).presentation
            Dict("generators"=>p.generators,"relators"=>p.relators,"metadata"=>p.metadata)
        end
        base["hc_hom_counts"] = Dict(String(g)=>catalog_attempt(()->string(CAT_H.hc_hom_count(hc,g)))
                                       for g in (:Z2,:Z3,:Z5,:S3))
        candidates = try
            forgot=CAT_H.forget_reference(CAT_H.realize_hc(hc))
            [r.id for r in forgot.regions if r.joining_disk]
        catch e
            base["conversion"] = catalog_error(e); emit(base); continue
        end
        base["joining_candidates"] = candidates
        if isempty(candidates)
            base["conversion"] = Dict("status"=>"not_computed","reason"=>"no_joining_disk_candidate")
            emit(base); continue
        end
        accepted = Pair{Int,CAT_H.HCClosedTriangulation}[]
        for disk in candidates
            row=deepcopy(base); row["joining_disk"]=disk
            tri = try
                CAT_H.hc_triangulation(hc; joining_disk=disk,signature=false)
            catch e
                row["conversion"]=catalog_error(e); emit(row); continue
            end
            row["conversion"]=Dict("status"=>"ok")
            row["triangulation"]=Dict(string(k)=>v for (k,v) in pairs(tri.topology))
            row["triangulation"]["tetrahedra"]=tri.tetrahedra
            filename="input-$(lpad(index,3,'0'))-disk-$disk.json"
            write(joinpath(output,"triangulations",filename),CAT_JSON.write(CAT_H.hc_triangulation_data(tri)))
            row["triangulation_file"]="triangulations/"*filename
            push!(accepted,disk=>tri)
            if backend_available
                row["snappy"] = catalog_attempt(() -> CAT_H.hc_snappy(tri; runner=runner))
                if row["snappy"]["status"]=="ok"
                    s=row["snappy"]["value"]
                    if get(s["homology"],"status","")=="ok" && base["hc_homology"]["status"]=="ok"
                        a,b=s["homology"]["value"],base["hc_homology"]["value"]
                        row["h1_match"]=a["free_rank"]==b["free_rank"] && sort(a["torsion"])==sort(b["torsion"])
                    end
                    for g in ("Z2","Z3","Z5","S3")
                        a,b=s["hom_counts"][g],base["hc_hom_counts"][g]
                        if a["status"]=="ok" && b["status"]=="ok"
                            row["hom_matches"][g]=a["value"]==b["value"]
                        end
                    end
                end
            else
                row["snappy"]=Dict("status"=>"not_computed","reason"=>"backend_unavailable")
            end
            emit(row)
        end
        for i in eachindex(candidates), j in (i+1):length(candidates)
            a,b=findfirst(p->p.first==candidates[i],accepted),findfirst(p->p.first==candidates[j],accepted)
            result = if a===nothing || b===nothing
                Dict("status"=>"not_computed","reason"=>"triangulation_unavailable")
            elseif !backend_available
                Dict("status"=>"not_computed","reason"=>"backend_unavailable")
            else
                catalog_attempt(()->CAT_H.hc_snappy_compare(accepted[a].second,accepted[b].second;runner=runner))
            end
            push!(comparisons,Dict("input"=>label,"left_disk"=>candidates[i],"right_disk"=>candidates[j],"result"=>result))
        end
    end
    write(joinpath(output,"joining-comparisons.json"),CAT_JSON.write(comparisons))
    mismatch = [r for r in rows if r["h1_match"]===false || any(v===false for v in values(r["hom_matches"]))]
    summary=Dict("inputs"=>length(manifest["inputs"]),"rows"=>length(rows),"mode"=>mode,
        "converted"=>count(r->r["conversion"]["status"]=="ok",rows),
        "mismatches"=>length(mismatch),"h1_comparisons_completed"=>count(r->r["h1_match"]!==nothing,rows),
        "joining_comparisons"=>length(comparisons),
        "snappy_completed"=>count(r->r["snappy"]!==nothing && get(r["snappy"],"status","")=="ok",rows))
    write(joinpath(output,"summary.json"),CAT_JSON.write(summary))
    write(joinpath(output,"mismatches.json"),CAT_JSON.write(mismatch))
    println(CAT_JSON.write(summary))
    summary
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    root=normpath(joinpath(@__DIR__,"..","..","..",".."))
    run_hc_triangulation_catalog(root)
end
