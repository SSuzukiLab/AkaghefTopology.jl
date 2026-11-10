# Bounded observation of the existing router. No router or core-data changes.
include(joinpath(@__DIR__,"..","..","src","heegaard","Heegaard.jl"))
using .Heegaard, JSON3, SHA
const H=Heegaard
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUT=joinpath(ROOT,"artifacts","heegaard","full_curves","literature")
const INPUT=joinpath(OUT,"klebanov-hempel-candidate.hgd.json")
mkpath(OUT)

function main()
    d=H.certify_heegaard(H.load_heegaard(INPUT))
    hash=H.canonical_hash(d)
    cut=H.cut_along_alpha(d;diagram_hash=hash)
    disks=H.default_puncture_placement(cut)
    ports=H.place_ports(cut,disks)
    wiring=H.build_wiring_problem(d,cut,ports)
    terminal=H._terminal_ribbon_report(d,cut,ports,wiring)
    println("Hempel input hash=",hash," n=",length(d.code.code.intersections),
        " g=",H.surface_genus(d)," terminal planar=",terminal.planar,
        " faces=",terminal.faces)
    flush(stdout)
    summary=Dict{String,Any}("input"=>basename(INPUT),"input_sha256"=>bytes2hex(sha256(read(INPUT))),
        "diagram_hash"=>hash,"intersections"=>length(d.code.code.intersections),
        "genus"=>H.surface_genus(d),"terminal"=>terminal,"core_valid"=>H.validate(d).valid,
        "attempts"=>Any[],"accepted"=>false,
        "scope"=>"existing greedy finite visibility graph; no backtracking; unchanged core and fixed 2-by-g disks",
        "source_status"=>"source correspondence inherited from the literature task; region decoration remains candidate")
    save()=open(io->JSON3.write(io,summary),joinpath(OUT,"hempel-routing-report.json"),"w")
    save()
    terminal.planar || return
    segmentids=Dict{Int,Int}(); perbeta=Dict{Int,Int}()
    for (k,w) in enumerate(wiring.wires)
        segmentids[k]=get(perbeta,w.beta_id,0)+1
        perbeta[w.beta_id]=segmentids[k]
    end
    order0=H._router_wire_order(wiring,ports)
    start=time()
    # For 13 wires, the default order/weight sequence repeats after 52 attempts.
    configs=[(step=.14,margin=1.5,attempts=52),
             (step=.18,margin=1.5,attempts=16),
             (step=.18,margin=2.2,attempts=12)]
    for config in configs
        nodes=H._router_grid_points(disks,ports;grid_step=config.step,margin=config.margin)
        for attempt in 1:config.attempts
            if time()-start>600
                summary["stop_reason"]="600-second bounded investigation budget"
                summary["elapsed_seconds"]=time()-start; save(); return
            end
            t0=time()
            order=circshift(order0,mod(fld(attempt-1,2),length(order0)))
            iseven(attempt) && reverse!(order)
            bend=.01*((attempt-1)%4)
            routes=H.Route[]; failed=nothing
            for wireindex in order
                w=wiring.wires[wireindex]
                path=H._router_candidate(nodes,w,ports,disks,routes;
                    bend_weight=bend,grid_step=config.step,max_expansions=25000)
                if path===nothing
                    failed=Dict("wire_index"=>wireindex,"beta_id"=>w.beta_id,
                        "segment_id"=>segmentids[wireindex],"from_port"=>w.from_port,"to_port"=>w.to_port)
                    break
                end
                push!(routes,H.Route(:beta,w.beta_id,segmentids[wireindex],w.from_port,w.to_port,path))
            end
            row=Dict{String,Any}("grid_step"=>config.step,"margin"=>config.margin,"attempt"=>attempt,
                "wire_order"=>order,"bend_weight"=>bend,"grid_nodes"=>length(nodes),
                "routed_prefix"=>length(routes),"failed_wire"=>failed,"seconds"=>time()-t0)
            push!(summary["attempts"],row)
            println("step=",config.step," margin=",config.margin," attempt=",attempt,
                    " prefix=",length(routes),"/",length(wiring.wires)," seconds=",round(time()-t0;digits=2))
            flush(stdout)
            if failed===nothing
                trace=Dict{String,Any}("backend"=>"pure_julia","catlab_used"=>false,
                    "post_router"=>"bounded-fixed-endpoint-feasible",
                    "router_limits"=>config,"selected_attempt"=>attempt,
                    "provenance"=>Dict("input"=>"klebanov-hempel-candidate",
                        "region_decoration_status"=>"candidate_not_source_transcribed",
                        "source_correspondence"=>"inherited from literature task; not re-audited here"))
                provisional=H.HeegaardLayout(hash,:alpha_cut,Dict{String,Any}(
                    "chart"=>"sphere-at-infinity","outer_halo"=>true),disks,ports,routes,
                    H._route_crossings(routes),trace,H.LayoutValidationReport(false,H.ValidationIssue[],Dict{Symbol,Bool}()))
                report=H.validate_layout(d,provisional;expected_hash=hash)
                row["valid"]=report.valid
                row["issues"]=[string(i.code) for i in report.issues]
                if report.valid
                    seed=H.HeegaardLayout(hash,:alpha_cut,provisional.chart,disks,ports,routes,
                        provisional.crossings,trace,report)
                    H.save_layout(joinpath(OUT,"Hempel-seed.hgd.layout.json"),seed)
                    H.render(seed;path=joinpath(OUT,"Hempel-seed.svg"))
                    summary["accepted"]=true; summary["stop_reason"]="verified feasible seed"; save()
                    fixed=H.plot(d,seed)
                    H.save_layout(joinpath(OUT,"Hempel.hgd.layout.json"),fixed)
                    H.render(fixed;path=joinpath(OUT,"Hempel-alpha-cut.svg"))
                    atlas=H.surface_chart(d,fixed)
                    svg=H.render(atlas;path=joinpath(OUT,"Hempel-surface.svg"))
                    data=Dict("id"=>"Hempel_candidate","title"=>"Hempel · reconstruction candidate",
                        "genus"=>2,"intersections"=>13,"beta"=>2,
                        "source"=>"Literature-task transcription; singleton region reconstruction candidate",
                        "seed"=>"Bounded existing-router verified seed",
                        "svg"=>svg,"base_svg"=>H.render(fixed),"chart"=>H._surface_chart_json(atlas),
                        "diagram"=>H._json_dict(d))
                    open(io->JSON3.write(io,data),joinpath(dirname(OUT),"Hempel_candidate.json"),"w")
                    # The dedicated URL works without silently promoting it into the baseline seven.
                    summary["fixed100"]=true; summary["atlas"]=true
                    summary["page"]="http://localhost:8741/surface_Hempel_candidate.html"
                    summary["elapsed_seconds"]=time()-start; save(); return
                end
            end
            save()
        end
    end
    summary["stop_reason"]="all recorded finite candidates failed; not a nonexistence result"
    summary["elapsed_seconds"]=time()-start
    save()
end
main()
