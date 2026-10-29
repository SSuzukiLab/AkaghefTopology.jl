# Rebuild both examples from punctured combinatorial inputs, without a layout fixture.
using Topology, JSON3
const H=Topology.Heegaard
include("surface_payload.jl")
include("koda_payload.jl")
const OUT=normpath(joinpath(@__DIR__,"../../artifacts/heegaard/full_curves"))
for (base,path,title) in (("Abalone","abalone/S3.punctured.hgd.json","S³ · abalone + D"),
                          ("PoincareD","PoincareD.punctured.hgd.json","Poincaré + D₀"),
                          ("KodaGenus3","KodaGenus3.punctured.hgd.json","Koda · genus 3 + D"))
    targets=filter(x->!startswith(x,"--"),ARGS)
    !isempty(targets) && !(base in targets) && continue
    object=H.load_punctured(joinpath(OUT,path))
    @info "combinatorial input loaded" base
    koda=koda_payload(object)
    write(joinpath(OUT,"$base.koda.json"),JSON3.write(koda))
    if isfile(joinpath(OUT,"$base.json")) # hand-placed baseline, when one exists
        baseline=JSON3.read(read(joinpath(OUT,"$base.json"),String),Dict{String,Any})
        baseline["koda"]=koda
        write(joinpath(OUT,"$base.json"),JSON3.write(baseline))
    end
    "--conversion-only" in ARGS && continue
    before=H.layout_heegaard(object.diagram;backend=:pure_julia,router=:planar)
    H.save_layout(joinpath(OUT,"$base-auto-before.hgd.layout.json"),before)
    @info "fixture-free planar layout validated" base
    after=H.plot(object.diagram,before)
    id=base*"Auto"
    H.save_layout(joinpath(OUT,"$id.hgd.layout.json"),after)
    data=surface_payload(object,after;id,title=title*" · 自動配置",
        source="同じ保存済み数学データとD",
        seed="組合せ平面埋込み → 極座標ポート → 周囲同位による2×g配置 → 共通固定100描画（配置fixtureなし）")
    data["koda"]=koda
    write(joinpath(OUT,"$id.json"),JSON3.write(data))
    write(joinpath(OUT,"$id.svg"),data["svg"])
    entries=JSON3.read(read(joinpath(OUT,"manifest.json"),String),Vector{Dict{String,Any}})
    filter!(x->x["id"]!=id,entries)
    push!(entries,Dict(k=>data[k] for k in ("id","title","genus","intersections","beta","source","seed")))
    write(joinpath(OUT,"manifest.json"),JSON3.write(entries))
    println("GENERATED ",id)
end
