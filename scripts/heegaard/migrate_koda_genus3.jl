# One-shot migration: supplied genus-3 marked cellulation -> ordinary punctured input.
# The written punctured JSON becomes the source for the unmarked 2D pipeline; this
# script is the only step that reads the cellulated document.
using Topology, JSON3, SHA
const H=Topology.Heegaard
root=normpath(joinpath(@__DIR__,"../.."))
src=joinpath(root,"vendor/koda_genus3/data/figure6.heegaard.marked.json")
doc=JSON3.read(read(src,String),Dict{String,Any})
projected=H.cellulated_punctured(doc)
ph=projected.object
code=ph.diagram.code
meta=Dict{String,Any}(
    "example"=>"Koda genus-3 figure with joining disk D",
    "surface_orientation"=>"dx wedge dy",
    "migrated_from"=>Dict("path"=>"vendor/koda_genus3/data/figure6.heegaard.marked.json",
        "file_sha256"=>bytes2hex(open(sha256,src)),
        "source_representation_sha256"=>projected.provenance["source_representation_sha256"],
        "disk_id"=>projected.provenance["disk_id"]),
    "equivalence"=>"labelled diagram up to orientation-preserving homeomorphism (S2 ambient isotopy of the alpha-cut picture); marking not retained")
h=H.certify_heegaard(H.decorate(code,code.regions.boundary_cycle_blocks;metadata=meta,cellular=true))
region=ph.joining_disk.region
out=H.make_punctured(h;joining_disk=H.joining_disk(h,region))
path=joinpath(root,"artifacts/heegaard/full_curves/KodaGenus3.punctured.hgd.json")
H.save_punctured(path,out)
back=H.load_punctured(path)
H._labelled_hash(back.diagram)==H._labelled_hash(h) || error("round trip changed the labelled diagram")
bundle=H.to_odata(back)
og=H._conversion_bundle_json(bundle)["ograph"]
got=JSON3.read(JSON3.write([og["gauss"],og["orientation"]]))
target=JSON3.read(read(joinpath(root,"artifacts/heegaard/cellulated/KodaGenus3.repository.odata.json"),String))
JSON3.write(got)==JSON3.write(target) || error("o-data differs from the verified cellulated conversion")
println("MIGRATED ",path," region=",region)
