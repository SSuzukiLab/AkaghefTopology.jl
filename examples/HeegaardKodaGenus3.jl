"""Owner-supplied genus-3 marked Heegaard diagram, through the shared library API.

The original ZIP's o-data uses slope traversal (D on right). The repository uses
the induced boundary direction (D on left): reverse the complete Gauss word,
not individual visits or vertex types. Retain the full marked input as source.
"""
module HeegaardKodaGenus3Example
using Topology, JSON3
const H=Topology.Heegaard
function fixture()
    path=joinpath(@__DIR__,"../vendor/koda_genus3/data/figure6.heegaard.marked.json")
    doc=JSON3.read(read(path,String),Dict{String,Any})
    p=H.cellulated_punctured(doc)
    (document=doc,object=p.object,projection=p.provenance)
end
function odata()
    H.to_odata(fixture().object).ograph
end
end
