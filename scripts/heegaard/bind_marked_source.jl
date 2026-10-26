# Attach a full source to an existing view after checking all saved incidence.
using Topology, JSON3, SHA
length(ARGS)==2 || error("usage: bind_marked_source.jl VIEW.json SOURCE.json")
viewpath,sourcepath=ARGS
view=JSON3.read(read(viewpath,String),Dict{String,Any})
bytes2hex(sha256(read(sourcepath)))==view["source_sha256"] || error("Source file hash mismatch")
view["source_document"]=JSON3.read(read(sourcepath,String),Dict{String,Any})
Topology.Heegaard.marked_display_source(view)
write(viewpath,JSON3.write(view))
println("Full source bound; curve words, crossings, cut pairing and oriented D boundary checked; geometry unchanged")
