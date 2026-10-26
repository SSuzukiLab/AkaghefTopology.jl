using Topology, JSON3
length(ARGS)==2 || error("usage: analyze_marked.jl INPUT.json REPORT.json")
doc=JSON3.read(read(ARGS[1],String),Dict{String,Any})
report=Topology.Heegaard.marked_topology_report(doc)
write(ARGS[2],JSON3.write(report))
c=report["compression"]
println(JSON3.write((surface=report["surface"],alpha_boundary_genera=c.alpha_boundary_genera,
    beta_boundary_genera=c.beta_boundary_genera,closed=c.closed,
    joining_disks=report["joining_disks"])))
