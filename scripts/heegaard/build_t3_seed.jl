using Topology, LinearAlgebra
const H=Topology.Heegaard
BLAS.set_num_threads(1)
include("t3_twist_seed.jl")
root=normpath(joinpath(@__DIR__,"..",".."))
diagram=H.certify_heegaard(H.load_heegaard(joinpath(root,"test","fixtures","t3_octahedron_v1.hgd.json")))
@time seed=octahedron_seed(diagram)
H.save_layout(joinpath(root,"test","fixtures","t3_octahedron_v1.hgd.layout.json"),seed)
H.render(seed;backend=:svg,path=joinpath(root,"artifacts","heegaard","app","T3-before.svg"))
println("T3 seed ready")
