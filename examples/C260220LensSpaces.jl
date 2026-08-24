using Topology

"""Real Gauss presentation used for the lens-space sweep in C260220LensSpaces.mlx."""
function lens_space_link(s::Integer,t::Integer)
    s>=2 || throw(ArgumentError("s must be at least 2"))
    gcd(s,t)==1 || throw(ArgumentError("s and t must be coprime"))
    code=vcat(-mod.((s:-1:1).*t.-1,s).-1,1:s)
    link=VirtualLink(); set_data!(link;rgauss=[code],orientation=ones(Int,s))
    link
end

"""Reproduce the coprime sweep and disk-edge sums in C260220LensSpaces.mlx."""
function run_C260220_lens_spaces(;maximum_s=10)
    maximum_s>=2 || throw(ArgumentError("maximum_s must be at least 2"))
    sweep=NamedTuple[]
    for s in 2:maximum_s, t in 1:s-1
        gcd(s,t)==1 || continue
        link=lens_space_link(s,t)
        disks=disk_table(link)
        push!(sweep,(s=s,t=t,path_lengths=disks.path_lengths,
                     delta=disks.delta,cp2=disks.cp2))
    end
    lens21=lens_space_link(2,1); base,kernel,hperp=calculate_weight(lens21;assign=true)
    edges=real_edge_table(lens21)
    edge_weights=[edge.weight for edge in edges]
    disk_sums=disk_table(lens21).delta*edge_weights
    (sweep=sweep,lens21=lens21,weight=(base=base,kernel=kernel,hperp=hperp),
     disk_edge_sums=disk_sums)
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C260220_lens_spaces())
