using Topology

"""Port the diagram-construction cells for p-twisted unknot complements."""
function run_C260319_unknots(max_twist=3)
    knot=VirtualLink(); set_data!(knot;rgauss=[[1,-1]],orientation=[1])
    knot_complement!(knot)
    base_code,base_orientation=real_gauss_code(knot)
    links=VirtualLink[]
    for twist in 1:max_twist
        first_component=reduce(vcat,([sign(x)*(abs(x)+4offset) for x in base_code[1]]
                                     for offset in 0:twist-1);init=Int[])
        second_component=reduce(vcat,([sign(x)*(abs(x)+4offset) for x in base_code[2]]
                                      for offset in twist-1:-1:0);init=Int[])
        link=VirtualLink()
        set_data!(link;rgauss=[first_component,second_component],
                  orientation=repeat(base_orientation,twist))
        push!(links,link)
    end
    links
end

"""Evaluate both Hopf algebras used by `C260319MSTinv_unknot.mlx`."""
function run_C260319_invariants(max_twist=3)
    links=run_C260319_unknots(max_twist)
    cyclic=cyclic_group_algebra(3)
    kac_paljutkin=kac_paljutkin_algebra()
    (
        links=links,
        cyclic=[hopf_invariant(link,cyclic) for link in links],
        kac_paljutkin=[hopf_invariant(link,kac_paljutkin) for link in links],
        cyclic_checks=verify_hopf(cyclic),
        kac_paljutkin_checks=verify_hopf(kac_paljutkin),
    )
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    display(run_C260319_invariants())
end
