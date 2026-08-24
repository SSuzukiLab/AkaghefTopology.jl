using Topology

"""Reproduce the numerical topology cells of `C250818MSTinv_uqsl2.mlx`."""
function run_C250818_uqsl2()
    algebra=uqsl2_borel_small(2//5;style=:L)

    lens=VirtualLink(); set_data!(lens;gauss=[[1,-3,2,-2,3,-1]],orientation=[1,1,0])
    lens_base,lens_kernel,_=calculate_weight(lens)
    lens_values=ComplexF64[]
    for weight in (lens_base,lens_base+3vec(lens_kernel[1,:]))
        set_weight!(lens,[weight]); push!(lens_values,hopf_invariant(lens,algebra))
    end

    abalone=VirtualLink(); set_data!(abalone;rgauss=[[1,-1]],orientation=[1])
    calculate_weight(abalone;assign=true)

    s2xs1=VirtualLink()
    set_data!(s2xs1;gauss=[[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],
              orientation=[1,1,1,1,0,0])
    s2_base,_,s2_hperp=calculate_weight(s2xs1)
    hperp=vec(s2_hperp)
    s2_values=ComplexF64[]
    for weight in (s2_base,hperp,s2_base+hperp,s2_base+2hperp,s2_base+3hperp)
        set_weight!(s2xs1,[weight]); push!(s2_values,hopf_invariant(s2xs1,algebra))
    end

    koda=VirtualLink()
    set_data!(koda;rgauss=[[2,-1,-2,1,4,-3,-4,3]],orientation=[-1,1,-1,1])
    calculate_weight(koda;assign=true)

    connected=deepcopy(abalone)
    disjoint_sum!(connected,abalone); connected_sum!(connected); calculate_weight(connected;assign=true)

    complement=deepcopy(abalone); knot_complement!(complement)
    set_weight_by_crossing!(complement,[2 -2;-3 1;3 4;4 -4],[1,-1,-1,1])
    complement2=deepcopy(complement)
    set_weight_by_crossing!(complement2,[2 -2;-3 1;3 4;4 -4;-2 -3;-4 -1],
                            [1,-1,-1,1,2,-2])
    (
        hopf_checks=verify_hopf(algebra),
        lens_values=lens_values,
        abalone_value=hopf_invariant(abalone,algebra),
        s2xs1_values=s2_values,
        koda_value=hopf_invariant(koda,algebra),
        connected_values=[hopf_invariant(abalone,algebra),hopf_invariant(connected,algebra)],
        knot_complement_values=[hopf_invariant(complement,algebra),
                                hopf_invariant(complement2,algebra)],
    )
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C250818_uqsl2())
