using Topology

"""Reproduce every invariant family in `C250823MSTinv_extalg.mlx`."""
function run_C250823_extalg()
    algebra=exterior_algebra(2)
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
    s2_values=ComplexF64[]
    for multiple in 0:3
        set_weight!(s2xs1,[s2_base+multiple*vec(s2_hperp)])
        push!(s2_values,hopf_invariant(s2xs1,algebra))
    end

    koda=VirtualLink()
    set_data!(koda;rgauss=[[2,-1,-2,1,4,-3,-4,3]],orientation=[-1,1,-1,1])
    calculate_weight(koda;assign=true)
    (
        hopf_checks=verify_hopf(algebra;degrees=count_ones.(0:3)),
        lens_values=lens_values,
        abalone_value=hopf_invariant(abalone,algebra),
        s2xs1_values=s2_values,
        koda_value=hopf_invariant(koda,algebra),
    )
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C250823_extalg())
