using Topology

"""Shared Julia execution of the topology cells in the invariant Live Scripts."""
function run_invariant_workflows()
    lens=VirtualLink(); set_data!(lens;gauss=[[1,-3,2,-2,3,-1]],orientation=[1,1,0])
    calculate_weight(lens;assign=true)
    s2xs1=VirtualLink()
    set_data!(s2xs1;gauss=[[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],
              orientation=[1,1,1,1,0,0])
    calculate_weight(s2xs1;assign=true)
    koda=VirtualLink()
    set_data!(koda;rgauss=[[2,-1,-2,1,4,-3,-4,3]],orientation=[-1,1,-1,1])
    calculate_weight(koda;assign=true)
    abalone=VirtualLink(); set_data!(abalone;rgauss=[[1,-1]],orientation=[1])
    calculate_weight(abalone;assign=true)

    z2=uqsl2_borel_small(1//2;style=:L)
    z5=uqsl2_borel_small(2//5;style=:L)
    sweedler=sweedler_algebra()
    exterior=exterior_algebra(2)
    Dict(
        "lens_uqsl2_5"=>hopf_invariant(lens,z5),
        "s2xs1_uqsl2_2"=>hopf_invariant(s2xs1,z2),
        "abalone_uqsl2_2"=>hopf_invariant(abalone,z2),
        "abalone_sweedler"=>hopf_invariant(abalone,sweedler),
        "koda_sweedler"=>hopf_invariant(koda,sweedler),
        "lens_exterior_2"=>hopf_invariant(lens,exterior),
        "s2xs1_exterior_2"=>hopf_invariant(s2xs1,exterior),
        "abalone_exterior_2"=>hopf_invariant(abalone,exterior),
        "koda_exterior_2"=>hopf_invariant(koda,exterior),
        "koda_ur_prime"=>ur_prime_invariant(koda),
    )
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    display(run_invariant_workflows())
end
