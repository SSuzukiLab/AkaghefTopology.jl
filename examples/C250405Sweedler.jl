using Topology

"""Replay `C250405Sweedler.mlx`, making its hidden initial `vl` dependency explicit."""
function run_C250405_sweedler(;initial_link=nothing)
    sweedler=sweedler_algebra()
    if initial_link===nothing
        initial_link=VirtualLink()
        set_data!(initial_link;rgauss=[[1,-1]],orientation=[1])
    end
    first=deepcopy(initial_link)
    set_data!(first;rgauss=[[1,-1]],orientation=[1])
    calculate_weight(first;assign=true)

    koda=VirtualLink()
    set_data!(koda;rgauss=[[2,-1,-2,1,4,-3,-4,3]],orientation=[-1,1,-1,1])
    calculate_weight(koda;assign=true)
    (
        hopf_checks=verify_hopf(sweedler),
        abalone_invariant=hopf_invariant(first,sweedler),
        koda_invariant=hopf_invariant(koda,sweedler),
        abalone=first,
        koda=koda,
    )
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C250405_sweedler())
