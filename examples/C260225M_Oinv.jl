using Topology
using LinearAlgebra

"""Reproduce the executable topology/algebra cells in `C260225M_Oinv.mlx`."""
function run_C260225_m_oinv()
    link=VirtualLink()
    set_data!(link;rgauss=[[-1,-2,3,4,-4,-3,2,1]],orientation=[1,1,1,1])
    base,kernel,hperp=calculate_weight(link;assign=true)

    uqsl2=uqsl2_borel_small(1//3;style=:L)
    antipode=zeros(ComplexF64,uqsl2.dimension,uqsl2.dimension)
    for ((output,input),value) in uqsl2.antipode.terms
        antipode[output,input]=value
    end
    cyclic=cyclic_group_algebra(3)
    (
        link=link,
        weight=(base=base,kernel=kernel,hperp=hperp),
        antipode_squared=antipode^2,
        cyclic_invariant=hopf_invariant(link,cyclic),
        uqsl2_checks=verify_hopf(uqsl2),
        cyclic_checks=verify_hopf(cyclic),
    )
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display(run_C260225_m_oinv())
