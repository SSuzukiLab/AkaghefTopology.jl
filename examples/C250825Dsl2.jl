using Topology

"""Replay the finite Hopf, dual, Drinfeld-double and Heisenberg-double cells."""
function run_C250825_dsl2(;order=3)
    algebra=uqsl2_borel_small(1//order;style=:L)
    dual=dual_algebra(algebra)
    drinfeld=drinfeld_double(algebra)
    heisenberg=heisenberg_double(algebra)
    K=basis_element(algebra,2)
    E=basis_element(algebra,order+1)
    pK=basis_element(dual,2)
    pE=basis_element(dual,order+1)
    (;algebra,dual,drinfeld,heisenberg,K,E,pK,pE,
      hopf_checks=verify_hopf(algebra),dual_checks=verify_hopf(dual),
      pairing=hopf_pairing(pE,E*K),
      commutator=basis_element(heisenberg,1+2order^2)*basis_element(heisenberg,1+order)-
                 basis_element(heisenberg,1+order)*basis_element(heisenberg,1+2order^2))
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    result=run_C250825_dsl2()
    display((result.pairing,result.commutator.coefficients))
end
