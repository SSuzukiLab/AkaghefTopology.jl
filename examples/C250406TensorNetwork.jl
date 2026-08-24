using Topology

"""Replay the three stored Sweedler tetrahedron contractions from C250406."""
function run_C250406_tensor_network()
    algebra=sweedler_algebra()
    tensors=legacy_hopf_tensors(algebra)
    expression1=evaluate_tensor_expression(
        "Mu{1,2,3}Il{3}An{4,2}De{7,5,6}De{6,4,8}(An^-2){8,9}" *
        "Mu{9,12,10}Mu{10,5,11}De{13,14,12}Cl{13}(An){14,15}",
        tensors,[11,1,7,15];dimension=4)
    expression2=evaluate_tensor_expression(
        "Mu{5,6,1}Mu{7,8,5}An{9,8}De{20,9,10}De{3,6,11}" *
        "De{11,12,13}An{12,14}Mu{2,14,15}Mu{15,16,17}(An^2){17,18}" *
        "De{18,16,19}(An^-2){19,20}(An^2){10,4}(An^2){13,7}",
        tensors,1:4;dimension=4)
    expression3=evaluate_tensor_expression(
        "Tp{3,5,6,1}Tm{11,5,9,2}Tm{10,7,6,12}Il{9}Cr{8}" *
        "(An^-2){8,7}(An^2){10,4}(An^-2){11,12}",
        tensors,1:4;dimension=4)
    dense_tensor(expression1,4),dense_tensor(expression2,4),dense_tensor(expression3,4)
end

"""The introductory Sweedler contractions preceding the tetrahedron cell."""
function run_C250406_exploratory_cells()
    tensors=legacy_hopf_tensors(sweedler_algebra())
    evaluate(expression,free)=dense_tensor(evaluate_tensor_expression(expression,tensors,free;dimension=4),4)
    (
        integral_cointegral=evaluate("Mu{1,2,3}Cl{3}(An^-2){1,2}",Int[]),
        coproduct_product=evaluate("Mu{1,2,3}(An^2){3,4}De{4,2,5}",[1,5]),
        integral_pair=evaluate("Il{2}Cr{1}",[2,1]),
        tp_twist=evaluate("Tp{0,1,2,3}(An^2){3,0}",[2,1]),
        tm_twist=evaluate("Tm{0,1,2,3}(An^2){0,3}",[1,2]),
        counit_unit=evaluate("Ep{1}Et{2}",[1,2]),
        mixed=evaluate("Mu{6,1,3}(An^2){3,4}De{4,5,2}An{5,6}",[1,2]),
    )
end

"""Sweedler framing comparison cells from the `問1` section."""
function run_C250406_sweedler_framing()
    tensors=legacy_hopf_tensors(sweedler_algebra())
    evaluate(expression)=dense_tensor(evaluate_tensor_expression(expression,tensors,1:4;dimension=4),4)
    full=evaluate("Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}Tp{16,9,10,15}" *
        "(An^4){7,6}(An^-2){9,8}(An^-4){11,10}(An^2){1,3}(An^-2){14,2}" *
        "(An^2){15,16}(An^-2){13,12}(An^2){4,14}")
    s4_identity=evaluate("Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}Tp{16,9,10,15}" *
        "(Un){7,6}(An^-2){9,8}(Un){11,10}(An^2){1,3}(An^-2){14,2}" *
        "(An^2){15,16}(An^-2){13,12}(An^2){4,14}")
    changed=evaluate("Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}Tp{16,9,10,15}" *
        "(Un){7,6}(An^-2){9,8}(Un){11,10}(An^2){1,3}(Un){4,2}(An^2){15,16}")
    (full=full,s4_identity=s4_identity,changed=changed)
end

"""Four distinct Sweedler twist cells from the `問1` exploration."""
function run_C250406_sweedler_twists()
    tensors=legacy_hopf_tensors(sweedler_algebra())
    evaluate(expression)=dense_tensor(evaluate_tensor_expression(expression,tensors,1:4;dimension=4),4)
    (
        twist1=evaluate("Tp{3,5,6,1}Tm{11,5,9,2}Tm{10,7,6,12}Tp{15,8,9,13}" *
            "(An^2){13,15}(An^-2){8,7}(An^2){10,4}(An^-2){11,12}"),
        twist2=evaluate("Tm{1,6,5,3}Tp{2,9,5,11}Tp{12,6,7,10}Tm{13,9,8,15}" *
            "(An^2){13,15}(An^-2){8,7}(An^2){10,4}(An^-2){11,12}"),
        twist3=evaluate("Tp{3,5,1,6}Tm{5,11,2,9}Tm{7,10,12,6}Tp{8,15,13,9}" *
            "(An^2){15,13}(An^-2){7,8}(An^2){4,10}(An^-2){12,11}"),
        twist4=evaluate("Tm{3,5,1,6}Tp{9,2,11,5}Tp{6,12,10,7}Tm{9,13,15,8}" *
            "(An^2){15,13}(An^-2){7,8}(An^2){4,10}(An^-2){12,11}"),
    )
end

"""Ring-addition cells from C250406 for Sweedler, small Borel, and cyclic Hopf algebras."""
function run_C250406_ring_cells()
    function prepare(algebra)
        tensors=legacy_hopf_tensors(algebra)
        tensors["ring"]=evaluate_tensor_expression(
            "Tp{1,2,3,4}(An^-2){1,2}",tensors,[3,4];dimension=algebra.dimension)
        tensors
    end
    evaluate(expression,tensors,dimension)=dense_tensor(
        evaluate_tensor_expression(expression,tensors,1:4;dimension=dimension),dimension)

    sweedler=prepare(sweedler_algebra())
    sweedler_with_ring=evaluate(
        "Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}Tp{16,9,10,15}" *
        "(An^4){7,6}(An^-2){9,8}(An^-4){11,10}(An^2){1,3}(An^-2){14,2}" *
        "(An^2){15,16}(An^-2){17,12}ring{13,17}(An^2){4,14}",sweedler,4)

    small_borel_3=prepare(uqsl2_borel_small(1//3;style=:L))
    small_borel_5=prepare(uqsl2_borel_small(1//5;style=:L))
    cyclic=prepare(cyclic_group_algebra(5))
    cyclic_with_ring=evaluate(
        "Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}(An^-2){9,8}Tp{16,9,10,15}" *
        "(An^4){7,6}(An^-4){11,10}(An^2){1,3}(An^-2){14,2}(An^2){15,16}" *
        "ring{13,17}(An^-2){17,12}(An^2){4,14}",cyclic,5)
    cyclic_without_ring=evaluate(
        "Tp{5,3,6,1}Tm{5,13,2,11}Tm{8,4,12,7}(An^-2){9,8}Tp{16,9,10,15}" *
        "(An^4){7,6}(An^-4){11,10}(An^2){1,3}(An^-2){14,2}(An^2){15,16}" *
        "Un{13,17}(An^-2){17,12}(An^2){4,14}",cyclic,5)
    (
        small_borel_3_checks=verify_hopf(uqsl2_borel_small(1//3;style=:L)),
        cyclic_checks=verify_hopf(cyclic_group_algebra(5)),
        sweedler_ring=dense_tensor(sweedler["ring"],4),
        sweedler_with_ring=sweedler_with_ring,
        small_borel_3_ring=dense_tensor(small_borel_3["ring"],9),
        small_borel_5_ring=dense_tensor(small_borel_5["ring"],25),
        cyclic_ring=dense_tensor(cyclic["ring"],5),
        cyclic_with_ring=cyclic_with_ring,
        cyclic_without_ring=cyclic_without_ring,
    )
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    (tetrahedron=first(run_C250406_tensor_network()),exploratory=run_C250406_exploratory_cells(),
     framing=run_C250406_sweedler_framing(),twists=run_C250406_sweedler_twists(),
     rings=run_C250406_ring_cells())
end
