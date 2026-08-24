using Topology

"""Direct non-weighted MP cells from `C250811MSTinvUR.mlx`.

The source resets the same Koda o-graph before its A2/B2/B4/A4 trials, then
continues the A4 result with C2.  This reproduces those executable move cells
without pretending to implement MATLAB's separate weighted/inverse branch.
"""
function run_C250811_mp_moves()
    seed()=begin
        link=VirtualLink()
        set_data!(link;rgauss=[[-1,-2,3,4,-4,-3,2,1]],orientation=[-1,1,-1,1])
        link
    end
    trials=Dict{String,VirtualLink}()
    for (parameter,vertices) in (("A2",[2,1]),("B2",[1,2]),("B4",[2,1]),("A4",[1,2]))
        link=seed(); mp_move!(link,parameter,vertices); trials[parameter]=link
    end
    chained=deepcopy(trials["A4"])
    mp_move!(chained,"C2",[3,5])
    (trials=trials,chained=chained)
end

"""Replay the three explicit URDiagram reduction sequences in C250811 cells 13, 15, and 16."""
function run_C250811_ur_reductions()
    cell13=URDiagram(); set_data_cvw!(cell13,1,[-5,1,-4,4,3,2,-3,5,-2,-1],[-1,-1,-1,1,1])
    for pair in ((-5,1),(5,6),(-6,-2),(3,2),(1,-5),(5,7),(-2,-7)); swap!(cell13,pair); end
    add_edges!(cell13,(6,7)); delete_edge!(cell13,6); dilation!(cell13,3,5); delete_edge!(cell13,5)
    add_edges!(cell13,(2,3)); dilation!(cell13,2,-1); swap!(cell13,(-1,2)); add_edges!(cell13,(1,5))
    swap!(cell13,(-1,-2)); dilation!(cell13,4,1); delete_edge!(cell13,1); swap!(cell13,(2,-2))
    dilation!(cell13,4,-2); add_edges!(cell13,(2,4))

    cell15=URDiagram(); set_data_cvw!(cell15,1,[1,-3,3,-2,2,-1],[1,1,1])
    dilation!(cell15,3,-2); add_edges!(cell15,(2,3)); dilation!(cell15,2,-1); delete_edge!(cell15,1)

    cell16=URDiagram(); set_data_cvw!(cell16,1,[-5,1,-4,4,3,2,-3,5,-2,-1],[-1,-1,-1,1,1])
    swap!(cell16,(2,3)); dilation!(cell16,3,5); swap!(cell16,(2,5)); add_edges!(cell16,(2,3))
    swap!(cell16,(-5,1)); swap!(cell16,(5,6)); add_edges!(cell16,(4,5)); add_edges!(cell16,(2,6))
    dilation!(cell16,4,1); add_edges!(cell16,(1,2)); swap!(cell16,(-1,1)); compose_dilation!(cell16,(1,4))
    (cell13=cell13,cell15=cell15,cell16=cell16)
end

"""Exact common-epsilon cyclic-shift sweep from C250811 cell 9."""
function run_C250811_epsilon_sweep()
    epsilon=laurent_variable(:epsilon)
    code=[-1,-2,3,4,-4,-3,2,1]
    values=Rational{BigInt}[]
    rational_forms=LaurentRational[]
    for _ in eachindex(code)
        link=VirtualLink()
        set_data!(link;rgauss=[code],orientation=[-1,1,-1,1])
        diagram=vl_to_urdiagram(link)
        # MATLAB substitutes every generated eps_i with the single `ep`.
        diagram.W=Any[orientation+epsilon for orientation in [-1,1,-1,1]]
        simplify2!(diagram)
        value=-(diagram.C*epsilon/diagram.W[1])
        push!(rational_forms,value); push!(values,limit_zero(value))
        code=circshift(code,-1)
    end
    (values=values,rational_forms=rational_forms)
end

"""Conversion and trace cells 1--8 from C250811, including multi-rank reductions."""
function run_C250811_conversion_cells()
    cases=(
        hopf=([[1,-2],[-1,2]],[1,1]),
        knot=([[1,2,-2,-3],[-4,-1,3,4]],[-1,1,-1,1]),
        three_component=([[1,2,-2,-3,-5],[-4,-1,3,4],[5]],[-1,1,-1,1,1]),
        lens=([[1,-3,2,-2,3,-1]],[1,1,0]),
        s2xs1=([[-1,-2,3,4,-4,-3,2,1]],[1,1,1,1]),
        two_zero=([[-1,5,-2,3,-6,4,-4,6,-3,2,-5,1]],[1,1,1,1,0,0]),
        virtual=([[2,-1,-2,1,4,-3,-4,3]],[-1,1,-1,1]),
    )
    result=Dict{Symbol,NamedTuple}()
    for (name,(gauss,orientation)) in pairs(cases)
        link=VirtualLink(); set_data!(link;gauss=gauss,orientation=orientation)
        diagram=vl_to_urdiagram(link)
        trace=trace_value(diagram)
        result[name]=(link=link,diagram=diagram,trace=trace)
    end
    result
end

"""Cell 11 weighted o-graph invariants, using the pinned q=exp(4πi/5) small Borel algebra."""
function run_C250811_weighted_invariant()
    link=VirtualLink()
    set_data!(link;gauss=[[4,-5,3,-4,2,-3,1,-2,5,-1]],orientation=ones(Int,5))
    base,kernel,orthogonal=calculate_weight(link;assign=true)
    algebra=uqsl2_borel_small(2//5;style=:L)
    (link=link,base=base,kernel=kernel,orthogonal=orthogonal,
     ur=ur_invariant(link),hopf=hopf_invariant(link,algebra))
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && display((moves=run_C250811_mp_moves(),reductions=run_C250811_ur_reductions(),epsilon_sweep=run_C250811_epsilon_sweep(),conversions=run_C250811_conversion_cells(),weighted=run_C250811_weighted_invariant()))
