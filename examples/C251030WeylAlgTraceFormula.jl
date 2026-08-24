using Topology

"""Julia port of the topology/trace-formula cells in C251030WeylAlgTraceFormula.mlx."""
function run_C251030_trace_formula(matrices)
    [(diagram=begin
          d=URDiagram(); set_data_cm!(d,1,matrix); d
      end,
      trace=standard_matrix_trace(1,matrix)) for matrix in matrices]
end

"""Symbolic standard-matrix cells for N=1 through `maximum_order`.

Each entry is the exact counterpart of the Live Script's `formula`, `B`, and
`dets` arrays: negative reciprocal of det(M+B_N), the cyclic correction matrix, and its
determinant.  Independent Laurent variables avoid a symbolic-CAS dependency.
"""
function run_C251030_standard_forms(;maximum_order=4)
    maximum_order>=1 || throw(ArgumentError("maximum_order must be positive"))
    [begin
        variables=laurent_variables(order^2;prefix="a$(order)_")
        matrix=reshape(variables,order,order)
        adjustment=standard_matrix_adjustment(matrix)
        determinant=Topology._det(matrix+adjustment)
        (order=order,matrix=matrix,adjustment=adjustment,
         determinant=determinant,trace=standard_matrix_trace(one(variables[1]),matrix))
    end for order in 1:maximum_order]
end

"""The three short symbolic reduction cells 5--7 in C251030."""
function run_C251030_symbolic_reductions()
    t=laurent_variable(:t)
    cases=(
        ([-2,1,-1,2],[t,-1]),
        ([-1,1,2,-2],[one(t)-t,t-1]),
    )
    result=URDiagram[]
    for (vertices,weights) in cases
        diagram=URDiagram(); set_data_cvw!(diagram,1,vertices,weights); simplify2!(diagram)
        push!(result,diagram)
    end
    # Cell 7 has independent s and t, so it uses one common bivariate Laurent
    # ring rather than incompatible one-variable coefficient rings.
    s,two_t=laurent_variables(2;prefix="t")
    cell7=URDiagram()
    set_data_cvw!(cell7,1,[-2,-1,1,2],[one(s)-s,two_t-one(s)])
    simplify2!(cell7); push!(result,cell7)
    result
end

"""Knot-variant trace cells 8 in C251030, including its singular second variant."""
function run_C251030_knot_trace_variants()
    t=laurent_variable(:t)
    first=URDiagram(); set_data_cvw!(first,1,
        [-7,-1,2,5,-5,-3,1,4,-4,-2,3,6,-6,7],
        [one(t)-t,one(t)-t,one(t)-t,t-one(t),t-one(t),t-one(t),-1])
    direct=trace3(deepcopy(first)); simplify2!(first); reduced_trace=trace2(first)
    second=URDiagram(); set_data_cvw!(second,1,
        [-1,2,5,-5,-3,1,4,-4,-2,3,6,-6],
        [one(t)-t,one(t)-t,one(t)-t,t-one(t),t-one(t),t-one(t)])
    simplify2!(second)
    singular=try
        trace2(second); nothing
    catch error
        typeof(error)
    end
    (first_direct=direct,first_reduced=reduced_trace,first_diagram=first,
     second_diagram=second,second_trace_failure=singular)
end

"""Run C251030's final matrix-to-URDiagram reduction at numeric coefficients.

This is a structural replay, not a substitution test for the symbolic
determinant formula: zero entries follow a different numeric reduction path.
The returned state is checked against the MATLAB source implementation.
"""
function run_C251030_numeric_reduction_probe(matrix::AbstractMatrix{<:Number})
    rows,columns=size(matrix); rows==columns || throw(ArgumentError("matrix must be square"))
    diagram=URDiagram(); set_data_cm!(diagram,1,matrix)
    filter!(v -> v!=0,diagram.V)
    simplify2!(diagram)
    (C=diagram.C,V=copy(diagram.V),W=copy(diagram.W),inverse_trace=inv(trace2(diagram)))
end

"""Port the final N=4 symbolic verification cell in C251030.

The MATLAB cell creates a generic square matrix `A`, forms the cyclic
correction `I-circshift(I,1)`, and verifies
`1/D.trace + det(A + correction) == 0`.  The Julia result retains all three
terms so the identity is an executable assertion rather than a prose claim.
"""
function run_C251030_final_matrix_workflow(;order::Integer=4)
    order>=1 || throw(ArgumentError("order must be positive"))
    variables=laurent_variables(order^2;prefix="a$(order)_")
    matrix=reshape(variables,order,order)
    correction=standard_matrix_adjustment(matrix)
    trace=standard_matrix_trace(one(variables[1]),matrix)
    observed=inv(trace)
    expected=-Topology._det(matrix+correction)
    (order=order,matrix=matrix,correction=correction,trace=trace,
     observed=observed,expected=expected,residual=observed-expected)
end

if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    matrices=[reshape([2.0],1,1),[2.0 1.0; 3.0 4.0],
              [2.0 1.0 0.0; 0.0 3.0 1.0; 1.0 0.0 4.0]]
    (numeric=run_C251030_trace_formula(matrices),symbolic=run_C251030_standard_forms(),
     reductions=run_C251030_symbolic_reductions(),variants=run_C251030_knot_trace_variants(),
     final_matrix=run_C251030_final_matrix_workflow())
end
