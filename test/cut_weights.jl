using Test
using Topology

@testset "cut weight equations use explicit relative boundary targets" begin
    closed=lens_space_link(2,1)
    calculate_weight(closed;assign=true)
    closed_data=boundary_weight_data(closed)
    @test all(iszero,closed_data.closed_residual)
    @test closed_data.measured+2 .* BigInt.(closed_data.isCycle)==closed_data.closed_residual
    @test validate_weight(closed;framed=true)
    @test weight_equation(closed).boundary_targets==zeros(BigInt,length(closed_data.isCycle))

    cut=VirtualLink()
    set_data!(cut;rgauss=[[1,-1]],orientation=[1],cut=[true])
    @test length(real_edge_table(cut))==3
    set_weight!(cut,zeros(length(real_edge_table(cut))))
    interval=VirtualLink()
    set_data!(interval;gauss=[Int[]],orientation=Int[],cut=[true])
    @test length(real_edge_table(interval))==1

    cut=lens_space_link(2,1)
    cut.gauss[1]=[cut.gauss[1]...,0]
    set_weight!(cut,zeros(length(real_edge_table(cut))))
    data=boundary_weight_data(cut)
    closed_rows=findall(identity,data.isCycle)
    closed_solution=solve_integer_weight_system(
        2 .* data.delta[closed_rows,:],-data.cp2[closed_rows])
    @test closed_solution.status==:integer_solution
    set_weight!(cut,Float64.(closed_solution.particular))
    data=boundary_weight_data(cut)
    @test data.measured==2 .* data.delta*data.edge_weights-data.Ndots
    open_rows=findall(value -> !value,data.isCycle)
    @test !isempty(open_rows)
    targets=Dict((data.boundary_start[row],data.boundary_end[row])=>data.measured[row]
                 for row in open_rows)

    @test_throws ArgumentError weight_equation(cut)
    @test_throws ArgumentError validate_weight(cut;framed=true)
    @test validate_weight(cut)
    equation=weight_equation(cut;boundary=targets)
    @test all(iszero,equation.residual)
    @test weight_residual(cut;boundary=targets)==zeros(BigInt,length(data.isCycle))
    @test validate_weight(cut;framed=true,boundary=targets)
    @test !is_framed_weight(cut;boundary=targets)

    solved=solve_integer_weights(cut;boundary=targets)
    @test solved.status==:integer_solution
    @test solved.A*solved.particular==solved.b
    lifted=solve_integer_weight_lift(cut,mod.(solved.particular,2);boundary=targets)
    @test lifted.status==:integer_solution
    @test lifted.particular==solved.particular
    @test verify_integer_weight_solution(lifted)

    assigned=deepcopy(cut)
    assigned_result=solve_integer_weights(assigned;boundary=targets,assign=true)
    @test assigned_result.assignment_status==:assigned
    @test validate_weight(assigned;framed=true,boundary=targets)
end
