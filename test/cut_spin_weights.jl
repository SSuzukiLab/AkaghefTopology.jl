using Test
using Topology

function _spin_cut_fixture(;assigned=true)
    v=VirtualLink(;weight_mode=:z2)
    set_data!(v;rgauss=[-2,-1,1,2,0],orientation=[1,1])
    assigned && set_weight!(v,zeros(length(real_edge_table(v)));mode=:z2)
    v
end

@testset "cut spin weights use explicit Z4 boundary targets" begin
    cut=_spin_cut_fixture()
    data=spin_boundary_data(cut)
    open_rows=findall(!,data.isCycle)
    @test open_rows==[1,3,4]
    @test all(value -> value in 0:3,data.measured4)
    @test data.boundary_identities[1]==(data.boundary_start[1],data.boundary_end[1])

    # The target map is keyed by the directed endpoint pair and preserves the
    # independently measured relative values of this presentation.
    targets=Dict((data.boundary_start[row],data.boundary_end[row])=>data.measured4[row]
                 for row in open_rows)
    equation=spin_weight_equation(cut;boundary=targets)
    @test equation.modulus==4
    @test equation.boundary_targets==BigInt[0,2,2,0]
    @test equation.measured4==data.measured4
    @test all(iszero,equation.residual)
    @test validate_spin_weight(cut;boundary=targets)
    @test is_spin_weight(cut;boundary=targets)
    @test !is_spin_weight(cut)

    # A row-aligned spelling is accepted, while closed rows retain fixed
    # target 2 and open targets remain explicit.
    row_targets=Any[data.measured4[1],2,data.measured4[3],data.measured4[4]]
    @test spin_weight_equation(cut;boundary=row_targets).boundary_targets==equation.boundary_targets
    @test_throws ArgumentError spin_weight_equation(cut)

    # Z4 parity is checked before division.  This target would require
    # dividing an odd value by two, so it is an explicit obstruction.
    bad=copy(data.measured4)
    bad[1]=1
    bad_targets=Dict((data.boundary_start[row],data.boundary_end[row])=>bad[row]
                     for row in open_rows)
    @test_throws ArgumentError solve_spin_weights(_spin_cut_fixture(;assigned=false);
                                                   boundary=bad_targets)
end

@testset "cut spin weights solve the divided GF2 system" begin
    source=_spin_cut_fixture()
    data=spin_boundary_data(source)
    open_rows=findall(!,data.isCycle)
    targets=Dict((data.boundary_start[row],data.boundary_end[row])=>data.measured4[row]
                 for row in open_rows)
    candidate=_spin_cut_fixture(;assigned=false)
    solved=solve_spin_weights(candidate;boundary=targets,assign=true)
    @test solved.status==:spin_solution
    @test solved.assignment_status==:assigned
    @test solved.modulus==2
    @test solved.boundary_targets==BigInt[0,2,2,0]
    @test validate_spin_weight(candidate;boundary=targets)
    @test spin_weight_equation(candidate;boundary=targets).residual==zeros(BigInt,4)
end

@testset "cut spin boundary key and parity failures are explicit" begin
    odd=VirtualLink(;weight_mode=:z2)
    set_data!(odd;rgauss=[1,-1,0],orientation=[1])
    set_weight!(odd,zeros(length(real_edge_table(odd)));mode=:z2)
    data=spin_boundary_data(odd)
    open_rows=findall(!,data.isCycle)
    @test data.Ndots[1]==1
    valid=Dict((data.boundary_start[row],data.boundary_end[row])=>data.measured4[row]
               for row in open_rows)
    @test spin_weight_equation(odd;boundary=valid).modulus==4
    odd_parity=copy(data.measured4)
    odd_parity[1]=0 # 0 + Ndots[1] is odd, so division by 2 is undefined.
    badparity=Dict((data.boundary_start[row],data.boundary_end[row])=>odd_parity[row]
                   for row in open_rows)
    @test_throws ArgumentError solve_spin_weights(odd;boundary=badparity)

    missing=Dict((data.boundary_start[row],data.boundary_end[row])=>valid[(data.boundary_start[row],data.boundary_end[row])]
                 for row in open_rows[1:end-1])
    @test_throws ArgumentError spin_weight_equation(odd;boundary=missing)
    unknown=Dict(((strand=99,side=:in,lane=1),(strand=99,side=:out,lane=1))=>0)
    @test_throws ArgumentError spin_weight_equation(odd;boundary=unknown)
    @test_throws BoundsError spin_weight_equation(odd;boundary=Dict(99=>0))

    duplicate=Dict{Any,Any}()
    first_row=first(open_rows)
    first_key=(data.boundary_start[first_row],data.boundary_end[first_row])
    duplicate[first_key]=valid[first_key]
    duplicate[first_row]=valid[first_key]
    @test_throws ArgumentError spin_weight_equation(odd;boundary=duplicate)
end

@testset "closed spin equation remains modulo 2" begin
    closed=VirtualLink(;weight_mode=:z2)
    set_data!(closed;rgauss=[1,2,-5,3,-2,-3,-4,4,-1,5,6,-7,7,-6],
              orientation=[1,1,1,-1,-1,-1,1])
    result=solve_spin_weights(closed;assign=true)
    equation=spin_weight_equation(closed)
    @test result.status==:spin_solution
    @test equation.modulus==2
    @test equation.residual==zeros(Int,length(equation.residual))
    @test validate_spin_weight(closed)
end
