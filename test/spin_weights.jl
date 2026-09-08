using Test, Topology, LinearAlgebra

@testset "Z2 storage and domain boundary" begin
    v=VirtualLink(); set_data!(v;rgauss=[1,-1],orientation=[1])
    @test weight_mode(v)==:integer
    set_weight!(v,[3,-2]); set_weight_mode!(v,:z2)
    @test v.weights==[[1,0]]
    @test weight_mode(deepcopy(v))==:z2
    set_weight!(v,[big(2)^100+1,-big(2)^100-1])
    @test v.weights==[[1,1]]
    pair=reshape(collect(first(real_edge_table(v)).crossing),1,2)
    set_weight_by_crossing!(v,pair,[big(2)^100+1])
    @test first(real_edge_table(v)).weight==1
    prior=deepcopy(v.weights)
    @test_throws ArgumentError set_weight!(v,[0,1//2])
    @test_throws ArgumentError set_weight_by_crossing!(v,pair,[1//2])
    @test v.weights==prior
    @test validate_weight(v)
    @test_throws ArgumentError validate_weight(v;framed=true)
    @test_throws ArgumentError weight_equation(v)
    @test_throws ArgumentError calculate_weight(v)
    @test_throws ArgumentError solve_integer_weights(v)
    @test_throws ArgumentError vl_to_urdiagram(v)
    @test_throws ArgumentError set_weight_mode!(v,:integer)
    set_weight!(v,[3,-2];mode=:integer)
    @test weight_mode(v)==:integer && v.weights==[[3,-2]]
    @test_throws ArgumentError VirtualLink(;weight_mode=:invalid)
end

@testset "Exact F2 affine solve versus exhaustive assignments" begin
    # Includes empty systems, repeated rows, a nontrivial free coordinate,
    # and an inconsistent system with an independently checked dual witness.
    systems=[(zeros(Int,0,2),Int[]),(zeros(Int,2,0),[0,1]),
             ([1 1 0;0 1 1],[1,0]),([1 1;1 1],[0,1]),
             ([2 3 -1;4 -1 1],[2,3])]
    for (A,b) in systems
        result=solve_spin_weight_system(A,b)
        n=size(A,2)
        allbits(k)=[Int[(i >> (j-1)) & 1 for j in 1:k] for i in 0:(2^k-1)]
        oracle=Set(Tuple(x) for x in allbits(n) if all(iszero,mod.(A*x-b,2)))
        @test (result.status==:spin_solution)==!isempty(oracle)
        @test mod.(result.certificate.row_transform*hcat(mod.(A,2),mod.(b,2)),2)==result.certificate.rref
        if isempty(oracle)
            row=result.obstruction.row_combination
            @test all(iszero,mod.(transpose(row)*A,2))
            @test mod(dot(row,b),2)==1
        else
            actual=Set(Tuple(mod.(result.particular+result.kernel*z,2)) for z in allbits(size(result.kernel,2)))
            @test actual==oracle
        end
    end
    @test_throws ArgumentError solve_spin_weight_system(reshape([0.5],1,1),[0])
end

@testset "Spin S equation and integer non-solvability are distinct" begin
    # The concrete CP RHS with integral obstruction 3*S=1 is Z2 solvable.
    v=VirtualLink(;weight_mode=:z2)
    set_data!(v;rgauss=[1,2,-5,3,-2,-3,-4,4,-1,5,6,-7,7,-6],
              orientation=[1,1,1,-1,-1,-1,1])
    solution=solve_spin_weights(v;assign=true)
    @test solution.status==:spin_solution
    @test solution.assignment_status==:assigned
    @test is_spin_weight(v)
    @test !is_framed_weight(v)
    disks=disk_table(v); w=Int[e.weight for e in real_edge_table(v)]
    # Independent form from Mih23 S: unsigned occurrence sum = Ndots/2+1.
    for row in axes(disks.delta,1)
        @test mod(sum(abs.(disks.delta[row,:]).*w),2)==mod(div(disks.Ndots[row],2)+1,2)
    end
    integer=deepcopy(v); set_weight!(integer,w;mode=:integer)
    @test solve_integer_weights(integer).status==:rational_only
    cut=VirtualLink(;weight_mode=:z2)
    set_data!(cut;rgauss=[1,-1,0],orientation=[1]); set_weight!(cut,[0,0,0])
    @test !is_spin_weight(cut)
    @test_throws ArgumentError solve_spin_weights(cut)
end
