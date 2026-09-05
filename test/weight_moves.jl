using Test
using Topology

@testset "integral weight equation and MATLAB weighted moves" begin
    weighted=VirtualLink()
    set_data!(weighted;
              gauss=[[1,-1,2,-3],[-2,3]],
              orientation=[1,1,0])
    base,kernel,hperp=calculate_weight(weighted;assign=true)

    @test weight_residual(weighted)==[0,0,0,0]
    @test validate_weight(weighted)

    lens=lens_space_link(2,1)
    _,_,lens_hperp=calculate_weight(lens)
    @test first_betti_number(lens)==size(lens_hperp,1)

    move!(weighted,"H";param=-2,vertices=[2])
    @test validate_weight(weighted)
    @test weighted.weights==[[-3.0,1.0,1.0],[-1.0]]

    move!(weighted,"MP";param="C3",vertices=[1,2])
    @test validate_weight(weighted)
    @test weighted.gauss==[[3,1,2,-3],[-1,-2]]
    @test weighted.orientation==[-1,1,1]
    @test weighted.weights==[[-3.0,1.0,0.0,0.0],[-1.0,0.0]]
end

@testset "integral 0-2 inverse weight transport" begin
    weighted=VirtualLink()
    set_data!(weighted;rgauss=[[-1,-2,2,1]],orientation=[1,-1])
    set_weight!(weighted,[[4,7,11,13]])

    @test_throws ArgumentError move!(weighted,"02";param="inv",vertices=[1,2])

    @test weight_residual(weighted) isa Vector{Int}
end

@testset "weights follow real arcs through virtual crossings" begin
    weighted=VirtualLink()
    set_data!(weighted;gauss=[[1,2,-1,-2]],orientation=[1,0])
    set_weight!(weighted,[[5,7]])

    @test [edge.weight for edge in real_edge_table(weighted)]==[5.0,7.0]
end

@testset "GaussCode lens weights preserve the disk equation" begin
    fixtures=((2,[-1,-1,1,0],[8,3,1],
               [2 0 2 0;-1 1 0 1;0 0 -1 0],[0,0,2]),
              (3,[-2//3,-2//3,-2//3,1,0,0],[10,3,4,1],
               [2 1 0 2 0 1;-1 1 0 0 0 1;0 -1 1 0 1 -1;0 0 0 -1 0 0],
               [0,0,0,2]),
              (4,[-1//2,-1//2,-1//2,-1//2,1,0,0,0],[12,3,4,4,1],
               [2 1 1 0 2 0 1 1;-1 1 0 0 0 0 0 1;0 -1 1 0 0 0 1 -1;
                0 0 -1 1 0 1 -1 0;0 0 0 0 -1 0 0 0],[0,0,0,0,2]))
    for (n,expected,expected_paths,expected_delta,expected_cp2) in fixtures
        arr=collect(1:2n-1)
        orientation=Int.((1 .- (-1).^arr)./2)
        signed=arr.*(-1).^arr
        lens=VirtualLink()
        set_data!(lens;gauss=[vcat(-signed,reverse(signed))],orientation=orientation)
        disks=disk_table(lens)
        base,_,_=calculate_weight(lens;assign=true)
        @test disks.path_lengths==expected_paths
        @test disks.delta==expected_delta
        @test disks.cp2==expected_cp2
        @test base==expected
        @test 2*disks.delta*base+disks.cp2==zeros(Rational{BigInt},size(disks.delta,1))
        @test lens.weights==[Float64.(base)]
    end
end

@testset "MATLAB-faithful L(2,1) disk table" begin
    disks=disk_table(lens_space_link(2,1))

    @test disks.ID==[1,2,3]
    @test disks.paths==[[1,3,4,10,8,6,11,7],[2,5,9],[12]]
    @test disks.edges_s==[[3,3,4,1,-4,-2,1,2],[-3,4,2],[-1]]
    @test disks.Ndots==[2,2,0]
    @test disks.isCycle==[true,true,true]
    @test disks.delta==[2 0 2 0;0 1 -1 1;-1 0 0 0]
    @test disks.cp==[0//1,0//1,1//1]
    @test disks.cp2==[0,0,2]
    @test disks.path_lengths==[8,3,1]

    @test [sum(sign == edge ? 1 : sign == -edge ? -1 : 0
               for sign in disk_edges)
           for disk_edges in disks.edges_s, edge in 1:size(disks.delta,2)] ==
          disks.delta
end

@testset "L(2,1) H move accumulates repeated edge incidences" begin
    base=[1,-1,-1,0]
    lens=lens_space_link(2,1)
    disks=disk_table(lens)
    hmoves=hmove_matrix(lens)
    cocycle=[0,1,0,-1]

    @test disks.delta*cocycle==zeros(Int,size(disks.delta,1))
    @test hmoves==hcat(-cocycle,cocycle)

    for vertex in 1:2
        moved=lens_space_link(2,1)
        set_weight!(moved,base)
        move!(moved,"H";param=1,vertices=[vertex])
        moved_weights=Int.(round.([edge.weight for edge in real_edge_table(moved)]))

        @test moved_weights-base==hmoves[:,vertex]
        @test 2*disks.delta*moved_weights+disks.cp2==zeros(Int,size(disks.delta,1))
        @test validate_weight(moved)
    end
end

@testset "integer weight solutions take precedence over rational pivots" begin
    A=[2 0 1]
    b=[1]
    rational=Topology._rational_particular(A,b)
    integer=Topology._integer_particular(A,b)

    @test rational==[1//2,0,0]
    @test integer isa Vector{Int}
    @test A*integer==b
    @test integer != Int[0,0,0]
end

@testset "Hopf weight exponents require integral weights" begin
    fractional=VirtualLink()
    set_data!(fractional;rgauss=[[1,-1]],orientation=[1])
    set_weight!(fractional,[[0.5,0.0]])

    @test_throws ArgumentError hopf_invariant(fractional,exterior_algebra(2))
end
