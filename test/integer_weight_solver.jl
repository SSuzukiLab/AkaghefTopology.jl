using Test
using Topology
using LinearAlgebra

_integer_solver_raw(v)=(deepcopy(v.gauss),copy(v.orientation),deepcopy(v.weights),
                       v.real_only,v.is_weighted,copy(v.strand_types))

@testset "exact integer feasibility and full kernel lattice" begin
    result=solve_integer_weight_system([2 1 1],[0])
    @test result.status==:integer_solution
    @test result.particular==zeros(BigInt,3)
    @test size(result.kernel)==(3,2)
    @test verify_integer_weight_solution(result)
    # Independent saturation oracle: in rank two the gcd of the maximal
    # minors is the index in the primitive lattice in its rational span.
    K=result.kernel
    minors=BigInt[K[i,1]*K[j,2]-K[i,2]*K[j,1] for (i,j) in ((1,2),(1,3),(2,3))]
    @test gcd(minors...)==1
    missing=BigInt[-1,1,1]
    # Cramer's rule on rows 1,3 expresses a vector omitted by the old
    # independently denominator-cleared basis [-1,2,0],[-1,0,2].
    d=K[1,1]*K[3,2]-K[1,2]*K[3,1]
    coefficients=[(missing[1]*K[3,2]-K[1,2]*missing[3])//d,
                  (K[1,1]*missing[3]-missing[1]*K[3,1])//d]
    @test all(denominator(x)==1 for x in coefficients)
    @test K*BigInt.(numerator.(coefficients))==missing

    # A nontrivial inhomogeneous system, including every solution in a small
    # independently enumerated integer box, must use integer kernel coordinates.
    A=BigInt[2 4 6;0 3 9]; b=A*BigInt[-1,2,0]
    r=solve_integer_weight_system(A,b)
    inverseV=inv(Rational{BigInt}.(r.certificate.V))
    count_solutions=0
    for x1 in -3:3,x2 in -3:3,x3 in -3:3
        x=BigInt[x1,x2,x3]
        A*x==b || continue
        count_solutions+=1
        coordinates=inverseV*(x-r.particular)
        @test all(iszero,coordinates[1:r.certificate.rank])
        @test all(denominator(c)==1 for c in coordinates)
        @test r.kernel*BigInt.(numerator.(coordinates[r.certificate.rank+1:end]))==x-r.particular
    end
    @test count_solutions>0

    # Known Smith factors survive integer elementary changes of coordinates.
    D=BigInt[2 0 0 0;0 6 0 0;0 0 0 0]
    L=BigInt[1 1 0;0 1 1;0 0 1]
    R=BigInt[1 2 0 0;0 1 1 0;0 0 1 1;0 0 0 1]
    A=L*D*R
    for (rhs,status) in (([4,12,0],:integer_solution),([1,0,0],:rational_only),
                         ([4,12,1],:inconsistent))
        r=solve_integer_weight_system(A,L*rhs)
        @test r.status==status
        @test diag(r.certificate.S)[1:2]==[2,6]
        @test verify_integer_weight_solution(r)
    end
end

@testset "recheckable obstructions, huge integers, and empty systems" begin
    rational=solve_integer_weight_system(reshape(BigInt[2],1,1),[1])
    @test rational.status==:rational_only
    @test rational.particular===nothing
    @test rational.rational_particular==[1//big(2)]
    @test rational.obstruction.divisor==2
    @test rational.obstruction.rhs_value==1
    @test rational.obstruction.remainder==1
    inconsistent=solve_integer_weight_system([2 0;0 0],[1,1])
    @test inconsistent.status==:inconsistent
    @test inconsistent.rational_particular===nothing
    o=inconsistent.obstruction
    @test transpose(o.row_combination)*inconsistent.A==zeros(BigInt,1,2)
    @test dot(o.row_combination,inconsistent.b)!=0
    @test o.divisor==0

    huge=big(2)^200+1
    r=solve_integer_weight_system(reshape(BigInt[1],1,1),[huge])
    @test r.particular==[huge]
    @test eltype(r.particular)==BigInt
    obstruction=solve_integer_weight_system(reshape(BigInt[big(2)^100],1,1),[1])
    @test obstruction.status==:rational_only
    @test obstruction.obstruction.divisor==big(2)^100
    @test verify_integer_weight_solution(obstruction)

    for (m,n,b,status) in ((0,0,Int[],:integer_solution),(0,3,Int[],:integer_solution),
                           (2,0,[0,0],:integer_solution),(2,0,[0,1],:inconsistent),
                           (2,3,[0,0],:integer_solution),(2,3,[1,0],:inconsistent))
        r=solve_integer_weight_system(zeros(Int,m,n),b)
        @test r.status==status
        @test size(r.kernel)==(n,n)
        @test verify_integer_weight_solution(r)
    end
    @test_throws ArgumentError solve_integer_weight_system([1.0 2.0],[1])
    @test_throws ArgumentError solve_integer_weight_system([1 2],[1.0])
    @test_throws DimensionMismatch solve_integer_weight_system([1 2],[1,2])
end

@testset "certificate tampering does not masquerade as saturation" begin
    r=solve_integer_weight_system([2 1 1],[0])
    # Same kernel equation and rational rank, but index four in the full lattice.
    badV=copy(r.certificate.V); badV[:,2:end].*=2
    bad=merge(r,(kernel=2 .* r.kernel,certificate=merge(r.certificate,(V=badV,))))
    @test bad.A*bad.kernel==zeros(BigInt,1,2)
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    badU=2 .* r.certificate.U
    bad=merge(r,(certificate=merge(r.certificate,(U=badU,S=2 .* r.certificate.S,
                         transformed_rhs=2 .* r.certificate.transformed_rhs)),))
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    bad=deepcopy(r);bad.certificate.S[1,1]+=1
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    bad=merge(r,(particular=BigInt[1,0,0],))
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    bad=merge(r,(status=:rational_only,))
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    @test_throws ArgumentError verify_integer_weight_solution(merge(r,(kernel=Float64.(r.kernel),)))
    r=solve_integer_weight_system(reshape([2],1,1),[1])
    bad=merge(r,(obstruction=merge(r.obstruction,(remainder=big(0),)),))
    @test_throws ArgumentError verify_integer_weight_solution(bad)
end

@testset "current presentation disk equations and initial assignment" begin
    for s in 2:4
        v=lens_space_link(s,1)
        prior=_integer_solver_raw(v)
        r=solve_integer_weights(v)
        @test _integer_solver_raw(v)==prior
        @test r.presentation.closed_normal
        @test r.A==2 .* BigInt.(disk_table(v).delta)
        @test r.b==-BigInt.(disk_table(v).cp2)
        @test r.status==(s==2 ? :integer_solution : :rational_only)
        @test verify_integer_weight_solution(r)
        assigned=solve_integer_weights(v;assign=true)
        if s==2
            @test assigned.assignment_status==:assigned
            @test v.is_weighted
            @test all(iszero,weight_residual(v))
            @test reduce(vcat,v.weights)==Float64.(assigned.particular)
        else
            @test assigned.assignment_status==:not_integer_feasible
            @test _integer_solver_raw(v)==prior
            @test assigned.rational_particular!==nothing
        end
        # An explicit 02 geometric relation changes this presentation's matrix.
        # Recheck the new matrix; never reinterpret a failed solve as a claim
        # about all manifold presentations. No weighted transport is requested.
        plain=lens_space_link(s,1)
        first_pair=plain.gauss[1][1:2]
        changed=move(plain,"02";param="dir",vertices=first_pair)
        changed_result=solve_integer_weights(changed)
        @test changed_result.A!=r.A
        @test changed_result.status==r.status
        @test changed_result.presentation.closed_normal
    end
    for s in 2:4
        arr=collect(1:2s-1)
        orientation=Int.((1 .- (-1).^arr)./2)
        signed=arr.*(-1).^arr
        v=VirtualLink();set_data!(v;gauss=[vcat(-signed,reverse(signed))],orientation=orientation)
        r=solve_integer_weights(v)
        @test r.status==(s==2 ? :integer_solution : :rational_only)
        @test !r.presentation.closed_normal
        @test size(r.A,2)==length(real_edge_table(v))==2s
        @test verify_integer_weight_solution(r)
    end
    v=VirtualLink();set_data!(v;gauss=[[1,-1,2,-3],[-2,3]],orientation=[1,1,0])
    before=_integer_solver_raw(v)
    r=solve_integer_weights(v;assign=true)
    @test r.assignment_status==:assigned
    @test v.gauss==before[1] && v.orientation==before[2] && v.real_only==before[4]
    @test [length(w) for w in v.weights]==[3,1]
    @test all(iszero,weight_residual(v))

    v=lens_space_link(2,1);set_weight!(v,[[9,8,7,6]])
    r=solve_integer_weights(v)
    huge=merge(r,(particular=r.particular+big(2)^100 .* r.kernel[:,1],))
    @test verify_integer_weight_solution(huge)
    prior=_integer_solver_raw(v)
    @test Topology._assign_integer_weight_particular!(v,huge)==:outside_safe_integer_range
    @test _integer_solver_raw(v)==prior
    @test_throws ArgumentError Topology._assign_integer_weight_particular!(v,
        solve_integer_weight_system(reshape([1],1,1),[1]))
    @test _integer_solver_raw(v)==prior
    cut=deepcopy(v);push!(cut.gauss[1],0)
    @test_throws ArgumentError solve_integer_weights(cut;assign=true)
    @test_throws ArgumentError solve_integer_weights(v;assign=1)
    impossible=lens_space_link(3,1);set_weight!(impossible,[fill(0.5,6)])
    prior=_integer_solver_raw(impossible)
    @test solve_integer_weights(impossible;assign=true).assignment_status==:not_integer_feasible
    @test _integer_solver_raw(impossible)==prior
    malformed=VirtualLink();set_data!(malformed;rgauss=[[1,1]],orientation=[1])
    prior=_integer_solver_raw(malformed)
    @test_throws ArgumentError solve_integer_weights(malformed;assign=true)
    @test _integer_solver_raw(malformed)==prior
end

@testset "specified parity can obstruct an otherwise solvable weight equation" begin
    v=lens_space_link(2,1)
    ordinary=solve_integer_weights(v)
    @test ordinary.status==:integer_solution
    z=BigInt[1,0,0,0]
    disks=disk_table(v)
    @test all(iseven,disks.cp2)
    # This parity satisfies the mod-2 cochain equation. The obstruction is to
    # an integer lift, not just failure of a mod-2 precondition.
    @test all(iszero,mod.(BigInt.(disks.delta)*z+div.(BigInt.(disks.cp2),2),2))
    prior=_integer_solver_raw(v)
    blocked=solve_integer_weights(v;parity=z,assign=true)
    @test blocked.status==:rational_only
    @test blocked.assignment_status==:not_integer_feasible
    @test blocked.particular===nothing
    @test blocked.rational_particular!==nothing
    @test _integer_solver_raw(v)==prior
    @test verify_integer_weight_solution(blocked)
    @test blocked.t_solution.A==4 .* BigInt.(disks.delta)
    @test blocked.t_solution.b==-BigInt.(disks.cp2)-2 .* BigInt.(disks.delta)*z

    z=mod.(ordinary.particular,2)
    lifted=solve_integer_weight_lift(v,z;assign=true)
    @test lifted.status==:integer_solution
    @test lifted.assignment_status==:assigned
    @test mod.(lifted.particular,2)==z
    @test lifted.kernel==2 .* lifted.t_solution.kernel
    @test all(iseven,lifted.kernel)
    @test all(iszero,weight_residual(v))
    @test verify_integer_weight_solution(lifted)
    bad=deepcopy(lifted);bad.parity[1]=1-bad.parity[1]
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    bad=merge(lifted,(kernel=lifted.t_solution.kernel,))
    @test_throws ArgumentError verify_integer_weight_solution(bad)
    @test_throws DimensionMismatch solve_integer_weight_lift(v,[0])
    @test_throws ArgumentError solve_integer_weight_lift(v,[0,0,0,2])
    @test_throws ArgumentError solve_integer_weight_lift(v,[0.0,0.0,0.0,0.0])
    empty_lift=solve_integer_weight_lift(zeros(Int,0,2),Int[],[1,0])
    @test empty_lift.particular==[1,0]
    @test empty_lift.kernel==2 .* Matrix{BigInt}(I,2,2)
end
