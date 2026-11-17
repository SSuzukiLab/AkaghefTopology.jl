using Test, Topology
isdefined(@__MODULE__,:vl_example) || include(joinpath(@__DIR__,"../examples/VLExample.jl"))

# cf:lib/topology.jl/src/hopf_invariant.jl @ af90b62, _hopf_invariant_factors — sparse reference; working-tree split: src/hopf_link.jl.
function group_sparse_reference(v,H)
    M=parentmodule(typeof(invariant_expression(v,H)))
    M._group_sparse_reference(v,H)
end

@testset "integer symbolic group MST" begin
    H=cyclic_group_algebra(5)
    @test all(c->c isa BigInt,values(H.multiplication.terms))
    M=parentmodule(typeof(invariant_expression(vl_example("lens_2_1"),H)))
    E=M.GroupConstraintExpression
    # Exhaustive oracle for small generic ternary networks, including repeated
    # variables, negative vertices, noncommutativity and inconsistent domains.
    table=[1 2 3 4 5 6;2 3 1 6 4 5;3 1 2 5 6 4;4 5 6 1 2 3;5 6 4 3 1 2;6 4 5 2 3 1]
    for K in (cyclic_group_algebra(2),cyclic_group_algebra(3),group_algebra(table))
        data=M._group_data(K); d=K.dimension
        for constraints in ([(1,2,3,1)],[(1,2,3,-1)],[(1,2,1,1)],[(1,1,2,1)],
                            [(1,1,1,-1)],[(1,2,3,1),(2,1,3,-1)],[(1,2,3,1),(2,3,1,1)])
            n=maximum(maximum(c[1:3]) for c in constraints)
            domains=[trues(d) for _ in 1:n]
            for limited in (false,true)
                limited && (domains[1][2:end].=false)
                expected=BigInt(0)
                for assignment in Iterators.product((1:d for _ in 1:n)...)
                    all(domains[i][assignment[i]] for i in 1:n) || continue
                    if all(assignment[c]==data[1][assignment[a],s==1 ? assignment[b] : data[2][assignment[b]]]
                           for (a,b,c,s) in constraints)
                        expected+=1
                    end
                end
                expr=E(data...,constraints,deepcopy(domains),true)
                before=deepcopy(expr.domains)
                @test evaluate_invariant(expr)==expected
                @test evaluate_invariant(simplify_invariant(expr))==expected
                @test evaluate_invariant(simplify_invariant(expr;strategy=:none))==expected
                @test expr.domains==before
            end
        end
    end
    noncomm=M._group_data(group_algebra(table))
    for (a,b,c,sign,expected) in ((2,4,6,1,1),(2,4,5,1,0),
                                (4,2,5,1,1),(4,2,5,-1,0),(4,2,6,-1,1))
        domains=[falses(6) for _ in 1:3]
        domains[1][a]=true;domains[2][b]=true;domains[3][c]=true
        @test evaluate_invariant(E(noncomm...,[(1,2,3,sign)],domains,true))==expected
    end
    data=M._group_data(H)
    free=E(data...,NTuple{4,Int}[],[trues(5) for _ in 1:100],true)
    @test evaluate_invariant(free)==big(5)^100
    @test evaluate_invariant(free;details=true).branches==0
    empty=E(data...,NTuple{4,Int}[],[falses(5)],true)
    @test evaluate_invariant(empty)==0
    impossible=E(data...,[(1,1,1,1)],[BitVector([false,true,false,false,false])],true)
    @test evaluate_invariant(impossible)==0
    @test !simplify_invariant(impossible).consistent
    for K in (cyclic_group_algebra(2),H,group_algebra(table)),
        name in ("s3_abalone","lens_2_1","lens_3_1_mst2_example23","quaternion","poincare","s2xs1_2vertex")
        v=vl_example(name)
        @test Z_MST(v,K) isa BigInt
        @test Z_MST(v,K)==group_sparse_reference(v,K)
        @test evaluate_invariant(simplify_invariant(invariant_expression(v,K)))==Z_MST(v,K)
    end
    for p in 2:5,K in (H,group_algebra(table))
        v=framed_lens_space_link(p,1)
        @test Z_MST(v,K)==group_sparse_reference(v,K)
        before=deepcopy(v.weights)
        @test Z_MST(v,K;strategy=:none)==Z_MST(v,K)
        @test v.weights==before
    end
    # Relabel the group identity; change input labels without changing its value.
    perm=[4,1,6,2,5,3]; backward=invperm(perm)
    renamed=[backward[table[perm[a],perm[b]]] for a in 1:6,b in 1:6]
    v=vl_example("quaternion")
    @test Z_MST(v,group_algebra(renamed))==Z_MST(v,group_algebra(table))
    modified=deepcopy(H); modified.positive_vertex.terms[first(keys(modified.positive_vertex.terms))]=2
    @test M._group_data(modified)===nothing
    @test !(invariant_expression(v,modified) isa E)
    @test Z_MST(v,modified)==M._group_sparse_reference(v,modified)
    @test_throws ArgumentError simplify_invariant(invariant_expression(v,H);strategy=:unknown)
    opened=vl_example("lens_2_1"); set_cuts!(opened,[true])
    @test_throws ArgumentError invariant_expression(opened,H)
end
