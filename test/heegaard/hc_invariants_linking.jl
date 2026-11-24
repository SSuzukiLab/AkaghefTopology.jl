# cf: .research-scratch/PJ260922/astra-lanes/theory/REPORT-2.md §独立比較 — MOO (1.1), CFH Proposition 1.2, registered forward inputs; primary Work registration pending, unverified
using Test, Topology, LinearAlgebra
import AkaghefAlgebra

const HCL_H = Topology.Heegaard
const HCL_TE = AkaghefAlgebra.TaftExact
const HCL_ROOT = normpath(joinpath(@__DIR__,"..","..","..",".."))
const HCL_INPUTS = [
    ("L(2,1)","threads/PJ260922/data/261112-example-table/A2-1/transport.hc",0,[2]),
    ("L(5,2)","threads/PJ260922/data/261112-example-table/A2-2/transport.hc",0,[5]),
    ("L(7,3)","threads/PJ260922/data/261112-example-table/A2-3/transport.hc",0,[7]),
    ("L(3,1)","threads/PJ260922/data/261112-example-table/E1-summand-L3/turn6-candidate/diagram.hc",0,[3]),
    ("Weeks",".research-scratch/261111-original-annotations-rebuild/PJ260922_source_transcriptions/cases/Weeks_fig1_inner_m/Weeks_fig1_inner_m.hc",0,[5,5]),
    ("T3",".research-scratch/261111-original-annotations-rebuild/PJ260922_source_transcriptions/cases/T3_fig3_curve_layer/T3_fig3_curve_layer.hc",3,Int[]),
    ("S3",".research-scratch/PJ260922/astra-lanes/hc-moves/turn3-snappy-baseline/input-001-all_family_twist-1-back.hc",0,Int[]),
    ("S2xS1","threads/PJ260922/data/261112-example-table/A3-registered/transport.hc",1,Int[]),
]

function hcl_embed(x,c)
    x isa Union{Integer,Rational} && return HCL_TE.Cyclotomic(c,x)
    a = x isa AkaghefAlgebra.Cyclotomic ? AkaghefAlgebra.canonical_value(x) : x
    n = Int(typeof(a).parameters[1])
    c%n==0 || throw(ArgumentError("incompatible embedding"))
    z = HCL_TE.zeta(c)
    sum((co*z^((j-1)*(c÷n)) for (j,co) in enumerate(HCL_TE.coefficients(a)));init=zero(z))
end
hcl_order(x::AkaghefAlgebra.Cyclotomic)=AkaghefAlgebra.root_spec(x).N
hcl_order(x::HCL_TE.Cyclotomic)=Int(typeof(x).parameters[1])
hcl_order(x::Union{Integer,Rational})=1
function hcl_equal(x,y)
    c = lcm(2,hcl_order(x),hcl_order(y))
    hcl_embed(x,c)==hcl_embed(y,c)
end
function hcl_conjugate(x)
    a = AkaghefAlgebra.canonical_value(x)
    c = hcl_order(x); z = HCL_TE.zeta(c)
    sum((co*z^mod(-(j-1),c) for (j,co) in enumerate(HCL_TE.coefficients(a)));init=zero(z))
end

# Independent original surgery sum. All four matrices below are positive
# definite, so signature=size(A,1) and normalization is G_N^(-size(A,1)).
function hcl_surgery_moo(A,N)
    n = size(A,1); z=HCL_TE.zeta(2N)
    G = sum((z^(k*k) for k in 0:N-1);init=zero(z))
    S = sum((z^Int(mod(dot(collect(x),A*collect(x)),2N))
        for x in Iterators.product(ntuple(_ -> 0:N-1,n)...));init=zero(z))
    S/(G^n)
end

function hcl_reverse_surface(hc)
    l=hc.lc
    lc=HCL_H.LickorishCode(l.genus,Dict(c=>-w for (c,w) in l.words),
        copy(l.places),deepcopy(l.disks),copy(l.points),copy(l.weights),l.source_path)
    HCL_H.HeegaardCode(lc,copy(hc.roles),deepcopy(hc.disks),copy(hc.points),copy(hc.weights),hc.hc_path)
end

@testset "HC linking and even MOO" begin
    forms=Dict{String,HCL_H.HCLinkingForm}()
    inputs=Dict{String,HCL_H.HeegaardCode}()
    @testset "frozen identified L(5,2), available without local corpus" begin
        # Exact LC and roles from identified forward row 368, A2-2/transport.*.
        # Joining metadata is copied; the invariant never chooses a disk.
        words=Dict(
            "l1"=>[-9,-11,-13,-1,-10,-12],
            "m1"=>[14,-15,1,17,16],
            "alpha1"=>[-7,6,4,2,5,3,8,15],
            "beta1"=>[-8,10,-2,11,-16,-3,12,-4,13,-17,-5,9,-6,7,-14])
        lc=HCL_H.LickorishCode(1,words,Dict{String,String}(),
            Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
        hc=HCL_H.resolve_hc(lc,Dict("alpha1"=>:alpha,"beta1"=>:beta);
            disks=Dict("joining"=>["F13"]))
        F=HCL_H.hc_linking_form(hc)
        forms["L(5,2)"]=F;inputs["L(5,2)"]=hc
        @test F.free_rank==0 && F.torsion==[5]
        @test HCL_H.Z_MOO(hc)==-1
        # Evaluate the named m₁ class, independently of the SNF's generator.
        coordinates=F.certificate.left*BigInt[0,1]
        indices=[j for j in 1:2 if F.certificate.smith[j,j]>1]
        @test HCL_H.linking_pairing(F,coordinates[indices],coordinates[indices])==big(2)//5
    end
    @testset "registered identified HC inputs" begin
        for (name,path,b1,torsion) in HCL_INPUTS
            file=joinpath(HCL_ROOT,path)
            if !isfile(file)
                @test_skip isfile(file) # Local acceptance corpus is not shipped with the package.
                continue
            end
            hc=HCL_H.read_hc(file); F=HCL_H.hc_linking_form(hc)
            forms[name]=F; inputs[name]=hc; c=F.certificate
            @test F.free_rank==b1 && F.torsion==torsion
            @test c.left*c.presentation*c.right==c.smith
            @test c.alpha_parts+c.beta_parts==F.generators*Diagonal(F.torsion)
            @test F.matrix==transpose(F.matrix)
            @test all(x -> 0<=x<1,F.matrix)
            @test all(i -> all(j -> denominator(F.torsion[i]*F.matrix[i,j])==1,
                eachindex(F.torsion)),eachindex(F.torsion))
            opposite=HCL_H.hc_linking_form(hc;boundary_side=:beta)
            @test opposite.matrix==mod.(-F.matrix,1)
        end
    end

    @testset "lens forms up to unoriented generator isometry" begin
        for (name,p,q) in (("L(2,1)",2,1),("L(5,2)",5,2),("L(7,3)",7,3),("L(3,1)",3,1))
            haskey(forms,name) || continue
            F=forms[name]
            @test length(F.torsion)==1
            @test any(mod(F.matrix[1,1]-s*(big(q)//p)*u^2,1)==0
                for s in (-1,1),u in 1:p-1 if gcd(u,p)==1)
            @test HCL_H.linking_pairing(F,[p],[1])==0
        end
    end

    @testset "original MOO surgery normalization and sign" begin
        i=HCL_TE.zeta(4)
        # Locks the sign independently of the HC and linking implementation.
        @test hcl_equal(hcl_surgery_moo(reshape([3],1,1),2),-i)
        @test hcl_equal(hcl_surgery_moo([3 1;1 2],2),-1)
        @test hcl_equal(hcl_surgery_moo([3 1 0;1 2 1;0 1 2],2),i)
        matrices=Dict("L(2,1)"=>reshape([2],1,1),"L(5,2)"=>[3 1;1 2],
            "L(7,3)"=>[3 1 0;1 2 1;0 1 2],"L(3,1)"=>reshape([3],1,1))
        for (name,A) in matrices
            haskey(forms,name) || continue
            # The source table identifies unoriented manifolds. One consistent
            # orientation (direct or conjugate) must work at every level.
            pairs=[(HCL_H.Z_MOO(forms[name],N),hcl_surgery_moo(A,N)) for N in (2,4,6,8,10,12)]
            direct=all(hcl_equal(x,y) for (x,y) in pairs)
            reverse=all(hcl_equal(hcl_conjugate(x),y) for (x,y) in pairs)
            @test direct || reverse
        end
        if haskey(forms,"L(5,2)")
            @test HCL_H.Z_MOO(inputs["L(5,2)"],2)==-1
            @test HCL_H.Z_MOO(forms["L(5,2)"],4)==1
            @test HCL_H.Z_MOO(forms["L(5,2)"],6)==1
        end
    end

    @testset "free homology and Weeks" begin
        z8=HCL_TE.zeta(8); sqrt2=z8+z8^7
        z24=HCL_TE.zeta(24); sqrt6=(z24^3+z24^21)*(z24^2+z24^22)
        roots=Dict(2=>sqrt2,4=>2,6=>sqrt6)
        for (name,b1) in (("S3",0),("S2xS1",1),("T3",3))
            haskey(forms,name) || continue
            for N in (2,4,6)
                expected=b1==0 ? 1 : b1==1 ? roots[N] : N*roots[N]
                @test hcl_equal(HCL_H.Z_MOO(forms[name],N),expected)
            end
        end
        if haskey(forms,"Weeks")
            @test HCL_H.Z_MOO(forms["Weeks"],10)==5
        end
        if haskey(forms,"S2xS1")
            @test HCL_H.Z_MOO(forms["S2xS1"],big(2)^200)==big(2)^100
            @test_throws ArgumentError HCL_H.Z_MOO(forms["S2xS1"],2_000_000_014;max_conductor=32)
        end
    end

    @testset "orientation reversal" begin
        if haskey(inputs,"L(3,1)")
            hc=inputs["L(3,1)"]; F=forms["L(3,1)"]
            reversed=HCL_H.hc_linking_form(hcl_reverse_surface(hc))
            @test reversed.torsion==F.torsion
            @test reversed.generators==F.generators
            @test reversed.matrix==mod.(-F.matrix,1)
            for N in (2,4,6)
                @test hcl_equal(HCL_H.Z_MOO(reversed,N),hcl_conjugate(HCL_H.Z_MOO(F,N)))
            end
        end
        if haskey(forms,"L(3,1)") && haskey(forms,"L(7,3)")
            # Both stored orientations give i. Equality of public invariant
            # values must survive their different torsion presentations.
            @test HCL_H.Z_MOO(forms["L(3,1)"],2)==HCL_H.Z_MOO(forms["L(7,3)"],2)
        end
    end

    @testset "mixed free and torsion, and two-primary forms" begin
        J=BigInt[0 0 1 0;0 0 0 1;-1 0 0 0;0 -1 0 0]
        A=BigInt[0 0;0 0;1 0;0 1]
        # Block sum of the [3] lens presentation and a free S²×S¹ summand.
        mixed=HCL_H._hcl_lattice_form(A,BigInt[3 0;0 0;-1 0;0 1],J)
        @test mixed.free_rank==1 && mixed.torsion==[3]
        z8=HCL_TE.zeta(8); sqrt2=z8+z8^7
        @test hcl_equal(HCL_H.Z_MOO(mixed,2),-z8^2*sqrt2)
        @test hcl_equal(HCL_H.Z_MOO(mixed,6),3sqrt2)
        hyperbolic=HCL_H._hcl_lattice_form(A,BigInt[2 0;0 2;0 -1;-1 0],J)
        diagonal=HCL_H._hcl_lattice_form(A,BigInt[2 0;0 2;-1 0;0 -1],J)
        @test hyperbolic.torsion==diagonal.torsion==[2,2]
        @test HCL_H.Z_MOO(hyperbolic,2)==2
        @test HCL_H.Z_MOO(diagonal,2)==0
        @test HCL_H.Z_MOO(diagonal,4)==2
        # Both columns are individually primitive; their full lattice is not.
        @test_throws HCL_H.HCLinkingInputError HCL_H._hcl_lattice_form(
            BigInt[1 1;1 -1;0 0;0 0],A,J)
        if haskey(inputs,"Weeks")
            hc=inputs["Weeks"]; roles=copy(hc.roles)
            roles["alpha1"]=:beta; roles["beta1"]=:alpha
            bad=HCL_H.HeegaardCode(hc.lc,roles,deepcopy(hc.disks),copy(hc.points),copy(hc.weights),hc.hc_path)
            err=try HCL_H.hc_linking_form(bad);nothing catch e;e end
            @test err isa HCL_H.HCLinkingInputError && err.reason==:family_intersection
        end
    end

    @testset "empty torsion, standard aliases, and input rejection" begin
        J=BigInt[0 1;-1 0]
        F=HCL_H._hcl_lattice_form(reshape(BigInt[0,1],2,1),reshape(BigInt[3,-1],2,1),J)
        @test F.torsion==[3] && F.matrix==reshape([big(2)//3],1,1)
        @test hcl_equal(HCL_H.Z_MOO(F,2),-HCL_TE.zeta(4))
        @test_throws HCL_H.HCLinkingInputError HCL_H._hcl_lattice_form(
            reshape(BigInt[0,2],2,1),reshape(BigInt[3,-1],2,1),J)
        @test_throws ArgumentError HCL_H.Z_MOO(F,1)
        @test_throws ArgumentError HCL_H.Z_MOO(F,0)
        @test_throws ArgumentError HCL_H.Z_MOO(F,-2)
        @test_throws ArgumentError HCL_H.Z_MOO(F,true)
        @test_throws ArgumentError HCL_H.Z_MOO(F,2;max_terms=2)
        @test_throws ArgumentError HCL_H.Z_MOO(F,2;max_conductor=2)
        @test_throws DimensionMismatch HCL_H.linking_pairing(F,Int[],Int[])
        @test_throws ArgumentError HCL_H._hcl_lattice_form(reshape(BigInt[0,1],2,1),
            reshape(BigInt[3,-1],2,1),J;boundary_side=:unknown)
        # Genus zero is the empty complete diagram for S³.
        lc=HCL_H.LickorishCode(0,Dict{String,Vector{Int}}(),Dict{String,String}(),
            Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
        hc=HCL_H.HeegaardCode(lc,Dict{String,Symbol}(),Dict{String,Vector{String}}(),
            Dict{String,String}(),Dict{String,Int}(),nothing)
        @test HCL_H.Z_MOO(hc,2)==1
        @test HCL_H.Z_MOO(hc,4)==HCL_H.Z_MOO(forms["L(5,2)"],4)
        if haskey(inputs,"L(2,1)")
            # Standard can bind a role to a reference K curve itself.
            source=inputs["L(2,1)"]
            removed=Set(abs.(source.lc.words["alpha1"]))
            words=Dict(c=>filter(x -> !(abs(x) in removed),w)
                for (c,w) in source.lc.words if c!="alpha1")
            reduced=HCL_H.LickorishCode(1,words,Dict{String,String}(),
                Dict{String,Vector{String}}(),Dict{String,String}(),Dict{String,Int}(),nothing)
            standard=HCL_H.resolve_hc(reduced,Dict("beta1"=>:beta);standard=true)
            @test isempty(HCL_H.hc_linking_form(standard).torsion)
            other=HCL_H.resolve_hc(reduced,Dict("beta1"=>:beta,"l1"=>:alpha))
            @test HCL_H.hc_linking_form(other).torsion==[2]
        end
    end
end
