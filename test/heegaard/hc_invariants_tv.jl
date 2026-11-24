# cf: Regina 7.4.1 Triangulation3.turaevViro — backend/API integration.
# cf: .research-scratch/PJ260922/astra-lanes/theory/REPORT-3.md §原典との独立比較 — TV (1992), Theorem 8.1.B (24)–(27); source Work unresolved, unverified registration.
using Test, Topology, JSON3, SHA
println("TV_TEST_SHA256 ",bytes2hex(sha256(read(@__FILE__))));flush(stdout)
import AkaghefAlgebra
const HCTV_H=Topology.Heegaard
const HCTV_TE=AkaghefAlgebra.TaftExact
const HCTV_ROOT=normpath(joinpath(@__DIR__,"..","..","..",".."))

function hctv_original(name,r,parity)
    c=isodd(r) && !parity ? r : 2r
    q=HCTV_TE.zeta(c)
    name=="S3" && return -(q-q^-1)^2/(2r)
    name=="S2xS1" && return one(q)
    name=="L(2,1)" && return (-q)^r == -one(q) ? (one(q)-q)*(one(q)-q^-1)/r : zero(q)
    k=fld(r-2,3)+1
    name=="L(3,1)" && return -(q^k-q^-k)^2/(2r)
    error("no original closed formula for $name")
end

const HCTV_CASES=[
    ("L(2,1)","threads/PJ260922/data/261112-example-table/A2-1/transport.hc"),
    ("L(3,1)","threads/PJ260922/data/261112-example-table/E1-summand-L3/turn6-candidate/diagram.hc"),
    ("S3",".research-scratch/PJ260922/astra-lanes/hc-moves/turn3-snappy-baseline/input-001-all_family_twist-1-back.hc"),
    ("S2xS1","threads/PJ260922/data/261112-example-table/A3-registered/transport.hc")]

const HCTV_STANDARD_CASES=[
    (name,"threads/PJ260922/data/261112-example-table/$folder/turn9-public-standard-$basis/diagram.hc")
    for (name,folder) in (("L(2,1)","A2-1"),("L(3,1)","E1-summand-L3"),("S3","A1"),("S2xS1","A3-registered"))
    for basis in ("l","m")]

@testset "Regina exact Turaev–Viro adapter" begin
    # Invalid parameters are rejected before launching the optional backend.
    t=HCTV_H.HCClosedTriangulation(0,NamedTuple[],Int[],nothing,nothing,Dict{String,Any}(),(;))
    @test_throws ArgumentError HCTV_H.Z_TuraevViro(t,2)
    @test_throws ArgumentError HCTV_H.Z_TuraevViro(t,3;timeout=Inf)
    @test_throws ArgumentError HCTV_H.Z_TuraevViro(t,3;runner="python")
    @test_throws ArgumentError HCTV_H.Z_TuraevViro(t,Tuple{Int,Bool}[])
    @test_throws ArgumentError HCTV_H.Z_TuraevViro([t],[(2,true)])
    available=try
        cmd=HCTV_H.hc_regina_runner()
        success(pipeline(`$cmd -c "import regina"`,stdout=devnull,stderr=devnull))
    catch
        false
    end
    if !available
        @info "Regina 7.4.1 exact tests require an already provisioned offline Python runner"
        @test_skip available
    else
        @test_throws HCTV_H.ReginaTVError HCTV_H.Z_TuraevViro(t)
        loaded=HCTV_H.HeegaardCode[];tris=HCTV_H.HCClosedTriangulation[];labels=Tuple{String,String}[]
        for (name,path) in vcat(HCTV_CASES,HCTV_STANDARD_CASES)
            full=joinpath(HCTV_ROOT,path)
            if !isfile(full)
                @test_skip isfile(full)
                continue
            end
            hc=HCTV_H.read_hc(full)
            seedfile=joinpath(HCTV_ROOT,"threads/PJ260922/data/261112-forward-sweep/epoch-2/triangulations/0019-001.json")
            seed=(name=="S3" && path==HCTV_CASES[3][2]) ? Int(JSON3.read(read(seedfile,String))["joining_disk"]) : nothing
            push!(tris,HCTV_H.hc_triangulation(hc;joining_disk=seed,signature=false))
            push!(loaded,hc);push!(labels,(name,path))
        end
        if !isempty(loaded)
            requests=[(r,p) for r in 3:6 for p in (true,false)]
            result=HCTV_H.Z_TuraevViro(tris,requests;certificate=true)
            for ((name,path),hc,values,job) in zip(labels,loaded,result.values,result.certificate["jobs"])
                for ((r,parity),value,cert) in zip(requests,values,job["results"])
                    @test AkaghefAlgebra.canonical_value(value)==hctv_original(name,r,parity)
                    @test cert["input_tetrahedra"]>0
                    @test cert["regina_distribution"]=="7.4.1"
                end
                for (i,p) in enumerate((true,false))
                    direct=HCTV_H.Z_TuraevViro_skeleton(hc;parity=p,certificate=true)
                    c=AkaghefAlgebra.root_spec(values[i]).N
                    @test HCTV_TE.Cyclotomic(c,direct.value)==AkaghefAlgebra.canonical_value(values[i])
                    @test direct.certificate.skeleton_chi==2
                end
                @test_throws ArgumentError HCTV_H.Z_TuraevViro_skeleton(hc,4)
            end
            # Scalar entry retains its own native file-argument path.
            @test AkaghefAlgebra.canonical_value(HCTV_H.Z_TuraevViro(first(loaded)))==hctv_original(first(labels)[1],3,true)
        end
        # These two original figure HCs have no stored disk mark. Reproduce the
        # explicit selections used by identified forward rows 2 and 1.
        # Independent Example3.weeks()/threeTorus() reference values were checked
        # at r=3,4,5, both parities; these rational constants do not use HC output.
        for (name,folder,index) in (("Weeks","Weeks_fig1_inner_m",2),("T3","T3_fig3_curve_layer",1),("S3","",19))
            rel=name=="S3" ? HCTV_CASES[3][2] : ".research-scratch/261111-original-annotations-rebuild/PJ260922_source_transcriptions/cases/$folder/$folder.hc"
            file=joinpath(HCTV_ROOT,rel)
            seedfile=joinpath(HCTV_ROOT,"threads/PJ260922/data/261112-forward-sweep/epoch-2/triangulations/$(lpad(string(index),4,'0'))-001.json")
            if !isfile(file) || !isfile(seedfile)
                @test_skip isfile(file) && isfile(seedfile)
                continue
            end
            hc=HCTV_H.read_hc(file)
            seed=Int(JSON3.read(read(seedfile,String))["joining_disk"])
            requests=[(r,p) for r in (3,4,5) for p in (true,false)]
            answer=HCTV_H.Z_TuraevViro(hc,requests;joining_disk=seed,certificate=true)
            recorded=[]
            for ((r,p),z,cert) in zip(requests,answer.values,answer.certificate["results"])
                c=AkaghefAlgebra.root_spec(z).N
                expected=name=="S3" ? hctv_original(name,r,p) : HCTV_TE.Cyclotomic(c,name=="Weeks" ? (r==4 ? big(1)//4 : big(1)//2) : big((r-1)^2)//1)
                @test AkaghefAlgebra.canonical_value(z)==expected
                push!(recorded,(r=r,parity=p,value=(field=c,coefficients=string.(HCTV_TE.coefficients(AkaghefAlgebra.canonical_value(z)))),certificate=cert))
            end
            direct=[]
            for (i,p) in enumerate((true,false))
                d=HCTV_H.Z_TuraevViro_skeleton(hc;parity=p,certificate=true)
                c=AkaghefAlgebra.root_spec(answer.values[i]).N
                matching=AkaghefAlgebra.canonical_value(answer.values[i])==HCTV_TE.Cyclotomic(c,d.value)
                @test matching
                push!(direct,(parity=p,value=string(d.value),regina_match=matching,certificate=d.certificate))
            end
            lf=HCTV_H.hc_linking_form(hc);moo=HCTV_H.Z_MOO(lf,2)
            receipt=(name=name,path=rel,style="general",status="ok",joining_disk=seed,
                selection_source=relpath(seedfile,HCTV_ROOT),hc_sha256=bytes2hex(sha256(read(file))),
                values=recorded,direct_skeleton=direct,public_hc_match=true,
                free_rank=lf.free_rank,torsion=string.(lf.torsion),
                moo2=(field=AkaghefAlgebra.root_spec(moo).N,coefficients=string.(HCTV_TE.coefficients(AkaghefAlgebra.canonical_value(moo)))))
            println("TV_SELECTED_HC ",JSON3.write(receipt));flush(stdout)
        end
    end
end
