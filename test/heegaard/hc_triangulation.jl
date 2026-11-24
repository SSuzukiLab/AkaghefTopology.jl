# cf: lib/topology.jl/docs/specs/heegaard-snappy.md — finite quotient and optional backend contract
using Test
using Topology
const HC_TRI_H = Topology.Heegaard

# cf: parent turn 10 — genus-one sums and fresh geometric handle cancellation.
@testset "HC surface operations" begin
    root=normpath(joinpath(@__DIR__,"..","..","..",".."))
    data=joinpath(root,"threads","PJ260922","data","261112-example-table")
    for basis in (:l,:m)
        lens=HC_TRI_H.read_hc(joinpath(data,"A2-1","turn9-public-standard-"*String(basis),"diagram.hc"))
        sum=HC_TRI_H.hc_connected_sum(lens,lens)
        @test sum.lc.genus==2
        @test count(==(:alpha),values(sum.hc.roles))==2
        @test count(==(:beta),values(sum.hc.roles))==2
        @test isempty(sum.hc.disks)
        @test sum.certificate["status"]=="joining-disk-pending"
        hom=HC_TRI_H.hc_first_homology(sum.hc)
        @test (hom.free_rank,hom.torsion)==(0,[2,2])
        @test HC_TRI_H.hc_hom_count(sum.hc,:S3)==16
        forgot=HC_TRI_H.forget_reference(HC_TRI_H.realize_hc(sum.hc))
        @test count(r->(r.chi,r.genus,r.boundary_components)==(0,0,2),forgot.regions)==1
        @test isempty([r for r in forgot.regions if r.joining_disk])

        l3=HC_TRI_H.read_hc(joinpath(data,"E1-summand-L3","turn9-public-standard-"*String(basis),"diagram.hc"))
        @test_throws HC_TRI_H.LCMoveError HC_TRI_H.hc_connected_sum(lens,l3)
        @test_throws HC_TRI_H.LCMoveError HC_TRI_H.hc_destabilize(l3;alpha=1,beta=1,alpha_basis=basis,disk_name="joining")
        reduced=HC_TRI_H.hc_destabilize(l3;alpha=1,beta=2,alpha_basis=basis,disk_name="joining")
        @test reduced.lc.genus==1
        hom=HC_TRI_H.hc_first_homology(reduced.hc)
        @test (hom.free_rank,hom.torsion)==(0,[3])
        @test HC_TRI_H.hc_hom_count(reduced.hc,:S3)==3
        tri=HC_TRI_H.hc_triangulation(reduced.hc;disk_name="joining",signature=false)
        @test tri.topology.vertices==1
        @test tri.topology.euler_characteristic==0

        torus=HC_TRI_H.read_hc(joinpath(data,"F1-2","turn9-public-standard-"*String(basis),"diagram.hc"))
        first=HC_TRI_H.hc_destabilize(torus;alpha=2,beta=1,alpha_basis=basis,disk_name="joining")
        second=HC_TRI_H.hc_destabilize(first.hc;alpha=3,beta=3,alpha_basis=basis,disk_name="joining")
        @test (first.lc.genus,second.lc.genus)==(4,3)
        hom=HC_TRI_H.hc_first_homology(second.hc)
        @test (hom.free_rank,hom.torsion)==(3,Int[])
        @test HC_TRI_H.hc_hom_count(second.hc,:S3)==48
        tri=HC_TRI_H.hc_triangulation(second.hc;disk_name="joining",signature=false)
        @test tri.topology.vertices==1
        @test tri.topology.euler_characteristic==0
    end
end

@testset "HC closed triangulation" begin
    # Registered o-graph inputs, not reconstructed examples.
    registry = Module(:HCTriangulationRegistry)
    Base.include(registry, abspath(joinpath(@__DIR__, "..", "..", "examples", "VLExample.jl")))
    for name in (:s3_abalone, :lens_2_1, :lens_3_1_mst2_example23, :poincare, :s2xs1_2vertex)
        v = registry.vl_example(String(name))
        t = HC_TRI_H.hc_triangulation(v; signature=false)
        @test t.tetrahedra == length(Topology.real_gauss_code(v)[2])
        @test t.topology.vertices == 1
        @test t.topology.euler_characteristic == 0
        @test t.isosig === nothing
        @test length(t.gluings) == 2t.tetrahedra
        faces = vcat([[(g.tet,g.omitted),(g.to_tet,g.to_omitted)] for g in t.gluings]...)
        @test length(unique(faces)) == 4t.tetrahedra
        reverse = HC_TRI_H.hc_triangulation_to_odata(t; orders=t.branch_orders)
        @test reverse.correspondence["exact_oriented_face_correspondence"]
        @test reverse.value.orientation == v.orientation
        @test reverse.correspondence["input_local_edge_to_disk"] == t.correspondence["tetrahedron_edges_to_disks"]
        @test HC_TRI_H.hc_triangulation(reverse.virtual_link;signature=false).topology == t.topology
        found = HC_TRI_H.hc_find_branching(t)
        @test length(found.orders) == t.tetrahedra
        @test HC_TRI_H.hc_find_branching(t;orders=found.orders).orders == found.orders
        @test_throws HC_TRI_H.HCTriangulationError HC_TRI_H.hc_find_branching(t;
            orders=[zeros(Int,4) for _ in 1:t.tetrahedra])
        @test_throws HC_TRI_H.HCTriangulationError HC_TRI_H._hc_triangulation_topology(
            t.tetrahedra, t.gluings[1:end-1])
        @test_throws HC_TRI_H.HCTriangulationError HC_TRI_H._hc_triangulation_topology(
            t.tetrahedra, vcat(t.gluings,t.gluings[1:1]))
    end
end

@testset "HC inverse certificates and explicit unresolved boundary" begin
    registry=Module(:HCInverseRegistry)
    Base.include(registry,abspath(joinpath(@__DIR__,"..","..","examples","VLExample.jl")))
    for name in ("s3_abalone","lens_2_1","lens_3_1_mst2_example23","poincare","s2xs1_2vertex")
        tri=HC_TRI_H.hc_triangulation(registry.vl_example(name);signature=false)
        enumeration=HC_TRI_H.hc_each_branching(tri)
        @test enumeration.complete
        @test !isempty(enumeration.solutions)
        @test all(o->HC_TRI_H.hc_find_branching(tri;orders=o).orders==o,enumeration.solutions)
        limited=HC_TRI_H.hc_each_branching(tri;max_solutions=1)
        @test !limited.complete
        @test limited.stop==:branching_solution_budget
        @test length(limited.solutions)==1
        f=HC_TRI_H.hc_flow_certificate(tri;orders=tri.branch_orders)
        @test f.certificate["boundary_chi"]==2
        @test f.certificate["bw_chi"]==[1,1]
        @test f.certificate["koda_roundtrip"]["equal_up_to_metadata"]
        cut=HC_TRI_H.hc_surface_cut_system(f.diagram)
        @test length(cut["residual_edges"])==2*cut["genus"]
        @test length(cut["polygon_word"])==4*cut["genus"]
        @test cut["pl_identification"]===nothing
        r=HC_TRI_H.hc_prepare_inverse(tri;max_depth=0,max_states=1)
        @test r["status"]=="koda_only"
        @test r["hc"]===nothing
        @test r["stop"]=="polygon_normalization_and_reference_overlay_required"
    end
end

const HC_CUT_TEST_PAIRS = []
@testset "HC genus-one surface transport" begin
    for (p,q) in ((2,1),(5,2),(7,3))
        source=HC_TRI_H.hc_triangulation(Topology.lens_space_link(p,q);signature=false)
        flow=HC_TRI_H.hc_flow_certificate(source;orders=source.branch_orders)
        cut=HC_TRI_H.hc_surface_cut_system(flow.diagram)
        result=HC_TRI_H.hc_cut_to_resolved(flow.diagram,cut;joining_region=flow.region)
        @test result.hc isa Topology.HeegaardDiagram
        @test HC_TRI_H.validate_diagram(result.hc)
        @test HC_TRI_H._lc_surface_stats(result.lc.words) |> x->x[1]-x[2]+x[3]==0
        @test result.certificate["joining_region_source"]==flow.region
        @test all(c["preserved"] for c in result.certificate["vertex_rotation_checks"])
        forgot=HC_TRI_H.forget_reference(HC_TRI_H.realize_hc(result.hc))
        @test Set(only(r for r in forgot.regions if r.id==result.joining_region_id).face_roots)==Set(result.joining_faces)
        h1=HC_TRI_H.hc_first_homology(result.hc)
        @test h1.free_rank==0 && h1.torsion==[p]
        mktempdir() do dir
            HC_TRI_H.write_lc(joinpath(dir,"diagram.lc"),result.lc)
            HC_TRI_H.write_hc(joinpath(dir,"diagram.hc"),result.hc;lc_filename="diagram.lc")
            reread=HC_TRI_H.read_hc(joinpath(dir,"diagram.hc"))
            @test reread.lc.words==result.lc.words
            @test reread.disks==result.hc.disks
            target=HC_TRI_H.hc_triangulation(reread;joining_disk=result.joining_region_id,signature=false)
            @test target.topology.vertices==1
            @test target.topology.euler_characteristic==0
            used=vcat([[(f.tet,f.omitted),(f.to_tet,f.to_omitted)] for f in target.gluings]...)
            @test length(used)==length(unique(used))==4target.tetrahedra
            push!(HC_CUT_TEST_PAIRS,(source,target,reread))
        end
        bad=deepcopy(cut);reverse!(bad["disk_boundary_darts"])
        @test_throws HC_TRI_H.HCTriangulationError HC_TRI_H.hc_cut_to_resolved(flow.diagram,bad;joining_region=flow.region)
    end
end

@testset "HC standard genus-two reference overlay" begin
    registry=Module(:HCStandardOctagonRegistry)
    Base.include(registry,abspath(joinpath(@__DIR__,"..","..","examples","VLExample.jl")))
    source=HC_TRI_H.hc_triangulation(registry.vl_example("lens_3_1_mst2_example23");signature=false)
    flow=HC_TRI_H.hc_flow_certificate(source;orders=source.branch_orders)
    cut=HC_TRI_H.hc_surface_cut_system(flow.diagram)
    @test cut["genus"]==2
    result=HC_TRI_H.hc_cut_to_resolved(flow.diagram,cut;joining_region=flow.region)
    @test result.hc isa Topology.HeegaardDiagram
    @test HC_TRI_H.validate_diagram(result.hc)
    v,e,f=HC_TRI_H._lc_surface_stats(result.lc.words)
    @test v-e+f == -2
    constants=Dict("l1"=>[-1,-3],"l2"=>[-2,4],"m1"=>[1],"m2"=>[2],"n1"=>[3,-4])
    for (name,w) in constants
        @test [x for x in result.lc.words[name] if abs(x)<=4]==w
    end
    @test length(result.certificate["reference_arcs"]["n1"])==2
    @test all(c["preserved"] for c in result.certificate["vertex_rotation_checks"])
    @test length(result.certificate["face_corner_checks"])==4length(result.certificate["old_vertex_to_lc_label"])
    @test result.certificate["joining_region_source"]==flow.region
    @test result.certificate["grid_trials"]<=result.certificate["grid_bound"]
    mktempdir() do dir
        HC_TRI_H.write_lc(joinpath(dir,"diagram.lc"),result.lc)
        HC_TRI_H.write_hc(joinpath(dir,"diagram.hc"),result.hc;lc_filename="diagram.lc")
        hc=HC_TRI_H.read_hc(joinpath(dir,"diagram.hc"))
        @test hc.lc.words==result.lc.words
        @test hc.disks==result.hc.disks
        forgot=HC_TRI_H.forget_reference(HC_TRI_H.realize_hc(hc))
        selected=only(r for r in forgot.regions if r.id==result.joining_region_id)
        @test Set(selected.face_roots)==Set(result.joining_faces)
        @test (selected.chi,selected.genus,selected.boundary_components)==(1,0,1)
        @test length(selected.boundary_ports)==4
        @test selected.joining_disk
        real=forgot.realization
        allfaces=Set(f.id for f in real.faces);allnames=Set(keys(hc.lc.words))
        closed=HC_TRI_H._region_cell_data(hc.lc.words,hc.roles,allfaces,real,allnames)
        @test closed[1:3]==(-2,2,0)
        @test isempty(closed[4])
        onecut=HC_TRI_H._region_cell_data(hc.lc.words,hc.roles,allfaces,real,setdiff(allnames,Set(["alpha1"])))
        @test onecut[1:3]==(-2,1,2)
        bothcut=HC_TRI_H._region_cell_data(hc.lc.words,hc.roles,allfaces,real,setdiff(allnames,Set(["alpha1","alpha2"])))
        @test bothcut[1:3]==(-2,0,4)
        forward=HC_TRI_H.hc_triangulation(hc;joining_disk=result.joining_region_id,signature=false)
        @test forward.topology.euler_characteristic==0
        @test forward.topology.vertices==1
        used=vcat([[(f.tet,f.omitted),(f.to_tet,f.to_omitted)] for f in forward.gluings]...)
        @test length(used)==length(unique(used))==4forward.tetrahedra
        push!(HC_CUT_TEST_PAIRS,(source,forward,hc))
    end
end

const HC_STANDARD_TEST_PAIRS = []
@testset "HC alpha-adapted standard public inverse" begin
    registry=Module(:HCStandardRegistry)
    Base.include(registry,abspath(joinpath(@__DIR__,"..","..","examples","VLExample.jl")))
    for name in ("lens_2_1","lens_3_1_mst2_example23")
        source=HC_TRI_H.hc_triangulation(registry.vl_example(name);signature=false)
        orders=source.branch_orders
        for basis in (:l,:m)
            result=HC_TRI_H.hc_triangulation_to_standard(source;alpha_basis=basis,orders=orders)
            @test HC_TRI_H.validate_diagram(result.hc)
            @test Set(HC_TRI_H.hc_curve_names(result.hc,:alpha))==Set("$(basis)$i" for i in 1:result.lc.genus)
            @test all(c["preserved"] for c in result.certificate["source_curve_checks"])
            @test result.certificate["reference_words_checked"]
            forgotten=HC_TRI_H.forget_reference(HC_TRI_H.realize_hc(result.hc))
            @test all((r.chi,r.genus,r.boundary_components)==(1,0,1) for r in forgotten.regions)
            selected=only(r for r in forgotten.regions if r.id==result.joining_region_id)
            @test Set(selected.face_roots)==Set(result.joining_faces)
            @test selected.joining_disk
            target=HC_TRI_H.hc_triangulation(result.hc;disk_name="joining",signature=false)
            @test target.topology.vertices==1
            @test target.topology.euler_characteristic==0
            @test length(unique(vcat([[(f.tet,f.omitted),(f.to_tet,f.to_omitted)] for f in target.gluings]...)))==4target.tetrahedra
            push!(HC_STANDARD_TEST_PAIRS,(source,target))
        end
        @test_throws HC_TRI_H.HCTriangulationError HC_TRI_H.hc_triangulation_to_standard(source;alpha_basis=:invalid,orders=orders)
    end
end


@testset "HC optional SnapPy" begin
    runner = try
        HC_TRI_H.hc_snappy_runner()
    catch e
        e isa HC_TRI_H.SnapPyUnavailable || rethrow()
        nothing
    end
    available = runner !== nothing && success(pipeline(
        `$runner -c "import snappy; from snappy.snap.t3mlite import Mcomplex"`, stdout=devnull, stderr=devnull))
    if !available
        @info "SKIP SnapPy checks: uv and an offline SnapPy 3.3.2 environment are required"
        @test_skip available
    else
        for (source,target,hc) in HC_CUT_TEST_PAIRS
            comparison=HC_TRI_H.hc_snappy_compare(source,target;runner=runner)
            @test comparison["same_manifold"] === true
            a=HC_TRI_H.hc_snappy(source;runner=runner)
            b=HC_TRI_H.hc_snappy(target;runner=runner)
            @test a["homology"]["value"]==b["homology"]["value"]
            for group in (:Z2,:Z3,:Z5,:S3)
                expected=string(HC_TRI_H.hc_hom_count(hc,group))
                @test a["hom_counts"][String(group)]["value"]==b["hom_counts"][String(group)]["value"]==expected
            end
        end
        root = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
        hc = HC_TRI_H.read_hc(joinpath(root,"threads","PJ260922","data",
            "261107-lickorish-data","group-18","case-01.hc"))
        forgot = HC_TRI_H.forget_reference(HC_TRI_H.realize_hc(hc))
        candidates = [r.id for r in forgot.regions if r.joining_disk]
        @test !isempty(candidates)
        for disk in candidates
            result = HC_TRI_H.hc_snappy(hc; joining_disk=disk,
                references=["m003(-3,1)"], runner=runner)
            @test result["all_cusps_filled"]
            h1 = HC_TRI_H.hc_first_homology(hc)
            @test result["homology"]["value"]["free_rank"] == h1.free_rank
            @test result["homology"]["value"]["torsion"] == h1.torsion
            for group in (:Z2,:Z3,:Z5,:S3)
                @test result["hom_counts"][String(group)]["value"] == string(HC_TRI_H.hc_hom_count(hc,group))
            end
            @test result["references"]["m003(-3,1)"]["value"] === true
            @test !isempty(result["isosig"])
        end
        reverse_result = HC_TRI_H.hc_snappy_catalog("m003(-3,1)"; runner=runner)
        @test reverse_result["topology"]["euler_characteristic"] == 0
        @test all(==(0),reverse_result["topology"]["vertex_link_genera"])
        @test !isempty(reverse_result["dual_spine"])
        @test reverse_result["branching"]["status"] == "not_computed"
        neighbors=HC_TRI_H.hc_pachner_neighbors(reverse_result["triangulation"];runner=runner)
        @test !isempty(neighbors["states"])
        @test all(s->s["topology"]["vertices"]==reverse_result["topology"]["vertices"],neighbors["states"])
        @test all(s->s["triangulation"]["tetrahedra"]==reverse_result["triangulation"]["tetrahedra"]+1,neighbors["states"])
        @test all(s->s["move"]["inverse_restores_isosig"],neighbors["states"])
        @test all(s->s["move"]["before_isosig"]==reverse_result["isosig"],neighbors["states"])
    end
end

@testset "HC standard independent SnapPy comparison" begin
    runner=try HC_TRI_H.hc_snappy_runner() catch; nothing end
    available=runner!==nothing && success(pipeline(`$runner -c "import snappy"`,stdout=devnull,stderr=devnull))
    if !available
        @info "SKIP standard SnapPy checks: uv and offline SnapPy 3.3.2 are required"
        @test_skip available
    else
        for (source,target) in HC_STANDARD_TEST_PAIRS
            comparison=HC_TRI_H.hc_snappy_compare(source,target;runner=runner)
            @test comparison["same_manifold"]===true
            a=HC_TRI_H.hc_snappy(source;runner=runner)
            b=HC_TRI_H.hc_snappy(target;runner=runner)
            @test a["homology"]["value"]==b["homology"]["value"]
            for group in ("Z2","Z3","Z5","S3")
                @test a["hom_counts"][group]==b["hom_counts"][group]
            end
        end
    end
end
