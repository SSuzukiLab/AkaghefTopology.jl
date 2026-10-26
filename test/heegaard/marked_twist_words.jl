module MarkedTwistWordCoreTests
using Test,Topology,JSON3
const H=Topology.Heegaard
const M=H.MarkedSurfaceDiagrams
const Q=M.Rat
const E=[[Q(i==j ? 1 : 0) for j in 1:3] for i in 1:3]
const CENTER=fill(Q(1//3),3)
rows_json(rows)=[M.qstr.(r) for r in rows]
function map_atom(source,target,s,t)
    Dict("source"=>collect(source),"target"=>collect(target),
        "source_weights"=>rows_json(s),"target_weights"=>rows_json(t))
end
function identity_atoms(mesh)
    [map_atom((f,i-1),(f,i-1),E,E) for (f,i) in sort(collect(keys(mesh.triangles)))]
end
function bump(mesh,target_center)
    result=identity_atoms(mesh);p=Tuple(first(result)["source"]);deleteat!(result,1)
    for (a,b) in ((1,2),(2,3),(3,1))
        push!(result,map_atom(p,p,[E[a],E[b],CENTER],[E[a],E[b],target_center]))
    end
    result,p
end

# Source: direct affine evaluation on one source triangle. Expected values in
# the test below are independently derived from 3/4*center + 1/4*corner.
function evaluate(atoms,mesh,parent,point)
    results=Set()
    for a in H._mtw_atoms(atoms,mesh)
        a.source==parent || continue
        triangle=[(w[2],w[3]) for w in a.source_weights];x,y,z=triangle
        den=M.orient(x,y,z)
        w=[M.orient(point,y,z)/den,M.orient(x,point,z)/den,M.orient(x,y,point)/den]
        all(>=(0),w) || continue
        target=[sum(w[j]*a.target_weights[j][k] for j in 1:3) for k in 1:3]
        push!(results,(a.target,Tuple(target)))
    end
    only(results)
end

@testset "Private exact RIGHT composition on one fixed prepared mesh" begin
    doc=H.MarkedDiagramBuilders.grid_document(2);mesh=H._marked_states(doc).geometry
    F,parent=bump(mesh,Q[1//4,1//2,1//4])
    G,_=bump(mesh,Q[1//4,1//4,1//2])
    snapshot=deepcopy((F,G,doc))
    FG=H._mtw_compose_maps(F,G,mesh)
    GF=H._mtw_compose_maps(G,F,mesh)
    @test (F,G,doc)==snapshot
    # F sends center to (1/2,1/4), G to (1/4,1/2). The latter is
    # 3/4*center + 1/4*C, so F(G(center))=(3/8,7/16), and vice versa.
    fg=evaluate(FG.atoms,mesh,parent,(Q(1//3),Q(1//3)))
    gf=evaluate(GF.atoms,mesh,parent,(Q(1//3),Q(1//3)))
    @test fg==(parent,(Q(3//16),Q(3//8),Q(7//16)))
    @test gf==(parent,(Q(3//16),Q(7//16),Q(3//8)))
    @test fg!=gf
    @test !isempty(FG.correspondence["pieces"])
    @test all(begin
        li=p["left_atom"];ri=p["right_atom"]
        p["middle_parent"]==F[li]["source"]==G[ri]["target"]
    end for p in FG.correspondence["pieces"])
    @test all(r["area_ratio"]=="1" for r in FG.correspondence["result"]["source_partition"])
    @test all(r["area_ratio"]=="1" for r in FG.correspondence["result"]["target_partition"])
    @test FG.correspondence["result"]["forward_continuity"]["paired_edges"]>0
    @test FG.correspondence["result"]["inverse_continuity"]["paired_edges"]>0
    serialized=JSON3.read(JSON3.write(FG.atoms),Vector{Any})
    @test H._mtw_certify(serialized,mesh).correspondence==FG.correspondence["result"]

    # The inverse here swaps the ACTUAL paired map rows; this is a private
    # affine-map check and is not a public native word inverse claim.
    invF=[Dict("source"=>copy(a["target"]),"target"=>copy(a["source"]),
        "source_weights"=>deepcopy(a["target_weights"]),"target_weights"=>deepcopy(a["source_weights"])) for a in F]
    identity=H._mtw_compose_maps(F,invF,mesh)
    @test evaluate(identity.atoms,mesh,parent,(Q(1//3),Q(1//3)))==(parent,(Q(1//3),Q(1//3),Q(1//3)))
    @test evaluate(identity.atoms,mesh,parent,(Q(1//5),Q(2//5)))==(parent,(Q(2//5),Q(1//5),Q(2//5)))

    overlap=[deepcopy(F);[deepcopy(first(F))]]
    @test_throws ArgumentError H._mtw_certify(overlap,mesh)
    missing=deepcopy(F);pop!(missing)
    @test_throws ArgumentError H._mtw_certify(missing,mesh)
    folded=deepcopy(F);folded[end]["target_weights"]=reverse(folded[end]["target_weights"])
    @test_throws ArgumentError H._mtw_certify(folded,mesh)
    # Full source and target partitions can each be valid while their affine
    # pairing is discontinuous. A corner cycle in only one old triangle must
    # fail across the actual old edges. Parallel carriers are checked below
    # independently of this vertex-image failure.
    broken=identity_atoms(mesh);broken[1]["target_weights"]=rows_json(E[[2,3,1]])
    @test_throws ArgumentError H._mtw_certify(broken,mesh)
    # A one-sided old-edge midpoint is a T junction even when the two local
    # area sums are correct. It must not be silently inserted by the validator.
    nonconforming=identity_atoms(mesh);p=Tuple(first(nonconforming)["source"]);deleteat!(nonconforming,1)
    mid=(E[1]+E[2])/2
    push!(nonconforming,map_atom(p,p,[E[1],mid,E[3]],[E[1],mid,E[3]]))
    push!(nonconforming,map_atom(p,p,[mid,E[2],E[3]],[mid,E[2],E[3]]))
    @test_throws ArgumentError H._mtw_certify(nonconforming,mesh)

    # On this torus several distinct fine edges have the SAME quotient
    # endpoints. Endpoint roots alone would identify them incorrectly.
    pairs=M.fine_edge_occurrences(mesh);parallel=Dict()
    for (o,p) in pairs
        isless(o,p) || continue
        endpoints=minmax(M.root!(mesh.vertex_dsu,(o[1],o[3])),M.root!(mesh.vertex_dsu,(o[1],o[4])))
        push!(get!(parallel,endpoints,Any[]),o)
    end
    candidates=sort([sort(es) for es in values(parallel) if length(es)>=2];by=first)
    @test !isempty(candidates)
    o,p=first(candidates)[1:2]
    function carrier(occ)
        t=mesh.triangles[occ[1:2]]
        u=only(findall(==(occ[3]),t));v=only(findall(==(occ[4]),t))
        first(H._mtw_carrier((occ[1],occ[2]-1),E[u],E[v],mesh,pairs))
    end
    @test carrier(o)!=carrier(p)
    @test minmax(M.root!(mesh.vertex_dsu,(o[1],o[3])),M.root!(mesh.vertex_dsu,(o[1],o[4])))==
        minmax(M.root!(mesh.vertex_dsu,(p[1],p[3])),M.root!(mesh.vertex_dsu,(p[1],p[4])))
end

@testset "Paired polygon rows and prestep RIGHT transvections" begin
    points=[(Q(0),Q(0)),(Q(1//2),Q(0)),(Q(1//2),Q(1//2)),(Q(0),Q(1//2))]
    source=[[1-x-y,x,y] for (x,y) in points]
    target=[[1-(x/2+1//8)-(y/2+1//16),x/2+1//8,y/2+1//16] for (x,y) in points]
    triangles=H._mtw_paired_triangles(points,source,target)
    @test length(triangles)==2
    @test all(begin
        r.source_weights==source[r.indices] && r.target_weights==target[r.indices] &&
        all(t[2]==s[2]/2+1//8 && t[3]==s[3]/2+1//16 for (s,t) in zip(r.source_weights,r.target_weights))
    end for r in triangles)
    @test sum(H._mt_ratio(r.source_weights) for r in triangles)==1//2
    @test sum(H._mt_ratio(r.target_weights) for r in triangles)==1//8

    original=Dict("a"=>[1,0],"b"=>[0,1]);snapshot=deepcopy(original)
    pairing=Dict(("a","a")=>0,("a","b")=>1,("b","a")=>-1,("b","b")=>0)
    A=H._mtw_right_update(original,"a",1,pairing)
    AB=H._mtw_right_update(A,"b",1,pairing)
    B=H._mtw_right_update(original,"b",1,pairing)
    BA=H._mtw_right_update(B,"a",1,pairing)
    @test original==snapshot
    @test hcat(A["a"],A["b"])==[1 1;0 1]
    @test hcat(B["a"],B["b"])==[1 0;-1 1]
    @test hcat(AB["a"],AB["b"])==[0 1;-1 1]
    @test hcat(BA["a"],BA["b"])==[1 1;-1 0]
    @test AB!=BA
    @test H._mtw_right_update(A,"a",-1,pairing)==original
    huge=Dict("a"=>[typemax(Int),0],"b"=>[1,1])
    exact=H._mtw_right_update(huge,"a",1,pairing)
    @test exact["b"]==[BigInt(typemax(Int))+1,BigInt(1)]
    @test eltype(exact["b"])==BigInt
    large_count=copy(pairing);large_count[("a","b")]=typemax(Int)
    @test H._mtw_right_update(huge,"a",1,large_count)["b"][1]==BigInt(typemax(Int))^2+1
    @test_throws ArgumentError H._mtw_right_update(original,"a",true,pairing)
    absent=copy(pairing);delete!(absent,("a","b"))
    @test_throws ArgumentError H._mtw_right_update(original,"a",1,absent)
end

@testset "Private original-core products share one prepared mesh pair" begin
    doc=H.MarkedDiagramBuilders.grid_document(2)
    H.MarkedDiagramBuilders.add_curve!(doc,"a","alpha",true,[[(0,1),(2,1)]])
    H.MarkedDiagramBuilders.add_curve!(doc,"b","beta",true,[[(1,0),(1,2)]])
    H.MarkedDiagramBuilders.declare_crossing!(doc,"x",1,1,"a","b",1)
    requested=["a","b","a"];snapshot=deepcopy((doc,requested))
    prepared=H._mtw_prepare_original_products(doc,requested)
    @test (doc,requested)==snapshot
    @test Set(keys(prepared.products))==Set(["a","b"])
    @test prepared.correspondence["requested_cores"]==["a","b","a"]
    @test prepared.correspondence["product_cores"]==["a","b"]
    @test prepared.products["a"]===prepared.preparation.product
    @test prepared.supports["a"]==prepared.preparation.support
    @test prepared.supports["a"]!=prepared.supports["b"]
    @test length(prepared.fine.triangles)==6length(prepared.coarse.triangles)
    @test length(prepared.preparation.history)==2
    @test Set(keys(prepared.fine.curves))==Set(["a","b"])
    @test prepared.correspondence["coarse_source_hash"]==H._mm_hash(prepared.preparation.coarse)
    @test prepared.correspondence["fine_source_hash"]==H._mm_hash(prepared.preparation.fine)
    @test prepared.basis.presentation.metadata["source_hash"]==M.representation_digest(prepared.preparation.fine)
    @test prepared.basis.presentation.metadata["source_marking_sha256"]==prepared.correspondence["fine_marking_hash"]
    for cid in ("a","b")
        product=prepared.products[cid]
        @test product.support_cells==prepared.supports[cid]
        @test Set(Tuple(row.id) for row in product.triangles)==prepared.supports[cid]
        @test all(haskey(prepared.fine.triangles,(row.id[1],row.id[2]+1)) for row in product.triangles)
        @test all(M.orient(row.points...)>0 for row in product.triangles)
        @test sum(M.orient(row.points...)/2 for row in product.triangles)==2
        @test product.core_order==[[d.face,d.u,d.v] for d in prepared.coarse.curves[cid][2]]
        @test product.certificate["core"]==cid && product.certificate["core_winding"]=="1"
        @test length(prepared.fine.curves[cid][2])==2length(prepared.coarse.curves[cid][2])
    end
    @test_throws ArgumentError H._mtw_prepare_original_products(doc,["a","missing"])
    @test_throws ArgumentError H._mtw_prepare_original_products(doc,String[])
    @test_throws ArgumentError H._mtw_prepare_original_products(doc,["a",17])
end


# DRAFT ONLY: these helpers require the shared extraction patch. They do not
# run against the currently frozen application source and are not acceptance
# of an ordered native word or its specified-factor inverse.
function total_candidate_spec(prepared,result)
    fine=prepared.preparation.fine;raw=result.raw
    spec=Dict{String,Any}("target_hash"=>H._mm_hash(result.document),
        "fine_source_mesh"=>Dict("surface"=>deepcopy(fine["surface"]),"patches"=>deepcopy(fine["patches"])),
        "fine_source_hash"=>H._mm_hash(fine),"child_cells"=>raw.records,
        "new_vertices"=>raw.vertex_records,"curve_images"=>raw.image_records,
        "region_seeds"=>H._mm_inverse_seeds(fine),
        "proof_crossing_locations"=>Dict(x["id"]=>x["at"] for x in raw.complex["crossings"]),
        "proof_region_seeds"=>H._mm_inverse_seeds(raw.complex),
        "rechart_inverse_specification"=>result.rechart===nothing ? nothing : result.rechart.inverse_specification)
    H._mm_infinity(fine)!==nothing && (spec["infinity_triangle"]=deepcopy(H._mm_infinity(fine)))
    H._mm_infinity(raw.complex)!==nothing && (spec["proof_infinity_triangle"]=deepcopy(H._mm_infinity(raw.complex)))
    spec
end

@testset "DRAFT private total-map assembly and actual prepared-source pullback" begin
    # Source: the unchanged embedded-D input in marked_twists.jl. This small
    # boundary-fixed bump is a DIFFERENT private map test; passing it cannot
    # replace the required whole-Dehn-action D constructor/inverse gate.
    doc=H.MarkedDiagramBuilders.grid_document(4)
    H.MarkedDiagramBuilders.add_curve!(doc,"a","alpha",true,[[(1,1),(3,1),(3,3),(1,3),(1,1)]])
    H.MarkedDiagramBuilders.add_curve!(doc,"b","beta",true,[[(0,2),(4,2)]])
    H.MarkedDiagramBuilders.add_curve!(doc,"cut","aux",true,[[(2,0),(2,4)]])
    H.MarkedDiagramBuilders.declare_crossing!(doc,"ab_right",3,2,"a","b",-1)
    H.MarkedDiagramBuilders.declare_crossing!(doc,"ab_left",1,2,"a","b",1)
    H.MarkedDiagramBuilders.declare_crossing!(doc,"ac_bottom",2,1,"a","cut",1)
    H.MarkedDiagramBuilders.declare_crossing!(doc,"ac_top",2,3,"a","cut",-1)
    H.MarkedDiagramBuilders.declare_crossing!(doc,"bc",2,2,"b","cut",1)
    doc["region_marks"]=[Dict("id"=>"D","seed"=>["P",10],"claim"=>"embedded_closed_disk")]
    doc["metadata"]=Dict("display"=>Dict("infinity_triangle"=>["P",0]))
    prepared=H._mtw_prepare_original_products(doc,["cut","a"])
    atoms,parent=bump(prepared.fine,Q[1//4,1//2,1//4])
    snapshot=deepcopy((doc,atoms))
    actual=H._mtw_realize_total_map(prepared,atoms)
    result=actual.result
    @test (doc,atoms)==snapshot
    @test H.validate(result.document;level=:surface).valid
    @test result.rechart===nothing
    @test length(result.raw.vertex_records)==1
    @test any(r["source_weights"]!=r["target_weights"] for r in result.records)
    @test result.document["curves"]==prepared.preparation.fine["curves"]
    @test result.document["crossings"]==prepared.preparation.fine["crossings"]
    @test result.document["region_marks"][1]["claim"]=="embedded_closed_disk"
    # The map changes an actual interior point, although this boundary-fixed
    # bump intentionally leaves the saved edge curves fixed. It is not a
    # substitute for the later noncommuting Dehn-word acceptance case.
    center=evaluate(atoms,prepared.fine,parent,(Q(1//3),Q(1//3)))
    @test center==(parent,(Q(1//4),Q(1//2),Q(1//4)))
    spec=JSON3.read(JSON3.write(total_candidate_spec(prepared,result)),Dict{String,Any})
    restored=H._mtw_restore_total_source_candidate(result.document,spec)
    @test restored.candidate==prepared.preparation.fine
    @test restored.forward_word_regeneration_required
    @test length(restored.regions)==length(result.region_map)
    @test length(restored.geometric_regions)==length(result.geometry_regions)
    wrong=deepcopy(spec);wrong["child_cells"][1]["source_weights"][1]=["0","0","0"]
    @test_throws ArgumentError H._mtw_restore_total_source_candidate(result.document,wrong)
    target=deepcopy(result.document);target["metadata"]["changed"]=true
    @test_throws ArgumentError H._mtw_restore_total_source_candidate(target,spec)
end
end
