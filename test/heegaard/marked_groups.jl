module MarkedGroupTests
using Test, Topology, JSON3
const H=Topology.Heegaard
const B=H.MarkedDiagramBuilders
loadcase(name)=JSON3.read(read(joinpath(@__DIR__,"../../vendor/marked_surface_diagrams/examples",name*".json"),String),Dict{String,Any})
@testset "marked cellular groups and curve words" begin
    for g in 1:4
        d=H.lickorish_reference(g);p=H.surface_group_presentation(d)
        @test length(p.generators)==2g
        @test length(p.relators)==1
        @test all(iszero,H.abelianization_matrix(p))
        @test H.first_homology(p).free_rank==2g
        @test length(H.marked_curve_words(d).curve_words)==3g-1
    end
    for (surface,rank) in ((:disk,0),(:annulus,1))
        p=H.surface_group_presentation(B.grid_document(4;surface))
        @test isempty(p.relators)
        @test H.first_homology(p).free_rank==rank
    end
    @test H.first_homology(loadcase("sphere_empty")).free_rank==0
    @test H.first_homology(loadcase("torus_mer_lon")).free_rank==0
    parallel=loadcase("torus_parallel_zero_crossings")
    @test H.first_homology(parallel).free_rank==1
    @test isempty(H.first_homology(parallel).torsion)
    @test H.intersection_matrix(parallel)==zeros(Int,1,1)
    @test H.intersection_matrix(loadcase("torus_mer_lon"))==[1;;]
    # A repeated meridian is an extra relation, never an extra handle generator.
    B.add_curve!(parallel,"a2","alpha",true,[[(0,2),(4,2)]])
    @test H.first_homology(parallel).free_rank==1
    @test !H.compression_boundary_report(parallel).alpha_complete
end

# Source: marked_groups.jl::marked_curve_words fixed-source edge interface and
# vendor builders.jl::waypoint_word. The diagonal has independently specified
# positive horizontal/vertical winding (1,1); no moved target or recomputed
# target presentation is used. Disk triangle identities test the boundary case.
@testset "fixed source edge-groupoid words" begin
    M=H.MarkedSurfaceDiagrams
    torus=loadcase("torus_mer_lon");st=M.validate(torus);data=H.marked_curve_words(torus)
    function edge_image(word)
        letters=Int[]
        for d in M.parse_word(word)
            i=first(st.incidence[d.face][M.edge(d.u,d.v)])[1]
            append!(letters,[g*s for (g,s) in data.edge_words[(d.face,i-1,d.u,d.v)]])
        end
        [(abs(x),sign(x)) for x in H._marked_reduce(letters)]
    end
    function exponents(word)
        result=zeros(Int,length(data.presentation.generators))
        for (g,s) in word;result[g]+=s;end
        result
    end
    a=only(c for c in torus["curves"] if c["role"]=="alpha")
    b=only(c for c in torus["curves"] if c["role"]=="beta")
    av=exponents(edge_image(a["word"]));bv=exponents(edge_image(b["word"]))
    @test abs(av[1]*bv[2]-av[2]*bv[1])==1
    diagonal=B.waypoint_word([[(0,1),(3,4)],[(3,0),(4,1)]])
    @test exponents(edge_image(diagonal))==av+bv
    @test all(data.edge_words[(f,k,v,u)]==[(g,-s) for (g,s) in reverse(w)]
        for ((f,k,u,v),w) in data.edge_words)
    @test all(data.edge_words[(s["source_triangle"]...,s["oriented_vertices"]...)]==[(s["generator"],1)]
        for s in data.presentation.metadata["generator_edges"])
    @test length(data.edge_correspondence)==6sum(length(p["triangles"]) for p in torus["patches"])
    @test count(r->r["retained_relator"]!==nothing,data.face_correspondence)==1
    retained=only(filter(r->r["retained_relator"]!==nothing,data.face_correspondence))
    @test retained["word"]==only(data.presentation.relators)
    @test sort([g*s for (g,s) in retained["word"]])==[-2,-1,1,2]
    @test all(r["retained_relator"]!==nothing || isempty(r["word"]) for r in data.face_correspondence)
    @test all(data.face_words[Tuple(r["source_triangle"])]==r["word"] for r in data.face_correspondence)
    edges=Dict(s["edge"]=>s for s in data.presentation.metadata["cellular_edge_sources"])
    function based_path_valid(path)
        current=data.presentation.metadata["base_vertex"]
        for e in path["signed_cellular_edges"]
            source=edges[abs(e)]
            start=e>0 ? source["start_vertex"] : source["end_vertex"]
            stop=e>0 ? source["end_vertex"] : source["start_vertex"]
            start==current || return false
            current=stop
        end
        current==path["vertex"]
    end
    @test all(based_path_valid,data.presentation.metadata["vertex_tree_paths"])
    disk=H.marked_curve_words(B.grid_document(3;surface=:disk))
    @test isempty(disk.presentation.generators) && all(isempty,values(disk.edge_words))
    @test all(isempty,values(disk.face_words)) && all(r["retained_relator"]===nothing for r in disk.face_correspondence)
    # Tuple-key maps stay Julia-side; these source records are the JSON surface.
    wire=JSON3.read(JSON3.write(data.edge_correspondence),Vector{Any})
    @test length(wire)==length(data.edge_words)
end

# Source: marked_groups.jl::marked_computation_report and the existing vendor
# fixtures/builders. Expected torus and boundary H1 values are independent
# topological examples; source maps are checked back against saved elements.
@testset "marked computation source and applicability" begin
    M=H.MarkedSurfaceDiagrams
    torus=loadcase("torus_mer_lon")
    torus["metadata"]=Dict("source"=>Dict("figure"=>"test locator"))
    before=deepcopy(torus);r=H.marked_computation_report(torus)
    @test r["diagram_valid"] && r["surface"]["applicable"] && r["manifold"]["applicable"]
    @test torus==before
    @test r["source"]["representation_sha256"]==M.representation_digest(torus)
    @test r["source"]["marking_sha256"]==M.marking_digest(torus)
    @test r["source"]["manifest"]==H.source_manifest(torus)
    r["source"]["manifest"]["source"]["figure"]="changed report"
    @test torus["metadata"]==before["metadata"]
    inter=r["intersections"]
    @test inter["matrix"]==[1;;]
    @test inter["row_curves"]==[c["id"] for c in torus["curves"] if c["role"]=="alpha"]
    @test inter["column_curves"]==[c["id"] for c in torus["curves"] if c["role"]=="beta"]
    @test only(inter["contributions"])["source_crossing"]==only(torus["crossings"])
    @test only(inter["contributions"])["contribution"]==1
    region_triangles=[Tuple(t) for region in r["regions"] for t in region["triangles"]]
    source_triangles=[(patch["face"],i-1) for patch in torus["patches"] for i in eachindex(patch["triangles"])]
    @test sort(region_triangles)==sort(source_triangles)
    reversed_order=deepcopy(torus)
    reverse!(reversed_order["crossings"][1]["curves"])
    reversed_order["crossings"][1]["sign"]*=-1
    @test H.marked_computation_report(reversed_order)["intersections"]["matrix"]==[1;;]
    sp=r["surface"]["presentation"];mp=r["manifold"]["presentation"]
    @test r["surface"]["first_homology"]["free_rank"]==2
    @test r["manifold"]["first_homology"]["free_rank"]==0
    @test !r["manifold"]["first_homology"]["element_coordinates_available"]
    @test length(mp["metadata"]["relator_sources"])==length(mp["relators"])
    for source in mp["metadata"]["relator_sources"]
        if source["kind"]=="attached_curve"
            @test mp["relators"][source["relator"]]==r["curves"][source["curve_id"]]["word"]
        end
    end
    for source in sp["metadata"]["generator_edges"]
        face,index=source["source_triangle"]
        patch=only(filter(p->p["face"]==face,torus["patches"]))
        @test all(v->v in patch["triangles"][index+1],source["oriented_vertices"])
    end
    @test all(e->any(s->s["edge"]==e,sp["metadata"]["cellular_edge_sources"]),sp["metadata"]["primal_tree_edges"])
    @test H.surface_group_presentation(torus).metadata["source_hash"]==M.representation_digest(torus)

    # Auxiliary reference curves remain source elements/words, never attachments.
    reference=H.lickorish_reference(2);rr=H.marked_computation_report(reference)
    @test Set(keys(rr["curves"]))==Set(c["id"] for c in reference["curves"])
    @test rr["source"]["representation_sha256"]==M.representation_digest(reference)
    @test Set(c["role"] for c in values(rr["curves"]))==Set(c["role"] for c in reference["curves"])
    @test rr["surface"]["first_homology"]["free_rank"]==4

    parallel=loadcase("torus_parallel_zero_crossings")
    B.add_curve!(parallel,"a2","alpha",true,[[(0,2),(4,2)]])
    pr=H.marked_computation_report(parallel)
    @test pr["manifold"]["applicable"] && pr["manifold"]["first_homology"]["free_rank"]==1
    @test size(pr["intersections"]["matrix"])==(2,1)
    @test pr["intersections"]["row_curves"]==sort([c["id"] for c in parallel["curves"] if c["role"]=="alpha"])
    @test isempty(pr["intersections"]["contributions"])

    for (surface,rank) in ((:disk,0),(:annulus,1))
        boundary=H.marked_computation_report(B.grid_document(4;surface))
        @test boundary["surface"]["first_homology"]["free_rank"]==rank
        @test !boundary["manifold"]["applicable"] && !isempty(boundary["manifold"]["reasons"])
        @test !haskey(boundary["manifold"],"presentation")
    end
    arc=B.grid_document(4;surface=:disk)
    B.add_curve!(arc,"a","alpha",false,[[(2,0),(2,4)]])
    ar=H.marked_computation_report(arc);path=ar["curves"]["a"]
    @test ar["surface"]["applicable"] && !ar["manifold"]["applicable"]
    @test path["kind"]=="unbased_proper_path" && !path["closed"]
    @test path["endpoint_anchors"]==M.endpoint_anchors(M.validate(arc),"a")
    @test length(path["endpoint_anchors"])==2
    @test !haskey(ar["manifold"],"first_homology")

    incomplete=B.grid_document(6)
    B.add_curve!(incomplete,"a","alpha",true,[[(1,1),(5,1),(5,5),(1,5),(1,1)]])
    ir=H.marked_computation_report(incomplete)
    @test ir["surface"]["applicable"] && !ir["manifold"]["applicable"]
    bad=deepcopy(torus);bad["crossings"][1]["sign"]*=-1
    invalid=H.marked_computation_report(bad)
    @test !invalid["diagram_valid"] && !invalid["surface"]["applicable"] && !invalid["manifold"]["applicable"]
    @test only(invalid["issues"])["code"]=="CROSSING_SIGN"
    # The report is consumable by the app's JSON path, including empty matrices.
    wire=JSON3.read(JSON3.write(ar),Dict{String,Any})
    @test wire["curves"]["a"]["kind"]=="unbased_proper_path"
    @test wire["intersections"]["matrix_shape"]==[1,0]
    @test wire["intersections"]["matrix_rows"]==[[]]
end
end # module MarkedGroupTests
