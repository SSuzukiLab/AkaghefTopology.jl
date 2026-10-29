# Source: sol61_build_periodic_g2_UNRUN.jl's actual native region/chain check;
# MarkedDisplayFinitePoints.compile_finite_base/old_point/point_stage and production
# build_rendered_mesh. Include inside Topology.Heegaard after marked_display_finite_points.jl.
# This intentionally rejects cap/support points; no general analytic map added.
function _md_has_periodic(input)
    input isa AbstractDict || return false
    haskey(input,"periodic_annotations") && return true
    haskey(get(input,"metadata",Dict()),"periodic_domain") && return true
    haskey(input,"source_document") && return _md_has_periodic(input["source_document"])
    false
end

function _md_periodic_preflight(input)
    _md_has_periodic(input) || return nothing
    document=get(input,"source_document",input)
    decoration=get(get(document,"metadata",Dict()),"periodic_domain",nothing)
    _mm_need(decoration isa AbstractDict,"Periodic render data requires native periodic_domain metadata")
    bound=deepcopy(document)
    pop!(bound["metadata"]["periodic_domain"],"source_binding_sha256",nothing)
    _mm_need(get(decoration,"source_binding_sha256",nothing)==_mm_hash(bound),"stale annotation source binding")
    nothing
end

# Source: native strictly interior barycentric points have unique parent cells;
# rational spelling is not identity. Called after each point's typed validation.
_md_periodic_points_distinct(a,b)=_mm_cell(a["triangle"])!=_mm_cell(b["triangle"]) ||
    _mt_live_bary(a["barycentric"])!=_mt_live_bary(b["barycentric"])

# Source: actual kernel regions plus effective recipe infinity; this helper is
# shared by compiler and narrow negative controls, not a caller waiver API.
function _md_periodic_require_bounded_annuli(regions,owners,chains,weights,infinity)
    for (i,c) in enumerate(chains)
        r=regions[owners[i]]
        if c["coefficient"]!=0
            _mm_need(r["genus"]==0 && r["boundary_components"]==2,"periodic annotations: nonzero region is not a genus-zero annulus")
        end
    end
    _mm_need(get(weights,_mm_cell(infinity),nothing)==0,"periodic annotations: nonzero region contains effective display infinity")
    nothing
end

function _md_periodic_annotations(recipe,document)
    P=MarkedDisplayFinitePoints; M=MarkedSurfaceDiagrams
    need(ok,msg)=_mm_need(ok,"periodic annotations: "*msg)
    decoration=get(get(document,"metadata",Dict()),"periodic_domain",nothing)
    decoration===nothing && return nothing
    need(get(document,"reference_context",nothing)===nothing && get(get(document,"metadata",Dict()),"reference_context",nothing)===nothing,"combined reference and annotations unsupported")
    need(decoration["schema"]=="scratch.periodic-domain/1","unknown input schema")
    # Source-bound authoring contract: exclude ONLY the hash field itself.
    # Copying old annotation metadata onto any changed source must fail.
    _md_periodic_preflight(document)
    need(recipe["source_document_sha256"]==_mm_hash(document),"wrong native source")
    st=_marked_states(document).mathematics
    regions=M.regions(st);chains=decoration["region_chains"]
    cells(r)=Set(_mm_cell(c) for c in r["triangles"])
    need(length(chains)==length(regions)==3,"requires all three actual regions")
    owners=Int[];weights=Dict{Tuple{String,Int},Int}()
    for c in chains
        n=c["coefficient"];need(n isa Integer && !(n isa Bool),"noninteger coefficient")
        hits=findall(r->cells(r)==cells(c),regions);need(length(hits)==1,"chain is not an actual full region")
        push!(owners,only(hits))
        for cell in cells(c);need(!haskey(weights,cell),"region chain overlap");weights[cell]=n;end
    end
    need(length(Set(owners))==3 && Set(keys(weights))==Set((f,k-1) for (f,k) in keys(st.triangles)),"chain omits native cells")
    # Source: compound even-odd fill represents bounded planar annuli only.
    # Check actual regions/effective infinity, not caller topology declarations.
    infinity=_mm_cell(recipe["infinity_extension"]["source_triangle"])
    _md_periodic_require_bounded_annuli(regions,owners,chains,weights,infinity)
    pair=M.fine_edge_occurrences(st);canonical(o)=min(o,pair[o])
    function term!(out,o,n)
        k=canonical(o);out[k]=get(out,k,0)+(o==k ? n : -n)
    end
    boundary=Dict();expected=Dict()
    for ((f,k),t) in st.triangles,(u,v) in M.tedges(t);term!(boundary,(f,k,u,v),weights[(f,k-1)]);end
    for (cid,n) in decoration["boundary"]
        need(n isa Integer && !(n isa Bool),"noninteger boundary")
        for d in st.curves[cid][2]
            k,u,v=first(st.incidence[d.face][M.edge(d.u,d.v)])
            term!(expected,(d.face,k,u,v),n*(u==d.u ? 1 : -1))
        end
    end
    filter!(p->last(p)!=0,boundary);filter!(p->last(p)!=0,expected)
    need(boundary==expected,"actual 2-chain boundary differs from declared curve chain")
    base=P.compile_finite_base(recipe,document;rendered_svg=true)
    # Private ordinary-call contract: fresh Python attachment geometry and
    # same-process _mdsm_verify_recipe_cells MUST precede this compiler. They
    # establish same T, S=600V, R/A boundary/topology/global-motion hypotheses.
    # No serialized checks/pass flag is consumed as authorization here.
    T=[Int.(t) for t in recipe["port_map"]["triangles"]]
    anchored=last(base.stages);R=Dict{Int,Tuple{M.Rat,M.Rat}}();A=Dict{Int,Tuple{M.Rat,M.Rat}}()
    for (t,atom) in zip(T,anchored.atoms)
        for j in 1:3
            v=t[j]
            need(!haskey(R,v) || (R[v]==atom.source[j] && A[v]==atom.target[j]),"inverse anchor mesh vertex ties disagree")
            R[v]=atom.source[j];A[v]=atom.target[j]
        end
    end
    # Exact closed AABBs only discard strict outside points. Validate every
    # target triangle first; preserve original increasing IDs and all ties.
    need(all(M.orient((A[v] for v in t)...)>0 for t in T),"inverse A triangle not positive")
    boxes=[ntuple(k->(minimum(A[v][k] for v in t),maximum(A[v][k] for v in t)),2) for t in T]
    ratstr(v)=string(numerator(v))*"/"*string(denominator(v))
    meshpayload(points)=Dict("triangles"=>T,"positions"=>Dict(string(v)=>ratstr.(collect(q)) for (v,q) in points))
    inverse_mesh=Dict("stage"=>"reverse certified anchor mesh A to R","R_sha256"=>_mm_hash(meshpayload(R)),"A_sha256"=>_mm_hash(meshpayload(A)),"triangles"=>length(T))
    region_weights_before=copy(weights)
    function image(point)
        # Julia nested assignments otherwise capture same-named enclosing locals.
        # These are point-local temporaries; outer chain/mesh data are immutable.
        local parent,w,q,mapped,exact,candidates,inverse,global_triangle,t,a,b,c,den,inverse_weights,inverse_point
        parent=_mm_cell(point["triangle"]);w=_mt_live_bary(point["barycentric"])
        need(length(w)==3 && all(>(0),w) && sum(w)==1,"point is not strict native triangle interior")
        q=try
            mapped=P.old_point(base,parent,P.xy(w));need(mapped.finite,"infinity-cap point unsupported")
            # Source: build_fixed_quotient F1 and final_graph; K0 is not part
            # of ordinary native-to-current endpoint transport.
            P.require_final_identity_attachment(base,mapped.point)
        catch err
            # Preserve exact failed source placement. Never search/relocate here.
            throw(ArgumentError("periodic source point "*JSON3.write(point)*" rejected: "*sprint(showerror,err)))
        end
        exact=[ratstr(v) for v in q]
        candidates=[i for i in eachindex(T) if all(boxes[i][k][1]<=q[k]<=boxes[i][k][2] for k in 1:2)]
        need(!isempty(candidates),"inverse point outside actual A mesh")
        inverse=_mdsm_map_point(q,A,R,T[candidates]);global_triangle=candidates[inverse.triangle]
        t=T[global_triangle];a,b,c=(A[v] for v in t);den=M.orient(a,b,c)
        inverse_weights=[M.orient(q,b,c)/den,M.orient(a,q,c)/den,M.orient(a,b,q)/den]
        need(sum(inverse_weights)==1 && minimum(inverse_weights)>=0,"inverse anchor point barycentric inverse_weights invalid")
        need(ntuple(k->sum(inverse_weights[j]*A[t[j]][k] for j in 1:3),2)==q &&
             ntuple(k->sum(inverse_weights[j]*R[t[j]][k] for j in 1:3),2)==inverse.point,"inverse anchor point coordinate identity failed")
        need(weights==region_weights_before,"inverse point computation changed native region coefficients")
        inverse_point=Dict("triangle_index"=>global_triangle,"triangle_vertices"=>t,"lambda"=>ratstr.(inverse_weights),"qA"=>exact,"qR"=>ratstr.(collect(inverse.point)))
        Dict("triangle"=>point["triangle"],"barycentric"=>point["barycentric"],"point_svg_exact"=>exact,"inverse_anchor_point"=>inverse_point,"attachment"=>"native-to-final-PL exact image; final analytic support identity; emitted marker requires exact same-region proof")
    end
    points=Any[]
    need(length(decoration["basepoints"])==2 && Set(p["id"] for p in decoration["basepoints"])==Set(["w","z"]),"requires distinct named w/z")
    for p in decoration["basepoints"]
        need(weights[_mm_cell(p["triangle"])]==0,"basepoint outside zero-coefficient region")
        push!(points,merge(image(p),Dict("name"=>p["id"],"kind"=>"basepoint","coefficient"=>0)))
    end
    need(_md_periodic_points_distinct(points[1],points[2]),"coincident basepoints")
    rendered_regions=Any[]
    for (i,c) in enumerate(chains)
        r=regions[owners[i]];cycles=Any[]
        for cycle in r["boundary_cycles"]
            ids=Set(only(st.used[(f,M.edge(u,v))]) for (f,_,u,v) in cycle)
            need(length(ids)==1,"first-case boundary must be one entire native curve")
            cid=only(ids);native=st.curves[cid][2]
            need(length(cycle)==length(native),"partial curve boundary unsupported")
            need(Set((f,M.edge(u,v)) for (f,_,u,v) in cycle)==Set((d.face,M.edge(d.u,d.v)) for d in native),"boundary cycle differs from full native curve")
            push!(cycles,cid)
        end
        # Explicit native centroid label; never relocate after failure.
        label=get(c,"label_point",Dict("triangle"=>first(c["triangles"]),"barycentric"=>["1/3","1/3","1/3"]))
        need(_mm_cell(label["triangle"]) in cells(c),"label outside owning source region")
        push!(points,merge(image(label),Dict("name"=>string(c["coefficient"]),"kind"=>"coefficient","coefficient"=>c["coefficient"],"region"=>i)))
        push!(rendered_regions,Dict("coefficient"=>c["coefficient"],"boundary_curve_cycles"=>cycles,"triangles"=>c["triangles"]))
    end
    need(weights==region_weights_before,"native region coefficients changed during point compilation")
    Dict("schema"=>"marked-periodic-annotations/1","inverse_anchor_mesh"=>inverse_mesh,"source_binding_sha256"=>decoration["source_binding_sha256"],"source_document_sha256"=>_mm_hash(document),"source_representation_sha256"=>recipe["source_representation_sha256"],"regions"=>rendered_regions,"points"=>points,"boundary_checked"=>true,"operation_transport"=>"unsupported; reject annotated operations")
end

# Source: operation entry contract. Call BEFORE any operation/preparation write.
function _md_guard_periodic_operation(document)
    _md_has_periodic(document) && throw(ArgumentError("Periodic annotation transport is unsupported for operations; remove/re-author explicitly"))
    nothing
end

# Source: ordinary marked_display_plot.jl temporary proof/envelope flow. The
# native compiler runs in this same Julia process; Python freshly certifies
# graph/point transport and renders the certified endpoints.
function _md_finish_periodic(result,recipe;root,uv,dir)
    compiled=_md_periodic_annotations(recipe,result["source_document"])
    compiled===nothing && return result
    source_before=deepcopy(result["source_document"])
    display_before=deepcopy(result["display_request"])
    input=joinpath(dir,"periodic-view.json");proof=joinpath(dir,"periodic-base-proof.json")
    points=joinpath(dir,"periodic-compiled.json");output=joinpath(dir,"periodic-output.json")
    finished=joinpath(dir,"periodic-output-proof.json")
    write(input,JSON3.write(result));write(proof,JSON3.write(recipe));write(points,JSON3.write(compiled))
    script=joinpath(root,"scripts/heegaard/finish_periodic_annotations.py")
    command=`$uv run --no-project --offline --python 3.12 --with numpy==2.5.4 --with scipy --with shapely==2.1.2 --with networkx --with triangle==20250106 python $script $input $proof $points $output $finished`
    diagnostic=IOBuffer()
    process=run(pipeline(ignorestatus(addenv(command,"OPENBLAS_NUM_THREADS"=>"1")),stdout=stdout,stderr=diagnostic))
    success(process) || error("Periodic drawing rejected: "*last(String(take!(diagnostic)),3000))
    after=JSON3.read(read(output,String),Dict{String,Any})
    0<filesize(finished)<=32*1024*1024 || error("Invalid annotated source-map decoded size")
    bytes=read(finished);envelope=after["source_map"]
    length(bytes)==envelope["uncompressed_bytes"] && bytes2hex(SHA.sha256(bytes))==envelope["sha256"] || error("Annotated proof differs from envelope")
    final=JSON3.read(String(bytes),Dict{String,Any})
    after["source_document"]==source_before && after["display_request"]==display_before || error("Annotation renderer changed source or marking")
    final["source_document_sha256"]==recipe["source_document_sha256"] &&
        final["source_representation_sha256"]==recipe["source_representation_sha256"] || error("Annotated source binding changed")
    annotation=deepcopy(after["periodic_annotations"])
    transport=pop!(annotation,"point_transport",nothing)
    transport isa AbstractDict && get(transport,"schema",nothing)=="periodic-moving-point-transport/1" || error("Missing freshly regenerated moving-point transport")
    transport["source_document_sha256"]==compiled["source_document_sha256"] &&
        transport["source_representation_sha256"]==compiled["source_representation_sha256"] || error("Moving point source binding differs")
    transport["base_geometry_sha256"]==recipe["output_geometry_sha256"] || error("Moving point geometry binding differs")
    transport["literal_F1_point_track_claim"]===false && transport["source_points_relocated"]===false || error("Wrong point-map semantics")
    tracks=transport["tracks"];endpoints=transport["points"]
    length(tracks)==length(endpoints)==length(compiled["points"])==5 || error("Moving transport omits authored points")
    times=_mdsm_q.(tracks[1]["times"])
    first(times)==0 && last(times)==1 && all(diff(times).>0) || error("Moving transport time coverage incomplete")
    steps=length(times)-1;steps in (16,32,64) || error("Unknown bounded track subdivision")
    coverage=transport["coverage"]
    coverage["slabs"]==steps && coverage["point_pairs"]==10*steps && coverage["hole_segments"]==20*steps || error("Moving transport coverage differs")
    expected=deepcopy(compiled)
    for i in 1:5
        target=expected["points"][i];track=tracks[i];endpoint=endpoints[i]
        _mdsm_q.(track["times"])==times && length(track["points"])==steps+1 || error("Unsynchronized point tracks")
        all(length(q)==2 for q in track["points"]) || error("Bad track coordinate dimension")
        _mdsm_q.(first(track["points"]))==_mdsm_q.(target["inverse_anchor_point"]["qR"]) || error("Moving track changed actual inverse start")
        endpoint["point_index"]==i-1 && endpoint["point_svg_exact"]==target["point_svg_exact"] || error("Moving endpoint changed source qA")
        endpoint["point_svg_moved_exact"]==last(track["points"]) || error("Moving endpoint not consumed from certified track")
        moved=_mdsm_q.(endpoint["point_svg_moved_exact"])
        decoded=_mdsm_q.(endpoint["point_svg_decoded"])
        decoded==_mdsm_q.(Float64.(moved))==_mdsm_q.(endpoint["point_svg"]) || error("Moving endpoint Float64 decoder differs")
        for key in ("point_svg_moved_exact","point_svg_decoded","point_svg")
            target[key]=endpoint[key]
        end
        transport["inverse_anchor"]["points"][i]["barycentric"]==target["inverse_anchor_point"] || error("Inverse barycentric provenance changed")
    end
    all(_mdsm_q.(expected["points"][i]["point_svg_decoded"])!=_mdsm_q.(expected["points"][j]["point_svg_decoded"]) for i in 1:5 for j in 1:i-1) || error("Decoded marked points collide")
    annotation==expected && after["periodic_annotations"]==final["periodic_annotations"] || error("Renderer changed source/track endpoint annotation data")
    bytes2hex(SHA.sha256(Vector{UInt8}(codeunits(after["svg"]))))==final["output_svg_sha256"] || error("Annotated SVG digest differs")
    for (k,v) in recipe
        k in ("periodic_annotations","checks","output_geometry_sha256","output_svg_sha256","transport_normalization") && continue
        get(final,k,nothing)==v || error("Annotation renderer changed validated source-map component: "*k)
    end
    for (k,v) in recipe["checks"]
        k in ("numeric","periodic_annotations") && continue
        get(final["checks"],k,nothing)==v || error("Annotation renderer changed existing source-map check: "*k)
    end
    # Only the consumed annotation payload and final digest/check fields change.
    for (k,v) in result
        k in ("svg","source_map","periodic_annotations") && continue
        get(after,k,nothing)==v || error("Annotation renderer changed existing geometry field: "*k)
    end
    after
end
