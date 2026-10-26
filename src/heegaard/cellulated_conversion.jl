"""Project a validated cellular marked diagram to (H,D), retaining its marking.

No layout or example-specific corner list is consulted. Each fine region boundary
is collapsed along its saved curve edges to a ribbon boundary; all regions must
match bijectively before the existing Koda construction may consume the result.
The projection is for conversion only, never a replacement for the marked input.
"""
function cellulated_punctured(doc::AbstractDict; disk_id=nothing)
    M=MarkedSurfaceDiagrams
    states=_marked_states(doc)
    st=states.mathematics
    # Source: marked_topology.jl::_marked_conversion_applicability. Reject the
    # actual applicability failure before legacy projection's narrower syntax.
    all(c["closed"] for (c,_) in values(st.curves)) ||
        throw(ArgumentError("cellulated Koda conversion requires closed alpha/beta curves"))
    all(M.family_cut_report(st,f)["cut_system"] for f in ("alpha","beta")) ||
        throw(ArgumentError("cellulated Koda conversion requires both complete cut systems"))
    rs=M.regions(st)
    all(r["normalized_disk"] for r in rs) ||
        throw(ArgumentError("non-cellular marked conversion is not supported"))
    core=M.legacy_projection(st) # Auxiliary cuts remain in the full source.
    marks=get(st.doc,"region_marks",[])
    selected=disk_id===nothing ? marks : [m for m in marks if m["id"]==disk_id]
    length(selected)==1 || throw(ArgumentError("select exactly one stored joining disk"))
    mark=only(selected)
    M.check_joining_disk(st,mark["seed"])["koda_joining"] ||
        throw(ArgumentError("selected disk fails the stored-orientation Koda conditions"))
    labels=sort([x["id"] for x in core["intersections"]])
    pids=Dict(x=>i for (i,x) in enumerate(labels))
    curves=Dict(c["id"] => (Symbol(f),i) for f in ("alpha","beta") for (i,c) in enumerate(core[f]))
    specs(f)=[Curve(Symbol(f),i,Int[pids[x] for x in c["cycle"]];label=c["id"]) for (i,c) in enumerate(core[f])]
    xs=[Intersection(pids[x["id"]],curves[x["alpha"]][2],curves[x["beta"]][2],x["sign"];source_label=x["id"]) for x in core["intersections"]]
    h=certify_heegaard(decorate(heegaard_code(specs("alpha"),specs("beta"),xs);cellular=true))
    surface_genus(h)==core["genus"] || error("cellulated projection changed genus")
    graph=ribbon_graph(h); rb=boundary_cycles(h)
    at=Dict(M.root!(st.vertex_dsu,Tuple(x["at"]))=>pids[x["id"]] for x in st.doc["crossings"])
    region_map=Int[]; witnesses=Any[]
    for r in rs
        events=[]
        for (f,_,u,v) in only(r["boundary_cycles"])
            p=get(at,M.root!(st.vertex_dsu,(f,u)),nothing)
            p===nothing && continue
            cid=only(st.used[(f,M.edge(u,v))])
            saved=only(d for d in st.curves[cid][2] if d.face==f && M.edge(d.u,d.v)==M.edge(u,v))
            push!(events,(p,curves[cid],(saved.u,saved.v)==(u,v)))
        end
        isempty(events) && error("region without intersection events is unsupported")
        darts=Int[]
        for i in eachindex(events)
            (p,key,forward)=events[i];q=events[mod1(i+1,length(events))][1]
            candidates=[d.id for d in graph.darts if d.origin==p && d.target==q &&
                (graph.edges[d.edge].family,graph.edges[d.edge].curve_id)==key && isodd(d.id)==forward]
            push!(darts,only(candidates))
        end
        j=only(findall(c->_cyclic_equal(c,darts),rb))
        push!(region_map,j)
        push!(witnesses,Dict("seed"=>r["seed"],"ribbon_region"=>j,"boundary_darts"=>darts))
    end
    sort(region_map)==collect(eachindex(rb)) || error("marked/ribbon regions are not bijective")
    j=only(findall(r->mark["seed"] in r["triangles"],rs))
    ph=make_punctured(h;joining_disk=joining_disk(h,region_map[j]),source_identity=:labelled)
    provenance=Dict("source_representation_sha256"=>M.representation_digest(states.document),
        "source_marking_sha256"=>M.marking_digest(states.document),"disk_id"=>mark["id"],
        "mathematical_projection_sha256"=>M.representation_digest(st.doc),
        "auxiliary_curves"=>sort([c["id"] for c in states.document["curves"] if c["role"]=="aux"]),
        "intersection_labels"=>Dict(string(i)=>x for (i,x) in enumerate(labels)),
        "region_map"=>witnesses,"policy"=>"marked input retained; cellular projection for conversion only; no rerouting")
    (object=ph,provenance=provenance)
end

"""Run the existing cellular Koda construction on a full marked document.

The entire marked source is retained for inspection, not used as an inverse.
The region/dart correspondence is checked before conversion; unsupported
non-cellular or non-minimal systems are rejected explicitly by the adapter.
"""
function to_odata(doc::AbstractDict;disk_id=nothing)
    if get(doc,"schema","")=="marked-cut-view/1"
        bundle=to_odata(marked_display_source(doc);disk_id=disk_id)
        bundle.ograph.metadata["display_source_incidence_checked"]=true
        return bundle
    end
    projected=cellulated_punctured(doc;disk_id=disk_id)
    bundle=to_odata(projected.object)
    bundle.ograph.metadata["marked_projection"]=projected.provenance
    bundle.ograph.metadata["marked_source"]=_plain_value(doc)
    bundle.certificate.label_map["marked_projection"]=deepcopy(projected.provenance)
    bundle.comparison["marked_source_retained"]=true
    bundle.comparison["inverse_constructed"]=false
    bundle
end
