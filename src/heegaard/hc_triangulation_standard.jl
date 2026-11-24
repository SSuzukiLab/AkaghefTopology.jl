# cf: docs/specs/heegaard-snappy.md — alpha-adapted standard inverse and source-cell certificate.
# cf: .research-scratch/PJ260922/astra-lanes/snappy-bridge/answer-9-1.md §§実装・必要なAPI前提 — signed dual routes and source-face ancestry; reconstructed using existing band-boundary primitives.
# cf: lc_moves_ribbon.jl::_move_certify_route, _move_barbell_boundary, _move_compile_boundary — ordered banks and exact word edits.

function _hc_standard_duals(diagram)
    bad(msg)=throw(HCTriangulationError(:standard_reference_construction,msg))
    rep=compression_boundary_report(diagram)
    rep.alpha_complete || bad("alpha not complete")
    code=_undecorated(diagram);g=rep.surface_genus;limit=3g-2
    oldlabels=sort([p.id for p in code.intersections])
    labels=Dict(v=>limit+i for (i,v) in enumerate(oldlabels))
    eps=Dict(p.id=>p.sign for p in code.intersections)
    words=Dict{String,Vector{Int}}()
    # forget_reference stores ε = +1 when the plus passage is on α.
    # Writing α = −ε would reverse every inherited passage and flip the linking form.
    for (i,c) in enumerate(code.alpha)
        words["l$i"]=[eps[v]*labels[v] for v in c.visits]
    end
    for (i,c) in enumerate(code.beta)
        words["beta$i"]=[-eps[v]*labels[v] for v in c.visits]
    end
    base=deepcopy(words);fresh=maximum(values(labels));stages=Any[]
    for i in 1:g
        faces,right=_lc_faces(words);occ=_lc_rotation_occ(words)
        tau(c)=_lc_tau(words,occ,c)
        cw(x)=(abs(x),x>0 ? "+" : "-","out")
        l="l$i";seed=cw(first(words[l]));start=right[tau(seed)];goal=right[seed]
        adj=Dict(f.id=>Any[] for f in faces)
        for j in 1:length(code.beta),x in words["beta$j"]
            c=cw(x);d=tau(c);a=right[c];b=right[d]
            push!(adj[a],(next=b,dart=c));push!(adj[b],(next=a,dart=d))
        end
        parent=Dict{String,Any}(start=>nothing);queue=[start]
        for face in queue
            haskey(parent,goal)&&break
            for edge in adj[face]
                haskey(parent,edge.next)&&continue
                parent[edge.next]=(face=face,dart=edge.dart);push!(queue,edge.next)
            end
        end
        haskey(parent,goal)||bad("no dual path at handle $i")
        path=Any[];face=goal
        while face!=start
            step=parent[face];push!(path,step.dart);face=step.face
        end
        reverse!(path)
        insertions=Dict{Tuple{String,Int},Vector{Int}}();mw=Int[i]
        for c in path
            name,pos=occ[c[1]][c[2]]
            edge=c[3]=="out" ? pos+1 : mod1(pos,length(words[name]))
            # _lc_rotation_occ stores zero-based occurrence positions.
            fresh+=1;token=(c[3]=="out" ? 1 : -1)*fresh
            push!(mw,token);push!(get!(insertions,(name,edge),Int[]),-token)
        end
        for (name,w) in collect(words)
            next=Int[]
            for k in eachindex(w)
                push!(next,w[k]);append!(next,get(insertions,(name,k),Int[]))
                name==l&&k==1&&push!(next,-i)
            end
            words[name]=next
        end
        words["m$i"]=mw
        nextfaces,_=_lc_faces(words)
        labs=Set(abs(x) for w in values(words) for x in w)
        chi=length(nextfaces)-length(labs)
        chi==2-2g||bad("dual overlay changed genus at handle $i: chi=$chi")
        push!(stages,Dict("handle"=>i,"path"=>path,"chi"=>chi,"words"=>deepcopy(words)))
    end
    (words=words,source_words=base,labels=labels,stages=stages)
end

function _hc_standard_reference(diagram;alpha_basis::Symbol)
    bad(msg)=throw(HCTriangulationError(:standard_reference_construction,msg))
    dual=_hc_standard_duals(diagram);words=deepcopy(dual.words)
    alpha_basis in (:l,:m)||bad("alpha_basis must be :l or :m")
    if alpha_basis==:m
        reversed=Set(abs(x) for (n,w) in words if startswith(n,"m") for x in w)
        rename(n)=startswith(n,"l") ? "m"*n[2:end] : startswith(n,"m") ? "l"*n[2:end] : n
        words=Dict(rename(n)=>(startswith(n,"m") ? reverse([abs(x) in reversed ? -x : x for x in w]) :
                  [abs(x) in reversed ? -x : x for x in w]) for (n,w) in words)
    end
    g=compression_boundary_report(diagram).surface_genus
    g>=1||bad("positive genus required")
    for j in 1:g-1;words["n$j"]=Int[];end
    routes=Any[]
    make_lc(w)=LickorishCode(g,w,Dict{String,String}(),Dict{String,Vector{String}}(),
                                     Dict{String,String}(),Dict{String,Int}(),nothing)
    for i in 1:g-1
        lc=make_lc(words);active=_lc_nonempty_words(words)
        faces,right=_lc_faces(active);occ=_lc_rotation_occ(active)
        tau(c)=_lc_tau(active,occ,c)
        start=(i,"+","out");finish=(i+1,"+","in")
        adj=Dict(f.id=>Any[] for f in faces)
        for j in 1:length(_undecorated(diagram).beta),x in words["beta$j"]
            c=(abs(x),x>0 ? "+" : "-","out");d=tau(c)
            push!(adj[right[c]],(next=right[d],dart=c));push!(adj[right[d]],(next=right[c],dart=d))
        end
        initial=right[start];goal=right[finish]
        parent=Dict{String,Any}(initial=>nothing);queue=[initial]
        for face in queue
            haskey(parent,goal)&&break
            for e in adj[face]
                haskey(parent,e.next)&&continue
                parent[e.next]=(face=face,dart=e.dart);push!(queue,e.next)
            end
        end
        haskey(parent,goal)||bad("connector band unavailable at handle $i")
        path=Any[];face=goal
        while face!=initial
            step=parent[face];push!(path,step.dart);face=step.face
        end
        reverse!(path)
        point(c)=Dict("dart"=>collect(c),"u"=>"1/2")
        band=Dict("start"=>point(start),"finish"=>point(finish),"crossings"=>point.(path))
        avoid=[n for n in keys(words) if !startswith(n,"beta")]
        route=_move_certify_route(lc,band;avoid=avoid)
        si,sj=-_move_dir(start),-_move_dir(finish)
        cycle=_move_barbell_boundary(si,sj);runs=Any[]
        for d in cycle
            if abs(d)==1
                push!(runs,(kind=:curve,component="m$i",edge=route.start.edge,side=si,direction=sign(d),survives=false))
            elseif abs(d)==2
                push!(runs,(kind=:curve,component="m$(i+1)",edge=route.finish.edge,side=sj,direction=sign(d),survives=false))
            else
                push!(runs,(kind=:bank,route=route,side=si*sign(d),direction=sign(d)))
            end
        end
        row=Dict("id"=>"reference-construction","inverse_id"=>nothing,"digest"=>"",
                 "relation"=>"insertion","layer"=>"LC")
        witness=_move_compile_boundary(lc,row,"n$i",runs,Dict("construction_only"=>true))
        projected,_,_,_,_,_,_=_kernel_apply(witness.names,lc.words,witness.before_frags,
            witness.after_frags,witness.locations,witness.vertices,witness.capacity;
            label_map=witness.extra["label_map"])
        words=Dict(witness.names[j]=>projected[j] for j in eachindex(projected))
        links=Dict{Int,Int}()
        for j in (i,i+1)
            labs=intersect(Set(abs.(words["l$j"])),Set(abs.(words["n$i"])))
            length(labs)==1||bad("connector $i has wrong intersection count with l$j")
            links[only(labs)]=g+2i-1+(j-i)
        end
        words=Dict(n=>[sign(x)*get(links,abs(x),abs(x)) for x in w] for (n,w) in words)
        push!(routes,Dict("handle"=>i,"cycle"=>cycle,"route"=>route))
    end
    full=make_lc(words);validate_lc_input(full)
    limit=3g-2
    projection=make_lc(Dict(n=>filter(x->abs(x)<=limit,words[n]) for n in _lc_k_names(g)))
    _lc_validate_k_constants!(projection)

    (lc=full,words=words,labels=dual.labels,dual=dual,routes=routes,alpha_basis=alpha_basis)
end

function _hc_standard_finalize(diagram,lc,labels;joining_region::Integer,alpha_basis::Symbol,extra=Dict{String,Any}())
    bad(msg)=throw(HCTriangulationError(:standard_source_transport,msg))
    code=_undecorated(diagram);graph=ribbon_graph(diagram)
    eps_source=Dict(v.id=>v.sign for v in code.intersections)
    inherited=Set(values(labels));curve_checks=Any[]
    for (family,curves) in ((:alpha,code.alpha),(:beta,code.beta)),(i,c) in enumerate(curves)
        name=family==:alpha ? "$(alpha_basis)$i" : "beta$i"
        expected=[(family==:alpha ? 1 : -1)*eps_source[v]*labels[v] for v in c.visits]
        observed=filter(x->abs(x) in inherited,lc.words[name])
        _lc_canon_cyclic(observed)==_lc_canon_cyclic(expected)||bad("original curve changed: $name")
        push!(curve_checks,Dict("family"=>String(family),"source_id"=>c.id,"lc_curve"=>name,"preserved"=>true))
    end
    alpha_basis in (:l,:m)||bad("unsupported basis")
    roles=Dict{String,Symbol}("$(alpha_basis)$i"=>:alpha for i in 1:lc.genus)
    merge!(roles,Dict("beta$i"=>:beta for i in eachindex(code.beta)))
    hc=resolve_hc(lc,roles);realized=realize_hc(hc)
    forgot=forget_reference(realized);realization=realized.value["realization"]
    lc_face_region=Dict(f=>r.id for r in forgot.regions for f in r.face_roots)
    sourcefaces=boundary_cycles(graph);face_of=Dict(d=>i for (i,f) in enumerate(sourcefaces) for d in f)
    outgoing=Dict((d.origin,d.origin_port)=>d.id for d in graph.darts)
    eps=Dict(v.id=>v.sign for v in code.intersections)
    all((r.chi,r.genus,r.boundary_components)==(1,0,1) for r in forgot.regions)||
        bad("a forgotten source face is not a disk")
    source_of=Dict{Int,Int}();checks=Any[]
    for v in graph.vertices,(sg,way) in _LC_ROT_SEQ
        fam=sg==_lc_sign_char(eps[v]) ? :alpha : :beta
        dart=outgoing[(v,_port(fam,Symbol(way)))]
        # α carries the plus passage when ε = +1, so this corner bounds the
        # outgoing dart's face. The opposite dart is the α = −ε correspondence.
        sf=face_of[dart];label=labels[v]
        region=lc_face_region[realization.corner_face[(label,sg,way)]]
        haskey(source_of,region)&&source_of[region]!=sf&&bad("source faces merged")
        source_of[region]=sf;push!(checks,Dict("source_vertex"=>v,"corner"=>(sg,way),"source_face"=>sf,"forgotten_region"=>region))
    end
    length(source_of)==length(forgot.regions)==length(sourcefaces)||bad("region count changed")
    Set(values(source_of))==Set(eachindex(sourcefaces))||bad("face transport is not bijective")
    matches=[r for r in forgot.regions if source_of[r.id]==joining_region]
    length(matches)==1||bad("physical source disk not mapped uniquely")
    region=only(matches);region.joining_disk||bad("mapped physical region is not joining")
    hc=resolve_hc(lc,roles;disks=Dict("joining"=>region.face_roots))
    (lc=lc,hc=hc,joining_region_id=region.id,joining_faces=region.face_roots,
     certificate=merge(Dict("alpha_basis"=>String(alpha_basis),"source_joining_region"=>joining_region,
        "joining_region_id"=>region.id,"corner_checks"=>checks,"source_of_forgotten_region"=>source_of,
        "source_labels"=>labels,"source_curve_checks"=>curve_checks,"all_source_regions_disks"=>true),extra))
end


"""Construct reference-bound alpha and free beta words on the original Koda surface.

`alpha_basis` is required and accepts `:l` or `:m`. The source must be cellular,
closed, positive genus, with complete minimal alpha and beta systems. The given
joining region is transported as the whole union of its child LC faces.
"""
function hc_koda_to_standard(diagram;joining_region::Integer,alpha_basis::Symbol)
    alpha_basis in (:l,:m)||throw(HCTriangulationError(:standard_alpha_basis,"choose :l or :m explicitly"))
    is_cellular(diagram)||throw(HCTriangulationError(:standard_source_not_cellular,"disk regions required"))
    report=compression_boundary_report(diagram)
    report.surface_genus>0||throw(HCTriangulationError(:standard_positive_genus,"positive genus required"))
    report.alpha_complete||throw(HCTriangulationError(:standard_alpha_not_complete,"alpha cut is not one sphere with 2g holes"))
    report.beta_complete||throw(HCTriangulationError(:standard_beta_not_complete,"closed minimal beta system required"))
    joining_disk(diagram,joining_region)
    x=_hc_standard_reference(diagram;alpha_basis)
    result=_hc_standard_finalize(diagram,x.lc,x.labels;joining_region,alpha_basis,
        extra=Dict("alpha_cut"=>report.alpha,"dual_stages"=>x.dual.stages,
                   "connector_routes"=>x.routes,"reference_words_checked"=>true))
    result
end

"""Actual triangulation → Koda → reference-bound HC, with an explicit alpha basis."""
function hc_triangulation_to_standard(input;alpha_basis::Symbol,orders=nothing,max_nodes::Int=100_000)
    flow=hc_flow_certificate(input;orders,max_nodes)
    result=hc_koda_to_standard(flow.diagram;joining_region=flow.region,alpha_basis)
    result.certificate["flow_certificate"]=flow.certificate
    result
end

export hc_koda_to_standard, hc_triangulation_to_standard
