# cf: docs/specs/heegaard-snappy.md — fresh one-intersection cancellation and marked-region transport.
# cf: answer-10-2.md §§1–5 — fresh geometric cancellation via row/column bands and punctured-torus excision; reconstructed.
function _hc_destabilize_ribbon(diagram;alpha::Int,beta::Int,joining_region::Int)
    bad(msg)=throw(LCMoveError(:destabilization_support,msg,Dict{String,Any}()))
    code=_undecorated(diagram);g=surface_genus(diagram);g>1||bad("positive output genus required")
    input_report=compression_boundary_report(diagram)
    input_report.alpha_complete&&input_report.beta_complete||bad("input systems are not complete")
    joining_disk(diagram,joining_region)
    a0=code.alpha[alpha].visits;b0=code.beta[beta].visits
    pivots=intersect(a0,b0);length(pivots)==1||bad("selected pair must meet exactly once")
    pivot=only(pivots);eps=Dict(p.id=>p.sign for p in code.intersections)
    after(w,p)=begin i=findfirst(==(p),w);[w[mod1(i+k,length(w))] for k in 1:length(w)-1] end
    aa=after(a0,pivot);bb=after(b0,pivot);removed=union(Set(a0),Set(b0))
    fresh=maximum(keys(eps));grid=Dict{Tuple{Int,Int},Int}();newsign=Dict{Int,Int}()
    for q in bb,r in aa
        fresh+=1;grid[(q,r)]=fresh;newsign[fresh]=-eps[q]*eps[r]*eps[pivot]
    end
    av=Vector{Int}[];bv=Vector{Int}[];aorder=setdiff(eachindex(code.alpha),[alpha]);border=setdiff(eachindex(code.beta),[beta])
    for i in aorder
        word=Int[]
        for q in code.alpha[i].visits
            if q in bb
                order=-eps[q]*eps[pivot]==1 ? aa : reverse(aa)
                append!(word,[grid[(q,r)] for r in order])
            else
                q in removed&&bad("unexpected removed alpha crossing")
                push!(word,q)
            end
        end
        isempty(word)&&bad("empty remaining alpha needs a noncellular chart")
        push!(av,word)
    end
    for i in border
        word=Int[]
        for r in code.beta[i].visits
            if r in aa
                order=-eps[r]*eps[pivot]==1 ? bb : reverse(bb)
                append!(word,[grid[(q,r)] for q in order])
            else
                r in removed&&bad("unexpected removed beta crossing")
                push!(word,r)
            end
        end
        isempty(word)&&bad("empty remaining beta needs a noncellular chart")
        push!(bv,word)
    end
    labels=sort(unique(reduce(vcat,av;init=Int[])))
    points=[Intersection(v,only(i for i in eachindex(av) if v in av[i]),only(i for i in eachindex(bv) if v in bv[i]),get(newsign,v,get(eps,v,0))) for v in labels]
    core=heegaard_code([Curve(:alpha,i,w) for (i,w) in enumerate(av)],[Curve(:beta,i,w) for (i,w) in enumerate(bv)],points)
    out=certify_heegaard(decorate(core;cellular=true))
    actual=surface_genus(out)
    actual==g-1||bad("surface genus $actual, expected $(g-1)")
    rep=compression_boundary_report(out);rep.alpha_complete&&rep.beta_complete||bad("remaining systems not complete")
    oldgraph=ribbon_graph(diagram);newgraph=ribbon_graph(out)
    newfaces=boundary_cycles(newgraph);faceof=Dict(d=>i for (i,f) in enumerate(newfaces) for d in f)
    outgoing=Dict((d.origin,d.origin_port)=>d.id for d in newgraph.darts)
    anchors=Any[];images=Int[]
    for did in boundary_cycles(oldgraph)[joining_region]
        d=oldgraph.darts[did];d.origin in removed&&continue
        nd=outgoing[(d.origin,d.origin_port)];push!(images,faceof[nd]);push!(anchors,Dict("vertex"=>d.origin,"port"=>string(d.origin_port),"old_dart"=>did,"new_dart"=>nd,"new_face"=>faceof[nd]))
    end
    length(unique(images))==1||bad("marked disk has $(length(unique(images))) images via retained corners")
    region=only(unique(images));joining_disk(out,region)
    cert=Dict("pivot"=>pivot,"alpha_index"=>alpha,"beta_index"=>beta,"pivot_sign"=>eps[pivot],"removed_labels"=>sort(collect(removed)),"grid"=>[Dict("q"=>q,"r"=>r,"label"=>v,"sign"=>newsign[v]) for ((q,r),v) in sort(collect(grid))],"alpha_source_indices"=>aorder,"beta_source_indices"=>border,"alpha_visits"=>av,"beta_visits"=>bv,"genus_before"=>g,"genus_after"=>actual,"joining_region_before"=>joining_region,"joining_region_after"=>region,"joining_anchors"=>anchors)
    (diagram=out,joining_region=region,certificate=cert)
end

"""Cancel a selected once-intersecting alpha/beta pair on a cellular HC.
`alpha` and `beta` index the forgotten diagram orders; `alpha_basis` is required.
The exact named joining region is transported through the surface surgery.
Returns the existing standard-constructor result with cancellation correspondence.
"""
# cf: koda_forward.jl::_hc_disk_selection/_resolve_hc_named_disk; preserve the exact named LC support.
"""
    hc_destabilize(hc; alpha::Int, beta::Int, alpha_basis::Symbol, joining_disk=nothing, disk_name=nothing, disk_id=nothing) -> NamedTuple

Return a standard HC after a supported one-intersection handle cancellation.

# Arguments
`hc` is a validated HC in the domain of
[HC surface operations](../specs/heegaard-snappy.md#hc-の曲面操作).
`alpha` and `beta` select the forgotten diagram's curve orders; in-range
indices are the caller's responsibility (`BoundsError` otherwise).
`alpha_basis` is required; values other than `:l` or `:m` raise
`HCTriangulationError` naming the standard-basis choice.
Supply an explicit joining-disk selection through the named keywords.
Unsupported marks, noncellular input, an unavailable disk transport, failed
one-intersection support or incomplete cut systems raise `LCMoveError`
with the scope or support reason. Input/selection errors from HC and Koda
validation propagate. Output certification errors also propagate.

# Returns
The named-tuple result of `hc_koda_to_standard`, including `lc`, `hc` and
`certificate`. The certificate adds `destabilization` (curve, crossing,
genus and disk-corner correspondence), `input_joining_region` and
`input_joining_faces`.

# Guarantees
1. A successful result has a certified cellular legacy diagram of genus
   one less than the input, with complete alpha and beta systems and the
   transported joining disk. [proved/agent: by construction]
2. The resulting HC is intended to describe the same closed 3-manifold
   through the specified cancellation. [unverified: a general geometric identification of the emitted standard HC with that surgery is not certified here]
3. The input HC is not modified. [proved/agent: by construction]

# Not guaranteed
The result does not prove minimal genus or independently identify the output
3-manifold. The general manifold-preservation statement in guarantee 2 is
not a consequence of homology agreement alone.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> root = normpath(joinpath(dirname(pathof(Topology)), "..", "..", ".."));

julia> h = read_hc(joinpath(root, "threads", "PJ260922", "data", "261112-example-table", "E1-summand-L3", "turn9-public-standard-l", "diagram.hc"));

julia> try hc_destabilize(h; alpha=1, beta=1, alpha_basis=:l, disk_name="joining"); false catch e; e isa LCMoveError end
true
```

See also [`hc_koda_to_standard`](@ref), [`hc_first_homology`](@ref).

Role: move
"""
function hc_destabilize(hc;alpha::Int,beta::Int,alpha_basis::Symbol,joining_disk=nothing,disk_name=nothing,disk_id=nothing)
    validate_hc_input(hc)
    isempty(hc.lc.disks)&&isempty(hc.lc.points)&&isempty(hc.lc.weights)&&length(hc.disks)<=1 ||
        throw(LCMoveError(:destabilization_scope,"additional mark transport is unsupported",Dict{String,Any}()))
    isempty(hc.points)&&isempty(hc.weights)||throw(LCMoveError(:destabilization_scope,"destabilization point/weight transport unsupported",Dict{String,Any}()))
    forgot=forget_reference(realize_hc(hc));forgot.diagram===nothing&&throw(LCMoveError(:destabilization_scope,"missing source ribbon diagram",Dict{String,Any}()))
    all(r->(r.chi,r.genus,r.boundary_components)==(1,0,1),forgot.regions)||throw(LCMoveError(:destabilization_scope,"destabilization requires cellular source complement",Dict{String,Any}()))
    selection=_hc_disk_selection(hc;joining_disk=joining_disk,disk_name=disk_name,disk_id=disk_id)
    selection===nothing&&throw(LCMoveError(:destabilization_scope,"explicit joining disk required",Dict{String,Any}()))
    selection=_resolve_hc_named_disk(hc,forgot,selection);region=_resolve_joining_region(forgot,selection)
    region.joining_disk||throw(LCMoveError(:destabilization_scope,"selected mark is not a joining disk",Dict{String,Any}()))
    disk=_find_ribbon_joining_disk(forgot.diagram,region);disk===nothing&&throw(LCMoveError(:destabilization_scope,"marked region has no ribbon boundary transport",Dict{String,Any}()))
    reduced=_hc_destabilize_ribbon(forgot.diagram;alpha=alpha,beta=beta,joining_region=disk.region)
    output=hc_koda_to_standard(reduced.diagram;joining_region=reduced.joining_region,alpha_basis=alpha_basis)
    output.certificate["destabilization"]=reduced.certificate
    output.certificate["input_joining_region"]=region.id;output.certificate["input_joining_faces"]=region.face_roots
    output
end
