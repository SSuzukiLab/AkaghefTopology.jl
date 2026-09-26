"""Construct the sectors of K08 §5, B = S ∪ Dα ∪ Dβ, then delete D.

Each compression disk inherits the orientation whose boundary agrees with
the disk-left boundary of D. Removing a D-edge glues its two remaining sheets.
No layout, example name, or target o-data is consulted.
"""
function _koda_complex(H, D, edges, directions)
    graph = ribbon_graph(H)
    twin(d) = isodd(d) ? d+1 : d-1
    cycles = boundary_cycles(H)
    cells = Tuple{Vector{String},Vector{Int}}[]
    for (region, block) in enumerate(H.code.regions.boundary_cycle_blocks)
        region == D.region && continue
        push!(cells, (["surface:$region"], copy(cycles[only(block)])))
    end
    for c in _curve_vector(H.code.code)
        ids = [e.id for e in graph.edges if (e.family,e.curve_id)==(c.family,c.id)]
        darts = directions[(c.family,c.id)] == 1 ? 2 .* ids .- 1 : 2 .* reverse(ids)
        push!(cells, (["$(c.family):$(c.id)"], darts))
    end
    gluings = Dict{String,Any}[]
    for d in D.boundary_darts
        occurrences = [(i,j) for (i,(_,word)) in enumerate(cells)
                        for j in eachindex(word) if word[j] in (d,twin(d))]
        length(occurrences)==2 || throw(ArgumentError("D-edge has other than two remaining sheets"))
        (i,a),(j,b) = occurrences
        i != j || throw(ArgumentError("removal does not merge two disk sectors; standardization required"))
        ci,wi=cells[i]; cj,wj=cells[j]
        wi[a] == twin(wj[b]) || throw(ArgumentError("incoherent sector gluing orientation"))
        tail(w,k) = [w[mod1(k+t,length(w))] for t in 1:length(w)-1]
        merged = vcat(tail(wi,a),tail(wj,b))
        push!(gluings,Dict("joining_dart"=>d,"removed_edge"=>graph.darts[d].edge,
            "left_source_cells"=>copy(ci),"right_source_cells"=>copy(cj),
            "result_source_darts"=>copy(merged)))
        cells[i]=(vcat(ci,cj),merged)
        deleteat!(cells,j)
    end
    paths = Dict{Int,Vector{Int}}()
    start_at = Dict{Int,Int}()
    for (i,e) in enumerate(edges)
        paths[i]=e.source_darts
        paths[-i]=twin.(reverse(e.source_darts))
        start_at[first(paths[i])]=i
        start_at[first(paths[-i])]=-i
    end
    regions=SpineRegion[]
    for (origins,word) in cells
        k=findfirst(d->haskey(start_at,d),word)
        k === nothing && throw(ArgumentError("sector has no surviving singular vertex"))
        rotated=[word[mod1(k+t,length(word))] for t in 0:length(word)-1]
        boundary=Int[]; pos=1
        while pos<=length(rotated)
            e=get(start_at,rotated[pos],0)
            e!=0 || throw(ArgumentError("sector does not follow a smoothed edge"))
            path=paths[e]; stop=pos+length(path)-1
            stop<=length(rotated) && rotated[pos:stop]==path ||
                throw(ArgumentError("sector does not contain the complete singular-edge path"))
            push!(boundary,e); pos=stop+1
        end
        push!(regions,SpineRegion(sort(origins),rotated,boundary))
    end
    incidence=[Tuple{Int,Int}[] for _ in edges]
    for (r,region) in enumerate(regions), e in region.boundary
        push!(incidence[abs(e)],(r,sign(e)))
    end
    closed=all(regions) do r
        all(eachindex(r.boundary)) do i
            e,f=r.boundary[i],r.boundary[mod1(i+1,length(r.boundary))]
            target=e>0 ? edges[e].target_vertex : edges[-e].source_vertex
            origin=f>0 ? edges[f].source_vertex : edges[-f].target_vertex
            target==origin
        end
    end
    triple=all(length(x)==3 for x in incidence)
    branching=all(sort(last.(x))==[-1,1,1] for x in incidence)
    vertices=Set(v for e in edges for v in (e.source_vertex,e.target_vertex))
    chi=length(vertices)-length(edges)+length(regions)
    closed && triple && branching && chi==1 ||
        throw(ArgumentError("Koda sector incidence check failed: closed=$closed, triple=$triple, branching=$branching, chi=$chi"))
    BranchedSpineComplex(regions,gluings,incidence,Dict{String,Any}(
        "closed_attaching_words"=>closed,"three_sheet_incidence"=>triple,
        "oriented_branching"=>branching,"euler_characteristic"=>chi,
        "disk_sectors"=>true,"scope"=>"finite incidence and disk-gluing checks; not an independent manifold recognition"))
end

function _complex_json(c::BranchedSpineComplex)
    Dict("regions"=>[Dict("id"=>i,"source_cells"=>r.source_cells,
        "source_darts"=>r.source_darts,"boundary"=>r.boundary) for (i,r) in enumerate(c.regions)],
        "gluings"=>c.gluings,"edge_incidence"=>[[collect(x) for x in inc] for inc in c.edge_incidence],
        "checks"=>c.checks)
end
