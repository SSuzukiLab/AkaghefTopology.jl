# Source: docs/specs/heegaard-koda-conversion.md §7.4, §7.3;
# conversion_spine.jl::to_ograph (inverse wiring only).
"""Decode o-data to a combinatorial branched spine; encode spine to o-data."""

struct DecodedSpine
    genus_hint::Int
    vertex_signs::Dict{Int, Int}
    edges::Vector{SpineEdge}
    gauss::Vector{Vector{Int}}
    metadata::Dict{String, Any}
end

function _odata_vertex_labels(G::OData)
    n = length(G.orientation)
    n == 0 && throw(KodaInputError(:odata_to_spine, "empty orientation vector", Dict()))
    1:n
end

function odata_to_spine(G::OData; tables=load_koda_local_tables())
    T = parentmodule(@__MODULE__)
    link = T.VirtualLink()
    T.set_data!(link; rgauss=deepcopy(G.gauss), orientation=copy(G.orientation))
    closed, cond = T.is_closed(link)
    closed || throw(KodaInputError(:odata_to_spine, "o-data must satisfy C1–C3", Dict("conditions" => cond)))
    labels = _odata_vertex_labels(G)
    sigma = Dict(i => G.orientation[i] for i in labels)
    table_spine = try
        koda_table_decode(G.gauss, sigma; tables=tables)
    catch err
        err isa KodaDomainExcluded && rethrow()
        err isa KodaConversionFailure && rethrow()
        throw(KodaImplementationError(:odata_to_spine, "table decode failed", Dict("error" => sprint(showerror, err), "input" => G.gauss)))
    end
    label_map = Dict(p => i for (i, p) in enumerate(labels))
    next = Dict{Tuple{Int, Symbol}, Tuple{Int, Symbol}}()
    for word in G.gauss
        prev_state = nothing
        start_state = nothing
        for (pos, p) in enumerate(word)
            v = labels[abs(p)]
            fam = p > 0 ? :beta : :alpha
            state = (v, fam)
            prev_state !== nothing && (next[prev_state] = state)
            start_state === nothing && (start_state = state)
            prev_state = state
            pos == length(word) && (next[state] = start_state)
        end
    end
    edges = SpineEdge[]
    for (src_v, src_f) in keys(next)
        tgt_v, tgt_f = next[(src_v, src_f)]
        push!(edges, SpineEdge(src_v, src_f, tgt_v, tgt_f, Int[]))
    end
    meta = Dict{String, Any}(
        "decoded_from" => "odata",
        "convention" => get(tables.orientation, "convention", "candidate_A"),
        "vertex_to_label" => Dict(string(p) => label_map[p] for p in labels),
        "table_spine" => table_spine,
        "status" => "constructed",
    )
    DecodedSpine(get(table_spine, "genus", 0), sigma, edges, deepcopy(G.gauss), meta)
end

function spine_to_odata(spine::DecodedSpine)
    OData(deepcopy(spine.gauss), collect(spine.vertex_signs[i] for i in sort(collect(keys(spine.vertex_signs)))), Dict{String, Any}(
        "representation" => "normal-o-graph",
        "convention" => "candidate_A",
        "source_status" => "constructed",
    ))
end

function _non_disk_sector(spine::FlowSpine)
    sectors = get(spine.complex.checks, "sectors", Dict{String, Any}[])
    for sec in sectors
        get(sec, "genus", 0) == 0 && get(sec, "boundary_count", 1) == 1 && continue
        return sec
    end
    nothing
end

# cf: lib/topology.jl/src/virtual_link.jl set_data! — store the o-graph Gauss code and vertex signs on a VirtualLink
function odata_to_virtual_link(G::OData)
    T = parentmodule(@__MODULE__)
    link = T.VirtualLink()
    T.set_data!(link; rgauss=deepcopy(G.gauss), orientation=copy(G.orientation))
    link
end

# cf: lib/topology.jl/src/virtual_link.jl VirtualLink — gauss and orientation are the stored o-data presentation
function virtual_link_to_odata(link)
    T = parentmodule(@__MODULE__)
    link isa T.VirtualLink || throw(KodaInputError(
        :virtual_link_to_odata,
        "expected a VirtualLink",
        Dict(),
    ))
    OData(deepcopy(link.gauss), copy(link.orientation), Dict{String, Any}(
        "representation" => "normal-o-graph",
        "source" => "VirtualLink",
    ))
end

# cf: [[K08#^k08-sec-0073]] §3 — beta isotopy; source-cell annulus cut is reconstructed.
# cf: lib/topology.jl/src/heegaard/lc_moves_support.jl::_move_left_boundary — convert ribbon/right banks to the existing left-face digon table.
function _koda_annulus_regularize(spine::FlowSpine)
    fail(reason, detail=Dict{String,Any}()) = koda_not_computed(:spine_to_odata,reason,detail)
    spine.source isa PuncturedHeegaardDiagram || fail(:regularization_source_required)
    regs=get(spine.complex.checks,"source_surface_regions",nothing)
    regs===nothing && fail(:regularization_surface_cells_required)
    annuli=[r for r in regs if r["genus"]==0 && r["boundary_count"]==2 && r["chi"]==0]
    length(annuli)==1 || fail(:regularization_requires_one_annulus)
    ann=only(annuli)
    all(r===ann || (r["genus"]==0 && r["boundary_count"]==1 && r["chi"]==1) for r in regs) ||
        fail(:regularization_other_surface_sector)
    ph=spine.source; h=ph.diagram; graph=ribbon_graph(h); old=_undecorated(h)
    length(old.intersections)-length(graph.edges)+sum(r["chi"] for r in regs)==2-2*spine.genus ||
        fail(:regularization_surface_euler_mismatch)
    banks=ann["boundary_darts"];length(banks)==2 || fail(:regularization_annulus_banks)
    removed=Set(graph.darts[d].edge for d in ph.joining_disk.boundary_darts)
    pairs=Tuple{Int,Int}[]
    for side in 1:2, p in banks[side], c in banks[3-side]
        pe=graph.edges[graph.darts[p].edge]; ce=graph.edges[graph.darts[c].edge]
        pe.family==:beta && ce.family==:alpha || continue
        pe.id in removed && continue
        ce.id in removed && continue
        push!(pairs,(p,c))
    end
    isempty(pairs) && fail(:regularization_no_relative_finger)
    # One topology-selected route; never choose another pair after a failed check.
    pd,cd=first(sort!(pairs)); pe=graph.edges[graph.darts[pd].edge];ce=graph.edges[graph.darts[cd].edge]
    cp,kp=_derived_digon(!isodd(pd),!isodd(cd))
    n=maximum(x.id for x in old.intersections)
    function refine(curve)
        edge=curve.family==:beta ? pe : ce; pat=curve.family==:beta ? cp : kp
        word=copy(curve.visits)
        if curve.family==edge.family && curve.id==edge.curve_id
            pos=only(i for i in eachindex(word) if word[i]==edge.from && word[mod1(i+1,length(word))]==edge.to)
            word=vcat(word[1:pos],n .+ abs.(pat),word[pos+1:end])
        end
        Curve(curve.family,curve.id,word;label=curve.label)
    end
    signs=Dict(abs(x)=>sign(x) for x in kp)
    xs=vcat(old.intersections,[Intersection(n+i,ce.curve_id,pe.curve_id,signs[i]) for i in 1:2])
    nh=certify_heegaard(decorate(heegaard_code(refine.(old.alpha),refine.(old.beta),xs);cellular=true))
    surface_genus(nh)==spine.genus || fail(:regularization_changed_surface_genus,
        Dict("before"=>spine.genus,"after"=>surface_genus(nh)))
    ng=ribbon_graph(nh);cycles=boundary_cycles(ng)
    length(cycles)==length(regs)+1 || fail(:regularization_cell_count_mismatch)
    key(g,d)=begin x=g.darts[d];e=g.edges[x.edge];(e.family,e.curve_id,x.origin,x.target,isodd(d)) end
    jkeys=Set(key(graph,d) for d in ph.joining_disk.boundary_darts)
    hits=[i for (i,c) in enumerate(cycles) if length(c)==length(jkeys) && Set(key(ng,d) for d in c)==jkeys]
    length(hits)==1 || fail(:regularization_changed_joining_boundary)
    jd=joining_disk(nh,only(hits));qph=make_punctured_combinatorial(nh,jd)
    # Deleting precisely the generated crossings recovers every oriented curve.
    for (before,after) in zip(vcat(old.alpha,old.beta),vcat(_undecorated(nh).alpha,_undecorated(nh).beta))
        filter(x->x<=n,after.visits)==before.visits || fail(:regularization_inverse_word_mismatch)
    end
    result=_koda_combinatorial_spine_from_punctured(qph;genus=spine.genus)
    G=to_ograph(result.value);v=odata_to_virtual_link(G)
    closed,conditions=parentmodule(@__MODULE__).is_closed(v)
    closed || fail(:regularization_not_closed,Dict("conditions"=>conditions))
    G.metadata["regularization"]=Dict{String,Any}(
        "kind"=>"annulus_beta_finger","surface_region"=>ann["id"],
        "opposite_banks"=>banks,"push_dart"=>pd,"cross_dart"=>cd,
        "new_crossings"=>[n+1,n+2],"new_crossing_signs"=>[signs[1],signs[2]],
        "source_diagram_hash"=>_labelled_hash(h),
        "joining_boundary_before"=>[collect(key(graph,d)) for d in ph.joining_disk.boundary_darts],
        "joining_boundary_after"=>[collect(key(ng,d)) for d in jd.boundary_darts],
        "joining_boundary_unchanged"=>true,"inverse_delete_crossings"=>[n+1,n+2],
        "surface_genus_before"=>spine.genus,"surface_genus_after"=>surface_genus(nh),
        "surface_disk_regions_after"=>length(cycles),"closed_conditions"=>collect(conditions))
    G
end

function spine_to_odata(spine::FlowSpine; regularize::Bool=false)
    if regularize && (get(spine.complex.checks, "vertex_free", false) || _non_disk_sector(spine) !== nothing)
        return _koda_annulus_regularize(spine)
    end
    if get(spine.complex.checks, "vertex_free", false) || isempty(spine.vertex_signs)
        koda_not_computed(:spine_to_odata, :vertex_free_odata_not_encoded, Dict{String, Any}(
            "spec" => "heegaard-koda-conversion.md §7.0 / §13.3",
            "kind" => "not_computed",
        ))
    end
    bad = _non_disk_sector(spine)
    if bad !== nothing
        koda_not_computed(:spine_to_odata, :non_disk_sector, Dict{String, Any}(
            "spec" => "heegaard-koda-conversion.md §6.3",
            "kind" => "not_computed",
            "sector" => get(bad, "id", nothing),
            "genus" => get(bad, "genus", nothing),
            "boundary_count" => get(bad, "boundary_count", nothing),
        ))
    end
    to_ograph(spine)
end



function _odata_word_rotations(word::Vector{Int})
    n = length(word)
    [vcat(word[i:n], word[1:(i - 1)]) for i in 1:n]
end

function _odata_all_perms(n::Int)
    n <= 1 && return [Int[1]]
    n > 8 && throw(ArgumentError("correspondence search too large"))
    out = Vector{Int}[]
    used = fill(false, n)
    cur = Int[]
    function rec!()
        length(cur) == n && (push!(out, copy(cur)); return)
        for j in 1:n
            used[j] && continue
            used[j] = true
            push!(cur, j)
            rec!()
            pop!(cur)
            used[j] = false
        end
    end
    rec!()
    out
end

"""All label bijections preserving oriented component cycles and passages."""
function _odata_oriented_correspondences(
    source::Vector{Vector{Int}},
    source_orient::Vector{Int},
    target::Vector{Vector{Int}},
    target_orient::Vector{Int},
)
    source_labels = sort(collect(Set(abs(x) for w in source for x in w)))
    target_labels = sort(collect(Set(abs(x) for w in target for x in w)))
    (length(source_labels) != length(target_labels) || length(source) != length(target)) && return Dict{String, Any}[]
    label_to_src = Dict(l => i for (i, l) in enumerate(source_labels))
    found = Dict{String, Any}[]
    seen = Set{Tuple{Vararg{Int}}}()
    for component_perm in _odata_all_perms(length(target))
        any(length(source[i]) != length(target[component_perm[i]]) for i in 1:length(source)) && continue
        rots = [_odata_word_rotations(target[component_perm[i]]) for i in 1:length(source)]
        lim = [length(r) for r in rots]
        idx = ones(Int, length(source))
        while true
            aligned = [rots[i][idx[i]] for i in 1:length(source)]
            mapping = Dict{Int, Int}()
            inverse = Dict{Int, Int}()
            ok = true
            for (sw, tw) in zip(source, aligned)
                for (s, t) in zip(sw, tw)
                    (s > 0) != (t > 0) && (ok = false; break)
                    a, b = abs(s), abs(t)
                    (get(mapping, a, b) != b || get(inverse, b, a) != a) && (ok = false; break)
                    mapping[a], inverse[b] = b, a
                end
                ok || break
            end
            if ok
                key = Tuple(mapping[l] for l in source_labels)
                if !(key in seen)
                    push!(seen, key)
                    orient_ok = all(source_orient[label_to_src[l]] == target_orient[mapping[l]] for l in source_labels)
                    push!(found, Dict(
                        "map" => Dict(string(k) => v for (k, v) in mapping),
                        "orientation_matches" => orient_ok,
                    ))
                end
            end
            pos = length(source)
            while pos >= 1
                idx[pos] < lim[pos] && (idx[pos] += 1; break)
                idx[pos] = 1
                pos -= 1
            end
            pos == 0 && break
        end
    end
    found
end

function compare_odata_words(a::OData, b::OData)
    corrs = _odata_oriented_correspondences(a.gauss, a.orientation, b.gauss, b.orientation)
    orient_ok = any(c -> get(c, "orientation_matches", false), corrs)
    Dict(
        "gauss_equal" => a.gauss == b.gauss,
        "orientation_equal" => a.orientation == b.orientation,
        "equal_up_to_metadata" => !isempty(corrs) && orient_ok,
        "correspondence_count" => length(corrs),
        "compared_by" => "directed_word_vertex_sign_up_to_relabel",
    )
end
