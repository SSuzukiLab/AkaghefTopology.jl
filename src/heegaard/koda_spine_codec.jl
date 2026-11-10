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

function spine_to_odata(spine::FlowSpine)
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
