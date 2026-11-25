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

"""
    odata_to_spine(G::OData; tables=load_koda_local_tables()) -> DecodedSpine

Decode a closed normal o-graph into a finite Koda spine and boundary record.

# Arguments
`G` uses [normal o-data, §7.1](../specs/heegaard-koda-conversion.md#71-保持するデータ).
The existing `VirtualLink` checks of C1–C3 are required; failure raises
`KodaInputError` with the condition results. An empty orientation vector
also raises `KodaInputError`. Errors while forming the `VirtualLink` can
propagate before those checks. `tables` must have the installed local-table
schema. Koda decoding errors propagate; other decoder failures are wrapped
in `KodaImplementationError` with the input and diagnostic.

# Returns
`DecodedSpine` with `genus_hint`, `vertex_signs`, directed `edges`, copied
`gauss` words, and `metadata` containing `table_spine`, label correspondence,
the convention and `status == "constructed"`. The genus hint is not a
certified Heegaard genus.

# Guarantees
1. The record preserves the input Gauss words and vertex signs and contains
   the finite local-table interpretation of their spine and boundary.
   [proved/agent: by construction]
2. Its intended geometric interpretation is the spine/thickening of §7,
   with the limits in §13.3. [unverified: the general source identification and P1 discrepancy are not resolved by this decoder]
3. The input and tables are not modified; the returned Gauss words are
   independent copies. [proved/agent: by construction]

# Not guaranteed
The input metadata is not preserved wholesale. Constructed status is not a
proof of the source identification, and the genus hint is not inferred from
a complete Heegaard construction.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§7–8, 13.3 — finite decoding and source-alignment scope.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); s = odata_to_spine(g);

julia> (s.gauss == [[1,-1]], s.vertex_signs[1])
(true, 1)
```

See also [`spine_to_heegaard_koda`](@ref), [`spine_to_odata`](@ref).

Role: convert
"""
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

"""
    spine_to_odata(spine::DecodedSpine) -> OData

Copy the stored signed o-data presentation of a decoded spine.

# Arguments
`spine` should come from `odata_to_spine`. For hand-built records, the caller
ensures that its Gauss labels and sign keys have the correspondence of
[Koda encoding §7.1](../specs/heegaard-koda-conversion.md#71-保持するデータ).
This method does not validate edges, closedness or that correspondence.

# Returns
An `OData` with copied Gauss words, signs in sorted vertex-key order, and
metadata for the normal-o-graph representation, candidate-A convention and
constructed source status.

# Guarantees
1. For a decoder-produced spine, the result retains its stored signed
   presentation. [proved/agent: by construction]
2. The input is not modified; word vectors and the sign vector are copied.
   [proved/agent: by construction]

# Not guaranteed
This does not reconstruct Gauss words from edited spine edges or preserve
arbitrary metadata. Constructed status is not a new geometric certificate.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §7.1 — data retained by the spine/o-data correspondence.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); b = spine_to_odata(odata_to_spine(g)); (b.gauss == g.gauss, b.orientation == g.orientation)
(true, true)
```

See also [`odata_to_virtual_link`](@ref).

Role: convert
"""
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
"""
    odata_to_virtual_link(G::OData) -> VirtualLink

Copy o-data words and vertex signs into the library's VirtualLink representation.

# Arguments
`G` must satisfy the real-Gauss input requirements of `Topology.set_data!`,
using [the o-data convention](../specs/heegaard-koda-conversion.md#71-保持するデータ).
Construction and consistency errors from that setter propagate. The caller
must separately establish the closed-normal conditions when a 3-manifold
interpretation is required; this wrapper does not call `is_closed`.

# Returns
A new `VirtualLink` populated from the real Gauss words and orientation vector.

# Guarantees
1. The stored real-Gauss input represents the same signed o-data, within
   the setter's domain and conventions. [proved/agent: by construction]
2. The input is not modified; its words and sign vector are copied before
   being supplied to the result. [proved/agent: by construction]

# Not guaranteed
O-data metadata is not transported. A successful copy is not a closedness,
manifold-identification or planar-drawing certificate.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}());

julia> Topology.real_gauss_code(odata_to_virtual_link(g))
([[1, -1]], [1])
```

See also [`virtual_link_to_odata`](@ref), [`odata_to_spine`](@ref).

Role: convert
"""
function odata_to_virtual_link(G::OData)
    T = parentmodule(@__MODULE__)
    link = T.VirtualLink()
    T.set_data!(link; rgauss=deepcopy(G.gauss), orientation=copy(G.orientation))
    link
end

# cf: lib/topology.jl/src/virtual_link.jl VirtualLink — gauss and orientation are the stored o-data presentation
"""
    virtual_link_to_odata(link) -> OData

Copy a VirtualLink's stored Gauss words and orientation signs into o-data.

# Arguments
`link` must be a `VirtualLink`; other inputs raise `KodaInputError` at stage
`:virtual_link_to_odata` saying that a VirtualLink was expected. The caller
establishes [the closed-normal o-data domain](../specs/heegaard-koda-conversion.md#71-保持するデータ)
before using the result for a 3-manifold construction; this copy does not
check that domain.

# Returns
An `OData` with copied Gauss words and signs, and metadata identifying the
normal-o-graph representation and VirtualLink source.

# Guarantees
1. The returned words and orientation equal the link's stored values.
   [proved/agent: by construction]
2. The input is not modified; its word vectors and sign vector are copied.
   [proved/agent: by construction]

# Not guaranteed
The metadata's representation label is not a closed-normal certificate.
Other VirtualLink data and arbitrary metadata are not retained.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); b = virtual_link_to_odata(odata_to_virtual_link(g)); (b.gauss == g.gauss, b.orientation == g.orientation)
(true, true)
```

See also [`odata_to_spine`](@ref).

Role: convert
"""
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

"""
    spine_to_odata(spine::FlowSpine; regularize::Bool=false) -> OData

Encode a supported flow spine as signed o-data, optionally regularizing its sectors.

# Arguments
`spine` follows [the forward encoding domain](../specs/heegaard-koda-conversion.md#7-branched-spine-と-o-data-の符号化復号).
Without regularization, vertex-free data raise `KodaNotComputed` with reason
`:vertex_free_odata_not_encoded`, and recorded nondisk sectors raise it with
reason `:non_disk_sector`. `regularize=true` requests the supported annulus
regularization; its domain and conversion failures propagate. Encoding
failures from `to_ograph` also propagate. Sector records are supplied data;
their completeness is the caller's responsibility.

# Returns
An `OData` encoding the supported spine, or the regularized result with its
regularization metadata.

# Guarantees
1. Without a required regularization, the result is the encoding returned
   by `to_ograph`. [proved/agent: by construction]
2. The input is not modified. [proved/agent: by construction]

# Not guaranteed
Regularization is not implemented for every singular spine. Identification
with the thickened 3-manifold inherits
[the Koda conversion limits](../specs/heegaard-koda-conversion.md#133-独立な-c2-と閉入力の追加照合);
this wrapper supplies no independent proof.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§6.3, 7.0, 13.3 — sector domain and encoding limitations.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}()); ph = odata_to_heegaard_koda(g).value;

julia> compare_odata_words(g, spine_to_odata(heegaard_to_spine(ph).value))["equal_up_to_metadata"]
true
```

See also [`to_ograph`](@ref).

Role: convert
"""
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

"""
    compare_odata_words(a::OData, b::OData) -> Dict{String,Any}

Compare two directed o-data word systems up to relabeling while checking vertex signs.

# Arguments
`a` and `b` use the [o-data encoding](../specs/heegaard-koda-conversion.md#74-巡回語と復号).
The caller must supply valid passage labels and vertex-sign vectors, and
each stored component word must be nonempty; malformed input can raise
`BoundsError` or lookup errors instead of a verdict. Matching component
counts above eight throw `ArgumentError` stating that the correspondence
search is too large. Metadata are ignored.

# Returns
A dictionary with `gauss_equal` and `orientation_equal` for literal field
equality, `equal_up_to_metadata` for a sign-preserving word correspondence,
`correspondence_count` for distinct passage-preserving label maps (including
those that fail vertex-sign agreement), and a `compared_by` description.

# Guarantees
1. On valid inputs within the search domain, `equal_up_to_metadata` is true
   exactly when a label bijection, component permutation and cyclic shifts
   preserve directed passages and vertex signs. Neither reversal nor a
   diagram move is included. [proved/agent: by construction]
2. Both inputs are unchanged. [proved/agent: by construction]

# Not guaranteed
This is representation comparison, not a test for spine moves or manifold
homeomorphism. `correspondence_count` counts label maps, not all choices of
component permutations and starting positions. Literal sign-vector equality
alone need not be sign agreement under a nontrivial correspondence.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-koda-conversion.md §§7.4, 10.1 — directed cyclic words and equality of their encoded data.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes

julia> g = OData([[1,-1]], [1], Dict{String,Any}());

julia> compare_odata_words(g, g)["equal_up_to_metadata"]
true
```

See also [`spine_to_odata`](@ref), [`odata_to_spine`](@ref).

Role: condition
"""
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
