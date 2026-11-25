"""
    ograph_path_endpoint_permutations(target, reference) -> Vector{Vector{Int}}

Return the cyclic-root and crossing-label identifications of two closed one-circuit presentations.

# Arguments
`target` and `reference` must be coherent `VirtualLink`-compatible presentations.
Cuts throw `ArgumentError` naming the closed-endpoint restriction. Multiple
components or unequal word lengths give an empty result. The caller supplies
valid passage incidence and vertex signs; these are not independently
validated and malformed inputs can give errors or meaningless identifications.

# Returns
Vectors `p` of target edge-occurrence indices in reference order, as specified
by the [endpoint-coordinate specification](../specs/ograph-search-db.md#affine-caches-and-endpoint-coordinates).

# Guarantees
1. For valid nonempty one-circuit inputs, every returned permutation preserves passage and vertex signs under a bijection of crossing labels. [proved/agent: by construction]
2. All matching cyclic roots are returned in increasing starting-index order. [proved/agent: by construction]
3. Assigned weights are not compared and neither input is modified. [proved/agent: by construction]

# Not guaranteed
No move-equivalence or manifold-homeomorphism search is performed. The result
does not give an induced homology marking. Two empty words return no
permutations, so an empty result is not always proof of nonisomorphism.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — explicit endpoint permutations and the homology-marking boundary.

# Examples
```jldoctest
julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4));

julia> code, signs = real_gauss_code(v); iso = VirtualLink(); set_data!(iso; rgauss=[circshift(only(code),-3)], orientation=signs);

julia> [6,7,8,1,2,3,4,5] in ograph_path_endpoint_permutations(iso, v)
true
```

Role: condition
"""
function ograph_path_endpoint_permutations(target,reference)
    any(_is_cut_component,target.gauss) && throw(ArgumentError("closed endpoints only"))
    any(_is_cut_component,reference.gauss) && throw(ArgumentError("closed endpoints only"))
    a,sa=real_gauss_code(target);b,sb=real_gauss_code(reference)
    length(a)==length(b)==1 || return Vector{Int}[]
    a=only(a);b=only(b);length(a)==length(b) || return Vector{Int}[]
    result=Vector{Int}[]
    for start in eachindex(a)
        indices=[mod1(start+i-1,length(a)) for i in eachindex(b)]
        mapping=Dict{Int,Int}();good=true
        for (i,j) in enumerate(indices)
            sign(a[j])==sign(b[i]) || (good=false;break)
            old=abs(a[j]);new=abs(b[i])
            haskey(mapping,old) && mapping[old]!=new && (good=false;break)
            mapping[old]=new
            sa[old]==sb[new] || (good=false;break)
        end
        good && length(unique(values(mapping)))==length(sa) && push!(result,indices)
    end
    result
end

"""
    lift_ograph_search_path!(store, run, path, new_source_state_id; source_permutation=nothing) -> Vector{NamedTuple}

Return the stored operation rows of a geometric move word lifted to an isomorphic unweighted source presentation.

# Arguments
- `store`, `run`: an open migrated search store and a registered run, as for
  [`expand_ograph_search_state!`](@ref).
- `path`: stored contiguous path rows returned by the search queries. The
  caller supplies coherent row data; it is not authenticated against the
  database as a complete path. Checked: moves other than MP or literal 02
  throw `ArgumentError`. An empty path returns immediately.
- `new_source_state_id`: a stored presentation isomorphic to the path source
  in the sense of [`ograph_path_endpoint_permutations`](@ref).
- `source_permutation`: one of those endpoint permutations, or `nothing` to
  choose the first. `ArgumentError` identifies nonisomorphic sources, an invalid
  permutation, a lost vertex marking or a nonunique lifted candidate.
  State loading, candidate enumeration and malformed-row errors propagate.

# Returns
The lifted operation rows in order, with the same field meanings as
[`query_ograph_search_paths`](@ref); an empty path gives an empty vector.

# Guarantees
1. On successful return, each row is a stored successful unweighted move with the original move kind, parameter and direction, transported using a consistent vertex marking. [proved/agent: by construction]
2. The first row starts at `new_source_state_id`, and subsequent rows start at the preceding lifted target. [proved/agent: by construction]
3. The database and run associations may change; the input path is not modified. [proved/agent: by construction]

# Not guaranteed
This is a lift of presentation isomorphisms, not an arbitrary manifold
homeomorphism or a weighted move word. A later failure can leave earlier
expansions persisted. Successful return does not certify framing transport.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — marked geometric path lifting and separate weighted replay.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_path_constraints.jl"); label="affine relation smoke");

julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); sid = put_ograph_search_state!(store, v).state_id;

julia> rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, sid; target_state_id=String(first(rows).target_state_id), kinds=("02",), max_depth=1);

julia> code, signs = real_gauss_code(v); iso = VirtualLink(); set_data!(iso; rgauss=[circshift(only(code),-3)], orientation=signs); iso_id = put_ograph_search_state!(store, iso).state_id;

julia> lifted = lift_ograph_search_path!(store, run, path, iso_id; source_permutation=[6,7,8,1,2,3,4,5]);

julia> only(lifted).source_state_id == iso_id
true

julia> close_store!(store);
```

Role: other
"""
function lift_ograph_search_path!(store,run,path,new_source_state_id;
        source_permutation=nothing)
    steps=collect(path);isempty(steps) && return NamedTuple[]
    all(e->e.move_kind in ("MP","02"),steps) ||
        throw(ArgumentError("path lifting currently supports MP and literal 02"))
    current=load_ograph_search_state(store,String(new_source_state_id))
    reference=load_ograph_search_state(store,first(steps).source_state_id)
    permutations=ograph_path_endpoint_permutations(current,reference)
    isempty(permutations) && throw(ArgumentError("word sources are not isomorphic"))
    p=source_permutation===nothing ? first(permutations) : Int.(source_permutation)
    p in permutations || throw(ArgumentError("source permutation is not an isomorphism"))
    cg,_=real_gauss_code(current);rg,_=real_gauss_code(reference)
    mapping=Dict(abs(only(rg)[i])=>abs(only(cg)[p[i]]) for i in eachindex(p))
    lifted=NamedTuple[]
    for e in steps
        reference=load_ograph_search_state(store,e.source_state_id)
        rg,rs=real_gauss_code(reference);cg,cs=real_gauss_code(current)
        matches=[q for q in ograph_path_endpoint_permutations(current,reference) if
            all(abs(only(cg)[q[i]])==mapping[abs(only(rg)[i])] for i in eachindex(q))]
        isempty(matches) && throw(ArgumentError("lift vertex marking was lost"))
        p=first(matches);data=e.candidate
        vertices=[sign(v)*mapping[abs(v)] for v in Int.(data["vertices"])]
        locations=[[1,p[Int(loc[2])]] for loc in data["locations"]]
        selector=e.move_kind=="MP" ? Dict("vertex_arities"=>[length(vertices)]) :
            Dict("direction"=>e.direction)
        rows=expand_ograph_search_state!(store,run,current,e.move_kind;selector=selector)
        ids=String[String(row.operation_id) for row in rows if String(row.outcome)=="success"]
        options=_ograph_search_path_rows(store,ids)
        matching=[a for a in options if a.candidate["parameter"]==data["parameter"] &&
            Int.(a.candidate["vertices"])==vertices && a.candidate["locations"]==locations &&
            a.direction==e.direction]
        length(matching)==1 || throw(ArgumentError("lift candidate is not unique"))
        a=only(matching)
        rc=_normalized_02(_opath_candidate(reference,data))
        cc=_normalized_02(_opath_candidate(current,a.candidate))
        rm,ra,_=_vertex_maps(rs,rc.vertices,rc.before,rc.after;role_ordered=e.move_kind=="MP")
        cm,ca,_=_vertex_maps(cs,cc.vertices,cc.before,cc.after;role_ordered=e.move_kind=="MP")
        next_map=Dict{Int,Int}()
        for (old,new) in rm
            haskey(cm,mapping[old]) && (next_map[new]=cm[mapping[old]])
        end
        for i in eachindex(ra);next_map[ra[i]]=ca[i];end
        mapping=next_map;push!(lifted,a)
        current=load_ograph_search_state(store,a.target_state_id)
    end
    lifted
end
