"""All explicit cyclic-root / label identifications of two closed one-circuit geometries.

Each permutation p indexes target edge occurrences in reference order. Weights
are not compared, and no manifold homology marking is silently inferred.
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

"""Lift a DB-selected MP/02 word to an isomorphic geometric source.

All original vertex roles are tracked through every move. Intermediate graph
automorphisms are therefore not chosen independently. Every lifted move is
recorded by expand_ograph_search_state!, with weights still excluded.
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
