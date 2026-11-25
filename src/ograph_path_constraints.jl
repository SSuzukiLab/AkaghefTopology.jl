# Exact affine transport on a FIXED geometric MP/literal-02 path. Geometry
# enumeration and state identity remain unweighted. H variables enter only
# after the geometric path has been selected.
const OGRAPH_PATH_CONSTRAINT_VERSION = "ograph-path-constraints/1"

_opath_flat(v) = BigInt[round(BigInt,x) for w in v.weights for x in w]

function _opath_set!(v, values)
    code,_ = real_gauss_code(v)
    cursor=0
    weights=Vector{Vector{BigInt}}()
    for word in code
        push!(weights,BigInt.(values[cursor+1:cursor+length(word)]))
        cursor+=length(word)
    end
    cursor==length(values) || throw(DimensionMismatch("path weight dimensions"))
    set_weight!(v,weights)
    v
end

"H columns in edge-occurrence order, with every incidence added."
function _opath_H(v)
    code,signs=real_gauss_code(v)
    B=zeros(BigInt,sum(length,code),length(signs)); cursor=0
    for word in code
        any(iszero,word) && throw(ArgumentError("closed paths only"))
        for (i,visit) in enumerate(word)
            B[cursor+i,abs(visit)]-=1
            B[cursor+mod1(i+1,length(word)),abs(visit)]+=1
        end
        cursor+=length(word)
    end
    B
end

function _opath_candidate(v,data)
    candidates=movable(v,String(data["type"]);
        param=String(data["parameter"]),vertices=Int.(data["vertices"]),strict=true)
    matching=filter(c -> _ograph_search_match_candidate(c,data),candidates)
    length(matching)==1 || throw(ArgumentError("path candidate is not authentic"))
    only(matching)
end

"""
    ograph_move_weight_relation(store, edge) -> NamedTuple

Return the integer affine weight relation of one stored geometric move.

# Arguments
- `store`: an open writable search store. Existing cached relations must be
  authentic; cache hits are decoded without fresh move authentication.
- `edge`: a successful closed MP or literal-02 path row. Other kinds throw
  `ArgumentError`. On a cache miss, the candidate must match the source and
  target; `ArgumentError` identifies unauthentic candidates, incompatible
  fixed coordinates, open fragments or target geometry mismatch. Database,
  malformed-row and numeric conversion errors propagate.
  See [exact path-weight scope](../specs/ograph-search-db.md#exact-h-insertion-feasibility)
  and [cache/coordinate requirements](../specs/ograph-search-db.md#affine-caches-and-endpoint-coordinates).

# Returns
Big-integer matrices `C`, `T` and vectors `d`, `b`. `C*x == d` describes the
fixed-coordinate admissibility conditions and `T*x + b` gives the destination
weights in edge-occurrence order. Matrix columns refer to source edges;
rows of `T` and entries of `b` refer to target edges.

# Guarantees
1. For authentic rows and cache entries, the result describes the integer affine extension of the move's weight transport under the linked path specification. [proved/agent: by construction]
2. Repeated calls can reuse the stored relation for the operation and algorithm version. [proved/agent: by construction]
3. The database may gain the relation table and a cached row; `edge` is not modified. [proved/agent: by construction]

# Not guaranteed
The relation alone does not verify disk equations, a framing, or a geometric
path's endpoint identification. Arbitrarily large integer solutions need not
be faithfully represented by the `Float64` weights of a `VirtualLink`;
replay through the public weighted API is required before claiming a move word.
Manually altered cached rows are not independently certified here.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — affine integer transport and required weighted replay.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_path_constraints.jl"); label="affine relation smoke");

julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); sid = put_ograph_search_state!(store, v).state_id;

julia> rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> edge = only(query_ograph_search_path(store, sid; target_state_id=String(first(rows).target_state_id), kinds=("02",), max_depth=1));

julia> relation = ograph_move_weight_relation(store, edge); relation == ograph_move_weight_relation(store, edge)
true

julia> close_store!(store);
```

Role: access
"""
function ograph_move_weight_relation(store,edge)
    edge.move_kind in ("MP","02") || throw(ArgumentError("MP and literal 02 only"))
    MoveGraphStore._dbexec(store.db,"""
        CREATE TABLE IF NOT EXISTS ograph_search_edge_weight_relations (
          operation_id TEXT NOT NULL REFERENCES ograph_search_operations(operation_id),
          algorithm_version TEXT NOT NULL, relation_json TEXT NOT NULL,
          PRIMARY KEY(operation_id,algorithm_version))
    """)
    saved=MoveGraphStore._one(store.db,"""
        SELECT relation_json FROM ograph_search_edge_weight_relations
        WHERE operation_id=? AND algorithm_version=?
    """,(String(edge.operation_id),OGRAPH_PATH_CONSTRAINT_VERSION))
    if saved!==nothing
        data=JSON.parse(String(saved.relation_json))
        function matrix(key,m,n)
            a=zeros(BigInt,m,n)
            for i in 1:m; a[i,:]=parse.(BigInt,String.(data[key][i])); end
            a
        end
        n=Int(data["source_edges"]);m=length(data["b"]);r=length(data["d"])
        return (C=matrix("C",r,n),d=parse.(BigInt,String.(data["d"])),
                T=matrix("T",m,n),b=parse.(BigInt,String.(data["b"])))
    end
    source=load_ograph_search_state(store,edge.source_state_id)
    target=load_ograph_search_state(store,edge.target_state_id)
    candidate=_opath_candidate(source,edge.candidate)
    normalized=_normalized_02(candidate)
    code,_=real_gauss_code(source)
    offsets=cumsum([0;length.(code)])
    n=sum(length,code)
    fixed=Dict{Int,BigInt}()
    function fix!(index,value)
        haskey(fixed,index) && fixed[index]!=value &&
            throw(ArgumentError("inconsistent internal source rules"))
        fixed[index]=value
    end
    common=min(length(normalized.before.gauss),length(normalized.after.gauss))
    for (i,fragment) in enumerate(normalized.before.gauss)
        component,start=normalized.locations[i]
        any(iszero,code[component]) && throw(ArgumentError("closed paths only"))
        rules=BigInt.(_move_integer.(normalized.before.weights[i]))
        for k in 1:length(fragment)-1
            fix!(offsets[component]+mod1(start+k,length(code[component])),rules[k+1])
        end
        if i>common
            length(fragment)==length(code[component]) ||
                throw(ArgumentError("source-only fragment is not a component"))
            fix!(offsets[component]+start,first(rules)+last(rules))
        end
    end
    origin=zeros(BigInt,n)
    for (i,value) in fixed; origin[i]=value; end
    function apply(values)
        v=_opath_set!(deepcopy(source),values)
        _apply_candidate!(v,candidate)
        real_gauss_code(v)==real_gauss_code(target) ||
            throw(ArgumentError("affine target geometry mismatch"))
        _opath_flat(v)
    end
    b=apply(origin)
    T=zeros(BigInt,length(b),n)
    for j in 1:n
        haskey(fixed,j) && continue
        basis=copy(origin); basis[j]+=1
        T[:,j]=apply(basis)-b
    end
    indices=sort!(collect(keys(fixed)))
    C=zeros(BigInt,length(indices),n)
    for (i,j) in enumerate(indices); C[i,j]=1; end
    d=BigInt[fixed[j] for j in indices]
    # Check several mixed coordinates against the actual transport. This is
    # a regression check; exactness follows from the arithmetic above.
    mixed=copy(origin)
    for j in 1:n; haskey(fixed,j) || (mixed[j]=mod(j*7,11)-5); end
    apply(mixed)==T*mixed+b || error("non-affine move arithmetic")
    relation=Dict("source_edges"=>n,"C"=>[string.(collect(row)) for row in eachrow(C)],
        "d"=>string.(d),"T"=>[string.(collect(row)) for row in eachrow(T)],"b"=>string.(b))
    MoveGraphStore._dbexec(store.db,"""
        INSERT OR IGNORE INTO ograph_search_edge_weight_relations VALUES (?,?,?)
    """,(String(edge.operation_id),OGRAPH_PATH_CONSTRAINT_VERSION,_ograph_canonical_json(relation)))
    (C=C,d=d,T=T,b=b)
end

"""
    solve_ograph_path_weights(store, weighted_source, path; target_weights=nothing, target_generator=nothing, target_permutation=nothing, relation_cache=Dict()) -> NamedTuple

Solve the integer H-insertion constraints along a fixed geometric o-graph path.

# Arguments
`store` is an open search store; `path` is a nonempty composable sequence of
stored step records starting at the geometry of `weighted_source`.
Empty paths, mismatched source geometry and noncomposable paths throw
`ArgumentError`. Supply a coherent closed `VirtualLink` with finite integral
weights in the [integer framing domain](../specs/ograph-search-db.md#exact-h-insertion-feasibility).
This function does not enforce that domain: source weights are rounded to
integers. `target_weights` defaults to those source weights.
`target_generator`, when present, supplies an integral endpoint shift direction.
Wrong target or generator length throws `DimensionMismatch`.
`target_permutation` must be the caller's valid endpoint coordinate permutation;
bijectivity is not checked. A supplied cache must contain the correct affine
relations for its operation IDs. Conversion, indexing, database and relation
errors propagate; an H coefficient outside `Int` can throw `InexactError`.

# Returns
A named tuple with `solution`, `h_moves`, `shift`, `shift_lattice`,
`target_weights`, `target_generator` and `target_permutation`.
The integer-system result determines feasibility. Successful solutions include
nonzero H insertions; a requested shift additionally has its chosen value and
integer lattice. Without an integer solution the insertions are empty and
both shift fields are `nothing`.

# Guarantees
1. Under the stated domain and relation preconditions, the equations express
   the linked exact H-insertion problem on this path, including insertions
   at its source and target. [proved/agent: by construction]
2. For an integer solution with a shift direction, `shift_lattice` describes
   all shifts admitted by that integer system; step zero means a singleton.
   [proved/agent: by construction]
3. The source and path are unchanged. Relation caching can write to the store
   and to `relation_cache`; returned generator and permutation values can
   alias the supplied values. [proved/agent: by construction]

# Not guaranteed
Feasibility concerns this path, not all move paths or a manifold classification.
The result is not independently replayed through the public weighted moves.
Use `record_weighted_ograph_path_trial!` to check that replay. Incorrect cached
relations, nonintegral source weights, invalid endpoint correspondences and
loss of integer precision in floating-point weight storage are outside the
mathematical guarantee.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — the fixed-path integer problem and its separate replay check.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); set_weight!(v, [1,-1,-1,0,1,0,-1,0]); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_path_constraints.jl"); label="affine relation smoke"); edges = expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, s.state_id; target_state_id=String(first(edges).target_state_id), kinds=("02",), max_depth=1); target = replay_weighted_ograph_path(store, v, path).target; move!(target, "H"; vertices=[1], param=3);

julia> solve_ograph_path_weights(store, v, path; target_weights=BigInt.(vcat(target.weights...))).solution.status
:integer_solution

julia> close_store!(store);
```

See also [`replay_weighted_ograph_path`](@ref), [`ograph_move_weight_relation`](@ref).

Role: convert
"""
function solve_ograph_path_weights(store,weighted_source,path;
        target_weights=nothing,target_generator=nothing,target_permutation=nothing,
        relation_cache=Dict())
    steps=collect(path)
    isempty(steps) && throw(ArgumentError("nonempty geometric path required"))
    ograph_search_state(weighted_source).state_id==first(steps).source_state_id ||
        throw(ArgumentError("source geometry mismatch"))
    geometries=[load_ograph_search_state(store,first(steps).source_state_id);
                [load_ograph_search_state(store,e.target_state_id) for e in steps]]
    blocks=[_opath_H(v) for v in geometries]
    counts=size.(blocks,2)
    offsets=cumsum([0;counts])
    variable_shift=target_generator!==nothing
    nvar=last(offsets)+Int(variable_shift)
    M=zeros(BigInt,length(_opath_flat(weighted_source)),nvar)
    M[:,1:counts[1]]=blocks[1]
    constant=_opath_flat(weighted_source)
    rows=Matrix{BigInt}[]; rhs=Vector{BigInt}[]
    for (i,edge) in enumerate(steps)
        edge.source_state_id==ograph_search_state(geometries[i]).state_id ||
            throw(ArgumentError("path is not composable"))
        relation=get!(relation_cache,edge.operation_id) do
            ograph_move_weight_relation(store,edge)
        end
        push!(rows,relation.C*M)
        push!(rhs,relation.d-relation.C*constant)
        M=relation.T*M
        constant=relation.T*constant+relation.b
        M[:,offsets[i+1]+1:offsets[i+2]]+=blocks[i+1]
    end
    if target_permutation!==nothing
        M=M[Int.(target_permutation),:]
        constant=constant[Int.(target_permutation)]
    end
    destination=target_weights===nothing ? _opath_flat(weighted_source) : BigInt.(target_weights)
    length(destination)==length(constant) || throw(DimensionMismatch("target weight length"))
    if variable_shift
        length(target_generator)==length(constant) || throw(DimensionMismatch("generator length"))
        M[:,end]-=BigInt.(target_generator)
    end
    push!(rows,M); push!(rhs,destination-constant)
    A=reduce(vcat,rows); b=reduce(vcat,rhs)
    solution=solve_integer_weight_system(A,b)
    h_moves=NamedTuple[]
    shift=nothing; lattice=nothing
    if solution.status==:integer_solution
        for stage in eachindex(blocks),vertex in 1:counts[stage]
            amount=solution.particular[offsets[stage]+vertex]
            iszero(amount) || push!(h_moves,(after_step=stage-1,vertex=vertex,amount=Int(amount)))
        end
        if variable_shift
            shift=solution.particular[end]
            period=foldl(gcd,abs.(solution.kernel[end,:]);init=big(0))
            lattice=(offset=shift,step=period)
        end
    end
    (solution=solution,h_moves=h_moves,shift=shift,shift_lattice=lattice,
     target_weights=destination,target_generator=target_generator,
     target_permutation=target_permutation)
end

"""
    record_ograph_path_constraint_trial!(store, run_id, source, path, result) -> String

Store an exact H-insertion feasibility result for a geometric move path and return its trial identifier.

# Arguments
`store` must be an open search store and `run_id` an existing run. `source`
is the weighted `VirtualLink` and `path` the ordered step records used to
compute `result` with [`solve_ograph_path_weights`](@ref). Their mutual
consistency is the caller's responsibility: this function does not solve or
replay them. Missing result fields, nonintegral source weights, serialization
errors and database constraints propagate errors.

# Returns
A `String` trial identifier, derived from the run, ordered operation IDs,
source/target settings and current constraint-algorithm version.

# Guarantees
1. Successful return stores the result's integer system, status and solution
   certificate or obstruction, together with its path and endpoint settings.
   Source, path and result values are unchanged. [proved/agent: by construction]
2. Repeating the same identity replaces the trial payload. The stored claim
   concerns [exact H-insertion feasibility](../specs/ograph-search-db.md#exact-h-insertion-feasibility)
   on that geometric path, with the scope of the supplied solver result.
   [proved/agent: by construction]

# Not guaranteed
The supplied result is not checked against the equations, path or source.
A stored positive solution is not a replayed weighted move word; use
`record_weighted_ograph_path_trial!` for that independent check. A negative
result on one path does not exclude other geometric paths.

# Sources
cf: lib/topology.jl/docs/specs/ograph-search-db.md §Exact H-insertion feasibility — persisted equations, witnesses and obstructions, followed by weighted replay.

# Examples
```jldoctest
julia> Base.include(@__MODULE__, joinpath(pkgdir(Topology), "examples", "VLExample.jl"));

julia> store = open_store(); migrate_ograph_search_store!(store); v = vl_example("s2xs1"); set_weight!(v, [1,-1,-1,0,1,0,-1,0]); s = put_ograph_search_state!(store, v);

julia> r = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_path_constraints.jl"); label="affine relation smoke"); edges = expand_ograph_search_state!(store, r, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, s.state_id; target_state_id=String(first(edges).target_state_id), kinds=("02",), max_depth=1); target = replay_weighted_ograph_path(store, v, path).target;

julia> move!(target, "H"; vertices=[1], param=3); weights = BigInt[e.weight for e in real_edge_table(target)]; result = solve_ograph_path_weights(store, v, path; target_weights=weights);

julia> record_ograph_path_constraint_trial!(store, r, v, path, result);

julia> length(query_ograph_path_constraint_trials(store; source_state_id=s.state_id, target_weights=weights, status=:integer_solution))
1

julia> close_store!(store);
```

See also [`solve_ograph_path_weights`](@ref), [`query_ograph_path_constraint_trials`](@ref).

Role: other
"""
function record_ograph_path_constraint_trial!(store,run_id,source,path,result)
    MoveGraphStore._dbexec(store.db,"""
        CREATE TABLE IF NOT EXISTS ograph_search_constraint_trials (
          trial_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES ograph_search_runs(run_id),
          path_fingerprint TEXT NOT NULL, operation_ids_json TEXT NOT NULL,
          algorithm_version TEXT NOT NULL, result_json TEXT NOT NULL)
    """)
    ids=String[e.operation_id for e in path]
    serialized=_ograph_canonical_json(ids)
    fingerprint=_ograph_sha256(serialized)
    sol=result.solution
    payload=Dict("status"=>string(sol.status),"source_weights"=>string.(_opath_flat(source)),
        "A"=>[string.(collect(row)) for row in eachrow(sol.A)],"b"=>string.(sol.b),
        "h_moves"=>result.h_moves,"shift_lattice"=>result.shift_lattice,
        "target_weights"=>string.(result.target_weights),
        "target_generator"=>result.target_generator===nothing ? nothing : string.(result.target_generator),
        "target_permutation"=>result.target_permutation)
    if sol.status==:integer_solution
        payload["particular"]=string.(sol.particular)
        payload["kernel"]=[string.(collect(row)) for row in eachrow(sol.kernel)]
        payload["smith_certificate"]=Dict("rank"=>sol.certificate.rank,
            "S"=>[string.(collect(row)) for row in eachrow(sol.certificate.S)],
            "U"=>[string.(collect(row)) for row in eachrow(sol.certificate.U)],
            "V"=>[string.(collect(row)) for row in eachrow(sol.certificate.V)])
    else
        obs=sol.obstruction
        payload["obstruction"]=Dict("kind"=>string(obs.kind),
            "row_combination"=>string.(obs.row_combination),
            "lhs_coefficients"=>string.(obs.lhs_coefficients),
            "divisor"=>string(obs.divisor),"rhs_value"=>string(obs.rhs_value))
    end
    settings=Dict("source_weights"=>payload["source_weights"],
        "target_weights"=>payload["target_weights"],
        "target_generator"=>payload["target_generator"],
        "target_permutation"=>payload["target_permutation"],
        "version"=>OGRAPH_PATH_CONSTRAINT_VERSION)
    id="ograph-constraint-trial-"*_ograph_sha256(String(run_id)*fingerprint*
        _ograph_canonical_json(settings))
    MoveGraphStore._dbexec(store.db,"""
        INSERT OR REPLACE INTO ograph_search_constraint_trials
        VALUES (?,?,?,?,?,?)
    """,(id,String(run_id),fingerprint,serialized,OGRAPH_PATH_CONSTRAINT_VERSION,
           _ograph_canonical_json(payload)))
    id
end

"""
    query_ograph_path_constraint_trials(store; source_state_id=nothing, target_weights=nothing, status=nothing, limit::Integer=100) -> Vector{NamedTuple}

Read stored integer path-weight solutions or obstructions and their separate replay outcomes.

# Arguments
- `store`: an open store containing constraint-trial and weighted-trial tables;
  absent tables or malformed stored data propagate database/decoding errors.
- `source_state_id`: optional exact source filter.
- `target_weights`: optional fixed integer target; conversion to `BigInt` must
  succeed. Variable-shift trials are excluded when this filter is supplied.
- `status`: `nothing`, `:integer_solution`, `:rational_only` or `:inconsistent`
  (the corresponding strings are also accepted). Unknown status throws
  `ArgumentError`. This is an equation-system status.
- `limit`: positive and representable as `Int`; nonpositive limits throw
  `ArgumentError`. Stored source weights and H data must also convert to `Int`.

# Returns
Rows ordered by trial ID, with `trial_id`, `run_id`, decoded `path`, `h_moves`,
stored solver `result`, and optional `replay_outcome` for the same run, path,
source weights and H recipe.

# Guarantees
1. Returned records satisfy the supplied filters and their count is at most `limit`. [proved/agent: by construction]
2. Results retain the fixed-path integer-equation meaning of the linked specification; composed walks are not excluded merely for repeated states. [proved/agent: by construction]
3. The database and arguments are not modified. [proved/agent: by construction]

# Not guaranteed
Reading a stored solution does not recompute its certificate or replay its
weighted move word. A missing replay outcome is not a failed replay, and a
stored equation obstruction does not exclude other geometric paths.

# Sources
cf:docs/specs/ograph-search-db.md#exact-h-insertion-feasibility — trial retrieval, composed paths and independent weighted replay.

# Examples
```jldoctest
julia> store = open_store(); migrate_ograph_search_store!(store);

julia> v = VirtualLink(); set_data!(v; rgauss=[[-1,-2,3,4,-4,-3,2,1]], orientation=ones(Int,4)); set_weight!(v, [1,-1,-1,0,1,0,-1,0]); sid = put_ograph_search_state!(store, v).state_id;

julia> run = begin_ograph_search_run!(store, joinpath(pkgdir(Topology), "test", "ograph_path_constraints.jl"); label="affine relation smoke"); rows = expand_ograph_search_state!(store, run, v, "02"; selector=Dict("direction" => "forward"));

julia> path = query_ograph_search_path(store, sid; target_state_id=String(first(rows).target_state_id), kinds=("02",), max_depth=1); target = replay_weighted_ograph_path(store, v, path).target;

julia> move!(target, "H"; vertices=[1], param=3); weights = BigInt[e.weight for e in real_edge_table(target)]; solution = solve_ograph_path_weights(store, v, path; target_weights=weights);

julia> record_ograph_path_constraint_trial!(store, run, v, path, solution); record_weighted_ograph_path_trial!(store, run, v, path; h_moves=solution.h_moves);

julia> saved = query_ograph_path_constraint_trials(store; source_state_id=sid, target_weights=weights, status=:integer_solution); (length(saved), only(saved).replay_outcome)
(1, "admissible")

julia> close_store!(store);
```

Role: access
"""
function query_ograph_path_constraint_trials(store;source_state_id=nothing,
        target_weights=nothing,status=nothing,limit::Integer=100)
    limit>0 || throw(ArgumentError("trial query limit must be positive"))
    clauses=String[];parameters=Any[]
    if source_state_id!==nothing
        push!(clauses,"first_operation.source_state_id=?");push!(parameters,String(source_state_id))
    end
    if status!==nothing
        string(status) in ("integer_solution","rational_only","inconsistent") ||
            throw(ArgumentError("unknown constraint status"))
        push!(clauses,raw"json_extract(trial.result_json,'$.status')=?")
        push!(parameters,string(status))
    end
    if target_weights!==nothing
        push!(clauses,raw"json_extract(trial.result_json,'$.target_generator') IS NULL")
        push!(clauses,raw"json_extract(trial.result_json,'$.target_weights')=?")
        push!(parameters,_ograph_canonical_json(string.(BigInt.(target_weights))))
    end
    sql=raw"""
        SELECT trial.* FROM ograph_search_constraint_trials AS trial
        JOIN ograph_search_operations AS first_operation
          ON first_operation.operation_id=json_extract(trial.operation_ids_json,'$[0]')
    """
    isempty(clauses) || (sql*=" WHERE "*join(clauses," AND "))
    sql*=" ORDER BY trial.trial_id LIMIT ?";push!(parameters,Int(limit))
    result=NamedTuple[]
    for row in MoveGraphStore._rows(store.db,sql,Tuple(parameters))
        data=JSON.parse(String(row.result_json))
        path=_ograph_search_path_rows(store,JSON.parse(String(row.operation_ids_json)))
        geometry=load_ograph_search_state(store,first(path).source_state_id)
        values=parse.(BigInt,String.(data["source_weights"]))
        code,_=real_gauss_code(geometry);cursor=0;weights=Vector{Vector{Int}}()
        for word in code
            push!(weights,Int.(values[cursor+1:cursor+length(word)]));cursor+=length(word)
        end
        replay=MoveGraphStore._one(store.db,"""
            SELECT outcome FROM ograph_search_weighted_trials
            WHERE run_id=? AND path_fingerprint=? AND source_weights_json=? AND h_moves_json=?
        """,(String(row.run_id),String(row.path_fingerprint),_ograph_canonical_json(weights),
              _ograph_canonical_json(data["h_moves"])))
        h_moves=[(after_step=Int(h["after_step"]),vertex=Int(h["vertex"]),amount=Int(h["amount"]))
                 for h in data["h_moves"]]
        push!(result,(trial_id=String(row.trial_id),run_id=String(row.run_id),path=path,
            h_moves=h_moves,result=data,replay_outcome=replay===nothing ? nothing : String(replay.outcome)))
    end
    result
end
