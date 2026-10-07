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

"""Extract C*x=d and y=T*x+b for an authenticated closed MP/02 edge.

The source code uses only addition/subtraction of edge portions after checking
fixed internal coordinates. Evaluating one admissible origin and each FREE
coordinate basis vector therefore extracts the exact affine map, not an
empirical guess from arbitrary samples. Fixed-coordinate columns are zero;
their effect is included in b. Public weighted replay checks solved witnesses.
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

"""Solve all integer H insertions on a DB-selected path at once.

H variables are present at EVERY stage, including the source and the target.
If target_generator is supplied, solve for a shift s as an additional integer
variable instead of fixing the final weight. shift_lattice=(offset,step)
describes every attainable s for this path and the given source weights.
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
