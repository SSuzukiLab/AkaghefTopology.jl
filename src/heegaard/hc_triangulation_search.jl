# cf: docs/specs/heegaard-snappy.md — fixed bounded inverse search and evidence boundaries.

"""
    hc_pachner_neighbors(input; kwargs...) -> Dict{String,Any}

Return backend-certified positive Pachner neighbors of a finite triangulation.

# Arguments
`input` is an `HCClosedTriangulation` or explicit face-pairing dictionary in
[the bridge domain](../specs/heegaard-snappy.md). Keywords include an existing
Python `runner::Cmd` and positive `timeout`; invalid runner or timeout raises
`ArgumentError`. Backend availability, execution and malformed replies raise
`SnapPyUnavailable` or `SnapPyExecutionError`.

# Returns
The bridge reply dictionary, including source `isosig` and neighbor `states`
with triangulations, topology data and inverse-replay move certificates.

# Guarantees
1. A successful return is a successful bridge reply for the unchanged
   input's native positive 2→3 moves. [proved/agent: by construction]
2. Registered backend fixtures have one additional tetrahedron per neighbor,
   unchanged vertex count and inverse replay restoring the source isoSig.
   [verified-computationally: test/heegaard/hc_triangulation.jl, native Pachner-neighbor checks]
3. The input is not modified; temporary files and a child process are used.
   [proved/agent: by construction]

# Not guaranteed
Returned identities are unoriented. This entry does not enumerate all
Pachner move kinds, construct a branching or produce a resolved HC. An
unavailable backend is not evidence that no neighbor exists.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — native Pachner moves and inverse-replay boundary.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.hc_pachner_neighbors(Dict{String,Any}(); timeout=0) catch e; e isa ArgumentError end
true
```

Role: move
"""
hc_pachner_neighbors(input;kwargs...) = _hc_snappy_call(
    input isa HCClosedTriangulation ? hc_triangulation_data(input) : input;
    extra_args=["--pachner-neighbors"],kwargs...)

"""
    hc_prepare_inverse(input; max_depth::Int=2, max_states::Int=64, max_extra_tetrahedra::Int=2, max_nodes::Int=100_000, max_solutions::Int=1024, runner=nothing, timeout::Real=120) -> Dict{String,Any}

Search a bounded set of triangulations and branchings for a checked Koda reconstruction.

# Arguments
`input` is an `HCClosedTriangulation` or explicit face-pairing dictionary in
[the inverse-search domain](../specs/heegaard-snappy.md).
Nonnegative depth and extra-tetrahedron limits and positive state count are
checked with `ArgumentError`. Branching budgets, runner and timeout have the
domains of `hc_find_branching` and `hc_pachner_neighbors`. Initial malformed
triangulation and backend errors can propagate.

# Returns
A dictionary with `status`, `stop`, `hc`, `limits`, `attempts`, `moves` and
current `triangulation`; a successful Koda stage adds flow, o-data, branching,
joining-region and surface-cut certificates.

# Guarantees
1. A `koda_only` result contains the checked flow reconstruction and surface
   cut of its reported triangulation. `hc` remains `nothing`.
   [proved/agent: by construction]
2. Multiple finite vertices stop with `multiple_finite_vertices_unreduced`.
   Search stops distinguish resource bounds and certificate failure.
   [proved/agent: by construction]
3. The input is not modified; backend calls can launch child processes and
   use temporary files. [proved/agent: by construction]

# Not guaranteed
No resolved LC-bound HC or general PL identification is produced. A bounded
search stop is not nonexistence of a suitable diagram. The timeout applies
to backend calls, not the entire search, and returned diagnostic data may
share the original input dictionary.

# Sources
cf: lib/topology.jl/docs/specs/heegaard-snappy.md — bounded inverse search and unresolved reference overlay.

# Examples
```jldoctest
julia> using Topology.HeegaardCodes; H = Topology.HeegaardCodes;

julia> try H.hc_prepare_inverse(Dict{String,Any}(); max_depth=-1) catch e; e isa ArgumentError end
true
```

Role: convert
"""
function hc_prepare_inverse(input;max_depth::Int=2,max_states::Int=64,max_extra_tetrahedra::Int=2,
                            max_nodes::Int=100_000,max_solutions::Int=1024,
                            runner=nothing,timeout::Real=120)
    max_depth>=0 && max_states>0 && max_extra_tetrahedra>=0 || throw(ArgumentError("invalid move limits"))
    original=input isa HCClosedTriangulation ? hc_triangulation_data(input) : input
    t=_hc_explicit_triangulation(original)
    limits=Dict("max_depth"=>max_depth,"max_states"=>max_states,"max_extra_tetrahedra"=>max_extra_tetrahedra,
        "branching_nodes"=>max_nodes,"branching_solutions"=>max_solutions,"backend_timeout_seconds"=>timeout)
    result=Dict{String,Any}("status"=>"stopped","stop"=>nothing,"hc"=>nothing,"limits"=>limits,
        "attempts"=>Any[],"moves"=>Any[],"triangulation"=>original)
    if t.topology.vertices!=1
        result["stop"]="multiple_finite_vertices_unreduced"
        result["detail"]="safe native 2->3 preserves finite vertices; a certified 4->1 reducer is not available"
        return result
    end
    frontier=Any[(data=original,depth=0,moves=Any[])]
    seen=Set{String}(); processed=0; truncated=false; any_branch_budget=false
    while !isempty(frontier) && processed<max_states
        state=popfirst!(frontier); processed+=1
        accepted=Ref{Any}(nothing); failures=Any[]; enumeration_stop=nothing
        try
            hc_find_branching(state.data;max_nodes=max_nodes,max_solutions=max_solutions,accept=orders->begin
                try
                    accepted[]=hc_flow_certificate(state.data;orders=orders,max_nodes=max_nodes)
                    true
                catch e
                    if e isa HCTriangulationError && e.reason in
                        (:not_normal_ograph,:not_one_suture,:bw_not_two_disks,:suture_not_simple)
                        push!(failures,Dict("reason"=>string(e.reason),"detail"=>e.detail));false
                    else
                        rethrow()
                    end
                end
            end)
        catch e
            if e isa HCTriangulationError && e.reason in
                (:not_branchable_on_input,:branching_filter_exhausted,:branching_search_budget,:branching_solution_budget)
                enumeration_stop=string(e.reason)
                any_branch_budget |= e.reason in (:branching_search_budget,:branching_solution_budget)
            else
                result["stop"]="inverse_certificate_failure"
                result["detail"]=sprint(showerror,e)
                result["triangulation"]=state.data;result["moves"]=state.moves
                return result
            end
        end
        push!(result["attempts"],Dict("state"=>processed,"depth"=>state.depth,
            "tetrahedra"=>state.data["tetrahedra"],"branching_stop"=>enumeration_stop,"rejected_branchings"=>failures))
        if accepted[]!==nothing
            f=accepted[]
            result["status"]="koda_only";result["stop"]="polygon_normalization_and_reference_overlay_required"
            result["triangulation"]=state.data;result["moves"]=state.moves
            result["flow_certificate"]=f.certificate;result["koda_boundary_code"]=f.raw
            result["odata"]=Dict("gauss"=>f.encoded.value.gauss,"orientation"=>f.encoded.value.orientation)
            result["branching"]=f.encoded.branching
            result["correspondence"]=f.encoded.correspondence
            result["joining_region"]=f.region
            result["surface_cut"]=hc_surface_cut_system(f.diagram)
            result["processed_states"]=processed
            return result
        end
        state.depth>=max_depth && continue
        Int(state.data["tetrahedra"])>=t.n+max_extra_tetrahedra && continue
        neighbors=hc_pachner_neighbors(state.data;runner=runner,timeout=timeout)
        push!(seen,neighbors["isosig"])
        for child in neighbors["states"]
            child["isosig"] in seen && continue
            if length(frontier)+processed>=max_states
                truncated=true;break
            end
            push!(seen,child["isosig"])
            push!(frontier,(data=child["triangulation"],depth=state.depth+1,moves=vcat(state.moves,Any[child["move"]])))
        end
    end
    result["processed_states"]=processed
    result["stop"]=truncated ? "pachner_state_budget" : any_branch_budget ? "branching_budget" : "no_flow_branching_within_bounds"
    result
end

export hc_prepare_inverse, hc_pachner_neighbors
